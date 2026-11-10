/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Adrian Kummerlaender
 *
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#ifndef COMMUNICATION_OVERLAP_BUFFER_LAYOUT_H
#define COMMUNICATION_OVERLAP_BUFFER_LAYOUT_H

#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <typeindex>
#include <utility>
#include <vector>

namespace olb {

namespace overlap {

/// Per-field, per-cell component masks negotiated for one halo exchange
/**
 * Ordered list of (field, per-cell component bitmasks). Each entry's mask
 * vector is index-aligned with the negotiated cell list; a set bit selects a
 * streaming direction (propagatable fields) resp. a field component. Order is
 * the field-request order and is load-bearing for the wire format. Empty means
 * a dense exchange (no field masked).
 **/
using FieldMasks = std::vector<std::pair<std::type_index, std::vector<std::uint64_t>>>;

/// Partition of a packed multi-field overlap buffer around one typed field
/**
 * Wire contract (shared by every platform's gather/scatter kernels and the
 * CPU-side MultiConcreteCommunicatable): fields are packed sequentially in
 * request order, each field storing its components component-major over the
 * communicated cells.
 *
 * This partition lets a task route one compile-time-known field (POPULATION)
 * through the typed gather/scatter kernels — which read fresh, launch-time
 * device pointers and thus need no per-step registry refresh — while all
 * remaining fields stay on the type-erased any-fields path. The byte layout
 * is unchanged: the typed segment sits at exactly the offset the any-fields
 * kernel would have used.
 **/
struct TypedFieldPartition {
  /// Position of the typed field in the request order, -1 if absent
  int index = -1;
  /// Byte offset of the typed field's segment
  std::size_t begin = 0;
  /// Byte offset just past the typed field's segment
  std::size_t end = 0;
  /// Fields preceding resp. following the typed segment, in request order
  std::vector<std::type_index> before;
  std::vector<std::type_index> after;
};

/// Compute the partition of `fields` around `typedField`
/**
 * sizer(field) must return the byte size of the field's segment for the
 * communicated cell set (i.e. Communicatable::size(cells)).
 *
 * The packed buffer gives no alignment guarantees between segments (the
 * CPU-side MultiConcreteCommunicatable packs tightly, the any-fields
 * kernels copy byte-wise) while the typed kernels issue element-typed
 * loads/stores. If the typed field's segment does not start at a multiple
 * of `alignment` (its element size), the partition folds back to the pure
 * type-erased layout: index = -1 with all fields in `before`.
 **/
template <typename SIZER>
TypedFieldPartition partitionAround(std::type_index typedField,
                                    const std::vector<std::type_index>& fields,
                                    std::size_t alignment,
                                    SIZER&& sizer)
{
  TypedFieldPartition partition;
  std::size_t offset = 0;
  for (unsigned iField=0; iField < fields.size(); ++iField) {
    const std::size_t size = sizer(fields[iField]);
    if (fields[iField] == typedField && offset % alignment == 0) {
      partition.index = static_cast<int>(iField);
      partition.begin = offset;
      partition.end   = offset + size;
    } else if (partition.index == -1) {
      partition.before.emplace_back(fields[iField]);
    } else {
      partition.after.emplace_back(fields[iField]);
    }
    offset += size;
  }
  return partition;
}

/// One field's segment in a packed multi-field overlap buffer
/**
 * Generalises TypedFieldPartition to N interleaved typed segments: every
 * field of a multi-field exchange becomes one segment in request order,
 * classified by how the GPU tasks route it. The byte layout is identical
 * to the CPU-side MultiConcreteCommunicatable / MaskedMultiCommunicatable
 * packing, so heterogeneous CPU/GPU pairs interoperate by construction.
 **/
struct FieldSegment {
  enum class Kind {
    MaskedPropagatable, ///< masked group layout, field-generic masked kernel
    TypedDense,         ///< dense typed kernel (fresh launch-time pointers)
    AnyField,           ///< type-erased any-fields kernel (per-step refresh)
  };
  Kind kind;
  std::type_index field;
  /// Byte offset of the field's segment
  std::size_t begin;
};

/// Ordered typed segments of a packed multi-field overlap buffer
struct MultiFieldPartition {
  std::vector<FieldSegment> segments;
  /// Total byte size of the packed buffer
  std::size_t size = 0;
};

/// Partition `fields` into request-order typed segments
/**
 * classify(field) selects the kind; sizer(field) returns the field's byte
 * size for the communicated cell set. A typed segment whose byte offset is
 * not a multiple of `alignment` (its element size) degrades to AnyField, so
 * the element-typed masked/dense kernels never issue misaligned accesses.
 **/
template <typename CLASSIFIER, typename SIZER>
MultiFieldPartition partitionMultiField(const std::vector<std::type_index>& fields,
                                        std::size_t alignment,
                                        CLASSIFIER&& classify,
                                        SIZER&& sizer)
{
  MultiFieldPartition partition;
  std::size_t offset = 0;
  for (std::type_index field : fields) {
    FieldSegment::Kind kind = classify(field);
    if (kind != FieldSegment::Kind::AnyField && offset % alignment != 0) {
      kind = FieldSegment::Kind::AnyField;
    }
    partition.segments.push_back(FieldSegment{kind, field, offset});
    offset += sizer(field);
  }
  partition.size = offset;
  return partition;
}

/// Buffer layout of a component-masked field exchange
/**
 * Single source of truth for the byte layout of a selector-restricted
 * (e.g. directional streaming) field segment, shared by the CPU-side
 * serialization and the GPU gather/scatter kernels. Both ends of an
 * exchange compute it from the identical negotiated per-cell masks, so
 * the layout agrees by construction — including across heterogeneous
 * CPU / GPU pairs.
 *
 * Cells are bucketed by distinct mask value into groups ordered by
 * ascending mask; cells keep their negotiated relative order within a
 * group. Each group stores its components component-major over its
 * cells. Cells with empty mask (e.g. corner cells of the directional
 * streaming exchange) are dropped entirely. With all-dense masks this
 * reproduces the dense component-major segment byte for byte.
 **/
struct MaskGroupLayout {
  struct Group {
    /// Common component mask of the group's cells
    std::uint64_t mask;
    /// Component indices decoded from mask, ascending
    std::vector<std::uint8_t> components;
    /// Indices into the negotiated cell (and mask) arrays
    std::vector<std::uint32_t> cellIndices;
    /// Byte offset of the group's segment
    std::size_t begin;
  };

  std::vector<Group> groups;
  /// Total byte size of the masked segment
  std::size_t size;
  /// Field components carried per mask bit (1 scalar, DESCRIPTOR::n vector)
  unsigned componentsPerDirection = 1;

  /// Compute the layout from negotiated per-cell masks
  /**
   * Each mask bit is a streaming direction; componentsPerDirection is the number of
   * field components it carries (1 for a scalar lattice, DESCRIPTOR::n for a
   * vector lattice). A set bit k decodes to the direction's contiguous
   * component block [k*n, k*n+n), matching descriptors::q(k,·). With
   * componentsPerDirection == 1 this reproduces the scalar per-population layout.
   **/
  static MaskGroupLayout compute(const std::vector<std::uint64_t>& masks,
                                 std::size_t elementSize,
                                 unsigned componentsPerDirection = 1)
  {
    MaskGroupLayout layout;
    layout.size = 0;
    layout.componentsPerDirection = componentsPerDirection;
    std::vector<std::uint64_t> distinct(masks);
    std::sort(distinct.begin(), distinct.end());
    distinct.erase(std::unique(distinct.begin(), distinct.end()), distinct.end());
    for (std::uint64_t mask : distinct) {
      if (mask == 0) {
        continue; // fully-masked-out cells are not transferred
      }
      Group group;
      group.mask = mask;
      for (std::uint8_t i=0; i < 64; ++i) {
        if (mask & (std::uint64_t{1} << i)) {
          for (unsigned k=0; k < componentsPerDirection; ++k) {
            group.components.emplace_back(i*componentsPerDirection + k);
          }
        }
      }
      for (std::uint32_t iCell=0; iCell < masks.size(); ++iCell) {
        if (masks[iCell] == mask) {
          group.cellIndices.emplace_back(iCell);
        }
      }
      group.begin = layout.size;
      layout.size += group.components.size() * group.cellIndices.size() * elementSize;
      layout.groups.emplace_back(std::move(group));
    }
    return layout;
  }
};

}

}

#endif
