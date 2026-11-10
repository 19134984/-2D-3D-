/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef COMMUNICATABLE_H
#define COMMUNICATABLE_H

#include <vector>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <typeindex>
#include <algorithm>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string_view>

#include "core/blockStructure.h"
#include "communication/overlapBufferLayout.h"
// NOTE: communication/directionalOverlap.h is included further down (after the
// ConcreteCommunicatable primary template and ConstSpan are defined) because it
// pulls descriptor.h, whose fields.h -> platform column.h chain re-enters this
// header and instantiates the GPU column.h ConcreteCommunicatable<Column>
// specializations, which need those definitions to already be visible.

namespace olb {

template <typename T>
class ConstSpan {
private:
  const T* const _base;
  const std::size_t _size;

public:
  ConstSpan(const T* base, std::size_t size):
    _base(base), _size(size) { }

  ConstSpan(const std::vector<T>& v):
    ConstSpan{v.data(), v.size()} { }

  const T* data() const {
    return _base;
  }

  T operator[](std::size_t i) const {
    return _base[i];
  }

  const T* begin() const {
    return _base;
  }
  const T* end() const {
    return _base + _size;
  }

  std::size_t size() const {
    return _size;
  }

};

}

#include "communication/canary.h"

namespace olb {

struct Communicatable {
  virtual ~Communicatable() { };

  virtual std::size_t size(ConstSpan<CellID> indices) const = 0;
  /// Serialize data at locations `indices` to `buffer`
  /**
   * Used for preparing inter-block communication messages
   * \returns size of serialized data in bytes
   **/
  virtual std::size_t serialize(ConstSpan<CellID> indices,
                                std::uint8_t* buffer) const = 0;
  /// Deserialize data at locations `indices` to `buffer`
  /**
   * Used for applying inter-block communication messages
   * \returns size of serialized data in bytes
   **/
  virtual std::size_t deserialize(ConstSpan<CellID> indices,
                                  const std::uint8_t* buffer) = 0;

  /// Ensure per-cell storage for the given cells exists before an exchange
  /**
   * No-op for dense fields. Sparse fields must have their inbound halo tiles
   * materialised host-side before a device scatter, which never allocates.
   **/
  virtual void ensureAllocated(ConstSpan<CellID> indices) { }

  /// Serialize a component-masked segment (see overlap::MaskGroupLayout)
  /**
   * Cross-platform access dispatches virtually so the implementation is
   * instantiated where the concrete columns are complete; supported only
   * by multi-component field communicatables.
   **/
  virtual std::size_t serializeMasked(ConstSpan<CellID> indices,
                                      const overlap::MaskGroupLayout& layout,
                                      std::uint8_t* buffer) const {
    throw std::runtime_error("Masked serialization is not supported by this Communicatable");
  }
  /// Deserialize a component-masked segment (see overlap::MaskGroupLayout)
  /**
   * masks must be index-aligned with indices; they are consumed by the
   * directional canary (see overlap::canaryEnabled).
   **/
  virtual std::size_t deserializeMasked(ConstSpan<CellID> indices,
                                        const overlap::MaskGroupLayout& layout,
                                        ConstSpan<std::uint64_t> masks,
                                        const std::uint8_t* buffer) {
    throw std::runtime_error("Masked deserialization is not supported by this Communicatable");
  }
};

namespace overlap {

/// Serialize a component-masked field segment (CPU reference implementation)
/**
 * fieldArray[iD][cellID] must yield component iD of the cell (e.g. a
 * block's FieldArrayD of POPULATION, whose cyclic columns resolve the
 * propagation shift inside operator[]).
 **/
template <typename FIELD_ARRAY>
std::size_t serializeMasked(FIELD_ARRAY& fieldArray,
                            ConstSpan<CellID> cells,
                            const MaskGroupLayout& layout,
                            std::uint8_t* buffer)
{
  using value_t = std::remove_reference_t<decltype(fieldArray[0][0])>;
  for (const auto& group : layout.groups) {
    std::uint8_t* target = buffer + group.begin;
    for (std::uint8_t component : group.components) {
      const auto& column = fieldArray[component];
      for (std::uint32_t iCell : group.cellIndices) {
        std::memcpy(target, reinterpret_cast<const void*>(&column[cells[iCell]]), sizeof(value_t));
        target += sizeof(value_t);
      }
    }
  }
  return layout.size;
}

/// Deserialize a component-masked field segment (CPU reference implementation)
template <typename FIELD_ARRAY>
std::size_t deserializeMasked(FIELD_ARRAY& fieldArray,
                              ConstSpan<CellID> cells,
                              const MaskGroupLayout& layout,
                              const std::uint8_t* buffer)
{
  using value_t = std::remove_reference_t<decltype(fieldArray[0][0])>;
  for (const auto& group : layout.groups) {
    const std::uint8_t* source = buffer + group.begin;
    for (std::uint8_t component : group.components) {
      auto& column = fieldArray[component];
      for (std::uint32_t iCell : group.cellIndices) {
        std::memcpy(reinterpret_cast<void*>(&column[cells[iCell]]), source, sizeof(value_t));
        source += sizeof(value_t);
      }
    }
  }
  return layout.size;
}

}

template <typename COMMUNICATEE>
class ConcreteCommunicatable final : public Communicatable {
private:
  COMMUNICATEE& _communicatee;

public:
  ConcreteCommunicatable(COMMUNICATEE& communicatee):
    _communicatee{communicatee} { }

  /// Get serialized size for data at locations `indices`
  std::size_t size(ConstSpan<CellID> indices) const override
  {
    return indices.size() * sizeof(typename COMMUNICATEE::value_type);
  }

  /// Serialize data at locations `indices` to `buffer`
  std::size_t serialize(ConstSpan<CellID> indices,
                        std::uint8_t* buffer) const override
  {
    if (meta::is_aligned<typename COMMUNICATEE::value_type>(buffer)) {
      auto* target = reinterpret_cast<typename COMMUNICATEE::value_type*>(buffer);
      for (CellID index : indices) {
        *(target++) = _communicatee[index];
      }
    } else {
      std::uint8_t* target = buffer;
      for (CellID index : indices) {
        std::memcpy(target,
                    reinterpret_cast<const void*>(&_communicatee[index]),
                    sizeof(typename COMMUNICATEE::value_type));
        target += sizeof(typename COMMUNICATEE::value_type);
      }
    }
    return indices.size() * sizeof(typename COMMUNICATEE::value_type);
  }

  /// Deserialize data at locations `indices` to `buffer`
  std::size_t deserialize(ConstSpan<CellID> indices,
                          const std::uint8_t* buffer) override
  {
    if (meta::is_aligned<typename COMMUNICATEE::value_type>(buffer)) {
      const auto* source = reinterpret_cast<const typename COMMUNICATEE::value_type*>(buffer);
      for (CellID index : indices) {
        _communicatee[index] = *(source++);
      }
    } else {
      const std::uint8_t* source = buffer;
      for (CellID index : indices) {
        std::memcpy(reinterpret_cast<void*>(&_communicatee[index]),
                    source,
                    sizeof(typename COMMUNICATEE::value_type));
        source += sizeof(typename COMMUNICATEE::value_type);
      }
    }
    return indices.size() * sizeof(typename COMMUNICATEE::value_type);
  }

  //// ADDITIONAL NON OVERWITTEN CALLS: Removing the necessity to provide indices.
  ///- Here, intendet to be used for std::array (for which std::tuple_size works),
  ///  for other COMMUNICATEE types, a different static size retrieval should
  ///  be implemented.
  ///- Actually, vector might also work here, however is unfortunately rather
  ///  passed to ConcreteCommunicatable<std::vector<COLUMN>>
  ///- TODO: Include properly into complete framework, if found usefull.

  /// Get serialized size for complete data
  std::size_t size() const
  {
    return std::tuple_size<COMMUNICATEE>::value * sizeof(typename COMMUNICATEE::value_type);
  }

  /// Serialize complete data
  std::size_t serialize(std::uint8_t* buffer) const
  {
    std::size_t noI = std::tuple_size<COMMUNICATEE>::value;
    if (meta::is_aligned<typename COMMUNICATEE::value_type>(buffer)) {
      auto* target = reinterpret_cast<typename COMMUNICATEE::value_type*>(buffer);
      for (CellID index=0; index<noI; ++index) {
        *(target++) = _communicatee[index];
      }
    } else {
      std::uint8_t* target = buffer;
      for (CellID index=0; index<noI; ++index) {
        std::memcpy(target,
                    reinterpret_cast<const void*>(&_communicatee[index]),
                    sizeof(typename COMMUNICATEE::value_type));
        target += sizeof(typename COMMUNICATEE::value_type);
      }
    }
    return noI * sizeof(typename COMMUNICATEE::value_type);
  }

  /// Deserialize complete data
  std::size_t deserialize(const std::uint8_t* buffer)
  {
    std::size_t noI = std::tuple_size<COMMUNICATEE>::value;
    if (meta::is_aligned<typename COMMUNICATEE::value_type>(buffer)) {
      const auto* source = reinterpret_cast<const typename COMMUNICATEE::value_type*>(buffer);
      for (CellID index=0; index<noI; ++index) {
        _communicatee[index] = *(source++);
      }
    } else {
      const std::uint8_t* source = buffer;
      for (CellID index=0; index<noI; ++index) {
        std::memcpy(reinterpret_cast<void*>(&_communicatee[index]),
                    source,
                    sizeof(typename COMMUNICATEE::value_type));
        source += sizeof(typename COMMUNICATEE::value_type);
      }
    }
    return noI * sizeof(typename COMMUNICATEE::value_type);
  }

};

} // namespace olb

// Deferred here (see note at the top include block): pulls descriptor.h and the
// platform column.h chain, which is only safe once ConstSpan and the
// ConcreteCommunicatable primary template above are defined.
#include "communication/directionalOverlap.h"

namespace olb {

template <typename COLUMN>
class ConcreteCommunicatable<std::vector<COLUMN>> final : public Communicatable {
private:
  std::vector<COLUMN>& _vector;

public:
  ConcreteCommunicatable(std::vector<COLUMN>& vector):
    _vector{vector} { }

  /// Get serialized size for data at locations `indices`
  std::size_t size(ConstSpan<CellID> indices) const override
  {
    std::size_t size = 0;
    for (unsigned iD=0; iD < _vector.size(); ++iD) {
      size += ConcreteCommunicatable<COLUMN>(_vector[iD]).size(indices);
    }
    return size;
  }

  /// Serialize data at locations `indices` to `buffer`
  std::size_t serialize(ConstSpan<CellID> indices,
                        std::uint8_t* buffer) const override
  {
    std::size_t size = ConcreteCommunicatable<COLUMN>(_vector[0]).size(indices);
    std::uint8_t* curr = buffer;
    #ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for schedule(static,1)
    #endif
    for (unsigned iD=0; iD < _vector.size(); ++iD) {
      ConcreteCommunicatable<COLUMN>(_vector[iD]).serialize(indices, curr + iD*size);
    }
    return _vector.size() * size;
  }

  /// Deserialize data at locations `indices` to `buffer`
  std::size_t deserialize(ConstSpan<CellID> indices,
                          const std::uint8_t* buffer) override
  {
    std::size_t size = ConcreteCommunicatable<COLUMN>(_vector[0]).size(indices);
    const std::uint8_t* curr = buffer;
    #ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for schedule(static,1)
    #endif
    for (unsigned iD=0; iD < _vector.size(); ++iD) {
      ConcreteCommunicatable<COLUMN>(_vector[iD]).deserialize(indices, curr + iD*size);
    }
    return _vector.size() * size;
  }

  std::size_t serializeMasked(ConstSpan<CellID> indices,
                              const overlap::MaskGroupLayout& layout,
                              std::uint8_t* buffer) const override
  {
    return overlap::serializeMasked(_vector, indices, layout, buffer);
  }

  std::size_t deserializeMasked(ConstSpan<CellID> indices,
                                const overlap::MaskGroupLayout& layout,
                                ConstSpan<std::uint64_t> masks,
                                const std::uint8_t* buffer) override
  {
    const std::size_t size = overlap::deserializeMasked(_vector, indices, layout, buffer);
#ifdef FEATURE_DIRECTIONAL_CANARY
    if (overlap::canaryEnabled()) {
      util::writeDirectionalCanary(_vector, _vector.size(), indices, masks, layout.componentsPerDirection);
    }
#endif
    return size;
  }

};

template <typename COMMUNICATEE>
class MultiConcreteCommunicatable final : public Communicatable {
private:
  COMMUNICATEE& _communicatee;
  const std::vector<std::type_index>& _fields;

public:
  MultiConcreteCommunicatable(COMMUNICATEE& communicatee,
                              const std::vector<std::type_index>& fields):
    _communicatee{communicatee},
    _fields{fields} { }

  /// Get serialized size for data at locations `indices`
  std::size_t size(ConstSpan<CellID> indices) const override
  {
    std::size_t size = 0;
    for (auto& field : _fields) {
      size += _communicatee.getCommunicatable(field).size(indices);
    }
    return size;
  }

  /// Serialize data at locations `indices` to `buffer`
  std::size_t serialize(ConstSpan<CellID> indices,
                        std::uint8_t* buffer) const override
  {
    std::uint8_t* curr = buffer;
    for (auto& field : _fields) {
      curr += _communicatee.getCommunicatable(field).serialize(indices, curr);
    }
    return curr - buffer;
  }

  /// Deserialize data at locations `indices` to `buffer`
  std::size_t deserialize(ConstSpan<CellID> indices,
                          const std::uint8_t* buffer) override
  {
    const std::uint8_t* curr = buffer;
    for (auto& field : _fields) {
      curr += _communicatee.getCommunicatable(field).deserialize(indices, curr);
    }
    return curr - buffer;
  }
};

namespace overlap {

/// Multi-field Communicatable with component-masked segments
/**
 * Drop-in alternative to MultiConcreteCommunicatable for exchanges with
 * negotiated per-cell masks: fields serialize sequentially in request
 * order exactly as there, except that every masked field uses its own
 * MaskGroupLayout segment (see overlap::maskedFieldTraits for the per-field
 * mask-bit granularity). The masks must be index-aligned with every indices
 * argument passed to size/serialize/deserialize (i.e. the negotiated cell
 * list). All field access dispatches through the communicatee's
 * Communicatable interface and is therefore safe across platform
 * translation units.
 **/
template <typename COMMUNICATEE>
class MaskedMultiCommunicatable final : public Communicatable {
private:
  /// Per-field masks and the buffer layout they induce
  struct MaskedField {
    std::vector<std::uint64_t> masks;
    MaskGroupLayout            layout;
  };

  COMMUNICATEE& _communicatee;
  const std::vector<std::type_index> _fields;
  std::map<std::type_index, MaskedField> _masked;

  const MaskedField* maskedFor(std::type_index field) const {
    auto it = _masked.find(field);
    return it != _masked.end() ? &it->second : nullptr;
  }

public:
  MaskedMultiCommunicatable(COMMUNICATEE& communicatee,
                            const std::vector<std::type_index>& fields,
                            const overlap::FieldMasks& fieldMasks):
    _communicatee{communicatee},
    _fields(fields)
  {
    // Each masked segment packs at its field's stored scalar (half under
    // store_precision), matching serializeMasked -- not sizeof(value_t).
    using T = typename COMMUNICATEE::value_t;
    using DESCRIPTOR = typename COMMUNICATEE::descriptor_t;
    for (const auto& [field, masks] : fieldMasks) {
      if (const auto traits = maskedFieldTraits<T,DESCRIPTOR>(field)) {
        _masked.emplace(field, MaskedField{
          masks,
          MaskGroupLayout::compute(masks, traits->elementSize, traits->componentsPerDirection)
        });
      }
    }
  }

  std::size_t size(ConstSpan<CellID> indices) const override
  {
    std::size_t size = 0;
    for (auto& field : _fields) {
      if (const auto* m = maskedFor(field)) {
        size += m->layout.size;
      } else {
        size += _communicatee.getCommunicatable(field).size(indices);
      }
    }
    return size;
  }

  std::size_t serialize(ConstSpan<CellID> indices,
                        std::uint8_t* buffer) const override
  {
    std::uint8_t* curr = buffer;
    for (auto& field : _fields) {
      if (const auto* m = maskedFor(field)) {
        curr += _communicatee.getCommunicatable(field).serializeMasked(indices, m->layout, curr);
      } else {
        curr += _communicatee.getCommunicatable(field).serialize(indices, curr);
      }
    }
    return curr - buffer;
  }

  std::size_t deserialize(ConstSpan<CellID> indices,
                          const std::uint8_t* buffer) override
  {
    const std::uint8_t* curr = buffer;
    for (auto& field : _fields) {
      if (const auto* m = maskedFor(field)) {
        curr += _communicatee.getCommunicatable(field).deserializeMasked(
          indices, m->layout, ConstSpan<std::uint64_t>(m->masks), curr);
      } else {
        curr += _communicatee.getCommunicatable(field).deserialize(indices, curr);
      }
    }
    return curr - buffer;
  }

};

/// Construct the matching multi-field Communicatable for an exchange
/**
 * Masked when per-field masks were negotiated (lattice exchanges only),
 * the plain dense MultiConcreteCommunicatable otherwise.
 **/
template <typename BLOCK>
std::unique_ptr<Communicatable> makeMaskedMultiCommunicatable(
  BLOCK& block,
  const std::vector<std::type_index>& fields,
  const overlap::FieldMasks& fieldMasks)
{
  if constexpr (requires { typename BLOCK::value_t; typename BLOCK::descriptor_t; }) {
    const bool anyMasked = std::any_of(fieldMasks.begin(), fieldMasks.end(), [&](const auto& fm) {
      return std::find(fields.begin(), fields.end(), fm.first) != fields.end();
    });
    if (anyMasked) {
      return std::make_unique<MaskedMultiCommunicatable<BLOCK>>(block, fields, fieldMasks);
    }
  }
  return std::make_unique<MultiConcreteCommunicatable<BLOCK>>(block, fields);
}

}

}

#endif
