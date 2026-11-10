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

#ifndef PLATFORM_GPU_SYCL_CANARY_H
#define PLATFORM_GPU_SYCL_CANARY_H

#ifdef FEATURE_DIRECTIONAL_CANARY

/// Directional-canary test scaffolding (see overlap::canaryEnabled)
/**
 * Included by communicator.hh after its kernel and dispatch machinery;
 * not a standalone header.
 **/

namespace olb {

namespace gpu {

namespace sycl {

/// Poison untransmitted slots of one masked FIELD with quiet NaN
/**
 * Overwrites the mask-complement slots of every exchanged cell so implicit
 * consumers of the dense halo fail loudly (see overlap::canaryEnabled).
 * Ordered after the scatter/copy that wrote the transmitted slots.
 **/
template <typename FIELD, typename BLOCK>
::sycl::event async_poison_untransmitted(
  ::sycl::queue* stream,
  BLOCK& block,
  gpu::sycl::Column<CellID>& cells,
  gpu::sycl::Column<std::uint64_t>& masks,
  unsigned componentsPerDirection,
  const std::vector<::sycl::event>& dependencies)
{
  using T = typename BLOCK::value_t;
  using S = ::olb::detail::field_scalar_t<T, typename BLOCK::descriptor_t, FIELD>;
  if constexpr (std::is_integral_v<S>) {
    return ::sycl::event{};
  } else {
    constexpr unsigned nComponents = BLOCK::descriptor_t::template size<FIELD>();
    auto field = propagatableAccessors<FIELD>(block);
    const CellID* cell_ptr = cells.deviceData();
    const std::uint64_t* mask_ptr = masks.deviceData();
    const std::size_t nCells = cells.size();
    return stream->submit([&](::sycl::handler& cgh) {
      cgh.depends_on(dependencies);
      cgh.parallel_for(::sycl::range<1>(nCells), [=](::sycl::id<1> idx) {
        std::size_t iIndex = idx[0];
        const CellID cell = cell_ptr[iIndex];
        const std::uint64_t mask = mask_ptr[iIndex];
        for (unsigned iPop=0; iPop < nComponents; ++iPop) {
          if (!(mask & (std::uint64_t{1} << (iPop / componentsPerDirection)))) {
            field[iPop][cell] = std::numeric_limits<S>::quiet_NaN();
          }
        }
      });
    });
  }
}

/// Per-field device copies of the negotiated masks (populated only when enabled)
using FieldCanaryMasks = std::map<std::type_index, gpu::sycl::Column<std::uint64_t>>;

inline FieldCanaryMasks canaryMasksFor(const std::vector<std::type_index>& fields,
                                       const overlap::FieldMasks& masks)
{
  FieldCanaryMasks deviceMasks;
  if (overlap::canaryEnabled()) {
    auto& q = device::getQueue();
    for (const auto& [field, mask] : masks) {
      if (std::find(fields.begin(), fields.end(), field) != fields.end() && !mask.empty()) {
        auto it = deviceMasks.emplace(std::piecewise_construct,
                                      std::forward_as_tuple(field),
                                      std::forward_as_tuple(mask.size())).first;
        q.memcpy(it->second.deviceData(), mask.data(), mask.size()*sizeof(std::uint64_t));
      }
    }
    q.wait();
  }
  return deviceMasks;
}

/// Poison the untransmitted halo of every masked propagatable segment
/**
 * Ordered after `dependencies` (the scatter resp. copy events); the returned
 * events transitively cover them, so they replace the completion token.
 **/
template <typename BLOCK>
std::vector<::sycl::event> async_poison_segments(
  ::sycl::queue* stream, BLOCK& block,
  const overlap::MultiFieldPartition& partition,
  gpu::sycl::Column<CellID>& cells,
  FieldCanaryMasks& canaryMasks,
  unsigned componentsPerDirection,
  const std::vector<::sycl::event>& dependencies)
{
  std::vector<::sycl::event> events;
  for (const auto& seg : partition.segments) {
    if (seg.kind == overlap::FieldSegment::Kind::MaskedPropagatable) {
      dispatchMaskedField<BLOCK>(seg.field, [&](auto id) {
        events.push_back(async_poison_untransmitted<typename decltype(id)::type>(
          stream, block, cells, canaryMasks.at(seg.field), componentsPerDirection, dependencies));
      });
    }
  }
  return events;
}

}

}

}

#endif

#endif
