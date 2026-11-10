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

#ifndef PLATFORM_GPU_HIP_CANARY_H
#define PLATFORM_GPU_HIP_CANARY_H

#ifdef FEATURE_DIRECTIONAL_CANARY

/// Directional-canary test scaffolding (see overlap::canaryEnabled)
/**
 * Included by communicator.hh after its kernel and dispatch machinery;
 * not a standalone header.
 **/

namespace olb {

namespace gpu {

namespace hip {

namespace kernel {

/// HIP kernel poisoning untransmitted population slots (directional canary)
/**
 * Overwrites the mask-complement slots of every exchanged cell with quiet
 * NaN so implicit consumers of the dense halo fail loudly (see
 * overlap::canaryEnabled). Fully dropped cells (mask 0) poison all q.
 **/
template <typename ACCESSORS>
void poison_untransmitted(ACCESSORS field, unsigned nComponents,
                          unsigned componentsPerDirection,
                          const CellID* cells,
                          const std::uint64_t* masks, std::size_t nCells) __global__ {
  const std::size_t iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nCells)) {
    return;
  }
  const CellID cell = cells[iIndex];
  const std::uint64_t mask = masks[iIndex];
  for (unsigned iPop=0; iPop < nComponents; ++iPop) {
    if (!(mask & (std::uint64_t{1} << (iPop / componentsPerDirection)))) {
      field[iPop][cell] = nan("");
    }
  }
}

}

/// Launch the directional-canary poison kernel for one masked FIELD
template <typename FIELD, typename BLOCK>
void async_poison_untransmitted(hipStream_t stream,
                                BLOCK& block,
                                const thrust::device_vector<CellID>& cells,
                                const thrust::device_vector<std::uint64_t>& masks,
                                unsigned componentsPerDirection) {
  using S = ::olb::detail::field_scalar_t<typename BLOCK::value_t, typename BLOCK::descriptor_t, FIELD>;
  if constexpr (!std::is_integral_v<S>) {
    const auto block_size = DEFAULT_HIP_WARP_SIZE;
    const auto block_count = (cells.size() + block_size - 1) / block_size;
    kernel::poison_untransmitted<<<block_count,block_size,0,stream>>>(
      propagatableAccessors<FIELD>(block),
      BLOCK::descriptor_t::template size<FIELD>(),
      componentsPerDirection,
      cells.data().get(), masks.data().get(), cells.size());
    device::check();
  }
}

/// Per-field device copies of the negotiated masks (populated only when enabled)
using FieldCanaryMasks = std::map<std::type_index, thrust::device_vector<std::uint64_t>>;

inline FieldCanaryMasks canaryMasksFor(const std::vector<std::type_index>& fields,
                                       const overlap::FieldMasks& masks)
{
  FieldCanaryMasks deviceMasks;
  if (overlap::canaryEnabled()) {
    for (const auto& [field, mask] : masks) {
      if (std::find(fields.begin(), fields.end(), field) != fields.end()) {
        deviceMasks.emplace(field, mask);
      }
    }
  }
  return deviceMasks;
}

/// Poison the untransmitted halo of every masked propagatable segment
template <typename BLOCK>
void async_poison_segments(hipStream_t stream, BLOCK& block,
                           const overlap::MultiFieldPartition& partition,
                           const thrust::device_vector<CellID>& cells,
                           const FieldCanaryMasks& masks,
                           unsigned componentsPerDirection) {
  for (const auto& seg : partition.segments) {
    if (seg.kind == overlap::FieldSegment::Kind::MaskedPropagatable) {
      dispatchMaskedField<BLOCK>(seg.field, [&](auto id) {
        async_poison_untransmitted<typename decltype(id)::type>(stream, block, cells, masks.at(seg.field), componentsPerDirection);
      });
    }
  }
}

}

}

}

#endif

#endif
