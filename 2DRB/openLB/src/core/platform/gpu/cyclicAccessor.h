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

#ifndef GPU_CYCLIC_ACCESSOR_H
#define GPU_CYCLIC_ACCESSOR_H

#include <cstddef>

#include "core/platform/platform.h"

namespace olb {

namespace gpu {

/// Addressing strategy for the cyclic population store (the streaming ring)
/**
 * - VirtualMemory: a single contiguous pointer into a double-mapped virtual
 *   address range. Wrap-around is handled by the mapping, so reads are
 *   branch-free. Requires the backend VMM API (cu/hipMemAddressReserve etc.).
 * - Branching: two pointers (head and wrapped tail) plus the split index; a
 *   single predicated select per access. Needs only a plain device allocation.
 *
 * Selected per backend by gpu::cyclicStrategyFor(Platform); override globally
 * with the FEATURE_GPU_CYCLIC_{BRANCHING,VMEM} build features (see config.git.mk).
 **/
enum struct CyclicStrategy {
  VirtualMemory,
  Branching
};

/// Cyclic population store strategy selected for a given GPU backend.
/**
 * Default: VirtualMemory on CUDA, Branching on HIP (which avoids the VMM driver
 * dependency). Override globally at build time via FEATURE_GPU_CYCLIC_BRANCHING
 * or FEATURE_GPU_CYCLIC_VMEM. Consumed as a constexpr template argument for
 * CyclicAccessor<T,STRATEGY> and as an `if constexpr` condition in the backend
 * column storage, so no preprocessor token plumbing is needed.
 **/
constexpr CyclicStrategy cyclicStrategyFor(Platform platform) {
#if defined(FEATURE_GPU_CYCLIC_BRANCHING)
  return CyclicStrategy::Branching;
#elif defined(FEATURE_GPU_CYCLIC_VMEM)
  return CyclicStrategy::VirtualMemory;
#else
  return platform == Platform::GPU_HIP ? CyclicStrategy::Branching
                                       : CyclicStrategy::VirtualMemory;
#endif
}

template <typename T, CyclicStrategy STRATEGY>
struct CyclicAccessor;

/// Branch-free accessor over a double-mapped virtual-memory ring
template <typename T>
struct CyclicAccessor<T,CyclicStrategy::VirtualMemory> {
  T* ptr;  ///< base + shift, valid for [0, count) thanks to the second mapping

  any_platform T& operator[](std::size_t i) const {
    return ptr[i];
  }
  any_platform T get(std::size_t i) const {
    return ptr[i];
  }
};

/// Two-pointer accessor wrapping at `remainder` without a virtual mapping
template <typename T>
struct CyclicAccessor<T,CyclicStrategy::Branching> {
  T* start0;                ///< base + shift (the un-wrapped head)
  T* start1;                ///< base - (count - shift) (the wrapped tail)
  std::ptrdiff_t remainder; ///< last index served from start0

  any_platform T& operator[](std::size_t i) const {
    return (static_cast<std::ptrdiff_t>(i) > remainder) ? start1[i] : start0[i];
  }
  any_platform T get(std::size_t i) const {
    return (static_cast<std::ptrdiff_t>(i) > remainder) ? start1[i] : start0[i];
  }
};

}

}

#endif
