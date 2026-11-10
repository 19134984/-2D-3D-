/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Adrian Kummerlaender
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

#ifndef PLATFORM_DEVICE_MEMORY_H
#define PLATFORM_DEVICE_MEMORY_H

#include <cstddef>

#include "platform.h"

#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

namespace olb {

// Forward declarations of the per-backend total-memory queries (defined in the
// respective device.hh). Declared here unconditionally under their platform
// guard so that the constexpr dispatch below can name them.
#ifdef PLATFORM_GPU_CUDA
namespace gpu { namespace cuda { namespace device { std::size_t getDeviceMemorySize(); } } }
#endif
#ifdef PLATFORM_GPU_HIP
namespace gpu { namespace hip { namespace device { std::size_t getDeviceMemorySize(); } } }
#endif
#ifdef PLATFORM_GPU_SYCL
namespace gpu { namespace sycl { namespace device { std::size_t getDeviceMemorySize(); } } }
#endif

/// Total physical memory of the host in bytes (0 if it cannot be determined)
inline std::size_t getHostTotalMemory()
{
#if defined(_SC_PHYS_PAGES) && defined(_SC_PAGE_SIZE)
  const long pages    = sysconf(_SC_PHYS_PAGES);
  const long pageSize = sysconf(_SC_PAGE_SIZE);
  if (pages > 0 && pageSize > 0) {
    return static_cast<std::size_t>(pages) * static_cast<std::size_t>(pageSize);
  }
#endif
  return 0;
}

/// Total physical memory of the device backing PLATFORM, in bytes
/**
 * Returns host RAM for CPU platforms and the device global memory size for GPU
 * platforms. Used to express per-block field storage as a fraction of the RAM
 * it lives on. Returns 0 when the quantity is unavailable.
 **/
template <Platform PLATFORM>
inline std::size_t getDeviceTotalMemory()
{
  if constexpr (isPlatformCPU(PLATFORM)) {
    return getHostTotalMemory();
  }
#ifdef PLATFORM_GPU_CUDA
  else if constexpr (PLATFORM == Platform::GPU_CUDA) {
    return gpu::cuda::device::getDeviceMemorySize();
  }
#endif
#ifdef PLATFORM_GPU_HIP
  else if constexpr (PLATFORM == Platform::GPU_HIP) {
    return gpu::hip::device::getDeviceMemorySize();
  }
#endif
#ifdef PLATFORM_GPU_SYCL
  else if constexpr (PLATFORM == Platform::GPU_SYCL) {
    return gpu::sycl::device::getDeviceMemorySize();
  }
#endif
  else {
    return 0;
  }
}

}

#endif
