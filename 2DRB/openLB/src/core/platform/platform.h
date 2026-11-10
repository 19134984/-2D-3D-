/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef PLATFORM_H
#define PLATFORM_H

#include <string>
#include <cstdint>
#include <stdexcept>

/// Top level namespace for all of OpenLB
namespace olb {

/// OpenLB execution targets
enum struct Platform : std::uint8_t {
  CPU_SISD, /// Basic scalar CPU
  CPU_SIMD, /// Vector CPU (AVX2 / AVX-512 collision)
  GPU_CUDA, /// GPU using CUDA
  GPU_HIP,  /// GPU using HIP
  GPU_SYCL  /// GPU using SYCL
};

/// Verifies requirements for using PLATFORM
template <Platform PLATFORM>
void checkPlatform();

/// OpenLB processing contexts
/**
 * Currently of no relevance for CPU_SISD and CPU_SIMD target platforms
 *
 * Used to control synchronization between mirrored device and host data
 * for non-host processed block lattices.
 *
 * Preliminary for first GPU release.
 **/
enum struct ProcessingContext {
  Evaluation, /// Data available on host for e.g. functor evaluation
  Simulation  /// Data available on device for evolving the simulation
};

namespace stage {

template <ProcessingContext CONTEXT>
struct PreContextSwitchTo { };

}

/// Define preprocessor macros for device-side functions, constant storage
// platform_constant_definition is inline as it is used for explicit variable
// template specializations in headers, which have external linkage and would
// otherwise violate the ODR once included in more than one translation unit.
#if defined(__CUDACC__) || defined(__HIPCC__)
  #define any_platform __device__ __host__
  #ifdef __CUDA_ARCH__
    #define platform_constant constexpr __constant__
    #define platform_constant_definition inline constexpr __constant__
  #else
    #define platform_constant constexpr
    #define platform_constant_definition inline constexpr
  #endif
#else
  #define any_platform
  #define platform_constant constexpr
  #define platform_constant_definition inline constexpr
#endif

/// Returns true if platform is equal to Platform::CPU_*
constexpr bool isPlatformCPU(Platform platform) {
  switch (platform) {
#ifdef PLATFORM_CPU_SISD
  case Platform::CPU_SISD:
    return true;
#endif
#ifdef PLATFORM_CPU_SIMD
  case Platform::CPU_SIMD:
    return true;
#endif
#ifdef PLATFORM_GPU_CUDA
  case Platform::GPU_CUDA:
    return false;
#endif
#ifdef PLATFORM_GPU_HIP
  case Platform::GPU_HIP:
    return false;
#endif
#ifdef PLATFORM_GPU_SYCL
  case Platform::GPU_SYCL:
    return false;
#endif
  default:
    throw std::invalid_argument("Invalid PLATFORM");
  }
}

/// Returns name of platform for profiling
constexpr const char* getName(Platform platform) {
  switch (platform) {
    case Platform::CPU_SISD: return "CPU_SISD";
    case Platform::CPU_SIMD: return "CPU_SIMD";
    case Platform::GPU_CUDA: return "GPU_CUDA";
    case Platform::GPU_HIP:  return "GPU_HIP";
    case Platform::GPU_SYCL: return "GPU_SYCL";
    default: return "UNKNOWN";
  }
}

/// CSE code-generation target a collision body is specialized for
/**
 * Spans the moments-storage strategy (whether rho/u are stored as fields versus
 * accessed via virtual dispatch / recomputed) and, in the future, hardware-specific
 * tuning. Selected per Platform by codegenTargetFor.
 **/
enum struct CodegenTarget {
  Generic,       /// Portable body; moments accessed via virtual dispatch / recompute
  StoredMoments  /// Stores rho/u as fields; for platforms without cheap device virtual dispatch
  // future: SyclSubgroup, CudaWarp, ...
};

/// Returns the CSE code-generation target appropriate for PLATFORM
/**
 * SYCL has no device virtual dispatch, so moments are always stored as fields.
 * Every other platform can recompute moments via virtual dispatch (Generic) but
 * may opt into stored moments through a per-platform FEATURE_<PLATFORM>_STORED_MOMENTS
 * build feature, useful e.g. for whole-lattice coupling readers. See config.git.mk.
 **/
constexpr CodegenTarget codegenTargetFor(Platform platform) {
  switch (platform) {
  case Platform::GPU_SYCL:
    return CodegenTarget::StoredMoments;
  case Platform::GPU_HIP:
#ifdef FEATURE_HIP_STORED_MOMENTS
    return CodegenTarget::StoredMoments;
#else
    return CodegenTarget::Generic;
#endif
  case Platform::GPU_CUDA:
#ifdef FEATURE_CUDA_STORED_MOMENTS
    return CodegenTarget::StoredMoments;
#else
    return CodegenTarget::Generic;
#endif
  case Platform::CPU_SIMD:
#ifdef FEATURE_CPU_SIMD_STORED_MOMENTS
    return CodegenTarget::StoredMoments;
#else
    return CodegenTarget::Generic;
#endif
  case Platform::CPU_SISD:
#ifdef FEATURE_CPU_SISD_STORED_MOMENTS
    return CodegenTarget::StoredMoments;
#else
    return CodegenTarget::Generic;
#endif
  default:
    return CodegenTarget::Generic;
  }
}

/// True if any platform enabled in this build resolves to the StoredMoments codegen target.
/** Gates compilation of moments-refresh machinery (RefreshMomentsO providers): stored-moment
 *  fields only exist on StoredMoments platforms, so on Generic-only builds (CPU / CUDA) the
 *  providers are dead weight and instantiating them just bloats compile time. **/
constexpr bool anyEnabledPlatformUsesStoredMoments() {
  bool result = false;
#ifdef PLATFORM_CPU_SISD
  result = result || codegenTargetFor(Platform::CPU_SISD) == CodegenTarget::StoredMoments;
#endif
#ifdef PLATFORM_CPU_SIMD
  result = result || codegenTargetFor(Platform::CPU_SIMD) == CodegenTarget::StoredMoments;
#endif
#ifdef PLATFORM_GPU_CUDA
  result = result || codegenTargetFor(Platform::GPU_CUDA) == CodegenTarget::StoredMoments;
#endif
#ifdef PLATFORM_GPU_HIP
  result = result || codegenTargetFor(Platform::GPU_HIP) == CodegenTarget::StoredMoments;
#endif
#ifdef PLATFORM_GPU_SYCL
  result = result || codegenTargetFor(Platform::GPU_SYCL) == CodegenTarget::StoredMoments;
#endif
  return result;
}

}

#endif
