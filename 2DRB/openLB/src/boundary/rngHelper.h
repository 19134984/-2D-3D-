/*  Lattice Boltzmann sample, written in C++, using the OpenLB library
 *
 *  Copyright (C) 2026 Leonardo Dorneles
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

#ifndef RNG_HELPER_H
#define RNG_HELPER_H

#include <cstdint>
#include <random>

#ifdef PLATFORM_GPU_HIP
#include <hiprand/hiprand_kernel.h>
#endif

#if defined(OLB_PRECOMPILING_SYCL)
#include <oneapi/mkl/rng/device.hpp>
#endif

#if defined(__CUDACC__)
#include <curand_kernel.h>
#endif

namespace olb {

template <typename T>
struct CpuRNG {
    std::mt19937 engine;

    CpuRNG(std::uint64_t seed, std::size_t id)
        : engine(seed + id) {}

    T  uniform()  { return std::uniform_real_distribution<T> (0,1)(engine); }
    T uniform_size()   { return std::uniform_real_distribution<T>(0.5,1.5)(engine); }
};

#if defined(OLB_PRECOMPILING_SYCL)

template <typename T>
struct SyclRNG {
    oneapi::mkl::rng::device::philox4x32x10<1> engine;
    oneapi::mkl::rng::device::uniform<T>  distr     {0.0, 1.0};
    oneapi::mkl::rng::device::uniform<T>  sizeDistr {0.5, 1.5};

    SyclRNG(std::uint64_t seed, std::size_t id) : engine(seed, id) {}

    T uniform() { return oneapi::mkl::rng::device::generate(distr,     engine); }
    T uniform_size()   { return oneapi::mkl::rng::device::generate(sizeDistr, engine); }
};

template <typename T>
using PlatformRNG = SyclRNG<T>;

#elif defined(__HIPCC__) || defined(__CUDACC__)

#if defined(__HIPCC__)
using philox_state_t = hiprandStatePhilox4_32_10_t;
#define platform_rand_init          hiprand_init
#define platform_rand_uniform       hiprand_uniform
#define platform_rand_uniform_double hiprand_uniform_double
#else
using philox_state_t = curandStatePhilox4_32_10_t;
#define platform_rand_init          curand_init
#define platform_rand_uniform       curand_uniform
#define platform_rand_uniform_double curand_uniform_double
#endif

template <typename T>
struct GpuRNG {
    philox_state_t state;

    __device__ GpuRNG(std::uint64_t seed, std::size_t id) {
        platform_rand_init(seed, id, 0, &state);
    }

    __device__ float  uniform()  {
      if constexpr (std::is_same_v<T, double>) {
          return platform_rand_uniform_double(&state);
        } else {
          return platform_rand_uniform(&state);
        }
      }
    __device__ T uniform_size()   { return 0.5 + uniform(); }
};

#if defined(__HIP_DEVICE_COMPILE__) || defined(__CUDA_ARCH__)
template <typename T>
using PlatformRNG = GpuRNG<T>;
#else
template <typename T>
using PlatformRNG = CpuRNG<T>;
#endif
#else
template <typename T>
using PlatformRNG = CpuRNG<T>;
#endif

}

#endif
