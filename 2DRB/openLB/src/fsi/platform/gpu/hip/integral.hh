/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
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

#ifndef FSI_PLATFORM_GPU_HIP_INTEGRAL_HH
#define FSI_PLATFORM_GPU_HIP_INTEGRAL_HH

#include "integral.h"

#include <thrust/partition.h>
#include <thrust/sort.h>
#include <thrust/reduce.h>
#include <thrust/execution_policy.h>

namespace olb {

namespace gpu {

namespace hip {

template <typename F, std::size_t... INDICES>
auto make_thrust_tuple_f(F&& f, std::index_sequence<INDICES...>) {
  return thrust::make_tuple(f(INDICES)...);
}

template <unsigned D, typename F>
auto make_thrust_tuple_f(F&& f) {
  return make_thrust_tuple_f<F>(std::forward<F&&>(f), std::make_index_sequence<D>{});
}

template <typename A, std::size_t... AIs>
auto thrust_tuple_cat(A& a, std::index_sequence<AIs...>) {
  return thrust::make_tuple(thrust::get<AIs>(a)...);
}

template <typename A, typename B, std::size_t... AIs, std::size_t... BIs>
auto thrust_tuple_cat(A& a, B& b, std::index_sequence<AIs...>, std::index_sequence<BIs...>) {
  return thrust::make_tuple(thrust::get<AIs>(a)...,
                            thrust::get<BIs>(b)...);
}

template <typename A, typename B, typename C, std::size_t... AIs, std::size_t... BIs, std::size_t... CIs>
auto thrust_tuple_cat(A& a, B& b, C& c, std::index_sequence<AIs...>, std::index_sequence<BIs...>, std::index_sequence<CIs...>) {
  return thrust::make_tuple(thrust::get<AIs>(a)...,
                            thrust::get<BIs>(b)...,
                            thrust::get<CIs>(c)...);
}

template <typename... ARGS>
auto thrust_tuple_cat(ARGS&&... args) {
  return thrust_tuple_cat(args...,
                          std::make_index_sequence<thrust::tuple_size<ARGS>::value>{}...);
}

template <typename T, unsigned D>
auto make_thrust_tuple_of(Vector<T,D>& f) {
  return gpu::hip::make_thrust_tuple_f<D>([&](unsigned iDim) {
    return f[iDim];
  });
}

template <typename T>
auto make_thrust_tuple_of(T& f) {
  return thrust::make_tuple(f);
}

template <typename VECTOR>
auto make_thrust_tuple_of_device_data(VECTOR& f) {
  return gpu::hip::make_thrust_tuple_f<std::decay_t<VECTOR>::d>([&](unsigned iDim) {
    return f[iDim].deviceData();
  });
}

struct first_tuple_element_non_zero {
  template <typename TUPLE>
  bool operator()(const TUPLE& t) __host__ __device__ {
    return thrust::get<0>(t) != 0;
  }
};

struct sum_of_thrust_tuple {
  template <typename TUPLE, std::size_t... INDICES>
  auto operator()(const TUPLE& t1, const TUPLE& t2,
                  std::index_sequence<INDICES...>) __host__ __device__ {
    return thrust::make_tuple((thrust::get<INDICES>(t1) + thrust::get<INDICES>(t2))...);
  }

  template <typename TUPLE>
  auto operator()(const TUPLE& t1, const TUPLE& t2) __host__ __device__ {
    return operator()(t1, t2, std::make_index_sequence<thrust::tuple_size<TUPLE>::value>{});
  }
};

template <typename T, typename OUT_PTRS, typename... FIELDS>
__global__ void zero_elements_kernel(int* out_tags, OUT_PTRS out_ptrs, std::size_t nElements) {
  std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= nElements) return;

  out_tags[idx] = 0;

  auto zero_field = [&]<std::size_t I>() {
    auto out_f_ptrs = thrust::get<I>(out_ptrs);
    auto zero_dim = [&]<std::size_t D>() {
      thrust::get<D>(out_f_ptrs)[idx] = T{};
    };
    constexpr std::size_t out_f_ptrs_size = thrust::tuple_size<decltype(out_f_ptrs)>::value;
    [&]<std::size_t... Ds>(std::index_sequence<Ds...>) {
      (zero_dim.template operator()<Ds>(), ...);
    }(std::make_index_sequence<out_f_ptrs_size>{});
  };

  [&]<std::size_t... Is>(std::index_sequence<Is...>) {
    (zero_field.template operator()<Is>(), ...);
  }(std::make_index_sequence<sizeof...(FIELDS)>{});
}

template <typename T, typename IN_PTRS, typename OUT_PTRS, typename... FIELDS>
__global__ void apply_elements_kernel(int* t_data, int* out_tags, IN_PTRS in_ptrs, OUT_PTRS out_ptrs, std::size_t nCells) {
  std::size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= nCells) return;

  int tag = t_data[idx];
  if (tag > 0) {
    int out_idx = tag - 1;
    out_tags[out_idx] = tag;

    auto apply_field = [&]<std::size_t I>() {
      auto in_f_ptrs = thrust::get<I>(in_ptrs);
      auto out_f_ptrs = thrust::get<I>(out_ptrs);

      auto apply_dim = [&]<std::size_t D>() {
        atomicAdd(&(thrust::get<D>(out_f_ptrs)[out_idx]), thrust::get<D>(in_f_ptrs)[idx]);
      };

      constexpr std::size_t in_f_ptrs_size = thrust::tuple_size<decltype(in_f_ptrs)>::value;
      [&]<std::size_t... Ds>(std::index_sequence<Ds...>) {
        (apply_dim.template operator()<Ds>(), ...);
      }(std::make_index_sequence<in_f_ptrs_size>{});
    };

    [&]<std::size_t... Is>(std::index_sequence<Is...>) {
      (apply_field.template operator()<Is>(), ...);
    }(std::make_index_sequence<sizeof...(FIELDS)>{});

    t_data[idx] = 0;
  }
}


}

}

template <typename... FIELDS>
template <typename T, typename DESCRIPTOR>
void IntegratePorousElementFieldsO<FIELDS...>::type<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_HIP>>::setup(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_HIP>& blockLattice)
{
  blockLattice.template getData<OperatorParameters<IntegratePorousElementFieldsO>>();
}

template <typename... FIELDS>
template <typename T, typename DESCRIPTOR>
void IntegratePorousElementFieldsO<FIELDS...>::type<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_HIP>>::apply(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_HIP>& blockLattice){
  // Reduction runs on the default stream; wait for the collection that fed us.
  gpu::hip::device::synchronizePostProcessor();
// Two reduction strategies sum the per-cell porous-element fields by element tag:
//  - default: sort_by_key + reduce_by_key
//  - FEATURE_FSI_ATOMIC_REDUCTION: faster atomicAdd scatter
// Neither is bit-reproducible on GPU (unstable sort, order-dependent float sums).
// See config.git.mk.
#ifdef FEATURE_FSI_ATOMIC_REDUCTION
  const std::size_t nCells = blockLattice.getNcells();

  auto t_data = blockLattice.template getField<fields::fsi::REDUCED_ELEMENT_TAG>()[0].deviceData();

  auto& parameters = blockLattice.template getData<OperatorParameters<IntegratePorousElementFieldsO<FIELDS...>>>().parameters;
  auto out_tags = parameters.template get<fields::array_of<fields::fsi::REDUCED_ELEMENT_TAG>>();
  std::size_t nElements = parameters.template get<fields::fsi::REDUCED_ELEMENTS_COUNT>();

  auto in_ptrs = thrust::make_tuple(
    gpu::hip::make_thrust_tuple_of_device_data(blockLattice.template getField<FIELDS>())...
  );

  auto out_ptrs = thrust::make_tuple(
    gpu::hip::make_thrust_tuple_of(parameters.template get<fields::array_of<FIELDS>>())...
  );

  constexpr int blockSize = 256;

  if (nElements > 0){
    int numBlocksZero = (nElements + blockSize - 1)/blockSize;

    (gpu::hip::zero_elements_kernel<T, decltype(out_ptrs), FIELDS...>)<<<dim3(numBlocksZero), dim3(blockSize)>>>(
      out_tags, out_ptrs, nElements
    );
  }

  HIP_CHECK(hipDeviceSynchronize());

  if (nCells > 0) {
      int numBlocksApply = (nCells + blockSize - 1) / blockSize;
          (gpu::hip::apply_elements_kernel<T, decltype(in_ptrs), decltype(out_ptrs), FIELDS...>)<<<dim3(numBlocksApply), dim3(blockSize)>>>(
          t_data, out_tags, in_ptrs, out_ptrs, nCells
      );
  }

  HIP_CHECK(hipDeviceSynchronize());
#else
  auto& tag = blockLattice.template getField<fields::fsi::REDUCED_ELEMENT_TAG>();

  auto begin = thrust::make_zip_iterator(gpu::hip::thrust_tuple_cat(
                 gpu::hip::make_thrust_tuple_of_device_data(tag),
                 gpu::hip::make_thrust_tuple_of_device_data(
                   blockLattice.template getField<FIELDS>())...
               ));
  auto end = thrust::partition(thrust::device,
                               begin,
                               begin + blockLattice.getNcells(),
                               gpu::hip::first_tuple_element_non_zero{});
  std::size_t nActiveCells = end - begin;

  if (nActiveCells > 0) {
    // Group per-cell field values by FSI tag
    thrust::sort_by_key(thrust::device,
                        tag[0].deviceData(),
                        tag[0].deviceData() + nActiveCells,
                        thrust::make_zip_iterator(gpu::hip::thrust_tuple_cat(
                          gpu::hip::make_thrust_tuple_of_device_data(
                            blockLattice.template getField<FIELDS>())...
                        )));
  }

  // Sum up grouped force values
  auto& parameters = blockLattice.template getData<OperatorParameters<IntegratePorousElementFieldsO>>().parameters;
  if (nActiveCells > 0) {
    auto reducedElementTag = parameters.template get<fields::array_of<fields::fsi::REDUCED_ELEMENT_TAG>>();
    const std::size_t nKeys = thrust::reduce_by_key(
      thrust::device,
      tag[0].deviceData(),
      tag[0].deviceData() + nActiveCells,
      thrust::make_zip_iterator(gpu::hip::thrust_tuple_cat(
        gpu::hip::make_thrust_tuple_of_device_data(
          blockLattice.template getField<FIELDS>())...
      )),
      reducedElementTag,
      thrust::make_zip_iterator(gpu::hip::thrust_tuple_cat(
        gpu::hip::make_thrust_tuple_of(parameters.template get<fields::array_of<FIELDS>>())...
      )),
      thrust::equal_to<T>{},
      gpu::hip::sum_of_thrust_tuple{}
    ).first - reducedElementTag;
    parameters.template set<fields::fsi::REDUCED_ELEMENTS_COUNT>(nKeys);
  } else {
    parameters.template set<fields::fsi::REDUCED_ELEMENTS_COUNT>(0);
  }

  if (nActiveCells > 0) {
    thrust::fill(thrust::device,
                 tag[0].deviceData(),
                 tag[0].deviceData() + nActiveCells,
                 0);
  }
#endif
}

}

#endif
