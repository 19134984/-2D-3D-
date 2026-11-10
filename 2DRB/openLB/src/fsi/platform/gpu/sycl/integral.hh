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

#ifndef FSI_PLATFORM_GPU_SYCL_INTEGRAL_HH
#define FSI_PLATFORM_GPU_SYCL_INTEGRAL_HH

#include "integral.h"

#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/numeric>
#include <oneapi/dpl/iterator>

#include <tuple>
#include "core/platform/gpu/sycl/device.h"

namespace olb {

namespace gpu {

namespace sycl {

template <typename F, std::size_t... INDICES>
auto make_dpl_tuple_f(F&& f, std::index_sequence<INDICES...>) {
  return std::make_tuple(f(INDICES)...);
}

template <unsigned D, typename F>
auto make_dpl_tuple_f(F&& f) {
  return make_dpl_tuple_f<F>(std::forward<F&&>(f), std::make_index_sequence<D>{});
}

template <typename A, std::size_t... AIs>
auto dpl_tuple_cat(A& a, std::index_sequence<AIs...>) {
  return std::make_tuple(std::get<AIs>(a)...);
}

template <typename A, typename B, std::size_t... AIs, std::size_t... BIs>
auto dpl_tuple_cat(A& a, B& b, std::index_sequence<AIs...>, std::index_sequence<BIs...>) {
  return std::make_tuple(std::get<AIs>(a)...,
                         std::get<BIs>(b)...);
}

template <typename A, typename B, typename C, std::size_t... AIs, std::size_t... BIs, std::size_t... CIs>
auto dpl_tuple_cat(A& a, B& b, C& c, std::index_sequence<AIs...>, std::index_sequence<BIs...>, std::index_sequence<CIs...>) {
  return std::make_tuple(std::get<AIs>(a)...,
                         std::get<BIs>(b)...,
                         std::get<CIs>(c)...);
}

template <typename... ARGS>
auto dpl_tuple_cat(ARGS&&... args) {
  return dpl_tuple_cat(args...,
                       std::make_index_sequence<std::tuple_size_v<std::remove_reference_t<ARGS>>>{}...);
}

template <typename T, unsigned D>
auto make_dpl_tuple_of(Vector<T,D>& f) {
  return gpu::sycl::make_dpl_tuple_f<D>([&](unsigned iDim) {
    return f[iDim];
  });
}

template <typename T>
auto make_dpl_tuple_of(T& f) {
  return std::make_tuple(f);
}

template <typename VECTOR>
auto make_dpl_tuple_of_device_data(VECTOR& f) {
  return gpu::sycl::make_dpl_tuple_f<VECTOR::d>([&](unsigned iDim) {
    return f[iDim].deviceData();
  });
}

struct first_tuple_element_non_zero {
  template <typename TUPLE>
  bool operator()(const TUPLE& t) const {
    return std::get<0>(t) != 0;
  }
};

struct sum_of_dpl_tuple {
  template <typename TUPLE1, typename TUPLE2, std::size_t... INDICES>
  auto operator()(const TUPLE1& t1, const TUPLE2& t2,
                  std::index_sequence<INDICES...>) const {
    return std::make_tuple((std::get<INDICES>(t1) + std::get<INDICES>(t2))...);
  }

  template <typename TUPLE1, typename TUPLE2>
  auto operator()(const TUPLE1& t1, const TUPLE2& t2) const {
    return operator()(t1, t2, std::make_index_sequence<std::tuple_size_v<std::remove_reference_t<TUPLE1>>>{});
  }
};

}

}

template <typename... FIELDS>
template <typename T, typename DESCRIPTOR>
void IntegratePorousElementFieldsO<FIELDS...>::type<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>::setup(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& blockLattice)
{
  blockLattice.template getData<OperatorParameters<IntegratePorousElementFieldsO>>();
}

template <typename... FIELDS>
template <typename T, typename DESCRIPTOR>
void IntegratePorousElementFieldsO<FIELDS...>::type<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>::apply(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& blockLattice)
{
  const std::size_t nCells = blockLattice.getNcells();
  auto queue = gpu::sycl::device::getQueue();

  auto t_data = blockLattice.template getField<fields::fsi::REDUCED_ELEMENT_TAG>()[0].deviceData();

  auto& parameters = blockLattice.template getData<OperatorParameters<IntegratePorousElementFieldsO<FIELDS...>>>().parameters;
  auto out_tags = parameters.template get<fields::array_of<fields::fsi::REDUCED_ELEMENT_TAG>>();
  std::size_t nElements = parameters.template get<fields::fsi::REDUCED_ELEMENTS_COUNT>();

  auto in_ptrs = std::make_tuple(
    gpu::sycl::make_dpl_tuple_of_device_data(blockLattice.template getField<FIELDS>())...
  );

  auto out_ptrs = std::make_tuple(
    gpu::sycl::make_dpl_tuple_of(parameters.template get<fields::array_of<FIELDS>>())...
  );

  auto zero_event = queue.submit([&](::sycl::handler& cgh) {
    cgh.parallel_for(::sycl::range<1>(nElements), [=](::sycl::id<1> idx) {
      out_tags[idx] = 0;

      auto zero_field = [&]<std::size_t I>() {
        auto out_f_ptrs = std::get<I>(out_ptrs);
        auto zero_dim = [&]<std::size_t D>() {
          std::get<D>(out_f_ptrs)[idx] = T{};
        };
        constexpr std::size_t out_f_ptrs_size = std::tuple_size_v<decltype(out_f_ptrs)>;
        [&]<std::size_t... Ds>(std::index_sequence<Ds...>) {
          (zero_dim.template operator()<Ds>(), ...);
        }(std::make_index_sequence<out_f_ptrs_size>{});
      };

      [&]<std::size_t... Is>(std::index_sequence<Is...>) {
        (zero_field.template operator()<Is>(), ...);
      }(std::make_index_sequence<sizeof...(FIELDS)>{});
    });
  });

  queue.submit([&](::sycl::handler& cgh) {
    cgh.depends_on(zero_event);
    cgh.parallel_for(::sycl::range<1>(nCells), [=](::sycl::id<1> idx) {
      int tag = t_data[idx];
      if (tag > 0) {
        int out_idx = tag - 1;
        out_tags[out_idx] = tag;

        auto apply_field = [&]<std::size_t I>() {
          auto in_f_ptrs = std::get<I>(in_ptrs);
          auto out_f_ptrs = std::get<I>(out_ptrs);

          auto apply_dim = [&]<std::size_t D>() {
            ::sycl::atomic_ref<T, ::sycl::memory_order::relaxed
                                , ::sycl::memory_scope::device
                                , ::sycl::access::address_space::global_space> ref(std::get<D>(out_f_ptrs)[out_idx]);
            ref += std::get<D>(in_f_ptrs)[idx];
          };

          constexpr std::size_t in_f_ptrs_size = std::tuple_size_v<decltype(in_f_ptrs)>;
          [&]<std::size_t... Ds>(std::index_sequence<Ds...>) {
            (apply_dim.template operator()<Ds>(), ...);
          }(std::make_index_sequence<in_f_ptrs_size>{});
        };

        [&]<std::size_t... Is>(std::index_sequence<Is...>) {
          (apply_field.template operator()<Is>(), ...);
        }(std::make_index_sequence<sizeof...(FIELDS)>{});

        t_data[idx] = 0;
      }
    });
  }).wait();
}

}

#endif
