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

#ifndef GPU_SYCL_STATISTICS_HH
#define GPU_SYCL_STATISTICS_HH

#include "statistics.h"
#include "device.h"
#include "core/fieldStatisticO.h"

#include <sycl/sycl.hpp>
#include <limits>

namespace olb {
namespace gpu {
namespace sycl {

template <typename T, unsigned D>
struct FieldStatReduce {
  const T*    values;
  std::size_t proj[D];
  int         core[D];
  int         padding;
  template <typename RC, typename RMin, typename RMax, typename RSum>
  void operator()(::sycl::item<1> it, RC& r_count, RMin& r_min, RMax& r_max, RSum& r_sum) const {
    const std::size_t idx = it.get_linear_id();
    int coord[D];
    if constexpr (D == 3) {
      coord[0] = idx / proj[0];
      const std::size_t rem = idx % proj[0];
      coord[1] = rem / proj[1];
      coord[2] = rem % proj[1];
    } else if constexpr (D == 2) {
      coord[0] = idx / proj[0];
      coord[1] = idx % proj[0];
    } else {
      coord[0] = idx;
    }
    for (unsigned i = 0; i < D; ++i) {
      if (coord[i] - padding < 0 || coord[i] - padding >= core[i]) {
        return;
      }
    }
    const T v = values[idx];
    r_count.combine(std::size_t{1});
    r_min.combine(v);
    r_max.combine(v);
    r_sum.combine(v);
  }
};

}
}
}

namespace olb {

template <typename T, typename DESCRIPTOR>
void StatisticsPostProcessor::type<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>::apply(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& blockLattice)
{
  if (!blockLattice.statisticsEnabled()) {
    return;
  }

  auto& statisticGenerated = blockLattice.template getField<descriptors::STATISTIC_GENERATED>();
  auto& statistic = blockLattice.template getField<descriptors::STATISTIC>();

  std::size_t nCells = blockLattice.getNcells();

  auto* d_generated = statisticGenerated[0].deviceData();
  T* d_rho = statistic[0].deviceData();
  T* d_uSqr = statistic[1].deviceData();

  std::size_t nCellsReduced = 0;
  T rhoSum  = 0;
  T uSqrSum = 0;
  T maxU    = std::numeric_limits<T>::lowest();
  T minRho  = std::numeric_limits<T>::max();
  T maxRho  = std::numeric_limits<T>::lowest();
  T minUSqr = std::numeric_limits<T>::max();
  T sumU    = 0;

  auto& queue = gpu::sycl::device::getQueue();

  {
    OLB_PROFILE_T("GPU_SYCL:Stats:Reduce", DESCRIPTOR);
    // Masked min/max/mean of rho and u via built-in SYCL reductions, split
    // into two kernels of four reductions each.
    // First kernel: count, sum(rho), sum(uSqr), max(uSqr).
    ::sycl::buffer<std::size_t, 1> buffer_nCells(&nCellsReduced, 1);
    ::sycl::buffer<T, 1> buffer_rho(&rhoSum, 1);
    ::sycl::buffer<T, 1> buffer_uSqr(&uSqrSum, 1);
    ::sycl::buffer<T, 1> buffer_maxU(&maxU, 1);
    queue.submit([&](::sycl::handler& cgh) {
      auto reduction_nCells = ::sycl::reduction(buffer_nCells, cgh, ::sycl::plus<std::size_t>());
      auto reduction_rho    = ::sycl::reduction(buffer_rho, cgh, ::sycl::plus<T>());
      auto reduction_uSqr   = ::sycl::reduction(buffer_uSqr, cgh, ::sycl::plus<T>());
      auto reduction_maxU   = ::sycl::reduction(buffer_maxU, cgh, ::sycl::maximum<T>());
      cgh.parallel_for(::sycl::range<1>(nCells),
                       reduction_nCells, reduction_rho, reduction_uSqr, reduction_maxU,
        [=](::sycl::item<1> it, auto& r_nCells, auto& r_rho, auto& r_uSqr, auto& r_maxU) {
          std::size_t idx = it.get_linear_id();
          if (d_generated[idx]) {
            r_nCells.combine(std::size_t{1});
            r_rho.combine(d_rho[idx]);
            r_uSqr.combine(d_uSqr[idx]);
            r_maxU.combine(d_uSqr[idx]);
          }
      });
    });
    // Second kernel: min(rho), max(rho), min(uSqr), sum(|u|)
    ::sycl::buffer<T, 1> buffer_minRho(&minRho, 1);
    ::sycl::buffer<T, 1> buffer_maxRho(&maxRho, 1);
    ::sycl::buffer<T, 1> buffer_minUSqr(&minUSqr, 1);
    ::sycl::buffer<T, 1> buffer_sumU(&sumU, 1);
    queue.submit([&](::sycl::handler& cgh) {
      auto reduction_minRho  = ::sycl::reduction(buffer_minRho, cgh, ::sycl::minimum<T>());
      auto reduction_maxRho  = ::sycl::reduction(buffer_maxRho, cgh, ::sycl::maximum<T>());
      auto reduction_minUSqr = ::sycl::reduction(buffer_minUSqr, cgh, ::sycl::minimum<T>());
      auto reduction_sumU    = ::sycl::reduction(buffer_sumU, cgh, ::sycl::plus<T>());
      cgh.parallel_for(::sycl::range<1>(nCells),
                       reduction_minRho, reduction_maxRho, reduction_minUSqr, reduction_sumU,
        [=](::sycl::item<1> it, auto& r_minRho, auto& r_maxRho, auto& r_minUSqr, auto& r_sumU) {
          std::size_t idx = it.get_linear_id();
          if (d_generated[idx]) {
            r_minRho.combine(d_rho[idx]);
            r_maxRho.combine(d_rho[idx]);
            r_minUSqr.combine(d_uSqr[idx]);
            r_sumU.combine(::sycl::sqrt(d_uSqr[idx]));
          }
      });
    });
  }

  if (nCellsReduced == 0) {
    maxU = 0;
  }

  typename LatticeStatistics<T,DESCRIPTOR>::Aggregatable statistics{
    .nCells   = nCellsReduced,
    .avRho    = rhoSum,
    .avEnergy = uSqrSum,
    .maxU     = maxU,
    .minRho   = minRho,
    .maxRho   = maxRho,
    .minUSqr  = minUSqr,
    .sumU     = sumU
  };
  blockLattice.getStatistics().incrementStats(statistics);

  blockLattice.getStatistics().reset();
}

/// Device reduction of a registered FIELD over the core cells (masked by isCore)
template <typename FIELD>
template <typename T, typename DESCRIPTOR>
void FieldStatisticO<FIELD>::type<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>::apply(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& block)
{
    if (!block.statisticsEnabled()) {
      return;
    }
    constexpr unsigned dim = DESCRIPTOR::d;
    const unsigned    d = DESCRIPTOR::template size<FIELD>();
    const std::size_t n = block.getNcells();
    auto& field = block.template getField<FIELD>();
    BlockStructureD<dim> structure = block;
    auto& queue = gpu::sycl::device::getQueue();

    gpu::sycl::FieldStatReduce<T,dim> reducer{};
    reducer.padding = structure.getPadding();
    auto coreExtent = structure.getCore();
    for (unsigned i = 0; i < dim; ++i) {
      reducer.core[i] = coreExtent[i];
    }
    if constexpr (dim == 3) {
      const std::size_t s1 = std::size_t(coreExtent[1]) + 2*reducer.padding;
      const std::size_t s2 = std::size_t(coreExtent[2]) + 2*reducer.padding;
      reducer.proj[0] = s1 * s2;
      reducer.proj[1] = s2;
    } else if constexpr (dim == 2) {
      reducer.proj[0] = std::size_t(coreExtent[1]) + 2*reducer.padding;
    }

    typename LatticeStatistics<T,DESCRIPTOR>::FieldResult r;
    r.components = d;
    r.min.resize(d);
    r.max.resize(d);
    r.sum.resize(d);
    r.count = 0;

    for (unsigned iD = 0; iD < d; ++iD) {
      reducer.values = field[iD].deviceData();
      std::size_t count = 0;
      T mn = std::numeric_limits<T>::max();
      T mx = std::numeric_limits<T>::lowest();
      T sm = 0;
      {
        ::sycl::buffer<std::size_t, 1> buffer_count(&count, 1);
        ::sycl::buffer<T, 1> buffer_min(&mn, 1);
        ::sycl::buffer<T, 1> buffer_max(&mx, 1);
        ::sycl::buffer<T, 1> buffer_sum(&sm, 1);
        queue.submit([&](::sycl::handler& cgh) {
          auto reduction_count = ::sycl::reduction(buffer_count, cgh, ::sycl::plus<std::size_t>());
          auto reduction_min   = ::sycl::reduction(buffer_min, cgh, ::sycl::minimum<T>());
          auto reduction_max   = ::sycl::reduction(buffer_max, cgh, ::sycl::maximum<T>());
          auto reduction_sum   = ::sycl::reduction(buffer_sum, cgh, ::sycl::plus<T>());
          cgh.parallel_for(::sycl::range<1>(n),
                           reduction_count, reduction_min, reduction_max, reduction_sum,
            reducer);
        });
      }
      r.min[iD] = mn;
      r.max[iD] = mx;
      r.sum[iD] = sm;
      if (iD == 0) {
        r.count = count;
      }
    }

    block.getStatistics().template setFieldResult<FIELD>(r);
}

}

#endif
