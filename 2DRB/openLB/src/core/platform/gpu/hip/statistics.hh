/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022, 2025 Adrian Kummerlaender, Leonardo Dorneles
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

#ifndef GPU_HIP_STATISTICS_HH
#define GPU_HIP_STATISTICS_HH

#include "statistics.h"

#include "utilities/performanceProfiler.h"

#include <thrust/reduce.h>
#include <thrust/transform_reduce.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/tuple.h>
#include <thrust/execution_policy.h>

#include "core/fieldStatisticO.h"

namespace olb {

namespace gpu {

namespace hip {

/// Plain pair type with single-value constructor for use in gpu::cuda::maximum_and_plus
template <typename T, typename U>
struct pair {
  T first;
  U second;

  pair() any_platform:
    first{},
    second{}
  { }

  pair(T init) any_platform:
    first{init},
    second{init}
  { }

  pair(T a, U b) any_platform:
    first{a},
    second{b}
  { }

};

/// Function object for simulateneously computing maximum and sum in a single thrust::reduce
template <typename T>
struct maximum_and_plus {
  using first_argument_type  = pair<T,T>;
  using second_argument_type = pair<T,T>;
  using result_type          = pair<T,T>;

  __device__ result_type operator()(const first_argument_type&  lhs,
                                    const second_argument_type& rhs) const {
    return {lhs.first < rhs.first ? rhs.first : lhs.first,
            lhs.second + rhs.second};
  }
};

/// Masked accumulator for the fused rho/u statistics reduction
template <typename T>
struct StatAccum {
  std::size_t count;
  T rhoSum, rhoMin, rhoMax;
  T uSqrSum, uSqrMin, uSum, uSqrMax;
};

/// Maps a (generated, rho, uSqr) cell tuple to a masked StatAccum
template <typename T>
struct statTransform {
  template <typename Tuple>
  __device__ StatAccum<T> operator()(const Tuple& t) const {
    if (thrust::get<0>(t) != 0) {
      const T rho  = thrust::get<1>(t);
      const T uSqr = thrust::get<2>(t);
      return StatAccum<T>{1, rho, rho, rho, uSqr, uSqr, util::sqrt(uSqr), uSqr};
    }
    return StatAccum<T>{0, T{}, std::numeric_limits<T>::max(), std::numeric_limits<T>::lowest(),
                        T{}, std::numeric_limits<T>::max(), T{}, std::numeric_limits<T>::lowest()};
  }
};

/// Combines two StatAccum (sum counts and sums, min/max the extremes)
template <typename T>
struct statCombine {
  __device__ StatAccum<T> operator()(const StatAccum<T>& a, const StatAccum<T>& b) const {
    return StatAccum<T>{
      a.count + b.count,
      a.rhoSum + b.rhoSum,
      a.rhoMin < b.rhoMin ? a.rhoMin : b.rhoMin,
      a.rhoMax > b.rhoMax ? a.rhoMax : b.rhoMax,
      a.uSqrSum + b.uSqrSum,
      a.uSqrMin < b.uSqrMin ? a.uSqrMin : b.uSqrMin,
      a.uSum + b.uSum,
      a.uSqrMax > b.uSqrMax ? a.uSqrMax : b.uSqrMax
    };
  }
};

}

}

template <typename T, typename DESCRIPTOR>
void StatisticsPostProcessor::type<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_HIP>>::apply(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_HIP>& blockLattice)
{
  OLB_PROFILE_T("GPU_HIP:Stats:Reduce", DESCRIPTOR);
  if (!blockLattice.statisticsEnabled()) {
    return;
  }

  // Single fused reduction of rho and u (min/max/mean), masked by the
  // per-cell STATISTIC_GENERATED flag so disabled cells do not skew min/max.
  auto& statisticGenerated = blockLattice.template getField<descriptors::STATISTIC_GENERATED>();
  auto& statistic          = blockLattice.template getField<descriptors::STATISTIC>();
  const std::size_t n = blockLattice.getNcells();

  auto first = thrust::make_zip_iterator(thrust::make_tuple(
    statisticGenerated[0].deviceData(), statistic[0].deviceData(), statistic[1].deviceData()));

  using Accum = gpu::hip::StatAccum<T>;
  const Accum init{0, T{}, std::numeric_limits<T>::max(), std::numeric_limits<T>::lowest(),
                   T{}, std::numeric_limits<T>::max(), T{}, std::numeric_limits<T>::lowest()};
  const Accum r = thrust::transform_reduce(thrust::device, first, first + n,
                                           gpu::hip::statTransform<T>{}, init,
                                           gpu::hip::statCombine<T>{});

  typename LatticeStatistics<T,DESCRIPTOR>::Aggregatable statistics{
    .nCells   = r.count,
    .avRho    = r.rhoSum,
    .avEnergy = r.uSqrSum,
    .maxU     = r.uSqrMax,
    .minRho   = r.rhoMin,
    .maxRho   = r.rhoMax,
    .minUSqr  = r.uSqrMin,
    .sumU     = r.uSum
  };
  blockLattice.getStatistics().incrementStats(statistics);

  blockLattice.getStatistics().reset();
}


namespace gpu {
namespace hip {

/// Accumulator for a registered field reduction (single component)
template <typename T>
struct FieldStatAccum {
  T count, min, max, sum;
};

/// Maps a cell index to its masked (core-only) FieldStatAccum
template <typename T, unsigned D>
struct fieldStatTransform {
  BlockStructureD<D> structure;
  const T*           values;
  __device__ FieldStatAccum<T> operator()(std::size_t iCell) const {
    if (structure.isInsideCore(structure.getLatticeR(iCell))) {
      const T v = values[iCell];
      return FieldStatAccum<T>{T(1), v, v, v};
    }
    return FieldStatAccum<T>{T(0), std::numeric_limits<T>::max(),
                             std::numeric_limits<T>::lowest(), T(0)};
  }
};

template <typename T>
struct fieldStatCombine {
  __device__ FieldStatAccum<T> operator()(const FieldStatAccum<T>& a, const FieldStatAccum<T>& b) const {
    return FieldStatAccum<T>{
      a.count + b.count,
      a.min < b.min ? a.min : b.min,
      a.max > b.max ? a.max : b.max,
      a.sum + b.sum
    };
  }
};

}
}

/// Device reduction of a registered FIELD over the core cells (masked by isCore)
template <typename FIELD>
template <typename T, typename DESCRIPTOR>
struct FieldStatisticO<FIELD>::type<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_HIP>> {
  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_HIP>&) { }

  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_HIP>& block) {
    if (!block.statisticsEnabled()) {
      return;
    }
    const unsigned    d = DESCRIPTOR::template size<FIELD>();
    const std::size_t n = block.getNcells();
    auto& field = block.template getField<FIELD>();
    const BlockStructureD<DESCRIPTOR::d> structure = block;

    typename LatticeStatistics<T,DESCRIPTOR>::FieldResult r;
    r.components = d;
    r.min.resize(d);
    r.max.resize(d);
    r.sum.resize(d);
    r.count = 0;

    using Accum = gpu::hip::FieldStatAccum<T>;
    const Accum init{T(0), std::numeric_limits<T>::max(), std::numeric_limits<T>::lowest(), T(0)};
    auto index = thrust::make_counting_iterator<std::size_t>(0);
    for (unsigned iD = 0; iD < d; ++iD) {
      const Accum res = thrust::transform_reduce(thrust::device, index, index + n,
        gpu::hip::fieldStatTransform<T,DESCRIPTOR::d>{structure, field[iD].deviceData()},
        init, gpu::hip::fieldStatCombine<T>{});
      r.min[iD] = res.min;
      r.max[iD] = res.max;
      r.sum[iD] = res.sum;
      if (iD == 0) {
        r.count = static_cast<std::size_t>(res.count);
      }
    }

    block.getStatistics().template setFieldResult<FIELD>(r);
  }
};


}

#endif
