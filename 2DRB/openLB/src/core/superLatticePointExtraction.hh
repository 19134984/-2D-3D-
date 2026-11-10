/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Adrian Kummerlaender
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

#ifndef CORE_SUPER_LATTICE_POINT_EXTRACTION_HH
#define CORE_SUPER_LATTICE_POINT_EXTRACTION_HH

#include "superLatticePointExtraction.h"
#include <functional>

namespace olb {

template <typename T, typename DESCRIPTOR, typename FUNCTOR>
void SuperLatticePointExtraction<T,DESCRIPTOR,FUNCTOR>::updateExtractionResultD()
{
  if (!_localResultD) {
    return;
  }

  LoadBalancer<T>& load = _sLattice.getLoadBalancer();
  std::map<int,int> rankLocalPointsIndex;
  std::size_t iLocal = 0;

  for (const auto& [iC, iExtracted, physR] : _rankLocalPoints) {
    auto& block = _sampleD.getBlock(load.loc(iC));
    auto iPoint = rankLocalPointsIndex[iC]++;
    auto point = block.get(iPoint);
    _localResultD->setField(iLocal++, point.template getField<descriptors::packed_field_t<typename FUNCTOR::result_t>>());
  }
}

template <typename T, typename DESCRIPTOR, typename FUNCTOR>
SuperLatticePointExtraction<T,DESCRIPTOR,FUNCTOR>::SuperLatticePointExtraction(
  SuperLattice<T,DESCRIPTOR>& sLattice,
  const PointGeneratorF& pointGenerator)
  : _sLattice(sLattice),
    _sampleD(sLattice.getLoadBalancer()),
    _couplingO(PointExtractionO<FUNCTOR>{},
               names::Lattice1{}, sLattice,
               names::Points{}, _sampleD)
{
  build([&](IndexedConsumerF emit) {
    std::size_t iG = 0;
    pointGenerator([&](const Vector<T,DESCRIPTOR::d>& physR) {
      emit(iG++, physR);
    });
  });
}

template <typename T, typename DESCRIPTOR, typename FUNCTOR>
SuperLatticePointExtraction<T,DESCRIPTOR,FUNCTOR>::SuperLatticePointExtraction(
  SuperLattice<T,DESCRIPTOR>& sLattice,
  const IndexedGeneratorF& pointGenerator)
  : _sLattice(sLattice),
    _sampleD(sLattice.getLoadBalancer()),
    _couplingO(PointExtractionO<FUNCTOR>{},
               names::Lattice1{}, sLattice,
               names::Points{}, _sampleD)
{
  build([&](IndexedConsumerF emit) {
    pointGenerator(emit);
  });
}

template <typename T, typename DESCRIPTOR, typename FUNCTOR>
SuperLatticePointExtraction<T,DESCRIPTOR,FUNCTOR>::SuperLatticePointExtraction(
  SuperLattice<T,DESCRIPTOR>& sLattice,
  const std::vector<Vector<T,DESCRIPTOR::d>>& points)
  : SuperLatticePointExtraction(sLattice,
      PointGeneratorF([&points](PointConsumerF emit) {
        for (const auto& p : points) {
          emit(p);
        }
      }))
{ }

template <typename T, typename DESCRIPTOR, typename FUNCTOR>
void SuperLatticePointExtraction<T,DESCRIPTOR,FUNCTOR>::build(
  const IndexedGeneratorF& generator)
{
  const auto& cDecomposition  = _sLattice.getCuboidDecomposition();
  LoadBalancer<T>& load = _sLattice.getLoadBalancer();

  for (int iC=0; iC < cDecomposition.size(); ++iC) {
    _rankLocalPointsSize[iC] = 0;
  }

  std::vector<std::pair<Vector<T,DESCRIPTOR::d>, Vector<T,DESCRIPTOR::d>>> localBounds;
  const T dx = _sLattice.getUnitConverter().getPhysDeltaX();
  const Vector<T,DESCRIPTOR::d> halfCell(T{0.5} * dx);
  for (int iC = 0; iC < load.size(); ++iC) {
    auto cuboid = cDecomposition.get(load.glob(iC));
    localBounds.push_back({cuboid.getOrigin() - halfCell,
                           cuboid.getOrigin() + cuboid.getExtent() * dx - halfCell});
  }

  std::size_t numLocal = 0;

  generator([&](std::size_t iG, const Vector<T,DESCRIPTOR::d>& physR) {
    bool mightBeLocal = false;
    for (const auto& bounds : localBounds) {
      bool inside = true;
      for (unsigned d=0; d<DESCRIPTOR::d; ++d) {
        if (physR[d] < bounds.first[d] || physR[d] > bounds.second[d]) {
          inside = false;
          break;
        }
      }
      if (inside) {
        mightBeLocal = true;
        break;
      }
    }

    if (mightBeLocal) {
      if (auto iC = cDecomposition.getC(physR)) {
        if (load.isLocal(*iC)) {
          _rankLocalPoints.emplace_back(*iC, iG, physR);
          _rankLocalPointsSize[*iC]++;
          numLocal++;
        }
      }
    }
  });

  _globalSize = numLocal;
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(_globalSize, MPI_SUM);
#endif

  if (numLocal > 0) {
    _globalIndices.resize(numLocal);
    _localPointD = std::make_unique<FieldArrayD<T,DESCRIPTOR,Platform::CPU_SISD,fields::PHYS_R>>(numLocal);
    _localResultD = std::make_unique<FieldArrayD<T,DESCRIPTOR,Platform::CPU_SISD,descriptors::packed_field_t<typename FUNCTOR::result_t>>>(numLocal);
  }

  std::map<int,int> rankLocalPointsIndex;
  for (auto [iC, size] : _rankLocalPointsSize) {
    if (load.isLocal(iC)) {
      Vector<int,DESCRIPTOR::d> blockSize(1);
      blockSize[0] = size;
      _sampleD.getBlock(load.loc(iC)).resize(blockSize);
      rankLocalPointsIndex[iC] = 0;
    }
  }

  std::size_t iLocal = 0;
  for (auto& [iC, globalIdx, physR] : _rankLocalPoints) {
    auto& block = _sampleD.getBlock(load.loc(iC));
    auto iPoint = rankLocalPointsIndex[iC]++;
    auto point = block.get(iPoint);
    point.template setField<fields::PHYS_R>(physR);

    _globalIndices[iLocal] = globalIdx;
    _localPointD->setField(iLocal, physR);
    iLocal++;
  }

  _sampleD.setProcessingContext(ProcessingContext::Simulation);
}

template <typename T, typename DESCRIPTOR, typename FUNCTOR>
void SuperLatticePointExtraction<T,DESCRIPTOR,FUNCTOR>::update()
{
  _couplingO.apply();
  _sampleD.setProcessingContext(ProcessingContext::Evaluation);
  updateExtractionResultD();
}

}

#endif
