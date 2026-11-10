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

#ifndef FSI_SUPER_POROUS_ELEMENT_REDUCTION_O_H
#define FSI_SUPER_POROUS_ELEMENT_REDUCTION_O_H

#include "fields.h"
#include "operators.h"
#include "integral.h"

#include "core/superD.h"

namespace olb {

namespace stage {

namespace fsi {

struct CollectFromFluid { };

}

}

template <typename T, typename DESCRIPTOR, typename... FIELDS>
class SuperPorousElementReductionO final {
private:
  /// Fluid lattice
  SuperLattice<T,DESCRIPTOR>& _sLattice;
  /// FSI region indicator
  FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>> _indicatorF;
  /// Reduced element data obtained from the fluid lattice
  std::unique_ptr<SuperD<T,descriptors::fsi::REDUCED_ELEMENTS<DESCRIPTOR::d>>> _superReducedElementsD;

  /// Reduced per-element fields
  std::vector<StaticParametersD<T,DESCRIPTOR,FIELDS...>> _fields;
  std::size_t _nElements;

  bool _rankDoesFSI;

  #ifdef PARALLEL_MODE_MPI
  MPI_Comm _mpiCommunicator;
  std::vector<T> _localReductionBuffer;
  std::vector<T> _globalReductionBuffer;
  static constexpr unsigned _totalFieldSize = (DESCRIPTOR::template size<FIELDS>() + ...);
  #endif

public:
  SuperPorousElementReductionO(SuperLattice<T,DESCRIPTOR>& sLattice,
                               FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicatorF)
    : _sLattice(sLattice)
    , _indicatorF{std::move(indicatorF)}
    , _superReducedElementsD(new SuperD<
        T,descriptors::fsi::REDUCED_ELEMENTS<DESCRIPTOR::d>>(sLattice.getLoadBalancer()))
    , _nElements{0}
    , _rankDoesFSI{singleton::mpi().isMainProcessor()}
  {
    OstreamManager clout(std::cout, "SuperPorousElementReductionO");
    auto& load = _sLattice.getLoadBalancer();

    for (int iC=0; iC < _sLattice.getLoadBalancer().size(); ++iC) {
      auto& block = _sLattice.getBlock(iC);
      auto& blockF = _indicatorF->getBlockIndicatorF(iC);
      if (!blockF.isEmpty()) {
        block.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> latticeR) {
          if (blockF(latticeR)) {
            if (block.isPadding(latticeR, 2) || block.isInsideCore(latticeR)) {
              block.addPostProcessor(typeid(stage::PostStream),
                                     latticeR,
                                     meta::id<GrowPaddingLayerO>{});
            }
          }
        });
      }
    }
    sLattice.template addPostProcessor<stage::fsi::CollectFromFluid>(
      std::forward<decltype(indicatorF)>(indicatorF),
      meta::id<IntegratePorousElementFieldsO<FIELDS...>>{});

    #ifdef PARALLEL_MODE_MPI
    for (int iC=0; iC < load.size(); ++iC) {
      _rankDoesFSI |= sLattice.getBlock(iC).hasPostProcessor(
        typeid(stage::fsi::CollectFromFluid),
        meta::id<IntegratePorousElementFieldsO<FIELDS...>>{});
    }
    MPI_Comm_split(MPI_COMM_WORLD, _rankDoesFSI ? 0 : MPI_UNDEFINED, singleton::mpi().getRank(), &_mpiCommunicator);
    #endif

    _superReducedElementsD->setProcessingContext(ProcessingContext::Simulation);

    for (int iC=0; iC < load.size(); ++iC) {
      auto& block = _sLattice.getBlock(iC);
      auto& elementsBlock = _superReducedElementsD->getBlock(iC);
      block.template setParameter<fields::array_of<fields::fsi::REDUCED_ELEMENT_TAG>>(
        elementsBlock.template getField<fields::fsi::ELEMENT_TAG>());
      (block.template setParameter<fields::array_of<FIELDS>>(
        elementsBlock.template getField<FIELDS>()), ...);
    }
  }

  void resize(std::size_t nElements) {
    _nElements = nElements;
    _fields.resize(nElements);

    #ifdef PARALLEL_MODE_MPI
    _localReductionBuffer.resize(nElements * _totalFieldSize, T{0});
    _globalReductionBuffer.resize(nElements * _totalFieldSize, T{0});
    #endif

    for (int iC=0; iC < _sLattice.getLoadBalancer().size(); ++iC) {
      auto& elementsBlock = _superReducedElementsD->getBlock(iC);
      elementsBlock.template getField<fields::fsi::ELEMENT_TAG>().resize(nElements);
      (elementsBlock.template getField<FIELDS>().resize(nElements), ...);
    }

    // Trigger device memory allocation (deviceData pointers potentially change here)
    _superReducedElementsD->setProcessingContext(ProcessingContext::Simulation);

    for (int iC=0; iC < _sLattice.getLoadBalancer().size(); ++iC) {
      auto& elementsBlock = _superReducedElementsD->getBlock(iC);
      auto& block = _sLattice.getBlock(iC);

      block.template getData<OperatorParameters<IntegratePorousElementFieldsO<FIELDS...>>>()
           .template set<fields::fsi::REDUCED_ELEMENTS_COUNT>(nElements);

      block.template setParameter<fields::array_of<fields::fsi::REDUCED_ELEMENT_TAG>>(
        elementsBlock.template getField<fields::fsi::ELEMENT_TAG>());
      (block.template setParameter<fields::array_of<FIELDS>>(
        elementsBlock.template getField<FIELDS>()), ...);
    }
  }

  void addCollectionO(PostProcessorPromise<T,DESCRIPTOR>&& collectionO) {
    for (int iC=0; iC < _sLattice.getLoadBalancer().size(); ++iC) {
      auto& block = _sLattice.getBlock(iC);
      auto& blockF = _indicatorF->getBlockIndicatorF(iC);
      if (!blockF.isEmpty()) {
        block.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> latticeR) {
          if (blockF(latticeR)) {
            if (block.isInsideCore(latticeR)) {
              block.addPostProcessor(typeid(stage::fsi::CollectFromFluid),
                                     latticeR,
                                     std::forward<decltype(collectionO)>(collectionO));
            }
          }
        });
      }
    }

    _sLattice.template addPostProcessor<stage::Evaluation>(
      std::forward<decltype(_indicatorF)>(_indicatorF),
      std::forward<decltype(collectionO)>(collectionO));
  }

  bool rankDoesFSI() const
  {
    return _rankDoesFSI;
  }

  void apply()
  {
    _sLattice.executePostProcessors(stage::fsi::CollectFromFluid{});

    if (!_rankDoesFSI) {
      return;
    }

    for (int iC=0; iC < _sLattice.getLoadBalancer().size(); ++iC) {
      auto& block = _superReducedElementsD->getBlock(iC);
      block.template getField<fields::fsi::ELEMENT_TAG>()
           .setProcessingContext(ProcessingContext::Evaluation);
      ((block.template getField<FIELDS>()
             .setProcessingContext(ProcessingContext::Evaluation))
       , ...);
    }

    #ifndef PARALLEL_MODE_MPI
    for (auto& f : _fields) {
      f = StaticParametersD<T,DESCRIPTOR,FIELDS...>{};
    }

    for (int iC=0; iC < _sLattice.getLoadBalancer().size(); ++iC) {
      auto& block = _superReducedElementsD->getBlock(iC);
      auto accumulate = [&]<typename FIELD>(meta::id<FIELD>, std::size_t iElement) {
        auto& fieldArray = block.template getField<FIELD>();
        FieldD<T,DESCRIPTOR,FIELD> f([&](std::size_t iD) {
          return fieldArray[iD][iElement];
        });
        _fields[iElement].template set<FIELD>(_fields[iElement].template get<FIELD>() + f);
      };
      for (std::size_t iElement=0; iElement < _nElements; ++iElement) {
        ((accumulate(meta::id<FIELDS>{}, iElement)), ...);
      }
    }
    #else
    std::fill(_localReductionBuffer.begin(), _localReductionBuffer.end(), T{0});

    for (int iC=0; iC < _sLattice.getLoadBalancer().size(); ++iC) {
      auto& block = _superReducedElementsD->getBlock(iC);
      unsigned fieldOffset = 0;

      auto accumulateToBuffer = [&]<typename FIELD>(meta::id<FIELD>) {
        const unsigned fSize = DESCRIPTOR::template size<FIELD>();
        auto& fieldArray = block.template getField<FIELD>();
        for (std::size_t iElement=0; iElement < _nElements; ++iElement) {
          for (std::size_t iD=0; iD < fSize; ++iD) {
            _localReductionBuffer[iElement * _totalFieldSize + fieldOffset + iD] += fieldArray[iD][iElement];
          }
        }
        fieldOffset += fSize;
      };

      (accumulateToBuffer(meta::id<FIELDS>{}), ...);
    }

    singleton::mpi().allreduce(_localReductionBuffer.data(),
                               _globalReductionBuffer.data(),
                               _globalReductionBuffer.size(),
                               MPI_SUM,
                               _mpiCommunicator);

    unsigned unpackOffset = 0;
    auto unpackFromBuffer = [&]<typename FIELD>(meta::id<FIELD>) {
      const unsigned fSize = DESCRIPTOR::template size<FIELD>();
      for (std::size_t iElement=0; iElement < _nElements; ++iElement) {
        FieldD<T,DESCRIPTOR,FIELD> f([&](std::size_t iD) {
          return _globalReductionBuffer[iElement * _totalFieldSize + unpackOffset + iD];
        });
        _fields[iElement].template set<FIELD>(f);
      }
      unpackOffset += fSize;
    };

    (unpackFromBuffer(meta::id<FIELDS>{}), ...);
    #endif
  }

  std::size_t getElementCount() const {
    return _nElements;
  }

  template <concepts::Field FIELD>
  auto getField(unsigned iElement) const {
    return _fields[iElement - 1].template get<FIELD>();
  }

};

}

#endif
