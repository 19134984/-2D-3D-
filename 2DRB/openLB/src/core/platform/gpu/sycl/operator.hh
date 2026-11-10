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

#ifndef GPU_SYCL_OPERATOR_HH
#define GPU_SYCL_OPERATOR_HH

#include "operator.h"

#include "context.hh"
#include "dynamics.hh"

#include "core/introspection.h"

#include "solver/names.h"

namespace olb {

/// Mask of non-overlap block subdomain
struct CollisionSubdomainMask;

namespace gpu {

namespace sycl {

/// Masked application of DYNAMICS::collide
template <typename T, typename DESCRIPTOR, typename DYNAMICS>
class MaskedCollision {
private:
  ParametersOfOperatorD<T,DESCRIPTOR,DYNAMICS> _parameters;
  bool* _mask;

  CellStatistic<T> apply(DeviceContext<T,DESCRIPTOR> lattice, CellID iCell) const {
    DataOnlyCell<T,DESCRIPTOR> cell(lattice, iCell);
    return DYNAMICS().collide(cell, _parameters);
  }

public:
  MaskedCollision(ParametersOfOperatorD<T,DESCRIPTOR,DYNAMICS>& parameters, bool* mask) :
    _parameters{parameters},
    _mask{mask}
  { }

  bool operator()(DeviceContext<T,DESCRIPTOR> lattice, CellID iCell) const {
    if (_mask[iCell]) {
      apply(lattice, iCell);
      return true;
    }
    return false;
  }

  bool operator()(DeviceContext<T,DESCRIPTOR> lattice, CellID iCell, CellStatistic<T>& statistic) const {
    if (_mask[iCell]) {
      statistic = apply(lattice, iCell);
      return true;
    }
    return false;
  }

};

/// List-based application of DYNAMICS::collide
template <typename T, typename DESCRIPTOR, typename DYNAMICS>
class ListedCollision {
private:
  ParametersOfOperatorD<T,DESCRIPTOR,DYNAMICS> _parameters;

public:
  ListedCollision(ParametersOfOperatorD<T,DESCRIPTOR,DYNAMICS>& parameters) :
    _parameters{parameters}
  { }

  bool operator()(DeviceContext<T,DESCRIPTOR> lattice, CellID iCell) const {
    DataOnlyCell<T,DESCRIPTOR> cell(lattice, iCell);
    DYNAMICS().collide(cell, _parameters);
    return true;
  }

  bool operator()(DeviceContext<T,DESCRIPTOR> lattice, CellID iCell, CellStatistic<T>& statistic) const {
    DataOnlyCell<T,DESCRIPTOR> cell(lattice, iCell);
    statistic = DYNAMICS().collide(cell, _parameters);
    return true;
  }

};


/// Masked application of OPERATOR::apply
template <typename OPERATOR>
class MaskedPostProcessor {
private:
  bool* _mask;

public:
  MaskedPostProcessor(bool* mask) :
    _mask{mask}
  { }

  template <typename T, typename DESCRIPTOR>
  bool operator()(DeviceBlockLattice<T,DESCRIPTOR> lattice, CellID iCell) const {
    if (_mask[iCell]) {
      Cell<T,DESCRIPTOR> cell(lattice, iCell);
      OPERATOR().apply(cell);
      return true;
    }
    return false;
  }

};

/// List-based application of OPERATOR::apply
template <typename OPERATOR>
struct ListedPostProcessor {
  template <typename T, typename DESCRIPTOR>
  bool operator()(DeviceBlockLattice<T,DESCRIPTOR> lattice, CellID iCell) const {
    Cell<T,DESCRIPTOR> cell(lattice, iCell);
    OPERATOR().apply(cell);
    return true;
  }

};

/// List-based application of OPERATOR::apply with parameters
template <typename T, typename DESCRIPTOR, typename OPERATOR>
class ListedPostProcessorWithParameters {
private:
  ParametersOfOperatorD<T,DESCRIPTOR,OPERATOR> _parameters;

public:
  ListedPostProcessorWithParameters(ParametersOfOperatorD<T,DESCRIPTOR,OPERATOR>& parameters) :
    _parameters{parameters}
  { }

  bool operator()(DeviceBlockLattice<T,DESCRIPTOR> lattice, CellID iCell) const {
    Cell<T,DESCRIPTOR> cell(lattice, iCell);
    OPERATOR().apply(cell, _parameters);
    return true;
  }

};

/// Unrestricted application of COUPLING::apply
template <typename COUPLER>
struct UnmaskedCoupling {
  template <typename CONTEXT>
  bool operator()(CONTEXT lattices, CellID iCell) const {
    auto cells = lattices.exchange_values([&](auto name) -> auto {
      return Cell{lattices.get(name), iCell};
    });
    COUPLER().apply(cells);
    return true;
  }
};

/// Unrestricted application of COUPLING::apply with parameters
template <typename COUPLER, typename COUPLEES>
class UnmaskedCouplingWithParameters {
private:
  using PresentParametersD = typename COUPLER::parameters::template decompose_into<
    AbstractCouplingO<COUPLEES>::ParametersD::template include_fields
  >;
  PresentParametersD _parameters;

public:
  UnmaskedCouplingWithParameters(PresentParametersD& parameters) :
    _parameters{parameters}
  { }

  template <typename CONTEXT>
  bool operator()(CONTEXT lattices, CellID iCell) const {
    auto cells = lattices.exchange_values([&](auto name) -> auto {
      return Cell{lattices.get(name), iCell};
    });
    COUPLER().apply(cells, _parameters);
    return true;
  }

};

template <typename COUPLER, typename COUPLEES>
class ParticleCouplingWithParameters {
private:
  typename AbstractCouplingO<COUPLEES>::ParametersD::template include<
    typename meta::merge<
      typename COUPLER::parameters,
      meta::list<
        fields::converter::PHYS_DELTA_X,
        fields::BLOCK_LOWER
      >
    >
  > _parameters;

public:
  ParticleCouplingWithParameters(
    const typename AbstractCouplingO<COUPLEES>::ParametersD::template include<
      typename meta::merge<
        typename COUPLER::parameters,
        meta::list<
          fields::converter::PHYS_DELTA_X,
          fields::BLOCK_LOWER
        >
      >
    >& parameters):
    _parameters{parameters}
  { }

  ParticleCouplingWithParameters(const ParticleCouplingWithParameters&) = default;

  template <typename CONTEXT>
  bool operator()(CONTEXT lattices, CellID iParticle) const {
    using V = typename CONTEXT::template value_t<names::Points>::value_t;
    using DESCRIPTOR = typename CONTEXT::template value_t<names::Points>::descriptor_t;
    const auto originR = _parameters.template get<fields::BLOCK_LOWER>();
    const auto deltaX = _parameters.template get<fields::converter::PHYS_DELTA_X>();
    const auto physRs = lattices.get(meta::id<names::Points>{})
                                .template getField<fields::PHYS_R>();
    const PhysR<V,DESCRIPTOR::d> physR([&](unsigned iD) {
      return physRs[iD][iParticle];
    });
    const LatticeR<DESCRIPTOR::d> latticeR = util::floor((physR - originR) / deltaX + V{0.5});
    auto cells = lattices.exchange_values([&](auto name) -> auto {
      if constexpr (std::is_same_v<typename decltype(name)::type, names::Points>) {
        return DataOnlyCell{lattices.get(name), iParticle};
      } else {
        const auto iCell = lattices.get(name).getCellId(latticeR);
        return Cell{lattices.get(name), iCell};
      }
    });
    COUPLER().apply(cells, _parameters);
    return true;
  }

};


template <typename CONTEXT, typename... ARGS>
void call_operators(CONTEXT& lattice, bool* subdomain, ARGS&&... args) {
  OLB_PROFILE_T("GPU_SYCL:Knl:Sync", ARGS..., CONTEXT);
  auto& queue = device::getQueue();
  std::size_t nCells = lattice.getNcells();
  queue.submit([&](::sycl::handler& cgh) {
    cgh.parallel_for(::sycl::range<1>(nCells), [=](::sycl::id<1> idx) {
      CellID iCell = idx[0];
      if (!subdomain[iCell]) return;
      (args(lattice, iCell) || ... );
    });
  }).wait_and_throw();
}

template <typename CONTEXT, typename... ARGS>
void async_call_operators(::sycl::queue* stream, CONTEXT& lattice, bool* subdomain, ARGS&&... args) {
  OLB_PROFILE_T("GPU_SYCL:Knl:Submit", ARGS..., CONTEXT);
  std::size_t nCells = lattice.getNcells();
  stream->submit([&](::sycl::handler& cgh)
    [[intel::kernel_args_restrict]]
  {
    cgh.parallel_for(::sycl::range<1>(nCells), [=](::sycl::id<1> idx) {
      CellID iCell = idx[0];
      if (!subdomain[iCell]) return;
      (args(lattice, iCell) || ... );
    });
  });
}

template <typename CONTEXT, typename... ARGS>
void call_operators_with_statistics(CONTEXT& lattice, bool* subdomain, ARGS&&... args) {
  OLB_PROFILE_T("GPU_SYCL:Knl:Stats:Sync", ARGS..., CONTEXT);
  auto& queue = device::getQueue();
  std::size_t nCells = lattice.getNcells();
  queue.submit([&](::sycl::handler& cgh) {
    cgh.parallel_for(::sycl::range<1>(nCells), [=](::sycl::id<1> idx) {
      CellID iCell = idx[0];
      if (!subdomain[iCell]) return;

      auto statistic = lattice.template getField<descriptors::STATISTIC>();
      auto statisticGenerated = lattice.template getField<descriptors::STATISTIC_GENERATED>()[0];
      CellStatistic<typename CONTEXT::value_t> cellStatistic{-1, -1};

      if ((args(lattice, iCell, cellStatistic) || ... )) {
        if (cellStatistic) {
          statisticGenerated[iCell] = 1;
          statistic[0][iCell] = cellStatistic.rho;
          statistic[1][iCell] = cellStatistic.uSqr;
        } else {
          statisticGenerated[iCell] = 0;
          statistic[0][iCell] = 0;
          statistic[1][iCell] = 0;
        }
      }
    });
  }).wait_and_throw();
}

template <typename CONTEXT, typename... ARGS>
::sycl::event async_call_operators_with_statistics(::sycl::queue* stream,
                                                   CONTEXT& lattice,
                                                   bool* subdomain,
                                                   ARGS&&... args) {
  OLB_PROFILE_T("GPU_SYCL:Knl:Stats:Submit", ARGS..., CONTEXT);
  std::size_t nCells = lattice.getNcells();
  return stream->submit([&](::sycl::handler& cgh) {
    cgh.parallel_for(::sycl::range<1>(nCells), [=](::sycl::id<1> idx)
      [[intel::kernel_args_restrict]]
    {
      CellID iCell = idx[0];
      if (!subdomain[iCell]) return;

      auto statistic = lattice.template getField<descriptors::STATISTIC>();
      auto statisticGenerated = lattice.template getField<descriptors::STATISTIC_GENERATED>()[0];
      CellStatistic<typename CONTEXT::value_t> cellStatistic{-1, -1};

      if ((args(lattice, iCell, cellStatistic) || ... )) {
        if (cellStatistic) {
          statisticGenerated[iCell] = 1;
          statistic[0][iCell] = cellStatistic.rho;
          statistic[1][iCell] = cellStatistic.uSqr;
        } else {
          statisticGenerated[iCell] = 0;
          statistic[0][iCell] = 0;
          statistic[1][iCell] = 0;
        }
      }
    });
  });
}

template <typename CONTEXT, typename... ARGS>
void call_list_operators(CONTEXT& lattice,
                         const gpu::sycl::Column<CellID>& cells,
                         ARGS&&... args) {
  OLB_PROFILE_T("GPU_SYCL:KnlList:Sync", ARGS..., CONTEXT);
  auto& queue = device::getQueue();
  std::size_t nIndices = cells.size();
  const CellID* indices = cells.deviceData();

  queue.submit([&](::sycl::handler& cgh) {
    cgh.parallel_for(::sycl::range<1>(nIndices), [=](::sycl::id<1> idx) {
      CellID iCell = indices[idx[0]];
      (args(lattice, iCell) || ... );
    });
  }).wait_and_throw();
}

template <typename CONTEXT, typename... ARGS>
::sycl::event async_call_list_collision_operators(::sycl::queue* stream,
                                                  const std::vector<::sycl::event>& dependencies,
                                                  CONTEXT& lattice,
                                                  const gpu::sycl::Column<CellID>& cells,
                                                  ARGS&&... args) {
  OLB_PROFILE_T("GPU_SYCL:KnlList:Coll:Submit", ARGS..., CONTEXT);
  std::size_t nIndices = cells.size();
  const CellID* indices = cells.deviceData();

#ifdef OLB_SYCL_EXPENSIVE_OPERATORS
  std::size_t localSize = 256;
  std::size_t globalSize = ((nIndices + localSize - 1) / localSize) * localSize;

  return stream->submit([&](::sycl::handler& cgh) {
    cgh.depends_on(dependencies);
    cgh.parallel_for(::sycl::nd_range<1>(globalSize, localSize), [=](::sycl::nd_item<1> item)
      [[intel::kernel_args_restrict]]
      [[sycl::reqd_sub_group_size(16)]]
    {
      std::size_t idx = item.get_global_id(0);
      if (idx < nIndices) {
        CellID iCell = indices[idx];
        (args(lattice, iCell) || ... );
      }
    });
  });
#else
  return stream->submit([&](::sycl::handler& cgh) {
    cgh.depends_on(dependencies);
    cgh.parallel_for(::sycl::range<1>(nIndices), [=](::sycl::id<1> idx)
      [[intel::kernel_args_restrict]]
    {
      CellID iCell = indices[idx[0]];
      (args(lattice, iCell) || ... );
    });
  });
#endif
}

template <typename CONTEXT, typename... ARGS>
::sycl::event async_call_list_operators(::sycl::queue* stream,
                                        const std::vector<::sycl::event>& dependencies,
                                        CONTEXT& lattice,
                                        const gpu::sycl::Column<CellID>& cells,
                                        ARGS&&... args) {
  OLB_PROFILE_T("GPU_SYCL:KnlList:Submit", ARGS..., CONTEXT);
  std::size_t nIndices = cells.size();
  const CellID* indices = cells.deviceData();

  return stream->submit([&](::sycl::handler& cgh) {
    cgh.depends_on(dependencies);
    cgh.parallel_for(::sycl::range<1>(nIndices), [=](::sycl::id<1> idx)
      [[intel::kernel_args_restrict]]
    {
      CellID iCell = indices[idx[0]];
      (args(lattice, iCell) || ... );
    });
  });
}

template <typename CONTEXT, typename... ARGS>
::sycl::event async_call_list_operators_with_statistics(::sycl::queue* stream,
                                                        const std::vector<::sycl::event>& dependencies,
                                                        CONTEXT& lattice,
                                                        const gpu::sycl::Column<CellID>& cells,
                                                        ARGS&&... args) {
  OLB_PROFILE_T("GPU_SYCL:KnlList:Stats:Submit", ARGS..., CONTEXT);
  std::size_t nIndices = cells.size();
  const CellID* indices = cells.deviceData();

#ifdef OLB_SYCL_EXPENSIVE_OPERATORS
  std::size_t localSize = 256;
  std::size_t globalSize = ((nIndices + localSize - 1) / localSize) * localSize;

  return stream->submit([&](::sycl::handler& cgh) {
    cgh.depends_on(dependencies);
    cgh.parallel_for(::sycl::nd_range<1>(globalSize, localSize), [=](::sycl::nd_item<1> item)
      [[intel::kernel_args_restrict]]
      [[sycl::reqd_sub_group_size(16)]]
    {
      std::size_t idx = item.get_global_id(0);
      if (idx < nIndices) {
        CellID iCell = indices[idx];
        auto statistic = lattice.template getField<descriptors::STATISTIC>();
        auto statisticGenerated = lattice.template getField<descriptors::STATISTIC_GENERATED>()[0];
        CellStatistic<typename CONTEXT::value_t> cellStatistic{-1, -1};

        if ((args(lattice, iCell, cellStatistic) || ... )) {
          if (cellStatistic) {
            statisticGenerated[iCell] = 1;
            statistic[0][iCell] = cellStatistic.rho;
            statistic[1][iCell] = cellStatistic.uSqr;
          } else {
            statisticGenerated[iCell] = 0;
            statistic[0][iCell] = 0;
            statistic[1][iCell] = 0;
          }
        }
      }
    });
  });
#else
  return stream->submit([&](::sycl::handler& cgh) {
    cgh.depends_on(dependencies);
    cgh.parallel_for(::sycl::range<1>(nIndices), [=](::sycl::id<1> idx)
      [[intel::kernel_args_restrict]]
    {
      CellID iCell = indices[idx[0]];
      auto statistic = lattice.template getField<descriptors::STATISTIC>();
      auto statisticGenerated = lattice.template getField<descriptors::STATISTIC_GENERATED>()[0];
      CellStatistic<typename CONTEXT::value_t> cellStatistic{-1, -1};

      if ((args(lattice, iCell, cellStatistic) || ... )) {
        if (cellStatistic) {
          statisticGenerated[iCell] = 1;
          statistic[0][iCell] = cellStatistic.rho;
          statistic[1][iCell] = cellStatistic.uSqr;
        } else {
          statisticGenerated[iCell] = 0;
          statistic[0][iCell] = 0;
          statistic[1][iCell] = 0;
        }
      }
    });
  });
#endif
}

template <typename CONTEXT, typename... ARGS>
void call_coupling_operators(CONTEXT& lattices, bool* subdomain, ARGS&&... args) {
  OLB_PROFILE_T("GPU_SYCL:Coupler:Sync", ARGS..., CONTEXT);
  auto& queue = device::getQueue();
  std::size_t nCells = lattices.template get<0>().getNcells();
  queue.submit([&](::sycl::handler& cgh) {
    cgh.parallel_for(::sycl::range<1>(nCells), [=](::sycl::id<1> idx) {
      CellID iCell = idx[0];
      if (!subdomain[iCell]) return;
      (args(lattices, iCell) || ... );
    });
  }).wait_and_throw();
}

template <typename CONTEXT, typename... ARGS>
void call_particle_coupling_operators(CONTEXT& lattices, ARGS&&... args) {
  OLB_PROFILE_T("GPU_SYCL:Coupler:Particle:Sync", ARGS..., CONTEXT);
  auto& queue = device::getQueue();
  std::size_t nParticles = lattices.get(meta::id<names::Points>{}).getNcells();
  queue.submit([&](::sycl::handler& cgh) {
    cgh.parallel_for(::sycl::range<1>(nParticles), [=](::sycl::id<1> idx) {
      CellID iParticle = idx[0];
      (args(lattices, iParticle) || ... );
    });
  }).wait_and_throw();
}

template <typename T, typename DESCRIPTOR, typename... DYNAMICS>
std::function<void(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>&)>
getFusedCollisionO() {
  return [](ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& block) {
    bool* subdomain = block.template getData<CollisionSubdomainMask>().deviceData();
    DeviceContext<T,DESCRIPTOR> lattice(block);
    if (block.statisticsEnabled()) {
      call_operators_with_statistics(
        lattice,
        subdomain,
        MaskedCollision<T,DESCRIPTOR,DYNAMICS>{
          block.template getData<OperatorParameters<DYNAMICS>>().parameters,
          block.template getData<DynamicsMask<DYNAMICS>>().deviceData()
        }...,
        DynamicDispatchCollision{});
    } else {
      call_operators(
        lattice,
        subdomain,
        MaskedCollision<T,DESCRIPTOR,DYNAMICS>{
          block.template getData<OperatorParameters<DYNAMICS>>().parameters,
          block.template getData<DynamicsMask<DYNAMICS>>().deviceData()
        }...,
        DynamicDispatchCollision{});
    }
  };
}

namespace kernel {
template <typename T, typename DESCRIPTOR, typename DYNAMICS, typename PARAMETERS=typename DYNAMICS::ParametersD>
void construct_dynamics(void* target, PARAMETERS* parameters)  {
  new (target) ConcreteDynamics<T,DESCRIPTOR,DYNAMICS>(parameters);
}
}

}

}

template <typename T, typename DESCRIPTOR, typename DYNAMICS>
ConcreteBlockCollisionO<T,DESCRIPTOR,Platform::GPU_SYCL,DYNAMICS>::ConcreteBlockCollisionO():
  _dynamics(new DYNAMICS()),
  _parameters(nullptr),
  _mask(nullptr),
  _outerCells(0),
  _innerCells(0),
  _modified(true)
{ }

template <typename T, typename DESCRIPTOR, typename DYNAMICS>
void ConcreteBlockCollisionO<T,DESCRIPTOR,Platform::GPU_SYCL,DYNAMICS>::applyDominant(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& block,
  ConcreteBlockMask<T,Platform::GPU_SYCL>&               subdomain)
{
  throw std::runtime_error("Dominant strategy not supported for SYCL");
}

template <typename T, typename DESCRIPTOR, typename DYNAMICS>
void ConcreteBlockCollisionO<T,DESCRIPTOR,Platform::GPU_SYCL,DYNAMICS>::applyIndividual(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& block,
  ConcreteBlockMask<T,Platform::GPU_SYCL>&               subdomain)
{
  OLB_PROFILE_T("GPU_SYCL:Coll:Apply", DYNAMICS, DESCRIPTOR);
  using namespace gpu::sycl;

  DeviceContext<T,DESCRIPTOR> lattice(block);

  if (_modified) {
    _outerCells.clear();
    _innerCells.clear();
    for (CellID iCell=0; iCell  < block.getNcells(); ++iCell) {
      if (_mask->operator[](iCell)) {
        if (block.isPaddingCore(iCell, 1)) {
          _outerCells.push_back(iCell);
        } else {
          _innerCells.push_back(iCell);
        }
      }
    }

    _outerCells.setProcessingContext(ProcessingContext::Simulation);
    _innerCells.setProcessingContext(ProcessingContext::Simulation);
    _modified = false;
  }

  std::vector<::sycl::event> outerDependencies;

  if (block.statisticsEnabled()) {
    if (_outerCells.size() > 0) {
      OLB_PROFILE_T("GPU_SYCL:Coll:OuterStats", DYNAMICS, DESCRIPTOR);
      outerDependencies.emplace_back(async_call_list_operators_with_statistics(
        &device::getQueue<queues::OuterCollision>(),
        {},
        lattice,
        _outerCells,
        ListedCollision<T,DESCRIPTOR,DYNAMICS>{_parameters->parameters}));
    }
    if (_innerCells.size() > 0) {
      OLB_PROFILE_T("GPU_SYCL:Coll:InnerStats", DYNAMICS, DESCRIPTOR);
      async_call_list_operators_with_statistics(
        &device::getQueue<queues::InnerCollision>(),
        outerDependencies,
        lattice,
        _innerCells,
        ListedCollision<T,DESCRIPTOR,DYNAMICS>{_parameters->parameters});
    }
  } else {
    if (_outerCells.size() > 0) {
      OLB_PROFILE_T("GPU_SYCL:Coll:Outer", DYNAMICS, DESCRIPTOR);
      outerDependencies.emplace_back(async_call_list_collision_operators(
        &device::getQueue<queues::OuterCollision>(),
        {},
        lattice,
        _outerCells,
        ListedCollision<T,DESCRIPTOR,DYNAMICS>{_parameters->parameters}));
    }
    if (_innerCells.size() > 0) {
      OLB_PROFILE_T("GPU_SYCL:Coll:Inner", DYNAMICS, DESCRIPTOR);
      async_call_list_collision_operators(
        &device::getQueue<queues::InnerCollision>(),
        outerDependencies,
        lattice,
        _innerCells,
        ListedCollision<T,DESCRIPTOR,DYNAMICS>{_parameters->parameters});
    }
  }

  // Record outer collision completion for overlap-aware halo communication to
  // order against (the transfer barriers on this per-block event, not the
  // shared OuterCollision queue). Recorded on the OuterCollision queue, so it
  // captures the boundary collision without waiting on the inner bulk.
  block.outerCollisionDone().record(queues::OuterCollision{});
}

template <typename T, typename DESCRIPTOR, typename DYNAMICS>
void ConcreteBlockCollisionO<T,DESCRIPTOR,Platform::GPU_SYCL,DYNAMICS>::setup(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& block)
{
  // Fetch pointers to DYNAMICS-specific parameter and mask data
  _parameters = &block.template getData<OperatorParameters<DYNAMICS>>();
  _mask = &block.template getData<DynamicsMask<DYNAMICS>>();
}

template <typename T, typename DESCRIPTOR, concepts::CellOperator OPERATOR>
ConcreteBlockO<T,DESCRIPTOR,Platform::GPU_SYCL,OPERATOR,OperatorScope::PerCell>::ConcreteBlockO():
  _cells(0),
  _modified{false}
{ }

template <typename T, typename DESCRIPTOR, concepts::CellOperator OPERATOR>
void ConcreteBlockO<T,DESCRIPTOR,Platform::GPU_SYCL,OPERATOR,OperatorScope::PerCell>::apply(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& block)
{
  if (_cells.size() > 0) {
    if (_modified) {
      OLB_PROFILE_T("GPU_SYCL:Post:Deduplicate", OPERATOR, DESCRIPTOR);
      _cells.deduplicate();
      _cells.setProcessingContext(ProcessingContext::Simulation);
      _modified = false;
    }
    OLB_PROFILE_T("GPU_SYCL:Post:Apply", OPERATOR, DESCRIPTOR);
    gpu::sycl::DeviceBlockLattice<T,DESCRIPTOR> lattice(block);
    gpu::sycl::async_call_list_operators(&gpu::sycl::device::getQueue(),
                                         {},
                                         lattice,
                                         _cells,
                                         gpu::sycl::ListedPostProcessor<OPERATOR>{});
  }
}

template <typename T, typename DESCRIPTOR, concepts::CellOperator OPERATOR>
ConcreteBlockO<T,DESCRIPTOR,Platform::GPU_SYCL,OPERATOR,OperatorScope::PerCellWithParameters>::ConcreteBlockO():
  _cells(0),
  _modified{false}
{ }

template <typename T, typename DESCRIPTOR, concepts::CellOperator OPERATOR>
void ConcreteBlockO<T,DESCRIPTOR,Platform::GPU_SYCL,OPERATOR,OperatorScope::PerCellWithParameters>::apply(
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& block)
{
  if (_cells.size() > 0) {
    if (_modified) {
      OLB_PROFILE_T("GPU_SYCL:Post:Deduplicate", OPERATOR, DESCRIPTOR);
      _cells.deduplicate();
      _cells.setProcessingContext(ProcessingContext::Simulation);
      _modified = false;
    }
    OLB_PROFILE_T("GPU_SYCL:Post:Apply", OPERATOR, DESCRIPTOR);
    using namespace gpu::sycl;
    DeviceBlockLattice<T,DESCRIPTOR> lattice(block);
    async_call_list_operators(&gpu::sycl::device::getQueue(),
                              {},
                              lattice,
                              _cells,
                              ListedPostProcessorWithParameters<T,DESCRIPTOR,OPERATOR>{_parameters->parameters});
  }
}


template <typename COUPLER, typename COUPLEES>
template <typename LATTICES>
ConcreteBlockCouplingO<COUPLEES,Platform::GPU_SYCL,COUPLER,OperatorScope::PerCell>
::ConcreteBlockCouplingO(LATTICES&& lattices)
  : _lattices{lattices}
{
  // If the lattices are introspectable, allocate fields accessed by the COUPLER
  bool allCoupledLatticesAreIntrospectable = true;
  COUPLEES::keys_t::for_each([&](auto id) {
    using name_t = typename decltype(id)::type;
    auto& lattice = *(_lattices.template get<name_t>());
    if (!lattice.isIntrospectable()) {
      allCoupledLatticesAreIntrospectable = false;
    }
  });

  if (allCoupledLatticesAreIntrospectable) {
    auto accessedFields = introspection::getFieldsAccessedByCoupling<COUPLER,COUPLEES>();
    COUPLEES::keys_t::for_each([&](auto id) {
      using name_t = typename decltype(id)::type;
      auto& lattice = *(_lattices.template get<name_t>());
      for (auto field : accessedFields.template get<name_t>()) {
        field.ensureAvailabilityIn(lattice);
      }
    });
  }
}

template <typename COUPLER, typename COUPLEES>
void ConcreteBlockCouplingO<COUPLEES,Platform::GPU_SYCL,COUPLER,OperatorScope::PerCell>::set(
  CellID iCell, bool state)
{
  if (!_mask) {
    _mask = std::make_unique<ConcreteBlockMask<typename COUPLEES::values_t::template get<0>::value_t,
                                               Platform::GPU_SYCL>>(
      _lattices.template get<0>()->template getData<CollisionSubdomainMask>()
    );
  }
  _mask->set(iCell, state);
}

template <typename COUPLER, typename COUPLEES>
void ConcreteBlockCouplingO<COUPLEES,Platform::GPU_SYCL,COUPLER,OperatorScope::PerCell>::apply()
{
  auto deviceLattice = _lattices.exchange_values([&](auto name) -> auto {
    return gpu::sycl::DeviceBlockLattice{*_lattices.get(name)};
  });
  if (_mask) {
    _mask->setProcessingContext(ProcessingContext::Simulation);
    gpu::sycl::call_coupling_operators(
      deviceLattice, _mask->deviceData(),
      gpu::sycl::UnmaskedCoupling<COUPLER>{});
  } else {
    auto& mask = _lattices.template get<0>()->template getData<CollisionSubdomainMask>();
    gpu::sycl::call_coupling_operators(
      deviceLattice, mask.deviceData(),
      gpu::sycl::UnmaskedCoupling<COUPLER>{});
  }
}

template <typename COUPLER, typename COUPLEES>
template <typename LATTICES>
ConcreteBlockCouplingO<COUPLEES,Platform::GPU_SYCL,COUPLER,OperatorScope::PerCellWithParameters>
::ConcreteBlockCouplingO(LATTICES&& lattices)
  : _lattices{lattices}
{
  // If the lattices are introspectable, allocate fields accessed by the COUPLER
  bool allCoupledLatticesAreIntrospectable = true;
  COUPLEES::keys_t::for_each([&](auto id) {
    using name_t = typename decltype(id)::type;
    auto& lattice = *(_lattices.template get<name_t>());
    if (!lattice.isIntrospectable()) {
      allCoupledLatticesAreIntrospectable = false;
    }
  });

  if (allCoupledLatticesAreIntrospectable) {
    auto accessedFields = introspection::getFieldsAccessedByCoupling<COUPLER,COUPLEES>();
    COUPLEES::keys_t::for_each([&](auto id) {
      using name_t = typename decltype(id)::type;
      auto& lattice = *(_lattices.template get<name_t>());
      for (auto field : accessedFields.template get<name_t>()) {
        field.ensureAvailabilityIn(lattice);
      }
    });
  }
}

template <typename COUPLER, typename COUPLEES>
void ConcreteBlockCouplingO<COUPLEES,Platform::GPU_SYCL,COUPLER,OperatorScope::PerCellWithParameters>::set(
  CellID iCell, bool state)
{
  if (!_mask) {
    _mask = std::make_unique<ConcreteBlockMask<typename COUPLEES::values_t::template get<0>::value_t,
                                               Platform::GPU_SYCL>>(
      _lattices.template get<0>()->template getData<CollisionSubdomainMask>()
    );
  }
  _mask->set(iCell, state);
}

template <typename COUPLER, typename COUPLEES>
void ConcreteBlockCouplingO<COUPLEES,Platform::GPU_SYCL,COUPLER,OperatorScope::PerCellWithParameters>::apply()
{
  auto deviceLattice = _lattices.exchange_values([&](auto name) -> auto {
    return gpu::sycl::DeviceBlockLattice{*_lattices.get(name)};
  });
  if (_mask) {
    _mask->setProcessingContext(ProcessingContext::Simulation);
    gpu::sycl::call_coupling_operators(
      deviceLattice, _mask->deviceData(),
      gpu::sycl::UnmaskedCouplingWithParameters<COUPLER,COUPLEES>{_parameters});
  } else {
    auto& mask = _lattices.template get<0>()->template getData<CollisionSubdomainMask>();
    gpu::sycl::call_coupling_operators(
      deviceLattice, mask.deviceData(),
      gpu::sycl::UnmaskedCouplingWithParameters<COUPLER,COUPLEES>{_parameters});
  }
}

template <typename COUPLER, typename COUPLEES>
void ConcreteBlockPointCouplingO<COUPLEES,Platform::GPU_SYCL,COUPLER,OperatorScope::PerCellWithParameters>::set(
  CellID iCell, bool state)
{
  throw std::bad_function_call();
}

template <typename COUPLER, typename COUPLEES>
void ConcreteBlockPointCouplingO<COUPLEES,Platform::GPU_SYCL,COUPLER,OperatorScope::PerCellWithParameters>::apply()
{
  if (_lattices.get(meta::id<names::Points>{})->getNcells() == 0) {
    return;
  }

  auto* particleLattice = _lattices.get(meta::id<names::Points>{});
  auto deviceLattice = _lattices.exchange_values([&](auto name) -> auto {
    if constexpr (std::is_same_v<typename decltype(name)::type, names::Points>) {
      return gpu::sycl::DeviceContext{*_lattices.get(name)};
    } else {
      return gpu::sycl::DeviceBlockLattice{*_lattices.get(name)};
    }
  });

  gpu::sycl::call_particle_coupling_operators(
    deviceLattice,
    gpu::sycl::ParticleCouplingWithParameters<COUPLER,COUPLEES>{_parameters});
}

}

/// State guarantees that internal structs are copyable to device
namespace sycl {

template <typename T, typename DESCRIPTOR, typename DYNAMICS>
struct is_device_copyable<olb::gpu::sycl::MaskedCollision<T, DESCRIPTOR, DYNAMICS>> : std::true_type { };

template <typename T, typename DESCRIPTOR, typename DYNAMICS>
struct is_device_copyable<olb::gpu::sycl::ListedCollision<T, DESCRIPTOR, DYNAMICS>> : std::true_type { };

template <typename T, typename DESCRIPTOR, typename OPERATOR>
struct is_device_copyable<olb::gpu::sycl::ListedPostProcessorWithParameters<T, DESCRIPTOR, OPERATOR>> : std::true_type { };

template <typename COUPLER, typename COUPLEES>
struct is_device_copyable<olb::gpu::sycl::UnmaskedCouplingWithParameters<COUPLER, COUPLEES>> : std::true_type { };

template <typename COUPLER, typename COUPLEES>
struct is_device_copyable<olb::gpu::sycl::ParticleCouplingWithParameters<COUPLER, COUPLEES>> : std::true_type { };

template <typename KEYS, typename VALUES>
struct is_device_copyable<olb::utilities::TypeIndexedTuple<olb::meta::plain_map<KEYS, VALUES>>> : std::true_type { };

}

#endif
