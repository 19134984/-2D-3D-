/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007-2024 Mathias J. Krause,
 *                Adrian Kummerlaender,
 *                Dennis Teutscher
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

#ifndef SUPER_LATTICE_HH
#define SUPER_LATTICE_HH

#include "utilities/log/channels.h"
#include "superLattice.h"

#include <cstdlib>
#include <string_view>
#include <array>
#include <vector>
#include <unordered_set>

#include "communication/mpiManager.h"
#include "communication/directionalOverlap.h"
#include "cell.h"
#include "io/base64.h"
#include "io/input.h"
#include "io/textTable.h"

#include "functors/analytical/analyticalF.h"

#include "functors/lattice/superBaseF2D.h"
#include "functors/lattice/superBaseF3D.h"
#include "functors/lattice/indicator/superIndicatorF2D.h"
#include "functors/lattice/indicator/superIndicatorF3D.h"

#include "io/serializerIO.h"

#include "geometry/cuboidDecomposition.h"
#include "geometry/superGeometry.hh"

#include "communication/loadBalancer.h"

#include "io/tupleParser.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::collectStatistics()
{
  _statisticsUpdateNeeded = false;

  executePostProcessors(stage::Statistics());

  using Result = typename LatticeStatistics<T,DESCRIPTOR>::Result;

  T sum_weight = 0;
  std::size_t sum_cells = 0;   // true summed cell count (not the volume weight)
  T average_rho = 0;
  T average_energy = 0;
  T mean_u = 0;
  T min_rho = std::numeric_limits<T>::max();
  T max_rho = std::numeric_limits<T>::lowest();
  T min_u   = std::numeric_limits<T>::max();
  T max_u   = 0;

  // Averages are volume-weighted (weight = cells * dx^3), min/max are unweighted
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    const T delta = this->_cuboidDecomposition.get(this->_loadBalancer.glob(iC)).getDeltaR();
    const Result r = _block[iC]->getStatistics().getResult();
    if (r.nCells == 0) {
      continue;
    }
    const T weight = r.nCells * delta * delta * delta;
    sum_weight     += weight;
    sum_cells      += r.nCells;
    average_rho    += r.avRho    * weight;
    average_energy += r.avEnergy * weight;
    mean_u         += r.meanU    * weight;
    min_rho = std::min(min_rho, r.minRho);
    max_rho = std::max(max_rho, r.maxRho);
    min_u   = std::min(min_u,   r.minU);
    max_u   = std::max(max_u,   r.maxU);
  }

#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(sum_weight, MPI_SUM);
  singleton::mpi().reduceAndBcast(sum_cells, MPI_SUM);
  singleton::mpi().reduceAndBcast(average_rho, MPI_SUM);
  singleton::mpi().reduceAndBcast(average_energy, MPI_SUM);
  singleton::mpi().reduceAndBcast(mean_u, MPI_SUM);
  singleton::mpi().reduceAndBcast(min_rho, MPI_MIN);
  singleton::mpi().reduceAndBcast(max_rho, MPI_MAX);
  singleton::mpi().reduceAndBcast(min_u, MPI_MIN);
  singleton::mpi().reduceAndBcast(max_u, MPI_MAX);
#endif

  if (sum_weight > 0) {
    average_rho    /= sum_weight;
    average_energy /= sum_weight;
    mean_u         /= sum_weight;
  }

  const Result global{average_rho, average_energy, max_u,
                      min_rho, max_rho, min_u, mean_u, sum_cells};
  getStatistics().setResult(global);
  getStatistics().incrementTime();

  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->getStatistics().setResult(global);
    _block[iC]->getStatistics().incrementTime();
  }

  // Combine the per-block registered-field reductions (filled above by the
  // FieldStatisticO operators) across blocks and ranks into the super statistics
  using FieldResult = typename LatticeStatistics<T,DESCRIPTOR>::FieldResult;
  for (const auto& rf : getStatistics().getRegisteredFields()) {
    FieldResult combined;
    combined.components = rf.components;
    combined.count = 0;
    combined.min.assign(rf.components, std::numeric_limits<T>::max());
    combined.max.assign(rf.components, std::numeric_limits<T>::lowest());
    combined.sum.assign(rf.components, T{});

    for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
      auto& results = _block[iC]->getStatistics().fieldResults();
      const auto it = results.find(rf.id);
      if (it == results.end() || it->second.count == 0) {
        continue;
      }
      const auto& b = it->second;
      combined.count += b.count;
      for (unsigned c = 0; c < rf.components; ++c) {
        combined.min[c] = std::min(combined.min[c], b.min[c]);
        combined.max[c] = std::max(combined.max[c], b.max[c]);
        combined.sum[c] += b.sum[c];
      }
    }

#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(combined.count, MPI_SUM);
    for (unsigned c = 0; c < rf.components; ++c) {
      singleton::mpi().reduceAndBcast(combined.min[c], MPI_MIN);
      singleton::mpi().reduceAndBcast(combined.max[c], MPI_MAX);
      singleton::mpi().reduceAndBcast(combined.sum[c], MPI_SUM);
    }
#endif

    getStatistics().fieldResults()[rf.id] = combined;
  }
}

template<typename T, typename DESCRIPTOR>
SuperLattice<T,DESCRIPTOR>::SuperLattice(CuboidDecomposition<T,DESCRIPTOR::d>& decomposition,
                                         LoadBalancer<T>& loadBalancer,
                                         unsigned overlap)
  : SuperStructure<T,DESCRIPTOR::d>(decomposition,
                                    loadBalancer,
                                    overlap)
  , _statistics()
{
  using namespace stage;

  auto& load = this->getLoadBalancer();

  OLB_LOG(logging::Setup, Debug, "constructing blocks");

  for (int iC = 0; iC < load.size(); ++iC) {
    auto& cuboid = this->_cuboidDecomposition.get(load.glob(iC));
    #ifdef PLATFORM_GPU_CUDA
    if (load.platform(iC) == Platform::GPU_CUDA) {
      if (gpu::cuda::device::getCount() == 0) {
        throw std::runtime_error("Load balancer requested GPU processing for cuboid "
                               + std::to_string(load.glob(iC))
                               + " but no CUDA device was found on rank "
                               + std::to_string(singleton::mpi().getRank()));

      }
    }
    #endif

    #ifdef PLATFORM_GPU_HIP
    if (load.platform(iC) == Platform::GPU_HIP) {
      if (gpu::hip::device::getCount() == 0) {
        throw std::runtime_error("Load balancer requested GPU processing for cuboid "
                               + std::to_string(load.glob(iC))
                               + " but no HIP device was found on rank "
                               + std::to_string(singleton::mpi().getRank()));

      }
    }
    #endif

    //#ifdef PLATFORM_GPU_SYCL
    //if (load.platform(iC) == Platform::GPU_SYCL) {
    //  if (gpu::sycl::device::getCount() == 0) {
    //    throw std::runtime_error("Load balancer requested GPU processing for cuboid "
    //                           + std::to_string(load.glob(iC))
    //                           + " but no SYCL device was found on rank "
    //                           + std::to_string(singleton::mpi().getRank()));

    //  }
    //}
    //#endif

    _block.emplace_back(constructUsingConcretePlatform<ConcretizableBlockLattice<T,DESCRIPTOR>>(
      load.platform(iC), cuboid, this->getOverlap()));
  }

  OLB_LOG(logging::Setup, Debug, "blocks constructed");

#ifdef PARALLEL_MODE_MPI
  OLB_LOG(logging::Setup, Debug, "MPI barrier");
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  {
    OLB_LOG(logging::Setup, Debug, "setting up PostCollide communicator");
    auto& communicator = getCommunicator(PostCollide());
    // The directional streaming exchange transmits only the populations
    // actually consumed by inter-block propagation (~4x less data on faces,
    // see docs/design/directional-granular-halo.md). Default for descriptors
    // whose per-cell population mask fits a std::uint64 (q <= 64); disable at
    // compile time via FEATURES=DISABLE_DIRECTIONAL_STREAMING_COMM.
#ifndef FEATURE_DISABLE_DIRECTIONAL_STREAMING_COMM
    if constexpr (DESCRIPTOR::q <= 64) {
      _directionalStreaming = true;
      DESCRIPTOR::template filter<descriptors::is_propagatable_field>
                ::for_each([&](auto field) {
        communicator.template requestField<typename decltype(field)::type>(
          overlap::streamingMask<DESCRIPTOR>());
      });
    } else
#endif
    {
      DESCRIPTOR::template filter<descriptors::is_propagatable_field>
                ::for_each([&](auto field) {
        communicator.template requestField<typename decltype(field)::type>();
      });
    }
    communicator.requestOverlap(1); // Required for inter-block propagation
    // Negotiation is deferred to initialize() so that setStreamingReversible
    // can still upgrade the exchanged population sets
  }
  OLB_LOG(logging::Setup, Debug, "PostCollide communicator ready");

  {
    OLB_LOG(logging::Setup, Debug, "exchanging Full requests");
    auto& communicator = getCommunicator(Full());
    DESCRIPTOR::fields_t::for_each([&](auto field) {
      communicator.template requestField<typename decltype(field)::type>();
    });
    // VTK output includes overlap of 1, some implicit dependencies for overlap of 2 exist
    communicator.requestOverlap(std::min(2, this->getOverlap()));
    communicator.exchangeRequests();
    OLB_LOG(logging::Setup, Debug, "Full requests exchanged");
  }

  _statisticsEnabled = true;
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
SuperLattice<T,DESCRIPTOR>::SuperLattice(UnitConverter<T,DESCRIPTOR> converter,
                                         CuboidDecomposition<T,DESCRIPTOR::d>& decomposition,
                                         LoadBalancer<T>& loadBalancer,
                                         unsigned overlap)
  : SuperLattice(decomposition,
                 loadBalancer,
                 overlap)
{
  setUnitConverter<DESCRIPTOR>(converter);
}


template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::setProcessingContext(ProcessingContext context, bool communicate)
{
  // Join all device work before host-side copies observe it
  #ifdef PLATFORM_GPU_CUDA
  gpu::cuda::device::synchronize();
  #endif
  #ifdef PLATFORM_GPU_HIP
  gpu::hip::device::synchronize();
  #endif
  #ifdef PLATFORM_GPU_SYCL
  gpu::sycl::device::synchronize();
  #endif

  // Communicate overlap prior to evaluation
  if (context == ProcessingContext::Evaluation) {
    if (communicate) {
      this->communicate();
    }
    waitForBackgroundTasks<stage::PreContextSwitchTo<ProcessingContext::Evaluation>>();
  }
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->setProcessingContext(context);
  }
}

template<typename T, typename DESCRIPTOR>
template<typename FIELD_TYPE>
void SuperLattice<T,DESCRIPTOR>::setProcessingContext(ProcessingContext context) {
  // Join all device work before host-side copies observe it
  #ifdef PLATFORM_GPU_CUDA
  gpu::cuda::device::synchronize();
  #endif
  #ifdef PLATFORM_GPU_HIP
  gpu::hip::device::synchronize();
  #endif
  #ifdef PLATFORM_GPU_SYCL
  gpu::sycl::device::synchronize();
  #endif

  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    if (_block[iC]->template hasData<FIELD_TYPE>()) {
      _block[iC]->template getData<FIELD_TYPE>().setProcessingContext(context);
    }
  }
}

template<typename T, typename DESCRIPTOR>
Cell<T,DESCRIPTOR> SuperLattice<T,DESCRIPTOR>::get(LatticeR<DESCRIPTOR::d+1> latticeR)
{
#ifdef PARALLEL_MODE_MPI
  if (this->_loadBalancer.isLocal(latticeR[0])) {
    if constexpr (DESCRIPTOR::d == 3) {
      return _block[this->_loadBalancer.loc(latticeR[0])]->get(latticeR[1], latticeR[2], latticeR[3]);
    } else {
      return _block[this->_loadBalancer.loc(latticeR[0])]->get(latticeR[1], latticeR[2]);
    }
  }
  else {
    throw std::domain_error("Cuboid iC must be locally available");
  }
#else
  if constexpr (DESCRIPTOR::d == 3) {
    return _block[this->_loadBalancer.loc(latticeR[0])]->get(latticeR[1], latticeR[2], latticeR[3]);
  } else {
    return _block[this->_loadBalancer.loc(latticeR[0])]->get(latticeR[1], latticeR[2]);
  }
#endif
}

template<typename T, typename DESCRIPTOR>
template<typename... R>
std::enable_if_t<sizeof...(R) == DESCRIPTOR::d+1, Cell<T,DESCRIPTOR>>
SuperLattice<T,DESCRIPTOR>::get(R... latticeR)
{
  return get({latticeR...});
}

template<typename T, typename DESCRIPTOR>
template <Platform PLATFORM, typename F>
void SuperLattice<T,DESCRIPTOR>::forBlocksOnPlatform(F f) {
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    if (_block[iC]->getPlatform() == PLATFORM) {
      f(static_cast<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>&>(*_block[iC]));
    }
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::initialize()
{
  if (!_initialized) {
    if (!_converter) {
      throw std::runtime_error("Each SuperLattice must be assigned a UnitConverter");
    }

    getCommunicator(stage::PostCollide()).exchangeRequests();

    // Resolve the per-cuboid collision execution policy from the load balancer.
    // The overlap axis (outer/inner split) is resolved here where the
    // decomposition topology is available. The dispatch axes stay Auto and are
    // resolved per operator from mask density at first apply.
    const auto& decomposition = this->getCuboidDecomposition();
    for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
      auto policy = this->_loadBalancer.collisionPolicy(iC);
      policy.overlap = this->_loadBalancer.splitsCollision(iC, decomposition, this->getOverlap())
                       ? CollisionOverlap::Split
                       : CollisionOverlap::Fused;
      _block[iC]->setCollisionPolicy(policy);
    }

    for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
      _block[iC]->initialize();
    }

    // Globally reduce whether any rank consumes live statistics
    int consumesLiveStatistics = 0;
    for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
      consumesLiveStatistics |= _block[iC]->consumesLiveStatistics();
    }
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(consumesLiveStatistics, MPI_MAX);
#endif
    _consumesLiveStatistics = consumesLiveStatistics != 0;

    // Install the registered-field reduction operators now that the blocks
    // (and their statistics) are initialized
    if (!_statisticsFieldsInstalled) {
      _statisticsFieldsInstalled = true;
      for (const auto& rf : getStatistics().getRegisteredFields()) {
        rf.addOperator(*this);
      }
    }

    // Write the per-component introspection summary unless suppressed centrally
    if (!_name.empty()
        && !singleton::input().resolveBuiltin<parameters::no_introspection>(false)) {
      writeSummary(io::text::filenameSafe(_name));
    }

    _communicationNeeded = true;
    _initialized = true;
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::printSummary() const
{
  _block[0]->writeDescription(this->clout);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::writeSummary(std::string fileName) const
{
  const auto& cGeometry = this->getCuboidDecomposition();

  // Map each descriptor field's type_index to its component count and byte size,
  // so the communication manifest can derive exchanged bytes. Field names come
  // from the communicator itself (SuperCommunicator::getFieldName), which also
  // covers non-descriptor fields a coupling might request.
  struct FieldMeta { unsigned components; std::size_t componentBytes; };
  std::map<std::type_index, FieldMeta> fieldInfo;
  DESCRIPTOR::fields_t::for_each([&](auto id) {
    using FIELD = typename decltype(id)::type;
    fieldInfo[typeid(FIELD)] = FieldMeta{ DESCRIPTOR::template size<FIELD>(),
                                          sizeof(typename FIELD::template value_type<T>) };
  });

  // Demangled names for the communication stages (type_index keys of _communicator).
  const std::map<std::type_index, std::string_view> stageNames = {
    {typeid(stage::PreCollide),           meta::nice_name<stage::PreCollide>()},
    {typeid(stage::PostCollide),          meta::nice_name<stage::PostCollide>()},
    {typeid(stage::PostStream),           meta::nice_name<stage::PostStream>()},
    {typeid(stage::PostPostProcess),      meta::nice_name<stage::PostPostProcess>()},
    {typeid(stage::IterativePostProcess), meta::nice_name<stage::IterativePostProcess>()},
    {typeid(stage::PreCoupling),          meta::nice_name<stage::PreCoupling>()},
    {typeid(stage::PostCoupling),         meta::nice_name<stage::PostCoupling>()},
    {typeid(stage::Full),                 meta::nice_name<stage::Full>()},
  };

  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    const int globC = this->_loadBalancer.glob(iC);
    {
      std::string dataFile = singleton::directories().getIntrospectionOutDir()
                           + fileName + "_globC_" + std::to_string(globC)
                           + ".log";
      std::ofstream fout(dataFile.c_str(), std::ios::trunc);
      if (!fout) {
        this->clout << "writeSummary error: can not open std::ofstream" << std::endl;
      } else {
        cGeometry.get(globC).write(fout);
        fout << std::endl;
        _block[iC]->writeDescription(fout);
        fout.close();
      }
    }
    {
      std::string dataFile = singleton::directories().getIntrospectionOutDir()
                           + fileName + "_globC_" + std::to_string(globC)
                           + ".dynamics";
      std::ofstream fout(dataFile.c_str(), std::ios::trunc);
      if (!fout) {
        this->clout << "writeSummary error: can not open std::ofstream" << std::endl;
      } else {
        _block[iC]->writeDynamicsAsCSV(fout);
        fout.close();
      }
    }
    {
      std::string dataFile = singleton::directories().getIntrospectionOutDir()
                           + fileName + "_globC_" + std::to_string(globC)
                           + ".operator";
      std::ofstream fout(dataFile.c_str(), std::ios::trunc);
      if (!fout) {
        this->clout << "writeSummary error: can not open std::ofstream" << std::endl;
      } else {
        _block[iC]->writeOperatorAsCSV(fout);
        fout.close();
      }
    }
    {
      std::string dataFile = singleton::directories().getIntrospectionOutDir()
                           + fileName + "_globC_" + std::to_string(globC)
                           + ".fields";
      std::ofstream fout(dataFile.c_str(), std::ios::trunc);
      if (!fout) {
        this->clout << "writeSummary error: can not open std::ofstream" << std::endl;
      } else {
        _block[iC]->writeFieldsAsCSV(fout);
        fout.close();
      }
    }
    {
      std::string dataFile = singleton::directories().getIntrospectionOutDir()
                           + fileName + "_globC_" + std::to_string(globC)
                           + ".communication";
      std::ofstream fout(dataFile.c_str(), std::ios::trunc);
      if (!fout) {
        this->clout << "writeSummary error: can not open std::ofstream" << std::endl;
      } else {
        fout << "stage; remoteCuboid; remoteRank; field; requestMode; effectiveCells; "
                "carriedPopulations; bytes; haloTiles" << std::endl;
        for (const auto& [stageId, communicator] : _communicator) {
          const auto stageIt = stageNames.find(stageId);
          const std::string_view stageName = stageIt != stageNames.end() ? stageIt->second
                                                                         : std::string_view(stageId.name());
          const auto& neighborhood = communicator->getBlockNeighborhood(iC);
          neighborhood.forNeighbors([&](int remoteC) {
            const int remoteRank = this->_loadBalancer.rank(remoteC);
            const auto& inboundCells = neighborhood.getCellsInboundFrom(remoteC);
            neighborhood.forFieldCommSummary(remoteC, [&](std::type_index field, const auto& summary) {
              const bool dense = summary.carriedPopulations == 0;
              const std::string_view fieldName = communicator->getFieldName(field);
              const auto info = fieldInfo.find(field);
              std::size_t bytes = 0;
              if (info != fieldInfo.end()) {
                bytes = dense ? summary.effectiveCells * info->second.components * info->second.componentBytes
                              : summary.carriedPopulations * info->second.componentBytes;
              }
              // Halo tiles a sparse field materialises for this neighbor's overlap cells (-1 if dense)
              const long haloTiles = _block[iC]->countHaloTilesForField(field, inboundCells);
              fout << stageName << "; " << remoteC << "; " << remoteRank << "; "
                   << fieldName << "; " << (dense ? "dense" : "selector") << "; "
                   << summary.effectiveCells << "; "
                   << (dense ? -1 : static_cast<long>(summary.carriedPopulations)) << "; "
                   << bytes << "; " << haloTiles << std::endl;
            });
          });
        }
        fout.close();
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
template<typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::defineDynamics()
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->template defineDynamics<DYNAMICS>();
  }
}

template<typename T, typename DESCRIPTOR>
template<template<typename...> typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::defineDynamics()
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->template defineDynamics<DYNAMICS<T,DESCRIPTOR>>();
  }
}

template<typename T, typename DESCRIPTOR>
template<typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::defineDynamics(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->template defineDynamics<DYNAMICS>(indicator->getBlockIndicatorF(iC));
  }
}

template<typename T, typename DESCRIPTOR>
template<typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::defineDynamics(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material)
{
  defineDynamics<DYNAMICS>(sGeometry.getMaterialIndicator(material));
}

template<typename T, typename DESCRIPTOR>
template<template<typename...> typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::defineDynamics(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->template defineDynamics<DYNAMICS<T,DESCRIPTOR>>(indicator->getBlockIndicatorF(iC));
  }
}

template<typename T, typename DESCRIPTOR>
template<template<typename...> typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::defineDynamics(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material)
{
  defineDynamics<DYNAMICS>(sGeometry.getMaterialIndicator(material));
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineDynamics(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                DynamicsPromise<T,DESCRIPTOR>&& promise)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->defineDynamics(indicator->getBlockIndicatorF(iC), std::forward<DynamicsPromise<T,DESCRIPTOR>&&>(promise));
  }
}


template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineRho(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                           AnalyticalF<DESCRIPTOR::d,T,T>& rho)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->defineRho(indicator->getBlockIndicatorF(iC),
                         rho);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineRho(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                           AnalyticalF<DESCRIPTOR::d,T,T>& rho)
{
  defineRho(sGeometry.getMaterialIndicator(material), rho);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineU(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                         AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  #ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for schedule(dynamic,1)
  #endif
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->defineU(indicator->getBlockIndicatorF(iC),
                        u);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineU(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                         AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  defineU(sGeometry.getMaterialIndicator(material), u);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineRhoU(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                            AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                            AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  #ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for schedule(dynamic,1)
  #endif
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->defineRhoU(indicator->getBlockIndicatorF(iC),
                          rho, u);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineRhoU(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                            AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                            AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  defineRhoU(sGeometry.getMaterialIndicator(material), rho, u);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::definePopulations(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                   AnalyticalF<DESCRIPTOR::d,T,T>& Pop)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->definePopulations(indicator->getBlockIndicatorF(iC), Pop);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::definePopulations(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                                   AnalyticalF<DESCRIPTOR::d,T,T>& Pop)
{
  definePopulations(sGeometry.getMaterialIndicator(material), Pop);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::definePopulations(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                   SuperF<DESCRIPTOR::d,T,T>& Pop)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->definePopulations(indicator->getBlockIndicatorF(iC), Pop.getBlockF(iC));
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::definePopulations(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                                   SuperF<DESCRIPTOR::d,T,T>& Pop)
{
  definePopulations(sGeometry.getMaterialIndicator(material), Pop);
}


template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::defineField(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                             FunctorPtr<SuperF<DESCRIPTOR::d,T,T>>&& field)
{
  if (field->getBlockFSize() == this->_loadBalancer.size()) {
    for (int iC=0; iC < this->_loadBalancer.size(); ++iC) {
      _block[iC]->template defineField<FIELD>(indicator->getBlockIndicatorF(iC),
                                              field->getBlockF(iC));
    }
  }
  else {
    FieldD<T,DESCRIPTOR,FIELD> fieldTmp;
    int coords[DESCRIPTOR::d+1];
    for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
      _block[iC]->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
        coords[0] = iC;
        for (int j = 0; j < DESCRIPTOR::d; ++j) {
          coords[j+1] = loc[j];
        }
        if (indicator(coords)) {
          field(fieldTmp.data(), coords);
          _block[iC]->get(loc).template setField<FIELD>(fieldTmp);
        }
      });
    }
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::defineField(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                             AnalyticalF<DESCRIPTOR::d,T,T>& field)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->template defineField<FIELD>(indicator->getBlockIndicatorF(iC),
                                            field);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::defineField(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                             FunctorPtr<SuperF<DESCRIPTOR::d,T,T>>&& field)
{
  defineField<FIELD>(sGeometry.getMaterialIndicator(material), std::forward<decltype(field)>(field));
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::defineField(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                             AnalyticalF<DESCRIPTOR::d,T,T>& field)

{
  defineField<FIELD>(sGeometry.getMaterialIndicator(material), field);
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::defineField(SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
                                             IndicatorF<T,DESCRIPTOR::d>& indicator,
                                             AnalyticalF<DESCRIPTOR::d,T,T>& field)
{
  SuperIndicatorFfromIndicatorF<T,DESCRIPTOR::d> indicatorF(indicator, sGeometry);
  defineField<FIELD>(indicatorF, field);
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::defineField(SuperGeometry<T,DESCRIPTOR::d>& sGeometry,
                                             IndicatorF<T,DESCRIPTOR::d>& indicator,
                                             FunctorPtr<SuperF<DESCRIPTOR::d,T,T>>&& field)
{
  SuperIndicatorFfromIndicatorF<T,DESCRIPTOR::d> indicatorF(indicator, sGeometry);
  defineField<FIELD>(indicatorF, std::forward<decltype(field)>(field));
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::allocateField(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator)
{
  int coords[DESCRIPTOR::d+1];
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    auto& field = _block[iC]->template getField<FIELD>();
    _block[iC]->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      coords[0] = iC;
      for (int j = 0; j < DESCRIPTOR::d; ++j) {
        coords[j+1] = loc[j];
      }
      if (indicator(coords)) {
        field.allocateRow(_block[iC]->getCellId(loc));
      }
    });
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::allocateField(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material)
{
  allocateField<FIELD>(sGeometry.getMaterialIndicator(material));
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::allocateField()
{
  static_assert(storageOf<DESCRIPTOR,FIELD>() != FieldStorage::Sparse,
                "sparse FIELD reserves no tiles here, so device writes are discarded; "
                "use defineField(region, value) to reserve the write region");
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->template getField<FIELD>();
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::setFieldStorage(FieldStorage mode)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->template getField<FIELD>().setStorageMode(mode);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::setTileShape(sparse::TileShape shape, unsigned capacity)
{
  [[maybe_unused]] const unsigned product = shape.tx * shape.ty * shape.tz;
  // Runtime capacity (cells per tile) equals the shape product and is bounded by the
  // field's compile-time tile capacity. A given `capacity` only cross-checks that product.
  OLB_PRECONDITION(capacity == 0 || capacity == product);
  OLB_PRECONDITION(product <= tileCapacityOf<FIELD>());
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    const auto extent = _block[iC]->getExtent();
    const int  pad    = _block[iC]->getPadding();
    const unsigned sx = extent[0] + 2*pad;
    const unsigned sy = DESCRIPTOR::d > 1 ? extent[1] + 2*pad : 1;
    const unsigned sz = DESCRIPTOR::d > 2 ? extent[2] + 2*pad : 1;
    std::size_t p0 = 1;
    std::size_t p1 = 1;
    if constexpr (DESCRIPTOR::d == 3) {
      p0 = std::size_t(sy)*sz;
      p1 = sz;
    } else if constexpr (DESCRIPTOR::d == 2) {
      p0 = sy;
    }
    _block[iC]->template getField<FIELD>().setTileShape(shape, p0, p1, sx, sy, sz, pad);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::setFieldStorage(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& region,
                                                 double denseThreshold)
{
  // Candidate (capacity x shape) pairs for auto-selection (3D only). Each candidate's tile
  // capacity is its shape product (tx*ty*tz); capacities span a bounded power-of-two set.
  // Coalescing floor: no capacity below 32, since a sub-warp tile wastes coalesced device
  // bandwidth, so smaller tiles are never offered even if they would pack marginally tighter.
  static const std::vector<sparse::TileShape> candidates{
    // capacity 32
    {4,4,2}, {8,4,1}, {2,2,8}, {1,1,32},
    // capacity 64
    {1,1,64}, {4,4,4}, {8,8,1}, {2,2,16},
    // capacity 128
    {4,4,8}, {8,8,2}, {2,4,16},
    // capacity 256
    {8,8,4}, {4,8,8}, {2,8,16},
  };
  // Only pairs whose capacity fits the field's compile-time tile capacity are eligible
  constexpr unsigned maxCapacity = tileCapacityOf<FIELD>();
  const std::size_t nCand = candidates.size();
  std::vector<std::vector<std::unordered_set<std::size_t>>> tilesPerCandidate;
  std::vector<std::size_t> gridPerCandidate;
  if constexpr (DESCRIPTOR::d == 3) {
    tilesPerCandidate.resize(nCand);
    gridPerCandidate.assign(nCand, 0);
    for (auto& perBlock : tilesPerCandidate) {
      perBlock.resize(this->_loadBalancer.size());
    }
  }

  int coords[DESCRIPTOR::d+1];
  std::size_t inRegion = 0;
  std::size_t total    = 0;
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    const auto extent = _block[iC]->getExtent();
    const int  pad    = _block[iC]->getPadding();
    const unsigned sx = extent[0] + 2*pad;
    const unsigned sy = DESCRIPTOR::d > 1 ? extent[1] + 2*pad : 1;
    const unsigned sz = DESCRIPTOR::d > 2 ? extent[2] + 2*pad : 1;
    std::vector<sparse::TileData> geos;
    if constexpr (DESCRIPTOR::d == 3) {
      geos.reserve(nCand);
      for (std::size_t s=0; s < nCand; ++s) {
        geos.push_back(sparse::makeTileData(candidates[s], std::size_t(sy)*sz, sz,
                                              sx, sy, sz, maxCapacity, unsigned(pad)));
        const unsigned cap = candidates[s].tx * candidates[s].ty * candidates[s].tz;
        // Directory spans the whole padded tile grid. The linear scheme keeps ceil(cells/cap)
        gridPerCandidate[s] += geos.back().linear
          ? (std::size_t(sx)*sy*sz + cap - 1) / cap
          : geos.back().gridTileCount();
      }
    }
    _block[iC]->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      // Only core cells count toward the fraction. Halo cells would inflate the total
      if (!_block[iC]->isInsideCore(loc)) {
        return;
      }
      coords[0] = iC;
      for (int j = 0; j < DESCRIPTOR::d; ++j) {
        coords[j+1] = loc[j];
      }
      ++total;
      if (region(coords)) {
        ++inRegion;
        if constexpr (DESCRIPTOR::d == 3) {
          const std::size_t cell = _block[iC]->getCellId(loc);
          for (std::size_t s=0; s < nCand; ++s) {
            tilesPerCandidate[s][iC].insert(geos[s].tileId(cell, maxCapacity));
          }
        }
      }
    });
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(inRegion, MPI_SUM);
  singleton::mpi().reduceAndBcast(total,    MPI_SUM);
#endif
  const bool dense = total > 0 && double(inRegion) > denseThreshold * double(total);
  if (dense) {
    setFieldStorage<FIELD>(FieldStorage::Dense);
    return;
  }
  if constexpr (DESCRIPTOR::d == 3) {
    // Pick the (capacity, shape) pair minimising predicted device memory: payload cells
    // (occupied tiles * capacity) plus the tile directory (grid tiles). Capacity is the
    // shape product. Pairs exceeding the field's compile-time capacity are skipped.
    std::vector<std::size_t> tileCount(nCand, 0);
    for (std::size_t s=0; s < nCand; ++s) {
      for (const auto& perBlock : tilesPerCandidate[s]) {
        tileCount[s] += perBlock.size();
      }
    }
#ifdef PARALLEL_MODE_MPI
    for (std::size_t s=0; s < nCand; ++s) {
      singleton::mpi().reduceAndBcast(tileCount[s],     MPI_SUM);
      singleton::mpi().reduceAndBcast(gridPerCandidate[s], MPI_SUM);
    }
#endif
    std::size_t best = 0;
    double bestBytes = -1;
    for (std::size_t s=0; s < nCand; ++s) {
      const unsigned cap = candidates[s].tx * candidates[s].ty * candidates[s].tz;
      if (cap > maxCapacity) {
        continue;
      }
      const double bytes = double(tileCount[s]) * cap * sizeof(T)
                         + double(gridPerCandidate[s]) * sizeof(std::int32_t);
      if (bestBytes < 0 || bytes < bestBytes) {
        bestBytes = bytes;
        best = s;
      }
    }
    setTileShape<FIELD>(candidates[best]);
  }
  setFieldStorage<FIELD>(FieldStorage::Sparse);
  allocateField<FIELD>(std::move(region));
}

template<typename T, typename DESCRIPTOR>
template <typename PARAMETER>
void SuperLattice<T,DESCRIPTOR>::setParameter(FieldD<T,DESCRIPTOR,PARAMETER> field)
{
  if (!PARAMETER::template isValid<T,DESCRIPTOR,PARAMETER>(field)) {
    const std::string error = "Value of parameter " + fields::name<PARAMETER>() + " is invalid";
    throw std::invalid_argument(error);
  }
  for (int iC=0; iC < this->getLoadBalancer().size(); ++iC) {
    _block[iC]->template setParameter<PARAMETER>(field);
  }
}

template<typename T, typename DESCRIPTOR>
template <typename PARAMETER, typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::setParameterOfDynamics(FieldD<T,DESCRIPTOR,PARAMETER>&& field)
{
  for (int iC=0; iC < this->getLoadBalancer().size(); ++iC) {
    _block[iC]->template getData<OperatorParameters<DYNAMICS>>().template set<PARAMETER>(
       std::forward<decltype(field)>(field));
  }
}

template<typename T, typename DESCRIPTOR>
template <typename PARAMETER, template<typename...> typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::setParameterOfDynamics(FieldD<T,DESCRIPTOR,PARAMETER>&& field)
{
  setParameterOfDynamics<PARAMETER,DYNAMICS<T,DESCRIPTOR>>(std::forward<decltype(field)>(field));
}

template<typename T, typename DESCRIPTOR>
template <typename PARAMETER, typename _DESCRIPTOR, typename FIELD>
void SuperLattice<T,DESCRIPTOR>::setParameter(AbstractFieldArrayD<T,_DESCRIPTOR,FIELD>& abstractFieldArray)
{
  static_assert(DESCRIPTOR::template size<PARAMETER>() == DESCRIPTOR::template size<FIELD>(),
                "PARAMETER size must equal FIELD size");
  static_assert(std::is_same_v<typename PARAMETER::template value_type<T>,
                               typename FIELD::template value_type<T>*>,
                "PARAMETER must store pointers to FIELD components");
  // auto test = PARAMETER::template getMaxValue<T,DESCRIPTOR>();
  // std::cout << test<< std::endl;
  for (int iC=0; iC < this->getLoadBalancer().size(); ++iC) {
    _block[iC]->template setParameter<PARAMETER>(abstractFieldArray);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::iniEquilibrium(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->iniEquilibrium(indicator->getBlockIndicatorF(iC), rho, u);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::iniEquilibrium(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  iniEquilibrium(sGeometry.getMaterialIndicator(material), rho, u);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::iniEquilibrium(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                SuperF<DESCRIPTOR::d,T,T>& u)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->iniEquilibrium(indicator->getBlockIndicatorF(iC), rho, u.getBlockF(iC));
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::iniEquilibrium(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                SuperF<DESCRIPTOR::d,T,T>& u)
{
  iniEquilibrium(sGeometry.getMaterialIndicator(material), rho, u);
}


template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::iniRegularized(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& u,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& pi)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->iniRegularized(indicator->getBlockIndicatorF(iC), rho, u, pi);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::iniRegularized(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& u,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& pi)
{
  iniRegularized(sGeometry.getMaterialIndicator(material), rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::collideAndStream()
{
  OLB_PROFILE_STEP();

  collide();
  stream();
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::collide()
{
  if (!_initialized) [[unlikely]] {
    initialize();
  }

  using namespace stage;

  waitForBackgroundTasks(PreCollide());
  auto& load = this->_loadBalancer;

  if (_consumesLiveStatistics) {
    auto& stats = getStatistics();
    statistics::live_parameters::for_each([&](auto id) {
      using PARAM = typename decltype(id)::type;
      static_assert(concepts::StatisticsField<PARAM,T,DESCRIPTOR>);
      setParameter<PARAM>(PARAM::compute(stats));
    });
  }

  // Optional pre processing stage
  executePostProcessors(PreCollide());

  // Execute custom tasks (arbitrary callables)
  // (used for multi-stage models such as bubble model)
  executeCustomTasks(PreCollide());

  // Start the streaming halo gather, hidden behind the bulk collision
  getCommunicator(PostCollide()).prepareCommunicate();

  #ifdef PARALLEL_MODE_OMP
  #pragma omp taskloop
  #endif
  for (int iC = 0; iC < load.size(); ++iC) {
    _block[iC]->collide();
  }

  // Complete the streaming halo exchange
  getCommunicator(PostCollide()).continueCommunicate();

  executePostProcessors(PostCollide());
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::stream()
{
  using namespace stage;
  auto& load = this->_loadBalancer;

  getCommunicator(PreStream()).communicate();

  // Join outer collision before the population views are rotated
  #ifdef PLATFORM_GPU_CUDA
  gpu::cuda::device::synchronizeOuterCollision();
  #endif
  #ifdef PLATFORM_GPU_HIP
  gpu::hip::device::synchronizeOuterCollision();
  #endif

  // Block-local propagation
  for (int iC = 0; iC < load.size(); ++iC) {
    _block[iC]->stream();
  }

  // Communicate (default) post processor neighborhood and apply them
  executePostProcessors(PostStream());

  // Execute custom tasks (arbitrary callables)
  // (used for multi-stage models such as free surface)
  executeCustomTasks(PostStream());

  // Final communication stage (e.g. for external coupling)
  getCommunicator(PostPostProcess()).communicate();

  if (_statisticsEnabled) {
    _statisticsUpdateNeeded = true;
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::requireDenseStreamingExchange() {
  if (_directionalStreaming) {
    OLB_LOG(logging::Setup, Debug, "directional streaming disabled at runtime; using dense halo exchange");
    auto& communicator = getCommunicator(stage::PostCollide());
    DESCRIPTOR::template filter<descriptors::is_propagatable_field>
              ::for_each([&](auto field) {
      communicator.template requestField<typename decltype(field)::type>();
    });
    if (_initialized) {
      communicator.exchangeRequests();
    }
    _directionalStreaming = false;
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::revertStreaming() {
  // Reverted propagation consumes the opposite population sets. Dense
  // covers both directions, making any further reverts free
  requireDenseStreamingExchange();

  auto& load = this->_loadBalancer;

  for (int iC = 0; iC < load.size(); ++iC) {
    _block[iC]->revertStreaming();
  }
}

template<typename T, typename DESCRIPTOR>
template<typename STAGE>
void SuperLattice<T,DESCRIPTOR>::executePostProcessors(STAGE stage)
{
  // Only sync global queue for non-collide stages
  if constexpr (!std::is_same_v<STAGE,stage::PostCollide>) {
  #ifdef PLATFORM_GPU_CUDA
  gpu::cuda::device::synchronizeOuterCollision();
  #endif

  #ifdef PLATFORM_GPU_HIP
  gpu::hip::device::synchronizeOuterCollision();
  #endif

  #ifdef PLATFORM_GPU_SYCL
    gpu::sycl::device::synchronize();
  #endif
  }
  // Join inner collision before the communicate: stage communicators may
  // gather bulk data (e.g. interpolated BC neighborhoods of overlap >= 2)
  // which the transfer streams do not order against otherwise. Post
  // processors read bulk cells as well, so this is needed every stage.
  #ifdef PLATFORM_GPU_CUDA
  gpu::cuda::device::synchronizeInnerCollision();
  #endif
  #ifdef PLATFORM_GPU_HIP
  gpu::hip::device::synchronizeInnerCollision();
  #endif
  if constexpr (!std::is_same_v<STAGE,stage::PostCollide>) {
    getCommunicator(stage).communicate();
  }

  if constexpr (std::is_same_v<STAGE,stage::PostCollide>) {
  #ifdef PLATFORM_GPU_SYCL
    gpu::sycl::device::synchronizeInnerCollision();
  #endif
  }
  auto& load = this->_loadBalancer;

  #ifdef PARALLEL_MODE_OMP
  #pragma omp taskloop
  #endif
  for (int iC = 0; iC < load.size(); ++iC) {
    _block[iC]->template postProcess<STAGE>();
  }

  if constexpr (std::is_same_v<STAGE, stage::PostCollide>) {
  #ifdef PLATFORM_GPU_CUDA
    gpu::cuda::device::synchronizePostProcessor();
  #endif
  #ifdef PLATFORM_GPU_HIP
    gpu::hip::device::synchronizePostProcessor();
  #endif
  }
}


template<typename T, typename DESCRIPTOR>
template<typename STAGE>
SuperCommunicator<T,SuperLattice<T,DESCRIPTOR>>& SuperLattice<T,DESCRIPTOR>::getCommunicator(STAGE stage)
{
  auto iter = _communicator.find(typeid(STAGE));
  if (iter == _communicator.end()) {
    iter = std::get<0>(_communicator.emplace(typeid(STAGE),
                                             std::make_unique<SuperCommunicator<T,SuperLattice>>(*this)));
  }
  return *std::get<1>(*iter);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::communicate()
{
  if (_communicationNeeded) {
    getCommunicator(stage::Full()).communicate();
    _communicationNeeded = false;
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::stripeOffDensityOffset(T offset)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->stripeOffDensityOffset(offset);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T,DESCRIPTOR>& SuperLattice<T,DESCRIPTOR>::getStatistics()
{
  if (_statisticsUpdateNeeded) {
    collectStatistics();
  }
  return _statistics;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::setStatisticsOn()
{
  _statisticsEnabled = true;
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->setStatisticsEnabled(true);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::setStatisticsOff()
{
  _statisticsEnabled = false;
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->setStatisticsEnabled(false);
  }
}

template<typename T, typename DESCRIPTOR>
const LatticeStatistics<T,DESCRIPTOR>& SuperLattice<T,DESCRIPTOR>::getStatistics() const
{
  if (_statisticsUpdateNeeded) {
    const_cast<SuperLattice<T,DESCRIPTOR>*>(this)->collectStatistics();
  }
  return _statistics;
}

template<typename T, typename DESCRIPTOR>
template<typename STAGE>
void SuperLattice<T,DESCRIPTOR>::addPostProcessor(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                  PostProcessorPromise<T,DESCRIPTOR>&& promise)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    getBlock(iC).addPostProcessor(typeid(STAGE),
                                  indicator->getBlockIndicatorF(iC),
                                  std::forward<decltype(promise)>(promise));
  }
}

template<typename T, typename DESCRIPTOR>
template<typename STAGE>
void SuperLattice<T,DESCRIPTOR>::addPostProcessor(PostProcessorPromise<T,DESCRIPTOR>&& promise)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    getBlock(iC).addPostProcessor(typeid(STAGE),
                                  std::forward<decltype(promise)>(promise));
  }
}

template<typename T, typename DESCRIPTOR>
template<typename STAGE>
bool SuperLattice<T,DESCRIPTOR>::removePostProcessor(PostProcessorPromise<T,DESCRIPTOR>&& promise)
{
  bool removed = false;
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    removed |= getBlock(iC).removePostProcessor(typeid(STAGE),
                                                std::forward<decltype(promise)>(promise));
  }
  return removed;
}

template<typename T, typename DESCRIPTOR>
template<typename STAGE>
void SuperLattice<T,DESCRIPTOR>::executeCustomTasks(STAGE stage)
{
  for (auto& f : _customTasks[typeid(STAGE)]) {
    f();
  }
}

template<typename T, typename DESCRIPTOR>
SuperLattice<T,DESCRIPTOR>::~SuperLattice()
{
  // Background output tasks capture this lattice; join any still pending
  // before its blocks are destroyed to avoid a use-after-free on the pool thread.
  for (auto& [stage, tasks] : _backgroundTasks) {
    singleton::pool().waitFor(tasks);
  }
}

template<typename T, typename DESCRIPTOR>
template <typename STAGE, typename F>
void SuperLattice<T,DESCRIPTOR>::scheduleBackgroundTask(F&& f)
{
  _backgroundTasks[typeid(STAGE)].emplace_back(singleton::pool().schedule(f));
}

template<typename T, typename DESCRIPTOR>
template <typename F>
void SuperLattice<T,DESCRIPTOR>::scheduleBackgroundOutput(F&& f)
{
  auto& load = this->getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    if (isPlatformCPU(load.platform(iC))) {
      scheduleBackgroundTask<stage::PreCollide>(std::bind(f, iC));
    }
  }
  for (int iC=0; iC < load.size(); ++iC) {
    if (   load.platform(iC) == Platform::GPU_CUDA
        || load.platform(iC) == Platform::GPU_HIP
        || load.platform(iC) == Platform::GPU_SYCL) {
      scheduleBackgroundTask<stage::PreContextSwitchTo<ProcessingContext::Evaluation>>(std::bind(f, iC));
    }
  }
}

template<typename T, typename DESCRIPTOR>
template <typename CONTEXT>
void SuperLattice<T,DESCRIPTOR>::scheduleBackgroundOutputVTK(CONTEXT&& vtkContext)
{
  vtkContext([](auto& writer, std::size_t iT) {
    writer.writePVD(iT);
  });
  scheduleBackgroundOutput([vtkContext](int iC) {
    vtkContext([&](auto& writer, std::size_t iT) {
      writer.writeVTI(iT, iC);
    });
  });
}

template<typename T, typename DESCRIPTOR>
template <typename STAGE>
void SuperLattice<T,DESCRIPTOR>::waitForBackgroundTasks(STAGE)
{
  if (!_backgroundTasks[typeid(STAGE)].empty()) {
    singleton::pool().waitFor(_backgroundTasks[typeid(STAGE)]);
    _backgroundTasks[typeid(STAGE)].clear();
  }
}

template<typename T, typename DESCRIPTOR>
std::size_t SuperLattice<T,DESCRIPTOR>::getNblock() const
{
  return std::accumulate(_block.begin(), _block.end(), size_t(0), [](std::size_t sum, auto& b) -> std::size_t {
    return sum + b->getNblock();
  });
}


template<typename T, typename DESCRIPTOR>
std::size_t SuperLattice<T,DESCRIPTOR>::getSerializableSize() const
{
  return std::accumulate(_block.begin(), _block.end(), size_t(0), [](std::size_t sum, auto& b) -> std::size_t {
    return sum + b->getSerializableSize();
  });
}

template<typename T, typename DESCRIPTOR>
bool* SuperLattice<T,DESCRIPTOR>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  for (std::size_t iC=0; iC < _block.size(); ++iC) {
    registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, getBlock(iC), loadingMode);
  }

  return dataPtr;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::postLoad()
{
  for (int iC=0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->postLoad();
  }
}

}

#endif
