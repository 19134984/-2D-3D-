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

#ifndef PLATFORM_GPU_SYCL_COMMUNICATOR_HH
#define PLATFORM_GPU_SYCL_COMMUNICATOR_HH

#include "communication/mpiRequest.h"
#include "communication/communicatable.h"
#include "communication/overlapBufferLayout.h"
#include "communication/directionalOverlap.h"
#include "communication/superCommunicationTagCoordinator.h"

#include "registry.h"
#include "context.hh"
#include "column.h"

#ifdef PARALLEL_MODE_MPI
#include "mpi.h"
#if defined(OPEN_MPI) && OPEN_MPI
#include <mpi-ext.h>
#endif
#endif

#include <limits>
#include <thread>
#include <vector>

namespace olb {

template <>
void checkPlatform<Platform::GPU_SYCL>()
{
  OstreamManager clout(std::cout, "GPU_SYCL");

  int nDevices = gpu::sycl::device::getCount();

#ifdef OLB_DEBUG
  clout.setMultiOutput(true);
  if (nDevices < 1) {
    clout << "No SYCL device found" << std::endl;
  } else if (nDevices > 1) {
    clout << "Found " << nDevices << " SYCL devices but only one can be used per MPI process." << std::endl;
  }
  auto& queue = gpu::sycl::device::getQueue();
  clout << queue.get_device().get_info<sycl::info::device::name>() << " visible" << std::endl;
  clout.setMultiOutput(false);
#endif

#ifdef OLB_PARALLEL_MODE_MPI_INJECT_ENTROPY
  clout << "MPI Entropy Injection (Incast Mitigation): ACTIVE" << std::endl;
#endif
}

namespace gpu {

namespace sycl {

#ifdef OLB_GPU_SYCL_SPARSE_CANARY
/// Read and reset the scratch-write counter and fail loudly if any inbound halo cell was
/// scattered into the discard scratch tile (ensureAllocated missed a tile). Called from the
/// scatter wrapper so no device symbol crosses the codegen.
inline void checkSparseScratchCanary() {
  auto& queue = device::getQueue();
  unsigned long long n = 0;
  queue.copy(sparse_scratch_writes, &n).wait();
  if (n > 0) {
    const unsigned long long zero = 0;
    queue.copy(&zero, sparse_scratch_writes).wait();
    throw std::runtime_error("sparse halo scatter hit the scratch tile "
                             + std::to_string(n) + " time(s): an inbound halo tile was not pre-allocated");
  }
}
#endif

/// Blocking gather of FIELD at given indices into buffer
template <typename FIELD, typename CONTEXT>
void gather_field(CONTEXT& lattice, gpu::sycl::Column<CellID>& indices, std::uint8_t* buffer) {
  auto& queue = device::getQueue();
  std::size_t nIndices = indices.size();
  const CellID* idx_ptr = indices.deviceData();
  auto* buf_ptr = reinterpret_cast<typename FIELD::template value_type<typename CONTEXT::value_t>*>(buffer);

  queue.submit([&](::sycl::handler& cgh) {
    cgh.parallel_for(::sycl::range<1>(nIndices), [=](::sycl::id<1> idx) {
      std::size_t iIndex = idx[0];
      auto* field = lattice.template getField<FIELD>();
      for (unsigned iD=0; iD < CONTEXT::descriptor_t::template size<FIELD>(); ++iD) {
        buf_ptr[iD*nIndices+iIndex] = field[iD][idx_ptr[iIndex]];
      }
    });
  }).wait_and_throw();
}

/// Producers a halo gather or source-side copy of `source` must order after
/**
 * Two sources, uniform for every field: the outer collision (recorded on its
 * own OuterCollision queue, so captured explicitly) and all prior work on the
 * default queue (post-processors, couplings, defineField). The default-lane
 * barrier is REQUIRED because the SYCL queues are out of order (see device.hh),
 * so a bare submit onto the shared comm queue would race its producers. The
 * whole-queue barrier is coarse, but a single one per send/copy keeps it cheap.
 **/
template <typename BLOCK>
std::vector<::sycl::event> gatherDependencies(BLOCK& source) {
  std::vector<::sycl::event> dependencies;
  if (source.outerCollisionDone().impl().event) {
    dependencies.push_back(*source.outerCollisionDone().impl().event);
  }
  dependencies.push_back(device::getQueue().ext_oneapi_submit_barrier());
  return dependencies;
}

/// Non-blocking gather of FIELD at given indices into buffer
template <typename FIELD, typename CONTEXT>
::sycl::event async_gather_field(::sycl::queue* stream,
                                 CONTEXT& lattice,
                                 gpu::sycl::Column<CellID>& indices,
                                 std::uint8_t* buffer,
                                 const std::vector<::sycl::event>& dependencies) {
  std::size_t nIndices = indices.size();
  const CellID* idx_ptr = indices.deviceData();
  auto* buf_ptr = reinterpret_cast<typename FIELD::template value_type<typename CONTEXT::value_t>*>(buffer);
  // Order the gather after the field's producers (outer collision + prior
  // default-queue work), uniformly for any field. Only the read accessor
  // differs: a propagatable field is read shift-aware via the cyclic population
  // accessor, other fields through the dense accessor.
  if constexpr (descriptors::is_propagatable_field<FIELD>::value) {
    return stream->submit([&](::sycl::handler& cgh) {
      cgh.depends_on(dependencies);
      cgh.parallel_for(::sycl::range<1>(nIndices), [=](::sycl::id<1> idx) {
        std::size_t iIndex = idx[0];
        for (unsigned iPop=0; iPop < CONTEXT::descriptor_t::template size<FIELD>(); ++iPop) {
          buf_ptr[iPop*nIndices+iIndex] = lattice.getPopulation(iPop)[idx_ptr[iIndex]];
        }
      });
    });
  } else {
    return stream->submit([&](::sycl::handler& cgh) {
      cgh.depends_on(dependencies);
      cgh.parallel_for(::sycl::range<1>(nIndices), [=](::sycl::id<1> idx) {
        std::size_t iIndex = idx[0];
        auto* field = lattice.template getField<FIELD>();
        for (unsigned iD=0; iD < CONTEXT::descriptor_t::template size<FIELD>(); ++iD) {
          buf_ptr[iD*nIndices+iIndex] = field[iD][idx_ptr[iIndex]];
        }
      });
    });
  }
}

/// Non-blocking copy of FIELD at given indices from sourceLattice to targetLattice
template <typename FIELD, typename SOURCE, typename TARGET>
::sycl::event async_copy_field(::sycl::queue* stream,
                               SOURCE& sourceLattice,
                               TARGET& targetLattice,
                               gpu::sycl::Column<CellID>& sourceIndices,
                               gpu::sycl::Column<CellID>& targetIndices,
                               const std::vector<::sycl::event>& dependencies) {
  std::size_t nIndices = sourceIndices.size();
  const CellID* src_idx_ptr = sourceIndices.deviceData();
  const CellID* dst_idx_ptr = targetIndices.deviceData();

  // Order the copy after the source's producers, uniformly for any field. Only
  // the read accessor differs: a propagatable field is read shift-aware via the
  // cyclic population accessor, other fields through the dense accessor.
  if constexpr (std::is_same_v<descriptors::POPULATION,FIELD>) {
    return stream->submit([&](::sycl::handler& cgh) {
      cgh.depends_on(dependencies);
      cgh.parallel_for(::sycl::range<1>(nIndices), [=](::sycl::id<1> idx) {
        std::size_t iIndex = idx[0];
        for (unsigned iPop=0; iPop < SOURCE::descriptor_t::template size<FIELD>(); ++iPop) {
          targetLattice.getPopulation(iPop)[dst_idx_ptr[iIndex]] = sourceLattice.getPopulation(iPop)[src_idx_ptr[iIndex]];
        }
      });
    });
  } else {
    return stream->submit([&](::sycl::handler& cgh) {
      cgh.depends_on(dependencies);
      cgh.parallel_for(::sycl::range<1>(nIndices), [=](::sycl::id<1> idx) {
        std::size_t iIndex = idx[0];
        auto* source = sourceLattice.template getField<FIELD>();
        auto* target = targetLattice.template getField<FIELD>();
        for (unsigned iD=0; iD < SOURCE::descriptor_t::template size<FIELD>(); ++iD) {
          target[iD][dst_idx_ptr[iIndex]] = source[iD][src_idx_ptr[iIndex]];
        }
      });
    });
  }
}


/// Non-blocking gather of fields at given indices into buffer
::sycl::event async_gather_any_fields(::sycl::queue* stream,
                                      gpu::sycl::Column<AnyDeviceFieldArrayD>& fields,
                                      gpu::sycl::Column<CellID>& indices,
                                      std::uint8_t* buffer,
                                      const std::vector<::sycl::event>& dependencies) {
  std::size_t nIndices = indices.size();
  std::size_t nFields = fields.size();
  const CellID* idx_ptr = indices.deviceData();
  const AnyDeviceFieldArrayD* fields_ptr = fields.deviceData();

  return stream->submit([&](::sycl::handler& cgh) {
    cgh.depends_on(dependencies);
    cgh.parallel_for(::sycl::range<1>(nIndices), [=](::sycl::id<1> idx) {
      std::size_t iIndex = idx[0];
      std::uint8_t* buf = buffer; // Local copy to increment
      for (unsigned iField=0; iField < nFields; ++iField) {
        auto& field = fields_ptr[iField];
        for (unsigned iD=0; iD < field.column_count; ++iD) {
          field.gather(iD, idx_ptr[iIndex],
                       buf + (iD*nIndices + iIndex)*field.element_size);
        }
        buf += nIndices*field.column_count*field.element_size;
      }
    });
  });
}

/// Non-blocking copy of fields at given indices from sourceIndices to targetIndices
::sycl::event async_copy_any_fields(::sycl::queue* stream,
                                    gpu::sycl::Column<AnyDeviceFieldArrayD>& sourceFields,
                                    gpu::sycl::Column<AnyDeviceFieldArrayD>& targetFields,
                                    gpu::sycl::Column<CellID>& sourceIndices,
                                    gpu::sycl::Column<CellID>& targetIndices,
                                    const std::vector<::sycl::event>& dependencies) {
  std::size_t nIndices = sourceIndices.size();
  std::size_t nFields = sourceFields.size();
  const CellID* src_idx = sourceIndices.deviceData();
  const CellID* dst_idx = targetIndices.deviceData();
  const AnyDeviceFieldArrayD* src_fields = sourceFields.deviceData();
  const AnyDeviceFieldArrayD* dst_fields = targetFields.deviceData();

  return stream->submit([&](::sycl::handler& cgh) {
    cgh.depends_on(dependencies);
    cgh.parallel_for(::sycl::range<1>(nIndices), [=](::sycl::id<1> idx) {
      std::size_t iIndex = idx[0];
      for (unsigned iField=0; iField < nFields; ++iField) {
        auto& sourceField = src_fields[iField];
        auto& targetField = dst_fields[iField];
        for (unsigned iD=0; iD < sourceField.column_count; ++iD) {
          std::uint8_t scratch[sizeof(double)*8];
          sourceField.gather(iD, src_idx[iIndex], scratch);
          targetField.scatter(iD, dst_idx[iIndex], scratch);
        }
      }
    });
  });
}

/// Non-blocking scatter of FIELD data in buffer to given indices
template <typename FIELD, typename CONTEXT>
::sycl::event async_scatter_field(::sycl::queue* stream,
                                  CONTEXT& lattice,
                                  gpu::sycl::Column<CellID>& indices,
                                  std::uint8_t* buffer) {
  std::size_t nIndices = indices.size();
  const CellID* idx_ptr = indices.deviceData();
  auto* buf_ptr = reinterpret_cast<typename FIELD::template value_type<typename CONTEXT::value_t>*>(buffer);

  // No producer dependency: the scatter writes disjoint padding cells and its
  // ordering against the received payload is enforced host-side (the receive
  // completes before poll() submits the scatter, and communicate() joins on the
  // scatter event before returning).
  if constexpr (std::is_same_v<descriptors::POPULATION,FIELD>) {
    return stream->submit([&](::sycl::handler& cgh) {
      cgh.parallel_for(::sycl::range<1>(nIndices), [=](::sycl::id<1> idx) {
        std::size_t iIndex = idx[0];
        for (unsigned iPop=0; iPop < CONTEXT::descriptor_t::template size<FIELD>(); ++iPop) {
          lattice.getPopulation(iPop)[idx_ptr[iIndex]] = buf_ptr[iPop*nIndices+iIndex];
        }
      });
    });
  } else {
    return stream->submit([&](::sycl::handler& cgh) {
      cgh.parallel_for(::sycl::range<1>(nIndices), [=](::sycl::id<1> idx) {
        std::size_t iIndex = idx[0];
        auto* field = lattice.template getField<FIELD>();
        for (unsigned iD=0; iD < CONTEXT::descriptor_t::template size<FIELD>(); ++iD) {
          field[iD][idx_ptr[iIndex]] = buf_ptr[iD*nIndices+iIndex];
        }
      });
    });
  }
}

/// Non-blocking scatter of fields data in buffer to given indices
::sycl::event async_scatter_any_fields(::sycl::queue* stream,
                                       gpu::sycl::Column<AnyDeviceFieldArrayD>& fields,
                                       gpu::sycl::Column<CellID>& indices,
                                       std::uint8_t* buffer) {
  std::size_t nIndices = indices.size();
  std::size_t nFields = fields.size();
  const CellID* idx_ptr = indices.deviceData();
  const AnyDeviceFieldArrayD* fields_ptr = fields.deviceData();

  auto event = stream->submit([&](::sycl::handler& cgh) {
    cgh.parallel_for(::sycl::range<1>(nIndices), [=](::sycl::id<1> idx) {
      std::size_t iIndex = idx[0];
      std::uint8_t* buf = buffer;
      for (unsigned iField=0; iField < nFields; ++iField) {
        auto& field = fields_ptr[iField];
        for (unsigned iD=0; iD < field.column_count; ++iD) {
          field.scatter(iD, idx_ptr[iIndex],
                        buf + (iD*nIndices + iIndex)*field.element_size);
        }
        buf += nIndices*field.column_count*field.element_size;
      }
    });
  });
#ifdef OLB_GPU_SYCL_SPARSE_CANARY
  // Serialize on the canary build so the scratch-write counter reflects this scatter
  event.wait();
  checkSparseScratchCanary();
#endif
  return event;
}

/// Device-resident group tables of an overlap::MaskGroupLayout
/**
 * SYCL counterpart of gpu::cuda::MaskedSegmentTables: flattened cell-index
 * and component lists over gpu::sycl::Column with per-group host-side
 * offsets driving one kernel submission per group.
 **/
struct MaskedSegmentTables {
  gpu::sycl::Column<std::uint32_t> cellIndices;
  gpu::sycl::Column<std::uint8_t>  components;
  std::vector<std::size_t> cellBegin;
  std::vector<std::size_t> componentBegin;

  static std::size_t totalCells(const overlap::MaskGroupLayout& layout) {
    std::size_t n = 0;
    for (const auto& group : layout.groups) {
      n += group.cellIndices.size();
    }
    return n;
  }
  static std::size_t totalComponents(const overlap::MaskGroupLayout& layout) {
    std::size_t n = 0;
    for (const auto& group : layout.groups) {
      n += group.components.size();
    }
    return n;
  }

  explicit MaskedSegmentTables(const overlap::MaskGroupLayout& layout):
    cellIndices(totalCells(layout)),
    components(totalComponents(layout))
  {
    std::vector<std::uint32_t> flatCells;
    std::vector<std::uint8_t>  flatComponents;
    for (const auto& group : layout.groups) {
      cellBegin.emplace_back(flatCells.size());
      componentBegin.emplace_back(flatComponents.size());
      flatCells.insert(flatCells.end(), group.cellIndices.begin(), group.cellIndices.end());
      flatComponents.insert(flatComponents.end(), group.components.begin(), group.components.end());
    }
    auto& q = device::getQueue();
    if (!flatCells.empty()) {
      q.memcpy(cellIndices.deviceData(), flatCells.data(), flatCells.size()*sizeof(std::uint32_t));
    }
    if (!flatComponents.empty()) {
      q.memcpy(components.deviceData(), flatComponents.data(), flatComponents.size()*sizeof(std::uint8_t));
    }
    q.wait();
  }
};

/// Per-field MaskGroupLayouts resp. device group tables of one exchange
using FieldMaskLayouts = std::map<std::type_index, overlap::MaskGroupLayout>;
using FieldMaskTables  = std::map<std::type_index, MaskedSegmentTables>;

/// Compute the MaskGroupLayout of every masked field participating in `fields`
template <typename T, typename DESCRIPTOR>
FieldMaskLayouts maskLayoutsFor(const std::vector<std::type_index>& fields,
                                const overlap::FieldMasks& masks)
{
  FieldMaskLayouts layouts;
  for (const auto& [field, mask] : masks) {
    if (std::find(fields.begin(), fields.end(), field) != fields.end()) {
      if (const auto traits = overlap::maskedFieldTraits<T,DESCRIPTOR>(field)) {
        layouts.emplace(field, overlap::MaskGroupLayout::compute(
          mask, traits->elementSize, traits->componentsPerDirection));
      }
    }
  }
  return layouts;
}

/// Build the device group tables of every masked field
inline FieldMaskTables maskTablesFor(const FieldMaskLayouts& layouts)
{
  FieldMaskTables tables;
  for (const auto& [field, layout] : layouts) {
    tables.emplace(std::piecewise_construct,
                   std::forward_as_tuple(field),
                   std::forward_as_tuple(layout));
  }
  return tables;
}

/// Build the launch-time device accessor bundle of FIELD's columns
template <typename FIELD, typename BLOCK>
auto propagatableAccessors(BLOCK& block) {
  constexpr unsigned Q = BLOCK::descriptor_t::template size<FIELD>();
  using acc_t = std::decay_t<decltype(block.template getField<FIELD>()[0].getDeviceAccessor())>;
  return Vector<acc_t, Q>([&block](unsigned iD) {
    return block.template getField<FIELD>()[iD].getDeviceAccessor();
  });
}

/// Per-group masked gather of one FIELD (see MaskGroupLayout)
/**
 * The field's components are read through the shift-aware cyclic accessor
 * bundle (built launch-time from getField<FIELD>()[iD].getDeviceAccessor()),
 * so POPULATION and any further propagatable field share one code path.
 **/
template <typename FIELD, typename BLOCK>
std::vector<::sycl::event> async_gather_propagatable_masked(
  ::sycl::queue* stream,
  BLOCK& block,
  gpu::sycl::Column<CellID>& cells,
  const overlap::MaskGroupLayout& layout,
  const MaskedSegmentTables& tables,
  std::uint8_t* buffer,
  const std::vector<::sycl::event>& dependencies)
{
  using T = typename BLOCK::value_t;
  using S = ::olb::detail::field_scalar_t<T, typename BLOCK::descriptor_t, FIELD>;
  auto field = propagatableAccessors<FIELD>(block);
  const CellID* cell_ptr = cells.deviceData();
  std::vector<::sycl::event> events;
  for (unsigned iG=0; iG < layout.groups.size(); ++iG) {
    const auto& group = layout.groups[iG];
    const std::size_t nCells = group.cellIndices.size();
    const unsigned nComponents = group.components.size();
    const std::uint32_t* cellIndices = tables.cellIndices.deviceData() + tables.cellBegin[iG];
    const std::uint8_t*  components  = tables.components.deviceData() + tables.componentBegin[iG];
    // Buffer packs at the field's stored scalar, not sizeof(T)
    auto* buf_ptr = reinterpret_cast<S*>(buffer + group.begin);
    events.push_back(stream->submit([&](::sycl::handler& cgh) {
      cgh.depends_on(dependencies);
      cgh.parallel_for(::sycl::range<1>(nCells), [=](::sycl::id<1> idx) {
        std::size_t iIndex = idx[0];
        const CellID cell = cell_ptr[cellIndices[iIndex]];
        for (unsigned iC=0; iC < nComponents; ++iC) {
          buf_ptr[iC*nCells + iIndex] = field[components[iC]][cell];
        }
      });
    }));
  }
  return events;
}

/// Per-group masked scatter of one FIELD (see MaskGroupLayout)
template <typename FIELD, typename BLOCK>
std::vector<::sycl::event> async_scatter_propagatable_masked(
  ::sycl::queue* stream,
  BLOCK& block,
  gpu::sycl::Column<CellID>& cells,
  const overlap::MaskGroupLayout& layout,
  const MaskedSegmentTables& tables,
  const std::uint8_t* buffer)
{
  using T = typename BLOCK::value_t;
  using S = ::olb::detail::field_scalar_t<T, typename BLOCK::descriptor_t, FIELD>;
  auto field = propagatableAccessors<FIELD>(block);
  // No producer dependency: scatter ordering is host-enforced (see async_scatter_field).
  const CellID* cell_ptr = cells.deviceData();
  std::vector<::sycl::event> events;
  for (unsigned iG=0; iG < layout.groups.size(); ++iG) {
    const auto& group = layout.groups[iG];
    const std::size_t nCells = group.cellIndices.size();
    const unsigned nComponents = group.components.size();
    const std::uint32_t* cellIndices = tables.cellIndices.deviceData() + tables.cellBegin[iG];
    const std::uint8_t*  components  = tables.components.deviceData() + tables.componentBegin[iG];
    auto* buf_ptr = reinterpret_cast<const S*>(buffer + group.begin);
    events.push_back(stream->submit([&](::sycl::handler& cgh) {
      cgh.parallel_for(::sycl::range<1>(nCells), [=](::sycl::id<1> idx) {
        std::size_t iIndex = idx[0];
        const CellID cell = cell_ptr[cellIndices[iIndex]];
        for (unsigned iC=0; iC < nComponents; ++iC) {
          field[components[iC]][cell] = buf_ptr[iC*nCells + iIndex];
        }
      });
    }));
  }
  return events;
}

/// Per-group masked lattice-to-lattice copy of one propagatable FIELD
template <typename FIELD, typename BLOCK>
std::vector<::sycl::event> async_copy_propagatable_masked(
  ::sycl::queue* stream,
  BLOCK& source, BLOCK& target,
  gpu::sycl::Column<CellID>& sourceCells,
  gpu::sycl::Column<CellID>& targetCells,
  const overlap::MaskGroupLayout& layout,
  const MaskedSegmentTables& tables,
  const std::vector<::sycl::event>& dependencies)
{
  auto sourceField = propagatableAccessors<FIELD>(source);
  auto targetField = propagatableAccessors<FIELD>(target);
  const CellID* src_ptr = sourceCells.deviceData();
  const CellID* dst_ptr = targetCells.deviceData();
  std::vector<::sycl::event> events;
  for (unsigned iG=0; iG < layout.groups.size(); ++iG) {
    const auto& group = layout.groups[iG];
    const std::size_t nCells = group.cellIndices.size();
    const unsigned nComponents = group.components.size();
    const std::uint32_t* cellIndices = tables.cellIndices.deviceData() + tables.cellBegin[iG];
    const std::uint8_t*  components  = tables.components.deviceData() + tables.componentBegin[iG];
    events.push_back(stream->submit([&](::sycl::handler& cgh) {
      cgh.depends_on(dependencies);
      cgh.parallel_for(::sycl::range<1>(nCells), [=](::sycl::id<1> idx) {
        std::size_t iIndex = idx[0];
        const std::uint32_t i = cellIndices[iIndex];
        for (unsigned iC=0; iC < nComponents; ++iC) {
          targetField[components[iC]][dst_ptr[i]]
            = sourceField[components[iC]][src_ptr[i]];
        }
      });
    }));
  }
  return events;
}


/// Fold a set of per-group events into a single completion token
inline ::sycl::event foldEvents(::sycl::queue& queue, const std::vector<::sycl::event>& events)
{
  if (events.empty()) {
    return ::sycl::event{};
  }
  if (events.size() == 1) {
    return events.front();
  }
  return queue.ext_oneapi_submit_barrier(events);
}

/// Per-segment dense field arrays over gpu::sycl::Column, empty otherwise
using SegmentDenseArrays = std::vector<gpu::sycl::Column<AnyDeviceFieldArrayD>>;

/// Device field arrays for `block`'s AnyField segments, empty otherwise
/**
 * Only the type-erased (AnyField) segments capture the registry's device
 * pointers. The typed masked and dense segments build fresh launch-time
 * accessors. SYCL captures the arrays once (no per-step refresh).
 **/
template <typename BLOCK>
SegmentDenseArrays segmentDenseArrays(BLOCK& block,
                                      const overlap::MultiFieldPartition& partition)
{
  auto& q = device::getQueue();
  SegmentDenseArrays denseArrays;
  denseArrays.reserve(partition.segments.size());
  for (const auto& seg : partition.segments) {
    if (seg.kind == overlap::FieldSegment::Kind::AnyField) {
      auto hostArrays = block.getDataRegistry().deviceFieldArrays({seg.field});
      denseArrays.emplace_back(hostArrays.size());
      if (!hostArrays.empty()) {
        // Blocking: hostArrays is a loop-local, so the copy must finish
        // before it is destroyed at the end of this iteration.
        q.memcpy(denseArrays.back().deviceData(), hostArrays.data(),
                 hostArrays.size()*sizeof(AnyDeviceFieldArrayD)).wait();
      }
    } else {
      denseArrays.emplace_back(std::size_t{0});
    }
  }
  return denseArrays;
}

/// Build the packed multi-field partition and per-segment dense field arrays
/**
 * Segments are classified in request order: a masked propagatable field
 * takes its own MaskGroupLayout segment, POPULATION in a dense exchange
 * takes the typed dense kernel, every other field the type-erased path.
 **/
template <typename BLOCK>
overlap::MultiFieldPartition makeSegmentPartition(
  BLOCK& block,
  const std::vector<std::type_index>& fields,
  const std::vector<CellID>& cells,
  const FieldMaskLayouts& layouts,
  SegmentDenseArrays& denseArrays)
{
  using T = typename BLOCK::value_t;
  auto partition = overlap::partitionMultiField(fields, alignof(T),
    [&](std::type_index f) {
      if (layouts.find(f) != layouts.end()) {
        return overlap::FieldSegment::Kind::MaskedPropagatable;
      }
      if (f == std::type_index(typeid(descriptors::POPULATION))) {
        return overlap::FieldSegment::Kind::TypedDense;
      }
      return overlap::FieldSegment::Kind::AnyField;
    },
    [&](std::type_index f) -> std::size_t {
      if (auto layout = layouts.find(f); layout != layouts.end()) {
        return layout->second.size;
      }
      return block.getCommunicatable(f).size(cells);
    });
  denseArrays = segmentDenseArrays(block, partition);
  return partition;
}

/// Dispatch the matching FIELD of a masked segment to `launch`
template <typename BLOCK, typename LAUNCH>
void dispatchMaskedField(std::type_index field, LAUNCH&& launch) {
  BLOCK::descriptor_t::fields_t::for_each([&](auto f) {
    using FIELD = typename decltype(f)::type;
    constexpr bool propagatable = descriptors::is_propagatable_field<FIELD>::value;
    if constexpr (propagatable || BLOCK::descriptor_t::template size<FIELD>() <= 64) {
      if (field == std::type_index(typeid(FIELD))) {
        launch(f);
      }
    }
  });
}

/// Gather every segment of a packed multi-field buffer (see makeSegmentPartition)
template <typename BLOCK>
std::vector<::sycl::event> async_gather_segments(
  ::sycl::queue* stream, BLOCK& block,
  const overlap::MultiFieldPartition& partition,
  SegmentDenseArrays& denseArrays,
  gpu::sycl::Column<CellID>& cells,
  const FieldMaskLayouts& layouts,
  const FieldMaskTables& tables,
  std::uint8_t* buffer,
  const std::vector<::sycl::event>& dependencies)
{
  using T = typename BLOCK::value_t;
  using DESCRIPTOR = typename BLOCK::descriptor_t;
  std::vector<::sycl::event> events;
  for (std::size_t i=0; i < partition.segments.size(); ++i) {
    const auto& seg = partition.segments[i];
    std::uint8_t* at = buffer + seg.begin;
    switch (seg.kind) {
    case overlap::FieldSegment::Kind::MaskedPropagatable:
      dispatchMaskedField<BLOCK>(seg.field, [&](auto id) {
        auto ge = async_gather_propagatable_masked<typename decltype(id)::type>(
          stream, block, cells, layouts.at(seg.field), tables.at(seg.field), at, dependencies);
        events.insert(events.end(), ge.begin(), ge.end());
      });
      break;
    case overlap::FieldSegment::Kind::TypedDense: {
      DeviceContext<T,DESCRIPTOR> lattice(block);
      events.push_back(async_gather_field<descriptors::POPULATION>(stream, lattice, cells, at, dependencies));
      break;
    }
    case overlap::FieldSegment::Kind::AnyField:
      events.push_back(async_gather_any_fields(stream, denseArrays[i], cells, at, dependencies));
      break;
    }
  }
  return events;
}

/// Scatter every segment of a packed multi-field buffer (see makeSegmentPartition)
template <typename BLOCK>
std::vector<::sycl::event> async_scatter_segments(
  ::sycl::queue* stream, BLOCK& block,
  const overlap::MultiFieldPartition& partition,
  SegmentDenseArrays& denseArrays,
  gpu::sycl::Column<CellID>& cells,
  const FieldMaskLayouts& layouts,
  const FieldMaskTables& tables,
  std::uint8_t* buffer)
{
  using T = typename BLOCK::value_t;
  using DESCRIPTOR = typename BLOCK::descriptor_t;
  std::vector<::sycl::event> events;
  for (std::size_t i=0; i < partition.segments.size(); ++i) {
    const auto& seg = partition.segments[i];
    std::uint8_t* at = buffer + seg.begin;
    switch (seg.kind) {
    case overlap::FieldSegment::Kind::MaskedPropagatable:
      dispatchMaskedField<BLOCK>(seg.field, [&](auto id) {
        auto groupEvents = async_scatter_propagatable_masked<typename decltype(id)::type>(
          stream, block, cells, layouts.at(seg.field), tables.at(seg.field), at);
        events.insert(events.end(), groupEvents.begin(), groupEvents.end());
      });
      break;
    case overlap::FieldSegment::Kind::TypedDense: {
      DeviceContext<T,DESCRIPTOR> lattice(block);
      events.push_back(async_scatter_field<descriptors::POPULATION>(stream, lattice, cells, at));
      break;
    }
    case overlap::FieldSegment::Kind::AnyField:
      events.push_back(async_scatter_any_fields(stream, denseArrays[i], cells, at));
      break;
    }
  }
  return events;
}

/// Copy every segment of a multi-field exchange source->target (no packed buffer)
template <typename BLOCK>
std::vector<::sycl::event> async_copy_segments(
  ::sycl::queue* stream, BLOCK& source, BLOCK& target,
  const overlap::MultiFieldPartition& partition,
  SegmentDenseArrays& sourceDense, SegmentDenseArrays& targetDense,
  gpu::sycl::Column<CellID>& sourceCells, gpu::sycl::Column<CellID>& targetCells,
  const FieldMaskLayouts& layouts,
  const FieldMaskTables& tables,
  const std::vector<::sycl::event>& dependencies)
{
  using T = typename BLOCK::value_t;
  using DESCRIPTOR = typename BLOCK::descriptor_t;
  std::vector<::sycl::event> events;
  for (std::size_t i=0; i < partition.segments.size(); ++i) {
    const auto& seg = partition.segments[i];
    switch (seg.kind) {
    case overlap::FieldSegment::Kind::MaskedPropagatable:
      dispatchMaskedField<BLOCK>(seg.field, [&](auto id) {
        auto groupEvents = async_copy_propagatable_masked<typename decltype(id)::type>(
          stream, source, target, sourceCells, targetCells, layouts.at(seg.field), tables.at(seg.field), dependencies);
        events.insert(events.end(), groupEvents.begin(), groupEvents.end());
      });
      break;
    case overlap::FieldSegment::Kind::TypedDense: {
      DeviceContext<T,DESCRIPTOR> sourceLattice(source);
      DeviceContext<T,DESCRIPTOR> targetLattice(target);
      events.push_back(async_copy_field<descriptors::POPULATION>(
        stream, sourceLattice, targetLattice, sourceCells, targetCells, dependencies));
      break;
    }
    case overlap::FieldSegment::Kind::AnyField:
      events.push_back(async_copy_any_fields(stream, sourceDense[i], targetDense[i],
                                             sourceCells, targetCells, dependencies));
      break;
    }
  }
  return events;
}

}

}

}

#ifdef FEATURE_DIRECTIONAL_CANARY
#include "canary.h"
#endif

namespace olb {

#ifdef PARALLEL_MODE_MPI

template <typename T, typename DESCRIPTOR>
class ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>::SendTask {
private:
  const bool _onlyPopulationField;

  /// Negotiated per-field, per-cell component masks restrict each masked
  /// field's segment to its overlap::MaskGroupLayout (e.g. directional
  /// streaming). Empty masks keep the dense layout.
  const gpu::sycl::FieldMaskLayouts _maskLayouts;
  const gpu::sycl::FieldMaskTables _maskTables;
  const bool _masked;

  /// Every field becomes one packed segment in request order (see
  /// makeSegmentPartition): masked propagatable fields take their own
  /// MaskGroupLayout, POPULATION in a dense exchange the typed kernel, the
  /// rest the type-erased path. The byte layout matches the CPU-side
  /// (Masked)MultiConcreteCommunicatable.
  gpu::sycl::SegmentDenseArrays _denseArrays;
  const overlap::MultiFieldPartition _segments;

  gpu::sycl::Column<CellID> _cells;
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& _source;
  gpu::sycl::device::unique_ptr<std::uint8_t> _buffer;

  std::unique_ptr<MpiSendRequest> _request;

  bool _payloadSent;

public:
  SendTask(MPI_Comm comm, int tag, int rank,
           const std::vector<std::type_index>& fields,
           const std::vector<CellID>& cells,
           const overlap::FieldMasks& masks,
           ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& block):
    _onlyPopulationField(fields.size() == 1 && fields[0] == typeid(descriptors::POPULATION)),
    _maskLayouts(gpu::sycl::maskLayoutsFor<T,DESCRIPTOR>(fields, masks)),
    _maskTables(gpu::sycl::maskTablesFor(_maskLayouts)),
    _masked(!_maskLayouts.empty()),
    _segments(gpu::sycl::makeSegmentPartition(block, fields, cells, _maskLayouts, _denseArrays)),
    _cells(cells.size()),
    _source(block),
    _payloadSent(false)
  {
    auto& q = gpu::sycl::device::getQueue();
    q.memcpy(_cells.deviceData(), cells.data(), cells.size()*sizeof(CellID));
    _buffer = gpu::sycl::device::malloc<std::uint8_t>(_segments.size);
    _request = std::make_unique<MpiSendRequest>(_buffer.get(), _segments.size, rank, tag, comm);
    q.wait();
  }

  ~SendTask()
  {
    wait();
  }

  void prepare()
  {
    OLB_PROFILE_T("GPU_SYCL:Comm:SendTask:prepare", DESCRIPTOR);
    _payloadSent = false;
  }

  ::sycl::event submitGather()
  {
    OLB_PROFILE_T("GPU_SYCL:Comm:SendTask:submitGather", DESCRIPTOR);
    auto& q = gpu::sycl::device::getQueue();
    auto dependencies = gpu::sycl::gatherDependencies(_source);
    if (_onlyPopulationField && !_masked) {
      gpu::sycl::DeviceContext<T,DESCRIPTOR> lattice(_source);
      return gpu::sycl::async_gather_field<descriptors::POPULATION>(
        &q, lattice, _cells, _buffer.get(), dependencies);
    }
    auto events = gpu::sycl::async_gather_segments(
      &q, _source, _segments, _denseArrays, _cells, _maskLayouts, _maskTables, _buffer.get(), dependencies);
    return gpu::sycl::foldEvents(q, events);
  }

  bool send()
  {
    if (_payloadSent) {
      return true;
    }
    OLB_PROFILE_T("GPU_SYCL:Comm:SendTask:send", DESCRIPTOR);
    _request->start();
    _payloadSent = true;
    return true;
  }

  void wait()
  {
    OLB_PROFILE_T("GPU_SYCL:Comm:SendTask:wait", DESCRIPTOR);
    _request->wait();
  }
};

template <typename T, typename DESCRIPTOR>
class ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>::RecvTask {
private:
  const int _tag;
  const int _rank;

  const bool _onlyPopulationField;

  /// Same per-field masked handling as SendTask (see there)
  const gpu::sycl::FieldMaskLayouts _maskLayouts;
  const gpu::sycl::FieldMaskTables _maskTables;
  const bool _masked;
#ifdef FEATURE_DIRECTIONAL_CANARY
  /// Per-field, per-cell masks for the directional canary (populated only when enabled)
  gpu::sycl::FieldCanaryMasks _canaryMasks;
#endif

  /// Same packed multi-segment routing as SendTask (see there)
  gpu::sycl::SegmentDenseArrays _denseArrays;
  const overlap::MultiFieldPartition _segments;

  gpu::sycl::Column<CellID> _cells;
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& _target;
  gpu::sycl::device::unique_ptr<std::uint8_t> _buffer;
  std::unique_ptr<MpiRecvRequest> _request;
  ::sycl::event _scatterEvent;
  bool _scatterSubmitted;

public:
  RecvTask(MPI_Comm comm, int tag, int rank,
           const std::vector<std::type_index>& fields,
           const std::vector<CellID>& cells,
           const overlap::FieldMasks& masks,
           ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& block):
    _tag(tag),
    _rank(rank),
    _onlyPopulationField(fields.size() == 1 && fields[0] == typeid(descriptors::POPULATION)),
    _maskLayouts(gpu::sycl::maskLayoutsFor<T,DESCRIPTOR>(fields, masks)),
    _maskTables(gpu::sycl::maskTablesFor(_maskLayouts)),
    _masked(!_maskLayouts.empty()),
#ifdef FEATURE_DIRECTIONAL_CANARY
    _canaryMasks(gpu::sycl::canaryMasksFor(fields, masks)),
#endif
    _segments(gpu::sycl::makeSegmentPartition(block, fields, cells, _maskLayouts, _denseArrays)),
    _cells(cells.size()),
    _target(block),
    _scatterSubmitted(false)
  {
    auto& q = gpu::sycl::device::getQueue();
    q.memcpy(_cells.deviceData(), cells.data(), cells.size()*sizeof(CellID));
    _buffer = gpu::sycl::device::malloc<std::uint8_t>(_segments.size);
    _request = std::make_unique<MpiRecvRequest>(_buffer.get(), _segments.size, _rank, _tag, comm);
    q.wait();
  }

  ~RecvTask()
  {
    wait();
  }

  void receive()
  {
    OLB_PROFILE_T("GPU_SYCL:Comm:RecvTask:receive", DESCRIPTOR);
    _scatterSubmitted = false;
    _request->start();
  }

  bool poll()
  {
    if (_scatterSubmitted) {
      return false;
    }

    if (_request->isDone()) {
      OLB_PROFILE_T("GPU_SYCL:Comm:RecvTask:pollDone", DESCRIPTOR);
      auto& q = gpu::sycl::device::getQueue();
      if (_onlyPopulationField && !_masked) {
        gpu::sycl::DeviceContext<T,DESCRIPTOR> lattice(_target);
        _scatterEvent = gpu::sycl::async_scatter_field<descriptors::POPULATION>(
          &q, lattice, _cells, _buffer.get());
      } else {
        auto events = gpu::sycl::async_scatter_segments(
          &q, _target, _segments, _denseArrays, _cells, _maskLayouts, _maskTables, _buffer.get());
#ifdef FEATURE_DIRECTIONAL_CANARY
        if (!_canaryMasks.empty()) {
          events = gpu::sycl::async_poison_segments(
            &q, _target, _segments, _cells, _canaryMasks,
            overlap::componentsPerDirection<DESCRIPTOR>(), events);
        }
#endif
        _scatterEvent = gpu::sycl::foldEvents(q, events);
      }
      _scatterSubmitted = true;
      return false;
    }
    return true;
  }

  void wait()
  {
    OLB_PROFILE_T("GPU_SYCL:Comm:RecvTask:wait", DESCRIPTOR);
    if (_scatterSubmitted) {
      _scatterEvent.wait();
    }
  }
};

#endif // PARALLEL_MODE_MPI

template <typename T, typename DESCRIPTOR>
class ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>::HomogeneousCopyTask
  : public ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>::CopyTask {
private:
  const bool _onlyPopulationField;

  /// Same per-field masked handling as SendTask (see there)
  const gpu::sycl::FieldMaskLayouts _maskLayouts;
  const gpu::sycl::FieldMaskTables _maskTables;
  const bool _masked;
#ifdef FEATURE_DIRECTIONAL_CANARY
  /// Per-field, per-cell masks for the directional canary (populated only when enabled)
  gpu::sycl::FieldCanaryMasks _canaryMasks;
#endif

  /// Same packed multi-segment routing as SendTask. Direct copies have no
  /// packed buffer, so the segment byte offsets are unused
  gpu::sycl::SegmentDenseArrays _sourceDense;
  const overlap::MultiFieldPartition _segments;
  gpu::sycl::SegmentDenseArrays _targetDense;

  gpu::sycl::Column<CellID> _targetCells;
  gpu::sycl::Column<CellID> _sourceCells;

  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& _target;
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& _source;

  ::sycl::event _copyEvent;
  bool _copySubmitted;

public:
  HomogeneousCopyTask(
    const std::vector<std::type_index>& fields,
    const std::vector<CellID>& targetCells, ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& target,
    const std::vector<CellID>& sourceCells, ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& source,
    const overlap::FieldMasks& masks = {}):
    _onlyPopulationField(fields.size() == 1 && fields[0] == typeid(descriptors::POPULATION)),
    _maskLayouts(gpu::sycl::maskLayoutsFor<T,DESCRIPTOR>(fields, masks)),
    _maskTables(gpu::sycl::maskTablesFor(_maskLayouts)),
    _masked(!_maskLayouts.empty()),
#ifdef FEATURE_DIRECTIONAL_CANARY
    _canaryMasks(gpu::sycl::canaryMasksFor(fields, masks)),
#endif
    _segments(gpu::sycl::makeSegmentPartition(source, fields, sourceCells, _maskLayouts, _sourceDense)),
    _targetDense(gpu::sycl::segmentDenseArrays(target, _segments)),
    _targetCells(targetCells.size()),
    _sourceCells(sourceCells.size()),
    _target(target),
    _source(source),
    _copySubmitted(false)
  {
    OLB_ASSERT(sourceCells.size() == targetCells.size(),
               "Source cell count must match target cell count");

    auto& q = gpu::sycl::device::getQueue();
    q.memcpy(_sourceCells.deviceData(), sourceCells.data(), sourceCells.size()*sizeof(CellID));
    q.memcpy(_targetCells.deviceData(), targetCells.data(), targetCells.size()*sizeof(CellID));

    q.wait();
  }

  ~HomogeneousCopyTask()
  {
    wait();
  }

  void copy() override
  {
    OLB_PROFILE_T("GPU_SYCL:Comm:CopyTask:copy", DESCRIPTOR);
    auto& q = gpu::sycl::device::getQueue();
    auto dependencies = gpu::sycl::gatherDependencies(_source);
    if (_onlyPopulationField && !_masked) {
      gpu::sycl::DeviceContext<T,DESCRIPTOR> sourceLattice(_source);
      gpu::sycl::DeviceContext<T,DESCRIPTOR> targetLattice(_target);
      _copyEvent = gpu::sycl::async_copy_field<descriptors::POPULATION>(
        &q, sourceLattice, targetLattice, _sourceCells, _targetCells, dependencies);
    } else {
      auto events = gpu::sycl::async_copy_segments(
        &q, _source, _target, _segments, _sourceDense, _targetDense,
        _sourceCells, _targetCells, _maskLayouts, _maskTables, dependencies);
#ifdef FEATURE_DIRECTIONAL_CANARY
      if (!_canaryMasks.empty()) {
        events = gpu::sycl::async_poison_segments(
          &q, _target, _segments, _targetCells, _canaryMasks,
          overlap::componentsPerDirection<DESCRIPTOR>(), events);
      }
#endif
      _copyEvent = gpu::sycl::foldEvents(q, events);
    }
    _copySubmitted = true;
  }

  void wait() override
  {
    OLB_PROFILE_T("GPU_SYCL:Comm:CopyTask:wait", DESCRIPTOR);
    if (_copySubmitted) {
      _copyEvent.wait();
      _copySubmitted = false;
    }
  }

};

/// Private implementation of heterogeneous copy task between CPU_* source and GPU_SYCL target
template <typename T, typename DESCRIPTOR, Platform SOURCE>
class HeterogeneousCopyTaskDataForGpuTarget : public ConcreteHeterogeneousCopyTask {
private:
  const bool _onlyPopulationField;

  /// Same per-field masked handling as SendTask (see there). The CPU side
  /// serializes the identical mask-group layouts via MaskedMultiCommunicatable
  const gpu::sycl::FieldMaskLayouts _maskLayouts;
  const gpu::sycl::FieldMaskTables _maskTables;
  const bool _masked;
#ifdef FEATURE_DIRECTIONAL_CANARY
  gpu::sycl::FieldCanaryMasks _canaryMasks;
#endif

  /// Same packed multi-segment routing as SendTask. The CPU-side serialize
  /// packs fields sequentially in request order, matching the segment offsets
  gpu::sycl::SegmentDenseArrays _denseArrays;
  const overlap::MultiFieldPartition _segments;

  gpu::sycl::Column<CellID>  _targetCells;
  const std::vector<CellID>& _sourceCells;

  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& _target;

  std::unique_ptr<Communicatable> _source;

  gpu::sycl::Column<std::uint8_t> _buffer;
  ::sycl::event _scatterEvent;
  bool _scatterSubmitted;

public:
  HeterogeneousCopyTaskDataForGpuTarget(
    const std::vector<std::type_index>& fields,
    const std::vector<CellID>& targetCells, ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& target,
    const std::vector<CellID>& sourceCells, ConcreteBlockLattice<T,DESCRIPTOR,SOURCE>&             source,
    const overlap::FieldMasks& masks = {}):
    _onlyPopulationField(fields.size() == 1 && fields[0] == typeid(descriptors::POPULATION)),
    _maskLayouts(gpu::sycl::maskLayoutsFor<T,DESCRIPTOR>(fields, masks)),
    _maskTables(gpu::sycl::maskTablesFor(_maskLayouts)),
    _masked(!_maskLayouts.empty()),
#ifdef FEATURE_DIRECTIONAL_CANARY
    _canaryMasks(gpu::sycl::canaryMasksFor(fields, masks)),
#endif
    _segments(gpu::sycl::makeSegmentPartition(target, fields, targetCells, _maskLayouts, _denseArrays)),
    _targetCells(targetCells.size()),
    _sourceCells(sourceCells),
    _target(target),
    _source(overlap::makeMaskedMultiCommunicatable(source, fields, masks)),
    _buffer(_source->size(_sourceCells)),
    _scatterSubmitted(false)
  {
    auto& q = gpu::sycl::device::getQueue();
    q.memcpy(_targetCells.deviceData(), targetCells.data(), targetCells.size()*sizeof(CellID));
    q.wait();
  }

  void copy() override {
    auto& q = gpu::sycl::device::getQueue();
    _source->serialize(_sourceCells, _buffer.data());
    _buffer.setProcessingContext(ProcessingContext::Simulation);

    if (_onlyPopulationField && !_masked) {
      gpu::sycl::DeviceContext<T,DESCRIPTOR> lattice(_target);
      _scatterEvent = gpu::sycl::async_scatter_field<descriptors::POPULATION>(
        &q, lattice, _targetCells, _buffer.deviceData());
    } else {
      auto events = gpu::sycl::async_scatter_segments(
        &q, _target, _segments, _denseArrays, _targetCells, _maskLayouts, _maskTables, _buffer.deviceData());
#ifdef FEATURE_DIRECTIONAL_CANARY
      if (!_canaryMasks.empty()) {
        events = gpu::sycl::async_poison_segments(
          &q, _target, _segments, _targetCells, _canaryMasks,
          overlap::componentsPerDirection<DESCRIPTOR>(), events);
      }
#endif
      _scatterEvent = gpu::sycl::foldEvents(q, events);
    }
    _scatterSubmitted = true;
  }

  void wait() override {
    if (_scatterSubmitted) {
      _scatterEvent.wait();
      _scatterSubmitted = false;
    }
  }

};

/// Private implementation of heterogeneous copy task between GPU_SYCL source and CPU_* target
template <typename T, typename DESCRIPTOR, Platform TARGET>
class HeterogeneousCopyTaskDataForGpuSource : public ConcreteHeterogeneousCopyTask {
private:
  const bool _onlyPopulationField;

  /// Same per-field masked handling as SendTask (see there). The CPU side
  /// deserializes the identical mask-group layouts via MaskedMultiCommunicatable
  const gpu::sycl::FieldMaskLayouts _maskLayouts;
  const gpu::sycl::FieldMaskTables _maskTables;
  const bool _masked;

  /// Same packed multi-segment routing as SendTask (see there)
  gpu::sycl::SegmentDenseArrays _denseArrays;
  const overlap::MultiFieldPartition _segments;

  const std::vector<CellID>& _targetCells;
  gpu::sycl::Column<CellID>  _sourceCells;

  std::unique_ptr<Communicatable> _target;
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& _source;

  gpu::sycl::Column<std::uint8_t> _buffer;
  ::sycl::event _gatherEvent;
  bool _gatherSubmitted;

public:
  HeterogeneousCopyTaskDataForGpuSource(
    const std::vector<std::type_index>& fields,
    const std::vector<CellID>& targetCells, ConcreteBlockLattice<T,DESCRIPTOR,TARGET>&             target,
    const std::vector<CellID>& sourceCells, ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& source,
    const overlap::FieldMasks& masks = {}):
    _onlyPopulationField(fields.size() == 1 && fields[0] == typeid(descriptors::POPULATION)),
    _maskLayouts(gpu::sycl::maskLayoutsFor<T,DESCRIPTOR>(fields, masks)),
    _maskTables(gpu::sycl::maskTablesFor(_maskLayouts)),
    _masked(!_maskLayouts.empty()),
    _segments(gpu::sycl::makeSegmentPartition(source, fields, sourceCells, _maskLayouts, _denseArrays)),
    _targetCells(targetCells),
    _sourceCells(sourceCells.size()),
    _target(overlap::makeMaskedMultiCommunicatable(target, fields, masks)),
    _source(source),
    _buffer(_target->size(_targetCells)),
    _gatherSubmitted(false)
  {
    auto& q = gpu::sycl::device::getQueue();
    q.memcpy(_sourceCells.deviceData(), sourceCells.data(), sourceCells.size()*sizeof(CellID));
    q.wait();
  }

  void copy() override {
    auto& q = gpu::sycl::device::getQueue();
    auto dependencies = gpu::sycl::gatherDependencies(_source);
    if (_onlyPopulationField && !_masked) {
      gpu::sycl::DeviceContext<T,DESCRIPTOR> lattice(_source);
      _gatherEvent = gpu::sycl::async_gather_field<descriptors::POPULATION>(
        &q, lattice, _sourceCells, _buffer.deviceData(), dependencies);
    } else {
      auto events = gpu::sycl::async_gather_segments(
        &q, _source, _segments, _denseArrays, _sourceCells, _maskLayouts, _maskTables, _buffer.deviceData(), dependencies);
      _gatherEvent = gpu::sycl::foldEvents(q, events);
    }
    _gatherSubmitted = true;
  }

  void wait() override {
    if (_gatherSubmitted) {
      _gatherEvent.wait();
      _gatherSubmitted = false;
      _buffer.setProcessingContext(ProcessingContext::Evaluation);
      _target->deserialize(_targetCells, _buffer.data());
    }
  }

};

template <typename T, typename DESCRIPTOR, Platform SOURCE>
HeterogeneousCopyTask<T,DESCRIPTOR,SOURCE,Platform::GPU_SYCL>::HeterogeneousCopyTask(
  const std::vector<std::type_index>& fields,
  const std::vector<CellID>& targetCells, ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& target,
  const std::vector<CellID>& sourceCells, ConcreteBlockLattice<T,DESCRIPTOR,SOURCE>&             source,
  const overlap::FieldMasks& masks):
  _impl(new HeterogeneousCopyTaskDataForGpuTarget<T,DESCRIPTOR,SOURCE>(fields,
                                                                       targetCells, target,
                                                                       sourceCells, source,
                                                                       masks))
{
  OLB_ASSERT(sourceCells.size() == targetCells.size(),
             "Source cell count must match target cell count");
}

template <typename T, typename DESCRIPTOR, Platform SOURCE>
void HeterogeneousCopyTask<T,DESCRIPTOR,SOURCE,Platform::GPU_SYCL>::copy()
{
  _impl->copy();
}

template <typename T, typename DESCRIPTOR, Platform SOURCE>
void HeterogeneousCopyTask<T,DESCRIPTOR,SOURCE,Platform::GPU_SYCL>::wait()
{
  _impl->wait();
}

template <typename T, typename DESCRIPTOR, Platform TARGET>
HeterogeneousCopyTask<T,DESCRIPTOR,Platform::GPU_SYCL,TARGET>::HeterogeneousCopyTask(
  const std::vector<std::type_index>& fields,
  const std::vector<CellID>& targetCells, ConcreteBlockLattice<T,DESCRIPTOR,TARGET>&             target,
  const std::vector<CellID>& sourceCells, ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& source,
  const overlap::FieldMasks& masks):
  _impl(new HeterogeneousCopyTaskDataForGpuSource<T,DESCRIPTOR,TARGET>(fields,
                                                                       targetCells, target,
                                                                       sourceCells, source,
                                                                       masks))
{
  OLB_ASSERT(sourceCells.size() == targetCells.size(),
             "Source cell count must match target cell count");
}

template <typename T, typename DESCRIPTOR, Platform TARGET>
void HeterogeneousCopyTask<T,DESCRIPTOR,Platform::GPU_SYCL,TARGET>::copy()
{
  _impl->copy();
}

template <typename T, typename DESCRIPTOR, Platform TARGET>
void HeterogeneousCopyTask<T,DESCRIPTOR,Platform::GPU_SYCL,TARGET>::wait()
{
  _impl->wait();
}

template <typename T, typename DESCRIPTOR>
ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>::ConcreteBlockCommunicator(
  SuperLattice<T,DESCRIPTOR>& super,
  LoadBalancer<T>& loadBalancer,
#ifdef PARALLEL_MODE_MPI
  SuperCommunicationTagCoordinator<T>& tagCoordinator,
  MPI_Comm comm,
#endif
  int iC,
  const BlockCommunicationNeighborhood<T,DESCRIPTOR::d>& neighborhood):
  _iC(iC)
#ifdef PARALLEL_MODE_MPI
, _mpiCommunicator(comm)
#endif
{
  // Materialise inbound halo tiles of any sparse field host-side before the next
  // Simulation sync. The device scatter never allocates (no-op for dense fields).
  // SYCL local pairs are homogeneous, so every neighbour uses the sparse scatter path.
  {
    auto& block = super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>(_iC);
    neighborhood.forNeighbors([&](int remoteC) {
      const auto& inbound = neighborhood.getCellsInboundFrom(remoteC);
      if (!inbound.empty()) {
        for (std::type_index field : neighborhood.getFieldsCommonWith(remoteC)) {
          block.getCommunicatable(field).ensureAllocated(inbound);
        }
      }
    });
  }
#ifdef PARALLEL_MODE_MPI
  neighborhood.forNeighbors([&](int remoteC) {
    if (loadBalancer.isLocal(remoteC)) {
      const Platform remotePlatform = loadBalancer.platform(loadBalancer.loc(remoteC));
      if (!neighborhood.getCellsInboundFrom(remoteC).empty()) {
        switch (remotePlatform) {
          case Platform::GPU_SYCL:
            _copyTasks.emplace_back(new HomogeneousCopyTask(
              neighborhood.getFieldsCommonWith(remoteC),
              neighborhood.getCellsInboundFrom(remoteC),   super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>(_iC),
              neighborhood.getCellsRequestedFrom(remoteC), super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>(loadBalancer.loc(remoteC)),
              neighborhood.getMasksInboundFrom(remoteC)));
            break;
#ifdef PLATFORM_CPU_SIMD
          case Platform::CPU_SIMD:
            // Use manual copy for local GPU-CPU communication due to better performance
            _copyTasks.emplace_back(new HeterogeneousCopyTask<T,DESCRIPTOR,Platform::CPU_SIMD,Platform::GPU_SYCL>(
              neighborhood.getFieldsCommonWith(remoteC),
              neighborhood.getCellsInboundFrom(remoteC),   super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>(_iC),
              neighborhood.getCellsRequestedFrom(remoteC), super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SIMD>>(loadBalancer.loc(remoteC)),
              neighborhood.getMasksInboundFrom(remoteC)));
            break;
#endif
#ifdef PLATFORM_CPU_SISD
          case Platform::CPU_SISD:
            // Use manual copy for local GPU-CPU communication due to better performance
            _copyTasks.emplace_back(new HeterogeneousCopyTask<T,DESCRIPTOR,Platform::CPU_SISD,Platform::GPU_SYCL>(
              neighborhood.getFieldsCommonWith(remoteC),
              neighborhood.getCellsInboundFrom(remoteC),   super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>(_iC),
              neighborhood.getCellsRequestedFrom(remoteC), super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>>(loadBalancer.loc(remoteC)),
              neighborhood.getMasksInboundFrom(remoteC)));
            break;
#endif
          default:
            throw std::runtime_error("Invalid remote PLATFORM");
        }
      }
    } else {
      if (!neighborhood.getCellsOutboundTo(remoteC).empty()) {
        _sendTasks.emplace_back(std::make_unique<SendTask>(
          _mpiCommunicator, tagCoordinator.get(loadBalancer.glob(_iC), remoteC),
          loadBalancer.rank(remoteC),
          neighborhood.getFieldsCommonWith(remoteC),
          neighborhood.getCellsOutboundTo(remoteC),
          neighborhood.getMasksOutboundTo(remoteC),
          super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>(_iC)));
      }
      if (!neighborhood.getCellsInboundFrom(remoteC).empty()) {
        _recvTasks.emplace_back(std::make_unique<RecvTask>(
          _mpiCommunicator, tagCoordinator.get(remoteC, loadBalancer.glob(_iC)),
          loadBalancer.rank(remoteC),
          neighborhood.getFieldsCommonWith(remoteC),
          neighborhood.getCellsInboundFrom(remoteC),
          neighborhood.getMasksInboundFrom(remoteC),
          super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>(_iC)));
      }
    }
  });

#else // not using PARALLEL_MODE_MPI
  neighborhood.forNeighbors([&](int localC) {
    if (!neighborhood.getCellsInboundFrom(localC).empty()) {
      switch (loadBalancer.platform(loadBalancer.loc(localC))) {
        case Platform::GPU_SYCL:
          _copyTasks.emplace_back(new HomogeneousCopyTask(
            neighborhood.getFieldsCommonWith(localC),
            neighborhood.getCellsInboundFrom(localC),   super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>(_iC),
            neighborhood.getCellsRequestedFrom(localC), super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>(loadBalancer.loc(localC)),
            neighborhood.getMasksInboundFrom(localC)));
          break;
#ifdef PLATFORM_CPU_SIMD
        case Platform::CPU_SIMD:
          _copyTasks.emplace_back(new HeterogeneousCopyTask<T,DESCRIPTOR,Platform::CPU_SIMD,Platform::GPU_SYCL>(
            neighborhood.getFieldsCommonWith(localC),
            neighborhood.getCellsInboundFrom(localC),   super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>(_iC),
            neighborhood.getCellsRequestedFrom(localC), super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SIMD>>(loadBalancer.loc(localC)),
            neighborhood.getMasksInboundFrom(localC)));
          break;
#endif
#ifdef PLATFORM_CPU_SISD
        case Platform::CPU_SISD:
          _copyTasks.emplace_back(new HeterogeneousCopyTask<T,DESCRIPTOR,Platform::CPU_SISD,Platform::GPU_SYCL>(
            neighborhood.getFieldsCommonWith(localC),
            neighborhood.getCellsInboundFrom(localC),   super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>(_iC),
            neighborhood.getCellsRequestedFrom(localC), super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>>(loadBalancer.loc(localC)),
            neighborhood.getMasksInboundFrom(localC)));
          break;
#endif
        default:
          throw std::runtime_error("Invalid remote PLATFORM");
      }
    }
  });
#endif
}

template <typename T, typename DESCRIPTOR>
ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>::~ConcreteBlockCommunicator()
{ }

#ifdef PARALLEL_MODE_MPI

template <typename T, typename DESCRIPTOR>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>::receive()
{
  OLB_PROFILE_T("GPU_SYCL:Comm:receive", DESCRIPTOR);
  for (auto& task : _sendTasks) {
    task->prepare();
  }
  for (auto& task : _recvTasks) {
    task->receive();
  }
}

template <typename T, typename DESCRIPTOR>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>::send()
{
  OLB_PROFILE_T("GPU_SYCL:Comm:send", DESCRIPTOR);
  std::vector<::sycl::event> gatherEvents;
  gatherEvents.reserve(_sendTasks.size());

  for (auto& task : _sendTasks) {
    gatherEvents.push_back(task->submitGather());
  }

  ::sycl::event::wait(gatherEvents);

  for (auto& task : _copyTasks) {
    task->copy();
  }

  bool pending = true;
  while (pending) {
    pending = false;
    for (auto& task : _sendTasks) {
      if (!task->send()) {
        pending = true;
      }
    }
    if (pending) {
      std::this_thread::yield();
    }
  }
}

template <typename T, typename DESCRIPTOR>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>::unpack()
{
  OLB_PROFILE_T("GPU_SYCL:Comm:unpack", DESCRIPTOR);
  bool pending = true;
  while (pending) {
    pending = false;
    for (auto& task : _recvTasks) {
      if (task->poll()) {
        pending = true;
      }
    }
    if (pending) {
      std::this_thread::yield();
    }
  }
}

template <typename T, typename DESCRIPTOR>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>::wait()
{
  OLB_PROFILE_T("GPU_SYCL:Comm:wait", DESCRIPTOR);
  for (auto& task : _copyTasks) {
    task->wait();
  }
  for (auto& task : _recvTasks) {
    task->wait();
  }
  for (auto& task : _sendTasks) {
    task->wait();
  }
}

#else // not using PARALLEL_MODE_MPI

template <typename T, typename DESCRIPTOR>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>>::copy()
{
  for (auto& task : _copyTasks) {
    task->copy();
  }
  // Copies run asynchronously on the shared queue which no consumer waits
  // for implicitly. Join them before any subsequent kernel dispatch reads
  // the updated overlap data (mirrors the MPI path's RecvTask::wait loop).
  for (auto& task : _copyTasks) {
    task->wait();
  }
}

#endif

}

#endif
