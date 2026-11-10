/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Adrian Kummerlaender
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

#ifndef PLATFORM_GPU_CUDA_COMMUNICATOR_HH
#define PLATFORM_GPU_CUDA_COMMUNICATOR_HH

#include "communication/mpiRequest.h"
#include "communication/communicatable.h"
#include "communication/overlapBufferLayout.h"
#include "communication/directionalOverlap.h"
#include "communication/superCommunicationTagCoordinator.h"

#include "registry.h"
#include "context.hh"

#include "utilities/performanceProfiler.h"

#include <thrust/device_vector.h>

#ifdef PARALLEL_MODE_MPI
#include "mpi.h"
#if defined(OPEN_MPI) && OPEN_MPI
#include <mpi-ext.h>
#endif
#endif

namespace olb {

template <>
void checkPlatform<Platform::GPU_CUDA>()
{
  OstreamManager clout(std::cout, "GPU_CUDA");

  int nDevices{};
  cudaGetDeviceCount(&nDevices);

#ifdef OLB_DEBUG
  clout.setMultiOutput(true);
#endif
  if (nDevices < 1) {
    clout << "No CUDA device found" << std::endl;
  } else if (nDevices > 1) {
    clout << "Found " << nDevices << " CUDA devices but only one can be used per MPI process." << std::endl;
  } else {
    clout << "Found a single CUDA device." << std::endl;
  }
  for (int deviceI=0; deviceI < nDevices; ++deviceI) {
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, deviceI);
    clout << "Detected " << deviceProp.name << "." << std::endl;
  }
#ifdef OLB_DEBUG
  clout.setMultiOutput(false);
#endif

#ifdef PARALLEL_MODE_MPI
#if defined(MPIX_CUDA_AWARE_SUPPORT) && MPIX_CUDA_AWARE_SUPPORT
  if (!MPIX_Query_cuda_support()) {
    clout << "The used MPI Library is not CUDA-aware. Multi-GPU execution will fail." << std::endl;
  }
#endif
#if defined(MPIX_CUDA_AWARE_SUPPORT) && !MPIX_CUDA_AWARE_SUPPORT
  clout << "The used MPI Library is not CUDA-aware. Multi-GPU execution will fail." << std::endl;
#endif
#if !defined(MPIX_CUDA_AWARE_SUPPORT)
  clout << "Unable to check for CUDA-aware MPI support. Multi-GPU execution may fail." << std::endl;
#endif
#endif // PARALLEL_MODE_MPI
}

namespace gpu {

namespace cuda {

#ifdef FEATURE_SPARSE_SCRATCH_CANARY
/// Counts sparse halo scatters that fell into the scratch tile (unallocated inbound cell)
__device__ unsigned long long sparse_scratch_writes = 0;

/// Read and reset the scratch-write counter and fail loudly if any inbound halo cell was
/// scattered into the discard scratch tile (ensureAllocated missed a tile). Called from the
/// nvcc-compiled scatter wrappers so no host-callable device symbol crosses the codegen.
inline void checkSparseScratchCanary() {
  cudaDeviceSynchronize();
  unsigned long long n = 0;
  cudaMemcpyFromSymbol(&n, sparse_scratch_writes, sizeof(n));
  if (n > 0) {
    const unsigned long long zero = 0;
    cudaMemcpyToSymbol(sparse_scratch_writes, &zero, sizeof(zero));
    throw std::runtime_error("sparse halo scatter hit the scratch tile "
                             + std::to_string(n) + " time(s): an inbound halo tile was not pre-allocated");
  }
}
#endif

namespace kernel {

/// CUDA kernel for gathering FIELD data of lattice at indices into buffer
template <typename CONTEXT, typename FIELD>
void gather_field(CONTEXT lattice,
                  const CellID* indices, std::size_t nIndices,
                  ::olb::detail::field_scalar_t<typename CONTEXT::value_t, typename CONTEXT::descriptor_t, FIELD>* buffer) __global__ {
  const CellID iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nIndices)) {
    return;
  }
  if constexpr (descriptors::is_propagatable_field<FIELD>::value
                && cyclicStrategyFor(Platform::GPU_CUDA) != CyclicStrategy::VirtualMemory) {
    for (unsigned iD=0; iD < CONTEXT::descriptor_t::template size<FIELD>(); ++iD) {
      buffer[iD*nIndices+iIndex] = lattice.getPopulation(iD)[indices[iIndex]];
    }
  } else {
    auto* field = lattice.template getField<FIELD>();
    for (unsigned iD=0; iD < CONTEXT::descriptor_t::template size<FIELD>(); ++iD) {
      buffer[iD*nIndices+iIndex] = field[iD].get(indices[iIndex]);
    }
  }
}

/// CUDA kernel for copying FIELD data of sourceLattice at sourceIndices into targetLattice at targetIndices
template <typename SOURCE, typename TARGET, typename FIELD>
void copy_field(SOURCE sourceLattice, TARGET targetLattice,
                const CellID* sourceIndices,
                const CellID* targetIndices,
                std::size_t nIndices) __global__ {
  const CellID iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nIndices)) {
    return;
  }
  if constexpr (descriptors::is_propagatable_field<FIELD>::value
                && cyclicStrategyFor(Platform::GPU_CUDA) != CyclicStrategy::VirtualMemory) {
    for (unsigned iD=0; iD < SOURCE::descriptor_t::template size<FIELD>(); ++iD) {
      targetLattice.getPopulation(iD)[targetIndices[iIndex]]
        = sourceLattice.getPopulation(iD)[sourceIndices[iIndex]];
    }
  } else {
    auto* source = sourceLattice.template getField<FIELD>();
    auto* target = targetLattice.template getField<FIELD>();
    for (unsigned iD=0; iD < SOURCE::descriptor_t::template size<FIELD>(); ++iD) {
      target[iD][targetIndices[iIndex]] = source[iD].get(sourceIndices[iIndex]);
    }
  }
}

/// CUDA kernel for gathering fields at indices into buffer
__global__ void gather_any_fields(AnyDeviceFieldArrayD* fields, std::size_t nFields,
                                  const CellID* indices, std::size_t nIndices,
                                  std::uint8_t* buffer) {
  const CellID iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nIndices)) {
    return;
  }
  for (unsigned iField=0; iField < nFields; ++iField) {
    auto& field = fields[iField];
    for (unsigned iD=0; iD < field.column_count; ++iD) {
      field.gather(iD, indices[iIndex],
                   buffer + (iD*nIndices + iIndex)*field.element_size);
    }
    buffer += nIndices*field.column_count*field.element_size;
  }
}

/// CUDA kernel for copying sourceFields at sourceIndices to targetFields at targetIndices
/**
 * source and target fields may be of different block lattices but must represent the
 * same field types in the same sequence
 **/
__global__ void copy_any_fields(AnyDeviceFieldArrayD* sourceFields,
                                AnyDeviceFieldArrayD* targetFields,
                                std::size_t nFields,
                                const CellID* sourceIndices,
                                const CellID* targetIndices,
                                std::size_t nIndices) {
  const CellID iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nIndices)) {
    return;
  }
  for (unsigned iField=0; iField < nFields; ++iField) {
    auto& sourceField = sourceFields[iField];
    auto& targetField = targetFields[iField];
    for (unsigned iD=0; iD < sourceField.column_count; ++iD) {
      std::uint8_t scratch[sizeof(double)*8];
      sourceField.gather(iD, sourceIndices[iIndex], scratch);
      targetField.scatter(iD, targetIndices[iIndex], scratch);
    }
  }
}

/// CUDA kernel for scattering FIELD data in buffer to indices in lattice
template <typename CONTEXT, typename FIELD>
void scatter_field(CONTEXT lattice,
                   const CellID* indices, std::size_t nIndices,
                   ::olb::detail::field_scalar_t<typename CONTEXT::value_t, typename CONTEXT::descriptor_t, FIELD>* buffer) __global__ {
  const CellID iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nIndices)) {
    return;
  }
  if constexpr (descriptors::is_propagatable_field<FIELD>::value
                && cyclicStrategyFor(Platform::GPU_CUDA) != CyclicStrategy::VirtualMemory) {
    for (unsigned iD=0; iD < CONTEXT::descriptor_t::template size<FIELD>(); ++iD) {
      lattice.getPopulation(iD)[indices[iIndex]] = buffer[iD*nIndices+iIndex];
    }
  } else {
    auto* field = lattice.template getField<FIELD>();
    for (unsigned iD=0; iD < CONTEXT::descriptor_t::template size<FIELD>(); ++iD) {
      field[iD][indices[iIndex]] = buffer[iD*nIndices+iIndex];
    }
  }
}

/// CUDA kernel for scattering fields in buffer to indices in lattice
__global__ void scatter_any_fields(AnyDeviceFieldArrayD* fields, std::size_t nFields,
                                   const CellID* indices, std::size_t nIndices,
                                   std::uint8_t* buffer) {
  const CellID iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nIndices)) {
    return;
  }
  for (unsigned iField=0; iField < nFields; ++iField) {
    auto& field = fields[iField];
    for (unsigned iD=0; iD < field.column_count; ++iD) {
      field.scatter(iD, indices[iIndex],
                    buffer + (iD*nIndices + iIndex)*field.element_size);
    }
    buffer += nIndices*field.column_count*field.element_size;
  }
}

}

/// Blocking gather of FIELD at given indices into buffer
template <typename FIELD, typename CONTEXT>
void gather_field(CONTEXT& lattice, const thrust::device_vector<CellID>& indices, std::uint8_t* buffer) {
  const auto block_size = 32;
  const auto block_count = (indices.size() + block_size - 1) / block_size;
  kernel::gather_field<CONTEXT,FIELD><<<block_count,block_size>>>(
    lattice,
    indices.data().get(), indices.size(),
    reinterpret_cast<::olb::detail::field_scalar_t<typename CONTEXT::value_t, typename CONTEXT::descriptor_t, FIELD>*>(buffer));
  device::check();
}

/// Non-blocking gather of FIELD at given indices into buffer
template <typename FIELD, typename CONTEXT>
void async_gather_field(cudaStream_t stream,
                        CONTEXT& lattice,
                        const thrust::device_vector<CellID>& indices,
                        std::uint8_t* buffer) {
  const auto block_size = 32;
  const auto block_count = (indices.size() + block_size - 1) / block_size;
  kernel::gather_field<CONTEXT,FIELD><<<block_count,block_size,0,stream>>>(
    lattice,
    indices.data().get(), indices.size(),
    reinterpret_cast<::olb::detail::field_scalar_t<typename CONTEXT::value_t, typename CONTEXT::descriptor_t, FIELD>*>(buffer));
  device::check();
}

/// Non-blocking copy of FIELD at given indices from sourceLattice to targetLattice
template <typename FIELD, typename SOURCE, typename TARGET>
void async_copy_field(cudaStream_t stream,
                      SOURCE& sourceLattice,
                      TARGET& targetLattice,
                      const thrust::device_vector<CellID>& sourceIndices,
                      const thrust::device_vector<CellID>& targetIndices) {
  const auto block_size = 32;
  const auto block_count = (sourceIndices.size() + block_size - 1) / block_size;
  kernel::copy_field<SOURCE,TARGET,FIELD><<<block_count,block_size,0,stream>>>(
    sourceLattice,
    targetLattice,
    sourceIndices.data().get(),
    targetIndices.data().get(),
    sourceIndices.size());
  device::check();
}


/// Blocking gather of fields at given indices into buffer
void gather_any_fields(thrust::device_vector<AnyDeviceFieldArrayD>& fields,
                       const thrust::device_vector<CellID>& indices,
                       std::uint8_t* buffer) {
  const auto block_size = 32;
  const auto block_count = (indices.size() + block_size - 1) / block_size;
  kernel::gather_any_fields<<<block_count,block_size>>>(
    fields.data().get(), fields.size(),
    indices.data().get(), indices.size(),
    buffer);
  device::check();
}

/// Non-blocking gather of fields at given indices into buffer
void async_gather_any_fields(cudaStream_t stream,
                             thrust::device_vector<AnyDeviceFieldArrayD>& fields,
                             const thrust::device_vector<CellID>& indices,
                             std::uint8_t* buffer) {
  const auto block_size = 32;
  const auto block_count = (indices.size() + block_size - 1) / block_size;
  kernel::gather_any_fields<<<block_count,block_size,0,stream>>>(
    fields.data().get(), fields.size(),
    indices.data().get(), indices.size(),
    buffer);
  device::check();
}

/// Non-blocking copy of fields at given indices from sourceIndices to targetIndices
void async_copy_any_fields(cudaStream_t stream,
                           thrust::device_vector<AnyDeviceFieldArrayD>& sourceFields,
                           thrust::device_vector<AnyDeviceFieldArrayD>& targetFields,
                           const thrust::device_vector<CellID>& sourceIndices,
                           const thrust::device_vector<CellID>& targetIndices) {
  const auto block_size = 32;
  const auto block_count = (sourceIndices.size() + block_size - 1) / block_size;
  kernel::copy_any_fields<<<block_count,block_size,0,stream>>>(
    sourceFields.data().get(),  targetFields.data().get(),  sourceFields.size(),
    sourceIndices.data().get(), targetIndices.data().get(), sourceIndices.size());
  device::check();
}

/// Blocking scatter of FIELD data in buffer to given indices
template <typename FIELD, typename CONTEXT>
void scatter_field(CONTEXT& lattice, const thrust::device_vector<CellID>& indices, std::uint8_t* buffer) {
  const auto block_size = 32;
  const auto block_count = (indices.size() + block_size - 1) / block_size;
  kernel::scatter_field<CONTEXT,FIELD><<<block_count,block_size>>>(
    lattice,
    indices.data().get(), indices.size(),
    reinterpret_cast<::olb::detail::field_scalar_t<typename CONTEXT::value_t, typename CONTEXT::descriptor_t, FIELD>*>(buffer));
  device::check();
}

/// Non-blocking scatter of FIELD data in buffer to given indices
template <typename FIELD, typename CONTEXT>
void async_scatter_field(cudaStream_t stream,
                         CONTEXT& lattice,
                         const thrust::device_vector<CellID>& indices,
                         std::uint8_t* buffer) {
  const auto block_size = 32;
  const auto block_count = (indices.size() + block_size - 1) / block_size;
  kernel::scatter_field<CONTEXT,FIELD><<<block_count,block_size,0,stream>>>(
    lattice,
    indices.data().get(), indices.size(),
    reinterpret_cast<::olb::detail::field_scalar_t<typename CONTEXT::value_t, typename CONTEXT::descriptor_t, FIELD>*>(buffer));
  device::check();
}

/// Blocking scatter of fields data in buffer to given indices
void scatter_any_fields(thrust::device_vector<AnyDeviceFieldArrayD>& fields,
                        const thrust::device_vector<CellID>& indices,
                        std::uint8_t* buffer) {
  const auto block_size = 32;
  const auto block_count = (indices.size() + block_size - 1) / block_size;
  kernel::scatter_any_fields<<<block_count,block_size>>>(
    fields.data().get(), fields.size(),
    indices.data().get(), indices.size(),
    buffer);
  device::check();
#ifdef FEATURE_SPARSE_SCRATCH_CANARY
  checkSparseScratchCanary();
#endif
}

/// Non-blocking scatter of fields data in buffer to given indices
void async_scatter_any_fields(cudaStream_t stream,
                              thrust::device_vector<AnyDeviceFieldArrayD>& fields,
                              const thrust::device_vector<CellID>& indices,
                              std::uint8_t* buffer) {
  const auto block_size = 32;
  const auto block_count = (indices.size() + block_size - 1) / block_size;
  kernel::scatter_any_fields<<<block_count,block_size,0,stream>>>(
    fields.data().get(), fields.size(),
    indices.data().get(), indices.size(),
    buffer);
  device::check();
#ifdef FEATURE_SPARSE_SCRATCH_CANARY
  checkSparseScratchCanary();
#endif
}

namespace kernel {

/// CUDA kernel gathering one mask group of a masked propagatable segment
/**
 * Implements one group of overlap::MaskGroupLayout: components of the
 * group's mask, component-major over the group's cells. The field's
 * components are read through the shift-aware cyclic accessor bundle
 * (built launch-time from getField<FIELD>()[iD].getDeviceAccessor()), so
 * POPULATION and any further propagatable field share one code path.
 **/
template <typename ACCESSORS, typename T>
void gather_propagatable_masked(ACCESSORS field,
                              const CellID* cells,
                              const std::uint32_t* cellIndices, std::size_t nCells,
                              const std::uint8_t* components, unsigned nComponents,
                              T* buffer) __global__ {
  const std::size_t iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nCells)) {
    return;
  }
  const CellID cell = cells[cellIndices[iIndex]];
  for (unsigned iC=0; iC < nComponents; ++iC) {
    buffer[iC*nCells + iIndex] = field[components[iC]][cell];
  }
}

/// CUDA kernel scattering one mask group of a masked propagatable segment
template <typename ACCESSORS, typename T>
void scatter_propagatable_masked(ACCESSORS field,
                               const CellID* cells,
                               const std::uint32_t* cellIndices, std::size_t nCells,
                               const std::uint8_t* components, unsigned nComponents,
                               const T* buffer) __global__ {
  const std::size_t iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nCells)) {
    return;
  }
  const CellID cell = cells[cellIndices[iIndex]];
  for (unsigned iC=0; iC < nComponents; ++iC) {
    field[components[iC]][cell] = buffer[iC*nCells + iIndex];
  }
}

}

/// Build the launch-time device accessor bundle of FIELD's columns
template <typename FIELD, typename BLOCK>
auto propagatableAccessors(BLOCK& block) {
  constexpr unsigned Q = BLOCK::descriptor_t::template size<FIELD>();
  // accessor element follows the field's stored scalar (T, or reduced precision) so
  // a store_precision field packs at its storage precision, halving halo traffic too
  using acc_t = std::decay_t<decltype(block.template getField<FIELD>()[0].getDeviceAccessor())>;
  return Vector<acc_t, Q>([&block](unsigned iD) {
    return block.template getField<FIELD>()[iD].getDeviceAccessor();
  });
}


namespace kernel {

/// CUDA kernel copying one mask group of a masked POPULATION exchange
/**
 * Direct lattice-to-lattice variant: cellIndices index both the source
 * and target cell arrays (the negotiated lists are index-aligned).
 **/
template <typename SOURCE_ACC, typename TARGET_ACC>
void copy_propagatable_masked(SOURCE_ACC sourceField, TARGET_ACC targetField,
                            const CellID* sourceCells, const CellID* targetCells,
                            const std::uint32_t* cellIndices, std::size_t nCells,
                            const std::uint8_t* components, unsigned nComponents) __global__ {
  const std::size_t iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nCells)) {
    return;
  }
  const std::uint32_t i = cellIndices[iIndex];
  for (unsigned iC=0; iC < nComponents; ++iC) {
    targetField[components[iC]][targetCells[i]]
      = sourceField[components[iC]][sourceCells[i]];
  }
}

}

/// Device-resident group tables of an overlap::MaskGroupLayout
/**
 * Flattened cell-index and component lists with per-group host-side
 * offsets driving one kernel launch per group.
 **/
struct MaskedSegmentTables {
  thrust::device_vector<std::uint32_t> cellIndices;
  thrust::device_vector<std::uint8_t>  components;
  std::vector<std::size_t> cellBegin;
  std::vector<std::size_t> componentBegin;

  MaskedSegmentTables() = default;

  explicit MaskedSegmentTables(const overlap::MaskGroupLayout& layout) {
    std::vector<std::uint32_t> flatCells;
    std::vector<std::uint8_t>  flatComponents;
    for (const auto& group : layout.groups) {
      cellBegin.emplace_back(flatCells.size());
      componentBegin.emplace_back(flatComponents.size());
      flatCells.insert(flatCells.end(), group.cellIndices.begin(), group.cellIndices.end());
      flatComponents.insert(flatComponents.end(), group.components.begin(), group.components.end());
    }
    cellIndices = flatCells;
    components  = flatComponents;
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
    tables.emplace(field, MaskedSegmentTables(layout));
  }
  return tables;
}


/// Launch the per-group masked gather kernels for one FIELD
template <typename FIELD, typename BLOCK>
void async_gather_propagatable_masked(cudaStream_t stream,
                                    BLOCK& block,
                                    const thrust::device_vector<CellID>& cells,
                                    const overlap::MaskGroupLayout& layout,
                                    const MaskedSegmentTables& tables,
                                    std::uint8_t* buffer) {
  using T = typename BLOCK::value_t;
  auto field = propagatableAccessors<FIELD>(block);
  const auto block_size = 32;
  for (unsigned iG=0; iG < layout.groups.size(); ++iG) {
    const auto& group = layout.groups[iG];
    const auto block_count = (group.cellIndices.size() + block_size - 1) / block_size;
    kernel::gather_propagatable_masked<<<block_count,block_size,0,stream>>>(
      field,
      cells.data().get(),
      tables.cellIndices.data().get() + tables.cellBegin[iG], group.cellIndices.size(),
      tables.components.data().get() + tables.componentBegin[iG], group.components.size(),
      // Buffer packs at the field's stored scalar (half under store_precision), not sizeof(T)
      reinterpret_cast<::olb::detail::field_scalar_t<T,typename BLOCK::descriptor_t,FIELD>*>(buffer + group.begin));
  }
  device::check();
}

/// Launch the per-group masked scatter kernels for one FIELD
template <typename FIELD, typename BLOCK>
void async_scatter_propagatable_masked(cudaStream_t stream,
                                     BLOCK& block,
                                     const thrust::device_vector<CellID>& cells,
                                     const overlap::MaskGroupLayout& layout,
                                     const MaskedSegmentTables& tables,
                                     const std::uint8_t* buffer) {
  using T = typename BLOCK::value_t;
  auto field = propagatableAccessors<FIELD>(block);
  const auto block_size = 32;
  for (unsigned iG=0; iG < layout.groups.size(); ++iG) {
    const auto& group = layout.groups[iG];
    const auto block_count = (group.cellIndices.size() + block_size - 1) / block_size;
    kernel::scatter_propagatable_masked<<<block_count,block_size,0,stream>>>(
      field,
      cells.data().get(),
      tables.cellIndices.data().get() + tables.cellBegin[iG], group.cellIndices.size(),
      tables.components.data().get() + tables.componentBegin[iG], group.components.size(),
      reinterpret_cast<const ::olb::detail::field_scalar_t<T,typename BLOCK::descriptor_t,FIELD>*>(buffer + group.begin));
  }
  device::check();
}

/// Launch the per-group masked copy kernels for one propagatable FIELD
template <typename FIELD, typename BLOCK>
void async_copy_propagatable_masked(cudaStream_t stream,
                                  BLOCK& source, BLOCK& target,
                                  const thrust::device_vector<CellID>& sourceCells,
                                  const thrust::device_vector<CellID>& targetCells,
                                  const overlap::MaskGroupLayout& layout,
                                  const MaskedSegmentTables& tables) {
  auto sourceField = propagatableAccessors<FIELD>(source);
  auto targetField = propagatableAccessors<FIELD>(target);
  const auto block_size = 32;
  for (unsigned iG=0; iG < layout.groups.size(); ++iG) {
    const auto& group = layout.groups[iG];
    const auto block_count = (group.cellIndices.size() + block_size - 1) / block_size;
    kernel::copy_propagatable_masked<<<block_count,block_size,0,stream>>>(
      sourceField, targetField,
      sourceCells.data().get(), targetCells.data().get(),
      tables.cellIndices.data().get() + tables.cellBegin[iG], group.cellIndices.size(),
      tables.components.data().get() + tables.componentBegin[iG], group.components.size());
  }
  device::check();
}

/// Device-side field-array lists per AnyField segment, empty otherwise
using SegmentDenseArrays = std::vector<thrust::device_vector<AnyDeviceFieldArrayD>>;

/// Device field arrays for `block`'s AnyField segments. Requests their refresh
/**
 * Only the type-erased (AnyField) segments read the registry's device
 * pointers per step and thus need the per-step refresh. The typed masked
 * and dense segments build fresh launch-time accessors.
 **/
template <typename BLOCK>
SegmentDenseArrays segmentDenseArrays(BLOCK& block,
                                      const overlap::MultiFieldPartition& partition)
{
  SegmentDenseArrays denseArrays(partition.segments.size());
  for (std::size_t i=0; i < partition.segments.size(); ++i) {
    if (partition.segments[i].kind == overlap::FieldSegment::Kind::AnyField) {
      denseArrays[i] = block.getDataRegistry().deviceFieldArrays({partition.segments[i].field});
      block.getDataRegistry().requestPerStepRefresh(partition.segments[i].field);
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
void async_gather_segments(cudaStream_t stream, BLOCK& block,
                           const overlap::MultiFieldPartition& partition,
                           SegmentDenseArrays& denseArrays,
                           const thrust::device_vector<CellID>& cells,
                           const FieldMaskLayouts& layouts,
                           const FieldMaskTables& tables,
                           std::uint8_t* buffer) {
  using T = typename BLOCK::value_t;
  using DESCRIPTOR = typename BLOCK::descriptor_t;
  for (std::size_t i=0; i < partition.segments.size(); ++i) {
    const auto& seg = partition.segments[i];
    std::uint8_t* at = buffer + seg.begin;
    switch (seg.kind) {
    case overlap::FieldSegment::Kind::MaskedPropagatable:
      dispatchMaskedField<BLOCK>(seg.field, [&](auto id) {
        async_gather_propagatable_masked<typename decltype(id)::type>(
          stream, block, cells, layouts.at(seg.field), tables.at(seg.field), at);
      });
      break;
    case overlap::FieldSegment::Kind::TypedDense: {
      DeviceContext<T,DESCRIPTOR> lattice(block);
      async_gather_field<descriptors::POPULATION>(stream, lattice, cells, at);
      break;
    }
    case overlap::FieldSegment::Kind::AnyField:
      async_gather_any_fields(stream, denseArrays[i], cells, at);
      break;
    }
  }
}

/// Scatter every segment of a packed multi-field buffer (see makeSegmentPartition)
template <typename BLOCK>
void async_scatter_segments(cudaStream_t stream, BLOCK& block,
                            const overlap::MultiFieldPartition& partition,
                            SegmentDenseArrays& denseArrays,
                            const thrust::device_vector<CellID>& cells,
                            const FieldMaskLayouts& layouts,
                            const FieldMaskTables& tables,
                            std::uint8_t* buffer) {
  using T = typename BLOCK::value_t;
  using DESCRIPTOR = typename BLOCK::descriptor_t;
  for (std::size_t i=0; i < partition.segments.size(); ++i) {
    const auto& seg = partition.segments[i];
    std::uint8_t* at = buffer + seg.begin;
    switch (seg.kind) {
    case overlap::FieldSegment::Kind::MaskedPropagatable:
      dispatchMaskedField<BLOCK>(seg.field, [&](auto id) {
        async_scatter_propagatable_masked<typename decltype(id)::type>(
          stream, block, cells, layouts.at(seg.field), tables.at(seg.field), at);
      });
      break;
    case overlap::FieldSegment::Kind::TypedDense: {
      DeviceContext<T,DESCRIPTOR> lattice(block);
      async_scatter_field<descriptors::POPULATION>(stream, lattice, cells, at);
      break;
    }
    case overlap::FieldSegment::Kind::AnyField:
      async_scatter_any_fields(stream, denseArrays[i], cells, at);
      break;
    }
  }
}

/// Copy every segment of a multi-field exchange source->target (no packed buffer)
template <typename BLOCK>
void async_copy_segments(cudaStream_t stream, BLOCK& source, BLOCK& target,
                         const overlap::MultiFieldPartition& partition,
                         SegmentDenseArrays& sourceDense,
                         SegmentDenseArrays& targetDense,
                         const thrust::device_vector<CellID>& sourceCells,
                         const thrust::device_vector<CellID>& targetCells,
                         const FieldMaskLayouts& layouts,
                         const FieldMaskTables& tables) {
  using T = typename BLOCK::value_t;
  using DESCRIPTOR = typename BLOCK::descriptor_t;
  for (std::size_t i=0; i < partition.segments.size(); ++i) {
    const auto& seg = partition.segments[i];
    switch (seg.kind) {
    case overlap::FieldSegment::Kind::MaskedPropagatable:
      dispatchMaskedField<BLOCK>(seg.field, [&](auto id) {
        async_copy_propagatable_masked<typename decltype(id)::type>(
          stream, source, target, sourceCells, targetCells, layouts.at(seg.field), tables.at(seg.field));
      });
      break;
    case overlap::FieldSegment::Kind::TypedDense: {
      DeviceContext<T,DESCRIPTOR> sourceLattice(source);
      DeviceContext<T,DESCRIPTOR> targetLattice(target);
      async_copy_field<descriptors::POPULATION>(stream, sourceLattice, targetLattice,
        sourceCells, targetCells);
      break;
    }
    case overlap::FieldSegment::Kind::AnyField:
      async_copy_any_fields(stream, sourceDense[i], targetDense[i], sourceCells, targetCells);
      break;
    }
  }
}

}

}

}

#ifdef FEATURE_DIRECTIONAL_CANARY
#include "canary.h"
#endif

namespace olb {

#ifdef PARALLEL_MODE_MPI

/// Wrapper for a non-blocking block propagation send request
template <typename T, typename DESCRIPTOR>
class ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::SendTask {
private:
  const bool _onlyPopulationField;

  /// Negotiated per-field, per-cell component masks restrict each masked
  /// field's segment to its overlap::MaskGroupLayout (e.g. directional
  /// streaming). Empty masks keep the dense layout.
  const gpu::cuda::FieldMaskLayouts _maskLayouts;
  const gpu::cuda::FieldMaskTables _maskTables;
  const bool _masked;

  /// Every field becomes one packed segment in request order (see
  /// makeSegmentPartition): masked propagatable fields take the shared
  /// MaskGroupLayout, POPULATION in a dense exchange the typed kernel, the
  /// rest the type-erased path (which alone needs the per-step refresh).
  /// The byte layout matches the CPU-side (Masked)MultiConcreteCommunicatable.
  gpu::cuda::SegmentDenseArrays _denseArrays;
  const overlap::MultiFieldPartition _segments;

  const thrust::device_vector<CellID> _cells;

  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& _source;

  gpu::cuda::device::unique_ptr<std::uint8_t> _buffer;
  std::unique_ptr<MpiSendRequest> _request;

  std::unique_ptr<gpu::cuda::device::Stream> _stream;
  Event<Platform::GPU_CUDA> _computeDone;

public:
  SendTask(MPI_Comm comm, int tag, int rank,
           const std::vector<std::type_index>& fields,
           const std::vector<CellID>& cells,
           const overlap::FieldMasks& masks,
           ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block):
    _onlyPopulationField(fields.size() == 1 && fields[0] == typeid(descriptors::POPULATION)),
    _maskLayouts(gpu::cuda::maskLayoutsFor<T,DESCRIPTOR>(fields, masks)),
    _maskTables(gpu::cuda::maskTablesFor(_maskLayouts)),
    _masked(!_maskLayouts.empty()),
    _segments(gpu::cuda::makeSegmentPartition(block, fields, cells, _maskLayouts, _denseArrays)),
    _cells(cells),
    _source(block),
    _stream(std::make_unique<gpu::cuda::device::Stream>(cudaStreamNonBlocking, []{ int l=0,g=0; cudaDeviceGetStreamPriorityRange(&l,&g); return g; }()))
  {
    _buffer = gpu::cuda::device::malloc<std::uint8_t>(_segments.size);
    _request = std::make_unique<MpiSendRequest>(
      _buffer.get(), _segments.size, rank, tag, comm);
  }

  ~SendTask()
  {
    _stream->synchronize();
    wait();
  }

  void prepare()
  {
    OLB_PROFILE_T("GPU_CUDA:Comm:SendTask:prepare", DESCRIPTOR);
    // Order the gather after this block's outer collision and post
    // processors (named non-blocking streams) as well as all legacy
    // default-stream work (couplings, PerBlock operators). The
    // non-blocking transfer stream waits for none of them implicitly.
    _source.outerCollisionDone().barrierOn(*_stream);
    _source.postProcessorsDone().barrierOn(*_stream);
    _computeDone.record(cudaStreamLegacy);
    _computeDone.barrierOn(*_stream);
    if (_onlyPopulationField && !_masked) {
      gpu::cuda::DeviceContext<T,DESCRIPTOR> lattice(_source);
      gpu::cuda::async_gather_field<descriptors::POPULATION>(_stream->get(), lattice, _cells, _buffer.get());
    } else {
      gpu::cuda::async_gather_segments(_stream->get(), _source, _segments, _denseArrays,
                                       _cells, _maskLayouts, _maskTables, _buffer.get());
    }
  }

  void send()
  {
    OLB_PROFILE_T("GPU_CUDA:Comm:SendTask:send", DESCRIPTOR);
    _stream->synchronize();
    _request->start();
  }

  void wait()
  {
    OLB_PROFILE_T("GPU_CUDA:Comm:SendTask:wait", DESCRIPTOR);
    _request->wait();
  }
};

/// Wrapper for a non-blocking block propagation receive request
template <typename T, typename DESCRIPTOR>
class ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::RecvTask {
private:
  const int _tag;
  const int _rank;

  const bool _onlyPopulationField;

  /// Same per-field masked handling as SendTask (see there)
  const gpu::cuda::FieldMaskLayouts _maskLayouts;
  const gpu::cuda::FieldMaskTables _maskTables;
  const bool _masked;
#ifdef FEATURE_DIRECTIONAL_CANARY
  /// Per-field, per-cell masks for the directional canary (populated only when enabled)
  const gpu::cuda::FieldCanaryMasks _canaryMasks;
#endif

  /// Same packed multi-segment routing as SendTask (see there)
  gpu::cuda::SegmentDenseArrays _denseArrays;
  const overlap::MultiFieldPartition _segments;

  const thrust::device_vector<CellID> _cells;

  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& _target;

  gpu::cuda::device::unique_ptr<std::uint8_t> _buffer;
  std::unique_ptr<MpiRecvRequest> _request;

  std::unique_ptr<gpu::cuda::device::Stream> _stream;

public:
  /// Manual replacement for std::reference_wrapper<RecvTask>
  /**
   * Used to track pending receive requests in std::set.
   *
   * This is a workaround for problematic external definition of
   * dependently-typed comparision operators for nested classes.
   * Reconsider as soon as depending on C++17 is allowed.
   **/
  class ref {
  private:
    RecvTask& _task;
  public:
    ref(std::unique_ptr<RecvTask>& task): _task(*task) { };

    RecvTask* operator->() const
    {
      return &_task;
    }

    bool operator <(const ref& rhs) const
    {
      return _task < rhs._task;
    }
  };

  RecvTask(MPI_Comm comm, int tag, int rank,
           const std::vector<std::type_index>& fields,
           const std::vector<CellID>& cells,
           const overlap::FieldMasks& masks,
           ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block):
    _tag(tag),
    _rank(rank),
    _onlyPopulationField(fields.size() == 1 && fields[0] == typeid(descriptors::POPULATION)),
    _maskLayouts(gpu::cuda::maskLayoutsFor<T,DESCRIPTOR>(fields, masks)),
    _maskTables(gpu::cuda::maskTablesFor(_maskLayouts)),
    _masked(!_maskLayouts.empty()),
#ifdef FEATURE_DIRECTIONAL_CANARY
    _canaryMasks(gpu::cuda::canaryMasksFor(fields, masks)),
#endif
    _segments(gpu::cuda::makeSegmentPartition(block, fields, cells, _maskLayouts, _denseArrays)),
    _cells(cells),
    _target(block),
    _stream(std::make_unique<gpu::cuda::device::Stream>(cudaStreamNonBlocking, []{ int l=0,g=0; cudaDeviceGetStreamPriorityRange(&l,&g); return g; }()))
  {
    _buffer = gpu::cuda::device::malloc<std::uint8_t>(_segments.size);
    _request = std::make_unique<MpiRecvRequest>(
      _buffer.get(), _segments.size, _rank, _tag, comm);
  }

  ~RecvTask()
  {
    wait();
  }

  bool operator<(const RecvTask& rhs) const
  {
    return  _rank  < rhs._rank
        || (_rank == rhs._rank && _tag < rhs._tag);
  }

  void receive()
  {
    OLB_PROFILE_T("GPU_CUDA:Comm:RecvTask:receive", DESCRIPTOR);
    _request->start();
  };

  bool isDone()
  {
    return _request->isDone();
  }

  void unpack()
  {
    if (_onlyPopulationField && !_masked) {
      gpu::cuda::DeviceContext<T,DESCRIPTOR> lattice(_target);
      gpu::cuda::async_scatter_field<descriptors::POPULATION>(_stream->get(), lattice, _cells, _buffer.get());
    } else {
      gpu::cuda::async_scatter_segments(_stream->get(), _target, _segments, _denseArrays,
                                        _cells, _maskLayouts, _maskTables, _buffer.get());
#ifdef FEATURE_DIRECTIONAL_CANARY
      if (!_canaryMasks.empty()) {
        gpu::cuda::async_poison_segments(_stream->get(), _target, _segments, _cells, _canaryMasks, overlap::componentsPerDirection<DESCRIPTOR>());
      }
#endif
    }
  }

  void wait()
  {
    OLB_PROFILE_T("GPU_CUDA:Comm:RecvTask:wait", DESCRIPTOR);
    _stream->synchronize();
  }

};

#endif // PARALLEL_MODE_MPI

/// Wrapper for a local plain-copy block communication request
template <typename T, typename DESCRIPTOR>
struct ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::CopyTask {
  virtual ~CopyTask() { }

  virtual void copy() = 0;
  virtual void wait() = 0;
};

/// Wrapper for a local plain-copy block communication request
template <typename T, typename DESCRIPTOR>
class ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::HomogeneousCopyTask
  : public ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::CopyTask {
private:
  const bool _onlyPopulationField;

  /// Same per-field masked handling as SendTask (see there)
  const gpu::cuda::FieldMaskLayouts _maskLayouts;
  const gpu::cuda::FieldMaskTables _maskTables;
  const bool _masked;
#ifdef FEATURE_DIRECTIONAL_CANARY
  /// Per-field, per-cell masks for the directional canary (populated only when enabled)
  const gpu::cuda::FieldCanaryMasks _canaryMasks;
#endif

  /// Same packed multi-segment routing as SendTask. Direct copies have no
  /// packed buffer, so the segment byte offsets are unused
  gpu::cuda::SegmentDenseArrays _sourceDense;
  const overlap::MultiFieldPartition _segments;
  gpu::cuda::SegmentDenseArrays _targetDense;

  const thrust::device_vector<CellID> _targetCells;
  const thrust::device_vector<CellID> _sourceCells;

  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& _target;
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& _source;

  std::unique_ptr<gpu::cuda::device::Stream> _stream;
  Event<Platform::GPU_CUDA> _computeDone;

public:
  HomogeneousCopyTask(
    const std::vector<std::type_index>& fields,
    const std::vector<CellID>& targetCells, ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& target,
    const std::vector<CellID>& sourceCells, ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& source,
    const overlap::FieldMasks& masks = {}):
    _onlyPopulationField(fields.size() == 1 && fields[0] == typeid(descriptors::POPULATION)),
    _maskLayouts(gpu::cuda::maskLayoutsFor<T,DESCRIPTOR>(fields, masks)),
    _maskTables(gpu::cuda::maskTablesFor(_maskLayouts)),
    _masked(!_maskLayouts.empty()),
#ifdef FEATURE_DIRECTIONAL_CANARY
    _canaryMasks(gpu::cuda::canaryMasksFor(fields, masks)),
#endif
    _segments(gpu::cuda::makeSegmentPartition(source, fields, sourceCells, _maskLayouts, _sourceDense)),
    _targetDense(gpu::cuda::segmentDenseArrays(target, _segments)),
    _targetCells(targetCells),
    _sourceCells(sourceCells),
    _target(target),
    _source(source),
    _stream(std::make_unique<gpu::cuda::device::Stream>(cudaStreamNonBlocking, []{ int l=0,g=0; cudaDeviceGetStreamPriorityRange(&l,&g); return g; }()))
  {
    OLB_ASSERT(_sourceCells.size() == _targetCells.size(),
               "Source cell count must match target cell count");
  }

  ~HomogeneousCopyTask()
  {
    wait();
  }

  void copy() override
  {
    OLB_PROFILE_T("GPU_CUDA:Comm:CopyTask:copy", DESCRIPTOR);
    // Order the copy after all compute dispatched so far (see SendTask::prepare)
    _source.outerCollisionDone().barrierOn(*_stream);
    _source.postProcessorsDone().barrierOn(*_stream);
    _computeDone.record(cudaStreamLegacy);
    _computeDone.barrierOn(*_stream);
    if (_onlyPopulationField && !_masked) {
      gpu::cuda::DeviceContext<T,DESCRIPTOR> sourceLattice(_source);
      gpu::cuda::DeviceContext<T,DESCRIPTOR> targetLattice(_target);
      gpu::cuda::async_copy_field<descriptors::POPULATION>(_stream->get(), sourceLattice, targetLattice, _sourceCells, _targetCells);
    } else {
      gpu::cuda::async_copy_segments(_stream->get(), _source, _target, _segments,
                                     _sourceDense, _targetDense,
                                     _sourceCells, _targetCells, _maskLayouts, _maskTables);
#ifdef FEATURE_DIRECTIONAL_CANARY
      if (!_canaryMasks.empty()) {
        gpu::cuda::async_poison_segments(_stream->get(), _target, _segments, _targetCells, _canaryMasks, overlap::componentsPerDirection<DESCRIPTOR>());
      }
#endif
    }
  }

  void wait() override
  {
    OLB_PROFILE_T("GPU_CUDA:Comm:CopyTask:wait", DESCRIPTOR);
    _stream->synchronize();
  }

};


/// Private implementation of heterogeneous copy task between CPU_* source and GPU_CUDA target
template <typename T, typename DESCRIPTOR, Platform SOURCE>
class HeterogeneousCopyTaskDataForGpuTarget : public ConcreteHeterogeneousCopyTask {
private:
  const bool _onlyPopulationField;

  /// Same per-field masked handling as SendTask (see there). The CPU side
  /// serializes the identical mask-group layouts via MaskedMultiCommunicatable
  const gpu::cuda::FieldMaskLayouts _maskLayouts;
  const gpu::cuda::FieldMaskTables _maskTables;
  const bool _masked;
#ifdef FEATURE_DIRECTIONAL_CANARY
  const gpu::cuda::FieldCanaryMasks _canaryMasks;
#endif

  /// Same packed multi-segment routing as SendTask. The CPU-side serialize
  /// packs fields sequentially in request order, matching the segment offsets
  gpu::cuda::SegmentDenseArrays _denseArrays;
  const overlap::MultiFieldPartition _segments;

  const thrust::device_vector<CellID> _targetCells;
  const std::vector<CellID>&          _sourceCells;

  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& _target;

  std::unique_ptr<Communicatable> _source;

  std::unique_ptr<gpu::cuda::device::Stream> _stream;

  gpu::cuda::Column<std::uint8_t> _buffer;

public:
  HeterogeneousCopyTaskDataForGpuTarget(
    const std::vector<std::type_index>& fields,
    const std::vector<CellID>& targetCells, ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& target,
    const std::vector<CellID>& sourceCells, ConcreteBlockLattice<T,DESCRIPTOR,SOURCE>&             source,
    const overlap::FieldMasks& masks = {}):
    _onlyPopulationField(fields.size() == 1 && fields[0] == typeid(descriptors::POPULATION)),
    _maskLayouts(gpu::cuda::maskLayoutsFor<T,DESCRIPTOR>(fields, masks)),
    _maskTables(gpu::cuda::maskTablesFor(_maskLayouts)),
    _masked(!_maskLayouts.empty()),
#ifdef FEATURE_DIRECTIONAL_CANARY
    _canaryMasks(gpu::cuda::canaryMasksFor(fields, masks)),
#endif
    _segments(gpu::cuda::makeSegmentPartition(target, fields, targetCells, _maskLayouts, _denseArrays)),
    _targetCells(targetCells),
    _sourceCells(sourceCells),
    _target(target),
    _source(overlap::makeMaskedMultiCommunicatable(source, fields, masks)),
    _stream(std::make_unique<gpu::cuda::device::Stream>(cudaStreamNonBlocking, []{ int l=0,g=0; cudaDeviceGetStreamPriorityRange(&l,&g); return g; }())),
    _buffer(_source->size(_sourceCells))
  { }

  void copy() override {
    _source->serialize(_sourceCells, _buffer.data());
    _buffer.setProcessingContext(ProcessingContext::Simulation, *_stream);

    if (_onlyPopulationField && !_masked) {
      gpu::cuda::DeviceContext<T,DESCRIPTOR> lattice(_target);
      gpu::cuda::async_scatter_field<descriptors::POPULATION>(_stream->get(), lattice, _targetCells, _buffer.deviceData());
    } else {
      gpu::cuda::async_scatter_segments(_stream->get(), _target, _segments, _denseArrays,
                                        _targetCells, _maskLayouts, _maskTables, _buffer.deviceData());
#ifdef FEATURE_DIRECTIONAL_CANARY
      if (!_canaryMasks.empty()) {
        gpu::cuda::async_poison_segments(_stream->get(), _target, _segments, _targetCells, _canaryMasks, overlap::componentsPerDirection<DESCRIPTOR>());
      }
#endif
    }
  }

  void wait() override {
    _stream->synchronize();
  }

};

/// Private implementation of heterogeneous copy task between GPU_CUDA source and CPU_* target
template <typename T, typename DESCRIPTOR, Platform TARGET>
class HeterogeneousCopyTaskDataForGpuSource : public ConcreteHeterogeneousCopyTask {
private:
  const bool _onlyPopulationField;

  /// Same per-field masked handling as SendTask (see there). The CPU side
  /// deserializes the identical mask-group layouts via MaskedMultiCommunicatable
  const gpu::cuda::FieldMaskLayouts _maskLayouts;
  const gpu::cuda::FieldMaskTables _maskTables;
  const bool _masked;

  /// Same packed multi-segment routing as SendTask (see there)
  gpu::cuda::SegmentDenseArrays _denseArrays;
  const overlap::MultiFieldPartition _segments;

  const std::vector<CellID>&          _targetCells;
  const thrust::device_vector<CellID> _sourceCells;

  std::unique_ptr<Communicatable> _target;
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& _source;

  std::unique_ptr<gpu::cuda::device::Stream> _stream;
  Event<Platform::GPU_CUDA> _computeDone;

  gpu::cuda::Column<std::uint8_t> _buffer;

public:
  HeterogeneousCopyTaskDataForGpuSource(
    const std::vector<std::type_index>& fields,
    const std::vector<CellID>& targetCells, ConcreteBlockLattice<T,DESCRIPTOR,TARGET>&             target,
    const std::vector<CellID>& sourceCells, ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& source,
    const overlap::FieldMasks& masks = {}):
    _onlyPopulationField(fields.size() == 1 && fields[0] == typeid(descriptors::POPULATION)),
    _maskLayouts(gpu::cuda::maskLayoutsFor<T,DESCRIPTOR>(fields, masks)),
    _maskTables(gpu::cuda::maskTablesFor(_maskLayouts)),
    _masked(!_maskLayouts.empty()),
    _segments(gpu::cuda::makeSegmentPartition(source, fields, sourceCells, _maskLayouts, _denseArrays)),
    _targetCells(targetCells),
    _sourceCells(sourceCells),
    _target(overlap::makeMaskedMultiCommunicatable(target, fields, masks)),
    _source(source),
    _stream(std::make_unique<gpu::cuda::device::Stream>(cudaStreamNonBlocking, []{ int l=0,g=0; cudaDeviceGetStreamPriorityRange(&l,&g); return g; }())),
    _buffer(_target->size(_targetCells))
  { }

  void copy() override {
    // Order the gather after all compute dispatched so far (see SendTask::prepare)
    _source.outerCollisionDone().barrierOn(*_stream);
    _source.postProcessorsDone().barrierOn(*_stream);
    _computeDone.record(cudaStreamLegacy);
    _computeDone.barrierOn(*_stream);
    if (_onlyPopulationField && !_masked) {
      gpu::cuda::DeviceContext<T,DESCRIPTOR> lattice(_source);
      gpu::cuda::async_gather_field<descriptors::POPULATION>(_stream->get(), lattice, _sourceCells, _buffer.deviceData());
    } else {
      gpu::cuda::async_gather_segments(_stream->get(), _source, _segments, _denseArrays,
                                       _sourceCells, _maskLayouts, _maskTables, _buffer.deviceData());
    }

    _buffer.setProcessingContext(ProcessingContext::Evaluation, *_stream);
  }

  void wait() override {
    _stream->synchronize();
    _target->deserialize(_targetCells, _buffer.data());
  }

};

template <typename T, typename DESCRIPTOR, Platform SOURCE>
HeterogeneousCopyTask<T,DESCRIPTOR,SOURCE,Platform::GPU_CUDA>::HeterogeneousCopyTask(
  const std::vector<std::type_index>& fields,
  const std::vector<CellID>& targetCells, ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& target,
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
void HeterogeneousCopyTask<T,DESCRIPTOR,SOURCE,Platform::GPU_CUDA>::copy()
{
  _impl->copy();
}

template <typename T, typename DESCRIPTOR, Platform SOURCE>
void HeterogeneousCopyTask<T,DESCRIPTOR,SOURCE,Platform::GPU_CUDA>::wait()
{
  _impl->wait();
}

template <typename T, typename DESCRIPTOR, Platform TARGET>
HeterogeneousCopyTask<T,DESCRIPTOR,Platform::GPU_CUDA,TARGET>::HeterogeneousCopyTask(
  const std::vector<std::type_index>& fields,
  const std::vector<CellID>& targetCells, ConcreteBlockLattice<T,DESCRIPTOR,TARGET>&             target,
  const std::vector<CellID>& sourceCells, ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& source,
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
void HeterogeneousCopyTask<T,DESCRIPTOR,Platform::GPU_CUDA,TARGET>::copy()
{
  _impl->copy();
}

template <typename T, typename DESCRIPTOR, Platform TARGET>
void HeterogeneousCopyTask<T,DESCRIPTOR,Platform::GPU_CUDA,TARGET>::wait()
{
  _impl->wait();
}


template <typename T, typename DESCRIPTOR>
ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::ConcreteBlockCommunicator(
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
  // Simulation sync. The device scatter never allocates (no-op for dense fields)
  {
    auto& block = super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>(_iC);
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
          case Platform::GPU_CUDA:
            // Use manual copy for local GPU-GPU communication due to better performance
            _copyTasks.emplace_back(new HomogeneousCopyTask(
              neighborhood.getFieldsCommonWith(remoteC),
              neighborhood.getCellsInboundFrom(remoteC),   super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>(_iC),
              neighborhood.getCellsRequestedFrom(remoteC), super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>(loadBalancer.loc(remoteC)),
              neighborhood.getMasksInboundFrom(remoteC)));
            break;
          case Platform::CPU_SIMD:
            // Use manual copy for local GPU-CPU communication due to better performance
            _copyTasks.emplace_back(new HeterogeneousCopyTask<T,DESCRIPTOR,Platform::CPU_SIMD,Platform::GPU_CUDA>(
              neighborhood.getFieldsCommonWith(remoteC),
              neighborhood.getCellsInboundFrom(remoteC),   super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>(_iC),
              neighborhood.getCellsRequestedFrom(remoteC), super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SIMD>>(loadBalancer.loc(remoteC)),
              neighborhood.getMasksInboundFrom(remoteC)));
            break;
          case Platform::CPU_SISD:
            // Use manual copy for local GPU-CPU communication due to better performance
            _copyTasks.emplace_back(new HeterogeneousCopyTask<T,DESCRIPTOR,Platform::CPU_SISD,Platform::GPU_CUDA>(
              neighborhood.getFieldsCommonWith(remoteC),
              neighborhood.getCellsInboundFrom(remoteC),   super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>(_iC),
              neighborhood.getCellsRequestedFrom(remoteC), super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>>(loadBalancer.loc(remoteC)),
              neighborhood.getMasksInboundFrom(remoteC)));
            break;
          default:
            throw std::runtime_error("Invalid remote PLATFORM");
        }
      }
    } else {
      // Handling of non-local GPU-GPU and GPU-^GPU communication in general
      if (!neighborhood.getCellsOutboundTo(remoteC).empty()) {
        _sendTasks.emplace_back(std::make_unique<SendTask>(
          _mpiCommunicator, tagCoordinator.get(loadBalancer.glob(_iC), remoteC),
          loadBalancer.rank(remoteC),
          neighborhood.getFieldsCommonWith(remoteC),
          neighborhood.getCellsOutboundTo(remoteC),
          neighborhood.getMasksOutboundTo(remoteC),
          super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>(_iC)));
      }
      if (!neighborhood.getCellsInboundFrom(remoteC).empty()) {
        _recvTasks.emplace_back(std::make_unique<RecvTask>(
          _mpiCommunicator, tagCoordinator.get(remoteC, loadBalancer.glob(_iC)),
          loadBalancer.rank(remoteC),
          neighborhood.getFieldsCommonWith(remoteC),
          neighborhood.getCellsInboundFrom(remoteC),
          neighborhood.getMasksInboundFrom(remoteC),
          super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>(_iC)));
      }
    }
  });

#else // not using PARALLEL_MODE_MPI
  neighborhood.forNeighbors([&](int localC) {
    if (!neighborhood.getCellsInboundFrom(localC).empty()) {
      switch (loadBalancer.platform(loadBalancer.loc(localC))) {
        case Platform::GPU_CUDA:
          _copyTasks.emplace_back(new HomogeneousCopyTask(
            neighborhood.getFieldsCommonWith(localC),
            neighborhood.getCellsInboundFrom(localC),   super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>(_iC),
            neighborhood.getCellsRequestedFrom(localC), super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>(loadBalancer.loc(localC)),
            neighborhood.getMasksInboundFrom(localC)));
          break;
        case Platform::CPU_SIMD:
          _copyTasks.emplace_back(new HeterogeneousCopyTask<T,DESCRIPTOR,Platform::CPU_SIMD,Platform::GPU_CUDA>(
            neighborhood.getFieldsCommonWith(localC),
            neighborhood.getCellsInboundFrom(localC),   super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>(_iC),
            neighborhood.getCellsRequestedFrom(localC), super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SIMD>>(loadBalancer.loc(localC)),
            neighborhood.getMasksInboundFrom(localC)));
          break;
        case Platform::CPU_SISD:
          _copyTasks.emplace_back(new HeterogeneousCopyTask<T,DESCRIPTOR,Platform::CPU_SISD,Platform::GPU_CUDA>(
            neighborhood.getFieldsCommonWith(localC),
            neighborhood.getCellsInboundFrom(localC),   super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>(_iC),
            neighborhood.getCellsRequestedFrom(localC), super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>>(loadBalancer.loc(localC)),
            neighborhood.getMasksInboundFrom(localC)));
          break;
        default:
          throw std::runtime_error("Invalid remote PLATFORM");
      }
    }
  });
#endif
}

template <typename T, typename DESCRIPTOR>
ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::~ConcreteBlockCommunicator()
{ }

#ifdef PARALLEL_MODE_MPI

template <typename T, typename DESCRIPTOR>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::receive()
{
  OLB_PROFILE_T("GPU_CUDA:Comm:receive", DESCRIPTOR);
  for (auto& task : _recvTasks) {
    task->receive();
  }
}

template <typename T, typename DESCRIPTOR>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::send()
{
  OLB_PROFILE_T("GPU_CUDA:Comm:send", DESCRIPTOR);
  for (auto& task : _sendTasks) {
    task->prepare();
  }
  for (auto& task : _sendTasks) {
    task->send();
  }
  for (auto& task : _copyTasks) {
    task->copy();
  }
}

template <typename T, typename DESCRIPTOR>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::unpack()
{
  OLB_PROFILE_T("GPU_CUDA:Comm:unpack", DESCRIPTOR);
  std::set<typename RecvTask::ref> pending(_recvTasks.begin(), _recvTasks.end());
  while (!pending.empty()) {
    auto task_iterator = pending.begin();
    while (task_iterator != pending.end()) {
      auto& task = *task_iterator;
      if (task->isDone()) {
        task->unpack();
        task_iterator = pending.erase(task_iterator);
      }
      else {
        ++task_iterator;
      }
    }
  }
}

template <typename T, typename DESCRIPTOR>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::wait()
{
  OLB_PROFILE_T("GPU_CUDA:Comm:wait", DESCRIPTOR);
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
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::copy()
{
  for (auto& task : _copyTasks) {
    task->copy();
  }
  // Scatters run on non-blocking streams which no consumer waits for
  // implicitly. Join them before any subsequent kernel dispatch reads
  // the updated overlap data (mirrors the MPI path's RecvTask::wait loop).
  for (auto& task : _copyTasks) {
    task->wait();
  }
}

#endif

}

#endif
