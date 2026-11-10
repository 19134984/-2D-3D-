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

#ifndef GPU_CUDA_REGISTRY_HH
#define GPU_CUDA_REGISTRY_HH

#include "registry.h"

#include <functional>

namespace olb {

namespace gpu {

namespace cuda {

/// Mapping of TYPE in CONTEXT to runtime-fixed index
/**
 * Used for dynamic access to field arrays
 **/
template <typename CONTEXT, typename TYPE>
__constant__ std::size_t field_type_index;

/// Host-side version of gpu::cuda::AnyDeviceFieldArrayD
struct FieldArrayPointer {
  gpu::cuda::device::unique_ptr<void*> data;
  const unsigned column_count;
  const unsigned element_size;
  const unsigned tile;
  const unsigned accessor_stride;
};

}

}

template <typename T, typename DESCRIPTOR>
struct FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::Data {
  Data():
    indexOnHost(2, nullptr),
    dynamicsTableHost(1, nullptr),
    dynamicsTableModified(true)
  { }

  std::map<std::type_index,gpu::cuda::FieldArrayPointer> fieldArrayPointers;
  /// Fields whose rotated device pointers must be refreshed each step
  /// (populated by the communicator for any-fields exchanges)
  std::set<std::type_index> perStepRefresh;
  /// Typed refreshers rebuilding each sparse field's stored device accessor after a sync
  std::vector<std::function<void()>> sparseRefreshers;

  std::vector<void**>                   indexOnHost;
  gpu::cuda::device::unique_ptr<void**> indexOnDevice;

  /// Distinct dynamics type -> stable non-zero dispatch-table index
  std::map<std::type_index,std::uint16_t> dynamicsIndex;
  /// Host mirror of the per-block dynamics dispatch table (index 0 == nullptr)
  std::vector<void*>                      dynamicsTableHost;
  gpu::cuda::device::unique_ptr<void*>    dynamicsTableDevice;
  bool                                    dynamicsTableModified;

};

template <typename T, typename DESCRIPTOR>
FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::FieldTypeRegistry():
  _data(new Data()),
  _modified{true} { }

template <typename T, typename DESCRIPTOR>
FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::~FieldTypeRegistry()
{ }

template <typename T, typename DESCRIPTOR>
template <typename FIELD_TYPE>
void FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::track(
  AnyFieldTypeD<T,DESCRIPTOR,Platform::GPU_CUDA>* fieldType)
{
  _index.template set<FIELD_TYPE>(fieldType);

  // Copy host-side field index to device-side field index
  std::size_t index = _index.template index<FIELD_TYPE>();
  cudaMemcpyToSymbol(gpu::cuda::field_type_index<FieldTypeRegistry,FIELD_TYPE>, &index, sizeof(std::size_t));
  gpu::cuda::device::check();
  if (index >= _data->indexOnHost.size()) {
    _data->indexOnHost.resize(2*index);
  }

  // Update device-side pointers to Array<FIELD> field types
  using ConcreteFieldType = typename FIELD_TYPE::template type<T,DESCRIPTOR,Platform::GPU_CUDA>;
  if constexpr (std::is_base_of_v<ColumnVectorBase, ConcreteFieldType>) {
    using field_t = typename ConcreteFieldType::field_t;
    constexpr unsigned dim = DESCRIPTOR::template size<field_t>();
    auto& fieldArray = *fieldType->template as<FIELD_TYPE>();

    if constexpr (detail::isSparseColumnBacked<T,DESCRIPTOR,Platform::GPU_CUDA,field_t>()) {
      // A sparse field stores its tile-directory/payload accessor per dimension. The device pointers
      // change on each sync, so it is refreshed afterwards (setProcessingContext -> refreshSparseAccessors)
      using Accessor = typename ConcreteFieldType::column_type::Accessor;
      std::array<Accessor,dim> accessors;
      for (unsigned iD=0; iD < dim; ++iD) {
        accessors[iD] = fieldArray[iD].getDeviceAccessor();
      }
      auto accessorsOnDevice = gpu::cuda::device::malloc<Accessor>(dim);
      gpu::cuda::device::copyToDevice(accessors.data(), accessorsOnDevice, dim*sizeof(Accessor));
      _data->indexOnHost[index] = reinterpret_cast<void**>(accessorsOnDevice);
      _data->fieldArrayPointers.emplace(typeid(field_t), gpu::cuda::FieldArrayPointer{
        gpu::cuda::device::unique_ptr<void*>(reinterpret_cast<void**>(accessorsOnDevice)),
        dim,
        sizeof(precision_scalar_t<typename field_t::template value_type<T>, precisionOf<DESCRIPTOR,field_t>()>),
        Accessor::tile_size,
        sizeof(Accessor)
      });
      _data->perStepRefresh.insert(typeid(field_t));
      _data->sparseRefreshers.push_back([this, &fieldArray]() {
        this->template refreshDeviceFieldArray<field_t>(fieldArray);
      });
    }
    else {
      // A dense field stores raw device pointers, layout-identical to its DeviceColumnAccessor
      std::array<void*,dim> componentPointers;
      for (unsigned iD=0; iD < dim; ++iD) {
        componentPointers[iD] = fieldArray[iD].deviceData();
      }
      auto componentPointersOnDevice = gpu::cuda::device::malloc<void*>(dim);
      gpu::cuda::device::copyToDevice(componentPointers.data(),
                                      componentPointersOnDevice,
                                      dim*sizeof(void*));
      _data->indexOnHost[index] = componentPointersOnDevice;
      _data->fieldArrayPointers.emplace(typeid(field_t), gpu::cuda::FieldArrayPointer{
        gpu::cuda::device::unique_ptr<void*>(componentPointersOnDevice),
        dim,
        sizeof(precision_scalar_t<typename field_t::template value_type<T>, precisionOf<DESCRIPTOR,field_t>()>),
        0,
        0
      });
    }
  }

  _modified = true;
}

template<typename T, typename DESCRIPTOR>
void*** FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::deviceData()
{
  if (_modified) {
    _data->indexOnDevice = gpu::cuda::device::malloc<void**>(_data->indexOnHost.size());
    gpu::cuda::device::copyToDevice(_data->indexOnHost.data(),
                                    _data->indexOnDevice.get(),
                                    _data->indexOnHost.size()*sizeof(void**));
    _modified = false;
  }
  return _data->indexOnDevice.get();
}

template<typename T, typename DESCRIPTOR>
std::uint16_t FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::registerDynamics(
  std::type_index type, void* deviceDynamics)
{
  auto it = _data->dynamicsIndex.find(type);
  if (it != _data->dynamicsIndex.end()) {
    return it->second;
  }
  std::uint16_t index = _data->dynamicsTableHost.size();
  _data->dynamicsTableHost.push_back(deviceDynamics);
  _data->dynamicsIndex.emplace(type, index);
  _data->dynamicsTableModified = true;
  return index;
}

template<typename T, typename DESCRIPTOR>
gpu::cuda::Dynamics<T,DESCRIPTOR>**
FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::dynamicsTable()
{
  if (_data->dynamicsTableModified) {
    _data->dynamicsTableDevice = gpu::cuda::device::malloc<void*>(_data->dynamicsTableHost.size());
    gpu::cuda::device::copyToDevice(_data->dynamicsTableHost.data(),
                                    _data->dynamicsTableDevice.get(),
                                    _data->dynamicsTableHost.size()*sizeof(void*));
    _data->dynamicsTableModified = false;
  }
  return reinterpret_cast<gpu::cuda::Dynamics<T,DESCRIPTOR>**>(_data->dynamicsTableDevice.get());
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
void FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::refreshDeviceFieldArray(
  FieldArrayD<T,DESCRIPTOR,Platform::GPU_CUDA,FIELD>& fieldArray)
{
  // Skip the (synchronous) pointer upload unless a consumer reads the registry
  // pointers per step. The typed collision / specialized single-field halo use
  // fresh accessors and never read them, so the common propagation path never
  // pays this refresh.
  if (_data->perStepRefresh.find(typeid(FIELD)) == _data->perStepRefresh.end()) {
    return;
  }
  constexpr unsigned dim = DESCRIPTOR::template size<FIELD>();
  auto& ptr = _data->fieldArrayPointers.at(typeid(FIELD));
  if constexpr (detail::isSparseColumnBacked<T,DESCRIPTOR,Platform::GPU_CUDA,FIELD>()) {
    // Re-capture the sparse tile directory/payload pointers into the stored accessor
    using Accessor = typename FieldArrayD<T,DESCRIPTOR,Platform::GPU_CUDA,FIELD>::column_type::Accessor;
    std::array<Accessor,dim> accessors;
    for (unsigned iD=0; iD < dim; ++iD) {
      accessors[iD] = fieldArray[iD].getDeviceAccessor();
    }
    gpu::cuda::device::copyToDevice(accessors.data(),
                                    reinterpret_cast<Accessor*>(ptr.data.get()),
                                    dim*sizeof(Accessor));
  }
  else {
    std::array<void*,dim> componentPointers;
    for (unsigned iD=0; iD < dim; ++iD) {
      componentPointers[iD] = fieldArray[iD].deviceData();
    }
    gpu::cuda::device::copyToDevice(componentPointers.data(),
                                    ptr.data.get(),
                                    dim*sizeof(void*));
  }
}

template <typename T, typename DESCRIPTOR>
void FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::refreshSparseAccessors()
{
  for (auto& refresh : _data->sparseRefreshers) {
    refresh();
  }
}

template <typename T, typename DESCRIPTOR>
void FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::requestPerStepRefresh(std::type_index field)
{
  _data->perStepRefresh.insert(field);
}

template <typename T, typename DESCRIPTOR>
bool FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::needsPerStepRefresh(std::type_index field) const
{
  return _data->perStepRefresh.find(field) != _data->perStepRefresh.end();
}

template <typename T, typename DESCRIPTOR>
gpu::cuda::AnyDeviceFieldArrayD
FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::deviceFieldArray(
  std::type_index field)
{
  auto& ptr = _data->fieldArrayPointers.at(field);
  return gpu::cuda::AnyDeviceFieldArrayD{
    ptr.data.get(),
    ptr.column_count,
    ptr.element_size,
    ptr.tile,
    ptr.accessor_stride
  };
}

template <typename T, typename DESCRIPTOR>
std::vector<gpu::cuda::AnyDeviceFieldArrayD>
FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>::deviceFieldArrays(
  const std::vector<std::type_index>& fields)
{
  std::vector<gpu::cuda::AnyDeviceFieldArrayD> deviceFields;
  deviceFields.reserve(fields.size());
  for (std::type_index field : fields) {
    deviceFields.emplace_back(deviceFieldArray(field));
  }
  return deviceFields;
}

}

#endif
