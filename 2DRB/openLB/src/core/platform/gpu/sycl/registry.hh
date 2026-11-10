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

#ifndef GPU_SYCL_REGISTRY_HH
#define GPU_SYCL_REGISTRY_HH

#include "registry.h"
#include "device.h"

#include <sycl/sycl.hpp>

#include <set>
#include <functional>

namespace olb {

namespace gpu {

namespace sycl {

template <typename CONTEXT, typename TYPE>
::sycl::ext::oneapi::experimental::device_global<std::size_t> field_type_index;

struct FieldArrayPointer {
  gpu::sycl::device::unique_ptr<void*> data;
  const unsigned column_count;
  const unsigned element_size;
  const unsigned tile;
  const unsigned accessor_stride;
};

}

}

template <typename T, typename DESCRIPTOR>
struct FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_SYCL>::Data {
  Data():
    indexOnHost(2, nullptr)
  { }

  std::map<std::type_index, gpu::sycl::FieldArrayPointer> fieldArrayPointers;
  /// Fields whose device pointers must be refreshed after a sync
  std::set<std::type_index> perStepRefresh;
  /// Typed refreshers rebuilding each sparse field's stored device accessor after a sync
  std::vector<std::function<void()>> sparseRefreshers;

  std::vector<void**>                     indexOnHost;
  gpu::sycl::device::unique_ptr<void**>   indexOnDevice;
};

template <typename T, typename DESCRIPTOR>
FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_SYCL>::FieldTypeRegistry():
  _data(new Data()),
  _modified{true}
{ }

template <typename T, typename DESCRIPTOR>
FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_SYCL>::~FieldTypeRegistry()
{ }

template <typename T, typename DESCRIPTOR>
template <typename FIELD_TYPE>
void FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_SYCL>::track(
  AnyFieldTypeD<T,DESCRIPTOR,Platform::GPU_SYCL>* fieldType)
{
  _index.template set<FIELD_TYPE>(fieldType);

  std::size_t index = _index.template index<FIELD_TYPE>();

  gpu::sycl::device::getQueue().copy(
    &index,
    gpu::sycl::field_type_index<FieldTypeRegistry,FIELD_TYPE>).wait();

  if (index >= _data->indexOnHost.size()) {
    _data->indexOnHost.resize(2*index);
  }

  using ConcreteFieldType = typename FIELD_TYPE::template type<T,DESCRIPTOR,Platform::GPU_SYCL>;

  if constexpr (std::is_base_of_v<ColumnVectorBase, ConcreteFieldType>) {
    using field_t = typename ConcreteFieldType::field_t;
    constexpr unsigned dim = DESCRIPTOR::template size<field_t>();
    auto& fieldArray = *fieldType->template as<FIELD_TYPE>();

    if constexpr (detail::isSparseColumnBacked<T,DESCRIPTOR,Platform::GPU_SYCL,field_t>()) {
      // A sparse field stores its tile-directory/payload accessor per dimension. The device pointers
      // change on each sync, so it is refreshed afterwards (setProcessingContext -> refreshSparseAccessors)
      using Accessor = typename ConcreteFieldType::column_type::Accessor;
      std::array<Accessor,dim> accessors;
      for (unsigned iD=0; iD < dim; ++iD) {
        accessors[iD] = fieldArray[iD].getDeviceAccessor();
      }
      auto accessorsOnDevice = gpu::sycl::device::malloc<Accessor>(dim);
      gpu::sycl::device::copyToDevice(accessors.data(), accessorsOnDevice, dim*sizeof(Accessor));
      _data->indexOnHost[index] = reinterpret_cast<void**>(accessorsOnDevice);
      _data->fieldArrayPointers.emplace(typeid(field_t), gpu::sycl::FieldArrayPointer{
        gpu::sycl::device::unique_ptr<void*>(reinterpret_cast<void**>(accessorsOnDevice)),
        dim,
        (unsigned)sizeof(precision_scalar_t<typename field_t::template value_type<T>, precisionOf<DESCRIPTOR,field_t>()>),
        Accessor::tile_size,
        (unsigned)sizeof(Accessor)
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
      auto componentPointersOnDevice = gpu::sycl::device::malloc<void*>(dim);
      gpu::sycl::device::copyToDevice(componentPointers.data(),
                                      componentPointersOnDevice,
                                      dim*sizeof(void*));
      _data->indexOnHost[index] = componentPointersOnDevice;
      _data->fieldArrayPointers.emplace(typeid(field_t), gpu::sycl::FieldArrayPointer{
        gpu::sycl::device::unique_ptr<void*>(componentPointersOnDevice),
        dim,
        (unsigned)sizeof(precision_scalar_t<typename field_t::template value_type<T>, precisionOf<DESCRIPTOR,field_t>()>),
        0,
        0
      });
    }
  }

  _modified = true;
}

template<typename T, typename DESCRIPTOR>
void*** FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_SYCL>::deviceData()
{
  if (_modified) {
    _data->indexOnDevice = gpu::sycl::device::malloc<void**>(_data->indexOnHost.size());
    gpu::sycl::device::copyToDevice(_data->indexOnHost.data(),
                                    _data->indexOnDevice.get(),
                                    _data->indexOnHost.size()*sizeof(void**));
    _modified = false;
  }
  return _data->indexOnDevice.get();
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
void FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_SYCL>::refreshDeviceFieldArray(
  FieldArrayD<T,DESCRIPTOR,Platform::GPU_SYCL,FIELD>& fieldArray)
{
  if (_data->perStepRefresh.find(typeid(FIELD)) == _data->perStepRefresh.end()) {
    return;
  }
  constexpr unsigned dim = DESCRIPTOR::template size<FIELD>();
  auto& ptr = _data->fieldArrayPointers.at(typeid(FIELD));
  if constexpr (detail::isSparseColumnBacked<T,DESCRIPTOR,Platform::GPU_SYCL,FIELD>()) {
    // Re-capture the sparse tile directory/payload pointers into the stored accessor
    using Accessor = typename FieldArrayD<T,DESCRIPTOR,Platform::GPU_SYCL,FIELD>::column_type::Accessor;
    std::array<Accessor,dim> accessors;
    for (unsigned iD=0; iD < dim; ++iD) {
      accessors[iD] = fieldArray[iD].getDeviceAccessor();
    }
    gpu::sycl::device::copyToDevice(accessors.data(),
                                    reinterpret_cast<Accessor*>(ptr.data.get()),
                                    dim*sizeof(Accessor));
  }
  else {
    std::array<void*,dim> componentPointers;
    for (unsigned iD=0; iD < dim; ++iD) {
      componentPointers[iD] = fieldArray[iD].deviceData();
    }
    gpu::sycl::device::copyToDevice(componentPointers.data(),
                                    ptr.data.get(),
                                    dim*sizeof(void*));
  }
}

template <typename T, typename DESCRIPTOR>
void FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_SYCL>::refreshSparseAccessors()
{
  for (auto& refresh : _data->sparseRefreshers) {
    refresh();
  }
}

template <typename T, typename DESCRIPTOR>
void FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_SYCL>::requestPerStepRefresh(std::type_index field)
{
  _data->perStepRefresh.insert(field);
}

template <typename T, typename DESCRIPTOR>
bool FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_SYCL>::needsPerStepRefresh(std::type_index field) const
{
  return _data->perStepRefresh.find(field) != _data->perStepRefresh.end();
}

template <typename T, typename DESCRIPTOR>
gpu::sycl::AnyDeviceFieldArrayD
FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_SYCL>::deviceFieldArray(
  std::type_index field)
{
  auto& ptr = _data->fieldArrayPointers.at(field);
  return gpu::sycl::AnyDeviceFieldArrayD{
    ptr.data.get(),
    ptr.column_count,
    ptr.element_size,
    ptr.tile,
    ptr.accessor_stride
  };
}

template <typename T, typename DESCRIPTOR>
std::vector<gpu::sycl::AnyDeviceFieldArrayD>
FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_SYCL>::deviceFieldArrays(
  const std::vector<std::type_index>& fields)
{
  std::vector<gpu::sycl::AnyDeviceFieldArrayD> deviceFields;
  deviceFields.reserve(fields.size());
  for (std::type_index field : fields) {
    deviceFields.emplace_back(deviceFieldArray(field));
  }
  return deviceFields;
}

}

#endif
