/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022, 2025 Adrian Kummerlaender, Leonardo Dorneles
 *
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

#ifndef GPU_HIP_COLUMN_HH
#define GPU_HIP_COLUMN_HH

#include "column.h"
#include "device.h"

#include <thrust/gather.h>
#include <thrust/scatter.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/transform_iterator.h>

#include <cstring>




namespace olb {

struct WrapIndex {
  std::size_t shift;
  std::size_t count;

  __host__ __device__
  std::size_t operator()(const CellID& idx) const {
    std::size_t i = static_cast<std::size_t>(idx) + shift;
    i -= (i >= count) * count;
    return i;
  }
};


namespace gpu {

namespace hip {


template<typename T>
struct Column<T>::Data {
  thrust::host_vector<T>   host;
  thrust::device_vector<T> device;

  Data(std::size_t count):
    host(count),
    device(count)
  { }

  void resize(std::size_t newCount)
  {
    host.resize(newCount);
    device.resize(newCount);
  }

};

template<typename T>
Column<T>::Column(std::size_t count):
  _count(count),
  _data(new Data(count))
{ }

template<typename T>
Column<T>::Column(Column<T>&& rhs):
  _count(rhs._count),
  _data(rhs._data.release())
{ }

template<typename T>
Column<T>::Column(const Column<T>& rhs):
  _count(rhs._count),
  _data(new Data(*rhs._data))
{ }

template<typename T>
Column<T>::~Column()
{ }

template<typename T>
void Column<T>::clear()
{
  _data->device.clear();
  _data->host.clear();
  _count = _data->host.size();
}

template<typename T>
void Column<T>::resize(std::size_t newCount)
{
  _data->resize(newCount);
  _count = newCount;
}

template<typename T>
void Column<T>::push_back(T value)
{
  _data->host.push_back(value);
  _count = _data->host.size();
}

template<typename T>
void Column<T>::deduplicate()
{
  std::sort(_data->host.begin(), _data->host.end());
  _data->host.erase(std::unique(_data->host.begin(), _data->host.end()), _data->host.end());
}

template<typename T>
const T& Column<T>::operator[](std::size_t i) const
{
  return _data->host[i];
}

template<typename T>
T& Column<T>::operator[](std::size_t i)
{
  return _data->host[i];
}

template<typename T>
std::size_t Column<T>::size() const
{
  return _count;
}

template<typename T>
const T* Column<T>::data() const
{
  return _data->host.data();
}

template<typename T>
T* Column<T>::data()
{
  return _data->host.data();
}

template<typename T>
const T* Column<T>::deviceData() const
{
  return _data->device.data().get();
}

template<typename T>
T* Column<T>::deviceData()
{
  return _data->device.data().get();
}

template<typename T>
void Column<T>::setProcessingContext(ProcessingContext context)
{
  if (_count != _data->device.size()) {
    _data->device.resize(_count);
  }
  switch (context) {
  case ProcessingContext::Evaluation:
    device::copyToHost(_data->device.data().get(), _data->host.data(), size()*sizeof(T));
    return;
  case ProcessingContext::Simulation:
    device::copyToDevice(_data->host.data(), _data->device.data().get(), size()*sizeof(T));
    return;
  }
}

template<typename T>
void Column<T>::setProcessingContext(ProcessingContext context, device::Stream& stream)
{
  if (_count != _data->device.size()) {
    _data->device.resize(_count);
  }
  switch (context) {
  case ProcessingContext::Evaluation:
    device::asyncCopyToHost(stream, _data->device.data().get(), _data->host.data(), size()*sizeof(T));
    return;
  case ProcessingContext::Simulation:
    device::asyncCopyToDevice(stream, _data->host.data(), _data->device.data().get(), size()*sizeof(T));
    return;
  }
}

template<typename T>
std::size_t Column<T>::getNblock() const
{
  return 2;
}

template<typename T>
std::size_t Column<T>::getSerializableSize() const
{
  return _count * sizeof(T) + sizeof(std::size_t);
}

template<typename T>
bool* Column<T>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _count);
  if (loadingMode && iBlock == 1) {
    resize(_count);
  }
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, *data(), _count);

  return dataPtr;
}


// Branching ring: a single plain device allocation, two-pointer wrap on access.
template<typename T>
struct CyclicColumnData<T,gpu::CyclicStrategy::Branching> {
  std::unique_ptr<T[]> host;
  const std::size_t size;
  hipDeviceptr_t ptr;

  CyclicColumnData(std::size_t count):
    host(new T[count] { }),
    size(count * sizeof(T)),
    ptr(0)  // hipDeviceptr_t is an integer on NVIDIA, so 0 not nullptr
  {
    HIP_CHECK(hipMalloc(reinterpret_cast<void**>(&ptr), size));
    // hipDeviceptr_t is void* on AMD but an integer (CUdeviceptr) on the NVIDIA
    // platform, so cast explicitly for the void* hipMemset parameter.
    HIP_CHECK(hipMemset(reinterpret_cast<void*>(ptr), 0, size));
  }

  ~CyclicColumnData() {
    HIP_CHECK(hipFree(reinterpret_cast<void*>(ptr)));
  }

  static std::size_t alignedCount(std::size_t count) { return count; }
  T* base() const { return reinterpret_cast<T*>(ptr); }
};

// Virtual-memory ring: a double-mapped address range, branch-free on access.
template<typename T>
struct CyclicColumnData<T,gpu::CyclicStrategy::VirtualMemory> {
  std::unique_ptr<T[]> host;
  const std::size_t size;

  hipMemGenericAllocationHandle_t handle;
  hipMemAllocationProp prop{};
  hipMemAccessDesc access{};
  hipDeviceptr_t ptr;

  CyclicColumnData(std::size_t count):
    host(new T[count] { }),
    size(count * sizeof(T))
  {
    const int device = device::get();

    prop.type = hipMemAllocationTypePinned;
    prop.location.type = hipMemLocationTypeDevice;
    prop.location.id = device;
    HIP_CHECK(hipMemAddressReserve(&ptr, 2 * size, 0, 0, 0));

    HIP_CHECK(hipMemCreate(&handle, size, &prop, 0));
    HIP_CHECK(hipMemMap(ptr, size, 0, handle, 0));
    // hipDeviceptr_t is void* on AMD but an integer (CUdeviceptr) on NVIDIA; round-trip
    // through a byte pointer so the byte offset is portable across both platforms.
    HIP_CHECK(hipMemMap(reinterpret_cast<hipDeviceptr_t>(reinterpret_cast<std::uint8_t*>(ptr) + size),
                        size, 0, handle, 0));

    access.location.type = hipMemLocationTypeDevice;
    access.location.id = device;
    access.flags = hipMemAccessFlagsProtReadWrite;
    HIP_CHECK(hipMemSetAccess(ptr, 2 * size, &access, 1));
  }

  ~CyclicColumnData() {
    HIP_CHECK(hipMemUnmap(ptr, 2 * size));
    HIP_CHECK(hipMemRelease(handle));
    HIP_CHECK(hipMemAddressFree(ptr, 2 * size));
  }

  static std::size_t alignedCount(std::size_t count) { return device::getPageAlignedCount<T>(count); }
  T* base() const { return reinterpret_cast<T*>(ptr); }
};

template<typename T>
CyclicColumn<T>::CyclicColumn(std::size_t count):
  _count(Data::alignedCount(count)),
  _size(_count * sizeof(T)),
  _data(new Data(_count)),
  _shift(0)
{
  _hostData = _data->host.get();
  _deviceBase = _data->base();
  _devicePopulation = _deviceBase;
}

template<typename T>
CyclicColumn<T>::~CyclicColumn()
{ }

template<typename T>
const T& CyclicColumn<T>::operator[](std::size_t i) const
{
  return const_cast<CyclicColumn<T>*>(this)->getHostAccessor()[i];
}

template<typename T>
T& CyclicColumn<T>::operator[](std::size_t i)
{
  return getHostAccessor()[i];
}

template<typename T>
typename CyclicColumn<T>::Accessor CyclicColumn<T>::getHostAccessor()
{
  if constexpr (cyclicStrategy == gpu::CyclicStrategy::Branching) {
    return Accessor{
      .start0 = _data->host.get() + _shift,
      .start1 = _data->host.get() - (_count - _shift),
      .remainder = static_cast<std::ptrdiff_t>(_count - _shift - 1)
    };
  } else {
    return Accessor{ .ptr = _data->host.get() };
  }
}

template<typename T>
typename CyclicColumn<T>::Accessor CyclicColumn<T>::getDeviceAccessor()
{
  if constexpr (cyclicStrategy == gpu::CyclicStrategy::Branching) {
    return Accessor{
      .start0 = deviceData() + _shift,
      .start1 = deviceData() - (_count - _shift),
      .remainder = static_cast<std::ptrdiff_t>(_count - _shift - 1)
    };
  } else {
    return Accessor{ .ptr = _devicePopulation };
  }
}

template<typename T>
void CyclicColumn<T>::setProcessingContext(ProcessingContext context)
{
  switch (context) {
  case ProcessingContext::Evaluation:
    device::copyToHost(deviceData(), _data->host.get(), _size);
    return;
  case ProcessingContext::Simulation:
    device::copyToDevice(_data->host.get(), deviceData(), _size);
    return;
  }
}

template<typename T>
std::size_t CyclicColumn<T>::getNblock() const
{
  return 3;
}

template<typename T>
std::size_t CyclicColumn<T>::getSerializableSize() const
{
  return _count * sizeof(T) + sizeof(std::ptrdiff_t) + sizeof(std::size_t);
}

template<typename T>
bool* CyclicColumn<T>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _shift);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _count);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, *_data->host.get(), _count);

  return dataPtr;
}


}

}


template <typename T>
std::size_t ConcreteCommunicatable<gpu::hip::Column<T>>::serialize(ConstSpan<CellID> indices,
                                                                    std::uint8_t* buffer) const
{
  thrust::gather(thrust::device,
                 thrust::device_pointer_cast(indices.begin()),
                 thrust::device_pointer_cast(indices.end()),
                 thrust::device_pointer_cast(_column.deviceData()),
                 thrust::device_pointer_cast(reinterpret_cast<T*>(buffer)));
  return indices.size() * sizeof(T);
}

template <typename T>
std::size_t ConcreteCommunicatable<gpu::hip::Column<T>>::deserialize(ConstSpan<CellID> indices,
                                                                      const std::uint8_t* buffer)
{
  thrust::scatter(thrust::device,
                  thrust::device_pointer_cast(reinterpret_cast<const T*>(buffer)),
                  thrust::device_pointer_cast(reinterpret_cast<const T*>(buffer) + indices.size()),
                  thrust::device_pointer_cast(indices.begin()),
                  thrust::device_pointer_cast(_column.deviceData()));
  return indices.size() * sizeof(T);
}

template <typename T>
std::size_t ConcreteCommunicatable<gpu::hip::CyclicColumn<T>>::serialize(
  ConstSpan<CellID> indices,
  std::uint8_t* buffer) const
{
  thrust::gather(thrust::device,
                 thrust::device_pointer_cast(indices.begin()),
                 thrust::device_pointer_cast(indices.end()),
                 thrust::device_pointer_cast(_column.deviceData()),
                 thrust::device_pointer_cast(reinterpret_cast<T*>(buffer)));
  return indices.size() * sizeof(T);
}

template <typename T>
std::size_t ConcreteCommunicatable<gpu::hip::CyclicColumn<T>>::deserialize(
  ConstSpan<CellID> indices,
  const std::uint8_t* buffer)
{
  thrust::scatter(thrust::device,
                  thrust::device_pointer_cast(reinterpret_cast<const T*>(buffer)),
                  thrust::device_pointer_cast(reinterpret_cast<const T*>(buffer) + indices.size()),
                  thrust::device_pointer_cast(indices.begin()),
                  thrust::device_pointer_cast(_column.deviceData()));
  return indices.size() * sizeof(T);
}


}

#include "sparseColumn.hh"

#endif
