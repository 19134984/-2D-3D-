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
#ifndef GPU_SYCL_COLUMN_HH
#define GPU_SYCL_COLUMN_HH

#include "column.h"
#include "device.h"

#include <algorithm>
#include <memory>
#include <optional>
#include <string>
#include <stdexcept>

#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <asm/unistd_64.h>

#if __has_include(<sycl/ext/oneapi/virtual_mem.hpp>)
#include <sycl/ext/oneapi/virtual_mem.hpp>
#endif

#if defined(SYCL_EXT_ONEAPI_VIRTUAL_MEM)
#include <sycl/ext/oneapi/backend/level_zero.hpp>
#include <level_zero/ze_api.h>

namespace oneapi_vm = ::sycl::ext::oneapi::experimental;
#endif

namespace olb {

namespace gpu {

namespace sycl {

template<typename T>
struct Column<T>::Data {
  gpu::sycl::device::unique_ptr<T> deviceBuffer;
  std::unique_ptr<T[]> hostBuffer;
  std::size_t hostCapacity = 0;
  std::size_t deviceCapacity = 0;
};

template<typename T>
Column<T>::Column(std::size_t count):
  _count(count),
  _data(std::make_unique<Data>())
{
  if (count > 0) {
    _data->deviceBuffer = gpu::sycl::device::malloc<T>(count);
    _data->hostBuffer = std::make_unique<T[]>(count);
    _data->deviceCapacity = count;
    _data->hostCapacity = count;
  }
}

template<typename T>
Column<T>::Column(Column<T>&& rhs):
  _count(rhs._count),
  _data(std::move(rhs._data))
{
  rhs._count = 0;
}

template<typename T>
Column<T>::Column(const Column<T>& rhs):
  _count(rhs._count),
  _data(std::make_unique<Data>())
{
  if (_count > 0) {
    _data->deviceBuffer = gpu::sycl::device::malloc<T>(_count);
    _data->hostBuffer = std::make_unique<T[]>(_count);
    _data->deviceCapacity = _count;
    _data->hostCapacity = _count;
    std::copy(rhs._data->hostBuffer.get(), rhs._data->hostBuffer.get() + _count, _data->hostBuffer.get());
    gpu::sycl::device::copyToDevice(rhs.deviceData(), deviceData(), _count * sizeof(T));
  }
}

template<typename T>
Column<T>::~Column()
{ }

template<typename T>
const T& Column<T>::operator[](std::size_t i) const
{
  return _data->hostBuffer[i];
}

template<typename T>
T& Column<T>::operator[](std::size_t i)
{
  return _data->hostBuffer[i];
}

template<typename T>
std::size_t Column<T>::size() const
{
  return _count;
}

template<typename T>
const T* Column<T>::data() const
{
  return _data->hostBuffer.get();
}

template<typename T>
T* Column<T>::data()
{
  return _data->hostBuffer.get();
}

template<typename T>
const T* Column<T>::deviceData() const
{
  return _data->deviceBuffer.get();
}

template<typename T>
T* Column<T>::deviceData()
{
  return _data->deviceBuffer.get();
}

template<typename T>
void Column<T>::clear()
{
  resize(0);
}

template<typename T>
void Column<T>::resize(std::size_t newCount)
{
  if (newCount == _count) {
    return;
  }
  if (newCount > _data->hostCapacity) {
    auto newBuffer = std::make_unique<T[]>(newCount);
    if (_count > 0) {
      std::copy(_data->hostBuffer.get(), _data->hostBuffer.get() + _count, newBuffer.get());
    }
    _data->hostBuffer = std::move(newBuffer);
    _data->hostCapacity = newCount;
  }
  _count = newCount;
}

template<typename T>
void Column<T>::push_back(T value)
{
  if (_count == _data->hostCapacity) {
    std::size_t newCapacity = _data->hostCapacity == 0 ? 1 : _data->hostCapacity * 2;
    auto newBuffer = std::make_unique<T[]>(newCapacity);
    if (_count > 0) {
      std::copy(_data->hostBuffer.get(), _data->hostBuffer.get() + _count, newBuffer.get());
    }
    _data->hostBuffer = std::move(newBuffer);
    _data->hostCapacity = newCapacity;
  }
  _data->hostBuffer[_count] = std::move(value);
  _count++;
}

template<typename T>
void Column<T>::deduplicate()
{
  if (_count <= 1) {
    return;
  }
  std::sort(_data->hostBuffer.get(), _data->hostBuffer.get() + _count);
  auto newEnd = std::unique(_data->hostBuffer.get(), _data->hostBuffer.get() + _count);
  _count = std::distance(_data->hostBuffer.get(), newEnd);
}

template<typename T>
void Column<T>::setProcessingContext(ProcessingContext context)
{
  if (_count != _data->deviceCapacity) {
    auto newBuffer = gpu::sycl::device::malloc<T>(_count);
    if (_count > 0 && _data->deviceCapacity > 0) {
      std::size_t copyCount = std::min(_data->deviceCapacity, _count);
      gpu::sycl::device::copyToDevice(_data->deviceBuffer.get(), newBuffer, copyCount * sizeof(T));
    }
    _data->deviceBuffer = std::move(newBuffer);
    _data->deviceCapacity = _count;
  }

  if (_count == 0) {
    return;
  }

  switch (context) {
  case ProcessingContext::Evaluation:
    gpu::sycl::device::copyToHost(deviceData(), _data->hostBuffer.get(), _count * sizeof(T));
    return;
  case ProcessingContext::Simulation:
    gpu::sycl::device::copyToDevice(_data->hostBuffer.get(), deviceData(), _count * sizeof(T));
    return;
  }
}

template<typename T>
void Column<T>::setProcessingContext(ProcessingContext context, device::Stream&)
{
  setProcessingContext(context);
}

template<typename T>
std::size_t Column<T>::getNblock() const
{
  return 1;
}

template<typename T>
std::size_t Column<T>::getSerializableSize() const
{
  return _count * sizeof(T);
}

template<typename T>
bool* Column<T>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  return nullptr;
}

template<typename T>
struct CyclicColumn<T>::Data {
  std::unique_ptr<T[]> hostBuffer;

#if defined(SYCL_EXT_ONEAPI_VIRTUAL_MEM)
  //std::optional<oneapi_vm::physical_mem> physMem;
  //std::uintptr_t vPtr = 0;
  //std::size_t vSize = 0;

  gpu::sycl::device::unique_ptr<T> deviceBuffer;
#endif

  //int shm_file = -1;
  //std::uint8_t* cpuVmmBuffer = nullptr;

  //bool oneapiVmmEnabled = false;
};

template<typename T>
CyclicColumn<T>::CyclicColumn(std::size_t count):
  //_count(device::getPageAlignedCount<T>(count)),
  _count(count),
  _size(_count * sizeof(T)),
  _data(std::make_unique<Data>()),
  _shift(0)
{
  _data->hostBuffer = std::make_unique<T[]>(_count);
  _data->deviceBuffer = gpu::sycl::device::malloc<T>(_count);

//#if defined(SYCL_EXT_ONEAPI_VIRTUAL_MEM)
  _deviceBase = _data->deviceBuffer.get();
  _devicePopulation = _deviceBase;

  //_data->oneapiVmmEnabled = true;

  //// In the future once SYCL allows resp. hardware supports mapping the same physical buffer multiple times
  //try {
  //  auto& queue = gpu::sycl::device::getQueue();
  //  const auto& context = queue.get_context();
  //  const auto& device = queue.get_device();

  //  _data->vPtr = oneapi_vm::reserve_virtual_mem(2 * _size, context);

  //  oneapi_vm::physical_mem physMem(device, context, _size);
  //  physMem.map(_data->vPtr,         _size, oneapi_vm::address_access_mode::read_write);
  //  physMem.map(_data->vPtr + _size, _size, oneapi_vm::address_access_mode::read_write);

  //  _data->physMem = std::move(physMem);
  //  _data->oneapiVmmEnabled = true;

  //  _deviceBase = reinterpret_cast<T*>(_data->vPtr);
  //  _devicePopulation = _deviceBase;

  //  auto zeContext = ::sycl::get_native<::sycl::backend::ext_oneapi_level_zero>(context);
  //  auto zeDevice  = ::sycl::get_native<::sycl::backend::ext_oneapi_level_zero>(device);

  //  ze_result_t status = zeContextMakeMemoryResident(
  //    zeContext,
  //    zeDevice,
  //    reinterpret_cast<void*>(_data->vPtr),
  //    _size
  //  );

  //  if (status != ZE_RESULT_SUCCESS) {
  //    throw std::runtime_error("Level Zero failed to make VMM resident!");
  //  }

  //  status = zeContextMakeMemoryResident(
  //    zeContext,
  //    zeDevice,
  //    reinterpret_cast<void*>(_data->vPtr + _size),
  //    _size
  //  );

  //  if (status != ZE_RESULT_SUCCESS) {
  //    throw std::runtime_error("Level Zero failed to make VMM resident!");
  //  }

  //  auto ze_driver = ::sycl::get_native<::sycl::backend::ext_oneapi_level_zero>(device.get_platform());
  //  typedef ze_result_t (*zex_driver_import_external_pointer_t)(ze_driver_handle_t, void*, size_t);
  //  zex_driver_import_external_pointer_t zexDriverImportExternalPointer = nullptr;
  //  zeDriverGetExtensionFunctionAddress(ze_driver,
  //                                      "zexDriverImportExternalPointer",
  //                                      reinterpret_cast<void**>(&zexDriverImportExternalPointer));
  //  if (zexDriverImportExternalPointer) {
  //    ze_result_t res1 = zexDriverImportExternalPointer(ze_driver, _deviceBase, _size);
  //    if (res1 != ZE_RESULT_SUCCESS) {
  //         std::cerr << "Warning: Failed to import first VMM range to driver tracking" << std::endl;
  //    }
  //    ze_result_t res2 = zexDriverImportExternalPointer(ze_driver, _deviceBase + _count, _size);
  //    if (res2 != ZE_RESULT_SUCCESS) {
  //         std::cerr << "Warning: Failed to import second VMM range to driver tracking" << std::endl;
  //    }
  //  }

  //  {
  //    T* base = _deviceBase;
  //    auto baseCount = _count;
  //    queue.submit([=](::sycl::handler& h) {
  //      h.single_task([=]() {
  //        T test1 = base[0]; // Access the first float of the VMM range
  //        T test2 = base[baseCount]; // Access the first float of the VMM range
  //        T test3 = base[int(1.5*baseCount)]; // Access the first float of the VMM range
  //      });
  //    }).wait_and_throw();
  //    std::cout << "VMM Hardware PDE Link: SUCCESS" << std::endl;
  //  }
  //} catch (const ::sycl::exception& e) {
  //  std::cerr << "VMM Exception: " << e.what() << std::endl;
  //  _data->oneapiVmmEnabled = false;
  //}
//#endif

//  if (!_data->oneapiVmmEnabled) {
//#ifdef __NR_memfd_create
//    _data->shm_file = syscall(__NR_memfd_create, "openlb", MFD_CLOEXEC);
//#else
//    std::string shm_path = "/openlb_block_XXXXXX";
//    const int shm_name = mkstemp(const_cast<char*>(shm_path.data()));
//    if (shm_name != -1) {
//      throw std::runtime_error("Could not generate unique shared memory object name");
//    }
//    _data->shm_file = shm_open(shm_path.c_str(), O_CREAT | O_RDWR | O_EXCL | O_CLOEXEC, S_IRUSR | S_IWUSR);
//    shm_unlink(shm_path.c_str());
//#endif
//    if (_data->shm_file == -1) {
//      throw std::runtime_error("Failed to create shared memory object");
//    }
//
//    if (ftruncate(_data->shm_file, _size) == -1) {
//      throw std::runtime_error("Failed to resize shared memory object");
//    }
//
//    _data->cpuVmmBuffer = static_cast<std::uint8_t*>(
//      //mmap(NULL, 2 * _size, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0));
//      mmap(NULL, _size, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0));
//
//    const int PROT_RW = PROT_READ | PROT_WRITE;
//    mmap(_data->cpuVmmBuffer,         _size, PROT_RW, MAP_SHARED | MAP_FIXED, _data->shm_file, 0);
//    //mmap(_data->cpuVmmBuffer + _size, _size, PROT_RW, MAP_SHARED | MAP_FIXED, _data->shm_file, 0);
//
//    _deviceBase = reinterpret_cast<T*>(_data->cpuVmmBuffer);
//    _devicePopulation = _deviceBase;
//  }
}

template<typename T>
CyclicColumn<T>::~CyclicColumn()
{
//#if defined(SYCL_EXT_ONEAPI_VIRTUAL_MEM)
//  if (_data->oneapiVmmEnabled) {
//    //const auto& context = gpu::sycl::device::getQueue().get_context();
//    //oneapi_vm::free_virtual_mem(_data->vPtr, 2 * _size, context);
//    return;
//  }
//#endif
//  if (_data->cpuVmmBuffer) {
//    munmap(_data->cpuVmmBuffer, 2 * _size);
//    close(_data->shm_file);
//  }
}

template<typename T>
void CyclicColumn<T>::rotate(std::ptrdiff_t offset)
{
  _shift -= offset;
  if (_shift >= _count) {
    _shift -= _count;
  }
  else if (_shift < 0) {
    _shift += _count;
  }
  refresh();
}

template<typename T>
CyclicColumn<T>::Accessor CyclicColumn<T>::getHostAccessor()
{
  return Accessor{
    .start0    = _data->hostBuffer.get() + _shift,
    .start1    = _data->hostBuffer.get() - (_count - _shift),
    .remainder = _count - _shift - 1
  };
}

template<typename T>
CyclicColumn<T>::Accessor CyclicColumn<T>::getDeviceAccessor()
{
  return Accessor{
    .start0    = deviceData() + _shift,
    .start1    = deviceData() - (_count - _shift),
    .remainder = _count - _shift - 1
  };
}

template<typename T>
void CyclicColumn<T>::refresh()
{
  //_devicePopulation = _deviceBase + _shift;
}

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
void CyclicColumn<T>::setProcessingContext(ProcessingContext context)
{
  switch (context) {
  case ProcessingContext::Evaluation:
    gpu::sycl::device::copyToHost(_devicePopulation, _data->hostBuffer.get(), _count * sizeof(T));
    return;
  case ProcessingContext::Simulation:
    gpu::sycl::device::copyToDevice(_data->hostBuffer.get(), _devicePopulation, _count * sizeof(T));
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
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, *_data->hostBuffer.get(), _count);

  return dataPtr;
}

} // namespace sycl
} // namespace gpu

template <typename T>
std::size_t ConcreteCommunicatable<gpu::sycl::Column<T>>::serialize(ConstSpan<CellID> indices,
                                                                    std::uint8_t* buffer) const
{
  const auto* deviceD = _column.deviceData();
  auto* bufferD = reinterpret_cast<T*>(buffer);
  const auto* indicesD = indices.data();
  gpu::sycl::device::getQueue().parallel_for(sycl::range<1>(indices.size()), [=](sycl::id<1> i) {
    bufferD[i] = deviceD[indicesD[i]];
  }).wait();
  return indices.size() * sizeof(T);
}

template <typename T>
std::size_t ConcreteCommunicatable<gpu::sycl::Column<T>>::deserialize(ConstSpan<CellID> indices,
                                                                      const std::uint8_t* buffer)
{
  auto* deviceD = _column.deviceData();
  const auto* bufferD = reinterpret_cast<const T*>(buffer);
  const auto* indicesD = indices.data();
  gpu::sycl::device::getQueue().parallel_for(sycl::range<1>(indices.size()), [=](sycl::id<1> i) {
    deviceD[indicesD[i]] = bufferD[i];
  }).wait();
  return indices.size() * sizeof(T);
}

template <typename T>
std::size_t ConcreteCommunicatable<gpu::sycl::CyclicColumn<T>>::serialize(ConstSpan<CellID> indices,
                                                                          std::uint8_t* buffer) const
{
  // Index through the shift-aware accessor: unlike the CUDA backend, SYCL
  // does not fold the cyclic shift into deviceData() (refresh() is a no-op).
  const auto column = _column.getDeviceAccessor();
  auto* bufferD = reinterpret_cast<T*>(buffer);
  const auto* indicesD = indices.data();
  gpu::sycl::device::getQueue().parallel_for(sycl::range<1>(indices.size()), [=](sycl::id<1> i) {
    bufferD[i] = column[indicesD[i]];
  }).wait();
  return indices.size() * sizeof(T);
}

template <typename T>
std::size_t ConcreteCommunicatable<gpu::sycl::CyclicColumn<T>>::deserialize(ConstSpan<CellID> indices,
                                                                            const std::uint8_t* buffer)
{
  const auto column = _column.getDeviceAccessor();
  const auto* bufferD = reinterpret_cast<const T*>(buffer);
  const auto* indicesD = indices.data();
  gpu::sycl::device::getQueue().parallel_for(sycl::range<1>(indices.size()), [=](sycl::id<1> i) {
    column[indicesD[i]] = bufferD[i];
  }).wait();
  return indices.size() * sizeof(T);
}

} // namespace olb

#include "sparseColumn.hh"

#endif
