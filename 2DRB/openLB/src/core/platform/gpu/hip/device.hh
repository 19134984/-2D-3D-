/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022, 2025 Adrian Kummerlaender, Leonardo Dorneles
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

#ifndef GPU_HIP_DEVICE_HH
#define GPU_HIP_DEVICE_HH

#include "device.h"

#include "core/concepts.h"
#include "utilities/performanceProfiler.h"

#include <stdexcept>

namespace olb {

namespace gpu {

namespace hip {

namespace device {

int getCount() {
  int devices{};
  HIP_CHECK(hipGetDeviceCount(&devices));
  return devices;
}

void check() {
  hipError_t error = hipGetLastError();
  if (error != hipSuccess) {
    throw std::runtime_error(hipGetErrorString(error));
  }
}

void synchronize() {
  if (getCount() > 0) {
    OLB_PROFILE("GPU_HIP:Dev:Sync");
    HIP_CHECK(hipDeviceSynchronize());
  }
}

int get() {
  int device{};
  HIP_CHECK(hipGetDevice(&device));
  return device;
}

void copyToHost(void* src, void* dst, std::size_t count) {
  OLB_PROFILE("GPU_HIP:Mem:DtoH");
  HIP_CHECK(hipMemcpy(dst, src, count, hipMemcpyDeviceToHost));
}

void copyToDevice(void* src, void* dst, std::size_t count) {
  OLB_PROFILE("GPU_HIP:Mem:HtoD");
  HIP_CHECK(hipMemcpy(dst, src, count, hipMemcpyHostToDevice));
}

template <typename T>
T* malloc(std::size_t size) {
  T* ptr{};
  HIP_CHECK(hipMalloc(&ptr, size*sizeof(T)));
  return ptr;
}

std::size_t getDevicePageSize() {
  std::size_t granularity = 0;
  hipMemAllocationProp prop = {};
  prop.type = hipMemAllocationTypePinned;
  prop.location.type = hipMemLocationTypeDevice;
  prop.location.id = device::get();
  HIP_CHECK(hipMemGetAllocationGranularity(&granularity, &prop, hipMemAllocationGranularityMinimum));
  return granularity;
}

std::size_t getDeviceMemorySize() {
  std::size_t freeBytes = 0;
  std::size_t totalBytes = 0;
  HIP_CHECK(hipMemGetInfo(&freeBytes, &totalBytes));
  return totalBytes;
}


template <typename T>
unique_ptr<T>::~unique_ptr() {
  if (_ptr != nullptr) {
    HIP_CHECK(hipFree(_ptr));
  }
}

template <typename T>
unique_ptr<T>& unique_ptr<T>::operator=(unique_ptr<T>&& rhs) {
  if (_ptr != nullptr) {
    HIP_CHECK(hipFree(_ptr));
  }
  _ptr = rhs.release();
  return *this;
}


Stream::Stream(unsigned int flags) {
  HIP_CHECK(hipStreamCreateWithFlags(&_stream, flags));
}

Stream::Stream(unsigned int flags, int priority) {
  HIP_CHECK(hipStreamCreateWithPriority(&_stream, flags, priority));
}

Stream::~Stream() {
  HIP_CHECK(hipStreamDestroy(_stream));

}

void Stream::synchronize() {
  HIP_CHECK(hipStreamSynchronize(_stream));
}

void asyncCopyToHost(Stream& stream, void* src, void* dst, std::size_t count) {
  HIP_CHECK(hipMemcpyAsync(dst, src, count, hipMemcpyDeviceToHost, stream.get()));
}

void asyncCopyToDevice(Stream& stream, void* src, void* dst, std::size_t count) {
  HIP_CHECK(hipMemcpyAsync(dst, src, count, hipMemcpyHostToDevice, stream.get()));
}

template <>
inline Stream& getStream<streams::OuterCollision>() {
  int leastPriority = 0, greatestPriority = 0;
  HIP_CHECK(hipDeviceGetStreamPriorityRange(&leastPriority, &greatestPriority));

  //int mediumPriority = (greatestPriority+leastPriority)/2;
  static Stream stream(hipStreamNonBlocking, greatestPriority);
  return stream;
}

template <>
inline Stream& getStream<streams::InnerCollision>() {
  int leastPriority = 0, greatestPriority = 0;
  HIP_CHECK(hipDeviceGetStreamPriorityRange(&leastPriority, &greatestPriority));
  static Stream stream(hipStreamNonBlocking, leastPriority);
  return stream;
}

template <>
inline Stream& getStream<streams::PostProcessor>() {
  int leastPriority = 0, greatestPriority = 0;
  HIP_CHECK(hipDeviceGetStreamPriorityRange(&leastPriority, &greatestPriority));
  static Stream stream(hipStreamNonBlocking, greatestPriority);
  return stream;
}



inline void synchronizeInnerCollision(){
  HIP_CHECK(hipStreamSynchronize(gpu::hip::device::getStream<streams::InnerCollision>().get()));
}

inline void synchronizeOuterCollision(){
  HIP_CHECK(hipStreamSynchronize(gpu::hip::device::getStream<streams::OuterCollision>().get()));
}

inline void synchronizePostProcessor(){
   HIP_CHECK(hipStreamSynchronize(gpu::hip::device::getStream<streams::PostProcessor>().get()));
}


} // device

}

}

Event<Platform::GPU_HIP>::Event() {
  HIP_CHECK(hipEventCreateWithFlags(&_event, hipEventDisableTiming));
}

Event<Platform::GPU_HIP>::~Event() {
  HIP_CHECK(hipEventDestroy(_event));
}

void Event<Platform::GPU_HIP>::record(gpu::hip::device::Stream& stream) {
  record(stream.get());
}

void Event<Platform::GPU_HIP>::barrierOn(gpu::hip::device::Stream& stream) const {
  barrierOn(stream.get());
}

void Event<Platform::GPU_HIP>::record(hipStream_t stream) {
  HIP_CHECK(hipEventRecord(_event, stream));
}

void Event<Platform::GPU_HIP>::barrierOn(hipStream_t stream) const {
  HIP_CHECK(hipStreamWaitEvent(stream, _event, 0));
}

void Event<Platform::GPU_HIP>::wait() const {
  HIP_CHECK(hipEventSynchronize(_event));
}

static_assert(concepts::Completable<Event<Platform::GPU_HIP>>,
              "Event<GPU_HIP> must model a completion token");
static_assert(concepts::DeviceOrderable<Event<Platform::GPU_HIP>,
                                        gpu::hip::device::Stream>,
              "Event<GPU_HIP> must order device work w.r.t. streams");

}

#endif
