/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Adrian Kummerlaender
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

#ifndef GPU_CUDA_DEVICE_HH
#define GPU_CUDA_DEVICE_HH

#include "device.h"

#include "core/concepts.h"
#include "utilities/performanceProfiler.h"

#include <stdexcept>

#include <cuda.h>

namespace olb {

namespace gpu {

namespace cuda {

namespace device {

int getCount() {
  int devices{};
  cudaGetDeviceCount(&devices);
  return devices;
}

void check() {
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess) {
    throw std::runtime_error(cudaGetErrorString(error));
  }
}

/// Check CUDA driver errors
void check(CUresult result) {
  if (result != CUDA_SUCCESS) {
    const char* description{};
    cuGetErrorString(result, &description);
    throw std::runtime_error(std::string(description));
  }
}

void synchronize() {
  if (getCount() > 0) {
    OLB_PROFILE("GPU_CUDA:Dev:Sync");
    cudaDeviceSynchronize();
    check();
  }
}

int get() {
  int device{};
  cudaGetDevice(&device);
  check();
  return device;
}

void copyToHost(void* src, void* dst, std::size_t count) {
  OLB_PROFILE("GPU_CUDA:Mem:DtoH");
  cudaMemcpy(dst, src, count, cudaMemcpyDeviceToHost);
  check();
}

void copyToDevice(void* src, void* dst, std::size_t count) {
  OLB_PROFILE("GPU_CUDA:Mem:HtoD");
  cudaMemcpy(dst, src, count, cudaMemcpyHostToDevice);
  check();
}

template <typename T>
T* malloc(std::size_t size) {
  T* ptr{};
  cudaMalloc(&ptr, size*sizeof(T));
  check();
  return ptr;
}

std::size_t getDevicePageSize() {
  std::size_t granularity = 0;
  CUmemAllocationProp prop = {};
  prop.type = CU_MEM_ALLOCATION_TYPE_PINNED;
  prop.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
  prop.location.id = device::get();
  check(cuMemGetAllocationGranularity(&granularity, &prop, CU_MEM_ALLOC_GRANULARITY_MINIMUM));
  return granularity;
}

std::size_t getDeviceMemorySize() {
  std::size_t freeBytes = 0;
  std::size_t totalBytes = 0;
  cudaMemGetInfo(&freeBytes, &totalBytes);
  check();
  return totalBytes;
}


template <typename T>
unique_ptr<T>::~unique_ptr() {
  if (_ptr != nullptr) {
    cudaFree(_ptr);
    check();
  }
}

template <typename T>
unique_ptr<T>& unique_ptr<T>::operator=(unique_ptr<T>&& rhs) {
  if (_ptr != nullptr) {
    cudaFree(_ptr);
    check();
  }
  _ptr = rhs.release();
  return *this;
}


Stream::Stream(unsigned int flags) {
  cudaStreamCreateWithFlags(&_stream, flags);
  check();
}

Stream::Stream(unsigned int flags, int priority) {
  cudaStreamCreateWithPriority(&_stream, flags, priority);
  check();
}

Stream::~Stream() {
  cudaStreamDestroy(_stream);
  check();
}

void Stream::synchronize() {
  cudaStreamSynchronize(_stream);
  check();
}

void asyncCopyToHost(Stream& stream, void* src, void* dst, std::size_t count) {
  cudaMemcpyAsync(dst, src, count, cudaMemcpyDeviceToHost, stream.get());
}

void asyncCopyToDevice(Stream& stream, void* src, void* dst, std::size_t count) {
  cudaMemcpyAsync(dst, src, count, cudaMemcpyHostToDevice, stream.get());
}

template <>
Stream& getStream<streams::OuterCollision>() {
  int leastPriority = 0, greatestPriority = 0;
  cudaDeviceGetStreamPriorityRange(&leastPriority, &greatestPriority);
  static Stream stream(cudaStreamNonBlocking, greatestPriority);
  return stream;
}

template <>
Stream& getStream<streams::InnerCollision>() {
  int leastPriority = 0, greatestPriority = 0;
  cudaDeviceGetStreamPriorityRange(&leastPriority, &greatestPriority);
  static Stream stream(cudaStreamNonBlocking, leastPriority);
  return stream;
}

template <>
Stream& getStream<streams::PostProcessor>() {
  int leastPriority = 0, greatestPriority = 0;
  cudaDeviceGetStreamPriorityRange(&leastPriority, &greatestPriority);
  static Stream stream(cudaStreamNonBlocking, greatestPriority);
  return stream;
}

void synchronizeOuterCollision() {
  OLB_PROFILE("GPU_CUDA:Dev:SyncOuter");
  getStream<streams::OuterCollision>().synchronize();
}

void synchronizeInnerCollision() {
  OLB_PROFILE("GPU_CUDA:Dev:SyncInner");
  getStream<streams::InnerCollision>().synchronize();
}

void synchronizePostProcessor() {
  OLB_PROFILE("GPU_CUDA:Dev:SyncPP");
  getStream<streams::PostProcessor>().synchronize();
}

}

}

}

Event<Platform::GPU_CUDA>::Event() {
  cudaEventCreateWithFlags(&_event, cudaEventDisableTiming);
  gpu::cuda::device::check();
}

Event<Platform::GPU_CUDA>::~Event() {
  cudaEventDestroy(_event);
  gpu::cuda::device::check();
}

void Event<Platform::GPU_CUDA>::record(gpu::cuda::device::Stream& stream) {
  record(stream.get());
}

void Event<Platform::GPU_CUDA>::barrierOn(gpu::cuda::device::Stream& stream) const {
  barrierOn(stream.get());
}

void Event<Platform::GPU_CUDA>::record(cudaStream_t stream) {
  cudaEventRecord(_event, stream);
  gpu::cuda::device::check();
}

void Event<Platform::GPU_CUDA>::barrierOn(cudaStream_t stream) const {
  cudaStreamWaitEvent(stream, _event, 0);
  gpu::cuda::device::check();
}

void Event<Platform::GPU_CUDA>::wait() const {
  cudaEventSynchronize(_event);
  gpu::cuda::device::check();
}

static_assert(concepts::Completable<Event<Platform::GPU_CUDA>>,
              "Event<GPU_CUDA> must model a completion token");
static_assert(concepts::DeviceOrderable<Event<Platform::GPU_CUDA>,
                                        gpu::cuda::device::Stream>,
              "Event<GPU_CUDA> must order device work w.r.t. streams");

}

#endif
