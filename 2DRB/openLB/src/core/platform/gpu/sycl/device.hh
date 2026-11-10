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

#ifndef GPU_SYCL_DEVICE_HH
#define GPU_SYCL_DEVICE_HH

#include "device.h"

#include "utilities/performanceProfiler.h"
#include <stdexcept>
#include <iostream>
#include <optional>

#include <sycl/sycl.hpp>

#include "core/concepts.h"

#if defined(SYCL_EXT_ONEAPI_VIRTUAL_MEM)
namespace oneapi_vm = ::sycl::ext::oneapi::experimental;
#endif

namespace olb {

namespace gpu {

namespace sycl {

namespace queues {

struct OuterCollision { };
struct InnerCollision { };

}

namespace device {

inline ::sycl::device& getDevice() {
  static ::sycl::device device(::sycl::default_selector_v);
  return device;
}

inline ::sycl::context& getContext() {
  static ::sycl::context context(getDevice());
  return context;
}

// Singleton queue for global operations
inline ::sycl::queue& getQueue() {
  static ::sycl::queue queue(getContext(), getDevice());
  return queue;
}

// Singleton queue for subsets of global operations
template <typename QUEUE>
inline ::sycl::queue& getQueue() {
  static ::sycl::queue queue(getContext(), getDevice(),
                             ::sycl::property_list{::sycl::ext::oneapi::property::queue::priority_normal()});
  return queue;
}

template <>
inline ::sycl::queue& getQueue<queues::OuterCollision>() {
  static ::sycl::queue queue(getContext(), getDevice(),
                             ::sycl::property_list{::sycl::ext::oneapi::property::queue::priority_high()});
  return queue;
}

template <>
inline ::sycl::queue& getQueue<queues::InnerCollision>() {
  static ::sycl::queue queue(getContext(), getDevice(),
                             ::sycl::property_list{::sycl::ext::oneapi::property::queue::priority_low()});
  return queue;
}

inline int getCount() {
  try {
    return ::sycl::device::get_devices(::sycl::info::device_type::gpu).size();
  } catch (...) {
    return 0;
  }
}

inline void check() {
  getQueue().throw_asynchronous();
}

void synchronize() {
  OLB_PROFILE("GPU_SYCL:Dev:SyncAll");
  getQueue().wait_and_throw();
}

void synchronizeInnerCollision() {
  OLB_PROFILE("GPU_SYCL:Dev:SyncIn");
  getQueue<queues::InnerCollision>().wait_and_throw();
}

inline int get() {
  return 0; // Default
}

inline void copyToHost(const void* src, void* dst, std::size_t count) {
  OLB_PROFILE("GPU_SYCL:Mem:DtoH");
  getQueue().memcpy(dst, src, count).wait();
}

inline void copyToDevice(const void* src, void* dst, std::size_t count) {
  OLB_PROFILE("GPU_SYCL:Mem:HtoD");
  getQueue().memcpy(dst, src, count).wait();
}

template <typename T>
T* malloc(std::size_t size) {
  return static_cast<T*>(::sycl::malloc_device(size * sizeof(T), getQueue()));
}

inline std::size_t getDevicePageSize() {
//#if defined(SYCL_EXT_ONEAPI_VIRTUAL_MEM)
//  auto& queue = getQueue();
//  const auto& context = queue.get_context();
//  const auto& device = queue.get_device();
//  std::size_t granularity = oneapi_vm::get_mem_granularity(device, context);
//  return std::max<std::size_t>(granularity, sysconf(_SC_PAGESIZE));
//#else
//  return sysconf(_SC_PAGESIZE);
//#endif
  return 2097152;
}

std::size_t getDeviceMemorySize() {
  return getDevice().get_info<::sycl::info::device::global_mem_size>();
}

template <typename T>
unique_ptr<T>::~unique_ptr() {
  if (_ptr != nullptr) {
    ::sycl::free(_ptr, getQueue());
  }
}

template <typename T>
unique_ptr<T>& unique_ptr<T>::operator=(unique_ptr<T>&& rhs) {
  if (_ptr != nullptr) {
    ::sycl::free(_ptr, getQueue());
  }
  _ptr = rhs.release();
  return *this;
}

struct Stream::Impl {
  ::sycl::queue queue;

  Impl(unsigned int)
  : queue(::sycl::default_selector_v) {}
};

Stream::Stream(unsigned int flags):
  _impl(new Impl(flags))
{ }

Stream::~Stream() {
  delete _impl;
}

void Stream::synchronize() {
  _impl->queue.wait();
}

void* Stream::getNativeHandle() {
  return &_impl->queue;
}

inline void asyncCopyToHost(Stream& stream, void* src, void* dst, std::size_t count) {
  static_cast<::sycl::queue*>(stream.getNativeHandle())->memcpy(dst, src, count);
}

inline void asyncCopyToDevice(Stream& stream, void* src, void* dst, std::size_t count) {
  static_cast<::sycl::queue*>(stream.getNativeHandle())->memcpy(dst, src, count);
}

}

}

}

struct Event<Platform::GPU_SYCL>::Impl {
  std::optional<::sycl::event> event;
};

Event<Platform::GPU_SYCL>::Event():
  _impl{new Impl{}}
{ }

Event<Platform::GPU_SYCL>::~Event() {
  delete _impl;
}

template <typename QUEUE>
void Event<Platform::GPU_SYCL>::record(QUEUE) {
  _impl->event = gpu::sycl::device::getQueue<QUEUE>().ext_oneapi_submit_barrier();
}

template <typename QUEUE>
void Event<Platform::GPU_SYCL>::barrierOn(QUEUE) const {
  if (_impl->event) {
    gpu::sycl::device::getQueue<QUEUE>().ext_oneapi_submit_barrier({*_impl->event});
  }
}

void Event<Platform::GPU_SYCL>::wait() const {
  if (_impl->event) {
    _impl->event->wait_and_throw();
  }
}

namespace gpu {

namespace sycl {

namespace device {

/// Register the last recorded event as dependency of a native command group
void dependsOn(::sycl::handler& cgh, const Event<Platform::GPU_SYCL>& event) {
  if (event.impl().event) {
    cgh.depends_on(*event.impl().event);
  }
}

}

}

}

static_assert(concepts::Completable<Event<Platform::GPU_SYCL>>,
              "Event<GPU_SYCL> must model a completion token");
static_assert(concepts::DeviceOrderable<Event<Platform::GPU_SYCL>,
                                        gpu::sycl::queues::OuterCollision>,
              "Event<GPU_SYCL> must order device work w.r.t. queue tags");
static_assert(concepts::Completable<::sycl::event>,
              "native SYCL events share the host-side waiting vocabulary");

}

#endif
