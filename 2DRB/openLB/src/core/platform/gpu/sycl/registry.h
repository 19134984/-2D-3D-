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

#ifndef GPU_SYCL_REGISTRY_H
#define GPU_SYCL_REGISTRY_H

#include "utilities/typeIndexedContainers.h"
#include "core/platform/gpu/sparseAccessor.h"

#include <map>
#include <memory>
#include <cstdint>
#include <cstddef>
#include <string_view>

// The canary is device-side state. Enable it only in the precompiled SYCL translation unit
// (compiled with -fsycl -> SYCL_LANGUAGE_VERSION). App TUs include this header without -fsycl.
#if defined(FEATURE_SPARSE_SCRATCH_CANARY) && defined(SYCL_LANGUAGE_VERSION)
#include <sycl/sycl.hpp>
#define OLB_GPU_SYCL_SPARSE_CANARY
#endif

namespace olb {

namespace gpu {

namespace sycl {

#ifdef OLB_GPU_SYCL_SPARSE_CANARY
/// Counts sparse halo scatters that fell into the discard scratch tile (unallocated inbound cell)
inline ::sycl::ext::oneapi::experimental::device_global<unsigned long long> sparse_scratch_writes;
#endif

/// Type-erased pointer to FieldArrayD device data
struct AnyDeviceFieldArrayD {
  void**   data;
  unsigned column_count;
  unsigned element_size;
  /// 0 for a dense field, else the sparse TILE size
  unsigned tile;
  /// byte stride between per-dimension sparse accessors in `data` (sizeof(DeviceSparseColumn))
  unsigned accessor_stride;

  std::uint8_t* operator[](unsigned iD) const {
    return reinterpret_cast<std::uint8_t*>(data[iD]);
  }

  /// Gather component iD of cell into dst (element_size bytes)
  void gather(unsigned iD, std::size_t cell, std::uint8_t* dst) const {
    if (tile == 0) {
      const std::uint8_t* src = reinterpret_cast<std::uint8_t*>(data[iD]) + cell*element_size;
      for (unsigned k=0; k < element_size; ++k) {
        dst[k] = src[k];
      }
    } else {
      const std::uint8_t* acc = reinterpret_cast<const std::uint8_t*>(data) + iD*accessor_stride;
      const std::int32_t* dir     = *reinterpret_cast<const std::int32_t* const*>(acc);
      const std::uint8_t* payload = *reinterpret_cast<std::uint8_t* const*>(acc + 8);
      if (dir == nullptr) { // dense mode: payload is a flat per-cell buffer
        const std::uint8_t* src = payload + cell*element_size;
        for (unsigned k=0; k < element_size; ++k) {
          dst[k] = src[k];
        }
        return;
      }
      const sparse::TileData* geo = reinterpret_cast<const sparse::TileData*>(acc + 24);
      const unsigned cap = geo->capacity ? geo->capacity : tile;
      const std::int32_t  slot = dir[geo->tileId(cell, tile)];
      const std::uint8_t* src = (slot < 0)
        ? acc + 16
        : payload + (std::size_t(slot)*cap + geo->offset(cell, tile))*element_size;
      for (unsigned k=0; k < element_size; ++k) {
        dst[k] = src[k];
      }
    }
  }

  /// Scatter element_size bytes from src into component iD of cell
  void scatter(unsigned iD, std::size_t cell, const std::uint8_t* src) const {
    if (tile == 0) {
      std::uint8_t* dst = reinterpret_cast<std::uint8_t*>(data[iD]) + cell*element_size;
      for (unsigned k=0; k < element_size; ++k) {
        dst[k] = src[k];
      }
    } else {
      const std::uint8_t* acc = reinterpret_cast<const std::uint8_t*>(data) + iD*accessor_stride;
      const std::int32_t* dir     = *reinterpret_cast<const std::int32_t* const*>(acc);
      std::uint8_t*       payload = *reinterpret_cast<std::uint8_t* const*>(acc + 8);
      if (dir == nullptr) { // dense mode: payload is a flat per-cell buffer
        std::uint8_t* dst = payload + cell*element_size;
        for (unsigned k=0; k < element_size; ++k) {
          dst[k] = src[k];
        }
        return;
      }
      const sparse::TileData* geo = reinterpret_cast<const sparse::TileData*>(acc + 24);
      const unsigned cap = geo->capacity ? geo->capacity : tile;
      const std::int32_t  slot = dir[geo->tileId(cell, tile)];
#ifdef OLB_GPU_SYCL_SPARSE_CANARY
      if (slot < 0) {
        ::sycl::atomic_ref<unsigned long long,
                           ::sycl::memory_order::relaxed,
                           ::sycl::memory_scope::device,
                           ::sycl::access::address_space::global_space> ref(sparse_scratch_writes.get());
        ref.fetch_add(1ULL);
      }
#endif
      const std::size_t base = (slot < 0 ? std::size_t(0) : std::size_t(slot));
      std::uint8_t* dst = payload + (base*cap + geo->offset(cell, tile))*element_size;
      for (unsigned k=0; k < element_size; ++k) {
        dst[k] = src[k];
      }
    }
  }
};

}

}

/// Maintain on-device structure for dynamic field access
template<typename T, typename DESCRIPTOR>
class FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_SYCL> {
private:
  utilities::TypeIndexedMap<AnyFieldTypeD<T,DESCRIPTOR,Platform::GPU_SYCL>*,FieldTypeRegistry> _index;

  struct Data;
  std::unique_ptr<Data> _data;

  bool _modified;

public:
  FieldTypeRegistry();
  ~FieldTypeRegistry();

  template <typename FIELD_TYPE>
  bool provides() const {
    return _index.template provides<FIELD_TYPE>();
  }

  template <typename FIELD_TYPE>
  const AnyFieldTypeD<T,DESCRIPTOR,Platform::GPU_SYCL>* get() const {
    return _index.template get<FIELD_TYPE>();
  }
  template <typename FIELD_TYPE>
  AnyFieldTypeD<T,DESCRIPTOR,Platform::GPU_SYCL>* get() {
    return _index.template get<FIELD_TYPE>();
  }

  template <typename FIELD_TYPE>
  void track(AnyFieldTypeD<T,DESCRIPTOR,Platform::GPU_SYCL>* fieldType);

  void*** deviceData();

  template <typename FIELD>
  void refreshDeviceFieldArray(FieldArrayD<T,DESCRIPTOR,Platform::GPU_SYCL,FIELD>& fieldArray);

  /// Mark FIELD as needing per-step device pointer refresh
  void requestPerStepRefresh(std::type_index field);
  bool needsPerStepRefresh(std::type_index field) const;

  /// Rebuild every sparse field's stored device accessor from its current (post-sync) tile buffers
  void refreshSparseAccessors();

  gpu::sycl::AnyDeviceFieldArrayD deviceFieldArray(std::type_index field);
  std::vector<gpu::sycl::AnyDeviceFieldArrayD> deviceFieldArrays(
    const std::vector<std::type_index>& fields);

};

}

#endif
