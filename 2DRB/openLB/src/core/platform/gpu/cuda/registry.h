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

#ifndef GPU_CUDA_REGISTRY_H
#define GPU_CUDA_REGISTRY_H

#include "utilities/typeIndexedContainers.h"
#include "core/platform/column.h"

#include <map>
#include <memory>
#include <cstdint>
#include <cstddef>

namespace olb {

namespace gpu {

namespace cuda {

/// Virtual interface for device-side dynamically-dispatched dynamics access
template <typename T, typename DESCRIPTOR>
struct Dynamics;

#ifdef FEATURE_SPARSE_SCRATCH_CANARY
extern __device__ unsigned long long sparse_scratch_writes;
#endif

/// Type-erased pointer to FieldArrayD device data
/**
 * Used for dynamic communication buffer (de)serialization
 **/
struct AnyDeviceFieldArrayD {
  void**   data;
  unsigned column_count;
  unsigned element_size;
  /// 0 for a dense field, else the sparse TILE size
  unsigned tile;
  /// byte stride between per-dimension sparse accessors in `data` (sizeof(DeviceSparseColumn))
  unsigned accessor_stride;

  std::uint8_t* operator[](unsigned iD) __device__ {
    return reinterpret_cast<std::uint8_t*>(data[iD]);
  }

  /// Gather component iD of cell into dst (element_size bytes)
  void gather(unsigned iD, std::size_t cell, std::uint8_t* dst) const __device__ {
    if (tile == 0) {
      memcpy(dst,
             reinterpret_cast<std::uint8_t*>(data[iD]) + cell*element_size,
             element_size);
    } else {
      const std::uint8_t* acc = reinterpret_cast<const std::uint8_t*>(data) + iD*accessor_stride;
      const std::int32_t* dir     = *reinterpret_cast<const std::int32_t* const*>(acc);
      const std::uint8_t* payload = *reinterpret_cast<std::uint8_t* const*>(acc + 8);
      if (dir == nullptr) { // dense mode: payload is a flat per-cell buffer
        memcpy(dst, payload + cell*element_size, element_size);
        return;
      }
      const sparse::TileData* geo = reinterpret_cast<const sparse::TileData*>(acc + 24);
      const unsigned cap = geo->capacity ? geo->capacity : tile;
      const std::int32_t  slot    = dir[geo->tileId(cell, tile)];
      if (slot < 0) {
        memcpy(dst, acc + 16, element_size);
      } else {
        memcpy(dst, payload + (std::size_t(slot)*cap + geo->offset(cell, tile))*element_size, element_size);
      }
    }
  }

  /// Scatter element_size bytes from src into component iD of cell
  void scatter(unsigned iD, std::size_t cell, const std::uint8_t* src) const __device__ {
    if (tile == 0) {
      memcpy(reinterpret_cast<std::uint8_t*>(data[iD]) + cell*element_size,
             src,
             element_size);
    } else {
      const std::uint8_t* acc = reinterpret_cast<const std::uint8_t*>(data) + iD*accessor_stride;
      const std::int32_t* dir     = *reinterpret_cast<const std::int32_t* const*>(acc);
      std::uint8_t*       payload = *reinterpret_cast<std::uint8_t* const*>(acc + 8);
      if (dir == nullptr) { // dense mode: payload is a flat per-cell buffer
        memcpy(payload + cell*element_size, src, element_size);
        return;
      }
      const sparse::TileData* geo = reinterpret_cast<const sparse::TileData*>(acc + 24);
      const unsigned cap = geo->capacity ? geo->capacity : tile;
      const std::int32_t  slot    = dir[geo->tileId(cell, tile)];
#if defined(FEATURE_SPARSE_SCRATCH_CANARY) && defined(__CUDA_ARCH__)
      if (slot < 0) {
        atomicAdd(&sparse_scratch_writes, 1ULL);
      }
#endif
      const std::size_t   base    = (slot < 0 ? std::size_t(0) : std::size_t(slot));
      memcpy(payload + (base*cap + geo->offset(cell, tile))*element_size, src, element_size);
    }
  }
};

}

}

/// Maintain on-device structure for dynamic field access
/**
 * Used to enable access to fields that are dynamically allocated
 * instead of being declared by the descriptor.
 **/
template<typename T, typename DESCRIPTOR>
class FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA> {
private:
  /// Host-side field type index
  utilities::TypeIndexedMap<AnyFieldTypeD<T,DESCRIPTOR,Platform::GPU_CUDA>*,FieldTypeRegistry> _index;

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
  const AnyFieldTypeD<T,DESCRIPTOR,Platform::GPU_CUDA>* get() const {
    return _index.template get<FIELD_TYPE>();
  }
  template <typename FIELD_TYPE>
  AnyFieldTypeD<T,DESCRIPTOR,Platform::GPU_CUDA>* get() {
    return _index.template get<FIELD_TYPE>();
  }

  template <typename FIELD_TYPE>
  void track(AnyFieldTypeD<T,DESCRIPTOR,Platform::GPU_CUDA>* fieldType);

  void*** deviceData();

  /// Register a distinct device dynamics object, returning its stable non-zero table index
  std::uint16_t registerDynamics(std::type_index type, void* deviceDynamics);
  /// Device pointer to the per-block dynamics dispatch table (index 0 reserved for nullptr)
  gpu::cuda::Dynamics<T,DESCRIPTOR>** dynamicsTable();

  template <typename FIELD>
  void refreshDeviceFieldArray(FieldArrayD<T,DESCRIPTOR,Platform::GPU_CUDA,FIELD>& fieldArray);

  /// Mark FIELD as needing per-step device pointer refresh in stream()
  /**
   * Only fields read via the type-erased any-fields communication path need
   * their rotated device pointers re-uploaded each step. The common single-
   * field POPULATION propagation uses the specialized typed path (fresh
   * accessor) and never reads these pointers, so its refresh is skipped.
   **/
  void requestPerStepRefresh(std::type_index field);
  bool needsPerStepRefresh(std::type_index field) const;

  /// Rebuild every sparse field's stored device accessor from its current (post-sync) tile buffers
  void refreshSparseAccessors();

  gpu::cuda::AnyDeviceFieldArrayD deviceFieldArray(std::type_index field);
  std::vector<gpu::cuda::AnyDeviceFieldArrayD> deviceFieldArrays(
    const std::vector<std::type_index>& fields);

};

}

#endif
