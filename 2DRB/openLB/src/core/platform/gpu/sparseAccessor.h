/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Adrian Kummerlaender
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

#ifndef GPU_SPARSE_ACCESSOR_H
#define GPU_SPARSE_ACCESSOR_H

#include <cstdint>
#include <cstddef>
#include <cassert>

#include "core/platform/platform.h"
#include "core/platform/sparseTiles.h"

namespace olb {

namespace gpu {

/// Device-side view of a SparseColumn passed into kernels. Absent tiles read the default
/**
 * A tile is materialised (host side) only where the field is written. The directory maps a cell's
 * tile (i / TILE) to its payload slot, or -1 when the tile is absent. Absent cells read the default.
 * Per the device convention the kernel never allocates, so writes target a pre-allocated tile
 * (reserved host-side via defineField over the write region before launch). Slot 0 is a discard scratch tile: a stray
 * write to an unallocated tile lands there harmlessly instead of out of bounds.
 *
 * Backend-agnostic (CUDA/HIP/SYCL): the accessor holds only raw pointers and POD geometry.
 **/
template<typename T, unsigned TILE>
struct DeviceSparseColumn {
  static constexpr unsigned tile_size = TILE;

  const std::int32_t* dir;
  T*                  payload;
  T                   defaultValue;
  /// Addressing data placed after defaultValue so the fixed dir, payload and default offsets hold
  sparse::TileData  geo;

  any_platform T get(std::size_t i) const {
    if (dir == nullptr) { // dense mode: payload is a flat per-cell buffer
      return payload[i];
    }
    const unsigned cap = geo.capacity ? geo.capacity : TILE;
    const std::int32_t slot = dir[geo.tileId(i, TILE)];
    return slot < 0 ? defaultValue : payload[std::size_t(slot)*cap + geo.offset(i, TILE)];
  }

  any_platform T& operator[](std::size_t i) const {
    if (dir == nullptr) { // dense mode: payload is a flat per-cell buffer
      return payload[i];
    }
    const unsigned cap = geo.capacity ? geo.capacity : TILE;
    const std::int32_t slot = dir[geo.tileId(i, TILE)];
#ifdef FEATURE_SPARSE_CANARY
    // Trap device writes to an unallocated tile; reserve tiles with defineField first
    assert(slot >= 0 && "sparse write to an unallocated tile: defineField over the write region first");
#endif
    // Slot 0 is the discard scratch tile: an unallocated write lands there harmlessly.
    return payload[std::size_t(slot < 0 ? 0 : slot)*cap + geo.offset(i, TILE)];
  }
};

// The type-erased any-fields comm addresses a DeviceSparseColumn through these fixed offsets
using DeviceSparseColumnAccessorF = DeviceSparseColumn<float,64>;
using DeviceSparseColumnAccessorD = DeviceSparseColumn<double,64>;
using DeviceSparseColumnAccessorI = DeviceSparseColumn<int,64>;
static_assert(offsetof(DeviceSparseColumnAccessorF, dir) == 0
           && offsetof(DeviceSparseColumnAccessorF, payload) == 8
           && offsetof(DeviceSparseColumnAccessorF, defaultValue) == 16
           && offsetof(DeviceSparseColumnAccessorF, geo) == 24);
static_assert(offsetof(DeviceSparseColumnAccessorD, dir) == 0
           && offsetof(DeviceSparseColumnAccessorD, payload) == 8
           && offsetof(DeviceSparseColumnAccessorD, defaultValue) == 16
           && offsetof(DeviceSparseColumnAccessorD, geo) == 24);
static_assert(offsetof(DeviceSparseColumnAccessorI, dir) == 0
           && offsetof(DeviceSparseColumnAccessorI, payload) == 8
           && offsetof(DeviceSparseColumnAccessorI, defaultValue) == 16
           && offsetof(DeviceSparseColumnAccessorI, geo) == 24);

}

}

#endif
