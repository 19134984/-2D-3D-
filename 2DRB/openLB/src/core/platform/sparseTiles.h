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

#ifndef SPARSE_TILES_H
#define SPARSE_TILES_H

#include "platform.h"

#include <cstdint>
#include <cstddef>
#include <memory>
#include <vector>
#include <stdexcept>
#include <unordered_map>
#include <set>

namespace olb {

namespace sparse {

/// Runtime tile shape. The extent product must equal the column TILE capacity and each extent must
/// be a power of two. The default 1x1x64 groups 64 consecutive cells along the fastest (z) axis.
struct TileShape {
  unsigned tx;
  unsigned ty;
  unsigned tz;
};

/// Addressing data mapping a flat cell index to its tile id and intra tile offset.
/**
 * POD memcpy'd verbatim into the device accessor, so the field order and sizes are part of the ABI
 * (see the byte offsets in the per backend registry). The linear fast path handles the default
 * 1x1xTILE shape with plain division. The shaped path decomposes the index into (x,y,z) via the
 * padded projection strides and packs it into a power of two tile grid.
 * Internal device-layout POD: the field order/sizes are ABI, not a public API.
 **/
struct alignas(8) TileData {
  /// Cells per x slab of the padded block (stride to recover x from a flat index)
  std::uint64_t p0{0};
  /// Cells per y row of the padded block (stride to recover y within an x slab)
  std::uint64_t p1{0};
  /// Tile grid extent along x
  std::uint32_t gx{1};
  /// Tile grid extent along y
  std::uint32_t gy{1};
  /// Tile grid extent along z
  std::uint32_t gz{1};
  /// log2 of the tile extent along x (shift a cell coord to its tile coord)
  std::uint32_t lx{0};
  /// log2 of the tile extent along y
  std::uint32_t ly{0};
  /// log2 of the tile extent along z
  std::uint32_t lz{0};
  /// Bitmask (tile extent minus one) selecting the intra tile x coord
  std::uint32_t maskx{0};
  /// Bitmask for the intra tile y coord
  std::uint32_t masky{0};
  /// Bitmask for the intra tile z coord
  std::uint32_t maskz{0};
  /// Tile extent along y (tx is recoverable as maskx plus one)
  std::uint32_t tileExtentY{1};
  /// Tile extent along z
  std::uint32_t tileExtentZ{1};
  /// 1 selects the linear 1x1xTILE fast path, 0 the shaped path
  std::uint32_t linear{1};
  /// Additive shift aligning the core origin to a tile corner along x (see makeTileData)
  std::uint32_t addx{0};
  /// Alignment shift along y
  std::uint32_t addy{0};
  /// Alignment shift along z
  std::uint32_t addz{0};
  /// Runtime cells per tile (tx*ty*tz). 0 means use the compile time TILE
  std::uint32_t capacity{0};

  any_platform std::size_t tileId(std::size_t i, unsigned tile) const {
    if (linear) {
      return i / tile;
    }
    const std::size_t x = i / p0 + addx;
    const std::size_t r = i % p0;
    const std::size_t y = r / p1 + addy;
    const std::size_t z = r % p1 + addz;
    return ((x >> lx)*gy + (y >> ly))*gz + (z >> lz);
  }

  any_platform unsigned offset(std::size_t i, unsigned tile) const {
    if (linear) {
      return i % tile;
    }
    const std::size_t x = i / p0 + addx;
    const std::size_t r = i % p0;
    const std::size_t y = r / p1 + addy;
    const std::size_t z = r % p1 + addz;
    return (unsigned(x & maskx)*tileExtentY + unsigned(y & masky))*tileExtentZ + unsigned(z & maskz);
  }

  std::size_t gridTileCount() const {
    return std::size_t(gx)*gy*gz;
  }
};

/// Build the tile addressing data for a shape given the block's padded projection and padded size
TileData makeTileData(TileShape shape,
                      std::size_t p0, std::size_t p1,
                      unsigned sx, unsigned sy, unsigned sz,
                      unsigned tile, unsigned pad = 0);

/// Platform agnostic host tile storage shared by every backend's sparse column.
/**
 * Owns the materialised tile set, the addressing data, the storage mode and the finalization guard.
 * Backends derive from this and add only their per platform dense buffer and device mirror. The flat
 * buffer conversion and the device pack/unpack live here so every backend reuses identical logic.
 * TILE is the compile time per field tile capacity (power of two).
 **/
template<typename T, unsigned TILE = 64>
class HostSparseTiles {
  static_assert(TILE > 0 && (TILE & (TILE-1)) == 0 && TILE <= 256,
                "sparse tile capacity must be a power of two and at most 256");
protected:
  std::size_t _count;
  T           _default;
  FieldStorage _mode{FieldStorage::Sparse};
  std::unordered_map<std::size_t,std::unique_ptr<T[]>> _tiles;
  TileData _geo;
  /// Set on the first Simulation sync. Storage mode and shape are then locked
  bool _finalized{false};

  T* tileFor(std::size_t i, bool create);
  void copyTilesFrom(const HostSparseTiles& rhs);
  void fillDenseFromTiles(T* dense) const;

public:
  HostSparseTiles(std::size_t count):
    _count(count),
    _default{}
  { }

  HostSparseTiles(HostSparseTiles&& rhs):
    _count(rhs._count),
    _default(rhs._default),
    _mode(rhs._mode),
    _tiles(std::move(rhs._tiles)),
    _geo(rhs._geo)
  { }

  HostSparseTiles(const HostSparseTiles& rhs):
    _count(rhs._count),
    _default(rhs._default),
    _mode(rhs._mode),
    _geo(rhs._geo)
  {
    copyTilesFrom(rhs);
  }

  /// Value returned for cells in unallocated tiles. Set before any write
  void setDefaultValue(T value) { _default = value; }

  FieldStorage getStorageMode() const { return _mode; }
  std::size_t size() const { return _count; }

  std::size_t getAllocatedTiles() const { return _tiles.size(); }
  unsigned getTileCapacity() const { return _geo.capacity ? _geo.capacity : TILE; }
  TileShape getTileShape() const {
    return _geo.linear ? TileShape{1,1,TILE} : TileShape{_geo.maskx+1, _geo.tileExtentY, _geo.tileExtentZ};
  }
  /// Directory length (one slot per grid tile), the int32 tile to payload map size
  std::size_t getDirectorySlots() const {
    return _geo.linear ? (_count + TILE - 1) / TILE : _geo.gridTileCount();
  }

  /// Number of distinct tiles the given cells map into (for halo materialisation accounting)
  template <typename CELLS>
  std::size_t countTilesFor(const CELLS& cells) const {
    std::set<std::size_t> tiles;
    for (const auto c : cells) {
      tiles.insert(_geo.tileId(std::size_t(c), TILE));
    }
    return tiles.size();
  }

  /// Set the runtime tile shape from the block's padded projection and size, clearing any tiles.
  /// Call at setup before writes. The default 1x1x64 keeps the plain linear scheme.
  void setTileShape(TileShape shape, std::size_t p0, std::size_t p1,
                    unsigned sx, unsigned sy, unsigned sz, unsigned pad = 0);

  /// Materialise the tile owning cell i (default filled) without changing any value
  void ensureAllocated(std::size_t i);

  /// Value of cell i in Sparse mode (default if its tile is absent)
  const T& getSparse(std::size_t i) const;

  /// Reference to cell i in Sparse mode, materialising its tile
  T& sparseRef(std::size_t i);

  /// Lock the storage mode and shape on the first Simulation sync
  void finalize(ProcessingContext context);

  /// Pack the materialised tiles into a device directory and payload. Slot 0 is the discard scratch
  /// tile, materialised tiles take slots 1..N. The payload packs at the runtime capacity
  void packTiles(std::vector<std::int32_t>& hostDir, std::vector<T>& hostPayload) const;

  /// Read the materialised tiles back from a packed device directory and payload
  void unpackTiles(const std::int32_t* hostDir, const T* hostPayload);
};

}

/// Sparse column for a platform, defaulting to the dense column for platforms without one.
/// TILE is the compile time per field tile capacity (power of two)
template <typename T, Platform PLATFORM, unsigned TILE = 64>
struct SparseImplementationOf {
  using type = typename ImplementationOf<AbstractColumn<T>,PLATFORM>::type;
};

/// Whether PLATFORM provides a genuine sparse column implementation
constexpr bool platformHasSparseColumn(Platform platform) {
  return platform == Platform::CPU_SISD
      || platform == Platform::GPU_CUDA
      || platform == Platform::GPU_HIP
      || platform == Platform::GPU_SYCL;
}

}

#endif

#include "sparseTiles.hh"
