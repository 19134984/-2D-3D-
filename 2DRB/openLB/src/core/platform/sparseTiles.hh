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

#ifndef SPARSE_TILES_HH
#define SPARSE_TILES_HH

#include "sparseTiles.h"

#include <algorithm>

namespace olb {

namespace sparse {

inline TileData makeTileData(TileShape shape,
                             std::size_t p0, std::size_t p1,
                             unsigned sx, unsigned sy, unsigned sz,
                             unsigned tile, unsigned pad)
{
  TileData g;
  g.capacity = shape.tx * shape.ty * shape.tz;
  if (shape.tx == 1 && shape.ty == 1 && shape.tz == tile) {
    return g; // default 1x1xTILE keeps the linear scheme (bit identical)
  }
  auto log2u = [](unsigned v) { unsigned l=0; while ((1u<<l) < v) ++l; return l; };
  // Align the tile grid to the core origin so surface tiles pack tightly
  auto align = [&](unsigned t) { return ((pad + t - 1) / t) * t - pad; };
  g.linear = 0;
  g.p0 = p0;
  g.p1 = p1;
  g.lx = log2u(shape.tx);
  g.ly = log2u(shape.ty);
  g.lz = log2u(shape.tz);
  g.maskx = shape.tx - 1;
  g.masky = shape.ty - 1;
  g.maskz = shape.tz - 1;
  g.tileExtentY = shape.ty;
  g.tileExtentZ = shape.tz;
  g.addx = align(shape.tx);
  g.addy = align(shape.ty);
  g.addz = align(shape.tz);
  g.gx = ((g.addx + sx - 1) >> g.lx) + 1;
  g.gy = ((g.addy + sy - 1) >> g.ly) + 1;
  g.gz = ((g.addz + sz - 1) >> g.lz) + 1;
  return g;
}

template<typename T, unsigned TILE>
T* HostSparseTiles<T,TILE>::tileFor(std::size_t i, bool create)
{
  const std::size_t key = _geo.tileId(i, TILE);
  T* result = nullptr;
  #ifdef PARALLEL_MODE_OMP
  #pragma omp critical(olb_sparse_tileFor)
  #endif
  {
    auto it = _tiles.find(key);
    if (it == _tiles.end() && create) {
      auto tile = std::unique_ptr<T[]>(new T[TILE]);
      std::fill(tile.get(), tile.get() + TILE, _default);
      it = _tiles.emplace(key, std::move(tile)).first;
    }
    if (it != _tiles.end()) {
      result = it->second.get();
    }
  }
  return result;
}

template<typename T, unsigned TILE>
void HostSparseTiles<T,TILE>::copyTilesFrom(const HostSparseTiles& rhs)
{
  for (const auto& [key,tile] : rhs._tiles) {
    auto copy = std::unique_ptr<T[]>(new T[TILE]);
    std::copy(tile.get(), tile.get() + TILE, copy.get());
    _tiles.emplace(key, std::move(copy));
  }
}

template<typename T, unsigned TILE>
void HostSparseTiles<T,TILE>::fillDenseFromTiles(T* dense) const
{
  if (_geo.linear) {
    for (const auto& [key,tile] : _tiles) {
      const std::size_t base = key*TILE;
      for (std::size_t j=0; j < TILE && base+j < _count; ++j) {
        dense[base+j] = tile[j];
      }
    }
  } else {
    for (std::size_t i=0; i < _count; ++i) {
      const auto it = _tiles.find(_geo.tileId(i, TILE));
      if (it != _tiles.end()) {
        dense[i] = it->second[_geo.offset(i, TILE)];
      }
    }
  }
}

template<typename T, unsigned TILE>
void HostSparseTiles<T,TILE>::setTileShape(TileShape shape, std::size_t p0, std::size_t p1,
                                           unsigned sx, unsigned sy, unsigned sz, unsigned pad)
{
  if (_finalized) {
    throw std::runtime_error("sparse field tile shape must be set before the first "
                             "ProcessingContext::Simulation. Runtime reconfiguration is unsupported");
  }
  _tiles.clear();
  _geo = makeTileData(shape, p0, p1, sx, sy, sz, TILE, pad);
}

template<typename T, unsigned TILE>
void HostSparseTiles<T,TILE>::ensureAllocated(std::size_t i)
{
  if (_mode == FieldStorage::Sparse) {
    tileFor(i, true);
  }
}

template<typename T, unsigned TILE>
const T& HostSparseTiles<T,TILE>::getSparse(std::size_t i) const
{
  const auto it = _tiles.find(_geo.tileId(i, TILE));
  return it == _tiles.end() ? _default : it->second[_geo.offset(i, TILE)];
}

template<typename T, unsigned TILE>
T& HostSparseTiles<T,TILE>::sparseRef(std::size_t i)
{
  return tileFor(i, true)[_geo.offset(i, TILE)];
}

template<typename T, unsigned TILE>
void HostSparseTiles<T,TILE>::finalize(ProcessingContext context)
{
  if (context == ProcessingContext::Simulation) {
    _finalized = true;
  }
}

template<typename T, unsigned TILE>
void HostSparseTiles<T,TILE>::packTiles(std::vector<std::int32_t>& hostDir,
                                        std::vector<T>& hostPayload) const
{
  const std::size_t cap = _geo.capacity ? _geo.capacity : TILE;
  const std::size_t tileGrid = _geo.linear ? (_count + TILE - 1) / TILE
                                           : _geo.gridTileCount();
  hostDir.assign(tileGrid, -1);
  hostPayload.assign((_tiles.size() + 1) * cap, _default);
  std::int32_t slot = 1;
  for (const auto& [key,tile] : _tiles) {
    hostDir[key] = slot;
    std::copy(tile.get(), tile.get() + cap, hostPayload.begin() + std::size_t(slot)*cap);
    ++slot;
  }
}

template<typename T, unsigned TILE>
void HostSparseTiles<T,TILE>::unpackTiles(const std::int32_t* hostDir, const T* hostPayload)
{
  const std::size_t cap = _geo.capacity ? _geo.capacity : TILE;
  for (const auto& [key,tile] : _tiles) {
    const std::int32_t slot = hostDir[key];
    if (slot >= 1) {
      std::copy(hostPayload + std::size_t(slot)*cap,
                hostPayload + std::size_t(slot+1)*cap,
                tile.get());
    }
  }
}

}

}

#endif
