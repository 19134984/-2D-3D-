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

#ifndef CPU_SISD_SPARSE_COLUMN_H
#define CPU_SISD_SPARSE_COLUMN_H

#include <memory>
#include <stdexcept>
#include <unordered_map>

namespace olb {

namespace cpu {

namespace sisd {

/// Tile-backed column materialising storage only for touched cells (drop-in AbstractColumn)
/**
 * Thin wrapper over HostSparseTiles: adds only the SISD dense buffer and the AbstractColumn /
 * Serializable interface. All tile logic lives in the shared host base.
 **/
template<typename T, unsigned TILE = 64>
class SparseColumn final : public AbstractColumn<T>
                         , public Serializable
                         , public sparse::HostSparseTiles<T,TILE> {
private:
  /// Flat contiguous storage, allocated only in Dense mode (identical layout to Column)
  std::unique_ptr<T[]> _dense;

public:
  using value_t = T;

  SparseColumn(std::size_t count):
    sparse::HostSparseTiles<T,TILE>(count)
  { }

  SparseColumn():
    SparseColumn(0)
  { }

  SparseColumn(SparseColumn<T,TILE>&& rhs):
    sparse::HostSparseTiles<T,TILE>(std::move(rhs)),
    _dense(rhs._dense.release())
  { }

  SparseColumn(const SparseColumn<T,TILE>& rhs):
    sparse::HostSparseTiles<T,TILE>(rhs)
  {
    if (rhs._dense) {
      _dense = std::unique_ptr<T[]>(new T[this->_count]);
      std::copy(rhs._dense.get(), rhs._dense.get() + this->_count, _dense.get());
    }
  }

  virtual ~SparseColumn() = default;

  /// Select storage mode. Dense materialises a flat per-cell buffer equivalent to a
  /// dense Column, Sparse keeps tile-backed storage. Set at setup, before writes.
  void setStorageMode(FieldStorage mode)
  {
    if (this->_finalized) {
      throw std::runtime_error("sparse field storage mode must be set before the first "
                               "ProcessingContext::Simulation; runtime re-configuration is unsupported");
    }
    if (mode == this->_mode) {
      return;
    }
    if (mode == FieldStorage::Dense) {
      _dense = std::unique_ptr<T[]>(new T[this->_count]);
      std::fill(_dense.get(), _dense.get() + this->_count, this->_default);
      this->fillDenseFromTiles(_dense.get());
      this->_tiles.clear();
    } else {
      _dense.reset();
    }
    this->_mode = mode;
  }

  void resize(std::size_t count)
  {
    if (this->_mode == FieldStorage::Dense) {
      auto data = std::unique_ptr<T[]>(new T[count]);
      std::fill(data.get(), data.get() + count, this->_default);
      std::copy(_dense.get(), _dense.get() + std::min(this->_count, count), data.get());
      _dense.swap(data);
    }
    this->_count = count;
  }

  const T& operator[](std::size_t i) const override
  {
    if (this->_mode == FieldStorage::Dense) {
      return _dense[i];
    }
    return this->getSparse(i);
  }

  T& operator[](std::size_t i) override
  {
    if (this->_mode == FieldStorage::Dense) {
      return _dense[i];
    }
    return this->sparseRef(i);
  }

  std::size_t size() const override
  {
    return this->_count;
  }

  void setProcessingContext(ProcessingContext context) {
    this->finalize(context);
  };

  // Serialisation of the variable tile set is a later phase, opt out for now
  std::size_t getNblock() const override { return 0; }
  std::size_t getSerializableSize() const override { return 0; }
  bool* getBlock(std::size_t, std::size_t&, bool) override { return nullptr; }

};

}

}

/// Declare cpu::sisd::SparseColumn as the sparse column for CPU SISD targets
template <typename T, unsigned TILE>
struct SparseImplementationOf<T,Platform::CPU_SISD,TILE> {
  using type = cpu::sisd::SparseColumn<T,TILE>;
};

}

#endif
