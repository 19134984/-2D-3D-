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

#ifndef GPU_SYCL_SPARSE_COLUMN_H
#define GPU_SYCL_SPARSE_COLUMN_H

#include <memory>
#include <cstdint>
#include <vector>

#include "core/platform/gpu/sparseAccessor.h"

namespace olb {

namespace gpu {

namespace sycl {

/// Tile-backed column for SYCL GPU targets. setProcessingContext packs tiles into a device directory + payload
/**
 * Thin wrapper over HostSparseTiles: adds only the host Dense buffer, the USM device mirror and the
 * device accessor. All tile logic lives in the shared host base.
 **/
template<typename T, unsigned TILE = 64>
class SparseColumn final : public AbstractColumn<T>
                         , public Serializable
                         , public sparse::HostSparseTiles<T,TILE> {
private:
  /// Host flat storage, used only in Dense mode
  std::vector<T> _dense;

  struct DeviceData;
  std::unique_ptr<DeviceData> _device;

public:
  using value_t = T;
  using Accessor = gpu::DeviceSparseColumn<T,TILE>;

  SparseColumn(std::size_t count);
  SparseColumn(SparseColumn<T,TILE>&& rhs);
  SparseColumn(const SparseColumn<T,TILE>& rhs);
  ~SparseColumn();

  void setStorageMode(FieldStorage mode);

  const T& operator[](std::size_t i) const override;
  T& operator[](std::size_t i) override;

  std::size_t size() const override { return this->_count; }
  void resize(std::size_t newCount);

  /// Device view over the currently packed directory + payload
  Accessor getDeviceAccessor();

  void setProcessingContext(ProcessingContext);

  std::size_t getNblock() const override { return 0; }
  std::size_t getSerializableSize() const override { return 0; }
  bool* getBlock(std::size_t, std::size_t&, bool) override { return nullptr; }

};

}

}

/// Declare gpu::sycl::SparseColumn as the sparse column for GPU SYCL targets
template <typename T, unsigned TILE>
struct SparseImplementationOf<T,Platform::GPU_SYCL,TILE> {
  using type = gpu::sycl::SparseColumn<T,TILE>;
};

}

#endif
