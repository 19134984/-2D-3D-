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

#ifndef GPU_SYCL_SPARSE_COLUMN_HH
#define GPU_SYCL_SPARSE_COLUMN_HH

#include "device.h"

namespace olb {

namespace gpu {

namespace sycl {

template<typename T, unsigned TILE>
struct SparseColumn<T,TILE>::DeviceData {
  gpu::sycl::device::unique_ptr<std::int32_t> dir;
  gpu::sycl::device::unique_ptr<T>            payload;
  /// Flat per-cell device buffer, used only in Dense mode
  gpu::sycl::device::unique_ptr<T>            dense;
  std::size_t dirCapacity     = 0;
  std::size_t payloadCapacity = 0;
  std::size_t denseCapacity   = 0;
  /// Packed sizes of the last Simulation sync (needed to read back on Evaluation)
  std::size_t dirCount        = 0;
  std::size_t payloadCount    = 0;
};

template<typename T, unsigned TILE>
SparseColumn<T,TILE>::SparseColumn(std::size_t count):
  sparse::HostSparseTiles<T,TILE>(count),
  _device(new DeviceData())
{ }

template<typename T, unsigned TILE>
SparseColumn<T,TILE>::SparseColumn(SparseColumn<T,TILE>&& rhs):
  sparse::HostSparseTiles<T,TILE>(std::move(rhs)),
  _dense(std::move(rhs._dense)),
  _device(rhs._device.release())
{ }

template<typename T, unsigned TILE>
SparseColumn<T,TILE>::SparseColumn(const SparseColumn<T,TILE>& rhs):
  sparse::HostSparseTiles<T,TILE>(rhs),
  _dense(rhs._dense),
  _device(new DeviceData())
{ }

template<typename T, unsigned TILE>
SparseColumn<T,TILE>::~SparseColumn()
{ }

template<typename T, unsigned TILE>
void SparseColumn<T,TILE>::setStorageMode(FieldStorage mode)
{
  if (this->_finalized) {
    throw std::runtime_error("sparse field storage mode must be set before the first "
                             "ProcessingContext::Simulation; runtime re-configuration is unsupported");
  }
  if (mode == this->_mode) {
    return;
  }
  if (mode == FieldStorage::Dense) {
    _dense.assign(this->_count, this->_default);
    this->fillDenseFromTiles(_dense.data());
    this->_tiles.clear();
  } else {
    _dense.clear();
    _dense.shrink_to_fit();
  }
  this->_mode = mode;
}

template<typename T, unsigned TILE>
void SparseColumn<T,TILE>::resize(std::size_t newCount)
{
  if (this->_mode == FieldStorage::Dense) {
    _dense.resize(newCount, this->_default);
  }
  this->_count = newCount;
}

template<typename T, unsigned TILE>
const T& SparseColumn<T,TILE>::operator[](std::size_t i) const
{
  if (this->_mode == FieldStorage::Dense) {
    return _dense[i];
  }
  return this->getSparse(i);
}

template<typename T, unsigned TILE>
T& SparseColumn<T,TILE>::operator[](std::size_t i)
{
  if (this->_mode == FieldStorage::Dense) {
    return _dense[i];
  }
  return this->sparseRef(i);
}

template<typename T, unsigned TILE>
void SparseColumn<T,TILE>::setProcessingContext(ProcessingContext context)
{
  this->finalize(context);
  if (this->_mode == FieldStorage::Dense) {
    // Flat host <-> device buffer, identical to a dense Column sync
    if (this->_count == 0) {
      return;
    }
    if (_device->denseCapacity < this->_count) {
      _device->dense = gpu::sycl::device::malloc<T>(this->_count);
      _device->denseCapacity = this->_count;
    }
    if (context == ProcessingContext::Simulation) {
      gpu::sycl::device::copyToDevice(_dense.data(), _device->dense.get(), this->_count*sizeof(T));
    } else {
      gpu::sycl::device::copyToHost(_device->dense.get(), _dense.data(), this->_count*sizeof(T));
    }
    return;
  }
  if (context == ProcessingContext::Simulation) {
    std::vector<std::int32_t> hostDir;
    std::vector<T> hostPayload;
    this->packTiles(hostDir, hostPayload);
    if (_device->dirCapacity < hostDir.size()) {
      _device->dir = gpu::sycl::device::malloc<std::int32_t>(hostDir.size());
      _device->dirCapacity = hostDir.size();
    }
    if (_device->payloadCapacity < hostPayload.size()) {
      _device->payload = gpu::sycl::device::malloc<T>(hostPayload.size());
      _device->payloadCapacity = hostPayload.size();
    }
    _device->dirCount     = hostDir.size();
    _device->payloadCount = hostPayload.size();
    if (!hostDir.empty()) {
      gpu::sycl::device::copyToDevice(hostDir.data(), _device->dir.get(), hostDir.size()*sizeof(std::int32_t));
    }
    if (!hostPayload.empty()) {
      gpu::sycl::device::copyToDevice(hostPayload.data(), _device->payload.get(), hostPayload.size()*sizeof(T));
    }
  }
  else {
    std::vector<std::int32_t> hostDir(_device->dirCount);
    std::vector<T> hostPayload(_device->payloadCount);
    if (_device->dirCount > 0) {
      gpu::sycl::device::copyToHost(_device->dir.get(), hostDir.data(), hostDir.size()*sizeof(std::int32_t));
    }
    if (_device->payloadCount > 0) {
      gpu::sycl::device::copyToHost(_device->payload.get(), hostPayload.data(), hostPayload.size()*sizeof(T));
    }
    this->unpackTiles(hostDir.data(), hostPayload.data());
  }
}

template<typename T, unsigned TILE>
typename SparseColumn<T,TILE>::Accessor SparseColumn<T,TILE>::getDeviceAccessor()
{
  if (this->_mode == FieldStorage::Dense) {
    // dir == nullptr signals a flat per-cell payload to the device accessor
    return Accessor{
      nullptr,
      _device->dense.get(),
      this->_default,
      this->_geo
    };
  }
  return Accessor{
    _device->dir.get(),
    _device->payload.get(),
    this->_default,
    this->_geo
  };
}

}

}

}

#endif
