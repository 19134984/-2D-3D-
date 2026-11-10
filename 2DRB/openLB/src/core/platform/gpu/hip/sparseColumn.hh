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

#ifndef GPU_HIP_SPARSE_COLUMN_HH
#define GPU_HIP_SPARSE_COLUMN_HH

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

namespace olb {

namespace gpu {

namespace hip {

template<typename T, unsigned TILE>
struct SparseColumn<T,TILE>::DeviceData {
  thrust::device_vector<std::int32_t> dir;
  thrust::device_vector<T>            payload;
  /// Flat per-cell device buffer, used only in Dense mode
  thrust::device_vector<T>            dense;
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
    if (context == ProcessingContext::Simulation) {
      _device->dense = _dense;
    } else {
      thrust::host_vector<T> hostDense = _device->dense;
      std::copy(hostDense.begin(), hostDense.end(), _dense.begin());
    }
    return;
  }
  if (context == ProcessingContext::Simulation) {
    std::vector<std::int32_t> hostDir;
    std::vector<T> hostPayload;
    this->packTiles(hostDir, hostPayload);
    _device->dir = hostDir;
    _device->payload = hostPayload;
  }
  else {
    const thrust::host_vector<std::int32_t> hostDir = _device->dir;
    const thrust::host_vector<T> hostPayload = _device->payload;
    this->unpackTiles(thrust::raw_pointer_cast(hostDir.data()),
                      thrust::raw_pointer_cast(hostPayload.data()));
  }
}

template<typename T, unsigned TILE>
typename SparseColumn<T,TILE>::Accessor SparseColumn<T,TILE>::getDeviceAccessor()
{
  if (this->_mode == FieldStorage::Dense) {
    // dir == nullptr signals a flat per-cell payload to the device accessor
    return Accessor{
      nullptr,
      thrust::raw_pointer_cast(_device->dense.data()),
      this->_default,
      this->_geo
    };
  }
  return Accessor{
    thrust::raw_pointer_cast(_device->dir.data()),
    thrust::raw_pointer_cast(_device->payload.data()),
    this->_default,
    this->_geo
  };
}

}

}

}

#endif
