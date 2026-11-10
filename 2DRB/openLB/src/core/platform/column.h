/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef COLUMN_H
#define COLUMN_H

#include "platform.h"

#include <cstdint>
#include <cstddef>

namespace olb {

/// Storage class of a field: dense per-cell columns or sparse tile-backed.
/// A sparse-capable field's column can switch between these at runtime.
enum class FieldStorage {
  Dense,
  Sparse,
};

/// Storage precision of a field, orthogonal to FieldStorage. Inherit keeps the
/// lattice's compute type T; F16/F32/F64 store the field in that precision while
/// all cell arithmetic still happens in T (conversion at the cell boundary).
enum class FieldPrecision {
  Inherit,
  F16,
  F32,
  F64,
};

/// Abstract declarator of Column-like storage
template <typename T>
struct AbstractColumn {
  using value_type = T;

  virtual std::size_t size() const = 0;

  virtual const T& operator[](std::size_t i) const = 0;
  virtual       T& operator[](std::size_t i)       = 0;
};

/// Abstract declarator of cyclic Column-like storage
template <typename T>
struct AbstractCyclicColumn {
  using value_type = T;

  virtual std::size_t size() const = 0;

  virtual const T& operator[](std::size_t i) const = 0;
  virtual       T& operator[](std::size_t i)       = 0;
};

/// Specializable declarator for concrete implementations of abstract storage types
template <typename ABSTRACT, Platform PLATFORM>
struct ImplementationOf;

}

#include "sparseTiles.h"

#endif

#include "cpu/sisd/column.h"

#ifdef PLATFORM_CPU_SIMD
#include "cpu/simd/column.h"
#endif

#ifdef PLATFORM_GPU_CUDA
#include "gpu/cuda/column.h"
#endif

#ifdef PLATFORM_GPU_HIP
#include "gpu/hip/column.h"
#endif

#ifdef PLATFORM_GPU_SYCL
#include "gpu/sycl/column.h"
#endif
