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

#ifndef OLB_CORE_HALF_H
#define OLB_CORE_HALF_H

#include <cstdint>

#include "core/platform/platform.h"

#ifdef __CUDACC__
#include <cuda_fp16.h>
#endif

namespace olb {

namespace detail {

/// IEEE binary16 conversion, platform-portable (software). Stage 2 replaces the
/// device path with hardware cvt (__half); host stays software here.
any_platform inline std::uint16_t floatToHalfBits(float f) {
  const std::uint32_t x    = __builtin_bit_cast(std::uint32_t, f);
  const std::uint32_t sign = (x >> 16) & 0x8000u;
  const std::int32_t  e    = static_cast<std::int32_t>((x >> 23) & 0xffu) - 112; // 127-15
  const std::uint32_t mant = x & 0x7fffffu;
  if (e >= 0x1f) {                                    // overflow / inf / nan
    if (((x >> 23) & 0xffu) == 0xffu) {
      return static_cast<std::uint16_t>(sign | 0x7c00u | (mant ? 0x200u : 0u));
    }
    return static_cast<std::uint16_t>(sign | 0x7c00u);
  } else if (e <= 0) {                                // subnormal / underflow
    if (e < -10) {
      return static_cast<std::uint16_t>(sign);
    }
    std::uint32_t m = (mant | 0x800000u) >> (1 - e);
    if (m & 0x1000u) { m += 0x2000u; }               // round to nearest
    return static_cast<std::uint16_t>(sign | (m >> 13));
  } else {                                           // normal
    std::uint32_t h = sign | (static_cast<std::uint32_t>(e) << 10) | (mant >> 13);
    if (mant & 0x1000u) { h += 1; }                  // round to nearest
    return static_cast<std::uint16_t>(h);
  }
}

any_platform inline float halfBitsToFloat(std::uint16_t h) {
  const std::uint32_t sign = static_cast<std::uint32_t>(h & 0x8000u) << 16;
  const std::uint32_t exp  = (h >> 10) & 0x1fu;
  const std::uint32_t mant = h & 0x3ffu;
  std::uint32_t out;
  if (exp == 0) {
    if (mant == 0) {
      out = sign;
    } else {                                         // subnormal
      std::int32_t  ee = -1;
      std::uint32_t mm = mant;
      do { ee++; mm <<= 1; } while ((mm & 0x400u) == 0);
      out = sign | (static_cast<std::uint32_t>(127 - 15 - ee) << 23) | ((mm & 0x3ffu) << 13);
    }
  } else if (exp == 0x1f) {                          // inf / nan
    out = sign | 0x7f800000u | (mant << 13);
  } else {                                           // normal
    out = sign | (static_cast<std::uint32_t>(exp + 112) << 23) | (mant << 13);
  }
  return __builtin_bit_cast(float, out);
}

}

/// IEEE binary16 storage scalar: a 2-byte POD carrying ONLY conversions to/from
/// float. Deliberately no arithmetic operators, so half arithmetic cannot happen
/// implicitly -- values promote to the compute type before any math. Selected as
/// a per-field storage type via FieldPrecision::F16; all lattice arithmetic stays
/// in the compute type through the cell-boundary conversion proxy.
struct half {
  std::uint16_t bits;

  half() = default;

  any_platform half(float f) {
#ifdef __CUDA_ARCH__
    bits = __half_as_ushort(__float2half_rn(f));   // hardware cvt on device
#else
    bits = detail::floatToHalfBits(f);             // portable software cvt on host
#endif
  }

  any_platform operator float() const {
#ifdef __CUDA_ARCH__
    return __half2float(__ushort_as_half(bits));   // hardware cvt on device
#else
    return detail::halfBitsToFloat(bits);
#endif
  }
};

}

#endif
