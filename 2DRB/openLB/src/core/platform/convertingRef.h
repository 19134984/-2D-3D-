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

#ifndef OLB_CORE_PLATFORM_CONVERTING_REF_H
#define OLB_CORE_PLATFORM_CONVERTING_REF_H

#include <type_traits>

#include "core/platform/platform.h"

namespace olb {

/// Conversion proxy bridging a stored element (precision S) and the compute type V.
/// Reads promote to V before any arithmetic; writes convert back to S. Deliberately
/// exposes only value access and assignment -- no operator that would let half math
/// happen implicitly. Used at the cell population seam when a field's storage
/// precision differs from the lattice compute type; when they match the seams return
/// a plain V& instead (see makeConvertingRef), keeping untagged fields zero-overhead.
template <typename V, typename S>
class ConvertingRef {
private:
  S* _p;

public:
  any_platform ConvertingRef(S* p): _p(p) { }

  any_platform operator V() const { return static_cast<V>(*_p); }

  any_platform ConvertingRef& operator=(V v)  { *_p = S(v);                    return *this; }
  any_platform ConvertingRef& operator=(const ConvertingRef& rhs) { *_p = S(static_cast<V>(rhs)); return *this; }
  any_platform ConvertingRef& operator+=(V v) { *_p = S(static_cast<V>(*_p) + v); return *this; }
  any_platform ConvertingRef& operator-=(V v) { *_p = S(static_cast<V>(*_p) - v); return *this; }
  any_platform ConvertingRef& operator*=(V v) { *_p = S(static_cast<V>(*_p) * v); return *this; }
  any_platform ConvertingRef& operator/=(V v) { *_p = S(static_cast<V>(*_p) / v); return *this; }
};

/// Wrap a stored element as either a plain reference (when its type already is the
/// compute type V -- the untagged/inherit default, bit-identical) or a converting
/// ConvertingRef (when storage precision differs).
template <typename V, typename S>
any_platform decltype(auto) makeConvertingRef(S& s) {
  if constexpr (std::is_same_v<V,S>) {
    return (S&)s;
  } else {
    return ConvertingRef<V,S>(&s);
  }
}

}

#endif
