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

#ifndef OLB_CORE_MOMENTS_DEMAND_H
#define OLB_CORE_MOMENTS_DEMAND_H

#include <vector>
#include <concepts>

#include "vector.h"

namespace olb {

enum class MomentSet {
  U       = 0b0010,
  RhoU    = 0b0011,
  RhoUPi  = 0b0111,
  RhoUJ   = 0b1011,
  RhoUPiJ = 0b1111,
};

constexpr MomentSet operator|(MomentSet a, MomentSet b) {
  return MomentSet(static_cast<int>(a) | static_cast<int>(b));
}

enum class Freshness {
  Fresh,
  AsStored,
};

template <unsigned D>
struct MomentsDemand {
  MomentSet                           moments;
  Freshness                           freshness = Freshness::Fresh;
  std::vector<Vector<std::int32_t,D>> stencil;
};

namespace concepts {

template <typename OPERATOR>
concept DeclaresMomentsDemand = requires {
  OPERATOR::getMomentsDemand();
};

}

template <unsigned D, typename OPERATOR>
std::vector<MomentsDemand<D>> getMomentsDemandOf() {
  if constexpr (concepts::DeclaresMomentsDemand<OPERATOR>) {
    return OPERATOR::getMomentsDemand();
  } else {
    return {};
  }
}

}

#endif
