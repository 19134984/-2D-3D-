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

#ifndef COMMUNICATION_CANARY_H
#define COMMUNICATION_CANARY_H

#include <cstdint>
#include <cstdlib>
#include <limits>
#include <string_view>
#include <type_traits>

#include "core/blockStructure.h"

namespace olb {

template <typename T>
class ConstSpan;

namespace overlap {

/// Returns true iff untransmitted population slots are to be poisoned
/**
 * The directional exchange leaves padding populations outside R(p)
 * undefined by contract. In builds with FEATURES=DIRECTIONAL_CANARY and
 * OLB_DIRECTIONAL_CANARY=1 at runtime, every masked exchange overwrites
 * the untransmitted slots of its cells with quiet NaN so that any
 * consumer implicitly relying on the full dense halo fails loudly
 * instead of silently reading stale data. Stripped entirely from
 * regular builds.
 **/
#ifdef FEATURE_DIRECTIONAL_CANARY
inline bool canaryEnabled()
{
  static const bool enabled = []() {
    const char* canary = std::getenv("OLB_DIRECTIONAL_CANARY");
    return canary && std::string_view(canary) == "1";
  }();
  return enabled;
}
#else
constexpr bool canaryEnabled() { return false; }
#endif

}

namespace util {

/// Poison the mask-complement slots of the given cells (directional canary)
/**
 * The masks must be index-aligned with cells; fully dropped cells (mask 0)
 * poison all components. Each mask bit is a streaming direction carrying
 * componentsPerDirection components, so component iD maps to direction bit
 * iD / componentsPerDirection (1:1 for a scalar lattice).
 **/
template <typename FIELD_ARRAY>
void writeDirectionalCanary(FIELD_ARRAY& fieldArray,
                            unsigned components,
                            ConstSpan<CellID> cells,
                            ConstSpan<std::uint64_t> masks,
                            unsigned componentsPerDirection = 1)
{
  using value_t = std::remove_reference_t<decltype(fieldArray[0][0])>;
  if constexpr (std::is_floating_point_v<value_t>) {
    const value_t poison = std::numeric_limits<value_t>::quiet_NaN();
    for (std::size_t iCell=0; iCell < masks.size(); ++iCell) {
      for (unsigned iD=0; iD < components; ++iD) {
        if (!(masks[iCell] & (std::uint64_t{1} << (iD / componentsPerDirection)))) {
          fieldArray[iD][cells[iCell]] = poison;
        }
      }
    }
  }
}

}

}

#endif
