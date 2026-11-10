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

#ifndef OLB_CORE_MOMENTS_PROVIDER_H
#define OLB_CORE_MOMENTS_PROVIDER_H

#include "momentsDemand.h"

namespace olb {

template <typename DYNAMICS, typename WRITE_SET, typename MOMENTS_STORE = void>
struct RefreshMomentsO {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return -10;
  }

  template <concepts::DynamicCell CELL, typename V = typename CELL::value_t>
  void apply(CELL& cell) any_platform {
    using DESCRIPTOR = typename CELL::descriptor_t;
    constexpr bool writeRho = WRITE_SET::template contains<fields::moments::RHO>();
    constexpr bool writeU   = WRITE_SET::template contains<fields::moments::U>();
    constexpr bool writePi  = WRITE_SET::template contains<fields::moments::PI>();
    constexpr bool writeJ   = WRITE_SET::template contains<fields::moments::J>();

    V rho{};
    V u[DESCRIPTOR::d] { };
    V pi[util::TensorVal<DESCRIPTOR>::n] { };
    V j[DESCRIPTOR::d] { };

    if constexpr (writePi) {
      typename DYNAMICS::MomentaF().computeAllMomenta(cell, rho, u, pi);
    } else if constexpr (writeU) {
      typename DYNAMICS::MomentaF().computeRhoU(cell, rho, u);
    } else if constexpr (writeRho) {
      rho = typename DYNAMICS::MomentaF().computeRho(cell);
    }
    if constexpr (writeJ) {
      typename DYNAMICS::MomentaF().computeJ(cell, j);
    }

    if constexpr (writeRho) {
      cell.template getFieldPointer<fields::moments::RHO>()[0] = rho;
    }
    if constexpr (writeU) {
      for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
        cell.template getFieldPointer<fields::moments::U>()[iD] = u[iD];
      }
    }
    if constexpr (writePi) {
      for (unsigned iPi=0; iPi < util::TensorVal<DESCRIPTOR>::n; ++iPi) {
        cell.template getFieldPointer<fields::moments::PI>()[iPi] = pi[iPi];
      }
    }
    if constexpr (writeJ) {
      for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
        cell.template getFieldPointer<fields::moments::J>()[iD] = j[iD];
      }
    }
  }
};

template <MomentSet SET>
using moment_write_set_t = std::conditional_t<
  SET == MomentSet::U,      meta::list<fields::moments::U>,
  std::conditional_t<
  SET == MomentSet::RhoU,   meta::list<fields::moments::RHO, fields::moments::U>,
  std::conditional_t<
  SET == MomentSet::RhoUPi, meta::list<fields::moments::RHO, fields::moments::U, fields::moments::PI>,
  std::conditional_t<
  SET == MomentSet::RhoUJ,  meta::list<fields::moments::RHO, fields::moments::U, fields::moments::J>,
                            meta::list<fields::moments::RHO, fields::moments::U, fields::moments::PI, fields::moments::J>>>>>;

}

#endif
