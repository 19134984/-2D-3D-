/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2026 Shota Ito
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

#ifndef SMOOTH_INFLOW_O_H
#define SMOOTH_INFLOW_O_H

#include <core/operator.h>
#include "../singleLatticeO.h"

namespace olb {

namespace operators {

struct SmoothConstInflowO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct INFLOW_VELOCITY : public descriptors::FIELD_BASE<0,1> {};

  using parameters = meta::list<INFLOW_VELOCITY>;

  int getPriority() const {
   return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    auto u = parameters.template get<INFLOW_VELOCITY>();
    cell.template setField<momenta::FixedVelocityMomentumGeneric::VELOCITY>(u);
  }
};

}

}
#endif
