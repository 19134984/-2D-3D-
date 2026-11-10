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

#ifndef TRACK_AVERAGE_FUNCTOR_O_H
#define TRACK_AVERAGE_FUNCTOR_O_H

#include <core/operator.h>
#include "../singleLatticeO.h"

namespace olb {

namespace operators {

template <concepts::PrimitiveFunctor FUNCTOR, typename AVERAGE_FIELD>
struct TrackAverageFunctorO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = typename FUNCTOR::parameters_t::template include<descriptors::LATTICE_TIME_T>;

  int getPriority() const {
   return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using DESCRIPTOR = typename CELL::descriptor_t;
    static_assert(FUNCTOR::result_t::template size<DESCRIPTOR>() == AVERAGE_FIELD::template size<DESCRIPTOR>(),
                  "TrackAverageFunctorO: Mismatch of dimension for FUNCTOR and AVERAGE_FIELD");

    auto average = cell.template getField<AVERAGE_FIELD>();
    auto iT = parameters.template get<descriptors::LATTICE_TIME_T>();
    auto value = FUNCTOR{}.compute(cell, parameters);
    cell.template setField<AVERAGE_FIELD>((average * (iT-1) + value) / iT);
  }
};

}

}
#endif
