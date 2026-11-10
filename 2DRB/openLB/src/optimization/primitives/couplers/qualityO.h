/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Shota Ito
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

#ifndef QUALITY_O_H
#define QUALITY_O_H

#include <solver/names.h>
#include <core/superLatticeCoupling.h>

namespace olb {

namespace couplers {

// Operator which computes the adjoint-based quality of interest
struct QualityOfInterestO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform {
    using V = typename CELLS::template value_t<names::Primal>::value_t;
    using DESCRIPTOR = typename CELLS::template value_t<names::Primal>::descriptor_t;
    auto& primal_cell = cells.template get<names::Primal>();
    auto& dual_cell = cells.template get<names::Dual>();

    auto knudsenMetricPerPopulation = functors::KnudsenMetricPerPopulationF{}.compute(primal_cell, parameters);

    FieldD<V,DESCRIPTOR,descriptors::QUALITY_OF_INTEREST> quality {};
    for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      quality += knudsenMetricPerPopulation[iPop]*dual_cell[iPop];
    }
    dual_cell.template setField<descriptors::QUALITY_OF_INTEREST>(quality);
  }
};

}

}

#endif
