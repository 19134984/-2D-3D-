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

#ifndef QUALITY_F_H
#define QUALITY_F_H

namespace olb {

namespace functors {

/// Computes the local Knudsen metric
/// (cf. https://arxiv.org/pdf/1507.06767)
struct KnudsenMetricF {
  using parameters_t = meta::list<>;

  using result_t = descriptors::KNUDSEN_METRIC;
  using fields_t = meta::list<>;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;

    V rho {};
    V u[DESCRIPTOR::d] {};
    cell.computeRhoU(rho, u);
    V uSqr {};
    for (unsigned D=0; D<DESCRIPTOR::d; ++D) {
      uSqr += u[D]*u[D];
    }
    V sum = {};
    for (unsigned iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
      const V fEq = equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr);
      sum += util::abs((cell[iPop] - fEq) / (fEq + descriptors::t<V,DESCRIPTOR>(iPop)));
    }

    FieldD<V,DESCRIPTOR,typename KnudsenMetricF::result_t> metric = sum / DESCRIPTOR::q;
    return metric;
  }
};

/// Computes the local Knudsen metric for each population
/// (cf. https://arxiv.org/pdf/1507.06767)
struct KnudsenMetricPerPopulationF {
  using parameters_t = meta::list<>;

  using result_t = descriptors::KNUDSEN_METRIC_PER_POPULATION;
  using fields_t = meta::list<>;

  template <typename CELL, typename PARAMETERS>
  auto compute(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;

    V rho {};
    V u[DESCRIPTOR::d] {};
    cell.computeRhoU(rho, u);
    V uSqr {};
    for (unsigned D=0; D<DESCRIPTOR::d; ++D) {
      uSqr += u[D]*u[D];
    }
    V sum = {};
    FieldD<V,DESCRIPTOR,KnudsenMetricPerPopulationF::result_t> metricPerPopulation {};
    for (unsigned iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
      const V fEq = equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr);
      metricPerPopulation[iPop] = util::abs((cell[iPop] - fEq) / (fEq + descriptors::t<V,DESCRIPTOR>(iPop)));
    }
    return metricPerPopulation / DESCRIPTOR::q;
  }
};

}

}

#endif
