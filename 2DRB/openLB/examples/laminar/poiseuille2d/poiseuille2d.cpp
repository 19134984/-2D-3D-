/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2007, 2012 Jonas Latt, Mathias J. Krause
 *  Vojtech Cvrcek, Peter Weisbrod
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

/* poiseuille2d.cpp:
 * This example examines a 2D Poiseuille flow
 * It illustrates the usage of different flow setups, boundary conditions
 * and computation of error norms.
 * The following simulation options can be combined freely:
 * - if the compiler flag ENABLE_MRT is set, mrt collision operators are used
 * - forced/ nonForced flow
 * - different boundary conditions
 * - simulation only or eoc convergence analysis
 */

// the main code of the simulation is in poiseuille2d.h as it is also used by the
// example poiseuille2dEoc

// set flag in order to use mrt collision operators instead of bgk
//#define ENABLE_MRT

#include "case.h"

int main( int argc, char* argv[] )
{
  OstreamManager clout(std::cout,"simulatePoiseuille");

  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION>(50);
    myCaseParameters.set<DOMAIN_EXTENT>({2., 1.});
    myCaseParameters.set<PHYS_CHAR_VELOCITY>(1.0);
    myCaseParameters.set<REYNOLDS>(10);
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.8);
    myCaseParameters.set<PHYS_MAX_T>(30.);
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1.0);
    myCaseParameters.set<TUNER_PARAM>(0.5);

    myCaseParameters.set<FLOW_TYPE>(FlowType::nonForced);
    myCaseParameters.set<BOUNDARY_TYPE>(BoundaryType::interpolated);

    myCaseParameters.set<CONVERGENCE_PRECISION>(1e-9);
    myCaseParameters.set<CONVERGED>(false);
    myCaseParameters.set<CONV_ITER>(0.25);
    myCaseParameters.set<VTK_ENABLED>(true);
    myCaseParameters.set<GNUPLOT_ENABLED>(true);
    myCaseParameters.set<COMPUTE_ERROR>(true);
    myCaseParameters.set<EOC>(false);

    myCaseParameters.set<PHYS_VTK_ITER_T>([] (MyCase::ParametersD& p) {
      return p.get<PHYS_MAX_T>()/20.;
    });

    myCaseParameters.set<PHYS_STAT_ITER_T>([] (MyCase::ParametersD& p) {
      return p.get<PHYS_MAX_T>()/20.;
    });

  }

  if (myCaseParameters.get<parameters::EOC>()) {
    // Convergence study: one simulation per resolution, the error plot
    // tmp/gnuplotData/error.png collects the errors of all of them
    myCaseParameters.set<parameters::VTK_ENABLED>(false);
    for (int resolution : {30, 50, 70, 100}) {
      myCaseParameters.set<parameters::RESOLUTION>(resolution);
      myCaseParameters.set<parameters::CONVERGED>(false);
      run(myCaseParameters);
    }
  }
  else {
    run(myCaseParameters);
  }
}
