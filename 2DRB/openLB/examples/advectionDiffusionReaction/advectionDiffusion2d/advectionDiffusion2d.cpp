
/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2020 Louis Kronberg, Stephan Simonis
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

/*  adePeriodic2d:
 *  The solution to a linear, scalar, two-dimensional advection-diffusion
 *  equation is approximated.
 *  The numerical setup and the analytical solution are taken from
 *  Simonis, S., Frank, M., and Krause M. J. Applied Mathematics Letters
 *  (2023) 135:108484, DOI: https://doi.org/10.1016/j.aml.2022.108484.
 *  Error norms are calculated for three subsequent resolutions and stored
 *  in the respective /tmp folders. A python script is provided to calculate
 *  the experimental order of convergence towards the analytical solution.
 */

#include "case.h"

int main(int argc, char *argv[])
{
  OstreamManager clout(std::cout,"main");
  initialize(&argc, &argv);

  singleton::directories().setOutputDir("./tmp/");

  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION>(50);
    myCaseParameters.set<LATTICE_STAT_ITER_T>(50);
    myCaseParameters.set<PHYS_CHAR_LENGTH>(2.);
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1.);
    myCaseParameters.set<PHYS_DIFFUSIVITY>(0.05);
    myCaseParameters.set<PECLET>(100.); // Peclet number (Pe = u*L/mue)
    myCaseParameters.set<PULSE_DIFF_BOUND>(1e-1);
    myCaseParameters.set<RESOLUTION_MAX>(myCaseParameters.get<RESOLUTION>());
  }

  Mesh mesh = createMesh(myCaseParameters);

  MyCase myCase(myCaseParameters, mesh);

  prepareGeometry(myCase);

  prepareLattice(myCase);

  setInitialValues(myCase);

  simulate(myCase);
}
