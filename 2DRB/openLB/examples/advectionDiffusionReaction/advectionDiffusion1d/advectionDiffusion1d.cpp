
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

/*  advectionDiffusion1d:
 *  The solution to a linear, scalar, one-dimensional advection-diffusion
 *  equation is approximated. The computation takes place in 2D and is
 *  projected to 1D by slicing the domain along a centerline.
 *  The numerical setup and the analytical solution are taken from
 *  [Simonis, S., Frank, M., and Krause, M. J. 2020. Phil. Trans. R. Soc. A378:
 *  20190400. DOI: 10.1098/rsta.2019.0400]
 *  Error norms are calculated for three subsequent resolutions and stored
 *  in the respective /tmp folders. A python script is provided to calculate
 *  the experimental order of convergence towards the analytical solution.
 */

#include "case.h"

using namespace olb;
using namespace olb::graphics;
using namespace olb::names;


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
    myCaseParameters.set<PHYS_DIFFUSIVITY>(1.5);
    myCaseParameters.set<PECLET>(40./3.);
    myCaseParameters.set<PULSE_DIFF_BOUND>(1e-3);
  }

  Mesh mesh = createMesh(myCaseParameters);

  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase);

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);
}
