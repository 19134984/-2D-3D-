/*  Lattice Boltzmann sample, written in C++, using the OpenLB library
 *
 *  Copyright (C) 2006-2019 Jonas Latt, Mathias J. Krause,
 *  Vojtech Cvrcek, Peter Weisbrod, Adrian Kummerländer
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

/* cylinder2d.h:
 * This example examines a steady flow past a cylinder placed in a channel.
 * The cylinder is offset somewhat from the center of the flow to make the
 * steady-state symmetrical flow unstable. At the inlet, a Poiseuille profile is
 * imposed on the velocity, whereas the outlet implements a Dirichlet pressure
 * condition set by p = 0.
 * Inspired by "Benchmark Computations of Laminar Flow Around
 * a Cylinder" by M.Schäfer and S.Turek. For high resolution, low
 * latticeU, and enough time to converge, the results for pressure drop, drag
 * and lift lie within the estimated intervals for the exact results.
 * An unsteady flow with Karman vortex street can be created by changing the
 * Reynolds number to Re=100.
 */

#include <olb.h>
#include "case.h"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::names;

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION             >(10);
    myCaseParameters.set<REYNOLDS               >(20.);
    myCaseParameters.set<PHYS_MAX_T             >(16.);
    myCaseParameters.set<DOMAIN_EXTENT          >({2.2, .41});
    myCaseParameters.set<PHYS_CHAR_VELOCITY     >(0.2);
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.56);
    myCaseParameters.set<PHYS_CHAR_DENSITY      >(1.0);
    myCaseParameters.set<RADIUS_CYLINDER        >(0.05);
    myCaseParameters.set<PHYS_CHAR_LENGTH       >([](MyCase::ParametersD& p){
      return 2 * p.get<RADIUS_CYLINDER>();
    });
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>([](MyCase::ParametersD& p){
      return p.get<PHYS_CHAR_LENGTH>() *
             p.get<parameters::PHYS_CHAR_VELOCITY>() /
             p.get<parameters::REYNOLDS>();
    });
    myCaseParameters.set<CENTER_CYLINDER>({0.2, 0.2});
    myCaseParameters.set<PHYS_VTK_ITER_T>(0.3);
    myCaseParameters.set<PHYS_STAT_ITER_T>(0.1);
    myCaseParameters.set<RAMP_UP_T>(0.01);
    myCaseParameters.set<RAMP_UP_END_FRACTION>(0.4);
  }

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
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
