/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Tim Bingert, Michael Rennick
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

/** flatInterface2d.cpp
 * In this example a flat interface test is performed for several multi-
 * phase models such as incompressible Allen-Cahn or well-balanced Cahn-
 * Hilliard. A rectangular domain of fluid 2 is immersed in fluid 1.
 * A diffusive interface forms with the profile of a hyperbolic tangent
 * whose accuracy is measured for multiple resolutions in order to test
 * the models experimental order of convergence. The equilibrium pressure
 * is also investigated in a similar way. The grid refinement can be
 * performed with either a constant or a decreasing Cahn number. Multiple
 * different error norms can then be written into an output EOC.dat file.
 *
 * This example shows the simplest application of the hybrid (and local)
 * phase field Allen-Cahn model with periodic boundaries, based on:
 *
 * Liu, Xi, Zhenhua Chai, and Baochang Shi. "Improved hybrid Allen-Cahn
 * phase-field-based lattice Boltzmann method for incompressible two-phase
 * flows." Physical Review E 107.3 (2023): 035308.
 *
 * It also shows the same application of the well-balanced Cahn-Hilliard
 * model that has no spurious currents (machine accuracy), from:
 *
 * Ju, Long, et al. "A well-balanced lattice Boltzmann model for binary
 * fluids based on the incompressible phase-field theory." arXiv preprint
 * arXiv:2311.10827 (2023).
 */

// Choose model from loading three different cases
#include "flatInterfaceConservativePhaseField2d.h"
// #include "flatInterfaceHybridAllenCahn2d.h"
// #include "flatInterfaceCahnHilliard2d.h"

const bool Cahn_const = true;  // scaling with const Cahn number?
const int Ny = 100;            // initial resolution
const double w = 5.;           // initial interface thickness
const double sigma = 0.01;     // initial surface tension

int main( int argc, char *argv[] )
{
  initialize( &argc, &argv );

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<DOMAIN_EXTENT          >({2e-6, 100e-6}); // domain size [physical units]
    myCaseParameters.set<RESOLUTION             >(Ny);             // resolution y [lattice units]
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.8);       // tau liquid [lattice units]
    myCaseParameters.set<RELAXATION_TIME_PF     >(0.8);            // tau mobility [lattice units]
    myCaseParameters.set<REYNOLDS               >(0.);             // Reynolds number []
    myCaseParameters.set<RHO_VAPOR              >(1.);             // physDensity gas/vapor [physical units]
    myCaseParameters.set<RHO_LIQUID             >(1000.);          // physDensity liquid [physical units]
    myCaseParameters.set<NU_VAPOR               >(9e-7);           // physViscosity gas/vapor [physical units]
    myCaseParameters.set<NU_LIQUID              >(9e-6);           // physViscosity liquid [physical units]
    myCaseParameters.set<PHYS_CHAR_PRESSURE     >(1e5);            // physPressure [physical units]
    myCaseParameters.set<SURFACE_TENSION        >(sigma);          // surface tension [lattice units]
    myCaseParameters.set<parameters::INTERFACE_WIDTH>(w);     // interface width [lattice units]
    myCaseParameters.set<C_RHO             >(1000.);          // conversion factor density [physical units]
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
  simulate(myCase, Cahn_const);
}
