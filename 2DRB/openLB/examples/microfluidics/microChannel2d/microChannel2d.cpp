/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Fedor Bukreev, Adrian Kummerländer,
 *  Mathias J. Krause
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
// This is an 2D example for a laminar flow at higher Knudsen numbers
// with Knudsen slip condition and thermal creep along the wall

#include "microChannel2d.h"

int main(int argc, char* argv[]) {
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<KNUDSEN>(0.45);
    myCaseParameters.set<B_COEFF>(-1);
    myCaseParameters.set<ACCOMODATION_COEFF>(1);
    myCaseParameters.set<RESOLUTION>(21);
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(1.13);
    myCaseParameters.set<PHYS_MAX_T>(2.e-6);
    myCaseParameters.set<PHYS_START_T>([] (MyCase::ParametersD& p) {
      return 0.05*p.get<parameters::PHYS_MAX_T>();
    });
    myCaseParameters.set<KINETIC_DIAMETER>(364e-12); //kinetic diamter of nitrogen [m]
    myCaseParameters.set<AVERAGE_PRESSURE>(1000);
    myCaseParameters.set<PRESSURE_DROP>(0.5);
    myCaseParameters.set<MOLAR_MASS>(28.96e-3);
    myCaseParameters.set<R_SPECIFIC>(296.8);
    myCaseParameters.set<parameters::TEMPERATURE>(300.);
    myCaseParameters.set<parameters::TEMPERATURE_LEFT>(298.);
    myCaseParameters.set<parameters::TEMPERATURE_RIGHT>(302.);
    myCaseParameters.set<PHYS_CHAR_DENSITY>([] (MyCase::ParametersD& p) {
      return util::idealGasDensity( p.get<MOLAR_MASS>(), p.get<AVERAGE_PRESSURE>(), p.get<parameters::TEMPERATURE>());
    });
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>([] (MyCase::ParametersD& p) {
      return util::gasDynamicViscosity( p.get<MOLAR_MASS>(), p.get<AVERAGE_PRESSURE>(), p.get<parameters::TEMPERATURE>(), p.get<KINETIC_DIAMETER>())/p.get<PHYS_CHAR_DENSITY>();
    });
    myCaseParameters.set<MFP>([] (MyCase::ParametersD& p) {
      return util::meanFreePath(p.get<MOLAR_MASS>(), p.get<AVERAGE_PRESSURE>(), p.get<parameters::TEMPERATURE>(), p.get<PHYS_CHAR_VISCOSITY>()*p.get<PHYS_CHAR_DENSITY>());
    });
    myCaseParameters.set<DOMAIN_EXTENT>([] (MyCase::ParametersD& p) {
      return Vector{10.*p.get<MFP>()/p.get<KNUDSEN>(), p.get<MFP>()/p.get<KNUDSEN>()};
    });
    myCaseParameters.set<SLIPCOEFF>([] (MyCase::ParametersD& p) {
      return (2. - p.get<ACCOMODATION_COEFF>())/p.get<ACCOMODATION_COEFF>()/(1. - p.get<B_COEFF>()*p.get<KNUDSEN>())*p.get<MFP>();
    });
    myCaseParameters.set<CREEPCOEFF>([] (MyCase::ParametersD& p) {
      return 3./4.*p.get<PHYS_CHAR_VISCOSITY>()*p.get<PHYS_CHAR_DENSITY>()*p.get<R_SPECIFIC>()/p.get<AVERAGE_PRESSURE>();
    });
    myCaseParameters.set<PHYS_CHAR_VELOCITY>([] (MyCase::ParametersD& p) {
      return p.get<DOMAIN_EXTENT>()[1]*p.get<DOMAIN_EXTENT>()[1]/2./p.get<PHYS_CHAR_VISCOSITY>()/p.get<PHYS_CHAR_DENSITY>()
      * (-p.get<PRESSURE_DROP>())/p.get<DOMAIN_EXTENT>()[0]
      * (0.25 - 0.5 - p.get<SLIPCOEFF>()/p.get<DOMAIN_EXTENT>()[1])
      + p.get<CREEPCOEFF>()*(p.get<parameters::TEMPERATURE_RIGHT>()
      - p.get<parameters::TEMPERATURE_LEFT>())/p.get<DOMAIN_EXTENT>()[0];
    });
    myCaseParameters.set<PHYS_INTERVAL_T>(2.e-6/50.);
    myCaseParameters.set<RESIDUUM>(1e-6);
    myCaseParameters.set<CONVERGED>(false);
    myCaseParameters.set<ERROR_NORMS>({0,0,0,0,0,0});
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
