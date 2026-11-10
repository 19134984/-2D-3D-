/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Shota Ito
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

/* 3D force-driven Poiseuille flow with different non-Newtonian viscosity models
 * Implemented options are Newtonian, PowerLaw, Casson, and Carreau-Yasuda models.
 * The first three models are validated for the same pressure drop and dynamic
 * viscosity used for the unit-conversion. The pressure drop is computed according
 * the analytical solution for the Newtonian case.
 * Carreau-Yasuda model is provided as a popular alternative model but still
 * remains untested yet.
 *
 * Pass parameter VISCOSITY_MODEL to choose (default=0)
 * 0 -> Newtonian
 * 1 -> Power Law
 * 2 -> Casson
 * 3 -> Carreau Yasuda
*/

#include <olb.h>
#include "case.h"

template <typename MODEL>
void setupAndSimulate( MyCase::ParametersD& myCaseParameters, std::string viscosityModel ) {
  std::string outDir = "./tmp/" + viscosityModel + "/";
  singleton::directories().setOutputDir(outDir);

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Lattice ===
  prepareLattice<MODEL>(myCase);

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues<MODEL>(myCase);

  /// === Step 8: Simulate ===
  simulatePoiseuilleWith<MODEL,MyCase::value_t>( myCase );

}

int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  initialize( &argc, &argv );

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  setGetParameters(myCaseParameters, argc, argv);

  std::string viscosityModel = "Newtonian";
  switch ( myCaseParameters.get<parameters::VISCOSITY_MODEL>() ) {
    case ViscosityModel::NEWTONIAN:
    default:
      setupAndSimulate<Newtonian>(myCaseParameters, viscosityModel);
      break;
    case ViscosityModel::POWER_LAW:
      viscosityModel = "PowerLaw";
      setupAndSimulate<PowerLaw>(myCaseParameters, viscosityModel);
      break;
    case ViscosityModel::CASSON:
      viscosityModel = "Casson";
      setupAndSimulate<Casson>(myCaseParameters, viscosityModel);
      break;
    case ViscosityModel::CARREAU_YASUDA:
      viscosityModel = "CarreauYasuda";
      setupAndSimulate<CarreauYasuda>( myCaseParameters, viscosityModel );
      break;
  }
}
