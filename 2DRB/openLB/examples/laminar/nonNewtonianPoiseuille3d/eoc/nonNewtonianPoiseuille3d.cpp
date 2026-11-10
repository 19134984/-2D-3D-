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
};
*/

#include <olb.h>

#include "../case.h"

using namespace olb;
using namespace olb::names;

template <typename MODEL, typename T>
void setupAndEOC( MyCase::ParametersD& myCaseParameters,
                  std::string viscosityModel,
                  std::vector<Vector<T,3>>& errors ) {
  std::string outDir = "./tmp/" + viscosityModel + "_eoc/";
  singleton::directories().setOutputDir(outDir);

  OstreamManager clout( std::cout,"eocPoiseuille" );
  Vector res = myCaseParameters.get<parameters::EOC_RESOLUTIONS>();

  for (int iter = 0; iter < (int) res.size(); ++iter) {
    /// === Step 2a: Increase resolution ===
    myCaseParameters.set<parameters::RESOLUTION>(res[iter]);

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
    errors.push_back(simulatePoiseuilleWith<MODEL,T>( myCase ));
  }
}

int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  initialize( &argc, &argv );
  OstreamManager clout( std::cout,"main" );

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  using T = MyCase::value_t;
  setGetParameters( myCaseParameters, argc, argv );

  std::vector<Vector<T,3>> errors;
  switch ( myCaseParameters.get<parameters::VISCOSITY_MODEL>() ) {
    case ViscosityModel::NEWTONIAN:
    default:
      setupAndEOC<Newtonian,T>( myCaseParameters, "Newtonian", errors );
      break;
    case ViscosityModel::POWER_LAW:
      setupAndEOC<PowerLaw,T>( myCaseParameters, "PowerLaw", errors );
      break;
    case ViscosityModel::CASSON:
      setupAndEOC<Casson,T>( myCaseParameters, "Casson", errors );
      break;
    case ViscosityModel::CARREAU_YASUDA:
      setupAndEOC<CarreauYasuda,T>( myCaseParameters, "CarreauYasuda", errors );
      break;
  }

  Vector res = myCaseParameters.get<parameters::EOC_RESOLUTIONS>();
  std::vector<T> eoc_L1, eoc_L2, eoc_Linf;
  for (int i = 0; i < (int) res.size() - 1; ++i) {
    eoc_L1.push_back( (util::log(errors.at(i)[0]) - util::log(errors.at(i+1)[0]))/(util::log(res[i]) - util::log(res[i+1])) );
    eoc_L2.push_back( (util::log(errors.at(i)[1]) - util::log(errors.at(i+1)[1]))/(util::log(res[i]) - util::log(res[i+1])) );
    eoc_Linf.push_back( (util::log(errors.at(i)[2]) - util::log(errors.at(i+1)[2]))/(util::log(res[i]) - util::log(res[i+1])) );
  }

  for (int i = 0; i < (int) res.size() -1; ++i) {
    clout << "Error L1: " << errors.at(i)[0] << std::endl;
    clout << "Error L2: " << errors.at(i)[1] << std::endl;
    clout << "Error Linf: " << errors.at(i)[2] << std::endl;
  }
  clout << "EOC with N: " << res << std::endl;
  clout << "EOC L1: " << eoc_L1 << std::endl;
  clout << "EOC L2: " << eoc_L2 << std::endl;
  clout << "EOC Linf: " << eoc_Linf << std::endl;
}
