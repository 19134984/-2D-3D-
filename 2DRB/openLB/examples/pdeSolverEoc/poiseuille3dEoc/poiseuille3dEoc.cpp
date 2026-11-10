/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2018 Marc Haußmann, Mathias J. Krause, Jonathan Jeppener-Haltenhoff
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

/* poiseuille3d.cpp:
 * This example examines a 3D Poseuille flow
 * It illustrates the usage of different flow setups, boundary conditions
 * and computation of error norms.
 * The following simulation options can be combined freely:
 * - if the compiler flag ENABLE_MRT is set, mrt collision operators are used
 * - forced/ nonForced flow
 * - different boundary conditions
 * - simulation only or eoc convergence analysis
 */

// the main code of the simulation is in poiseuille3d.h as it is also used by the
// example ../../laminar/poiseuille3d

// set flag in order to use mrt collision operators instead of bgk
//#define ENABLE_MRT

#include "../../laminar/poiseuille3d/case.h"

void simulatePoiseuilleForEOC(MyCase::ParametersD& parameters, Gnuplot<MyCase::value_t>& gplot) {
  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(parameters);

  /// === Step 4: Create Case ===
  MyCase myCase(parameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase);

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);

  const FlowType      flowType        = parameters.get<parameters::FLOW_TYPE>();
  const BoundaryType  boundaryType    = parameters.get<parameters::BOUNDARY_TYPE>();
  const bool          noslipBoundary  = ((boundaryType != BoundaryType::FREE_SLIP) && (boundaryType != BoundaryType::PARTIAL_SLIP));
  const std::size_t   res             = parameters.get<parameters::RESOLUTION>();

  if (noslipBoundary) {
    if (flowType == FlowType::NON_FORCED) {
      gplot.setData(
        MyCase::value_t(res),
        { parameters.get<parameters::VELOCITY_L1_ABS_ERROR>(),
          parameters.get<parameters::VELOCITY_L2_ABS_ERROR>(),
          parameters.get<parameters::VELOCITY_LINF_ABS_ERROR>(),
          parameters.get<parameters::STRAIN_RATE_L1_ABS_ERROR>(),
          parameters.get<parameters::STRAIN_RATE_L2_ABS_ERROR>(),
          parameters.get<parameters::STRAIN_RATE_LINF_ABS_ERROR>(),
          parameters.get<parameters::WSS_L1_ABS_ERROR>(),
          parameters.get<parameters::WSS_L2_ABS_ERROR>(),
          parameters.get<parameters::WSS_LINF_ABS_ERROR>(),
          parameters.get<parameters::PRESSURE_L1_ABS_ERROR>(),
          parameters.get<parameters::PRESSURE_L2_ABS_ERROR>(),
          parameters.get<parameters::PRESSURE_LINF_ABS_ERROR>() },
        { "velocity L1 abs Error","velocity L2 abs Error",
          "velocity Linf abs error",
          "strain rate L1 abs error", "strain rate L2 abs error",
          "strain rate Linf abs error",
          "wall shear stress L1 abs error", "wall shear stress L2 abs error",
          "wall shear stress Linf abs error",
          "pressure L1 abs error", "pressure L2 abs error",
          "pressure Linf abs error" },
        "top right",
        { 'p','p','p','p','p','p','p','p','p','p','p','p' } );
    } else {
      // same as above, but without pressure computation
      gplot.setData (
        MyCase::value_t(res),
        { parameters.get<parameters::VELOCITY_L1_ABS_ERROR>(),
          parameters.get<parameters::VELOCITY_L2_ABS_ERROR>(),
          parameters.get<parameters::VELOCITY_LINF_ABS_ERROR>(),
          parameters.get<parameters::STRAIN_RATE_L1_ABS_ERROR>(),
          parameters.get<parameters::STRAIN_RATE_L2_ABS_ERROR>(),
          parameters.get<parameters::STRAIN_RATE_LINF_ABS_ERROR>(),
          parameters.get<parameters::WSS_L1_ABS_ERROR>(),
          parameters.get<parameters::WSS_L2_ABS_ERROR>(),
          parameters.get<parameters::WSS_LINF_ABS_ERROR>() },
        { "velocity L1 abs Error","velocity L2 abs Error",
          "velocity Linf abs error",
          "strain rate L1 abs error", "strain rate L2 abs error",
          "strain rate Linf abs error",
          "wall shear stress L1 abs error", "wall shear stress L2 abs error",
          "wall shear stress Linf abs error",},
        "top right",
        { 'p','p','p','p','p','p','p','p', 'p' } );
    }
  }
}

int main( int argc, char* argv[] )
{
  OstreamManager clout( std::cout,"main" );

  // === 1st Step: Initialization ===
  initialize( &argc, &argv );

  MyCase::ParametersD myCaseParameters;
  setDefaultParameters(myCaseParameters);
  myCaseParameters.set<parameters::RESOLUTION_START>(21);
  myCaseParameters.set<parameters::RESOLUTION_MAX>(52);
  myCaseParameters.set<parameters::RESOLUTION_STEP>(10);
  myCaseParameters.set<parameters::EOC>(true);

  BoundaryType boundaryType = myCaseParameters.get<parameters::BOUNDARY_TYPE>();
  bool forbiddenEOCCombination = (boundaryType == BoundaryType::FREE_SLIP) || (boundaryType == BoundaryType::PARTIAL_SLIP);
  if (forbiddenEOCCombination) {
    std::runtime_error("eoc computation is currently not supported for slip boundary conditions");
  }

  std::string bcName = "bc" + std::to_string(int(myCaseParameters.get<parameters::BOUNDARY_TYPE>()));
  std::string runName = (myCaseParameters.get<parameters::FLOW_TYPE>() == FlowType::FORCED) ? bcName + "_force" : bcName + "_nonForce";
  singleton::directories().setOutputDir( "./tmp/" + runName + "/" );

  // Initialize gnuplot
  Gnuplot<MyCase::value_t> gplot(
    "Velocity_and_StrainRate_eoc",
    false,
    "set terminal png size 720, 720 font 'Arial,10'",
    Gnuplot<MyCase::value_t>::LOGLOG,
    Gnuplot<MyCase::value_t>::LINREG);

  // set the labels for the plot
  gplot.setLabel("Resolution test", "average Error");

  std::size_t startN = myCaseParameters.get<parameters::RESOLUTION_START>();
  std::size_t maxN   = myCaseParameters.get<parameters::RESOLUTION_MAX>();
  std::size_t stepN  = myCaseParameters.get<parameters::RESOLUTION_STEP>();

  for(size_t simuN = startN; simuN < maxN; simuN += stepN){
    clout << "Starting next simulation with N = " << simuN << std::endl;
    myCaseParameters.set<parameters::RESOLUTION>(simuN);
    simulatePoiseuilleForEOC(myCaseParameters, gplot);
  }

  gplot.writePNG();
  return 0;
}
