/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Florian Kaiser
 *                2006, 2007, 2012 Jonas Latt, Mathias J. Krause
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

/* bstep3d.cpp:
 * The implementation of a backward facing step. It is furthermore
 * shown how to use checkpointing to save the state of the
 * simulation regularly.
 */

#include <olb.h>
using namespace olb;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes,
  Lattice<float,
          descriptors::D3Q19<>
  >
>;

namespace olb::parameters {

struct PHYS_STEP_EXTENT : public descriptors::FIELD_BASE<0,1> {};
struct PHYS_STEP_ORIGIN : public descriptors::FIELD_BASE<0,1> {};

} // namespace olb::parameters

/// @brief Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters)
{
  using T = MyCase::value_t;

  // setup channel
  const Vector extendChannel = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector originChannel = parameters.get<parameters::DOMAIN_ORIGIN>();
  IndicatorCuboid3D<T> channel(extendChannel, originChannel);

  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();
  Mesh<T, MyCase::d> mesh(channel, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

/// @brief Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
/// @note The material numbers will be used to assign physics to lattice nodes
void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;
  using T          = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry  = myCase.getGeometry();

  geometry.rename( 0, 2 );

  geometry.rename( 2, 1, {1, 1, 1} );

  const Vector channelExtent  = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector stepExtent     = parameters.get<parameters::PHYS_STEP_EXTENT>();
  const Vector stepOrigin     = parameters.get<parameters::PHYS_STEP_ORIGIN>();
  const T physDeltaX          = parameters.get<parameters::PHYS_DELTA_X>();

  IndicatorCuboid3D<T> step( stepExtent, stepOrigin );

  geometry.rename( 1, 2, step );

  // Set material number for inflow
  const Vector inflowExtent{physDeltaX, channelExtent[1], channelExtent[2]};
  const Vector inflowOrigin{-physDeltaX/2., 0., 0.};
  IndicatorCuboid3D<T> inflow( inflowExtent, inflowOrigin );
  geometry.rename( 2, 3, 1, inflow );

  // Set material number for outflow
  const Vector outflowExtent{2 * physDeltaX, channelExtent[1], channelExtent[2]};
  const Vector outflowOrigin{channelExtent[0] - 3. / 2. * physDeltaX, 0., 0.};
  IndicatorCuboid3D<T> outflow( outflowExtent, outflowOrigin );
  geometry.rename( 2, 4, 1, outflow );

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  // Removes all not needed boundary voxels inside the surface
  geometry.innerClean();
  geometry.checkForErrors();
  geometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

/// @brief Set lattice dynamics
/// @param myCase The Case instance which keeps the simulation data
void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare lattice ..." << std::endl;

  using T = MyCase::value_t;

  auto& geometry   = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes{});
  using DESCRIPTOR = MyCase::descriptor_t;

  lattice.setUnitConverter<UnitConverter<T,DESCRIPTOR>>(
    parameters.get<parameters::PHYS_DELTA_X>(),
    parameters.get<parameters::PHYS_DELTA_T>(),
    parameters.get<parameters::PHYS_CHAR_LENGTH>(),
    parameters.get<parameters::PHYS_CHAR_VELOCITY>(),
    parameters.get<parameters::PHYS_CHAR_VISCOSITY>(),
    parameters.get<parameters::PHYS_CHAR_DENSITY>()
  );
  const auto& converter = lattice.getUnitConverter();

  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("bstep3d");

  auto bulkIndicator = geometry.getMaterialIndicator({1, 3, 4});

  // Material=1 -->bulk dynamics
  // Material=3 -->bulk dynamics (inflow)
  // Material=4 -->bulk dynamics (outflow)
  dynamics::set<BGKdynamics>(lattice, bulkIndicator);
  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(lattice, geometry, 2);

  //if boundary conditions are chosen to be local
  boundary::set<boundary::LocalVelocity>(lattice, geometry, 3);
  boundary::set<boundary::LocalPressure>(lattice, geometry, 4);

  //if boundary conditions are chosen to be interpolated
  // boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 3);
  // boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 4);

  clout << "Prepare lattice ... OK" << std::endl;
}

/// Set initial condition for primal variables (velocity and density)
/// @param myCase The Case instance which keeps the simulation data
/// @note Be careful: initial values have to be set using lattice units
void setInitialValues(MyCase& myCase)
{
  auto& lattice  = myCase.getLattice(NavierStokes{});
  lattice.setParameter<descriptors::OMEGA>(
    lattice.getUnitConverter().getLatticeRelaxationFrequency());
  lattice.initialize();
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{
  OstreamManager clout(std::cout, "setTemporalValues");
  using T = MyCase::value_t;

  auto&             lattice    = myCase.getLattice(NavierStokes {});
  auto&             converter  = lattice.getUnitConverter();
  auto&             geometry   = myCase.getGeometry();
  auto&             parameters = myCase.getParameters();
  const std::size_t iTmaxStart = converter.getLatticeTime(parameters.get<parameters::PHYS_START_T>());
  const std::size_t iTUpdate   = converter.getLatticeTime(parameters.get<parameters::PHYS_BOUNDARY_VALUE_UPDATE_T>());

  if ( iT % iTUpdate == 0 && iT <= iTmaxStart ) {
    T frac{};
    PolynomialStartScale<T,std::size_t>(iTmaxStart, T(1))(&frac, &iT);
    std::vector<T> maxVelocity(3, 0);
    maxVelocity[0] = 2.25 * frac * converter.getCharPhysVelocity();

    T distance2Wall = converter.getPhysDeltaX() / 2.;
    RectanglePoiseuille3D<T> poiseuilleU( geometry, 3, maxVelocity, distance2Wall, distance2Wall, distance2Wall );
    momenta::setVelocity(lattice, geometry.getMaterialIndicator(3), poiseuilleU);

    if (iT % (10 * iTUpdate) == 0 && iT <= iTmaxStart) {
      clout << "step=" << iT << "; maxVel=" << maxVelocity[0] << std::endl;
    }

    lattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

void getResults(MyCase& myCase,
                std::size_t iT,
                util::Timer<MyCase::value_t> &timer)
{
  OstreamManager clout(std::cout, "getResults");
  using   T             = MyCase::value_t;
  using   DESCRIPTOR    = MyCase::descriptor_t;
  auto&   parameters    = myCase.getParameters();
  auto&   lattice       = myCase.getLattice(NavierStokes{});
  auto&   converter     = lattice.getUnitConverter();
  auto&   geometry      = myCase.getGeometry();

  SuperVTMwriter3D<T> vtmWriter( "bstep3d" );

  const std::size_t vtkIter  = converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T>());
  const std::size_t statIter = converter.getLatticeTime(parameters.get<parameters::PHYS_STAT_ITER_T>());

  if ( iT == 0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( lattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( lattice );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Writes the ppm files
  if ( iT % vtkIter == 0 ) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( lattice, converter );
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( lattice, converter );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );

    vtmWriter.write( iT );
    SuperEuklidNorm3D<T> normVel(velocity);
    BlockReduction3D2D<T> planeReduction(normVel,
                                         Hyperplane3D<T>().centeredIn(geometry.getCuboidDecomposition().getMotherCuboid()).normalTo({0,0,1}),
                                         600,
                                         BlockDataSyncMode::ReduceOnly);
    // write output as JPEG
    heatmap::plotParam<T> jpeg_Param;
    jpeg_Param.maxValue = converter.getCharPhysVelocity() * 3./2.;
    jpeg_Param.minValue = 0.0;
    jpeg_Param.fullScreenPlot = true;
    heatmap::write(planeReduction, iT, jpeg_Param);
  }

  // Writes output on the console
  if ( iT % statIter == 0 && iT >= 0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    lattice.getStatistics().print( iT, converter.getPhysTime( iT ) );
  }
}

void simulate(MyCase& myCase)
{
  OstreamManager clout(std::cout, "simulate");

  using T           = MyCase::value_t;
  auto& parameters  = myCase.getParameters();
  auto& lattice     = myCase.getLattice(NavierStokes{});

  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(
    parameters.get<parameters::PHYS_MAX_T>()
  );

  util::Timer<T> timer(iTmax,
                       myCase.getGeometry().getStatistics().getNvoxel());

  clout << "Starting simulation..." << std::endl;
  timer.start();

  for (std::size_t iT = 0; iT < iTmax; ++iT) {

    setTemporalValues(myCase, iT);

    lattice.collideAndStream();

    getResults(myCase, iT, timer);
  }

  clout << "Simulation finished." << std::endl;
  lattice.setProcessingContext(ProcessingContext::Evaluation);
  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);
  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION>(20);
    myCaseParameters.set<TIME_RESOLUTION>(25);
    myCaseParameters.set<PHYS_CHAR_LENGTH>(1.);
    myCaseParameters.set<PHYS_DELTA_X>([] (MyCase::ParametersD& p) {return
      p.get<PHYS_CHAR_LENGTH>() / p.get<RESOLUTION>();
    });
    myCaseParameters.set<PHYS_DELTA_T>([] (MyCase::ParametersD& p) {return
      p.get<PHYS_DELTA_X>() / p.get<TIME_RESOLUTION>();
    });
    myCaseParameters.set<PHYS_MAX_T>(40.0);
    myCaseParameters.set<PHYS_CHAR_VELOCITY>(1.0);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(1.0 / 100.);
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1.0);
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.518);

    // Startup
    myCaseParameters.set<PHYS_START_T>([] (MyCase::ParametersD& p) {return
      0.2 * p.get<PHYS_MAX_T>();
    });
    myCaseParameters.set<PHYS_BOUNDARY_VALUE_UPDATE_T>(0.01);

    // Domain
    myCaseParameters.set<DOMAIN_EXTENT>({18.0, 1.5, 1.5});
    myCaseParameters.set<DOMAIN_ORIGIN>({0., 0., 0.});
    myCaseParameters.set<PHYS_STEP_EXTENT>({5.0, 0.75, 1.5});
    myCaseParameters.set<PHYS_STEP_ORIGIN>({0., 0., 0.});

    // Output
    myCaseParameters.set<PHYS_STAT_ITER_T>(0.2);
    myCaseParameters.set<PHYS_VTK_ITER_T>(0.2);
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

  return 0;
}
