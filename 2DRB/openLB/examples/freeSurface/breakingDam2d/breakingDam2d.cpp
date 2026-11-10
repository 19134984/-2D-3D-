/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Claudius Holeksa, Maximilian Schecher, Adrian Kummerlaender
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net
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

#include <olb.h>

using namespace olb;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D2Q9<descriptors::FORCE,
  FreeSurface::MASS,
  FreeSurface::EPSILON,
  FreeSurface::CELL_TYPE,
  FreeSurface::CELL_FLAGS,
  FreeSurface::TEMP_MASS_EXCHANGE,
  FreeSurface::PREVIOUS_VELOCITY>>
>;

/// @brief Create a simulation mesh, based on user-specified geometry
/// @return An instance of a mesh with the relevant information
Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;

  const Vector<T,2> extend( parameters.get<parameters::DOMAIN_EXTENT>()[0], parameters.get<parameters::DOMAIN_EXTENT>()[1] );
  const Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid( extend, origin );

  T characteristic_length = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T physDeltaX = characteristic_length / parameters.get<parameters::RESOLUTION>();

  Mesh<T,MyCase::d> mesh = Mesh<MyCase::value_t, MyCase::d>(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(3);

  return mesh;
}

/// @brief Define which cells should be filled with liquid at the beginning
template<typename T, typename DESCRIPTOR>
class FreeSurfaceBreakingDam2D final : public AnalyticalF2D<T,T> {
private:
  T lattice_size;
  std::array<T,3> cell_values;
  std::array<T, 2> area;
public:
  FreeSurfaceBreakingDam2D(T lattice_size_, const std::array<T,3>& cell_vals, const std::array<T,2>& area_):AnalyticalF2D<T,T>{1}, lattice_size{lattice_size_}, cell_values{cell_vals}, area{area_}{}

  bool operator()(T output[], const T x[]) override {
    output[0] = cell_values[0];

    if(x[1] <= area[1] * 0.6 && x[0] <= area[0] * 0.5){
      output[0] = cell_values[2];
    }else if(x[1] - lattice_size * 1.1 <= area[1] * 0.6 && x[0] - lattice_size * 1.1 <= area[0] * 0.5){
      output[0] = cell_values[1];
    }

    return true;
  }
};

/// @brief Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
/// @note The material numbers will be used to assign physics to lattice nodes
void prepareGeometry( MyCase& myCase ) {
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  auto& geometry = myCase.getGeometry();

  geometry.rename( 0,2 );
  geometry.rename( 2,1,{1,1} );

  geometry.clean();
  geometry.innerClean();
  geometry.checkForErrors();

  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

/// @brief Set lattice dynamics
/// @param myCase The Case instance which keeps the simulation data
void prepareLattice( MyCase& myCase ) {
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  const int resolution = parameters.get<parameters::RESOLUTION>();
  const T latticeRelaxationTime = parameters.get<parameters::LATTICE_RELAXATION_TIME>();
  const T char_phys_length = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T char_phys_vel = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T viscosity = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T density = parameters.get<parameters::PHYS_CHAR_DENSITY>();

  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR>>(
    int {resolution},            // resolution: number of voxels per charPhysL
    (T) latticeRelaxationTime,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T) char_phys_length,        // charPhysLength: reference length of simulation geometry
    (T) char_phys_vel,           // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T) viscosity,               // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T) density                  // physDensity: physical density in __kg / m^3__
  );

  dynamics::set<SmagorinskyForcedBGKdynamics>(lattice, geometry.getMaterialIndicator({1}));
  boundary::set<boundary::BounceBack>(lattice, geometry, 2);
  //setSlipBoundary<T,DESCRIPTOR>(lattice, geometry, 2);

  lattice.setParameter<descriptors::OMEGA>(lattice.getUnitConverter().getLatticeRelaxationFrequency());
  lattice.setParameter<collision::LES::SMAGORINSKY>(T(0.2));

  std::array<T,2> area;
  area[0] = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  area[1] = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  FreeSurfaceBreakingDam2D<T,DESCRIPTOR> cells_analytical{ lattice.getUnitConverter().getPhysDeltaX(), {0., 1., 2.}, area};
  FreeSurfaceBreakingDam2D<T,DESCRIPTOR> mass_analytical{ lattice.getUnitConverter().getPhysDeltaX(), {0., 0.5, 1.}, area};

  AnalyticalConst2D<T,T> force_zero{0., 0.};

  fields::set<FreeSurface::MASS>(lattice, geometry.getMaterialIndicator({0,1,2}), 0.);
  fields::set<FreeSurface::EPSILON>(lattice, geometry.getMaterialIndicator({0,1,2}), 0.);
  fields::set<FreeSurface::CELL_TYPE>(lattice, geometry.getMaterialIndicator({0,1,2}), 0.);
  fields::set<FreeSurface::CELL_FLAGS>(lattice, geometry.getMaterialIndicator({0,1,2}), 0.);
  fields::set<descriptors::FORCE>(lattice, geometry.getMaterialIndicator({0,1,2}), force_zero);

  fields::set<FreeSurface::CELL_TYPE>(lattice, geometry.getMaterialIndicator({1}), cells_analytical);
  fields::set<FreeSurface::MASS>(lattice, geometry.getMaterialIndicator({1}), mass_analytical);
  fields::set<FreeSurface::EPSILON>(lattice, geometry.getMaterialIndicator({1}), mass_analytical);

  //fields::set<FreeSurface::MASS>(lattice, geometry.getMaterialIndicator({0,2}), 1.);
  fields::set<FreeSurface::EPSILON>(lattice, geometry.getMaterialIndicator({0,2}), 1.);
  fields::set<FreeSurface::CELL_TYPE>(lattice, geometry.getMaterialIndicator({0,2}), 4.);


  static auto fluidIndicator = geometry.getMaterialIndicator({1});

  std::array<T,2> gravity;
  gravity[0] = parameters.get<parameters::GRAVITY>()[0];
  gravity[1] = parameters.get<parameters::GRAVITY>()[1];

  T force_factor = T(1) / lattice.getUnitConverter().getConversionFactorForce() * lattice.getUnitConverter().getConversionFactorMass();
  AnalyticalConst2D<T,T> force_a{gravity[0] * force_factor, gravity[1] * force_factor};
  fields::set<descriptors::FORCE>(lattice, fluidIndicator, force_a);

  T surface_tension_coefficient_factor = std::pow(lattice.getUnitConverter().getConversionFactorTime(), 2) / (density * std::pow(lattice.getUnitConverter().getPhysDeltaX(),3));

  static FreeSurface2DSetup<T,DESCRIPTOR> free_surface_setup{lattice};
  free_surface_setup.addPostProcessor();

  const bool drop_isolated_cells = parameters.get<FreeSurface::DROP_ISOLATED_CELLS>();
  const bool has_surface_tension = parameters.get<FreeSurface::HAS_SURFACE_TENSION>();
  const double surface_tension_coefficient = parameters.get<FreeSurface::SURFACE_TENSION_PARAMETER>();
  const double transitionThreshold = parameters.get<FreeSurface::TRANSITION>();
  const double lonelyThreshold = parameters.get<FreeSurface::LONELY_THRESHOLD>();

  // Set variables from freeSurfaceHelpers.h
  lattice.setParameter<FreeSurface::DROP_ISOLATED_CELLS>(drop_isolated_cells);
  lattice.setParameter<FreeSurface::TRANSITION>(transitionThreshold);
  lattice.setParameter<FreeSurface::LONELY_THRESHOLD>(lonelyThreshold);
  lattice.setParameter<FreeSurface::HAS_SURFACE_TENSION>(has_surface_tension);
  lattice.setParameter<FreeSurface::SURFACE_TENSION_PARAMETER>(surface_tension_coefficient_factor * surface_tension_coefficient);
  //lattice.setParameter<FreeSurface::FORCE_DENSITY>({gravity[0] * force_factor, gravity[1] * force_factor});

  clout << "Prepare Lattice ... OK" << std::endl;
}

/// Set initial condition for primal variables (velocity and density)
/// @param myCase The Case instance which keeps the simulation data
/// @note Be careful: initial values have to be set using lattice units
void setInitialValues( MyCase& myCase ){
  OstreamManager clout( std::cout,"setInitialValues" );

  auto& lattice = myCase.getLattice(NavierStokes{});

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();

  AnalyticalConst2D<T,T> u{0., 0.};
  AnalyticalConst2D<T,T> one(1.);

  lattice.defineRhoU( geometry.getMaterialIndicator({0,1,2}), one, u );
  for (int i: {0,1,2}) {
    lattice.iniEquilibrium( geometry, i, one, u );
  }

  // Set up free surface communicator stages
  FreeSurface::initialize(lattice);
  // Make the lattice ready for simulation
  lattice.initialize();
}

/// Update boundary values at times (and external fields, if they exist)
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
/// @note Be careful: boundary values have to be set using lattice units
void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{
  // Nothing to do here, because simulation does not depend on time
}


/// Compute simulation results at times
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step, timer
void getResults (MyCase& myCase, int iT, util::Timer<MyCase::value_t>& timer)
{
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();

  const int vtmIter  = 100;
  const int statIter = 100;

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperVTMwriter2D<T> vtmWriter( "breakingDam2d" );
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid(lattice);
    SuperLatticeRank2D<T, DESCRIPTOR> rank(lattice);
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }

  // Writes output for ParaView
  if ( iT%vtmIter==0 ) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperVTMwriter2D<T> vtmWriter( "breakingDam2d" );
    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity( lattice, converter );
    SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure( lattice, converter );
    SuperLatticeField2D<T, DESCRIPTOR, FreeSurface::EPSILON> epsilon( lattice );
    SuperLatticeField2D<T, DESCRIPTOR, FreeSurface::CELL_TYPE> cells( lattice );
    SuperLatticeField2D<T, DESCRIPTOR, FreeSurface::MASS> mass( lattice );
    epsilon.getName() = "epsilon";
    cells.getName() = "cell_type";
    mass.getName() = "mass";
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    vtmWriter.addFunctor( epsilon );
    vtmWriter.addFunctor( cells );
    vtmWriter.addFunctor( mass );

    vtmWriter.write(iT);
  }

  // Writes console outputs
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    lattice.getStatistics().print( iT, converter.getPhysTime( iT ) );
  }
}


/// @brief Execute simulation: set initial values and run time loop
/// @param myCase The Case instance which keeps the simulation data
void simulate( MyCase& myCase ){

  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  // Main Loop with Timer
  OstreamManager clout( std::cout,"starting simulation..." );

  T physTime = parameters.get<parameters::PHYS_MAX_T>();
  util::Timer<MyCase::value_t> timer( lattice.getUnitConverter().getLatticeTime( physTime ), geometry.getStatistics().getNvoxel() );
  timer.start();

  for ( std::size_t iT = 0; iT < lattice.getUnitConverter().getLatticeTime( physTime ); ++iT ) {
    getResults( myCase, iT, timer );
    lattice.collideAndStream();
  }

  timer.stop();
  timer.printSummary();

}

/// Setup and run a simulation
int main(int argc, char **argv)
{
  initialize(&argc, &argv, false, false);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION>(500);
    myCaseParameters.set<PHYS_CHAR_VELOCITY>(0.1);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(1e-4);
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1e3);
    myCaseParameters.set<PHYS_MAX_T>(30.);
    myCaseParameters.set<DOMAIN_EXTENT>({7.31, 0.42});
    myCaseParameters.set<GRAVITY>({0., -9.81});
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(.501);
    myCaseParameters.set<FreeSurface::DROP_ISOLATED_CELLS>(true);
    myCaseParameters.set<FreeSurface::HAS_SURFACE_TENSION>(true);
    myCaseParameters.set<FreeSurface::SURFACE_TENSION_PARAMETER>(0.05);
    myCaseParameters.set<FreeSurface::TRANSITION>(1e-3);
    myCaseParameters.set<FreeSurface::LONELY_THRESHOLD>(1.0);
  }

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry( myCase );

  /// === Step 6: Prepare Lattice ===
  prepareLattice( myCase );

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);

}