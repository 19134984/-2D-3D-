/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Claudius Holeksa
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
  NavierStokes, Lattice<double, descriptors::D3Q27<descriptors::FORCE,
  FreeSurface::MASS,
  FreeSurface::EPSILON,
  FreeSurface::CELL_TYPE,
  FreeSurface::CELL_FLAGS,
  FreeSurface::TEMP_MASS_EXCHANGE,
  FreeSurface::PREVIOUS_VELOCITY,
  FreeSurface::HAS_INTERFACE_NBRS>>
>;

/// @brief Create a simulation mesh, based on user-specified geometry
/// @return An instance of a mesh with the relevant information
Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;

  Vector<T,3> extend(parameters.get<parameters::DOMAIN_EXTENT>()[0], parameters.get<parameters::DOMAIN_EXTENT>()[1], parameters.get<parameters::DOMAIN_EXTENT>()[2]);
  Vector<T,3> origin;
  IndicatorCuboid3D<T> cuboid(extend, origin);

  T characteristic_length = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T physDeltaX = characteristic_length / parameters.get<parameters::RESOLUTION>();

  Mesh<T,MyCase::d> mesh = Mesh<MyCase::value_t, MyCase::d>(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(3);

  return mesh;
}


template <typename T, typename DESCRIPTOR>
class FreeSurfaceRayleighInstability3D final : public AnalyticalF<DESCRIPTOR::d, T,T> {
private:
  T lattice_size;
  std::array<T, DESCRIPTOR::d> area;
  // This value doesn't depend on the dimension of the descriptor. It's always 3
  std::array<T, 3> cell_values;
  bool interpolate;
public:
  FreeSurfaceRayleighInstability3D(T lattice_size_, const std::array<T,DESCRIPTOR::d>& area_, const std::array<T,3>& cell_vals, bool interpolate_):AnalyticalF<DESCRIPTOR::d,T,T>{1}, lattice_size{lattice_size_}, area{area_}, cell_values{cell_vals}, interpolate{interpolate_}{}

  bool operator()(T output[], const T x[]) override {
    output[0] = cell_values[0];

    T middle = area[1] * 0.5;

    std::array<T,DESCRIPTOR::d-1> diff {
      x[1] - middle,
      x[2] - middle
    };

    T length_2 = diff[0]*diff[0] + diff[1] * diff[1];

    // Number of bumps in the water jet
    T number = 7.;
    // Optimal decay values for $number of bumps
    T rad = area[0] / (number * 9.);
    T frequency = 2. * 3.141529 / (rad * 9.);

    T rad_2 = rad*rad;

    // Set the sinusoidal instabilities
    if(length_2 < rad_2){
      output[0] = cell_values[2];
    }else{
      for(int i = -1; i <= 1; ++i){
        for(int j = -1; j <= 1; ++j){
          for(int k = -1; k <= 1; ++k){
            std::array<T,DESCRIPTOR::d-1> s_diff{
              diff[0]+j*lattice_size,
              diff[1]+k*lattice_size
            };

            T s_len_2 = s_diff[0]*s_diff[0] + s_diff[1] * s_diff[1];

            if(s_len_2 < rad_2){
              if(interpolate){
                T inter = 0.5 + 0.45*std::sin(frequency * x[0]);

                output[0] = inter;

                return true;
              }else{
                output[0] = cell_values[1];
                return true;
              }
            }
          }
        }
      }
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
  geometry.rename( 2,1,{1,1,1} );

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

  auto& converter = lattice.getUnitConverter();

  T lattice_size = parameters.get<parameters::DOMAIN_EXTENT>()[0] / parameters.get<parameters::RESOLUTION>();
  T force_factor = T(1) / converter.getConversionFactorForce() * converter.getConversionFactorMass();

  // Convert kg / s^2
  // Basically it is multiplied with s^2 / kg = s^2 * m^3 / (kg * m^2 * m) = 1. / (velocity_factor^2 * density * length_factor)
  T surface_tension_coefficient_factor = std::pow(converter.getConversionFactorTime(),2)/ (density * std::pow(converter.getPhysDeltaX(),3));
  T surface_tension_coefficient = parameters.get<FreeSurface::SURFACE_TENSION_PARAMETER>();
  clout<<"Surface: "<<surface_tension_coefficient_factor * surface_tension_coefficient <<std::endl;
  clout<<"Lattice Size: "<<converter.getPhysDeltaX()<<std::endl;


  dynamics::set<SmagorinskyForcedBGKdynamics>(lattice, geometry.getMaterialIndicator({1}));
  boundary::set<boundary::BounceBack>(lattice, geometry, 2);
  //setSlipBoundary<T,DESCRIPTOR>(lattice, geometry, 2);

  lattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  lattice.setParameter<collision::LES::SMAGORINSKY>(T(0.2));

  // prepareRayleighInstability
  std::array<T,3> area;
  area[0] = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  area[1] = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  area[2] = parameters.get<parameters::DOMAIN_EXTENT>()[2];
  FreeSurfaceRayleighInstability3D<T,DESCRIPTOR> cells_analytical{ lattice_size, area, {0., 1., 2.}, false};
  FreeSurfaceRayleighInstability3D<T,DESCRIPTOR> mass_analytical{ lattice_size, area, {0., 0.5, 1.}, true};

  AnalyticalConst3D<T,T> force_zero{0., 0., 0.};

  // Set border values
  fields::set<FreeSurface::MASS>(lattice, geometry.getMaterialIndicator({0,2}), 0.);
  fields::set<FreeSurface::EPSILON>(lattice, geometry.getMaterialIndicator({0,2}), 0.);
  fields::set<FreeSurface::CELL_TYPE>(lattice, geometry.getMaterialIndicator({0,2}), 4.);

  // Set simulation area values
  fields::set<FreeSurface::MASS>(lattice, geometry.getMaterialIndicator({1}), mass_analytical);
  fields::set<FreeSurface::EPSILON>(lattice, geometry.getMaterialIndicator({1}), mass_analytical);
  fields::set<FreeSurface::CELL_TYPE>(lattice, geometry.getMaterialIndicator({1}), cells_analytical);

  // Define to zero
  fields::set<FreeSurface::CELL_FLAGS>(lattice, geometry.getMaterialIndicator({0,1,2}), 0.);
  fields::set<descriptors::FORCE>(lattice, geometry.getMaterialIndicator({0,1,2}), force_zero);


  // Needs to be set to a calculated velocity of the current velocity
  fields::set<FreeSurface::PREVIOUS_VELOCITY>(lattice, geometry.getMaterialIndicator({0,1,2}), force_zero);

  static FreeSurface3DSetup<T,DESCRIPTOR> free_surface_setup{lattice};
  free_surface_setup.addPostProcessor();

  bool drop_isolated_cells = parameters.get<FreeSurface::DROP_ISOLATED_CELLS>();
  bool has_surface_tension = parameters.get<FreeSurface::HAS_SURFACE_TENSION>();
  T transitionThreshold = parameters.get<FreeSurface::TRANSITION>();
  T lonelyThreshold = parameters.get<FreeSurface::LONELY_THRESHOLD>();
  std::array<T,2> gravity;
  gravity[0] = parameters.get<parameters::GRAVITY>()[0];
  gravity[1] = parameters.get<parameters::GRAVITY>()[1];
  gravity[2] = parameters.get<parameters::GRAVITY>()[2];

  // Set variables from freeSurfaceHelpers.h
  lattice.setParameter<FreeSurface::DROP_ISOLATED_CELLS>(drop_isolated_cells);
  lattice.setParameter<FreeSurface::TRANSITION>(transitionThreshold);
  lattice.setParameter<FreeSurface::LONELY_THRESHOLD>(lonelyThreshold);
  lattice.setParameter<FreeSurface::HAS_SURFACE_TENSION>(has_surface_tension);
  lattice.setParameter<FreeSurface::SURFACE_TENSION_PARAMETER>(surface_tension_coefficient_factor * surface_tension_coefficient);
  lattice.setParameter<FreeSurface::FORCE_DENSITY>({gravity[0] * force_factor, gravity[1] * force_factor, gravity[2] * force_factor});

  {
    auto& communicator = lattice.getCommunicator(stage::Full());
    communicator.requestField<FreeSurface::MASS>();
    communicator.requestField<FreeSurface::EPSILON>();
    communicator.requestField<FreeSurface::CELL_TYPE>();
    communicator.requestOverlap(2);
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}


/// Set initial condition for primal variables (velocity and density)
/// @param myCase The Case instance which keeps the simulation data
/// @note Be careful: initial values have to be set using lattice units
void setInitialValues( MyCase& myCase ){
  OstreamManager clout( std::cout,"setInitialValues" );

  auto& lattice = myCase.getLattice(NavierStokes{});

  // Set up free surface communicator stages
  FreeSurface::initialize(lattice);
  // Make the lattice ready for simulation
  lattice.initialize();
}


/// Update boundary values at times (and external fields, if they exist)
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
/// @note Be careful: boundary values have to be set using lattice units
void setTemporalValues(MyCase& myCase, std::size_t iT)
{
  // Nothing to do here, because simulation does not depend on time
}


/// Compute simulation results at times
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step, timer
void getResults(MyCase& myCase, int iT, util::Timer<MyCase::value_t>& timer)
{
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& converter = lattice.getUnitConverter();

  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter3D<T> vtmWriter( "rayleighInstability3d" );

  const int vtmIter  = 100;//converter.getLatticeTime( maxPhysT/2000. );
  const int statIter =  50;//converter.getLatticeTime( maxPhysT/2000. );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( lattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( lattice );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }

  // Writes the vtm files and profile text file
  if ( iT%vtmIter==0 ) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( lattice, converter );
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( lattice, converter );
    SuperLatticeField3D<T, DESCRIPTOR, FreeSurface::EPSILON> epsilon( lattice );
    SuperLatticeField3D<T, DESCRIPTOR, FreeSurface::CELL_TYPE> cells( lattice );
    SuperLatticeField3D<T, DESCRIPTOR, FreeSurface::MASS> mass( lattice );
    epsilon.getName() = "epsilon";
    cells.getName() = "cell_type";
    mass.getName() = "mass";
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    vtmWriter.addFunctor( epsilon );
    vtmWriter.addFunctor( cells );
    vtmWriter.addFunctor( mass );

    vtmWriter.write( iT );
  }

  // Writes output on the console
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    lattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }
}


/// @brief Execute simulation: set initial values and run time loop
/// @param myCase The Case instance which keeps the simulation data
void simulate(MyCase& myCase){

  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& converter = lattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  T physTime = parameters.get<parameters::PHYS_MAX_T>();

  // Main Loop with Timer
  std::cout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( physTime ), geometry.getStatistics().getNvoxel() );
  timer.start();

  for ( std::size_t iT = 0; iT < converter.getLatticeTime( physTime ); ++iT ) {
    getResults(myCase, iT, timer );
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
    myCaseParameters.set<RESOLUTION>(512);
    myCaseParameters.set<PHYS_CHAR_VELOCITY>(0.1);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(1./6.);
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1e0);
    myCaseParameters.set<PHYS_MAX_T>(500.);
    myCaseParameters.set<DOMAIN_EXTENT>({512., 50., 50.});
    myCaseParameters.set<GRAVITY>({0., 0., 0.});
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(.6);
    myCaseParameters.set<FreeSurface::DROP_ISOLATED_CELLS>(true);
    myCaseParameters.set<FreeSurface::HAS_SURFACE_TENSION>(true);
    myCaseParameters.set<FreeSurface::SURFACE_TENSION_PARAMETER>(1.0);
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
