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
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <memory>
#include <fstream>
#include <sstream>

using namespace olb;
using namespace olb::names;
using namespace olb::descriptors;

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

  Vector<T,3> extend(parameters.get<parameters::DOMAIN_EXTENT>()[0], parameters.get<parameters::DOMAIN_EXTENT>()[1],parameters.get<parameters::DOMAIN_EXTENT>()[2]);
  Vector<T,3> origin;
  IndicatorCuboid3D<T> cuboid (extend, origin);

  T characteristic_length = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  T physDeltaX = characteristic_length / parameters.get<parameters::RESOLUTION>();

  Mesh<T, MyCase::d> mesh = Mesh<MyCase::value_t, MyCase::d>(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(3);

  return mesh;
}

/// @brief Define which cells should be filled with liquid at the beginning
template <typename T, typename DESCRIPTOR>
class FreeSurfaceFallingDrop3D final : public AnalyticalF3D<T, T> {
  private:
    T lattice_size;
    // This value doesn't depend on the dimension of the descriptor. It's always 3
    std::array<T, 3> cell_values;
  public:
    FreeSurfaceFallingDrop3D(T lattice_size_, const std::array<T, 3>& cell_vals) :AnalyticalF3D<T, T>{ 1 }, lattice_size{ lattice_size_ }, cell_values{ cell_vals }{}

    bool operator()(T output[], const T x[]) override {
      output[0] = cell_values[0];
      T radius = 0.00155;

      if (x[2] <= radius) {
        output[0] = cell_values[2];
        }
      else if (x[2] <= radius + lattice_size * 1.5) {
        output[0] = cell_values[1];
        }

      std::array<T, DESCRIPTOR::d> point = { 0.015, 0.015, 2 * radius + lattice_size * 4 };
      std::array<T, DESCRIPTOR::d> diff = { x[0] - point[0], x[1] - point[1], x[2] - point[2] };

      if ((diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]) <= radius * radius) {
        output[0] = cell_values[2];
        }
        else {
            for (int i = -1; i <= 1; ++i) {
                for (int j = -1; j <= 1; ++j) {
                    for (int k = -1; k <= 1; ++k) {
                        std::array<T, DESCRIPTOR::d> shifted_diff = { diff[0] + i * lattice_size * T{1.1}, diff[1] + j * lattice_size * T{1.1}, diff[2] + k * lattice_size * T{1.1} };
                        if ((shifted_diff[0] * shifted_diff[0] + shifted_diff[1] * shifted_diff[1] + shifted_diff[2] * shifted_diff[2]) <= radius * radius) {
                            output[0] = cell_values[1];
                            return true;
                        }
                    }
                }
            }
        }

      return true;
    }
};

/// @brief Define the initial state of the falling drop
template <typename T, typename DESCRIPTOR>
class FallingDropVel3D final : public AnalyticalF<3, T, T> {
private:
    T lattice_size;
    std::array<T, DESCRIPTOR::d> lattice_speed;
public:
    FallingDropVel3D(T lattice_size_, const std::array<T, DESCRIPTOR::d>& lattice_speed_) :AnalyticalF<3, T, T>{ 3 }, lattice_size{ lattice_size_ }, lattice_speed{ lattice_speed_ }{}

    bool operator()(T output[], const T x[]) override {

        T radius = 0.00155;
        std::array<T, DESCRIPTOR::d> point = { 0.015, 0.015, 2 * radius + lattice_size * 2 };
        std::array<T, DESCRIPTOR::d> diff = { x[0] - point[0], x[1] - point[1], x[2] - point[2] };

        output[0] = 0.;
        output[1] = 0.;
        output[2] = 0.;
        for (int i = -1; i <= 1; ++i) {
            for (int j = -1; j <= 1; ++j) {
                for (int k = -1; k <= 1; ++k) {
                    std::array<T, DESCRIPTOR::d> shifted_diff = { diff[0] + i * lattice_size * T{1.1}, diff[1] + j * lattice_size * T{1.1},diff[2] + k * lattice_size * T{1.1} };
                    if ((shifted_diff[0] * shifted_diff[0] + shifted_diff[1] * shifted_diff[1] + shifted_diff[2] * shifted_diff[2]) <= radius * radius) {
                        output[0] = lattice_speed[0];
                        output[1] = lattice_speed[1];
                        output[2] = lattice_speed[2];
                        return true;
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
void prepareGeometry(MyCase& myCase) {

    OstreamManager clout(std::cout, "prepareGeometry");
    clout << "Prepare Geometry ..." << std::endl;

    auto& geometry = myCase.getGeometry();

    geometry.rename(0, 2);
    geometry.rename(2, 1, { 1,1,1 });

    geometry.clean();
    geometry.innerClean();
    geometry.checkForErrors();

    geometry.print();

    clout << "Prepare Geometry ... OK" << std::endl;
}

/// @brief Set up the falling drop with the freeSurface information
void prepareFallingDrop(MyCase& myCase)
{
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  const T char_phys_length = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const int resolution = parameters.get<parameters::RESOLUTION>();
  std::array<T,3> gravity_force;
  gravity_force[0] = parameters.get<parameters::GRAVITY>()[0];
  gravity_force[1] = parameters.get<parameters::GRAVITY>()[1];
  gravity_force[2] = parameters.get<parameters::GRAVITY>()[2];
  auto& converter = lattice.getUnitConverter();

  const T lattice_size = char_phys_length / resolution;

  FreeSurfaceFallingDrop3D<T, DESCRIPTOR> cells_analytical{ lattice_size, {0., 1., 2.} };
  FreeSurfaceFallingDrop3D<T, DESCRIPTOR> mass_analytical{ lattice_size, {0., 0.5, 1.} };

  AnalyticalConst3D<T, T> force_zero{ 0., 0., 0. };

  fields::set<FreeSurface::MASS>(lattice, geometry.getMaterialIndicator({0,1,2}), 0.);
  fields::set<FreeSurface::EPSILON>(lattice, geometry.getMaterialIndicator({0,1,2}), 0.);
  fields::set<FreeSurface::CELL_TYPE>(lattice, geometry.getMaterialIndicator({0,1,2}), 0.);
  fields::set<FreeSurface::CELL_FLAGS>(lattice, geometry.getMaterialIndicator({0,1,2}), 0.);
  fields::set<descriptors::FORCE>(lattice, geometry.getMaterialIndicator({0,1,2}), force_zero);

  fields::set<FreeSurface::CELL_TYPE>(lattice, geometry.getMaterialIndicator({1}), cells_analytical);
  fields::set<FreeSurface::MASS>(lattice, geometry.getMaterialIndicator({1}), mass_analytical);
  fields::set<FreeSurface::EPSILON>(lattice, geometry.getMaterialIndicator({1}), mass_analytical);

  fields::set<FreeSurface::MASS>(lattice, geometry.getMaterialIndicator({0,2}), 1.);
  fields::set<FreeSurface::EPSILON>(lattice, geometry.getMaterialIndicator({0,2}), 1.);
  fields::set<FreeSurface::CELL_TYPE>(lattice, geometry.getMaterialIndicator({0,2}), 4.);

  T force_factor = T(1) / converter.getConversionFactorForce() * converter.getConversionFactorMass();
  AnalyticalConst3D<T, T> force_a{ gravity_force[0] * force_factor, gravity_force[1] * force_factor, gravity_force[2] * force_factor };

  fields::set<descriptors::FORCE>(lattice, geometry.getMaterialIndicator({1}), force_a);
}

/// @brief Set lattice dynamics
/// @param myCase The Case instance which keeps the simulation data
void prepareLattice(MyCase& myCase) {
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const int resolution = parameters.get<parameters::RESOLUTION>();
  const T latticeRelaxationTime = parameters.get<parameters::LATTICE_RELAXATION_TIME>();
  const T char_phys_length = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T char_phys_vel = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T viscosity = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T density = parameters.get<parameters::PHYS_CHAR_DENSITY>();

  lattice.setUnitConverter < UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR>>(
    int{ resolution },          // resolution: number of voxels per charPhysL
    (T)latticeRelaxationTime,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)char_phys_length,        // charPhysLength: reference length of simulation geometry
    (T)char_phys_vel,           // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)viscosity,               // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)density                  // physDensity: physical density in __kg / m^3__
  );
  auto& converter = lattice.getUnitConverter();

  dynamics::set<SmagorinskyForcedBGKdynamics>(lattice, geometry.getMaterialIndicator({1}));
  boundary::set<boundary::BounceBack>(lattice, geometry, 2);
  //setSlipBoundary<T,DESCRIPTOR>(lattice, geometry, 2);

  lattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  lattice.setParameter<collision::LES::SMAGORINSKY>(T(0.2));

  prepareFallingDrop(myCase);

  const bool drop_isolated_cells = parameters.get<FreeSurface::DROP_ISOLATED_CELLS>();
  const bool has_surface_tension = parameters.get<FreeSurface::HAS_SURFACE_TENSION>();
  const T surface_tension_coefficient = parameters.get<FreeSurface::SURFACE_TENSION_PARAMETER>();
  const T transitionThreshold = parameters.get<FreeSurface::TRANSITION>();
  const T lonelyThreshold = parameters.get<FreeSurface::LONELY_THRESHOLD>();

  // Compute the scaling factor for surface tension
  T surface_tension_factor = std::pow(converter.getConversionFactorTime(), 2) / (density * std::pow(converter.getPhysDeltaX(), 3));

  static FreeSurface3DSetup<T, DESCRIPTOR> free_surface_setup{ lattice };
  free_surface_setup.addPostProcessor();

  lattice.setParameter<FreeSurface::DROP_ISOLATED_CELLS>(drop_isolated_cells);
  lattice.setParameter<FreeSurface::TRANSITION>(transitionThreshold);
  lattice.setParameter<FreeSurface::LONELY_THRESHOLD>(lonelyThreshold);
  lattice.setParameter<FreeSurface::HAS_SURFACE_TENSION>(has_surface_tension);
  lattice.setParameter<FreeSurface::SURFACE_TENSION_PARAMETER>(surface_tension_factor * surface_tension_coefficient);
  clout << "Prepare Lattice ... OK" << std::endl;
}


/// Set initial condition for primal variables (velocity and density)
/// @param myCase The Case instance which keeps the simulation data
/// @note Be careful: initial values have to be set using lattice units
void setInitialValues(MyCase& myCase){
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  auto& converter = lattice.getUnitConverter();

  OstreamManager clout( std::cout,"setInitialValues" );

  std::array<T,3> initial_falling_speed;
  initial_falling_speed[0] = parameters.get<parameters::INITIAL_FALLING_SPEED>()[0];
  initial_falling_speed[1] = parameters.get<parameters::INITIAL_FALLING_SPEED>()[1];
  initial_falling_speed[2] = parameters.get<parameters::INITIAL_FALLING_SPEED>()[2];

  std::array<T, DESCRIPTOR::d> lattice_speed;
  for (size_t i = 0; i < DESCRIPTOR::d; ++i) {
      lattice_speed[i] = initial_falling_speed[i] * converter.getPhysDeltaT() / converter.getPhysDeltaX();
  }

  T characteristic_length = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  T lattice_size = characteristic_length / parameters.get<parameters::RESOLUTION>();

  FallingDropVel3D<T, DESCRIPTOR> u{ lattice_size, lattice_speed };
  AnalyticalConst3D<T, T> one(1.);

  lattice.defineRhoU(geometry.getMaterialIndicator({ 0,1,2 }), one, u);
  for (int i : {0, 1, 2}) {
        lattice.iniEquilibrium(geometry, i, one, u);
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
void getResults(MyCase& myCase, int iT, util::Timer<MyCase::value_t>& timer)
{
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();
  auto& parameters = myCase.getParameters();
  const T physTime = parameters.get<parameters::PHYS_MAX_T>();
  OstreamManager clout(std::cout, "getResults");

    SuperVTMwriter3D<T> vtmWriter("fallingDrop3d");
    const int vtmIter = converter.getLatticeTime(physTime / 50.);
    const int statIter = converter.getLatticeTime(physTime / 100.);

    if (iT == 0) {
        // Writes the geometry, cuboid no. and rank no. as vti file for visualization
        SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(lattice);
        SuperLatticeRank3D<T, DESCRIPTOR> rank(lattice);
        vtmWriter.write(cuboid);
        vtmWriter.write(rank);

        vtmWriter.createMasterFile();
    }

    // Writes the vtm files and profile text file
    if (iT % vtmIter == 0) {
        lattice.setProcessingContext(ProcessingContext::Evaluation);

        SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(lattice, converter);
        SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(lattice, converter);
        SuperLatticeField3D<T, DESCRIPTOR, FreeSurface::EPSILON> epsilon(lattice);
        SuperLatticeField3D<T, DESCRIPTOR, FreeSurface::CELL_TYPE> cells(lattice);
        SuperLatticeField3D<T, DESCRIPTOR, FreeSurface::MASS> mass(lattice);
        epsilon.getName() = "epsilon";
        cells.getName() = "cell_type";
        mass.getName() = "mass";
        vtmWriter.addFunctor(velocity);
        vtmWriter.addFunctor(pressure);
        vtmWriter.addFunctor(epsilon);
        vtmWriter.addFunctor(cells);
        vtmWriter.addFunctor(mass);

        vtmWriter.write(iT);
    }

    // Writes output on the console
    if (iT % statIter == 0) {
        // Timer console output
        timer.update(iT);
        timer.printStep();

        // Lattice statistics console output
        lattice.getStatistics().print(iT, converter.getPhysTime(iT));
    }
}

/// @brief Execute simulation: set initial values and run time loop
/// @param myCase The Case instance which keeps the simulation data
void simulate( MyCase& myCase ){
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& converter = lattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  // Main Loop with Timer
  OstreamManager clout( std::cout,"starting simulation..." );

  T physTime = parameters.get<parameters::PHYS_MAX_T>();
    util::Timer<T> timer(converter.getLatticeTime(physTime), geometry.getStatistics().getNvoxel());
    timer.start();
    setInitialValues(myCase);

    for (std::size_t iT = 0; iT < converter.getLatticeTime(physTime); ++iT) {
        getResults(myCase, iT, timer);
        lattice.collideAndStream();
    }

    timer.stop();
    timer.printSummary();

}

/// Setup and run a simulation
int main(int argc, char** argv)
{
  initialize(&argc, &argv, false, false);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(1e-4); //viscosity
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1e3); //density
    myCaseParameters.set<PHYS_MAX_T>(0.01); //physTime
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(.516); //latticeRelaxationTime
    myCaseParameters.set<RESOLUTION>(100); // resolution: number of voxels per charPhysL
    myCaseParameters.set<DOMAIN_EXTENT>({0.03, 0.03, 0.025}); //formerly area
    myCaseParameters.set<GRAVITY>({0., 0., -9.81});
    myCaseParameters.set<PHYS_CHAR_VELOCITY>(0.1); //char_phys_vel
    myCaseParameters.set<INITIAL_FALLING_SPEED>({0., 0., -6.03});

    myCaseParameters.set<FreeSurface::DROP_ISOLATED_CELLS>(true);
    myCaseParameters.set<FreeSurface::HAS_SURFACE_TENSION>(true);
    myCaseParameters.set<FreeSurface::SURFACE_TENSION_PARAMETER>(0.0661);
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
