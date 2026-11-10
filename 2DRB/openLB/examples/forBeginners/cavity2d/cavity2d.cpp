/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Shota Ito, Julius Jessberger, Adrian Kummerländer, Mathias J. Krause
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

#include <olb.h>

using namespace olb;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D2Q9<>>
>;

/// @brief Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector origin{0, 0};
  IndicatorCuboid2D<T> cuboid(extent, origin);

  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

/// @brief Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
/// @note The material numbers will be used to assign physics to lattice nodes
void prepareGeometry(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();

  const T physLength = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();

  /// Set material numbers
  geometry.rename(0, 2);
  geometry.rename(2, 1, {1,1});

  const Vector lidExtend {physLength + 2*physDeltaX, 2*physDeltaX};
  const Vector lidOrigin {-physDeltaX, physLength - physDeltaX};
  IndicatorCuboid2D lid (lidExtend, lidOrigin);
  geometry.rename(2,3,1,lid);

  geometry.getStatistics().print();
}

/// @brief Set lattice dynamics
/// @param myCase The Case instance which keeps the simulation data
void prepareLattice(MyCase& myCase) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});

  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();
  const T physDeltaT = parameters.get<parameters::PHYS_DELTA_T>();
  const T physLength = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T physCharVelocity = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T physCharViscosity = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T physCharDensity = parameters.get<parameters::PHYS_CHAR_DENSITY>();

  // Set up a unit converter with the characteristic physical units
  myCase.getLattice(NavierStokes{}).setUnitConverter<UnitConverter<T,DESCRIPTOR>>(
    physDeltaX,        // physDeltaX: spacing between two lattice cells in [m]
    physDeltaT,        // physDeltaT: time step in [s]
    physLength,        // charPhysLength: reference length of simulation geometry in [m]
    physCharVelocity,  // charPhysVelocity: highest expected velocity during simulation in [m/s]
    physCharViscosity, // physViscosity: physical kinematic viscosity in [m^2/s]
    physCharDensity    // physDensity: physical density [kg/m^3]
  );
  lattice.getUnitConverter().print();

  /// Material=1 -->bulk dynamics
  lattice.defineDynamics<ConstRhoBGKdynamics>(myCase.getGeometry(), 1);

  /// Material=2,3 -->bulk dynamics, velocity boundary
  boundary::set<boundary::InterpolatedVelocity, ConstRhoBGKdynamics>(lattice, myCase.getGeometry(), 2);
  boundary::set<boundary::InterpolatedVelocity, ConstRhoBGKdynamics>(lattice, myCase.getGeometry(), 3);

  lattice.setParameter<descriptors::OMEGA>(lattice.getUnitConverter().getLatticeRelaxationFrequency());
}

/// Set initial condition for primal variables (velocity and density)
/// @param myCase The Case instance which keeps the simulation data
/// @note Be careful: initial values have to be set using lattice units
void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});

  /// Initialize density to one everywhere
  /// Initialize velocity to be tangential at the lid and zero in the bulk and at the walls
  AnalyticalConst2D<T,T> rho(1);
  AnalyticalConst2D<T,T> uLid(lattice.getUnitConverter().getCharLatticeVelocity(), 0);
  AnalyticalConst2D<T,T> uRest(0, 0);

  /// Initialize populations to equilibrium state
  auto domain = myCase.getGeometry().getMaterialIndicator({1,2,3});
  lattice.iniEquilibrium(domain, rho, uRest);
  lattice.defineRhoU(domain, rho, uRest);

  /// Assign the non-zero velocity to the lid
  lattice.defineU(myCase.getGeometry(), 3, uLid);
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
/// @param iT The time step
void getResults(MyCase& myCase, std::size_t iT)
{
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();

  const T physMaxT = parameters.get<parameters::PHYS_MAX_T>();
  const std::size_t iTvtk = converter.getLatticeTime(physMaxT/100.);

  SuperVTMwriter2D<T> vtmWriter("cavity2d");
  SuperLatticePhysVelocity2D velocity(lattice, converter);
  SuperLatticePhysPressure2D pressure(lattice, converter);
  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(pressure);

  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  if (iT%iTvtk == 0 && iT > 0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }

  myCase.report(iT);
}


/// @brief Execute simulation: set initial values and run time loop
/// @param myCase The Case instance which keeps the simulation data
void simulate(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  const T physMaxT = parameters.get<parameters::PHYS_MAX_T>();

  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(physMaxT);

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    myCase.getLattice(NavierStokes{}).collideAndStream();

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, iT);
  }

  myCase.finalize(iTmax);
}

/// Setup and run a simulation
int main(int argc, char* argv[]) {
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<PHYS_DELTA_X       >( 0.0078125);
    myCaseParameters.set<PHYS_DELTA_T       >(0.00078125);
    myCaseParameters.set<PHYS_CHAR_VELOCITY >(       1.0);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(     0.001);
    myCaseParameters.set<PHYS_CHAR_DENSITY  >(       1.0);
    myCaseParameters.set<DOMAIN_EXTENT      >({1.0, 1.0});
    myCaseParameters.set<PHYS_MAX_T         >(       30.);
    myCaseParameters["reporter"]["statistics"].set<PHYS_ITER_T>(1.5);
    myCaseParameters["reporter"]["timer"].set<PHYS_ITER_T>(1.5);
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
