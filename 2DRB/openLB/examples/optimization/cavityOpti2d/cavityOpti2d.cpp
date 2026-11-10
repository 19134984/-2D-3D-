/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Shota Ito, Julius Jessberger
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


/** @file cavity2d.cpp
 * @brief In this example, a parameter identification problem is solved via numerical optimization.
 *
 * In a lid-driven cavity simulation, the goal is do find the unknown viscosity, s.t. experimental
 * results from Ghia et al. are matched.
 * The state variable (f) are the lattice populations.
 * The control variable (alpha) is 1.e+3 * viscosity.
 * The side condition (G(f, alpha)=0) is the BGK-Boltzmann equation plus boundary conditions,
 * it is implemented in cavity2d.h. For given alpha, there is a corresponding state f=f(alpha),
 * s.t. the side condition is fulfilled.
 * The objective functional is the relative error of the x-component of the velocity, i.e.,
 * J(f) = sum_{i=0}^{2} |u(f(p_i)) - u*(p_i)|^2 / sum_{i=0}^{2} |u*(p_i)|^2,
 * where p_i are probe points, u is the x-component of the simulated velocity, u* is the x-component
 * of the simulated velocity.
 *
 * As every OpenLB simulation, it consists of eight steps:
 * @li Step 1: Declarations
 * @li Step 2: Initialization
 * @li Step 3: Create mesh
 * @li Step 4: Create case
 * @li Step 5: Prepare geometry
 * @li Step 6: Prepare lattice
 * @li Step 7: Set initial values
 * @li Step 8: Simulate
 *
 * Here, it shows the basic structure of numerical optimization in OpenLB and should be suitable for beginners.
 *
 * As every OpenLB optimization example, it consists of additional six steps:
 * @li Step A: Declarations
 * @li Step B: Initialization (G)
 * @li Step C: Set initial control value (alpha)
 * @li Step D: Define objective routine (J)
 * @li Step E: Create optimizer
 * @li Step F: Optimize
 *
 * For more information, we refer to the OpenLB user guide.
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::opti;

/// @brief Step 1: Declare simulation structure.
/// Model name and lattice type are collected in a Case class
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D2Q9<>>
>;
/// @brief Step A: Declarations
/// State the type(s) of the used simulation cases and give a name to each case
using MyOptiCase = OptiCaseCDQ<
  Controlled, MyCase
>;

/// @brief Step 3: Create a simulation mesh, based on user-specific geometry
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

/// @brief Step 5: Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
/// @note The material numbers are used to assign physics to lattice nodes
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

/// @brief Step 6: Set lattice dynamics and boundary conditions
/// @param myCase The Case instance which keeps the simulation data
void prepareLattice(MyCase& myCase)
{
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

  /// Set up a unit converter with the characteristic physical units
  lattice.setUnitConverter<UnitConverter<T,DESCRIPTOR>>(
    physDeltaX,        // physDeltaX: spacing between two lattice cells in [m]
    physDeltaT,        // physDeltaT: time step in [s]
    physLength,        // charPhysLength: reference length of simulation geometry in [m]
    physCharVelocity,  // charPhysVelocity: highest expected velocity during simulation in [m/s]
    physCharViscosity, // physViscosity: physical kinematic viscosity in [m^2/s]
    physCharDensity    // physDensity: physical density [kg/m^3]
  );
  lattice.getUnitConverter().print();

  /// Material=1 --> bulk dynamics
  dynamics::set<ConstRhoBGKdynamics>(lattice, myCase.getGeometry().getMaterialIndicator({1}));

  /// Material=2,3 --> bulk dynamics, velocity boundary
  boundary::set<boundary::InterpolatedVelocity, ConstRhoBGKdynamics>(lattice, myCase.getGeometry(), 2);
  boundary::set<boundary::InterpolatedVelocity, ConstRhoBGKdynamics>(lattice, myCase.getGeometry(), 3);

  lattice.setParameter<descriptors::OMEGA>(lattice.getUnitConverter().getLatticeRelaxationFrequency());
}

/// Step 7: Set initial values for primal variables (e.g. velocity, density) and additional fields
/// @param myCase The Case instance which keeps the simulation data
/// @note Initial values have to be set using lattice units
void setInitialValues(MyCase& myCase)
{
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});

  /// Initialize density to one everywhere
  /// Initialize velocity to be tangential at the lid and zero in the bulk and at the walls
  AnalyticalConst2D<T,T> uLid(lattice.getUnitConverter().getCharPhysVelocity(), 0);

  /// Assign the non-zero velocity to the lid
  momenta::setVelocity(lattice, myCase.getGeometry().getMaterialIndicator({3}), uLid);
  lattice.initialize();
}

/// Step 8.1: Update boundary values at times (and additional fields, if needed)
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
/// @note Boundary values have to be set using lattice units
void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{
  // Nothing to do here, because simulation does not depend on time
}

/// Step 8.3: Compute simulation results at times
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  /// Write vtk plots every 0.3 seconds (of phys. simulation time)
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();

  const T physMaxT = parameters.get<parameters::PHYS_MAX_T>();
  const std::size_t iTlog = converter.getLatticeTime(physMaxT/20.);
  const std::size_t iTvtk = converter.getLatticeTime(physMaxT/100.);

  SuperVTMwriter2D<T> vtmWriter("cavity2d");
  SuperLatticePhysVelocity2D velocity(lattice, converter);
  SuperLatticePhysPressure2D pressure(lattice, converter);
  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(pressure);

  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  // Writes the VTK files
  if (iT%iTvtk == 0 && iT > 0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }

  /// Print some (numerical and computational) statistics
  if (iT%iTlog == 0) {
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
    timer.print(iT);
  }
}

/// Step 8: Execute simulation
/// @param myCase The Case instance which keeps the simulation data
/// Run time loop
void simulate(MyCase& myCase)
{
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  const T physMaxT = parameters.get<parameters::PHYS_MAX_T>();

  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(physMaxT);
  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT)
  {
    /// @li Step 8.1: Update the Boundary Values and Fields at Times
    setTemporalValues(myCase, iT);

    /// @li Step 8.2: Collide and Stream Execution
    myCase.getLattice(NavierStokes{}).collideAndStream();

    /// @li Step 8.3: Computation and Output of the Results
    getResults(myCase, timer, iT);
  }

  timer.stop();
  timer.printSummary();
}

/// Step C: Set initial control value
/// @param optiCase The optiCase instance which keeps the optimization data
void setInitialControl(MyOptiCase& optiCase) {
  using T = MyOptiCase::value_t;
  const T physViscosity = optiCase.getCase(Controlled{}).getParameters().get<parameters::PHYS_CHAR_VISCOSITY>();
  std::vector<T> initialControls({physViscosity*1000});
  optiCase.getController().set(initialControls);
}

/// Step D.1: Apply the control values in the simulation
/// @param optiCase The optiCase instance which keeps the optimization data
void applyControl(MyOptiCase& optiCase) {
  auto& control = optiCase.getController();
  /// @li Employ the control values in the simulation
  optiCase.getCase(Controlled{}).getParameters().set<parameters::PHYS_CHAR_VISCOSITY>(control[0]/1000);
}

MyOptiCase::value_t objectiveF(MyOptiCase& optiCase) {
  using T = MyOptiCase::value_t;
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& lattice = controlledCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();
  lattice.setProcessingContext(ProcessingContext::Evaluation);

  /// @li Step D.2: Compute and return the objective value on the basis of the simulation results
  /// @li Load experimental reference data
  const Vector<T,3> y_coord({122, 64, 13});
  const Vector<T,3> vel_ghia_RE100({0.68717, -0.20581, -0.06434});

  /// @li Get point probes from the simulation
  Vector<T,3> vel_simulation;
  SuperLatticePhysVelocity2D velocityF(lattice, converter);
  AnalyticalFfromSuperF2D<T> interpolationF(velocityF, true, 1);
  for (std::size_t y=0; y < y_coord.size(); ++y) {
    T position[2] = {0.5, y_coord[y]/T(128)};
    T velocity[2] = {T(), T()};
    interpolationF(velocity, position);
    vel_simulation[y] = velocity[0];
  }

  /// @li Compute relative, squared norm error
  const T relError = norm(vel_simulation - vel_ghia_RE100) / norm(vel_ghia_RE100);
  return util::pow(relError, 2.0);
}

/// Step D: Define objective routine
/// @param optiCase The optiCase instance which keeps the optimization data
/// @return The objective value/ value of the goal functional
/// Define the routine to get the objective from the updated controls
/// In this case, we compare the simulated velocity with experimental results from Ghia et al.
MyOptiCase::value_t computeObjective(MyOptiCase& optiCase) {
  /// @par Contents
  /// @li Store references for simplified access
  auto& controlledCase = optiCase.getCase(Controlled{});

  /// @li Clear old simulation data
  controlledCase.resetLattices();

  /// @li Step D.1: Communicate control values to the simulation setup
  applyControl(optiCase);

  /// @li Prepare and run simulation (as in non-optimization examples)
  prepareGeometry(controlledCase);
  prepareLattice(controlledCase);
  setInitialValues(controlledCase);
  simulate(controlledCase);

  // Evaluate objective functor to compute objective value
  return objectiveF(optiCase);
}

/// Steps B-G: Set up and run the optimization routine
int main(int argc, char* argv[])
{
  // Step 2: Initialization
  initialize(&argc, &argv);

  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<PHYS_DELTA_X       >(  0.015625);
    myCaseParameters.set<PHYS_DELTA_T       >( 0.0015625);
    myCaseParameters.set<PHYS_CHAR_VELOCITY >(       1.0);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(     0.002);
    myCaseParameters.set<PHYS_CHAR_DENSITY  >(       1.0);
    myCaseParameters.set<DOMAIN_EXTENT      >({1.0, 1.0});
    myCaseParameters.set<PHYS_MAX_T         >(       30.);
  }

  OptimizerParametersD<MyCase::value_t> optiParameters;
  {
    using namespace parameters;
    optiParameters.set<optimizer::EPS>(1e-7);
    optiParameters.set<optimizer::MAX_OPTIMIZATION_STEPS>(20);
    optiParameters.set<linesearch::MAX_STEP_ATTEMPTS>(10);
    optiParameters.set<linesearch::STEP_CONDITION>("StrongWolfe");
    optiParameters.set<optimizer::VERBOSE>(true);
  }

  // Step 3: Create Mesh
  Mesh mesh = createMesh(myCaseParameters);

  // Step 4: Create Case
  MyCase myCase(myCaseParameters, mesh);

  // ==== BELOW HERE OPTIMIZATION SPECIFIC ====
  /// @li Step B: Create OptiCase
  MyOptiCase optiCase;
  optiCase.setCase<Controlled>(myCase);

  /// @li Step C: Set initial control value
  setInitialControl(optiCase);

  /// @li Step D: Define objective routine
  optiCase.setObjective(computeObjective);

  /// @li Step E: Create an Optimizer
  OptimizerLBFGS<MyCase::value_t,std::vector<MyCase::value_t>> optimizer(
    optiCase.getController().size(), optiParameters,
    {OptimizerLogType::objective, OptimizerLogType::control, OptimizerLogType::derivative});

  /// @li Step F: Optimize
  optimizer.optimize(optiCase);

  OstreamManager clout(std::cout, "cavity2dOpti");
  clout << "Reynolds number in experiments = 100" << std::endl;
  clout << "Reynolds number found by optimization algorithm = "
        << optiCase.getCase(Controlled{}).getLattice(NavierStokes{}).getUnitConverter().getReynoldsNumber()
        << std::endl;
}
