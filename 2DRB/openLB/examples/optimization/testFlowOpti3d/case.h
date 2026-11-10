/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006-2025 Mathias J. Krause, Fabian Klemens,
 *  Julius Jessberger, Shota Ito
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
#include "analyticalSolutionTestFlow3D.h"

using namespace olb;
using namespace olb::names;
using namespace olb::opti;

/// @brief Step 1: Declare simulation structure.
/// Model name and lattice type are collected in a Case class
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q19<descriptors::FORCE>>
>;
/// State the type(s) of the used simulation cases and give a name to each case
using MyOptiCase = OptiCaseAdjoint<
  Controlled, MyCase,
  Adjoint, MyCase,
  Reference, MyCase
>;

using BulkDynamics = ForcedBGKdynamics<MyCase::value_t,
                                       MyCase::descriptor_t_of<NavierStokes>>;
using ObjectiveF = functors::L2DistanceF<functors::VelocityF>;
using ControlledField = descriptors::FORCE;

enum class AdjointType : int {
  discrete = 0,
  continuous = 1
};

namespace olb::parameters {
  struct ADJOINT_TYPE : public descriptors::TYPED_FIELD_BASE<AdjointType,1> { };
} // namespace olb::parameters

/// @brief Step 3: Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector origin{0, 0, 0};
  IndicatorCuboid3D<T> cuboid(extent, origin);

  const T physDeltaX = extent[0] / parameters.get<parameters::RESOLUTION>();
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

/// @brief Step 5: Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
void prepareGeometry(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;
  auto& geometry = myCase.getGeometry();

  /// @li Set material numbers
  geometry.rename(0, 2);
  geometry.rename(2, 1, {1,1,1});
  geometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

/// @brief Step 6: Set lattice dynamics and boundary conditions
/// @param myCase The Case instance which keeps the simulation data
void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& parameters = myCase.getParameters();

  /// @li Set up a unit converter with the characteristic physical units
  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<MyCase::value_t,
                                                                         MyCase::descriptor_t_of<NavierStokes>>>(
    parameters.get<parameters::RESOLUTION>(),               // resolution
    parameters.get<parameters::LATTICE_RELAXATION_TIME>(),  // latticeRelaxationTime
    parameters.get<parameters::DOMAIN_EXTENT>()[0],         // charPhysLength: reference length of simulation geometry in [m]
    parameters.get<parameters::PHYS_CHAR_VELOCITY>(),       // charPhysVelocity: highest expected velocity during simulation in [m/s]
    parameters.get<parameters::PHYS_CHAR_VISCOSITY>(),      // physViscosity: physical kinematic viscosity in [m^2/s]
    parameters.get<parameters::PHYS_CHAR_DENSITY>()         // physDensity: physical density [kg/m^3]
  );
  auto& converter = lattice.getUnitConverter();
  converter.print();

  /// @li Material=1 --> bulk dynamics
  dynamics::set<ForcedBGKdynamics>(lattice, geometry.getMaterialIndicator({1}));
  /// @li Material=2,3 --> velocity boundary
  boundary::set<boundary::LocalVelocity>(lattice, geometry.getMaterialIndicator({2}));
  /// @li Set lattice relaxation frequency
  lattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
}

// Step 7: Set initial values for primal variables (e.g. velocity, density) and additional fields
/// @param myCase The Case instance which keeps the simulation data
/// @note Initial values have to be set using lattice units
void setInitialValues(MyCase& myCase)
{
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& converter = lattice.getUnitConverter();

  /// @li Set force field
  ForceTestFlow3D<T,T,MyCase::descriptor_t> forceF(converter);
  const T latticeScaling(converter.getConversionFactorMass() / converter.getConversionFactorForce());
  AnalyticalScaled3D<T,T> scaledForceF(forceF, latticeScaling);
  fields::set<descriptors::FORCE>(lattice, geometry.getMaterialIndicator({1}), scaledForceF);
  lattice.initialize();
}

/// Step 8.1: Update boundary values at times (and additional fields, if needed)
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
/// In order to start the simulation slowly, we start with zero boundary velocity and gradually increase it
/// @note Boundary values have to be set using lattice units
void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& converter = lattice.getUnitConverter();
  const T physStartT = myCase.getParameters().get<parameters::PHYS_MAX_T>() * (2.0 / 3.0);

  const std::size_t itStart = converter.getLatticeTime(physStartT);
  if (iT <= itStart) {
    /// @li Compute scaling factor
    PolynomialStartScale<T,T> startScaleF(itStart, T(1));
    const T iTvec[1] = {T(iT)};
    T frac[1] = {};
    startScaleF(frac, iTvec);

    /// @li Take analytical velocity solution, scale it to lattice units, set the boundary data
    VelocityTestFlow3D<T,T,MyCase::descriptor_t> velocityF(converter);
    AnalyticalScaled3D<T,T> uBoundaryStartF(velocityF, frac[0]);
    momenta::setVelocity(lattice, geometry.getMaterialIndicator({2}), uBoundaryStartF);
    lattice.template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }
}

/// Step 8.3: Compute simulation results at times
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();
  const T physMaxT = myCase.getParameters().get<parameters::PHYS_MAX_T>();
  const std::size_t iTlog = converter.getLatticeTime(physMaxT/5);

  /// @li Print some (numerical and computational) statistics
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
  const T physMaxT = myCase.getParameters().get<parameters::PHYS_MAX_T>();
  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(physMaxT);
  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT)
  {
    /// @li Step 8.1: Update the Boundary Values and Fields at Times
    if (!myCase.hasName(Adjoint{})) {
      setTemporalValues(myCase, iT);
    }
    /// @li Step 8.2: Collide and Stream Execution
    myCase.getLattice(NavierStokes{}).collideAndStream();
    /// @li Stripe off density offset due to Dirichlet boundary conditions
    myCase.getLattice(NavierStokes{}).stripeOffDensityOffset(
      myCase.getLattice(NavierStokes{}).getStatistics().getAverageRho() - T{1});
    /// @li Step 8.3: Computation and Output of the Results
    getResults(myCase, timer, iT);
  }
  timer.stop();
  timer.printSummary();
  myCase.getLattice(NavierStokes{}).setProcessingContext(ProcessingContext::Evaluation);
}

void prepareAdjointLattice(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  auto& adjointLattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& geometry = optiCase.getCase(Adjoint{}).getGeometry();
  auto& parameters = optiCase.getCase(Adjoint{}).getParameters();

  // Set adjoint unit converter as the controlled case
  adjointLattice.setUnitConverter(controlledLattice.getUnitConverter());

  // Define dual physics
  switch (parameters.template get<parameters::ADJOINT_TYPE>()) {
    case AdjointType::discrete:
      dynamics::set<olb::Dual<BulkDynamics>>(adjointLattice, geometry.getMaterialIndicator({1}));
      break;
    case AdjointType::continuous:
      dynamics::set<DualForcedBGKDynamics<T,MyCase::descriptor_t>>(adjointLattice, geometry.getMaterialIndicator({1}));
      break;
    default:
      throw std::runtime_error("Unknown adjoint type");
  }
  boundary::set<boundary::BounceBack>(adjointLattice, geometry.getMaterialIndicator({2}));
  adjointLattice.template setParameter<descriptors::OMEGA>(adjointLattice.getUnitConverter().getLatticeRelaxationFrequency());
  // Revert the direction of streaming for adjoint backwards propagation
  adjointLattice.revertStreaming();
}

void setAdjointInitialValues(MyOptiCase& optiCase) {
  auto& adjointLattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& referenceLattice = optiCase.getCase(Reference{}).getLattice(NavierStokes{});
  auto& converter = adjointLattice.getUnitConverter();
  auto& geometry = optiCase.getCase(Adjoint{}).getGeometry();
  auto& parameters = optiCase.getCase(Adjoint{}).getParameters();
  auto bulkIndicator = geometry.getMaterialIndicator({1});

  // Initialize dual problem
  // This needs to be before copying the fields as otherwise the the copied fields will be overwritten
  adjointLattice.initialize();

  // Compute source term for the dual simulation
  auto velocityO = makeWriteFunctorO<functors::VelocityF,descriptors::VELOCITY>(referenceLattice);
  velocityO->restrictTo(bulkIndicator);
  velocityO->template setParameter<parameters::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
  velocityO->apply();

  auto objectiveDerivativeO = makeWriteFunctorO<functors::DerivativeF<ObjectiveF,descriptors::POPULATION,BulkDynamics>,
                                                opti::DJDF>(controlledLattice);
  objectiveDerivativeO->restrictTo(bulkIndicator);
  objectiveDerivativeO->template setParameter<parameters::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
  objectiveDerivativeO->template setParameter<parameters::NORMALIZE>(computeL2Norm<descriptors::VELOCITY>(referenceLattice, bulkIndicator, converter.getPhysDeltaX()));
  objectiveDerivativeO->template setParameter<parameters::FACTOR>(util::pow(converter.getPhysDeltaX(), 3));
  objectiveDerivativeO->apply();

  // Write primal population in dedicated fields
  switch (parameters.template get<parameters::ADJOINT_TYPE>()) {
    case AdjointType::discrete:
      copyFields<descriptors::POPULATION,opti::F>(controlledLattice, adjointLattice);
      break;
    case AdjointType::continuous:
      writeFunctorTo<functors::PopulationF,opti::F>(controlledLattice, bulkIndicator);
      copyFields<opti::F,opti::F>(controlledLattice, adjointLattice);
      break;
    default:
      throw std::runtime_error("Unknown adjoint type");
  }

  // Provide fields required by the dual collision operator
  copyFields<ControlledField,ControlledField>(controlledLattice, adjointLattice);
  copyFields<opti::DJDF,opti::DJDF>(controlledLattice, adjointLattice);
}

void setInitialControl(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& lattice = controlledCase.getLattice(NavierStokes{});
  auto& geometry = controlledCase.getGeometry();
  auto& control = optiCase.getController();

  // Intialize controlled field in superLattice
  fields::set<ControlledField>(lattice, geometry.getMaterialIndicator(1), Vector<T,3>{0.0, 0.0, 0.0});
  control.template setProjection<projection::ForceFactor<T>>(lattice.getUnitConverter());
  control.template set<ControlledField>(geometry, 1, lattice);
}

void applyControl(MyOptiCase& optiCase) {
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& lattice = controlledCase.getLattice(NavierStokes{});

  optiCase.getController().template setUpdatedControlsOnField<ControlledField>(lattice);
  lattice.template setProcessingContext<Array<ControlledField>>(ProcessingContext::Simulation);
}

MyCase::value_t objectiveF(MyOptiCase& optiCase) {
  auto& referenceLattice = optiCase.getCase(Reference{}).getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& converter = controlledLattice.getUnitConverter();
  auto& objectiveDomain = optiCase.getController().getDesignDomain<MyCase::descriptor_t>();

  // Evaluate functor for objective computation
  auto objectiveO = makeWriteFunctorO<ObjectiveF,opti::J>(controlledLattice);
  objectiveO->restrictTo(objectiveDomain);

  // Write velocity field in the reference case
  auto velocityO = makeWriteFunctorO<functors::VelocityF,ObjectiveF::Reference>(referenceLattice);
  velocityO->template setParameter<parameters::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
  velocityO->restrictTo(objectiveDomain);
  velocityO->apply();

  // Get solution from the reference simulation for the inverse problem
  copyFields<ObjectiveF::Reference,ObjectiveF::Reference>(referenceLattice, controlledLattice);
  objectiveO->template setParameter<parameters::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
  objectiveO->template setParameter<parameters::NORMALIZE>(computeL2Norm<ObjectiveF::Reference>(referenceLattice, objectiveDomain, converter.getPhysDeltaX()));
  objectiveO->apply();
  return integrateField<opti::J>(controlledLattice, objectiveDomain, converter.getPhysDeltaX())[0];
}

std::vector<MyCase::value_t> derivativeF(MyOptiCase& optiCase) {
  auto& adjointLattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& converter = controlledLattice.getUnitConverter();

  // Evaluate optimality condition
  auto optimalityO = makeWriteFunctorO<functors::OptimalityF<BulkDynamics,ControlledField>,
                                       opti::SENSITIVITY<ControlledField>>(adjointLattice);
  optiCase.getController().template setUpdatedProjectionDerivativesOnField<opti::DPROJECTIONDALPHA<ControlledField>>(adjointLattice);
  adjointLattice.template setProcessingContext<Array<opti::DPROJECTIONDALPHA<ControlledField>>>(ProcessingContext::Simulation);
  // Compute jacobian of collision operator regarding control variable
  auto dCDalphaO = makeWriteFunctorO<functors::DerivativeF<functors::CollisionF<BulkDynamics>,ControlledField,BulkDynamics>,
                                     opti::DCDALPHA<ControlledField>>(controlledLattice);
  dCDalphaO->template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  dCDalphaO->template setParameter<parameters::FACTOR>(1.0);
  dCDalphaO->apply();
  // Jacobian is computed on primal lattice as jacobian is evaluated for primal populations
  copyFields<opti::DCDALPHA<ControlledField>,opti::DCDALPHA<ControlledField>>(controlledLattice, adjointLattice);
  optimalityO->apply();

  // Return serial vector containing total derivatives of objective regarding controls
  adjointLattice.setProcessingContext(ProcessingContext::Evaluation);
  return getSerializedFromField<opti::SENSITIVITY<ControlledField>>(adjointLattice, optiCase.getController().getDesignDomain<MyCase::descriptor_t>());
}

void getOptiResults(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& converter = controlledLattice.getUnitConverter();
  std::size_t iT = optiCase.getOptimizationStep();

  SuperVTMwriter3D<T> vtmWriter("testFlowOpti3d");
  SuperLatticePhysVelocity3D velocityF(controlledLattice, converter);
  SuperLatticePhysPressure3D pressureF(controlledLattice, converter);
  vtmWriter.addFunctor(velocityF);
  vtmWriter.addFunctor(pressureF);

  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  // Writes the VTK files
  if (iT > 0) {
    controlledLattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }
}


MyCase::value_t computeObjective(MyOptiCase& optiCase) {
  auto& controlledCase = optiCase.getCase(Controlled{});
  // Reset prior simulation lattice
  controlledCase.resetLattices();

  // Prepare new lattices
  prepareLattice(controlledCase);
  setInitialValues(controlledCase);

  // Set updated controls in the simulation
  applyControl(optiCase);

  // Execute simulation
  simulate(controlledCase);

  // Optimization related IO
  getOptiResults(optiCase);

  // Compute Objective value from simulation results
  return objectiveF(optiCase);
}

std::vector<MyCase::value_t> computeDerivative(MyOptiCase& optiCase) {
  auto& adjointCase = optiCase.getCase(Adjoint{});
  // Reset prior simulation lattice
  adjointCase.resetLattices();

  // Prepare new lattices
  prepareAdjointLattice(optiCase);
  setAdjointInitialValues(optiCase);

  // Execute simulation
  simulate(adjointCase);

  // Compute Derivative value from simulation results
  return derivativeF(optiCase);
}
