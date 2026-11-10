/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Shota Ito
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

// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::opti;

using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q19<>>
>;
using MyOptiCase = OptiCaseAdjoint<
  Controlled, MyCase,
  Adjoint, MyCase,
  Reference, MyCase
>;

using BulkDynamics = PorousBGKdynamics<MyCase::value_t,
                                       MyCase::descriptor_t_of<NavierStokes>>;
using ObjectiveF = functors::L2DistanceF<functors::VelocityF>;
using ControlledField = descriptors::POROSITY;

// Switch between different combination of boundary conditions in primal/adjoint simulation
enum class BoundaryType : int {
  PRIMAL_LOCAL_ADJOINT_BOUNCE_BACK = 0,
  PRIMAL_LOCAL_ADJOINT_ZERO = 1,
  PRIMAL_LOCAL_ADJOINT_LOCAL = 2, // adjoint pressure BC leads to infinite increase in sensitivity
  PRIMAL_ZOUHE_ADJOINT_ZOUHE = 3  // same here
};

namespace olb::parameters{

struct DESIGN_DOMAIN_EXTENT : public descriptors::FIELD_BASE<0,1>{ };
struct DESIGN_DOMAIN_ORIGIN : public descriptors::FIELD_BASE<0,1>{ };
struct REFERENCE_OBJECT_EXTENT : public descriptors::FIELD_BASE<0,1>{ };
struct REFERENCE_OBJECT_ORIGIN : public descriptors::FIELD_BASE<0,1>{ };
struct BOUNDARY_TYPE : public descriptors::TYPED_FIELD_BASE<BoundaryType,1>{ };

}

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

void prepareGeometry(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector originDesign = parameters.get<parameters::DESIGN_DOMAIN_ORIGIN>();
  const Vector extentDesign = parameters.get<parameters::DESIGN_DOMAIN_EXTENT>();
  const T physDeltaX = extent[0] / parameters.get<parameters::RESOLUTION>();
  auto& superGeometry = myCase.getGeometry();

  superGeometry.rename(0, 2);
  superGeometry.rename(2, 1, {1,1,1});

  Vector<T,3> origin{0.,0.,0.};
  origin[0] += physDeltaX/2.;
  origin[1] -= physDeltaX/2.;
  origin[2] += physDeltaX/2.;

  Vector<T,3> extentInflow{extent[0], 0., extent[2]};
  extentInflow[0] -= 2.*physDeltaX/2.;
  extentInflow[1] +=    physDeltaX/2.;
  extentInflow[2] -= 2.*physDeltaX/2.;

  IndicatorCuboid3D<T> inflow(extentInflow, origin);
  superGeometry.rename(2, 3, inflow);

  origin = {0., extent[1], 0.};
  origin[0] += physDeltaX/2.;
  origin[1] -= physDeltaX/2.;
  origin[2] += physDeltaX/2.;

  Vector<T,3> extentOutflow{extent[0], extent[1], extent[2]};
  extentOutflow[0] -= 2.*physDeltaX/2.;
  extentOutflow[1] +=    physDeltaX/2.;
  extentOutflow[2] -= 2.*physDeltaX/2.;

  IndicatorCuboid3D<T> outflow(extentOutflow, origin);
  superGeometry.rename(2, 4, outflow);

  // Indicators for material 6 (designDomain)
  IndicatorCuboid3D<T> designDomain(extentDesign, originDesign);
  superGeometry.rename(1, 6, designDomain);

  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
}

void prepareLattice(MyCase& myCase) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  sLattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    parameters.get<parameters::RESOLUTION>(),               // resolution
    parameters.get<parameters::LATTICE_RELAXATION_TIME>(),  // relaxationTime
    parameters.get<parameters::DOMAIN_EXTENT>()[0],         // charPhysLength: reference length of simulation geometry in [m]
    parameters.get<parameters::PHYS_CHAR_VELOCITY>(),       // charPhysVelocity: highest expected velocity during simulation in [m/s]
    parameters.get<parameters::PHYS_CHAR_VISCOSITY>(),      // physViscosity: physical kinematic viscosity in [m^2/s]
    parameters.get<parameters::PHYS_CHAR_DENSITY>()         // physDensity: physical density [kg/m^3]
  );
  auto& converter = sLattice.getUnitConverter();
  converter.print();

  auto bulkIndicator = superGeometry.getMaterialIndicator({1,3,4,6});
  dynamics::set<PorousBGKdynamics>(sLattice, bulkIndicator);

  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);
  switch (parameters.get<parameters::BOUNDARY_TYPE>()) {
    case BoundaryType::PRIMAL_ZOUHE_ADJOINT_ZOUHE: {
      boundary::set<boundary::ZouHeVelocity>(sLattice, superGeometry.getMaterialIndicator({3}));
      boundary::set<boundary::ZouHePressure>(sLattice, superGeometry.getMaterialIndicator({4}));
      break;
    }
    default: {
      boundary::set<boundary::LocalVelocity>(sLattice, superGeometry.getMaterialIndicator({3}));
      boundary::set<boundary::LocalPressure>(sLattice, superGeometry.getMaterialIndicator({4}));
      break;
    }
  }
  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  const Vector originRefObj = parameters.get<parameters::REFERENCE_OBJECT_ORIGIN>();
  const Vector extentRefObj = parameters.get<parameters::REFERENCE_OBJECT_EXTENT>();
  auto bulkIndicator = superGeometry.getMaterialIndicator({1,2,3,4,6});

  AnalyticalConst3D<T,T> one(1.);
  AnalyticalConst3D<T,T> zero(0.);

  IndicatorCuboid3D<T> referenceObject(extentRefObj, originRefObj);
  SuperIndicatorFfromIndicatorF<T,DESCRIPTOR::d> indicatorF(referenceObject, superGeometry);
  fields::set<descriptors::POROSITY>(sLattice, bulkIndicator, one);
  fields::set<descriptors::POROSITY>(sLattice, indicatorF, zero);

  sLattice.initialize();
}

void setTemporalValues(MyCase& myCase, std::size_t iT) {
  using T = MyCase::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& converter = sLattice.getUnitConverter();
  auto& superGeometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  const T physStartT = parameters.get<parameters::PHYS_START_T>();
  const T physBoundaryT = parameters.get<parameters::PHYS_BOUNDARY_VALUE_UPDATE_T>();

  const std::size_t itStart = converter.getLatticeTime(physStartT);
  const std::size_t itUpdate = converter.getLatticeTime(physBoundaryT);
  if (iT <= itStart && iT % itUpdate == 0) {
    PolynomialStartScale<T,std::size_t> StartScale( itStart, T( 1 ) );
    std::size_t iTvec[1] = {iT};
    T frac[1] = {};
    StartScale( frac,iTvec );

    AnalyticalConst3D<T,T> uF(0., frac[0] * converter.getCharPhysVelocity(), 0.);
    momenta::setVelocity(sLattice, superGeometry.getMaterialIndicator({3}), uF);
    sLattice.template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }
}

void getResults(MyCase& myCase, std::size_t iT, util::Timer<MyCase::value_t> timer) {
  using T = MyCase::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& parameters = myCase.getParameters();
  const T physMaxT = parameters.get<parameters::PHYS_MAX_T>();
  auto& converter = sLattice.getUnitConverter();
  const std::size_t iTlog = converter.getLatticeTime(physMaxT/5);

  // Get statistics
  if (iT%iTlog == 0) {
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    timer.print(iT);
  }
}

void simulate(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  auto& converter = sLattice.getUnitConverter();
  const T physMaxT = parameters.get<parameters::PHYS_MAX_T>();

  // === 4th Step: Main Loop with Timer ===
  const std::size_t iTmax = converter.getLatticeTime(physMaxT);
  util::Timer<T> timer(iTmax, superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    if (!myCase.hasName(Adjoint{})) {
      setTemporalValues(myCase, iT);
    }
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    if (myCase.hasName(Adjoint{})) {
      myCase.getLattice(NavierStokes{}).stripeOffDensityOffset(
        myCase.getLattice(NavierStokes{}).getStatistics().getAverageRho() - T{1});
    }
    // === 7th Step: Computation and Output of the Results ===
    getResults(myCase, iT, timer);
  }
  timer.stop();
  timer.printSummary();
  myCase.getLattice(NavierStokes{}).setProcessingContext(ProcessingContext::Evaluation);
}

void prepareAdjointLattice(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& adjointLattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& geometry = optiCase.getCase(Adjoint{}).getGeometry();
  auto& parameters = optiCase.getCase(Adjoint{}).getParameters();

  // Set adjoint unit converter as the controlled case
  adjointLattice.setUnitConverter(controlledLattice.getUnitConverter());

  // Define dual physics
  dynamics::set<olb::Dual<BulkDynamics>>(adjointLattice, geometry.getMaterialIndicator({1,6}));
  boundary::set<boundary::BounceBack>(adjointLattice, geometry.getMaterialIndicator({2}));

  switch (parameters.get<parameters::BOUNDARY_TYPE>()) {
    case BoundaryType::PRIMAL_LOCAL_ADJOINT_BOUNCE_BACK: {
      boundary::set<boundary::BounceBack>(adjointLattice, geometry.getMaterialIndicator({3,4}));
      break;
    }
    case BoundaryType::PRIMAL_LOCAL_ADJOINT_ZERO: {
      boundary::set<boundary::ZeroPopulation>(adjointLattice, geometry.getMaterialIndicator({3,4}));
      break;
    }
    case BoundaryType::PRIMAL_LOCAL_ADJOINT_LOCAL: {
      dynamics::set<DualLocalVelocityBoundary<T,DESCRIPTOR,1,-1>>(adjointLattice, geometry.getMaterialIndicator({3}));
      dynamics::set<DualLocalPressureBoundary<T,DESCRIPTOR,1,1>>(adjointLattice, geometry.getMaterialIndicator({4}));
      break;
    }
    case BoundaryType::PRIMAL_ZOUHE_ADJOINT_ZOUHE: {
      dynamics::set<DualZouHeVelocityBoundary<T,DESCRIPTOR,1,-1>>(adjointLattice, geometry.getMaterialIndicator({3}));
      dynamics::set<DualZouHePressureBoundary<T,DESCRIPTOR,1,1>>(adjointLattice, geometry.getMaterialIndicator({4}));
      break;
    }
  }
  adjointLattice.template setParameter<descriptors::OMEGA>(adjointLattice.getUnitConverter().getLatticeRelaxationFrequency());

  adjointLattice.revertStreaming();
}

void setAdjointInitialValues(MyOptiCase& optiCase) {
  auto& adjointLattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& referenceLattice = optiCase.getCase(Reference{}).getLattice(NavierStokes{});
  auto& converter = adjointLattice.getUnitConverter();
  auto& geometry = optiCase.getCase(Adjoint{}).getGeometry();
  auto bulkIndicator = geometry.getMaterialIndicator({1,2,3,4,6});
  auto objectiveDomain = geometry.getMaterialIndicator({1,6});

  adjointLattice.initialize();

  // Initialize dual problem
  auto velocityO = makeWriteFunctorO<functors::VelocityF,descriptors::VELOCITY>(referenceLattice);
  velocityO->template setParameter<parameters::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
  velocityO->restrictTo(objectiveDomain);
  velocityO->apply();
  // Compute source term for the dual simulation
  auto objectiveDerivativeO = makeWriteFunctorO<functors::DerivativeF<ObjectiveF,descriptors::POPULATION,BulkDynamics>,
                                                opti::DJDF>(controlledLattice);
  objectiveDerivativeO->restrictTo(objectiveDomain);
  objectiveDerivativeO->template setParameter<parameters::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
  objectiveDerivativeO->template setParameter<parameters::NORMALIZE>(computeL2Norm<descriptors::VELOCITY>(referenceLattice, objectiveDomain, converter.getPhysDeltaX()));
  objectiveDerivativeO->template setParameter<parameters::FACTOR>(util::pow(converter.getPhysDeltaX(),3));
  objectiveDerivativeO->apply();

  // Write primal population in dedicated fields
  writeFunctorTo<functors::PopulationF,opti::F>(controlledLattice, bulkIndicator);

  // Provide fields required by the dual collision operator
  copyFields<ControlledField,ControlledField>(controlledLattice, adjointLattice);
  copyFields<opti::F,opti::F>(controlledLattice, adjointLattice);
  copyFields<opti::DJDF,opti::DJDF>(controlledLattice, adjointLattice);
}

void setInitialControl(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& lattice = controlledCase.getLattice(NavierStokes{});
  auto& geometry = controlledCase.getGeometry();
  auto& control = optiCase.getController();
  auto& converter = lattice.getUnitConverter();
  auto& parameters = controlledCase.getParameters();

  // Intialize controlled field in superLattice
  T porosity = projection::permeabilityToPorosity(parameters.get<parameters::INITIAL_CONTROL_SCALAR>(), converter);
  std::cout << "POROSITY: " << porosity << std::endl;
  AnalyticalConst3D<T,T> controls(porosity);
  fields::set<ControlledField>(lattice, geometry.getMaterialIndicator({6}), controls);
  control.template setProjection<projection::Sigmoid<T>>();
  control.template set<ControlledField>(geometry, 6, lattice);
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
  auto& geometry = optiCase.getCase(Controlled{}).getGeometry();
  auto objectiveDomain = geometry.getMaterialIndicator({1,6});

  // Evaluate functor for objective computation
  auto objectiveO = makeWriteFunctorO<ObjectiveF,opti::J>(controlledLattice);
  objectiveO->restrictTo(objectiveDomain);

  // Write velocity field in the reference case
  auto velocityO = makeWriteFunctorO<functors::VelocityF,ObjectiveF::Reference>(referenceLattice);
  velocityO->restrictTo(objectiveDomain);
  velocityO->template setParameter<parameters::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
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

  SuperVTMwriter3D<T> vtmWriter("domainId3d");
  SuperLatticePhysVelocity3D velocityF(controlledLattice, converter);
  SuperLatticePhysPressure3D pressureF(controlledLattice, converter);
  SuperLatticeField3D<T,MyCase::descriptor_t,ControlledField> porosityF(controlledLattice);
  vtmWriter.addFunctor(velocityF);
  vtmWriter.addFunctor(pressureF);
  vtmWriter.addFunctor(porosityF);

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

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);

  MyCase::ParametersD myCaseParameters;
  {
    using namespace parameters;
    myCaseParameters.set<DOMAIN_EXTENT               >({1.0, 1.0, 1.0});
    myCaseParameters.set<PHYS_CHAR_VELOCITY          >(            1.0);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY         >(            0.1);
    myCaseParameters.set<PHYS_CHAR_DENSITY           >(            1.0);
    myCaseParameters.set<PHYS_MAX_T                  >(            6.0);
    myCaseParameters.set<RESOLUTION                  >(             20);
    myCaseParameters.set<LATTICE_RELAXATION_TIME     >(            0.8);
    myCaseParameters.set<PHYS_START_T                >(            1.0);
    myCaseParameters.set<PHYS_BOUNDARY_VALUE_UPDATE_T>(           0.02);
    myCaseParameters.set<DESIGN_DOMAIN_EXTENT        >({0.4, 0.4, 0.4});
    myCaseParameters.set<DESIGN_DOMAIN_ORIGIN        >({0.3, 0.3, 0.3});
    myCaseParameters.set<REFERENCE_OBJECT_EXTENT     >({0.2, 0.2, 0.2});
    myCaseParameters.set<REFERENCE_OBJECT_ORIGIN     >({0.4, 0.4, 0.4});
    myCaseParameters.set<INITIAL_CONTROL_SCALAR      >(         1.e-2);
    myCaseParameters.set<BOUNDARY_TYPE               >(BoundaryType::PRIMAL_LOCAL_ADJOINT_BOUNCE_BACK);
  }

  OptimizerParametersD<MyCase::value_t> optiParameters;
  {
    using namespace parameters;
    optiParameters.set<optimizer::EPS>(5e-6);
    optiParameters.set<optimizer::MAX_OPTIMIZATION_STEPS>(10);
    optiParameters.set<linesearch::MAX_STEP_ATTEMPTS>(20);
    optiParameters.set<linesearch::STEP_CONDITION>("Wolfe");
  }

  /// Step 3: Create Mesh
  Mesh mesh = createMesh(myCaseParameters);

  /// Step B: Create Cases and prepare
  // Prepare controlled case
  MyCase myCase(myCaseParameters, mesh);
  prepareGeometry(myCase);
  prepareLattice(myCase);

  // Create adjoint case for gradient computation
  MyCase adjointCase(myCaseParameters, mesh);
  prepareGeometry(adjointCase);

  // Compute solution for the objective functional
  MyCase referenceCase(myCaseParameters, mesh);
  prepareGeometry(referenceCase);
  prepareLattice(referenceCase);
  setInitialValues(referenceCase);
  simulate(referenceCase);

  /// Step C: Create OptiCase and set cases
  MyOptiCase optiCase;
  optiCase.setCase<Controlled>(myCase);
  optiCase.setCase<Adjoint>(adjointCase);
  optiCase.setCase<Reference>(referenceCase);

  /// Step D: Set initial control
  setInitialControl(optiCase);

  /// Step E: Define objective routine
  optiCase.setObjective(computeObjective);

  /// Step F: Define derivative routine
  optiCase.setDerivative(computeDerivative);

  OptimizerLBFGS<MyCase::value_t,std::vector<MyCase::value_t>> optimizer(
    optiCase.getController().size(), optiParameters
  );
  optimizer.optimize(optiCase);
}
