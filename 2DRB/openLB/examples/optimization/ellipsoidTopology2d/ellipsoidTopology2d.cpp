/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Jiaqi Lan, Shota Ito
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

namespace olb::names {
// Optimality system with dissipation objective
struct DissipationOS {};
// Optimality system with pressuredrop objective
struct PressureDropOS {};
}

// Define the structure of the case
using MyCase = Case< // olb\src\case\case.h
  NavierStokes, Lattice<double,descriptors::D2Q9<descriptors::POROSITY,descriptors::MATERIAL>>
>;

// Define the optimization problem to specify the optimality system
using DissipationEllipsoid = SteadyOptimalitySystem<
  functors::AddF<functors::AddF<functors::DissipationF,functors::PorousDissipationF>, // Objective functional
  functors::TikhonovRegularizationF<descriptors::POROSITY>>,
  PorousBGKdynamics<MyCase::value_t,MyCase::descriptor_t>,                            // Primal dynamics
  descriptors::POROSITY                                                               // Controlled field
>;

using PressureDropEllipsoid = SteadyOptimalitySystem<
  functors::AddF<functors::PressureDropF<7,8>,                                        // Objective functional
  functors::TikhonovRegularizationF<descriptors::POROSITY>>,
  PorousBGKdynamics<MyCase::value_t,MyCase::descriptor_t>,                            // Primal dynamics
  descriptors::POROSITY                                                               // Controlled field
>;

using MyOptiCase = OptiCaseAdjoint<
  Controlled, MyCase,
  Adjoint, MyCase,
  DissipationOS, DissipationEllipsoid,
  PressureDropOS, PressureDropEllipsoid
>;

enum class ObjectiveType : int {
  DISSIPATION = 0,
  PRESSUREDROP = 1,
};

namespace olb::parameters {
  struct OBJECTIVE_TYPE: public descriptors::TYPED_FIELD_BASE<ObjectiveType,1> {};
}

auto createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const T L = parameters.get<parameters::PHYS_CHAR_LENGTH>() / parameters.get<parameters::RESOLUTION>();
  const T lengthX = 3*parameters.get<parameters::PHYS_CHAR_LENGTH>();
  const T lengthY = parameters.get<parameters::PHYS_CHAR_LENGTH>();
  const Vector<T,2> origin(0., 0.);
  const Vector<T,2> extend(lengthX, lengthY);
  IndicatorCuboid2D<T> channel(extend, origin);
  Mesh<T,MyCase::d> mesh{channel, L, singleton::mpi().getSize()};
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({false, true});
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& superGeometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  const T L = parameters.get<parameters::PHYS_CHAR_LENGTH>() / parameters.get<parameters::RESOLUTION>();
  const T lengthX = 3*parameters.get<parameters::PHYS_CHAR_LENGTH>();
  const T lengthY = parameters.get<parameters::PHYS_CHAR_LENGTH>();

  // Flow domain
  Vector<T,2> origin(0., 0.);
  Vector<T,2> extend(lengthX, lengthY);
  IndicatorCuboid2D<T> channel(extend, origin);

  // Edge point
  auto min = channel.getMin();
  auto max = channel.getMax();

  // Inflow regions
  origin = Vector<T,2>(min[0] - 0.5*L, min[1] - 0.5*L);
  extend = Vector<T,2>(L, max[1] + L);
  IndicatorCuboid2D<T> inflow(extend, origin);
  origin[0] = min[0] + 1.5*L;
  IndicatorCuboid2D<T> objectiveDomainInflow(extend, origin);

  // Outflow regions
  origin = Vector<T,2>(max[0] - 0.5*L, min[1] - 0.5*L);
  extend = Vector<T,2>(L, max[1] + L);
  IndicatorCuboid2D<T> outflow(extend, origin);
  origin[0] = max[0] - 2.5*L;
  IndicatorCuboid2D<T> objectiveDomainOutflow(extend, origin);

  // Design regions
  origin = Vector<T,2>(min[0] + max[1] - 0.5*L, min[1] + 0.25*max[1] - 0.5*L);
  extend = Vector<T,2>(max[1] + L, 0.5*max[1] + L);
  IndicatorCuboid2D<T> designDomain(extend, origin);

  // Material numbers
  // 1: fluid, 3: inflow, 4: outflow, 6: designDomain
  superGeometry.rename(0, 1, channel);
  superGeometry.rename(1, 3, inflow);
  superGeometry.rename(1 ,4, outflow);
  superGeometry.clean();
  superGeometry.rename(1 ,6, designDomain);
  superGeometry.rename(1 ,7, objectiveDomainInflow);
  superGeometry.rename(1 ,8, objectiveDomainOutflow);
}

void prepareLattice(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& superGeometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& parameters = myCase.getParameters();

  // define UnitConverter
  sLattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,MyCase::descriptor_t>>(
    (int) parameters.get<parameters::RESOLUTION>(),
    ( T ) parameters.get<parameters::LATTICE_RELAXATION_TIME>(),
    ( T ) parameters.get<parameters::PHYS_CHAR_LENGTH>(),
    ( T ) parameters.get<parameters::PHYS_CHAR_VELOCITY>(),
    ( T ) parameters.get<parameters::PHYS_CHAR_VISCOSITY>(),
    ( T ) parameters.get<parameters::PHYS_CHAR_DENSITY>()
  );
  auto& converter = sLattice.getUnitConverter();

  // define dynamics and bc
  auto bulkIndicator = superGeometry.getMaterialIndicator({1,3,4,6,7,8});
  dynamics::set<DissipationEllipsoid::dynamics_t>(sLattice, bulkIndicator);
  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);
  // set parameters
  sLattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& superGeometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes{});

  // Set porosity field
  auto bulkIndicator = superGeometry.getMaterialIndicator({1,3,4,6,7,8});
  auto solidIndicator = superGeometry.getMaterialIndicator({0});
  AnalyticalConst2D<T,T> zero(0.);
  AnalyticalConst2D<T,T> one(1.);
  AnalyticalConst2D<T,T> mat6(6);
  AnalyticalConst2D<T,T> mat7(7);
  AnalyticalConst2D<T,T> mat8(8);

  fields::set<descriptors::POROSITY>(sLattice, bulkIndicator, one);
  fields::set<descriptors::POROSITY>(sLattice, solidIndicator, zero);
  fields::set<descriptors::MATERIAL>(sLattice, superGeometry.getMaterialIndicator({6}), mat6);
  fields::set<descriptors::MATERIAL>(sLattice, superGeometry.getMaterialIndicator({7}), mat7);
  fields::set<descriptors::MATERIAL>(sLattice, superGeometry.getMaterialIndicator({8}), mat8);

  sLattice.initialize();
}

void setTemporalValues(MyCase& myCase, std::size_t iT) {
  using T = MyCase::value_t;
  auto& superGeometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& converter = sLattice.getUnitConverter();
  auto& parameters = myCase.getParameters();
  const T physStartUpT = parameters.get<parameters::PHYS_MAX_T>() / 4;
  const T physBoundaryValueUpdateT = parameters.get<parameters::PHYS_BOUNDARY_VALUE_UPDATE_T>();

  std::size_t iTmaxStart = converter.getLatticeTime(physStartUpT);
  std::size_t iTupdate = converter.getLatticeTime(physBoundaryValueUpdateT);
  if ( iT%iTupdate==0 && iT<= iTmaxStart ) {
    PolynomialStartScale<T,T> StartScale( iTmaxStart, T( 1 ) );
    T iTvec[1] = {T( iT )}; T frac[1] = {};
    StartScale( frac,iTvec );
    Vector<T,2> u_free {frac[0]*converter.getCharPhysVelocity(), 0.0};
    AnalyticalConst2D<T,T> u(u_free);
    momenta::setVelocity(sLattice, superGeometry.getMaterialIndicator({3}), u);
    sLattice.template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

void getResults(MyCase& myCase, util::Timer<MyCase::value_t>& timer, std::size_t iT) {
  using T = MyCase::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& converter = sLattice.getUnitConverter();
  auto& parameters = myCase.getParameters();
  const T physMaxT = parameters.get<parameters::PHYS_MAX_T>();

  const std::size_t iTlog = converter.getLatticeTime(physMaxT/5.);

  if (iT%iTlog == 0) {
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    timer.print(iT);
  }
}

void simulate(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& superGeometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& converter = sLattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  const std::size_t iTmax = converter.getLatticeTime(parameters.get<parameters::PHYS_MAX_T>());
  util::Timer<T> timer(iTmax, superGeometry.getStatistics().getNvoxel());
  timer.start();

  // main loop
  for (std::size_t iT=0; iT < iTmax; ++iT)
  {
    if (myCase.hasName(Controlled{})) {
      setTemporalValues(myCase, iT);
    }
    myCase.getLattice(NavierStokes{}).collideAndStream();
    getResults(myCase, timer, iT);
  }

  timer.stop();
  timer.printSummary();
  myCase.getLattice(NavierStokes{}).setProcessingContext(ProcessingContext::Evaluation);
}

void setInitialControl(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& lattice = controlledCase.getLattice(NavierStokes{});
  auto& geometry = controlledCase.getGeometry();
  auto& converter = lattice.getUnitConverter();
  auto& control = optiCase.getController();
  auto& parameters = controlledCase.getParameters();

  // Intialize controlled field in superLattice
  OstreamManager clout(std::cout, "setInitialControl");
  T porosity = projection::permeabilityToPorosity(parameters.get<parameters::INITIAL_CONTROL_SCALAR>(), converter);
  clout << "Initial value of the controlled porosity: " << porosity << std::endl;
  fields::set<DissipationEllipsoid::controls_t>(lattice, geometry.getMaterialIndicator({6}), porosity);
  control.template setProjection<projection::SigmoidWithGridTerm<T>>(converter);
  control.template set<DissipationEllipsoid::controls_t>(geometry, 6, lattice);
}

void prepareAdjointLattice(MyOptiCase& optiCase) {
  auto& adjointLattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& geometry = optiCase.getCase(Adjoint{}).getGeometry();

  adjointLattice.setUnitConverter(controlledLattice.getUnitConverter());

  // Define dual dynamics
  dynamics::set<olb::Dual<DissipationEllipsoid::dynamics_t>>(adjointLattice, geometry.getMaterialIndicator({1,6,7,8}));
  boundary::set<boundary::BounceBack>(adjointLattice, geometry.getMaterialIndicator({3,4}));
  adjointLattice.template setParameter<descriptors::OMEGA>(adjointLattice.getUnitConverter().getLatticeRelaxationFrequency());
  adjointLattice.revertStreaming();
}

void setAdjointInitialValues(MyOptiCase& optiCase) {
  auto& adjointCase = optiCase.getCase(Adjoint{});
  auto& adjointLattice = adjointCase.getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& parameters = optiCase.getCase(Controlled{}).getParameters();

  setInitialValues(adjointCase);

  switch (parameters.template get<parameters::OBJECTIVE_TYPE>()) {
    case ObjectiveType::DISSIPATION:
    default: {
      auto& optimalitySystem = optiCase.get(DissipationOS{});
      auto objectiveDomain = optiCase.getCase(Controlled{}).getGeometry().getMaterialIndicator({6});
      optimalitySystem.initializeAdjointProblem(adjointLattice, controlledLattice, objectiveDomain);
    break;
    }
    case ObjectiveType::PRESSUREDROP: {
      auto& optimalitySystem = optiCase.get(PressureDropOS{});
      auto objectiveDomain = optiCase.getCase(Controlled{}).getGeometry().getMaterialIndicator({6,7,8});
      optimalitySystem.initializeAdjointProblem(adjointLattice, controlledLattice, objectiveDomain);
    break;
    }
  }
}

void applyControl(MyOptiCase& optiCase) {
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& lattice = controlledCase.getLattice(NavierStokes{});

  optiCase.getController().template setUpdatedControlsOnField<DissipationEllipsoid::controls_t>(lattice);
  lattice.template setProcessingContext<Array<DissipationEllipsoid::controls_t>>(ProcessingContext::Simulation);
}

MyCase::value_t computeObjective(MyOptiCase& optiCase) {
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& converter = controlledLattice.getUnitConverter();
  auto& parameters = optiCase.getCase(Controlled{}).getParameters();
  MyCase::value_t objective = 0;

  switch (parameters.template get<parameters::OBJECTIVE_TYPE>()) {
    case ObjectiveType::DISSIPATION:
    default: {
      auto& optimalitySystem = optiCase.get(DissipationOS{});
      auto objectiveDomain = optiCase.getCase(Controlled{}).getGeometry().getMaterialIndicator({6});
      optimalitySystem.setObjectiveParameters(controlledLattice,
                                              parameters::OMEGA{}, converter.getLatticeRelaxationFrequency(),
                                              parameters::DT{}, converter.getPhysDeltaT(),
                                              parameters::PHYS_CHAR_VISCOSITY{}, converter.getPhysViscosity(),
                                              parameters::DX{}, converter.getPhysDeltaX(),
                                              parameters::LATTICE_VISCOSITY{}, converter.getLatticeViscosity(),
                                              parameters::CONVERSION_VELOCITY{}, converter.getConversionFactorVelocity(),
                                              parameters::REG_ALPHA{}, parameters.get<parameters::REGULARIZATION_FACTOR>());
      objective = optimalitySystem.evaluateObjectiveF(controlledLattice, objectiveDomain);
      break;
    }

    case ObjectiveType::PRESSUREDROP: {
      auto& optimalitySystem = optiCase.get(PressureDropOS{});
      auto objectiveDomain = optiCase.getCase(Controlled{}).getGeometry().getMaterialIndicator({7,8});
      optimalitySystem.setObjectiveParameters(controlledLattice,
                                              parameters::CONVERSION_PRESSURE{}, converter.getConversionFactorPressure(),
                                              parameters::REG_ALPHA{}, parameters.get<parameters::REGULARIZATION_FACTOR>());
      objective = optimalitySystem.evaluateObjectiveF(controlledLattice, objectiveDomain);
      break;
    }
  }

  return objective;
}

std::vector<MyCase::value_t> computeDerivative(MyOptiCase& optiCase) {
  auto& adjointLattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& converter = controlledLattice.getUnitConverter();
  auto& parameters = optiCase.getCase(Controlled{}).getParameters();
  auto objectiveDomain = optiCase.getCase(Controlled{}).getGeometry().getMaterialIndicator({6});
  std::vector<MyCase::value_t> derivatives;

  optiCase.getController().template setUpdatedProjectionDerivativesOnField<opti::DPROJECTIONDALPHA<DissipationEllipsoid::controls_t>>(adjointLattice);
  adjointLattice.template setProcessingContext<Array<opti::DPROJECTIONDALPHA<DissipationEllipsoid::controls_t>>>(ProcessingContext::Simulation);

  switch (parameters.template get<parameters::OBJECTIVE_TYPE>()) {
    case ObjectiveType::DISSIPATION:
    default: {
      auto& optimalitySystem = optiCase.get(DissipationOS{});
      // Set parameters which are required in the Primal collision
      optimalitySystem.setPrimalCollisionParameters(controlledLattice,
                                                    descriptors::OMEGA{}, converter.getLatticeRelaxationFrequency());
      derivatives = optimalitySystem.evaluateOptimalityCondition(adjointLattice, controlledLattice, objectiveDomain);
      break;
    }

    case ObjectiveType::PRESSUREDROP: {
      auto& optimalitySystem = optiCase.get(PressureDropOS{});

      // Set parameters which are required in the Primal collision
      optimalitySystem.setPrimalCollisionParameters(controlledLattice,
                                                    descriptors::OMEGA{}, converter.getLatticeRelaxationFrequency());
      derivatives = optimalitySystem.evaluateOptimalityCondition(adjointLattice, controlledLattice, objectiveDomain);
      break;
    }
  }

  // Return serial vector containing total derivatives of objective regarding controls
  return derivatives;
}

void getOptiResults(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& adjointLattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& converter = controlledLattice.getUnitConverter();
  auto objectiveDomain = optiCase.getCase(Controlled{}).getGeometry().getMaterialIndicator({6});
  auto controls = optiCase.getControlVector();

  std::size_t iT = optiCase.getOptimizationStep();

  SuperVTMwriter2D<T> vtmWriter("ellipsoidTopology2d");
  SuperLatticePhysVelocity2D velocityF(controlledLattice, converter);
  SuperLatticePhysPressure2D pressureF(controlledLattice, converter);
  SuperLatticeField2D<T,DESCRIPTOR,descriptors::POROSITY> porosityF(controlledLattice);
  SuperLatticeField2D<T,DESCRIPTOR,opti::J> jF(controlledLattice);
  SuperLatticeField2D<T,DESCRIPTOR,opti::DJDF> djdfF(adjointLattice);
  vtmWriter.addFunctor(velocityF);
  vtmWriter.addFunctor(pressureF);
  vtmWriter.addFunctor(porosityF);
  vtmWriter.addFunctor(jF);
  vtmWriter.addFunctor(djdfF);

  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  // Writes the VTK files
  if (iT >= 0) {
    controlledLattice.setProcessingContext(ProcessingContext::Evaluation);
    adjointLattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }

  // Print current solid-fraction
  projection::SigmoidWithGridTerm<T>sigmoid(converter);
  T sumPorosity = 0.0;
  for (auto e : controls) {sumPorosity += sigmoid.project(e);}
  T solidFraction = 1. - sumPorosity / controls.size();
  OstreamManager clout(std::cout, "getOptiResult");
  clout << "Current solid fraction: " << solidFraction << std::endl;
}

MyCase::value_t solvePrimalAndComputeObjective(MyOptiCase& optiCase) {
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

  // Compute Objective value from simulation results
  return computeObjective(optiCase);
}

std::vector<MyCase::value_t> solveAdjointAndComputeDerivative(MyOptiCase& optiCase) {
  auto& adjointCase = optiCase.getCase(Adjoint{});

  // Reset prior simulation lattice
  adjointCase.resetLattices();

  // Prepare new lattices
  prepareAdjointLattice(optiCase);
  setAdjointInitialValues(optiCase);

  // Execute simulation
  simulate(adjointCase);

  // Optimization related IO
  getOptiResults(optiCase);

  // Compute Derivative value from simulation results
  return computeDerivative(optiCase);
}

int main(int argc, char **argv) {
  initialize(&argc, &argv);

  MyCase::ParametersD myCaseParametersD;
  {
    using namespace parameters;
    myCaseParametersD.set<RESOLUTION                  >(    50);
    myCaseParametersD.set<LATTICE_RELAXATION_TIME     >(  0.51);
    myCaseParametersD.set<PHYS_CHAR_LENGTH            >(   0.4);
    myCaseParametersD.set<PHYS_CHAR_VISCOSITY         >( 0.001);
    myCaseParametersD.set<REYNOLDS                    >(    20);
    myCaseParametersD.set<PHYS_CHAR_VELOCITY          >([](MyCase::ParametersD& p){
      return p.get<PHYS_CHAR_VISCOSITY>() *
             p.get<REYNOLDS>() /
             p.get<PHYS_CHAR_LENGTH>();
    });
    myCaseParametersD.set<PHYS_CHAR_DENSITY           >(   1.0);
    myCaseParametersD.set<PHYS_MAX_T                  >(   30.);
    myCaseParametersD.set<PHYS_BOUNDARY_VALUE_UPDATE_T>(   0.1);
    myCaseParametersD.set<INITIAL_CONTROL_SCALAR      >(  5e-4);
    myCaseParametersD.set<REGULARIZATION_FACTOR       >(   0.2);
    myCaseParametersD.set<OBJECTIVE_TYPE              >(ObjectiveType::PRESSUREDROP);
  }

  OptimizerParametersD<MyCase::value_t> optiParameters;
  {
    using namespace parameters;
    optiParameters.set<optimizer::EPS>(1e-10);
    optiParameters.set<optimizer::MAX_OPTIMIZATION_STEPS>(8);
    optiParameters.set<linesearch::MAX_STEP_ATTEMPTS>(20);
    optiParameters.set<linesearch::STEP_CONDITION>("Wolfe");
  }

  Mesh mesh = createMesh(myCaseParametersD);
  MyCase myCase(myCaseParametersD, mesh);
  prepareGeometry(myCase);
  prepareLattice(myCase);

  MyCase adjointCase(myCaseParametersD, mesh);
  prepareGeometry(adjointCase);

  MyOptiCase optiCase;
  optiCase.setCase<Controlled>(myCase);
  optiCase.setCase<Adjoint>(adjointCase);

  DissipationEllipsoid dissipationOS;
  optiCase.set<DissipationOS>(dissipationOS);
  PressureDropEllipsoid pressureDropOS;
  optiCase.set<PressureDropOS>(pressureDropOS);

  setInitialControl(optiCase);
  optiCase.setObjective(solvePrimalAndComputeObjective);
  optiCase.setDerivative(solveAdjointAndComputeDerivative);

  OptimizerLBFGS<MyCase::value_t,std::vector<MyCase::value_t>> optimizer(
    optiCase.getController().size(), optiParameters);

  optimizer.optimize(optiCase);
}
