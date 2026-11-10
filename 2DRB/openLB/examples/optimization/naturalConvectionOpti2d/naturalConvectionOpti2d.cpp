/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2026 Ecem Öz, Shota Ito
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
  NavierStokes, Lattice<double, descriptors::D2Q9<descriptors::FORCE, descriptors::SCALAR, opti::DCnseDFxphi, opti::DCadeDFxlambda, opti::F, opti::DJDF>>,
  Temperature,  Lattice<double, descriptors::D2Q5<descriptors::VELOCITY, descriptors::SCALAR, opti::DCnseDGxphi, opti::DCadeDGxlambda, opti::G, opti::DJDG>>
>;

using MyOptiCase = OptiCaseAdjoint<
  Controlled, MyCase,
  Adjoint, MyCase,
  Reference, MyCase
>;

using NSEBulkDynamics = ForcedBGKdynamics<MyCase::value_t,
                                          MyCase::descriptor_t_of<NavierStokes>>;

using ADEBulkDynamics = AdvectionDiffusionBGKdynamics<MyCase::value_t,
                                                      MyCase::descriptor_t_of<Temperature>>;

using ObjectiveF = functors::L2DistanceF<functors::VelocityF>;

namespace olb::parameters{
struct ADJOINT_GRADIENTS  : public descriptors::FIELD_BASE<1> { };
struct FDQ_GRADIENTS      : public descriptors::FIELD_BASE<1> { };
struct ADJOINT_MAX_T : public descriptors::FIELD_BASE<1> { };
}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t_of<NavierStokes>;
  const T dx = parameters.get<parameters::PHYS_DELTA_X>();
  const Vector domainExtend = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector extend = {domainExtend[0] + dx, domainExtend[1] + dx};
  const Vector domainOrigin = parameters.get<parameters::ORIGIN>();
  const Vector origin = {domainOrigin[0] - (T) dx / 2, domainOrigin[1] - (T) dx / 2};
  IndicatorCuboid2D<T> cuboid(extend, origin);

  Mesh<T,MyCase::d> mesh(cuboid, dx, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({ false, false });
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  using T = MyCase::value_t_of<NavierStokes>;

  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  const T dx                = parameters.get<parameters::PHYS_DELTA_X>();
  const Vector extend       = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector origin       = parameters.get<parameters::ORIGIN>();

  geometry.rename(0, 4, {1, 1});

  IndicatorCuboid2D<T> bulkCuboid(extend, origin);
  geometry.rename(4, 1, bulkCuboid);

  Vector extendwallleft{dx, extend[1] + dx};
  Vector originwallleft{origin[0] - (T) dx / 2, origin[1] - (T) dx / 2};
  IndicatorCuboid2D<T> wallleft(extendwallleft, originwallleft);
  geometry.rename( 4, 2, 1, wallleft  );

  Vector extendwallright{dx, extend[1] + dx};
  Vector originwallright{origin[0] + extend[0], origin[1] - (T) dx / 2};
  IndicatorCuboid2D<T> wallright(extendwallright, originwallright);

  geometry.rename( 4, 3, 1, wallright );

  geometry.clean();
  geometry.checkForErrors();
  geometry.print();
}

void prepareLattice(MyCase& myCase) {
  using T             = MyCase::value_t_of<NavierStokes>;
  using NSEDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ADEDESCRIPTOR = MyCase::descriptor_t_of<Temperature>;

  auto& geometry    = myCase.getGeometry();
  auto& NSElattice  = myCase.getLattice(NavierStokes{});
  auto& ADElattice  = myCase.getLattice(Temperature{});
  auto& parameters  = myCase.getParameters();

  NSElattice.setUnitConverter<ThermalUnitConverter<T,NSEDESCRIPTOR,ADEDESCRIPTOR>>(
    (T) parameters.get<parameters::PHYS_DELTA_X>(),
    (T) (parameters.get<parameters::LATTICE_RELAXATION_TIME>() - (T) 1 / 2) /
          descriptors::invCs2<T,NSEDESCRIPTOR>() *
          parameters.get<parameters::PHYS_DELTA_X>() *
          parameters.get<parameters::PHYS_DELTA_X>() /
          parameters.get<parameters::PHYS_KINEMATIC_VISCOSITY>(),
    (T) parameters.get<parameters::PHYS_CHAR_LENGTH>(),
    (T) parameters.get<parameters::PHYS_CHAR_VELOCITY>(),
    (T) parameters.get<parameters::PHYS_KINEMATIC_VISCOSITY>(),
    (T) parameters.get<parameters::PHYS_CHAR_DENSITY>(),
    (T) parameters.get<parameters::PHYS_THERMAL_CONDUCTIVITY>(),
    (T) parameters.get<parameters::PHYS_HEAT_CAPACITY>(),
    (T) parameters.get<parameters::PHYS_THERMAL_EXPANSION>(),
    (T) parameters.get<parameters::T_COLD>(),
    (T) parameters.get<parameters::T_HOT>()
  );
  const auto& converter = NSElattice.getUnitConverter();
  ADElattice.setUnitConverter(converter);
  converter.print();

  // Dynamics
  dynamics::set<ForcedBGKdynamics>(NSElattice, geometry.getMaterialIndicator({1}));
  boundary::set<boundary::BounceBack>(NSElattice, geometry.getMaterialIndicator({2,3,4}));
  dynamics::set<AdvectionDiffusionBGKdynamics>(ADElattice, geometry.getMaterialIndicator({1}));
  boundary::set<boundary::BounceBack>(ADElattice, geometry.getMaterialIndicator(4));
  boundary::set<boundary::AdvectionDiffusionDirichlet>(ADElattice, geometry.getMaterialIndicator({ 2, 3 }));

  // Coupling
  T boussinesqForcePrefactor =   parameters.get<parameters::GRAVITATIONAL_ACC>()
                                 * converter.getConversionFactorTime()
                                 * converter.getCharPhysTemperatureDifference()
                                 * converter.getPhysThermalExpansionCoefficient()
                                 / converter.getConversionFactorVelocity();

  auto& coupling = myCase.setCouplingOperator(
    "Boussinesq",
    NavierStokesAdvectionDiffusionCoupling{},
    names::NavierStokes{}, NSElattice,
    names::Temperature{},  ADElattice
  );
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::T0>(converter.getLatticeTemperature(parameters.get<parameters::T_COLD>()));
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::FORCE_PREFACTOR>(boussinesqForcePrefactor * Vector<T,2>{0.0,1.0});
  coupling.restrictTo(geometry.getMaterialIndicator({1}));
}

void setInitialValues(MyCase& myCase) {
  using T               = MyCase::value_t_of<NavierStokes>;
  auto& geometry        = myCase.getGeometry();
  auto& NSElattice      = myCase.getLattice(NavierStokes{});
  auto& ADElattice      = myCase.getLattice(Temperature{});
  const auto& converter = NSElattice.getUnitConverter();

  const T Tcold         = converter.getCharPhysLowTemperature();
  const T Thot          = converter.getCharPhysHighTemperature();
  const T Tmean         = (Thot + Tcold) / 2.;

  /// define initial conditions
  momenta::setTemperature(ADElattice, geometry.getMaterialIndicator(1), Tmean);
  momenta::setTemperature(ADElattice, geometry.getMaterialIndicator(2), Thot);
  momenta::setTemperature(ADElattice, geometry.getMaterialIndicator(3), Tcold);

  NSElattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  ADElattice.setParameter<descriptors::OMEGA>(converter.getLatticeThermalRelaxationFrequency());

  /// Make the lattice ready for simulation
  NSElattice.initialize();
  ADElattice.initialize();
}

void getResults(MyCase& myCase, std::size_t iT, util::Timer<MyCase::value_t> timer, std::string vtkName) {
  using T = MyCase::value_t_of<NavierStokes>;
  using NSEDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ADEDESCRIPTOR = MyCase::descriptor_t_of<Temperature>;

  auto& NSElattice        = myCase.getLattice(NavierStokes{});
  auto& ADElattice        = myCase.getLattice(Temperature{});
  const auto& converter   = NSElattice.getUnitConverter();
  auto& parameters        = myCase.getParameters();

  const int statIter      = converter.getLatticeTime(parameters.get<parameters::PHYS_STAT_ITER_T>());
  const int vtkIter       = converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T>());
  SuperVTMwriter2D<T> vtkWriter(vtkName);

  if (iT == 0) {
    vtkWriter.createMasterFile();
  }

  if (iT % vtkIter == 0 && iT > 0) {
    NSElattice.setProcessingContext(ProcessingContext::Evaluation);
    ADElattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticePhysVelocity2D<T, NSEDESCRIPTOR>                    velocity(NSElattice, converter);
    SuperLatticePhysPressure2D<T, NSEDESCRIPTOR>                    pressure(NSElattice, converter);
    SuperLatticePhysTemperature2D<T, NSEDESCRIPTOR, ADEDESCRIPTOR>  temperature(ADElattice, converter);

    vtkWriter.addFunctor(velocity);
    vtkWriter.addFunctor(pressure);
    vtkWriter.addFunctor(temperature);

    vtkWriter.write(iT);
  }

  if (iT % statIter == 0) {
    timer.printStep();
    NSElattice.getStatistics().print(iT,converter.getPhysTime(iT));
    ADElattice.getStatistics().print(iT,converter.getPhysTime(iT));
  }
}

void getAdjointResults(MyOptiCase& optiCase, util::Timer<MyCase::value_t> timer, std::size_t iT) {
  using T = MyCase::value_t;
  using NSEDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ADEDESCRIPTOR = MyCase::descriptor_t_of<Temperature>;

  auto& adjointCase = optiCase.getCase(Adjoint{});
  auto& parameters  = adjointCase.getParameters();
  auto& adjointNSE  = adjointCase.getLattice(NavierStokes{});
  auto& adjointADE  = adjointCase.getLattice(Temperature{});
  auto& converter   = adjointNSE.getUnitConverter();
  const int statIter      = converter.getLatticeTime(parameters.get<parameters::PHYS_STAT_ITER_T>());
  const int vtkIter       = converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T>());

  SuperVTMwriter2D<T> vtkWriter("adjointCavity2d");
  SuperLatticeField2D<T, NSEDESCRIPTOR, opti::DJDF> sourceDJDF(adjointNSE);
  SuperLatticeField2D<T, ADEDESCRIPTOR, opti::DJDG> sourceDJDG(adjointADE);
  SuperLatticeField2D<T, NSEDESCRIPTOR, opti::DCnseDFxphi> dCnsedFxphi(adjointNSE);
  SuperLatticeField2D<T, NSEDESCRIPTOR, opti::DCadeDFxlambda> dCadedFxlambda(adjointNSE);
  SuperLatticeField2D<T, ADEDESCRIPTOR, opti::DCnseDGxphi> dCnsedGxphi(adjointADE);
  SuperLatticeField2D<T, ADEDESCRIPTOR, opti::DCadeDGxlambda> dCadedGxlambda(adjointADE);
  SuperLatticeField2D<T, NSEDESCRIPTOR, opti::F> primalF(adjointNSE);
  SuperLatticeField2D<T, ADEDESCRIPTOR, opti::G> primalG(adjointADE);
  SuperLatticeField2D<T, NSEDESCRIPTOR, descriptors::POPULATION> adjointPhi(adjointNSE);
  adjointPhi.getName() = "adjointPhi";
  SuperLatticeField2D<T, ADEDESCRIPTOR, descriptors::POPULATION> adjointLambda(adjointADE);
  adjointLambda.getName() = "adjointLambda";

  vtkWriter.addFunctor(sourceDJDF);
  vtkWriter.addFunctor(sourceDJDG);
  vtkWriter.addFunctor(dCnsedFxphi);
  vtkWriter.addFunctor(dCadedFxlambda);
  vtkWriter.addFunctor(dCnsedGxphi);
  vtkWriter.addFunctor(dCadedGxlambda);
  vtkWriter.addFunctor(primalF);
  vtkWriter.addFunctor(primalG);
  vtkWriter.addFunctor(adjointPhi);
  vtkWriter.addFunctor(adjointLambda);

  if (iT == 0) {
    vtkWriter.createMasterFile();
  }

  if (iT % statIter == 0) {
    timer.printStep();
    adjointNSE.getStatistics().print(iT,converter.getPhysTime(iT));
    adjointADE.getStatistics().print(iT,converter.getPhysTime(iT));
  }

  if (iT % vtkIter == 0 ) {
    adjointNSE.setProcessingContext(ProcessingContext::Evaluation);
    adjointADE.setProcessingContext(ProcessingContext::Evaluation);
    vtkWriter.write(iT);
  }
}

void simulate(MyCase& myCase, const std::string& vtkName) {
  OstreamManager clout(std::cout,"PrimalSimulation");
  clout << "Starting Simulation ..." << std::endl;
  using T = MyCase::value_t_of<NavierStokes>;
  auto& NSElattice = myCase.getLattice(NavierStokes{});
  auto& ADElattice = myCase.getLattice(Temperature{});
  const auto& converter = NSElattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  const std::size_t iTmax = converter.getLatticeTime(parameters.get<parameters::PHYS_MAX_T>());
  const std::size_t convIter = parameters.get<parameters::CONV_ITER>();
  util::ValueTracer<T> converge(6, parameters.get<parameters::CONVERGENCE_PRECISION>());

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();
  for (std::size_t iT = 0; iT < iTmax; ++iT)
  {
    if (converge.hasConverged()) {
      parameters.set<parameters::ADJOINT_MAX_T>(iT);
      clout << "Simulation converged." << std::endl;
      break;
    }

    NSElattice.collideAndStream();
    ADElattice.collideAndStream();
    myCase.getOperator("Boussinesq").apply();

    getResults(myCase, iT, timer, vtkName);

    if(iT % convIter == 0) {
      converge.takeValue(NSElattice.getStatistics().getAverageRho(), true);
    }

    timer.update(iT);
  }

  timer.stop();
  timer.printSummary();
}


void simulateAdjoint(MyOptiCase& optiCase) {
  OstreamManager clout(std::cout,"AdjointSimulation");
  clout << "Starting Simulation ..." << std::endl;
  using T = MyCase::value_t_of<NavierStokes>;
  auto& adjointCase = optiCase.getCase(Adjoint{});
  auto& adjointNSE  = adjointCase.getLattice(NavierStokes{});
  auto& adjointADE  = adjointCase.getLattice(Temperature{});

  auto& parameters = adjointCase.getParameters();
  const std::size_t iTmax = parameters.get<parameters::ADJOINT_MAX_T>();

  util::Timer<T> timer(iTmax, adjointCase.getGeometry().getStatistics().getNvoxel());
  timer.start();
  for (std::size_t iT = 0; iT < iTmax; ++iT) {
    // reverse order compared to primal
    adjointCase.getOperator("AdjointBoussinesq").apply();
    adjointADE.collideAndStream();
    adjointNSE.collideAndStream();

    getAdjointResults(optiCase, timer, iT);

    timer.update(iT);
  }
  timer.stop();
  timer.printSummary();
}


void prepareAdjointLattice(MyOptiCase& optiCase) {
  using T             = MyCase::value_t_of<NavierStokes>;
  auto& adjointCase = optiCase.getCase(Adjoint{});
  auto& adjointNSELattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& adjointADELattice = optiCase.getCase(Adjoint{}).getLattice(Temperature{});
  auto& controlledNSELattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& geometry = optiCase.getCase(Adjoint{}).getGeometry();
  auto& parameters = optiCase.getCase(Controlled{}).getParameters();

  // Set adjoint unit converter as the controlled case
  adjointNSELattice.setUnitConverter(controlledNSELattice.getUnitConverter());
  auto& converter = adjointNSELattice.getUnitConverter();

  adjointADELattice.setUnitConverter(converter);

  // Define dual physics
  dynamics::set<DualNavierStokesBGKDynamics>(adjointNSELattice, geometry.getMaterialIndicator({1}));
  dynamics::set<DualAdvectionDiffusionBGKDynamics>(adjointADELattice, geometry.getMaterialIndicator({1}));

  // Dual boundary conditions
  boundary::set<boundary::BounceBack>(adjointNSELattice, geometry.getMaterialIndicator({2,3,4}));
  boundary::set<boundary::BounceBack>(adjointADELattice, geometry.getMaterialIndicator({4}));
  boundary::set<boundary::AdvectionDiffusionDirichlet>(adjointADELattice, geometry.getMaterialIndicator({2,3}));

  adjointNSELattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  adjointADELattice.template setParameter<descriptors::OMEGA>(converter.getLatticeThermalRelaxationFrequency());

  // Dual coupling operator
  T boussinesqForcePrefactor = (parameters.get<parameters::GRAVITATIONAL_ACC>())         // gravity g
                                 * converter.getConversionFactorTime()                   // unit conversion
                                 * converter.getCharPhysTemperatureDifference()          // T_hot - T_cold
                                 * parameters.get<parameters::PHYS_THERMAL_EXPANSION>()  // alpha; this is our control parameter
                                 / converter.getConversionFactorVelocity();              // unit conversion

  auto& adjointCoupling = adjointCase.setCouplingOperator(
    "AdjointBoussinesq",
    AdjointBoussinesqCoupling{},
    names::NavierStokes{}, adjointNSELattice,
    names::Temperature{}, adjointADELattice
  );
  adjointCoupling.setParameter<AdjointBoussinesqCoupling::T0>(converter.getLatticeTemperature(parameters.get<parameters::T_COLD>()));
  adjointCoupling.setParameter<AdjointBoussinesqCoupling::FORCE_PREFACTOR>(boussinesqForcePrefactor * Vector<T,2>{0.0,1.0});
  adjointCoupling.setParameter<AdjointBoussinesqCoupling::OMEGA_NSE>(converter.getLatticeRelaxationFrequency());
  adjointCoupling.setParameter<AdjointBoussinesqCoupling::OMEGA_ADE>(converter.getLatticeThermalRelaxationFrequency());
  adjointCoupling.restrictTo(geometry.getMaterialIndicator({1}));

  // Dual streaming operation
  adjointNSELattice.revertStreaming();
  adjointADELattice.revertStreaming();
}


void setAdjointInitialValues(MyOptiCase& optiCase) {
  using T = MyCase::value_t_of<NavierStokes>;
  auto& adjointNSELattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& adjointADELattice = optiCase.getCase(Adjoint{}).getLattice(Temperature{});
  auto& controlledNSELattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& controlledADELattice = optiCase.getCase(Controlled{}).getLattice(Temperature{});
  auto& geometry = optiCase.getCase(Adjoint{}).getGeometry();
  auto& referenceNSELattice = optiCase.getCase(Reference{}).getLattice(NavierStokes{});
  auto& converter = adjointNSELattice.getUnitConverter();
  auto& parameters = optiCase.getCase(Controlled{}).getParameters();
  auto objectiveDomain = geometry.getMaterialIndicator({1});

  adjointNSELattice.initialize();
  adjointADELattice.initialize();

  // Initialize dual problem
  auto djdf = makeWriteFunctorO<functors::DerivativeF<ObjectiveF,descriptors::POPULATION,NSEBulkDynamics>,
                                opti::DJDF>(controlledNSELattice);
  djdf->template setParameter<parameters::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
  djdf->template setParameter<parameters::NORMALIZE>(computeL2Norm<ObjectiveF::Reference>(referenceNSELattice, objectiveDomain, converter.getPhysDeltaX()));
  djdf->template setParameter<parameters::FACTOR>(util::pow(converter.getPhysDeltaX(),2));
  djdf->restrictTo(objectiveDomain);
  djdf->apply();

  // Provide fields required by the dual collision operator
  copyFields<descriptors::POPULATION,opti::F>(controlledNSELattice, adjointNSELattice);
  copyFields<descriptors::POPULATION,opti::G>(controlledADELattice, adjointADELattice);
  copyFields<opti::DJDF,opti::DJDF>(controlledNSELattice, adjointNSELattice);

  // Compute DJDG; use coupling as both adjoint lattices are required
  // The objective depends on g via the velocity correction due to forcing
  SuperLatticeCoupling djdg(
    AdjointBoussinesqObjectiveDerivativeO{},
    names::NavierStokes{}, adjointNSELattice,
    names::Temperature{}, adjointADELattice
  );

  // Parameterize coupling
  const T boussinesqForcePrefactor = (parameters.get<parameters::GRAVITATIONAL_ACC>())        // gravity g
                                       * converter.getConversionFactorTime()                  // unit conversion
                                       * converter.getCharPhysTemperatureDifference()         // t_hot - t_cold
                                       * parameters.get<parameters::PHYS_THERMAL_EXPANSION>() // alpha; this is our control parameter
                                       / converter.getConversionFactorVelocity();             // unit conversion
  // Pass reference velocity for objective evaluation
  copyFields<ObjectiveF::Reference,AdjointBoussinesqObjectiveDerivativeO::REFERENCE_VELOCITY>(referenceNSELattice, adjointNSELattice);
  djdg.setParameter<AdjointBoussinesqObjectiveDerivativeO::T0>(converter.getLatticeTemperature(parameters.get<parameters::T_COLD>()));
  djdg.setParameter<AdjointBoussinesqObjectiveDerivativeO::FORCE_PREFACTOR>(boussinesqForcePrefactor * Vector<T,2>{0.0,1.0});
  djdg.setParameter<AdjointBoussinesqObjectiveDerivativeO::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
  djdg.setParameter<AdjointBoussinesqObjectiveDerivativeO::NORMALIZE>(computeL2Norm<ObjectiveF::Reference>(referenceNSELattice,
    objectiveDomain, converter.getPhysDeltaX()));
  djdg.setParameter<AdjointBoussinesqObjectiveDerivativeO::FACTOR>(util::pow(converter.getPhysDeltaX(),2));
  djdg.restrictTo(objectiveDomain);
  djdg.apply();
}

void setInitialControl(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  const T alpha = optiCase.getCase(Controlled{}).getParameters().get<parameters::PHYS_THERMAL_EXPANSION>();
  optiCase.getController().set({1.1 * alpha}); // Perturbate the thermal expansion coefficient by 10%
}

void applyControl(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& parameters = controlledCase.getParameters();
  const T thermalExpansion = optiCase.getController().get()[0];
  parameters.set<parameters::PHYS_THERMAL_EXPANSION>(thermalExpansion);
}


MyCase::value_t objectiveF(MyOptiCase& optiCase) {
  auto& referenceNSELattice = optiCase.getCase(Reference{}).getLattice(NavierStokes{});
  auto& controlledNSELattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& converter = controlledNSELattice.getUnitConverter();
  auto& geometry = optiCase.getCase(Controlled{}).getGeometry();
  auto objectiveDomain = geometry.getMaterialIndicator({1});

  // Evaluate functor for objective computation
  auto objectiveO = makeWriteFunctorO<ObjectiveF,opti::J>(controlledNSELattice);
  objectiveO->restrictTo(objectiveDomain);

  // Write velocity field in the reference case
  auto velocityO = makeWriteFunctorO<functors::VelocityF,ObjectiveF::Reference>(referenceNSELattice);
  velocityO->restrictTo(objectiveDomain);
  velocityO->template setParameter<parameters::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
  velocityO->apply();

  // Get solution from the reference simulation for the inverse problem
  copyFields<ObjectiveF::Reference,ObjectiveF::Reference>(referenceNSELattice, controlledNSELattice);
  objectiveO->template setParameter<parameters::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
  objectiveO->template setParameter<parameters::NORMALIZE>(computeL2Norm<ObjectiveF::Reference>(referenceNSELattice, objectiveDomain, converter.getPhysDeltaX()));
  objectiveO->apply();
  return integrateField<opti::J>(controlledNSELattice, objectiveDomain, converter.getPhysDeltaX())[0];
}


std::vector<MyCase::value_t> derivativeF(MyOptiCase& optiCase) {
  using T = MyCase::value_t_of<NavierStokes>;
  auto& adjointNSELattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& adjointADELattice = optiCase.getCase(Adjoint{}).getLattice(Temperature{});
  auto& parameters = optiCase.getCase(Controlled{}).getParameters();
  auto& geometry = optiCase.getCase(Controlled{}).getGeometry();
  auto& referenceNSELattice = optiCase.getCase(Reference{}).getLattice(NavierStokes{});
  auto& converter = adjointNSELattice.getUnitConverter();
  auto objectiveDomain = geometry.getMaterialIndicator({1});

  adjointNSELattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  adjointADELattice.template setParameter<descriptors::OMEGA>(converter.getLatticeThermalRelaxationFrequency());

  // Total derivative of objective regarding the controls
  SuperLatticeCoupling djdalpha(
    AdjointBoussinesqOptimalityO{},
    names::NavierStokes{}, adjointNSELattice,
    names::Temperature{}, adjointADELattice
  );

  const T DboussinesqForcePrefactorDALPHA = (parameters.get<parameters::GRAVITATIONAL_ACC>())   // gravity g
                                 * converter.getConversionFactorTime()                          // unit conversion
                                 * converter.getCharPhysTemperatureDifference()                 // t_hot - t_cold
                                 / converter.getConversionFactorVelocity();                     // unit conversion

  const T boussinesqForcePrefactor = DboussinesqForcePrefactorDALPHA                            // gravity g
                             * parameters.get<parameters::PHYS_THERMAL_EXPANSION>();            // alpha; this is our control parameter

  copyFields<ObjectiveF::Reference,AdjointBoussinesqOptimalityO::REFERENCE_VELOCITY>(referenceNSELattice, adjointNSELattice);
  djdalpha.setParameter<AdjointBoussinesqOptimalityO::T0>(converter.getLatticeTemperature(parameters.get<parameters::T_COLD>()));
  djdalpha.setParameter<AdjointBoussinesqOptimalityO::OMEGA_NSE>(converter.getLatticeRelaxationFrequency());
  djdalpha.setParameter<AdjointBoussinesqOptimalityO::OMEGA_ADE>(converter.getLatticeThermalRelaxationFrequency());
  djdalpha.setParameter<AdjointBoussinesqOptimalityO::D_FORCE_PREFACTOR_D_ALPHA>(DboussinesqForcePrefactorDALPHA * Vector<T,2>{0.0,1.0});
  djdalpha.setParameter<AdjointBoussinesqOptimalityO::FORCE_PREFACTOR>(boussinesqForcePrefactor * Vector<T,2>{0.0,1.0});
  djdalpha.setParameter<AdjointBoussinesqOptimalityO::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
  djdalpha.setParameter<AdjointBoussinesqOptimalityO::NORMALIZE>(computeL2Norm<ObjectiveF::Reference>(
    referenceNSELattice, objectiveDomain, converter.getPhysDeltaX()));
  djdalpha.setParameter<AdjointBoussinesqOptimalityO::FACTOR>(util::pow(converter.getPhysDeltaX(),2));
  djdalpha.apply();

  adjointNSELattice.setProcessingContext(ProcessingContext::Evaluation);
  auto cellwiseContribution = getSerializedFromField<AdjointBoussinesqOptimalityO::SENSITIVITY>(adjointNSELattice, geometry.getMaterialIndicator({1}));
  std::vector<T> gradient;
  gradient.push_back(util::sumElements(cellwiseContribution));
  return gradient;
}

MyCase::value_t computeObjective(MyOptiCase& optiCase) {
  auto& controlledCase = optiCase.getCase(Controlled{});
  controlledCase.resetLattices();

  applyControl(optiCase);
  prepareLattice(controlledCase);
  setInitialValues(controlledCase);
  simulate(controlledCase, "controlledCavity2d");
  return objectiveF(optiCase);
}

std::vector<MyCase::value_t> computeDerivative(MyOptiCase& optiCase) {
  auto& adjointCase = optiCase.getCase(Adjoint{});
  adjointCase.resetLattices();

  applyControl(optiCase);
  prepareAdjointLattice(optiCase);
  setAdjointInitialValues(optiCase);
  simulateAdjoint(optiCase);
  return derivativeF(optiCase);
}

void computeDerivativesWithFDQ(MyOptiCase& optiCase, MyCase::value_t objective, MyCase::value_t shift = 1e-8) {
  using T = MyCase::value_t;
  Vector<T,1> gradients(0.0);

  auto control = optiCase.getController().get();
  auto shiftedControl = control;
  shiftedControl[0] += shift;

  auto shiftedObjective = optiCase.computeObjective(shiftedControl, 0);
  gradients[0] = (shiftedObjective - objective) / (shift);

  optiCase.getCase(Controlled{}).getParameters().set<parameters::FDQ_GRADIENTS>(gradients[0]);
  optiCase.getController().set(control);
}

auto computeDerivativesWithAdjoints(MyCase::ParametersD& parameters, MyOptiCase& optiCase){
  auto objective = optiCase.computeObjective(optiCase.getController().get());
  auto derivative = optiCase.computeDerivatives(optiCase.getController().get());
  parameters.set<parameters::ADJOINT_GRADIENTS>({derivative[0]});
  return objective;
}

int main(int argc, char* argv[])
{
  using T = MyCase::value_t_of<NavierStokes>;
  initialize(&argc, &argv);
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<PHYS_MAX_T                 >(1000);
    myCaseParameters.set<RESOLUTION                 >(64);
    myCaseParameters.set<GRAVITATIONAL_ACC          >(9.81);
    myCaseParameters.set<ORIGIN                     >({0., 0.});
    myCaseParameters.set<LATTICE_RELAXATION_TIME    >(0.65);
    myCaseParameters.set<PHYS_THERMAL_EXPANSION     >(3.41e-3);
    myCaseParameters.set<PHYS_THERMAL_CONDUCTIVITY  >(25.684e-3);
    myCaseParameters.set<PHYS_HEAT_CAPACITY         >(1.01309e3);
    myCaseParameters.set<PHYS_CHAR_DENSITY          >(1.19);
    myCaseParameters.set<PHYS_KINEMATIC_VISCOSITY   >(1.5126e-5);
    myCaseParameters.set<PRANDTL                    >(0.71);
    myCaseParameters.set<RAYLEIGH                   >(1e4);
    myCaseParameters.set<T_HOT                      >(285.15);
    myCaseParameters.set<T_COLD                     >(275.15);
    myCaseParameters.set<CONV_ITER                  >(1000);
    myCaseParameters.set<PHYS_VTK_ITER_T            >(0.1);
    myCaseParameters.set<PHYS_STAT_ITER_T           >(1.);
    myCaseParameters.set<CONVERGENCE_PRECISION      >(1e-9);
  }
  const int N              = myCaseParameters.get<parameters::RESOLUTION               >();
  const std::size_t Ra     = myCaseParameters.get<parameters::RAYLEIGH                 >();
  const T physViscosity    = myCaseParameters.get<parameters::PHYS_KINEMATIC_VISCOSITY >();
  const T Pr               = myCaseParameters.get<parameters::PRANDTL                  >();
  const T g                = myCaseParameters.get<parameters::GRAVITATIONAL_ACC        >();
  const T Thot             = myCaseParameters.get<parameters::T_HOT                    >();
  const T Tcold            = myCaseParameters.get<parameters::T_COLD                   >();
  const T thermalExpansion = myCaseParameters.get<parameters::PHYS_THERMAL_EXPANSION   >();
  // lx from Rayleigh number
  const T lx               = util::pow(Ra * physViscosity * physViscosity / (Pr * g * (Thot - Tcold) * thermalExpansion),
                               (T) 1 / 3);  // length of the square
  myCaseParameters.set<parameters::PHYS_CHAR_LENGTH>(lx);
  myCaseParameters.set<parameters::PHYS_DELTA_X>(lx / N);
  myCaseParameters.set<parameters::DOMAIN_EXTENT>([](MyCase::ParametersD& p) -> Vector<T, MyCase::d> {
    return {p.get<parameters::PHYS_CHAR_LENGTH>(), p.get<parameters::PHYS_CHAR_LENGTH>()};
  });
  const T charU = physViscosity / lx * sqrt((T) Ra / Pr);
  myCaseParameters.set<parameters::PHYS_CHAR_VELOCITY>(charU);

  // OptiCase used for adjoint gradient computation
  MyOptiCase optiCase;
  Mesh mesh = createMesh(myCaseParameters);

  MyCase myCase(myCaseParameters, mesh);
  prepareGeometry(myCase);
  prepareLattice(myCase);

  MyCase adjointCase(myCaseParameters, mesh);
  prepareGeometry(adjointCase);

  MyCase referenceCase(myCaseParameters, mesh);
  prepareGeometry(referenceCase);
  prepareLattice(referenceCase);
  setInitialValues(referenceCase);
  simulate(referenceCase, "referenceCavity2d");

  optiCase.setCase<Controlled>(myCase);
  optiCase.setCase<Adjoint>(adjointCase);
  optiCase.setCase<Reference>(referenceCase);

  setInitialControl(optiCase);

  optiCase.setObjective(computeObjective);
  optiCase.setDerivative(computeDerivative);

  // Evaluate Adjoint-based gradients
  auto objective = computeDerivativesWithAdjoints(myCaseParameters, optiCase);

  // Evaluate FDQ-based gradients
  computeDerivativesWithFDQ(optiCase, objective);

  // Compare gradients
  OstreamManager clout( std::cout,"main" );
  auto adjoint = myCaseParameters.template get<parameters::ADJOINT_GRADIENTS>();
  auto fdq = myCaseParameters.template get<parameters::FDQ_GRADIENTS>();
  auto error = std::abs(adjoint - fdq);
  clout << " : Adjoint = " << adjoint
        << " , FDQ = " << fdq
        << " , Error = " << error
        << std::endl;
}
