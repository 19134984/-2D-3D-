/*  Lattice Boltzmann sample, written in C++, using the OpenLBlibrary
 *
 *  Copyright (C) 2026 Shota Ito
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

// This app is based on the work in doi: 10.1016/j.camwa.2024.05.026.
// The simulation setup for a single component is based on Krüger et al. 2017 page 324 (doi: 10.1007/978-3-319-44649-3).

// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::opti;

namespace olb::parameters {

struct HEIGHT            : public descriptors::FIELD_BASE<1> { };

struct REACTION_CONSTANT : public descriptors::FIELD_BASE<1> { };
struct STOICHIOMETRY     : public descriptors::FIELD_BASE<3> { };
struct REACTION_ORDER    : public descriptors::FIELD_BASE<3> { };

struct REACTION_UPDATE_T : public descriptors::FIELD_BASE<1> { };
}

using MyCase = Case<
  Concentration0, Lattice<double, descriptors::D3Q7<descriptors::VELOCITY>>,
  Concentration1, Lattice<double, descriptors::D3Q7<descriptors::VELOCITY>>,
  Concentration2, Lattice<double, descriptors::D3Q7<descriptors::VELOCITY>>
>;
using MyADfCase = MyCase::exchange_value_t<util::ADf<double, 3>>;
using MyOptiCase = OptiCaseADf<
  Controlled, MyCase,
  Derivatives, MyADfCase,
  Reference, MyCase
>;

using OBJECTIVE = functors::L2DistanceF<functors::ConcentrationF>;

template <typename PARAMETERS>
auto createMesh(PARAMETERS& parameters) {
  using T = PARAMETERS::value_t;
  T height = parameters.template get<parameters::DOMAIN_EXTENT>()[2];
  T width = parameters.template get<parameters::DOMAIN_EXTENT>()[1];
  T length = parameters.template get<parameters::DOMAIN_EXTENT>()[0];
  const Vector<T,3> origin;
  const Vector<T,3> extend(length, height, width);
  IndicatorCuboid3D<T> cuboid(extend, origin);

  Mesh<T,3> mesh(cuboid, parameters.template get<parameters::PHYS_DELTA_X>(), singleton::mpi().getSize());

  mesh.setOverlap(parameters.template get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({false, false, true});

  return mesh;
 }

// prepare numerical domain
template <typename CASE>
void prepareGeometry(CASE& myCase){
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = CASE::value_t;
  auto& parameters = myCase.getParameters();
  auto& superGeometry   = myCase.getGeometry();
  superGeometry.setWriteIncrementalVTK(true);

  // === MATERIAL NUMBERS === //
  // solid: MN = 0            //
  // bulk: MN = 1             //
  // inlet: MN = 2            //
  // outlet: MN = 3           //
  // topWall: MN = 4          //
  // bottomWall: MN = 5       //
  // ======================== //

  superGeometry.rename(0, 1);
  T eps = parameters.template get<parameters::PHYS_DELTA_X>() * T(0.5);
  T height = parameters.template get<parameters::DOMAIN_EXTENT>()[2];
  T width = parameters.template get<parameters::DOMAIN_EXTENT>()[1];
  T length = parameters.template get<parameters::DOMAIN_EXTENT>()[0];
  T dx = parameters.template get<parameters::PHYS_DELTA_X>();

  {
    Vector<T,3>origin(-eps, -eps, -eps - dx);
    Vector<T,3>extend(eps * T(2.), height + T(2.) * eps, width + T(2.) * eps + T(2.) * dx);
    IndicatorCuboid3D<T> inflow(extend, origin);
    superGeometry.rename(1, 2, inflow);

    origin[0] += length;
    IndicatorCuboid3D<T> outflow(extend, origin);
    superGeometry.rename(1, 3, outflow);
  }
  {
    Vector<T,3>origin(eps, height - eps, -eps - dx);
    Vector<T,3>extend(length - eps * T(2.), eps * T(2.), width + T(2.) * eps + T(2.) * dx);
    IndicatorCuboid3D<T> topWall(extend, origin);
    superGeometry.rename(1, 4, topWall);

    origin[1] -= height;
    IndicatorCuboid3D<T> bottomWall(extend, origin);
    superGeometry.rename(1, 5, bottomWall);
  }
  // Clean gemetry
  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

template <typename CASE>
void prepareLattice(CASE& myCase){
  using T = CASE::value_t;
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;
  auto& superGeometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  auto& ALattice = myCase.getLattice(Concentration0{});
  auto& BLattice = myCase.getLattice(Concentration1{});
  auto& CLattice = myCase.getLattice(Concentration2{});

  ALattice.template setUnitConverter<AdeUnitConverter<T,MyCase::descriptor_t_of<Concentration0>>>(
    (T)   parameters.template get<parameters::PHYS_DELTA_X>(),                                                                               // physDeltaX
    (T)   util::pow(parameters.template get<parameters::HEIGHT>() / parameters.template get<parameters::RESOLUTION>(), 2),                            // physDeltaY
    (T)   parameters.template get<parameters::HEIGHT>(),                                                                                     // charPhysLength: reference length of simulation geometry
    (T)   parameters.template get<parameters::PECLET>() * parameters.template get<parameters::PHYS_DIFFUSIVITY>() / parameters.template get<parameters::HEIGHT>(),  // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   parameters.template get<parameters::PHYS_DIFFUSIVITY>(),                                                                                // physDiffusiviy
    (T)   parameters.template get<parameters::PHYS_CHAR_DENSITY>()                                                                           // physDensity: physical density in __kg / m^3__
  );
  ALattice.getUnitConverter().print();
  const auto& converter = ALattice.getUnitConverter();
  BLattice.setUnitConverter(converter);
  CLattice.setUnitConverter(converter);

  // define dynamics and boundary conditions
  const auto bulkIndicator = superGeometry.getMaterialIndicator({1, 2, 3, 4, 5});
  dynamics::set<SourcedAdvectionDiffusionBGKdynamics>(ALattice, bulkIndicator);
  dynamics::set<SourcedAdvectionDiffusionBGKdynamics>(BLattice, bulkIndicator);
  dynamics::set<SourcedAdvectionDiffusionBGKdynamics>(CLattice, bulkIndicator);

  // Component A
  boundary::set<boundary::AdvectionDiffusionDirichlet>(ALattice, superGeometry.getMaterialIndicator({2}));
  setZeroGradientBoundary<T, MyCase::descriptor_t_of<Concentration0>>(ALattice, superGeometry.getMaterialIndicator({3}));
  setZeroGradientBoundary<T, MyCase::descriptor_t_of<Concentration0>>(ALattice, superGeometry.getMaterialIndicator({4}));
  boundary::set<boundary::AdvectionDiffusionDirichlet>(ALattice, superGeometry.getMaterialIndicator({5}));

  // Component B
  boundary::set<boundary::AdvectionDiffusionDirichlet>(BLattice, superGeometry.getMaterialIndicator({2}));
  setZeroGradientBoundary<T, MyCase::descriptor_t_of<Concentration1>>(BLattice, superGeometry.getMaterialIndicator({3}));
  boundary::set<boundary::AdvectionDiffusionDirichlet>(BLattice, superGeometry.getMaterialIndicator({4}));
  setZeroGradientBoundary<T, MyCase::descriptor_t_of<Concentration1>>(BLattice, superGeometry.getMaterialIndicator({5}));

  // Component C
  boundary::set<boundary::AdvectionDiffusionDirichlet>(CLattice, superGeometry.getMaterialIndicator({2}));
  setZeroGradientBoundary<T, MyCase::descriptor_t_of<Concentration2>>(CLattice, superGeometry.getMaterialIndicator({3}));
  setZeroGradientBoundary<T, MyCase::descriptor_t_of<Concentration2>>(CLattice, superGeometry.getMaterialIndicator({4}));
  setZeroGradientBoundary<T, MyCase::descriptor_t_of<Concentration2>>(CLattice, superGeometry.getMaterialIndicator({5}));

  auto& coupling = myCase.setCouplingOperator(
    "HomogeneousReaction",
    ReactionCoupling{},
    Concentration0{}, ALattice,
    Concentration1{}, BLattice,
    Concentration2{}, CLattice);
  coupling.template setParameter<ReactionCoupling::LATTICE_REACTION_COEFF>(
    converter.getConversionFactorTime() * parameters.template get<parameters::REACTION_CONSTANT>());
  coupling.template setParameter<ReactionCoupling::STOCH_COEFF>(
    parameters.template get<parameters::STOICHIOMETRY>());
  coupling.template setParameter<ReactionCoupling::REACTION_ORDER>(
    parameters.template get<parameters::REACTION_ORDER>());

  clout << "Prepare Lattice ... OK" << std::endl;
}

template <typename CASE>
void setInitialValues(CASE& myCase) {
  using T = CASE::value_t;
  OstreamManager clout(std::cout, "setInitialValues");
  clout << "Setting initial values ..." << std::endl;

  auto& superGeometry = myCase.getGeometry();
  auto& ALattice = myCase.getLattice(Concentration0{});
  auto& BLattice = myCase.getLattice(Concentration1{});
  auto& CLattice = myCase.getLattice(Concentration2{});
  const auto& converter = ALattice.getUnitConverter();

  // set the convective velocity field
  Vector<T,3> u(ALattice.getUnitConverter().getCharPhysVelocity(), 0.0, 0.0);
  fields::setVelocity<descriptors::VELOCITY>(ALattice,
                                             superGeometry.getMaterialIndicator({1, 2, 3, 4, 5}),
                                             u);
  fields::setVelocity<descriptors::VELOCITY>(BLattice,
                                             superGeometry.getMaterialIndicator({1, 2, 3, 4, 5}),
                                             u);
  fields::setVelocity<descriptors::VELOCITY>(CLattice,
                                             superGeometry.getMaterialIndicator({1, 2, 3, 4, 5}),
                                             u);

  momenta::setDensity(ALattice, superGeometry.getMaterialIndicator({1,2,3,4}), 1e-8);
  momenta::setDensity(ALattice, superGeometry.getMaterialIndicator(5), 10.);
  momenta::setDensity(BLattice, superGeometry.getMaterialIndicator({1,2,3,5}), 1e-8);
  momenta::setDensity(BLattice, superGeometry.getMaterialIndicator(4), 10.);
  momenta::setDensity(CLattice, superGeometry.getMaterialIndicator({1,2,3,4,5}), 1e-8);

  ALattice.template setParameter<descriptors::OMEGA>(
    converter.getLatticeRelaxationFrequency());
  BLattice.template setParameter<descriptors::OMEGA>(
    converter.getLatticeRelaxationFrequency());
  CLattice.template setParameter<descriptors::OMEGA>(
    converter.getLatticeRelaxationFrequency());

  ALattice.initialize();
  BLattice.initialize();
  CLattice.initialize();

  clout << "Setting initial values OK" << std::endl;
}

// compute results and generate output
template <typename CASE, typename TIMER, typename CONVERGE>
void getResults(CASE& myCase, TIMER& timer, CONVERGE& converge, std::size_t iT)
{
  using T = CASE::value_t;
  auto& ALattice = myCase.getLattice(Concentration0{});
  auto& BLattice = myCase.getLattice(Concentration1{});
  auto& CLattice = myCase.getLattice(Concentration2{});
  auto& converter = ALattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  // generate vtm files for visualisation
  SuperVTMwriter3D<T> vtmWriter("convectedPlate");
  SuperLatticeDensity3D<T,MyCase::descriptor_t_of<Concentration0>> concentration0(ALattice);
  SuperLatticePhysField3D<T,MyCase::descriptor_t_of<Concentration0>,descriptors::VELOCITY> velocity(ALattice, converter.getConversionFactorVelocity());
  SuperLatticeDensity3D<T,MyCase::descriptor_t_of<Concentration1>> concentration1(BLattice);
  SuperLatticeDensity3D<T,MyCase::descriptor_t_of<Concentration2>> concentration2(CLattice);

  SuperLatticeField3D<T,MyCase::descriptor_t_of<Concentration0>,descriptors::SOURCE> productionA(ALattice);
  SuperLatticeField3D<T,MyCase::descriptor_t_of<Concentration1>,descriptors::SOURCE> productionB(BLattice);
  SuperLatticeField3D<T,MyCase::descriptor_t_of<Concentration2>,descriptors::SOURCE> productionC(CLattice);

  concentration0.getName() = "Concentration0";
  concentration1.getName() = "Concentration1";
  concentration2.getName() = "Concentration2";
  productionA.getName() = "Production0";
  productionB.getName() = "Production1";
  productionC.getName() = "Production2";
  vtmWriter.addFunctor(concentration0);
  vtmWriter.addFunctor(concentration1);
  vtmWriter.addFunctor(concentration2);
  vtmWriter.addFunctor(productionA);
  vtmWriter.addFunctor(productionB);
  vtmWriter.addFunctor(productionC);
  vtmWriter.addFunctor(velocity);

  const std::size_t vtkIter = converter.getLatticeTime(parameters.template get<parameters::PHYS_MAX_T>() / 100);
  const std::size_t statIter = converter.getLatticeTime(parameters.template get<parameters::PHYS_MAX_T>() / 100);

  if (iT == 0) {
    // write the geometry
    vtmWriter.createMasterFile();

    // write inital field
    vtmWriter.write(iT);
  }

  if (iT%statIter == 0) {
    ALattice.setProcessingContext(ProcessingContext::Evaluation);

    // Timer console output
    timer.update(iT);
    timer.printStep();

    // Lattice statistics console output
    ALattice.getStatistics().print(iT, converter.getPhysTime(iT));
    converge.takeValue(ALattice.getStatistics().getAverageRho(), false);
  }
  // write the vtk files
  if (iT % vtkIter == 0 && iT > 0) {
    vtmWriter.write(iT);
  }
}

template <typename CASE>
void simulate(CASE& myCase){
  OstreamManager clout(std::cout,"simulate");
  using T = CASE::value_t;
  auto& parameters = myCase.getParameters();
  auto& ALattice = myCase.getLattice(Concentration0{});
  auto& BLattice = myCase.getLattice(Concentration1{});
  auto& CLattice = myCase.getLattice(Concentration2{});
  auto& converter = ALattice.getUnitConverter();
  auto& superGeometry = myCase.getGeometry();

  util::Timer<T> timer(converter.getLatticeTime(parameters.template get<parameters::PHYS_MAX_T>()), superGeometry.getStatistics().getNvoxel());

  // Convergence check
  util::ValueTracer<T> converge(converter.getLatticeTime(parameters.template get<parameters::PHYS_CONVERGENCE_CHECK_T>()), parameters.template get<parameters::RESIDUUM>());
  timer.start();

  const std::size_t reactionUpdateT = converter.getLatticeTime(parameters.template get<parameters::REACTION_UPDATE_T>());
  for (std::size_t iT = 0; iT < converter.getLatticeTime(parameters.template get<parameters::PHYS_MAX_T>()); ++iT) {

    // Reaction coupling A + B -> C
    if (iT%reactionUpdateT == 0) {
      myCase.getOperator("HomogeneousReaction").apply();
    }

    // === Collide and stream ===
    ALattice.collideAndStream();
    BLattice.collideAndStream();
    CLattice.collideAndStream();

    // === Computation and output of the results ===
    getResults(myCase, timer, converge, iT);
    if (converge.hasConverged()) {
      clout << "Simulation stops." << std::endl;
      break;
    }
  }

  timer.stop();
  timer.printSummary();
}

void setInitialControl(MyOptiCase& optiCase) {
  // control[0] : reaction constant / 1e-6
  // control[1] : reaction order comp. A
  // control[2] : reaction order comp. B
  optiCase.getController().set({1.0, 1.0, 1.0});
}

template <typename CASE, typename CONTROLLER>
void applyControl(CASE& myCase, CONTROLLER& controller) {
  auto& parameters = myCase.getParameters();
  auto controls = controller.get();
  parameters.template set<parameters::REACTION_CONSTANT>(controls[0] * 1e-6);
  parameters.template set<parameters::REACTION_ORDER>({controls[1], controls[2], 0});
}

template <typename NAME>
auto computeL2Distance(MyOptiCase& optiCase, NAME) {
  auto& referenceCase = optiCase.getCase(Reference{});
  auto& referenceALattice = referenceCase.getLattice(Concentration0{});
  auto& myCase = optiCase.getCase(NAME{});
  auto& ALattice = myCase.getLattice(Concentration0{});
  auto& refConverter = referenceALattice.getUnitConverter();
  auto& converter = ALattice.getUnitConverter();
  auto refObjectiveDomain = referenceCase.getGeometry().getMaterialIndicator({1});
  auto objectiveDomain = myCase.getGeometry().getMaterialIndicator({1});

  auto concentrationAO = makeWriteFunctorO<functors::ConcentrationF, OBJECTIVE::Reference>(referenceALattice);
  concentrationAO->apply();
  copyFields<OBJECTIVE::Reference,OBJECTIVE::Reference>(referenceALattice, ALattice);

  auto objectiveO = makeWriteFunctorO<OBJECTIVE,opti::J>(ALattice);
  objectiveO->template setParameter<parameters::NORMALIZE>(computeL2Norm<OBJECTIVE::Reference>(referenceALattice, refObjectiveDomain, refConverter.getPhysDeltaX()));
  objectiveO->restrictTo(objectiveDomain);
  objectiveO->apply();

  return integrateField<opti::J>(ALattice, objectiveDomain, converter.getPhysDeltaX())[0];
}

template <typename NAME>
auto computeObjective(MyOptiCase& optiCase, NAME) {
  using CASE = MyOptiCase::case_t<NAME>;
  using T = CASE::value_t;
  auto& myCase = optiCase.getCase(NAME{});

  myCase.resetLattices();
  applyControl(myCase, optiCase.getController<T>());
  prepareLattice(myCase);
  setInitialValues(myCase);
  simulate(myCase);
  return computeL2Distance(optiCase, NAME{});
}

int main(int argc, char* argv[]) {
    /// === 1st Step: Initialization ===
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION       >(       40);
    myCaseParameters.set<HEIGHT           >(      160);
    myCaseParameters.set<BULK_VELOCITY    >(     0.05);
    myCaseParameters.set<PHYS_CHAR_DENSITY>(      1.0);
    myCaseParameters.set<PECLET           >(  108.005);
    myCaseParameters.set<PHYS_MAX_T       >( 100000.0);
    myCaseParameters.set<RESIDUUM         >(    1e-10);
    myCaseParameters.set<PHYS_DIFFUSIVITY >(  0.07407);
    myCaseParameters.set<PHYS_DELTA_X     >([](MyCase::ParametersD& p) {
      return p.get<HEIGHT>()/p.get<RESOLUTION>();
    });
    myCaseParameters.set<DOMAIN_EXTENT>([](MyCase::ParametersD& p) {
      return Vector{p.get<RESOLUTION>()*p.get<PHYS_DELTA_X>()*10, p.get<PHYS_DELTA_X>(), p.get<HEIGHT>()};
    });
    myCaseParameters.set<PHYS_CONVERGENCE_CHECK_T>([](MyCase::ParametersD& p) {
      return p.get<PHYS_MAX_T>()*0.001;
    });
    myCaseParameters.set<REACTION_CONSTANT    >(      1.1e-6);
    myCaseParameters.set<STOICHIOMETRY        >({-1.,-1.,1.});
    myCaseParameters.set<REACTION_ORDER       >({1.3,1.2,0.});
    myCaseParameters.set<REACTION_UPDATE_T    >(      10000.);
  }
  auto myADfParameters = myCaseParameters.copyAs<MyADfCase::value_t>();

  auto mesh = createMesh(myCaseParameters);
  auto aDfMesh = createMesh(myADfParameters);

  MyCase myCase(myCaseParameters, mesh);
  prepareGeometry(myCase);

  MyADfCase myADfCase(myADfParameters, aDfMesh);
  prepareGeometry(myADfCase);

  MyCase referenceCase(myCaseParameters, mesh);
  prepareGeometry(referenceCase);
  prepareLattice(referenceCase);
  setInitialValues(referenceCase);
  simulate(referenceCase);

  MyOptiCase optiCase;
  optiCase.setCase<Controlled>(myCase);
  optiCase.setCase<Derivatives>(myADfCase);
  optiCase.setCase<Reference>(referenceCase);

  setInitialControl(optiCase);

  optiCase.setObjective([](auto& optiCase, auto name) {
    return computeObjective(optiCase, name);
  });

  OptimizerParametersD<double> optiParameters;
  {
    using namespace parameters;
    optiParameters.set<optimizer::EPS>(1.e-8);
    optiParameters.set<optimizer::MAX_OPTIMIZATION_STEPS>(30);
    optiParameters.set<linesearch::INITIAL_STEP_SIZE>(1e-2);
    optiParameters.set<linesearch::MAX_STEP_ATTEMPTS>(30);
    optiParameters.set<linesearch::STEP_CONDITION>("Wolfe");
    optiParameters.set<optimizer::LOG_FILE_NAME>("log");
    optiParameters.set<optimizer::WITH_UPPER_BOUND>(true);
    optiParameters.set<optimizer::UPPER_BOUND>(5.0);
    optiParameters.set<optimizer::WITH_LOWER_BOUND>(true);
    optiParameters.set<optimizer::CONTROL_EPS>(0.);
  }

  OptimizerLBFGS<MyCase::value_t,std::vector<MyCase::value_t>> optimizer(
    optiCase.getController().size(), optiParameters,
    {OptimizerLogType::objective, OptimizerLogType::control, OptimizerLogType::derivative});

  /// @li Step F: Optimize
  optimizer.optimize(optiCase);
}
