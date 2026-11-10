/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2026 Johannes Grafen
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

/* Two-parameter force identification on the periodic plate with the adjoint LBM.
 * The analytical body force F* is scaled component-wise by two scalar controls,
 * F = alpha o F*, realized with the optional field opti::CONTROL_VECTOR
 * of the StationaryLinearElasticityMRT dynamics (default (1,1)). The cell-wise
 * sensitivities of the adjoint optimality condition are summed to the gradient of
 * the two controls.
 *
 *  The corresponding forward problem is the example
 *  `structuralmechanics/linearElasticity/quasistationary/periodicPlate`.
 */

/**
 * Grafen et al. "Benchmark cases for gradient-based optimization in linear elastic solid
 * mechanics using lattice Boltzmann methods" arXiv preprint arXiv:2609.31347 (2026).
*/

#pragma once
#include <olb.h>
#include "analyticalSolution.h"

using namespace olb;
using namespace olb::names;
using namespace olb::opti;

namespace olb::parameters {
struct CONVERGENCE_PRECISION_ADJOINT : public descriptors::FIELD_BASE<1> {};
struct INITIAL_CONTROL : public descriptors::FIELD_BASE<0,1> {};
struct CONVERGED_IT : public descriptors::TYPED_FIELD_BASE<std::size_t, 1> {};
struct ANALYTICAL_REFERENCE : public descriptors::TYPED_FIELD_BASE<bool, 1> {};
struct SET_ADJOINT_MAX_IT : public descriptors::TYPED_FIELD_BASE<bool, 1> {};
struct CONV_VALUE : public descriptors::FIELD_BASE<1> {};
} // namespace olb::parameters

// === Step 1: Declarations ===
using MyCase = Case<NavierCauchy,
    Lattice<
      double,
      descriptors::D2Q8<
        descriptors::FORCE,
        structuralmechanics::MOMENT_VECTOR,
        structuralmechanics::DISPLACEMENT,
        structuralmechanics::CAUCHY_STRESS
      >
    >
  >;

using PeriodicPlate = SteadyOptimalitySystem<
  functors::L2DistanceF<functors::DisplacementF>,                  //Objective function
  structuralmechanics::StationaryLinearElasticityMRT
    <MyCase::value_t, MyCase::descriptor_t>,                       //Primal Dynamics
  opti::CONTROL_VECTOR                                             //Controlled Field
>;

using MyOptiCase = OptiCaseAdjoint<
  Controlled, MyCase,
  Adjoint, MyCase,
  Reference, MyCase,
  Optimality, PeriodicPlate
>;

/// @brief Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
auto createMesh(MyCase::ParametersD& parameters)
{
  using T = MyCase::value_t;

  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector origin = parameters.get<parameters::ORIGIN>();
  const T      dx     = parameters.get<parameters::PHYS_DELTA_X>();

  IndicatorCuboid2D<T> cuboid(extent - Vector {dx, dx}, origin);

  Mesh<T, MyCase::d> mesh(cuboid, dx, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({true, true});
  return mesh;
}

/// @brief Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
/// @note The material numbers will be used to assign physics to lattice nodes
void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  auto& geometry = myCase.getGeometry();

  /// Set material numbers
  geometry.rename(0, 1);
  geometry.clean();
  geometry.checkForErrors();
  geometry.getStatistics().print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

/// @brief Set lattice dynamics
/// @param myCase The Case instance which keeps the simulation data
void prepareLattice(MyCase& myCase)
{
  using DESCRIPTOR = MyCase::descriptor_t;
  using T          = MyCase::value_t;
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(NavierCauchy {});
  auto& parameters = myCase.getParameters();

  dynamics::set<structuralmechanics::StationaryLinearElasticityMRT>(lattice, geometry.getMaterialIndicator({1}));

  const T physCharLength       = parameters.get<parameters::PHYS_CHAR_LENGTH>();
  const T dx                   = parameters.get<parameters::PHYS_DELTA_X>();
  const T physCharDisplacement = parameters.get<parameters::PHYS_CHAR_DISPLACEMENT>();
  const T physDensity          = parameters.get<parameters::PHYS_DENSITY>();
  const T physYoungsModulus    = parameters.get<parameters::YOUNGS_MODULUS>();
  const T physPoissonRatio     = parameters.get<parameters::POISSON_RATIO>();
  const T kappa                = parameters.get<parameters::DAMPING_FACTOR>();
  const std::string solidCase  = parameters.template get<parameters::SOLID_CASE>();

  const T dt                   = parameters.get<parameters::PHYS_DELTA_T>();

  // Set up a unit converter with the characteristic physical units
  lattice.template setUnitConverter<StationaryMRTUnitConverter<T, DESCRIPTOR>>(
    dx,
    dt,
    physCharLength,
    T(1),
    physCharDisplacement,
    physYoungsModulus,
    physPoissonRatio,
    physDensity,
    kappa,
    solidCase
  );
  const auto& converter = lattice.getUnitConverter();
  converter.print();

  const T theta = T(1) / descriptors::invCs2<T, DESCRIPTOR>();

  const T omega_11  = T(1.) / (converter.getLatticeLameMu() / theta  + T(0.5));
  const T omega_d   = T(1.) / (T(2.) *  converter.getLatticeLameMu()                                   / (T(1.) - theta) + T(0.5));
  const T omega_s   = T(1.) / (T(2.) * (converter.getLatticeLameMu() + converter.getLatticeLameLambda()) / (T(1.) + theta) + T(0.5));

  const T tau_12 = T(0.5);
  const T tau_22 = T(0.5);
  const T tau_21 = tau_12;

  const T omega_12  = T(1.) / (tau_12 + T(0.5));
  const T omega_21  = T(1.) / (tau_21 + T(0.5));
  const T omega_22  = T(1.) / (tau_22 + T(0.5));

  Vector<T,DESCRIPTOR::q> Omega =
  {
    0,
    0,
    omega_11,
    omega_s,
    omega_d,
    omega_12,
    omega_21,
    omega_22
  };
  lattice.template setParameter<structuralmechanics::OMEGA_SOLID>(Omega);
  parameters.template set<parameters::OMEGA_SOLID>(Omega);

  clout << "Prepare Lattice ... OK" << std::endl;
}

/// Set initial condition for primal variables (velocity and density)
/// @param myCase The Case instance which keeps the simulation data
/// @note Be careful: initial values have to be set using lattice units
void setInitialValues(MyCase& myCase)
{
  using T = MyCase::value_t;
  OstreamManager clout(std::cout, "setInitialValues");
  clout << "Set Initial Values ..." << std::endl;

  auto& lattice   = myCase.getLattice(names::NavierCauchy {});
  auto& converter = lattice.getUnitConverter();
  auto& geometry  = myCase.getGeometry();

  ForceField2D<T> force(converter.getLameMu(), converter.getLameLambda(), converter.getConversionFactorSolidForce());
  fields::set<descriptors::FORCE>(lattice, geometry.getMaterialIndicator({1}), force);

  lattice.initialize();

  clout << "Set Initial Values ... OK" << std::endl;
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
/// @param iT The time step
void getResults(MyCase& myCase, util::Timer<typename MyCase::value_t>& timer, std::size_t iT)
{
  OstreamManager clout(std::cout, "getResults");

  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;

  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(names::NavierCauchy {});
  auto& converter  = lattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  const std::size_t iTMax     = converter.getLatticeTime(parameters.template get<parameters::PHYS_MAX_T>());
  const std::size_t iTlog     = parameters.template get<parameters::LATTICE_LOG_ITER_T>();
  const std::size_t convIter  = parameters.template get<parameters::CONV_ITER>();
  const bool        converged = parameters.template get<parameters::CONVERGED>();

  if (iT % iTlog == 0 || iT % iTMax == 0 || converged) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    ManufacturedSolutionU2D<T>      anaPhysDisp(true);
    ManufacturedSolutionStress2D<T> anaPhysStress(converter.getLameMu(), converter.getLameLambda());

    const T charU = converter.getCharPhysDisplacement();
    const T dx    = converter.getPhysDeltaX();
    const T dt    = converter.getPhysDeltaT();
    const T kappa = converter.getDampingFactor();

    SuperLatticePhysField2D<T, DESCRIPTOR, olb::structuralmechanics::DISPLACEMENT> numPhysDisp(lattice, dx);
    SuperLatticePhysField2D<T, DESCRIPTOR, olb::structuralmechanics::CAUCHY_STRESS> numPhysStress(lattice, kappa * dx * charU / dt);

    T   l2UResult[2]        = {T(), T()};
    T   lInfUResult[2]      = {T(), T()};
    T   l2StressResult[2]   = {T(), T()};
    T   lInfStressResult[2] = {T(), T()};
    int tmp[]               = {int()};

    auto                            indicatorF = geometry.getMaterialIndicator({1});
    SuperRelativeErrorL2Norm2D<T>   relUErrorL2Norm(lattice, numPhysDisp, anaPhysDisp, indicatorF);
    SuperRelativeErrorLinfNorm2D<T> relUErrorLinfNorm(lattice, numPhysDisp, anaPhysDisp, indicatorF);
    SuperRelativeErrorL2Norm2D<T>   relStressErrorL2Norm(lattice, numPhysStress, anaPhysStress, indicatorF);
    SuperRelativeErrorLinfNorm2D<T> relStressErrorLinfNorm(lattice, numPhysStress, anaPhysStress, indicatorF);

    relUErrorL2Norm(l2UResult, tmp);
    relUErrorLinfNorm(lInfUResult, tmp);
    relStressErrorL2Norm(l2StressResult, tmp);
    relStressErrorLinfNorm(lInfStressResult, tmp);

    clout << "N\t" << "L2 U Error\t" << "LInf U Error\t" << "L2 Stress Error\t" << "LInf Stress Error" << std::endl;
    clout << converter.getResolution() << "\t" << l2UResult[0] << "\t" << lInfUResult[0] << "\t" << l2StressResult[0]
          << "\t\t" << lInfStressResult[0] << std::endl;

    timer.print(iT);

    lattice.setProcessingContext(ProcessingContext::Simulation);
  }

  if (iT % convIter == 0) {
    const T avgDisplacementSqr = lattice.getStatistics().getAverageEnergy();
    parameters.template set<parameters::CONV_VALUE>(avgDisplacementSqr);
  }

  if (converged) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    auto indicatorF = geometry.getMaterialIndicator({1});

    T   l2UResult[2]        = {T(), T()};
    T   lInfUResult[2]      = {T(), T()};
    T   l2StressResult[2]   = {T(), T()};
    T   lInfStressResult[2] = {T(), T()};
    int tmp[]               = {int()};

    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::DISPLACEMENT> dispNumPhys(lattice, converter.getConversionFactorDisplacement());
    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::CAUCHY_STRESS      > stressNumPhys(lattice, converter.getConversionFactorSolidStress());

    ManufacturedSolutionU2D<T>                    dispAnaPhys(true);

    ManufacturedSolutionStress2D<T>               stressAnaPhys(converter.getLameMu(), converter.getLameLambda());

    SuperRelativeErrorL2Norm2D<T>   relUErrorL2Norm(lattice, dispNumPhys, dispAnaPhys, indicatorF);
    SuperRelativeErrorLinfNorm2D<T> relUErrorLinfNorm(lattice, dispNumPhys, dispAnaPhys, indicatorF);
    relUErrorL2Norm(l2UResult, tmp);
    relUErrorLinfNorm(lInfUResult, tmp);

    SuperRelativeErrorL2Norm2D<T>   relStressErrorL2Norm(lattice, stressNumPhys, stressAnaPhys, indicatorF);
    SuperRelativeErrorLinfNorm2D<T> relStressErrorLinfNorm(lattice, stressNumPhys, stressAnaPhys, indicatorF);
    relStressErrorL2Norm(l2StressResult, tmp);
    relStressErrorLinfNorm(lInfStressResult, tmp);

    CSV<T> csvWriterU("disp_err", ',', {"N", "L2", "LInf"}, ".csv", true);
    CSV<T> csvWriterSigma("stress_err", ',', {"N", "L2", "LInf"}, ".csv", true);

    csvWriterU.writeDataFile(converter.getResolution(), {l2UResult[0], lInfUResult[0]});
    csvWriterSigma.writeDataFile(converter.getResolution(), {l2StressResult[0], lInfStressResult[0]});
  }
}

void getOptiResults(MyOptiCase& optiCase)
{
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  using ControlledField = PeriodicPlate::controls_t;
  OstreamManager clout(std::cout, "getResults");

  auto& lattice    = optiCase.getCase(Controlled {}).getLattice(NavierCauchy {});
  auto& adjointLattice    = optiCase.getCase(Adjoint {}).getLattice(NavierCauchy {});
  auto& converter  = lattice.getUnitConverter();

  const std::size_t                             iT = optiCase.getOptimizationStep();
  SuperVTMwriter2D<T>                           vtmWriter("periodicPlate");
  ManufacturedSolutionU2D<T>                    anaPhysDisp(true);
  SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> anaPhysDispOnLattice(anaPhysDisp, lattice);

  ManufacturedSolutionStress2D<T> anaPhysStress(converter.getLameMu(), converter.getLameLambda());
  SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> anaPhysStressOnLattice(anaPhysStress, lattice);

  const T charU = converter.getCharPhysDisplacement();
  const T dx    = converter.getPhysDeltaX();
  const T dt    = converter.getPhysDeltaT();
  const T kappa = converter.getDampingFactor();

  SuperLatticeField2D<T, DESCRIPTOR, olb::structuralmechanics::DISPLACEMENT>     numLatticeDisp(lattice);
  SuperLatticeField2D<T, DESCRIPTOR, olb::structuralmechanics::CAUCHY_STRESS>           numLatticeStress(lattice);
  SuperLatticePhysField2D<T, DESCRIPTOR, olb::structuralmechanics::DISPLACEMENT> numPhysDisp(lattice, dx);
  SuperLatticePhysField2D<T, DESCRIPTOR, olb::structuralmechanics::CAUCHY_STRESS> numPhysStress(lattice, kappa * dx * charU / dt);
  SuperLatticeField<T,DESCRIPTOR, ControlledField> controlOnLattice(lattice);
  SuperLatticeField<T,DESCRIPTOR, opti::J> objective(lattice);

  // opti::J is not part of the descriptor, so the Full communicator does not include it by default
  auto& controlledCommunicator = lattice.getCommunicator(stage::Full());
  controlledCommunicator.requestField<opti::J>();
  controlledCommunicator.requestOverlap(1);
  controlledCommunicator.exchangeRequests();
  controlledCommunicator.communicate();

  auto& communicator = adjointLattice.getCommunicator(stage::Full());
  communicator.requestField<opti::DJDALPHA<ControlledField>>();
  communicator.requestField<opti::DCDALPHA<ControlledField>>();
  communicator.requestField<opti::SENSITIVITY<ControlledField>>();
  communicator.requestField<opti::DJDF>();
  communicator.requestOverlap(1);
  communicator.exchangeRequests();
  communicator.communicate();

  SuperLatticeField<T,DESCRIPTOR, opti::SENSITIVITY<ControlledField>> Sensitivity(adjointLattice);
  SuperLatticeField<T,DESCRIPTOR, opti::DJDALPHA<ControlledField>> gradientObjective(adjointLattice);
  SuperLatticeField<T,DESCRIPTOR, opti::DCDALPHA<ControlledField>> jacobianCollision(adjointLattice);
  SuperLatticeField<T,DESCRIPTOR, opti::DJDF> djdf(adjointLattice);

  ForceField2D<T> anaPhysForce(converter.getLameMu(), converter.getLameLambda(), 1.);
  SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> anaPhysForceOnLattice(anaPhysForce, lattice);

  vtmWriter.addFunctor(numLatticeDisp, "num lattice disp");
  vtmWriter.addFunctor(numPhysDisp, "num phys disp");
  vtmWriter.addFunctor(anaPhysDispOnLattice, "ana phys disp");

  vtmWriter.addFunctor(numLatticeStress, "num lattice stress");
  vtmWriter.addFunctor(numPhysStress, "num phys stress");
  vtmWriter.addFunctor(anaPhysStressOnLattice, "ana phys stress");
  vtmWriter.addFunctor(anaPhysForceOnLattice, "ana phys force");

  vtmWriter.addFunctor(controlOnLattice, "control");
  vtmWriter.addFunctor(objective, "objective");
  vtmWriter.addFunctor(Sensitivity, "sensitivity");
  vtmWriter.addFunctor(gradientObjective, "DJDALPHA");
  vtmWriter.addFunctor(jacobianCollision, "DCDalpha");
  vtmWriter.addFunctor(djdf, "DJDF");

  if (iT == 0) {
    vtmWriter.createMasterFile();
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid(lattice);
    SuperLatticeRank2D<T, DESCRIPTOR>   rank(lattice);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.write(iT);
  }

  if (iT > 0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    adjointLattice.setProcessingContext(ProcessingContext::Evaluation);
    clout << "Writing regular results to VTK ..." << std::endl;
    vtmWriter.write(iT);
  }
}

/// @brief Execute simulation: set initial values and run time loop
/// @param myCase The Case instance which keeps the simulation data
std::size_t simulate(MyCase& myCase, MyCase::value_t const convergenceEPS)
{
  using T = MyCase::value_t;
  OstreamManager clout(std::cout, "simulate");
  clout << "Starting simulation ..." << std::endl;

  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(NavierCauchy {});
  auto& parameters = myCase.getParameters();

  const T           physMaxT = parameters.get<parameters::PHYS_MAX_T>();
  const std::size_t convIt   = parameters.get<parameters::CONV_ITER>();
  const std::size_t iTmax    = lattice.getUnitConverter().getLatticeTime(physMaxT);

  parameters.template set<parameters::CONVERGED>(false);
  util::ValueTracer<T> converge(convIt, convergenceEPS);
  std::size_t          convergedIt = iTmax;
  std::size_t          zeroSamples = 0;

  util::Timer<T> timer(iTmax, geometry.getStatistics().getNvoxel());
  timer.start();
  bool convCondition2 = parameters.get<parameters::SET_ADJOINT_MAX_IT>() ? (!myCase.hasName(Adjoint {})) : true;
  for (std::size_t iT = 0; iT < iTmax; ++iT) {

    if (parameters.get<parameters::SET_ADJOINT_MAX_IT>()
        && myCase.hasName(Adjoint {})
        && iT == myCase.getParameters().template get<parameters::CONVERGED_IT>()) {
      clout << "Adjoint simulation reached the number of iterations of the controlled case: " << iT << std::endl;
      break;
    }

    /// Collide and Stream Execution ===
    lattice.collideAndStream();

    /// Computation and Output of the Results ===
    getResults(myCase, timer, iT);

    if (convCondition2 && iT % convIt == 0 ) {
      const T energy = lattice.getStatistics().getAverageEnergy();
      converge.takeValue(energy, true);
      // Test for zero average displacement, to prevent ValueTracer from throwing "Simulation diverged" when the state is stationary at zero
      zeroSamples = (util::abs(energy) < std::numeric_limits<T>::min()) ? zeroSamples + 1 : 0;
    }

    // A state that stays zero (e.g. zero control) is stationary
    const bool zeroState = zeroSamples >= convIt;
    bool convCondition = parameters.get<parameters::SET_ADJOINT_MAX_IT>()
                       ? (!myCase.hasName(Adjoint {}) && (zeroState || converge.hasConverged()))
                       : (zeroState || converge.hasConverged());
    if (convCondition) {
      parameters.set<parameters::CONVERGED>(true);
      getResults(myCase, timer, iT);
      clout << "Convergence reached at time step" << iT << std::endl;
      convergedIt = iT;
      break;
    }
  }

  timer.stop();
  timer.printSummary();

  clout << "Simulation finished ..." << std::endl;
  return convergedIt;
}

void prepareAdjointLattice(MyOptiCase& optiCase)
{
  OstreamManager clout(std::cout, "pepareAdjointLattice");
  clout << "Prepare Adjoint Lattice ..." << std::endl;
  auto& adjointLattice    = optiCase.getCase(Adjoint {}).getLattice(NavierCauchy {});
  auto& controlledLattice = optiCase.getCase(Controlled {}).getLattice(NavierCauchy {});
  auto& geometry          = optiCase.getCase(Adjoint {}).getGeometry();
  auto& controlparameters = optiCase.getCase(Controlled {}).getParameters();

  // Set adjoint unit converter as the controlled case
  adjointLattice.setUnitConverter(controlledLattice.getUnitConverter());

  // Define dual physics
  dynamics::set<olb::Dual<PeriodicPlate::dynamics_t>>(adjointLattice, geometry.getMaterialIndicator({1}));

  adjointLattice.setParameter<structuralmechanics::OMEGA_SOLID>(controlparameters.template get<parameters::OMEGA_SOLID>());

  adjointLattice.revertStreaming();

  clout << "Prepare Adjoint Lattice ... OK" << std::endl;
}

void setAdjointInitialValues(MyOptiCase& optiCase)
{
  using T          = MyCase::value_t;

  OstreamManager clout(std::cout, "setAdjointInitialValues");
  clout << "Set Adjoint Initial Values ..." << std::endl;

  auto& adjointLattice    = optiCase.getCase(Adjoint {}).getLattice(NavierCauchy {});
  auto& controlledLattice = optiCase.getCase(Controlled {}).getLattice(NavierCauchy {});
  auto& referenceLattice  = optiCase.getCase(Reference {}).getLattice(NavierCauchy {});
  auto& converter         = adjointLattice.getUnitConverter();
  auto& geometry          = optiCase.getCase(Adjoint {}).getGeometry();
  auto  bulkIndicator     = geometry.getMaterialIndicator({1});
  auto& parameters        = optiCase.getCase(Adjoint {}).getParameters();
  auto& controlparameters = optiCase.getCase(Controlled {}).getParameters();
  auto& optimalitySystem  = optiCase.get(Optimality{});

  adjointLattice.initialize();

  if(parameters.get<parameters::ANALYTICAL_REFERENCE>()) {
    // For comparison with analytical solution, compute displacement field from functor and store in reference lattice
    ManufacturedSolutionU2D<T> AnalyticalDisplacement(true);
    fields::set<structuralmechanics::DISPLACEMENT>(referenceLattice,geometry.getMaterialIndicator({1}), AnalyticalDisplacement);
    referenceLattice.template setProcessingContext<Array<structuralmechanics::DISPLACEMENT>>(ProcessingContext::Simulation);
  } else {
    // For comparison with discrete solution, compute normalization factor for objective (derivative) using reference simulation
    auto displacementO = makeWriteFunctorO<functors::DisplacementF, structuralmechanics::DISPLACEMENT>(referenceLattice);
    displacementO->restrictTo(bulkIndicator);
    displacementO->template setParameter<parameters::CONVERSION_DISPLACEMENT>(converter.getCharPhysDisplacement());
    displacementO->template setParameter<parameters::OMEGA_SOLID>(controlparameters.template get<parameters::OMEGA_SOLID>());
    displacementO->apply();
  }

  //Compute DJDF
  optimalitySystem.initializeAdjointProblem(adjointLattice, controlledLattice, bulkIndicator);
}

/// @brief Set the two scalar controls, which scale the components of the analytical force
void setInitialControls(MyOptiCase& optiCase)
{
  using T               = MyCase::value_t;
  using ControlledField = PeriodicPlate::controls_t;

  OstreamManager clout(std::cout, "setInitialControls");
  clout << "Set Initial Control ..." << std::endl;
  auto& controller = optiCase.getController();
  auto& parameters = optiCase.getCase(Controlled {}).getParameters();
  auto& geometry   = optiCase.getCase(Controlled {}).getGeometry();
  auto& lattice    = optiCase.getCase(Controlled {}).getLattice(NavierCauchy {});

  controller.setProjection<projection::Identity<T>>();
  controller.set<ControlledField>(geometry, 1, lattice);

  const auto initialControl = parameters.get<parameters::INITIAL_CONTROL>();
  controller.set(std::vector<T>{initialControl[0], initialControl[1]});
}

/// @brief Write the two controls into the CONTROL_VECTOR field of every bulk cell
void applyControls(Controller<MyCase::value_t>& controller, MyCase& controlledCase)
{
  using T               = MyCase::value_t;
  using ControlledField = PeriodicPlate::controls_t;

  OstreamManager clout(std::cout, "applyControls");
  clout << "Applying control ..." << std::endl;

  auto& lattice  = controlledCase.getLattice(NavierCauchy {});
  auto& geometry = controlledCase.getGeometry();

  AnalyticalConst2D<T, T> control(controller.get()[0], controller.get()[1]);
  fields::set<ControlledField>(lattice, geometry.getMaterialIndicator({1}), control);
  lattice.template setProcessingContext<Array<ControlledField>>(ProcessingContext::Simulation);
}

MyCase::value_t objectiveF(MyOptiCase& optiCase)
{
  using ObjectiveF = PeriodicPlate::objective_t;
  auto& controlledCase    = optiCase.getCase(Controlled {});
  auto& controlledLattice = controlledCase.getLattice(NavierCauchy {});
  auto& controlledParameters = controlledCase.getParameters();
  auto  objectiveDomain   = controlledCase.getGeometry().getMaterialIndicator({1});
  auto& converter         = controlledLattice.getUnitConverter();
  auto& referenceLattice   = optiCase.getCase(Reference {}).getLattice(NavierCauchy {});
  auto  refObjectiveDomain = optiCase.getCase(Reference {}).getGeometry().getMaterialIndicator({1});
  auto& parameters         = controlledCase.getParameters();

  auto& optimalitySystem   = optiCase.get(Optimality {});

  // Write displacement field in the reference case
  if(parameters.get<parameters::ANALYTICAL_REFERENCE>()) {
    ManufacturedSolutionU2D<MyCase::value_t> anaPhysDisp(true);
    fields::set<ObjectiveF::Reference>(referenceLattice, refObjectiveDomain, anaPhysDisp);
    referenceLattice.template setProcessingContext<Array<ObjectiveF::Reference>>(ProcessingContext::Simulation);
  } else {
    auto displacementO = makeWriteFunctorO<functors::DisplacementF, ObjectiveF::Reference>(referenceLattice);
    displacementO->restrictTo(refObjectiveDomain);
    displacementO->template setParameter<parameters::CONVERSION_DISPLACEMENT>(converter.getCharPhysDisplacement());
    displacementO->template setParameter<parameters::OMEGA_SOLID>(controlledParameters.template get<parameters::OMEGA_SOLID>());
    displacementO->apply();
  }

  copyFields<ObjectiveF::Reference, ObjectiveF::Reference>(referenceLattice, controlledLattice);

  optimalitySystem.setObjectiveParameters(controlledLattice,
                                          parameters::CONVERSION_DISPLACEMENT{},
                                          converter.getCharPhysDisplacement(),
                                          parameters::OMEGA_SOLID{},
                                          controlledParameters.get<parameters::OMEGA_SOLID>(),
                                          parameters::NORMALIZE{},
                                          computeL2Norm<ObjectiveF::Reference>(referenceLattice, refObjectiveDomain, converter.getPhysDeltaX())
  );

  return optimalitySystem.evaluateObjectiveF(controlledLattice, objectiveDomain);
}

/// @brief Evaluate the optimality condition and sum the cell-wise sensitivities
///        to the derivatives with respect to the two scalar controls
std::vector<MyCase::value_t> derivativeF(MyOptiCase& optiCase)
{
  using T               = MyCase::value_t;
  using ControlledField = PeriodicPlate::controls_t;

  OstreamManager clout(std::cout, "derivativeF");
  auto& adjointLattice       = optiCase.getCase(Adjoint {}).getLattice(NavierCauchy {});
  auto& controlledLattice    = optiCase.getCase(Controlled {}).getLattice(NavierCauchy {});
  auto  objectiveDomain      = optiCase.getCase(Controlled {}).getGeometry().getMaterialIndicator({1});
  auto& optimalitySystem     = optiCase.get(Optimality{});
  auto& controlledParameters = optiCase.getCase(Controlled {}).getParameters();

  AnalyticalConst2D<T, T> dProjectionDalpha(T(1), T(1));
  fields::set<opti::DPROJECTIONDALPHA<ControlledField>>(adjointLattice,
    optiCase.getCase(Adjoint {}).getGeometry().getMaterialIndicator({1}), dProjectionDalpha);
  adjointLattice.template setProcessingContext<Array<opti::DPROJECTIONDALPHA<ControlledField>>>(
      ProcessingContext::Simulation);

  optimalitySystem.setPrimalCollisionParameters(controlledLattice,
                                                structuralmechanics::OMEGA_SOLID{},
                                                controlledParameters.template get<parameters::OMEGA_SOLID>()
  );

  auto sensitivities = optimalitySystem.evaluateOptimalityCondition(adjointLattice, controlledLattice, objectiveDomain);
  getOptiResults(optiCase);

  std::vector<T> sum(MyCase::d, T(0));
  for (std::size_t i = 0; i < sensitivities.size(); i += 2) {
    sum[0] += sensitivities[i];
    sum[1] += sensitivities[i + 1];
  }

  return sum;
}

MyCase::value_t computeObjective(MyOptiCase& optiCase)
{
  OstreamManager clout(std::cout, "computeObjective");
  clout << "Computing Objective ..." << std::endl;

  auto& controlledCase = optiCase.getCase(Controlled {});
  auto& parameters     = controlledCase.getParameters();

  controlledCase.resetLattices();
  prepareLattice(controlledCase);
  setInitialValues(controlledCase);
  applyControls(optiCase.getController(), controlledCase);

  std::size_t convergedIt = simulate(controlledCase, parameters.get<parameters::CONVERGENCE_PRECISION>());
  parameters.set<parameters::CONVERGED_IT>(convergedIt);

  return objectiveF(optiCase);
}

std::vector<MyCase::value_t> computeDerivative(MyOptiCase& optiCase)
{
  OstreamManager clout(std::cout, "computeDerivative");
  clout << "Computing Derivative ..." << std::endl;
  auto& adjointCase = optiCase.getCase(Adjoint {});
  auto& parameters   = adjointCase.getParameters();

  adjointCase.resetLattices();
  prepareAdjointLattice(optiCase);
  setAdjointInitialValues(optiCase);
  simulate(adjointCase, parameters.get<parameters::CONVERGENCE_PRECISION_ADJOINT>());

  return derivativeF(optiCase);
}
