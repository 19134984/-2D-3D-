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

/**
 *  The corresponding forward problem is the example
 *  `structuralmechanics/linearElasticity/quasistationary/periodicPlate`.
 *
 * Grafen et al. "Benchmark cases for gradient-based optimization in linear elastic solid
 * mechanics using lattice Boltzmann methods" arXiv preprint arXiv:2609.31347 (2026).
*/

#pragma once
// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include <olb.h>

#include "analyticalSolution.h"

using namespace olb;
using namespace olb::names;
using namespace olb::opti;

namespace olb::parameters {
  struct CONV_VALUE : public fields::SCALAR<10.> {};
  struct INITIAL_CONTROL : public descriptors::FIELD_BASE<0,1> {};
  struct ANALYTICAL_REFERENCE : public fields::TYPED_SCALAR<bool, false> {};
} // namespace olb::parameters

//2D ADf type
using ADf = util::ADf<double, 2>;

// === Step 1: Declarations ===
using MyCase = Case<names::NavierCauchy,
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
//using MyADfCase = Case<names::NavierCauchy, Lattice<ADf, descriptors::D2Q8<>>>;
using MyADfCase = MyCase::exchange_value_t<ADf>;
using MyOptiCase = OptiCaseADf<Controlled, MyCase, Derivatives, MyADfCase, Reference, MyCase>;

/// @brief Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
template <typename PARAMETERS>
auto createMesh(PARAMETERS& parameters)
{
  using T          = PARAMETERS::value_t;

  const Vector extent = parameters.template get<parameters::DOMAIN_EXTENT>();
  const Vector origin = parameters.template get<parameters::ORIGIN>();
  const T      dx     = parameters.template get<parameters::PHYS_DELTA_X>();

  IndicatorCuboid2D<T> cuboid(extent - Vector<T, MyCase::d> {dx, dx}, origin);

  Mesh<T, MyCase::d> mesh(cuboid, dx, singleton::mpi().getSize());
  mesh.setOverlap(parameters.template get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({true, true});
  return mesh;
}

/// @brief Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
/// @note The material numbers will be used to assign physics to lattice nodes
template <typename CASE>
void prepareGeometry(CASE& myCase)
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
template <typename CASE>
void prepareLattice(CASE& myCase)
{
  using T          = CASE::value_t;
  using DESCRIPTOR = CASE::descriptor_t;
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(NavierCauchy {});
  auto& parameters = myCase.getParameters();

  dynamics::set<structuralmechanics::StationaryLinearElasticityMRT>(lattice, geometry.getMaterialIndicator({1}));

  const T physCharLength        = parameters.template get<parameters::PHYS_CHAR_LENGTH       >();
  const T dx                    = parameters.template get<parameters::PHYS_DELTA_X           >();
  const T physCharDisplacement  = parameters.template get<parameters::PHYS_CHAR_DISPLACEMENT >();
  const T physDensity           = parameters.template get<parameters::PHYS_DENSITY           >();
  const T physYoungsModulus     = parameters.template get<parameters::YOUNGS_MODULUS         >();
  const T physPoissonRatio      = parameters.template get<parameters::POISSON_RATIO          >();
  const T kappa                 = parameters.template get<parameters::DAMPING_FACTOR         >();
  const std::string solidCase   = parameters.template get<parameters::SOLID_CASE             >();
  const T dt                    = parameters.template get<parameters::PHYS_DELTA_T>();

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
template <typename CASE>
void setInitialValues(CASE& myCase)
{
  using T = CASE::value_t;
  OstreamManager clout(std::cout, "setInitialValues");
  clout << "Set Initial Values ..." << std::endl;

  auto& lattice    = myCase.getLattice(names::NavierCauchy {});
  auto& geometry   = myCase.getGeometry();
  auto& converter  = lattice.getUnitConverter();

  ForceField2D<T> force(converter.getLameMu(), converter.getLameLambda(), converter.getConversionFactorSolidForce());
  fields::set<descriptors::FORCE>(lattice, geometry.getMaterialIndicator({1}), force);

  lattice.initialize();

  clout << "Set Initial Values ... OK" << std::endl;
}

/// Update boundary values at times (and external fields, if they exist)
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
/// @note Be careful: boundary values have to be set using lattice units
template <typename CASE>
void setTemporalValues(CASE& myCase, std::size_t iT)
{
  // Nothing to do here, because simulation does not depend on time
}

/// Compute simulation results at times
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
template <typename CASE>
void getResults(CASE& myCase, util::Timer<typename CASE::value_t>& timer, std::size_t iT)
{
  OstreamManager clout(std::cout, "getResults");

  using T          = CASE::value_t;
  using DESCRIPTOR = CASE::descriptor_t;

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

template <typename NAME>
void getOptiResults(MyOptiCase& optiCase)
{
  using CASE       = MyOptiCase::case_t<NAME>;
  using T          = CASE::value_t;
  using DESCRIPTOR = CASE::descriptor_t;

  OstreamManager clout(std::cout, "getResults");
  auto& lattice    = optiCase.getCase(NAME {}).getLattice(NavierCauchy {});
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
  SuperLatticePhysField2D<T, DESCRIPTOR, olb::structuralmechanics::CAUCHY_STRESS> numPhysStress(lattice, kappa * dx * dt / charU);

  ForceField2D<T> anaPhysForce(converter.getLameMu(), converter.getLameLambda(), 1.);
  SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> anaPhysForceOnLattice(anaPhysForce, lattice);

  vtmWriter.addFunctor(numLatticeDisp, "num lattice disp");
  vtmWriter.addFunctor(numPhysDisp, "num phys disp");
  vtmWriter.addFunctor(anaPhysDispOnLattice, "ana phys disp");

  vtmWriter.addFunctor(numLatticeStress, "num lattice stress");
  vtmWriter.addFunctor(numPhysStress, "num phys stress");
  vtmWriter.addFunctor(anaPhysStressOnLattice, "ana phys stress");
  vtmWriter.addFunctor(anaPhysForceOnLattice, "ana phys force");

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
    vtmWriter.write(iT);
  }
}

/// @brief Execute simulation: set initial values and run time loop
/// @param myCase The Case instance which keeps the simulation data
template <typename CASE>
void simulate(CASE& myCase)
{
  using T = CASE::value_t;
  OstreamManager clout(std::cout, "simulate");
  clout << "Starting simulation ..." << std::endl;

  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(names::NavierCauchy {});
  auto& parameters = myCase.getParameters();

  const T           physMaxT = parameters.template get<parameters::PHYS_MAX_T>();
  const std::size_t convIt   = parameters.template get<parameters::CONV_ITER>();
  const std::size_t iTmax    = lattice.getUnitConverter().getLatticeTime(physMaxT);

  parameters.template set<parameters::CONVERGED>(false);
  clout << "Max number of time steps: " << iTmax << std::endl;

  util::ValueTracer<T> converge(convIt, parameters.template get<parameters::CONVERGENCE_PRECISION>());

  util::Timer<T> timer(iTmax, geometry.getStatistics().getNvoxel());
  timer.start();
  for (std::size_t iT = 0; iT < iTmax; ++iT) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    //setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    lattice.collideAndStream();

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, timer, iT);

    if (iT % convIt == 0) {
      converge.takeValue(parameters.template get<parameters::CONV_VALUE>(), true);
    }

    if (converge.hasConverged()) {
      parameters.template set<parameters::CONVERGED>(true);
      getResults(myCase, timer, iT);
      clout << "Simulation converged." << std::endl;
      break;
    }
  }

  timer.stop();
  timer.printSummary();

  clout << "Simulation finished ..." << std::endl;
}

void setInitialControls(MyOptiCase& optiCase)
{
  OstreamManager clout(std::cout, "setInitialControls");
  clout << "Set Initial Control ..." << std::endl;
  auto& parameters          = optiCase.getCase(Controlled {}).getParameters();
  const auto initialControl = parameters.template get<parameters::INITIAL_CONTROL>();
  optiCase.getController().set(std::vector<MyCase::value_t>({initialControl[0], initialControl[1]}));
}

template <typename CASE, typename CONTROLLER>
void applyControls(CASE& controlledCase, const CONTROLLER& controller)
{
  using T          = CASE::value_t;
  auto& lattice    = controlledCase.getLattice(NavierCauchy {});
  auto& converter  = lattice.getUnitConverter();
  auto& geometry   = controlledCase.getGeometry();

  OstreamManager clout(std::cout, "applyControls");
  clout << "Applying control ..." << std::endl;

  clout << "Number of Controls: " << controller.get().size() << std::endl;
  clout << "Control [";
  for (std::size_t i = 0; i < controller.get().size(); i++) {
    clout << controller.get()[i] << ", ";
  }
  clout << "]" << std::endl;

  ForceField2D<T> force(converter.getLameMu(), converter.getLameLambda(), converter.getConversionFactorSolidForce());

  T c0 = controller.get()[0] + std::numeric_limits<double>::epsilon();
  T c1 = controller.get()[1] + std::numeric_limits<double>::epsilon();

  AnalyticalConst2D<T, T> control(c0, c1);
  clout << "target dims: " << force.getTargetDim() << " " << control.getTargetDim() << std::endl;

  fields::set<descriptors::FORCE>(lattice, geometry.getMaterialIndicator({1}), force * control);
  lattice.template setProcessingContext<Array<descriptors::FORCE>>(ProcessingContext::Simulation);
}

template <typename NAME>
auto displacementObjective(MyOptiCase& optiCase)
{
  using OBJECTIVE = functors::L2DistanceF<functors::DisplacementF>;

  auto& controlledCase    = optiCase.getCase(NAME{});
  auto& controlledLattice = controlledCase.getLattice(NavierCauchy {});
  auto& converter         = controlledLattice.getUnitConverter();
  auto  objectiveDomain   = controlledCase.getGeometry().getMaterialIndicator({1});
  auto& controlparameters = controlledCase.getParameters();

  auto& referenceLattice   = optiCase.getCase(Reference {}).getLattice(NavierCauchy {});
  auto& refConverter       = referenceLattice.getUnitConverter();
  auto  refObjectiveDomain = optiCase.getCase(Reference {}).getGeometry().getMaterialIndicator({1});
  auto& parameters         = optiCase.getCase(Reference {}).getParameters();

  // Evaluate functor for objective computation
  auto objectiveO = makeWriteFunctorO<OBJECTIVE, opti::J>(controlledLattice);
  objectiveO->restrictTo(objectiveDomain);

  // Write displacement field in the reference case
  if(parameters.get<parameters::ANALYTICAL_REFERENCE>()) {
    ManufacturedSolutionU2D<MyCase::value_t> anaPhysDisp(true);
    fields::set<OBJECTIVE::Reference>(referenceLattice, refObjectiveDomain, anaPhysDisp);
    referenceLattice.template setProcessingContext<Array<OBJECTIVE::Reference>>(ProcessingContext::Simulation);
  } else {
    auto displacementO = makeWriteFunctorO<functors::DisplacementF, OBJECTIVE::Reference>(referenceLattice);
    displacementO->restrictTo(refObjectiveDomain);
    displacementO->template setParameter<parameters::CONVERSION_DISPLACEMENT>(refConverter.getCharPhysDisplacement());
    displacementO->template setParameter<parameters::OMEGA_SOLID>(parameters.template get<parameters::OMEGA_SOLID>());
    displacementO->apply();
  }

  copyFields<OBJECTIVE::Reference, OBJECTIVE::Reference>(referenceLattice, controlledLattice);

  objectiveO->template setParameter<parameters::CONVERSION_DISPLACEMENT>(converter.getCharPhysDisplacement());
  objectiveO->template setParameter<parameters::OMEGA_SOLID>(controlparameters.template get<parameters::OMEGA_SOLID>());
  objectiveO->template setParameter<parameters::NORMALIZE>(
      computeL2Norm<OBJECTIVE::Reference>(referenceLattice, refObjectiveDomain, refConverter.getPhysDeltaX()));
  objectiveO->apply();

  return integrateField<opti::J>(controlledLattice, objectiveDomain, converter.getPhysDeltaX())[0];
}

template <typename NAME>
auto computeObjective(MyOptiCase& optiCase, NAME)
{
  using T = MyOptiCase::case_t<NAME>::value_t;
  OstreamManager clout(std::cout, "computeObjective");
  clout << "Computing Objective ..." << std::endl;

  auto& controlledCase = optiCase.getCase(NAME{});

  controlledCase.resetLattices();
  prepareLattice(controlledCase);
  setInitialValues(controlledCase);
  applyControls(controlledCase, optiCase.getController<T>());
  simulate(controlledCase);

  if constexpr (!util::is_adf_v<T>) {
    clout << "Getting Optimization Results ..." << std::endl;
    getOptiResults<NAME>(optiCase);
  }

  return displacementObjective<NAME>(optiCase);
}
