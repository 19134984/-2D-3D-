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

/** \file A two-dimensional parameter identification problem for quasistationary
 *  linear elasticity: Young's modulus of the combined ellipse domain is
 *  recovered from a reference displacement field.
 *
 *  The corresponding forward problem is the example
 *  `structuralmechanics/linearElasticity/quasistationary/combinedEllipse`.
 *
 *  Grafen et al. "Benchmark cases for gradient-based optimization in linear elastic solid
 *  mechanics using lattice Boltzmann methods" arXiv preprint arXiv:2609.31347 (2026).
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
struct CONV_VALUE            : public descriptors::FIELD_BASE<1> { };
struct INITIAL_CONTROL       : public descriptors::FIELD_BASE<1> { };
struct ANALYTICAL_REFERENCE  : public fields::SCALAR<false>      { };
struct TARGET_YOUNGS_MODULUS : public descriptors::FIELD_BASE<1> { };
struct OPTI_VTK_ENABLED      : public fields::SCALAR<true>       { };
struct FIXED_TIME_STEP       : public descriptors::TYPED_FIELD_BASE<bool,1> { };
struct FORWARD_ONLY          : public descriptors::TYPED_FIELD_BASE<bool,1> { };
} // namespace olb::parameters

/// One-dimensional ADf type: Young's modulus is the only control
using ADf = util::ADf<double, 1>;

// === Step 1: Declarations ===
using MyCase = Case<names::NavierCauchy,Lattice<double,descriptors::D2Q8<>>>;
using MyADfCase  = MyCase::exchange_value_t<ADf>;
using MyOptiCase = OptiCaseADf<
  Controlled,  MyCase,
  Derivatives, MyADfCase,
  Reference,   MyCase
>;

/// @brief Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
template <typename PARAMETERS>
auto createMesh(PARAMETERS& parameters)
{
  using T = PARAMETERS::value_t;
  constexpr unsigned D = PARAMETERS::descriptor_t::d;

  const Vector domainExtent = parameters.template get<parameters::DOMAIN_EXTENT>();
  const Vector origin       = parameters.template get<parameters::ORIGIN       >();
  const T      dx           = parameters.template get<parameters::PHYS_DELTA_X >();

  IndicatorCuboid2D<T> cuboid(domainExtent + Vector<T,D>{dx, dx}, origin);

  Mesh<T,D> mesh(cuboid, dx, singleton::mpi().getSize());
  mesh.setOverlap(parameters.template get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({false, false});
  return mesh;
}

/// @brief Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
/// @note The material numbers will be used to assign physics to lattice nodes
template <typename CASE>
void prepareGeometry(CASE& myCase)
{
  using T = CASE::value_t;
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  auto& geometry = myCase.getGeometry();

  ellipse2D<T> ellipseCase(myCase.getParameters().template get<parameters::PHYS_CHAR_LENGTH>());

  /// Set material numbers
  geometry.rename(0, 2);
  geometry.rename(2, 1, ellipseCase.ellipse1);
  geometry.rename(1, 3, ellipseCase.ellipse2);
  geometry.rename(1, 4, ellipseCase.ellipse3);

  geometry.innerClean();
  geometry.clean();
  geometry.checkForErrors();

  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

/// @brief Set lattice dynamics, boundaries and the unit converter
/// @param myCase The Case instance which keeps the simulation data
template <typename CASE>
void prepareLattice(CASE& myCase)
{
  using T          = CASE::value_t;
  using DESCRIPTOR = CASE::template descriptor_t_of<names::NavierCauchy>;

  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(names::NavierCauchy{});
  auto& parameters = myCase.getParameters();

  const T dx                   = parameters.template get<parameters::PHYS_DELTA_X          >();
  const T physCharLength       = parameters.template get<parameters::PHYS_CHAR_LENGTH      >();
  const T physCharDisplacement = parameters.template get<parameters::PHYS_CHAR_DISPLACEMENT>();
  const T physDensity          = parameters.template get<parameters::PHYS_DENSITY          >();
  const T youngsModulus        = parameters.template get<parameters::YOUNGS_MODULUS        >();
  const T poissonRatio         = parameters.template get<parameters::POISSON_RATIO         >();
  const T kappa                = parameters.template get<parameters::DAMPING_FACTOR        >();
  const T targetYoungsModulus  = parameters.template get<parameters::TARGET_YOUNGS_MODULUS >();
  const std::string solidCase  = parameters.template get<parameters::SOLID_CASE            >();

  ellipse2D<T> ellipseCase(physCharLength);

  const bool fixedTimeStep = parameters.template get<parameters::FIXED_TIME_STEP>();
  const T dt = T(0.1) * kappa * dx * dx / (fixedTimeStep ? targetYoungsModulus : youngsModulus);
  parameters.template set<parameters::PHYS_DELTA_T>(dt);

  // Set up a unit converter with the characteristic physical units
  lattice.template setUnitConverter<StationaryMRTUnitConverter<T, DESCRIPTOR>>(
    dx,
    dt,
    physCharLength,
    T(1),
    physCharDisplacement,
    youngsModulus,
    poissonRatio,
    physDensity,
    kappa,
    solidCase
  );

  const auto& converter = lattice.getUnitConverter();
  converter.print();

  // Converter of the sought-after material, used for the moduli of the manufactured
  // Neumann traction and force (the loads must not depend on the control).
  StationaryMRTUnitConverter<T, DESCRIPTOR> targetConverter(
    dx,
    dt,
    physCharLength,
    T(1),
    physCharDisplacement,
    targetYoungsModulus,
    poissonRatio,
    physDensity,
    kappa,
    solidCase
  );

  auto bulkIndicator = geometry.getMaterialIndicator(1);

  dynamics::set<structuralmechanics::StationaryLinearElasticityMRTBoundary>(lattice, bulkIndicator);

  setStructuralMechanicsBoundary<T, DESCRIPTOR, structuralmechanics::StationaryLinearElasticityMRTNeumannValidationFixedMaterial>(
    lattice, geometry.getMaterialIndicator(2), bulkIndicator, ellipseCase.ellipse1, true);
  setStructuralMechanicsBoundary<T, DESCRIPTOR, structuralmechanics::StationaryLinearElasticityMRTDirichletValidation>(
    lattice, geometry.getMaterialIndicator(3), bulkIndicator, ellipseCase.ellipse2, false);
  setStructuralMechanicsBoundary<T, DESCRIPTOR, structuralmechanics::StationaryLinearElasticityMRTDirichletValidation>(
    lattice, geometry.getMaterialIndicator(4), bulkIndicator, ellipseCase.ellipse3, false);

  lattice.template setParameter<descriptors::DX                    >(dx                               );
  lattice.template setParameter<descriptors::DT                    >(dt                               );
  lattice.template setParameter<structuralmechanics::DAMPING_FACTOR>(converter.getDampingFactor()     );
  lattice.template setParameter<structuralmechanics::LAME_MU       >(converter.getLatticeLameMu()     );
  lattice.template setParameter<structuralmechanics::LAME_LAMBDA   >(converter.getLatticeLameLambda() );
  lattice.template setParameter<structuralmechanics::BULK_MODULUS  >(converter.getLatticeBulkModulus());

  lattice.template setParameter<structuralmechanics::CHAR_U         >(converter.getConversionFactorDisplacement());

  lattice.template setParameter<structuralmechanics::REFERENCE_LAME_MU     >(targetConverter.getLatticeLameMu()     );
  lattice.template setParameter<structuralmechanics::REFERENCE_BULK_MODULUS>(targetConverter.getLatticeBulkModulus());

  const T theta = T(1) / descriptors::invCs2<T, DESCRIPTOR>();

  const T omega_11 = T(1.) / (        converter.getLatticeLameMu()                                    /         theta  + T(0.5));
  const T omega_d  = T(1.) / (T(2.) *  converter.getLatticeLameMu()                                   / (T(1.) - theta) + T(0.5));
  const T omega_s  = T(1.) / (T(2.) * (converter.getLatticeLameMu() + converter.getLatticeLameLambda()) / (T(1.) + theta) + T(0.5));

  const T tau_12 = T(0.5);
  const T tau_22 = T(0.5);
  const T tau_21 = tau_12;

  const T omega_12 = T(1.) / (tau_12 + T(0.5));
  const T omega_21 = T(1.) / (tau_21 + T(0.5));
  const T omega_22 = T(1.) / (tau_22 + T(0.5));

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

  Vector neumann_constants = {
    (T(2.) * (T(1.) - theta) * (converter.getLatticeBulkModulus() - converter.getLatticeLameMu()))
      / (theta * (T(1.) - theta - T(4.) * converter.getLatticeLameMu())),
     T(2.) * converter.getLatticeLameMu() / (theta - T(2.) * converter.getLatticeLameMu()),
     T(4.) * converter.getLatticeLameMu() / (T(1.) - theta - T(4.) * converter.getLatticeLameMu())
  };
  lattice.template setParameter<structuralmechanics::NEUMANN_SOLID_C>(neumann_constants);

  {
    auto& communicator = lattice.getCommunicator(stage::PostCollide());
    communicator.template requestField<descriptors::POPULATION>();
    communicator.template requestField<structuralmechanics::DISPLACEMENT>();
    communicator.template requestField<structuralmechanics::CAUCHY_STRESS>();
    communicator.template requestField<structuralmechanics::BARED_MOMENT_VECTOR>();
    communicator.template requestField<structuralmechanics::MOMENT_VECTOR>();
    communicator.template requestField<structuralmechanics::PREVIOUS_CELL>();
    communicator.requestOverlap(parameters.template get<parameters::OVERLAP>());
    communicator.exchangeRequests();
  }

  {
    auto& communicator = lattice.getCommunicator(stage::PreStream());
    communicator.template requestField<descriptors::POPULATION>();
    communicator.requestOverlap(parameters.template get<parameters::OVERLAP>());
    communicator.exchangeRequests();
  }

  clout << "  simulated E       = " << youngsModulus       << std::endl;
  clout << "  target E          = " << targetYoungsModulus  << std::endl;
  clout << "  simulated Poisson = " << poissonRatio         << std::endl;
  clout << "Prepare Lattice ... OK" << std::endl;
}

/// Set initial condition and the body force of the manufactured solution
/// @param myCase The Case instance which keeps the simulation data
/// @note The body force belongs to the sought-after (target) material, so that
///       the reference displacement is reproduced exactly at the target control.
///       It is converted to lattice units with the controlled time step.
template <typename CASE>
void setInitialValues(CASE& myCase)
{
  using T          = CASE::value_t;
  using DESCRIPTOR = CASE::template descriptor_t_of<names::NavierCauchy>;

  OstreamManager clout(std::cout, "setInitialValues");
  clout << "Set Initial Values ..." << std::endl;

  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(names::NavierCauchy{});
  auto& parameters = myCase.getParameters();

  const auto& converter = lattice.getUnitConverter();

  const T physCharLength       = parameters.template get<parameters::PHYS_CHAR_LENGTH      >();
  const T physCharDisplacement = parameters.template get<parameters::PHYS_CHAR_DISPLACEMENT>();
  const T physDensity          = parameters.template get<parameters::PHYS_DENSITY          >();
  const T poissonRatio         = parameters.template get<parameters::POISSON_RATIO         >();
  const T kappa                = parameters.template get<parameters::DAMPING_FACTOR        >();
  const T targetYoungsModulus  = parameters.template get<parameters::TARGET_YOUNGS_MODULUS >();
  const std::string solidCase  = parameters.template get<parameters::SOLID_CASE            >();

  StationaryMRTUnitConverter<T, DESCRIPTOR> targetConverter(
    converter.getPhysDeltaX(),
    converter.getPhysDeltaT(),
    physCharLength,
    T(1),
    physCharDisplacement,
    targetYoungsModulus,
    poissonRatio,
    physDensity,
    kappa,
    solidCase
  );

  lattice.initialize();

  ForceField2D<T> force(targetConverter.getLameMu(),
                        targetConverter.getLameLambda(),
                        converter.getConversionFactorSolidForce());
  fields::set<descriptors::FORCE>(lattice, geometry.getMaterialIndicator({1}), force);

  lattice.setProcessingContext(ProcessingContext::Simulation);

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
/// @param timer The timer of the running pseudo time loop
/// @param iT The time step
template <typename CASE>
void getResults(CASE& myCase, util::Timer<typename CASE::value_t>& timer, std::size_t iT)
{
  using T          = CASE::value_t;
  using DESCRIPTOR = CASE::template descriptor_t_of<names::NavierCauchy>;

  OstreamManager clout(std::cout, "getResults");

  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(names::NavierCauchy{});
  auto& parameters = myCase.getParameters();

  const auto& converter = lattice.getUnitConverter();

  const std::size_t iTlog     = parameters.template get<parameters::LATTICE_LOG_ITER_T>();
  const std::size_t iTMax     = parameters.template get<parameters::LATTICE_MAX_T>();
  const std::size_t convIter  = parameters.template get<parameters::CONV_ITER>();
  const bool        converged = parameters.template get<parameters::CONVERGED>();

  if (iT % iTlog == 0 || iT % iTMax == 0 || converged) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    ManufacturedSolutionU2D<T>      anaPhysDisp(true);
    ManufacturedSolutionStress2D<T> anaPhysStress(converter.getLameMu(), converter.getLameLambda());

    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::DISPLACEMENT> numPhysDisp(lattice, converter.getConversionFactorDisplacement());
    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::CAUCHY_STRESS      > numPhysStress(lattice, converter.getConversionFactorSolidStress());

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
    parameters.template set<parameters::CONV_VALUE>(lattice.getStatistics().getAverageEnergy());
  }
}

/// Write the relative errors of the final state against the manufactured solution
/// to disp_err.csv and stress_err.csv
template <typename CASE>
void writeErrors(CASE& myCase)
{
  using T          = CASE::value_t;
  using DESCRIPTOR = CASE::template descriptor_t_of<names::NavierCauchy>;

  auto& geometry  = myCase.getGeometry();
  auto& lattice   = myCase.getLattice(names::NavierCauchy{});
  auto& converter = lattice.getUnitConverter();

  lattice.setProcessingContext(ProcessingContext::Evaluation);
  auto indicatorF = geometry.getMaterialIndicator({1});

  T   l2UResult[2]        = {T(), T()};
  T   lInfUResult[2]      = {T(), T()};
  T   l2StressResult[2]   = {T(), T()};
  T   lInfStressResult[2] = {T(), T()};
  int tmp[]               = {int()};

  SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::DISPLACEMENT> dispNumPhys(lattice, converter.getConversionFactorDisplacement());
  SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::CAUCHY_STRESS      > stressNumPhys(lattice, converter.getConversionFactorSolidStress());

  ManufacturedSolutionU2D<T>      dispAnaPhys(true);
  ManufacturedSolutionStress2D<T> stressAnaPhys(converter.getLameMu(), converter.getLameLambda());

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

  const T resolution = myCase.getParameters().template get<parameters::RESOLUTION>();
  csvWriterU.writeDataFile(resolution, {l2UResult[0], lInfUResult[0]});
  csvWriterSigma.writeDataFile(resolution, {l2StressResult[0], lInfStressResult[0]});

  lattice.setProcessingContext(ProcessingContext::Simulation);
}

/// Write the state of one optimization step of the ADf-typed case to VTK
/// @param optiCase The OptiCase instance holding all cases
/// @note Every field is written by its primal value; the numerical fields and the
///       objective are written by their derivatives w.r.t. the controls as well
template <typename NAME>
void getOptiResults(MyOptiCase& optiCase)
{
  using CASE       = MyOptiCase::case_t<NAME>;
  using T          = CASE::value_t;
  using DESCRIPTOR = CASE::template descriptor_t_of<names::NavierCauchy>;
  static_assert(util::is_adf_v<T>, "getOptiResults expects the ADf-typed case");

  OstreamManager clout(std::cout, "getOptiResults");

  auto& myCase     = optiCase.getCase(NAME{});
  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(names::NavierCauchy{});

  const auto& converter = lattice.getUnitConverter();

  // Count the calls locally: the optimization step is not yet incremented between the
  // first two derivative evaluations, which would overwrite the output of the first one
  static std::size_t callCounter = 0;
  const std::size_t iT = callCounter++;

  SuperVTMwriter<T,CASE::d,float,BaseType<T>> vtmWriter("parameterIdentificationEllipseCombo");
  SuperGeometryF<T,CASE::d> materials(geometry);

  ManufacturedSolutionU2D<T>                    anaPhysDisp(true);
  SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> anaPhysDispOnLattice(anaPhysDisp, lattice);

  ManufacturedSolutionStress2D<T>               anaPhysStress(converter.getLameMu(), converter.getLameLambda());
  SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> anaPhysStressOnLattice(anaPhysStress, lattice);

  ForceField2D<T>                               anaPhysForce(converter.getLameMu(), converter.getLameLambda(), T(1));
  SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> anaPhysForceOnLattice(anaPhysForce, lattice);

  SuperLatticeField2D<    T, DESCRIPTOR, structuralmechanics::DISPLACEMENT> numLatticeDisp(lattice);
  SuperLatticeField2D<    T, DESCRIPTOR, structuralmechanics::CAUCHY_STRESS      > numLatticeStress(lattice);
  SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::DISPLACEMENT> numPhysDisp(lattice, converter.getConversionFactorDisplacement());
  SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::CAUCHY_STRESS      > numPhysStress(lattice, converter.getConversionFactorSolidStress());
  SuperLatticeField2D<    T, DESCRIPTOR, opti::J                  > objective(lattice);

  // Adapters have to outlive the writer calls below
  std::vector<std::unique_ptr<SuperADfLayerF<CASE::d,T>>> layers;
  auto addPrimal = [&](SuperF2D<T,T>& f, const std::string& name) {
    layers.emplace_back(std::make_unique<SuperADfLayerF<CASE::d,T>>(f));
    vtmWriter.addFunctor(*layers.back(), name);
  };
  auto addPrimalAndDerivatives = [&](SuperF2D<T,T>& f, const std::string& name) {
    addPrimal(f, name);
    for (unsigned iD = 0; iD < T::dim; ++iD) {
      layers.emplace_back(std::make_unique<SuperADfLayerF<CASE::d,T>>(f, iD));
      vtmWriter.addFunctor(*layers.back(), "d(" + name + ")/dc" + std::to_string(iD));
    }
  };

  addPrimal(materials, "material");
  addPrimalAndDerivatives(numLatticeDisp,   "num lattice disp");
  addPrimalAndDerivatives(numPhysDisp,      "num phys disp");
  addPrimal(anaPhysDispOnLattice,           "ana phys disp");
  addPrimalAndDerivatives(numLatticeStress, "num lattice stress");
  addPrimalAndDerivatives(numPhysStress,    "num phys stress");
  addPrimal(anaPhysStressOnLattice,         "ana phys stress");
  addPrimal(anaPhysForceOnLattice,          "ana phys force");
  addPrimalAndDerivatives(objective,        "objective");

  lattice.setProcessingContext(ProcessingContext::Evaluation);

  if (iT == 0) {
    vtmWriter.createMasterFile();
    using T_PRIMAL          = MyCase::value_t;
    using DESCRIPTOR_PRIMAL = MyCase::descriptor_t;
    auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(names::NavierCauchy{});

    SuperVTMwriter<T_PRIMAL,CASE::d>                  decompositionWriter("parameterIdentificationEllipseCombo");
    SuperLatticeCuboid2D<T_PRIMAL, DESCRIPTOR_PRIMAL> cuboid(controlledLattice);
    SuperLatticeRank2D<T_PRIMAL, DESCRIPTOR_PRIMAL>   rank(controlledLattice);
    decompositionWriter.write(cuboid);
    decompositionWriter.write(rank);
  }
  vtmWriter.write(iT);

  lattice.setProcessingContext(ProcessingContext::Simulation);
}

/// @brief Execute simulation: run the pseudo time loop until convergence
/// @param myCase The Case instance which keeps the simulation data
template <typename CASE>
void simulate(CASE& myCase)
{
  using T = CASE::value_t;
  OstreamManager clout(std::cout, "simulate");

  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(names::NavierCauchy{});
  auto& parameters = myCase.getParameters();

  parameters.template set<parameters::CONVERGED>(false);

  const std::size_t maxIt  = parameters.template get<parameters::LATTICE_MAX_T>();
  const std::size_t convIt = parameters.template get<parameters::CONV_ITER>();

  util::ValueTracer<T> converge(convIt, parameters.template get<parameters::CONVERGENCE_PRECISION>());

  clout << "Awaiting " << maxIt << " Time steps. Starting simulation..." << std::endl;
  util::Timer<T> timer(maxIt, geometry.getStatistics().getNvoxel());
  timer.start();
  for (std::size_t iT = 0; iT <= maxIt; ++iT) {
    timer.update(iT);

    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Computation and Output of the Results ===
    getResults(myCase, timer, iT);

    /// === Step 8.3: Collide and Stream Execution ===
    lattice.collideAndStream();

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
}

// ==== BELOW HERE OPTIMIZATION SPECIFIC ====

/// @brief Set the initial control from the INITIAL_CONTROL parameter
/// @note The control is the square root of Young's modulus, which keeps the
///       identified modulus positive throughout the optimization
void setInitialControls(MyOptiCase& optiCase)
{
  OstreamManager clout(std::cout, "setInitialControls");
  clout << "Set Initial Control ..." << std::endl;

  auto&      parameters     = optiCase.getCase(Controlled{}).getParameters();
  const auto initialControl = parameters.template get<parameters::INITIAL_CONTROL>();

  optiCase.getController().set(std::vector<MyCase::value_t>{util::sqrt(initialControl)});
}

/// @brief Write the current control into the parameters of the controlled case
template <typename CASE, typename CONTROLLER>
void applyControls(CASE& controlledCase, const CONTROLLER& controller)
{
  using T = CASE::value_t;

  OstreamManager clout(std::cout, "applyControls");
  clout << "Applying control ..." << std::endl;

  auto& parameters = controlledCase.getParameters();

  const T c0            = controller.get()[0];
  const T youngsModulus = c0 * c0;
  parameters.template set<parameters::YOUNGS_MODULUS>(youngsModulus);

  clout << "Number of controls: "       << controller.get().size() << std::endl;
  clout << "Young's modulus control = " << youngsModulus           << std::endl;
}

/// @brief Evaluate the L2 distance between controlled and reference displacement
template <typename NAME>
auto displacementObjective(MyOptiCase& optiCase)
{
  using OBJECTIVE = functors::L2DistanceF<functors::DisplacementF>;

  OstreamManager clout(std::cout, "displacementObjective");

  auto& controlledCase       = optiCase.getCase(NAME{});
  auto& controlledLattice    = controlledCase.getLattice(NavierCauchy{});
  auto& converter            = controlledLattice.getUnitConverter();
  auto  objectiveDomain      = controlledCase.getGeometry().getMaterialIndicator({1});
  auto& controlledParameters = controlledCase.getParameters();

  auto& referenceLattice   = optiCase.getCase(Reference{}).getLattice(NavierCauchy{});
  auto& refConverter       = referenceLattice.getUnitConverter();
  auto  refObjectiveDomain = optiCase.getCase(Reference{}).getGeometry().getMaterialIndicator({1});
  auto& parameters         = optiCase.getCase(Reference{}).getParameters();

  // Evaluate functor for objective computation
  auto objectiveO = makeWriteFunctorO<OBJECTIVE, opti::J>(controlledLattice);
  objectiveO->restrictTo(objectiveDomain);

  // Write displacement field in the reference case
  if (parameters.template get<parameters::ANALYTICAL_REFERENCE>()) {
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

  // Get solution from the reference simulation for the inverse problem
  copyFields<OBJECTIVE::Reference, OBJECTIVE::Reference>(referenceLattice, controlledLattice);

  objectiveO->template setParameter<parameters::CONVERSION_DISPLACEMENT>(converter.getCharPhysDisplacement());
  objectiveO->template setParameter<parameters::OMEGA_SOLID>(controlledParameters.template get<parameters::OMEGA_SOLID>());
  objectiveO->template setParameter<parameters::NORMALIZE>(
    computeL2Norm<OBJECTIVE::Reference>(referenceLattice, refObjectiveDomain, refConverter.getPhysDeltaX()));
  objectiveO->apply();

  const auto J = integrateField<opti::J>(controlledLattice, objectiveDomain, converter.getPhysDeltaX())[0];
  clout << "J = " << J << std::endl;
  return J;
}

/// @brief Run the forward problem for the current control and evaluate the objective
template <typename NAME>
auto computeObjective(MyOptiCase& optiCase, NAME)
{
  using T = MyOptiCase::case_t<NAME>::value_t;

  OstreamManager clout(std::cout, "computeObjective");
  clout << "Computing Objective ..." << std::endl;

  // decide whether we solve for value or derivatives
  auto& controlledCase = optiCase.getCase(NAME{});
  // Reset prior simulation lattice
  controlledCase.resetLattices();

  // Set updated controls in the simulation
  applyControls(controlledCase, optiCase.getController<T>());

  // Prepare new lattices
  prepareLattice(controlledCase);
  setInitialValues(controlledCase);

  // Execute simulation
  simulate(controlledCase);

  // Compute Objective value from simulation results
  const auto J = displacementObjective<NAME>(optiCase);

  // Only the ADf-typed case carries the derivatives next to the primal values; the
  // objective field exists only after the objective has been evaluated
  if constexpr (util::is_adf_v<T>) {
    if (controlledCase.getParameters().template get<parameters::OPTI_VTK_ENABLED>()) {
      clout << "Getting Optimization Results ..." << std::endl;
      getOptiResults<NAME>(optiCase);
    }
  }

  return J;
}
