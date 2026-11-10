/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Florian Kaiser, Stephan Simonis
 *                2026 Florian Kaiser
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

#include "../manufactured_solution_periodic.h"

// === Step 1: Declarations ===
using MyCase = Case<
  names::NavierCauchy,
  Lattice<
    FLOATING_POINT_TYPE,
    descriptors::D2Q8<
      descriptors::FORCE,
      structuralmechanics::MOMENT_VECTOR,
      structuralmechanics::DISPLACEMENT,
      structuralmechanics::CAUCHY_STRESS
    >
  >
>;

namespace olb::parameters {
  struct CONV_VALUE   : public descriptors::FIELD_BASE<1> { };
}

/// @brief Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<names::NavierCauchy>;

  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector origin = parameters.get<parameters::ORIGIN       >();
  const T      dx     = parameters.get<parameters::PHYS_DELTA_X >();

  IndicatorCuboid2D<T> cuboid(extent - Vector<T, DESCRIPTOR::d>{dx, dx}, origin);

  Mesh<T,MyCase::d> mesh(cuboid, dx, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity( { true, true } );
  return mesh;
}

/// @brief Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
/// @note The material numbers will be used to assign physics to lattice nodes
void prepareGeometry(MyCase& myCase) {
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
void prepareLattice(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T           = MyCase::value_t;
  using DESCRIPTOR  = MyCase::descriptor_t_of<names::NavierCauchy>;

  auto& geometry    = myCase.getGeometry();
  auto& lattice     = myCase.getLattice(names::NavierCauchy{});
  auto& parameters  = myCase.getParameters();

  lattice.defineDynamics<structuralmechanics::StationaryLinearElasticityMRT>(geometry, 1);

  const T physCharLength        = parameters.get<parameters::PHYS_CHAR_LENGTH       >();
  const T dx                    = parameters.get<parameters::PHYS_DELTA_X           >();
  const T physCharDisplacement  = parameters.get<parameters::PHYS_CHAR_DISPLACEMENT >();
  const T physDensity           = parameters.get<parameters::PHYS_DENSITY           >();
  const T physYoungsModulus     = parameters.get<parameters::YOUNGS_MODULUS         >();
  const T physPoissonRatio      = parameters.get<parameters::POISSON_RATIO          >();
  const T kappa                 = parameters.get<parameters::DAMPING_FACTOR         >();
  const std::string solidCase   = parameters.get<parameters::SOLID_CASE             >();

  const T dt = T(0.1) / physYoungsModulus * kappa * dx * dx;

  lattice.setUnitConverter<StationaryMRTUnitConverter<T, DESCRIPTOR>>(
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

  Vector<T,DESCRIPTOR::q> Omega = {
    0,
    0,
    omega_11,
    omega_s,
    omega_d,
    omega_12,
    omega_21,
    omega_22
  };
  lattice.setParameter<structuralmechanics::OMEGA_SOLID>(Omega);

  clout << "Prepare Lattice ... OK" << std::endl;
}

/// Set initial condition for primal variables (velocity and density)
/// @param myCase The Case instance which keeps the simulation data
/// @note Be careful: initial values have to be set using lattice units
void setInitialValues(MyCase& myCase) {
  OstreamManager clout(std::cout, "setInitialValues");
  clout << "Set Initial Values ..." << std::endl;

  using T = MyCase::value_t;

  auto& geometry  = myCase.getGeometry();
  auto& lattice   = myCase.getLattice(names::NavierCauchy{});
  auto& converter = lattice.getUnitConverter();

  lattice.initialize();

  ForceField2D<T> force(converter.getLameLambda(), converter.getLameMu(), converter.getConversionFactorSolidForce());
  lattice.defineField<descriptors::FORCE>(geometry, 1, force);

  lattice.setProcessingContext<Array<descriptors::FORCE>>(ProcessingContext::Simulation);

  clout << "Set Initial Values ... OK" << std::endl;
}

/// Update boundary values at times (and external fields, if they exist)
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
/// @note Be careful: boundary values have to be set using lattice units
void setTemporalValues(MyCase& myCase,
                       std::size_t iT) {
  // Nothing to do here, because simulation does not depend on time
}

/// Compute simulation results at times
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT) {
  OstreamManager clout(std::cout, "getResults");

  using T           = MyCase::value_t;
  using DESCRIPTOR  = MyCase::descriptor_t_of<names::NavierCauchy>;

  auto& geometry    = myCase.getGeometry();
  auto& lattice     = myCase.getLattice(names::NavierCauchy{});
  auto& converter   = lattice.getUnitConverter();
  auto& parameters  = myCase.getParameters();

  const std::size_t iTMax     = converter.getLatticeTime(parameters.get<parameters::PHYS_MAX_T>());
  const std::size_t iTlog     = parameters.get<parameters::LATTICE_LOG_ITER_T                 >();
  const std::size_t iTvtk     = parameters.get<parameters::LATTICE_VTK_ITER_T                 >();
  const std::size_t convIter  = parameters.get<parameters::CONV_ITER                          >();
  const bool        converged = parameters.get<parameters::CONVERGED                          >();
  const bool        writeVTK  = parameters.get<parameters::VTK_ENABLED                        >();

  if (iT == 0) {
    SuperVTMwriter2D<T> vtmWriter("periodicPlate");

    SuperLatticeCuboid2D<T, DESCRIPTOR>   cuboid(lattice);
    SuperLatticeRank2D<T, DESCRIPTOR>     rank(lattice);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();
  }

  if (writeVTK && (iT % iTvtk == 0 || iT % iTMax == 0 || converged)) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperVTMwriter2D<T> vtmWriter("periodicPlate");

    ManufacturedSolutionU2D<T>                    anaPhysDisp(true);
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> anaPhysDispOnLattice(anaPhysDisp, lattice);

    ManufacturedSolutionStress2D<T>               anaPhysStress(converter.getLameLambda(), converter.getLameMu());
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> anaPhysStressOnLattice(anaPhysStress, lattice);

    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::DISPLACEMENT  > numPhysDisp(lattice, converter.getConversionFactorDisplacement());
    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::CAUCHY_STRESS > numPhysStress(lattice, converter.getConversionFactorSolidStress());
    SuperLatticePhysField2D<T, DESCRIPTOR, descriptors::FORCE                 > numPhysForce(lattice, converter.getConversionFactorSolidForce());

    vtmWriter.addFunctor(numPhysDisp,            "num phys disp"  );
    vtmWriter.addFunctor(anaPhysDispOnLattice,   "ana phys disp"  );

    vtmWriter.addFunctor(numPhysStress,          "num phys stress");
    vtmWriter.addFunctor(anaPhysStressOnLattice, "ana phys stress");

    vtmWriter.addFunctor(numPhysForce,           "num phys force" );
    vtmWriter.write(iT);

    lattice.setProcessingContext(ProcessingContext::Simulation);
  }

  if (iT % iTlog == 0 || iT % iTMax == 0 || converged) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    ManufacturedSolutionU2D<T>      anaPhysDisp(true);
    ManufacturedSolutionStress2D<T> anaPhysStress(converter.getLameLambda(), converter.getLameMu());

    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::DISPLACEMENT  > numPhysDisp(lattice, converter.getConversionFactorDisplacement());
    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::CAUCHY_STRESS > numPhysStress(lattice, converter.getConversionFactorSolidStress());

    T   l2UResult[2]        = {T(), T()};
    T   lInfUResult[2]      = {T(), T()};
    T   l2StressResult[2]   = {T(), T()};
    T   lInfStressResult[2] = {T(), T()};
    int tmp[]               = {int()};

    auto indicatorF = geometry.getMaterialIndicator(1);
    SuperRelativeErrorL2Norm2D<T>   relUErrorL2Norm(lattice, numPhysDisp, anaPhysDisp, indicatorF);
    SuperRelativeErrorLinfNorm2D<T> relUErrorLinfNorm(lattice, numPhysDisp, anaPhysDisp, indicatorF);
    SuperRelativeErrorL2Norm2D<T>   relStressErrorL2Norm( lattice, numPhysStress, anaPhysStress, indicatorF );
    SuperRelativeErrorLinfNorm2D<T> relStressErrorLinfNorm( lattice, numPhysStress, anaPhysStress, indicatorF );

    relUErrorL2Norm(l2UResult, tmp);
    relUErrorLinfNorm(lInfUResult, tmp);
    relStressErrorL2Norm(l2StressResult, tmp);
    relStressErrorLinfNorm(lInfStressResult, tmp);

    clout << "N\t" << "L2 U Error\t" << "LInf U Error\t" << "L2 Stress Error\t" << "LInf Stress Error" << std::endl;
    clout << converter.getResolution() << "\t" << l2UResult[0] << "\t" << lInfUResult[0] << "\t" << l2StressResult[0] << "\t\t" << lInfStressResult[0] << std::endl;

    timer.print(iT);
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));

    lattice.setProcessingContext(ProcessingContext::Simulation);
  }


  if ( iT % convIter == 0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    T   lInfUResult[2]  = {T(), T()};
    int tmp[]           = {int()};

    auto bulk           = geometry.getMaterialIndicator(1);

    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::DISPLACEMENT> numPhysDisp(lattice, converter.getConversionFactorDisplacement());
    ManufacturedSolutionU2D<T> anaPhysDisp(true);

    SuperRelativeErrorLinfNorm2D<T> relUErrorLinfNorm(lattice, numPhysDisp, anaPhysDisp, bulk);
    relUErrorLinfNorm(lInfUResult, tmp);

    parameters.set<parameters::CONV_VALUE>(lInfUResult[0]);

    lattice.setProcessingContext(ProcessingContext::Simulation);
  }

  if (converged) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    auto indicatorF = geometry.getMaterialIndicator({1});

    T   l2UResult[2]        = {T(), T()};
    T   lInfUResult[2]      = {T(), T()};
    T   l2StressResult[2]   = {T(), T()};
    T   lInfStressResult[2] = {T(), T()};
    int tmp[]               = {int()};

    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::DISPLACEMENT  > dispNumPhys(lattice, converter.getConversionFactorDisplacement());
    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::CAUCHY_STRESS > stressNumPhys(lattice, converter.getConversionFactorSolidStress());

    ManufacturedSolutionU2D<T>      dispAnaPhys(true);

    ManufacturedSolutionStress2D<T> stressAnaPhys(converter.getLameLambda(), converter.getLameMu());

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

/// @brief Execute simulation: set initial values and run time loop
/// @param myCase The Case instance which keeps the simulation data
void simulate(MyCase& myCase) {
  OstreamManager clout(std::cout, "simulate");
  clout << "Starting simulation ..." << std::endl;

  using T = MyCase::value_t;

  auto& geometry            = myCase.getGeometry();
  auto& lattice             = myCase.getLattice(names::NavierCauchy{});
  auto& parameters          = myCase.getParameters();

  const T physMaxT          = parameters.get<parameters::PHYS_MAX_T>();
  const std::size_t convIt  = parameters.get<parameters::CONV_ITER>();
  const std::size_t iTmax   = lattice.getUnitConverter().getLatticeTime(physMaxT);

  util::ValueTracer<T> converge(convIt, parameters.get<parameters::CONVERGENCE_PRECISION>());

  util::Timer<T> timer(iTmax, geometry.getStatistics().getNvoxel());
  timer.start();
  for (std::size_t iT = 0; iT < iTmax; ++iT) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    lattice.collideAndStream();

    lattice.communicate();

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, timer, iT);

    if (iT % convIt == 0) {
      converge.takeValue(parameters.get<parameters::CONV_VALUE>(), true);
    }

    if (converge.hasConverged()) {
      parameters.set<parameters::CONVERGED>(true);
      getResults(myCase, timer, iT);
      clout << "Simulation converged." << std::endl;
      break;
    }

    if (iT == iTmax) {
      parameters.set<parameters::CONVERGED>(true);
      getResults(myCase, timer, iT);

      clout << "Simulation not converged. Writing error at last time step." << std::endl;
      break;
    }

  }

  timer.stop();
  timer.printSummary();

  clout << "Simulation finished ..." << std::endl;
}

/// Setup and run a simulation
int main(int argc, char* argv[]) {
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<SOLID_CASE               >( "continua"     );
    myCaseParameters.set<RESOLUTION               >(   40           );
    myCaseParameters.set<YOUNGS_MODULUS           >(    0.11        );
    myCaseParameters.set<PHYS_DENSITY             >(    1           );
    myCaseParameters.set<POISSON_RATIO            >(    0.8         );
    myCaseParameters.set<DAMPING_FACTOR           >(    1.          );
    myCaseParameters.set<PHYS_CHAR_LENGTH         >(    1.0         );
    myCaseParameters.set<DOMAIN_EXTENT            >(  { 1.0, 1.0 }  );
    myCaseParameters.set<ORIGIN                   >(  { 0.0, 0.0 }  );
    myCaseParameters.set<PHYS_MAX_T               >( 1000           );
    myCaseParameters.set<LATTICE_LOG_ITER_T       >( 1000           );
    myCaseParameters.set<LATTICE_VTK_ITER_T       >( 1000           );
    myCaseParameters.set<CONV_ITER                >(   10           );
    myCaseParameters.set<CONV_VALUE               >(    1           );
    myCaseParameters.set<CONVERGENCE_PRECISION    >(    1e-7        );
    myCaseParameters.set<VTK_ENABLED              >(  true          );
  }

  myCaseParameters.set<parameters::PHYS_DELTA_X>(
    myCaseParameters.get<parameters::DOMAIN_EXTENT>()[0] / myCaseParameters.get<parameters::RESOLUTION>()
  );
  myCaseParameters.set<parameters::PHYS_DELTA_T>(
    myCaseParameters.get<parameters::PHYS_DELTA_X>() * myCaseParameters.get<parameters::PHYS_DELTA_X>()
  );
  myCaseParameters.set<parameters::PHYS_CHAR_DISPLACEMENT>(
    myCaseParameters.get<parameters::PHYS_DELTA_X>()
  );

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
