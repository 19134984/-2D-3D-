/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Florian Kaiser
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
#include "olb.h"

using namespace olb;

#include "../manufactured_solution_boundaries.h"

using MyCase = Case<
  names::NavierCauchy,
  Lattice<
    FLOATING_POINT_TYPE,
    descriptors::D2Q8<
      descriptors::FORCE,
      structuralmechanics::PREVIOUS_CELL,
      structuralmechanics::BARED_MOMENT_VECTOR,
      structuralmechanics::DISPLACEMENT,
      structuralmechanics::CAUCHY_STRESS
    >
  >
>;

namespace olb::parameters {
  struct CONV_VALUE : public descriptors::FIELD_BASE<1> { };
}

class ellipse2D {
  Vector<MyCase::value_t,MyCase::d> center1;
  Vector<MyCase::value_t,MyCase::d> center2;
  Vector<MyCase::value_t,MyCase::d> center3;

  public:
    ellipse2D() :
      center1(Vector<MyCase::value_t,MyCase::d>(.75, .75)),
      center2(Vector<MyCase::value_t,MyCase::d>(.95, .6)),
      center3(Vector<MyCase::value_t,MyCase::d>(.5, .9)),
      ellipse1(IndicatorEllipse2D<MyCase::value_t>(center1, .693, .548, -20)),
      ellipse2(IndicatorEllipse2D<MyCase::value_t>(center2, .1, .134)),
      ellipse3(IndicatorEllipse2D<MyCase::value_t>(center3, .187, .1))
    {};

    IndicatorEllipse2D<MyCase::value_t> ellipse1;
    IndicatorEllipse2D<MyCase::value_t> ellipse2;
    IndicatorEllipse2D<MyCase::value_t> ellipse3;
};

/// @brief Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters)
{
  using T             = MyCase::value_t;

  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const T      dx     = parameters.get<parameters::PHYS_DELTA_X >();
  const Vector origin = parameters.get<parameters::ORIGIN       >();

  IndicatorCuboid2D<T> cuboid(extent, origin);

  Mesh<T,MyCase::d> mesh(cuboid, dx, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity( { false, false } );
  return mesh;
}

void prepareGeometry( MyCase& myCase,
                      ellipse2D& ellipseCase) {
  auto& geometry = myCase.getGeometry();
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  geometry.rename( 0, 2 );
  geometry.rename( 2, 1, ellipseCase.ellipse1);
  geometry.rename( 1, 3, ellipseCase.ellipse2);
  geometry.rename( 1, 4, ellipseCase.ellipse3);

  geometry.innerClean();
  geometry.clean();
  geometry.checkForErrors();

  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase,
                    ellipse2D& ellipseCase) {
  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T                       = MyCase::value_t;
  using DESCRIPTOR              = MyCase::descriptor_t_of<names::NavierCauchy>;

  auto& geometry                = myCase.getGeometry();
  auto& lattice                 = myCase.getLattice(names::NavierCauchy{});
  auto& parameters              = myCase.getParameters();

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

  auto bulkIndicator = geometry.getMaterialIndicator( 1 );

  lattice.defineDynamics<structuralmechanics::StationaryLinearElasticityMRTBoundary>( bulkIndicator );

  setStructuralMechanicsBoundary<T, DESCRIPTOR, structuralmechanics::StationaryLinearElasticityMRTDirichletValidation>(lattice, geometry.getMaterialIndicator( 2 ), bulkIndicator, ellipseCase.ellipse1);
  setStructuralMechanicsBoundary<T, DESCRIPTOR, structuralmechanics::StationaryLinearElasticityMRTDirichletValidation>(lattice, geometry.getMaterialIndicator( 3 ), bulkIndicator, ellipseCase.ellipse2);
  setStructuralMechanicsBoundary<T, DESCRIPTOR, structuralmechanics::StationaryLinearElasticityMRTDirichletValidation>(lattice, geometry.getMaterialIndicator( 4 ), bulkIndicator, ellipseCase.ellipse3);

  lattice.setParameter<structuralmechanics::CHAR_U      >(converter.getConversionFactorDisplacement());
  lattice.setParameter<structuralmechanics::LAME_MU     >(converter.getLatticeLameMu());
  lattice.setParameter<structuralmechanics::BULK_MODULUS>(converter.getLatticeBulkModulus());

  const T theta     = T(1) / descriptors::invCs2<T, DESCRIPTOR>();

  const T omega_11  = T(1.) / (converter.getLatticeLameMu() / theta  + T(0.5));
  const T omega_d   = T(1.) / (T(2.) *  converter.getLatticeLameMu() / (T(1.) - theta) + T(0.5));
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

  {
    auto& communicator = lattice.getCommunicator(stage::PostCollide());
    communicator.template requestField<structuralmechanics::DISPLACEMENT>();
    communicator.template requestField<structuralmechanics::CAUCHY_STRESS>();
    communicator.template requestField<structuralmechanics::PREVIOUS_CELL>();
    communicator.template requestField<structuralmechanics::BARED_MOMENT_VECTOR>();
    communicator.requestOverlap(parameters.get<parameters::OVERLAP>());
    communicator.exchangeRequests();
  }

  {
    auto& communicator = lattice.getCommunicator(stage::PreStream());
    communicator.template requestField<descriptors::POPULATION>();
    communicator.requestOverlap(parameters.get<parameters::OVERLAP>());
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  OstreamManager clout(std::cout, "setInitialValues");
  clout << "Set Initial Values ..." << std::endl;

  using T                       = MyCase::value_t;

  auto& geometry                = myCase.getGeometry();
  auto& lattice                 = myCase.getLattice(names::NavierCauchy{});
  const auto& converter         = lattice.getUnitConverter();

  lattice.initialize();

  ForceField2D<T> force(converter.getLameLambda(), converter.getLameMu(), converter.getConversionFactorSolidForce());
  lattice.defineField<descriptors::FORCE>(geometry, 1, force);

  lattice.setProcessingContext(ProcessingContext::Simulation);

  clout << "Set Initial Values ... OK" << std::endl;
}

/// Update boundary values at times (and external fields, if they exist)
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
/// @note Be careful: boundary values have to be set using lattice units
void setTemporalValues(MyCase& myCase,
                       std::size_t iT) {
  // Nothing to do here
}

/// Compute simulation results at times
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT) {
  OstreamManager clout(std::cout, "getResults");

  using T                     = MyCase::value_t;
  using DESCRIPTOR            = MyCase::descriptor_t_of<names::NavierCauchy>;

  auto& geometry              = myCase.getGeometry();
  auto& lattice               = myCase.getLattice(names::NavierCauchy{});
  const auto& converter       = lattice.getUnitConverter();
  auto& parameters            = myCase.getParameters();

  const std::size_t iTlog     = parameters.get<parameters::LATTICE_LOG_ITER_T                 >();
  const std::size_t iTvtk     = parameters.get<parameters::LATTICE_VTK_ITER_T                 >();
  const std::size_t iTMax     = converter.getLatticeTime(parameters.get<parameters::PHYS_MAX_T>());
  const std::size_t convIter  = parameters.get<parameters::CONV_ITER                          >();
  const bool        converged = parameters.get<parameters::CONVERGED                          >();
  const bool        writeVTK  = parameters.get<parameters::VTK_ENABLED                        >();

  if (iT == 0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperVTMwriter2D<T> vtmWriter("dirichletEllipse");

    SuperLatticeCuboid2D<T, DESCRIPTOR>   cuboid(lattice);
    SuperLatticeRank2D<T, DESCRIPTOR>     rank(lattice);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();

    lattice.setProcessingContext(ProcessingContext::Simulation);
  }

  if (writeVTK && (iT % iTvtk == 0 || iT % iTMax == 0 || converged)) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperVTMwriter2D<T> vtmWriter("dirichletEllipse");

    ManufacturedSolutionU2D<T>                    anaPhysDisp(true);
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> anaPhysDispOnLattice(anaPhysDisp, lattice);

    ManufacturedSolutionStress2D<T>               anaPhysStress(converter.getLameMu(), converter.getLameLambda());
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> anaPhysStressOnLattice(anaPhysStress, lattice);

    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::DISPLACEMENT  > numPhysDisp(lattice, converter.getConversionFactorDisplacement());
    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::CAUCHY_STRESS > numPhysStress(lattice, converter.getConversionFactorSolidStress());

    SuperLatticePhysField2D<T, DESCRIPTOR, descriptors::FORCE> numPhysForce(lattice, converter.getConversionFactorSolidForce());

    SuperLatticeField2D<T, DESCRIPTOR, structuralmechanics::SOLID_DISTANCE_FIELD> q_ij(lattice);
    SuperLatticeField2D<T, DESCRIPTOR, structuralmechanics::BOUNDARY_COORDS_X   > bcX(lattice);
    SuperLatticeField2D<T, DESCRIPTOR, structuralmechanics::BOUNDARY_COORDS_Y   > bcY(lattice);
    SuperLatticeField2D<T, DESCRIPTOR, structuralmechanics::PREVIOUS_CELL       > prevCell(lattice);
    SuperLatticeField2D<T, DESCRIPTOR, structuralmechanics::BARED_MOMENT_VECTOR > bared_moments(lattice);

    vtmWriter.addFunctor(numPhysDisp,              "num phys disp"  );
    vtmWriter.addFunctor(anaPhysDispOnLattice,     "ana phys disp"  );

    vtmWriter.addFunctor(numPhysStress,            "num phys stress");
    vtmWriter.addFunctor(anaPhysStressOnLattice,   "ana phys stress");
    vtmWriter.addFunctor(numPhysForce,             "num force"      );

    vtmWriter.addFunctor(q_ij,                     "q_ij"           );
    vtmWriter.addFunctor(bcX,                      "bcX"            );
    vtmWriter.addFunctor(bcY,                      "bcY"            );
    vtmWriter.addFunctor(prevCell,                 "prev cell"      );
    vtmWriter.addFunctor(bared_moments,            "bared moments"  );

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

    auto indicatorF = geometry.getMaterialIndicator({ 1 });
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
    T   lInfStressResult[2]  = {T(), T()};
    int tmp[]           = {int()};

    auto bulk           = geometry.getMaterialIndicator({ 1 });

    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::CAUCHY_STRESS>  numPhysStress(lattice, converter.getConversionFactorSolidStress());
    ManufacturedSolutionStress2D<T> anaPhysStress(converter.getLameLambda(), converter.getLameMu());

    SuperRelativeErrorLinfNorm2D<T> relStressErrorLinfNorm(lattice, numPhysStress, anaPhysStress, bulk);
    relStressErrorLinfNorm(lInfStressResult, tmp);

    parameters.set<parameters::CONV_VALUE>(lInfStressResult[0]);

    lattice.setProcessingContext(ProcessingContext::Simulation);
  }

  if (converged) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    auto indicatorF = geometry.getMaterialIndicator(1);

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

void simulate(MyCase& myCase) {
  OstreamManager clout(std::cout, "simulate");

  using T = MyCase::value_t;

  auto& geometry          = myCase.getGeometry();
  auto& lattice           = myCase.getLattice(names::NavierCauchy{});
  const auto& converter   = lattice.getUnitConverter();
  auto& parameters        = myCase.getParameters();

  const std::size_t maxIt   = converter.getLatticeTime(parameters.get<parameters::PHYS_MAX_T>());
  const std::size_t convIt  = parameters.get<parameters::CONV_ITER>();

  util::ValueTracer<T> converge(convIt, parameters.get<parameters::CONVERGENCE_PRECISION>());

  clout << "Awaiting " << maxIt << " Time steps. Starting simulation..." << std::endl;
  util::Timer<T> timer(maxIt,
                       geometry.getStatistics().getNvoxel());
  timer.start();
  for (std::size_t iT = 0; iT <= maxIt; ++iT) {
    timer.update(iT);

    setTemporalValues(myCase, iT);

    getResults(myCase, timer, iT);

    lattice.collideAndStream();

    lattice.communicate();

    if (iT % convIt == 0) {
      converge.takeValue(parameters.get<parameters::CONV_VALUE>(), true);
    }

    if (converge.hasConverged()) {
      parameters.set<parameters::CONVERGED>(true);
      getResults(myCase, timer, iT);

      clout << "Simulation converged." << std::endl;
      break;
    }

    if (iT == maxIt) {
      parameters.set<parameters::CONVERGED>(true);
      getResults(myCase, timer, iT);

      clout << "Simulation not converged. Writing error at last time step." << std::endl;
      break;
    }
  }

  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[]) {
  // === 1st Step: Initialization ===
  OstreamManager clout(std::cout, "main");

  initialize(&argc, &argv);

  singleton::directories().setOutputDir("./tmp/");

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<SOLID_CASE               >( "continua"     );
    myCaseParameters.set<RESOLUTION               >(    40          );
    myCaseParameters.set<PHYS_CHAR_LENGTH         >(     1.5        );
    myCaseParameters.set<DOMAIN_EXTENT            >(    {1.5, 1.5}  );
    myCaseParameters.set<YOUNGS_MODULUS           >(     0.1        );
    myCaseParameters.set<PHYS_DENSITY             >(     1          );
    myCaseParameters.set<POISSON_RATIO            >(     0.7        );
    myCaseParameters.set<DAMPING_FACTOR           >(     1.0        );
    myCaseParameters.set<PHYS_MAX_T               >(  1000          );
    myCaseParameters.set<LATTICE_LOG_ITER_T       >(   100          );
    myCaseParameters.set<LATTICE_VTK_ITER_T       >(  1000          );
    myCaseParameters.set<CONVERGENCE_PRECISION    >(     1e-7       );
    myCaseParameters.set<CONV_ITER                >(    10          );
    myCaseParameters.set<CONV_VALUE               >(     1          );
    myCaseParameters.set<VTK_ENABLED              >(  true          );
  }

  myCaseParameters.set<parameters::PHYS_DELTA_X>([](MyCase::ParametersD& p) {
    return p.get<parameters::DOMAIN_EXTENT>()[0] / p.get<parameters::RESOLUTION>();
  });
  myCaseParameters.set<parameters::PHYS_DELTA_T>(
    myCaseParameters.get<parameters::PHYS_DELTA_X>() * myCaseParameters.get<parameters::PHYS_DELTA_X>()
  );
  myCaseParameters.set<parameters::ORIGIN>([](MyCase::ParametersD& p) -> Vector<MyCase::value_t, MyCase::d>{
    return {p.get<parameters::PHYS_DELTA_X>(), p.get<parameters::PHYS_DELTA_X>()};
  });
  myCaseParameters.set<parameters::PHYS_CHAR_DISPLACEMENT>(
    myCaseParameters.get<parameters::PHYS_DELTA_X>()
  );


  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  ellipse2D ellipseCase;
  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase, ellipseCase);

  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase, ellipseCase);

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);
}
