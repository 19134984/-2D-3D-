/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2026 Florian Kaiser
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

#include "../../manufactured_solution_fields.h"

using MyCase = Case<
  names::NavierCauchy,
  Lattice<
    FLOATING_POINT_TYPE,
    descriptors::VECTOR<
      descriptors::D2Q4<
        descriptors::POROSITY,
        descriptors::FORCE,
        structuralmechanics::DISPLACEMENT,
        structuralmechanics::VELOCITY,
        structuralmechanics::CAUCHY_STRESS,
        structuralmechanics::U_ZERO_N<5>,
        structuralmechanics::VLBM_SOLID_DISTANCE_FIELD,
        structuralmechanics::BOUNDARY_COORDS_X,
        structuralmechanics::BOUNDARY_COORDS_Y,
        structuralmechanics::NEUMANN_NORMAL_X,
        structuralmechanics::NEUMANN_NORMAL_Y
      >, 5
    >
  >
>;

namespace olb::parameters {
  struct SCALING_FACTOR_DT : public descriptors::FIELD_BASE<1> { };
  struct MS_AMPLITUDE      : public descriptors::FIELD_BASE<1> { };
  struct MS_WAVENUMBER     : public descriptors::FIELD_BASE<1> { };
}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<names::NavierCauchy>;

  const Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector origin = parameters.get<parameters::ORIGIN       >();
  const T      dx     = parameters.get<parameters::PHYS_DELTA_X >();

  Vector<T, DESCRIPTOR::d> motherOrigin = origin - T(0.5) * dx;
  Vector<T, DESCRIPTOR::d> motherExtend = extend + T(1.5) * dx;

  IndicatorCuboid2D<T> cuboid(motherExtend, motherOrigin);

  Mesh<T,MyCase::d> mesh(cuboid, dx, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity( { false, false } );
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<names::NavierCauchy>;

  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  const Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector origin = parameters.get<parameters::ORIGIN       >();
  const T          dx = parameters.get<parameters::PHYS_DELTA_X >();

  IndicatorCuboid2D<T> cuboid(extend, origin);

  geometry.rename(0, 3);

  Vector<T, DESCRIPTOR::d> extendBulk = extend - dx;
  IndicatorCuboid2D<T>     bulk(extendBulk, origin);

  Vector<T, DESCRIPTOR::d> originBottom = origin - dx;

  Vector<T, DESCRIPTOR::d> extendBottom = Vector<T, DESCRIPTOR::d>{extend[0], origin[1]} + dx;
  IndicatorCuboid2D<T>     bottom(extendBottom, originBottom);

  geometry.rename(3, 1, bulk);
  geometry.rename(3, 2, bottom);

  geometry.clean();

  geometry.checkForErrors();
  geometry.getStatistics().print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using DESCRIPTOR = MyCase::descriptor_t_of<names::NavierCauchy>;
  using T = MyCase::value_t;

  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(names::NavierCauchy{});
  auto& parameters = myCase.getParameters();

  auto bulkIndicator = geometry.getMaterialIndicator( 1 );
  lattice.defineDynamics<structuralmechanics::ElastodynamicsVectorLBM>(bulkIndicator);

  const Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector origin = parameters.get<parameters::ORIGIN       >();
  const T          dx = parameters.get<parameters::PHYS_DELTA_X >();

  Vector<T, DESCRIPTOR::d> extendBulk = extend - dx;
  IndicatorCuboid2D<T>     bulk(extendBulk, origin);

  auto dirichletIndicator = geometry.getMaterialIndicator(2);
  auto neumannIndicator   = geometry.getMaterialIndicator(3);

  setStructuralMechanicsBoundary<T, DESCRIPTOR, structuralmechanics::ElastodynamicsVectorLBMNeumannValidation>(
    lattice, neumannIndicator, bulkIndicator, bulk
  );
  setStructuralMechanicsBoundary<T, DESCRIPTOR, structuralmechanics::ElastodynamicsVectorLBMDirichletValidation>(
    lattice, dirichletIndicator, bulkIndicator, bulk
  );

  const T charL               = parameters.get<parameters::PHYS_CHAR_LENGTH            >();
  const T charT               = parameters.get<parameters::PHYS_CHAR_TIME              >();
  const T charV               = parameters.get<parameters::PHYS_CHAR_VELOCITY          >();
  const T E                   = parameters.get<parameters::YOUNGS_MODULUS              >();
  const T nu                  = parameters.get<parameters::POISSON_RATIO               >();
  const T rho0                = parameters.get<parameters::PHYS_DENSITY                >();
  const T omega               = parameters.get<parameters::LATTICE_RELAXATION_FREQUENCY>();
  const std::string solidCase = parameters.get<parameters::SOLID_CASE                  >();

  const T c_mu     = E / (T(2.) * (T(1.) + nu));
  const T c_K      = E * (T(1.0) - nu) / (rho0 * (T(1.0) + nu) * (T(1.0) - T(2.0) * nu));

  const T dtFactor = parameters.get<parameters::SCALING_FACTOR_DT>();
  const T dt       = dtFactor * T(1) / util::sqrt(descriptors::invCs2<T, DESCRIPTOR>()) * dx / util::sqrt(c_mu * c_mu + c_K * c_K);

  lattice.setUnitConverter<ElastodynamicsVectorLBMUnitConverter<T, DESCRIPTOR>>(
    dx,
    dt,
    charL,
    charT,
    charV,
    E,
    nu,
    rho0,
    solidCase
  );
  const auto& converter = lattice.getUnitConverter();
  converter.print();

  lattice.setParameter<descriptors::DT                   >(converter.getConversionFactorTime()    );
  lattice.setParameter<descriptors::OMEGA                >(omega                                  );
  lattice.setParameter<structuralmechanics::C_K          >(converter.getLatticeCK()               );
  lattice.setParameter<structuralmechanics::C_MU         >(converter.getLatticeCMu()              );
  lattice.setParameter<structuralmechanics::LATTICE_SPEED>(converter.getLatticeSpeed()            );
  lattice.setParameter<structuralmechanics::LAME_MU      >(converter.getLatticeLameMu()           );
  lattice.setParameter<structuralmechanics::LAME_LAMBDA  >(converter.getLatticeLameLambda()       );
  lattice.setParameter<structuralmechanics::CHAR_V       >(converter.getConversionFactorVelocity());
  lattice.setParameter<structuralmechanics::MS_AMPLITUDE >(parameters.get<parameters::MS_AMPLITUDE>());
  lattice.setParameter<structuralmechanics::MS_WAVENUMBER>(parameters.get<parameters::MS_WAVENUMBER>());

  {
    auto& communicator = lattice.getCommunicator(stage::PostCollide());
    communicator.requestField<structuralmechanics::VELOCITY>();
    communicator.requestField<structuralmechanics::DISPLACEMENT>();
    communicator.requestField<structuralmechanics::CAUCHY_STRESS>();
    communicator.requestOverlap(parameters.get<parameters::OVERLAP>());
    communicator.exchangeRequests();
  }

  {
    auto& communicator = lattice.getCommunicator(stage::PreStream());
    communicator.requestField<descriptors::POPULATION>();
    communicator.requestOverlap(parameters.get<parameters::OVERLAP>());
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  OstreamManager clout(std::cout, "setInitialValues");
  clout << "Set Initial Values ..." << std::endl;

  using DESCRIPTOR = MyCase::descriptor_t_of<names::NavierCauchy>;
  using T = MyCase::value_t;

  auto& geometry        = myCase.getGeometry();
  auto& lattice         = myCase.getLattice(names::NavierCauchy{});
  const auto& converter = lattice.getUnitConverter();

  lattice.initialize();

  auto bulkIndicator = geometry.getMaterialIndicator( 1 );

  AnalyticalConst2D<T, T> poro(T(10));
  lattice.defineField<descriptors::POROSITY>(bulkIndicator, poro);

  AnalyticalVelocity<T>     v0(T(0), converter.getConversionFactorVelocity());
  AnalyticalDisplacement<T> u0(T(0), converter.getConversionFactorDisplacement());
  AnalyticalForce<T>        f0(T(0), converter.getPhysDensity(), converter.getLameLambda(), converter.getLameMu(), converter.getConversionFactorSolidForce());

  lattice.defineField<structuralmechanics::VELOCITY     >(bulkIndicator, v0);
  lattice.defineField<structuralmechanics::DISPLACEMENT >(bulkIndicator, u0);
  lattice.defineField<descriptors::FORCE                >(bulkIndicator, f0);

  lattice.setProcessingContext<Array<descriptors::POROSITY            >>(ProcessingContext::Simulation);
  lattice.setProcessingContext<Array<descriptors::FORCE               >>(ProcessingContext::Simulation);
  lattice.setProcessingContext<Array<structuralmechanics::DISPLACEMENT>>(ProcessingContext::Simulation);
  lattice.setProcessingContext<Array<structuralmechanics::VELOCITY    >>(ProcessingContext::Simulation);

  lattice.communicate();

  SuperLatticeCoupling ini1(
    structuralmechanics::IniElastodynamicsUZero<T, DESCRIPTOR>{},
    names::NavierStokes{},
    lattice
  );
  ini1.setParameter<structuralmechanics::IniElastodynamicsUZero<T, DESCRIPTOR>::C_K>(converter.getLatticeCK());
  ini1.setParameter<structuralmechanics::IniElastodynamicsUZero<T, DESCRIPTOR>::C_MU>(converter.getLatticeCMu());

  SuperLatticeCoupling ini2(
    structuralmechanics::IniPopulationsElastodynamics<T, DESCRIPTOR>{},
    names::NavierStokes{},
    lattice
  );
  ini2.setParameter<structuralmechanics::IniPopulationsElastodynamics<T, DESCRIPTOR>::SECOND_ORDER>(T(1));
  ini2.setParameter<structuralmechanics::IniPopulationsElastodynamics<T, DESCRIPTOR>::C_K>(converter.getLatticeCK());
  ini2.setParameter<structuralmechanics::IniPopulationsElastodynamics<T, DESCRIPTOR>::C_MU>(converter.getLatticeCMu());
  ini2.setParameter<structuralmechanics::IniPopulationsElastodynamics<T, DESCRIPTOR>::LATTICE_SPEED>(converter.getLatticeSpeed());

  struct IniExchange1 {};
  struct IniExchange2 {};
  struct IniExchange3 {};

  auto& iniComm1 = lattice.getCommunicator(IniExchange1{});
  iniComm1.requestField<descriptors::POROSITY>();
  iniComm1.requestField<structuralmechanics::DISPLACEMENT>();
  iniComm1.requestOverlap(lattice.getOverlap());
  iniComm1.exchangeRequests();

  auto& iniComm2 = lattice.getCommunicator(IniExchange2{});
  iniComm2.requestField<descriptors::POROSITY>();
  iniComm2.requestField<structuralmechanics::U_ZERO_N<DESCRIPTOR::n>>();
  iniComm2.requestOverlap(lattice.getOverlap());
  iniComm2.exchangeRequests();

  auto& iniComm3 = lattice.getCommunicator(IniExchange3{});
  iniComm3.requestField<structuralmechanics::DISPLACEMENT>();
  iniComm3.requestField<structuralmechanics::VELOCITY>();
  iniComm3.requestField<descriptors::POPULATION>();
  iniComm3.requestOverlap(lattice.getOverlap());
  iniComm3.exchangeRequests();

  // 5. Execute sequence
  iniComm1.communicate();
  ini1.apply();

  iniComm2.exchangeRequests();
  iniComm2.communicate();

  ini2.apply();

  iniComm3.communicate();

  clout << "Set Initial Values ... OK" << std::endl;
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT) {
  using T = MyCase::value_t;

  auto& geometry        = myCase.getGeometry();
  auto& lattice         = myCase.getLattice(names::NavierCauchy{});
  const auto& converter = lattice.getUnitConverter();

  AnalyticalForce<T> forceSol(
    converter.getPhysTime(iT),
    converter.getPhysDensity(),
    converter.getLameLambda(),
    converter.getLameMu(),
    converter.getConversionFactorSolidForce()
  );

  lattice.setProcessingContext(ProcessingContext::Evaluation);
  lattice.defineField<descriptors::FORCE>(geometry, 1, forceSol);
  lattice.setParameter<descriptors::MAGIC>(converter.getPhysTime(iT) + T(0.5) * converter.getConversionFactorTime());
  lattice.setProcessingContext(ProcessingContext::Simulation);
}

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

  const std::size_t iTMax     = converter.getLatticeTime(parameters.get<parameters::PHYS_MAX_T      >());
  const std::size_t iTlog     = converter.getLatticeTime(parameters.get<parameters::PHYS_STAT_ITER_T>());
  const std::size_t iTvtk     = converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T >());
  const bool        writeVTK  = parameters.get<parameters::VTK_ENABLED                              >();

  SuperVTMwriter2D<T> vtmWriter("combinedPlate");

  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  if (iT % iTlog == 0 || iT % iTMax == 0) {
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
    timer.update(iT);
    timer.printStep();
  }

  if (writeVTK && (iT % iTvtk == 0 || iT % iTMax == 0)) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperGeometryF<T,DESCRIPTOR::d> materials(geometry);

    SuperLatticeField2D<T, DESCRIPTOR, descriptors::POPULATION> population(lattice);
    SuperLatticePhysField2D<T, DESCRIPTOR, descriptors::FORCE > physForce(lattice, converter.getConversionFactorSolidForce());

    AnalyticalVelocity<T>     veloAnaPhys(converter.getPhysTime(iT));
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> veloAnaPhysOnLattice(veloAnaPhys, lattice);

    AnalyticalDisplacement<T> dispAnaPhys(converter.getPhysTime(iT));
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> dispAnaPhysOnLattice(dispAnaPhys, lattice);

    AnalyticalCauchyStress<T> stressAnaPhys(converter.getPhysTime(iT), converter.getLameLambda(), converter.getLameMu());
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> stressAnaPhysOnLattice(stressAnaPhys, lattice);

    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::VELOCITY      > veloNumPhys(lattice, converter.getConversionFactorVelocity());
    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::DISPLACEMENT  > dispNumPhys(lattice, converter.getConversionFactorDisplacement());
    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::CAUCHY_STRESS > stressNumPhys(lattice, converter.getConversionFactorSolidStress());

    auto indicatorF = geometry.getMaterialIndicator(1);

    vtmWriter.addFunctor(materials,               "materials");
    vtmWriter.addFunctor(population,              "population");
    vtmWriter.addFunctor(physForce,               "phys force");

    vtmWriter.addFunctor(veloNumPhys,             "num phys vel");
    vtmWriter.addFunctor(veloAnaPhysOnLattice,    "ana phys vel");

    vtmWriter.addFunctor(dispNumPhys,             "num phys disp");
    vtmWriter.addFunctor(dispAnaPhysOnLattice,    "ana phys disp");

    vtmWriter.addFunctor(stressNumPhys,           "num phys stress");
    vtmWriter.addFunctor(stressAnaPhysOnLattice,  "ana phys stress");

    vtmWriter.write(iT);

    T   l2VResult[2]        = {T(), T()};
    T   lInfVResult[2]      = {T(), T()};
    T   l2UResult[2]        = {T(), T()};
    T   lInfUResult[2]      = {T(), T()};
    T   l2StressResult[2]   = {T(), T()};
    T   lInfStressResult[2] = {T(), T()};
    int tmp[]               = {int()};

    SuperRelativeErrorL2Norm2D<T>   relVErrorL2Norm(lattice, veloNumPhys, veloAnaPhys, indicatorF);
    SuperRelativeErrorLinfNorm2D<T> relVErrorLinfNorm(lattice, veloNumPhys, veloAnaPhys, indicatorF);
    relVErrorL2Norm(l2VResult, tmp);
    relVErrorLinfNorm(lInfVResult, tmp);

    SuperRelativeErrorL2Norm2D<T>   relUErrorL2Norm(lattice, dispNumPhys, dispAnaPhys, indicatorF);
    SuperRelativeErrorLinfNorm2D<T> relUErrorLinfNorm(lattice, dispNumPhys, dispAnaPhys, indicatorF);
    relUErrorL2Norm(l2UResult, tmp);
    relUErrorLinfNorm(lInfUResult, tmp);

    SuperRelativeErrorL2Norm2D<T>   relStressErrorL2Norm(lattice, stressNumPhys, stressAnaPhys, indicatorF);
    SuperRelativeErrorLinfNorm2D<T> relStressErrorLinfNorm(lattice, stressNumPhys, stressAnaPhys, indicatorF);
    relStressErrorL2Norm(l2StressResult, tmp);
    relStressErrorLinfNorm(lInfStressResult, tmp);

    clout << "N\t"     << "L2 V Error\t" << "LInf V Error\t" << "L2 U Error\t" << "LInf U Error\t" << "L2 Stress Error\t" << "LInf Stress Error" <<std::endl;
    clout <<  converter.getResolution() <<"\t" <<  l2VResult[0]  << "\t" << lInfVResult[0] << "\t" << l2UResult[0] << "\t" << lInfUResult[0] << "\t" << l2StressResult[0] << "\t" << lInfStressResult[0] << std::endl;

    lattice.setProcessingContext(ProcessingContext::Simulation);
  }

  if (iT == iTMax) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    auto indicatorF = geometry.getMaterialIndicator({1});

    T   l2VResult[2]        = {T(), T()};
    T   lInfVResult[2]      = {T(), T()};
    T   l2UResult[2]        = {T(), T()};
    T   lInfUResult[2]      = {T(), T()};
    T   l2StressResult[2]   = {T(), T()};
    T   lInfStressResult[2] = {T(), T()};
    int tmp[]               = {int()};

    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::VELOCITY      > veloNumPhys(lattice, converter.getConversionFactorVelocity());
    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::DISPLACEMENT  > dispNumPhys(lattice, converter.getConversionFactorDisplacement());
    SuperLatticePhysField2D<T, DESCRIPTOR, structuralmechanics::CAUCHY_STRESS > stressNumPhys(lattice, converter.getConversionFactorSolidStress());

    AnalyticalVelocity<T>     veloAnaPhys(converter.getPhysTime(iT));

    AnalyticalDisplacement<T> dispAnaPhys(converter.getPhysTime(iT));

    AnalyticalCauchyStress<T> stressAnaPhys(converter.getPhysTime(iT), converter.getLameLambda(), converter.getLameMu());

    SuperRelativeErrorL2Norm2D<T>   relVErrorL2Norm(lattice, veloNumPhys, veloAnaPhys, indicatorF);
    SuperRelativeErrorLinfNorm2D<T> relVErrorLinfNorm(lattice, veloNumPhys, veloAnaPhys, indicatorF);
    relVErrorL2Norm(l2VResult, tmp);
    relVErrorLinfNorm(lInfVResult, tmp);

    SuperRelativeErrorL2Norm2D<T>   relUErrorL2Norm(lattice, dispNumPhys, dispAnaPhys, indicatorF);
    SuperRelativeErrorLinfNorm2D<T> relUErrorLinfNorm(lattice, dispNumPhys, dispAnaPhys, indicatorF);
    relUErrorL2Norm(l2UResult, tmp);
    relUErrorLinfNorm(lInfUResult, tmp);

    SuperRelativeErrorL2Norm2D<T>   relStressErrorL2Norm(lattice, stressNumPhys, stressAnaPhys, indicatorF);
    SuperRelativeErrorLinfNorm2D<T> relStressErrorLinfNorm(lattice, stressNumPhys, stressAnaPhys, indicatorF);
    relStressErrorL2Norm(l2StressResult, tmp);
    relStressErrorLinfNorm(lInfStressResult, tmp);

    CSV<T> csvWriterV("vel_err", ',', {"N", "L2", "LInf"}, ".csv", true);
    CSV<T> csvWriterU("disp_err", ',', {"N", "L2", "LInf"}, ".csv", true);
    CSV<T> csvWriterSigma("stress_err", ',', {"N", "L2", "LInf"}, ".csv", true);

    csvWriterV.writeDataFile(converter.getResolution(), {l2VResult[0], lInfVResult[0]});
    csvWriterU.writeDataFile(converter.getResolution(), {l2UResult[0], lInfUResult[0]});
    csvWriterSigma.writeDataFile(converter.getResolution(), {l2StressResult[0], lInfStressResult[0]});

    lattice.setProcessingContext(ProcessingContext::Simulation);
  }
}

void simulate(MyCase& myCase) {
  OstreamManager clout(std::cout, "simulate");
  clout << "Starting simulation ..." << std::endl;

  using T = MyCase::value_t;

  auto& geometry            = myCase.getGeometry();
  auto& lattice             = myCase.getLattice(names::NavierCauchy{});
  auto& parameters          = myCase.getParameters();

  const T physMaxT          = parameters.get<parameters::PHYS_MAX_T>();
  const std::size_t iTmax   = lattice.getUnitConverter().getLatticeTime(physMaxT);

  util::Timer<T> timer(iTmax, geometry.getStatistics().getNvoxel());
  timer.start();
  for (std::size_t iT = 0; iT <= iTmax; ++iT) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    lattice.collideAndStream();

    lattice.communicate();

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, timer, iT);
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

    myCaseParameters.set<RESOLUTION                   >(   40           );
    myCaseParameters.set<SOLID_CASE                   >( "continua"     );
    myCaseParameters.set<YOUNGS_MODULUS               >(    1           );
    myCaseParameters.set<PHYS_DENSITY                 >(    1           );
    myCaseParameters.set<POISSON_RATIO                >(    0.4         );
    myCaseParameters.set<PHYS_CHAR_LENGTH             >(    1.0         );
    myCaseParameters.set<PHYS_CHAR_TIME               >(    1.0         );
    myCaseParameters.set<PHYS_CHAR_VELOCITY           >(    1.0         );
    myCaseParameters.set<PHYS_CHAR_DISPLACEMENT       >(    1.0         );
    myCaseParameters.set<LATTICE_RELAXATION_FREQUENCY >(    2.0         );
    myCaseParameters.set<DOMAIN_EXTENT                >(  { 1.0, 1.0 }  );
    myCaseParameters.set<ORIGIN                       >(  { 0.0, 0.0 }  );
    myCaseParameters.set<PHYS_MAX_T                   >(    1.95        );
    myCaseParameters.set<PHYS_STAT_ITER_T             >(    0.1         );
    myCaseParameters.set<PHYS_VTK_ITER_T              >(    0.1         );
    myCaseParameters.set<SCALING_FACTOR_DT            >(    0.4         );
    myCaseParameters.set<MS_AMPLITUDE                 >(    1           );
    myCaseParameters.set<MS_WAVENUMBER                >(    2 * M_PI    );
    myCaseParameters.set<VTK_ENABLED                  >( true           );
  }

  myCaseParameters.set<parameters::PHYS_DELTA_X>(
    myCaseParameters.get<parameters::PHYS_CHAR_LENGTH>() / myCaseParameters.get<parameters::RESOLUTION>()
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
