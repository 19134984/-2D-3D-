/*  Lattice Boltzmann sample, written in C++, using the OpenLBlibrary
 *
 *  Copyright (C) 2023 Shota Ito
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

// This app is used for the validation of the ADE-LBM and of the Dirchlet and Zero-Gradient BC regarding the transported concentration.
// The simulation setup and analytic solution is taken from Krüger et al. 2017 page 324 (doi: https://doi.org/10.1007/978-3-319-44649-3).

#include <olb.h>
#include <cmath>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;

namespace olb::parameters {
  struct HEIGHT : public descriptors::FIELD_BASE<1> { };
}

using MyCase = Case<NavierStokes, Lattice<double, D3Q7<VELOCITY,OMEGA>>>;

Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  T height = parameters.get<parameters::DOMAIN_EXTENT>()[2];
  T width = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  T length = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const Vector<T,3> origin;
  const Vector<T,3> extend(length, height, width);
  IndicatorCuboid3D<T> cuboid(extend, origin);

  Mesh<T,MyCase::d> mesh(cuboid, parameters.get<parameters::PHYS_DELTA_X>(), singleton::mpi().getSize());

  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({false, false, true});

  return mesh;
 }

// prepare numerical domain
void prepareGeometry(MyCase& myCase){
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& superGeometry   = myCase.getGeometry();

  // === MATERIAL NUMBERS === //
  // solid: MN = 0            //
  // bulk: MN = 1             //
  // inlet: MN = 2            //
  // outlet: MN = 3           //
  // topWall: MN = 4          //
  // bottomWall: MN = 5       //
  // ======================== //

  superGeometry.rename(0, 1);
  T eps = parameters.get<parameters::PHYS_DELTA_X>() * T(0.5);
  T height = parameters.get<parameters::DOMAIN_EXTENT>()[2];
  T width = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  T length = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  T dx = parameters.get<parameters::PHYS_DELTA_X>();

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

void prepareLattice(MyCase& myCase){
  using T = MyCase::value_t;
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  auto& superGeometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  auto& adLattice = myCase.getLattice(NavierStokes{});

  adLattice.setUnitConverter<AdeUnitConverter<T, MyCase::descriptor_t_of<NavierStokes>>>(
    (T)   parameters.get<parameters::PHYS_DELTA_X>(),                                                                               // physDeltaX
    (T)   util::pow(parameters.get<parameters::HEIGHT>() / parameters.get<parameters::RESOLUTION>(), 2),                            // physDeltaY
    (T)   parameters.get<parameters::HEIGHT>(),                                                                                     // charPhysLength: reference length of simulation geometry
    (T)   parameters.get<parameters::PECLET>() * parameters.get<parameters::PHYS_DIFFUSIVITY>() / parameters.get<parameters::HEIGHT>(),  // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   parameters.get<parameters::PHYS_DIFFUSIVITY>(),                                                                                // physDiffusiviy
    (T)   parameters.get<parameters::PHYS_CHAR_DENSITY>()                                                                           // physDensity: physical density in __kg / m^3__
  );

  // define dynamics and boundary conditions
  const auto bulkIndicator = superGeometry.getMaterialIndicator({1, 2, 3, 4, 5});

  adLattice.template defineDynamics<NoDynamics>(superGeometry, 0);
  adLattice.template defineDynamics<AdvectionDiffusionBGKdynamics>(bulkIndicator);

  boundary::set<boundary::AdvectionDiffusionDirichlet>(adLattice, superGeometry.getMaterialIndicator({2}));
  setZeroGradientBoundary<T, MyCase::descriptor_t_of<NavierStokes>>(adLattice, superGeometry.getMaterialIndicator({3}));
  setZeroGradientBoundary<T, MyCase::descriptor_t_of<NavierStokes>>(adLattice, superGeometry.getMaterialIndicator({4}));
  boundary::set<boundary::AdvectionDiffusionDirichlet>(adLattice, superGeometry.getMaterialIndicator({5}));

  clout << "Prepare Lattice ... OK" << std::endl;
  adLattice.getUnitConverter().print();
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  OstreamManager clout(std::cout, "setInitialValues");
  clout << "Setting initial values ..." << std::endl;

  auto& superGeometry = myCase.getGeometry();
  auto& adLattice = myCase.getLattice(NavierStokes{});
  const auto& converter = adLattice.getUnitConverter();

  // set the convective velocity field
  Vector<T,3> u(adLattice.getUnitConverter().getCharPhysVelocity(), 0.0, 0.0);
  fields::setVelocity<descriptors::VELOCITY>(adLattice,
                                             superGeometry.getMaterialIndicator({1, 2, 3, 4, 5}),
                                             u);

  momenta::setDensity(adLattice, superGeometry.getMaterialIndicator({1,2,3,4}), 1e-8);
  momenta::setDensity(adLattice, superGeometry.getMaterialIndicator(5), 1);

  adLattice.template setParameter<descriptors::OMEGA>(
    converter.getLatticeRelaxationFrequency());

  adLattice.initialize();

  clout << "Setting initial values OK" << std::endl;
}

// compute error to analytical solution
void errorCalc(MyCase& myCase, MyCase::value_t& relError){
  using T = MyCase::value_t;
  OstreamManager clout(std::cout, "Error2Analytic");
  auto& superGeometry = myCase.getGeometry();
  auto& adLattice = myCase.getLattice(NavierStokes{});
  auto& converter = adLattice.getUnitConverter();


  T result[3] = {T(), T(), T()};
  int tmp[] = {int()};

  // Compute L2-norm error to the analytic concentration field
  SuperLatticeDensity3D<T,MyCase::descriptor_t_of<NavierStokes>> concentration(adLattice);
 AnalyticalFfromCallableF<MyCase::d, T, T> analyticSol([&] (Vector<T,3> input)->Vector<T,1>{

    if (input[0] >= 0.0 && input[1] >= 0.0) {
      return erfc(input[1] / sqrt(4 * converter.getPhysDiffusivity() * input[0] / converter.getCharPhysVelocity()));
    }
    else {
      return 0.0;
    }
  });

  auto indicatorF = superGeometry.getMaterialIndicator({1});

  SuperRelativeErrorL2Norm3D<T> relCErrorL2Norm(concentration, analyticSol, indicatorF);
  SuperAbsoluteErrorL2Norm3D<T> absCErrorL2Norm(concentration, analyticSol, indicatorF);

  relCErrorL2Norm(result, tmp);
  relError = result[0];
}

// compute results and generate output
void getResults(MyCase& myCase,
                std::size_t iT,
                util::Timer<MyCase::value_t>& timer,
                util::ValueTracer<MyCase::value_t> converge,
                MyCase::value_t& relError)
{
  using T = MyCase::value_t;
  OstreamManager clout(std::cout, "Results");
  auto& adLattice = myCase.getLattice(NavierStokes{});
  auto& converter = adLattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  // generate vtm files for visualisation
  SuperVTMwriter3D<T> vtmWriter("convectedPlate");
  SuperLatticeDensity3D<T,MyCase::descriptor_t_of<NavierStokes>> concentration(adLattice);
  SuperLatticePhysField3D<T,MyCase::descriptor_t_of<NavierStokes>,VELOCITY> velocity(adLattice, adLattice.getUnitConverter().getConversionFactorVelocity());

  AnalyticalFfromCallableF<MyCase::d, T, T> analyticSol([&] (Vector<T,3> input)->Vector<T,1>{
    if (input[0] >= 0.0 && input[1] >= 0.0) {
      return erfc(input[1] / sqrt(4 * converter.getPhysDiffusivity() * input[0] / converter.getCharPhysVelocity()));
    }
    else {
      return 0.0;
    }
  });

  SuperLatticeFfromAnalyticalF3D<T,MyCase::descriptor_t_of<NavierStokes>> analyticalC(analyticSol, adLattice);

  concentration.getName() = "C";
  velocity.getName() = "u0";
  analyticalC.getName() = "analyticC";
  vtmWriter.addFunctor(concentration);
  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(analyticalC);

  const std::size_t vtkIter = adLattice.getUnitConverter().getLatticeTime(parameters.get<parameters::PHYS_MAX_T>() / 100);
  const std::size_t statIter = adLattice.getUnitConverter().getLatticeTime(parameters.get<parameters::PHYS_MAX_T>() / 100);

  if (iT == 0) {
    // write the geometry
    vtmWriter.createMasterFile();

    // write inital field
    vtmWriter.write(iT);
  }

  if (iT%statIter == 0) {
    adLattice.setProcessingContext(ProcessingContext::Evaluation);

    // print L2 error
    errorCalc(myCase, relError);

    // Timer console output
    timer.update(iT);
    timer.printStep();

    // Lattice statistics console output
    adLattice.getStatistics().print(iT, adLattice.getUnitConverter().getPhysTime(iT));
    clout << "concentration-L2-relative-error = " << relError << std::endl;

    converge.takeValue( relError, false );
  }
  // write the vtk files
  if (iT % vtkIter == 0 && iT > 0) {
    vtmWriter.write(iT);
  }
}


void simulate(MyCase& myCase){
  OstreamManager clout(std::cout,"simulate");
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& adLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();

  util::Timer<T> timer(adLattice.getUnitConverter().getLatticeTime(parameters.get<parameters::PHYS_MAX_T>()), superGeometry.getStatistics().getNvoxel());

  // Convergence check
  util::ValueTracer<T> converge(adLattice.getUnitConverter().getLatticeTime(parameters.get<parameters::PHYS_CONVERGENCE_CHECK_T>()), parameters.get<parameters::RESIDUUM>());
  timer.start();
  T relError{};

  for (std::size_t iT = 0; iT < adLattice.getUnitConverter().getLatticeTime(parameters.get<parameters::PHYS_MAX_T>()); ++iT) {

    // === Collide and stream ===
    adLattice.collideAndStream();

    // === Computation and output of the results ===
    getResults(myCase, iT, timer,  converge, relError);
    if (converge.hasConverged()) {
      clout << "Simulation stops." << std::endl;
      break;
    }
  }

  timer.stop();
  timer.printSummary();


}

int main(int argc, char* argv[]) {
    /// === 1st Step: Initialization ===
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

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

    myCaseParameters.set<PHYS_DELTA_X>([] (MyCase::ParametersD& p) {
      return p.get<HEIGHT>()/p.get<RESOLUTION>();
    });
    myCaseParameters.set<DOMAIN_EXTENT>([] (MyCase::ParametersD& p) {
      return Vector{p.get<RESOLUTION>()*p.get<PHYS_DELTA_X>()*10, p.get<PHYS_DELTA_X>(), p.get<HEIGHT>()};
    });
    myCaseParameters.set<PHYS_CONVERGENCE_CHECK_T>([] (MyCase::ParametersD& p) {
      return p.get<PHYS_MAX_T>()*0.001;
    });

  }

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
