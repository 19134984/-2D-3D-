/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2016 Robin Trunk, Mathias J. Krause
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

/* bifurcation3d.cpp:
 * This example examines a steady particulate flow past a bifurcation. At the inlet,
 * an inflow condition with grad_n u = 0 and rho = 1 is implemented.
 * At both outlets, a Poiseuille profile is imposed on the velocity.
 * After a start time, particles are put into the bifurcation by imposing
 * a inflow condition with rho = 1 on the second euler phase at the inlet.
 * The particles are simulated as a continuum with a advection-diffusion equation
 * and experience a stokes drag force.
 *
 * A publication using the same geometry can be found here:
 * http://link.springer.com/chapter/10.1007/978-3-642-36961-2_5
 *  *
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;

using MyCase = Case<
  NavierStokes,   Lattice<double, descriptors::D3Q19<>>,
  Concentration0, Lattice<double, descriptors::D3Q7<descriptors::VELOCITY, descriptors::VELOCITY2>>
>;

// === Step 1: Declarations ===
namespace olb::parameters {

struct INLET_RADIUS : public descriptors::FIELD_BASE<1> {};
struct INLET_CENTER : public descriptors::FIELD_BASE<0, 1> {};
struct INLET_NORMAL : public descriptors::FIELD_BASE<0, 1> {};
struct OUTLET_RADIUS0 : public descriptors::FIELD_BASE<1> {};
struct OUTLET_CENTER0 : public descriptors::FIELD_BASE<0, 1> {};
struct OUTLET_NORMAL0 : public descriptors::FIELD_BASE<0, 1> {};
struct OUTLET_RADIUS1 : public descriptors::FIELD_BASE<1> {};
struct OUTLET_CENTER1 : public descriptors::FIELD_BASE<0, 1> {};
struct OUTLET_NORMAL1 : public descriptors::FIELD_BASE<0, 1> {};

}

Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters)
{
  using T                 = MyCase::value_t;
  const T      physDeltaX = 2 * parameters.get<parameters::INLET_RADIUS>() / parameters.get<parameters::RESOLUTION>();
  STLreader<T> stlReader("../bifurcation3d.stl", physDeltaX);
  IndicatorLayer3D<T> extendedDomain(stlReader, physDeltaX);
  std::size_t         noOfCuboids = util::max(20, singleton::mpi().getSize());
  Mesh<T, MyCase::d>  mesh(extendedDomain, physDeltaX, noOfCuboids);
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T          = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();

  const T      physDeltaX = 2 * parameters.get<parameters::INLET_RADIUS>() / parameters.get<parameters::RESOLUTION>();
  STLreader<T> stlReader("../bifurcation3d.stl", physDeltaX);
  IndicatorLayer3D<T> extendedDomain(stlReader, physDeltaX);
  geometry.rename(0, 2, extendedDomain);
  geometry.rename(2, 1, stlReader);
  geometry.clean();

  // rename the material at the inlet
  IndicatorCircle3D<T>   inletCircle(parameters.get<parameters::INLET_CENTER>(),
                                     parameters.get<parameters::INLET_NORMAL>(),
                                     parameters.get<parameters::INLET_RADIUS>());
  IndicatorCylinder3D<T> inlet(inletCircle, 2 * physDeltaX);

  geometry.rename(2, 3, 1, inlet);

  // rename the material at the outlet0
  IndicatorCircle3D<T>   outletCircle0(parameters.get<parameters::OUTLET_CENTER0>(),
                                       parameters.get<parameters::OUTLET_NORMAL0>(),
                                       0.95 * parameters.get<parameters::OUTLET_RADIUS0>());
  IndicatorCylinder3D<T> outlet0(outletCircle0, 4 * physDeltaX);
  geometry.rename(2, 4, outlet0);

  // rename the material at the outlet1
  IndicatorCircle3D<T>   outletCircle1(parameters.get<parameters::OUTLET_CENTER1>(),
                                       parameters.get<parameters::OUTLET_NORMAL1>(),
                                       0.95 * parameters.get<parameters::OUTLET_RADIUS1>());
  IndicatorCylinder3D<T> outlet1(outletCircle1, 4 * physDeltaX);
  geometry.rename(2, 5, outlet1);

  IndicatorCircle3D<T>   inletCircleExtended(parameters.get<parameters::INLET_CENTER>(),
                                             parameters.get<parameters::INLET_NORMAL>(),
                                             parameters.get<parameters::INLET_RADIUS>() + 2 * physDeltaX);
  IndicatorCylinder3D<T> inletExtended(inletCircleExtended, 2 * physDeltaX);
  geometry.rename(2, 6, inletExtended);

  geometry.clean();
  geometry.innerClean(true);
  geometry.checkForErrors();

  clout << "Prepare Geometry ... OK" << std::endl;
  return;
}

void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T            = MyCase::value_t;
  auto& parameters   = myCase.getParameters();
  auto& geometry     = myCase.getGeometry();
  auto& NSlattice    = myCase.getLattice(NavierStokes {});
  auto& ADlattice    = myCase.getLattice(Concentration0 {});
  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ADDESCRIPTOR = MyCase::descriptor_t_of<Concentration0>;

  NSlattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, NSDESCRIPTOR>>(
      int {parameters.get<parameters::RESOLUTION>()}, // resolution: number of voxels per charPhysL
      (T)0.557646, // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
      (T)parameters.get<parameters::INLET_RADIUS>() * 2., // charPhysLength: reference length of simulation geometry
      (T)parameters.get<parameters::REYNOLDS>() * 1.5e-5 /
          (parameters.get<parameters::INLET_RADIUS>() *
           2),   // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
      (T)1.5e-5, // physViscosity: physical kinematic viscosity in __m^2 / s__
      (T)1.225   // physDensity: physical density in __kg / m^3__
  );

  const auto& converter = NSlattice.getUnitConverter();
  converter.print();
  ADlattice.setUnitConverter(converter);

  // Material=1 --> bulk dynamics
  // Material=3 --> bulk dynamics (inflow)
  auto inflowIndicator = geometry.getMaterialIndicator({1, 3});
  dynamics::set<BGKdynamics>(NSlattice, inflowIndicator);
  dynamics::set<ParticleAdvectionDiffusionBGKdynamics>(ADlattice, inflowIndicator);

  // Material=2 --> bounce-back / zero distribution dynamics
  dynamics::set<BounceBack>(NSlattice, geometry.getMaterialIndicator({2}));
  dynamics::set<ZeroDistributionDynamics>(ADlattice, geometry.getMaterialIndicator({2}));

  // Material=4,5 -->bulk dynamics / do-nothing (outflow)
  auto outflowIndicator = geometry.getMaterialIndicator({4, 5});
  dynamics::set<BGKdynamics>(NSlattice, outflowIndicator);

  // Material=6 --> bounce-back / bounce-back
  dynamics::set<BounceBack>(NSlattice, geometry, 6);
  dynamics::set<BounceBack>(ADlattice, geometry, 6);

  // Setting of the boundary conditions
  boundary::set<boundary::InterpolatedPressure>(NSlattice, geometry, 3);
  boundary::set<boundary::InterpolatedVelocity>(NSlattice, outflowIndicator);
  boundary::set<boundary::ZeroDistribution>(ADlattice, geometry, 2);
  boundary::set<boundary::AdvectionDiffusionDirichlet>(ADlattice, geometry, 3);
  setZeroGradientBoundary<T, ADDESCRIPTOR>(ADlattice, outflowIndicator);
  boundary::set<boundary::ExternalField<T, ADDESCRIPTOR, descriptors::VELOCITY, descriptors::VELOCITY2>>(
      ADlattice, geometry.getMaterialIndicator({2, 3, 4, 5, 6}));

  auto& coupling =
      myCase.setCouplingOperator("DragForce", AdvectionDiffusionParticleCoupling3D<ade_forces::AdvDiffDragForce3D> {},
                                 names::NavierStokes {}, NSlattice, names::AdvectionDiffusion {}, ADlattice);
  coupling.restrictTo(geometry.getMaterialIndicator(1));

  {
    auto& communicator = ADlattice.getCommunicator(stage::PostCoupling());
    communicator.requestField<descriptors::VELOCITY>();
    communicator.requestField<descriptors::VELOCITY2>();
    communicator.requestOverlap(ADlattice.getOverlap());
    communicator.exchangeRequests();
  }
  clout << "Prepare Lattice ... OK" << std::endl;
  return;
}

void setInitialValues(MyCase& myCase)
{
  OstreamManager clout(std::cout, "setInitialValues");
  clout << "Set Initial Values ..." << std::endl;
  using T            = MyCase::value_t;
  auto& parameters   = myCase.getParameters();
  auto& geometry     = myCase.getGeometry();
  auto& NSlattice    = myCase.getLattice(NavierStokes {});
  auto& ADlattice    = myCase.getLattice(Concentration0 {});
  auto& converter   = NSlattice.getUnitConverter();
  using ADDESCRIPTOR = MyCase::descriptor_t_of<Concentration0>;

  // Initial conditions
  AnalyticalConst3D<T, T> rho0(1.e-8);
  std::vector<T>          velocity(3, T());
  AnalyticalConst3D<T, T> u0(velocity);

  // Initialize all values of distribution functions to their local equilibrium
  ADlattice.iniEquilibrium(geometry.getMaterialIndicator({1, 2, 4, 5, 6}), rho0, u0);
  const T omega = converter.getLatticeRelaxationFrequency();
  const T omegaAD =
      converter.getLatticeRelaxationFrequencyFromDiffusivity<ADDESCRIPTOR>(parameters.get<parameters::DIFFUSION>());
  NSlattice.setParameter<descriptors::OMEGA>(omega);
  ADlattice.setParameter<descriptors::OMEGA>(omegaAD);

  // Lattice initialize
  NSlattice.initialize();
  ADlattice.initialize();
  clout << "Set Initial Values ... OK" << std::endl;
}


void setBoundaryValues(MyCase& myCase, std::size_t iT)
{

  OstreamManager clout(std::cout, "setBoundaryValues");

  using T                = MyCase::value_t;
  auto&       parameters = myCase.getParameters();
  auto&       geometry   = myCase.getGeometry();
  auto&       NSlattice  = myCase.getLattice(NavierStokes {});
  const auto& converter  = NSlattice.getUnitConverter();

  // No of time steps for smooth start-up
  std::size_t iTmaxStart = converter.getLatticeTime(0.8 * parameters.get<parameters::PHYS_MAX_T>());
  // Set inflow velocity
  T maxVelocity = converter.getCharPhysVelocity() * 3. / 4. *
                  util::pow(parameters.get<parameters::INLET_RADIUS>(), 2) /
                  util::pow(parameters.get<parameters::OUTLET_RADIUS0>(), 2);
  if (iT % converter.getLatticeTime(parameters.get<parameters::PHYS_STAT_ITER_T>()) == 0) {
    if (iT <= iTmaxStart) {
      SinusStartScale<T, std::size_t> startScale(iTmaxStart, T(1));
      std::size_t                     iTvec[1] = {iT};
      T                               frac[1]  = {T(0)};
      startScale(frac, iTvec);
      maxVelocity *= frac[0];
    }

    Vector<T, 3>          outletCenter0 = parameters.get<parameters::OUTLET_CENTER0>();
    Vector<T, 3>          outletNormal0 = parameters.get<parameters::OUTLET_NORMAL0>();
    Vector<T, 3>          outletCenter1 = parameters.get<parameters::OUTLET_CENTER1>();
    Vector<T, 3>          outletNormal1 = parameters.get<parameters::OUTLET_NORMAL1>();
    CirclePoiseuille3D<T> poiseuilleU4(outletCenter0[0], outletCenter0[1], outletCenter0[2], outletNormal0[0],
                                       outletNormal0[1], outletNormal0[2],
                                       parameters.get<parameters::OUTLET_RADIUS0>() * 0.95, -maxVelocity);
    CirclePoiseuille3D<T> poiseuilleU5(outletCenter1[0], outletCenter1[1], outletCenter1[2], outletNormal1[0],
                                       outletNormal1[1], outletNormal1[2],
                                       parameters.get<parameters::OUTLET_RADIUS1>() * 0.95, -maxVelocity);
    momenta::setVelocity(NSlattice, geometry.getMaterialIndicator({4}), poiseuilleU4);
    momenta::setVelocity(NSlattice, geometry.getMaterialIndicator({5}), poiseuilleU5);

    NSlattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
        ProcessingContext::Simulation);
  }
}

void getResults(MyCase& myCase, std::size_t iT, util::Timer<MyCase::value_t>& timer)
{
  using T                = MyCase::value_t;
  auto&       parameters = myCase.getParameters();
  auto&       geometry   = myCase.getGeometry();
  auto&       NSlattice  = myCase.getLattice(NavierStokes {});
  auto&       ADlattice  = myCase.getLattice(Concentration0 {});
  const auto& converter  = NSlattice.getUnitConverter();
  using NSDESCRIPTOR     = MyCase::descriptor_t_of<NavierStokes>;
  using ADDESCRIPTOR     = MyCase::descriptor_t_of<Concentration0>;

  OstreamManager                              clout(std::cout, "getResults");
  SuperVTMwriter3D<T>                         vtmWriter("bifurcation3d_fluid");
  SuperVTMwriter3D<T>                         vtmWriterAD("bifurcation3d_particle");
  SuperLatticePhysVelocity3D<T, NSDESCRIPTOR> velocity(NSlattice, converter);
  SuperLatticeVelocity3D<T, NSDESCRIPTOR>     latticeVelocity(NSlattice);

  SuperLatticePhysPressure3D<T, NSDESCRIPTOR>                     pressure(NSlattice, converter);
  SuperLatticeDensity3D<T, ADDESCRIPTOR>                          particles(ADlattice);
  SuperLatticePhysField3D<T, ADDESCRIPTOR, descriptors::VELOCITY> extField(ADlattice,
                                                                           converter.getConversionFactorVelocity());

  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(pressure);
  vtmWriterAD.addFunctor(particles);
  vtmWriterAD.addFunctor(extField);

  if (iT == 0) {
    SuperLatticeCuboid3D<T, NSDESCRIPTOR> cuboid(NSlattice);
    SuperLatticeRank3D<T, NSDESCRIPTOR>   rank(NSlattice);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();
    vtmWriterAD.createMasterFile();

    // Print some output of the chosen simulation setup
    clout << "N=" << parameters.get<parameters::RESOLUTION>()
          << "; maxTimeSteps=" << parameters.get<parameters::PHYS_MAX_T>()
          << "; noOfCuboid=" << geometry.getCuboidDecomposition().size()
          << "; Re=" << parameters.get<parameters::REYNOLDS>() << "; St="
          << (2. * parameters.get<parameters::PART_RHO>() * parameters.get<parameters::PART_RADIUS>() *
              parameters.get<parameters::PART_RADIUS>() * converter.getCharPhysVelocity()) /
                 (9. * converter.getPhysViscosity() * converter.getPhysDensity() * converter.getCharPhysLength())
          << std::endl;
  }

  if (iT % converter.getLatticeTime(parameters.get<parameters::PHYS_STAT_ITER_T>()) == 0) {
    NSlattice.setProcessingContext(ProcessingContext::Evaluation);
    ADlattice.setProcessingContext(ProcessingContext::Evaluation);

    // Writes the vtk files
    vtmWriter.write(iT);
    vtmWriterAD.write(iT);

    // GIF Writer
    SuperEuklidNorm3D<T>   normVel(velocity);
    HyperplaneLattice3D<T> gifLattice(
        geometry.getCuboidDecomposition(),
        Hyperplane3D<T>().centeredIn(geometry.getCuboidDecomposition().getMotherCuboid()).normalTo({0, -1, 0}), 600);
    BlockReduction3D2D<T> planeReductionVelocity(normVel, gifLattice, BlockDataSyncMode::ReduceOnly);
    BlockReduction3D2D<T> planeReductionParticles(particles, gifLattice, BlockDataSyncMode::ReduceOnly);
    // write output as JPEG
    heatmap::write(planeReductionVelocity, iT);
    heatmap::write(planeReductionParticles, iT);

    // Writes output on the console
    timer.update(iT);
    timer.printStep();
    NSlattice.getStatistics().print(iT, iT * converter.getCharLatticeVelocity() / T(converter.getResolution()));

    // preparation for flux computations
    const std::vector<int> materials = {1, 3, 4, 5};
    IndicatorCircle3D<T>   inlet(parameters.get<parameters::INLET_CENTER>() +
                                     2. * converter.getPhysDeltaX() * parameters.get<parameters::INLET_NORMAL>(),
                                 parameters.get<parameters::INLET_NORMAL>(),
                                 parameters.get<parameters::INLET_RADIUS>() + 2. * converter.getPhysDeltaX());
    IndicatorCircle3D<T>   outlet0(parameters.get<parameters::OUTLET_CENTER0>() +
                                       2. * converter.getPhysDeltaX() * parameters.get<parameters::OUTLET_NORMAL0>(),
                                   parameters.get<parameters::OUTLET_NORMAL0>(),
                                   parameters.get<parameters::OUTLET_RADIUS0>() + 2. * converter.getPhysDeltaX());
    IndicatorCircle3D<T>   outlet1(parameters.get<parameters::OUTLET_CENTER1>() +
                                       2. * converter.getPhysDeltaX() * parameters.get<parameters::OUTLET_NORMAL1>(),
                                   parameters.get<parameters::OUTLET_NORMAL1>(),
                                   parameters.get<parameters::OUTLET_RADIUS1>() + 2. * converter.getPhysDeltaX());

    // Flux of the fluid at the inlet and outlet regions
    SuperPlaneIntegralFluxVelocity3D<T> vFluxInflow(NSlattice, converter, geometry, inlet, materials);
    vFluxInflow.print("inflow", "ml/s");
    SuperPlaneIntegralFluxPressure3D<T> pFluxInflow(NSlattice, converter, geometry, inlet, materials);
    pFluxInflow.print("inflow", "N", "Pa");
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOutflow0(NSlattice, converter, geometry, outlet0, materials);
    vFluxOutflow0.print("outflow0", "ml/s");
    SuperPlaneIntegralFluxPressure3D<T> pFluxOutflow0(NSlattice, converter, geometry, outlet0, materials);
    pFluxOutflow0.print("outflow0", "N", "Pa");
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOutflow1(NSlattice, converter, geometry, outlet1, materials);
    vFluxOutflow1.print("outflow1", "ml/s");
    SuperPlaneIntegralFluxPressure3D<T> pFluxOutflow1(NSlattice, converter, geometry, outlet1, materials);
    pFluxOutflow1.print("outflow1", "N", "Pa");

    int input[4] = {0};
    T   mFlux[5] = {0.}, mFlux0[5] = {0.}, mFlux1[5] = {0.};
    // Flux of particles at the inlet and outlet regions: Inflow, Outflow0 and Outlfow1
    SuperPlaneIntegralFluxMass3D<T> mFluxInflow(latticeVelocity, particles, geometry,
                                                converter.getConversionFactorMass(),
                                                converter.getConversionFactorTime(), inlet, materials);
    SuperPlaneIntegralFluxMass3D<T> mFluxOutflow0(latticeVelocity, particles, geometry,
                                                  converter.getConversionFactorMass(),
                                                  converter.getConversionFactorTime(), outlet0, materials);
    SuperPlaneIntegralFluxMass3D<T> mFluxOutflow1(latticeVelocity, particles, geometry,
                                                  converter.getConversionFactorMass(),
                                                  converter.getConversionFactorTime(), outlet1, materials);

    mFluxInflow(mFlux, input);
    mFluxOutflow0(mFlux0, input);
    mFluxOutflow1(mFlux1, input);

    // Since more diffusion is added to ensure stability the computed escaperate falls short of the real value,
    // therefore it is scaled by the factor 1.4 computed by a simulation without drag force. This value is computed
    // for this specific setup. For further information see R.Trunk, T.Henn, W.Dörfler, H.Nirschl, M.J.Krause,
    // "Inertial Dilute Particulate Fluid Flow Simulations with an Euler-Euler Lattice Boltzmann Method"
    T escr = -(mFlux0[0] + mFlux1[0]) / mFlux[0] * 1.4;
    clout << "escapeRate=" << escr << "; captureRate=" << 1 - escr << std::endl;
  }
}

void simulate(MyCase& myCase)
{
  using T                = MyCase::value_t;
  auto&       parameters = myCase.getParameters();
  auto&       geometry   = myCase.getGeometry();
  auto&       NSlattice  = myCase.getLattice(NavierStokes {});
  auto&       ADlattice  = myCase.getLattice(Concentration0 {});
  const auto& converter  = NSlattice.getUnitConverter();
  using NSDESCRIPTOR     = MyCase::descriptor_t_of<NavierStokes>;
  using ADDESCRIPTOR     = MyCase::descriptor_t_of<Concentration0>;
  auto& coupling         = dynamic_cast<
              SuperLatticeCoupling<AdvectionDiffusionParticleCoupling3D<ade_forces::AdvDiffDragForce3D>,
                                   meta::map<names::NavierStokes, descriptors::VALUED_DESCRIPTOR<T, NSDESCRIPTOR>,
                                             names::AdvectionDiffusion, descriptors::VALUED_DESCRIPTOR<T, ADDESCRIPTOR>>>&>
                                             ( myCase.getOperator("DragForce"));
  // Compute the drag force parameters
  ade_forces::AdvDiffDragForce3D::computeParametersFromRhoAndRadius<T>(
      parameters.get<parameters::PART_RHO>(), parameters.get<parameters::PART_RADIUS>(), coupling, converter);

  // ===Step 9: Main Loop with Timer ===
  util::Timer<T> timer(converter.getLatticeTime(parameters.get<parameters::PHYS_MAX_T>()),
                       geometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT <= converter.getLatticeTime(parameters.get<parameters::PHYS_MAX_T>()); ++iT) {
    ADlattice.setParameter<descriptors::LATTICE_TIME>(iT);
    coupling.setParameter<AdvectionDiffusionParticleCoupling3D<ade_forces::AdvDiffDragForce3D>::LATTICE_TIME>(iT);
    getResults(myCase, iT, timer);
    setBoundaryValues(myCase, iT);
    coupling.apply();
    ADlattice.getCommunicator(stage::PostCoupling()).communicate();
    NSlattice.collideAndStream();
    ADlattice.collideAndStream();
  }

  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[])
{

  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<REYNOLDS>(50);
    myCaseParameters.set<RESOLUTION>(19);
    myCaseParameters.set<PHYS_STAT_ITER_T>(0.1419);   // time interval in s after which boundary conditions are updated
    myCaseParameters.set<DIFFUSION>(1e-6);     // diffusion coefficient for advection-diffusion equation
    myCaseParameters.set<PART_RADIUS>(1.5e-4); // particles radius
    myCaseParameters.set<PART_RHO>(998.2);     // particles density
    myCaseParameters.set<PHYS_MAX_T>(10);      // max. simulation time in s, SI unit
    myCaseParameters.set<INLET_CENTER>({0., 0., 0.0786395});
    myCaseParameters.set<OUTLET_CENTER0>({-0.0235929682287551, -0.000052820468762797, -0.021445708949909});
    myCaseParameters.set<OUTLET_CENTER1>({0.0233643529416147, 0.00000212439067050152, -0.0211994104877918});
    myCaseParameters.set<INLET_RADIUS>(0.00999839);
    myCaseParameters.set<OUTLET_RADIUS0>(0.007927);
    myCaseParameters.set<OUTLET_RADIUS1>(0.00787134);
    myCaseParameters.set<INLET_NORMAL>({0., 0., -1.});
    myCaseParameters.set<OUTLET_NORMAL0>({0.505126, -0.04177, 0.862034});
    myCaseParameters.set<OUTLET_NORMAL1>({-0.483331, -0.0102764, 0.875377});
  }

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase);

  /// === Step 7: Set Initial Conditions ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);
}
