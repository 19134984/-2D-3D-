/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2022 Florian Raichle
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

/** @file micromixer3d
 * Simulating adsorption in a static mixing reactor using an Euler-Euler approach.
 * The model is based on the linear driving force model and uses advection
 * diffusion reaction lattices for particles, solute and particle loading.
 *
 * Different isotherms and mass transfer models can be used.
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;
//using namespace olb::graphics;

/// === Step 1: Declare Simulation Structure ===
using MyCase = Case<
  NavierStokes,       Lattice<double, descriptors::D3Q19<>>,
  Concentration0,     Lattice<double, descriptors::D3Q7<descriptors::VELOCITY>>,
  Concentration1,     Lattice<double, descriptors::D3Q7<descriptors::VELOCITY>>,
  Concentration2,     Lattice<double, descriptors::D3Q7<descriptors::VELOCITY>>
>;

namespace olb::parameters {

struct PARTICLE_RADIUS : public descriptors::FIELD_BASE<1> { };
struct PARTICLE_CONCENTRATION : public descriptors::FIELD_BASE<1> { };
struct PARTICLE_DENSITY : public descriptors::FIELD_BASE<1> { };

struct ISO_CONST_A : public descriptors::FIELD_BASE<1> { };
struct ISO_CONST_B : public descriptors::FIELD_BASE<1> { };
struct K_F : public descriptors::FIELD_BASE<1> { };
struct C_0 : public descriptors::FIELD_BASE<1> { };
struct D_S : public descriptors::FIELD_BASE<1> { };

}

/// === Step 3: Create Mesh ===
Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  const T dx = params.get<parameters::DOMAIN_L>() / params.get<parameters::RESOLUTION>();
  const int noOfCuboids = util::max(16, 4 * singleton::mpi().getSize());

  STLreader<T> stlReader("microMixer3d_small.stl",
                         dx, params.get<parameters::STL_SIZE>());
  IndicatorLayer3D<T> extendedDomain(stlReader, dx);

  CuboidDecomposition3D<T> cuboidDecomposition(extendedDomain, dx, noOfCuboids);

  Mesh<T,MyCase::d> mesh(extendedDomain, dx, noOfCuboids);
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  return mesh;
}

/// === Step 5: Prepare Geometry ===
void prepareGeometry(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;

  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  const T dx = params.get<parameters::DOMAIN_L>() / params.get<parameters::RESOLUTION>();
  STLreader<T> stlReader("microMixer3d_small.stl",
                         dx, params.get<parameters::STL_SIZE>());
  geometry.rename(0, 2, stlReader);
  geometry.rename(2, 1, 1);

// Returns the minimum phys position in each direction for material 2
  Vector<T, 3> minR = geometry.getStatistics().getMinPhysR(2);
  Vector<T, 3> maxR = geometry.getStatistics().getMaxPhysR(2);
  Vector<T, 3> centerR = geometry.getStatistics().getCenterPhysR(2);
  Vector<T, 3> extend = geometry.getStatistics().getPhysExtend(2);
  extend[2] = dx;

  // sets circle of both inflows and the outflow, with direction and radius
  IndicatorCircle3D<T> inflow1(minR[0], minR[1], centerR[2], 0., -1., 0.,
                               (maxR[0] - minR[0]) / T(3));
  IndicatorCircle3D<T> inflow2(maxR[0], minR[1], centerR[2], 0., -1., 0.,
                               (maxR[0] - minR[0]) / T(3));
  IndicatorCircle3D<T> outflow(minR[0], maxR[1], centerR[2],
                               0., 1., 0., (maxR[0] - minR[0]) / T(3));

  // sets cylinder on that in-/out-flow circles with length
  IndicatorCylinder3D<T> layerInflow1(inflow1, dx);
  IndicatorCylinder3D<T> layerInflow2(inflow2, dx);
  IndicatorCylinder3D<T> layerOutflow(outflow, dx);
  // renames all boundary voxels of material fromBcMat to toBcMat if two neighbour voxel
  // in the direction of the discrete normal are fluid voxel with material fluidM in the region
  // where the indicator function is fulfilled
  geometry.rename(2, 3, 1, layerInflow1); // layer of inflow1 gets mat = 3
  geometry.rename(2, 4, 1, layerInflow2); // layer of inflow2 gets mat = 4
  geometry.rename(2, 5, 1, layerOutflow); // layer of outflow gets mat = 5

  geometry.clean();
  geometry.innerClean();
  geometry.checkForErrors();

  geometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

template<typename T, typename ADDESCRIPTOR>
void prepareLatticeAD(
  MyCase& myCase,
  Lattice<T, ADDESCRIPTOR>& lattice,
  int inlet, int noInlet)
{
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  const T viscosity = params.get<parameters::PHYS_CHAR_VELOCITY>() * params.get<parameters::DOMAIN_L>() / params.get<parameters::REYNOLDS>();
  lattice.template setUnitConverter<AdsorptionConverter<T, ADDESCRIPTOR>>(
      (T)   params.get<parameters::DOMAIN_L>()/params.get<parameters::RESOLUTION>(),  //physDeltaX
      (T)   myCase.getLattice(NavierStokes{}).getUnitConverter().getPhysDeltaT(),         //physDeltaT
      (T)   params.get<parameters::DOMAIN_L>(),     //charPhysLength
      (T)   params.get<parameters::PHYS_CHAR_VELOCITY>(),    //charPhysVelocity
      (T)   viscosity / params.get<parameters::SCHMIDT>(),
      (T)   1.,
      (T)   1.,
      (T)   viscosity
  );
  const auto& converter = lattice.getUnitConverter();

  // dynamics for ADE
  dynamics::set<SourcedAdvectionDiffusionBGKdynamics>(lattice, geometry.getMaterialIndicator({1, 5}));
  boundary::set<boundary::BounceBack>(lattice, geometry, 2);
  boundary::set<boundary::BounceBack>(lattice, geometry, noInlet);

  // boundary for ADE
  dynamics::set<SourcedAdvectionDiffusionBGKdynamics>(lattice, geometry, inlet);
  boundary::set<boundary::AdvectionDiffusionDirichlet>(lattice, geometry, inlet);
  setZeroGradientBoundary<T,ADDESCRIPTOR>(lattice, geometry, 5);

  const T omegaAD = converter.getLatticeRelaxationFrequency();
  lattice.template setParameter<descriptors::OMEGA>(omegaAD);
  lattice.initialize();
}

void prepareLattice(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& NSlattice = myCase.getLattice(NavierStokes{});
  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

  const T viscosity = params.get<parameters::PHYS_CHAR_VELOCITY>() * params.get<parameters::DOMAIN_L>() / params.get<parameters::REYNOLDS>();
  NSlattice.template setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, NSDESCRIPTOR>> (
      (T) params.get<parameters::RESOLUTION>(),  // resolution: number of voxels per charPhysL
      (T) params.get<parameters::LATTICE_RELAXATION_TIME>(),         // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
      (T) params.get<parameters::DOMAIN_L>(),       // charPhysLength: reference channelLength of simulation geometry
      (T) params.get<parameters::PHYS_CHAR_VELOCITY>(),      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
      (T) viscosity,                             // physViscosity: physical kinematic viscosity in __m^2 / s__
      (T) 1000                                   // physDensity: physical density in __kg / m^3__
  );
  const auto& converter = NSlattice.getUnitConverter();
  converter.print();

  // dynamics for fluid
  dynamics::set<BGKdynamics>(NSlattice, geometry.getMaterialIndicator({1, 3, 4, 5}));

  // boundary conditions for fluid
  // wall
  boundary::set<boundary::BounceBack>(NSlattice, geometry, 2);
  // inlet
  boundary::set<boundary::InterpolatedPressure>(NSlattice, geometry, 3);
  boundary::set<boundary::InterpolatedPressure>(NSlattice, geometry, 4);
  // outlet
  boundary::set<boundary::InterpolatedVelocity>(NSlattice, geometry, 5);

  NSlattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  // converter, dynamics, boundary conditions for ADE lattices
  prepareLatticeAD(myCase, myCase.getLattice(Concentration0{}), 3, 4);
  prepareLatticeAD(myCase, myCase.getLattice(Concentration1{}), 4, 3);
  prepareLatticeAD(myCase, myCase.getLattice(Concentration2{}), 3, 4);
  myCase.getLattice(Concentration0{}).getUnitConverter().print();

  {
    auto &communicator = NSlattice.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(NSlattice.getOverlap());
    communicator.exchangeRequests();
  }

  auto& coupling = myCase.setCouplingOperator(
    "ParticleTransport",
    AdsorptionFullCoupling3D<AdsorptionReaction<Isotherm::LangmuirIsotherm>,
                             ade_forces::AdvDiffDragForce3D>{},
    names::NavierStokes{},   NSlattice,
    names::Concentration0{}, myCase.getLattice(Concentration0{}),
    names::Concentration1{}, myCase.getLattice(Concentration1{}),
    names::Concentration2{}, myCase.getLattice(Concentration2{}));

  // Setting Isotherm Parameters
  Isotherm::LangmuirIsotherm::setParameters<T>(
    params.get<parameters::ISO_CONST_A>(),
    params.get<parameters::ISO_CONST_B>(),
    coupling);

  // Setting Adsorption Reaction Parameters
  coupling.template setParameter<AdsorptionReaction<Isotherm::LangmuirIsotherm>::K_F>(params.get<parameters::K_F>());
  coupling.template setParameter<AdsorptionReaction<Isotherm::LangmuirIsotherm>::D_S>(params.get<parameters::D_S>());
  coupling.template setParameter<AdsorptionReaction<Isotherm::LangmuirIsotherm>::C_0>(params.get<parameters::C_0>());
  coupling.template setParameter<AdsorptionReaction<Isotherm::LangmuirIsotherm>::R_P>(params.get<parameters::PARTICLE_RADIUS>());

  // Compute the interaction parameters
  AdsorptionReaction<Isotherm::LangmuirIsotherm>::computeParameters<T>(
    coupling,
    myCase.getLattice(Concentration0{}).getUnitConverter());

  AdsorptionReaction<Isotherm::LangmuirIsotherm>::print<T>(
    clout,
    coupling,
    myCase.getLattice(Concentration0{}).getUnitConverter());

  // Compute the drag force parameters
  ade_forces::AdvDiffDragForce3D::computeParametersFromRhoAndRadius<T>(
    params.get<parameters::PARTICLE_DENSITY>(),
    params.get<parameters::PARTICLE_RADIUS>(),
    coupling,
    converter);

  clout << "Prepare Lattice ... OK" << std::endl;
}

template<typename T, typename ADDESCRIPTOR>
void setInitialValuesAD(
  MyCase& myCase,
  Lattice<T, ADDESCRIPTOR>& lattice,
  int inlet, int noInlet, T rhoInlet)
{
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  AnalyticalConst3D<T, T> rhoI(rhoInlet);
  AnalyticalConst3D<T, T> u0(0., 0., 0.);
  AnalyticalConst3D<T, T> rhoSmall(0);

  momenta::setDensity(lattice, geometry.getMaterialIndicator({1, 2, 5}), rhoSmall);
  momenta::setDensity(lattice, geometry.getMaterialIndicator({inlet}), rhoI);
  momenta::setDensity(lattice, geometry.getMaterialIndicator({noInlet}), rhoSmall);

  lattice.iniEquilibrium(geometry.getMaterialIndicator({1, 2, 5}), rhoSmall, u0);
  lattice.iniEquilibrium(geometry, inlet, rhoI, u0);
  lattice.iniEquilibrium(geometry, noInlet, rhoSmall, u0);

  lattice.initialize();
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& NSlattice = myCase.getLattice(NavierStokes{});

  // initialisation for fluid
  AnalyticalConst3D<T, T> rho1(1.);
  AnalyticalConst3D<T, T> u0(0., 0., 0.);

  NSlattice.initialize();

  // initial values for ADE lattices
  setInitialValuesAD(myCase, myCase.getLattice(Concentration0{}), 3, 4, 1.0);
  setInitialValuesAD(myCase, myCase.getLattice(Concentration1{}), 4, 3, 1.0);
  setInitialValuesAD(myCase, myCase.getLattice(Concentration2{}), 3, 4, 0.0);
}

/// === Step 7.2.1: Update the Boundary Values and Fields at Times ===
void setTemporalValues(MyCase& myCase,
                       std::size_t iT) {

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& lattice = myCase.getLattice(NavierStokes{});
  const auto& converter = lattice.getUnitConverter();

  std::vector<T> maxVelocity(3, T());
  const T distanceToBoundary = converter.getPhysDeltaX() / T(2);
  const T physVelNS = converter.getCharPhysVelocity();
  const size_t itStartTime = converter.getLatticeTime(params.get<parameters::PHYS_START_T>());

  if (iT <= itStartTime && iT % 50 == 0) {

    SinusStartScale<T, size_t> startScale(itStartTime, T(1)); //1+-amplitude
    size_t help[1] = { iT };
    T frac[3] = { T() };
    startScale(frac, help);

    // set velocity on boundary
    maxVelocity[1] = physVelNS * frac[0];

    RectanglePoiseuille3D<T> u5(geometry, 5, maxVelocity,
                                distanceToBoundary, distanceToBoundary, distanceToBoundary);
    momenta::setVelocity(lattice, geometry.getMaterialIndicator({5}), u5);
    lattice.setProcessingContext<olb::Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

void getResults(
  MyCase& myCase,
  util::Timer<double>& timer,
  size_t iT)
{
  OstreamManager clout(std::cout, "getResults");
  using T = MyCase::value_t;
  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ADDESCRIPTOR = MyCase::descriptor_t_of<Concentration0>;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();
  auto& NSlattice = myCase.getLattice(NavierStokes{});
  auto& conc0lattice = myCase.getLattice(Concentration0{});
  auto& conc1lattice = myCase.getLattice(Concentration1{});
  auto& conc2lattice = myCase.getLattice(Concentration2{});
  const auto& converter = NSlattice.getUnitConverter();

  const int itTotalTime = NSlattice.getUnitConverter().getLatticeTime(params.get<parameters::PHYS_MAX_T>());
  const int iTperiodConsole = itTotalTime / 10;   // Console output
  const int iTperiodVTK = itTotalTime / 10;   // Writes the vtk files

  SuperVTMwriter3D<T> vtmWriter( "microMixer3d" );

  if (iT == 0) {
    NSlattice.setProcessingContext(ProcessingContext::Evaluation);
    SuperLatticeCuboid3D<T, NSDESCRIPTOR> cuboid(NSlattice);
    SuperLatticeRank3D<T, NSDESCRIPTOR> rank(NSlattice);

    vtmWriter.write(cuboid);
    vtmWriter.write(rank);

    // have to be called before calling write(int iT=0), since it creates
    // the master pvd file, where all vti are linked!
    vtmWriter.createMasterFile();
  }


  // console output
  if (iT % iTperiodConsole == 0) {
    timer.update(iT);
    timer.printStep();

    // output for latticeStatistics
    NSlattice.getStatistics().print(iT, converter.getPhysTime(iT));
    conc0lattice.getStatistics().print(iT, converter.getPhysTime(iT));
  }

  // vtk output
  if (iT % iTperiodVTK == 0) {
    NSlattice.setProcessingContext(ProcessingContext::Evaluation);
    conc0lattice.setProcessingContext(ProcessingContext::Evaluation);
    conc1lattice.setProcessingContext(ProcessingContext::Evaluation);
    conc2lattice.setProcessingContext(ProcessingContext::Evaluation);
    SuperGeometryF<T,3> materials(geometry);
    SuperLatticePhysVelocity3D<T, NSDESCRIPTOR> velocityNS(NSlattice, converter);
    SuperLatticePhysPressure3D<T, NSDESCRIPTOR> pressure(NSlattice, converter);
    SuperLatticeDensity3D<T, ADDESCRIPTOR> adsorptive(conc0lattice);
    SuperLatticeDensity3D<T, ADDESCRIPTOR> loading(conc2lattice);
    SuperLatticeDensity3D<T, ADDESCRIPTOR> soluteConcentration(conc1lattice);
    vtmWriter.addFunctor(velocityNS);
    vtmWriter.addFunctor(pressure);
    vtmWriter.addFunctor(adsorptive, "particle concentration");
    vtmWriter.addFunctor(loading, "loading");
    vtmWriter.addFunctor(soluteConcentration, "solute concentration");
    vtmWriter.addFunctor(materials);
    vtmWriter.write(iT);
  }
}

/// === Step 7: Simulate ===
void simulate(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();

  /// === Step 7.1: Definition of Initial and Boundary Values and Fields ===
  setInitialValues(myCase);

  /// === Step 7.2: Main Loop with Timer ===
  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(
    parameters.get<parameters::PHYS_MAX_T>());

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    /// === Step 7.2.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    myCase.getOperator("ParticleTransport").apply();

    /// === Step 7.2.2: Computation and Output of the Results ===
    getResults(myCase, timer, iT);

    /// === Step 7.2.3: Collide and Stream Execution ===
    myCase.getLattice(Concentration0{}).collideAndStream();
    myCase.getLattice(Concentration1{}).collideAndStream();
    myCase.getLattice(Concentration2{}).collideAndStream();
    myCase.getLattice(NavierStokes{}).collideAndStream();
  }

  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[]) {
  /// === Step 2: Initialization ===
  initialize(&argc, &argv);

  /// === Step 2.1: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<PHYS_CHAR_VELOCITY     >(0.2);
    myCaseParameters.set<DOMAIN_L               >(0.0133);  // hydraulic diameter = 4 * surface / perimeter
    myCaseParameters.set<RESOLUTION             >(20);       // resolution of the hydraulic diameter // 20
    myCaseParameters.set<REYNOLDS               >(50);
    myCaseParameters.set<SCHMIDT                >(100);
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.6125);
    myCaseParameters.set<PARTICLE_RADIUS        >(5.e-5);
    myCaseParameters.set<PARTICLE_CONCENTRATION >(5);
    myCaseParameters.set<PARTICLE_DENSITY       >(1700);

    myCaseParameters.set<ISO_CONST_A>(45);
    myCaseParameters.set<ISO_CONST_B>(0.5);
    myCaseParameters.set<K_F        >(0);
    myCaseParameters.set<C_0        >(1);
    myCaseParameters.set<D_S        >(5.e-11);

    myCaseParameters.set<PHYS_MAX_T   >(5);    // time for fluid simulation
    myCaseParameters.set<PHYS_START_T >(0.6);  // time to start fluid pulsation
    myCaseParameters.set<STL_SIZE     >(10);
  }

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase);

  /// === Step 7: Simulate ===
  simulate(myCase);
}
