/* Lattice Boltzmann sample, written in C++, using the OpenLB
 * library

 * Copyright (C) 2019-2022 Mathias J. Krause, Julius Jeßberger,
 * E-mail contact: info@openlb.net
 * The most recent release of OpenLB can be downloaded at
 * <http:  //www.openlb.net/>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

/** @file microMixerOpti3d.cpp
 * In this example, the operation of a small mixer is optimized: pulsation
 * at the inlet is chosen, s.t. the segregation intensity at the outlet is
 * minimized.
 * Control variables: period length, amplitude and acuteness of the inflow.
 * Side condition: LB-Navier-Stokes equation for fluid, LB-Advection-Diffusion
 * equation for solute concentration, with one-way coupling.
 * Objective: Danckwerts' segregation intensity at the outlet, averaged over
 * one time period.
 * Cf. https://doi.org/10.3390/fluids7050144 for the full study.
 * Note: In this example, the particle diffusivity has been increased and the
 * resolution has been reduced, in order to allow for quick testing.
 * For the physically reasonable diffusivity D=1.e-9, the resolution N=96 yields
 * good results. The optimization runs stable already at N=36.
 */

#undef PLATFORM_GPU_CUDA

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;
using namespace olb::opti;

using MyCase = Case<
  NavierStokes,   Lattice<double,D3Q19<>>,
  Concentration0, Lattice<double,D3Q7<VELOCITY>>>;

using MyOptiCase = OptiCaseFDQ<Controlled, MyCase>;

namespace olb::parameters {

// control variables
struct PHYS_PERIOD  : public descriptors::FIELD_BASE<1> { };
struct AMPLITUDE_PHYS_PRESSURE : public descriptors::FIELD_BASE<1> { };
struct DIFFERENCE_PERIOD : public descriptors::FIELD_BASE<1> { };
}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  const int nC = util::max(16, 4 * singleton::mpi().getSize());
  const T dx = params.get<parameters::DX>();
  STLreader<T> stlReader("microMixer3d.stl", dx, T{1});
  IndicatorLayer3D<T> extendedDomain(stlReader, dx);
  Mesh<T,MyCase::d> mesh(extendedDomain, dx, nC);
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();
  const T dx = params.get<parameters::DX>();

  STLreader<T> stlReader("microMixer3d.stl", dx, T{1.});
  geometry.rename(0, 2, stlReader);
  geometry.rename(2, 1, {1, 1, 1});

  // Returns the minimum phys position in each direction for material 2
  Vector<T,3> minR = geometry.getStatistics().getMinPhysR(2);
  Vector<T,3> maxR = geometry.getStatistics().getMaxPhysR(2);
  Vector<T,3> centerR = geometry.getStatistics().getCenterPhysR(2);

  // sets circle of both inflows and the outflow, with direction and radius
  IndicatorCircle3D<T> inflow1(minR[0], minR[1], centerR[2], 0., -1., 0.,
                               (maxR[0] - minR[0]) / 3.);
  IndicatorCircle3D<T> inflow2(maxR[0], minR[1], centerR[2], 0., -1., 0.,
                               (maxR[0] - minR[0]) / 3.);
  IndicatorCircle3D<T> outflow((maxR[0] + minR[0]) / 2., maxR[1], centerR[2],
                               0., 1., 0., (maxR[0] - minR[0]) / 3.);

  // sets cylinder on that in-/out-flow circles with length
  IndicatorCylinder3D<T> layerInflow1(inflow1, dx);
  IndicatorCylinder3D<T> layerInflow2(inflow2, dx);
  IndicatorCylinder3D<T> layerOutflow(outflow, dx);

  geometry.rename(2, 3, 1, layerInflow1); // layer of inflow1 gets mat = 3
  geometry.rename(2, 4, 1, layerInflow2); // layer of inflow2 gets mat = 4
  geometry.rename(2, 5, 1, layerOutflow); // layer of outflow gets mat = 5

  geometry.clean(false);
  geometry.innerClean(false);
  geometry.checkForErrors(false);
}

void prepareLattice(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& params = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& sLatticeAD = myCase.getLattice(Concentration0{});

  sLattice.setUnitConverter<UnitConverterFromResolutionAndLatticeVelocity<T,D3Q19<>>>(
    (int)   params.get<parameters::RESOLUTION>(),        //resolution
    ( T )   params.get<parameters::LATTICE_CHAR_VELOCITY>(),         //charLatticeVelocity
    ( T )   params.get<parameters::PHYS_CHAR_LENGTH>(),  //charPhysLength
    ( T )   params.get<parameters::PHYS_CHAR_VELOCITY>(),  //charPhysVelocity
    ( T )   params.get<parameters::PHYS_CHAR_VISCOSITY>(), //physViscosity
    ( T )   params.get<parameters::PHYS_CHAR_DENSITY>()    //physDensity
  );
  sLatticeAD.setUnitConverter<UnitConverterFromResolutionAndLatticeVelocity<T,D3Q7<VELOCITY>>>(
    (int)   params.get<parameters::RESOLUTION>(),        //resolution
    ( T )   params.get<parameters::LATTICE_CHAR_VELOCITY>(),         //charLatticeVelocity
    ( T )   params.get<parameters::PHYS_CHAR_LENGTH>(),  //charPhysLength
    ( T )   params.get<parameters::PHYS_CHAR_VELOCITY>(),  //charPhysVelocity
    ( T )   params.get<parameters::PHYS_CHAR_VISCOSITY>(), //physViscosity
    ( T )   params.get<parameters::PHYS_CHAR_DENSITY>()    //physDensity
  );

  auto& converter = sLattice.getUnitConverter();

  auto bulkIndicator = geometry.getMaterialIndicator({1,3,4,5});
  // dynamics for fluid
  dynamics::set<BGKdynamics>(sLattice, bulkIndicator);

  // boundary conditions for fluid
  boundary::set<boundary::BounceBack>(sLattice, geometry, 2);
  boundary::set<boundary::InterpolatedPressure>(sLattice, geometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, geometry, 4);
  boundary::set<boundary::InterpolatedVelocity>(sLattice, geometry, 5);

  // dynamics for adsorptive
  auto bulkIndicatorAD = geometry.getMaterialIndicator({1,3,5});
  dynamics::set<ParticleAdvectionDiffusionBGKdynamics>(sLatticeAD, bulkIndicatorAD);

  // boundary conditions for adsorptive
  boundary::set<boundary::BounceBack>(sLatticeAD, geometry, 2);
  boundary::set<boundary::BounceBack>(sLatticeAD, geometry, 4);
  boundary::set<boundary::AdvectionDiffusionDirichlet>(sLatticeAD, geometry, 3);
  setZeroGradientBoundary<T,D3Q7<VELOCITY>>(sLatticeAD, geometry.getMaterialIndicator({5}));

  // set parameters
  sLattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  sLatticeAD.template setParameter<descriptors::OMEGA>(
    converter.getLatticeRelaxationFrequencyFromDiffusivity<D3Q7<VELOCITY>>(
      params.get<parameters::DIFFUSION>()));

  // define lattice coupling operator
  auto& coupling = myCase.setCouplingOperator(
    "NavierStokesAdvectionDiffusionCoupling",
    NavierStokesAdvectionDiffusionVelocityCoupling{},
    names::NavierStokes{},   sLattice,
    names::Concentration0{}, sLatticeAD
  );
  coupling.restrictTo(geometry.getMaterialIndicator({1}));
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& sLatticeAD = myCase.getLattice(Concentration0{});

  // initialisation for fluid
  AnalyticalConst3D<T, T> u0(0., 0., 0.);

  // initialisation for adsorptive
  auto initIndicatorAD = geometry.getMaterialIndicator({1,2,4,5});
  AnalyticalConst3D<T, T> rhoSmall(1.e-8);
  momenta::setDensity(sLatticeAD, initIndicatorAD, rhoSmall);
  sLatticeAD.iniEquilibrium(initIndicatorAD, rhoSmall, u0);

  sLattice.initialize();
  sLatticeAD.initialize();
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT) {
  using T = MyCase::value_t;
  auto& params = myCase.getParameters();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& converter = sLattice.getUnitConverter();

  std::vector<T> maxVelocity(3, T());
  const T distanceToBoundary = converter.getConversionFactorLength() / 2.;
  const T velNS = converter.getCharPhysVelocity();
  const std::size_t itStartTime = converter.getLatticeTime(params.get<parameters::PHYS_START_T>());

  if (iT <= itStartTime && iT % 50 == 0) {
    SinusStartScale<T,int> startScale(itStartTime, T(1));
    int help[1] = {(int) iT};
    T frac[3] = {T()};
    startScale(frac, help);

    // set lattice velocity on boundary
    maxVelocity[1] = velNS * frac[0];
    RectanglePoiseuille3D<T> u5(geometry, 5, maxVelocity,
                                distanceToBoundary, distanceToBoundary, distanceToBoundary);
    momenta::setVelocity(sLattice, geometry.getMaterialIndicator({5}), u5);
  }

  const int itStartPeriodTime = converter.getLatticeTime(params.get<parameters::PHYS_START_PERIOD>());
  const T amplitude = params.get<parameters::AMPLITUDE_PHYS_PRESSURE>() / converter.getConversionFactorPressure();

  T rho = 1.;
  if ((iT <= itStartTime + 0.5 * itStartPeriodTime) && (iT > itStartTime)) {
    Cosinus<T,T> cos(params.get<parameters::PHYS_START_PERIOD>(), T(0.5) * amplitude);
    T help[1] = {converter.getPhysTime(iT - itStartTime)};
    T frac[1] = {T()};
    cos(frac, help);

    rho = util::densityFromPressure<T,D3Q19<>>(T(-0.5) * amplitude + frac[0]);
    AnalyticalConst3D<T,T> rhovar(converter.getPhysDensity(rho));
    momenta::setDensity(sLattice, geometry.getMaterialIndicator({4}), rhovar);
  }

  if (iT > itStartTime + 0.5 * itStartPeriodTime)
  {
    CosinusComposite<T,T> cosComp(params.get<parameters::PHYS_PERIOD>(), amplitude, params.get<parameters::DIFFERENCE_PERIOD>());
    T help[1] = { converter.getPhysTime(iT - itStartTime - 0.5 * itStartPeriodTime) };
    T frac[1] = { T() };
    cosComp(frac, help);

    rho = util::densityFromPressure<T,D3Q19<>>(-frac[0]);
    AnalyticalConst3D<T,T> rhovar(converter.getPhysDensity(rho));
    momenta::setDensity(sLattice, geometry.getMaterialIndicator({4}), rhovar);
  }
}

void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT,
                util::TimeIntegrator<MyCase::value_t>& density,
                util::TimeIntegrator<MyCase::value_t>& variance)
{
  using T = MyCase::value_t;
  auto& params = myCase.getParameters();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& sLatticeAD = myCase.getLattice(Concentration0{});
  auto& geometry = myCase.getGeometry();
  auto& converter = sLattice.getUnitConverter();

  SuperLatticePhysVelocity3D<T,D3Q19<>> velocityNS(sLattice, converter);
  SuperLatticePhysPressure3D<T,D3Q19<>> pressureNS(sLattice, converter);
  SuperLatticeDensity3D<T,D3Q7<VELOCITY>> adsorptive(sLatticeAD);

  // activate this block if you want more simulation output (timer, stability, VTK)
  if constexpr (false) {
    const unsigned vtkSaveT = converter.getLatticeTime(1);
    const unsigned statIter = converter.getLatticeTime(1);
    SuperVTMwriter3D<T> vtmWriter("microMixer3d");
    vtmWriter.addFunctor(velocityNS);
    vtmWriter.addFunctor(adsorptive);
    if (iT == 0) {
      vtmWriter.createMasterFile();
    }

    if (iT % vtkSaveT == 0) {
      sLattice.setProcessingContext(ProcessingContext::Evaluation);
      sLatticeAD.setProcessingContext(ProcessingContext::Evaluation);
      vtmWriter.write(iT);
    }

    if (iT % statIter == 0) {
      timer.update(iT);
      timer.printStep();
      sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    }
  }

  const T time = converter.getPhysTime(iT);
  const T dt = converter.getPhysDeltaT();
  const T physMaxTime = params.get<parameters::PHYS_MAX_T>();
  const T physPeriod = params.get<parameters::PHYS_PERIOD>();

  if (iT == 0) {
    density.reset(physMaxTime - T(2) * physPeriod, physMaxTime - physPeriod, dt);
    variance.reset(physMaxTime - physPeriod, physMaxTime, dt);
  }

  if (time >= physMaxTime - 2.0 * physPeriod - dt) {

    int input[1] = { };
    T output[adsorptive.getTargetDim()+1];  // for concentration

    // average concentration at the outlet
    SuperAverage3D<T>(adsorptive, geometry, 5).operator()(output, input);
    density.takeValue(iT, output[0] / physPeriod);
  }

  if (time >= physMaxTime - 1.0 * physPeriod - dt) {

    int input[1] = { };
    T output1[1+1];  // for variance

    // variance of the density at the outlet
    const T mu = density.getResult();
    SuperConst3D<T> expectedValue(geometry, mu);
    SuperAverage3D<T>((adsorptive - expectedValue) * (adsorptive - expectedValue), geometry, 5).operator()(output1, input);
    variance.takeValue(iT, output1[0] / physPeriod);
  }
}

void simulate(MyCase& myCase,
  util::TimeIntegrator<MyCase::value_t>& density,
  util::TimeIntegrator<MyCase::value_t>& variance)
{
  using T = MyCase::value_t;
  OstreamManager clout(std::cout, "simulate");
  auto& params = myCase.getParameters();

  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(
    params.get<parameters::PHYS_MAX_T>());
  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (size_t iT = 0; iT < iTmax; ++iT) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    myCase.getLattice(NavierStokes{}).collideAndStream();
    myCase.getOperator("NavierStokesAdvectionDiffusionCoupling").apply();
    myCase.getLattice(Concentration0{}).collideAndStream();

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, timer, iT, density, variance);
  }

  timer.stop();
  timer.printSummary();
}

void setInitialControl(MyOptiCase& optiCase) {
  optiCase.getController().set({0.6, 3., 0.5});
}

void applyControl(MyOptiCase& optiCase) {
  auto& control = optiCase.getController();
  auto& params = optiCase.getCase(Controlled{}).getParameters();
  params.set<parameters::PHYS_PERIOD>(control[0]);
  params.set<parameters::AMPLITUDE_PHYS_PRESSURE>(control[1]);
  params.set<parameters::DIFFERENCE_PERIOD>(control[2]);

  OstreamManager clout(std::cout, "applyControl");
  clout << "Control: period length = " << std::setprecision (12) << control[0]
    << ", amplitude = " << control[1]
    << ", akuteness = " << control[2] << std::endl;
}

/// Return time-averaged std. deviation of the adsorptive density
MyCase::value_t objectiveF(MyOptiCase& optiCase,
  util::TimeIntegrator<MyCase::value_t>& density,
  util::TimeIntegrator<MyCase::value_t>& variance)
{
  using T = MyCase::value_t;
  OstreamManager clout (std::cout, "objectiveF");
  clout << "Optimize by segregation intensity." << std::endl;

  const T refVariance = density.getResult() * (T(1) - density.getResult());
  const T segIntensity = variance.getResult() / refVariance;
  clout << "Average density                  = " << std::setprecision (12) << T(density.getResult()) << std::endl;
  clout << "Average variance                 = " << std::setprecision (12) << T(variance.getResult()) << std::endl;
  clout << "Danckwerts reference variance    = " << std::setprecision (12) << refVariance << std::endl;
  clout << "Danckwerts segregation intensity = " << std::setprecision (12) << segIntensity << std::endl << std::endl;

  return segIntensity;
}

MyCase::value_t computeObjective(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  auto& controlledCase = optiCase.getCase(Controlled{});

  // Set control variables
  controlledCase.resetLattices();
  applyControl(optiCase);

  // initialize auxiliary variables for objective computation
  util::TimeIntegrator<T> density(0.0, 1.0, 0.1);
  util::TimeIntegrator<T> variance(0.0, 1.0, 0.1);

  // Prepare Case
  prepareGeometry(controlledCase);
  prepareLattice(controlledCase);
  setInitialValues(controlledCase);
  simulate(controlledCase, density, variance);

  // Evaluate objective functor to compute objective value
  return objectiveF(optiCase, density, variance);
}

int main(int argc, char* argv[])
{
  using T = MyCase::value_t;
  initialize(&argc, &argv);

  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<LATTICE_CHAR_VELOCITY>(0.1);
    myCaseParameters.set<PHYS_CHAR_VELOCITY   >(0.0332);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY  >(1.e-6);
    myCaseParameters.set<PHYS_CHAR_LENGTH     >(0.00133);
    myCaseParameters.set<PHYS_MAX_T           >(5.);   // time for fluid simulation
    myCaseParameters.set<PHYS_START_T         >(0.6);  // time to start fluid pulsation
    myCaseParameters.set<PHYS_START_PERIOD    >(0.4);  // time to start fluid pulsation
    myCaseParameters.set<PHYS_CHAR_DENSITY    >(1000);
    myCaseParameters.set<DIFFUSION            >(1.e-7);  // 1.e-9
    myCaseParameters.set<RESOLUTION           >(7);    // resolution of the hydraulic diameter  // 36
    myCaseParameters.set<parameters::DX>([](MyCase::ParametersD& p){
      return p.get<PHYS_CHAR_LENGTH>()
           / p.get<RESOLUTION>();
    });
  }

  Mesh mesh = createMesh(myCaseParameters);

  MyCase myCase(myCaseParameters, mesh);

  MyOptiCase optiCase;
  optiCase.setCase<Controlled>(myCase);

  prepareGeometry(myCase);
  prepareLattice(myCase);
  setInitialValues(myCase);
  setInitialControl(optiCase);
  // uncomment if a single simulation is desired
  // simulate(myCase);

  optiCase.setObjective(computeObjective);

  OptimizerParametersD<T> optiParameters;
  {
    using namespace parameters;
    optiParameters.set<optimizer::EPS>(1.e-7);
    optiParameters.set<optimizer::MAX_OPTIMIZATION_STEPS>(2);
    optiParameters.set<linesearch::INITIAL_STEP_SIZE>(0.5);
    optiParameters.set<linesearch::MAX_STEP_ATTEMPTS>(10);
    optiParameters.set<LBFGS::INITIAL_HESSIAN_SCALING>(1.e-6);
    optiParameters.set<optimizer::LOG_FILE_NAME>("log");
    optiParameters.set<optimizer::WITH_UPPER_BOUND>(true);
    optiParameters.set<optimizer::UPPER_BOUND>(2.0);
    optiParameters.set<optimizer::WITH_LOWER_BOUND>(true);
    optiParameters.set<optimizer::LOWER_BOUND>(0.1);
    optiParameters.set<optimizer::CONTROL_EPS>(0.);
  }

  OptimizerLBFGS<T,std::vector<T>> optimizer(
    optiCase.getController().size(), optiParameters,
    {OptimizerLogType::objective, OptimizerLogType::control, OptimizerLogType::derivative});

  optimizer.optimize(optiCase);
}
