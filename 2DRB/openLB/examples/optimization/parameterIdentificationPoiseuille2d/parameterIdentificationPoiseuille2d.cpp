/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2007, 2012, 2022 Jonas Latt, Mathias J. Krause
 *  Vojtech Cvrcek, Peter Weisbrod, Julius Jeßberger
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

/** \file A simple two-dimensional fluid flow optimization problem is solved.
 * The setup is a planar channel flow, similar to the example poiseuille2d.
 * For a given pressure drop, the steady flow is simulated and the mass flow
 * rate is computed. In the optimization problem, the inlet pressure is
 * determined s.t. a pre-defined mass flow rate is achieved.
 */

// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::opti;

namespace olb::parameters {

struct WANTED_MASS_FLOW : public descriptors::FIELD_BASE<1> { };
struct REFERENCE_INLET_PRESSURE : public descriptors::FIELD_BASE<1> { };
struct INLET_PRESSURE : public descriptors::FIELD_BASE<1> { };

}

using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D2Q9<>>
>;
using MyADfCase = MyCase::exchange_value_t<util::ADf<double,1>>;

using MyOptiCase = OptiCaseADf<
  Controlled, MyCase,
  Derivatives, MyADfCase
>;

template <typename PARAMETERS>
auto createMesh(PARAMETERS& parameters) {
  using T = PARAMETERS::value_t;
  const Vector extent = parameters.template get<parameters::DOMAIN_EXTENT>();
  const Vector origin {0., 0.};
  IndicatorCuboid2D<T> cuboid(extent, origin);

  Mesh<T,2> mesh(cuboid,
                 extent[1] / parameters.template get<parameters::RESOLUTION>(),
                 singleton::mpi().getSize());
  mesh.setOverlap(parameters.template get<parameters::OVERLAP>());
  return mesh;
}

template<typename CASE>
void prepareGeometry(CASE& myCase)
{
  using T = CASE::value_t;
  auto& parameters = myCase.getParameters();
  const Vector domain_extent = parameters.template get<parameters::DOMAIN_EXTENT>();
  const int N = parameters.template get<parameters::RESOLUTION>();
  auto& superGeometry = myCase.getGeometry();

  superGeometry.rename(0, 2);
  superGeometry.rename(2, 1, {1,1});

  const T physSpacing = domain_extent[1] / N;
  const Vector<T,2> extent    {physSpacing / T(2), domain_extent[1]};
  const Vector<T,2> originIn  {-physSpacing / T(4), 0};
  const Vector<T,2> originOut {domain_extent[0]-physSpacing / T(4), 0};

  IndicatorCuboid2D<T> inflow(extent, originIn);
  superGeometry.rename(2, 3, 1, inflow);

  IndicatorCuboid2D<T> outflow(extent, originOut);
  superGeometry.rename(2, 4, 1, outflow);

  superGeometry.clean(false);
  superGeometry.innerClean(false);
  superGeometry.checkForErrors(false);
}

template<typename CASE>
void prepareLattice(CASE& myCase)
{
  using T = CASE::value_t;
  using DESCRIPTOR = CASE::descriptor_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  const int N = parameters.template get<parameters::RESOLUTION>();
  const T Re = parameters.template get<parameters::REYNOLDS>();
  sLattice.template setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    int {N},     // resolution: number of voxels per charPhysL
    (T)   0.8,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   1,     // charPhysLength: reference length of simulation geometry
    (T)   1,     // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1./Re, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0    // physDensity: physical density in __kg / m^3__
  );
  auto& converter = sLattice.getUnitConverter();
  const T omega = converter.getLatticeRelaxationFrequency();

  dynamics::set<BGKdynamics>(sLattice, superGeometry.getMaterialIndicator({1}));
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);

  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);

  sLattice.template setParameter<descriptors::OMEGA>(omega);
}

template <typename CASE>
void setInitialValues(CASE& myCase) {
  using T = CASE::value_t;
  using DESCRIPTOR = CASE::descriptor_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  auto& converter = sLattice.getUnitConverter();
  auto& parameters = myCase.getParameters();
  const Vector domain_extent = parameters.template get<parameters::DOMAIN_EXTENT>();
  const T inletPressure = parameters.template get<parameters::INLET_PRESSURE>();
  std::cout << "starting simulation with pressure = " << inletPressure << "..." << std::endl;

  // Initial conditions
  AnalyticalLinear2D<T,T> rho(
    - inletPressure / domain_extent[0] * descriptors::invCs2<T,DESCRIPTOR>(),
    0,
    inletPressure * descriptors::invCs2<T,DESCRIPTOR>() + 1 );

  const T Lx = converter.getLatticeLength( domain_extent[0] ) - 1;
  const T Ly = converter.getLatticeLength( domain_extent[1] ) - 1;
  const T maxLatticeVelocity = inletPressure * Ly * Ly
    / (8.0 * converter.getLatticeViscosity() * Lx);
  const T maxVelocity = converter.getPhysVelocity(maxLatticeVelocity);

  const T radius = T(0.5) * (domain_extent[1] - converter.getPhysDeltaX());
  const std::vector<T> axisPoint { domain_extent[0]/T(2), domain_extent[1]/T(2) };
  const std::vector<T> axisDirection { 1, 0 };
  Poiseuille2D<T> u( axisPoint, axisDirection, maxVelocity, radius );

  const auto domain = superGeometry.getMaterialIndicator({1,2,3,4});
  //sLattice.defineRhoU( domain, rho, u );
  momenta::setVelocity(sLattice, domain, u);
  momenta::setDensity(sLattice, domain, rho);
  //sLattice.iniEquilibrium( domain, rho, u );

  sLattice.initialize();
}

template<typename CASE, typename TIMER>
void getResults(CASE& myCase, TIMER& timer, std::size_t iT, bool hasConverged) {
  using T = CASE::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& converter = sLattice.getUnitConverter();
  const T maxPhysT = myCase.getParameters().template get<parameters::PHYS_MAX_T>();

  SuperVTMwriter2D<T> vtmWriter( "poiseuille2d" );
  const bool lastTimeStep
   = ( hasConverged || (iT + 1 == converter.getLatticeTime( maxPhysT )) );
  const int vtmIter  = converter.getLatticeTime( maxPhysT/20. );

  if (iT == 0) {
      vtmWriter.createMasterFile();
  }

  if ( iT%vtmIter==0 || lastTimeStep )
  {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    timer.update( iT );
    timer.printStep();
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
    vtmWriter.write( iT );
  }
}

template <typename CASE>
void simulate(CASE& myCase) {
  using T = CASE::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  auto& converter = sLattice.getUnitConverter();
  const T maxPhysT = myCase.getParameters().template get<parameters::PHYS_MAX_T>();

  util::Timer<T> timer(converter.getLatticeTime(maxPhysT),
                       superGeometry.getStatistics().getNvoxel());
  util::ValueTracer<T> converge(converter.getLatticeTime(0.25), 1e-9);
  timer.start();
  for (std::size_t iT=0 ; iT < converter.getLatticeTime(maxPhysT); ++iT ) {
    if ( converge.hasConverged() ) {
      std::cout << "Simulation converged." << std::endl;
      getResults(myCase, timer, iT, converge.hasConverged());
      break;
    }
    sLattice.collideAndStream();
    converge.takeValue(sLattice.getStatistics().getMaxU(), false );
  }
  timer.stop();
  timer.printSummary();
}

void setInitialControls(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  std::vector<T> control({optiCase.getCase(Controlled{}).getParameters().template get<parameters::INLET_PRESSURE>()});
  optiCase.getController().set(control);
}

template <typename CASE, typename CONTROLLER>
void applyControls(CASE& myCase, CONTROLLER& controller) {
  auto control = controller.get();
  myCase.getParameters().template set<parameters::INLET_PRESSURE>(control[0]);
}

/// Compute squared error between simulated and wanted mass flow rate
template <typename CASE>
CASE::value_t massFlowError(CASE& myCase)
{
  using T = CASE::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  auto& converter = sLattice.getUnitConverter();
  auto& parameters = myCase.getParameters();
  const Vector domain_extent = parameters.template get<parameters::DOMAIN_EXTENT>();
  sLattice.setProcessingContext(ProcessingContext::Evaluation);

  SuperLatticeVelocity2D velocity(sLattice);
  SuperLatticeDensity2D density(sLattice);
  Vector<T,2> tmp{0, 1};
  SuperPlaneIntegralFluxMass2D<T> massFlowRate(
    velocity, density, superGeometry, converter.getConversionFactorMass(),
    converter.getConversionFactorTime(), Vector<T,2>({T(0.5)*domain_extent[0],
                                                      T(0.5)*domain_extent[1]}),
    tmp, BlockDataReductionMode::Analytical
  );
  const int input[3] = {0};
  T mFlow[4] = {0.};
  massFlowRate(mFlow, input);
  std::cout << "Mass flow rate = " << mFlow[0] << std::endl;
  const T res = mFlow[0];
  const T wantedRes = parameters.template get<parameters::WANTED_MASS_FLOW>();
  return 0.5 * (res - wantedRes) * (res - wantedRes);
}

template <typename NAME>
auto computeObjective(MyOptiCase& optiCase, NAME) {
  using T = MyOptiCase::case_t<NAME>::value_t;
  auto& myCase = optiCase.getCase(NAME{});
  auto& controller = optiCase.getController<T>();

  myCase.resetLattices();
  applyControls(myCase, controller);
  prepareGeometry(myCase);
  prepareLattice(myCase);
  setInitialValues(myCase);
  simulate(myCase);
  return massFlowError(myCase);
}

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  initialize( &argc, &argv );

  // === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParametersD;
  {
    using namespace parameters;
    myCaseParametersD.set<RESOLUTION              >(        50);
    myCaseParametersD.set<DOMAIN_EXTENT           >(  {2., 1.});
    myCaseParametersD.set<REYNOLDS                >(       10.);
    myCaseParametersD.set<PHYS_MAX_T              >(       30.);
    myCaseParametersD.set<WANTED_MASS_FLOW        >(0.00026159);
    myCaseParametersD.set<REFERENCE_INLET_PRESSURE>(  0.000659);
    myCaseParametersD.set<INLET_PRESSURE          >(    0.0001);
  }
  auto myADfCaseParametersD = myCaseParametersD.copyAs<MyADfCase::value_t>();

  // === Step 3: Create Mesh ===
  auto mesh = createMesh(myCaseParametersD);
  auto ADfmesh = createMesh(myADfCaseParametersD);

  // === Step 4: Create Case ===
  MyCase myCase(myCaseParametersD, mesh);
  MyADfCase myADfCase(myADfCaseParametersD, ADfmesh);

  // Create OptiCase
  MyOptiCase optiCase;

  // Set Case
  optiCase.setCase<Controlled>(myCase);
  optiCase.setCase<Derivatives>(myADfCase);

  // Define initial values for the control
  setInitialControls(optiCase);

  // Define objective computation routine
  optiCase.setObjective([](auto& optiCase, auto name) {
    return computeObjective(optiCase, name);
  });

  OptimizerParametersD<MyCase::value_t> optiParameters;
  {
    using namespace parameters;
    optiParameters.set<optimizer::EPS>(1.e-7);
    optiParameters.set<optimizer::MAX_OPTIMIZATION_STEPS>(20);
    optiParameters.set<linesearch::MAX_STEP_ATTEMPTS>(10);
    optiParameters.set<optimizer::LOG_FILE_NAME>("log");
    optiParameters.set<optimizer::WITH_UPPER_BOUND>(true);
    optiParameters.set<optimizer::UPPER_BOUND>(0.01);
    optiParameters.set<optimizer::WITH_LOWER_BOUND>(true);
    optiParameters.set<optimizer::CONTROL_EPS>(0.);
  }

  OptimizerLBFGS<MyCase::value_t,std::vector<MyCase::value_t>> optimizer(
    1, optiParameters,
    {OptimizerLogType::objective, OptimizerLogType::control, OptimizerLogType::derivative});
  optimizer.optimize(optiCase);
}
