/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Felix Schuhmann, Shota Ito
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

/* ellipsoidShape3d.cpp:
 * This example examines a steady flow past an ellipsoid placed inside a rectangular channel.
 * The half-axes are alligned along the x, y and z coordinate directions and periodic boundary treatment
 * is applied in y and z directions.
 *
 * Contains:
 * - Standard channel flow simulation
 * - Flow simulation while computing derivatives with respect to ellipsoid radii in y and z direction using ADf
 * - Shape optimization using LBFGS:
 *   radii in y and z direction are optimized to minimize L2-Norm of dissipation over the entire channel
 *   radius in x direction is adjusted to enforce fixed ellipsoid volume
 */

// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include "olb.h"

using namespace olb;
using namespace olb::names;
using namespace olb::opti;

using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q19<>>
>;
using MyADfCase = Case<
  NavierStokes, Lattice<util::ADf<double,2>, descriptors::D3Q19<>>
>;
using MyOptiCase = OptiCaseADf<
  Controlled, MyCase,
  Derivatives, MyADfCase
>;

namespace olb::parameters{

struct CONTROLS : public descriptors::FIELD_BASE<2>{ };
struct ELLIPSOID_VOLUME : public descriptors::FIELD_BASE<1>{ };
struct ELLIPSOID_POS : public descriptors::FIELD_BASE<0,1>{ };
struct SMOOTH_LAYER_THICKNESS : public descriptors::FIELD_BASE<1>{ };

}

template <typename PARAMETERS>
auto createMesh(PARAMETERS& parameters) {
  using T = PARAMETERS::value_t;
  const Vector extent = parameters.template get<parameters::DOMAIN_EXTENT>();
  const Vector origin{0., 0., 0.};
  IndicatorCuboid3D<T> channel( extent, origin );

  Mesh<T,3> mesh(channel,
                 extent[1] / parameters.template get<parameters::RESOLUTION>(),
         singleton::mpi().getSize());
  mesh.setOverlap(parameters.template get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({false, true, true});
  return mesh;
}

template<typename CASE>
void prepareGeometry(CASE& myCase)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;
  using T = CASE::value_t;
  auto& parameters = myCase.getParameters();
  auto& superGeometry = myCase.getGeometry();
  const T delta = parameters.template get<parameters::DOMAIN_EXTENT>()[1] /
                  parameters.template get<parameters::RESOLUTION>() / 2.0;
  const Vector extent = parameters.template get<parameters::DOMAIN_EXTENT>();
  const Vector origin_channel{0., 0., 0.};
  IndicatorCuboid3D<T> channel(extent, origin_channel);
  superGeometry.rename(0, 1, channel);

  auto min = channel.getMin();
  auto max = channel.getMax();
  Vector<T,3> origin { min[0]-delta, min[1]-delta, min[2]-delta};
  Vector<T,3> extend { delta * 2.0, max[1] + delta*2.0, max[2] + delta*2. };
  IndicatorCuboid3D<T> inflow(extend, origin);
  superGeometry.rename(1,3,inflow);

  Vector<T,3> origin2 { max[0]-delta, min[1]-delta, min[2]-delta};
  Vector<T,3> extend2 { delta * 2.0, max[1] + delta*2.0, max[2] + delta*2. };
  IndicatorCuboid3D<T> outflow(extend2, origin2);
  superGeometry.rename(1,4,outflow);

  superGeometry.clean();
  superGeometry.checkForErrors();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

template<typename CASE>
void prepareLattice(CASE& myCase) {
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;
  using T = CASE::value_t;
  using DESCRIPTOR = CASE::descriptor_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  sLattice.template setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    parameters.template get<parameters::RESOLUTION>(),                 // resolution: number of voxels per charPhysL
    parameters.template get<parameters::LATTICE_RELAXATION_TIME>(),    // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    parameters.template get<parameters::DOMAIN_EXTENT>()[1],           // charPhysLength: reference length of simulation geometry
    parameters.template get<parameters::PHYS_CHAR_VELOCITY>(),         // physVelocity: maximal/highest expected velocity during simulation in __m / s__
    parameters.template get<parameters::PHYS_CHAR_VISCOSITY>(),        // charPhysViscosity: physical kinematic viscosity in __m^2 / s__
    parameters.template get<parameters::PHYS_CHAR_DENSITY>()           // physDensity: physical density in __kg / m^3__
  );
  auto& converter = sLattice.getUnitConverter();
  converter.print();
  const T omega = converter.getLatticeRelaxationFrequency();

  // Bulk dynamics for all materials
  auto bulkIndicator = superGeometry.getMaterialIndicator({1});
  dynamics::set<PorousBGKdynamics>(sLattice, bulkIndicator);

  // boundary conditions
  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);

  sLattice.template setParameter<descriptors::OMEGA>(omega);
  clout << "Prepare Lattice ... OK" << std::endl;
}

template <typename CASE>
void setInitialValues(CASE& myCase) {
  using T = CASE::value_t;
  auto& parameters = myCase.getParameters();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto bulkIndicator = myCase.getGeometry().getMaterialIndicator( { 1,3,4,5 } );

  // Set required porosities
  Vector radiusYZ = parameters.template get<parameters::CONTROLS>();
  const T radiusX = 0.75*parameters.template get<parameters::ELLIPSOID_VOLUME>() /
        (M_PI*radiusYZ[0]*radiusYZ[1]);
  Vector radius{radiusX, radiusYZ[0], radiusYZ[1]};
  Vector radius_plus_eps = radius;
  const T eps = parameters.template get<parameters::SMOOTH_LAYER_THICKNESS>() *
            parameters.template get<parameters::DOMAIN_EXTENT>()[1] /
            parameters.template get<parameters::RESOLUTION>();
  for (int i=0; i < 3; ++i) {
    radius_plus_eps[i] += eps;
  }
  std::shared_ptr<IndicatorF3D<T>> ellipsoid =
    std::make_shared<IndicatorEllipsoid3D<T>>(parameters.template get<parameters::ELLIPSOID_POS>(), radius_plus_eps);

  // Smooth ellipsoid to set porosity field
  std::shared_ptr<AnalyticalF3D<T,T>> smoothEllipsoid =
    std::make_shared<SmoothIndicatorSigmoidEllipsoid3D<T,T,false>>(parameters.template get<parameters::ELLIPSOID_POS>(), radius, eps);
  std::shared_ptr<AnalyticalF3D<T,T>> one = std::make_shared<AnalyticalConst3D<T,T>>(1.);
  std::shared_ptr<AnalyticalF3D<T,T>> ellipsoidField = one - smoothEllipsoid;

  // Initialize porosity
  fields::set<descriptors::POROSITY>(sLattice, bulkIndicator, *ellipsoidField);

  sLattice.initialize();
}

template<typename CASE>
void setBoundaryValues(CASE& myCase, int iT)
{
  using T = CASE::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  auto& converter = sLattice.getUnitConverter();
  auto& parameters = myCase.getParameters();
  const T rampStartT = parameters.template get<parameters::PHYS_START_T>();
  const T rampUpdateT = parameters.template get<parameters::PHYS_BOUNDARY_VALUE_UPDATE_T>();
  int iTmaxStart = converter.getLatticeTime(rampStartT);
  int iTupdate = converter.getLatticeTime(rampUpdateT);

  if (iT%iTupdate == 0 && iT <= iTmaxStart) {
    PolynomialStartScale<T,int> startScale(iTmaxStart, T(1));
    int iTvec[1] = {iT};
    T frac[1] = {};
    startScale(frac, iTvec);
    Vector<T,3> u_free {frac[0]*converter.getCharPhysVelocity(), 0.0, 0.0};
    AnalyticalConst3D<T,T> u(u_free);
    momenta::setVelocity(sLattice, superGeometry.getMaterialIndicator({3}), u);

    sLattice.template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

template<typename CASE>
void getResults(CASE& myCase,
                int iT,
                util::Timer<typename CASE::value_t>& timer)
{
  OstreamManager clout(std::cout,"getResults");
  using T = CASE::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& converter = sLattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  const T maxPhysT = parameters.template get<parameters::PHYS_MAX_T>();
  const int statIter = converter.getLatticeTime( maxPhysT / 10 );

  // Writes output on the console
  if ( iT%statIter == 0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print(iT, converter.getPhysTime( iT ));
  }
}

template<typename CASE>
void getOptiResults(CASE& myCase,
                    int iT)
{
  using T = CASE::value_t;
  using DESCRIPTOR = CASE::descriptor_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& converter = sLattice.getUnitConverter();

  SuperVTMwriter3D<T> vtmWriter("ellipsoidOpti3d");
  SuperLatticePhysVelocity3D velocity(sLattice, converter);
  SuperLatticePhysPressure3D pressure(sLattice, converter);
  SuperLatticePorosity3D porosity(sLattice);
  SuperLatticePhysDissipation3D viscous_dissipation(sLattice, converter);
  SuperLatticeFfromCallableF porous_dissipation(sLattice, [&](T* output, auto cell){
    T uTemp[DESCRIPTOR::d];
    cell.computeU(uTemp);
    const T porosity = cell.template getField<descriptors::POROSITY>();

    const T invPermeability = projection::porosityToInvPermeability(porosity, converter);
    const T uNormSq = util::euklidN2(uTemp, DESCRIPTOR::d);
    output[0] = converter.getPhysViscosity() * invPermeability * uNormSq;
  });

  viscous_dissipation.getName() = "viscous dissipation";
  porous_dissipation.getName() = "porous dissipation";

  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );
  vtmWriter.addFunctor( porosity );
  vtmWriter.addFunctor( viscous_dissipation );
  vtmWriter.addFunctor( porous_dissipation );

  if ( iT==0 ) {
    vtmWriter.createMasterFile();
  }

  // Writes the vtk files, currently only for last simulation
  vtmWriter.write(iT);
}

template <typename CASE>
void simulate(CASE& myCase) {
  using T = CASE::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  auto& converter = sLattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  // === 4th Step: Main Loop with Timer ===
  const T maxPhysT = parameters.template get<parameters::PHYS_MAX_T>();
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  auto maxiT = converter.getLatticeTime( maxPhysT );
  timer.start();

  for (std::size_t iT = 0; iT < maxiT; ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions === can be skipped for forced flow type
    setBoundaryValues(myCase, iT);

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults(myCase, iT, timer );
  }
  timer.stop();
  timer.printSummary();
}

template <typename CASE>
void setInitialControls(MyOptiCase& optiCase) {
  Vector controls = optiCase.getCase(Controlled{}).getParameters().template get<parameters::CONTROLS>();
  optiCase.getController().set({controls[0], controls[1]});
}

template <typename CASE, typename CONTROLLER>
void applyControls(CASE& myCase, CONTROLLER& controller) {
  OstreamManager clout(std::cout, "applyControls");
  auto controls = controller.get();
  myCase.getParameters().template set<parameters::CONTROLS>({controls[0], controls[1]});
  clout << "Controls: " << myCase.getParameters().template get<parameters::CONTROLS>() << std::endl;
}

template<typename CASE>
CASE::value_t computeDissipation(CASE& myCase) {
  OstreamManager clout(std::cout, "dissipationF");
  using T = CASE::value_t;
  using DESCRIPTOR = CASE::descriptor_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& converter = sLattice.getUnitConverter();
  auto& superGeometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  sLattice.setProcessingContext(ProcessingContext::Evaluation);
  T result = 0.0;
  Vector<T,3> center = parameters.template get<parameters::ELLIPSOID_POS>();
  IndicatorCuboid3D<T> integrationDomain(1.0, 0.4, 0.4, center);
  SuperIndicatorFfromIndicatorF3D<T> iDomain(integrationDomain, superGeometry);

  // Viscous Dissipation
  SuperLatticePhysDissipation3D<T,DESCRIPTOR> viscous_dissipation( sLattice, converter );
  SuperIntegral3D<T,T> dissipationIntegral1(viscous_dissipation, iDomain);
  T vDissipation[1] = {0.};
  int input1[1];
  dissipationIntegral1( vDissipation, input1 );
  clout << "Viscous dissipation = " << vDissipation[0] << std::endl;
  result = vDissipation[0];

  // Porous Dissipation
  SuperLatticeFfromCallableF<T,DESCRIPTOR> porous_Dissipation(sLattice, [&](T* output, auto cell){
    T uTemp[DESCRIPTOR::d];
    cell.computeU(uTemp);
    const T porosity = cell.template getField<descriptors::POROSITY>();

    const T invPermeability = projection::porosityToInvPermeability(porosity, converter);
    const T uNormSq = util::euklidN2(uTemp, DESCRIPTOR::d);
    output[0] = converter.getPhysViscosity() * invPermeability * uNormSq;
  });
  SuperIntegral3D<T,T> dissipationIntegral2(porous_Dissipation, iDomain);
  T pDissipation[1] = {0.};
  dissipationIntegral2( pDissipation, input1 );
  clout << "Porous dissipation = " << pDissipation[0] << std::endl;
  result += pDissipation[0];
  return result*result;
}

template <typename NAME>
auto computeObjective(MyOptiCase& optiCase, NAME) {
  using T = MyOptiCase::case_t<NAME>::value_t;
  auto& myCase = optiCase.getCase(NAME{});
  myCase.resetLattices();
  applyControls(myCase, optiCase.getController<T>());
  prepareGeometry(myCase);
  prepareLattice(myCase);
  setInitialValues(myCase);
  simulate(myCase);
  getOptiResults(myCase, optiCase.getOptimizationStep());
  return computeDissipation(myCase);
}

int main( int argc, char* argv[] ) {
  initialize( &argc, &argv );
  MyCase::ParametersD myCaseParametersD;
  {
    using namespace parameters;
    myCaseParametersD.template set<DOMAIN_EXTENT               >({1.2, 0.4, 0.4});
    myCaseParametersD.template set<RESOLUTION                  >(             20);
    myCaseParametersD.template set<LATTICE_RELAXATION_TIME     >(           0.55);
    myCaseParametersD.template set<PHYS_CHAR_VISCOSITY         >(          0.001);
    myCaseParametersD.template set<REYNOLDS                    >(            50.);
    myCaseParametersD.template set<PHYS_MAX_T                  >(            50.);
    myCaseParametersD.template set<PHYS_START_T                >(            30.);
    myCaseParametersD.template set<PHYS_BOUNDARY_VALUE_UPDATE_T>(           0.01);
    myCaseParametersD.template set<PHYS_CHAR_VELOCITY          >([](MyCase::ParametersD& p) {
      return p.get<PHYS_CHAR_VISCOSITY>() *
             p.get<REYNOLDS>() /
         p.get<DOMAIN_EXTENT>()[1];
    });
    myCaseParametersD.template set<PHYS_CHAR_DENSITY           >(            1.0);
    myCaseParametersD.template set<CONTROLS                    >(   {0.08, 0.08});
    myCaseParametersD.template set<ELLIPSOID_VOLUME            >(          0.001);
    myCaseParametersD.template set<SMOOTH_LAYER_THICKNESS      >(             5.);
    myCaseParametersD.template set<ELLIPSOID_POS               >([](MyCase::ParametersD& p) {
      return FieldD<MyCase::value_t,MyCase::descriptor_t,ELLIPSOID_POS>{
        0.5,
        p.get<DOMAIN_EXTENT>()[1] / 2.0,
    p.get<DOMAIN_EXTENT>()[2] / 2.0};
    });
  }
  auto myADfCaseParametersD = myCaseParametersD.copyAs<MyADfCase::value_t>();

  OptimizerParametersD<MyCase::value_t> optiParameters;
  {
    using namespace parameters;
    optiParameters.set<optimizer::EPS>(1e-16);
    optiParameters.set<optimizer::MAX_OPTIMIZATION_STEPS>(10);
    optiParameters.set<linesearch::INITIAL_STEP_SIZE>(0.01);
    optiParameters.set<linesearch::MAX_STEP_ATTEMPTS>(10);
    optiParameters.set<linesearch::STEP_CONDITION>("StrongWolfe");
    optiParameters.set<optimizer::VERBOSE>(true);
    optiParameters.set<optimizer::LOG_FILE_NAME>("log");
    optiParameters.set<optimizer::WITH_UPPER_BOUND>(true);
    optiParameters.set<optimizer::UPPER_BOUND>(0.19);
    optiParameters.set<optimizer::WITH_LOWER_BOUND>(true);
    optiParameters.set<optimizer::LOWER_BOUND>(0.01);
    optiParameters.set<optimizer::CONTROL_EPS>(0.);
  }

  auto mesh = createMesh(myCaseParametersD);
  auto ADfmesh = createMesh(myADfCaseParametersD);
  MyCase myCase(myCaseParametersD, mesh);
  MyADfCase myADfCase(myADfCaseParametersD, ADfmesh);

  MyOptiCase optiCase;
  optiCase.setCase<Controlled>(myCase);
  optiCase.setCase<Derivatives>(myADfCase);
  setInitialControls<MyCase>(optiCase);

  optiCase.setObjective([](auto& optiCase, auto name) {
    return computeObjective(optiCase, name);
  });

  OptimizerLBFGS<MyCase::value_t,std::vector<MyCase::value_t>> optimizer(
    2, optiParameters,
    {OptimizerLogType::objective, OptimizerLogType::control, OptimizerLogType::derivative}
  );
  optimizer.optimize(optiCase);
}
