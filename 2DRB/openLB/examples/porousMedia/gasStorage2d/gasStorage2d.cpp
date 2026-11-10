/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Alice Raeli, Luiz Czelusniak, Tim Bingert
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

/**
 * This example shows a liquid (water) entering inside a porous media
 * filled with gas (compressed hydrogen at 100bar) and replacing the
 * gas, only to a certain extent as some of the gas remains trapped in
 * the porous domain which gives the effective porosity for such removal
 * of stored gas processes.
 *
 * To run this example, download the appropriate vti file
 *
 * curl "https://openlb.net/data/gas_storage/gasStorage2d.vti" -o gasStorage2d.vti
 *
 * and give the file name as the VTI_INPUT parameter.
 **/

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, D2Q9<RHO,NABLARHO,FORCE,EXTERNAL_FORCE,TAU_EFF,SCALAR>>,
  Component1,  Lattice<double, D2Q9<CONV_POPS,FORCE,VELOCITY,OLD_PHIU,STATISTIC,PSI,NORMGRADPSI,SCALAR,PSI0,THETA,BOUNDARY>>
>;

using NSBulkDynamics = MultiPhaseIncompressibleInterfaceTRTdynamics<MyCase::value_t,MyCase::descriptor_t>;
using ACBulkDynamics = ConservativePhaseFieldBGKdynamics<MyCase::value_t,MyCase::descriptor_t_of<Component1>>;
using Coupling = LiangPostProcessor;

namespace olb::parameters {

struct VTI_INPUT : public descriptors::TYPED_FIELD_BASE<std::string,1> { };
struct ARRAY_NAME : public descriptors::TYPED_FIELD_BASE<std::string,1> { };
struct LATTICE_RELAXATION_TIME_PF : public descriptors::FIELD_BASE<1> { };
struct INLET_LENGTH_F : public descriptors::FIELD_BASE<1> { };
struct OUTLET_LENGTH_F : public descriptors::FIELD_BASE<1> { };
struct SCALE : public descriptors::FIELD_BASE<1> { };
struct INLET_LENGTH : public descriptors::FIELD_BASE<1> { };
struct OUTLET_LENGTH : public descriptors::FIELD_BASE<1> { };
struct C_RHO : public descriptors::FIELD_BASE<1> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T inletLength = params.get<parameters::INLET_LENGTH>();
  const T outletLength = params.get<parameters::OUTLET_LENGTH>();

  // Indicator containing porous medium, inlet and outlet zones
  IndicatorCuboid2D<T> cuboid( extent[0] + inletLength + outletLength, extent[1],
                               { ( extent[0] + outletLength - inletLength )/2., extent[1]/2. }, 0);

  const T dx = extent[1] / params.get<parameters::RESOLUTION>();
  Mesh<T,MyCase::d> mesh(cuboid, dx, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry( MyCase& myCase )
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  const T dx = params.get<parameters::DOMAIN_EXTENT>()[1] / params.get<parameters::RESOLUTION>();
  const Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T deltaRsample = params.get<parameters::SCALE>();
  const T inletLength = params.get<parameters::INLET_LENGTH>();
  const T outletLength = params.get<parameters::OUTLET_LENGTH>();
  std::string vtiFile = params.get<parameters::VTI_INPUT>();
  std::string arrayName = params.get<parameters::ARRAY_NAME>();
  BlockVTIreader2D<T,T> vtiReader( vtiFile, arrayName );
  auto cuboidSample = vtiReader.getCuboid();
  BlockData<2,T,T>& block = vtiReader.getBlockData();

  Vector<int, 2> extentSample     = cuboidSample.getExtent();
  Vector<T, 2>   originSamplePhys = cuboidSample.getOrigin() * deltaRsample;

  Vector<T, 2> extentSamplePhys = {deltaRsample * T(extentSample[0]), deltaRsample * T(extentSample[1])};
  IndicatorBlockData2Dvti<T> indicator(block, extentSamplePhys, originSamplePhys, deltaRsample, false);

  geometry.rename(0, 2);
  geometry.rename(2, 1, {1, 1});

  // internal domain, solid rocks - not fluid
  geometry.rename(1, 5, indicator);

  // inlet and outlet
  IndicatorCuboid2D<T> inlet( dx, extent[1] - dx, { -inletLength, extent[1]/2. }, 0 );
  geometry.rename(2, 3, inlet);

  IndicatorCuboid2D<T> outlet( dx, extent[1] - dx, { outletLength + extent[0] - 0.5*dx, extent[1]/2. }, 0);
  geometry.rename(2, 4, outlet);

  geometry.clean();
  geometry.innerClean();
  geometry.outerClean();
  geometry.checkForErrors();
  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

template <typename STAGE>
void signedDistanceFunction( MyCase& myCase )
{
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();
  const T w = params.get<parameters::INTERFACE_WIDTH>();

  auto& sLatticeAC = myCase.getLattice(Component1{});
  using PFDESCRIPTOR = MyCase::descriptor_t_of<Component1>;

  T max = 1.;
  while ( max <= 1.5*w ) {
    int in[2];
    sLatticeAC.getCommunicator(STAGE{}).communicate();
    sLatticeAC.executePostProcessors(STAGE{});
    SuperLatticeField2D<T, PFDESCRIPTOR, PSI> psi( sLatticeAC );
    SuperMax2D<T,T> Max_psi_(psi, geometry, 1);
    Max_psi_(&max, in);
  }
}

template <typename T>
T helperConvectiveU( MyCase& myCase )
{
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();
  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T outletLength = params.get<parameters::OUTLET_LENGTH>();

  auto& sLatticeNS = myCase.getLattice(NavierStokes{});
  using NSDESCRIPTOR = MyCase::descriptor_t;

  const auto& converter = sLatticeNS.getUnitConverter();
  T dx = converter.getPhysDeltaX();
  IndicatorCuboid2D<T> beforeOutlet_( 1.1*dx, extent[1], { extent[0] + outletLength - 2.*dx, extent[1]/2. }, 0 );
  SuperIndicatorFfromIndicatorF2D<T> beforeOutlet(beforeOutlet_, geometry);
  int in[2];
  T uMax[2];
  SuperLatticeVelocity2D<T, NSDESCRIPTOR> u(sLatticeNS);
  SuperMax2D<T,T> uMax_(u, beforeOutlet);
  uMax_(uMax, in);
  return uMax[0];
}

void prepareLattice( MyCase& myCase )
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& sLatticeNS = myCase.getLattice(NavierStokes{});
  auto& sLatticeAC = myCase.getLattice(Component1{});

  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ACDESCRIPTOR = MyCase::descriptor_t_of<Component1>;

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const int N = params.get<parameters::RESOLUTION>();
  const T tau_mobil = params.get<parameters::LATTICE_RELAXATION_TIME_PF>();
  const T rho_g = params.get<parameters::RHO_VAPOR>();
  const T rho_l = params.get<parameters::RHO_LIQUID>();
  const T nu_g = params.get<parameters::NU_VAPOR>();
  const T nu_l = params.get<parameters::NU_LIQUID>();
  const T tau_l = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T tau_g = nu_g/nu_l*(tau_l-0.5)+0.5;
  const T sigma = params.get<parameters::SURFACE_TENSION>();
  const T w = params.get<parameters::INTERFACE_WIDTH>();
  const T C_rho = params.get<parameters::C_RHO>();

  sLatticeNS.setUnitConverter<MultiPhaseUnitConverterFromRelaxationTime<T,NSDESCRIPTOR>>(
    int   {N},           // resolution
    (T)   tau_l,         // lattice relaxation time
    (T)   rho_l/C_rho,   // lattice density
    (T)   extent[1],     // charPhysLength: reference length of simulation geometry
    (T)   nu_l,          // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   rho_l          // physDensity: physical density in __kg / m^3__
  );

  const auto& converter = sLatticeNS.getUnitConverter();
  converter.print();

  sLatticeAC.setUnitConverter(converter);

  // define lattice Dynamics
  dynamics::set<NSBulkDynamics>(sLatticeNS, geometry.getMaterialIndicator(1));
  dynamics::set<ACBulkDynamics>(sLatticeAC, geometry.getMaterialIndicator(1));

  auto bulk   = geometry.getMaterialIndicator(1);
  auto wall   = geometry.getMaterialIndicator({2, 5});
  auto inlet  = geometry.getMaterialIndicator(3);
  auto outlet = geometry.getMaterialIndicator(4);
  auto fluid  = geometry.getMaterialIndicator({1, 3});
  auto all    = geometry.getMaterialIndicator({0, 1, 2, 3, 4, 5});

  boundary::set<boundary::BounceBackIncompressible>(sLatticeNS, wall);
  boundary::set<boundary::PhaseFieldCurvedWall>(sLatticeAC, wall);
  boundary::set<boundary::RegularizedTemperature>(sLatticeAC, inlet);
  boundary::set<boundary::IncompressibleZouHePressure<T,NSDESCRIPTOR,NSBulkDynamics,momenta::IncRegularizedPressureBoundaryTuple>>(sLatticeNS, inlet);
  boundary::set<boundary::IncompressibleZouHePressure<T,NSDESCRIPTOR,NSBulkDynamics,momenta::IncRegularizedPressureBoundaryTuple>>(sLatticeNS, outlet);
  setConvectivePhaseFieldBoundary<T,ACDESCRIPTOR>(sLatticeAC, outlet);

  sLatticeAC.addPostProcessor<stage::InitOutlet>(outlet, meta::id<SetOutletCells2D<1, 0>>());
  sLatticeAC.addPostProcessor<stage::PreCoupling>(fluid,meta::id<RhoStatistics>());

  sLatticeAC.addPostProcessor<stage::PreCoupling>(meta::id<initialPsi>{});
  sLatticeAC.addPostProcessor<stage::IterativePostProcess>(bulk,meta::id<normGradPsi>{});
  sLatticeAC.addPostProcessor<stage::IterativePostProcess>(meta::id<psiEvolve>{});
  sLatticeAC.setParameter<psiEvolve::DELTAT>(0.7);
  setSignedDistanceBoundary<T,ACDESCRIPTOR>(sLatticeAC, wall);
  setSignedDistanceBoundary<T,ACDESCRIPTOR>(sLatticeAC, inlet);
  setSignedDistanceBoundary<T,ACDESCRIPTOR>(sLatticeAC, outlet);

  auto& coupling = myCase.setCouplingOperator(
      "Coupling",
      LiangPostProcessor {},
      names::NavierStokes {}, sLatticeNS,
      names::Component1 {}, sLatticeAC);
  coupling.restrictTo(geometry.getMaterialIndicator({1, 4}));

  auto& velocityCoupling = myCase.setCouplingOperator(
      "VeloCoupling",
      VelocityCoupling{},
      names::NavierStokes {}, sLatticeNS,
      names::Component1 {}, sLatticeAC);
  velocityCoupling.restrictTo(geometry.getMaterialIndicator({3}));

  coupling.setParameter<Coupling::SIGMA>( sigma / converter.getConversionFactorSurfaceTension() );
  coupling.setParameter<Coupling::W>( w );
  coupling.setParameter<Coupling::TAUS>({tau_g, tau_l});
  coupling.setParameter<Coupling::RHOS>( {rho_g/C_rho, rho_l/C_rho});
  coupling.template setParameter<Coupling::SWITCH>(1);

  sLatticeNS.setParameter<descriptors::OMEGA>( 1. / tau_l );
  T maxRhoGradient = (rho_l-rho_g)/converter.getConversionFactorDensity()/w;
  sLatticeNS.setParameter<collision::ITRT::TAU_MINUS>( T(1.5) );
  sLatticeNS.setParameter<collision::ITRT::MAXNABLARHO>( maxRhoGradient );
  sLatticeAC.setParameter<descriptors::OMEGA>( 1. / tau_mobil );
  sLatticeAC.addPostProcessor<stage::PhiLimiter>(fluid,meta::id<dispersionLimiter>{});
  sLatticeAC.setParameter<descriptors::EPSILON>( 1.5*w );
  sLatticeAC.setParameter<descriptors::INTERFACE_WIDTH>( w );

  T uMax = helperConvectiveU<T>( myCase );
  sLatticeAC.setParameter<descriptors::MAX_VELOCITY>(uMax);

  {
    auto& communicator = sLatticeAC.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(1);
    communicator.requestField<STATISTIC>();
    communicator.requestField<PSI>();
    communicator.exchangeRequests();
  }
  {
    auto& communicator =
        sLatticeAC.getCommunicator(stage::IterativePostProcess());
    communicator.requestOverlap(1);
    communicator.requestField<PSI>();
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  OstreamManager clout(std::cout, "setInitialValues");
  clout << "Set initial values ..." << std::endl;

  using T = MyCase::value_t;
  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ACDESCRIPTOR = MyCase::descriptor_t_of<Component1>;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& sLatticeNS = myCase.getLattice(NavierStokes{});
  auto& sLatticeAC = myCase.getLattice(Component1{});

  const auto& converter = sLatticeNS.getUnitConverter();

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T inletLength = params.get<parameters::INLET_LENGTH>();
  const T outletLength = params.get<parameters::OUTLET_LENGTH>();
  const int N = params.get<parameters::RESOLUTION>();
  const T dx = extent[1] / N;
  const T rho_g = params.get<parameters::RHO_VAPOR>();
  const T rho_l = params.get<parameters::RHO_LIQUID>();
  const T nu_g = params.get<parameters::NU_VAPOR>();
  const T nu_l = params.get<parameters::NU_LIQUID>();
  const T w = params.get<parameters::INTERFACE_WIDTH>();
  const T theta = params.get<parameters::THETA>();

  AnalyticalConst2D<T, T> rhov(rho_g/converter.getConversionFactorDensity());
  AnalyticalConst2D<T, T> rhol(rho_l/converter.getConversionFactorDensity());
  T tau_g = converter.computeRelaxationTimefromPhysViscosity( nu_g );
  AnalyticalConst2D<T, T> tauv( tau_g );
  T tau_l = converter.computeRelaxationTimefromPhysViscosity( nu_l );
  AnalyticalConst2D<T, T> taul( tau_l );
  AnalyticalConst2D<T, T> angleInside( theta );
  AnalyticalConst2D<T, T> angleOutside( M_PI/2. );
  AnalyticalConst2D<T, T> p0(0.);
  AnalyticalConst2D<T, T> u0(0, 0);

  auto bulk   = geometry.getMaterialIndicator(1);
  auto wall   = geometry.getMaterialIndicator({2, 5});
  auto inlet  = geometry.getMaterialIndicator(3);
  auto outlet = geometry.getMaterialIndicator(4);
  auto fluid  = geometry.getMaterialIndicator({1, 3});
  auto all    = geometry.getMaterialIndicator({0, 1, 2, 3, 4, 5});

  // water enters on the left, gas in the porous rock
  SmoothIndicatorFactoredCuboid2D<T, T> phi0( {-inletLength, extent[1]/2.},
                                              2.2*inletLength, 2.*extent[1],
                                              w*dx/2., 0, {0, 0}, 0, 1. );
  AnalyticalIdentity2D<T, T> rho0(rhov + (rhol - rhov) * phi0);
  AnalyticalIdentity2D<T, T> tau0(tauv + (taul - tauv) * phi0);

  SmoothIndicatorFactoredCuboid2D<T,T> fringe( {-inletLength, extent[1]/2.},
                                                2.*(inletLength+extent[0]+0.5*outletLength), 0,
                                                0.25*outletLength, 0, {0,0}, 0, 1. );

  fields::set<descriptors::SCALAR>(sLatticeNS, all, fringe);
  fields::set<descriptors::RHO>( sLatticeNS, all, rho0);
  fields::set<descriptors::TAU_EFF>( sLatticeNS, all, tau0);
  fields::set<descriptors::OLD_PHIU>( sLatticeAC, all, u0);
  fields::set<descriptors::BOUNDARY>( sLatticeAC, wall, 2);
  fields::set<descriptors::BOUNDARY>( sLatticeAC, geometry.getMaterialIndicator({1, 3, 4}), 1);
  fields::set<descriptors::THETA>( sLatticeAC, geometry.getMaterialIndicator({2}), angleOutside);
  fields::set<descriptors::THETA>( sLatticeAC, geometry.getMaterialIndicator({5}), angleInside);

  momenta::setDensity<T, ACDESCRIPTOR>(sLatticeAC, all, phi0);
  momenta::setVelocity<T, ACDESCRIPTOR>(sLatticeAC, all, u0);
  sLatticeAC.iniEquilibrium(all, phi0, u0);
  momenta::setDensity<T, NSDESCRIPTOR>(sLatticeNS, all, p0);
  momenta::setVelocity<T, NSDESCRIPTOR>(sLatticeNS, all, u0);
  sLatticeNS.iniEquilibrium(all, p0, u0);

  sLatticeAC.executePostProcessors(stage::InitOutlet());
  sLatticeAC.executePostProcessors(stage::PreCoupling());
  signedDistanceFunction<stage::IterativePostProcess>(myCase);

  sLatticeNS.initialize();
  sLatticeAC.initialize();

  sLatticeAC.getCommunicator(stage::PreCoupling()).communicate();

  clout << "Set initial values ... OK" << std::endl;
}

void setTemporalValues(MyCase& myCase, std::size_t iT)
{
  OstreamManager clout(std::cout, "setTemporalValues");

  using T = MyCase::value_t;
  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();
  auto& sLatticeNS = myCase.getLattice(NavierStokes{});
  const auto& converter = sLatticeNS.getUnitConverter();

  const T pressureDrop = params.get<parameters::PRESSURE_DROP>();

  std::size_t iTmaxStart = 10000;
  int  iTupdate   = 1;
  auto inlet      = geometry.getMaterialIndicator(3);

  // Ramp pressure gradient of domain by increasing inlet pressure
  if (iT % iTupdate == 0 && iT <= iTmaxStart) {
    // Smooth start curve, polynomial
    PolynomialStartScale<T, T> StartScale(iTmaxStart, T(1));

    T iTvec[1] = {T(iT)};
    T frac[1]  = {};
    StartScale(frac, iTvec);
    const T imposedPressure = pressureDrop / converter.getConversionFactorPressure() * frac[0];

    momenta::setDensity<T, NSDESCRIPTOR>(sLatticeNS, inlet, imposedPressure);
  }
}

void getResults( MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT )
{
  using T = MyCase::value_t;
  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ACDESCRIPTOR = MyCase::descriptor_t_of<Component1>;
  auto& params = myCase.getParameters();

  auto& sLatticeNS = myCase.getLattice(NavierStokes{});
  auto& sLatticeAC = myCase.getLattice(Component1{});

  const auto& converter = sLatticeNS.getUnitConverter();

  const T statIter = params.get<parameters::PHYS_STAT_ITER_T>();
  const T saveIter = params.get<parameters::PHYS_VTK_ITER_T>();

  SuperVTMwriter2D<T> vtkWriter("gasStorage2d");
  if (iT == 0) {
    SuperLatticeCuboid2D cuboid(sLatticeNS);
    SuperLatticeRank2D rank(sLatticeNS);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.createMasterFile();
  }

  if (iT % converter.getLatticeTime(statIter) == 0) {
    timer.update(iT);
    timer.printStep();
    sLatticeNS.getStatistics().print(iT, converter.getPhysTime(iT));
    sLatticeAC.getStatistics().print(iT, converter.getPhysTime(iT));
  }

  if (iT % converter.getLatticeTime(saveIter) == 0) {
    SuperLatticeDensity2D<T, NSDESCRIPTOR> p_total(sLatticeNS);
    p_total.getName() = "p_total";

    SuperLatticeField2D<T, NSDESCRIPTOR, RHO> rho_L(sLatticeNS);
    AnalyticalConst2D<T,T> ConversionDensity_(converter.getConversionFactorDensity());
    SuperLatticeFfromAnalyticalF2D<T, NSDESCRIPTOR> ConversionDensity(ConversionDensity_, sLatticeNS);
    SuperIdentity2D<T,T> rho(rho_L * ConversionDensity);
    rho.getName() = "rho";

    SuperLatticeVelocity2D<T, NSDESCRIPTOR> velocity(sLatticeNS);
    velocity.getName() = "u";

    SuperLatticeField2D<T, ACDESCRIPTOR, BOUNDARY> boundary(sLatticeAC);
    boundary.getName() = "boundary";

    vtkWriter.addFunctor(p_total);
    vtkWriter.addFunctor(rho);
    vtkWriter.addFunctor(velocity);
    vtkWriter.addFunctor(boundary); // can be used to cut the rock parts in i.e. Paraview
    vtkWriter.write(iT);
  }
}

void simulate(MyCase& myCase)
{
  OstreamManager clout( std::cout,"simulate" );

  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();

  auto& sLatticeNS = myCase.getLattice(NavierStokes{});
  auto& sLatticeAC = myCase.getLattice(Component1{});

  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(
    parameters.get<parameters::PHYS_MAX_T>());

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for ( std::size_t iT=0; iT<=iTmax; ++iT ) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    sLatticeNS.collideAndStream();
    sLatticeAC.collideAndStream();

    sLatticeAC.getCommunicator(stage::PreCoupling()).communicate();
    sLatticeAC.executePostProcessors(stage::PreCoupling());
    if ( iT%100==0 ) {
      signedDistanceFunction<stage::IterativePostProcess>(myCase);
      sLatticeAC.executePostProcessors(stage::PhiLimiter());
    }
    sLatticeAC.getCommunicator(stage::PreCoupling()).communicate();

    myCase.getOperator("Coupling").apply();

    myCase.getOperator("VeloCoupling").apply();

    T uMax = helperConvectiveU<T>( myCase );
    sLatticeAC.setParameter<descriptors::MAX_VELOCITY>(uMax);

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, timer, iT);

    if (std::isnan(sLatticeNS.getStatistics().getAverageEnergy())) {
      clout << "Code Diverged iteration: " << iT << std::endl;
      break;
    }
  }
  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);
  using T = MyCase::value_t;

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION                >(1000);               // resolution inlet [lattice units]
    myCaseParameters.set<VTI_INPUT                 >("gasStorage2d.vti"); // vti file name
    myCaseParameters.set<ARRAY_NAME                >("Tiff Scalars");     // data array name
    myCaseParameters.set<LATTICE_RELAXATION_TIME   >(0.52);    // tau liquid [lattice units]
    myCaseParameters.set<LATTICE_RELAXATION_TIME_PF>(0.8);     // tau mobility [lattice units]
    myCaseParameters.set<INLET_LENGTH_F            >(0.1);     // length of inlet region; % of total length []
    myCaseParameters.set<OUTLET_LENGTH_F           >(0.1);     // length of outlet region; % of total length []
    myCaseParameters.set<PRESSURE_DROP             >(50000.);  // pressure drop [physical units]
    myCaseParameters.set<RHO_LIQUID                >(992.);    // liquid density [physical units]
    myCaseParameters.set<RHO_VAPOR                 >(7.1);     // gas density [physical units]
    myCaseParameters.set<NU_LIQUID                 >(5.5e-7);  // liquid kinematic viscosity [physical units]
    myCaseParameters.set<NU_VAPOR                  >(1.34e-6); // gas kinematic viscosity [physical units]
    myCaseParameters.set<SURFACE_TENSION           >(0.072);   // surface tension [physical units]
    myCaseParameters.set<C_RHO                     >(500.);    // conversion factor density [physical units]
    myCaseParameters.set<parameters::THETA>(M_PI*40./180.); // contact angle [radians]
    myCaseParameters.set<PHYS_MAX_T      >(0.012);   // max simulation time [physical units]
    myCaseParameters.set<PHYS_VTK_ITER_T >(0.012/4000.);   // write simulation output [physical units]
    myCaseParameters.set<PHYS_STAT_ITER_T>(0.012/4000.);   // simulation statistics time [physical units]
    myCaseParameters.set<parameters::INTERFACE_WIDTH>(6.); // resolution inlet [lattice units]

    myCaseParameters.set<SCALE>(0.0000004); // scale down of vti to physical size []
    myCaseParameters.set<DOMAIN_EXTENT>([] (MyCase::ParametersD& p) {
      std::string vtiFile = p.get<VTI_INPUT>();
      std::string arrayName = p.get<ARRAY_NAME>();
      BlockVTIreader2D<T,T> vtiReader( vtiFile, arrayName );

      T deltaRsample = p.get<SCALE>();
      auto cuboidSample = vtiReader.getCuboid();

      Vector<int, 2> extentSample     = cuboidSample.getExtent();
      Vector<T, 2> extentSamplePhys = {deltaRsample * T(extentSample[0]), deltaRsample * T(extentSample[1])};
      return extentSamplePhys;
    });
    myCaseParameters.set<INLET_LENGTH>([] (MyCase::ParametersD& p) {
      Vector extent = p.get<DOMAIN_EXTENT>();
      return extent[0]*p.get<INLET_LENGTH_F>();
    });
    myCaseParameters.set<OUTLET_LENGTH>([] (MyCase::ParametersD& p) {
      Vector extent = p.get<DOMAIN_EXTENT>();
      return extent[0]*p.get<OUTLET_LENGTH_F>();
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
