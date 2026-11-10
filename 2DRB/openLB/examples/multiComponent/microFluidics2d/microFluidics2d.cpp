/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2019 Sam Avis
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

/*  microfluidics2d.cpp:
 *  This example shows a microfluidic channel creating droplets of
 *  two fluid components. Poiseuille velocity profiles are imposed
 *  at the various channel inlets, while a constant density outlet
 *  is imposed at the end of the channel to allow the droplets to
 *  exit the simulation.
 *
 *  This example demonstrates the use of three fluid components
 *  with the free energy model. It also shows the use of open
 *  boundary conditions, specifically velocity inlet and density
 *  outlet boundaries.
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::descriptors;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, D2Q9<CHEM_POTENTIAL,FORCE>>,
  Component1,   Lattice<double, D2Q9<CHEM_POTENTIAL,FORCE>>,
  Component2,   Lattice<double, D2Q9<CHEM_POTENTIAL,FORCE>>
>;

namespace olb::parameters {

struct TAU    : public descriptors::FIELD_BASE<1> { };

struct ALPHA  : public descriptors::FIELD_BASE<1> { };
struct KAPPA1 : public descriptors::FIELD_BASE<1> { };
struct KAPPA2 : public descriptors::FIELD_BASE<1> { };
struct KAPPA3 : public descriptors::FIELD_BASE<1> { };
struct GAMMA  : public descriptors::FIELD_BASE<1> { };
struct H1     : public descriptors::FIELD_BASE<1> { };
struct H2     : public descriptors::FIELD_BASE<1> { };
struct H3     : public descriptors::FIELD_BASE<1> { };

struct INLET_VELCOCITY1 : public descriptors::FIELD_BASE<1> { };
struct INLET_VELCOCITY2 : public descriptors::FIELD_BASE<1> { };
struct INLET_VELCOCITY3 : public descriptors::FIELD_BASE<1> { };

struct XL1 : public descriptors::FIELD_BASE<1> { };
struct XL2 : public descriptors::FIELD_BASE<1> { };
struct XL3 : public descriptors::FIELD_BASE<1> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extent, origin);

  const T physDeltaX = extent[1] / params.get<parameters::RESOLUTION>();
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  return mesh;
}

double helperFunction( double alpha, double kappa1, double kappa2,
                       double h1, double h2, int latticeNumber )
{
  double addend = 0;
  if (latticeNumber==1) {
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (h2/kappa2) );
  }
  else if (latticeNumber==2) {
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (-h2/kappa2) );
  }
  else if (latticeNumber==3) {
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (h2/kappa2) );
  }
  return addend;
}

double helperFunction( double alpha, double kappa1, double kappa2, double kappa3,
                       double h1, double h2, double h3, int latticeNumber )
{
  double addend = 0;
  if (latticeNumber==1) {
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (h2/kappa2) + (h3/kappa3) );
  }
  else if (latticeNumber==2) {
    addend = 1./(alpha*alpha) * ( (h1/kappa1) + (-h2/kappa2) );
  }
  else if (latticeNumber==3) {
    addend = 1./(alpha*alpha) * ( (h3/kappa3) );
  }
  return addend;
}

void prepareGeometry( MyCase& myCase )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  const T nx = params.get<parameters::DOMAIN_EXTENT>()[0];
  const T ny = params.get<parameters::DOMAIN_EXTENT>()[1];
  const T xl1 = params.get<parameters::XL1>();
  const T xl2 = params.get<parameters::XL2>();
  const T xl3 = params.get<parameters::XL3>();
  const T dx = ny / params.get<parameters::RESOLUTION>();

  const T yl1 = ny / 4.;
  const T yl2 = ny;
  const T yl3 = ny / 4.;
  const T xl4 = 175. / 7.;
  const T yl4 = ny;
  const T xl5 = nx - 175.;
  const T yl5 = ny / 2.;

  std::shared_ptr<IndicatorF2D<T>> section1 = std::make_shared<IndicatorCuboid2D<T>>( xl1, yl1, std::vector<T> {xl1/T(2), ny/T(2)} );
  std::shared_ptr<IndicatorF2D<T>> section2 = std::make_shared<IndicatorCuboid2D<T>>( xl2, yl2, std::vector<T> {xl1 + xl2/T(2), ny/T(2)} );
  std::shared_ptr<IndicatorF2D<T>> section3 = std::make_shared<IndicatorCuboid2D<T>>( xl3, yl3, std::vector<T> {xl1 + xl2 + xl3/T(2), ny/T(2)} );
  std::shared_ptr<IndicatorF2D<T>> section4 = std::make_shared<IndicatorCuboid2D<T>>( xl4, yl4, std::vector<T> {xl1 + xl2 + xl3 + xl4/T(2), ny/T(2)} );
  std::shared_ptr<IndicatorF2D<T>> section5 = std::make_shared<IndicatorCuboid2D<T>>( xl5, yl5, std::vector<T> {xl1 + xl2 + xl3 + xl4 + xl5/T(2), ny/T(2)} );
  IndicatorIdentity2D<T> channel( section1 + section2 + section3 + section4 + section5 );

  geometry.rename( 0, 2, channel );
  geometry.rename( 2,1,{1,1} );

  // Inlets and outlet
  IndicatorCuboid2D<T> inlet1 ( dx, yl1, {0., ny/T(2)} );
  IndicatorCuboid2D<T> inlet21( xl2 - dx, dx, {xl1 + xl2/T(2), 0.} );
  IndicatorCuboid2D<T> inlet22( xl2 - dx, dx, {xl1 + xl2/T(2), ny} );
  IndicatorCuboid2D<T> inlet31( xl4 - dx, dx, {xl1 + xl2 + xl3 + xl4/T(2), 0.} );
  IndicatorCuboid2D<T> inlet32( xl4 - dx, dx, {xl1 + xl2 + xl3 + xl4/T(2), ny} );
  IndicatorCuboid2D<T> outlet( dx, yl5, {nx, ny/T(2)} );
  geometry.rename( 2, 3, 1, inlet1 );
  geometry.rename( 2, 4, 1, inlet21 );
  geometry.rename( 2, 5, 1, inlet22 );
  geometry.rename( 2, 6, 1, inlet31 );
  geometry.rename( 2, 7, 1, inlet32 );
  geometry.rename( 2, 8, 1, outlet );

  geometry.innerClean();
  geometry.checkForErrors();
  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( MyCase& myCase )
{
  OstreamManager clout( std::cout,"prepareLattice" );

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& sLattice1 = myCase.getLattice(NavierStokes{});
  auto& sLattice2 = myCase.getLattice(Component1{});
  auto& sLattice3 = myCase.getLattice(Component2{});

  using DESCRIPTOR = MyCase::descriptor_t;

  // TODO for now, to be combined with unit converter refactor
  const T outletDensity = 1.;

  const int N = params.get<parameters::RESOLUTION>();
  const T ny = params.get<parameters::DOMAIN_EXTENT>()[1];

  const T tau = params.get<parameters::TAU>();
  const T alpha = params.get<parameters::ALPHA>();
  const T kappa1 = params.get<parameters::KAPPA1>();
  const T kappa2 = params.get<parameters::KAPPA2>();
  const T kappa3 = params.get<parameters::KAPPA3>();
  const T gama = params.get<parameters::GAMMA>();
  const T h1 = params.get<parameters::H1>();
  const T h2 = params.get<parameters::H2>();
  const T h3 = params.get<parameters::H3>();

  sLattice1.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    int   {N}, // resolution
    (T)   tau, // lattice relaxation time
    (T)   ny, // charPhysLength: reference length of simulation geometry
    (T)   1.e-6, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   0.1, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1. // physDensity: physical density in __kg / m^3__
  );

  const auto& converter = sLattice1.getUnitConverter();
  converter.print();

  sLattice2.setUnitConverter(converter);
  sLattice3.setUnitConverter(converter);

  // define lattice dynamics
  clout << "Prepare Lattice: Define lattice dynamics ..." << std::endl;
  dynamics::set<ForcedBGKdynamics>( sLattice1, geometry.getMaterialIndicator({1}));
  dynamics::set<FreeEnergyBGKdynamics>( sLattice2, geometry.getMaterialIndicator({1}));
  dynamics::set<FreeEnergyBGKdynamics>( sLattice3, geometry.getMaterialIndicator({1}));

  // Defining walls
  auto walls = geometry.getMaterialIndicator({2});

  // Compute Addends
  T addend1 = helperFunction( alpha, kappa1, kappa2, kappa3, h1, h2, h3, 1 );
  T addend2 = helperFunction( alpha, kappa1, kappa2, kappa3, h1, h2, h3, 2 );
  T addend3 = helperFunction( alpha, kappa1, kappa2, kappa3, h1, h2, h3, 3 );

  // Add wall boundary
  clout << "Prepare Lattice: Add wall boundary ..." << std::endl;
  boundary::set<boundary::FreeEnergyWallMomentum>(sLattice1, walls);
  sLattice1.setParameter<descriptors::ADDEND>( addend1 );
  boundary::set<boundary::FreeEnergyWallOrderParameter>(sLattice2, walls);
  sLattice2.setParameter<descriptors::ADDEND>( addend2 );
  boundary::set<boundary::FreeEnergyWallOrderParameter>(sLattice3, walls);
  sLattice3.setParameter<descriptors::ADDEND>( addend3 );

  // add inlet boundaries
  clout << "Prepare Lattice: Add inlet boundaries ..." << std::endl;
  T omega = converter.getLatticeRelaxationFrequency();
  auto inlet1Indicator = geometry.getMaterialIndicator(3);
  boundary::set<boundary::FreeEnergyVelocity>(sLattice1, inlet1Indicator);
  boundary::set<boundary::FreeEnergyOrderParameter>(sLattice2, inlet1Indicator);
  boundary::set<boundary::FreeEnergyOrderParameter>(sLattice3, inlet1Indicator);

  auto inlet2Indicator = geometry.getMaterialIndicator({4, 5});
  boundary::set<boundary::FreeEnergyVelocity>(sLattice1, inlet2Indicator);
  boundary::set<boundary::FreeEnergyOrderParameter>(sLattice2, inlet2Indicator);
  boundary::set<boundary::FreeEnergyOrderParameter>(sLattice3, inlet2Indicator);

  auto inlet3Indicator = geometry.getMaterialIndicator({6, 7});
  boundary::set<boundary::FreeEnergyVelocity>(sLattice1, inlet3Indicator);
  boundary::set<boundary::FreeEnergyOrderParameter>(sLattice2, inlet3Indicator);
  boundary::set<boundary::FreeEnergyOrderParameter>(sLattice3, inlet3Indicator);

  // add outlet boundary
  clout << "Prepare Lattice: Add outlet boundary ..." << std::endl;
  auto outletIndicator = geometry.getMaterialIndicator(8);
  boundary::set<boundary::FreeEnergyPressureConvective>(sLattice1, outletIndicator);
  boundary::set<boundary::FreeEnergyOrderParameterConvective>(sLattice2, outletIndicator);
  boundary::set<boundary::FreeEnergyOrderParameterConvective>(sLattice3, outletIndicator);

  sLattice1.setParameter<descriptors::OMEGA>(omega);
  sLattice2.setParameter<descriptors::OMEGA>(omega);
  sLattice2.setParameter<collision::FreeEnergy::GAMMA>(gama);
  sLattice3.setParameter<descriptors::OMEGA>(omega);
  sLattice3.setParameter<collision::FreeEnergy::GAMMA>(gama);

  {
    auto& communicator = sLattice1.getCommunicator(stage::PostPostProcess());
    communicator.requestField<POPULATION>();
    communicator.requestOverlap(sLattice1.getOverlap());
    communicator.exchangeRequests();
  }
  {
    auto& communicator = sLattice2.getCommunicator(stage::PostPostProcess());
    communicator.requestField<POPULATION>();
    communicator.requestOverlap(sLattice2.getOverlap());
    communicator.exchangeRequests();
  }
  {
    auto& communicator = sLattice3.getCommunicator(stage::PostPostProcess());
    communicator.requestField<POPULATION>();
    communicator.requestOverlap(sLattice3.getOverlap());
    communicator.exchangeRequests();
  }

  auto& coupling1 = myCase.setCouplingOperator(
    "Outlet_density",
    DensityOutletCoupling2D{},
    names::A{}, sLattice1,
    names::B{}, sLattice2,
    names::C{}, sLattice3
  );
  coupling1.setParameter<DensityOutletCoupling2D::RHO>(outletDensity);
  coupling1.restrictTo(geometry.getMaterialIndicator({8}));

  auto& coupling2 = myCase.setCouplingOperator(
    "Chemical_potential",
    ChemicalPotentialCoupling2D{},
    names::A{}, sLattice1,
    names::B{}, sLattice2,
    names::C{}, sLattice3
  );
  coupling2.setParameter<ChemicalPotentialCoupling2D::ALPHA>(alpha);
  coupling2.setParameter<ChemicalPotentialCoupling2D::KAPPA1>(kappa1);
  coupling2.setParameter<ChemicalPotentialCoupling2D::KAPPA2>(kappa2);
  coupling2.setParameter<ChemicalPotentialCoupling2D::KAPPA2>(kappa3);
  coupling2.restrictTo(geometry.getMaterialIndicator({1}));

  auto& coupling3 = myCase.setCouplingOperator(
    "Force",
    ForceCoupling2D{},
    names::A{}, sLattice2,
    names::B{}, sLattice1,
    names::C{}, sLattice3
  );
  coupling3.restrictTo(geometry.getMaterialIndicator({1}));

  auto& coupling4 = myCase.setCouplingOperator(
    "Inlet_outlet",
    InletOutletCoupling2D{},
    names::A{}, sLattice2,
    names::B{}, sLattice1,
    names::C{}, sLattice3
  );
  coupling4.restrictTo(geometry.getMaterialIndicator({3,4,5,6,7,8}));

  sLattice1.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());
  sLattice2.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());
  sLattice3.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());

  {
    auto& communicator = sLattice1.getCommunicator(stage::PostCoupling());
    communicator.requestField<CHEM_POTENTIAL>();
    communicator.requestField<RhoStatistics>();
    communicator.requestOverlap(sLattice1.getOverlap());
    communicator.exchangeRequests();
  }
  {
    auto& communicator = sLattice2.getCommunicator(stage::PreCoupling());
    communicator.requestField<CHEM_POTENTIAL>();
    communicator.requestField<RhoStatistics>();
    communicator.requestOverlap(sLattice2.getOverlap());
    communicator.exchangeRequests();
  }
  {
    auto& communicator = sLattice3.getCommunicator(stage::PreCoupling());
    communicator.requestField<CHEM_POTENTIAL>();
    communicator.requestField<RhoStatistics>();
    communicator.requestOverlap(sLattice3.getOverlap());
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  OstreamManager clout( std::cout,"initialValues" );
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& sLattice1 = myCase.getLattice(NavierStokes{});
  auto& sLattice2 = myCase.getLattice(Component1{});
  auto& sLattice3 = myCase.getLattice(Component2{});

  const T outletDensity = 1.;

  const int N = params.get<parameters::RESOLUTION>();
  const T ny = params.get<parameters::DOMAIN_EXTENT>()[1];
  const T dx = ny / N;
  const T xl1 = params.get<parameters::XL1>();
  const T xl2 = params.get<parameters::XL2>();
  const T xl3 = params.get<parameters::XL3>();

  const T inlet1Velocity = params.get<parameters::INLET_VELCOCITY1>();
  const T inlet2Velocity = params.get<parameters::INLET_VELCOCITY2>();
  const T inlet3Velocity = params.get<parameters::INLET_VELCOCITY3>();

  // bulk initial conditions
  clout << "Prepare Lattice: Bulk initial conditions ..." << std::endl;
  std::vector<T> v( 2,T() );
  AnalyticalConst2D<T,T> zeroVelocity( v );

  AnalyticalConst2D<T,T> zero ( 0. );
  AnalyticalConst2D<T,T> one ( 1. );
  IndicatorCuboid2D<T> ind1(xl1+dx, ny, {xl1/T(2), ny/T(2)});
  SmoothIndicatorCuboid2D<T,T> section1( ind1, 0. );
  IndicatorCuboid2D<T> ind2(xl2 + xl3, ny, {xl1 + (xl2 + xl3)/T(2), ny/T(2)});
  SmoothIndicatorCuboid2D<T,T> section2( ind2, 0. );

  AnalyticalIdentity2D<T,T> c1( section1 );
  AnalyticalIdentity2D<T,T> c2( section2 );
  AnalyticalIdentity2D<T,T> rho( one );
  AnalyticalIdentity2D<T,T> phi( c1 - c2 );
  AnalyticalIdentity2D<T,T> psi( rho - c1 - c2 );

  auto allIndicator = geometry.getMaterialIndicator({1, 2, 3, 4, 5, 6});
  auto inlet1Indicator = geometry.getMaterialIndicator(3);
  auto inlet2Indicator = geometry.getMaterialIndicator({4, 5});
  auto inlet3Indicator = geometry.getMaterialIndicator({6, 7});
  auto outletIndicator = geometry.getMaterialIndicator(8);

  sLattice1.iniEquilibrium( allIndicator, rho, zeroVelocity );
  sLattice2.iniEquilibrium( allIndicator, phi, zeroVelocity );
  sLattice3.iniEquilibrium( allIndicator, psi, zeroVelocity );

  // inlet boundary conditions
  clout << "Inlet boundary conditions ..." << std::endl;
  Poiseuille2D<T> inlet1U( geometry, 3, 1.5*inlet1Velocity, 0. );
  sLattice1.defineU( inlet1Indicator, inlet1U );
  sLattice2.defineRho( inlet1Indicator, phi );
  sLattice3.defineRho( inlet1Indicator, psi );

  Poiseuille2D<T> inlet21U( geometry, 4, 1.5*inlet2Velocity, 0. );
  Poiseuille2D<T> inlet22U( geometry, 5, 1.5*inlet2Velocity, 0. );
  sLattice1.defineU( geometry, 4, inlet21U );
  sLattice1.defineU( geometry, 5, inlet22U );
  sLattice2.defineRho( inlet2Indicator, phi );
  sLattice3.defineRho( inlet2Indicator, psi );

  Poiseuille2D<T> inlet31U( geometry, 6, 1.5*inlet3Velocity, 0. );
  Poiseuille2D<T> inlet32U( geometry, 7, 1.5*inlet3Velocity, 0. );
  sLattice1.defineU( geometry, 6, inlet31U );
  sLattice1.defineU( geometry, 7, inlet32U );
  sLattice2.defineRho( inlet3Indicator, phi );
  sLattice3.defineRho( inlet3Indicator, psi );

  // outlet initial / boundary conditions
  clout << "Outlet initial / Boundary conditions ..." << std::endl;
  AnalyticalConst2D<T,T> rhoOutlet( outletDensity );
  AnalyticalIdentity2D<T,T> phiOutlet( zero );
  AnalyticalIdentity2D<T,T> psiOutlet( rhoOutlet );
  sLattice1.defineRho( outletIndicator, rhoOutlet );
  sLattice2.defineRho( outletIndicator, phiOutlet );
  sLattice3.defineRho( outletIndicator, psiOutlet );

  clout << "Initialize lattices ..." << std::endl;
  sLattice1.initialize();
  sLattice2.initialize();
  sLattice3.initialize();

  sLattice1.communicate();
  sLattice2.communicate();
  sLattice3.communicate();
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{ }

void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& sLattice1 = myCase.getLattice(NavierStokes{});
  auto& sLattice2 = myCase.getLattice(Component1{});
  auto& sLattice3 = myCase.getLattice(Component2{});
  const auto& converter = sLattice1.getUnitConverter();

  SuperVTMwriter2D<T> vtkWriter( "microFluidics2d" );

  // TODO such times should also be parameters
  const int statIter = 2000;
  const int saveIter  = 1000;

  if ( iT==0 ) {
    SuperLatticeCuboid2D cuboid( sLattice1 );
    SuperLatticeRank2D rank( sLattice1 );
    vtkWriter.write( cuboid );
    vtkWriter.write( rank );
    vtkWriter.createMasterFile();
  }

  if ( iT%statIter==0 ) {
    timer.update( iT );
    timer.printStep();
    sLattice1.getStatistics().print( iT, converter.getPhysTime(iT) );
    sLattice2.getStatistics().print( iT, converter.getPhysTime(iT) );
    sLattice3.getStatistics().print( iT, converter.getPhysTime(iT) );
  }

  // Writes the VTK files
  if ( iT%saveIter==0 ) {
    sLattice1.setProcessingContext(ProcessingContext::Evaluation);
    sLattice2.setProcessingContext(ProcessingContext::Evaluation);
    sLattice3.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticeVelocity2D velocity( sLattice1 );
    SuperLatticeDensity2D density1( sLattice1 );
    density1.getName() = "rho";
    SuperLatticeDensity2D density2( sLattice2 );
    density2.getName() = "phi";
    SuperLatticeDensity2D density3( sLattice3 );
    density3.getName() = "density-fluid-3";

    AnalyticalConst2D<T,T> half_( 0.5 );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> half(half_, sLattice1);

    SuperIdentity2D<T,T> c1 (half*(density1+density2-density3));
    c1.getName() = "density-fluid-1";
    SuperIdentity2D<T,T> c2 (half*(density1-density2-density3));
    c2.getName() = "density-fluid-2";

    vtkWriter.addFunctor( velocity );
    vtkWriter.addFunctor( density1 );
    vtkWriter.addFunctor( density2 );
    vtkWriter.addFunctor( density3 );
    vtkWriter.addFunctor( c1 );
    vtkWriter.addFunctor( c2 );
    vtkWriter.write( iT );
  }
}

void simulate(MyCase& myCase) {
  using T = MyCase::value_t;

  const std::size_t iTmax = 1000000;

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for ( std::size_t iT=0; iT<iTmax; ++iT ) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, timer, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    myCase.getLattice(NavierStokes{}).collideAndStream();
    myCase.getLattice(Component1{}).collideAndStream();
    myCase.getLattice(Component2{}).collideAndStream();

    myCase.getLattice(NavierStokes{}).executePostProcessors(stage::PreCoupling());
    myCase.getLattice(Component1{}).executePostProcessors(stage::PreCoupling());
    myCase.getLattice(Component2{}).executePostProcessors(stage::PreCoupling());

    myCase.getLattice(NavierStokes{}).getCommunicator(stage::PreCoupling()).communicate();
    myCase.getLattice(Component1{}).getCommunicator(stage::PreCoupling()).communicate();
    myCase.getLattice(Component2{}).getCommunicator(stage::PreCoupling()).communicate();

    // Execute coupling between the two lattices
    myCase.getOperator("Outlet_density").apply();
    myCase.getOperator("Chemical_potential").apply();

    myCase.getLattice(NavierStokes{}).getCommunicator(stage::PostCoupling()).communicate();
    myCase.getLattice(NavierStokes{}).executePostProcessors(stage::PostCoupling());

    myCase.getOperator("Force").apply();
    myCase.getOperator("Inlet_outlet").apply();
  }

  timer.stop();
  timer.printSummary();
}

int main( int argc, char *argv[] ) {
  initialize( &argc, &argv );

  // === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<DOMAIN_EXTENT   >({800., 100.});// [lattice units]
    myCaseParameters.set<RESOLUTION      >(50);          // [lattice units]
    myCaseParameters.set<XL1             >(175.*2./7.);  // [lattice units]
    myCaseParameters.set<XL2             >(175./7.);     // [lattice units]
    myCaseParameters.set<XL3             >(175.*3./7.);  // [lattice units]
    myCaseParameters.set<parameters::TAU >(1.);          // [lattice units]
    myCaseParameters.set<INLET_VELCOCITY1>(0.00056);     // [lattice units]
    myCaseParameters.set<INLET_VELCOCITY2>(0.00055);     // [lattice units]
    myCaseParameters.set<INLET_VELCOCITY3>(0.0014);      // [lattice units]
    myCaseParameters.set<ALPHA           >(1.);          // Interfacial width          [lattice units]
    myCaseParameters.set<KAPPA1          >(0.0132);      // For surface tensions       [lattice units]
    myCaseParameters.set<KAPPA2          >(0.0012);      // For surface tensions       [lattice units]
    myCaseParameters.set<KAPPA3          >(0.0013);      // For surface tensions       [lattice units]
    myCaseParameters.set<parameters::GAMMA>(1.);         // For mobility of interfaces [lattice units]
    myCaseParameters.set<H1              >(0.);          // Contact angle 90 degrees   [lattice units]
    myCaseParameters.set<H2              >(0.);          // Contact angle 90 degrees   [lattice units]
    myCaseParameters.set<H3              >(0.);          // Contact angle 90 degrees   [lattice units]
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
