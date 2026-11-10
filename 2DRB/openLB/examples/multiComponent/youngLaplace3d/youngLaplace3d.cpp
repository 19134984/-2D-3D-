/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2018 Robin Trunk
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

/* youngLaplace3d.cpp
 * In this example a Young-Laplace test is performed. A spherical domain
 * of fluid 2 is immersed in fluid 1. A diffusive interface forms and the
 * surface tension can be calculated using the Laplace pressure relation.
 * The pressure difference is calculated between a point in the middle of
 * the circular domain and a point furthest away from it in the
 * computational domain (here left bottom corner).
 *
 * This example shows the simplest case for the free-energy model with two
 * fluid components.
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::descriptors;
using namespace olb::descriptors;
using namespace olb::graphics;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, D3Q19<CHEM_POTENTIAL,FORCE>>,
  Component1,   Lattice<double, D3Q19<CHEM_POTENTIAL,FORCE>>
>;

// Parameters for the simulation setup
namespace olb::parameters {

struct ALPHA  : public descriptors::FIELD_BASE<1> { };
struct KAPPA1 : public descriptors::FIELD_BASE<1> { };
struct KAPPA2 : public descriptors::FIELD_BASE<1> { };
struct GAMMA  : public descriptors::FIELD_BASE<1> { };
}



Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T nx = extent[0];
  const T physDeltaX = nx / params.get<parameters::RESOLUTION>();

  std::vector<T> extend = { nx, nx, nx };
  std::vector<T> origin = { 0, 0, 0 };
  IndicatorCuboid3D<T> cuboid(extend,origin);

  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());

  // set periodic boundaries to the domain
  mesh.getCuboidDecomposition().setPeriodicity({ true, true, true });
  return mesh;
}

void prepareGeometry( MyCase& myCase )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  auto& geometry = myCase.getGeometry();

  geometry.rename( 0,1 );

  geometry.innerClean();
  geometry.checkForErrors();
  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( MyCase& myCase )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;

  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& sLattice1 = myCase.getLattice(NavierStokes{});
  auto& sLattice2 = myCase.getLattice(Component1{});

  const int N = params.get<parameters::RESOLUTION>();
  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T nx = extent[0];
  const T gama = params.get<parameters::GAMMA>();
  const T alpha = params.get<parameters::ALPHA>();
  const T kappa1 = params.get<parameters::KAPPA1>();
  const T kappa2 = params.get<parameters::KAPPA2>();

  sLattice1.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    int   {N}, // resolution
    (T)   1., // lattice relaxation time
    (T)   nx, // charPhysLength: reference length of simulation geometry
    (T)   1.e-6, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   0.1, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1. // physDensity: physical density in __kg / m^3__
  );

  const auto& converter = sLattice1.getUnitConverter();
  converter.print();

  sLattice2.setUnitConverter(converter);
  dynamics::set<ForcedBGKdynamics>(sLattice1, geometry.getMaterialIndicator({1}));
  dynamics::set<FreeEnergyBGKdynamics>(sLattice2, geometry.getMaterialIndicator({1}));

  sLattice1.setParameter<descriptors::OMEGA>( sLattice1.getUnitConverter().getLatticeRelaxationFrequency() );
  sLattice2.setParameter<descriptors::OMEGA>( sLattice2.getUnitConverter().getLatticeRelaxationFrequency() );
  sLattice2.setParameter<collision::FreeEnergy::GAMMA>(gama);

  auto& coupling1 = myCase.setCouplingOperator(
  "Chemical_potential",
  ChemicalPotentialCoupling3D{},
  names::A{}, sLattice1,
  names::B{}, sLattice2);

  coupling1.template setParameter<ChemicalPotentialCoupling3D::ALPHA>(alpha);
  coupling1.template setParameter<ChemicalPotentialCoupling3D::KAPPA1>(kappa1);
  coupling1.template setParameter<ChemicalPotentialCoupling3D::KAPPA2>(kappa2);

  myCase.setCouplingOperator(
  "Force",
  ForceCoupling3D{},
  names::A{}, sLattice2,
  names::B{}, sLattice1);

  {
    auto& communicator = sLattice1.getCommunicator(stage::PostCoupling());
    communicator.requestField<CHEM_POTENTIAL>();
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

  sLattice1.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());
  sLattice2.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());

  sLattice1.initialize();
  sLattice2.initialize();


  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues( MyCase& myCase )
{
  OstreamManager clout( std::cout,"initialValues" );

  using T = MyCase::value_t;

  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();
  auto& sLattice1 = myCase.getLattice(NavierStokes{});
  auto& sLattice2 = myCase.getLattice(Component1{});

  const Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T nx = extent[0];
  const T radius = params.get<parameters::RADIUS>();
  const T alpha = params.get<parameters::ALPHA>();

  // bulk initial conditions
  // define circular domain for fluid 2
  std::vector<T> v( 2,T() );
  AnalyticalConst3D<T,T> zeroVelocity( v );

  AnalyticalConst3D<T,T> one ( 1. );
  IndicatorSphere3D<T> sphere( {nx/T(2), nx/T(2), nx/T(2)}, radius );
  SmoothIndicatorSphere3D<T,T> smoothSphere( sphere, 10.*alpha );

  AnalyticalIdentity3D<T,T> rho( one );
  AnalyticalIdentity3D<T,T> phi( one - smoothSphere - smoothSphere );

  sLattice1.iniEquilibrium( geometry, 1, rho, zeroVelocity );
  sLattice2.iniEquilibrium( geometry, 1, phi, zeroVelocity );


  sLattice1.initialize();
  sLattice2.initialize();

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
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{ }

void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  OstreamManager clout( std::cout,"getResults" );

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;

  auto& sLattice1 = myCase.getLattice(NavierStokes{});
  auto& sLattice2 = myCase.getLattice(Component1{});
  auto& params = myCase.getParameters();


  const int N = params.get<parameters::RESOLUTION>();
  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T nx = extent[0];
  const T radius = params.get<parameters::RADIUS>();
  const T alpha = params.get<parameters::ALPHA>();
  const T kappa1 = params.get<parameters::KAPPA1>();
  const T kappa2 = params.get<parameters::KAPPA2>();

  SuperVTMwriter3D<T> vtmWriter( "youngLaplace3d" );

  const int vtkIter  = 200;
  const int statIter = 200;

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice1 );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice1 );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Get statistics
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();
    sLattice1.getStatistics().print( iT, sLattice1.getUnitConverter().getPhysTime(iT) );
    sLattice2.getStatistics().print( iT, sLattice2.getUnitConverter().getPhysTime(iT) );
  }

  // Writes the VTK files
  if ( iT%vtkIter==0 ) {
    sLattice1.setProcessingContext(ProcessingContext::Evaluation);
    AnalyticalConst3D<T,T> half_( 0.5 );
    SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> half(half_, sLattice1);

    SuperLatticeDensity3D<T, DESCRIPTOR> density1( sLattice1 );
    density1.getName() = "rho";
    SuperLatticeDensity3D<T, DESCRIPTOR> density2( sLattice2 );
    density2.getName() = "phi";

    SuperIdentity3D<T,T> c1 (half*(density1+density2));
    c1.getName() = "density-fluid-1";
    SuperIdentity3D<T,T> c2 (half*(density1-density2));
    c2.getName() = "density-fluid-2";

    vtmWriter.addFunctor( density1 );
    vtmWriter.addFunctor( density2 );
    vtmWriter.addFunctor( c1 );
    vtmWriter.addFunctor( c2 );
    vtmWriter.write( iT );

    // calculate bulk pressure, pressure difference and surface tension
    if (iT%statIter==0) {
      AnalyticalConst3D<T,T> two_( 2. );
      AnalyticalConst3D<T,T> onefive_( 1.5 );
      AnalyticalConst3D<T,T> k1_( kappa1 );
      AnalyticalConst3D<T,T> k2_( kappa2 );
      AnalyticalConst3D<T,T> cs2_( 1./descriptors::invCs2<T,DESCRIPTOR>() );
      SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> two(two_, sLattice1);
      SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> onefive(onefive_, sLattice1);
      SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> k1(k1_, sLattice1);
      SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> k2(k2_, sLattice1);
      SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> cs2(cs2_, sLattice1);

      // Calculation of bulk pressure:
      SuperIdentity3D<T,T> bulkPressure ( density1*cs2
                                          + k1*( onefive*c1*c1*c1*c1 - two*c1*c1*c1 + half*c1*c1 )
                                          + k2*( onefive*c2*c2*c2*c2 - two*c2*c2*c2 + half*c2*c2 ) );

      AnalyticalFfromSuperF3D<T, T> interpolPressure( bulkPressure, true, 1);
      T position[3] = {nx/T(2), nx/T(2), nx/T(2)};
      T pressureIn = 0.;
      T pressureOut = 0.;
      interpolPressure(&pressureIn, position);
      position[0] = (T(N)/T(100))*sLattice1.getUnitConverter().getPhysDeltaX();
      position[1] = (T(N)/T(100))*sLattice1.getUnitConverter().getPhysDeltaX();
      position[2] = (T(N)/T(100))*sLattice1.getUnitConverter().getPhysDeltaX();
      interpolPressure(&pressureOut, position);

      clout << "Pressure Difference: " << pressureIn-pressureOut << "  ;  ";
      clout << "Surface Tension: " << radius*(pressureIn-pressureOut)/2 << std::endl;
      clout << "Analytical Pressure Difference: " << alpha/(3.*radius) * (kappa1 + kappa2) << "  ;  ";
      clout << "Analytical Surface Tension: " << alpha/6. * (kappa1 + kappa2) << std::endl;
    }
  }
}

void simulate(MyCase& myCase) {
  OstreamManager clout( std::cout,"simulate" );

  using T = MyCase::value_t;

  auto& geometry = myCase.getGeometry();
  auto& sLattice1 = myCase.getLattice(NavierStokes{});
  auto& sLattice2 = myCase.getLattice(Component1{});
  auto& coupling1 = myCase.getOperator("Chemical_potential");
  auto& coupling2 = myCase.getOperator("Force");

  const int maxIter  = 60000;

  int iT = 0;
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( maxIter, geometry.getStatistics().getNvoxel() );
  timer.start();


  for ( iT=0; iT<=maxIter; ++iT ) {
    // Computation and output of the results
    getResults(myCase, timer, iT);

    // Collide and stream execution
    sLattice1.collideAndStream();
    sLattice2.collideAndStream();

    // Execute coupling between the two lattices
    sLattice1.executePostProcessors(stage::PreCoupling());

    sLattice1.getCommunicator(stage::PreCoupling()).communicate();
    coupling1.apply();
    sLattice1.getCommunicator(stage::PostCoupling()).communicate();

    sLattice2.executePostProcessors(stage::PreCoupling());

    sLattice2.getCommunicator(stage::PreCoupling()).communicate();
    coupling2.apply();
    sLattice2.getCommunicator(stage::PostCoupling()).communicate();
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
    myCaseParameters.set<DOMAIN_EXTENT    >(100);         // [lattice units]
    myCaseParameters.set<RESOLUTION       >(100);         // [lattice units]
    myCaseParameters.set<ALPHA            >(1.5);         // Interfacial width          [lattice units]
    myCaseParameters.set<KAPPA1           >(0.0075);      // For surface tensions       [lattice units]
    myCaseParameters.set<KAPPA2           >(0.005);       // For surface tensions       [lattice units]
    myCaseParameters.set<parameters::GAMMA>(1.);          // For mobility of interfaces [lattice units]
    myCaseParameters.set<parameters::RADIUS>([] (MyCase::ParametersD& p) {
      return p.get<DOMAIN_EXTENT>()[0] * 0.25;
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
