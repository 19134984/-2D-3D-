/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2020 Johannes Nguyen, Stephan Simonis, Sam Avis
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

 /* binaryShearFlow2d.cpp
  * A circular domain of fluid phase II is immersed in a square
  * filled with fluid phase I. The top and bottom walls are
  * moving in opposite directions, such that the droplet consisting
  * of phase II is exposed to shear flow and deforms accordingly.
  * The default parameter setting is given in
  *  A.E. Komrakova et al. International Journal of Multiphase
  *  Flow 59 (2014) 24–43.
  * The present example uses the ternary free energy model by
  *  C Semprebon et al. Physical Review E 93.3 (2016) p. 033305.
  * Reference results are published in
  *  S. Simonis et al. Discrete and Continuous Dynamical
  *  Systems - Series S (2023), doi:10.3934/dcdss.2023069
  *
  * The droplet breakup case should be run on parallel mode
  * because of increased computation time.
  */

#include <olb.h>

#include <chrono>

using namespace olb;
using namespace olb::names;
using namespace olb::descriptors;
using namespace olb::graphics;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, D2Q9<CHEM_POTENTIAL,FORCE>>,
  Component1,   Lattice<double, D2Q9<CHEM_POTENTIAL,FORCE>>
>;

namespace olb::parameters {

struct NORM_VTK_ITER_T       : public descriptors::FIELD_BASE<1> { };
struct NORM_STAT_ITER_T       : public descriptors::FIELD_BASE<1> { };
struct MAX_NORM_T : public descriptors::FIELD_BASE<1> { };
struct CONT_VISC  : public descriptors::FIELD_BASE<1> { };
struct DROPLET_REYNOLDS : public descriptors::FIELD_BASE<1> { };
struct PHYS_INTERVAL : public descriptors::FIELD_BASE<1> { };
struct CAPILLARY : public descriptors::FIELD_BASE<1> { };
struct CAHN  : public descriptors::FIELD_BASE<1> { };
struct H1  : public descriptors::FIELD_BASE<1> { };
struct H2  : public descriptors::FIELD_BASE<1> { };
}

// // DEFAULT [lu] Setup Komrakova et al. 2013 "Lattice Boltzmann ... " Fig. 9
// // Simulated Domain LxHxW = 12*radius x 8*radius x 2*radius;
// // Resolution: radius = 30 [lu]
// // Re = 0.0625
// // Ca = 0.42 (breakup for Stokes flow)
// // Ch = 0.0379
// // Pe = 0.43
// // Breakup happens at normTime ~30

// Guideline parameters:
// Ca = 0.25 for steady state
// Ca = 1.0  for drop breakup
// Recommended xi = 1.14, min 1.0


Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T nx = extent[0];
  const T ny = extent[1];
  const T N = params.get<parameters::RESOLUTION>();
  const T radius = params.get<parameters::RADIUS>();
  std::vector<T> extend = { nx, ny };
  std::vector<T> origin = { 0, 0 };
  IndicatorCuboid2D<T> cuboid(extend,origin);

  T charL = 2*util::pow(radius, 2.)/ny;
  T physDeltaX = charL/N;

  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());

  // set periodic boundaries to the domain
  mesh.getCuboidDecomposition().setPeriodicity({ true, false });
  return mesh;
}

// Interpolation for large axis - smallest circle outside droplet
// This will only work for single ellipsoid droplets with stationary center
// Particularly in cases of droplet breakup, this will not work
std::vector<MyCase::value_t> circleApproximationBig(MyCase& myCase,
      MyCase::value_t minDistance,    // stop criterion: minimum distance for radii
      int points,       // number of points on the circle
      MyCase::value_t factor ){       // expansion and retraction of the circle with factor 1+factor or 1-factor
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;

  auto* sLatticeTest = &myCase.getLattice(Component1{});
  auto& params = myCase.getParameters();

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T nx = extent[0];
  const T ny = extent[1];
  const T radius = params.get<parameters::RADIUS>();

  SuperLatticeDensity2D<T, DESCRIPTOR> densityCircle( *sLatticeTest );

  AnalyticalFfromSuperF2D<T> aDensity(densityCircle, true, 1);

  T oldRadius = 0.;
  T newRadius = radius;

  T safety = 0.5;
  T coordinates[2] = {T(0),T(0)};
  T theta = 0.;

  while(std::abs(oldRadius-newRadius) > minDistance){
    T maxAngle = std::numbers::pi_v<T>/2.;      // third quadrant, as shear direction is known beforehand
    if( ny/2. + newRadius >= ny - safety) maxAngle = util::asin( ( ny/2. - safety) / newRadius );       // eliminate points in proximity to no-slip wall (phi=0 there)
    for(int i=0; i<=points-1; i++){
      coordinates[0] = nx/2. + util::cos(maxAngle * (T) i/ ((T) (points-1.)) ) * newRadius;
      coordinates[1] = ny/2. + util::sin(maxAngle * (T) i/ ((T) (points-1.)) ) * newRadius;
      T output[1];

      aDensity( output, coordinates );

      if(output[0] <= -0.0) {
        oldRadius = newRadius;
        newRadius *= 1. + factor;
        //theta is only assigned if the circle intersects the droplet
        theta = 90./std::numbers::pi_v<T>/2. * util::atan( (coordinates[1] - ny/2.) / (coordinates[0] - nx/2.) );
        break;
      } else if ( i == points - 1){
        oldRadius = newRadius;
        newRadius *= 1. - factor;
      }
    }
    factor *= 0.8;
  }
  return {oldRadius, theta};
}

// Interpolation for small axis - largest circle inside droplet
MyCase::value_t circleApproximationSmall(MyCase& myCase,
      MyCase::value_t minDistance,
      int points,
      MyCase::value_t factor){
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;

  auto& params = myCase.getParameters();

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T nx = extent[0];
  const T ny = extent[1];
  const T radius = params.get<parameters::RADIUS>();

  auto* sLatticeTest = &myCase.getLattice(Component1{});
  SuperLatticeDensity2D<T, DESCRIPTOR> densityCircle( *sLatticeTest );

  AnalyticalFfromSuperF2D<T> aDensity(densityCircle, true, 1);

  T oldRadius = 0.;
  T newRadius = radius;

  T coordinates[2] = {T(0),T(0)};

    while(std::abs(oldRadius-newRadius) > minDistance){
    T maxAngle = std::numbers::pi_v<T>;
    for(int i=0; i<=points-1; i++){
      coordinates[0] = nx/2. + util::cos(maxAngle * (T) i/ ((T) (points-1.)) ) * newRadius;
      coordinates[1] = ny/2. + util::sin(maxAngle * (T) i/ ((T) (points-1.)) ) * newRadius;
      T output[1];

      aDensity( output, coordinates );

      if(output[0] >= -0.0) {
        oldRadius = newRadius;
        newRadius *= 1. - factor;
        break;
      } else if ( i == points - 1){
        oldRadius = newRadius;
        newRadius *= 1. + factor;
      }
    }
    factor *= 0.8;
  }
  return oldRadius;
}

// The accuracy of this is dependent on the resolution
// Additionally droplet deformation is only valid for ellipsoidal droplets
MyCase::value_t measureDroplet(MyCase& myCase){
  OstreamManager clout( std::cout,"Measurement" );

  using T = MyCase::value_t;

  //axis approximation
  std::vector<T> circleMeasurementsBig = circleApproximationBig( myCase, 0.001, 5000, 0.8 );
  T circleMeasurementsSmall = circleApproximationSmall( myCase, 0.001, 5000, 0.8 );

  T deformation = (circleMeasurementsBig[0] - circleMeasurementsSmall) / (circleMeasurementsBig[0] + circleMeasurementsSmall);

  clout << "Circle Measurement:  " << " L: " << circleMeasurementsBig[0] << ", H: " << circleMeasurementsSmall << ", theta: " << circleMeasurementsBig[1] << ", deformation: " << deformation<< std::endl;
  return deformation;
}

void prepareGeometry( MyCase& myCase )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;
  using T = MyCase::value_t;

  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T nx = extent[0];
  const T ny = extent[1];
  const int N = params.get<parameters::RESOLUTION>();

  geometry.rename( 0, 2 );

  // bulk, MN=1
  geometry.rename( 2, 1, {0, 1} );

  T eps = nx/N;
  T edge = 0.5*eps;
  T safety = 2*eps;

  // top wall, MN=3
  std::vector<T> origin = {T(0) - safety, T(ny)-edge};
  std::vector<T> extend = {nx + 2*safety, edge+safety};
  IndicatorCuboid2D<T> top( extend, origin );
  geometry.rename(2, 3, top );

  // bottom wall, MN=4
  origin[1] = 0.0 - safety;
  IndicatorCuboid2D<T> bottom( extend, origin);
  geometry.rename(2, 4, bottom);

  // clean up
  geometry.clean();
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

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T ny = extent[1];
  const T radius = params.get<parameters::RADIUS>();
  const T shearRate = params.get<parameters::SHEAR_RATE>();
  const int N = params.get<parameters::RESOLUTION>();
  const T latticeRelaxationTime = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T contVisc = params.get<parameters::CONT_VISC>();
  const T physDensity = params.get<parameters::PHYS_CHAR_DENSITY>();
  const T capillaryNr = params.get<parameters::CAPILLARY>();
  const T cahnNr = params.get<parameters::CAHN>();
  const T pecletNr = params.get<parameters::PECLET>();
  const T dropRe = params.get<parameters::DROPLET_REYNOLDS>();
  const T physInterval = params.get<parameters::PHYS_INTERVAL>();
  const T maxPhysT = params.get<parameters::PHYS_MAX_T>();
  const T physStatIterT = params.get<parameters::PHYS_STAT_ITER_T>();
  const T physVtkIterT = params.get<parameters::PHYS_VTK_ITER_T>();
  const T residuum = params.get<parameters::CONVERGENCE_PRECISION>();
  const T maxNormTime = params.get<parameters::MAX_NORM_T>();

  T shearVelo = ny * shearRate / 2.;  // [pu] physical shear speed corr. to dropRe [m/s]
  T charL = 2*util::pow(radius, 2.)/ny;

  sLattice1.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    int   {N},                            // resolution
    (T)   latticeRelaxationTime,        // lattice relaxation time
    (T)   charL,   // charPhysLength: reference length of simulation geometry (dropRe!!)
    (T)   shearVelo,                    // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   contVisc,                     // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   physDensity                   // physDensity: physical density in __kg / m^3__
  );

  const auto& converter = sLattice1.getUnitConverter();
  converter.print();

  sLattice2.setUnitConverter(converter);

  // define lattice Dynamics
  dynamics::set<ForcedBGKdynamics>(sLattice1, geometry.getMaterialIndicator({1}));
  dynamics::set<FreeEnergyBGKdynamics>(sLattice2, geometry.getMaterialIndicator({1}));

  // moving walls (inlet boundary with tangetial velocity condition)
  T omega = converter.getLatticeRelaxationFrequency();

  auto topIndicator = geometry.getMaterialIndicator(3);
  boundary::set<boundary::FreeEnergyVelocity>(sLattice1, topIndicator);
  boundary::set<boundary::FreeEnergyOrderParameter>(sLattice2, topIndicator);

  auto bottomIndicator = geometry.getMaterialIndicator(4);
  boundary::set<boundary::FreeEnergyVelocity>(sLattice1, bottomIndicator);
  boundary::set<boundary::FreeEnergyOrderParameter>(sLattice2, bottomIndicator);

  sLattice1.setParameter<OMEGA>(omega);
  sLattice2.setParameter<OMEGA>(omega);

  // 1. compute surface tension from Ca
  const T surfTen = radius * shearRate * physDensity * contVisc / capillaryNr;  // [pu] physical surface tension
  const T surfTenLatt = surfTen / (converter.getConversionFactorPressure() * converter.getPhysDeltaX());  // [lu] lattice surface tension

  // 2. compute interface thickness from Ch
  const T xiThickness = cahnNr * converter.getLatticeLength(radius);  // [lu] characteristical thickness of interface (xi from Komrakova2013)

  // 3. compute kappa from xi (Komrakova2013 equation (12))
  const T kappaKomra = (3./4.) * surfTenLatt * xiThickness; // [lu] this is kappa from Komrakova2013

  // 4. compute A from kappa (Komrakova2013 equation (11), neglect minus sign, since solely multiplications here)
  const T aKomra = 2. * kappaKomra / util::pow(xiThickness, 2); // [lu] this is A from Komrakova2013

  // 5. compute M from Pe (Komrakova2013 equation (17), there's a typo in the paper, should be a minus on the right hand side... we dropped that anyway.)
  const T transportCoeff = shearRate * converter.getConversionFactorTime() * converter.getLatticeLength(radius) * xiThickness / (pecletNr * aKomra);  // [lu] this is M from Komrakova2013

  // 6. compute gamma from M (Komrakova2013 equation (8), caution: \Delta t is unity in lattice units!)
  const T gama = transportCoeff / (converter.getLatticeRelaxationTime() - .5);  // [lu] gama "Diffusivity of the interface" same as in Komrakova2013

  // Parameter fit from Komrakova to OpenLB (Semprebon2016):
  const T kappa1 = 4. * aKomra; // [lu] variable for interface tension tuning
  const T kappa2 = kappa1;  // [lu] variable for interface tension tuning
  const T alpha = xiThickness / 2.; // [lu] free energy alpha in surface tension

  const int maxIter = converter.getLatticeTime(maxPhysT);
  const int statIter = converter.getLatticeTime(physStatIterT);
  const int vtkIter = converter.getLatticeTime(physVtkIterT);

  // now print everything additionally to the unit converter (helps for double-checking)
  clout << "----------------------------------------------------------------------" << std::endl;
  clout << "----------------------------------------------------------------------" << std::endl;
  clout << "shear rate             " << shearRate      << "[pu], " << converter.getLatticeVelocity(shearRate) << "[lu]" << std::endl;
  clout << "shear rate             " << shearRate      << "[pu], " << shearRate * converter.getPhysDeltaT()   << "[lu]" << std::endl;
  clout << "Re_drop                " << dropRe         << std::endl;
  clout << "Ca                     " << capillaryNr    << std::endl;
  clout << "radius                 " << radius         << "[pu], " << converter.getLatticeLength(radius)      << "[lu]" << std::endl;
  clout << "normalized t_end       " << maxNormTime    << std::endl;
  clout << "surface tension        " << surfTen        << "[pu], " << surfTenLatt << "[lu]" << std::endl;
  clout << "Cahn Nr.               " << cahnNr         << std::endl;
  clout << "xi                     " << xiThickness    << "[lu]" << std::endl;
  clout << "alpha                  " << alpha          << "[lu]" << std::endl;
  clout << "kappa1                 " << kappa1         << "[lu]" << std::endl;
  clout << "kappa2                 " << kappa2         << "[lu]" << std::endl;
  clout << "Peclet Nr.             " << pecletNr       << std::endl;
  clout << "A                      " << aKomra         << "[lu]" << std::endl;
  clout << "Mobility M             " << transportCoeff << "[lu]" << std::endl;
  clout << "Gamma                  " << gama           << "[lu]" << std::endl;
  clout << "vtkIter                " << vtkIter        << "[lu]" << std::endl;
  clout << "statIter               " << statIter       << "[lu]" << std::endl;
  clout << "maxIter                " << maxIter        << "[lu]" << std::endl;
  clout << "physInterval convCheck " << physInterval   << "[pu]" << std::endl;
  clout << "residuum convCheck     " << residuum       << std::endl;
  clout << "----------------------------------------------------------------------" << std::endl;
  clout << "----------------------------------------------------------------------" << std::endl;

  sLattice2.setParameter<collision::FreeEnergy::GAMMA>(gama);

  // Add the lattice couplings
  // The chemical potential coupling must come before the force coupling
  auto& coupling1 = myCase.setCouplingOperator(
  "Chemical_potential",
  ChemicalPotentialCoupling2D{},
  names::A{}, sLattice1,
  names::B{}, sLattice2);

  coupling1.template setParameter<ChemicalPotentialCoupling2D::ALPHA>(alpha);
  coupling1.template setParameter<ChemicalPotentialCoupling2D::KAPPA1>(kappa1);
  coupling1.template setParameter<ChemicalPotentialCoupling2D::KAPPA2>(kappa2);

  auto& coupling2 = myCase.setCouplingOperator(
  "Force",
  ForceCoupling2D{},
  names::A{}, sLattice2,
  names::B{}, sLattice1);

  coupling1.restrictTo(geometry.getMaterialIndicator({1}));
  coupling2.restrictTo(geometry.getMaterialIndicator({1}));

  // Walls
  auto& coupling3 = myCase.setCouplingOperator(
  "Inlet_outlet",
  InletOutletCoupling2D{},
  names::A{}, sLattice2,
  names::B{}, sLattice1);

  coupling3.restrictTo(geometry.getMaterialIndicator({3,4}));

  sLattice1.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());
  sLattice2.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());

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

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues( MyCase& myCase )
{
  OstreamManager clout( std::cout,"initialValues" );

  using T = MyCase::value_t;

  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();
  auto& converter = myCase.getLattice(NavierStokes{}).getUnitConverter();

  auto& sLattice1 = myCase.getLattice(NavierStokes{});
  auto& sLattice2 = myCase.getLattice(Component1{});

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T nx = extent[0];
  const T ny = extent[1];
  const T radius = params.get<parameters::RADIUS>();
  const T cahnNr = params.get<parameters::CAHN>();

  // 1. compute interface thickness from Ch
  const T xiThickness = cahnNr * converter.getLatticeLength(radius);

  // 2. compute M from Pe (Komrakova2013 equation (17), there's a typo in the paper, should be a minus on the right hand side... we dropped that anyway.)

  // Parameter fit from Komrakova to OpenLB (Semprebon2016):
  const T alpha = xiThickness / 2.;

  // define circular domain for fluid 2
  std::vector<T> v( 2,T() );
  AnalyticalConst2D<T,T> zeroVelocity( v );

  AnalyticalConst2D<T,T> one ( 1. );

  SmoothIndicatorCircle2D<T,T> circle( {T(nx)/T(2), T(ny)/T(2)}, radius, converter.getPhysLength(alpha) );

  AnalyticalIdentity2D<T,T> rho( one );
  AnalyticalIdentity2D<T,T> phi( one - circle - circle );

  // shear velocity
  AnalyticalConst2D<T,T> uTop   ( converter.getCharPhysVelocity(), T( 0 ) );
  AnalyticalConst2D<T,T> uBottom( -1.*converter.getCharPhysVelocity(), T( 0 ) );

  // top
  auto topIndicator = geometry.getMaterialIndicator(3);
  momenta::setVelocity(sLattice1, topIndicator, uTop);

  // bottom
  auto bottomIndicator = geometry.getMaterialIndicator(4);
  momenta::setVelocity(sLattice1, bottomIndicator, uBottom);

  // equilibrium population initialization
  sLattice2.iniEquilibrium( geometry, 1, phi, zeroVelocity );

  sLattice1.initialize();
  sLattice2.initialize();

  sLattice1.communicate();
  sLattice2.communicate();
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{ }

void getResults(MyCase& myCase,
                std::size_t iT)
{
  OstreamManager clout( std::cout,"getResults" );

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;

  auto& sLattice1 = myCase.getLattice(NavierStokes{});
  auto& sLattice2 = myCase.getLattice(Component1{});
  auto& params = myCase.getParameters();
  auto& converter = myCase.getLattice(NavierStokes{}).getUnitConverter();

  const T physVtkIterT = params.get<parameters::PHYS_VTK_ITER_T>();
  const int vtkIter = converter.getLatticeTime(physVtkIterT);

  SuperVTMwriter2D<T> vtmWriter( "binaryShearFlow2d" );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice1 );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice1 );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Writes the VTK files
  if ( iT%vtkIter==0 ) {
    sLattice1.setProcessingContext(ProcessingContext::Evaluation);
    AnalyticalConst2D<T,T> half_( 0.5 );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> half(half_, sLattice1);

    SuperLatticeDensity2D<T, DESCRIPTOR> density1( sLattice1 );
    density1.getName() = "rho";
    SuperLatticeDensity2D<T, DESCRIPTOR> density2( sLattice2 );
    density2.getName() = "phi";

    SuperIdentity2D<T,T> c1 (half*(density1+density2));
    c1.getName() = "density-fluid-1";
    SuperIdentity2D<T,T> c2 (half*(density1-density2));
    c2.getName() = "density-fluid-2";

    // velocity (seperately and combined)
    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity1( sLattice1, converter );
    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity2( sLattice2, converter );
    SuperIdentity2D<T,T> velocityMittel (half*(velocity1+velocity2));
    velocityMittel.getName() = "velocity-Mittel";
    velocity1.getName() = "velocity1";
    velocity2.getName() = "velocity2";

    vtmWriter.addFunctor( density1 );
    vtmWriter.addFunctor( density2 );
    vtmWriter.addFunctor( c1 );
    vtmWriter.addFunctor( c2 );
    vtmWriter.addFunctor( velocityMittel );
    vtmWriter.write( iT );
  }

  myCase.report(iT);
}

void simulate(MyCase& myCase) {
  OstreamManager clout( std::cout,"simulate" );

  using T = MyCase::value_t;

  auto& converter = myCase.getLattice(NavierStokes{}).getUnitConverter();
  auto& params = myCase.getParameters();
  auto& sLattice1 = myCase.getLattice(NavierStokes{});
  auto& sLattice2 = myCase.getLattice(Component1{});
  auto& coupling1 = myCase.getOperator("Chemical_potential");
  auto& coupling2 = myCase.getOperator("Force");
  auto& coupling3 = myCase.getOperator("Inlet_outlet");

  const T physInterval = params.get<parameters::PHYS_INTERVAL>();
  const T maxPhysT = params.get<parameters::PHYS_MAX_T>();
  const T physStatIterT = params.get<parameters::PHYS_STAT_ITER_T>();
  const int maxIter = converter.getLatticeTime(maxPhysT);
  const int statIter = converter.getLatticeTime(physStatIterT);
  const T residuum = params.get<parameters::CONVERGENCE_PRECISION>();
  const T defResiduum = params.get<parameters::RESIDUUM>();

  int iT = 0;
  clout << "starting simulation..." << std::endl;
  util::ValueTracer<T> converge( converter.getLatticeTime( physInterval ), residuum );

  T old_deformation = 0;
  T deformation = 0;

  for ( iT=0; iT<=maxIter; ++iT ) {

    sLattice2.setProcessingContext(ProcessingContext::Evaluation);
    sLattice1.setProcessingContext(ProcessingContext::Evaluation);
    // Measure droplet
    if ( iT%statIter==0 ) {
      old_deformation = deformation;
      clout<<"--------------- t: "<< converter.getPhysTime(iT) <<" ----------------"<<std::endl;
      deformation = measureDroplet( myCase );
    }

    // Convergence check
    if ( converge.hasConverged() && util::abs(deformation - old_deformation) < defResiduum ) {
      clout << "Simulation converged." << std::endl;
      getResults( myCase, iT );
      break;
    }
    sLattice1.setProcessingContext(ProcessingContext::Simulation);
    sLattice2.setProcessingContext(ProcessingContext::Simulation);

    // Computation and output of the results
    getResults( myCase, iT );

    // Collide and stream execution
    sLattice1.collideAndStream();
    sLattice2.collideAndStream();

    // MPI communication for lattice data
    sLattice1.communicate();
    sLattice2.communicate();

    sLattice1.executePostProcessors(stage::PreCoupling());
    sLattice2.executePostProcessors(stage::PreCoupling());

    // Execute coupling between the two lattices
    coupling1.apply();
    coupling2.apply();
    coupling3.apply();

  converge.takeValue( sLattice2.getStatistics().getAverageRho(), true );
  }

  myCase.finalize(maxIter);
}
int main( int argc, char *argv[] )
{
  initialize( &argc, &argv );

  // === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<parameters::RADIUS     >(1.);    // [pu] physical droplet radius
    myCaseParameters.set<RESOLUTION             >(10);    // resolution N cells per charL = radius/2
    myCaseParameters.set<CONT_VISC              >(1.);    // [pu] physical(!) kinematic viscosity of continuous phase
    myCaseParameters.set<PHYS_CHAR_DENSITY      >(1.);    // [pu] physical density
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(1.);    // [lu] tau for first lattice (f)
    myCaseParameters.set<DROPLET_REYNOLDS       >(1.);    // definition: shearRate * radius^2 / continuousViscosity
    myCaseParameters.set<CAPILLARY              >(1.);    // capillary number for droplet breakup
    myCaseParameters.set<CAHN                   >(0.057); // Cahn Number = xi / radius;
    myCaseParameters.set<PECLET                 >(.43);   // Peclet Number = shearRate * radius * xi / (M A)
    myCaseParameters.set<MAX_NORM_T             >(30.);   // normalized time normT = T * shearRate [either in lu or pu]
    myCaseParameters.set<NORM_VTK_ITER_T        >(0.2);   // normalized vtk output (i.e. per normalized time step, one output!)
    myCaseParameters.set<NORM_STAT_ITER_T       >(0.2);   // normalized stat output (i.e. per normalized time step, one output!)
    myCaseParameters.set<PHYS_INTERVAL          >(0.1);   // interval for the convergence check in s
    myCaseParameters.set<CONVERGENCE_PRECISION  >(1e-5);  // residuum for the convergence check
    myCaseParameters.set<H1                     >(0.);    // Contact angle 90 degrees   [lattice units]
    myCaseParameters.set<H2                     >(0.);    // Contact angle 90 degrees   [lattice units]
    myCaseParameters.set<CONVERGENCE_PRECISION  >(0.0015);  // deformation residuum criterion
    myCaseParameters.set<parameters::SHEAR_RATE>([] (MyCase::ParametersD& p) {
      return p.get<CONT_VISC>() * p.get<DROPLET_REYNOLDS>() / util::pow(p.get<parameters::RADIUS>(), 2);
    }); // [pu] physical shear rate [1/s]
    myCaseParameters.set<parameters::PHYS_MAX_T>([] (MyCase::ParametersD& p) {
      return p.get<MAX_NORM_T>() / p.get<SHEAR_RATE>();
    });
    myCaseParameters.set<parameters::PHYS_VTK_ITER_T>([] (MyCase::ParametersD& p) {
      return p.get<NORM_VTK_ITER_T>() / p.get<SHEAR_RATE>();
    });
    myCaseParameters.set<parameters::PHYS_STAT_ITER_T>([] (MyCase::ParametersD& p) {
      return p.get<NORM_STAT_ITER_T>() / p.get<SHEAR_RATE>();
    });
    myCaseParameters["reporter"]["statistics"].set<parameters::PHYS_ITER_T>([] (MyCase::ParametersD& p) {
      return p.root().get<NORM_STAT_ITER_T>() / p.root().get<SHEAR_RATE>();
    });
    myCaseParameters["reporter"]["timer"].set<parameters::PHYS_ITER_T>([] (MyCase::ParametersD& p) {
      return p.root().get<NORM_STAT_ITER_T>() / p.root().get<SHEAR_RATE>();
    });
    myCaseParameters.set<parameters::DOMAIN_EXTENT>([] (MyCase::ParametersD& p) -> Vector<MyCase::value_t, 2> {
      return {20 * p.get<parameters::RADIUS>(), 4 * p.get<parameters::RADIUS>()};
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
