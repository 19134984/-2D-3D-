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

/* fourRollMill2d.cpp
 * A spherical domain of fluid phase II is immersed in a cuboid
 * filled with fluid phase I. The bottom left and top right
 * cylinders begin to spin in counterclock-wise direction.
 * The top left and bottom right cylinders spin in clock-wise direction-
 * The velocity field of extensional type deforms the droplet accordingly.
 * The flow configuration is taken from
 *  S. Simonis et al. Discrete and Continuous Dynamical
 *  Systems - Series S (2023), doi:10.3934/dcdss.2023069
 * The present example uses the ternary free energy model by
 *  C Semprebon et al. Physical Review E 93.3 (2016) p. 033305.
 *
 * The droplet breakup case should be run on parallel mode
 * because of increased computation time.
 */

#include <olb.h>

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

struct CONT_VISC  : public descriptors::FIELD_BASE<1> { };
struct DROPLET_REYNOLDS : public descriptors::FIELD_BASE<1> { };
struct CAPILLARY : public descriptors::FIELD_BASE<1> { };
struct CAHN  : public descriptors::FIELD_BASE<1> { };
struct H1  : public descriptors::FIELD_BASE<1> { };
struct H2  : public descriptors::FIELD_BASE<1> { };
}


Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T nx = extent[0];
  const T N = params.get<parameters::RESOLUTION>();
  const T radius = params.get<parameters::RADIUS>();

  std::vector<T> extend = { nx, nx };
  std::vector<T> origin = { 0, 0 };

  //Set physical input parameters
  const T cylinderOffset = nx / 2. / 3.;              // [pu] cylinder offset from x- and y- axis, symmetrical!
  const T radiusCylinder = cylinderOffset * 0.625;    // [pu] radius of the cylinder

  // Center of rotating cylinders
  const T cylinderOffsetMinus = (nx/2. - cylinderOffset);
  const T cylinderOffsetPlus = (nx/2. + cylinderOffset);
  std::vector<T> center1 = {cylinderOffsetMinus, cylinderOffsetMinus};
  std::vector<T> center2 = {cylinderOffsetPlus, cylinderOffsetMinus};
  std::vector<T> center3 = {cylinderOffsetPlus, cylinderOffsetPlus};
  std::vector<T> center4 = {cylinderOffsetMinus, cylinderOffsetPlus};

  std::shared_ptr<IndicatorF2D<T>> cuboid = std::make_shared<IndicatorCuboid2D<T>>( extend, origin );
  std::shared_ptr<IndicatorF2D<T>> cylind1 = std::make_shared<IndicatorCircle2D<T>>( center1, radiusCylinder );
  std::shared_ptr<IndicatorF2D<T>> cylind2 = std::make_shared<IndicatorCircle2D<T>>( center2, radiusCylinder );
  std::shared_ptr<IndicatorF2D<T>> cylind3 = std::make_shared<IndicatorCircle2D<T>>( center3, radiusCylinder );
  std::shared_ptr<IndicatorF2D<T>> cylind4 = std::make_shared<IndicatorCircle2D<T>>( center4, radiusCylinder );

  T charL = util::pow(radius, 2.);
  T physDeltaX = charL/N;

  Mesh<T,MyCase::d> mesh(*(cuboid-(cylind1+cylind2+cylind3+cylind4)), physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());

  // set periodic boundaries to the domain
  mesh.getCuboidDecomposition().setPeriodicity({ true, true });
  return mesh;
}

// Functors for rotating cylinders;
// cylinders rotate with tangential speed rotVel[lu]
template <typename T>
class RotatingMill2D : public AnalyticalF2D<T,T> {
  public:
    RotatingMill2D(std::vector<T> center, T radius, T rotVel, bool counterClock);
    bool operator()(T output[], const T input[]) override;
  private:
    T _center[2];
    T _radius;
    T _rotVel;
    bool _counterClock;
};

template <typename T>
RotatingMill2D<T>::RotatingMill2D(std::vector<T> center, T radius, T rotVel, bool counterClock) : AnalyticalF2D<T,T>(2){
  this->_center[0] = center[0];
  this->_center[1] = center[1];
  this->_radius = radius;
  this->_rotVel = rotVel;
  this->_counterClock = counterClock;
}

template <typename T>
bool RotatingMill2D<T>::operator()(T output[], const T input[]){
  output[0] = -(input[1]-_center[1])*_rotVel/_radius;
  output[1] =  (input[0]-_center[0])*_rotVel/_radius;
  if(!_counterClock){
    output[0] =  (input[1]-_center[1])*_rotVel/_radius;
    output[1] = -(input[0]-_center[0])*_rotVel/_radius;
  }
  return true;
}

// Measurement of droplet: interpolation in x-direction
// Simple interpolation on y = ny/2 with assumption of droplet center in domain center
MyCase::value_t horizontalMeasure(MyCase& myCase, int points){
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;

  auto* latticeTest = &myCase.getLattice(Component1{});
  auto& params = myCase.getParameters();

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T nx = extent[0];

  SuperLatticeDensity2D<T, DESCRIPTOR> density2( *latticeTest );

  AnalyticalFfromSuperF2D<T> aDensityH(density2, true, 1);

  T dropletLength = 0.;
  T coordinates[2] = {T(0),T(nx)/T(2)}; //interpolation coordinates, y fixed

  for(int i=0; i<=points-1; i++){
    coordinates[0] = nx/2. + ((nx/2.) * (T) i/ ((T) (points-1.)));
    T output[1];
    aDensityH( output, coordinates );

    if(output[0] >= 0.0){           //fluid equals continuous phase
      dropletLength = coordinates[0] - nx/2.;
      break;
    }
  }
  return dropletLength;
}

// Measurement of droplet: interpolation in y-direction
// Simple interpolation on x = nx/2 with assumption of droplet center in domain center
MyCase::value_t verticalMeasure(MyCase& myCase, int points){
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;

  auto* latticeTest = &myCase.getLattice(Component1{});
  auto& params = myCase.getParameters();

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T nx = extent[0];

  SuperLatticeDensity2D<T, DESCRIPTOR> density2( *latticeTest );

  AnalyticalFfromSuperF2D<T> aDensityV(density2, true, 1);

  T dropletLength = 0.;
  T coordinates[2] = {T(nx)/T(2),T(0)}; //interpolation coordinates, x fixed

  for(int i=0; i<=points-1; i++){
    coordinates[1] = nx/2. + ((nx/2.) * (T) i/ ((T) (points-1.)));
    T output[1];
    aDensityV( output, coordinates );

    if(output[0] >= 0.0){           //fluid equals continuous phase
      dropletLength = coordinates[1] - nx/2.;
      break;
    }
  }
  return dropletLength;
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

  //Set physical input parameters
  const T cylinderOffset = nx / 2. / 3.;              // [pu] cylinder offset from x- and y- axis, symmetrical!
  const T radiusCylinder = cylinderOffset * 0.625;    // [pu] radius of the cylinder

  // Center of rotating cylinders
  const T cylinderOffsetMinus = (nx/2. - cylinderOffset);
  const T cylinderOffsetPlus = (nx/2. + cylinderOffset);
  std::vector<T> center1 = {cylinderOffsetMinus, cylinderOffsetMinus};
  std::vector<T> center2 = {cylinderOffsetPlus, cylinderOffsetMinus};
  std::vector<T> center3 = {cylinderOffsetPlus, cylinderOffsetPlus};
  std::vector<T> center4 = {cylinderOffsetMinus, cylinderOffsetPlus};

  std::shared_ptr<IndicatorF2D<T>> cylind1 = std::make_shared<IndicatorCircle2D<T>>( center1, radiusCylinder );
  std::shared_ptr<IndicatorF2D<T>> cylind2 = std::make_shared<IndicatorCircle2D<T>>( center2, radiusCylinder );
  std::shared_ptr<IndicatorF2D<T>> cylind3 = std::make_shared<IndicatorCircle2D<T>>( center3, radiusCylinder );
  std::shared_ptr<IndicatorF2D<T>> cylind4 = std::make_shared<IndicatorCircle2D<T>>( center4, radiusCylinder );

  geometry.rename( 0,2 );
  geometry.rename( 2,1,{0,0} );
  geometry.rename( 1,3,cylind1 );    //Material number 3: cylinder bottom left
  geometry.rename( 1,4,cylind2 );    //Material number 4: cylinder bottom right
  geometry.rename( 1,5,cylind3 );    //Material number 5: cylinder top right
  geometry.rename( 1,6,cylind4 );    //Material number 6: cylinder top left

  // clean up
  geometry.clean();
  geometry.innerClean();
  geometry.checkForErrors();
  geometry.getStatistics().print();

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

  auto& lattice1 = myCase.getLattice(NavierStokes{});
  auto& lattice2 = myCase.getLattice(Component1{});

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T nx = extent[0];
  const T radius = params.get<parameters::RADIUS>();
  const T strainRate = params.get<parameters::STRAIN_RATE>();
  const int N = params.get<parameters::RESOLUTION>();
  const T latticeRelaxationTime = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T contVisc = params.get<parameters::CONT_VISC>();
  const T physDensity = params.get<parameters::PHYS_CHAR_DENSITY>();
  const T capillaryNr = params.get<parameters::CAPILLARY>();
  const T cahnNr = params.get<parameters::CAHN>();
  const T pecletNr = params.get<parameters::PECLET>();
  const T dropRe = params.get<parameters::DROPLET_REYNOLDS>();
  const T maxPhysT = params.get<parameters::PHYS_MAX_T>();
  const T physStatIterT = params["reporter"]["statistics"].get<parameters::PHYS_ITER_T>();
  const T physVtkIterT = params.get<parameters::PHYS_VTK_ITER_T>();

  // apparatus specific constant (~0.35) is only for this simulations setup with this cylinderRadius and offset
  // has to be calculated empirically: https://doi.org/10.1098/rspa.1934.0169
  T rotVelo = strainRate / 0.35;
  T charL = util::pow(radius, 2.);

  lattice1.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    (T)   N,                     // resolution
    (T)   latticeRelaxationTime, // lattice relaxation time
    (T)   charL,                 // charPhylength: reference length of simulation geometry (dropRe!!)
    (T)   rotVelo,               // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   contVisc,              // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   physDensity            // physDensity: physical density in __kg / m^3__
  );

  const auto& converter = lattice1.getUnitConverter();
  converter.print();

  lattice2.setUnitConverter(converter);

  //Set physical input parameters
  const T cylinderOffset = nx / 2. / 3.;              // [pu] cylinder offset from x- and y- axis, symmetrical!
  const T radiusCylinder = cylinderOffset * 0.625;    // [pu] radius of the cylinder

  const T deborahNumber = strainRate * converter.getLatticeRelaxationTime();

  // 1. compute surface tension from Ca
  const T surfTen = radius * strainRate * physDensity * contVisc / capillaryNr;

  const T surfTenLatt = surfTen / (converter.getConversionFactorPressure() * converter.getPhysDeltaX());

  // 2. compute interface thickness from Ch
  const T xiThickness = cahnNr * converter.getLatticeLength(radius);

  // 3. compute kappa from xi (Komrakova2013 equation (12))
  const T kappaKomra = (3./4.) * surfTenLatt * xiThickness;

  // 4. compute A from kappa (Komrakova2013 equation (11), neglect minus sign, since solely multiplications here)
  const T aKomra = 2. * kappaKomra / util::pow(xiThickness, 2);

  // 5. compute M from Pe (Komrakova2013 equation (17), there's a typo in the paper, should be a minus on the right hand side... we dropped that anyway.)
  const T transportCoeff = strainRate * converter.getConversionFactorTime() * converter.getLatticeLength(radius) * xiThickness / (pecletNr * aKomra);

  // 6. compute gamma from M (Komrakova2013 equation (8), caution: \Delta t is unity in lattice units!)
  const T gama = transportCoeff / (converter.getLatticeRelaxationTime() - .5);

  // Parameter fit from Komrakova to OpenLB (Semprebon2016):
  const T kappa1 = 4. * aKomra;
  const T kappa2 = kappa1;
  const T alpha = xiThickness / 2.;

  // transform normalized time to physical and then to lattice:
  const T maxIter = converter.getLatticeTime( maxPhysT );
  const T vtkIter = converter.getLatticeTime( physVtkIterT );
  const T statIter = converter.getLatticeTime( physStatIterT );;

  T omega = converter.getLatticeRelaxationFrequency();

  dynamics::set<ForcedBGKdynamics>( lattice1, geometry.getMaterialIndicator(1) );
  dynamics::set<FreeEnergyBGKdynamics>( lattice2, geometry.getMaterialIndicator(1) );

  auto cylinderIndicator1 = geometry.getMaterialIndicator(3);
  boundary::set<boundary::FreeEnergyVelocity>(lattice1, cylinderIndicator1);
  boundary::set<boundary::FreeEnergyOrderParameter>(lattice2, cylinderIndicator1);
  auto cylinderIndicator2 = geometry.getMaterialIndicator(4);
  boundary::set<boundary::FreeEnergyVelocity>(lattice1, cylinderIndicator2);
  boundary::set<boundary::FreeEnergyOrderParameter>(lattice2, cylinderIndicator2);
  auto cylinderIndicator3 = geometry.getMaterialIndicator(5);
  boundary::set<boundary::FreeEnergyVelocity>(lattice1, cylinderIndicator3);
  boundary::set<boundary::FreeEnergyOrderParameter>(lattice2, cylinderIndicator3);
  auto cylinderIndicator4 = geometry.getMaterialIndicator(6);
  boundary::set<boundary::FreeEnergyVelocity>(lattice1, cylinderIndicator4);
  boundary::set<boundary::FreeEnergyOrderParameter>(lattice2, cylinderIndicator4);

  lattice1.setParameter<OMEGA>(omega);
  lattice2.setParameter<OMEGA>(omega);
  lattice2.setParameter<collision::FreeEnergy::GAMMA>(gama);
  //
  // now print everything additionally to the unit converter (helps for double-checking)
  clout << "----------------------------------------------------------------------" << std::endl;
  clout << "----------------------------------------------------------------------" << std::endl;
  clout << "strain rate           " << strainRate     << "[pu], " << converter.getLatticeVelocity(strainRate) << "[lu]" << std::endl;
  clout << "Re_drop               " << dropRe         << std::endl;
  clout << "Deborah               " << deborahNumber  << std::endl;
  clout << "Ca                    " << capillaryNr    << std::endl;
  clout << "radius                " << radius         << "[pu], " << converter.getLatticeLength(radius) << "[lu]" << std::endl;
  clout << "surface tension       " << surfTen        << "[pu], " << surfTenLatt << "[lu]" << std::endl;
  clout << "Cahn Nr.              " << cahnNr         << std::endl;
  clout << "xi                    " << xiThickness    << "[lu]" << std::endl;
  clout << "kappa1                " << kappa1         << "[lu]" << std::endl;
  clout << "kappa2                " << kappa2         << "[lu]" << std::endl;
  clout << "Peclet Nr.            " << pecletNr       << std::endl;
  clout << "A                     " << aKomra         << "[lu]" << std::endl;
  clout << "Mobility M            " << transportCoeff << "[lu]" << std::endl;
  clout << "Gamma                 " << gama           << "[lu]" << std::endl;
  clout << "cylinderRadius        " << radiusCylinder << "[pu], " << converter.getLatticeLength(radiusCylinder) << "[lu]" << std::endl;
  clout << "vtkIter               " << vtkIter        << "[lu]" << std::endl;
  clout << "statIter              " << statIter       << "[lu]" << std::endl;
  clout << "maxIter               " << maxIter        << "[lu]" << std::endl;
  clout << "----------------------------------------------------------------------" << std::endl;
  clout << "----------------------------------------------------------------------" << std::endl;

  //prepareCoupling( lattice1, lattice2, geometry );
  auto& coupling1 = myCase.setCouplingOperator(
  "Chemical_potential",
  ChemicalPotentialCoupling2D{},
  names::A{}, lattice1,
  names::B{}, lattice2);

  coupling1.restrictTo(geometry.getMaterialIndicator({1}));

  coupling1.template setParameter<ChemicalPotentialCoupling2D::ALPHA>(alpha);
  coupling1.template setParameter<ChemicalPotentialCoupling2D::KAPPA1>(kappa1);
  coupling1.template setParameter<ChemicalPotentialCoupling2D::KAPPA2>(kappa2);

  auto& coupling2 = myCase.setCouplingOperator(
  "Force",
  ForceCoupling2D{},
  names::A{}, lattice2,
  names::B{}, lattice1);

  coupling2.restrictTo(geometry.getMaterialIndicator({1}));

  // Walls
  auto& coupling3 = myCase.setCouplingOperator(
  "Inlet_outlet",
  InletOutletCoupling2D{},
  names::A{}, lattice2,
  names::B{}, lattice1);

  coupling3.restrictTo(geometry.getMaterialIndicator({3,4,5,6}));

  lattice1.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());
  lattice2.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());

  {
    auto& communicator = lattice1.getCommunicator(stage::PostCoupling());
    communicator.requestField<CHEM_POTENTIAL>();
    communicator.requestOverlap(lattice1.getOverlap());
    communicator.exchangeRequests();
  }
  {
    auto& communicator = lattice2.getCommunicator(stage::PreCoupling());
    communicator.requestField<CHEM_POTENTIAL>();
    communicator.requestField<RhoStatistics>();
    communicator.requestOverlap(lattice2.getOverlap());
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

  auto& lattice1 = myCase.getLattice(NavierStokes{});
  auto& lattice2 = myCase.getLattice(Component1{});

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T nx = extent[0];
  const T radius = params.get<parameters::RADIUS>();
  const T cahnNr = params.get<parameters::CAHN>();

  const T xiThickness = cahnNr * converter.getLatticeLength(radius);
  const T alpha = xiThickness / 2.;

  //Set physical input parameters
  const T cylinderOffset = nx / 2. / 3.;              // [pu] cylinder offset from x- and y- axis, symmetrical!
  const T radiusCylinder = cylinderOffset * 0.625;    // [pu] radius of the cylinder

  // Center of rotating cylinders
  const T cylinderOffsetMinus = (nx/2. - cylinderOffset);
  const T cylinderOffsetPlus = (nx/2. + cylinderOffset);
  std::vector<T> center1 = {cylinderOffsetMinus, cylinderOffsetMinus};
  std::vector<T> center2 = {cylinderOffsetPlus, cylinderOffsetMinus};
  std::vector<T> center3 = {cylinderOffsetPlus, cylinderOffsetPlus};
  std::vector<T> center4 = {cylinderOffsetMinus, cylinderOffsetPlus};

  // bulk initial conditions
  // define circular domain for fluid 2
  std::vector<T> v( 2,T() );
  AnalyticalConst2D<T,T> zeroVelocity( v );

  AnalyticalConst2D<T,T> one ( 1. );

  SmoothIndicatorCircle2D<T,T> circle( {T(nx)/T(2), T(nx)/T(2)}, radius, converter.getPhysLength(alpha) );
  AnalyticalIdentity2D<T,T> rho( one );
  AnalyticalIdentity2D<T,T> phi( one - circle - circle );

  RotatingMill2D<T> rotator1 ( center1, radiusCylinder, converter.getCharPhysVelocity(), true);
  RotatingMill2D<T> rotator2 ( center2, radiusCylinder, converter.getCharPhysVelocity(), false);
  RotatingMill2D<T> rotator3 ( center3, radiusCylinder, converter.getCharPhysVelocity(), true);
  RotatingMill2D<T> rotator4 ( center4, radiusCylinder, converter.getCharPhysVelocity(), false);

  auto cylinderIndicator1 = geometry.getMaterialIndicator(3);
  auto cylinderIndicator2 = geometry.getMaterialIndicator(4);
  auto cylinderIndicator3 = geometry.getMaterialIndicator(5);
  auto cylinderIndicator4 = geometry.getMaterialIndicator(6);

  momenta::setVelocity( lattice1, cylinderIndicator1, rotator1 );
  momenta::setVelocity( lattice1, cylinderIndicator2, rotator2 );
  momenta::setVelocity( lattice1, cylinderIndicator3, rotator3 );
  momenta::setVelocity( lattice1, cylinderIndicator4, rotator4 );
  momenta::setDensity ( lattice1, cylinderIndicator1, phi );
  momenta::setDensity ( lattice1, cylinderIndicator2, phi );
  momenta::setDensity ( lattice1, cylinderIndicator3, phi );
  momenta::setDensity ( lattice1, cylinderIndicator4, phi );

  // equilibrium population initialization
  lattice2.iniEquilibrium( geometry, 1, phi, zeroVelocity );

  lattice1.initialize();
  lattice2.initialize();

  lattice1.communicate();
  lattice2.communicate();
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

  auto& lattice1 = myCase.getLattice(NavierStokes{});
  auto& lattice2 = myCase.getLattice(Component1{});
  auto& params = myCase.getParameters();
  auto& converter = myCase.getLattice(NavierStokes{}).getUnitConverter();

  const T physStatIterT = params["reporter"]["statistics"].get<parameters::PHYS_ITER_T>();
  const int statIter = converter.getLatticeTime(physStatIterT);

  SuperVTMwriter2D<T> vtmWriter( "fourRollMill2D" );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( lattice1 );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( lattice1 );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Writes the VTK files
   if ( iT%statIter==0 ) {
    lattice1.setProcessingContext(ProcessingContext::Evaluation);
    AnalyticalConst2D<T,T> half_( 0.5 );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> half(half_, lattice1);

    SuperLatticeDensity2D<T, DESCRIPTOR> density1( lattice1 );
    density1.getName() = "rho";
    SuperLatticeDensity2D<T, DESCRIPTOR> density2( lattice2 );
    density2.getName() = "phi";

    SuperIdentity2D<T,T> c1 (half*(density1+density2));
    c1.getName() = "density-fluid-1";
    SuperIdentity2D<T,T> c2 (half*(density1-density2));
    c2.getName() = "density-fluid-2";

    // velocity (seperately and combined)
    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity1( lattice1, converter );
    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity2( lattice2, converter );
    SuperIdentity2D<T,T> velocityMittel (half*(velocity1+velocity2));
    velocityMittel.getName() = "velocity-Mittel";

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
  auto& lattice1 = myCase.getLattice(NavierStokes{});
  auto& lattice2 = myCase.getLattice(Component1{});
  auto& coupling1 = myCase.getOperator("Chemical_potential");
  auto& coupling2 = myCase.getOperator("Force");
  auto& coupling3 = myCase.getOperator("Inlet_outlet");

  const int N = params.get<parameters::RESOLUTION>();
  const T radius = params.get<parameters::RADIUS>();
  const T physInterval = params.get<parameters::PHYS_INTERVAL_T>();
  const T maxPhysT = params.get<parameters::PHYS_MAX_T>();
  const T physStatIterT = params["reporter"]["statistics"].get<parameters::PHYS_ITER_T>();
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

    lattice2.setProcessingContext(ProcessingContext::Evaluation);
    lattice1.setProcessingContext(ProcessingContext::Evaluation);

    // doplet measurement
  if ( iT%statIter==0 ) {
      clout<<"--------------- t: "<< converter.getPhysTime(iT) <<" ----------------"<<std::endl;
    old_deformation = deformation;
      T lengthHorizontal = horizontalMeasure( myCase, N * 1000);
      T lengthVertical   = verticalMeasure(   myCase, N * 1000);
    deformation = (lengthHorizontal - lengthVertical) / (lengthHorizontal + lengthVertical);
      clout << "Length(horizontal): " << lengthHorizontal << std::endl;
      clout << "Length(vertical):   " << lengthVertical   << std::endl;
      clout << "Deformation:        " << deformation << std::endl;
      clout << "Length/radius       " << lengthHorizontal / radius << std::endl;
  }

    // Convergence check
    if ( converge.hasConverged() && util::abs(deformation - old_deformation) < defResiduum ) {
      clout << "Simulation converged." << std::endl;
      getResults( myCase, iT );
      break;
    }
    lattice1.setProcessingContext(ProcessingContext::Simulation);
    lattice2.setProcessingContext(ProcessingContext::Simulation);

    // Computation and output of the results
    getResults( myCase, iT );

    // Collide and stream execution
    lattice1.collideAndStream();
    lattice2.collideAndStream();

    lattice1.executePostProcessors(stage::PreCoupling());
    lattice2.executePostProcessors(stage::PreCoupling());

    // Execute coupling between the two lattices
    lattice1.getCommunicator(stage::PreCoupling()).communicate();
    coupling1.apply();
    lattice1.getCommunicator(stage::PostCoupling()).communicate();

    lattice2.executePostProcessors(stage::PreCoupling());

    lattice2.getCommunicator(stage::PreCoupling()).communicate();
    coupling2.apply();
    lattice2.getCommunicator(stage::PostCoupling()).communicate();

    coupling3.apply();

    converge.takeValue( lattice2.getStatistics().getAverageRho(), true );
  }

  myCase.finalize(maxIter);
}

int main( int argc, char *argv[] )
{
  // === 1st Step: Initialization ===
  initialize( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

// Guideline parameters:
// Ca = 0.1  for steady state
// Ca = 0.42 for drop breakup
// Recommended xi = 1.14, min 1.0

  // === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<parameters::RADIUS     >(1.);        // [pu] physical droplet radius
    myCaseParameters.set<DOMAIN_EXTENT          >({20, 20});  // [pu] physical domain size
    myCaseParameters.set<RESOLUTION             >(10);        // resolution N cells per charL = radius/2
    myCaseParameters.set<CONT_VISC              >(1.);        // [pu] physical(!) kinematic viscosity of continuous phase
    myCaseParameters.set<PHYS_CHAR_DENSITY      >(1.);        // [pu] physical density
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(1.);        // [lu] tau for first lattice (f)
    myCaseParameters.set<DROPLET_REYNOLDS       >(0.0625);    // definition: strainRate * radius^2 / continuousViscosity
    myCaseParameters.set<CAPILLARY              >(0.42);      // capillary number for droplet breakup
    myCaseParameters.set<CAHN                   >(0.114);     // Cahn Number = xi / radius;
    myCaseParameters.set<PECLET                 >(0.43);      // Peclet Number = strainRate * radius * xi / (M A)
    myCaseParameters.set<PHYS_MAX_T             >(160.);
    myCaseParameters.set<PHYS_VTK_ITER_T        >(1.);
    myCaseParameters["reporter"]["statistics"].set<PHYS_ITER_T>(1.);
    myCaseParameters["timer"     ].set<PHYS_TIMER_ITER_T>(1.);
    myCaseParameters.set<PHYS_INTERVAL_T      >(1.);        // interval for the convergence check in s
    myCaseParameters.set<CONVERGENCE_PRECISION>(7e-5);      // residuum for the convergence check
    myCaseParameters.set<H1                   >(0.);        // Contact angle 90 degrees   [lattice units]
    myCaseParameters.set<H2                   >(0.);        // Contact angle 90 degrees   [lattice units]
    myCaseParameters.set<CONVERGENCE_PRECISION>(0.001);  // deformation residuum criterion
    myCaseParameters.set<parameters::STRAIN_RATE>([] (MyCase::ParametersD& p) {
      return p.get<CONT_VISC>() * p.get<DROPLET_REYNOLDS>() * 0.35 / util::pow(p.get<parameters::RADIUS>(), 2);
    }); // [pu] physical strain rate [1/s]
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
