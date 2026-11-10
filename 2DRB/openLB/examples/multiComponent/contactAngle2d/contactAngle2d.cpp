/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Michael Rennick, Tim Bingert
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

/* contactAngle2d.cpp
 * This example shows a liquid droplet in contact with a solid wall
 * under a certain contact angle using the L. Ju et al. (2025) Well-Balanced
 * model with an energy contribution to the free energy for the order parameter
 * boundary condition. This model benefits from low spurious velocity and can
 * handle large density and viscosity ratios.
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, D2Q9<RHO,NABLARHO,FORCE,EXTERNAL_FORCE,TAU_EFF>>,
  Component1,  Lattice<double, D2Q9<FORCE,SOURCE,SOURCE_OLD,PHIWETTING,VELOCITY,OLD_PHIU,STATISTIC,CHEM_POTENTIAL,BOUNDARY>>
>;

using NSBulkDynamics = MultiPhaseIncompressibleBGKdynamics<MyCase::value_t, MyCase::descriptor_t>;
using CHBulkDynamics = WellBalancedCahnHilliardBGKdynamics<MyCase::value_t_of<Component1>, MyCase::descriptor_t_of<Component1>>;
using MixtureRules = LinearTauViscosity;
using Coupling = WellBalancedCahnHilliardPostProcessor<MixtureRules>;

namespace olb::parameters {

struct LATTICE_RELAXATION_TIME_PF : public descriptors::FIELD_BASE<1> { };
struct LATTICE_RELAXATION_TIME_2  : public descriptors::FIELD_BASE<1> { };
struct C_RHO                      : public descriptors::FIELD_BASE<1> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extent, origin);

  const T dx = 1.; // lattice units case
  Mesh<T,MyCase::d> mesh(cuboid, dx, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({true,false});
  return mesh;
}

// Labels for boundary and fluid locations
void prepareGeometry( MyCase& myCase )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;
  auto& geometry = myCase.getGeometry();
  // Fluid nodes labelled 2
  geometry.rename( 0,2 );
  // Label edges as 1
  geometry.rename( 2, 1, {0, 1} );

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
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& sLatticeNS = myCase.getLattice(NavierStokes{});
  auto& sLatticeCH = myCase.getLattice(Component1{});

  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

  const Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const int N = params.get<parameters::RESOLUTION>();
  const T tau_l = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T tau_g = params.get<parameters::LATTICE_RELAXATION_TIME_2>();
  const T tau_mobil = params.get<parameters::LATTICE_RELAXATION_TIME_PF>();
  const T L_char = params.get<parameters::PHYS_CHAR_LENGTH>();
  const T C_rho = params.get<parameters::C_RHO>();
  const T viscosityH2O = params.get<parameters::NU_LIQUID>();
  const T rho_l = params.get<parameters::RHO_LIQUID>();
  const T rho_g = params.get<parameters::RHO_VAPOR>();
  const T sigma = params.get<parameters::SURFACE_TENSION>();
  const T w = params.get<parameters::INTERFACE_WIDTH>();
  const T theta = params.get<parameters::THETA>();

  sLatticeNS.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,NSDESCRIPTOR>>(
    int   {N},          // resolution
    (T)   tau_l,        // lattice relaxation time
    (T)   L_char,       // charPhysLength: reference length of simulation geometry
    (T)   0,            // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   viscosityH2O, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   C_rho         // physDensity: physical density in __kg / m^3__
  );

  const auto& converter = sLatticeNS.getUnitConverter();
  converter.print();

  sLatticeCH.setUnitConverter(converter);

  // lattice dynamics for Navier-Stokes equation, only needed in fluid domain with mn 1
  dynamics::set<NSBulkDynamics>(sLatticeNS, geometry.getMaterialIndicator({1}));

  // lattice dynamics for Cahn-Hilliard equation, only needed in fluid domain with mn 1
  dynamics::set<CHBulkDynamics>(sLatticeCH, geometry.getMaterialIndicator({1}));

  // use interpolated bounce back conditions for the distribution function
  std::vector<T> origin = { -1.5, 0.5 };
  std::vector<T> extend = { extent[0]+2, extent[1]-2 };
  IndicatorCuboid2D<T> rectangle( extend, origin ); // rectangle spans from (0,0.5) to (Nx-1,Ny-1.5), in the x direction we add padding for the periodic boundaries
  setBouzidiBoundary( sLatticeNS, geometry, 2, rectangle );
  setBouzidiWellBalanced( sLatticeCH, geometry, 2, rectangle );

  auto bulk = geometry.getMaterialIndicator( 1 );

  // postprocessor to update rhowetting
  sLatticeCH.addPostProcessor<stage::PostStream>(bulk,meta::id<RhoWettingStatistics>());

  auto& coupling = myCase.setCouplingOperator(
    "Coupling",
    Coupling{},
    names::NavierStokes{}, sLatticeNS,
    names::Component1{}, sLatticeCH);
  coupling.restrictTo(geometry.getMaterialIndicator({1}));
  // give coupling the values of tau_l, tau_g, rho_l and rho_g so that it can update the viscosity and density
  coupling.setParameter<MixtureRules::TAUS>({tau_g,tau_l});
  coupling.setParameter<MixtureRules::RHOS>({rho_g,rho_l});

  // postprocessor to calculate the chemical potential
  sLatticeCH.addPostProcessor<stage::ChemPotCalc>(meta::id<ChemPotentialPhaseFieldProcessor>());

  // parameters needed by models
  sLatticeNS.setParameter<descriptors::OMEGA>( 1./tau_l );     // collision rate for Navier-Stokes lattice
  sLatticeCH.setParameter<descriptors::OMEGA>( 1./tau_mobil ); // collision rate for Cahn-Hilliard lattice
  sLatticeCH.setParameter<descriptors::THETA>( (M_PI-theta*M_PI/180.) );  // contact angle
  sLatticeCH.setParameter<descriptors::INTERFACE_WIDTH>( w );  // diffuse interface width
  sLatticeCH.setParameter<descriptors::SCALAR>( sigma );       // surface tension

  // define fields that must be communicated by mpi (needed for finite difference gradients)
  {
    auto& communicator = sLatticeCH.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(2);
    communicator.requestField<STATISTIC>();
    communicator.requestField<PHIWETTING>();
    communicator.requestField<CHEM_POTENTIAL>();
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& sLatticeNS = myCase.getLattice(NavierStokes{});
  auto& sLatticeCH = myCase.getLattice(Component1{});

  const Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const int diameter = params.get<parameters::RESOLUTION>();
  const T tau_l = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T tau_g = params.get<parameters::LATTICE_RELAXATION_TIME_2>();
  const T rho_l = params.get<parameters::RHO_LIQUID>();
  const T rho_g = params.get<parameters::RHO_VAPOR>();
  const T w = params.get<parameters::INTERFACE_WIDTH>();

  // set velocity to zero everywhere initially
  Vector<T,2> u(0., 0.);
  AnalyticalConst2D<T,T> zeroVelocity( u );

  // analytical constants to use in initialisation
  AnalyticalConst2D<T,T> one ( 1. );
  AnalyticalConst2D<T,T> two ( 2. );
  AnalyticalConst2D<T,T> zero ( 0. );
  AnalyticalConst2D<T,T> rhov ( rho_g );
  AnalyticalConst2D<T,T> rhol ( rho_l );
  AnalyticalConst2D<T,T> tauv ( tau_g );
  AnalyticalConst2D<T,T> taul ( tau_l );

  // regions with different material ids
  auto bulk = geometry.getMaterialIndicator( 1 );
  auto walls = geometry.getMaterialIndicator( 2 );
  auto all = geometry.getMaterialIndicator( {0,1,2} );

  // initialise phi with a circle with centre at (Nx/2,1) and radius (diameter/2), smoothed by w/2
  SmoothIndicatorFactoredCircle2D<T,T> circle( {extent[0]/2., 1.}, diameter/2., w/2, 0, {0,0}, 0, 1. );
  AnalyticalIdentity2D<T,T> phi( one - circle );

  // initial values for interpolated rho and tau across interfaces
  AnalyticalIdentity2D<T,T> rho( rhov + (rhol-rhov)*phi );
  AnalyticalIdentity2D<T,T> tau( tauv + (taul-tauv)*phi );

  // initial (hydrodynamic) pressure
  AnalyticalIdentity2D<T,T> pressure( zero );

  // initialise distribution functions in equilibrium with initial values of phi, velocity and pressure
  sLatticeCH.defineRhoU( all, phi, zeroVelocity );
  sLatticeCH.iniEquilibrium( all, phi, zeroVelocity );
  sLatticeNS.defineRhoU( all, pressure, zeroVelocity );
  sLatticeNS.iniEquilibrium( all, pressure, zeroVelocity );

  // set the initial values for the fields
  fields::set<descriptors::RHO>(sLatticeNS, all, rho );                     // density
  fields::set<descriptors::TAU_EFF>(sLatticeNS, geometry.getMaterialIndicator(1), tau ); // relaxation time for navier-stokes lattice
  fields::set<descriptors::SOURCE>(sLatticeCH, geometry.getMaterialIndicator(1), zero ); // source term for Cahn-Hilliard lattice
  fields::set<descriptors::PHIWETTING>(sLatticeCH, geometry.getMaterialIndicator(1), phi );      // fluid concentration used by wetting boundaries
  fields::set<descriptors::CHEM_POTENTIAL>(sLatticeCH, geometry.getMaterialIndicator(1), zero ); // chemical potential
  fields::set<descriptors::BOUNDARY>(sLatticeCH, walls, two );              // boundary field used by wetting boundaries
  fields::set<descriptors::BOUNDARY>(sLatticeCH, bulk, zero );              // boundary field used by wetting boundaries

  // some initial steps
  sLatticeCH.executePostProcessors(stage::PreCoupling());
  sLatticeCH.getCommunicator(stage::PreCoupling()).communicate();
  sLatticeCH.executePostProcessors(stage::ChemPotCalc());
  sLatticeCH.getCommunicator(stage::PreCoupling()).communicate();
  sLatticeNS.initialize();
  sLatticeCH.initialize();
  sLatticeCH.iniEquilibrium( all, phi, zeroVelocity );
  sLatticeCH.getCommunicator(stage::PreCoupling()).communicate();
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{ }

// for statistic output and vtk saving
template <typename T>
T getResults( MyCase& myCase,
              util::Timer<MyCase::value_t>& timer,
              std::size_t iT )
{
  OstreamManager clout( std::cout,"getResults" );
  auto& params = myCase.getParameters();

  auto& sLatticeNS = myCase.getLattice(NavierStokes{});
  auto& sLatticeCH = myCase.getLattice(Component1{});

  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using CHDESCRIPTOR = MyCase::descriptor_t_of<Component1>;
  const auto& converter = sLatticeNS.getUnitConverter();

  const Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const int statIter = params.get<parameters::LATTICE_STAT_ITER_T>();
  const int saveIter = params.get<parameters::LATTICE_VTK_ITER_T>();
  const T theta = params.get<parameters::THETA>();

  SuperVTMwriter2D<T> vtmWriter( "contactAngle2d" );
  if ( iT==0 ) {
    SuperLatticeCuboid2D cuboid( sLatticeNS );
    SuperLatticeRank2D rank( sLatticeNS );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  if ( iT%statIter==0 ) {
    timer.update( iT );
    timer.printStep();
    sLatticeNS.getStatistics().print( iT, converter.getPhysTime(iT) );
    sLatticeCH.getStatistics().print( iT, converter.getPhysTime(iT) );
  }

  T contact_angle = 0;
  if ( iT%saveIter==0 ) {

    SuperLatticeDensity2D<T, NSDESCRIPTOR> p_hydro( sLatticeNS );
    p_hydro.getName() = "p_hydro";

    SuperLatticeField2D<T, CHDESCRIPTOR, STATISTIC> phi( sLatticeCH );
    phi.getName() = "phi";

    SuperLatticeField2D<T, NSDESCRIPTOR, RHO> rho_L( sLatticeNS );
    AnalyticalConst2D<T,T> ConversionDensity_( converter.getConversionFactorDensity() );
    SuperLatticeFfromAnalyticalF2D<T, NSDESCRIPTOR> ConversionDensity(ConversionDensity_, sLatticeNS);
    SuperIdentity2D<T,T> rho ( rho_L * ConversionDensity );
    rho.getName() = "rho";

    SuperLatticeVelocity2D<T, NSDESCRIPTOR> velocity( sLatticeNS );
    velocity.getName() = "u";

    vtmWriter.addFunctor( p_hydro );
    vtmWriter.addFunctor( phi );
    vtmWriter.addFunctor( rho );
    vtmWriter.addFunctor( velocity );
    vtmWriter.write( iT );

    // function to interpolate phi at a given position
    AnalyticalFfromSuperF2D<T,T> interpolPhi( phi, true, 1 );
    T point1 = 0, point2 = 0, point3 = 0;

    // contact angle fitting, will be done by finding three points on the circlular droplet

    // point 1, x = ? y = 2
    T pos[2] = {0., 2.};
    for (int ix=0; ix<0.5*extent[0]; ix++) {
      T phi1, phi2;
      pos[0] = ix;
      interpolPhi( &phi1, pos );
      // if phi at ix is greater than 0.5, we have passed the interface of the liquid phase
      if (phi1 < 0.5) {
        pos[0] = (ix-1); // check the value at the previous point
        interpolPhi( &phi2, pos );
        point1 = ix - 1 + (0.5-phi2)/(phi1-phi2); // interpolate between these values
        break;
      }
    }

    // point 2, x = ? y = 2
    for (int ix=0.5*extent[0]; ix<extent[0]; ix++) {
      T phi1, phi2;
      pos[0] = ix;
      interpolPhi( &phi1, pos );
      if (phi1 > 0.5) {
        pos[0] = (ix-1);
        interpolPhi( &phi2, pos );
        point2 = ix - 1 + (0.5-phi2)/(phi1-phi2);
        break;
      }
    }

    // point 3, x = Nx/2 y = ?
    pos[0] = 0.5*(extent[0]);
    for (int iy=3; iy<extent[1]; iy++) {
      T phi1, phi2;
      pos[1] = iy;
      interpolPhi( &phi1, pos );
      if (phi1 > 0.5) {
        pos[0] = (iy-1);
        interpolPhi( &phi2, pos );
        point3 = iy - 1 + (0.5-phi2)/(phi1-phi2);
        break;
      }
    }

    T x1=point1;
    T y1=2;
    T x2=point2;
    T y2=2;
    T x3=0.5*(extent[0]);
    T y3=point3;

    // estimate centre and radius from three points, we must solve three simulatneous equations
    T s1 = x1*x1 + y1*y1;
    T s2 = x2*x2 + y1*y1;
    T s3 = x3*x3 + y3*y3;
    T M11 = x1*y2 + x2*y3 + x3*y1 - (x2*y1 + x3*y2 + x1*y3);
    T M12 = s1*y2 + s2*y3 + s3*y1 - (s2*y1 + s3*y2 + s1*y3);
    T M13 = s1*x2 + s2*x3 + s3*x1 - (s2*x1 + s3*x2 + s1*x3);
    T xc =  0.5*M12/M11;
    T yc = -0.5*M13/M11;
    T r = sqrt(pow(x2 - xc,2) + pow(y2 - yc,2));

    contact_angle = util::acos(-(yc-0.5) / r) * 180 / M_PI;

    clout << "Numerical Contact angle: " << contact_angle << std::endl;
    clout << "Analytical contact angle: " << theta <<  std::endl;
  }
  return contact_angle;
}

// run the simulation
void simulate( MyCase& myCase)
{
  OstreamManager clout( std::cout,"simulate" );
  using T = MyCase::value_t;
  auto& params = myCase.getParameters();

  auto& sLatticeNS = myCase.getLattice(NavierStokes{});
  auto& sLatticeCH = myCase.getLattice(Component1{});

  const std::size_t iTmax = params.get<parameters::LATTICE_MAX_T>();
  const int saveIter = params.get<parameters::LATTICE_VTK_ITER_T>();

  const auto& converter = sLatticeNS.getUnitConverter();

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  std::ofstream outfile;
  outfile.open ("contactAngleVsTime.dat");
  T old_angle = 1.;

  for (std::size_t iT=0; iT <= iTmax; ++iT) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    sLatticeNS.collideAndStream();
    sLatticeCH.collideAndStream();

    sLatticeCH.getCommunicator(stage::PreCoupling()).communicate();
    sLatticeCH.executePostProcessors(stage::PreCoupling());
    sLatticeCH.getCommunicator(stage::PreCoupling()).communicate();
    sLatticeCH.executePostProcessors(stage::ChemPotCalc());
    sLatticeCH.getCommunicator(stage::PreCoupling()).communicate();

    myCase.getOperator("Coupling").apply();

    /// === Step 8.3: Computation and Output of the Results ===
    T angle = getResults<T>( myCase, timer, iT );
    if ( iT%saveIter == 0 ) {
      outfile << iT*converter.getPhysDeltaT() << "," ;
      outfile << angle << "\n";
      if ( fabs(angle-old_angle)/fabs(old_angle) < 5e-5 ) {
        clout << "contact angle converged..." << std::endl;
        break;
      }
      old_angle = angle;
    }
    if ( std::isnan( sLatticeNS.getStatistics().getAverageEnergy() ) ) {
      break;
    }
  }
  outfile.close();
  timer.stop();
  timer.printSummary();
}

int main( int argc, char *argv[] )
{
  initialize( &argc, &argv );

  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<DOMAIN_EXTENT             >({160, 100}); // domain size [lattice units]
    myCaseParameters.set<RESOLUTION                >(70);         // diameter of droplet [lattice units]
    myCaseParameters.set<LATTICE_RELAXATION_TIME   >(1.);         // tau liquid [lattice units]
    myCaseParameters.set<LATTICE_RELAXATION_TIME_2 >(0.52);     // tau gas [lattice units]
    myCaseParameters.set<LATTICE_RELAXATION_TIME_PF>(1.);      // tau mobility [lattice units]
    myCaseParameters.set<LATTICE_MAX_T             >(500000);     // max iterations [lattice units]
    myCaseParameters.set<LATTICE_VTK_ITER_T        >(1000);       // vtk iterations [lattice units]
    myCaseParameters.set<LATTICE_STAT_ITER_T       >(1000);       // statistics iterations [lattice units]
    myCaseParameters.set<PHYS_CHAR_LENGTH          >(70e-6);        // charPhysLength [physical units]
    myCaseParameters.set<C_RHO                     >(1000.);        // conversion factor density [physical units]
    myCaseParameters.set<NU_LIQUID                 >(9e-7);         // physViscosity liquid [physical units]
    myCaseParameters.set<RHO_LIQUID                >(1.);           // lattice density liquid [lattice units]
    myCaseParameters.set<RHO_VAPOR                 >(1.);           // lattice density gas [lattice units]
    myCaseParameters.set<SURFACE_TENSION           >(0.01);         // lattice surface tension [lattice units]
    myCaseParameters.set<parameters::INTERFACE_WIDTH>(4.);// interface thickness [lattice units]
    myCaseParameters.set<parameters::THETA >(100.);       // contact angle in degrees [deg]
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
