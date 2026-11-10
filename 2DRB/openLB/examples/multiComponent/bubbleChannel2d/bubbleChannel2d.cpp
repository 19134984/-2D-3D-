/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Tim Bingert
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

/* bubbleChannel2d.cpp
 * This example shows an advecting bubble of air in a channel of
 * water. The bubble has the size of 40 micrometers and has the
 * exact surface tension of an air-water mixture. This high surface
 * tension creates the necessity for a trick at the outlet, using
 * a fringe zone. This setup achieves a nice shape of the bubble
 * at the artificial outlet, which is roughly five bubble diameters
 * before the actual domain outlet.
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;

using MyCase = Case<
  NavierStokes, Lattice<double, D2Q9<RHO,NABLARHO,FORCE,EXTERNAL_FORCE,TAU_EFF,STATISTIC,SCALAR>>,
  Component1,  Lattice<double, D2Q9<CONV_POPS,FORCE,SOURCE,VELOCITY,OLD_PHIU,STATISTIC>>
>;

using NSBulkDynamics = MultiPhaseIncompressibleTRTdynamics<MyCase::value_t,MyCase::descriptor_t>;
using ACBulkDynamics = ConservativePhaseFieldBGKdynamics<MyCase::value_t,MyCase::descriptor_t_of<Component1>>;

namespace olb::parameters {

struct AC_RELAXATION_TIME    : public descriptors::FIELD_BASE<1> { }; // AC lattice relaxation time
struct C_RHO                 : public descriptors::FIELD_BASE<1> { }; // Density conversion factor
struct DOMAIN_EXTENT_LATTICE : public descriptors::FIELD_BASE<0,1>{ }; // Domain in lattice units

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extent, origin);

  const T dx = params.get<parameters::PHYS_CHAR_LENGTH>() / params.get<parameters::RESOLUTION>();
  Mesh<T,MyCase::d> mesh(cuboid, dx, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();
  using T = MyCase::value_t;
  T Nx = params.get<parameters::DOMAIN_EXTENT_LATTICE>()[0];
  T Ny = params.get<parameters::DOMAIN_EXTENT_LATTICE>()[1];
  const T dx = params.get<parameters::PHYS_CHAR_LENGTH>() / params.get<parameters::RESOLUTION>();

  geometry.rename(0, 2);
  // bulk, MN=1
  geometry.rename(2, 1, {1, 1});

  // inlet, MN=3
  std::vector<T> origin = {-0.5*dx, dx};
  std::vector<T> extend = {1.0*dx, dx*(Ny-2.)};
  IndicatorCuboid2D<T> inlet( extend, origin );
  geometry.rename(2, 3, inlet);

  // outlet, MN=4
  origin[0] = dx*(Nx-0.5);
  IndicatorCuboid2D<T> outlet( extend, origin );
  geometry.rename( 2, 4, outlet );

  geometry.clean();
  geometry.innerClean();
  geometry.checkForErrors();
  geometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

MyCase::value_t helperConvectiveU(MyCase& myCase, SuperIndicatorF2D<MyCase::value_t>& indicator) {
  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using T = MyCase::value_t;

  auto& latticeNS = myCase.getLattice(NavierStokes{});

  int in[2];
  T uMax[2];
  latticeNS.setProcessingContext(ProcessingContext::Evaluation);
  SuperLatticeVelocity2D<T, NSDESCRIPTOR> u( latticeNS );
  SuperMax2D<T,T> uMaxFinder(u, indicator);
  uMaxFinder(uMax, in);
  return uMax[0];
}

void prepareLattice(MyCase& myCase)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice with conservative phase field model ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& latticeNS = myCase.getLattice(NavierStokes{});
  auto& latticeAC = myCase.getLattice(Component1{});

  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

  const int N       = params.get<parameters::RESOLUTION>();
  const T char_l    = params.get<parameters::PHYS_CHAR_LENGTH>();
  const T rho_g     = params.get<parameters::RHO_VAPOR>();
  const T rho_l     = params.get<parameters::RHO_LIQUID>();
  const T nu_g      = params.get<parameters::NU_VAPOR>();
  const T nu_l      = params.get<parameters::NU_LIQUID>();
  const T tau_l     = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T tau_g     = nu_g/nu_l*(tau_l-0.5)+0.5;
  const T sigma     = params.get<parameters::SURFACE_TENSION>();
  const T w         = params.get<parameters::INTERFACE_WIDTH>();
  const T C_rho     = params.get<parameters::C_RHO>();

  latticeNS.setUnitConverter<MultiPhaseUnitConverterFromRelaxationTime<T,NSDESCRIPTOR>>(
    int   {N},               // resolution
    (T)   tau_l,             // lattice relaxation time
    (T)   rho_l/C_rho,       // lattice density heavier phase
    (T)   char_l,            // charPhysLength: reference length of simulation geometry
    (T)   nu_l,              // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   rho_l              // physDensity: physical density in __kg / m^3__
  );

  const auto& converter = latticeNS.getUnitConverter();
  converter.print();
  const T sigma_Lattice = sigma / converter.getConversionFactorSurfaceTension();
  clout << "Lattice Surface Tension: " << sigma_Lattice << std::endl;

  latticeAC.setUnitConverter(converter);

  // define lattice Dynamics
  dynamics::set<NSBulkDynamics>(latticeNS, geometry, 1);
  dynamics::set<ACBulkDynamics>(latticeAC, geometry, 1);

  auto& coupling = myCase.setCouplingOperator(
    "Coupling",
    LiangPostProcessor{},
    names::NavierStokes{}, latticeNS,
    names::Component1{}, latticeAC
  );
  coupling.restrictTo(geometry.getMaterialIndicator({1,4}));

  auto& velocityCoupling = myCase.setCouplingOperator(
    "VelocityCoupling",
    VelocityCoupling{},
    names::NavierStokes{}, latticeNS,
    names::Component1{}, latticeAC
  );
  velocityCoupling.restrictTo(geometry.getMaterialIndicator({3}));

  coupling.setParameter<LiangPostProcessor::SIGMA>(sigma_Lattice);
  coupling.setParameter<LiangPostProcessor::W>(w);
  coupling.setParameter<LiangPostProcessor::TAUS>({tau_g, tau_l});
  coupling.setParameter<LiangPostProcessor::RHOS>({rho_g/C_rho,
                                                  rho_l/C_rho});
  coupling.setParameter<LiangPostProcessor::SWITCH>(1.);

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ACDESCRIPTOR = MyCase::descriptor_t_of<Component1>;
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();
  auto& latticeNS = myCase.getLattice(NavierStokes{});
  auto& latticeAC = myCase.getLattice(Component1{});
  const auto& converter = latticeNS.getUnitConverter();

  const T Re        = params.get<parameters::REYNOLDS>();
  const T tau_mobil = params.get<parameters::AC_RELAXATION_TIME>();
  const T rho_g     = params.get<parameters::RHO_VAPOR>();
  const T rho_l     = params.get<parameters::RHO_LIQUID>();
  const T nu_g      = params.get<parameters::NU_VAPOR>();
  const T nu_l      = params.get<parameters::NU_LIQUID>();
  const T tau_l     = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T tau_g     = nu_g/nu_l*(tau_l-0.5)+0.5;
  const T sigma     = params.get<parameters::SURFACE_TENSION>();
  const T w         = params.get<parameters::INTERFACE_WIDTH>();
  const T theta     = params.get<parameters::THETA>();
  const T C_rho     = params.get<parameters::C_RHO>();

  // conversion properties
  const T dx = converter.getPhysDeltaX();
  const T dt = converter.getPhysDeltaT();
  const T C_sigma = converter.getConversionFactorSurfaceTension();
  const T C_p = C_sigma/dx;
  const T diameter_lattice = params.get<parameters::PHYS_CHAR_LENGTH>() / dx;
  const T Nx = params.get<parameters::DOMAIN_EXTENT_LATTICE>()[0];
  const T Ny = params.get<parameters::DOMAIN_EXTENT_LATTICE>()[1];
  const T sigma_Lattice = sigma / C_sigma;

    // initial conditions
  AnalyticalConst2D<T,T> two ( 2. );
  AnalyticalConst2D<T,T> one ( 1. );
  AnalyticalConst2D<T,T> rhog ( rho_g/C_rho );
  AnalyticalConst2D<T,T> rhol ( rho_l/C_rho );
  AnalyticalConst2D<T,T> taug ( tau_g );
  AnalyticalConst2D<T,T> taul ( tau_l );
  AnalyticalConst2D<T,T> u0( 0,0 );

  auto bulk = geometry.getMaterialIndicator(1);
  auto wall = geometry.getMaterialIndicator(2);
  auto inlet = geometry.getMaterialIndicator(3);
  auto outlet = geometry.getMaterialIndicator(4);
  auto fluid = geometry.getMaterialIndicator({1,3});
  auto all = geometry.getMaterialIndicator({0,1,2,3,4});

  std::vector<T> pos = {dx*2.*diameter_lattice, dx*Ny/2.};
  CircularInterface2D<T> phi0(pos, dx*diameter_lattice/2., dx*w, 1., true);
  AnalyticalIdentity2D<T,T> rho0( rhog + (rhol-rhog)*phi0 );
  AnalyticalIdentity2D<T,T> tau0( taug + (taul-taug)*phi0 );
  std::shared_ptr<AnalyticalF2D<T,T>> bubblePressure(new LaplacePressure2D<T>( pos, dx*diameter_lattice/2., dx*w, sigma_Lattice*C_sigma/C_p ));
  std::shared_ptr<AnalyticalF2D<T,T>> halfYLPressure(new AnalyticalConst2D<T,T>( 2.*sigma_Lattice/diameter_lattice/2. ));
  std::shared_ptr<AnalyticalF2D<T,T>> ConvertP(new AnalyticalConst2D<T,T>( C_p ));
  auto p0 = bubblePressure + halfYLPressure;
  auto p0_phys = p0*ConvertP;
  SmoothIndicatorFactoredCuboid2D<T,T> fringe( {0., dx*Ny/2.}, dx*2.*(Nx-5.0*diameter_lattice), 0, dx*1.5*diameter_lattice, 0, {0,0}, 0, 1. );

  latticeNS.defineField<descriptors::RHO>( all, rho0 );
  latticeNS.defineField<descriptors::TAU_EFF>( all, tau0 );
  latticeNS.defineField<descriptors::SCALAR>( all, fringe );
  latticeNS.defineField<descriptors::SCALAR>( wall, two );
  latticeAC.defineField<descriptors::OLD_PHIU>( all, u0 );

  // walls
  std::vector<T> origin = {-2.*dx, dx*0.5};
  std::vector<T> extend = {dx*(Nx+4), dx*(Ny-2)};
  IndicatorCuboid2D<T> wallLocation( extend, origin );
  setBouzidiBoundary(latticeNS, geometry, 2, wallLocation);
  setBouzidiPhaseField<IsoPhaseFieldCurvedWall2D>(latticeAC, geometry, 2, wallLocation);

  // inlet and outlet
  boundary::set<boundary::IncompressibleZouHeVelocity<T,NSDESCRIPTOR,NSBulkDynamics,momenta::IncRegularizedVelocityBoundaryTuple>>(latticeNS, inlet);
  boundary::set<boundary::RegularizedTemperature>(latticeAC, inlet);
  boundary::set<boundary::IncompressibleZouHePressure<T,NSDESCRIPTOR,NSBulkDynamics,momenta::IncRegularizedPressureBoundaryTuple>>(latticeNS, outlet);
  setConvectivePhaseFieldBoundary<T,ACDESCRIPTOR>(latticeAC, outlet);

  const T maxVelocity_phys = Re/Ny*((tau_l-0.5)/3.)*dx/dt;
  const T radius = T(0.5)*dx*(Ny - 2);
  std::vector<T> axisPoint( 2,T() );
  axisPoint[0] = dx*Nx/2.;
  axisPoint[1] = dx*(Ny-1.)/2.;
  std::vector<T> axisDirection( 2,T() );
  axisDirection[0] = 1;
  axisDirection[1] = 0;
  Poiseuille2D<T> poiseuille( axisPoint, axisDirection, maxVelocity_phys/dx*dt, radius );
  Poiseuille2D<T> poiseuille_phys( axisPoint, axisDirection, maxVelocity_phys, radius );

  momenta::setOrderParameter(latticeAC, all, phi0);
  momenta::setVelocity(latticeAC, all, poiseuille_phys);

  momenta::setIncompressiblePressure(latticeNS, all, *p0_phys);
  momenta::setVelocity(latticeNS, all, poiseuille_phys);
  latticeNS.iniEquilibrium( all, *p0, poiseuille );

  latticeAC.addPostProcessor<stage::PreCoupling>(fluid, meta::id<RhoStatistics>());
  latticeAC.addPostProcessor<stage::InitOutlet>(outlet, meta::id<SetOutletCells2D<1,0>>());

  latticeNS.setParameter<descriptors::OMEGA>( 1./tau_l );
  latticeAC.setParameter<descriptors::OMEGA>( 1./tau_mobil );
  latticeAC.setParameter<descriptors::INTERFACE_WIDTH>( w );
  latticeAC.setParameter<descriptors::THETA>( theta );
  latticeAC.setParameter<descriptors::EPSILON>( 1.5*w );
  latticeNS.setParameter<collision::TRT::MAGIC>( (tau_l-T(0.5))*(1.6-0.5) );

  {
    auto& communicator = latticeAC.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(1);
    communicator.requestField<STATISTIC>();
    communicator.exchangeRequests();
  }

  latticeNS.initialize();
  latticeAC.initialize();

  latticeAC.iniEquilibrium( all, phi0, poiseuille );
  latticeNS.iniEquilibrium( all, *p0, poiseuille );

  latticeNS.setProcessingContext(ProcessingContext::Simulation);
  latticeAC.setProcessingContext(ProcessingContext::Simulation);

  latticeAC.executePostProcessors(stage::InitOutlet());
  latticeAC.getCommunicator(stage::PreCoupling()).communicate();
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{ }

void getResults(
  MyCase& myCase,
  util::Timer<MyCase::value_t>& timer,
  std::size_t iT)
{
  OstreamManager clout( std::cout,"getResults" );

  using T = MyCase::value_t;
  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ACDESCRIPTOR = MyCase::descriptor_t_of<Component1>;
  auto& params = myCase.getParameters();
  auto& latticeNS = myCase.getLattice(NavierStokes{});
  auto& latticeAC = myCase.getLattice(Component1{});
  const auto& converter = latticeNS.getUnitConverter();
  const int statIter = params.get<parameters::LATTICE_STAT_ITER_T>();
  const int vtkIter = params.get<parameters::LATTICE_VTK_ITER_T>();

  SuperVTMwriter2D<T> vtmWriter( "bubbleChannel2d" );
  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, NSDESCRIPTOR> cuboid( latticeNS );
    SuperLatticeRank2D<T, NSDESCRIPTOR> rank( latticeNS );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }
  // Get statistics
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();
    latticeNS.getStatistics().print( iT, converter.getPhysTime(iT), "NSlattice" );
    latticeAC.getStatistics().print( iT, converter.getPhysTime(iT), "AClattice" );
  }

  // Writes the VTK files
  if ( iT%vtkIter==0 ) {
    latticeNS.setProcessingContext(ProcessingContext::Evaluation);
    latticeAC.setProcessingContext(ProcessingContext::Evaluation);
    SuperLatticePhysIncPressure2D<T, NSDESCRIPTOR> p_total( latticeNS, converter );
    p_total.getName() = "p_total";
    SuperLatticePhysField2D<T, NSDESCRIPTOR, RHO> rho(latticeNS, converter.getConversionFactorDensity());
    rho.getName() = "rho";
    SuperLatticePhysVelocity2D<T, NSDESCRIPTOR> velocity( latticeNS, converter );
    velocity.getName() = "u";
    SuperLatticeField2D<T, ACDESCRIPTOR, STATISTIC> phi( latticeAC );
    phi.getName() = "phi";
    SuperLatticeField2D<T, NSDESCRIPTOR, SCALAR> scale( latticeNS );
    scale.getName() = "scale";
    vtmWriter.addFunctor( scale );

    vtmWriter.addFunctor( p_total );
    vtmWriter.addFunctor( rho );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( phi );
    vtmWriter.write( iT );
  }
}

void simulate(MyCase& myCase) {
  OstreamManager clout( std::cout,"simulate" );

  using T = MyCase::value_t;

  auto& latticeNS = myCase.getLattice(NavierStokes{});
  auto& latticeAC = myCase.getLattice(Component1{});
  auto& geometry = myCase.getGeometry();
  const auto& converter = latticeNS.getUnitConverter();
  auto& params = myCase.getParameters();

  const std::size_t maxIter = params.get<parameters::LATTICE_MAX_T>();

  const T dx = converter.getPhysDeltaX();
  const T Nx = params.get<parameters::DOMAIN_EXTENT_LATTICE>()[0];
  const T Ny = params.get<parameters::DOMAIN_EXTENT_LATTICE>()[1];

  Vector<T,2> origin2 = {dx*(Nx-2.), dx};
  Vector<T,2> extend2 = {1.5*dx, dx*(Ny-2.)};
  IndicatorCuboid2D<T> beforeOutlet_( extend2, origin2 );
  SuperIndicatorFfromIndicatorF2D<T> beforeOutlet(beforeOutlet_, geometry);

  util::Timer<T> timer(maxIter, myCase.getGeometry().getStatistics().getNvoxel());
  clout << "Starting simulation ";
  clout << "with " << maxIter << " timesteps " << std::endl;
  timer.start();

  for (std::size_t iT=0; iT<=maxIter; ++iT ) {
    // compute velocity for convective outlet
    if ( iT%500==0 ) {
      T uMax = helperConvectiveU(myCase, beforeOutlet);
      latticeAC.setParameter<descriptors::MAX_VELOCITY>( uMax );
    }

    // Collide and stream (and coupling) execution
    latticeNS.collideAndStream();
    latticeAC.collideAndStream();

    latticeAC.getCommunicator(stage::PreCoupling()).communicate();
    latticeAC.executePostProcessors(stage::PreCoupling());
    latticeAC.getCommunicator(stage::PreCoupling()).communicate();

    myCase.getOperator("Coupling").apply();
    myCase.getOperator("VelocityCoupling").apply();

    // Computation and output of the results
    getResults( myCase, timer, iT );
    if ( std::isnan( latticeNS.getStatistics().getAverageEnergy() ) ) {
      break;
    }

  }
  timer.stop();
  timer.printSummary();
}

int main( int argc, char *argv[] )
{
  /// === Step 1: Initialize olb ===
  initialize( &argc, &argv );

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<DOMAIN_EXTENT          >({64e-5, 16e-5}); // domain size [physical units]
    myCaseParameters.set<RESOLUTION             >(40);             // resolution y [lattice units]
    myCaseParameters.set<PHYS_CHAR_LENGTH       >(4e-5);           // characteristic length (bubble radius) [phys units]
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.53);      // tau liquid [lattice units]
    myCaseParameters.set<AC_RELAXATION_TIME     >(0.8);            // tau mobility [lattice units]
    myCaseParameters.set<REYNOLDS               >(100.);           // Reynolds number []
    myCaseParameters.set<RHO_VAPOR              >(1.2);            // physDensity gas/vapor [physical units]
    myCaseParameters.set<RHO_LIQUID             >(1000.);          // physDensity liquid [physical units]
    myCaseParameters.set<NU_VAPOR               >(117e-7);         // physViscosity gas/vapor [physical units]
    myCaseParameters.set<NU_LIQUID              >(9e-7);           // physViscosity liquid [physical units]
    myCaseParameters.set<SURFACE_TENSION        >(0.072);          // surface tension [physical units]
    myCaseParameters.set<parameters::THETA >(M_PI*40./180.);  // contact angle (<90 is wetting) [radians]
    myCaseParameters.set<parameters::INTERFACE_WIDTH>(6.);    // interface width [lattice units]
    myCaseParameters.set<C_RHO              >(100.);           // conversion factor density [physical units]
    myCaseParameters.set<LATTICE_VTK_ITER_T >(2000.);
    myCaseParameters.set<LATTICE_STAT_ITER_T>(2000.);
  }
  {
    using T = MyCase::value_t;
    myCaseParameters.set<parameters::LATTICE_MAX_T>([] (MyCase::ParametersD& p) {
      const T char_l = p.get<parameters::PHYS_CHAR_LENGTH>();
      const T N = p.get<parameters::RESOLUTION>();
      const T Nx = p.get<parameters::DOMAIN_EXTENT>()[0] / char_l * N;
      const T Ny = p.get<parameters::DOMAIN_EXTENT>()[1] / char_l * N;
      const T Re = p.get<parameters::REYNOLDS>();
      const T tau_l = p.get<parameters::LATTICE_RELAXATION_TIME>();
      const T maxIter = Nx/(Re/Ny*((tau_l-0.5)/3.))*1.1;
      return maxIter;
    });
    myCaseParameters.set<parameters::DOMAIN_EXTENT_LATTICE>([] (MyCase::ParametersD& p) {
      const T char_l = p.get<parameters::PHYS_CHAR_LENGTH>();
      const T N = p.get<parameters::RESOLUTION>();
      const T Nx = p.get<parameters::DOMAIN_EXTENT>()[0] / char_l * N;
      const T Ny = p.get<parameters::DOMAIN_EXTENT>()[1] / char_l * N;
      return Vector<T,2> {Nx, Ny};
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
