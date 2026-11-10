/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Tim Bingert
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

/*stratifiedHorizontal3d.cpp
 * Simulating stratified water-air flow of horizontal pipe in 3d,
 * using the Incompressible Conservative Phase Field model with an
 * outlet stabilizing fringe zone.
 */

//#define DISABLE_CSE
#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<float, D3Q19<RHO,NABLARHO,LAPLACERHO,FORCE,TAU_EFF,SCALAR,MOMENTUM>>,
  Component1,  Lattice<float, D3Q7<FORCE,VELOCITY,OLD_PHIU,STATISTIC,SCALAR>>
>;

using NSBulkDynamics = CorrectedMultiPhaseIncompressibleSmagorinskyRLBdynamics<MyCase::value_t,MyCase::descriptor_t>;
using ACBulkDynamics = ConservativePhaseFieldRLBdynamics<MyCase::value_t,MyCase::descriptor_t_of<Component1>>;
using MixtureRules = LinearTauViscosity;
using Coupling = CorrectedLiangPostProcessor<MixtureRules>;

namespace olb::parameters {

struct C_RHO : public descriptors::FIELD_BASE<1> { };
struct MAX_PHYS_START_T : public descriptors::FIELD_BASE<1> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  Vector extent = params.get<parameters::DOMAIN_EXTENT>();

  T radius = 0.5*extent[1];
  T length = extent[0];
  Vector<T, 3> center0(0, radius, radius);
  Vector<T, 3> center1(length, radius, radius);
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  const T dx = extent[1] / params.get<parameters::RESOLUTION>();
  IndicatorLayer3D<T> extendedDomain(pipe, dx);

  Mesh<T,MyCase::d> mesh(extendedDomain, dx, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry( MyCase& myCase )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& params = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T dx = extent[1] / params.get<parameters::RESOLUTION>();

  T radius = 0.5*extent[1];
  T length = extent[0];

  geometry.rename( 0,2 );
  // bulk, MN=1
  Vector<T, 3> center0(0*dx, radius, radius);
  Vector<T, 3> center1(length, radius, radius);
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  geometry.rename(2, 1, pipe);

  Vector<T, 3> origin(0, radius, radius);
  Vector<T, 3> extend = origin;

  // gas inlet, MN=3
  origin[0] = -dx;
  extend[0] = dx;
  IndicatorCylinder3D<T> inflow(origin, extend, radius);
  geometry.rename(2, 3, 1, inflow);

  // outlet, MN=4
  origin[0] = length - dx;
  extend[0] = length + dx;
  IndicatorCylinder3D<T> outflow(origin, extend, radius);
  geometry.rename(2, 4, 1, outflow);

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
  auto& sLatticeAC = myCase.getLattice(Component1{});

  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ACDESCRIPTOR = MyCase::descriptor_t_of<Component1>;

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const int N = params.get<parameters::RESOLUTION>();
  const T rho_g = params.get<parameters::RHO_VAPOR>();
  const T rho_l = params.get<parameters::RHO_LIQUID>();
  const T nu_g = params.get<parameters::NU_VAPOR>();
  const T nu_l = params.get<parameters::NU_LIQUID>();
  const T tau_l = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T tau_g = nu_g/nu_l*(tau_l-0.5)+0.5;
  const T sigma = params.get<parameters::SURFACE_TENSION>();
  const T w = params.get<parameters::INTERFACE_WIDTH>();
  const T theta = params.get<parameters::THETA>();
  const T C_rho = params.get<parameters::C_RHO>();
  Vector g = params.get<parameters::GRAVITY>();

  sLatticeNS.setUnitConverter<MultiPhaseUnitConverterFromRelaxationTime<T,NSDESCRIPTOR>>(
    int   {N},           // resolution
    (T)   tau_l,         // lattice relaxation time
    (T)   rho_l/C_rho,   // lattice density heavier phase
    (T)   extent[1],     // charPhysLength: reference length of simulation geometry
    (T)   nu_l,          // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   rho_l          // physDensity: physical density in __kg / m^3__
  );
  const auto& converter = sLatticeNS.getUnitConverter();
  converter.print();
  sLatticeAC.setUnitConverter(converter);

  // conversion properties
  T dx = converter.getPhysDeltaX();
  T dt = converter.getPhysDeltaT();
  T C_sigma = converter.getConversionFactorSurfaceTension();
  T sigma_L = sigma / C_sigma;
  clout << "lattice surface tension: " << sigma_L << std::endl;

  // define lattice Dynamics
  dynamics::set<NSBulkDynamics>(sLatticeNS, geometry, 1);
  dynamics::set<ACBulkDynamics>(sLatticeAC, geometry, 1);

  auto zero = geometry.getMaterialIndicator(0);
  auto bulk = geometry.getMaterialIndicator(1);
  auto wall = geometry.getMaterialIndicator(2);
  auto inlet = geometry.getMaterialIndicator(3);
  auto outlet = geometry.getMaterialIndicator(4);
  auto fluid = geometry.getMaterialIndicator({1,3});
  auto all = geometry.getMaterialIndicator({0,1,2,3,4});

  // walls
  T radius = 0.5*extent[1];
  T length = extent[0];
  Vector center0(-0.5*dx, radius, radius);
  Vector center1(length+0.5*dx, radius, radius);
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  setBouzidiBoundary(sLatticeNS, geometry, 2, pipe);
  setBouzidiPhaseField<IsoPhaseFieldCurvedWall3D>(sLatticeAC, geometry, 2, pipe);

  // inlet_gas, inlet_liquid and outlet
  boundary::set<boundary::
                IncompressibleZouHeVelocity<T,NSDESCRIPTOR,
                                            NSBulkDynamics,
                                            momenta::IncRegularizedVelocityBoundaryTuple>>
                (sLatticeNS, inlet);
  boundary::set<boundary::RegularizedPhaseFieldDirichlet>(sLatticeAC, inlet);
  boundary::set<boundary::
                IncompressibleLocalPressure<T,NSDESCRIPTOR,
                                            NSBulkDynamics,
                                            momenta::IncRegularizedPressureBoundaryTuple>>
                (sLatticeNS, outlet);
  setNeumannPhaseFieldBoundary<T,ACDESCRIPTOR>(sLatticeAC, outlet);

  sLatticeAC.addPostProcessor<stage::InitOutlet>(outlet,meta::id<SetOutletGhosts3D<1,0,0>>());
  sLatticeAC.addPostProcessor<stage::PreCoupling>(fluid,meta::id<PhiStatistics>());
  sLatticeNS.addPostProcessor<stage::PreCollide>(outlet,meta::id<MomentumExtrapolatePostProcessor3D<T,NSDESCRIPTOR,0,1>>());

  auto& coupling = myCase.setCouplingOperator(
      "Coupling",
      Coupling{},
      names::NavierStokes{}, sLatticeNS,
      names::Component1{}, sLatticeAC);
  coupling.restrictTo(geometry.getMaterialIndicator({1, 3, 4}));

  phasefield::init<Coupling,MixtureRules>(
    sLatticeNS, sLatticeAC,  // both lattices
    coupling, all,           // phase field coupling, all material numbers
    wall, zero,              // wall and zero material numbers
    sigma_L, w, theta,       // phase field lattice parameters
    tau_g, tau_l,            // relaxation times for mixture rule
    rho_g/C_rho, rho_l/C_rho // lattice densities for mixture rule
  );
  coupling.setParameter<Coupling::GAMMA>( T(-22.5) );
  coupling.setParameter<Coupling::GRAVITY>(g*dt*dt/dx);
  coupling.setParameter<Coupling::BODY_FORCE>({T(0),T(0),T(0)});

  sLatticeNS.setParameter<descriptors::OMEGA>( T(1)/tau_l );
  const T Pe = 6.0;
  const T tau_mobil = T(0.5) + (dx*dx*w)/(dt*Pe/T(3));
  clout << "lattice relaxation time: " << tau_mobil << std::endl;
  sLatticeAC.setParameter<descriptors::OMEGA>( T(1)/tau_mobil );
  sLatticeNS.setParameter<collision::LES::SMAGORINSKY>( T(0.15) );

  PlanarInterface3D<T> fringe( {T(0), extent[1]/T(2), extent[2]/T(2)}, {extent[0]-T(60)*dx, 0, 0}, dx*T(90), T(1) );
  phasefield::setFringeLattice<collision::LES::SMAGORINSKY_F>( sLatticeNS, all, wall, zero, fringe, T(10) );
  phasefield::setFringeLattice<descriptors::INTERFACE_WIDTH_F,descriptors::THETA_F>( sLatticeAC, all, wall, zero, fringe, T(2)*w, M_PI/T(2) );
  phasefield::setFringeCoupling<Coupling::SIGMA_F,Coupling::W_F,Coupling::TAU_F>( coupling, T(0), T(2)*w, T(1) );

  {
    auto& communicator = sLatticeAC.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(2);
    communicator.requestField<STATISTIC>();
    communicator.requestField<SCALAR>();
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  OstreamManager clout(std::cout, "setInitialValues");
  clout << "Set initial values ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& sLatticeNS = myCase.getLattice(NavierStokes{});
  auto& sLatticeAC = myCase.getLattice(Component1{});

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const int N = params.get<parameters::RESOLUTION>();
  const T dx = extent[1] / N;
  const T rho_l = params.get<parameters::RHO_LIQUID>();
  const T rho_g = params.get<parameters::RHO_VAPOR>();
  const T C_rho = params.get<parameters::C_RHO>();
  const T nu_g = params.get<parameters::NU_VAPOR>();
  const T nu_l = params.get<parameters::NU_LIQUID>();
  const T tau_l = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T tau_g = nu_g/nu_l*(tau_l-0.5)+0.5;
  const T w = params.get<parameters::INTERFACE_WIDTH>();

  // initial conditions
  AnalyticalConst3D<T,T> one ( T(1) );
  AnalyticalConst3D<T,T> zero ( T(0) );
  AnalyticalConst3D<T,T> rhol ( rho_l/C_rho );
  AnalyticalConst3D<T,T> rhog ( rho_g/C_rho );
  AnalyticalConst3D<T,T> taul ( tau_l );
  AnalyticalConst3D<T,T> taug ( tau_g );
  AnalyticalConst3D<T,T> u0( 0,0,0 );
  AnalyticalConst3D<T,T> stats0( 0,0 );

  auto all = geometry.getMaterialIndicator({0,1,2,3,4});
  T radius = 0.5*extent[1];
  T waterAngle = M_PI*T(3)/T(8);
  PlanarInterface3D<T> phi0( {T(0.5)*dx, T(0.5)*dx, radius},
                             {extent[1]/T(2)*util::sin<T>(waterAngle)*util::sin<T>(M_PI/T(2)-waterAngle),
                              extent[1]/T(2)*util::sin<T>(waterAngle)*util::cos<T>(M_PI/T(2)-waterAngle),
                              T(0)},
                              dx*w*util::sin<T>(waterAngle), T(1) );
  AnalyticalIdentity3D<T,T> rho0( rhog + (rhol-rhog)*phi0 );
  AnalyticalIdentity3D<T,T> tau0( taug + (taul-taug)*phi0 );

  sLatticeNS.defineField<descriptors::RHO>( all, rho0 );
  sLatticeNS.defineField<descriptors::TAU_EFF>( all, tau0 );
  sLatticeNS.defineField<descriptors::MOMENTUM>( all, u0 );
  sLatticeAC.defineField<descriptors::OLD_PHIU>( all, u0 );
  sLatticeAC.defineField<descriptors::STATISTIC>( geometry.getMaterialIndicator(0), stats0 );

  momenta::setOrderParameter(sLatticeAC, all, phi0);
  momenta::setVelocity(sLatticeAC, all, u0);

  momenta::setIncompressiblePressure(sLatticeNS, all, zero);
  momenta::setVelocity(sLatticeNS, all, u0);
  sLatticeNS.iniEquilibrium( all, zero, u0 );

  sLatticeNS.initialize();
  sLatticeAC.initialize();

  sLatticeAC.iniEquilibrium( all, phi0, u0 );
  sLatticeNS.iniEquilibrium( all, zero, u0 );
  sLatticeNS.setProcessingContext(ProcessingContext::Simulation);
  sLatticeAC.setProcessingContext(ProcessingContext::Simulation);

  sLatticeAC.executePostProcessors(stage::InitOutlet());
  sLatticeAC.getCommunicator(stage::PreCoupling()).communicate();

  clout << "Set initial values ... OK" << std::endl;
}

void setTemporalValues(MyCase& myCase, std::size_t iT)
{
  OstreamManager clout(std::cout, "setTemporalValues");

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();
  const T maxStartT = params.get<parameters::MAX_PHYS_START_T>();
  const int N = params.get<parameters::RESOLUTION>();
  const T Re = params.get<parameters::REYNOLDS>();
  const T nu_g = params.get<parameters::NU_VAPOR>();
  const T nu_l = params.get<parameters::NU_LIQUID>();
  const T tau_l = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T tau_g = nu_g/nu_l*(tau_l-0.5)+0.5;

  auto& sLatticeNS = myCase.getLattice(NavierStokes{});
  auto& sLatticeAC = myCase.getLattice(Component1{});

  const auto& converter = sLatticeNS.getUnitConverter();
  auto inlet = geometry.getMaterialIndicator(3);
  std::size_t iTmaxStart = converter.getLatticeTime(maxStartT);
  int iTupdate = converter.getLatticeTime(maxStartT/200);

  // Ramp inlet velocities
  if (iT % iTupdate == 0 && iT <= iTmaxStart) {
    // Smooth start curve, polynomial
    PolynomialStartScale<T, T> StartScale(iTmaxStart, T(1));

    T iTvec[1] = {T(iT)};
    T frac[1]  = {};
    StartScale(frac, iTvec);

    const T u_x_scalar = Re/N*((tau_g-T(0.5))/T(3)) * frac[0];
    AnalyticalConst3D<T,T> u_ramp( u_x_scalar,0,0 );

    sLatticeNS.defineU(inlet, u_ramp);
    sLatticeAC.defineU(inlet, u_ramp);
    sLatticeNS.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

void getResults( MyCase& myCase,
                 util::Timer<MyCase::value_t>& timer,
                 std::size_t iT )
{
  OstreamManager clout(std::cout, "getResults");

  using T = MyCase::value_t;
  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ACDESCRIPTOR = MyCase::descriptor_t_of<Component1>;
  auto& params = myCase.getParameters();
  auto& geometry = myCase.getGeometry();

  auto& sLatticeNS = myCase.getLattice(NavierStokes{});
  auto& sLatticeAC = myCase.getLattice(Component1{});

  const auto& converter = sLatticeNS.getUnitConverter();

  const T writeT = params.get<parameters::PHYS_STAT_ITER_T>();
  Vector extent = params.get<parameters::DOMAIN_EXTENT>();

  if ( iT==0 ) {
    SuperVTMwriter3D<T> vtmWriter( "horizontalGasLiquidPipe3d" );
    vtmWriter.createMasterFile();
  }

  if ( iT%converter.getLatticeTime(T(5)*writeT)==0 ){
    sLatticeAC.setProcessingContext(ProcessingContext::Evaluation);
    SuperLatticeDensity3D<T, ACDESCRIPTOR> phi( sLatticeAC );
    T dx = converter.getPhysDeltaX();
    T pos_yz = 0.5*extent[1];
    int in[3];
    T voidFrac[1];
    Vector<T, 3> center0(2./5.*extent[0], pos_yz, pos_yz);
    Vector<T, 3> center1(3./5.*extent[0], pos_yz, pos_yz);
    IndicatorCylinder3D<T> _sampleZone(center0, center1, pos_yz-dx);
    SuperIndicatorFfromIndicatorF3D<T> sampleZone(_sampleZone, geometry);
    SuperAverage3D<T,T> voidFracFinder(phi, sampleZone);
    voidFracFinder(voidFrac, in);
    clout << "gas void fraction in sample zone: " << T(1) - voidFrac[0] << std::endl;
  }

  if ( iT%converter.getLatticeTime(writeT)==0 ) {
    timer.update( iT );
    timer.printStep();
    sLatticeNS.getStatistics().print( iT, converter.getPhysTime(iT) );
    sLatticeAC.getStatistics().print( iT, converter.getPhysTime(iT) );
  }

  // Writes the VTK files
  if ( iT%converter.getLatticeTime(writeT)==0 || std::isnan( sLatticeNS.getStatistics().getAverageEnergy() ) ) {
    sLatticeNS.setProcessingContext(ProcessingContext::Evaluation);
    sLatticeNS.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtkWriter("horizontalGasLiquidPipe3d");

      SuperLatticePhysIncPressure3D<T, NSDESCRIPTOR> p_total( sLatticeNS, converter );
      p_total.getName() = "p_total";
      SuperLatticePhysField3D<T, NSDESCRIPTOR, RHO> rho(sLatticeNS, converter.getConversionFactorDensity());
      rho.getName() = "rho";
      SuperLatticeField3D<T, NSDESCRIPTOR, TAU_EFF> taueff( sLatticeNS );
      taueff.getName() = "taueff";
      SuperLatticePhysVelocity3D<T, NSDESCRIPTOR> velocity( sLatticeNS, converter );
      velocity.getName() = "u";
      SuperLatticeField3D<T, NSDESCRIPTOR, SCALAR> scale( sLatticeNS );
      scale.getName() = "scale";
      vtkWriter.addFunctor( scale );

      vtkWriter.addFunctor( p_total );
      vtkWriter.addFunctor( rho );
      vtkWriter.addFunctor( velocity );
      vtkWriter.addFunctor( taueff );

      task(vtkWriter, iT);
    });
  }
}

void simulate( MyCase& myCase )
{
  OstreamManager clout( std::cout,"simulate" );

  using T = MyCase::value_t;
  auto& params = myCase.getParameters();
  auto& sLatticeNS = myCase.getLattice(NavierStokes{});
  auto& sLatticeAC = myCase.getLattice(Component1{});

  int iB = 0;
  const std::size_t iTmax = sLatticeNS.getUnitConverter().getLatticeTime(
    params.get<parameters::PHYS_MAX_T>());

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();
  for ( std::size_t iT=0; iT<=iTmax; ++iT ) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    sLatticeNS.collideAndStream();
    sLatticeAC.collideAndStream();

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, timer, iT);

    sLatticeAC.getCommunicator(stage::PreCoupling()).communicate();
    sLatticeAC.executePostProcessors(stage::PreCoupling());
    sLatticeAC.getCommunicator(stage::PreCoupling()).communicate();

    myCase.getOperator("Coupling").apply();

    // Computation and output of the results
    if ( std::isnan( sLatticeNS.getStatistics().getAverageEnergy() ) ) {
      clout << "Simulation breaks at iT=" << iT << std::endl;
      if (iB > 10) break;
      iB++;
    }

  }
  timer.stop();
  timer.printSummary();
}

int main( int argc, char *argv[] )
{
  initialize( &argc, &argv );

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<DOMAIN_EXTENT          >({150e-3,35e-3,35e-3}); // domain sizes length and diameter [physical units]
    myCaseParameters.set<RESOLUTION             >(145);               // resolution inlet [lattice units]
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.5004); // tau liquid [lattice units]
    myCaseParameters.set<REYNOLDS               >(2000.);   // Reynolds gas phase [physical units]
    myCaseParameters.set<RHO_LIQUID             >(1000.);   // liquid density [physical units]
    myCaseParameters.set<RHO_VAPOR              >(1.2);     // gas density [physical units]
    myCaseParameters.set<NU_LIQUID              >(9e-7);    // liquid kinematic viscosity [physical units]
    myCaseParameters.set<NU_VAPOR               >(13.*9e-7);// gas kinematic viscosity [physical units]
    myCaseParameters.set<SURFACE_TENSION        >(0.072);   // surface tension [physical units]
    myCaseParameters.set<C_RHO                  >(100.);    // conversion factor density [physical units]
    myCaseParameters.set<parameters::THETA>(M_PI*45./180.); // contact angle [radians]
    myCaseParameters.set<PHYS_MAX_T      >(2.0);     // max simulation time [physical units]
    myCaseParameters.set<MAX_PHYS_START_T>(1.0);    // max ramp time [physical units]
    myCaseParameters.set<PHYS_STAT_ITER_T>(0.05);   // simulation statistics time [physical units]
    myCaseParameters.set<parameters::INTERFACE_WIDTH>(5.); // resolution inlet [lattice units]
    myCaseParameters.set<GRAVITY>({0.,-9.81,0.});      // physical acceleration gravity [physical units]
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