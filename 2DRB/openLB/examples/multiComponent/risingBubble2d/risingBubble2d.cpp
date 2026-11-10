/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2026 Tim Bingert
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

/* risingBubble2d.cpp
 * This example showcases a rising bubble under different conditions
 * of bubble Reynolds and Bond number with its current configuration
 * resulting in a cap-like shape.
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;

using MyCase = Case<
  NavierStokes, Lattice<float, D2Q9<RHO,NABLARHO,LAPLACERHO,FORCE,TAU_EFF,SCALAR>>,
  Component1,  Lattice<float, D2Q9<FORCE,VELOCITY,OLD_PHIU,STATISTIC>>
>;

using NSBulkDynamics = CorrectedMultiPhaseIncompressibleRLBdynamics<MyCase::value_t,MyCase::descriptor_t_of<NavierStokes>>;
using PFBulkDynamics = ConservativePhaseFieldBGKdynamics<MyCase::value_t,MyCase::descriptor_t_of<Component1>>;
using MixtureRules = LinearTauViscosity;
using Coupling = CorrectedLiangPostProcessor<MixtureRules>;

namespace olb::parameters {

  // Explanations for global variables are in main()
  struct C_RHO                      : public descriptors::FIELD_BASE<1> { };
  struct PF_LATTICE_RELAXATION_TIME : public descriptors::FIELD_BASE<1> { };
  struct TAU_VAPOR                  : public descriptors::FIELD_BASE<1> { };
  struct BOND                       : public descriptors::FIELD_BASE<1> { };
  struct RISING_VELOCITY            : public descriptors::FIELD_BASE<1> { };
  struct G_CONST                    : public descriptors::FIELD_BASE<1> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extent, origin);

  const T dx = params.get<parameters::PHYS_CHAR_LENGTH>() / params.get<parameters::RESOLUTION>();
  Mesh<T,MyCase::d> mesh(cuboid, dx, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({ true, false });

  return mesh;
}

void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  auto& geometry = myCase.getGeometry();

  geometry.rename( 0,2 );
  geometry.rename( 2, 1, {0, 1} );

  geometry.clean();
  geometry.innerClean();
  geometry.checkForErrors();
  geometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& latticeNS = myCase.getLattice(NavierStokes{});
  auto& latticePF = myCase.getLattice(Component1{});

  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

  const int N           = params.get<parameters::RESOLUTION>();
  const T char_l        = params.get<parameters::PHYS_CHAR_LENGTH>();
  const T rho_g         = params.get<parameters::RHO_VAPOR>();
  const T rho_l         = params.get<parameters::RHO_LIQUID>();
  const T nu_l          = params.get<parameters::NU_LIQUID>();
  const T tau_l         = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T tau_g         = params.get<parameters::TAU_VAPOR>();
  const T tau_mobil     = params.get<parameters::PF_LATTICE_RELAXATION_TIME>();
  const T sigma_lattice = params.get<parameters::SURFACE_TENSION>();
  const T w             = params.get<parameters::INTERFACE_WIDTH>();
  const T C_rho         = params.get<parameters::C_RHO>();
  const T g             = params.get<parameters::G_CONST>();
  const T theta         = params.get<parameters::THETA>();

  latticeNS.setUnitConverter<MultiPhaseUnitConverterFromRelaxationTime<T,NSDESCRIPTOR>>(
    int   {N},                  // resolution
    (T)   tau_l,                // lattice relaxation time
    (T)   rho_l/C_rho,          // lattice density heavier phase
    (T)   char_l,               // charPhysLength: reference length of simulation geometry in __m__
    (T)   nu_l,                 // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   rho_l                 // physDensity: physical density in __kg / m^3__
  );

  const auto& converter = latticeNS.getUnitConverter();
  converter.print();

  latticePF.setUnitConverter(converter);

  // conversion properties
  const T dx = converter.getPhysDeltaX();
  const T C_density = converter.getConversionFactorDensity();
  clout << "Lattice Surface Tension: " << sigma_lattice << " should be smaller than 0.1!" << std::endl;

  // define lattice Dynamics
  dynamics::set<NSBulkDynamics>(latticeNS, geometry, 1);
  dynamics::set<PFBulkDynamics>(latticePF, geometry, 1);

  // walls
  const T Nx = converter.getLatticeLength(params.get<parameters::DOMAIN_EXTENT>()[0]);
  const T Ny = converter.getLatticeLength(params.get<parameters::DOMAIN_EXTENT>()[1]);
  std::vector<T> origin = {-T(2)*dx, dx*T(0.5)};
  std::vector<T> extent = {dx*(Nx+T(4)), dx*(Ny-T(2))};
  IndicatorCuboid2D<T> wallLocation( extent, origin );
  setBouzidiBoundary(latticeNS, geometry, 2, wallLocation);
  setBouzidiPhaseField<IsoPhaseFieldCurvedWall2D>(latticePF, geometry, 2, wallLocation);

  auto zero = geometry.getMaterialIndicator(0);
  auto bulk = geometry.getMaterialIndicator(1);
  auto wall = geometry.getMaterialIndicator(2);
  auto all = geometry.getMaterialIndicator({0,1,2});
  latticePF.addPostProcessor<stage::PreCoupling>(bulk, meta::id<PhiStatistics>());

  auto& coupling = myCase.setCouplingOperator(
    "Coupling",
    Coupling{},
    names::NavierStokes{}, latticeNS,
    names::Component1{}, latticePF
  );
  coupling.restrictTo(geometry.getMaterialIndicator(1));

  phasefield::init<Coupling,MixtureRules>(
    latticeNS, latticePF,    // both lattices
    coupling, all,           // phase field coupling, all material numbers
    wall, zero,              // wall and zero material numbers
    sigma_lattice, w, theta, // phase field lattice parameters
    tau_g, tau_l,            // relaxation times for mixture rule
    rho_g/C_rho, rho_l/C_rho // lattice densities for mixture rule
  );

  coupling.setParameter<Coupling::GRAVITY>({T(0),-g});
  coupling.setParameter<Coupling::BODY_FORCE>({T(0),T(0)});

  latticeNS.setParameter<descriptors::OMEGA>( 1./tau_l );
  latticePF.setParameter<descriptors::OMEGA>( 1./tau_mobil );

  {
    auto& communicator = latticePF.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(1);
    communicator.requestField<STATISTIC>();
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  OstreamManager clout( std::cout,"setInitialValues" );
  clout << "Set Initial Values ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();
  auto& latticeNS = myCase.getLattice(NavierStokes{});
  auto& latticePF = myCase.getLattice(Component1{});
  const auto& converter = latticeNS.getUnitConverter();

  const T diameter_lattice = params.get<parameters::RESOLUTION>();
  const T rho_g            = params.get<parameters::RHO_VAPOR>();
  const T rho_l            = params.get<parameters::RHO_LIQUID>();
  const T tau_l            = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T tau_g            = params.get<parameters::TAU_VAPOR>();
  const T sigma_lattice    = params.get<parameters::SURFACE_TENSION>();
  const T w                = params.get<parameters::INTERFACE_WIDTH>();
  const T C_rho            = params.get<parameters::C_RHO>();
  const T Nx               = converter.getLatticeLength(params.get<parameters::DOMAIN_EXTENT>()[0]);
  const T dx               = converter.getPhysDeltaX();
  const T C_sigma          = converter.getConversionFactorSurfaceTension();
  const T C_p              = C_sigma/dx; // Conversion factor pressure
  const T rho_g_lattice    = rho_g / C_rho;
  const T rho_l_lattice    = rho_l / C_rho;
  const std::vector<T> u0 { 0,0 };

  auto bulk = geometry.getMaterialIndicator(1);
  auto wall = geometry.getMaterialIndicator(2);
  auto all = geometry.getMaterialIndicator({0,1,2});

  std::vector<T> pos = {dx*Nx/T(2), dx*diameter_lattice};
  std::shared_ptr<AnalyticalF2D<T,T>> phi0 = std::make_shared<CircularInterface2D<T>>(pos, dx*diameter_lattice/T(2), dx*w, T(1), true);
  std::shared_ptr<AnalyticalF2D<T,T>> rho0( rho_g_lattice + (rho_l_lattice - rho_g_lattice) * phi0 );
  std::shared_ptr<AnalyticalF2D<T,T>> tau0( tau_g + (tau_l-tau_g) * phi0 );
  std::shared_ptr<AnalyticalF2D<T,T>> bubblePressure = std::make_shared<LaplacePressure2D<T>>( pos, dx*diameter_lattice/T(2), dx*w, sigma_lattice*C_sigma/C_p );
  const T halfYLPressure = T(2)*sigma_lattice/diameter_lattice/T(2);

  std::shared_ptr<AnalyticalF2D<T,T>> p0 = bubblePressure + halfYLPressure;
  std::shared_ptr<AnalyticalF2D<T,T>> p0_phys = p0*C_p;

  fields::set<descriptors::RHO>( latticeNS, all, *rho0 );
  fields::set<descriptors::TAU_EFF>( latticeNS, all, *tau0 );
  fields::set<descriptors::SCALAR>( latticeNS, all, T(1) );
  fields::set<descriptors::SCALAR>( latticeNS, wall, T(2) );
  fields::set<descriptors::OLD_PHIU>( latticePF, all, u0 );

  momenta::setOrderParameter(latticePF, all, *phi0);
  momenta::setIncompressiblePressure(latticeNS, all, *p0_phys);

  AnalyticalConst2D<T,T> _u0(u0);
  latticeNS.iniEquilibrium( all, *p0, _u0 );

  latticeNS.initialize();
  latticePF.initialize();

  latticePF.iniEquilibrium( all, *phi0, _u0 );
  latticeNS.iniEquilibrium( all, *p0, _u0 );

  latticeNS.setProcessingContext(ProcessingContext::Simulation);
  latticePF.setProcessingContext(ProcessingContext::Simulation);

  latticePF.getCommunicator(stage::PreCoupling()).communicate();

  clout << "Set Initial Values ... OK" << std::endl;
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
  using PFDESCRIPTOR = MyCase::descriptor_t_of<Component1>;
  auto& params = myCase.getParameters();
  auto& latticeNS = myCase.getLattice(NavierStokes{});
  auto& latticePF = myCase.getLattice(Component1{});
  const auto& converter = latticeNS.getUnitConverter();

  const std::size_t vtkIter  = converter.getLatticeTime(params.get<parameters::PHYS_VTK_ITER_T>());
  const std::size_t statIter = converter.getLatticeTime(params.get<parameters::PHYS_STAT_ITER_T>());

  SuperVTMwriter2D<T> vtmWriter( "risingBubble2d" );
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
    latticeNS.getStatistics().print( iT, converter.getPhysTime(iT) );
    latticePF.getStatistics().print( iT, converter.getPhysTime(iT) );
  }

  // Writes the VTK files
  if ( iT%vtkIter==0 ) {
    latticeNS.setProcessingContext(ProcessingContext::Evaluation);
    latticePF.setProcessingContext(ProcessingContext::Evaluation);
    SuperLatticePhysIncPressure2D<T, NSDESCRIPTOR> p_total( latticeNS, converter );
    p_total.getName() = "p_total";
    SuperLatticePhysField2D<T, NSDESCRIPTOR, RHO> rho( latticeNS, converter.getConversionFactorDensity() );
    rho.getName() = "rho";
    SuperLatticePhysVelocity2D<T, NSDESCRIPTOR> velocity( latticeNS, converter );
    velocity.getName() = "u";
    SuperLatticeField2D<T, PFDESCRIPTOR, STATISTIC> phi( latticePF );
    phi.getName() = "phi";

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
  auto& latticePF = myCase.getLattice(Component1{});
  auto& params = myCase.getParameters();

  const T physMaxT = params.get<parameters::PHYS_MAX_T>();
  const std::size_t iTmax = latticeNS.getUnitConverter().getLatticeTime(physMaxT);
  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());

  clout << "Starting simulation ..." << std::endl;
  timer.start();

  for (std::size_t iT=0; iT<=iTmax; ++iT ) {
    // Collide and stream (and coupling) execution
    latticeNS.collideAndStream();
    latticePF.collideAndStream();

    latticePF.getCommunicator(stage::PreCoupling()).communicate();
    latticePF.executePostProcessors(stage::PreCoupling());
    latticePF.getCommunicator(stage::PreCoupling()).communicate();

    myCase.getOperator("Coupling").apply();

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
    using T = MyCase::value_t;
    myCaseParameters.set<PHYS_CHAR_LENGTH>(40e-6);          // bubble radius [phys units]
    myCaseParameters.set<RESOLUTION      >(40);             // bubble radius [lattice units]
    myCaseParameters.set<DOMAIN_EXTENT   >([](MyCase::ParametersD& p){
      return Vector<T,2> {4 * p.get<PHYS_CHAR_LENGTH>(), 8 * p.get<PHYS_CHAR_LENGTH>()};
    });
    myCaseParameters.set<C_RHO                     >(100.);           // conversion factor density [physical units]
    myCaseParameters.set<PF_LATTICE_RELAXATION_TIME>(0.8 ); // lattice relaxation time for interface mobility [lattice units]
    myCaseParameters.set<NU_LIQUID                 >(1e-6);           // physViscosity liquid [physical units]
    myCaseParameters.set<RHO_VAPOR                 >(1);              // physDensity gas/vapor [physical units]
    myCaseParameters.set<RHO_LIQUID                >(1000);          // physDensity liquid [physical units]
    myCaseParameters.set<REYNOLDS                  >(35);             // Reynolds number of rising bubble
    myCaseParameters.set<BOND                      >(100);            // Bond number of rising bubble
    myCaseParameters.set<RISING_VELOCITY           >(0.02);           // Bubble rising velocity [lattice units]
    myCaseParameters.set<G_CONST                   >([](MyCase::ParametersD& p) {
      return p.get<RISING_VELOCITY>() * p.get<RISING_VELOCITY>() / p.get<RESOLUTION>();
    });                                                     // gravity [lattice units]
    myCaseParameters.set<LATTICE_RELAXATION_TIME>([] (MyCase::ParametersD& p) {
      return 3.*(p.get<RISING_VELOCITY>() * p.get<RESOLUTION>() / p.get<REYNOLDS>()) + 0.5;
    });                                                     // lattice relaxation time H2O liquid [lattice units]
    myCaseParameters.set<TAU_VAPOR       >([](MyCase::ParametersD& p) {
      return 10.*(p.get<LATTICE_RELAXATION_TIME>() - 0.5) + 0.5;
    });                                                     // lattice relaxation time air gas [lattice units]
    myCaseParameters.set<parameters::SURFACE_TENSION >([](MyCase::ParametersD& p){
      const T rho_l = p.get<RHO_LIQUID>();
      const T C_rho = p.get<C_RHO>();
      const T g = p.get<G_CONST>();
      const T lat_diameter = p.get<RESOLUTION>();
      const T Bo = p.get<BOND>();
      return rho_l/C_rho*g*lat_diameter*lat_diameter/Bo;
    });                                                     // surface tension from Bond number [lattice units]
    myCaseParameters.set<parameters::THETA          >(M_PI*90./180.);  // contact angle (<90 is wetting) [radians]
    myCaseParameters.set<parameters::INTERFACE_WIDTH>(5.);  // interface width [lattice units]
    myCaseParameters.set<PHYS_MAX_T         >(0.0012);      // max physical time [physical units]
    myCaseParameters.set<PHYS_VTK_ITER_T    >(0.0012/30.);  // vtk output interval [physical units]
    myCaseParameters.set<PHYS_STAT_ITER_T   >(0.0012/30.);  // statistics output interval [physical units]
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
