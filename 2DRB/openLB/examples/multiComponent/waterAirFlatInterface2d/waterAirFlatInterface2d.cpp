/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2023 Tim Bingert
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

/* waterAirFlatInterface2d.cpp
 * Benchmark for thermodynamic consistency of the multi-component-
 * multi-phase Shan-Chen model at standard conditions.
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<
  Component1, Lattice<double, D2Q9<VELOCITY,FORCE,EXTERNAL_FORCE,STATISTIC,SCALAR,PSI>>,
  Component2, Lattice<double, D2Q9<VELOCITY,FORCE,EXTERNAL_FORCE,STATISTIC,SCALAR,PSI>>,
  Component3, Lattice<double, D2Q9<VELOCITY,FORCE,EXTERNAL_FORCE,STATISTIC,SCALAR,PSI>>
>;

/*==============================================================
 * go to src/dynamics/shanChenForcedPostProcessor.h and enable *
 * the THIRD_COMPONENT by uncommenting #define THIRD_COMPONENT *
 ==============================================================*/

using BulkDynamics = MultiComponentForcedBGKdynamics<MyCase::value_t, MyCase::descriptor_t>;
constexpr unsigned N_COMPONENTS = 3;
using COUPLING = MCMPForcedPostProcessor<N_COMPONENTS>;
using STATISTICS = RhoPsiStatistics<interaction::MCPRpseudoPotential<N_COMPONENTS>,N_COMPONENTS>;

namespace olb::parameters {

struct A0_LATTICE : public descriptors::FIELD_BASE<1> { };
struct FEED       : public descriptors::FIELD_BASE<N_COMPONENTS> { };
struct A          : public descriptors::FIELD_BASE<N_COMPONENTS> { };
struct B          : public descriptors::FIELD_BASE<N_COMPONENTS> { };
struct M          : public descriptors::FIELD_BASE<N_COMPONENTS> { };
struct T_C        : public descriptors::FIELD_BASE<N_COMPONENTS> { };
struct P_C        : public descriptors::FIELD_BASE<N_COMPONENTS> { };
struct OMEGA_PER_COMPONENT : public descriptors::FIELD_BASE<N_COMPONENTS> { };
struct DEVI       : public descriptors::FIELD_BASE<N_COMPONENTS> { };

struct RHO0_L      : public descriptors::FIELD_BASE<N_COMPONENTS> { };
struct RHO0_V      : public descriptors::FIELD_BASE<N_COMPONENTS> { };

struct ALPHA      : public descriptors::FIELD_BASE<N_COMPONENTS*N_COMPONENTS> { };
struct G_I        : public descriptors::FIELD_BASE<N_COMPONENTS*N_COMPONENTS> { };
struct G_II       : public descriptors::FIELD_BASE<N_COMPONENTS*N_COMPONENTS> { };

struct VLE        : public descriptors::FIELD_BASE<2*(N_COMPONENTS+1)> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extent, origin);

  const T dx = 1.; // lattice units case
  Mesh<T,MyCase::d> mesh(cuboid, dx, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({true,true});
  return mesh;
}

void prepareGeometry( MyCase& myCase )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;
  auto& geometry = myCase.getGeometry();
  geometry.rename( 0,1 );
  // bulk, MN=1

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

  auto& sLattice1 = myCase.getLattice(Component1{});
  auto& sLattice2 = myCase.getLattice(Component2{});
  auto& sLattice3 = myCase.getLattice(Component3{});

  using DESCRIPTOR = MyCase::descriptor_t;

  const int N = params.get<parameters::RESOLUTION>();
  const T tau_H2O = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T L_char = params.get<parameters::PHYS_CHAR_LENGTH>();
  const T viscosityH2O = params.get<parameters::NU_LIQUID>();
  const T viscosityAir = params.get<parameters::NU_VAPOR>();
  const T tau_Air = viscosityAir/viscosityH2O*(tau_H2O-0.5)+0.5;
  const T a_0L = params.get<parameters::A0_LATTICE>();
  const T eps = params.get<parameters::EPSILON>();
  const T sigma = params.get<parameters::SURFACE_TENSION>();
  const T pressure = params.get<parameters::PHYS_CHAR_PRESSURE>();
  const T temperature = params.get<parameters::PHYS_CHAR_TEMPERATURE>();

  const Vector<T,N_COMPONENTS> z = params.get<parameters::FEED>();
  const Vector<T,N_COMPONENTS> a = params.get<parameters::A>();
  const Vector<T,N_COMPONENTS> b = params.get<parameters::B>();
  const Vector<T,N_COMPONENTS> M = params.get<parameters::M>();
  const Vector<T,N_COMPONENTS> T_c = params.get<parameters::T_C>();
  const Vector<T,N_COMPONENTS> p_c = params.get<parameters::P_C>();
  const Vector<T,N_COMPONENTS> omega = params.get<parameters::OMEGA_PER_COMPONENT>();
  const Vector<T,N_COMPONENTS> devi = params.get<parameters::DEVI>();

  Vector<T,N_COMPONENTS> rho0L(0.,0.,0.);
  Vector<T,N_COMPONENTS> rho0V(0.,0.,0.);

  const Vector<T,N_COMPONENTS*N_COMPONENTS> alpha = params.get<parameters::ALPHA>();
  const Vector<T,N_COMPONENTS*N_COMPONENTS> g_I = params.get<parameters::G_I>();
  const Vector<T,N_COMPONENTS*N_COMPONENTS> g_II = params.get<parameters::G_II>();

  sLattice1.setUnitConverter<MultiPhaseUnitConverter<T,DESCRIPTOR>>(
    int   {N},         // resolution
    (T)   L_char,      // charPhysLength: reference length of simulation geometry
    (T)   0.,          // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   viscosityH2O,// physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   a[0],        // physEoSa: H2O energy parameter in __kg m^5 / mol^2 s^2__
    (T)   a_0L,        // latticeEoSa: first component's energy parameter in lattice units
    (T)   b[0],        // physEoSb: H2O co-volume parameter in __m^3 / mol__
    (T)   M[0],        // physMolarMass: H2O molar mass for EoS in __kg / mol__
    (T)   sigma,       // physSurfaceTension: physical surface tension of mixture in __kg / s^2__
    (T)   temperature, // charPhysTemperature: temperature of VLE in __K__
    (T)   pressure     // charPhysPressure: pressure of VLE in __kg / m s^2__
  );

  const auto& converter = sLattice1.getUnitConverter();
  converter.print();

  sLattice2.setUnitConverter(converter);
  sLattice3.setUnitConverter(converter);

  // define lattice Dynamics
  dynamics::set<BulkDynamics>(sLattice1, geometry.getMaterialIndicator({1}));
  dynamics::set<BulkDynamics>(sLattice2, geometry.getMaterialIndicator({1}));
  dynamics::set<BulkDynamics>(sLattice3, geometry.getMaterialIndicator({1}));

  //thermodynamic initial conditions in lattice units
  T p_L = pressure/converter.getConversionFactorPressure();
  T T_L = temperature/converter.getConversionFactorTemperature();
  std::vector<T> z_L(N_COMPONENTS), a_L(10), b_L(10), M_L(10), Tc_L(10), pc_L(10), omega_L(10), devi_L(10), alpha_L(100), gI_L(100), gII_L(100);
  for(std::size_t i = 0; i < N_COMPONENTS; i++){
    z_L[i]     = z[i];
    a_L[i]     = a[i]/converter.getConversionFactorEoSa();
    b_L[i]     = b[i]/converter.getConversionFactorEoSb();
    M_L[i]     = M[i]/converter.getConversionFactorMolarMass();
    Tc_L[i]    = T_c[i]/converter.getConversionFactorTemperature();
    pc_L[i]    = p_c[i]/converter.getConversionFactorPressure();
    omega_L[i] = omega[i];
    devi_L[i]  = devi[i];
  }
  T C_R = 8.314462618;
  T C_temp = converter.getConversionFactorTemperature();
  for(std::size_t i = 0; i < N_COMPONENTS*N_COMPONENTS; i++){
    alpha_L[i] = alpha[i];
    gI_L[i]    = g_I[i]/(C_R*C_temp);
    gII_L[i]   = g_II[i]/C_R;
  }

  MultiComponentPengRobinson VLE(p_L, T_L, z_L, a_L, b_L, M_L, Tc_L, pc_L, omega_L, devi_L, alpha_L, gI_L, gII_L);
  T beta0 = (z[1]+z[2]);
  std::vector<T> vxVLE = VLE.iterate_VLE(1e-11, beta0);   // Molar volumes and fractions of equilibrium phases [lattice units]
  clout <<"VLE: "<<vxVLE[0]<<", "<<vxVLE[1]<<", "<<vxVLE[2]<<", "<<vxVLE[3]<<", "<<vxVLE[4]<<", "<<vxVLE[5]<<", "<<vxVLE[6]<<", "<<vxVLE[7]<<std::endl;
  std::vector<T> chi = VLE.getChis(3);                    // Force split factor from VLE
  clout <<"Chis: "<<chi[0]<<", "<<chi[1]<<", "<<chi[2]<<std::endl;
  params.set<parameters::VLE>(vxVLE);

  for(std::size_t i=0; i<N_COMPONENTS; i++){
    rho0L[i] = 1./vxVLE[0]*vxVLE[i+2]*M_L[i]*converter.getConversionFactorDensity();
    rho0V[i] = 1./vxVLE[1]*vxVLE[i+5]*M_L[i]*converter.getConversionFactorDensity();
  }

  params.set<parameters::RHO0_L>(rho0L);
  params.set<parameters::RHO0_V>(rho0V);

  sLattice1.setParameter<descriptors::OMEGA>( 1./tau_H2O );
  sLattice2.setParameter<descriptors::OMEGA>( 1./tau_Air );
  sLattice3.setParameter<descriptors::OMEGA>( 1./tau_Air );

  auto& coupling = myCase.setCouplingOperator(
    "Coupling",
    COUPLING{},
    names::Component1{}, sLattice1,
    names::Component2{}, sLattice2,
    names::Component3{}, sLattice3);

  auto& statistics = myCase.setCouplingOperator(
    "Statistics",
    STATISTICS{},
    names::Component1{}, sLattice1,
    names::Component2{}, sLattice2,
    names::Component3{}, sLattice3);

  coupling.setParameter<COUPLING::CHI>(chi);
  coupling.setParameter<COUPLING::G>(-1.);
  coupling.setParameter<COUPLING::SIGMA>(sigma);
  coupling.setParameter<COUPLING::EPSILON>(eps);

  statistics.setParameter<STATISTICS::TEMPERATURE>(T_L);
  statistics.setParameter<STATISTICS::G>(-1.);
  statistics.setParameter<STATISTICS::K>(1.);
  statistics.setParameter<STATISTICS::A>(a_L);
  statistics.setParameter<STATISTICS::B>(b_L);
  statistics.setParameter<STATISTICS::MM>(M_L);
  statistics.setParameter<STATISTICS::TCRIT>(Tc_L);
  statistics.setParameter<STATISTICS::DEVI>(devi_L);
  statistics.setParameter<STATISTICS::ALPHA>(alpha_L);
  statistics.setParameter<STATISTICS::GI>(gI_L);
  statistics.setParameter<STATISTICS::GII>(gII_L);

  {
    auto& communicator = sLattice1.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(1);
    communicator.requestField<STATISTIC>();
    communicator.requestField<PSI>();
    communicator.exchangeRequests();
  }

  {
    auto& communicator = sLattice2.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(1);
    communicator.requestField<STATISTIC>();
    communicator.exchangeRequests();
  }

  {
    auto& communicator = sLattice3.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(1);
    communicator.requestField<STATISTIC>();
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& sLattice1 = myCase.getLattice(Component1{});
  auto& sLattice2 = myCase.getLattice(Component2{});
  auto& sLattice3 = myCase.getLattice(Component3{});

  const auto& converter = sLattice1.getUnitConverter();

  const Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const int N = params.get<parameters::RESOLUTION>();
  Vector g = params.get<parameters::GRAVITY>();

  const Vector<T,N_COMPONENTS> M = params.get<parameters::M>();
  const Vector<T,2*(N_COMPONENTS+1)> vle = params.get<parameters::VLE>();

  std::vector<T> M_L(10);
  for(std::size_t i = 0; i < N_COMPONENTS; i++){
    M_L[i] = M[i]/converter.getConversionFactorMolarMass();
  }

  std::vector<T> v = {0., 0.};
  AnalyticalConst2D<T,T> zeroVelocity( v );

  AnalyticalConst2D<T,T> liquidH2O ( 1./vle[0]*vle[2]*M_L[0] );
  AnalyticalConst2D<T,T> liquidN2  ( 1./vle[0]*vle[3]*M_L[1] );
  AnalyticalConst2D<T,T> liquidO2  ( 1./vle[0]*vle[4]*M_L[2] );

  SmoothIndicatorFactoredCuboid2D<T,T> vaporH2O( {extent[0]/2., extent[1]/2.}, 0, N, 4., 0, {0,0}, 0, (1./vle[1]*vle[5]*M_L[0] - 1./vle[0]*vle[2]*M_L[0]) );
  SmoothIndicatorFactoredCuboid2D<T,T> vaporN2 ( {extent[0]/2., extent[1]/2.}, 0, N, 4., 0, {0,0}, 0, (1./vle[1]*vle[6]*M_L[1] - 1./vle[0]*vle[3]*M_L[1]) );
  SmoothIndicatorFactoredCuboid2D<T,T> vaporO2 ( {extent[0]/2., extent[1]/2.}, 0, N, 4., 0, {0,0}, 0, (1./vle[1]*vle[7]*M_L[2] - 1./vle[0]*vle[4]*M_L[2]) );

  AnalyticalIdentity2D<T,T> rhoH2O( liquidH2O + vaporH2O );
  AnalyticalIdentity2D<T,T> rhoN2( liquidN2 + vaporN2 );
  AnalyticalIdentity2D<T,T> rhoO2( liquidO2 + vaporO2 );

  sLattice1.defineRhoU( geometry, 1, rhoH2O, zeroVelocity );
  sLattice2.defineRhoU( geometry, 1, rhoN2, zeroVelocity );
  sLattice3.defineRhoU( geometry, 1, rhoO2, zeroVelocity );

  sLattice1.iniEquilibrium( geometry, 1, rhoH2O, zeroVelocity );
  sLattice2.iniEquilibrium( geometry, 1, rhoN2, zeroVelocity );
  sLattice3.iniEquilibrium( geometry, 1, rhoO2, zeroVelocity );

  g = g/(converter.getPhysDeltaX()/converter.getPhysDeltaT()/converter.getPhysDeltaT());
  AnalyticalConst2D<T,T> f( g );
  fields::set<descriptors::EXTERNAL_FORCE>(sLattice1, geometry.getMaterialIndicator({1}), f);
  fields::set<descriptors::EXTERNAL_FORCE>(sLattice2, geometry.getMaterialIndicator({1}), f);
  fields::set<descriptors::EXTERNAL_FORCE>(sLattice3, geometry.getMaterialIndicator({1}), f);

  sLattice1.initialize();
  sLattice2.initialize();
  sLattice3.initialize();
  myCase.getOperator("Statistics").apply();
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{ }

template <typename T>
T getResults( MyCase& myCase,
              util::Timer<MyCase::value_t>& timer,
              std::size_t iT )
{
  OstreamManager clout( std::cout,"getResults" );
  auto& params = myCase.getParameters();

  auto& sLattice1 = myCase.getLattice(Component1{});
  auto& sLattice2 = myCase.getLattice(Component2{});
  auto& sLattice3 = myCase.getLattice(Component3{});

  using DESCRIPTOR = MyCase::descriptor_t;
  const auto& converter = sLattice1.getUnitConverter();

  const Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const Vector<T,N_COMPONENTS> rho0L = params.get<parameters::RHO0_L>();
  const Vector<T,N_COMPONENTS> rho0V = params.get<parameters::RHO0_V>();
  const int statIter = params.get<parameters::LATTICE_STAT_ITER_T>();
  const int saveIter = params.get<parameters::LATTICE_VTK_ITER_T>();

  SuperVTMwriter2D<T> vtkWriter( "waterAirFlatInterface2d" );
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

  T error = 0;
  if ( iT%saveIter==0 ) {
    sLattice1.setProcessingContext(ProcessingContext::Evaluation);
    sLattice2.setProcessingContext(ProcessingContext::Evaluation);
    sLattice3.setProcessingContext(ProcessingContext::Evaluation);

    //Factors for conversion back to physical units
    AnalyticalConst2D<T,T> _C_rho( converter.getConversionFactorDensity() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> C_rho(_C_rho, sLattice1);
    AnalyticalConst2D<T,T> _C_u( converter.getConversionFactorVelocity() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> C_u(_C_u, sLattice2);
    AnalyticalConst2D<T,T> _C_p( converter.getConversionFactorPressure() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> C_p(_C_p, sLattice1);

    SuperLatticeDensity2D<T, DESCRIPTOR> density1L( sLattice1 );
    SuperIdentity2D<T,T> density1( C_rho*density1L );
    density1.getName() = "rhoH2O";
    SuperLatticeDensity2D<T, DESCRIPTOR> density2L( sLattice2 );
    SuperIdentity2D<T,T> density2( C_rho*density2L );
    density2.getName() = "rhoN2";
    SuperLatticeDensity2D<T, DESCRIPTOR> density3L( sLattice3 );
    SuperIdentity2D<T,T> density3( C_rho*density3L );
    density3.getName() = "rhoO2";
    SuperIdentity2D<T,T> density( density1+density2+density3 );
    density.getName() = "rho";

    SuperLatticeVelocity2D<T, DESCRIPTOR> velocityL( sLattice1 );
    SuperIdentity2D<T,T> velocity( C_u*velocityL );
    velocity.getName() = "velocity";

    SuperLatticeField2D<T, DESCRIPTOR, SCALAR> bulkPressureL( sLattice1 );
    SuperIdentity2D<T,T> bulkPressure( C_p*bulkPressureL );
    bulkPressure.getName() = "bulkPressure";

    vtkWriter.addFunctor( density1 );
    vtkWriter.addFunctor( density2 );
    vtkWriter.addFunctor( density3 );
    vtkWriter.addFunctor( density );
    vtkWriter.addFunctor( velocity );
    vtkWriter.addFunctor( bulkPressure );
    vtkWriter.write( iT );

    // calculate total error of phase compositions
    AnalyticalFfromSuperF2D<T, T> interpolRho1( density1, true, 1);
    AnalyticalFfromSuperF2D<T, T> interpolRho2( density2, true, 1);
    AnalyticalFfromSuperF2D<T, T> interpolRho3( density3, true, 1);
    double position[2] = { 0.5*extent[0], 0.5*extent[1]+1 };
    double liquidDensities[3] = {0.,0.,0.};
    double vaporDensities[3] = {0.,0.,0.};
    interpolRho1(&vaporDensities[0], position);
    interpolRho2(&vaporDensities[1], position);
    interpolRho3(&vaporDensities[2], position);
    position[1] = 1.;
    interpolRho1(&liquidDensities[0], position);
    interpolRho2(&liquidDensities[1], position);
    interpolRho3(&liquidDensities[2], position);
    double liquidError[3] = {abs(liquidDensities[0]-rho0L[0])/rho0L[0],
                             abs(liquidDensities[1]-rho0L[1])/rho0L[1],
                             abs(liquidDensities[2]-rho0L[2])/rho0L[2]};
    double vaporError[3] = {abs(vaporDensities[0]-rho0V[0])/rho0V[0],
                            abs(vaporDensities[1]-rho0V[1])/rho0V[1],
                            abs(vaporDensities[2]-rho0V[2])/rho0V[2]};
    for(std::size_t i = 0; i < N_COMPONENTS; i++){
      error += vaporError[i]*vaporError[i] + liquidError[i]*liquidError[i];
    }
    error = util::pow(error/6.,0.5);
    clout << "Average error of all component densities: " << std::setprecision(6) << error << std::endl;
  }
  return error;
}

void simulate( MyCase& myCase )
{
  using T = MyCase::value_t;
  auto& params = myCase.getParameters();

  auto& sLattice1 = myCase.getLattice(Component1{});
  auto& sLattice2 = myCase.getLattice(Component2{});
  auto& sLattice3 = myCase.getLattice(Component3{});

  sLattice1.writeSummary("lattice1");
  sLattice2.writeSummary("lattice2");
  sLattice3.writeSummary("lattice3");

  const std::size_t iTmax = params.get<parameters::LATTICE_MAX_T>();

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT <= iTmax; ++iT) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    sLattice1.collideAndStream();
    sLattice2.collideAndStream();
    sLattice3.collideAndStream();

    myCase.getOperator("Statistics").apply();
    sLattice1.getCommunicator(stage::PreCoupling()).communicate();
    sLattice2.getCommunicator(stage::PreCoupling()).communicate();
    sLattice3.getCommunicator(stage::PreCoupling()).communicate();
    myCase.getOperator("Coupling").apply();

    /// === Step 8.3: Computation and Output of the Results ===
    getResults<T>(myCase, timer, iT);
    if ( std::isnan( sLattice1.getStatistics().getAverageEnergy() ) ) {
      break;
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
    myCaseParameters.set<DOMAIN_EXTENT          >({4, 100}); // domain size [lattice units]
    myCaseParameters.set<RESOLUTION             >(50);       // resolution [lattice units]
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.6);      // tau_H2O [lattice units]
    myCaseParameters.set<LATTICE_MAX_T          >(250000);   // max iterations [lattice units]
    myCaseParameters.set<LATTICE_VTK_ITER_T     >(1000);     // vtk iterations [lattice units]
    myCaseParameters.set<LATTICE_STAT_ITER_T    >(1000);     // statistics iterations [lattice units]
    myCaseParameters.set<PHYS_CHAR_LENGTH       >(1e-6);        // charPhysLength [physical units]
    myCaseParameters.set<NU_VAPOR               >(1.532e-5);    // physViscosity N2+O2 lattice [physical units]
    myCaseParameters.set<NU_LIQUID              >(1e-6);        // physViscosity H2O lattice [physical units]
    myCaseParameters.set<GRAVITY                >({0.,0.});     // gravitational acceleration [physical units]
    myCaseParameters.set<A0_LATTICE             >(3./245.);     // tune for stability/accuracy [lattice units]
    myCaseParameters.set<parameters::EPSILON>(2.2355);   // tune for chemical potential equilibrium in phases []
    myCaseParameters.set<SURFACE_TENSION      >(0.07);        // For unit conversion [physical units], not crucial for phase composition
    myCaseParameters.set<PHYS_CHAR_PRESSURE   >(1.013e5);// For Vapor-Liquid-Equilibrium [physical units]
    myCaseParameters.set<PHYS_CHAR_TEMPERATURE>(298.15); // For Vapor-Liquid-Equilibrium [physical units]
    //H2O, N2, O2
    myCaseParameters.set<FEED >({0.99, 0.0079, 0.0021});  // feed composition for Vapor-Liquid-Equilibrium as molar fraction
    myCaseParameters.set<parameters::A>({0.5995808, 0.1480650, 0.1506765});
    myCaseParameters.set<parameters::B>({1.8955853e-5, 2.4010114e-5, 1.9893672e-5});
    myCaseParameters.set<M                  >({0.01802, 0.02801, 0.03200});
    myCaseParameters.set<T_C                >({647.3, 126.2, 155.0});
    myCaseParameters.set<P_C                >({22089000, 3400000, 5040000});
    myCaseParameters.set<OMEGA_PER_COMPONENT>({0.34, 0.0377, 0.025});
    myCaseParameters.set<DEVI               >({0.867805648, 0.432399567, 0.41302780});
    //H2OH2O, H2ON2, H2OO2, N2H2O, N2N2, N2O2, O2H2O, O2N2, O2O2
    myCaseParameters.set<ALPHA>(
      {0., 0.199222317, 0.193233601, 0.199222317, 0., 0., 0.193233601, 0., 0.});
    myCaseParameters.set<G_I  >(
      {0., -4.29088111e3, -1.95640777e2, 3.02126911e4, 0., 5.65244934e2, 6.13396078e4, -5.01392189e2, 0.});
    myCaseParameters.set<G_II >(
      {0., 3.47847412e1, 2.10776021e1, -3.70834075e1, 0., 0., -1.22744109e2, 0., 0.});
    myCaseParameters.set<RHO0_L >({0., 0., 0.});
    myCaseParameters.set<RHO0_V >({0., 0., 0.});
    myCaseParameters.set<VLE    >({0., 0., 0., 0., 0., 0., 0., 0.});
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
