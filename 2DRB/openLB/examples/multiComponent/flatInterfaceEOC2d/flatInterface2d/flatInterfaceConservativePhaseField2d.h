/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Tim Bingert, Michael Rennick
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

/** flatInterfaceConservativePhaseField2d.h
 * Respective header for conseravtive phase field.
 */

#include <olb.h>
#include "../EOC/flatInterfaceErrors2d.h"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, D2Q9<RHO,NABLARHO,FORCE,EXTERNAL_FORCE,TAU_EFF,SCALAR>>,
  Component1,  Lattice<double, D2Q9<FORCE,VELOCITY,OLD_PHIU,STATISTIC,SCALAR>>
>;

using NSBulkDynamics = CorrectedMultiPhaseIncompressibleBGKdynamics<MyCase::value_t,MyCase::descriptor_t>;
using PFBulkDynamics = ConservativePhaseFieldBGKdynamics<MyCase::value_t,MyCase::descriptor_t_of<Component1>>;
using MixtureRules = LinearTauViscosity;
using Coupling = CorrectedLiangPostProcessor<MixtureRules>;

namespace olb::parameters {

struct RELAXATION_TIME_PF : public descriptors::FIELD_BASE<1> { };
struct C_RHO              : public descriptors::FIELD_BASE<1> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extent, origin);

  const T dx = extent[1] / params.get<parameters::RESOLUTION>();
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
  geometry.clean();
  geometry.innerClean();
  geometry.checkForErrors();
  geometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( MyCase& myCase )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice with conservative phase field model ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& sLatticeNS = myCase.getLattice(NavierStokes{});
  auto& sLatticePF = myCase.getLattice(Component1{});

  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const int N = params.get<parameters::RESOLUTION>();
  const T tau_mobil = params.get<parameters::RELAXATION_TIME_PF>();
  const T Re = params.get<parameters::REYNOLDS>();
  const T rho_g = params.get<parameters::RHO_VAPOR>();
  const T rho_l = params.get<parameters::RHO_LIQUID>();
  const T nu_g = params.get<parameters::NU_VAPOR>();
  const T nu_l = params.get<parameters::NU_LIQUID>();
  const T tau_l = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T tau_g = nu_g/nu_l*(tau_l-0.5)+0.5;
  const T sigma = params.get<parameters::SURFACE_TENSION>();
  const T w = params.get<parameters::INTERFACE_WIDTH>();
  const T C_rho = params.get<parameters::C_RHO>();
  const T theta = M_PI*45./180.;

  sLatticeNS.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,NSDESCRIPTOR>>(
    int   {N},               // resolution
    (T)   tau_l,             // lattice relaxation time
    (T)   extent[1],         // charPhysLength: reference length of simulation geometry
    (T)   Re/extent[1]*nu_l, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   nu_l,              // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   C_rho              // physDensity: physical density in __kg / m^3__
  );

  const auto& converter = sLatticeNS.getUnitConverter();
  converter.print();
  clout << "Physical Surface Tension: " << converter.getConversionFactorMass()/converter.getPhysDeltaT()/converter.getPhysDeltaT()*sigma << std::endl;

  sLatticePF.setUnitConverter(converter);

  // define lattice Dynamics
  dynamics::set<NSBulkDynamics>(sLatticeNS, geometry, 1);
  dynamics::set<PFBulkDynamics>(sLatticePF, geometry, 1);

  sLatticePF.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());

  auto& coupling = myCase.setCouplingOperator(
    "Coupling",
    Coupling{},
    names::NavierStokes{}, sLatticeNS,
    names::Component1{}, sLatticePF
  );

  auto zero = geometry.getMaterialIndicator(0);
  auto all  = geometry.getMaterialIndicator(1);
  auto wall = geometry.getMaterialIndicator(2);

  phasefield::init<Coupling,MixtureRules>(
    sLatticeNS, sLatticePF,  // both lattices
    coupling, all,           // phase field coupling, all material numbers
    wall, zero,              // wall and zero material numbers
    sigma, w, theta,         // phase field lattice parameters
    tau_g, tau_l,            // relaxation times for mixture rule
    rho_g/C_rho, rho_l/C_rho // lattice densities for mixture rule
  );
  coupling.setParameter<Coupling::GAMMA>( T(0) );
  coupling.setParameter<Coupling::GRAVITY>({T(0),T(0)});
  coupling.setParameter<Coupling::BODY_FORCE>({T(0),T(0)});

  sLatticeNS.setParameter<descriptors::OMEGA>( 1./tau_l );
  sLatticePF.setParameter<descriptors::OMEGA>( 1./tau_mobil );
  sLatticePF.setParameter<descriptors::INTERFACE_WIDTH>( w );

  {
    auto& communicator = sLatticePF.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(2);
    communicator.template requestField<STATISTIC>();
    communicator.template requestField<SCALAR>();
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  OstreamManager clout( std::cout,"setInitialValues" );
  clout << "Set Initial values with conservative phase field model ..." << std::endl;
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& sLatticeNS = myCase.getLattice(NavierStokes{});
  auto& sLatticePF = myCase.getLattice(Component1{});

  const auto& converter = sLatticeNS.getUnitConverter();

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const int N = params.get<parameters::RESOLUTION>();
  const int phaseLength = N/2;
  const T dx = extent[1] / N;
  const T Re = params.get<parameters::REYNOLDS>();
  const T rho_g = params.get<parameters::RHO_VAPOR>();
  const T rho_l = params.get<parameters::RHO_LIQUID>();
  const T nu_g = params.get<parameters::NU_VAPOR>();
  const T nu_l = params.get<parameters::NU_LIQUID>();
  const T tau_l = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T tau_g = nu_g/nu_l*(tau_l-0.5)+0.5;
  const T phys_pressure = params.get<parameters::PHYS_CHAR_PRESSURE>();
  const T w = params.get<parameters::INTERFACE_WIDTH>();
  const T C_rho = params.get<parameters::C_RHO>();

  Vector<T,2> u_l(0., Re/extent[1]*nu_l/converter.getConversionFactorVelocity());
  AnalyticalConst2D<T,T> zeroVelocity_l( u_l );
  Vector<T,2> u_p(0., Re/extent[1]*nu_l);
  AnalyticalConst2D<T,T> zeroVelocity_p( u_p );

  AnalyticalConst2D<T,T> one ( 1. );
  AnalyticalConst2D<T,T> zero ( 0. );
  T lattice_pressure = phys_pressure/converter.getConversionFactorPressure();
  AnalyticalConst2D<T,T> pressure_l ( lattice_pressure );
  AnalyticalConst2D<T,T> pressure_p ( phys_pressure );
  AnalyticalConst2D<T,T> rhov ( rho_g/C_rho );
  AnalyticalConst2D<T,T> rhol ( rho_l/C_rho );
  AnalyticalConst2D<T,T> tauv ( tau_g );
  AnalyticalConst2D<T,T> taul ( tau_l );

  SmoothIndicatorFactoredCuboid2D<T,T> interfaceAC( {extent[0]/2., (N/2.+0.5)*dx}, 0, phaseLength*dx, w*dx/2, 0, {0,0}, 0, -1. );
  AnalyticalIdentity2D<T,T> phi( one + interfaceAC );
  AnalyticalIdentity2D<T,T> rho( rhov + (rhol-rhov)*phi );
  AnalyticalIdentity2D<T,T> tau( tauv + (taul-tauv)*phi );

  sLatticeNS.defineField<descriptors::RHO>( geometry, 1, rho );
  sLatticeNS.defineField<descriptors::TAU_EFF>( geometry, 1, tau );
  sLatticeNS.defineField<descriptors::SCALAR>( geometry, 1, one );
  sLatticePF.defineField<descriptors::OLD_PHIU>( geometry, 1, zeroVelocity_l );

  auto all = geometry.getMaterialIndicator({0,1});

  momenta::setIncompressiblePressure(sLatticeNS, all, pressure_p);
  momenta::setVelocity(sLatticeNS, all, zeroVelocity_p);
  sLatticeNS.iniEquilibrium( geometry, 1, pressure_l, zeroVelocity_l );

  momenta::setOrderParameter(sLatticePF, all, phi);
  momenta::setVelocity(sLatticePF, all, zeroVelocity_p);

  sLatticeNS.initialize();
  sLatticePF.initialize();

  sLatticeNS.iniEquilibrium( geometry, 1, pressure_l, zeroVelocity_l );
  sLatticePF.iniEquilibrium( geometry, 1, phi, zeroVelocity_l );

  sLatticeNS.setProcessingContext(ProcessingContext::Simulation);
  sLatticePF.setProcessingContext(ProcessingContext::Simulation);

  sLatticePF.executePostProcessors(stage::PreCoupling());
  sLatticePF.getCommunicator(stage::PreCoupling()).communicate();
  clout << "Set Initial Velues ... OK" << std::endl;
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{ }

void getResults( MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT )
{
  using T = MyCase::value_t;
  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using PFDESCRIPTOR = MyCase::descriptor_t_of<Component1>;
  auto& params = myCase.getParameters();

  auto& sLatticeNS = myCase.getLattice(NavierStokes{});
  auto& sLatticePF = myCase.getLattice(Component1{});

  const auto& converter = sLatticeNS.getUnitConverter();

  const std::size_t statIter = converter.getLatticeTime(2.5e-4);
  const std::size_t saveIter = converter.getLatticeTime(2.5e-4);

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const int N = params.get<parameters::RESOLUTION>();
  const int phaseLength = N/2;
  const T dx = extent[1] / N;
  const T w = params.get<parameters::INTERFACE_WIDTH>();

  SuperVTMwriter2D<T> vtkWriter( "flatInterfaceConservativePhaseField2d" );

  SmoothIndicatorFactoredCuboid2D<T,T> interface_diffuse( {extent[0]/2., (N/2.+0.5)*dx}, 0, phaseLength*dx, w*dx/2, 0, {0,0}, 0, -1. );
  if ( iT==0 ) {
    SuperLatticeCuboid2D cuboid( sLatticeNS );
    SuperLatticeRank2D rank( sLatticeNS );
    vtkWriter.write( cuboid );
    vtkWriter.write( rank );
    vtkWriter.createMasterFile();
  }

  if ( iT%statIter==0 ) {
    timer.update( iT );
    timer.printStep();
    sLatticeNS.getStatistics().print( iT, converter.getPhysTime(iT) );
    sLatticePF.getStatistics().print( iT, converter.getPhysTime(iT) );
  }

  // Writes the VTK files
  if ( iT%saveIter==0 ) {
    sLatticeNS.setProcessingContext(ProcessingContext::Evaluation);
    sLatticePF.setProcessingContext(ProcessingContext::Evaluation);
    SuperLatticePhysIncPressure2D<T, NSDESCRIPTOR> p_hydro( sLatticeNS, converter );
    p_hydro.getName() = "p_hydro";
    SuperLatticeDensity2D<T, PFDESCRIPTOR> phi( sLatticePF );
    phi.getName() = "phi";
    SuperLatticeField2D<T, NSDESCRIPTOR, RHO> rho_L( sLatticeNS );
    AnalyticalConst2D<T,T> ConversionDensity_( converter.getConversionFactorDensity() );
    SuperLatticeFfromAnalyticalF2D<T, NSDESCRIPTOR> ConversionDensity(ConversionDensity_, sLatticeNS);
    SuperIdentity2D<T,T> rho( rho_L * ConversionDensity );
    rho.getName() = "rho";
    SuperLatticePhysVelocity2D<T, NSDESCRIPTOR> velocity( sLatticeNS, converter );
    velocity.getName() = "u";

    AnalyticalConst2D<T,T> one ( 1. );
    AnalyticalIdentity2D<T,T> phi_ana( one + interface_diffuse );
    SuperLatticeFfromAnalyticalF2D<T, PFDESCRIPTOR> phi_analytical_diff(phi_ana, sLatticePF);
    phi_analytical_diff.getName() = "phi_analytical_diff";
    vtkWriter.addFunctor( phi_analytical_diff );

    AnalyticalConst2D<T,T> zero ( 0. );
    Vector<T,2> extend( 2./100.*extent[1], 0.25*extent[1] );
    Vector<T,2> origin1( 0., 0.);
    Vector<T,2> origin2( 0., 3./4.*extent[1]);
    IndicatorCuboid2D<T> ind1( extend, origin1 );
    IndicatorCuboid2D<T> ind2( extend, origin2 );
    SmoothIndicatorCuboid2D<T,T> interface_sharp1( ind1, T(0) );
    SmoothIndicatorCuboid2D<T,T> interface_sharp2( ind2, T(0) );
    AnalyticalIdentity2D<T,T> phi_ana_sharp( zero + interface_sharp1 + interface_sharp2 );
    SuperLatticeFfromAnalyticalF2D<T, PFDESCRIPTOR> phi_analytical_sharp(phi_ana_sharp, sLatticePF);
    phi_analytical_sharp.getName() = "phi_analytical_sharp";
    vtkWriter.addFunctor( phi_analytical_sharp );

    vtkWriter.addFunctor( p_hydro );
    vtkWriter.addFunctor( phi );
    vtkWriter.addFunctor( rho );
    vtkWriter.addFunctor( velocity );
    vtkWriter.write( iT );
  }
}

void simulate(MyCase& myCase, bool Cahn_const) {
  OstreamManager clout( std::cout,"simulate" );

  using T = MyCase::value_t;

  auto& sLatticeNS = myCase.getLattice(NavierStokes{});
  auto& sLatticePF = myCase.getLattice(Component1{});

  const auto& converter = sLatticeNS.getUnitConverter();

  const T physCheckInterval = 2.5e-4;                                 // [s]
  const std::size_t checkIter = converter.getLatticeTime(physCheckInterval);
  const std::size_t iTmax     = converter.getLatticeTime(2e-2);

  const T tol           = 1e-5;
  const T latticeUFloor = 1e-8;   // below this, u is just round-off
  const int hitsNeeded  = 3;
  int hits = 0;

  T orderParameterL1RelError_old = 1.;

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  myCase.getOperator("Coupling").apply();
  for ( std::size_t iT=0; iT<=iTmax; ++iT ) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    sLatticeNS.collideAndStream();
    sLatticePF.collideAndStream();

    sLatticePF.executePostProcessors(stage::PreCoupling());
    sLatticePF.getCommunicator(stage::PreCoupling()).communicate();
    myCase.getOperator("Coupling").apply();

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, timer, iT);
    if ( std::isnan( sLatticeNS.getStatistics().getAverageEnergy() ) ) {
      break;
    }
    if (iT % checkIter == 0) {
      sLatticeNS.setProcessingContext(ProcessingContext::Evaluation);
      sLatticePF.setProcessingContext(ProcessingContext::Evaluation);
      errorFlatInterface(myCase, Cahn_const);

      if (iT > 0) {
        const T dPhi = util::abs(orderParameterL1RelError - orderParameterL1RelError_old)
                     / util::abs(orderParameterL1RelError);
        const T uLat = velocityLinfAbsError / converter.getConversionFactorVelocity();

        if (dPhi < tol && uLat < latticeUFloor) { ++hits; } else { hits = 0; }
        if (hits >= hitsNeeded) {
          clout << "Converged at physT=" << converter.getPhysTime(iT) << std::endl;
          break;
        }
      }
      orderParameterL1RelError_old = orderParameterL1RelError;
    }
  }
  timer.stop();
  timer.printSummary();
}