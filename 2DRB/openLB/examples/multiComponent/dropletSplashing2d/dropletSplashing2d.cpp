/*Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Luiz Eduardo Czelusniak
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

/* dropletSplashing2d.cpp
 * In this example a droplet falls in a liquid film.
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, D2Q9<FORCE,EXTERNAL_FORCE,STATISTIC,SCALAR>>
>;

using BulkDynamics = MultiphaseForcedBGKdynamics<MyCase::value_t, MyCase::descriptor_t>;
using COUPLING = PseudopotentialForcedPostProcessor<interaction::Polinomial>;

namespace olb::parameters {

struct U_DROPLET : public descriptors::FIELD_BASE<1> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  T radius = params.get<parameters::RADIUS>();
  Vector extent(8*radius, 4*radius);
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extent, origin);

  const T physDeltaX = extent[0] / params.get<parameters::RESOLUTION>();
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({true,false});
  return mesh;
}

/**
 * Correction for initial population
 * Take into account the force effect in initial population
 */
struct InitialPopulationCorrectionO {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform
  {
    using T = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    const auto force = cell.template getField<descriptors::FORCE>();
    const T rho = cell.computeRho();

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      T c_F{};
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        c_F += descriptors::c<DESCRIPTOR>(iPop,iD)*force[iD];
      }
      c_F *= descriptors::invCs2<T,DESCRIPTOR>();
      cell[iPop] -= descriptors::t<T,DESCRIPTOR>(iPop) * 0.5 * rho * c_F;
    }

  }
};

void prepareGeometry( MyCase& myCase )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;
  auto& geometry = myCase.getGeometry();
  geometry.rename( 0,2 );
  geometry.rename( 2,1,{0,1} );
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

  auto& sLattice = myCase.getLattice(NavierStokes{});

  // TODO for now, to be combined with unit converter refactor
  const int Nx = params.get<parameters::RESOLUTION>();
  const T tau = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T radius = params.get<parameters::RADIUS>();
  const T Lx = 8*radius;
  const T nu_vapor = params.get<parameters::NU_VAPOR>();
  const T nu_liquid = params.get<parameters::NU_LIQUID>();
  const T rho_vapor = params.get<parameters::RHO_VAPOR>();
  const T rho_liquid = params.get<parameters::RHO_LIQUID>();
  const T sigma = params.get<parameters::SURFACE_TENSION>();
  const T w = params.get<parameters::INTERFACE_WIDTH>();

  // Unit Converter
  sLattice.setUnitConverter<MultiPhaseUnitConverterFromRelaxationTime<T,DESCRIPTOR>>(
    (T)   Nx,                // resolution
    (T)   tau,               // lattice relaxation time
    (T)   rho_liquid/1000.,  // lattice density
    (T)   Lx,                // charPhysLength: reference length of simulation geometry
    (T)   nu_liquid,         // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   rho_liquid         // physDensity: physical density in __kg / m^3__
  );
  const auto& converter = sLattice.getUnitConverter();
  converter.print();

  // define lattice Dynamics
  dynamics::set<BulkDynamics>(sLattice, geometry, 1);
  boundary::set<boundary::BounceBack>(sLattice, geometry, 2);

  // global relaxation frequency (it can be initialized as one)
  AnalyticalConst2D<T,T> one( 1. );
  sLattice.setParameter<OMEGA>( 1. );
  sLattice.setParameter<multiphase::RHO_VAPOR>( rho_vapor/converter.getConversionFactorDensity() );
  sLattice.setParameter<multiphase::RHO_LIQUID>( rho_liquid/converter.getConversionFactorDensity() );

  sLattice.setParameter<multiphase::OMEGA_VAPOR>( 1./converter.computeRelaxationTimefromPhysViscosity( nu_vapor ) );
  sLattice.setParameter<multiphase::OMEGA_LIQUID>( 1./converter.computeRelaxationTimefromPhysViscosity( nu_liquid ) );

  sLattice.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());

  {
    auto& communicator = sLattice.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(1);
    communicator.requestField<STATISTIC>();
    communicator.exchangeRequests();
  }

  auto& couplingInteractionForce = myCase.setCouplingOperator(
    "Force_coupling",
    COUPLING{},
    names::Component1{}, sLattice
  );
  couplingInteractionForce.template setParameter<interaction::Polinomial::RHOV>(
                                                  rho_vapor/converter.getConversionFactorDensity());
  couplingInteractionForce.template setParameter<interaction::Polinomial::RHOL>(
                                                  rho_liquid/converter.getConversionFactorDensity());
  couplingInteractionForce.template setParameter<interaction::Polinomial::THICKNESS>(w);
  couplingInteractionForce.template setParameter<interaction::Polinomial::SURFTENSION>(
                                                  sigma/converter.getConversionFactorSurfaceTension());

  // Compute the interaction parameters
  interaction::Polinomial::computeParameters<T>(couplingInteractionForce);

  // Display the value of surface tension parameter
  // recommended kappaP not much larger than 1
  auto kappaP = couplingInteractionForce.template getParameter<interaction::Polinomial::KAPPAP>();
  clout << "Surface tension parameter: " << kappaP[0] << std::endl;

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& sLattice = myCase.getLattice(NavierStokes{});
  const auto& converter = sLattice.getUnitConverter();

  const T radius = params.get<parameters::RADIUS>();
  const T Lx = 8*radius;
  const T Ly = 4*radius;
  const T U_droplet = params.get<parameters::U_DROPLET>();
  const T rho_vapor = params.get<parameters::RHO_VAPOR>();
  const T rho_liquid = params.get<parameters::RHO_LIQUID>();
  const T w = params.get<parameters::INTERFACE_WIDTH>();

  std::vector<T> dropletVelocity_x = {0.};
  AnalyticalConst2D<T,T> _dropletVelocity_x( dropletVelocity_x );

  SmoothIndicatorFactoredCircle2D<T,T> dropletVelocity_y( {Lx/2., Ly/2.}, radius,
                                                           sqrt(2.)*w*converter.getConversionFactorLength(),
                                                           0, {0,0}, 0,
                                                           -U_droplet);
  AnalyticalIdentity2D<T,T> _dropletVelocity_y( dropletVelocity_y );

  std::shared_ptr<AnalyticalF2D<T,T>> fluidVelocity( new AnalyticalComposed2D<T,T>(_dropletVelocity_x,_dropletVelocity_y));


  std::shared_ptr<AnalyticalF2D<T,T>> vapor ( new AnalyticalConst2D<T,T>(rho_vapor));
  std::shared_ptr<AnalyticalF2D<T,T>> liquid( new SmoothIndicatorFactoredCircle2D<T,T>(
                                                 {Lx/2., Ly/2.}, radius,
                                                 sqrt(2.)*w*converter.getConversionFactorLength(),
                                                 0, {0,0}, 0,
                                                 ( rho_liquid - rho_vapor )));
  std::shared_ptr<AnalyticalF2D<T,T>> film ( new SmoothIndicatorFactoredCuboid2D<T,T>(
                                                 {Lx/2., 0.}, 2.*Lx, radius,
                                                 sqrt(2.)*w*converter.getConversionFactorLength(),
                                                 0, {0,0}, 0,
                                                 ( rho_liquid - rho_vapor )));

  std::shared_ptr<AnalyticalF2D<T,T>> fluidDensity( vapor + liquid + film );

  auto bulkIndicator = geometry.getMaterialIndicator({1,2});
  momenta::setVelocity(sLattice, bulkIndicator, *fluidVelocity);
  momenta::setDensity(sLattice, bulkIndicator, *fluidDensity);

  std::shared_ptr<AnalyticalF2D<T,T>> latticeFluidDensity( fluidDensity / converter.getConversionFactorDensity());
  std::shared_ptr<AnalyticalF2D<T,T>> latticeFluidVelocity( fluidVelocity / converter.getConversionFactorVelocity());
  sLattice.iniEquilibrium( bulkIndicator, *latticeFluidDensity, *latticeFluidVelocity );

  std::vector<T> fnull( 2,T() );
  fields::set<descriptors::EXTERNAL_FORCE>(sLattice, geometry.getMaterialIndicator(1), fnull);

  sLattice.initialize();
  sLattice.getCommunicator(stage::PreCoupling()).communicate();
  sLattice.executePostProcessors(stage::PreCoupling());

  // Compute initial force
  sLattice.executePostProcessors(stage::PreCoupling());
  sLattice.getCommunicator(stage::PreCoupling()).communicate();
  myCase.getOperator("Force_coupling").apply();

  // Correct velocity
  sLattice.addPostProcessor<stage::PostCoupling>(meta::id<InitialPopulationCorrectionO>());
  sLattice.executePostProcessors(stage::PostCoupling());
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{ }

//std::vector<T>
void getResults( MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT )
{
  OstreamManager clout( std::cout,"getResults" );

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;

  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& params = myCase.getParameters();
  const auto& converter = sLattice.getUnitConverter();

  SuperVTMwriter2D<T> vtkWriter( "dropletSplashing2d" );

  const int statIter = params.get<parameters::LATTICE_STAT_ITER_T>();
  const int saveIter = params.get<parameters::LATTICE_VTK_ITER_T>();

  if ( iT==0 ) {
    SuperLatticeCuboid2D cuboid( sLattice );
    SuperLatticeRank2D rank( sLattice );
    vtkWriter.write( cuboid );
    vtkWriter.write( rank );
    vtkWriter.createMasterFile();
  }

  if ( iT%statIter==0 ) {
    timer.update( iT );
    timer.printStep();
    sLattice.getStatistics().print( iT, converter.getPhysTime(iT) );
  }

  // Save vtk files
  if ( iT%saveIter==0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    // rho_lat -> density in lattice units
    SuperLatticeDensity2D rho_lat( sLattice );
    rho_lat.getName() = "rho_lat";

    // rho_phs -> density in physical units
    AnalyticalConst2D<T,T> _C_rho( converter.getConversionFactorDensity() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> __C_rho( _C_rho, sLattice );
    SuperIdentity2D<T,T> rho_phs( rho_lat * __C_rho );
    rho_phs.getName() = "rho_phs";

    // velocity_lat -> velocity in lattice units
    SuperLatticeVelocity2D<T, DESCRIPTOR> velocity_lat( sLattice );
    velocity_lat.getName() = "velocity_lat";

    // velocity_phs -> velocity in physical units
    AnalyticalConst2D<T,T> _C_U( converter.getConversionFactorVelocity() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> __C_U( _C_U, sLattice );
    SuperIdentity2D<T,T> velocity_phs( velocity_lat * __C_U );
    velocity_phs.getName() = "velocity_phs";

    // force_lat -> force in lattice units
    SuperLatticeField2D<T, DESCRIPTOR, FORCE> force_lat( sLattice );
    force_lat.getName() = "force_lat";

    // force_phs -> force in physical units
    AnalyticalConst2D<T,T> _C_F( converter.getConversionFactorForce() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> __C_F( _C_F, sLattice );
    SuperIdentity2D<T,T> force_phs( force_lat * __C_F );
    force_phs.getName() = "force_phs";

    // p_lat -> pressure in lattice units
    SuperLatticeField2D<T, DESCRIPTOR, SCALAR> p_lat( sLattice );
    p_lat.getName() = "p_lat";

    // p_phs -> physical in physical units
    AnalyticalConst2D<T,T> _C_P( converter.getConversionFactorPressure() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> __C_P( _C_P, sLattice );
    SuperIdentity2D<T,T> p_phs( p_lat * __C_P );
    p_phs.getName() = "p_phs";

    // omega
    SuperLatticeField2D<T, DESCRIPTOR, OMEGA> omega( sLattice );
    omega.getName() = "Omega";

    vtkWriter.addFunctor( rho_lat );
    vtkWriter.addFunctor( rho_phs );
    vtkWriter.addFunctor( velocity_lat );
    vtkWriter.addFunctor( velocity_phs );
    vtkWriter.addFunctor( force_lat );
    vtkWriter.addFunctor( force_phs );
    vtkWriter.addFunctor( p_lat );
    vtkWriter.addFunctor( p_phs );
    vtkWriter.addFunctor( omega );
    vtkWriter.write( iT );
  }
}

void simulate(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& params = myCase.getParameters();
  const std::size_t iTmax = params.get<parameters::LATTICE_MAX_T>();
  auto& sLattice = myCase.getLattice(NavierStokes{});

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for ( std::size_t iT=0; iT<=iTmax; ++iT ) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.3: Computation and Output of the Results ===
    getResults( myCase, timer, iT );

    /// === Step 8.2: Collide and Stream Execution ===
    sLattice.collideAndStream();

    sLattice.executePostProcessors(stage::PreCoupling());
    sLattice.getCommunicator(stage::PreCoupling()).communicate();
    myCase.getOperator("Force_coupling").apply();
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
    myCaseParameters.set<RESOLUTION             >(400);      // Nx [lattice units]
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.6);    // relaxation time [lattice units]
    myCaseParameters.set<parameters::RADIUS>(5e-6);     // radius of liquid phase [m]
    myCaseParameters.set<U_DROPLET         >(6.);       // droplet velocity [m.s-1]
  // Properties of R134 vapor and liquid
  // at saturation condition 30 degree Celsius
    myCaseParameters.set<NU_VAPOR       >(3.39e-7);  // vapor kinematic viscosity [m2.s-1]
    myCaseParameters.set<NU_LIQUID      >(1.58e-7);  // liquid kinematic viscosity [m2.s-1]
    myCaseParameters.set<RHO_VAPOR      >(37.5);     // vapor density [kg.m-3]
    myCaseParameters.set<RHO_LIQUID     >(1187);     // liquid density [kg.m-3]
    myCaseParameters.set<SURFACE_TENSION>(7.58e-3);// surface tension for water-air [N.m-1]
  // Lattice parameters
    myCaseParameters.set<parameters::INTERFACE_WIDTH>(1.);// interface thickness in grid nodes; physThickness = delta_x * thickness
    myCaseParameters.set<LATTICE_MAX_T      >(1000);   // max iterations [lattice units]
    myCaseParameters.set<LATTICE_STAT_ITER_T>(10);     // statistics iterations [lattice units]
    myCaseParameters.set<LATTICE_VTK_ITER_T >(10);     // vtk iterations density [lattice units]
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
