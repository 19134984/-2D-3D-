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

/* isothermalEvaporation2d.cpp
 * Flat interface evaporates due to a lower pressure assigned
 * to the boundary. This example is runned in lattice units, but
 * can be easily changed to physical units by prescribing the values
 * of the parameters in physical units.
 */

// Warning! This example is not numerically stable, i.e. it does not converge every time it is run.
// TODO: should be fixed in later versions
#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::descriptors;
using namespace olb::graphics;

using MyCase = Case<NavierStokes, Lattice<double, descriptors::D2Q9<STATISTIC,      // Store density field
                                                                    CHEM_POTENTIAL, // Chemical potential for multiphase system
                                                                    FORCE,          // Thermodynamic forced based on chemical potential
                                                                    EXTERNAL_FORCE, // Other forces such as gravity (optional in this case)
                                                                    OMEGA,          // Store density dependent relaxation frequency
                                                                    SCALAR>>>;      // Store laplacian of the density


using COUPLING1 = ChemicalPotentialPostProcessor<EOS::Landau>;
using COUPLING2 = FreeEnergyForcedPostProcessor;

/// @brief Create a simulation mesh, based on user-specified geometry
/// @return An instance of a mesh with the relevant information
Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;

  Vector<T,2> extend(parameters.get<parameters::DOMAIN_EXTENT>()[0], parameters.get<parameters::DOMAIN_EXTENT>()[1]);
  Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid(extend, origin);

  T physDeltaX = parameters.get<parameters::DOMAIN_EXTENT>()[0] / parameters.get<parameters::RESOLUTION>();

  Mesh<T,MyCase::d> mesh = Mesh<T, MyCase::d>(cuboid, physDeltaX , singleton::mpi().getSize());
  mesh.setOverlap(3);
  mesh.getCuboidDecomposition().setPeriodicity({false,true});
  return mesh;
}

/// @brief Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
void prepareGeometry(MyCase& myCase)
{
  auto& geometry = myCase.getGeometry();
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;
  geometry.rename( 0,2 );
  geometry.rename( 2,1,{1,0} );
  geometry.innerClean();
  geometry.checkForErrors();
  geometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

/**
 * Correction for lattice velocity
 * Lattice velocity Ulat in this model is not the real velocity
 * Ureal = Ulat + 0.5*F
 * We compute Ulat = Ureal - 0.5*F as initial condition
 */
struct initialVelocityCorrection  {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform
  {
    using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
    using T = MyCase::value_t;

    T U[descriptors::d<DESCRIPTOR>()];
    const auto force = cell.template getField<descriptors::FORCE>();

    cell.computeU(U);

    for(int iD = 0; iD < DESCRIPTOR::d; ++iD){
      U[iD] = U[iD] - 0.5*force[iD];
    }

    cell.defineU(U);

  }
};

/// @brief Set lattice dynamics
/// @param myCase The Case instance which keeps the simulation data
void prepareLattice(MyCase& myCase)
{
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const int Nx = parameters.get<parameters::RESOLUTION>();
  const T tau = parameters.get<parameters::LATTICE_RELAXATION_TIME>();
  const T rho_liquid = parameters.get<parameters::RHO_LIQUID>();
  const T Lx = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T nu = parameters.get<parameters::NU>();

  // Unit Converter
  lattice.setUnitConverter<MultiPhaseUnitConverterFromRelaxationTime<T,DESCRIPTOR>>(
    (T)   Nx,                        // resolution
    (T)   tau,                       // lattice relaxation time
    (T)   rho_liquid,                // lattice density
    (T)   Lx,                        // charPhysLength: reference length of simulation geometry
    (T)   nu,                        // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   rho_liquid                 // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  lattice.getUnitConverter().print();
  const auto& converter = lattice.getUnitConverter();

  // define lattice Dynamics
  dynamics::set<ForcedWagnerBGKdynamics>(lattice, geometry.getMaterialIndicator({1, 2}));

  clout << "Check Point" << std::endl;

  // Set boundary pressure condition
  boundary::set<boundary::LocalPressure<T, DESCRIPTOR, ForcedWagnerBGKdynamics<T, DESCRIPTOR>>>(lattice, geometry, 2);

  std::vector<T> zeroV = {0., 0.};
  AnalyticalConst2D<T,T> zeroVelocity( zeroV );

  // Computing the density
  const T rho_V = parameters.get<parameters::RHO_VAPOR_ANALYTICAL>();
  const T Ly = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const T length = parameters.get<parameters::LIQUID_PHASE_LENGTH>();
  const T thickness = parameters.get<parameters::THICKNESS>();
  const T rho_L = parameters.get<parameters::RHO_LIQUID_ANALYTICAL>();
  const T Ui_sol = parameters.get<parameters::INTERFACE_VELOCITY_ANALYTICAL>();

  std::shared_ptr<AnalyticalF2D<T,T>> liquid = std::make_shared<SmoothIndicatorFactoredCuboid2D<T,T>>(
                                                 Vector<T,2>(Lx/2., Ly/2.), length, Lx,
                                                 sqrt(2.)*thickness*converter.getConversionFactorLength(),
                                                 0, Vector<T,2>(0, 0), 0,
                                                 ( rho_L - rho_V ));
  std::shared_ptr<AnalyticalF2D<T,T>> rho(rho_V + liquid);
  std::shared_ptr<AnalyticalF2D<T,T>> latticeRho(rho / converter.getConversionFactorDensity());

  // sign function
  std::shared_ptr<AnalyticalF2D<T,T>> right = std::make_shared<SmoothIndicatorFactoredCuboid2D<T,T>>(
                                                 Vector<T,2>(Lx, Ly/2.), Lx, Lx,
                                                 sqrt(2.)*thickness*converter.getConversionFactorLength(),
                                                 0, Vector<T,2>(0, 0), 0,
                                                 2.);
  std::shared_ptr<AnalyticalF2D<T,T>> sign( -1. + right );

  // Computing the velocity
  std::shared_ptr<AnalyticalF2D<T,T>> solVelocity(sign * Ui_sol * rho_L / rho - (sign * Ui_sol));
  std::shared_ptr<AnalyticalF2D<T,T>> latticeSolVelocity(solVelocity / converter.getConversionFactorVelocity());

  momenta::setVelocity(lattice, geometry.getMaterialIndicator({1, 2}), *solVelocity);
  momenta::setDensity(lattice, geometry.getMaterialIndicator({1, 2}), *rho);
  lattice.iniEquilibrium(geometry.getMaterialIndicator({1, 2}), *latticeRho, *latticeSolVelocity);

  // Set Chemical Potential
  const T chemical = parameters.get<parameters::BOUNDARY_CHEMICAL_POTENTIAL>() / converter.getConversionFactorChemicalPotential();
  fields::set<descriptors::CHEM_POTENTIAL>(lattice, geometry.getMaterialIndicator(2), chemical);

  // Set OMEGA
  fields::set<descriptors::OMEGA>(lattice, geometry.getMaterialIndicator(1), 1./ converter.computeRelaxationTimefromPhysViscosity(nu));

  lattice.addPostProcessor<stage::PreCoupling>(meta::id<RhoStatistics>());

  {
    auto& communicator = lattice.getCommunicator(stage::PreCoupling());
    communicator.requestOverlap(1);
    communicator.requestField<STATISTIC>();
    communicator.requestField<CHEM_POTENTIAL>();
    communicator.exchangeRequests();
  }

  lattice.initialize();
  lattice.getCommunicator(stage::PreCoupling()).communicate();
  lattice.executePostProcessors(stage::PreCoupling());

  // Set chemical and force couplings
  auto& chemicalCoupling = myCase.setCouplingOperator(
    "Chemical",
    COUPLING1{},
    names::Component1{}, lattice
  );

  chemicalCoupling.restrictTo(geometry.getMaterialIndicator({1}));

  const T rho_vapor = parameters.get<parameters::RHO_VAPOR>();
  const T surfTension = parameters.get<parameters::SURFACE_TENSION>();
  chemicalCoupling.template setParameter<EOS::Landau::RHOV>(
                  rho_vapor/converter.getConversionFactorDensity());
  chemicalCoupling.template setParameter<EOS::Landau::RHOL>(
                  rho_liquid/converter.getConversionFactorDensity());
  chemicalCoupling.template setParameter<EOS::Landau::THICKNESS>(thickness);
  chemicalCoupling.template setParameter<EOS::Landau::SURFTENSION>(
                  surfTension/converter.getConversionFactorSurfaceTension());

  // Compute the EOS parameters
  EOS::Landau::computeParameters<T>(chemicalCoupling);

  auto& forceCoupling = myCase.setCouplingOperator(
    "Force",
    COUPLING2{},
    names::Component1{}, lattice
  );

  forceCoupling.restrictTo(geometry.getMaterialIndicator({1}));

  // Initial condition
  // Compute force
  lattice.executePostProcessors(stage::PreCoupling());
  lattice.getCommunicator(stage::PreCoupling()).communicate();
  chemicalCoupling.apply();
  lattice.executePostProcessors(stage::PreCoupling());
  lattice.getCommunicator(stage::PreCoupling()).communicate();
  forceCoupling.apply();

  // Correct velocity
  lattice.addPostProcessor<stage::PostCoupling>(meta::id<initialVelocityCorrection>());
  lattice.executePostProcessors(stage::PostCoupling());

  clout << "Prepare Lattice ... OK" << std::endl;
}


void boundaryCondition(MyCase& myCase)
{
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  const T rho_V = parameters.get<parameters::RHO_VAPOR_ANALYTICAL>();
  const T mu_b = parameters.get<parameters::BOUNDARY_CHEMICAL_POTENTIAL>();

  // Compute density at MN1
  SuperLatticeDensity2D<T, DESCRIPTOR> _density( lattice );
  AnalyticalFfromSuperF2D<T,T> interpolRHO( _density, true, 1 );
  T dx = lattice.getUnitConverter().getConversionFactorLength();
  T pos_MN1[2] = {dx, 0.};
  T lattice_rho_MN1;
  interpolRHO( &lattice_rho_MN1, pos_MN1 );

  // Set boundary density at MN2
  T rho_MN1 = lattice.getUnitConverter().getPhysDensity(lattice_rho_MN1);
  T rho_MN2 = (0.5 * rho_V) + (0.5 * rho_MN1) ;

  momenta::setDensity(lattice, geometry.getMaterialIndicator(2), rho_MN2);
  fields::set<descriptors::CHEM_POTENTIAL>(
    lattice, geometry.getMaterialIndicator(2), mu_b / lattice.getUnitConverter().getConversionFactorChemicalPotential()
  );
}

/**
 * Variable relaxation frequency increase stability and
 * decrease oscillations due to initialization
 * Bulk: omega = 1
 * Interface: omega = 1/tau
*/
void computeRelaxationFrequency(MyCase& myCase)
{
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  // Get density
  SuperLatticeDensity2D<T, DESCRIPTOR> _density( lattice );

  // Find interface position
  const T rho_V = parameters.get<parameters::RHO_VAPOR_ANALYTICAL>();
  const T rho_L = parameters.get<parameters::RHO_LIQUID_ANALYTICAL>();
  T rhom = 0.5*(rho_V+rho_L)/lattice.getUnitConverter().getConversionFactorDensity();
  T dx = lattice.getUnitConverter().getConversionFactorLength();
  AnalyticalFfromSuperF2D<T,T> interpolRHO( _density, true, 1 );
  T pos1[2] = {0., 0.};
  T pos2[2] = {0., 0.};
  T h1 = 0, h2 = 0;
  const T Nx = parameters.get<parameters::RESOLUTION>();
  for (int ix=1; ix<Nx-1; ix++) {
    pos1[0] = ix * dx;
    pos2[0] = ( ix + 1 ) * dx;
    T rho1, rho2;
    interpolRHO( &rho1, pos1 );
    interpolRHO( &rho2, pos2 );
    if ( rho1 < rhom && rho2 > rhom )
      h1 = ix * dx;
    if ( rho1 > rhom && rho2 < rhom )
      h2 = ix * dx;
  }

  // Computing indicator
  const T Lx = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T Ly = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const T thickness = parameters.get<parameters::THICKNESS>();
  const T tau = parameters.get<parameters::LATTICE_RELAXATION_TIME>();

  const T bulk = 1./tau;
  std::shared_ptr<AnalyticalF2D<T,T>> left = std::make_shared<SmoothIndicatorFactoredCuboid2D<T,T>>(
    Vector<T,2>(0., Ly/2.), 2.*h1-15.*thickness*dx,
    Lx, 1., 0,
    Vector<T,2>(0, 0), 0,
    1. - 1./tau
  );
  std::shared_ptr<AnalyticalF2D<T,T>> center = std::make_shared<SmoothIndicatorFactoredCuboid2D<T,T>>(
    Vector<T,2>(Lx/2., Ly/2.), h2-h1-15.*thickness*dx,
    Lx, 1., 0,
    Vector<T,2>(0, 0), 0,
    1. - 1./tau
  );
  std::shared_ptr<AnalyticalF2D<T,T>> right = std::make_shared<SmoothIndicatorFactoredCuboid2D<T,T>>(
    Vector<T,2>(Lx, Ly/2.), 2.*h1-15.*thickness*dx,
    Lx, 1., 0,
    Vector<T,2>(0, 0), 0,
    1. - 1./tau
  );
  std::shared_ptr<AnalyticalF2D<T,T>> _omega ( bulk + left + center + right );

  // Setting Omega field
  fields::set<descriptors::OMEGA>(lattice, geometry.getMaterialIndicator({1, 2}), *_omega);
}

/// Compute simulation results at times
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step, timer
void getResults(MyCase& myCase, int iT, util::Timer<MyCase::value_t>& timer)
{
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& parameters = myCase.getParameters();

  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter2D<T> vtmWriter( "isothermalEvaporation2D" );

  const int vtkIter = parameters.get<parameters::LATTICE_VTK_ITER_T>();
  const int statIter = parameters.get<parameters::LATTICE_STAT_ITER_T>();

  if ( iT==0 ) {
    // Writes the cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( lattice );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( lattice );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }
  // Get statistics
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();
    lattice.getStatistics().print( iT, iT );
  }

  //T delta_p = 0.;
  //T final_radius = 1.;
  // Writes the VTK files
  if ( iT%vtkIter==0 ) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticeDensity2D<T, DESCRIPTOR> density_lat( lattice );
    SuperIdentity2D<T,T> _density_lat( density_lat );
    _density_lat.getName() = "lattice density";

    AnalyticalConst2D<T,T> C_rho( lattice.getUnitConverter().getConversionFactorDensity() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> _C_rho( C_rho, lattice );
    SuperIdentity2D<T,T> _density_phs( density_lat * _C_rho );
    _density_phs.getName() = "physical density";

    SuperLatticeVelocity2D<T, DESCRIPTOR> velocityEq_lat( lattice );
    SuperIdentity2D<T,T> _velocityEq_lat( velocityEq_lat );
    _velocityEq_lat.getName() = "equibibrium distribution lattice velocity";

    SuperLatticeField2D<T, DESCRIPTOR, FORCE> force_lat( lattice );
    SuperIdentity2D<T,T> _force_lat( force_lat );
    _force_lat.getName() = "lattice force";

    AnalyticalConst2D<T,T> C_F( lattice.getUnitConverter().getConversionFactorForce() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> _C_F( C_F, lattice );
    SuperIdentity2D<T,T> _force_phs( force_lat * _C_F );
    _force_phs.getName() = "physical force";

    SuperLatticeField2D<T, DESCRIPTOR, CHEM_POTENTIAL> chemical_lat( lattice );
    SuperIdentity2D<T,T> _chemical_lat( chemical_lat );
    _chemical_lat.getName() = "lattice chemical potential";

    AnalyticalConst2D<T,T> C_MU( lattice.getUnitConverter().getConversionFactorChemicalPotential() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> _C_MU( C_MU, lattice );
    SuperIdentity2D<T,T> _chemical_phs( chemical_lat * _C_MU );
    _chemical_phs.getName() = "physical chemical potential";

    SuperLatticeField2D<T, DESCRIPTOR, OMEGA> omega( lattice );
    SuperIdentity2D<T,T> _omega( omega );
    _omega.getName() = "omega";

    AnalyticalConst2D<T,T> half( T(0.5) );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> _half( half, lattice );
    SuperIdentity2D<T,T> _velocity_lat( velocityEq_lat + _half * force_lat );
    _velocity_lat.getName() = "lattice velocity";

    AnalyticalConst2D<T,T> C_U( lattice.getUnitConverter().getConversionFactorVelocity() );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> _C_U( C_U, lattice );
    SuperIdentity2D<T,T> _velocity_phs( _velocity_lat * _C_U );
    _velocity_phs.getName() = "physical velocity";

    vtmWriter.addFunctor( _density_lat );
    vtmWriter.addFunctor( _density_phs );
    vtmWriter.addFunctor( _force_lat );
    vtmWriter.addFunctor( _force_phs );
    vtmWriter.addFunctor( _chemical_lat );
    vtmWriter.addFunctor( _chemical_phs );
    vtmWriter.addFunctor( _omega );
    vtmWriter.addFunctor( _velocity_lat );
    vtmWriter.addFunctor( _velocity_phs );
    vtmWriter.write( iT );

    // Compute error in Interface velocity
    // first: compute vapor phase velocity
    const T rho_V = parameters.get<parameters::RHO_VAPOR_ANALYTICAL>();
    const T rho_L = parameters.get<parameters::RHO_LIQUID_ANALYTICAL>();
    const T Ui_sol = parameters.get<parameters::INTERFACE_VELOCITY_ANALYTICAL>();
    AnalyticalFfromSuperF2D<T,T> interpolVelocity( _velocity_phs, true, 1 );
    T dx = lattice.getUnitConverter().getConversionFactorLength();
    T pos[2] = {dx, 0.};
    T Ui_num[2] = {0., 0.};
    interpolVelocity( &Ui_num[0], pos );
    // Compute interface velocity from vapor phase velocity
    Ui_num[0] = rho_V * Ui_num[0] / ( rho_V - rho_L );
    clout << "Numerical interface velocity [m/s]: " << Ui_num[0] << std::endl;
    clout << "Analytical interface velocity [m/s]: " << Ui_sol << std::endl;
    clout << "Relative error [%]: " << 100. * abs( Ui_num[0] - Ui_sol ) / abs( Ui_sol ) << std::endl;
  }
}

/// @brief Execute simulation: set initial values and run time loop
/// @param myCase The Case instance which keeps the simulation data
void simulate(MyCase& myCase)
{
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  //auto& converter = lattice.getUnitConverter();

  // Main Loop with Timer
  OstreamManager clout(std::cout, "simulate");
  int iT = 0;
  const size_t maxIter = parameters.get<parameters::LATTICE_MAX_T>();
  const T length = parameters.get<parameters::LIQUID_PHASE_LENGTH>();
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( maxIter, geometry.getStatistics().getNvoxel() );
  timer.start();
  std::vector<T> output = {1./length,0};
  for ( iT=0; iT<=maxIter; ++iT ) {

    // Apply boundary conditions
    boundaryCondition(myCase);

    // Compute relaxation frequency
    computeRelaxationFrequency(myCase);

    // Computation and output of the results
    getResults(myCase, iT, timer);
    if ( std::isnan(lattice.getStatistics().getAverageEnergy())) {
      clout << "unstable!" << std::endl;
      break;
    }

    // Collide and stream (and coupling) execution
    lattice.collideAndStream();

    // Compute chemical potential
    lattice.executePostProcessors(stage::PreCoupling());
    lattice.getCommunicator(stage::PreCoupling()).communicate();
    myCase.getOperator("Chemical").apply();

    // Compute thermodynamic force
    lattice.executePostProcessors(stage::PreCoupling());
    lattice.getCommunicator(stage::PreCoupling()).communicate();
    myCase.getOperator("Force").apply();

  }
  timer.stop();
  timer.printSummary();

  clout << "Finish" << std::endl;
}

/// Update boundary values at times (and external fields, if they exist)
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
/// @note Be careful: boundary values have to be set using lattice units
void setTemporalValues(MyCase& myCase, std::size_t iT)
{
  // Nothing to do here, because simulation does not depend on time
}

/// Set initial condition for primal variables
/// @param myCase The Case instance which keeps the simulation data
/// @note Be careful: initial values have to be set using lattice units
void setInitialValues(MyCase& myCase)
{
  // Nothing to do here
}

/// Setup and run a simulation
int main( int argc, char *argv[] )
{
  initialize( &argc, &argv );

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION>(400);
    myCaseParameters.set<LIQUID_PHASE_LENGTH>(300); // liquid phase length [m]
    myCaseParameters.set<parameters::NU>(0.1); // interface viscosity [m2.s-1]
    myCaseParameters.set<RHO_VAPOR>(1.); // vapor density [kg.m-3]
    myCaseParameters.set<RHO_LIQUID>(2.); // liquid density [kg.m-3]
    myCaseParameters.set<SURFACE_TENSION>(1.e-3); // surface tension [N.m-1]
    myCaseParameters.set<DOMAIN_EXTENT>({400, 10}); // domain length and height [m]

   // Boundary conditions in physical units
   /**
    * Boundary conditions for: p_sat/p_c = 0.99,
    * we use the analytical solution to initialize
    * the system close to the final solution
    */
    myCaseParameters.set<RHO_VAPOR_ANALYTICAL>(0.9809107952180959);
    myCaseParameters.set<RHO_LIQUID_ANALYTICAL>(1.9926932521825274);
    myCaseParameters.set<BOUNDARY_CHEMICAL_POTENTIAL>(-3.4274274827508194e-05);
    myCaseParameters.set<INTERFACE_VELOCITY_ANALYTICAL>(0.0021753492521563193);

    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.8); // Relaxation time in lattice units
    myCaseParameters.set<THICKNESS>(5.); // Interface thickness in lattice units; physThickness = delta_x * thickness

    // Time parameters in lattice units
    myCaseParameters.set<LATTICE_MAX_T>(30000);
    myCaseParameters.set<LATTICE_VTK_ITER_T>(100);
    myCaseParameters.set<LATTICE_STAT_ITER_T>(100);
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
