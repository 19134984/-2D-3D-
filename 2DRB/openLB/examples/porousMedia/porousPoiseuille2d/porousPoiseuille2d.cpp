/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2019 Fabian Klemens, Davide Dapelo, Mathias J. Krause
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

/* porousPoiseuille2d.cpp:
 * This example examines a 2D Poseuille flow with porous media.
 * Two porous media LB methods can be used here:
 * Spaid and Phelan (doi:10.1063/1.869392), or
 * Guo and Zhao (doi:10.1103/PhysRevE.66.036304)
 *
 * Case-specific arguments:
 * BOUNDARY_TYPE: 0=bounceBack, 1=local, 2=interpolated
 * POROSITY_TYPE: 0=BGK, 1=SpaidPhelan, 2=GuoZhao, 3=GuoZhaoSmagorinsky
 * PERMEABILIY: default 1e-2
 * Default: Resolution=50, Permeability=1e-2, BoundaryTpe=interpolated, PorosityType=GuoZhao
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;

enum class BoundaryType: int {
  bounceBack    = 0,
  local         = 1,
  interpolated  = 2
};

enum class PorosityType: int {
  BGK             = 0,
  SPAID_PHELAN    = 1,
  GUO_ZHAO        = 2,
  GUO_ZHAO_SMAGO  = 3
};

namespace olb::parameters {

struct BOUNDARY_TYPE  : public descriptors::TYPED_FIELD_BASE<BoundaryType,1> { };
struct POROSITY_TYPE  : public descriptors::TYPED_FIELD_BASE<PorosityType,1> { };
struct PERMEABILITY   : public descriptors::FIELD_BASE<1> { };
struct INITIAL_PRESSURE_L : public descriptors::FIELD_BASE<1> { };
struct CONVERGENCE_CHECK_T        : public descriptors::FIELD_BASE<1> { };

}

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D2Q9<>>
>;

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector origin{0, 0};
  IndicatorCuboid2D<T> cuboid(extent, origin);

  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

// Functor to convert physical to lattice velocity
template <typename T, typename DESCRIPTOR>
class PhysicalToLatticeVelocityF2D: public AnalyticalF2D<T,T> {
protected:
  AnalyticalF2D<T,T>* f;
  UnitConverter<T,DESCRIPTOR> converter;
public:
  PhysicalToLatticeVelocityF2D(AnalyticalF2D<T,T>* f_, UnitConverter<T,DESCRIPTOR> const& converter_)
    : AnalyticalF2D<T,T>(2), f(f_), converter(converter_) {};
  bool operator()(T output[], const T x[]) override
  {
    (*f)(output, x);
    for (int i=0; i<2; ++i) {
      output[i] = converter.getLatticeVelocity( output[i] );
    }
    return true;
  };
};

/// Velocity profile in a pipe filled with isotropic porous media
template <typename T>
class PorousPoiseuille2D : public AnalyticalF2D<T,T> {
protected:
  T Kin, dp, mu, epsilon, radius, wallOffset;
public:
  PorousPoiseuille2D( MyCase& myCase, T radius_, T wallOffset_=0.)
    : AnalyticalF2D<T,T>(2), radius(radius_), wallOffset(wallOffset_)
  {
    auto& parameters = myCase.getParameters();
    Kin       = parameters.get<parameters::PERMEABILITY>();
    dp        = parameters.get<parameters::PRESSURE_DROP>();
    mu        = parameters.get<parameters::PHYS_KINEMATIC_VISCOSITY>();
    epsilon   = parameters.get<parameters::EPSILON>();
  }

  bool operator()(T output[], const T x[]) override
  {
    T r = util::sqrt(epsilon/Kin);
    output[0] = dp / mu * Kin
     * (1. - util::cosh(r*(x[1] - radius)) / util::cosh(r*(radius-wallOffset)));
    output[1] = 0.;
    if ( x[1] < wallOffset || x[1] > 2.*radius - wallOffset ) {
      output[0] = 0.;
    }
    return true;
  }
};

void prepareGeometry( MyCase& myCase ) {
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  T physDeltaX2 = parameters.get<parameters::PHYS_DELTA_X>() / 2.;

  geometry.rename( 0, 2 );
  geometry.rename( 2, 1, {1, 1} );

  Vector<T,2> extentDomain = parameters.get<parameters::DOMAIN_EXTENT>();
  Vector<T,2> extent = extentDomain;
  extent[0] = 2 * physDeltaX2;
  Vector<T,2> origin(-physDeltaX2);

  // Set material number for inflow
  IndicatorCuboid2D<T> inflow( extent, origin );
  geometry.rename( 2, 3, 1, inflow );

  // Set material number for outflow
  origin[0] = extentDomain[0] - physDeltaX2;
  IndicatorCuboid2D<T> outflow( extent, origin );
  geometry.rename( 2, 4, 1, outflow );

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  // Removes all not needed boundary voxels inside the surface
  geometry.innerClean();
  geometry.checkForErrors();

  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( MyCase& myCase ) {
  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});
  size_t  N   = parameters.get<parameters::RESOLUTION>();
  T       Re  = parameters.get<parameters::REYNOLDS>();

  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR>>(
    size_t {N},  // resolution: number of voxels per charPhysL
    (T)   0.8,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   1,     // charPhysLength: reference length of simulation geometry
    (T)   1,     // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1./Re, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.     // physDensity: physical density in __kg / m^3__
  );
  auto& converter = lattice.getUnitConverter();
  converter.print();

  // Material=1,3,4 -->bulk dynamics
  switch ( parameters.get<parameters::POROSITY_TYPE>() ) {
  case PorosityType::BGK:
    dynamics::set<BGKdynamics>(lattice, geometry.getMaterialIndicator({1,3,4}));
    break;
  case PorosityType::SPAID_PHELAN:
    dynamics::set<PorousBGKdynamics>(lattice, geometry.getMaterialIndicator({1,3,4}));
    break;
  case PorosityType::GUO_ZHAO:
    dynamics::set<GuoZhaoBGKdynamics>(lattice, geometry.getMaterialIndicator({1,3,4}));
    break;
  case PorosityType::GUO_ZHAO_SMAGO:
  default:
    dynamics::set<SmagorinskyGuoZhaoBGKdynamics>(lattice, geometry.getMaterialIndicator({1,3,4}));
    break;
  }

  BoundaryType boundaryType = parameters.get<parameters::BOUNDARY_TYPE>();
  switch ( boundaryType ) {
    case BoundaryType::bounceBack:
      boundary::set<boundary::BounceBack>(lattice, geometry, 2);
      boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 3);
      boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 4);
      break;
    case BoundaryType::local:
      boundary::set<boundary::LocalVelocity>(lattice, geometry, 2);
      boundary::set<boundary::LocalVelocity>(lattice, geometry, 3);
      boundary::set<boundary::LocalPressure>(lattice, geometry, 4);
      break;
    case BoundaryType::interpolated:
    default:
      boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 2);
      boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 3);
      boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 4);
      break;
  }

  PorosityType porosityType = parameters.get<parameters::POROSITY_TYPE>();
  T Kin = parameters.get<parameters::PERMEABILITY>();
  switch ( porosityType ) {
    case PorosityType::BGK:
      break;
    case PorosityType::SPAID_PHELAN: {
      T tau = converter.getLatticeRelaxationTime();
      T nu = converter.getLatticeViscosity();
      T h = converter.getPhysDeltaX();
      T d = 1. - (h*h*nu*tau/Kin);
      clout << "Lattice Porosity: " << d << std::endl;
      clout << "Kmin: " << h*h*nu*tau << std::endl;
      if (Kin < h*h*nu*tau) {
        clout << "WARNING: Chosen K is too small!" << std::endl;
        exit(1);
      }
      AnalyticalConst2D<T,T> porosity( d );
      for (int i: {
            0,1,2,3,4
          }) {
        fields::set<descriptors::POROSITY>(lattice, geometry.getMaterialIndicator(i), porosity);
      }
    }
      break;
    case PorosityType::GUO_ZHAO_SMAGO:
      lattice.setParameter<collision::LES::SMAGORINSKY>(T(0.14));
    case PorosityType::GUO_ZHAO:
    default:
      AnalyticalConst2D<T,T> eps( parameters.get<parameters::EPSILON>() );
      AnalyticalConst2D<T,T> Nu( converter.getLatticeViscosity() );
      AnalyticalConst2D<T,T> k( Kin/util::pow(converter.getPhysDeltaX(), 2.) );
      for (int i: {0,1,2,3,4}) {
        fields::set<descriptors::K>(lattice, geometry.getMaterialIndicator(i), k);
        fields::set<descriptors::EPSILON>(lattice, geometry.getMaterialIndicator(i), eps);
        fields::set<descriptors::NU>(lattice, geometry.getMaterialIndicator(i), Nu);
      }
      break;
  }

  lattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  clout << "Prepare Lattice ... OK" << std::endl;
}

/// Set initial condition for primal variables (velocity and density)
void setInitialValues( MyCase& myCase ) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& lattice   = myCase.getLattice(NavierStokes{});
  auto& geometry  = myCase.getGeometry();
  auto& converter = lattice.getUnitConverter();
  auto& parameters= myCase.getParameters();
  T     lx        = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  T     ly        = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  BoundaryType boundaryType = parameters.get<parameters::BOUNDARY_TYPE>();

  // Initial conditions
  // Pressure for Poiseuille flow with maximum velocity of charU at K->infty
  T p0 = 8.*converter.getPhysViscosity()*converter.getCharPhysVelocity()*lx/( ly*ly );
  // Pressure for PorousPoiseuille with maximum velocity of charU for every permeability K
  //p0 = converter.getCharPhysVelocity() * converter.getPhysViscosity() * converter.getPhysDensity() / Kin * lx / (1. - 1./util::cosh(util::sqrt(1./Kin)*ly/2.));

  T p0L = converter.getLatticePressure(p0);

  AnalyticalLinear2D<T,T> rho( converter.getPhysDensity( - p0L/lx*descriptors::invCs2<T,DESCRIPTOR>()),
                                0,
                               converter.getPhysDensity( p0L*descriptors::invCs2<T,DESCRIPTOR>()+1 ));

  T dp = p0/lx;
  T mu = converter.getPhysViscosity()*converter.getPhysDensity();

  parameters.set<parameters::INITIAL_PRESSURE_L >(  p0L );
  parameters.set<parameters::PRESSURE_DROP  >(   dp );
  parameters.set<parameters::PHYS_KINEMATIC_VISCOSITY          >(   mu );

  const T wallOffset = ( boundaryType == BoundaryType::bounceBack ) ? 0.5 * converter.getPhysDeltaX() : T(0.);
  PorousPoiseuille2D<T> uSol( myCase, ly/2., wallOffset );
  PhysicalToLatticeVelocityF2D<T,DESCRIPTOR> u( &uSol, converter );

  // Initialize all values of distribution functions to their local equilibrium
  for ( int i: { 0,1,2,3,4 } ) {
    momenta::setVelocity(lattice, geometry.getMaterialIndicator({i}), uSol);
    momenta::setDensity(lattice, geometry.getMaterialIndicator({i}), rho);
    //lattice.iniEquilibrium( geometry, i, rho, u ); // gives problems with non-standard equilibria
  }

  lattice.initialize();
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{ }

void convergenceCheck(MyCase& myCase) {
  OstreamManager clout( std::cout,"error" );
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& parameters  = myCase.getParameters();
  auto& geometry    = myCase.getGeometry();
  auto& lattice     = myCase.getLattice(NavierStokes{});
  auto& converter   = lattice.getUnitConverter();
  const BoundaryType  boundaryType  = parameters.get<parameters::BOUNDARY_TYPE>();
  const T             lx            = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T             ly            = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const T             p0L           = parameters.get<parameters::INITIAL_PRESSURE_L>();

  int tmp[]= { };
  T result[2]= { };

  const T wallOffset = (boundaryType == BoundaryType::bounceBack) ? 0.5 * converter.getPhysDeltaX() : T(0.);
  PorousPoiseuille2D<T> uSol( myCase, ly/2., wallOffset );
  SuperLatticePhysVelocity2D<T,DESCRIPTOR> u( lattice,converter );
  auto indicatorF = geometry.getMaterialIndicator(1);

  // velocity error
  SuperAbsoluteErrorL1Norm2D<T> absVelocityErrorNormL1(u, uSol, indicatorF);
  absVelocityErrorNormL1(result, tmp);
  clout << "velocity-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm2D<T> relVelocityErrorNormL1(u, uSol, indicatorF);
  relVelocityErrorNormL1(result, tmp);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm2D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, tmp);
  clout << "velocity-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm2D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, tmp);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm2D<T> absVelocityErrorNormLinf(u, uSol, indicatorF);
  absVelocityErrorNormLinf(result, tmp);
  clout << "velocity-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm2D<T> relVelocityErrorNormLinf(u, uSol, indicatorF);
  relVelocityErrorNormLinf(result, tmp);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;

  // pressure error
  AnalyticalLinear2D<T,T> pressureSol( -converter.getPhysPressure( p0L )/lx, 0, converter.getPhysPressure( p0L ) );
  SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure( lattice,converter );

  SuperAbsoluteErrorL1Norm2D<T> absPressureErrorNormL1(pressure, pressureSol, indicatorF);
  absPressureErrorNormL1(result, tmp);
  clout << "pressure-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm2D<T> relPressureErrorNormL1(pressure, pressureSol, indicatorF);
  relPressureErrorNormL1(result, tmp);
  clout << "; pressure-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm2D<T> absPressureErrorNormL2(pressure, pressureSol, indicatorF);
  absPressureErrorNormL2(result, tmp);
  clout << "pressure-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm2D<T> relPressureErrorNormL2(pressure, pressureSol, indicatorF);
  relPressureErrorNormL2(result, tmp);
  clout << "; pressure-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm2D<T> absPressureErrorNormLinf(pressure, pressureSol, indicatorF);
  absPressureErrorNormLinf(result, tmp);
  clout << "pressure-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm2D<T> relPressureErrorNormLinf(pressure, pressureSol, indicatorF);
  relPressureErrorNormLinf(result, tmp);
  clout << "; pressure-Linf-error(rel)=" << result[0] << std::endl;
}

void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  OstreamManager clout( std::cout,"getResults" );
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& lattice   = myCase.getLattice(NavierStokes{});
  auto& geometry  = myCase.getGeometry();
  auto& converter = lattice.getUnitConverter();
  auto& parameters= myCase.getParameters();
  const BoundaryType  boundaryType  = parameters.get<parameters::BOUNDARY_TYPE>();
  const T             maxPhysT      = parameters.get<parameters::PHYS_MAX_T>();
  const T             lx            = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T             ly            = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  bool                hasConverged  = parameters.get<parameters::CONVERGED>();

  SuperVTMwriter2D<T> vtmWriter( "porousPoiseuille2d" );
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity( lattice, converter );
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure( lattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  const T wallOffset = ( boundaryType == BoundaryType::bounceBack ) ? 0.5 * converter.getPhysDeltaX() : T(0.);
  PorousPoiseuille2D<T> uSol( myCase, ly/2., wallOffset );
  SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR> uSolF( uSol, lattice);
  vtmWriter.addFunctor( uSolF );

  const int vtmIter  = converter.getLatticeTime( maxPhysT/20. );
  const int statIter = converter.getLatticeTime( maxPhysT/20. );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( lattice );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( lattice );
    geometry.rename( 0,2 );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }

  // Writes the vtm files and profile text file
  if ( iT%vtmIter==0 || hasConverged ) {
    vtmWriter.write( iT );

    SuperEuklidNorm2D<T, DESCRIPTOR> normVel( velocity );
    BlockReduction2D2D<T> planeReduction( normVel, 600, BlockDataSyncMode::ReduceOnly );
    // write output as JPEG
    heatmap::write(planeReduction, iT);
  }

  if ( hasConverged ) {
    Gnuplot<T> gplot( "centerVelocity" );
    T Ly = converter.getLatticeLength( ly );
    for ( int iY=0; iY<=Ly; ++iY ) {
      T dx = 1. / T(converter.getResolution());
      T point[2]= {T(),T()};
      point[0] = lx/2.;
      point[1] = ( T )iY/Ly;

      const T wallOffset = ( boundaryType == BoundaryType::bounceBack ) ? 0.5 * converter.getPhysDeltaX() : T(0.);
      PorousPoiseuille2D<T> uSol( myCase, ly/2., wallOffset );
      T analytical[2] = {T(),T()};
      uSol( analytical,point );
      SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity( lattice, converter );
      AnalyticalFfromSuperF2D<T> intpolateVelocity( velocity, true );
      T numerical[2] = {T(),T()};
      intpolateVelocity( numerical,point );
      gplot.setData( iY*dx, {analytical[0],numerical[0]}, {"analytical","numerical"} );
    }
    // Create PNG file
    gplot.writePNG();
  }

  // Writes output on the console
  if ( iT%statIter==0 || hasConverged ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    lattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Error norms
    convergenceCheck( myCase );
  }
}


void simulate(MyCase& myCase) {
  OstreamManager clout( std::cout, "simulation" );
  using T = MyCase::value_t;
  auto& parameters      = myCase.getParameters();
  auto& lattice         = myCase.getLattice(NavierStokes{});
  auto& converter       = lattice.getUnitConverter();
  const T physMaxT      = parameters.get<parameters::PHYS_MAX_T>();
  const T physInterval  = parameters.get<parameters::CONVERGENCE_CHECK_T>();
  const T residuum      = parameters.get<parameters::RESIDUUM>();
  const size_t iTmax    = lattice.getUnitConverter().getLatticeTime(physMaxT);

  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( iTmax, myCase.getGeometry().getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( converter.getLatticeTime( physInterval ), residuum );
  timer.start();

  for ( std::size_t iT = 0; iT < iTmax; ++iT ) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << std::endl;
      parameters.set<parameters::CONVERGED>( true );
      getResults( myCase, timer, iT );
      break;
    }
    lattice.collideAndStream();

    getResults( myCase, timer, iT );
    converge.takeValue( lattice.getStatistics().getAverageEnergy(), true );
  }

  timer.stop();
  timer.printSummary();
}

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<CONVERGENCE_CHECK_T>(      .25 );
    myCaseParameters.set<RESIDUUM           >(     1e-5 );
    myCaseParameters.set<DOMAIN_EXTENT      >( {2., 1.} );
    myCaseParameters.set<RESOLUTION         >(       50 );
    myCaseParameters.set<REYNOLDS           >(       1. );
    myCaseParameters.set<PERMEABILITY       >(     1e-2 );
    myCaseParameters.set<BOUNDARY_TYPE      >( BoundaryType::interpolated );
    myCaseParameters.set<POROSITY_TYPE      >( PorosityType::GUO_ZHAO );
    myCaseParameters.set<PHYS_MAX_T         >(       20 );
    myCaseParameters.set<CONVERGED          >(    false );
    myCaseParameters.set<EPSILON            >(       1. );  // Porosity (Spaid and Phelan can only handle epsilon=1)
    OLB_ASSERT( myCaseParameters.get<RESOLUTION>()    >= 1, "Fluid domain is too small" );
    OLB_ASSERT( myCaseParameters.get<PERMEABILITY>()  >= 0, "Permeability must be non-negative" );
  }
  if ( myCaseParameters.get<parameters::PHYS_DELTA_X>() == 0 ) {
    myCaseParameters.set<parameters::PHYS_DELTA_X>( 1. / myCaseParameters.get<parameters::RESOLUTION>() );
  } else {
    myCaseParameters.set<parameters::RESOLUTION>( 1. / myCaseParameters.get<parameters::PHYS_DELTA_X>() );
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
