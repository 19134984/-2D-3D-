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

/* porousPoiseuille3d.cpp:
 * This example examines a 3D Poseuille flow with porous media.
 * Two porous media LB methods can be used here:
 * Spaid and Phelan (doi:10.1063/1.869392), or
 * Guo and Zhao (doi:10.1103/PhysRevE.66.036304)
 *
 * Case-specific arguments:
 * BOUNDARY_TYPE: 0=bounceBack, 1=local, 2=interpolated
 * POROSITY_TYPE: 0=BGK, 1=SpaidPhelan, 2=GuoZhao
 * PERMEABILIY: default 1e-2
 * Default: Resolution=50, Permeability=1e-2, BoundaryTpe=interpolated, PorosityType=SpaidPhelan
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;

enum class PorosityType: int {
  BGK             = 0,
  SPAID_PHELAN    = 1,
  GUO_ZHAO        = 2
};

namespace olb::parameters {

struct POROSITY_TYPE  : public descriptors::TYPED_FIELD_BASE<PorosityType,1> { };
struct PERMEABILITY   : public descriptors::FIELD_BASE<1> { };
struct INITIAL_PRESSURE_L : public descriptors::FIELD_BASE<1> { };
struct CONVERGENCE_CHECK_T        : public descriptors::FIELD_BASE<1> { };

}

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q19<>>
>;

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const T length = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T radius = parameters.get<parameters::DOMAIN_EXTENT>()[1] / 2.;
  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();

  Vector<T, 3> center0(0, radius, radius);
  Vector<T, 3> center1(length + 0.5 * physDeltaX, radius, radius);
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  IndicatorLayer3D<T> extendedDomain(pipe, physDeltaX);

  Mesh<T,MyCase::d> mesh(extendedDomain, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

template <typename T, typename DESCRIPTOR>
class PhysicalToLatticeVelocityF3D: public AnalyticalF3D<T,T> {
protected:
  AnalyticalF3D<T,T>* f;
  UnitConverter<T,DESCRIPTOR> converter;
public:
  PhysicalToLatticeVelocityF3D(AnalyticalF3D<T,T>* f_, UnitConverter<T,DESCRIPTOR> const& converter_)
    : AnalyticalF3D<T,T>(3), f(f_), converter(converter_) {};
  bool operator()(T output[], const T x[]) override
  {
    (*f)(output, x);
    for (int i=0; i<3; ++i) {
      output[i] = converter.getLatticeVelocity( output[i] );
    }
    return true;
  };
};

// Approximation of the modified Bessel function (doi:10.1088/1742-6596/1043/1/012003)
template <typename T>
T besselApprox( T x )
{
  return util::cosh(x) / util::pow( 1 + 0.25*util::pow(x,2), 0.25 ) * ( 1 + 0.24273*util::pow(x,2) )/( 1 + 0.43023*util::pow(x,2) );
}

/// Functional to calculate velocity profile on pipe with porous media.
template <typename T>
class PorousPoiseuille3D : public AnalyticalF3D<T,T> {
protected:
  T Kin, mu, dp, epsilon, radius;
  bool trunc;
public:
  PorousPoiseuille3D( MyCase& myCase, T radius_ )
    : AnalyticalF3D<T,T>(3), radius(radius_)
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
    T dist = util::sqrt( util::pow(x[1]-radius, 2.) + util::pow(x[2]-radius, 2.) );
    output[0] = Kin / mu * dp * ( 1. - besselApprox(r*dist) / besselApprox(r*radius)  );
    output[1] = 0.;
    output[2] = 0.;
    if ( dist > radius ) { output[0] = 0.; }
    return true;
  };
};

void prepareGeometry( MyCase& myCase ) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry    = myCase.getGeometry();
  auto& parameters  = myCase.getParameters();
  T physDeltaX      = parameters.get<parameters::PHYS_DELTA_X>();
  T length          = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  T radius          = parameters.get<parameters::DOMAIN_EXTENT>()[1] / 2.;

  geometry.rename(0, 2);

  Vector<T,3> center0(-physDeltaX * 0.2, radius, radius);
  Vector<T,3> center1(length, radius, radius);
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  geometry.rename(2, 1, pipe);

  Vector<T, 3> origin(0, radius, radius);
  Vector<T, 3> extend = origin;

  // Set material number for inflow
  origin[0] = -physDeltaX * 2;
  extend[0] = physDeltaX * 2;
  IndicatorCylinder3D<T> inflow(origin, extend, radius);
  geometry.rename(2, 3, 1, inflow);

  // Set material number for outflow
  origin[0] = length - 2 * physDeltaX;
  extend[0] = length + 2 * physDeltaX;
  IndicatorCylinder3D<T> outflow(extend, origin, radius);
  geometry.rename(2, 4, 1, outflow);

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  // Removes all not needed boundary voxels inside the surface
  geometry.innerClean();

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
  T       physU     = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  T       physRho   = parameters.get<parameters::PHYS_CHAR_DENSITY>();
  T       tau       = parameters.get<parameters::LATTICE_RELAXATION_TIME>();
  T       length    = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  T       diameter  = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  T       radius    = diameter / 2.;

  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR>>(
    (size_t) N,                 // resolution: number of voxels per charPhysL
    (T)     tau,                // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)     diameter,           // charPhysLength: reference length of simulation geometry
    (T)     physU,              // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)     diameter*physU/Re,  // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)     physRho             // physDensity: physical density in __kg / m^3__
  );
  auto& converter = lattice.getUnitConverter();
  // Prints the converter log as console output
  //converter.print();
  // Writes the converter log in a file
  converter.write("porousPoiseuille3d");

  // Material=1,3,4 -->bulk dynamics
  switch ( parameters.get<parameters::POROSITY_TYPE>() ) {
  case PorosityType::BGK:
    dynamics::set<BGKdynamics>(lattice, geometry.getMaterialIndicator({1,3,4}));
    break;
  case PorosityType::SPAID_PHELAN:
    dynamics::set<PorousBGKdynamics>(lattice, geometry.getMaterialIndicator({1,3,4}));
    break;
  case PorosityType::GUO_ZHAO:
  default:
    dynamics::set<GuoZhaoBGKdynamics>(lattice, geometry.getMaterialIndicator({1,3,4}));
    break;
  }

  PorosityType porosityType = parameters.get<parameters::POROSITY_TYPE>();
  T Kin = parameters.get<parameters::PERMEABILITY>();
  T h = converter.getPhysDeltaX();
  switch ( porosityType ) {
    case PorosityType::BGK:
      break;
    case PorosityType::SPAID_PHELAN: {
      T tau = converter.getLatticeRelaxationTime();
      T nu = converter.getLatticeViscosity();
      T d = 1. - ( h*h*nu*tau/Kin );
      clout << "Lattice Porosity: " << d << std::endl;
      clout << "Kmin: " << h*h*nu*tau << std::endl;
      if (Kin < h*h*nu*tau) {
        clout << "WARNING: Chosen K is too small!" << std::endl;
        exit(1);
      }
      AnalyticalConst3D<T,T> porosity( d );
      for (int i: { 0,1,2,3,4 }) {
        fields::set<descriptors::POROSITY>(lattice, geometry.getMaterialIndicator({i}), porosity);
      }
    }
      break;
    case PorosityType::GUO_ZHAO:
    default:
      AnalyticalConst3D<T,T> eps( parameters.get<parameters::EPSILON>() );
      AnalyticalConst3D<T,T> Nu( converter.getLatticeViscosity() );
      AnalyticalConst3D<T,T> k( Kin / ( h*h ) );
      for (int i: {0,1,2,3,4}) {
        fields::set<descriptors::EPSILON>(lattice, geometry.getMaterialIndicator({i}), eps);
        fields::set<descriptors::NU>(lattice, geometry.getMaterialIndicator({i}), Nu);
        fields::set<descriptors::K>(lattice, geometry.getMaterialIndicator({i}), k);
      }
      break;
  }

  // Bouzidi
  Vector<T, 3> center0(0, radius, radius);
  Vector<T, 3> center1(length, radius, radius);
  center0[0] -= 0.5*converter.getPhysDeltaX();
  center1[0] += 0.5*converter.getPhysDeltaX();
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  setBouzidiBoundary(lattice, geometry, 2, pipe);

  // Interp
  boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 2);
  boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 3);
  boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 4);

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
  T     length        = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  T     radius        = parameters.get<parameters::DOMAIN_EXTENT>()[1] / 2.;

  // Initial conditions
  // Pressure for Poiseuille flow with maximum velocity of charU at K->infty
  T p0 = 4. * converter.getPhysViscosity() * converter.getCharPhysVelocity() * length / (radius * radius);
  T p0L = converter.getLatticePressure(p0);
  AnalyticalLinear3D<T, T> rho( converter.getPhysDensity(-p0L / length * descriptors::invCs2<T,DESCRIPTOR>()), 0, 0,
                                converter.getPhysDensity(p0L * descriptors::invCs2<T,DESCRIPTOR>() + 1));

  T dp = p0/length;
  T mu = converter.getPhysViscosity()*converter.getPhysDensity();

  parameters.set<parameters::INITIAL_PRESSURE_L >(  p0L );
  parameters.set<parameters::PRESSURE_DROP  >(   dp );
  parameters.set<parameters::PHYS_KINEMATIC_VISCOSITY          >(   mu );

  //CirclePoiseuille3D<T> uSol( {0., radius, radius}, {1, 0, 0}, converter.getCharPhysVelocity(), radius );
  PorousPoiseuille3D<T> uSol( myCase, radius );

  // Initialize all values of distribution functions to their local equilibrium
  for (int i: { 0,1,2,3,4 }) {
    momenta::setVelocity(lattice, geometry.getMaterialIndicator({i}), uSol);
    momenta::setDensity(lattice, geometry.getMaterialIndicator({i}), rho);
  }

  lattice.initialize();
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{ }

void convergenceCheck(MyCase& myCase) {
            // AnalyticalF3D<T,T>& uSol)

  OstreamManager clout( std::cout,"error" );
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& parameters  = myCase.getParameters();
  auto& geometry    = myCase.getGeometry();
  auto& lattice     = myCase.getLattice(NavierStokes{});
  auto& converter   = lattice.getUnitConverter();
  const T length    = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T radius    = parameters.get<parameters::DOMAIN_EXTENT>()[1] / 2.;
  const T p0L       = parameters.get<parameters::INITIAL_PRESSURE_L>();

  int tmp[]= { };
  T result[2]= { };

  PorousPoiseuille3D<T> uSol( myCase, radius );
  auto indicatorF = geometry.getMaterialIndicator(1);
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> u( lattice, converter );

  // Velocity error
  SuperAbsoluteErrorL1Norm3D<T> absVelocityErrorNormL1(u, uSol, indicatorF);
  absVelocityErrorNormL1(result, tmp);
  clout << "velocity-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm3D<T> relVelocityErrorNormL1(u, uSol, indicatorF);
  relVelocityErrorNormL1(result, tmp);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm3D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, tmp);
  clout << "velocity-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm3D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, tmp);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm3D<T> absVelocityErrorNormLinf(u, uSol, indicatorF);
  absVelocityErrorNormLinf(result, tmp);
  clout << "velocity-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm3D<T> relVelocityErrorNormLinf(u, uSol, indicatorF);
  relVelocityErrorNormLinf(result, tmp);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;

  // Pressure error
  AnalyticalLinear3D<T,T> pressureSol( -converter.getPhysPressure( p0L )/length, 0, 0, converter.getPhysPressure( p0L ) );
  SuperLatticePhysPressure3D<T,DESCRIPTOR> pressure( lattice,converter );

  SuperAbsoluteErrorL1Norm3D<T> absPressureErrorNormL1(pressure, pressureSol, indicatorF);
  absPressureErrorNormL1(result, tmp);
  clout << "pressure-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm3D<T> relPressureErrorNormL1(pressure, pressureSol, indicatorF);
  relPressureErrorNormL1(result, tmp);
  clout << "; pressure-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm3D<T> absPressureErrorNormL2(pressure, pressureSol, indicatorF);
  absPressureErrorNormL2(result, tmp);
  clout << "pressure-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm3D<T> relPressureErrorNormL2(pressure, pressureSol, indicatorF);
  relPressureErrorNormL2(result, tmp);
  clout << "; pressure-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm3D<T> absPressureErrorNormLinf(pressure, pressureSol, indicatorF);
  absPressureErrorNormLinf(result, tmp);
  clout << "pressure-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm3D<T> relPressureErrorNormLinf(pressure, pressureSol, indicatorF);
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
  const T             maxPhysT      = parameters.get<parameters::PHYS_MAX_T>();
  const T             length        = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T             diameter      = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const T             radius        = diameter / 2.;
  bool                hasConverged  = parameters.get<parameters::CONVERGED>();

  SuperVTMwriter3D<T> vtmWriter( "porousPoiseuille3d" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( lattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( lattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  //CirclePoiseuille3D<T> uSol( {0., radius, radius}, {1, 0, 0}, converter.getCharPhysVelocity(), radius );
  PorousPoiseuille3D<T> uSol( myCase, radius );
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> uSolF( uSol, lattice );
  vtmWriter.addFunctor( uSolF );

  const int vtmIter  = converter.getLatticeTime( maxPhysT/20. );
  const int statIter = converter.getLatticeTime( maxPhysT/20. );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( lattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( lattice );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }

  // Writes the vtm files and profile text file
  if ( iT%vtmIter==0 || hasConverged ) {
    vtmWriter.write( iT );
  }


  // Writes output on the console
  if ( iT%statIter==0 || hasConverged ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    lattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Calculate inflow and outflow flux
    std::vector<int> materials = { 1, 3, 4 };
    Vector<T,3> normal( 1, 0, 0 );
    auto mode = BlockDataReductionMode::Discrete;
    Vector<T,3> posInflow = geometry.getStatistics().getMinPhysR( 1 );
    Vector<T,3> posOutflow = geometry.getStatistics().getMaxPhysR( 1 );

    SuperPlaneIntegralFluxVelocity3D<T> vFluxIn( lattice, converter,
        geometry, posInflow, normal, materials, mode );
    SuperPlaneIntegralFluxPressure3D<T> pFluxIn( lattice, converter,
        geometry, posInflow, normal, materials, mode );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOut( lattice, converter,
        geometry, posOutflow, normal, materials, mode );
    SuperPlaneIntegralFluxPressure3D<T> pFluxOut( lattice, converter,
        geometry, posOutflow, normal, materials, mode );

    vFluxIn.print( "Inflow" );
    pFluxIn.print( "Inflow" );
    vFluxOut.print( "Outflow" );
    pFluxOut.print( "Outflow" );

    // Error norms
    AnalyticalFfromSuperF3D<T> intpolatePressure( pressure, true );

    T point1[3] = {0, radius, radius};
    T point2[3] = {0, radius, radius};

    point1[0] = length*0.5 - length*0.01;
    point2[0] = length*0.5 + length*0.01;

    T p1, p2;
    intpolatePressure( &p1,point1 );
    intpolatePressure( &p2,point2 );

    clout << "pressure1=" << p1;
    clout << "; pressure2=" << p2;

    T pressureDrop = p1-p2;
    clout << "; pressureDrop=" << pressureDrop << std::endl;

    convergenceCheck( myCase );

    // Gnuplot
    Gnuplot<T> gplot( "velocityProfile" );
    T uAnalytical[3] = {};
    T uNumerical[3] = {};
    AnalyticalFfromSuperF3D<T> intpolateVelocity( velocity, true );
    for (int i=0; i<101; i++) {
      T yInput = diameter*T(i)/T(100);
      T input[3] = {length*T(0.5), yInput, radius};
      uSol(uAnalytical, input);
      intpolateVelocity(uNumerical, input);
      gplot.setData( yInput, {uAnalytical[0], uNumerical[0]}, {"analytical","numerical"} );
    }

    // Create PNG file
    gplot.writePNG();
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
  const T residuum      = parameters.get<parameters::CONVERGENCE_PRECISION>();
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
  initialize( &argc, &argv );

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<CONVERGENCE_PRECISION  >(     1e-5 );
    myCaseParameters.set<DOMAIN_EXTENT          >( {2., 1., 1.} );  // length x diameter^2
    myCaseParameters.set<RESOLUTION             >(       21 );
    myCaseParameters.set<REYNOLDS               >(       1. );
    myCaseParameters.set<PERMEABILITY           >(     1e-2 );
    myCaseParameters.set<POROSITY_TYPE          >( PorosityType::SPAID_PHELAN );
    myCaseParameters.set<PHYS_MAX_T             >(       20 );
    myCaseParameters.set<CONVERGENCE_CHECK_T    >( 0.0125*myCaseParameters.get<PHYS_MAX_T>() );
    myCaseParameters.set<PHYS_CHAR_VELOCITY     >(       1. );
    myCaseParameters.set<PHYS_CHAR_DENSITY      >(       1. );
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(       .8 );
    myCaseParameters.set<CONVERGED              >(    false );
    myCaseParameters.set<EPSILON                >(       1. );  // Porosity (Spaid and Phelan can only handle epsilon=1)
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
