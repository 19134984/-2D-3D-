/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Shota Ito
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

/* 3D force-driven Poiseuille flow with different non-Newtonian viscosity models
 * Implemented options are Newtonian, PowerLaw, Casson, and Carreau-Yasuda models.
 * The first three models are validated for the same pressure drop and dynamic
 * viscosity used for the unit-conversion. The pressure drop is computed according
 * the analytical solution for the Newtonian case.
 * Carreau-Yasuda model is provided as a popular alternative model but still
 * remains untested yet.
*/

#include <olb.h>

using namespace olb;
using namespace olb::names;

using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q19<descriptors::FORCE,descriptors::OMEGA>>
>;

enum class ViscosityModel: int {
  NEWTONIAN       = 0,
  POWER_LAW       = 1,
  CASSON          = 2,
  CARREAU_YASUDA  = 3
};

namespace olb::parameters {

struct DYNAMIC_VISCOSITY  : public descriptors::FIELD_BASE<1> {};
struct POWER_LAW_EXPONENT : public descriptors::FIELD_BASE<1> {};
struct K0                 : public descriptors::FIELD_BASE<1> {};
struct K1                 : public descriptors::FIELD_BASE<1> {};
struct N_CY               : public descriptors::FIELD_BASE<1> {};
struct MODEL_CONSTANT_A   : public descriptors::FIELD_BASE<1> {};
struct CHAR_TIME_CONSTANT : public descriptors::FIELD_BASE<1> {};
struct MU_ZERO            : public descriptors::FIELD_BASE<1> {};
struct MU_INF             : public descriptors::FIELD_BASE<1> {};
struct N_PL               : public descriptors::FIELD_BASE<1> {};
struct CONSISTENCY_INDEX  : public descriptors::FIELD_BASE<1> {};
struct DIAMETER           : public descriptors::FIELD_BASE<1> {};
struct EOC                : public descriptors::TYPED_FIELD_BASE<bool,1> {};
struct VISCOSITY_MODEL    : public descriptors::TYPED_FIELD_BASE<ViscosityModel,1> {};
struct AXIS               : public descriptors::FIELD_BASE<0,1> {};
struct EOC_RESOLUTIONS    : public descriptors::TYPED_FIELD_BASE<int,4> {};

} // namespace olb::parameters

Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters)
{
  using T             = MyCase::value_t;
  const T physDeltaX  = parameters.get<parameters::PHYS_DELTA_X>();
  const T length      = parameters.get<parameters::DOMAIN_L>();
  const T radius      = parameters.get<parameters::RADIUS>();

  Vector<T, 3> center0(T(0), radius, radius);
  Vector<T, 3> center1(length + 0.5 * physDeltaX, radius, radius);
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  IndicatorLayer3D<T> extendedDomain(pipe, physDeltaX);

  Mesh<T, MyCase::d> mesh(extendedDomain, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({true, false, false});
  return mesh;
}

struct Newtonian
{
  using T           = MyCase::value_t;
  using DESCRIPTOR  = MyCase::descriptor_t_of<NavierStokes>;

  T charVelocity(MyCase& myCase) {
    auto& parameters  = myCase.getParameters();
    T Re              = parameters.get<parameters::REYNOLDS>();
    T viscosity       = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
    T diameter        = parameters.get<parameters::DIAMETER>();
    return Re * viscosity / diameter;
  }

  std::unique_ptr<AnalyticalF3D<T,T>> getLatticeVelocityProfile(MyCase& myCase) {
    auto& parameters  = myCase.getParameters();
    auto& converter   = myCase.getLattice(NavierStokes{}).getUnitConverter();
    Vector origin     = parameters.get<parameters::ORIGIN>();
    Vector axis       = parameters.get<parameters::AXIS>();
    T     radius      = parameters.get<parameters::RADIUS>();
    return std::make_unique<CirclePoiseuille3D<T>>(origin, axis, converter.getCharLatticeVelocity(), radius);
  }

  std::unique_ptr<AnalyticalF3D<T,T>> getPhysVelocityProfile(MyCase& myCase) {
    auto& parameters  = myCase.getParameters();
    auto& converter   = myCase.getLattice(NavierStokes{}).getUnitConverter();
    Vector origin     = parameters.get<parameters::ORIGIN>();
    Vector axis       = parameters.get<parameters::AXIS>();
    T     radius      = parameters.get<parameters::RADIUS>();
    return std::make_unique<CirclePoiseuille3D<T>>(origin, axis, converter.getCharPhysVelocity(), radius);
  }

  void setModelParameter(MyCase& myCase) {};
};

struct PowerLaw
{
  using T           = MyCase::value_t;
  using DESCRIPTOR  = MyCase::descriptor_t_of<NavierStokes>;

  T charVelocity(MyCase& myCase) {
    auto&   parameters  = myCase.getParameters();
    const T radius      = parameters.get<parameters::RADIUS>();
    const T n           = parameters.get<parameters::POWER_LAW_EXPONENT>();
    T     pressureDrop  = parameters.get<parameters::PRESSURE_DROP>();
    T consistencyIndex  = parameters.get<parameters::CONSISTENCY_INDEX>();
    return ((n)/(n+1))*(util::pow(pressureDrop/(2.*consistencyIndex),1./n))*util::pow(radius,(n+1)/n);
  }

  std::unique_ptr<AnalyticalF3D<T,T>> getLatticeVelocityProfile(MyCase& myCase) {
    auto&   parameters  = myCase.getParameters();
    auto&   converter   = myCase.getLattice(NavierStokes{}).getUnitConverter();
    Vector  origin      = parameters.get<parameters::ORIGIN>();
    Vector  axis        = parameters.get<parameters::AXIS>();
    const T radius      = parameters.get<parameters::RADIUS>();
    const T n           = parameters.get<parameters::POWER_LAW_EXPONENT>();
    return std::make_unique<CirclePowerLaw3D<T>>(origin, axis, converter.getCharLatticeVelocity(), radius, n);
  }

  std::unique_ptr<AnalyticalF3D<T,T>> getPhysVelocityProfile(MyCase& myCase) {
    auto&   parameters  = myCase.getParameters();
    auto&   converter   = myCase.getLattice(NavierStokes{}).getUnitConverter();
    Vector  origin      = parameters.get<parameters::ORIGIN>();
    Vector  axis        = parameters.get<parameters::AXIS>();
    const T radius      = parameters.get<parameters::RADIUS>();
    const T n           = parameters.get<parameters::POWER_LAW_EXPONENT>();
    return std::make_unique<CirclePowerLaw3D<T>>(origin, axis, converter.getCharPhysVelocity(), radius, n);
  }

  void setModelParameter(MyCase& myCase) {
    auto&   parameters        = myCase.getParameters();
    auto&   lattice           = myCase.getLattice(NavierStokes{});
    auto&   converter         = lattice.getUnitConverter();
    const T n                 = parameters.get<parameters::POWER_LAW_EXPONENT>();
    const T consistencyIndex  = parameters.get<parameters::CONSISTENCY_INDEX>();
    const T physDeltaT        = converter.getPhysDeltaT();
    const T physDeltaX        = converter.getPhysDeltaX();
    const T physRho           = parameters.get<parameters::PHYS_CHAR_DENSITY>();
    const T nuMin             = 2.9686e-3;
    const T nuMax             = 3.1667;
    T conversionConsistency   = util::pow(physDeltaT,n-2)*util::pow(physDeltaX,2.);
    lattice.setParameter<powerlaw::M>(consistencyIndex / (conversionConsistency * physRho));
    lattice.setParameter<powerlaw::N>(n);
    lattice.setParameter<powerlaw::OMEGA_MIN>(1./(nuMax*descriptors::invCs2<T,DESCRIPTOR>() + 0.5));
    lattice.setParameter<powerlaw::OMEGA_MAX>(1./(nuMin*descriptors::invCs2<T,DESCRIPTOR>() + 0.5));
  };
};

struct Casson
{
  using T           = MyCase::value_t;
  using DESCRIPTOR  = MyCase::descriptor_t_of<NavierStokes>;

  T charVelocity(MyCase& myCase) {
    auto&   parameters  = myCase.getParameters();
    const T Re          = parameters.get<parameters::REYNOLDS>();
    const T viscosity   = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
    const T diameter    = parameters.get<parameters::DIAMETER>();
    return  Re * viscosity / diameter;
  }

  std::unique_ptr<AnalyticalF3D<T,T>> getLatticeVelocityProfile(MyCase& myCase) {
    auto&   parameters  = myCase.getParameters();
    auto&   converter   = myCase.getLattice(NavierStokes{}).getUnitConverter();
    Vector  origin      = parameters.get<parameters::ORIGIN>();
    Vector  axis        = parameters.get<parameters::AXIS>();
    const T radius      = parameters.get<parameters::RADIUS>();
    const T k0          = parameters.get<parameters::K0>();
    const T k1          = parameters.get<parameters::K1>();
    const T pressureDrop = parameters.get<parameters::PRESSURE_DROP>();
    return std::make_unique<CircleCasson3D<T>>(origin, axis, radius, k1*k1, pressureDrop, k0*k0, 1.0/converter.getConversionFactorVelocity());
  }

  std::unique_ptr<AnalyticalF3D<T,T>> getPhysVelocityProfile(MyCase& myCase) {
    auto&   parameters  = myCase.getParameters();
    Vector  origin      = parameters.get<parameters::ORIGIN>();
    Vector  axis        = parameters.get<parameters::AXIS>();
    const T radius      = parameters.get<parameters::RADIUS>();
    const T k0          = parameters.get<parameters::K0>();
    const T k1          = parameters.get<parameters::K1>();
    const T pressureDrop = parameters.get<parameters::PRESSURE_DROP>();
    return std::make_unique<CircleCasson3D<T>>(origin, axis, radius, k1*k1, pressureDrop, k0*k0, 1.0);
  }

  void setModelParameter(MyCase& myCase) {
    auto&   lattice       = myCase.getLattice(NavierStokes {});
    auto&   converter     = lattice.getUnitConverter();
    auto&   parameters    = myCase.getParameters();
    const T k0            = parameters.get<parameters::K0>();
    const T k1            = parameters.get<parameters::K1>();
    const T physRho       = parameters.get<parameters::PHYS_CHAR_DENSITY>();
    const T conversionK1  = util::sqrt( converter.getConversionFactorViscosity());
    const T conversionK0  = conversionK1 / util::sqrt(converter.getConversionFactorTime());
    lattice.setParameter<visco::K_ZERO>(k0 / (util::sqrt(physRho) * conversionK0));
    lattice.setParameter<visco::K_ONE>(k1 / (util::sqrt(physRho) * conversionK1));
  };
};

struct CarreauYasuda
{
  using T           = MyCase::value_t;
  using DESCRIPTOR  = MyCase::descriptor_t_of<NavierStokes>;

  T charVelocity(MyCase& myCase) {
    auto&   parameters  = myCase.getParameters();
    const T Re          = parameters.get<parameters::REYNOLDS>();
    const T mu_zero     = parameters.get<parameters::MU_ZERO>();
    const T diameter    = parameters.get<parameters::DIAMETER>();
    const T physRho     = parameters.get<parameters::PHYS_CHAR_DENSITY>();
    return  Re * mu_zero / (diameter * physRho);
  }

  std::unique_ptr<AnalyticalF3D<T,T>> getLatticeVelocityProfile(MyCase& myCase) {
    auto&   parameters  = myCase.getParameters();
    auto&   converter   = myCase.getLattice(NavierStokes{}).getUnitConverter();
    Vector  origin      = parameters.get<parameters::ORIGIN>();
    Vector  axis        = parameters.get<parameters::AXIS>();
    const T radius      = parameters.get<parameters::RADIUS>();
    const T n_PL        = parameters.get<parameters::N_PL>();
    return std::make_unique<CirclePowerLaw3D<T>>(origin, axis, converter.getCharLatticeVelocity(), radius, n_PL);
  }

  std::unique_ptr<AnalyticalF3D<T,T>> getPhysVelocityProfile(MyCase& myCase) {
    auto&   parameters  = myCase.getParameters();
    auto&   converter   = myCase.getLattice(NavierStokes{}).getUnitConverter();
    Vector  origin      = parameters.get<parameters::ORIGIN>();
    Vector  axis        = parameters.get<parameters::AXIS>();
    const T radius      = parameters.get<parameters::RADIUS>();
    const T n_PL        = parameters.get<parameters::N_PL>();
    return std::make_unique<CirclePowerLaw3D<T>>(origin, axis, converter.getCharPhysVelocity(), radius, n_PL);
  }

  void setModelParameter(MyCase& myCase) {
    OstreamManager clout(std::cout, "setModelParameter");
    auto&   lattice       = myCase.getLattice(NavierStokes {});
    auto&   converter     = lattice.getUnitConverter();
    const T conversionMU  = converter.getConversionFactorViscosity();
    auto&   parameters    = myCase.getParameters();
    const T n_CY          = parameters.get<parameters::N_CY>();
    const T a             = parameters.get<parameters::MODEL_CONSTANT_A>();
    const T mu_zero       = parameters.get<parameters::MU_ZERO>();
    const T lambda        = parameters.get<parameters::CHAR_TIME_CONSTANT>();
    const T physRho       = parameters.get<parameters::PHYS_CHAR_DENSITY>();
    const T mu_inf        = parameters.get<parameters::MU_INF>();
    lattice.setParameter<visco::N>(n_CY);
    lattice.setParameter<visco::A>(a);
    lattice.setParameter<visco::LAMBDA>(lambda / converter.getConversionFactorTime());
    lattice.setParameter<visco::MU_ZERO>(mu_zero / (physRho * conversionMU));
    lattice.setParameter<visco::MU_INF>(mu_inf / (physRho * conversionMU));
  }
};

void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;
  using T             = MyCase::value_t;
  auto&   geometry    = myCase.getGeometry();
  auto&   parameters  = myCase.getParameters();
  const T physDeltaX  = parameters.get<parameters::PHYS_DELTA_X>();
  const T length      = parameters.get<parameters::DOMAIN_L>();
  const T radius      = parameters.get<parameters::RADIUS>();

  Vector<T, 3> center0(-physDeltaX * 0.2, radius, radius);
  Vector<T, 3> center1(length, radius, radius);
  center0[0] -= 3.*physDeltaX;
  center1[0] += 3.*physDeltaX;
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  geometry.rename(0, 2);
  geometry.rename(2, 1, pipe);
  geometry.clean();
  geometry.innerClean();
  geometry.checkForErrors();
  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

template<typename MODEL>
void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;
  using T           = MyCase::value_t;
  using DESCRIPTOR  = MyCase::descriptor_t_of<NavierStokes>;
  auto& geometry    = myCase.getGeometry();
  auto& lattice     = myCase.getLattice(NavierStokes {});
  auto& parameters  = myCase.getParameters();
  const T length    = parameters.get<parameters::DOMAIN_L>();
  const T radius    = parameters.get<parameters::RADIUS>();
  const T physRho   = parameters.get<parameters::PHYS_CHAR_DENSITY>();
  const T pressureDrop = parameters.get<parameters::PRESSURE_DROP>();
  const size_t N    = parameters.get<parameters::RESOLUTION>();
  const T tau       = parameters.get<parameters::LATTICE_RELAXATION_TIME>();
  const T diameter  = parameters.get<parameters::DIAMETER>();
  const T viscosity = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();

  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR>>(
    (T)   N,
    (T)   tau,
    (T)   diameter,
    (T)   MODEL{}.charVelocity( myCase ),
    (T)   viscosity,
    (T)   physRho
  );
  auto& converter   = lattice.getUnitConverter();
  converter.print();
  converter.write("nonNewtonianPoiseuille3d");

  switch ( parameters.get<parameters::VISCOSITY_MODEL>() ) {
    case ViscosityModel::NEWTONIAN:
    default:
      dynamics::set<ForcedBGKdynamics<T,DESCRIPTOR>>(lattice, geometry.getMaterialIndicator({1}));
      break;
    case ViscosityModel::POWER_LAW:
      dynamics::set<PowerLawForcedBGKdynamics<T,DESCRIPTOR>>(lattice, geometry.getMaterialIndicator({1}));
      break;
    case ViscosityModel::CASSON:
      dynamics::set<CassonForcedBGKdynamics<T,DESCRIPTOR>>(lattice, geometry.getMaterialIndicator({1}));
      break;
    case ViscosityModel::CARREAU_YASUDA:
      dynamics::set<CarreauYasudaForcedBGKdynamics<T,DESCRIPTOR>>(lattice, geometry.getMaterialIndicator({1}));
      break;
  }

  Vector<T, 3> center0(0, radius, radius);
  Vector<T, 3> center1(length, radius, radius);
  center0[0] -= 0.5*converter.getPhysDeltaX();
  center1[0] += 0.5*converter.getPhysDeltaX();
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  setBouzidiBoundary<T, DESCRIPTOR, BouzidiPostProcessor>(lattice, geometry, 2, pipe);

  Vector<T,3> forceTerm = {0,0,0};
  T conversionF_force = converter.getConversionFactorLength() / util::pow(converter.getConversionFactorTime(), 2.);
  forceTerm[0] = pressureDrop / physRho / conversionF_force;

  AnalyticalConst3D<T,T> force (forceTerm);
  fields::set<descriptors::FORCE>(lattice, geometry.getMaterialIndicator({1}), force);
  fields::set<descriptors::FORCE>(lattice, geometry.getMaterialIndicator({2}), force);

  clout << "Prepare Lattice ... OK" << std::endl;
}

template<typename MODEL>
void setInitialValues(MyCase& myCase) {
  OstreamManager clout(std::cout, "setInitialValues");
  clout << "Prepare setInitialValues ..." << std::endl;
  using T           = MyCase::value_t;
  auto& geometry    = myCase.getGeometry();
  auto& lattice     = myCase.getLattice(NavierStokes {});
  auto& converter   = lattice.getUnitConverter();

  std::unique_ptr<AnalyticalF3D<T,T>> profileU = MODEL{}.getPhysVelocityProfile( myCase );
  momenta::setVelocity(lattice, geometry.getMaterialIndicator({1}), *profileU);
  momenta::setVelocity(lattice, geometry.getMaterialIndicator({2}), *profileU);

  lattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  MODEL{}.setModelParameter( myCase );

  lattice.initialize();

  clout << "Prepare setInitialValues ... OK" << std::endl;
}

// Compute error norms
template<typename MODEL, typename T>
Vector<T,3> convergenceCheck( MyCase& myCase)
{
  OstreamManager clout( std::cout,"error" );
  auto& geometry    = myCase.getGeometry();
  auto& lattice     = myCase.getLattice(NavierStokes {});
  auto& converter   = lattice.getUnitConverter();
  using DESCRIPTOR  = MyCase::descriptor_t_of<NavierStokes>;
  Vector<T,3> errors;
  int tmp[]= { };
  T result[2]= { };

  lattice.setProcessingContext(ProcessingContext::Evaluation);
  std::unique_ptr<AnalyticalF3D<T,T>> uSol = MODEL{}.getPhysVelocityProfile( myCase );
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> u( lattice,converter );
  auto indicatorF = geometry.getMaterialIndicator(1);

  SuperAbsoluteErrorL1Norm3D<T> absVelocityErrorNormL1(u, *uSol, indicatorF);
  absVelocityErrorNormL1(result, tmp);
  clout << "velocity-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm3D<T> relVelocityErrorNormL1(u, *uSol, indicatorF);
  relVelocityErrorNormL1(result, tmp);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;
  errors[0] = result[0];

  SuperAbsoluteErrorL2Norm3D<T> absVelocityErrorNormL2(u, *uSol, indicatorF);
  absVelocityErrorNormL2(result, tmp);
  clout << "velocity-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm3D<T> relVelocityErrorNormL2(u, *uSol, indicatorF);
  relVelocityErrorNormL2(result, tmp);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;
  errors[1] = result[0];

  SuperAbsoluteErrorLinfNorm3D<T> absVelocityErrorNormLinf(u, *uSol, indicatorF);
  absVelocityErrorNormLinf(result, tmp);
  clout << "velocity-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm3D<T> relVelocityErrorNormLinf(u, *uSol, indicatorF);
  relVelocityErrorNormLinf(result, tmp);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;
  errors[2] = result[0];

  return errors;
}

template<typename MODEL>
void getResults(MyCase& myCase, std::size_t iT, util::Timer<MyCase::value_t>& timer, bool converged=false)
{
  OstreamManager clout( std::cout,"getResults" );
  using T               = MyCase::value_t;
  using DESCRIPTOR      = MyCase::descriptor_t_of<NavierStokes>;
  auto&       lattice   = myCase.getLattice(NavierStokes {});
  const auto& converter = lattice.getUnitConverter();
  auto&       parameters= myCase.getParameters();
  T           maxPhysT  = parameters.get<parameters::PHYS_MAX_T>();
  const bool lastTimeStep = (iT + 1 == converter.getLatticeTime( maxPhysT ));
  const int statIter = converter.getLatticeTime( maxPhysT / 10. );

  SuperVTMwriter3D<T> vtmWriter( "nonNewtonianPoiseuille3d" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( lattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( lattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  std::unique_ptr<AnalyticalF3D<T,T>> uSol = MODEL{}.getPhysVelocityProfile( myCase );
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> analyticalVelocityLattice(*uSol, lattice);
  analyticalVelocityLattice.getName() = "analytical solution";
  vtmWriter.addFunctor(analyticalVelocityLattice);

  SuperLatticeField3D<T,DESCRIPTOR,descriptors::OMEGA> omega(lattice);
  omega.getName() = "omega";
  vtmWriter.addFunctor(omega);

  const int vtmIter  = converter.getLatticeTime( maxPhysT/100. );

  if ( iT==0 ) {
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( lattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( lattice );

    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }

  // Writes the vtm files and profile text file
  if ( iT%vtmIter==0 || lastTimeStep ) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write( iT );
  }

  // Writes output on the console
  if ( iT%statIter==0 || lastTimeStep ) {
    timer.update( iT );
    timer.printStep();
    lattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
    convergenceCheck<MODEL,T>( myCase );
  }
}

template<typename MODEL, typename T>
Vector<T,3> simulatePoiseuilleWith( MyCase& myCase)
{
  OstreamManager clout( std::cout,"simulatePoiseuille" );
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  auto& converter = lattice.getUnitConverter();
  T maxPhysT = parameters.get<parameters::PHYS_MAX_T>();

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), geometry.getStatistics().getNvoxel() );
  timer.start();

  for ( std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {
    lattice.collideAndStream();
    getResults<MODEL>( myCase, iT, timer );
  }
  timer.stop();
  timer.printSummary();

  return convergenceCheck<MODEL,T>( myCase );
}

void setGetParameters( MyCase::ParametersD& myCaseParameters, int& argc, char** argv ) {

  using namespace olb::parameters;
  using T = MyCase::value_t;

  myCaseParameters.set<EOC>(false);
  myCaseParameters.set<VISCOSITY_MODEL>(ViscosityModel::NEWTONIAN);
  myCaseParameters.set<RESOLUTION>(31);
  myCaseParameters.set<EOC_RESOLUTIONS>({41, 61, 81, 101});
  myCaseParameters.set<DOMAIN_EXTENT>({3e-3,1e-3,1e-3});
  myCaseParameters.set<DOMAIN_L>([] (MyCase::ParametersD& p) {
    return p.get<DOMAIN_EXTENT>()[0];
  });
  myCaseParameters.set<DIAMETER>([] (MyCase::ParametersD& p) {
    return p.get<DOMAIN_EXTENT>()[1];
  });
  myCaseParameters.set<RADIUS>([] (MyCase::ParametersD& p) {
    return p.get<DIAMETER>() / 2.;
  });
  myCaseParameters.set<PHYS_CHAR_LENGTH>([] (MyCase::ParametersD& p) {
    return p.get<DIAMETER>();
  });
  myCaseParameters.set<PHYS_DELTA_X>([] (MyCase::ParametersD& p) {
    return p.get<PHYS_CHAR_LENGTH>() / p.get<RESOLUTION>();
  });
  myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.64265);
  myCaseParameters.set<PHYS_CHAR_DENSITY>(1060);
  myCaseParameters.set<PHYS_MAX_T>(.3);
  myCaseParameters.set<AXIS>({1.,0.,0.});
  myCaseParameters.set<ORIGIN>([] (MyCase::ParametersD& p) {
    T l = p.get<DOMAIN_L>();
    T r = p.get<RADIUS>();
    return Vector(l,r,r);
  });
  myCaseParameters.set<REYNOLDS>(100);
  myCaseParameters.set<DYNAMIC_VISCOSITY>(0.0040755);
  myCaseParameters.set<PHYS_CHAR_VISCOSITY>([] (MyCase::ParametersD& p) {
    return p.get<DYNAMIC_VISCOSITY>() / p.get<PHYS_CHAR_DENSITY>();
  });
  myCaseParameters.set<PRESSURE_DROP>([] (MyCase::ParametersD& p) {
    return (16.*util::pow(p.get<DYNAMIC_VISCOSITY>(),2.)*p.get<REYNOLDS>())/
      (p.get<PHYS_CHAR_DENSITY>()*util::pow(p.get<DOMAIN_EXTENT>()[1],3.));
  });

  // PowerLaw model parameter
  myCaseParameters.set<POWER_LAW_EXPONENT>(0.65);
  myCaseParameters.set<CONSISTENCY_INDEX>([] (MyCase::ParametersD& p) {
    T dynamicViscosity  = p.get<DYNAMIC_VISCOSITY>();
    T n                 = p.get<POWER_LAW_EXPONENT>();
    T pressureDrop      = p.get<PRESSURE_DROP>();
    T radius            = p.get<RADIUS>();
    return  util::pow(dynamicViscosity,n)
            * util::pow((radius*pressureDrop/8.)
                        * util::pow(n/(n+1),n)
                        , 1-n);
  });  // consistency index in Pa s^n, SI unit

  // Casson model parameter
  myCaseParameters.set<K0>(0.07);  // k0 constant (Pa)^0.5, k0^2 corresponds to yield stress
  myCaseParameters.set<K1>([] (MyCase::ParametersD& p) {
    return  util::sqrt(p.get<DYNAMIC_VISCOSITY>());
  });  // k1 constant (Pa*s)^0.5, k1^2 corresponds to Casson viscosity

  // CarreauYasuda model parameter
  myCaseParameters.set<N_CY>(0.9);  // CY fluid index
  myCaseParameters.set<MODEL_CONSTANT_A>(1.5);  // model constant
  myCaseParameters.set<CHAR_TIME_CONSTANT>(3.313);  // characteristic time constant lambda (s)
  myCaseParameters.set<MU_ZERO>([] (MyCase::ParametersD& p) {
    return p.get<DYNAMIC_VISCOSITY>()*1.45;
  });  // zero viscosity (Pa*s)
  myCaseParameters.set<MU_INF>([] (MyCase::ParametersD& p) {
    return p.get<DYNAMIC_VISCOSITY>()/10.;
  });  // infinity viscosity (Pa*s)
  myCaseParameters.set<N_PL>(0.708);  // PL fluid index

}
