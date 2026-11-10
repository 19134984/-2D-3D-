/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Fedor Bukreev, Adrian Kummerländer,
 *  Mathias J. Krause
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

#include <olb.h>

using namespace olb;
using namespace olb::names;

// Anayltical solution for velocity profile
template <typename T, typename S>
class VelocityKnudsenProfile3D : public AnalyticalF3D<T, S>
{
private:
    T height, length, mu, pIn, pOut, accomodC, Kn, slipCoeff, creepCoeff, tempLeft, tempRight;

public:
    VelocityKnudsenProfile3D(T height_, T length_, T mu_, T pIn_, T pOut_, T slipCoeff_, T creepCoeff_, T tempLeft_, T tempRight_) : AnalyticalF3D<T, S>(3),
           height(height_), length(length_), mu(mu_), pIn(pIn_), pOut(pOut_), slipCoeff(slipCoeff_), creepCoeff(creepCoeff_), tempLeft(tempLeft_), tempRight(tempRight_)
    {
        this->getName() = "VelocityKnudsenProfile3D";
    };

    bool operator()(T output[3], const S x[3])
    {
        T radius = util::sqrt((x[1]-height/2.)*(x[1]-height/2.) + (x[2]-height/2.)*(x[2]-height/2.));
        output[0] = height/2.*height/2./T(4)/mu*(pIn-pOut)/length*(1. - radius*radius/height/height*4. + 4.*slipCoeff/height)+creepCoeff*(tempRight-tempLeft)/length;
        output[1] = T(0);
        output[2] = T(0);

        return true;
    };
};

using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q19<descriptors::NORMAL,descriptors::VELOCITY,descriptors::TEMPERATURE,descriptors::TEMPGRADIENT>>
>;

namespace olb::parameters {

// Microfluidics
struct B_COEFF : public descriptors::FIELD_BASE<1> { };
struct ACCOMODATION_COEFF : public descriptors::FIELD_BASE<1> { };
struct KINETIC_DIAMETER : public descriptors::FIELD_BASE<1> { };
struct MOLAR_MASS : public descriptors::FIELD_BASE<1> { };
struct R_SPECIFIC : public descriptors::FIELD_BASE<1> { };
struct TEMPERATURE_LEFT : public descriptors::FIELD_BASE<1> { };
struct TEMPERATURE_RIGHT : public descriptors::FIELD_BASE<1> { };
struct SLIPCOEFF : public descriptors::FIELD_BASE<1> { };
struct CREEPCOEFF : public descriptors::FIELD_BASE<1> { };
struct MFP : public descriptors::FIELD_BASE<1> { };
struct ERROR_NORMS : public descriptors::FIELD_BASE<6> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const T physDeltaX = extend[1]/parameters.get<parameters::RESOLUTION>();

  Vector<T, 3> center0(0, extend[1]/2.-0.5 * physDeltaX, extend[1]/2.-0.5 * physDeltaX);
  Vector<T, 3> center1(extend[0] + 0.5 * physDeltaX, extend[1]/2.-0.5 * physDeltaX, extend[1]/2.-0.5 * physDeltaX);
  IndicatorCylinder3D<T> pipe(center0, center1, extend[1]/2.);
  IndicatorLayer3D<T> pipeLayer( pipe, 3.*physDeltaX );

  Mesh<T,MyCase::d> mesh(pipeLayer, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

// Stores geometry information in form of material numbers
void prepareGeometry(MyCase& myCase) {
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();

  const Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const T physDeltaX = extend[1]/parameters.get<parameters::RESOLUTION>();

  Vector<T, 3> center0(-physDeltaX * 0.2, extend[1]/2., extend[1]/2.);
  Vector<T, 3> center1(extend[0], extend[1]/2., extend[1]/2.);
  IndicatorCylinder3D<T> pipe(center0, center1, extend[1]/2.);

  geometry.rename(0, 2);

  geometry.rename(2, 1, pipe);

  geometry.clean();
  Vector<T, 3> origin(0, extend[1]/2., extend[1]/2.);
  Vector<T, 3> extendI = origin;

  // Set material number for inflow
  origin[0] = -physDeltaX * 2;
  extendI[0] = physDeltaX * 2;
  IndicatorCylinder3D<T> inflow(origin, extendI, extend[1]/2.);
  geometry.rename(2, 3, 1, inflow);

  // Set material number for outflow
  origin[0] = extend[0] - 2 * physDeltaX;
  extendI[0] = extend[0] + 2 * physDeltaX;
  IndicatorCylinder3D<T> outflow(extendI, origin, extend[1]/2.);
  geometry.rename(2, 4, 1, outflow);

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  // Removes all not needed boundary voxels inside the surface
  geometry.innerClean();
  geometry.checkForErrors();

  geometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
  geometry.communicate();
}

// Set up the geometry of the simulation
void prepareLattice(MyCase& myCase) {
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});

  const int N = parameters.get<parameters::RESOLUTION>();
  const T relaxationTime = parameters.get<parameters::LATTICE_RELAXATION_TIME>();
  Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const T physCharVelocity = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T physCharViscosity = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T physCharDensity = parameters.get<parameters::PHYS_CHAR_DENSITY>();
  const T tempLeft = parameters.get<parameters::TEMPERATURE_LEFT>();
  const T tempRight = parameters.get<parameters::TEMPERATURE_RIGHT>();
  const T outletPressure = parameters.get<parameters::AVERAGE_PRESSURE>()-0.5*parameters.get<parameters::PRESSURE_DROP>();
  const T slipCoeff = parameters.get<parameters::SLIPCOEFF>();
  const T creepCoeff = parameters.get<parameters::CREEPCOEFF>();

  // Set up a unit converter with the characteristic physical units
  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR>>(
    N,
    relaxationTime,
    extend[1],        // charPhysLength: reference length of simulation geometry in [m]
    physCharVelocity,  // charPhysVelocity: highest expected velocity during simulation in [m/s]
    physCharViscosity, // physViscosity: physical kinematic viscosity in [m^2/s]
    physCharDensity,    // physDensity: physical density [kg/m^3]
    outletPressure
  );
  lattice.getUnitConverter().print();

  using BulkDynamics = BGKdynamics<T,DESCRIPTOR>::template wrap_collision<collision::SaveVelocity>;
  dynamics::set<BulkDynamics>(lattice, geometry.getMaterialIndicator({1,3,4}));

  boundary::set<boundary::LocalVelocity<T,DESCRIPTOR,BulkDynamics>>(lattice, geometry, 3);
  boundary::set<boundary::LocalPressure<T,DESCRIPTOR,BulkDynamics>>(lattice, geometry, 4);

  AnalyticalLinear3D<T,T> temp((tempRight-tempLeft)/extend[0], 0, 0, tempLeft);
  fields::set<descriptors::TEMPERATURE>(lattice, geometry.getMaterialIndicator({0,2}), temp);

  Vector<T, 3> center02(-10.*lattice.getUnitConverter().getPhysDeltaX(), extend[1]/2., extend[1]/2.);
  Vector<T, 3> center12(extend[0] + 20 * lattice.getUnitConverter().getPhysDeltaX(), extend[1]/2., extend[1]/2.);
  IndicatorCylinder3D<T> pipe(center02, center12, extend[1]/2.);

  setBouzidiBoundary<T,DESCRIPTOR,KnudsenVelocityPostProcessor<true>>(lattice, geometry.getMaterialIndicator({2}), geometry.getMaterialIndicator({1,3,4}), pipe);
  setBouzidiKnudsenSlipVelocity<T,DESCRIPTOR,true>(lattice, geometry.getMaterialIndicator({2}), geometry.getMaterialIndicator({1,3,4}), pipe);

  lattice.template setParameter<KnudsenVelocityPostProcessor<true>::SLIPCOEFF>(slipCoeff/lattice.getUnitConverter().getPhysDeltaX());
  lattice.template setParameter<KnudsenVelocityPostProcessor<true>::CREEPCOEFF>(creepCoeff/lattice.getUnitConverter().getPhysDeltaX()/lattice.getUnitConverter().getConversionFactorVelocity());

  lattice.template setParameter<descriptors::OMEGA>(lattice.getUnitConverter().getLatticeRelaxationFrequency());

  {
    auto& communicator = lattice.getCommunicator(stage::PostStream());
    communicator.requestField<descriptors::POPULATION>();
    communicator.requestField<descriptors::VELOCITY>();
    communicator.requestOverlap(parameters.get<parameters::OVERLAP>());
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});
  const T outletPressure = parameters.get<parameters::AVERAGE_PRESSURE>()-0.5*parameters.get<parameters::PRESSURE_DROP>();
  // Initialize all values of distribution functions to their local equilibrium
  momenta::setPressure( lattice, geometry.getMaterialIndicator({0,1,2,3,4}), outletPressure);
  // Make the lattice ready for simulation
  lattice.initialize();
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{
  OstreamManager clout( std::cout,"setTemporalValues" );
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});

  // No of time steps for smooth start-up
  const T startTime = parameters.get<parameters::PHYS_START_T>();
  int iTmaxStart = lattice.getUnitConverter().getLatticeTime( startTime );

  if ( int(iT)<= iTmaxStart ) {
    PolynomialStartScale<T,T> StartScale( iTmaxStart, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    T iTvec[1] = {T( iT )};
    T frac[1] = {};
    StartScale( frac,iTvec );
    Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
    const T outletPressure = parameters.get<parameters::AVERAGE_PRESSURE>()-0.5*parameters.get<parameters::PRESSURE_DROP>();
    const T inletPressure = outletPressure + parameters.get<parameters::PRESSURE_DROP>();
    const T slipCoeff = parameters.get<parameters::SLIPCOEFF>();
    const T physCharViscosity = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
    const T physCharDensity = parameters.get<parameters::PHYS_CHAR_DENSITY>();
    const T creepCoeff = parameters.get<parameters::CREEPCOEFF>();
    const T tempLeft = parameters.get<parameters::TEMPERATURE_LEFT>();
    const T tempRight = parameters.get<parameters::TEMPERATURE_RIGHT>();
    T pLin = (inletPressure-outletPressure)*frac[0] + outletPressure;
    VelocityKnudsenProfile3D<T,T> uSol(extend[1], extend[0], physCharDensity*physCharViscosity, inletPressure, outletPressure, slipCoeff, creepCoeff, tempLeft, tempRight);

    momenta::setPressure( lattice, geometry.getMaterialIndicator(3), pLin );
    momenta::setVelocity( lattice, geometry.getMaterialIndicator(3), uSol );

    lattice.setProcessingContext<Array<momenta::FixedDensity::RHO>>(
      ProcessingContext::Simulation);
    lattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

// Compute error norms
void evaluateError(MyCase& myCase) {
  OstreamManager clout( std::cout,"error" );
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});
  lattice.setProcessingContext(ProcessingContext::Evaluation);

  int tmp[]= { };
  T result[2]= { };

  // velocity error
  const Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const T outletPressure = parameters.get<parameters::AVERAGE_PRESSURE>()-0.5*parameters.get<parameters::PRESSURE_DROP>();
  const T inletPressure = outletPressure + parameters.get<parameters::PRESSURE_DROP>();
  const T slipCoeff = parameters.get<parameters::SLIPCOEFF>();
  const T physCharViscosity = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T physCharDensity = parameters.get<parameters::PHYS_CHAR_DENSITY>();
  const T creepCoeff = parameters.get<parameters::CREEPCOEFF>();
  const T tempLeft = parameters.get<parameters::TEMPERATURE_LEFT>();
  const T tempRight = parameters.get<parameters::TEMPERATURE_RIGHT>();
  VelocityKnudsenProfile3D<T,T> uSol(extend[1], extend[0], physCharDensity*physCharViscosity, inletPressure, outletPressure, slipCoeff, creepCoeff, tempLeft, tempRight);
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> u( lattice,lattice.getUnitConverter() );
  auto indicatorF = geometry.getMaterialIndicator(1);

  Vector<T,6> errors;
  SuperAbsoluteErrorL1Norm3D<T> absVelocityErrorNormL1(u, uSol, indicatorF);
  absVelocityErrorNormL1(result, tmp);
  clout << "velocity-L1-error(abs)=" << result[0];
  errors[0] = result[0];
  SuperRelativeErrorL1Norm3D<T> relVelocityErrorNormL1(u, uSol, indicatorF);
  relVelocityErrorNormL1(result, tmp);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;
  errors[1] = result[0];

  SuperAbsoluteErrorL2Norm3D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, tmp);
  clout << "velocity-L2-error(abs)=" << result[0];
  errors[2] = result[0];
  SuperRelativeErrorL2Norm3D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, tmp);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;
  errors[3] = result[0];

  SuperAbsoluteErrorLinfNorm3D<T> absVelocityErrorNormLinf(u, uSol, indicatorF);
  absVelocityErrorNormLinf(result, tmp);
  clout << "velocity-Linf-error(abs)=" << result[0];
  errors[4] = result[0];
  SuperRelativeErrorLinfNorm3D<T> relVelocityErrorNormLinf(u, uSol, indicatorF);
  relVelocityErrorNormLinf(result, tmp);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;
  errors[5] = result[0];

  parameters.set<parameters::ERROR_NORMS>(errors);
}

// Output to console and files
void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  OstreamManager clout( std::cout,"getResults" );
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});
  const T maxPhysT = parameters.get<parameters::PHYS_MAX_T>();
  bool hasConverged = parameters.get<parameters::CONVERGED>();
  const bool lastTimeStep = ( hasConverged || (iT + 1 == lattice.getUnitConverter().getLatticeTime( maxPhysT )) );
  const int statIter = lattice.getUnitConverter().getLatticeTime( maxPhysT/100. );

  // VTK and image output only if no EOC analysis
  SuperVTMwriter3D<T> vtmWriter( "microChannel3d" );
  const int vtmIter  = lattice.getUnitConverter().getLatticeTime( maxPhysT/10. );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    vtmWriter.createMasterFile();
  }

  // Writes the vtm files and profile text file
  if ( iT%vtmIter==0 || lastTimeStep ) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    SuperGeometryF<T,3> geom(geometry);
    vtmWriter.addFunctor( geom );

    SuperLatticePhysVelocity3D<T,DESCRIPTOR> velocity( lattice, lattice.getUnitConverter() );
    vtmWriter.addFunctor( velocity );

    SuperLatticePhysPressure3D<T,DESCRIPTOR> pressure( lattice, lattice.getUnitConverter() );
    vtmWriter.addFunctor( pressure );

    const Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
    const T outletPressure = parameters.get<parameters::AVERAGE_PRESSURE>()-0.5*parameters.get<parameters::PRESSURE_DROP>();
    const T inletPressure = outletPressure + parameters.get<parameters::PRESSURE_DROP>();
    const T slipCoeff = parameters.get<parameters::SLIPCOEFF>();
    const T physCharViscosity = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
    const T physCharDensity = parameters.get<parameters::PHYS_CHAR_DENSITY>();
    const T creepCoeff = parameters.get<parameters::CREEPCOEFF>();
    const T tempLeft = parameters.get<parameters::TEMPERATURE_LEFT>();
    const T tempRight = parameters.get<parameters::TEMPERATURE_RIGHT>();
    VelocityKnudsenProfile3D<T,T> uSol(extend[1], extend[0], physCharDensity*physCharViscosity, inletPressure, outletPressure, slipCoeff, creepCoeff, tempLeft, tempRight);
    SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> analyticalVelocityLattice(uSol, lattice);
    lattice.communicate();
    analyticalVelocityLattice.getName() = "analytical solution";
    vtmWriter.addFunctor(analyticalVelocityLattice);

    SuperLatticeField3D<T,DESCRIPTOR,descriptors::TEMPERATURE> temperatureF( lattice );
    temperatureF.getName() = "temperatureF";
    vtmWriter.addFunctor( temperatureF );
    vtmWriter.write();
  }

  // Output on the console
  if ( iT%statIter==0 || lastTimeStep ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    lattice.getStatistics().print( iT,lattice.getUnitConverter().getPhysTime( iT ) );

    // Error norms
    if ( lastTimeStep ) {
      evaluateError(myCase);
    }
  }
}

void simulate(MyCase& myCase) {
  OstreamManager clout( std::cout,"simulate" );
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  const T maxPhysT = parameters.get<parameters::PHYS_MAX_T>();
  const T physInterval = parameters.get<parameters::PHYS_INTERVAL_T>();
  const T residuum = parameters.get<parameters::RESIDUUM>();

  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(maxPhysT);
  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();
  util::ValueTracer<T> converge( myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime( physInterval ), residuum );
  timer.start();

  parameters.set<parameters::CONVERGED>(false);
  for (std::size_t iT=0; iT < iTmax; ++iT) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << std::endl;
      parameters.set<parameters::CONVERGED>(true);
      getResults(myCase, timer, iT);

      break;
    }

    setTemporalValues(myCase, iT);

    myCase.getLattice(NavierStokes{}).collideAndStream();

    getResults(myCase, timer, iT);
    converge.takeValue( myCase.getLattice(NavierStokes{}).getStatistics().getMaxU(), false );
  }

  timer.stop();
  timer.printSummary();
  singleton::pool().wait();
}
