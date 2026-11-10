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
// Knudsen pump based on the temperature gradeint along the wall between two vessels
// The thermal creep flow leads to pressure gradient, which forces the fluid to stream backwards
// in the center of the pipe. The total MFR is therefore 0.
// Analytical solution is valid for b = 0 at Kn < 0.5

#include <olb.h>

using namespace olb;
using namespace olb::names;

// Anayltical solution for velocity profile
template <typename T, typename S>
class PressureProfile3D : public AnalyticalF3D<T, S>
{
private:
    T Rspec, radius, length, mu, pM, slipCoeff, xStart, xEnd, tempLeft, tempRight;

public:
    PressureProfile3D(T Rspec_, T radius_, T length_, T mu_, T pM_, T slipCoeff_, T xStart_, T xEnd_, T tempLeft_, T tempRight_) : AnalyticalF3D<T, S>(1),
           Rspec(Rspec_), radius(radius_), length(length_), mu(mu_), pM(pM_), slipCoeff(slipCoeff_), xStart(xStart_), xEnd(xEnd_), tempLeft(tempLeft_), tempRight(tempRight_)
    {
        this->getName() = "PressureProfile3D";
    };

    bool operator()(T output[1], const S x[3])
    {
        T dP = 6.*mu*mu*Rspec/ radius/ radius/ pM * (tempRight - tempLeft)/(1. + 8.*slipCoeff/(2*radius));
        if(x[0] < xStart) {
          output[0] = pM-0.5*dP;
        }
        else if(x[0] > xEnd) {
          output[0] = pM+0.5*dP;
        }
        else{
          output[0] = pM-0.5*dP + dP/(xEnd - xStart)*(x[0] - xStart);
        }
        return true;
    };
};

// Temperature profile
template <typename T, typename S>
class TemperatureProfile3D : public AnalyticalF3D<T, S>
{
private:
    T tempLeft, tempRight, xStart, xEnd;

public:
    TemperatureProfile3D(T tempLeft_, T tempRight_, T xStart_, T xEnd_) : AnalyticalF3D<T, S>(1),
           tempLeft(tempLeft_), tempRight(tempRight_), xStart(xStart_), xEnd(xEnd_)
    {
        this->getName() = "TemperatureProfile3D";
    };

    bool operator()(T output[1], const S x[3])
    {
        if(x[0] < xStart) {
          output[0] = tempLeft;
        }
        else if(x[0] > xEnd) {
          output[0] = tempRight;
        }
        else{
          output[0] = tempLeft + (tempRight - tempLeft)/(xEnd - xStart)*(x[0] - xStart);
        }

        return true;
    };
};

using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q19<descriptors::NORMAL,descriptors::VELOCITY,
                                                   descriptors::TEMPERATURE,descriptors::TEMPGRADIENT>>
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
struct ERROR_NORMS : public descriptors::FIELD_BASE<3> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const T physDeltaX = extend[1]/parameters.get<parameters::RESOLUTION>();

  Vector<T, 3> center0(0, extend[1]/2.-0.5 * physDeltaX, extend[1]/2.-0.5 * physDeltaX);
  Vector<T, 3> center1(extend[0], extend[1]/2.-0.5 * physDeltaX, extend[1]/2.-0.5 * physDeltaX);
  std::shared_ptr<IndicatorCylinder3D<T> > pipe( new IndicatorCylinder3D<T>(center0, center1, extend[1]/2.));
  Vector<T, 3> centerL(-extend[0]/2., extend[1]/2.-0.5 * physDeltaX, extend[1]/2.-0.5 * physDeltaX);
  std::shared_ptr<IndicatorCylinder3D<T> > tankL( new IndicatorCylinder3D<T>(centerL, center0, extend[1]));
  Vector<T, 3> centerR(1.5*extend[0], extend[1]/2.-0.5 * physDeltaX, extend[1]/2.-0.5 * physDeltaX);
  std::shared_ptr<IndicatorCylinder3D<T> > tankR( new IndicatorCylinder3D<T>(center1, centerR, extend[1]));
  auto knudsenPump = tankL + pipe + tankR;
  IndicatorLayer3D<T> knudsenPumpLayer( *knudsenPump, 2.*physDeltaX );

  Mesh<T,MyCase::d> mesh(knudsenPumpLayer, physDeltaX, singleton::mpi().getSize());
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
  Vector<T, 3> center0(0, extend[1]/2.-0.5 * physDeltaX, extend[1]/2.-0.5 * physDeltaX);
  Vector<T, 3> center1(extend[0], extend[1]/2.-0.5 * physDeltaX, extend[1]/2.-0.5 * physDeltaX);
  std::shared_ptr<IndicatorCylinder3D<T> > pipe( new IndicatorCylinder3D<T>(center0, center1, extend[1]/2.));
  Vector<T, 3> centerL(-extend[0]/2., extend[1]/2.-0.5 * physDeltaX, extend[1]/2.-0.5 * physDeltaX);
  std::shared_ptr<IndicatorCylinder3D<T> > tankL( new IndicatorCylinder3D<T>(centerL, center0, extend[1]));
  Vector<T, 3> centerR(1.5*extend[0], extend[1]/2.-0.5 * physDeltaX, extend[1]/2.-0.5 * physDeltaX);
  std::shared_ptr<IndicatorCylinder3D<T> > tankR( new IndicatorCylinder3D<T>(center1, centerR, extend[1]));
  auto knudsenPump = tankL + pipe + tankR;
  IndicatorLayer3D<T> knudsenPumpLayer( *knudsenPump, 2.*physDeltaX );

  geometry.rename(0, 2, knudsenPumpLayer);

  geometry.rename(2, 1, knudsenPump);

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

  TemperatureProfile3D<T,T> temp(tempLeft, tempRight, 0., extend[0]);
  fields::set<descriptors::TEMPERATURE>(lattice, geometry.getMaterialIndicator({0,2}), temp);

  T physDeltaX = lattice.getUnitConverter().getPhysDeltaX();
  Vector<T, 3> center0(0, extend[1]/2.-0.5 * physDeltaX, extend[1]/2.-0.5 * physDeltaX);
  Vector<T, 3> center1(extend[0], extend[1]/2.-0.5 * physDeltaX, extend[1]/2.-0.5 * physDeltaX);
  std::shared_ptr<IndicatorCylinder3D<T> > pipe( new IndicatorCylinder3D<T>(center0, center1, extend[1]/2.));
  Vector<T, 3> centerL(-extend[0]/2., extend[1]/2.-0.5 * physDeltaX, extend[1]/2.-0.5 * physDeltaX);
  std::shared_ptr<IndicatorCylinder3D<T> > tankL( new IndicatorCylinder3D<T>(centerL, center0, extend[1]));
  Vector<T, 3> centerR(1.5*extend[0], extend[1]/2.-0.5 * physDeltaX, extend[1]/2.-0.5 * physDeltaX);
  std::shared_ptr<IndicatorCylinder3D<T> > tankR( new IndicatorCylinder3D<T>(center1, centerR, extend[1]));
  auto knudsenPump = tankL + pipe + tankR;

  setBouzidiBoundary<T,DESCRIPTOR,KnudsenVelocityPostProcessor<true>>(lattice, geometry.getMaterialIndicator({2}), geometry.getMaterialIndicator({1,3,4}), *knudsenPump);
  setBouzidiKnudsenSlipVelocity<T,DESCRIPTOR,true>(lattice, geometry.getMaterialIndicator({2}), geometry.getMaterialIndicator({1,3,4}), *knudsenPump);

  lattice.setParameter<KnudsenVelocityPostProcessor<true>::SLIPCOEFF>(slipCoeff/lattice.getUnitConverter().getPhysDeltaX());
  lattice.setParameter<KnudsenVelocityPostProcessor<true>::CREEPCOEFF>(creepCoeff/lattice.getUnitConverter().getPhysDeltaX()/lattice.getUnitConverter().getConversionFactorVelocity());

  lattice.setParameter<descriptors::OMEGA>(lattice.getUnitConverter().getLatticeRelaxationFrequency());

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
  const T pM = parameters.get<parameters::AVERAGE_PRESSURE>();
  T pL0 = lattice.getUnitConverter().getLatticePressure(pM);
  AnalyticalConst3D<T,T> rho(pL0*descriptors::invCs2<T,DESCRIPTOR>()+T(1));

  momenta::setPressure(lattice, geometry.getMaterialIndicator({0,1,2}), pL0);

  lattice.initialize();
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
  const T avPressure = parameters.get<parameters::AVERAGE_PRESSURE>();
  const T slipCoeff = parameters.get<parameters::SLIPCOEFF>();
  const T physCharViscosity = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T physCharDensity = parameters.get<parameters::PHYS_CHAR_DENSITY>();
  const T Rspec = parameters.get<parameters::R_SPECIFIC>();
  const T tempLeft = parameters.get<parameters::TEMPERATURE_LEFT>();
  const T tempRight = parameters.get<parameters::TEMPERATURE_RIGHT>();
  PressureProfile3D<T,T> presSol(Rspec, extend[1]/2., extend[0], physCharDensity*physCharViscosity, avPressure, slipCoeff, 0., extend[0], tempLeft, tempRight);
  SuperLatticePhysPressure3D<T,DESCRIPTOR> pres( lattice,lattice.getUnitConverter() );
  auto indicatorF = geometry.getMaterialIndicator(1);

  Vector<T,3> errors;
  SuperRelativeErrorL1Norm3D<T> relPressureErrorNormL1(pres, presSol, indicatorF);
  relPressureErrorNormL1(result, tmp);
  clout << "; pressure-L1-error(rel)=" << result[0] << std::endl;
  errors[0] = result[0];

  SuperRelativeErrorL2Norm3D<T> relPressureErrorNormL2(pres, presSol, indicatorF);
  relPressureErrorNormL2(result, tmp);
  clout << "; pressure-L2-error(rel)=" << result[0] << std::endl;
  errors[1] = result[0];

  SuperRelativeErrorLinfNorm3D<T> relPressureErrorNormLinf(pres, presSol, indicatorF);
  relPressureErrorNormLinf(result, tmp);
  clout << "; pressure-Linf-error(rel)=" << result[0] << std::endl;
  errors[2] = result[0];
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
  const int statIter = lattice.getUnitConverter().getLatticeTime( maxPhysT/50. );

  // VTK and image output only if no EOC analysis
  SuperVTMwriter3D<T> vtmWriter( "knudsenPump3d" );
  const int vtmIter  = lattice.getUnitConverter().getLatticeTime( maxPhysT/50. );

  if ( iT==0 ) {
    // Writes the geometry, pipe no. and rank no. as vti file for visualization
    vtmWriter.createMasterFile();
  }

  // Writes the vtm files and profile text file
  if ( iT%vtmIter==0 || lastTimeStep ) {
    lattice.executePostProcessors(stage::Evaluation{});
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    lattice.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriter("knudsenPump3d");
      SuperLatticePhysVelocity3D velocity(lattice, lattice.getUnitConverter());
      SuperLatticePhysPressure3D pressure(lattice, lattice.getUnitConverter());
      SuperGeometryF<T,3> geom( geometry );
      vtmWriter.addFunctor(velocity);
      vtmWriter.addFunctor(pressure);
      vtmWriter.addFunctor(geom);
      const Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
      const T avPressure = parameters.get<parameters::AVERAGE_PRESSURE>();
      const T slipCoeff = parameters.get<parameters::SLIPCOEFF>();
      const T physCharViscosity = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
      const T physCharDensity = parameters.get<parameters::PHYS_CHAR_DENSITY>();
      const T Rspec = parameters.get<parameters::R_SPECIFIC>();
      const T tempLeft = parameters.get<parameters::TEMPERATURE_LEFT>();
      const T tempRight = parameters.get<parameters::TEMPERATURE_RIGHT>();
      PressureProfile3D<T,T> presSol(Rspec, extend[1]/2., extend[0], physCharDensity*physCharViscosity, avPressure, slipCoeff, 0., extend[0], tempLeft, tempRight);
      SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> analyticalPressureLattice(presSol, lattice);
      analyticalPressureLattice.getName() = "analytical pressure";
      vtmWriter.addFunctor(analyticalPressureLattice);
      SuperLatticeField3D<T,DESCRIPTOR,descriptors::TEMPERATURE> temperatureF( lattice );
      temperatureF.getName() = "temperatureF";
      vtmWriter.addFunctor( temperatureF );
      task(vtmWriter, iT);
    });
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

    myCase.getLattice(NavierStokes{}).collideAndStream();
    myCase.getLattice(NavierStokes{}).stripeOffDensityOffset(myCase.getLattice(NavierStokes{}).getStatistics().getAverageRho()-(T)1);

    getResults(myCase, timer, iT);
    converge.takeValue( myCase.getLattice(NavierStokes{}).getStatistics().getMaxU(), false );
  }

  timer.stop();
  timer.printSummary();
  singleton::pool().wait();
}

int main(int argc, char* argv[]) {
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<KNUDSEN>(0.45);
    myCaseParameters.set<B_COEFF>(-0.05);
    myCaseParameters.set<ACCOMODATION_COEFF>(1);
    myCaseParameters.set<RESOLUTION>(41);
    myCaseParameters.set<OVERLAP>(5);
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.9);
    myCaseParameters.set<PHYS_MAX_T>(2.e-6);
    myCaseParameters.set<PHYS_START_T>([] (MyCase::ParametersD& p) {
      return 0.05*p.get<parameters::PHYS_MAX_T>();
    });
    myCaseParameters.set<KINETIC_DIAMETER>(364e-12); //kinetic diamter of nitrogen [m]
    myCaseParameters.set<AVERAGE_PRESSURE>(1000);
    myCaseParameters.set<PRESSURE_DROP>(0.5);
    myCaseParameters.set<MOLAR_MASS>(28.96e-3);
    myCaseParameters.set<R_SPECIFIC>(296.8);
    myCaseParameters.set<parameters::TEMPERATURE>(300.);
    myCaseParameters.set<parameters::TEMPERATURE_LEFT>(298.);
    myCaseParameters.set<parameters::TEMPERATURE_RIGHT>(302.);
    myCaseParameters.set<PHYS_CHAR_DENSITY>([] (MyCase::ParametersD& p) {
      return util::idealGasDensity( p.get<MOLAR_MASS>(), p.get<AVERAGE_PRESSURE>(), p.get<parameters::TEMPERATURE>());
    });
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>([] (MyCase::ParametersD& p) {
      return util::gasDynamicViscosity( p.get<MOLAR_MASS>(), p.get<AVERAGE_PRESSURE>(), p.get<parameters::TEMPERATURE>(), p.get<KINETIC_DIAMETER>())/p.get<PHYS_CHAR_DENSITY>();
    });
    myCaseParameters.set<MFP>([] (MyCase::ParametersD& p) {
      return util::meanFreePath(p.get<MOLAR_MASS>(), p.get<AVERAGE_PRESSURE>(), p.get<parameters::TEMPERATURE>(), p.get<PHYS_CHAR_VISCOSITY>()*p.get<PHYS_CHAR_DENSITY>());
    });
    myCaseParameters.set<DOMAIN_EXTENT>([] (MyCase::ParametersD& p) {
      return Vector{4.*p.get<MFP>()/p.get<KNUDSEN>(), p.get<MFP>()/p.get<KNUDSEN>(), p.get<MFP>()/p.get<KNUDSEN>()};
    });
    myCaseParameters.set<SLIPCOEFF>([] (MyCase::ParametersD& p) {
      return (2. - p.get<ACCOMODATION_COEFF>())/p.get<ACCOMODATION_COEFF>()/(1. - p.get<B_COEFF>()*p.get<KNUDSEN>())*p.get<MFP>();
    });
    myCaseParameters.set<CREEPCOEFF>([] (MyCase::ParametersD& p) {
      return 3./4.*p.get<PHYS_CHAR_VISCOSITY>()*p.get<PHYS_CHAR_DENSITY>()*p.get<R_SPECIFIC>()/p.get<AVERAGE_PRESSURE>();
    });
    myCaseParameters.set<PHYS_CHAR_VELOCITY>([] (MyCase::ParametersD& p) {
      return p.get<DOMAIN_EXTENT>()[1]*p.get<DOMAIN_EXTENT>()[1]/2./p.get<PHYS_CHAR_VISCOSITY>()/p.get<PHYS_CHAR_DENSITY>()
      * (-p.get<PRESSURE_DROP>())/p.get<DOMAIN_EXTENT>()[0]
      * (0.25 - 0.5 - p.get<SLIPCOEFF>()/p.get<DOMAIN_EXTENT>()[1])
      + p.get<CREEPCOEFF>()*(p.get<parameters::TEMPERATURE_RIGHT>()
      - p.get<parameters::TEMPERATURE_LEFT>())/p.get<DOMAIN_EXTENT>()[0];
    });
    myCaseParameters.set<PHYS_INTERVAL_T>(2.e-6/50.);
    myCaseParameters.set<RESIDUUM>(1e-6);
    myCaseParameters.set<CONVERGED>(false);
    myCaseParameters.set<ERROR_NORMS>({0,0,0});
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
