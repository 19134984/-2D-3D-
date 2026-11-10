/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2025 Fedor Bukreev, Adrian Kummerländer,
 *  Shota Ito, Mathias J. Krause
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

 /* pipeWithValve3d.cpp:
 * This example represents flow in a pipe with rotated valve inside.
 * Pipe and valve geometries are read from STL files.
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;

using MyCase = Case<
  NavierStokes, Lattice<FLOATING_POINT_TYPE, descriptors::D3Q19<>>
>;

namespace olb::parameters {

struct ROTATION_POINT  : public descriptors::FIELD_BASE<0,1> { };
struct ROTATION_AXIS : public descriptors::FIELD_BASE<0,1> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  return Mesh<MyCase::value_t,MyCase::d>::fromSTL(parameters);
}

// Stores data from stl file in geometry in form of material numbers
void prepareGeometry(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();

  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();

  auto pipeI = myCase.getMesh().getIndicator("pipe.stl");
  IndicatorLayer3D<T> extendedDomain(pipeI, physDeltaX);

  geometry.rename(0,2, extendedDomain);
  geometry.rename(2,1, pipeI);
  geometry.clean();

  Vector<T,3> origin = geometry.getStatistics().getMinPhysR( 2 );
  origin[1] += physDeltaX/2.;
  origin[2] += physDeltaX/2.;

  Vector<T,3> extend = geometry.getStatistics().getMaxPhysR( 2 );
  extend[1] = extend[1]-origin[1]-physDeltaX/2.;
  extend[2] = extend[2]-origin[2]-physDeltaX/2.;

  // Set material number for inflow
  origin[0] = geometry.getStatistics().getMinPhysR( 2 )[0]-physDeltaX;
  extend[0] = 2*physDeltaX;
  IndicatorCuboid3D<T> inflow( extend,origin );
  geometry.rename( 2,3,inflow );

  // Set material number for outflow
  origin[0] = geometry.getStatistics().getMaxPhysR( 2 )[0]-physDeltaX;
  extend[0] = 2*physDeltaX;
  IndicatorCuboid3D<T> outflow(extend, origin);
  geometry.rename(2,4, outflow);

  // Rotate the valve and set material number for it
  STLreader<T> valveI("valve.stl", physDeltaX, 0.001);
  const auto rotationPoint = parameters.get<parameters::ROTATION_POINT>();
  const auto rotationAxis = parameters.get<parameters::ROTATION_AXIS>();
  const T angle = parameters.get<parameters::ANGLE>();
  T rotationAngle = std::numbers::pi_v<T> / T(180) *  angle;
  IndicatorRotate<T,3> valveRot(rotationPoint, rotationAxis, rotationAngle, valveI);
  geometry.rename(1,5, valveRot);

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  geometry.checkForErrors();
  geometry.getStatistics().print();
}

// Set up the lattice of the simulation
void prepareLattice(MyCase& myCase) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});

  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();
  const T physDeltaT = parameters.get<parameters::PHYS_DELTA_T>();
  const T physLength = parameters.get<parameters::PHYS_CHAR_LENGTH>();
  const T physCharVelocity = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T physCharViscosity = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T physCharDensity = parameters.get<parameters::PHYS_CHAR_DENSITY>();

  // Set up a unit converter with the characteristic physical units
  myCase.getLattice(NavierStokes{}).setUnitConverter<UnitConverter<T,DESCRIPTOR>>(
    physDeltaX,        // physDeltaX: spacing between two lattice cells in [m]
    physDeltaT,        // physDeltaT: time step in [s]
    physLength,        // charPhysLength: reference length of simulation geometry in [m]
    physCharVelocity,  // charPhysVelocity: highest expected velocity during simulation in [m/s]
    physCharViscosity, // physViscosity: physical kinematic viscosity in [m^2/s]
    physCharDensity    // physDensity: physical density [kg/m^3]
  );
  lattice.getUnitConverter().print();

  // Material=1 -->bulk dynamics
  lattice.defineDynamics<BGKdynamics>(geometry, 1);

  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(lattice, geometry, 2);

  // Material=3 -->fixed velocity
  boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 3);

  // Material=4 -->fixed pressure
  boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 4);

  // Material=5 -->bounce back
  boundary::set<boundary::BounceBack>(lattice, geometry, 5);

  lattice.setParameter<descriptors::OMEGA>(lattice.getUnitConverter().getLatticeRelaxationFrequency());
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  // Initial conditions
  AnalyticalConst3D<T,T> rhoF( 1 );
  AnalyticalConst3D<T,T> uF(T(0), T(0), T(0));

  // Initialize all values of distribution functions to their local equilibrium
  auto domain = myCase.getGeometry().getMaterialIndicator({1,2,3,4});
  lattice.defineRhoU( domain, rhoF, uF );
  lattice.iniEquilibrium( domain, rhoF, uF );

  // Make the lattice ready for simulation
  lattice.initialize();
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();

  // No of time steps for smooth start-up
  const T startTime = parameters.get<parameters::PHYS_START_T>();
  const T updateTime = parameters.get<parameters::PHYS_UPDATE_T>();
  int iTmaxStart = converter.getLatticeTime( startTime );
  int iTupdate = converter.getLatticeTime( updateTime );

  if ( int(iT)%iTupdate == 0 && int(iT) <= iTmaxStart ) {
    PolynomialStartScale<T,int> StartScale( iTmaxStart, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    int iTvec[1] = {int(iT)};
    T frac[1] = {};
    StartScale( frac,iTvec );
    std::vector<T> maxVelocity( 3,0 );
    maxVelocity[0] = T(2.25)*frac[0]*converter.getCharLatticeVelocity();
    T distance2Wall = converter.getPhysDeltaX()/T(2);
    CirclePoiseuille3D<T> poiseuilleU(myCase.getGeometry(), 3, maxVelocity[0], distance2Wall);
    lattice.defineU( myCase.getGeometry(), 3, poiseuilleU );

    // Update velocity on GPU
    lattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

// Computes the pressure drop between the voxels before and after the cylinder
void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  /// Write vtk plots every 0.3 seconds (of phys. simulation time)
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();

  const T physVtkIter = parameters.get<parameters::PHYS_VTK_ITER_T>();
  const T physStatIter = parameters.get<parameters::PHYS_STAT_ITER_T>();

  const std::size_t vtkIter  = converter.getLatticeTime( physVtkIter );
  const std::size_t statIter = converter.getLatticeTime( physStatIter );

  SuperVTMwriter3D<T> vtmWriter( "pipeWithValve3d" );
  SuperLatticePhysVelocity3D velocity( lattice, converter );
  SuperLatticePhysPressure3D pressure( lattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  if ( iT==0 ) {
    vtmWriter.createMasterFile();
  }

  // Writes the vtk files
  if (iT%vtkIter == 0) {
    // Send values from GPU to CPU for evaluation
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }

  // Writes output on the console
  if ( iT%statIter == 0 ) {
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
    timer.print(iT);
  }
}

void simulate(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  const T maxPhysT = parameters.get<parameters::PHYS_MAX_T>();

  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(maxPhysT);
  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    myCase.getLattice(NavierStokes{}).collideAndStream();

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, timer, iT);
  }

  timer.stop();
  timer.printSummary();
}

/// Setup and run a simulation
int main(int argc, char* argv[]) {
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;

    myCaseParameters.set<STL_PATH>("pipe.stl");
    myCaseParameters.set<STL_SCALING>(0.001);

    myCaseParameters.set<PHYS_DELTA_X       >(       0.00125);
    myCaseParameters.set<PHYS_DELTA_T       >(     0.0003125);
    myCaseParameters.set<PHYS_CHAR_LENGTH   >(           0.5);
    myCaseParameters.set<PHYS_CHAR_VELOCITY >(           0.2);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(        0.0005);
    myCaseParameters.set<PHYS_CHAR_DENSITY  >(        1000.0);
    myCaseParameters.set<PHYS_MAX_T         >(           16.);
    myCaseParameters.set<ROTATION_POINT     >({0.15, 0., 0.});
    myCaseParameters.set<ROTATION_AXIS      >(  {0., 0., 1.});
    myCaseParameters.set<ANGLE              >(           75.);
    myCaseParameters.set<PHYS_START_T       >(           6.4);
    myCaseParameters.set<PHYS_UPDATE_T      >(      0.009375);
    myCaseParameters.set<PHYS_VTK_ITER_T    >(           0.3);
    myCaseParameters.set<PHYS_STAT_ITER_T   >(           0.1);
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
