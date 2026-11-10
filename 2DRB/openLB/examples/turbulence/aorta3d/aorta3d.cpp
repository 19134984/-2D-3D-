/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2014 Mathias J. Krause
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

/* aorta3d.cpp:
 * In this example the fluid flow through a bifurcation is
 * simulated. The geometry is obtained from a mesh in stl-format.
 * With Bouzidi boundary conditions the curved boundary is
 * adequately mapped and initialized fully automatically. As
 * dynamics a Smagorinsky turbulent BGK model is used to stabilize
 * the simulation for low resolutions. As output the flux at the
 * inflow and outflow region is computed. The wall stress can be
 * visualized on the stl Mesh with the Mesh.pvd file in paraview.
 * The results has been validated by comparison with other results
 * obtained with FEM and FVM.
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::descriptors;
using namespace olb::graphics;

using MyCase = Case<
  NavierStokes, Lattice<double, D3Q19<>>
>;

namespace olb::parameters {

struct BOUZIDI_ENABLED : public descriptors::TYPED_FIELD_BASE<bool,1> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  return Mesh<MyCase::value_t,MyCase::d>::fromSTL(parameters);
}

// Stores data from stl file in geometry in form of material numbers
void prepareGeometry(MyCase& myCase)
{
  using T = MyCase::value_t;

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  const T dx = parameters.get<parameters::PHYS_CHAR_LENGTH>() / parameters.get<parameters::RESOLUTION>();
  auto stlI = myCase.getMesh().getIndicator(parameters.get<parameters::STL_PATH>());
  IndicatorLayer3D<T> indicator(*stlI, dx);

  geometry.rename(0,2, indicator);
  geometry.rename(2,1, *stlI);

  geometry.clean();

  // Set material number for inflow
  IndicatorCircle3D<T> inflow(  0.218125,0.249987,0.0234818, 0., 1.,0., 0.0112342 );
  IndicatorCylinder3D<T> layerInflow( inflow, 2.* dx);
  geometry.rename( 2,3,1,layerInflow );

  // Set material number for outflow0
  //IndicatorCircle3D<T> outflow0(0.2053696,0.0900099,0.0346537,  2.5522,5.0294,-1.5237, 0.0054686 );
  IndicatorCircle3D<T> outflow0( 0.2053696,0.0900099,0.0346537, 0.,-1.,0., 0.0054686 );
  IndicatorCylinder3D<T> layerOutflow0( outflow0, 2.*dx );
  geometry.rename( 2,4,1,layerOutflow0 );

  // Set material number for outflow1
  //IndicatorCircle3D<T> outflow1(0.2388403,0.0900099,0.0343228, -1.5129,5.1039,-2.8431, 0.0058006 );
  IndicatorCircle3D<T> outflow1( 0.2388403,0.0900099,0.0343228, 0.,-1.,0., 0.0058006 );
  IndicatorCylinder3D<T> layerOutflow1( outflow1, 2.*dx );
  geometry.rename( 2,5,1,layerOutflow1 );

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  // Removes all not needed boundary voxels inside the surface
  geometry.innerClean( 3 );
  geometry.checkForErrors();

  geometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(MyCase& myCase)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using BulkDynamics = SmagorinskyBGKdynamics<T,DESCRIPTOR>;

  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  auto& lattice = myCase.getLattice(NavierStokes{});

  const T physDeltaX = parameters.get<parameters::PHYS_CHAR_LENGTH>() / parameters.get<parameters::RESOLUTION>();
  const T physDeltaT = parameters.get<parameters::PHYS_CHAR_LENGTH>() / (parameters.get<parameters::TIME_RESOLUTION>() * parameters.get<parameters::RESOLUTION>());
  const T physLength = parameters.get<parameters::PHYS_CHAR_LENGTH>();
  const T physCharVelocity = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T physCharViscosity = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T physCharDensity = parameters.get<parameters::PHYS_CHAR_DENSITY>();
  const bool bouzidiOn = parameters.get<parameters::BOUZIDI_ENABLED>();

  myCase.getLattice(NavierStokes{}).setUnitConverter(
    physDeltaX,        // physDeltaX: spacing between two lattice cells in [m]
    physDeltaT,        // physDeltaT: time step in [s]
    physLength,        // charPhysLength: reference length of simulation geometry in [m]
    physCharVelocity,  // charPhysVelocity: highest expected velocity during simulation in [m/s]
    physCharViscosity, // physViscosity: physical kinematic viscosity in [m^2/s]
    physCharDensity    // physDensity: physical density [kg/m^3]
  );
  lattice.getUnitConverter().print();
  lattice.getUnitConverter().write("aorta3d");

  util::Timer<T> timer1( lattice.getUnitConverter().getLatticeTime(parameters.get<parameters::PHYS_MAX_T>()), geometry.getStatistics().getNvoxel() );
  timer1.start();

  auto stlI = myCase.getMesh().getIndicator(parameters.get<parameters::STL_PATH>());

  // material=1 --> bulk dynamics
  dynamics::set<BulkDynamics>(lattice, geometry.getMaterialIndicator({1}));

  if ( bouzidiOn ) {
    // material=2 --> no dynamics + bouzidi zero velocity
    setBouzidiBoundary<T,DESCRIPTOR>(lattice, geometry, 2, *stlI);
    // material=3 --> no dynamics + bouzidi velocity (inflow)
    setBouzidiBoundary<T,DESCRIPTOR,BouzidiVelocityPostProcessor>(lattice, geometry, 3, *stlI);
  }
  else {
    // material=2 --> bounceBack dynamics
    boundary::set<boundary::BounceBack>(lattice, geometry, 2);
    // material=3 --> bulk dynamics + velocity (inflow)
    dynamics::set<BulkDynamics>(lattice, geometry.getMaterialIndicator(3));
    boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 3);
  }

  // material=4,5 --> bulk dynamics + pressure (outflow)
  dynamics::set<BulkDynamics>(lattice, geometry.getMaterialIndicator({4, 5}));
  boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 4);
  boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 5);

  clout << "Prepare Lattice ... OK" << std::endl;

  timer1.stop();
  timer1.printSummary();
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;

  auto& lattice = myCase.getLattice(NavierStokes{});

  const T omega = lattice.getUnitConverter().getLatticeRelaxationFrequency();

  lattice.setParameter<descriptors::OMEGA>(omega);
  lattice.setParameter<collision::LES::SMAGORINSKY>(T(0.1));
  // Lattice initialize
  lattice.initialize();
}

// Generates a slowly increasing sinuidal inflow
void setBoundaryValues(MyCase& myCase, std::size_t iT)
{
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();

  // No of time steps for smooth start-up
  std::size_t iTperiod = lattice.getUnitConverter().getLatticeTime( 0.5 );
  std::size_t iTupdate = 50;

  if ( iT%iTupdate == 0 ) {
    // Smooth start curve, sinus
    SinusStartScale<T,std::size_t> nSinusStartScale( iTperiod, lattice.getUnitConverter().getCharPhysVelocity() );

    // Creates and sets the Poiseuille inflow profile using functors
    std::size_t iTvec[1]= {iT};
    T maxVelocity[1]= {T()};
    nSinusStartScale( maxVelocity,iTvec );

    if (parameters.get<parameters::BOUZIDI_ENABLED>()) {
      // setBouzidiVelocity expects the velocity in lattice units
      CirclePoiseuille3D<T> velocity( geometry,3,lattice.getUnitConverter().getLatticeVelocity(maxVelocity[0]), T() );
      setBouzidiVelocity(lattice, geometry, 3, velocity);
      lattice.setProcessingContext<Array<descriptors::BOUZIDI_VELOCITY>>(
        ProcessingContext::Simulation);
    }
    else {
      CirclePoiseuille3D<T> velocity( geometry,3,maxVelocity[0], T() );
      momenta::setVelocity(lattice, geometry.getMaterialIndicator(3), velocity);
      lattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
        ProcessingContext::Simulation);
    }
  }
}

// Computes flux at inflow and outflow
void getResults(MyCase& myCase, SuperVtuSurfaceWriter<MyCase::value_t>& vtuWriter, util::Timer<MyCase::value_t>& timer, std::size_t iT)
{
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& converter = lattice.getUnitConverter();
  auto& parameters = myCase.getParameters();
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

  const T physDeltaX = converter.getPhysDeltaX();
  const bool bouzidiOn = parameters.get<parameters::BOUZIDI_ENABLED>();

  auto stlI = myCase.getMesh().getSTL(parameters.get<parameters::STL_PATH>());

  OstreamManager clout( std::cout,"getResults" );

  const std::size_t vtkIter  = converter.getLatticeTime( .1 );
  const std::size_t statIter = converter.getLatticeTime( .1 );

  if ( iT==0 ) {
    SuperVTMwriter3D<T> vtmWriter("aorta3d");
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization

    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( lattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( lattice );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
    vtuWriter.createMasterFile();
  }

  // Writes the vtk files
  if ( iT%vtkIter==0 ) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    lattice.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriter("aorta3d");
      SuperLatticePhysVelocity3D velocity(lattice, converter);
      SuperLatticePhysPressure3D pressure(lattice, converter);
      vtmWriter.addFunctor(velocity);
      vtmWriter.addFunctor(pressure);
      task(vtmWriter, iT);
    });

    // Write interpolated wall shear stress on the STL surface
    {
      SuperLatticeDensity3D densityF(lattice);
      AnalyticalFfromSuperF3D smoothDensityF(densityF);

      SuperLatticeStress3D stressF(lattice);
      AnalyticalFfromSuperF3D smoothStressF(stressF);

      PhysWallShearStressOnSurface3D<T,DESCRIPTOR> interpolatedWssF(converter, smoothDensityF, smoothStressF, *stlI);
      interpolatedWssF.getName() = "interpolatedWss";
      vtuWriter.addFunctor(interpolatedWssF);

      vtuWriter.write(iT);
    }
  }

  // Writes output on the console
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    lattice.getStatistics().print( iT, converter.getPhysTime( iT ) );

    // Flux at the inflow and outflow region
    std::vector<int> materials = { 1, 3, 4, 5 };

    IndicatorCircle3D<T> inflow(  0.218125,0.249987-2.*physDeltaX,0.0234818, 0., -1.,0., 0.0112342+2*physDeltaX );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxInflow( lattice, converter, geometry, inflow, materials, BlockDataReductionMode::Discrete );
    vFluxInflow.print( "inflow","ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxInflow( lattice, converter, geometry, inflow, materials, BlockDataReductionMode::Discrete );
    pFluxInflow.print( "inflow","N","mmHg" );

    IndicatorCircle3D<T> outflow0( 0.2053696,0.0900099+2.*physDeltaX,0.0346537, 0.,1.,0., 0.0054686+2*physDeltaX );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOutflow0( lattice, converter, geometry, outflow0, materials, BlockDataReductionMode::Discrete );
    vFluxOutflow0.print( "outflow0","ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxOutflow0( lattice, converter, geometry, outflow0, materials, BlockDataReductionMode::Discrete );
    pFluxOutflow0.print( "outflow0","N","mmHg" );

    IndicatorCircle3D<T> outflow1( 0.2388403,0.0900099+2.*physDeltaX,0.0343228, 0.,1.,0., 0.0058006+2*physDeltaX );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOutflow1( lattice, converter, geometry, outflow1, materials, BlockDataReductionMode::Discrete );
    vFluxOutflow1.print( "outflow1","ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxOutflow1( lattice, converter, geometry, outflow1, materials, BlockDataReductionMode::Discrete );
    pFluxOutflow1.print( "outflow1","N","mmHg" );

    if ( bouzidiOn ) {
      SuperLatticeYplus3D<T, DESCRIPTOR> yPlus( lattice, converter, geometry, *stlI, 3 );
      SuperMax3D<T> yPlusMaxF( yPlus, geometry, 1 );
      int input[4]= {};
      T yPlusMax[1];
      yPlusMaxF( yPlusMax,input );
      clout << "yPlusMax=" << yPlusMax[0] << std::endl;
    }
  }

  if ( lattice.getStatistics().getMaxU() > 0.3 ) {
    clout << "PROBLEM uMax=" << lattice.getStatistics().getMaxU() << std::endl;
    std::exit(0);
  }
}

void simulate(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();

  const T maxPhysT = parameters.get<parameters::PHYS_MAX_T>();

  OstreamManager clout(std::cout, "simulate");

  auto stlI = myCase.getMesh().getSTL(parameters.get<parameters::STL_PATH>());
  SuperVtuSurfaceWriter<T> vtuWriter("surface", myCase.getMesh().getCuboidDecomposition(), myCase.getMesh().getLoadBalancer(), *stlI);

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), geometry.getStatistics().getNvoxel() );
  timer.start();

  for ( std::size_t iT = 0; iT <= converter.getLatticeTime( maxPhysT ); iT++ ) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(myCase, iT);

    // === 6th Step: Collide and Stream Execution ===
    lattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults(myCase, vtuWriter, timer, iT);
  }

  timer.stop();
  timer.printSummary();
}

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<STL_PATH>("aorta3d.stl");
    myCaseParameters.set<STL_SCALING>(0.001);
    myCaseParameters.set<STL_RAY_MODE>(RayMode::FastRayZ);
    myCaseParameters.set<DECOMPOSITION_STRATEGY>("volume");
    myCaseParameters.set<DECOMPOSITION_MULTIPLIER>(3);

    myCaseParameters.set<RESOLUTION         >(40);
    myCaseParameters.set<TIME_RESOLUTION    >(20);
    myCaseParameters.set<PHYS_CHAR_LENGTH   >(0.02246); // reference length of simulation geometry in m
    myCaseParameters.set<PHYS_CHAR_VELOCITY >(0.45);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(0.003/1055.);
    myCaseParameters.set<PHYS_CHAR_DENSITY  >(1055);
    myCaseParameters.set<PHYS_MAX_T         >(2.); // max. simulation time in s, SI unit
    myCaseParameters.set<BOUZIDI_ENABLED    >(true);
    myCaseParameters.set<PHYS_DELTA_X       >([] (MyCase::ParametersD& p) {
      return p.get<parameters::PHYS_CHAR_LENGTH>()
           / p.get<parameters::RESOLUTION>();
    });
  }

  Mesh mesh = createMesh(myCaseParameters);

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
