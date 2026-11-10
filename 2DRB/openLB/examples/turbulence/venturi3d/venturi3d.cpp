/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2014 Mathias J. Krause, Thomas Henn,
 *  Cyril Masquelier
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

/* venturi3d.cpp:
 * This example examines a steady flow in a venturi tube. At the
 * main inlet, a Poiseuille profile is imposed as Dirichlet velocity
 * boundary condition, whereas at the outlet and the minor inlet
 * a Dirichlet pressure condition is set by p=0 (i.e. rho=1).
 *
 * The example shows the usage of the Indicator functors to
 * build up a geometry and explains how to set boundary conditions
 * automatically.
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::descriptors;
using namespace olb::graphics;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<FLOATING_POINT_TYPE, descriptors::D3Q19<>>
>;

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  std::string fName("venturi3d.xml");
  XMLreader config(fName);

  std::shared_ptr<IndicatorF3D<T> > inflow = createIndicatorCylinder3D<T>(config["Geometry"]["Inflow"]["IndicatorCylinder3D"], false);
  std::shared_ptr<IndicatorF3D<T> > outflow0 = createIndicatorCylinder3D<T>(config["Geometry"]["Outflow0"]["IndicatorCylinder3D"], false);
  std::shared_ptr<IndicatorF3D<T> > outflow1 = createIndicatorCylinder3D<T>(config["Geometry"]["Outflow1"]["IndicatorCylinder3D"], false);

  std::shared_ptr<IndicatorF3D<T> > venturi = createIndicatorF3D<T>(config["Geometry"]["Venturi"], false);

  int N = config["Application"]["Discretization"]["Resolution"].get<int>();

  Mesh<T,MyCase::d> mesh(*venturi, 1./N, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();

  std::string fName("venturi3d.xml");
  XMLreader config(fName);

  std::shared_ptr<IndicatorF3D<T> > inflow = createIndicatorCylinder3D<T>(config["Geometry"]["Inflow"]["IndicatorCylinder3D"], false);
  std::shared_ptr<IndicatorF3D<T> > outflow0 = createIndicatorCylinder3D<T>(config["Geometry"]["Outflow0"]["IndicatorCylinder3D"], false);
  std::shared_ptr<IndicatorF3D<T> > outflow1 = createIndicatorCylinder3D<T>(config["Geometry"]["Outflow1"]["IndicatorCylinder3D"], false);

  std::shared_ptr<IndicatorF3D<T> > venturi = createIndicatorF3D<T>(config["Geometry"]["Venturi"], false);

  // Set boundary voxels by rename material numbers
  geometry.rename( 0,2, venturi );
  geometry.rename( 2,1,{1,1,1} );
  geometry.rename( 2,3,1, inflow );
  geometry.rename( 2,4,1, outflow0 );
  geometry.rename( 2,5,1, outflow1 );

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  // Removes all not needed boundary voxels inside the surface
  geometry.innerClean();
  geometry.checkForErrors();

  geometry.getStatistics().print();
  geometry.communicate();
}

void prepareLattice(MyCase& myCase) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;

  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();

  std::string fName("venturi3d.xml");
  XMLreader config(fName);

  UnitConverter<T, DESCRIPTOR>* converter = createUnitConverter<T, DESCRIPTOR>(config);
  myCase.getLattice(NavierStokes{}).setUnitConverter(*converter);
  lattice.getUnitConverter().print();

  // Material=1 -->bulk dynamics
  dynamics::set<RLBdynamics>(lattice, geometry, 1);

  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(lattice, geometry, 2);

  // Setting of the boundary conditions
  boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 3);
  boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 4);
  boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 5);

  lattice.setParameter<descriptors::OMEGA>(lattice.getUnitConverter().getLatticeRelaxationFrequency());

  lattice.initialize();

}

void setInitialValues(MyCase& myCase) {
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{
  using T = MyCase::value_t;

  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();

  const T maxPhysT = parameters.get<parameters::PHYS_MAX_T>();
  int iTmaxStart = lattice.getUnitConverter().getLatticeTime( maxPhysT*0.8 );
  int iTperiod = 50;

  if (iT % iTperiod == 0 && iT <= iTmaxStart) {
    PolynomialStartScale<T,std::size_t> startScale( iTmaxStart, T( 1 ) );
    T frac = T();
    startScale(&frac, &iT);

    // Creates and sets the Poiseuille inflow profile using functors
    CirclePoiseuille3D<T> poiseuilleU( geometry, 3, frac*lattice.getUnitConverter().getCharPhysVelocity(), T(), lattice.getUnitConverter().getPhysDeltaX() );
    momenta::setVelocity(lattice, geometry.getMaterialIndicator(3), poiseuilleU);

    lattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  /// Write vtk plots every 0.3 seconds (of phys. simulation time)
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;

  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();

  const std::size_t iTlog = converter.getLatticeTime(1.);
  const std::size_t iTvtk = converter.getLatticeTime(1.);

  SuperVTMwriter3D<T> vtmWriter( "venturi3d" );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( lattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( lattice );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Writes the vtm files
  if ( iT % iTvtk == 0 ) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    // Create the data-reading functors...
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( lattice, converter );
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( lattice, converter );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    vtmWriter.write( iT );

    SuperEuklidNorm3D<T> normVel( velocity );
    BlockReduction3D2D<T> planeReduction( normVel, {0, 0, 1} );

    // write output as JPEG
    heatmap::write(planeReduction, iT);

    // write output as JPEG and changing properties
    heatmap::plotParam<T> jpeg_Param;
    jpeg_Param.name             = "outflow";
    jpeg_Param.contourlevel     = 5;
    jpeg_Param.colour           = "blackbody";
    jpeg_Param.zoomOrigin       = {0.6, 0.3};
    jpeg_Param.zoomExtend       = {0.4, 0.7};
    heatmap::write(planeReduction, iT, jpeg_Param);
  }

  // Writes output on the console
  if ( iT % iTlog == 0 ) {
    timer.update( iT );
    timer.printStep();
    lattice.getStatistics().print( iT, lattice.getUnitConverter().getPhysTime( iT ) );
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


int main( int argc, char* argv[] )
{
  initialize( &argc, &argv );

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<PHYS_MAX_T>(200.);
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
