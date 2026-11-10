/*  Lattice Boltzmann sample, written in C++, using the OpenLB library
 *
 *  Copyright (C) 2006-2019 Jonas Latt, Mathias J. Krause,
 *  Vojtech Cvrcek, Peter Weisbrod, Adrian Kummerländer
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

/* cylinder2d.h:
 * This example examines a steady flow past a cylinder placed in a channel.
 * The cylinder is offset somewhat from the center of the flow to make the
 * steady-state symmetrical flow unstable. At the inlet, a Poiseuille profile is
 * imposed on the velocity, whereas the outlet implements a Dirichlet pressure
 * condition set by p = 0.
 * Inspired by "Benchmark Computations of Laminar Flow Around
 * a Cylinder" by M.Schäfer and S.Turek. For high resolution, low
 * latticeU, and enough time to converge, the results for pressure drop, drag
 * and lift lie within the estimated intervals for the exact results.
 * An unsteady flow with Karman vortex street can be created by changing the
 * Reynolds number to Re=100.
 */
#pragma once

#include <olb.h>

using namespace olb;
using namespace olb::graphics;
using namespace olb::names;

#define BOUZIDI

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<FLOATING_POINT_TYPE, descriptors::D2Q9<>>
>;

namespace olb::parameters {

struct RADIUS_CYLINDER            : public descriptors::FIELD_BASE<1> { };
struct CENTER_CYLINDER            : public descriptors::FIELD_BASE<0, 1> { };
struct RAMP_UP_END_FRACTION       : public descriptors::FIELD_BASE<1> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const T physLengthX = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T physDeltaX  = 2*parameters.get<parameters::RADIUS_CYLINDER>()/parameters.get<parameters::RESOLUTION>();
  const T physLengthY = parameters.get<parameters::DOMAIN_EXTENT>()[1] + physDeltaX;

  const Vector extent{physLengthX, physLengthY};
  const Vector origin{0, 0};

  IndicatorCuboid2D<T> cuboid(extent, origin);
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

// Stores geometry information in form of material numbers
void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry  = myCase.getGeometry();

  const T physLengthX = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T physDeltaX  = 2*parameters.get<parameters::RADIUS_CYLINDER>()
                      /   parameters.get<parameters::RESOLUTION>();
  const T physLengthY = parameters.get<parameters::DOMAIN_EXTENT>()[1] + physDeltaX;

  geometry.rename(0, 2);
  geometry.rename(2, 1, {1, 1});

  Vector<T, 2> extent(physLengthX, physLengthY);
  Vector<T, 2> origin;

  // Set material number for inflow
  extent[0] = 2.*physDeltaX;
  origin[0] = -physDeltaX;
  IndicatorCuboid2D<T> inflow(extent, origin );
  geometry.rename(2, 3, 1, inflow);

  // Set material number for outflow
  origin[0] = physLengthX-physDeltaX;
  IndicatorCuboid2D<T> outflow(extent, origin);
  geometry.rename(2, 4, 1, outflow);

  // Set material number for cylinder
  const T radiusCylinder  = parameters.get<parameters::RADIUS_CYLINDER>();

  Vector<T,2> center(parameters.get<parameters::CENTER_CYLINDER>()[0], parameters.get<parameters::CENTER_CYLINDER>()[1] + physDeltaX/2.);
  std::shared_ptr<IndicatorF2D<T>> circle = std::make_shared<IndicatorCircle2D<T>>(center, radiusCylinder);
  geometry.rename(1, 5, circle);

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
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
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes{});
  auto& geometry   = myCase.getGeometry();

  const T physDeltaX            = parameters.get<parameters::PHYS_CHAR_LENGTH>() / parameters.get<parameters::RESOLUTION>();
  lattice.setUnitConverter(parameters);
  lattice.getUnitConverter().print();

  // Material=1 -->bulk dynamics
  auto bulkIndicator = geometry.getMaterialIndicator({1});
  dynamics::set<BGKdynamics>(lattice, bulkIndicator);

  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(lattice, geometry, 2);

  //if boundary conditions are chosen to be interpolatedy, 3);
  boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 3);
  boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 4);

  // Material=5 -->bouzidi / bounce back
  const T centerCylinderX = parameters.get<parameters::CENTER_CYLINDER>()[0];
  const T centerCylinderY = parameters.get<parameters::CENTER_CYLINDER>()[1] + physDeltaX/2.;
  const T radiusCylinder  = parameters.get<parameters::RADIUS_CYLINDER>();

  Vector<T,2> center(centerCylinderX, centerCylinderY);
  std::shared_ptr<IndicatorF2D<T>> circle = std::make_shared<IndicatorCircle2D<T>>(center, radiusCylinder);

  #ifdef BOUZIDI
  setBouzidiBoundary(lattice, geometry, 5, *circle);
  #else
  boundary::set<boundary::BounceBack>(lattice, geometry, 5);
  #endif

  lattice.setParameter<descriptors::OMEGA>(lattice.getUnitConverter().getLatticeRelaxationFrequency());
  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  myCase.getLattice(NavierStokes()).initialize();
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setTemporalValues(MyCase& myCase,
                      std::size_t iT)
{
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes{});
  auto& geometry   = myCase.getGeometry();

  const T maxPhysT   = parameters.get<parameters::PHYS_MAX_T>();
  const T physDeltaX = 2*parameters.get<parameters::RADIUS_CYLINDER>()
                     /   parameters.get<parameters::RESOLUTION>();

  // No of time steps for smooth start-up
  std::size_t iTmaxStart = lattice.getUnitConverter().getLatticeTime(
    maxPhysT*parameters.get<parameters::RAMP_UP_END_FRACTION>());
  std::size_t iTupdate   = lattice.getUnitConverter().getLatticeTime(
    parameters.get<parameters::RAMP_UP_T>());

  if (iT % iTupdate == 0 && iT <= iTmaxStart) {
    T frac{};
    PolynomialStartScale<T,std::size_t>(iTmaxStart, T(1))(&frac, &iT);
    T maxVelocity   = 1.5*lattice.getUnitConverter().getCharPhysVelocity()*frac;
    T distance2Wall = 0.5*physDeltaX;
    Poiseuille2D<T> poiseuilleU(geometry, 3, maxVelocity, distance2Wall);
    momenta::setVelocity(lattice, geometry.getMaterialIndicator(3), poiseuilleU);
    lattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

// Computes the pressure drop between the voxels before and after the cylinder
void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  OstreamManager clout( std::cout,"getResults" );
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes{});
  auto& geometry   = myCase.getGeometry();

  const UnitConverter<T,DESCRIPTOR>& converter = lattice.getUnitConverter();

  // Gnuplot constructor (must be static!)
  // for real-time plotting: gplot("name", true) // experimental!
  static Gnuplot<T> gplot("drag");

  SuperVTMwriter2D<T> vtmWriter("cylinder2d");
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(lattice, converter);
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure(lattice, converter);

  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(pressure);

  const std::size_t vtkIter  = converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T>());
  const std::size_t statIter = converter.getLatticeTime(parameters.get<parameters::PHYS_STAT_ITER_T>());

  const T maxPhysT        = parameters.get<parameters::PHYS_MAX_T>();
  const T physDeltaX      = 2.*parameters.get<parameters::RADIUS_CYLINDER>()/parameters.get<parameters::RESOLUTION>();
  const T centerCylinderX = parameters.get<parameters::CENTER_CYLINDER>()[0];
  const T centerCylinderY = parameters.get<parameters::CENTER_CYLINDER>()[1] + physDeltaX/2.;
  const T radiusCylinder  = parameters.get<parameters::RADIUS_CYLINDER>();

  T point[2] = {};
  point[0] = centerCylinderX + 3*radiusCylinder;
  point[1] = centerCylinderY;
  AnalyticalFfromSuperF2D<T> intpolateP(pressure, true);
  T p;
  intpolateP(&p, point);

  if (iT == 0) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid(lattice);
    SuperLatticeRank2D<T, DESCRIPTOR> rank(lattice);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);

    vtmWriter.createMasterFile();
  }

  if (iT%statIter == 0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    // Timer console output
    timer.update(iT);
    timer.printStep();

    // Lattice statistics console output
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));

    // Drag, lift, pressure drop
    AnalyticalFfromSuperF2D<T> intpolatePressure(pressure, true);
    SuperLatticePhysDrag2D<T,DESCRIPTOR> drag(lattice, geometry, 5, converter );

    T point1[2] = {};
    T point2[2] = {};

    point1[0] = centerCylinderX - radiusCylinder;
    point1[1] = centerCylinderY;

    point2[0] = centerCylinderX + radiusCylinder;
    point2[1] = centerCylinderY;

    T p1, p2;
    intpolatePressure(&p1, point1);
    intpolatePressure(&p2, point2);

    clout << "pressure1=" << p1;
    clout << "; pressure2=" << p2;

    T pressureDrop = p1-p2;
    clout << "; pressureDrop=" << pressureDrop;

    int input[3] = {};
    T _drag[drag.getTargetDim()];
    drag(_drag, input);
    clout << "; drag=" << _drag[0] << "; lift=" << _drag[1] << std::endl;

    myCase.getParameters().set<parameters::COEFFICIENT_OF_DRAG>(_drag[0]);

    // set data for gnuplot: input={xValue, yValue(s), names (optional), position of key (optional)}
    gplot.setData(converter.getPhysTime(iT), {_drag[0], 5.58}, {"drag(openLB)", "drag(schaeferTurek)"}, "bottom right", {'l','l'});

    // every (iT%vtkIter) write an png of the plot
    if (iT % (vtkIter) == 0) {
      // writes pngs: input={name of the files (optional), x range for the plot (optional)}
      gplot.writePNG(iT, maxPhysT);
    }
  }

  // Writes the vtk files
  if (iT % vtkIter == 0 && iT > 0) {
    vtmWriter.write(iT);

    {
      SuperEuklidNorm2D<T, DESCRIPTOR> normVel(velocity);
      BlockReduction2D2D<T> planeReduction(normVel, 600, BlockDataSyncMode::ReduceOnly);
      // write output as JPEG
      heatmap::write(planeReduction, iT);
    }
  }

  // write pdf at last time step
  if (iT == converter.getLatticeTime(maxPhysT) - 1) {
    // writes pdf
    gplot.writePDF();
  }
}

void simulate(MyCase& myCase){
  OstreamManager clout(std::cout, "simulate");
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  const T physMaxT = parameters.get<parameters::PHYS_MAX_T>();
  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(physMaxT);

  clout << "starting simulation..." << std::endl;
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
