/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2013 Mathias J. Krause, Thomas Henn, Tim Dornieden
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

#pragma once

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::names;

#define BOUZIDI

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<float, descriptors::D3Q19<>>
>;

namespace olb::parameters {

struct RADIUS_CYLINDER            : public descriptors::FIELD_BASE<1> { };
struct RAMP_UP_END_FRACTION       : public descriptors::FIELD_BASE<1> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const T physDeltaX  = 2.*parameters.get<parameters::RADIUS_CYLINDER>()/parameters.get<parameters::RESOLUTION>();

  STLreader<T> stlReader("../../laminar/cylinder3d/cylinder3d.stl", physDeltaX, 0.001);
  IndicatorLayer3D<T> extendedDomain(stlReader, physDeltaX);

  #ifdef PARALLEL_MODE_MPI
    const int noOfCuboids = singleton::mpi().getSize();
  #else
    const int noOfCuboids = 7;
  #endif

  Mesh<T,MyCase::d> mesh(extendedDomain, physDeltaX, noOfCuboids);
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

// Stores data from stl file in geometry in form of material numbers
void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry  = myCase.getGeometry();

  const T physDeltaX  = 2.*parameters.get<parameters::RADIUS_CYLINDER>()/parameters.get<parameters::RESOLUTION>();
  STLreader<T> stlReader("../../laminar/cylinder3d/cylinder3d.stl", physDeltaX, 0.001);
  IndicatorLayer3D<T> extendedDomain(stlReader, physDeltaX);

  geometry.rename(0, 2, extendedDomain);
  geometry.rename(2, 1, stlReader);
  geometry.clean();

  Vector<T,3> origin = geometry.getStatistics().getMinPhysR(2);
  origin[1] += physDeltaX/2.;
  origin[2] += physDeltaX/2.;

  Vector<T,3> extend = geometry.getStatistics().getMaxPhysR(2);
  extend[1] = extend[1]-origin[1]-physDeltaX/2.;
  extend[2] = extend[2]-origin[2]-physDeltaX/2.;

  // Set material number for inflow
  origin[0] = geometry.getStatistics().getMinPhysR(2)[0]-physDeltaX;
  extend[0] = 2*physDeltaX;
  IndicatorCuboid3D<T> inflow(extend, origin);
  geometry.rename(2, 3, inflow);

  // Set material number for outflow
  origin[0] = geometry.getStatistics().getMaxPhysR(2)[0] - physDeltaX;
  extend[0] = 2*physDeltaX;
  IndicatorCuboid3D<T> outflow(extend, origin);
  geometry.rename(2, 4, outflow);

  // Set material number for cylinder
  origin[0] = geometry.getStatistics().getMinPhysR(2)[0] + physDeltaX;
  extend[0] = (geometry.getStatistics().getMaxPhysR(2)[0] - geometry.getStatistics().getMinPhysR(2)[0])/2.;
  std::shared_ptr<IndicatorF3D<T>> cylinder = std::make_shared<IndicatorCuboid3D<T>>(extend, origin);
  geometry.rename(2, 5, cylinder);

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  geometry.checkForErrors();
  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes{});
  auto& geometry   = myCase.getGeometry();

  const int resolution          = parameters.get<parameters::RESOLUTION>();
  const T physDeltaX            = 2.*parameters.get<parameters::RADIUS_CYLINDER>()/resolution;
  const T latticeRelaxationTime = parameters.get<parameters::LATTICE_RELAXATION_TIME>();
  const T diameterCylinder      = 2.*parameters.get<parameters::RADIUS_CYLINDER>();
  const T charPhysVelocity      = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T Re                    = parameters.get<parameters::REYNOLDS>();
  const T physDensity           = parameters.get<parameters::PHYS_CHAR_DENSITY>();

  STLreader<T> stlReader("../../laminar/cylinder3d/cylinder3d.stl", physDeltaX, 0.001);

  myCase.getLattice(NavierStokes{}).setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    resolution,                           // resolution: number of voxels per charPhysL
    latticeRelaxationTime,                // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    diameterCylinder,                     // charPhysLength: reference length of simulation geometry
    charPhysVelocity,                     // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    charPhysVelocity*diameterCylinder/Re, // physViscosity: physical kinematic viscosity in __m^2 / s__
    physDensity                           // physDensity: physical density in __kg / m^3__
  );
  lattice.getUnitConverter().print();

  // Material=1 -->bulk dynamics
  auto bulkIndicator = geometry.getMaterialIndicator({1});
  dynamics::set<BGKdynamics>(lattice, bulkIndicator);

  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(lattice, geometry, 2);

  // Setting of the boundary conditions

  // if local boundary conditions are chosen
  //boundary::set<boundary::LocalVelocity>(sLattice, geometry, 3);
  //boundary::set<boundary::LocalPressure>(sLattice, geometry, 4);

  //if interpolated boundary conditions are chosen
  boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 3);
  boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 4);

  // Material=5 -->bouzidi / bounce back
  #ifdef BOUZIDI
    setBouzidiBoundary<T,DESCRIPTOR>(lattice, geometry, 5, stlReader);
  #else
    boundary::set<boundary::BounceBack>(lattice, geometry, 5);
  #endif

  lattice.setParameter<descriptors::OMEGA>(lattice.getUnitConverter().getLatticeRelaxationFrequency());
  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  auto& lattice = myCase.getLattice(NavierStokes{});
  lattice.initialize();
}


// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setTemporalValues(MyCase& myCase,
                      std::size_t iT)
{
  OstreamManager clout(std::cout, "setTemporalValues");
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes{});
  auto& geometry   = myCase.getGeometry();

  const T maxPhysT   = parameters.get<parameters::PHYS_MAX_T>();
  const T physDeltaX = 2.*parameters.get<parameters::RADIUS_CYLINDER>()/parameters.get<parameters::RESOLUTION>();

  // No of time steps for smooth start-up
  std::size_t iTmaxStart = lattice.getUnitConverter().getLatticeTime(maxPhysT*parameters.get<parameters::RAMP_UP_END_FRACTION>());
  std::size_t iTupdate   = lattice.getUnitConverter().getLatticeTime(parameters.get<parameters::RAMP_UP_T>());

  if (iT%iTupdate == 0 && iT <= iTmaxStart) {
    T frac{};
    PolynomialStartScale<T,std::size_t>(iTmaxStart, T(1))(&frac, &iT);
    std::vector<T> maxVelocity(3, 0);
    maxVelocity[0] = 2.25*frac*lattice.getUnitConverter().getCharPhysVelocity();

    T distance2Wall = physDeltaX/2.;
    RectanglePoiseuille3D<T> poiseuilleU(geometry, 3, maxVelocity, distance2Wall, distance2Wall, distance2Wall);
    momenta::setVelocity(lattice, geometry.getMaterialIndicator(3), poiseuilleU);

    clout << "step=" << iT << "; maxVel=" << maxVelocity[0] << std::endl;

    lattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

// Computes the pressure drop between the voxels before and after the cylinder
void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  OstreamManager clout(std::cout, "getResults");
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes{});
  auto& geometry   = myCase.getGeometry();

  const auto& converter = lattice.getUnitConverter();
  //STLreader<T> stlReader("../../laminar/cylinder3d/cylinder3d.stl", converter.getPhysDeltaX(), 0.001);

  static Gnuplot<T> gplot("drag");

  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(lattice, converter);
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(lattice, converter);
  //SuperLatticeYplus3D<T, DESCRIPTOR> yPlus(lattice, converter, geometry, stlReader, 5);
  SuperLatticeRefinementMetricKnudsen3D<T, DESCRIPTOR> quality(lattice, converter);
  SuperRoundingF3D<T, T> roundedQuality(quality, RoundingMode::NearestInteger);
  SuperDiscretizationF3D<T> discretization(roundedQuality, 0., 2.);

  const std::size_t vtkIter  = converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T>());
  const std::size_t statIter = converter.getLatticeTime(parameters.get<parameters::PHYS_STAT_ITER_T>());
  const T maxPhysT           = parameters.get<parameters::PHYS_MAX_T>();

  if (parameters.get<parameters::VTK_ENABLED>()) {
    SuperVTMwriter3D<T> vtmWriter("cylinder3d");
    vtmWriter.addFunctor(quality);
    vtmWriter.addFunctor(velocity);
    vtmWriter.addFunctor(pressure);
    //vtmWriter.addFunctor(yPlus);
    if (iT == 0) {
      // Writes the geometry, cuboid no. and rank no. as vti file for visualization
      SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(lattice);
      SuperLatticeRank3D<T, DESCRIPTOR> rank(lattice);
      vtmWriter.write(cuboid);
      vtmWriter.write(rank);

      vtmWriter.createMasterFile();
    }

    // Writes the vtk files
    if (iT%vtkIter == 0) {
      vtmWriter.write(iT);

      {
        SuperEuklidNorm3D<T> normVel(velocity);
        BlockReduction3D2D<T> planeReduction(normVel, Vector<T,3>{0, 0, 1});
        // write output as JPEG
        heatmap::write(planeReduction, iT);
      }

      {
        BlockReduction3D2D<T> planeReduction(discretization, Vector<T,3>{0, 0, 1});
        heatmap::plotParam<T> jpeg_scale;
        jpeg_scale.colour = "blackbody";
        jpeg_scale.name = "quality";
        heatmap::write(planeReduction, iT, jpeg_scale);
      }
    }
  }

  // Writes output on the console
  if (iT%statIter == 0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    // Timer console output
    timer.update(iT);
    timer.printStep();

    // Lattice statistics console output
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));

    // Drag, lift, pressure drop
    AnalyticalFfromSuperF3D<T> intpolatePressure(pressure, true);
    SuperLatticePhysDrag3D<T,DESCRIPTOR> drag(lattice, geometry, 5, converter);

    olb::Vector<T, 3> point1V = geometry.getStatistics().getCenterPhysR(5);
    olb::Vector<T, 3> point2V = geometry.getStatistics().getCenterPhysR(5);
    T point1[3] = {};
    T point2[3] = {};
    for (int i = 0; i<3; i++) {
      point1[i] = point1V[i];
      point2[i] = point2V[i];
    }
    point1[0] = geometry.getStatistics().getMinPhysR(5)[0] - converter.getPhysDeltaX();
    point2[0] = geometry.getStatistics().getMaxPhysR(5)[0] + converter.getPhysDeltaX();

    T p1, p2;
    intpolatePressure(&p1, point1);
    intpolatePressure(&p2, point2);

    clout << "pressure1=" << p1;
    clout << "; pressure2=" << p2;

    T pressureDrop = p1-p2;
    clout << "; pressureDrop=" << pressureDrop;

    T dragA[3];
    int input1[0];
    drag(dragA, input1);
    clout << "; drag=" << dragA[0] << "; lift=" << dragA[1] << std::endl;

    myCase.getParameters().set<parameters::COEFFICIENT_OF_DRAG>(dragA[0]);

    //int input[4] = {};
    //SuperMax3D<T> yPlusMaxF(yPlus, geometry, 1);
    //T yPlusMax[1];
    //yPlusMaxF(yPlusMax, input);
    //clout << "yPlusMax=" << yPlusMax[0] << std::endl;

    if (parameters.get<parameters::GNUPLOT_ENABLED>()) {
      // set data for gnuplot: input={xValue, yValue(s), names (optional), position of key (optional)}
      gplot.setData(converter.getPhysTime(iT), {dragA[0], 5.58}, {"drag(openLB)", "drag(schaeferTurek)"}, "bottom right", {'l','l'});

      // every (iT%vtkIter) write an png of the plot
      if (iT%(vtkIter) == 0) {
        // writes pngs: input={name of the files (optional), x range for the plot (optional)}
        gplot.writePNG(iT, maxPhysT);
      }
    }
  }
}

void simulate(MyCase& myCase){
  OstreamManager clout(std::cout, "simulate");
  using T = MyCase::value_t;
  auto& parameters        = myCase.getParameters();
  const T physMaxT        = parameters.get<parameters::PHYS_MAX_T>();
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
