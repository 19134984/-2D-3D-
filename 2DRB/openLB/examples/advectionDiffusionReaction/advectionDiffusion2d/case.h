
/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2020 Louis Kronberg, Stephan Simonis
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

/// uncomment for console output of other norms
//#define allNorms

using namespace olb;
using namespace olb::graphics;
using namespace olb::names;

using MyCase = Case<
  AdvectionDiffusion, Lattice<double, descriptors::D2Q5<descriptors::VELOCITY>>
>;

namespace olb::parameters {
  struct PULSE_DIFF_BOUND : public descriptors::FIELD_BASE<1> {};
}

template <typename T>
class AdePhysTemp2D : public AnalyticalF2D<T,T> {

protected:
  MyCase& myCase;
  T t;
public:
  AdePhysTemp2D(T time, MyCase& _myCase) : AnalyticalF2D<T,T>(1), myCase(_myCase),
    t(time) {};

  bool operator()(T output[], const T input[]) override
  {
    auto& lattice = myCase.getLattice(AdvectionDiffusion{});
    const auto& converter = lattice.getUnitConverter();
    const T mue = converter.getPhysDiffusivity();
    const T uMag = converter.getCharPhysVelocity();
    const T res = converter.getResolution();

    const T x = input[0];
    const T y = input[1];
    const T gf = res/(res+1.);
    const T uMagX = uMag;
    const T uMagY = uMag;

    // initial condition (2D)
    output[0] = util::sin(M_PI*(x-uMagX*t)*gf);
    output[0] *= util::sin(M_PI*(y-uMagY*t)*gf);
    output[0] *= util::exp(-2*mue*util::pow(M_PI,2)*t*util::pow(gf,2));
    output[0] += 1.;

    return true;
  };
};

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  const T physCharLength = params.get<parameters::PHYS_CHAR_LENGTH>();
  Vector extent(physCharLength, physCharLength);
  Vector origin(-physCharLength/2, -physCharLength/2);

  IndicatorCuboid2D<T> cuboid(extent, origin);

#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif

  const T physDeltaX = physCharLength / params.get<parameters::RESOLUTION>();
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, noOfCuboids);
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({true,true});
  return mesh;
}

void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();

  geometry.rename(0,1);
  geometry.communicate();
  geometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  geometry.innerClean();
  geometry.checkForErrors();
  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase)
{
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<AdvectionDiffusion>;
  using namespace olb::parameters;

  auto& lattice = myCase.getLattice(AdvectionDiffusion{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const T physDeltaX = parameters.get<PHYS_CHAR_LENGTH>() / parameters.get<RESOLUTION>();
  const T physDeltaT = physDeltaX * physDeltaX;
  const T physLength = parameters.get<PHYS_CHAR_LENGTH>();
  const T mue = parameters.get<PHYS_DIFFUSIVITY>();
  const T physVelocity = parameters.get<PECLET>() * mue / physLength;
  const T physDensity = parameters.get<PHYS_CHAR_DENSITY>();

  lattice.setUnitConverter<AdeUnitConverter<T,DESCRIPTOR>>(
    physDeltaX,    // physDeltaX
    physDeltaT,    // physDeltaT (diffusive scaling)
    physLength,    // charPhysLength
    physVelocity,  // charPhysVelocity from Peclet
    mue,           // physDiffusivity
    physDensity    // physDensity,
  );

  const auto& converter = lattice.getUnitConverter();
  converter.print();

  dynamics::set<AdvectionDiffusionBGKdynamics>( lattice, geometry.getMaterialIndicator({1}) );
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(AdvectionDiffusion{});
  auto& geometry = myCase.getGeometry();
  const auto& converter = lattice.getUnitConverter();

  AnalyticalConst2D<T,T> u0( converter.getCharLatticeVelocity(), converter.getCharLatticeVelocity() );

  AdePhysTemp2D<T> Tinit( 0.0, myCase);

  auto bulkIndicator = geometry.getMaterialIndicator({0,1});

  fields::set<descriptors::VELOCITY>( lattice, bulkIndicator, u0 );
  momenta::setDensity( lattice, bulkIndicator, Tinit );

  lattice.setParameter<descriptors::OMEGA>( converter.getLatticeAdeRelaxationFrequency() );

  /// Make the lattice ready for simulation
  lattice.initialize();
}

// compute absolute and relative errors with different norms
MyCase::value_t errorEval( MyCase& myCase, const int iT)
{
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(AdvectionDiffusion{});
  auto& geometry = myCase.getGeometry();
  const auto& converter = lattice.getUnitConverter();

  using TDESCRIPTOR = MyCase::descriptor_t_of<AdvectionDiffusion>;

  OstreamManager clout(std::cout,"error");
  CSV<T> csvWriter("errors");

  T result[2] = {T(), T()};
  int tmp[] = {int()};

  SuperLatticeDensity2D<T, TDESCRIPTOR> temperature(lattice);
  AdePhysTemp2D<T> temperatureSol( converter.getPhysTime(iT), myCase);

  auto indicatorF = geometry.getMaterialIndicator({1});

  T l2rel_error;

  SuperRelativeErrorL2Norm2D<T> relTemperatureErrorL2Norm(temperature, temperatureSol, indicatorF);
  SuperAbsoluteErrorL2Norm2D<T> absTemperatureErrorL2Norm(temperature, temperatureSol, indicatorF);

  /// Console output of other error norms
#ifdef allNorms
  SuperAbsoluteErrorL1Norm2D<T> absTemperatureErrorL1Norm(temperature, temperatureSol, indicatorF);
  absTemperatureErrorL1Norm(result, tmp);
  clout << "temperature-L1-absolute-error=" << result[0] << std::endl;
  csvWriter.writeDataFile(iT, result[0], "l1-abs-error");

  SuperRelativeErrorL1Norm2D<T> relTemperatureErrorL1Norm(temperature, temperatureSol, indicatorF);
  relTemperatureErrorL1Norm(result, tmp);
  clout << "temperature-L1-relative-error=" << result[0] << std::endl;
  csvWriter.writeDataFile(iT, result[0],  "l1-rel-error");

  SuperAbsoluteErrorLinfNorm2D<T> absTemperatureErrorLinfNorm(temperature, temperatureSol, indicatorF);
  absTemperatureErrorLinfNorm(result, tmp);
  clout << "temperature-Linf-absolute-error=" << result[0] << std::endl;
  csvWriter.writeDataFile(iT, result[0],  "linf-abs-error");

  SuperRelativeErrorLinfNorm2D<T> relTemperatureErrorLinfNorm(temperature, temperatureSol, indicatorF);
  relTemperatureErrorLinfNorm(result, tmp);
  clout << "temperature-Linf-relative-error=" << result[0] << std::endl;
  csvWriter.writeDataFile(iT, result[0], "linf-rel-error");
#endif

  absTemperatureErrorL2Norm(result, tmp);
  clout << "temperature-L2-absolute-error=" << result[0] << std::endl;
  csvWriter.writeDataFile(iT, result[0],  "l2-abs-error");

  relTemperatureErrorL2Norm(result, tmp);
  clout << "temperature-L2-relative-error=" << result[0] << std::endl;
  csvWriter.writeDataFile(iT, result[0], "l2-rel-error");

  l2rel_error = result[0];

  return l2rel_error;
}

// will always return ZERO for iT == 0
MyCase::value_t getResults( MyCase& myCase,
              int iT,
              int statIter,
              util::Timer<MyCase::value_t>& timer )
{

  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(AdvectionDiffusion{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  const auto& converter = lattice.getUnitConverter();

  using TDESCRIPTOR = MyCase::descriptor_t_of<AdvectionDiffusion>;

  OstreamManager clout(std::cout,"getResults");
  SuperVTMwriter2D<T> vtkWriter("advectionDiffusion2d");

  AdePhysTemp2D<T> temperatureSol( converter.getPhysTime(iT), myCase);

  const T maxN = parameters.get<parameters::RESOLUTION_MAX>();

  T avg = 0;

  if (iT == 0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D<T, TDESCRIPTOR> cuboid(lattice);
    SuperLatticeRank2D<T, TDESCRIPTOR> rank(lattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();

    SuperLatticeDensity2D<T, TDESCRIPTOR> temperature(lattice);
    SuperLatticeFfromAnalyticalF2D<T, TDESCRIPTOR> solution(temperatureSol, lattice);

    vtkWriter.addFunctor( temperature );
    vtkWriter.addFunctor( solution );
    vtkWriter.write( iT );

    /// Creates contour plots as jpeg files for highest resolution run
    if (converter.getResolution() == maxN) {
      heatmap::plotParam<T> plotParam_temp0;
      plotParam_temp0.minValue = 0;
      plotParam_temp0.maxValue = 2;
      plotParam_temp0.contourlevel = 10;
      BlockReduction2D2D<T> planeReduction0( temperature, converter.getResolution(), BlockDataSyncMode::ReduceOnly );
      heatmap::write( planeReduction0, iT, plotParam_temp0 );
    }

    timer.update(iT);
    timer.printStep();
    lattice.getStatistics().print( iT, converter.getPhysTime(iT) );
  }
  else if (iT % statIter == 0 ) {
    /// Writes the VTK files
    SuperLatticeDensity2D<T, TDESCRIPTOR> temperature(lattice);
    SuperLatticeFfromAnalyticalF2D<T, TDESCRIPTOR> solution(temperatureSol, lattice);
    lattice.setProcessingContext(ProcessingContext::Evaluation);


    vtkWriter.addFunctor( temperature );
    vtkWriter.addFunctor( solution );
    vtkWriter.write( iT );

    /// Creates contour plots as jpeg files for highest resolution run
    if (converter.getResolution() == maxN) {
      heatmap::plotParam<T> plotParam_temp;
      plotParam_temp.minValue = 0;
      plotParam_temp.maxValue = 2;
      plotParam_temp.contourlevel = 10;
      BlockReduction2D2D<T> planeReduction( temperature, converter.getResolution(), BlockDataSyncMode::ReduceOnly );
      heatmap::write( planeReduction, iT, plotParam_temp );
    }

    /// lattice statistics console output
    timer.update(iT);
    timer.printStep();
    lattice.getStatistics().print(iT, converter.getPhysTime( iT ));

    // compute relative error
    avg = errorEval(myCase, iT);
    clout << "Relative L2-error norm: "  << avg << std::endl;
  }

  return avg;
}


void simulate(MyCase& myCase)
{
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(AdvectionDiffusion{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  const auto& converter = lattice.getUnitConverter();

  const T statIter = parameters.get<parameters::LATTICE_STAT_ITER_T>();
  const std::size_t N = std::size_t(parameters.get<parameters::RESOLUTION>());

  OstreamManager clout(std::cout,"simulate");
  clout << "Executing the simulation with N=" << std::to_string(N) << std::endl;

  // switch outdirectory only if there are multiple simulation runs
  if (parameters.get<parameters::NUMBER_OF_RUNS>() > 1) {
    singleton::directories().setOutputDir("./tmp/N" + std::to_string(N) + "/");
  }

  util::Timer<T> timer( geometry.getStatistics().getNvoxel() );
  timer.start();

  T pulseDiffBound = parameters.get<parameters::PULSE_DIFF_BOUND>();
  int timeCount = util::ceil( -1./converter.getPhysDeltaT() * util::log(pulseDiffBound)/(2.*converter.getPhysDiffusivity()*util::pow(M_PI,2)) /statIter );
  int iTmax = timeCount * statIter;

  int iT;
  T simulationAverage = .0;
  for (iT = 0; iT < iTmax; ++iT) {
    simulationAverage += getResults(myCase, iT, statIter, timer);

    lattice.collideAndStream();
  }

  simulationAverage /= timeCount;
  parameters.set<parameters::L2_ERROR>(simulationAverage);

  clout << "Simulation Average Relative L2 Error: " << simulationAverage << std::endl;

  lattice.setProcessingContext(ProcessingContext::Evaluation);
  timer.stop();
  timer.printSummary();
}
