
/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2021 Stephan Simonis, Davide Dapelo
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
using namespace olb::graphics;
using namespace olb::names;

using MyCase = Case<
  AdvectionDiffusion, Lattice<double, descriptors::D3Q7<descriptors::VELOCITY>>
>;

/// uncomment for console output of other norms
#define allNorms

namespace olb::parameters {
  struct PULSE_DIFF_BOUND : public descriptors::FIELD_BASE<1> {};
  struct NON_SMOOTH_ENABLED : public descriptors::TYPED_FIELD_BASE<bool,1> {};
}

template <typename T>
class AdePhysTemp3D : public AnalyticalF3D<T,T> {

protected:
  MyCase& myCase;
  T t;

public:
    AdePhysTemp3D(T time, MyCase& _myCase) :
      AnalyticalF3D<T, T>(1),
      myCase(_myCase),
      t(time)
      {};

    bool operator()(T output[], const T input[]) override {
      auto& lattice = myCase.getLattice(AdvectionDiffusion{});
      auto& parameters = myCase.getParameters();
      const auto& converter = lattice.getUnitConverter();
      const T mue = converter.getPhysDiffusivity();
      const T uMag = converter.getCharPhysVelocity();
      const T res = converter.getResolution();
      const T deltaT = converter.getPhysDeltaT();
      const bool nonsmooth = parameters.get<parameters::NON_SMOOTH_ENABLED>();

      const T x = input[0];
      const T y = input[1];
      const T z = input[2];
      T gf = res/(res+1.); // domain correction
      if (nonsmooth){
        gf = 1.;
      }
      const T uMagX = uMag;
      const T uMagY = uMag;
      const T uMagZ = uMag;

      if (nonsmooth) {
        /// dirac comb (see https://doi.org/10.1016/j.jocs.2021.101363)
        const T startingPt = 0.;
        const T tShift = deltaT;
        const int numSpikes = 10;
        output[0] = 0.;
        for (int i=-numSpikes; i<=numSpikes; ++i){
          output[0] += util::exp( -util::pow(gf*(x-uMagX*(t+tShift)-startingPt-i*(res+1)*util::sqrt(deltaT)),2) / (4*(t+tShift)*mue*util::pow(gf,2)) );
        }
        const T initialConcentr = 1./util::sqrt(4*M_PI*mue);
        output[0] *= initialConcentr * ( 1 / util::sqrt((t+tShift)*util::pow(gf,2)) );
        output[0] += 1.;
      }
      else {
        /// smooth modes (3D extrude of https://doi.org/10.1098/rsta.2019.0400)
        output[0] = util::sin(M_PI*(x-uMagX*t)*gf);
        output[0] *= util::sin(M_PI*(y-uMagY*t)*gf);
        output[0] *= util::sin(M_PI*(z-uMagZ*t)*gf);
        output[0] *= util::exp(-3*mue*util::pow(M_PI,2)*t*util::pow(gf,2));
        output[0] += 1.;
      }
      return true;
    };
};

template<typename T>
class AdePhysTemp3Dinit : public AnalyticalF3D<T, T> {

protected:
  MyCase& myCase;
  T t;

public:
  AdePhysTemp3Dinit(T time, MyCase& _myCase) :
    AnalyticalF3D<T, T>(1),
    myCase(_myCase),
    t(time)
    {};

  bool operator()(T output[], const T input[]) override {
    auto& lattice = myCase.getLattice(AdvectionDiffusion{});
    auto& parameters = myCase.getParameters();
    const auto& converter = lattice.getUnitConverter();
    const T mue = converter.getPhysDiffusivity();
    const T uMag = converter.getCharPhysVelocity();
    const T res = converter.getResolution();
    const T deltaT = converter.getPhysDeltaT();
    const bool nonsmooth = parameters.get<parameters::NON_SMOOTH_ENABLED>();

    if (nonsmooth) {
      const T x = input[0];
      const T gf = 1.;
      /// naive dirac comb (single spike)
      if (std::fabs(x) <= 1. / res) {
          output[0] = 1 / util::sqrt(deltaT * util::pow(gf, 2));
          output[0] += 1.;
      } else {
          output[0] = 1.;
      }
    }
    else {
      const T x = input[0];
      const T y = input[1];
      const T z = input[2];
      const T gf = res/(res+1.); // domain correction
      /// smooth initial condition (3D)
      const T uMagX = uMag;
      const T uMagY = uMag;
      const T uMagZ = uMag;
      output[0] = util::sin(M_PI*(x-uMagX*t)*gf);
      output[0] *= util::sin(M_PI*(y-uMagY*t)*gf);
      output[0] *= util::sin(M_PI*(z-uMagZ*t)*gf);
      output[0] *= util::exp(-3*mue*util::pow(M_PI,2)*t*util::pow(gf,2));
      output[0] += 1.;
    }

    return true;
  };
};

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  const T physCharLength = params.get<parameters::PHYS_CHAR_LENGTH>();

  Vector extent(physCharLength, physCharLength, physCharLength);
  Vector origin(-physCharLength/2, -physCharLength/2, -physCharLength/2);

  IndicatorCuboid3D <T> cuboid(extent, origin);

  /// Instantiation of an empty cuboidDecomposition
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif

  const T physDeltaX = physCharLength / params.get<parameters::RESOLUTION>();
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, noOfCuboids);
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({true,true,true});
  return mesh;
}

// geometry preparation for cube
void prepareGeometry(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();

  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  geometry.rename(0, 1);
  geometry.communicate();
  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}


void prepareLattice(MyCase& myCase) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<AdvectionDiffusion>;
  using namespace olb::parameters;

  auto& lattice = myCase.getLattice(AdvectionDiffusion{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const T physDeltaX = parameters.get<PHYS_CHAR_LENGTH>() / parameters.get<RESOLUTION>();
  const T physDeltaT = physDeltaX * physDeltaX;
  const T physLength = parameters.get<PHYS_CHAR_LENGTH>();
  const T mue = parameters.get<PHYS_CHAR_VELOCITY>() * physLength / parameters.get<PECLET>();
  const T physVelocity = parameters.get<PHYS_CHAR_VELOCITY>();
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

  dynamics::set<AdvectionDiffusionBGKdynamics>(lattice, geometry.getMaterialIndicator({1}));
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(AdvectionDiffusion{});
  auto& geometry = myCase.getGeometry();
  const auto& converter = lattice.getUnitConverter();
  auto& parameters = myCase.getParameters();
  const bool nonsmooth = parameters.get<parameters::NON_SMOOTH_ENABLED>();

  // Initial advection velocity
  AnalyticalConst3D<T,T> u0( converter.getCharLatticeVelocity(),
                             converter.getCharLatticeVelocity(),
                             converter.getCharLatticeVelocity() );
  if(nonsmooth){
      AnalyticalConst3D <T, T> u0(converter.getCharLatticeVelocity(), 0.0, 0.0);
  }

  AdePhysTemp3Dinit<T> Tinit( 0.0, myCase);

  auto bulkIndicator = geometry.getMaterialIndicator({0, 1});

  fields::set<descriptors::VELOCITY>( lattice, bulkIndicator, u0 );
  momenta::setDensity( lattice, bulkIndicator, Tinit );

  lattice.setParameter<descriptors::OMEGA>( converter.getLatticeAdeRelaxationFrequency() );
  lattice.initialize();
}

// compute absolute and relative errors with different norms
MyCase::value_t errorEval( MyCase& myCase, const int iT )
{
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<AdvectionDiffusion>;
  auto& lattice = myCase.getLattice(AdvectionDiffusion{});
  auto& geometry = myCase.getGeometry();
  const auto& converter = lattice.getUnitConverter();

  OstreamManager clout(std::cout, "error");
  CSV<T> csvWriter("errors");

  T result[3] = {T(), T(), T()};
  int tmp[] = {int()};

  SuperLatticeDensity3D <T, DESCRIPTOR> temperature(lattice);
  AdePhysTemp3D<T> temperatureSol(converter.getPhysTime(iT), myCase);

  auto indicatorF = geometry.getMaterialIndicator({1});

  T l2rel_error;

  SuperRelativeErrorL2Norm3D <T> relTemperatureErrorL2Norm(temperature, temperatureSol, indicatorF);
  SuperAbsoluteErrorL2Norm3D <T> absTemperatureErrorL2Norm(temperature, temperatureSol, indicatorF);

  /// Console output of other error norms
#ifdef allNorms
    SuperAbsoluteErrorL1Norm3D<T> absTemperatureErrorL1Norm(temperature, temperatureSol, indicatorF);
    absTemperatureErrorL1Norm(result, tmp);
    clout << "temperature-L1-absolute-error=" << result[0] << std::endl;
    csvWriter.writeDataFile(iT, result[0], "l1-abs-error");

    SuperRelativeErrorL1Norm3D<T> relTemperatureErrorL1Norm(temperature, temperatureSol, indicatorF);
    relTemperatureErrorL1Norm(result, tmp);
    clout << "temperature-L1-relative-error=" << result[0] << std::endl;
    csvWriter.writeDataFile(iT, result[0],  "l1-rel-error");

    SuperAbsoluteErrorLinfNorm3D<T> absTemperatureErrorLinfNorm(temperature, temperatureSol, indicatorF);
    absTemperatureErrorLinfNorm(result, tmp);
    clout << "temperature-Linf-absolute-error=" << result[0] << std::endl;
    csvWriter.writeDataFile(iT, result[0],  "linf-abs-error");

    SuperRelativeErrorLinfNorm3D<T> relTemperatureErrorLinfNorm(temperature, temperatureSol, indicatorF);
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

MyCase::value_t getResults( MyCase& myCase,
             int iT,
             int statIter,
             util::Timer <MyCase::value_t>& timer)
{
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<AdvectionDiffusion>;
  auto& lattice = myCase.getLattice(AdvectionDiffusion{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  const auto& converter = lattice.getUnitConverter();

  OstreamManager clout(std::cout, "getResults");
  SuperVTMwriter3D <T> vtkWriter("advectionDiffusion3d");

  AdePhysTemp3D<T> temperatureSol(converter.getPhysTime(iT), myCase);

  T avg = 0;

  if (iT == 0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D <T, DESCRIPTOR> cuboid(lattice);
    SuperLatticeRank3D <T, DESCRIPTOR> rank(lattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();

    SuperLatticeDensity3D <T, DESCRIPTOR> temperature(lattice);
    SuperLatticeFfromAnalyticalF3D <T, DESCRIPTOR> solution(temperatureSol, lattice);

    vtkWriter.addFunctor(temperature);
    vtkWriter.addFunctor(solution);
    vtkWriter.write(iT);

    timer.update(iT);
    timer.printStep();
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
  } else if (iT % statIter == 0 && iT != 0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    /// Writes the VTK files
    SuperLatticeDensity3D <T, DESCRIPTOR> temperature(lattice);
    SuperLatticeFfromAnalyticalF3D <T, DESCRIPTOR> solution(temperatureSol, lattice);

    vtkWriter.addFunctor(temperature);
    vtkWriter.addFunctor(solution);
    vtkWriter.write(iT);

    /// lattice statistics console output
    timer.update(iT);
    timer.printStep();
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));

    // compute relative error
    avg = errorEval(myCase, iT);
    clout << "Relative L2-error norm: " << avg << std::endl;
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
  const std::size_t runs = parameters.get<parameters::NUMBER_OF_RUNS>();
  const std::size_t N = parameters.get<parameters::RESOLUTION>();
  const std::size_t statIter = parameters.get<parameters::LATTICE_STAT_ITER_T>();
  const T peclet = parameters.get<parameters::PECLET>();
  const T physTmax = parameters.get<parameters::PHYS_MAX_T>();
  const T pulseDiffBound = parameters.get<parameters::PULSE_DIFF_BOUND>();
  const bool nonsmooth = parameters.get<parameters::NON_SMOOTH_ENABLED>();

  OstreamManager clout(std::cout, "simulate");
  clout << "Executing the simulation with N=" << std::to_string(N) << std::endl;

  // switch outdirectory if there multiple simulation runs
  // store results in subfolder corresponding to peclet number "p_XX"
  if (runs > 1) {
      singleton::directories().setOutputDir("./tmp/p_" + std::to_string((int)peclet) + "/N" + std::to_string(N) + "/");
  }

  int iTmax;
  int timeCount;

  if (nonsmooth){
    iTmax = converter.getLatticeTime(physTmax);
    timeCount = ceil(iTmax / statIter);
  } else {
    timeCount = ceil( -1./converter.getPhysDeltaT() * log(pulseDiffBound)/(3.*converter.getPhysDiffusivity()*util::pow(M_PI,2)) /statIter );
    iTmax = timeCount * statIter;
  }

  util::Timer<T> timer(iTmax, geometry.getStatistics().getNvoxel());
  timer.start();

  int iT;
  T simulationAverage = .0;
  for (iT = 0; iT < iTmax; ++iT) {
      simulationAverage += getResults(myCase, iT, statIter, timer);

      lattice.collideAndStream();
  }

  simulationAverage /= timeCount;
  parameters.set<parameters::L2_ERROR>(simulationAverage);

  clout << "Simulation Average Relative L2 Error: " << simulationAverage << std::endl;

  timer.stop();
  timer.printSummary();
}
