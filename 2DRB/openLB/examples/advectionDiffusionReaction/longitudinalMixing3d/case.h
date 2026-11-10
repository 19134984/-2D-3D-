/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2021 Stephan Simonis, Davide Dapelo,
 *                2024 Marc Heinzelmann
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

// Macro to switch between two approaches
#define BOUNDARY_CONDITION_APPROACH
//#define IN_BULK_APPROACH // IN_BULK_APPROACH only accurate for Peclet<0.1

using namespace olb;
using namespace olb::names;
using namespace olb::graphics;

using MyCase = Case<
  AdvectionDiffusion, Lattice<double, descriptors::D3Q7<descriptors::VELOCITY, descriptors::OMEGA>>
>;

namespace olb::parameters {

  struct ERROR_TIME : public descriptors::FIELD_BASE<1> {};
  struct C_EQ : public descriptors::FIELD_BASE<1> {};

}

/// functor for analytical solution
template<typename T, typename S>
class Solution : public AnalyticalF3D<T,S> {
private:
  MyCase& myCase;
  std::size_t iT;

public:

  Solution(MyCase& _myCase, std::size_t _iT) : AnalyticalF3D<T,S>(1), myCase(_myCase), iT(_iT) {};

  bool operator()(T output[], const S input[])
  {
    auto& lattice = myCase.getLattice(AdvectionDiffusion{});
    const auto& converter = lattice.getUnitConverter();

    const T t = converter.getPhysTime(iT);
    const T u = converter.getCharPhysVelocity();
    const T D = converter.getPhysDiffusivity();
    const T Ceq = myCase.getParameters().get<parameters::C_EQ>();

    const T x = input[0];
    output[0] = Ceq*( 0.5*erfc((x-u*t)/(2*sqrt(D*t)))
                        + sqrt(u*u*t/(M_PI*D)) * exp(-pow(x-u*t, 2)/(4*D*t))
                        - 0.5*(1+u*x/D+u*u*t/D) * exp(u*x/D) * erfc((x+u*t)/(2*sqrt(D*t))) );

    return true;
  }
};

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  const T physCharLength = params.get<parameters::PHYS_CHAR_LENGTH>();
  const T physDeltaX = physCharLength / params.get<parameters::RESOLUTION>();


  //minimal length in y,z direction
  Vector extent(physCharLength, physDeltaX, physDeltaX);
  Vector origin(0.0, 0.0, 0.0);

  IndicatorCuboid3D <T> cuboid(extent, origin);

  /// Instantiation of an empty cuboidDecomposition
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif

  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, noOfCuboids);
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({false,true,true});
  return mesh;
}


void prepareGeometry(MyCase& myCase) {
  // === MATERIAL NUMBERS === //
  // solid: MN = 0            //
  // bulk: MN = 1             //
  // wall: MN = 2             //
  // inlet: MN = 3            //
  // outlet: MN = 4           //
  // bulk reaction: MN = 5    //
  // ======================== //

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();


  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  const T delta = myCase.getParameters().get<parameters::PHYS_CHAR_LENGTH>() /
                  myCase.getParameters().get<parameters::RESOLUTION>();
  const T physLength = myCase.getParameters().get<parameters::PHYS_CHAR_LENGTH>();

  Vector<T,3> extend( 2*delta, physLength+2*delta, physLength+2*delta );
  Vector<T,3> origin(-delta, -delta, -delta);

  geometry.rename( 0,2 );
  geometry.rename( 2,1,{1,0,0} );

  // Set material number for inflow
  origin[0] = -delta;
  IndicatorCuboid3D<T> inflow( extend, origin );
  geometry.rename( 2,3,inflow );

#ifdef IN_BULK_APPROACH
  IndicatorCuboid3D<T> layer( extend, origin );
  geometry.rename( 1,5,layer );
#endif

  // Set material number for outflow
  origin[0] = physLength-delta;
  IndicatorCuboid3D<T> outflow( extend, origin );
  geometry.rename( 2,4,outflow );

  geometry.communicate();
  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}


void prepareLattice(MyCase& myCase)
{
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<AdvectionDiffusion>;
  using namespace olb::parameters;

  auto& ADlattice = myCase.getLattice(AdvectionDiffusion{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const T N = parameters.get<RESOLUTION>();
  const T relaxationTime = parameters.get<LATTICE_RELAXATION_TIME>();
  const T physLength = parameters.get<PHYS_CHAR_LENGTH>();
  const T Pe = parameters.get<PECLET>();
  const T physDiffusivity = parameters.get<PHYS_DIFFUSIVITY>();
  const T physDensity = parameters.get<PHYS_DENSITY>();
  const T Ceq = parameters.get<C_EQ>();

  ADlattice.setUnitConverter<AdeUnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    N,                              // Resolution
    relaxationTime,                 // Relaxation Time
    physLength,                     // charPhysLength
    Pe*physDiffusivity/physLength,  // charPhysVelocity from peclet
    physDiffusivity,                // physDiffusivity
    physDensity                     // physDensity
  );

  const auto& converter = ADlattice.getUnitConverter();
  converter.print();

  const T omega = converter.getLatticeRelaxationFrequency();

#ifdef BOUNDARY_CONDITION_APPROACH
  auto bulkIndicator = geometry.getMaterialIndicator({1,3,4});
  dynamics::set<AdvectionDiffusionBGKdynamics>( ADlattice, bulkIndicator );
  // reactive Robin boundary at left side MN=3
  boundary::set<boundary::Robin>(ADlattice, geometry.getMaterialIndicator({3}));
  T reactionRate = converter.getCharLatticeVelocity(); //determines speed of reaction, here equal to velocity
  AnalyticalConst3D <T,T> coefficients(reactionRate, -converter.getLatticeDiffusivity(), reactionRate * Ceq); //set robin coefficients as found in the paper
  fields::set<descriptors::G>( ADlattice, geometry.getMaterialIndicator({3}), (coefficients) ); //save coefficient in a field

  // ZeroGradient using RobinBoundary at right side MN=4
  boundary::set<boundary::Robin>(ADlattice, geometry.getMaterialIndicator({4}));
  AnalyticalConst3D <T,T> coefficients2(0., 1., 0.);
  fields::set<descriptors::G>( ADlattice, geometry.getMaterialIndicator({4}), (coefficients2) );
#endif


#ifdef IN_BULK_APPROACH
  auto bulkIndicator = geometry.getMaterialIndicator({1,3,4,5});
  dynamics::set<SourcedAdvectionDiffusionBGKdynamics>( ADlattice, bulkIndicator );
  // BounceBack for MN=3
  boundary::set<boundary::BounceBack>(ADlattice, geometry.getMaterialIndicator({3}));
  // Zero Gradient at right side MN=4
  setZeroGradientBoundary<T, DESCRIPTOR>(ADlattice, geometry.getMaterialIndicator({4}));
  // constants used in defineField
  AnalyticalConst3D<T,T> mat1(1.);
  AnalyticalConst3D<T,T> mat0(0.);
  // necessary for Zero Gradient Boundary
  fields::set<descriptors::SCALAR>( ADlattice, geometry.getMaterialIndicator({1,4}), mat1 );
  fields::set<descriptors::SCALAR>( ADlattice, geometry.getMaterialIndicator({0,3,5}), mat0 );
  // use field GAMMA to define where reaction is allowed to happen
  fields::set<descriptors::GAMMA>( ADlattice, geometry.getMaterialIndicator({5}), mat1 );
  fields::set<descriptors::GAMMA>( ADlattice, geometry.getMaterialIndicator({0,1,3,4}), mat0 );
#endif

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& ADlattice = myCase.getLattice(AdvectionDiffusion{});
  auto& geometry = myCase.getGeometry();
  const auto& converter = ADlattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  const T omega = converter.getLatticeRelaxationFrequency();

  #ifdef BOUNDARY_CONDITION_APPROACH
    const auto bulkIndicator = geometry.getMaterialIndicator({1,3,4});
  #endif

  #ifdef IN_BULK_APPROACH
    const auto bulkIndicator = geometry.getMaterialIndicator({1,3,4,5});
  #endif

  // Initial advection velocity and concentration
  AnalyticalConst3D <T, T> u0(converter.getCharLatticeVelocity(), 0.0, 0.0);
  AnalyticalConst3D<T,T> rho0(0.0);

  fields::set<descriptors::VELOCITY>( ADlattice, bulkIndicator, u0);
  momenta::setDensity( ADlattice, bulkIndicator, rho0);
  ADlattice.iniEquilibrium(bulkIndicator, rho0, u0);

  ADlattice.setParameter<descriptors::OMEGA>( omega );
  ADlattice.initialize();
}


// compute absolute and relative errors with different norms
void errorEval( MyCase& myCase, std::size_t iT)
{
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<AdvectionDiffusion>;
  auto& ADlattice = myCase.getLattice(AdvectionDiffusion{});
  auto& geometry = myCase.getGeometry();
  const auto& converter = ADlattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  OstreamManager clout(std::cout, "error");

  T result[3] = {T(), T(), T()};
  int tmp[] = {int()};

  SuperLatticeDensity3D <T, DESCRIPTOR> numericalConc(ADlattice);
  Solution<T,T> analyticalConc(myCase, iT);

#ifdef IN_BULK_APPROACH
  auto indicatorF = geometry.getMaterialIndicator({1,5}); // MN 3 and 4 should not be included for IN_BULK
#endif
#ifdef BOUNDARY_CONDITION_APPROACH
  auto indicatorF = geometry.getMaterialIndicator({1,3,4});
#endif

  SuperRelativeErrorL2Norm3D<T> relTemperatureErrorL2Norm(numericalConc, analyticalConc, indicatorF);
  SuperRelativeErrorL1Norm3D<T> relTemperatureErrorL1Norm(numericalConc, analyticalConc, indicatorF);
  SuperRelativeErrorLinfNorm3D<T> relTemperatureErrorLinfNorm(numericalConc, analyticalConc, indicatorF);

  /// Console output of every error norm
  relTemperatureErrorL2Norm(result, tmp);
  parameters.set<parameters::L2_ERROR>(result[0]);
  relTemperatureErrorL1Norm(result, tmp);
  parameters.set<parameters::L1_ERROR>(result[0]);
  relTemperatureErrorLinfNorm(result, tmp);
  parameters.set<parameters::LINF_ERROR>(result[0]);
}

void getResults(MyCase& myCase, const std::size_t iT,
                util::Timer <MyCase::value_t> &timer)
{
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<AdvectionDiffusion>;
  auto& ADlattice = myCase.getLattice(AdvectionDiffusion{});
  const auto& converter = ADlattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  const T statIter = parameters.get<parameters::PHYS_STAT_ITER_T>();
  const T errorTime = parameters.get<parameters::ERROR_TIME>();

  OstreamManager clout(std::cout, "getResults");
  SuperVTMwriter3D <T> vtkWriter("longitudinalMixing3d");

  SuperLatticeDensity3D <T, DESCRIPTOR> concentration(ADlattice);
  vtkWriter.addFunctor( concentration );
  concentration.getName() = "concentration";

  Solution<T,T> analyticalSolution(myCase, iT);
  SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> functor(analyticalSolution, ADlattice);
  functor.getName() = "analyticalSolution";
  vtkWriter.addFunctor( functor );

  if (iT == 0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization

    SuperLatticeCuboid3D <T, DESCRIPTOR> cuboid(ADlattice);
    SuperLatticeRank3D <T, DESCRIPTOR> rank(ADlattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
    vtkWriter.write(iT);

    timer.update(iT);
    timer.printStep();
    ADlattice.getStatistics().print(iT, converter.getPhysTime(iT));
  } else if (iT % converter.getLatticeTime(statIter) == 0) {
    ADlattice.setProcessingContext(ProcessingContext::Evaluation);

    /// ADLattice statistics console output
    timer.update(iT);
    timer.printStep();
    ADlattice.getStatistics().print(iT, converter.getPhysTime(iT));
    vtkWriter.write(iT);
  }
  if(iT == converter.getLatticeTime(errorTime)){
    errorEval(myCase, iT);
  }
}


void simulate(MyCase& myCase)
{
  using T = MyCase::value_t;
  auto& ADlattice = myCase.getLattice(AdvectionDiffusion{});
  auto& geometry = myCase.getGeometry();
  const auto& converter = ADlattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  const T N = converter.getResolution();
  const T Ceq = parameters.get<parameters::C_EQ>();
  const T maxPhysT = parameters.get<parameters::PHYS_MAX_T>();
  const T runs = parameters.get<parameters::NUMBER_OF_RUNS>();

  OstreamManager clout(std::cout, "simulate");
  clout << "Executing the simulation with N=" << std::to_string((std::size_t)N) << std::endl;

  // change directory for different runs
  if (runs > 1) {
    singleton::directories().setOutputDir("./tmp/N" + std::to_string((std::size_t)N) + "/");
  }

  util::Timer<T> timer(converter.getLatticeTime(maxPhysT), geometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {

#ifdef IN_BULK_APPROACH
    // Prepare reaction at left wall for IN_BULK_APPROACH
    SuperLatticeCoupling Reaction(
      LongitudinalMixingReactionCoupling<T>{},
      names::Concentration0{}, ADlattice);
    Reaction.template setParameter<typename LongitudinalMixingReactionCoupling<T>::
      REACTION_CONSTANT>({converter.getCharLatticeVelocity()}); //determines speed of reaction, here equal to velocity
    Reaction.template setParameter<typename LongitudinalMixingReactionCoupling<T>::EQUILIBRIUM>({Ceq}); //equilibrium concentration

    Reaction.apply();
#endif
    ADlattice.collideAndStream();

    getResults( myCase, iT, timer );
  }

  timer.stop();
  timer.printSummary();
}
