/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Mathias J. Krause
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

#include "analyticalSolutionTestFlow3D.h"
#include <vector>
#include <cmath>
#include <iostream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace std;

typedef double T;
#define DESCRIPTOR D3Q19<FORCE>

CuboidDecomposition3D<T>* cuboidDecomposition;
HeuristicLoadBalancer<T>* loadBalancer;

SuperGeometry<T,3> prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter, int geometryType)
{

  std::vector<T> origin(3,T());
  std::vector<T> extend(3,converter.getCharPhysLength());
  std::vector<T> center(origin);
  center[0]+=extend[0]/2.;
  center[1]+=extend[1]/2.;
  center[2]+=extend[2]/2.;
  IndicatorCuboid3D<T> cuboid(extend,origin);
  IndicatorSphere3D<T> sphere(center,converter.getCharPhysLength()/2.);
  IndicatorLayer3D<T> sphereLayer(sphere,converter.getPhysLength(1));

  cuboidDecomposition = new CuboidDecomposition3D<T>(sphereLayer, converter.getPhysLength(1), singleton::mpi().getSize() );
  loadBalancer = new HeuristicLoadBalancer<T>(*cuboidDecomposition);

  SuperGeometry<T,3> superGeometry(*cuboidDecomposition,*loadBalancer,2);

  superGeometry.rename(0,2,sphereLayer);
  superGeometry.rename(2,1,sphere);

  superGeometry.clean();
  superGeometry.checkForErrors();

  superGeometry.print();

  return superGeometry;
}

void prepareLattice(SuperLattice<T, DESCRIPTOR>& lattice,
                    UnitConverter<T,DESCRIPTOR> const& converter,
                    SuperGeometry<T,3>& superGeometry,
                    int N, int geometryType,
                    int bcType, int funcNum)
{
  OstreamManager cout(std::cout,"iniGeo");


  VelocityTestFlow3D<T,T,DESCRIPTOR> uSol(converter);
  AnalyticalScaled3D<T,T> uSolScaled(uSol,1./converter.getConversionFactorVelocity());


  std::vector<T> center(3,0.5);
  IndicatorSphere3D<T> sphere(center,converter.getCharPhysVelocity()/2.);

  ///material=0 -->do nothing
  lattice.defineDynamics<NoDynamics>(superGeometry, 0);

  if (bcType == 0) {
    ///material=1,2 -->bulk dynamics
    lattice.defineDynamics<ForcedBGKdynamics>(superGeometry.getMaterialIndicator({1, 2}));
    boundary::set<boundary::LocalVelocity>(lattice, superGeometry, 2);
  }
  else {
    ///material=1 -->bulk dynamics
    lattice.defineDynamics<ForcedBGKdynamics>(superGeometry, 1);
    setBouzidiBoundary<T,DESCRIPTOR,BouzidiVelocityPostProcessor>(lattice, superGeometry, 2, sphere);
  }

  // setting initial condition

  ForceTestFlow3D<T,T,DESCRIPTOR> forceSol(converter);

  AnalyticalConst3D<T,T> rhoF(1);
  std::vector<T> velocity(3,T());
  AnalyticalConst3D<T,T> uF(velocity);

  AnalyticalScaled3D<T,T> forceSolScaled(forceSol,1./converter.getConversionFactorForce()*converter.getConversionFactorMass());

  if (bcType == 0) {
    lattice.defineU(superGeometry, 2, uSolScaled);
  }
  else {
    setBouzidiVelocity(lattice, superGeometry, 2, uSolScaled);
  }
  lattice.defineField<descriptors::FORCE>(superGeometry.getMaterialIndicator({1, 2}), forceSolScaled);

  lattice.setParameter<OMEGA>(converter.getLatticeRelaxationFrequency());

  lattice.initialize();
}


void error(SuperGeometry<T,3>& superGeometry,
           SuperLattice<T, DESCRIPTOR>& lattice,
           UnitConverter<T,DESCRIPTOR> const& converter,
           int funcNum)
{
  OstreamManager clout(std::cout,"error");

  int tmp[1];
  T result[2];
  T result1[2];

  auto material = superGeometry.getMaterialIndicator(1);

  // velocity error
  VelocityTestFlow3D<T,T,DESCRIPTOR> uSol(converter);
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> u(lattice,converter);
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> uSolLattice(uSol,lattice);

  SuperL1Norm3D<T>   uSolL1Norm(uSolLattice,superGeometry.getMaterialIndicator(1));
  SuperL2Norm3D<T>   uSolL2Norm(uSolLattice,*material);
  SuperLinfNorm3D<T> uSolLinfNorm(uSolLattice,*material);

  SuperL1Norm3D<T>   uL1Norm(uSolLattice-u,*material);
  SuperL2Norm3D<T>   uL2Norm(uSolLattice-u,*material);
  SuperLinfNorm3D<T> uLinfNorm(uSolLattice-u,*material);

  uSolL1Norm(result1,tmp);
  uL1Norm(result,tmp);
  clout << "velocity-L1-error(util::abs)=" << result[0] << "; velocity-L1-error(rel)=" << result[0]/result1[0] << std::endl;

  uSolL2Norm(result1,tmp);
  uL2Norm(result,tmp);
  clout << "velocity-L2-error(util::abs)=" << result[0] << "; velocity-L2-error(rel)=" << result[0]/result1[0] << std::endl;

  uSolLinfNorm(result1,tmp);
  uLinfNorm(result,tmp);
  clout << "velocity-Linf-error(util::abs)=" << result[0] << "; velocity-Linf-error(rel)=" << result[0]/result1[0] << std::endl;

  // pressure error
  PressureTestFlow3D<T,T,DESCRIPTOR> pSol(converter);
  SuperLatticePhysPressure3D<T,DESCRIPTOR> p(lattice,converter);
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> pSolLattice(pSol,lattice);

  SuperL1Norm3D<T>   pSolL1Norm(pSolLattice,*material);
  SuperL2Norm3D<T>   pSolL2Norm(pSolLattice,*material);
  SuperLinfNorm3D<T> pSolLinfNorm(pSolLattice,*material);

  SuperL1Norm3D<T>   pL1Norm(pSolLattice-p,*material);
  SuperL2Norm3D<T>   pL2Norm(pSolLattice-p,*material);
  SuperLinfNorm3D<T> pLinfNorm(pSolLattice-p,*material);

  pSolL1Norm(result1,tmp);
  pL1Norm(result,tmp);
  clout << "pressure-L1-error(util::abs)=" << result[0] << "; pressure-L1-error(rel)=" << result[0]/result1[0] << std::endl;

  pSolL2Norm(result1,tmp);
  pL2Norm(result,tmp);
  clout << "pressure-L2-error(util::abs)=" << result[0] << "; pressure-L2-error(rel)=" << result[0]/result1[0] << std::endl;

  pSolLinfNorm(result1,tmp);
  pLinfNorm(result,tmp);
  clout << "pressure-Linf-error(util::abs)=" << result[0] << "; pressure-Linf-error(rel)=" << result[0]/result1[0] << std::endl;

  // strain rate error
  StrainRateTestFlow3D<T,T,DESCRIPTOR> sSol(converter);
  SuperLatticePhysStrainRate3D<T,DESCRIPTOR> s(lattice,converter);
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> sSolLattice(sSol,lattice);

  SuperL1Norm3D<T>   sSolL1Norm(sSolLattice,*material);
  SuperL2Norm3D<T>   sSolL2Norm(sSolLattice,*material);
  SuperLinfNorm3D<T> sSolLinfNorm(sSolLattice,*material);

  SuperL1Norm3D<T>   sL1Norm(sSolLattice-s,*material);
  SuperL2Norm3D<T>   sL2Norm(sSolLattice-s,*material);
  SuperLinfNorm3D<T> sLinfNorm(sSolLattice-s,*material);

  sSolL1Norm(result1,tmp);
  sL1Norm(result,tmp);
  clout << "strainRate-L1-error(util::abs)=" << result[0] << "; strainRate-L1-error(rel)=" << result[0]/result1[0] << std::endl;

  sSolL2Norm(result1,tmp);
  sL2Norm(result,tmp);
  clout << "strainRate-L2-error(util::abs)=" << result[0] << "; strainRate-L2-error(rel)=" << result[0]/result1[0] << std::endl;

  sSolLinfNorm(result1,tmp);
  sLinfNorm(result,tmp);
  clout << "strainRate-Linf-error(util::abs)=" << result[0] << "; strainRate-Linf-error(rel)=" << result[0]/result1[0] << std::endl;

  // stress error
  StressTensorTestFlow3D<T,T,DESCRIPTOR> stressSol(converter);
  SuperLatticePhysStressTensor3D<T,DESCRIPTOR> stress(lattice,converter);
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> stressSolLattice(stressSol,lattice);

  SuperL1Norm3D<T>   stressSolL1Norm(stressSolLattice,*material);
  SuperL2Norm3D<T>   stressSolL2Norm(stressSolLattice,*material);
  SuperLinfNorm3D<T> stressSolLinfNorm(stressSolLattice,*material);

  SuperL1Norm3D<T>   stressL1Norm(stressSolLattice-stress,*material);
  SuperL2Norm3D<T>   stressL2Norm(stressSolLattice-stress,*material);
  SuperLinfNorm3D<T> stressLinfNorm(stressSolLattice-stress,*material);

  stressSolL1Norm(result1,tmp);
  stressL1Norm(result,tmp);
  clout << "stress-L1-error(util::abs)=" << result[0] << "; stress-L1-error(rel)=" << result[0]/result1[0] << std::endl;

  stressSolL2Norm(result1,tmp);
  stressL2Norm(result,tmp);
  clout << "stress-L2-error(util::abs)=" << result[0] << "; stress-L2-error(rel)=" << result[0]/result1[0] << std::endl;

  stressSolLinfNorm(result1,tmp);
  stressLinfNorm(result,tmp);
  clout << "stress-Linf-error(util::abs)=" << result[0] << "; stress-Linf-error(rel)=" << result[0]/result1[0] << std::endl;

  // dissipation error
  DissipationTestFlow3D<T,T,DESCRIPTOR> dSol(converter);
  SuperLatticePhysDissipation3D<T,DESCRIPTOR> d(lattice,converter);
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> dSolLattice(dSol,lattice);

  SuperL1Norm3D<T>   dSolL1Norm(dSolLattice,*material);
  SuperL2Norm3D<T>   dSolL2Norm(dSolLattice,*material);
  SuperLinfNorm3D<T> dSolLinfNorm(dSolLattice,*material);

  SuperL1Norm3D<T>   dL1Norm(dSolLattice-d,*material);
  SuperL2Norm3D<T>   dL2Norm(dSolLattice-d,*material);
  SuperLinfNorm3D<T> dLinfNorm(dSolLattice-d,*material);

  dSolL1Norm(result1,tmp);
  dL1Norm(result,tmp);
  clout << "dissipation-L1-error(util::abs)=" << result[0] << "; dissipation-L1-error(rel)=" << result[0]/result1[0] << std::endl;

  dSolL2Norm(result1,tmp);
  dL2Norm(result,tmp);
  clout << "dissipation-L2-error(util::abs)=" << result[0] << "; dissipation-L2-error(rel)=" << result[0]/result1[0] << std::endl;

  dSolLinfNorm(result1,tmp);
  dLinfNorm(result,tmp);
  clout << "dissipation-Linf-error(util::abs)=" << result[0] << "; dissipation-Linf-error(rel)=" << result[0]/result1[0] << std::endl;
}





// VTK output
void writeVTK(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,3>& superGeometry,
              UnitConverter<T,DESCRIPTOR> const& converter, int iT)
{
  OstreamManager clout(std::cout,"getResults");
  SuperVTMwriter3D<T> vtkWriter("testFlow3d");

  if (iT==0) {
    /// Writes the converter log file
    // writeLogFile(converter, "testFlow3d");
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(sLattice);
    SuperLatticeRank3D<T, DESCRIPTOR> rank(sLattice);
    vtkWriter.write( cuboid );
    vtkWriter.write( rank );

    StressTensorTestFlow3D<T,T,DESCRIPTOR> stressTensorF(converter);
    SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> latticeStressTensorF(stressTensorF, sLattice);
    vtkWriter.write(latticeStressTensorF);

    vtkWriter.createMasterFile();
  }

  /// Writes the vtk files
  int iTperiodVTK = 10000;//converter.numTimeSteps(0.01);
  if (iT%iTperiodVTK==0) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
    vtkWriter.addFunctor( velocity );
    vtkWriter.addFunctor( pressure );
    vtkWriter.write(iT);
  }
}

void test(int TaskID,
          int const N,
          T const  epsilon,
          T const interval,
          int const geometryType,
          int const bcType,
          T const nu,
          T const uMax,
          T const maxPhysT,
          int const funcNum)
{
  OstreamManager cout(std::cout,"test");

  //params
  T charL = 1.;
  T charU = 1.;
  T charRho = 1200.;

  UnitConverterFromResolutionAndLatticeVelocity<T,DESCRIPTOR> const converter(
    int   {N},    // resolution: number of voxels per charPhysL
    (T)   uMax,   // charLatticeVelocity
    (T)   charL,  // charPhysLength: reference length of simulation geometry
    (T)   charU,  // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   nu,     // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   charRho // physDensity: physical density in __kg / m^3__
  );
  converter.print();

  int nT = converter.getLatticeTime(maxPhysT);

  // Prepare Geometry
  SuperGeometry<T,3> superGeometry(prepareGeometry(converter, geometryType) );

  // Prepare Lattice
  SuperLattice<T, DESCRIPTOR> sLattice(converter, superGeometry);

  prepareLattice(sLattice,
                 converter,
                 superGeometry,
                 N, geometryType,
                 bcType, funcNum);

  Timer<double> timer(nT, superGeometry.getStatistics().getNvoxel() );
  timer.start();

  T deltaT = converter.getLatticeTime(interval);
  util::ValueTracer<T> converge( deltaT, epsilon );

  int iT=0;
  for (iT=0; iT<=nT; ++iT) {
    //std::cout << iT << std::endl;
    if ( converge.hasConverged() ) {
      cout << "Simulation converged." << std::endl;
      sLattice.getStatistics().print(iT,converter.getPhysTime(iT));
      timer.print(iT);
      sLattice.setProcessingContext(ProcessingContext::Evaluation);
      error(superGeometry, sLattice, converter, funcNum);
      writeVTK(sLattice, superGeometry, converter, iT);
      break;
    }

    if (iT%1000==0) {
      sLattice.getStatistics().print(iT,converter.getPhysTime(iT));
      timer.print(iT);
      sLattice.setProcessingContext(ProcessingContext::Evaluation);
      error(superGeometry, sLattice, converter, funcNum);
    }

    sLattice.collideAndStream();

    //SuperLatticeDensity3D<T, DESCRIPTOR> rho(sLattice);
    //SuperAverage3D<T> avRho(rho, superGeometry, 1);

    sLattice.stripeOffDensityOffset(sLattice.getStatistics().getAverageRho() - T{1});

    //int tmp[4]; T output[2];
    //avRho(output,tmp);
    //sLattice.stripeOffDensityOffset ( output[0]-(T)1 );
    writeVTK(sLattice, superGeometry, converter, iT);
    converge.takeValue( sLattice.getStatistics().getAverageEnergy(), true);
  }
  sLattice.getStatistics().print(iT, converter.getPhysTime(iT));

  timer.stop();
  timer.printSummary();
}


/* main method
 * parameter list:
 * 1: N              default: 100
 * 2: epsilon        default: 1e-8
 * 3: interval       default: 1.0
 * 4: geometryType   default: (0) Sphere
 * 5: bcType         default: (1) bouzidi, others: (0) onLattice, (2) BB
 * 6: nu             default: 0.01 (Re = 1/nu)
 * 7: uMaxStart      default: 1/N
 * 8: maxPhysT       default: 15
 * 9: iterations     default: 1
 * 10: uFrac          default: 3/4
 * 11: funcNum       default: 1
 */

int main(int argc, char **argv)
{
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout,"main");

  int N = 100;
  T epsilon = 1e-8;
  T interval = 1.0;  // timespan considered by the value tracer
  int geometryType = 0;
  int bcType = 1;
  T nu = 0.01;
  T maxPhysT = 30;
  int iterations = 1;
  T uFrac = 0.75;
  int funcNum = 1;

  if (argc > 1) {
    N = atoi(argv[1]);
  }
  if (argc > 2) {
    epsilon = atof(argv[2]);
  }
  if (argc > 3) {
    interval = atof(argv[3]);
  }
  if (argc > 4) {
    geometryType = atoi(argv[4]);
  }
  if (argc > 5) {
    bcType = atoi(argv[5]);
  }
  if (argc > 6) {
    nu = atof(argv[6]);
  }
  T uMaxStart = (T) 1./N;
  if (argc > 7) {
    uMaxStart = atof(argv[7]);
  }
  if (argc > 8) {
    maxPhysT = atof(argv[8]);
  }
  if (argc > 9) {
    iterations = atoi(argv[9]);
  }
  if (argc > 10) {
    uFrac = atof(argv[10]);
  }
  if (argc > 11) {
    funcNum = atoi(argv[11]);
  }

  clout << "-------Parameters-------" << std::endl;
  clout << "N:           " << N << std::endl;
  clout << "epsilon:     " << epsilon << std::endl;
  clout << "interval:    " << interval << std::endl;
  clout << "geometryType:" << geometryType << std::endl;
  clout << "cells:       " << N*N*N << std::endl;
  clout << "bcType:      " << bcType << std::endl;
  clout << "nu:          " << nu << std::endl;
  clout << "uMaxStart:   " << uMaxStart << std::endl;
  clout << "maxPhysT:    " << maxPhysT << std::endl;
  clout << "iterations:  " << iterations << std::endl;
  clout << "uFrac:       " << uFrac << std::endl;
  clout << "funcNum:     " << funcNum << std::endl;

  int rank;
#ifdef PARALLEL_MODE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif

  for (int i = 0; i < iterations; i++) {
    clout << "Starting next simulation " << i << " with latticeU = " << uMaxStart << std::endl;
    test(rank, N, epsilon, interval, geometryType, bcType, nu, uMaxStart, maxPhysT, funcNum);
    uMaxStart *= uFrac;
  }

  delete cuboidDecomposition;
  delete loadBalancer;

  return 0;
}
