/*  This file is part of the OpenLB library
*
*  Copyright (C) 2021 Nicolas Hafen
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
#include <iostream>
#include <type_traits>
#include <gtest/gtest.h>
#include "helper/unitTestHelper.h"
#include "helper/framework.h"

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::particles;

typedef double T;

::testing::Environment* const env = ::testing::AddGlobalTestEnvironment(new OlbInitTestEnvironment);


#define VERBOSE_TESTS



TEST(TestMomentumExchangeFunctor, TestParticle2DCircleCentred)
{

  //Settings
  typedef D2Q9<> DESCRIPTOR;
  constexpr unsigned D = DESCRIPTOR::d;
  typedef ResolvedParticle2D PARTICLETYPE;
  PhysR<T,D> origin (0.0);
  const T baseSize = 1.0;
  T delX = 0.1;
  PhysR<T,D> particlePosition( .5*baseSize );
  T particleCircumRadius = 0.2*baseSize;
  Vector<T,D> fluidVelocity(2., 0.);

  //Set up environment
  Vector<T,D> extend( baseSize );
  IndicatorCuboid<T,D> cuboid( extend, origin );
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 4;
#endif
  CuboidDecomposition<T,D> cuboidDecomposition( cuboid, delX, noOfCuboids );
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,D> sGeometry( cuboidDecomposition, loadBalancer, 2 );
  sGeometry.rename(0, 1);
  sGeometry.clean();
  UnitConverter<T,DESCRIPTOR> converter( delX, delX, baseSize, 1., 1., 1.);
  SuperLattice<T, DESCRIPTOR> sLattice( converter, sGeometry );

  auto bulkIndicator = sGeometry.getMaterialIndicator(1);
  sLattice.defineDynamics<BGKdynamics>( bulkIndicator );

  //Initialize relaxation frequency
  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  //Set Velocity
  AnalyticalConst<D,T,T> density( 1.0 );
  AnalyticalConst<D,T,T> velocity( fluidVelocity );
  sLattice.defineRhoU( bulkIndicator, density, velocity );
  sLattice.iniEquilibrium( bulkIndicator, density, velocity );
  sLattice.initialize();

  //Create particle
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  T epsilon = .5*delX;
  T radius = particleCircumRadius-.5*epsilon;
  T rhoP = 1.0;
  creators::addResolvedCircle2D( particleSystem, particlePosition,
                                 radius, epsilon, rhoP );
  auto particle = particleSystem.get(0);

  //Retrieve hydrodynamic force and apply to particle
  SuperLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE> forceFunctor(
    sLattice, sGeometry, particleSystem, converter );
  particles::dynamics::applySerializableParticleForce<T,PARTICLETYPE>( forceFunctor, particleSystem );
  auto force = particle.getField<FORCING,FORCE>();

#ifdef VERBOSE_TESTS
  std::cout <<  "Force:     " << std::setprecision(16) << force
            << std::setprecision(6) << std::endl;
#endif

  //Test specifc
  EXPECT_NEAR(force[0], 5.329e-17, 1e-15);

//  SuperVTMwriter2D<T> vtkWriterLattice( "fluid" );
//  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocityF( sLattice, converter );
//  SuperGeometryF<T, DESCRIPTOR> geometryF( sGeometry );
//  SuperLatticeCuboid2D<T, DESCRIPTOR> cuboidF( sLattice );
//  vtkWriterLattice.write ( velocityF );
//  vtkWriterLattice.write ( geometryF );
//  vtkWriterLattice.write ( cuboidF );

}


TEST(TestMomentumExchangeFunctor, TestParticle3DSphereCentred)
{

  //Settings
  typedef D3Q19<> DESCRIPTOR;
  constexpr unsigned D = DESCRIPTOR::d;
  typedef ResolvedParticle3D PARTICLETYPE;
  PhysR<T,D> origin (0.0);
  const T baseSize = 1.0;
  T delX = 0.1;
  PhysR<T,D> particlePosition( .5*baseSize );
  T particleCircumRadius = 0.2*baseSize;
  Vector<T,D> fluidVelocity(2., 0., 0.);

  //Set up environment
  Vector<T,D> extend( baseSize );
  IndicatorCuboid<T,D> cuboid( extend, origin );
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 4;
#endif
  CuboidDecomposition<T,D> cuboidDecomposition( cuboid, delX, noOfCuboids );
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,D> sGeometry( cuboidDecomposition, loadBalancer, 2 );
  sGeometry.rename(0, 1);
  sGeometry.clean();
  UnitConverter<T,DESCRIPTOR> converter( delX, delX, baseSize, 1., 1., 1.);
  SuperLattice<T, DESCRIPTOR> sLattice( converter, sGeometry );
  auto bulkIndicator = sGeometry.getMaterialIndicator(1);
  sLattice.defineDynamics<BGKdynamics>( bulkIndicator );

  //Initialize relaxation frequency
  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  //Set Velocity
  AnalyticalConst<D,T,T> density( 1.0 );
  AnalyticalConst<D,T,T> velocity( fluidVelocity );
  sLattice.defineRhoU( bulkIndicator, density, velocity );
  sLattice.iniEquilibrium( bulkIndicator, density, velocity );
  sLattice.initialize();

  //Create particle
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  T epsilon = .5*delX;
  T radius = particleCircumRadius-.5*epsilon;
  T rhoP = 1.0;
  creators::addResolvedSphere3D( particleSystem, particlePosition,
                                 radius, epsilon, rhoP );
  auto particle = particleSystem.get(0);

  //Retrieve hydrodynamic force and apply to particle
  SuperLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE> forceFunctor(
    sLattice, sGeometry, particleSystem, converter );
  particles::dynamics::applySerializableParticleForce<T,PARTICLETYPE>( forceFunctor, particleSystem );
  auto force = particle.getField<FORCING,FORCE>();

#ifdef VERBOSE_TESTS
  std::cout <<  "Force:     " << std::setprecision(16) << force
            << std::setprecision(6) << std::endl;
#endif

  //Test specifc
  EXPECT_NEAR(force[0], 9.992e-18, 1e-15);

}

// Periodic boundary test - only checks the processed number of cells
TEST(TestMomentumExchangeFunctor, TestParticle3DCuboidPeriodicInX)
{
   //Settings
  typedef D3Q19<> DESCRIPTOR;
  constexpr unsigned D = DESCRIPTOR::d;
  typedef ResolvedParticle3D PARTICLETYPE;
  const PhysR<T,D> origin (0.0);
  const T baseSize = 1.0;
  const T delX = 0.1;
  const T particleCircumRadius = 0.2*baseSize;
  const T epsilon = .5*delX;
  const T radius = particleCircumRadius-.5*epsilon;
  const T rhoP = 1.0;
  const Vector<T,D> fluidVelocity(2., 0., 0.);
  const Vector<bool, D> periodicity(true, false, false);
  UnitConverter<T,DESCRIPTOR> converter( delX, delX, baseSize, 1., 1., 1.);

  // Setup Geometry
  Vector<T,D> extend( baseSize );
  IndicatorCuboid<T,D> cuboid( extend, origin );
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 4;
#endif
  CuboidDecomposition<T,D> cuboidDecomposition( cuboid, delX, noOfCuboids );
  cuboidDecomposition.setPeriodicity(periodicity);
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,D> sGeometry( cuboidDecomposition, loadBalancer, 2 );
  sGeometry.rename(0, 1);
  sGeometry.clean();

  // Setup Lattice
  SuperLattice<T, DESCRIPTOR> sLattice(converter, sGeometry );
  auto bulkIndicator = sGeometry.getMaterialIndicator(1);
  sLattice.defineDynamics<BGKdynamics>( bulkIndicator );
  //Initialize relaxation frequency
  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  //Set Velocity
  AnalyticalConst<D,T,T> density( 1.0 );
  AnalyticalConst<D,T,T> velocity( fluidVelocity );
  sLattice.defineRhoU( bulkIndicator, density, velocity );
  sLattice.iniEquilibrium( bulkIndicator, density, velocity );
  sLattice.initialize();

  // Do the whole procedure 3 times
  // 1. Get the number of cells without usage of the periodic boundary
  // 2. Get the number of cells with periodic boundary
  // 3. Get the number of cells with periodic boundary again, but on the other side of the domain
  int numVoxels[3];

  for(unsigned i = 0; i < 3; ++i) {
    PhysR<T,D> particlePosition( .5*baseSize );
    if(i == 1) {
      particlePosition[0] = baseSize - delX;
    }
    else if(i == 2) {
      particlePosition[0] = delX;
    }

    //Create particle
    ParticleSystem<T,PARTICLETYPE> particleSystem;
    creators::addResolvedCuboid3D( particleSystem, particlePosition,
                                   Vector<T,D>(2 * radius), epsilon, rhoP );

    //Retrieve hydrodynamic force and apply to particle
    SuperLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE> forceFunctor(
      sLattice, sGeometry, particleSystem, converter, periodicity );
    T output[forceFunctor.getTargetDim()];
    for(int j = 0; j < forceFunctor.getTargetDim(); ++j) {
      output[j] = T{0};
    }
    int input[1]{};
    forceFunctor(output, input);

    numVoxels[i] = output[forceFunctor.getTargetDim() - 1];
  }

  // Check that the number of cells are equal
  EXPECT_EQ(numVoxels[0], numVoxels[1]);
  EXPECT_EQ(numVoxels[0], numVoxels[2]);
}

TEST(TestMomentumExchangeFunctor, TestParticle3DCuboidPeriodicInY)
{
   //Settings
  typedef D3Q19<> DESCRIPTOR;
  constexpr unsigned D = DESCRIPTOR::d;
  typedef ResolvedParticle3D PARTICLETYPE;
  const PhysR<T,D> origin (0.0);
  const T baseSize = 1.0;
  const T delX = 0.1;
  const T particleCircumRadius = 0.2*baseSize;
  const T epsilon = .5*delX;
  const T radius = particleCircumRadius-.5*epsilon;
  const T rhoP = 1.0;
  const Vector<T,D> fluidVelocity(2., 0., 0.);
  const Vector<bool, D> periodicity(false, true, false);
  UnitConverter<T,DESCRIPTOR> converter( delX, delX, baseSize, 1., 1., 1.);

  // Setup Geometry
  Vector<T,D> extend( baseSize );
  IndicatorCuboid<T,D> cuboid( extend, origin );
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 4;
#endif
  CuboidDecomposition<T,D> cuboidDecomposition( cuboid, delX, noOfCuboids );
  cuboidDecomposition.setPeriodicity(periodicity);
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,D> sGeometry( cuboidDecomposition, loadBalancer, 2 );
  sGeometry.rename(0, 1);
  sGeometry.clean();

  // Setup Lattice
  SuperLattice<T, DESCRIPTOR> sLattice( converter,sGeometry );
  auto bulkIndicator = sGeometry.getMaterialIndicator(1);
  sLattice.defineDynamics<BGKdynamics>( bulkIndicator );
  //Initialize relaxation frequency
  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  //Set Velocity
  AnalyticalConst<D,T,T> density( 1.0 );
  AnalyticalConst<D,T,T> velocity( fluidVelocity );
  sLattice.defineRhoU( bulkIndicator, density, velocity );
  sLattice.iniEquilibrium( bulkIndicator, density, velocity );
  sLattice.initialize();

  // Do the whole procedure 3 times
  // 1. Get the number of cells without usage of the periodic boundary
  // 2. Get the number of cells with periodic boundary
  // 3. Get the number of cells with periodic boundary again, but on the other side of the domain
  int numVoxels[3];

  for(unsigned i = 0; i < 3; ++i) {
    PhysR<T,D> particlePosition( .5*baseSize );
    if(i == 1) {
      particlePosition[1] = baseSize - delX;
    }
    else if(i == 2) {
      particlePosition[1] = delX;
    }

    //Create particle
    ParticleSystem<T,PARTICLETYPE> particleSystem;
    creators::addResolvedCuboid3D( particleSystem, particlePosition,
                                   Vector<T,D>(2 * radius), epsilon, rhoP );

    //Retrieve hydrodynamic force and apply to particle
    SuperLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE> forceFunctor(
      sLattice, sGeometry, particleSystem, converter, periodicity );
    T output[forceFunctor.getTargetDim()];
    for(int j = 0; j < forceFunctor.getTargetDim(); ++j) {
      output[j] = T{0};
    }
    int input[1]{};
    forceFunctor(output, input);

    numVoxels[i] = output[forceFunctor.getTargetDim() - 1];
  }

  // Check that the number of cells are equal
  EXPECT_EQ(numVoxels[0], numVoxels[1]);
  EXPECT_EQ(numVoxels[0], numVoxels[2]);
}

TEST(TestMomentumExchangeFunctor, TestParticle3DCuboidPeriodicInZ)
{
   //Settings
  typedef D3Q19<> DESCRIPTOR;
  constexpr unsigned D = DESCRIPTOR::d;
  typedef ResolvedParticle3D PARTICLETYPE;
  const PhysR<T,D> origin (0.0);
  const T baseSize = 1.0;
  const T delX = 0.1;
  const T particleCircumRadius = 0.2*baseSize;
  const T epsilon = .5*delX;
  const T radius = particleCircumRadius-.5*epsilon;
  const T rhoP = 1.0;
  const Vector<T,D> fluidVelocity(2., 0., 0.);
  const Vector<bool, D> periodicity(false, false, true);
  UnitConverter<T,DESCRIPTOR> converter( delX, delX, baseSize, 1., 1., 1.);

  // Setup Geometry
  Vector<T,D> extend( baseSize );
  IndicatorCuboid<T,D> cuboid( extend, origin );
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 4;
#endif
  CuboidDecomposition<T,D> cuboidDecomposition( cuboid, delX, noOfCuboids );
  cuboidDecomposition.setPeriodicity(periodicity);
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,D> sGeometry( cuboidDecomposition, loadBalancer, 2 );
  sGeometry.rename(0, 1);
  sGeometry.clean();

  // Setup Lattice
  SuperLattice<T, DESCRIPTOR> sLattice( converter, sGeometry );
  auto bulkIndicator = sGeometry.getMaterialIndicator(1);
  sLattice.defineDynamics<BGKdynamics>( bulkIndicator );
  //Initialize relaxation frequency
  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  //Set Velocity
  AnalyticalConst<D,T,T> density( 1.0 );
  AnalyticalConst<D,T,T> velocity( fluidVelocity );
  sLattice.defineRhoU( bulkIndicator, density, velocity );
  sLattice.iniEquilibrium( bulkIndicator, density, velocity );
  sLattice.initialize();

  // Do the whole procedure 3 times
  // 1. Get the number of cells without usage of the periodic boundary
  // 2. Get the number of cells with periodic boundary
  // 3. Get the number of cells with periodic boundary again, but on the other side of the domain
  int numVoxels[3];

  for(unsigned i = 0; i < 3; ++i) {
    PhysR<T,D> particlePosition( .5*baseSize );
    if(i == 1) {
      particlePosition[2] = baseSize - delX;
    }
    else if(i == 2) {
      particlePosition[2] = delX;
    }

    //Create particle
    ParticleSystem<T,PARTICLETYPE> particleSystem;
    creators::addResolvedCuboid3D( particleSystem, particlePosition,
                                   Vector<T,D>(2 * radius), epsilon, rhoP );

    //Retrieve hydrodynamic force and apply to particle
    SuperLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE> forceFunctor(
      sLattice, sGeometry, particleSystem, converter, periodicity );
    T output[forceFunctor.getTargetDim()];
    for(int j = 0; j < forceFunctor.getTargetDim(); ++j) {
      output[j] = T{0};
    }
    int input[1]{};
    forceFunctor(output, input);

    numVoxels[i] = output[forceFunctor.getTargetDim() - 1];
  }

  // Check that the number of cells are equal
  EXPECT_EQ(numVoxels[0], numVoxels[1]);
  EXPECT_EQ(numVoxels[0], numVoxels[2]);
}

TEST(TestMomentumExchangeFunctor, TestParticle3DCuboidPeriodicInAll)
{
   //Settings
  typedef D3Q19<> DESCRIPTOR;
  constexpr unsigned D = DESCRIPTOR::d;
  typedef ResolvedParticle3D PARTICLETYPE;
  const PhysR<T,D> origin (0.0);
  const T baseSize = 1.0;
  const T delX = 0.1;
  const T particleCircumRadius = 0.2*baseSize;
  const T epsilon = .5*delX;
  const T radius = particleCircumRadius-.5*epsilon;
  const T rhoP = 1.0;
  const Vector<T,D> fluidVelocity(2., 0., 0.);
  const Vector<bool, D> periodicity(true, true, true);
  UnitConverter<T,DESCRIPTOR> converter( delX, delX, baseSize, 1., 1., 1.);

  // Setup Geometry
  Vector<T,D> extend( baseSize );
  IndicatorCuboid<T,D> cuboid( extend, origin );
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 4;
#endif
  CuboidDecomposition<T,D> cuboidDecomposition( cuboid, delX, noOfCuboids );
  cuboidDecomposition.setPeriodicity(periodicity);
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,D> sGeometry( cuboidDecomposition, loadBalancer, 2 );
  sGeometry.rename(0, 1);
  sGeometry.clean();

  // Setup Lattice
  SuperLattice<T, DESCRIPTOR> sLattice( converter,sGeometry );
  auto bulkIndicator = sGeometry.getMaterialIndicator(1);
  sLattice.defineDynamics<BGKdynamics>( bulkIndicator );
  //Initialize relaxation frequency
  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  //Set Velocity
  AnalyticalConst<D,T,T> density( 1.0 );
  AnalyticalConst<D,T,T> velocity( fluidVelocity );
  sLattice.defineRhoU( bulkIndicator, density, velocity );
  sLattice.iniEquilibrium( bulkIndicator, density, velocity );
  sLattice.initialize();

  // Do the whole procedure 3 times
  // 1. Get the number of cells without usage of the periodic boundary
  // 2. Get the number of cells with periodic boundary
  // 3. Get the number of cells with periodic boundary again, but on the other side of the domain
  int numVoxels[3];

  for(unsigned i = 0; i < 3; ++i) {
    PhysR<T,D> particlePosition( .5*baseSize );
    if(i == 1) {
      particlePosition = PhysR<T,D>(baseSize - delX);
    }
    else if(i == 2) {
      particlePosition = PhysR<T,D>(delX);
    }

    //Create particle
    ParticleSystem<T,PARTICLETYPE> particleSystem;
    creators::addResolvedCuboid3D( particleSystem, particlePosition,
                                   Vector<T,D>(2 * radius), epsilon, rhoP );

    //Retrieve hydrodynamic force and apply to particle
    SuperLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE> forceFunctor(
      sLattice, sGeometry, particleSystem, converter, periodicity );
    T output[forceFunctor.getTargetDim()];
    for(int j = 0; j < forceFunctor.getTargetDim(); ++j) {
      output[j] = T{0};
    }
    int input[1]{};
    forceFunctor(output, input);

    numVoxels[i] = output[forceFunctor.getTargetDim() - 1];
  }

  // Check that the number of cells are equal
  EXPECT_EQ(numVoxels[0], numVoxels[1]);
  EXPECT_EQ(numVoxels[0], numVoxels[2]);
}


TEST(TestMomentumExchangeFunctor, TestParticle2DCuboidPeriodicInX)
{
   //Settings
  typedef D2Q9<> DESCRIPTOR;
  constexpr unsigned D = DESCRIPTOR::d;
  typedef ResolvedParticle2D PARTICLETYPE;
  const PhysR<T,D> origin (0.0);
  const T baseSize = 1.0;
  const T delX = 0.1;
  const T particleCircumRadius = 0.2*baseSize;
  const T epsilon = .5*delX;
  const T radius = particleCircumRadius-.5*epsilon;
  const T rhoP = 1.0;
  const Vector<T,D> fluidVelocity(2., 0.);
  const Vector<bool, D> periodicity(true, false);
  UnitConverter<T,DESCRIPTOR> converter( delX, delX, baseSize, 1., 1., 1.);

  // Setup Geometry
  Vector<T,D> extend( baseSize );
  IndicatorCuboid<T,D> cuboid( extend, origin );
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 4;
#endif
  CuboidDecomposition<T,D> cuboidDecomposition( cuboid, delX, noOfCuboids );
  cuboidDecomposition.setPeriodicity(periodicity);
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,D> sGeometry( cuboidDecomposition, loadBalancer, 2 );
  sGeometry.rename(0, 1);
  sGeometry.clean();

  // Setup Lattice
  SuperLattice<T, DESCRIPTOR> sLattice( converter,sGeometry );
  auto bulkIndicator = sGeometry.getMaterialIndicator(1);
  sLattice.defineDynamics<BGKdynamics>( bulkIndicator );
  //Initialize relaxation frequency
  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  //Set Velocity
  AnalyticalConst<D,T,T> density( 1.0 );
  AnalyticalConst<D,T,T> velocity( fluidVelocity );
  sLattice.defineRhoU( bulkIndicator, density, velocity );
  sLattice.iniEquilibrium( bulkIndicator, density, velocity );
  sLattice.initialize();

  // Do the whole procedure 3 times
  // 1. Get the number of cells without usage of the periodic boundary
  // 2. Get the number of cells with periodic boundary
  // 3. Get the number of cells with periodic boundary again, but on the other side of the domain
  int numVoxels[3];

  for(unsigned i = 0; i < 3; ++i) {
    PhysR<T,D> particlePosition( .5*baseSize );
    if(i == 1) {
      particlePosition[0] = baseSize - delX;
    }
    else if(i == 2) {
      particlePosition[0] = delX;
    }

    //Create particle
    ParticleSystem<T,PARTICLETYPE> particleSystem;
    creators::addResolvedCuboid2D( particleSystem, particlePosition,
                                   Vector<T,D>(2 * radius), epsilon, rhoP );

    //Retrieve hydrodynamic force and apply to particle
    SuperLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE> forceFunctor(
      sLattice, sGeometry, particleSystem, converter, periodicity );
    T output[forceFunctor.getTargetDim()];
    for(int j = 0; j < forceFunctor.getTargetDim(); ++j) {
      output[j] = T{0};
    }
    int input[1]{};
    forceFunctor(output, input);

    numVoxels[i] = output[forceFunctor.getTargetDim() - 1];
  }

  // Check that the number of cells are equal
  EXPECT_EQ(numVoxels[0], numVoxels[1]);
  EXPECT_EQ(numVoxels[0], numVoxels[2]);
}

TEST(TestMomentumExchangeFunctor, TestParticle2DCuboidPeriodicInY)
{
   //Settings
  typedef D2Q9<> DESCRIPTOR;
  constexpr unsigned D = DESCRIPTOR::d;
  typedef ResolvedParticle2D PARTICLETYPE;
  const PhysR<T,D> origin (0.0);
  const T baseSize = 1.0;
  const T delX = 0.1;
  const T particleCircumRadius = 0.2*baseSize;
  const T epsilon = .5*delX;
  const T radius = particleCircumRadius-.5*epsilon;
  const T rhoP = 1.0;
  const Vector<T,D> fluidVelocity(2., 0.);
  const Vector<bool, D> periodicity(false, true);
  UnitConverter<T,DESCRIPTOR> converter( delX, delX, baseSize, 1., 1., 1.);

  // Setup Geometry
  Vector<T,D> extend( baseSize );
  IndicatorCuboid<T,D> cuboid( extend, origin );
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 4;
#endif
  CuboidDecomposition<T,D> cuboidDecomposition( cuboid, delX, noOfCuboids );
  cuboidDecomposition.setPeriodicity(periodicity);
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,D> sGeometry( cuboidDecomposition, loadBalancer, 2 );
  sGeometry.rename(0, 1);
  sGeometry.clean();

  // Setup Lattice
  SuperLattice<T, DESCRIPTOR> sLattice( converter,sGeometry );
  auto bulkIndicator = sGeometry.getMaterialIndicator(1);
  sLattice.defineDynamics<BGKdynamics>( bulkIndicator );
  //Initialize relaxation frequency
  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  //Set Velocity
  AnalyticalConst<D,T,T> density( 1.0 );
  AnalyticalConst<D,T,T> velocity( fluidVelocity );
  sLattice.defineRhoU( bulkIndicator, density, velocity );
  sLattice.iniEquilibrium( bulkIndicator, density, velocity );
  sLattice.initialize();

  // Do the whole procedure 3 times
  // 1. Get the number of cells without usage of the periodic boundary
  // 2. Get the number of cells with periodic boundary
  // 3. Get the number of cells with periodic boundary again, but on the other side of the domain
  int numVoxels[3];

  for(unsigned i = 0; i < 3; ++i) {
    PhysR<T,D> particlePosition( .5*baseSize );
    if(i == 1) {
      particlePosition[1] = baseSize - delX;
    }
    else if(i == 2) {
      particlePosition[1] = delX;
    }

    //Create particle
    ParticleSystem<T,PARTICLETYPE> particleSystem;
    creators::addResolvedCuboid2D( particleSystem, particlePosition,
                                   Vector<T,D>(2 * radius), epsilon, rhoP );

    //Retrieve hydrodynamic force and apply to particle
    SuperLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE> forceFunctor(
      sLattice, sGeometry, particleSystem, converter, periodicity );
    T output[forceFunctor.getTargetDim()];
    for(int j = 0; j < forceFunctor.getTargetDim(); ++j) {
      output[j] = T{0};
    }
    int input[1]{};
    forceFunctor(output, input);

    numVoxels[i] = output[forceFunctor.getTargetDim() - 1];
  }

  // Check that the number of cells are equal
  EXPECT_EQ(numVoxels[0], numVoxels[1]);
  EXPECT_EQ(numVoxels[0], numVoxels[2]);
}

TEST(TestMomentumExchangeFunctor, TestParticle2DCuboidPeriodicInAll)
{
   //Settings
  typedef D2Q9<> DESCRIPTOR;
  constexpr unsigned D = DESCRIPTOR::d;
  typedef ResolvedParticle2D PARTICLETYPE;
  const PhysR<T,D> origin (0.0);
  const T baseSize = 1.0;
  const T delX = 0.1;
  const T particleCircumRadius = 0.2*baseSize;
  const T epsilon = .5*delX;
  const T radius = particleCircumRadius-.5*epsilon;
  const T rhoP = 1.0;
  const Vector<T,D> fluidVelocity(2., 0.);
  const Vector<bool, D> periodicity(true, true);
  UnitConverter<T,DESCRIPTOR> converter( delX, delX, baseSize, 1., 1., 1.);

  // Setup Geometry
  Vector<T,D> extend( baseSize );
  IndicatorCuboid<T,D> cuboid( extend, origin );
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 4;
#endif
  CuboidDecomposition<T,D> cuboidDecomposition( cuboid, delX, noOfCuboids );
  cuboidDecomposition.setPeriodicity(periodicity);
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,D> sGeometry( cuboidDecomposition, loadBalancer, 2 );
  sGeometry.rename(0, 1);
  sGeometry.clean();

  // Setup Lattice
  SuperLattice<T, DESCRIPTOR> sLattice( converter, sGeometry );
  auto bulkIndicator = sGeometry.getMaterialIndicator(1);
  sLattice.defineDynamics<BGKdynamics>( bulkIndicator );
  //Initialize relaxation frequency
  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  //Set Velocity
  AnalyticalConst<D,T,T> density( 1.0 );
  AnalyticalConst<D,T,T> velocity( fluidVelocity );
  sLattice.defineRhoU( bulkIndicator, density, velocity );
  sLattice.iniEquilibrium( bulkIndicator, density, velocity );
  sLattice.initialize();

  // Do the whole procedure 3 times
  // 1. Get the number of cells without usage of the periodic boundary
  // 2. Get the number of cells with periodic boundary
  // 3. Get the number of cells with periodic boundary again, but on the other side of the domain
  int numVoxels[3];

  for(unsigned i = 0; i < 3; ++i) {
    PhysR<T,D> particlePosition( .5*baseSize );
    if(i == 1) {
      particlePosition = PhysR<T,D>(baseSize - delX);
    }
    else if(i == 2) {
      particlePosition = PhysR<T,D>(delX);
    }

    //Create particle
    ParticleSystem<T,PARTICLETYPE> particleSystem;
    creators::addResolvedCuboid2D( particleSystem, particlePosition,
                                   Vector<T,D>(2 * radius), epsilon, rhoP );

    //Retrieve hydrodynamic force and apply to particle
    SuperLatticeMomentumExchangeForce<T, DESCRIPTOR, PARTICLETYPE> forceFunctor(
      sLattice, sGeometry, particleSystem, converter, periodicity );
    T output[forceFunctor.getTargetDim()];
    for(int j = 0; j < forceFunctor.getTargetDim(); ++j) {
      output[j] = T{0};
    }
    int input[1]{};
    forceFunctor(output, input);

    numVoxels[i] = output[forceFunctor.getTargetDim() - 1];
  }

  // Check that the number of cells are equal
  EXPECT_EQ(numVoxels[0], numVoxels[1]);
  EXPECT_EQ(numVoxels[0], numVoxels[2]);
}
