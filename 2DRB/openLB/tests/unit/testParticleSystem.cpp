/*  This file is part of the OpenLB library
*
*  Copyright (C) 2022 Nicolas Hafen
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
using namespace olb::particles::dynamics;
//using namespace olb::particles::access;

using SubgridParticle3DsingleCuboid  = PARTICLE_DESCRIPTOR<3, GENERAL_TMP<3>,
      MOBILITY_VERLET_NO_ANGLE<3>, FORCING_SUBGRID<3>,
      PHYSPROPERTIES_SUBGRID<3>,
      DYNBEHAVIOUR_BASIC>;

typedef double T;
typedef ResolvedParticle3D            PARTICLETYPE;
typedef ResolvedDecomposedParticle3D    PARTICLETYPEPARA;
typedef SubgridParticle3DsingleCuboid PARTICLETYPESUB;

::testing::Environment* const env = ::testing::AddGlobalTestEnvironment(new OlbInitTestEnvironment);

//#define VERBOSE_TESTS


TEST(TestParticleSystem, ConstructorPlain)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  EXPECT_EQ(particleSystem.size(), 0);
}

TEST(TestParticleSystem, ConstructorAllocated)
{
  std::size_t noP = 10;
  ParticleSystem<T,PARTICLETYPE> particleSystem(noP);
  EXPECT_EQ(particleSystem.size(), noP);
}

TEST(TestParticleSystem, AddAndRetrieveDynamics)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
//  VerletParticleDynamics<T,PARTICLETYPE> dynamics;
  std::shared_ptr<ParticleDynamics<T,PARTICLETYPE>> dynamicsSPtr =
    std::make_shared<VerletParticleDynamics<T,PARTICLETYPE>>();
  particleSystem.addDynamics(dynamicsSPtr);
  auto& dynVec = particleSystem.getDynamicsVector();
  auto dynamicsPtr = particleSystem.getDynamics();
  std::string name = dynamicsPtr->getName();
  EXPECT_EQ(dynVec.size(), 1);
  EXPECT_EQ(name, "VerletParticleDynamics");
}

TEST(TestParticleSystem, ParticleInterface)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);
  auto particle = particleSystem.get(1);
  auto iParticle = particle.getId();
  EXPECT_EQ(iParticle, 1);
}

TEST(TestParticleSystem, ParticleInterfaceOperator)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem(5);
  auto particle = particleSystem[2];
  auto iParticle = particle.getId();
  EXPECT_EQ(iParticle, 2);
}

TEST(TestParticleSystem, ParticleInterfaceIteratorEnd)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem(4);
  //Retrieve iterator to first entry
  auto pIterator = particleSystem.begin();
  //Increment by postfix and prefix operator
  pIterator++;
  ++pIterator;
  //Retrieve particle by dereferencing
  auto particle = *pIterator;
  //Get id and compare
  auto iParticle = particle.getId();
  EXPECT_EQ(iParticle, 2);
}

TEST(TestParticleSystem, ParticleInterfaceIteratorForLoopA)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem(5);
  int idCount = 0;
  //Loop over particles
  for (auto particle : particleSystem){
    idCount += particle.getId();
  }
  EXPECT_EQ(idCount, 10);
}

TEST(TestParticleSystem, ParticleInterfaceIteratorForLoopB)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem(3);
  int idCount = 0;
  //Loop over particles
  for (auto it = particleSystem.begin(), end = particleSystem.end(); it != end; ++it) {
    auto particle = *it;
    idCount += particle.getId();
  }
  EXPECT_EQ(idCount, 3);
}

TEST(TestParticleSystem, AddAndRetrieveDynamicsViaParticle)
{
  //Add dynamics via particle interface
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);
  std::shared_ptr<ParticleDynamics<T,PARTICLETYPE>> dynamicsSPtr =
    std::make_shared<VerletParticleDynamics<T,PARTICLETYPE>>();
  auto particle = particleSystem.get(1);
  particle.addDynamics(dynamicsSPtr);
  //Retrieve dynamics via particle interface
  auto particleB = particleSystem.get(0);
  particleB.advanceId();
  auto dynamicsPtr = particle.getDynamics();
  std::string name = dynamicsPtr->getName();
  EXPECT_EQ(name, "VerletParticleDynamics");
}

TEST(TestParticleSystem, MultiDynamics)
{
  //Create Particlesystem
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);
  //Add dynamics A
  Vector<T,3> domainMin(0.);
  Vector<T,3> domainExt(1.);
  IndicatorCuboid3D<T> cuboid( domainExt, domainMin );
  std::unique_ptr<IndicatorF<T,3>> wallInd =
      std::make_unique<IndicInverse<T,3>>(cuboid);
  SolidBoundary<T,3> wall(std::move(wallInd), 2, 0);
  std::shared_ptr<ParticleDynamics<T,PARTICLETYPE>> dynamicsASPtr =
    std::make_shared<VerletParticleDynamicsVelocityWallReflection<T,PARTICLETYPE>>(
      wall);
  particleSystem.addDynamics(dynamicsASPtr);
  //Add dynamics B
  std::shared_ptr<ParticleDynamics<T,PARTICLETYPE>> dynamicsBSPtr =
    std::make_shared<VerletParticleDynamics<T,PARTICLETYPE>>();
  particleSystem.addDynamics(dynamicsBSPtr);
  //Retrieve dynamics via particle interface
  auto particle = particleSystem.get(0);
  auto dynamicsPtrA = particle.getDynamics(0);
  std::string nameA = dynamicsPtrA->getName();
  auto dynamicsPtrB = particle.getDynamics(1);
  std::string nameB = dynamicsPtrB->getName();
  EXPECT_EQ(nameA, "VerletParticleDynamicsVelocityWallReflection");
  EXPECT_EQ(nameB, "VerletParticleDynamics");
}

TEST(TestParticleSystem, Process)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);
  std::shared_ptr<ParticleDynamics<T,PARTICLETYPE>> dynamicsSPtr =
    std::make_shared<VerletParticleDynamics<T,PARTICLETYPE>>();
  particleSystem.addDynamics(dynamicsSPtr);
  //Particle interface
  auto particle = particleSystem.get(0);
  T timeStep = 0.01;
  particle.process(timeStep);
  particle.process(timeStep,0);
}

TEST(TestSuperParticleSystem, addAndRetrieveDynamics)
{
  Vector<T,3> suGEOo (0.);
  Vector<T,3> suGEOe (1.);
  IndicatorCuboid3D<T> cuboid( suGEOe, suGEOo );
  const int noOfCuboids = 2;
  CuboidDecomposition3D<T> cuboidDecomposition( cuboid, 0.1, noOfCuboids );
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,3> sGeometry( cuboidDecomposition, loadBalancer, 2 );
  SuperParticleSystem<T,PARTICLETYPEPARA> sParticleSystem(sGeometry);
  //Add dynamics
  std::shared_ptr<ParticleDynamics<T,PARTICLETYPEPARA>> dynamicsASPtr =
    std::make_shared<VerletParticleDynamics<T,PARTICLETYPEPARA>>();
  sParticleSystem.addDynamics(dynamicsASPtr);
  std::unique_ptr<IndicatorF<T,3>> wallInd =
      std::make_unique<IndicInverse<T,3>>(cuboid);
  SolidBoundary<T,3> wall(std::move(wallInd), 2, 0);
  std::shared_ptr<ParticleDynamics<T,PARTICLETYPEPARA>> dynamicsBSPtr =
    std::make_shared<VerletParticleDynamicsVelocityWallReflection<T,PARTICLETYPEPARA>>(
      wall);
  sParticleSystem.addDynamics(dynamicsBSPtr);
  //Retrieve particle systems
  auto& bParticleSystems = sParticleSystem.getBlockParticleSystems();
  int noSys = bParticleSystems.size();
  auto& particleSystemA = *bParticleSystems[0];
  auto& particleSystemB = *bParticleSystems[1];
  auto dynamicsPtrA = particleSystemA.getDynamics();
  std::string nameA = dynamicsPtrA->getName();
  auto dynamicsPtrB = particleSystemB.getDynamics(1);
  std::string nameB = dynamicsPtrB->getName();
  EXPECT_EQ(noSys, noOfCuboids);
  EXPECT_EQ(nameA, "VerletParticleDynamics");
  EXPECT_EQ(nameB, "VerletParticleDynamicsVelocityWallReflection");
}

TEST(TestSuperParticleSystem, addAndRetrieveMultiDynamics)
{
  Vector<T,3> suGEOo (0.);
  Vector<T,3> suGEOe (1.);
  IndicatorCuboid3D<T> cuboid( suGEOe, suGEOo );
  const int noOfCuboids = 2;
  CuboidDecomposition3D<T> cuboidDecomposition( cuboid, 0.1, noOfCuboids );
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,3> sGeometry( cuboidDecomposition, loadBalancer, 2 );
  SuperParticleSystem<T,PARTICLETYPEPARA> sParticleSystem(sGeometry);
  //Add dynamics
  std::shared_ptr<ParticleDynamics<T,PARTICLETYPEPARA>> dynamicsSPtr =
    std::make_shared<VerletParticleDynamics<T,PARTICLETYPEPARA>>();
  sParticleSystem.addDynamics(dynamicsSPtr);
  //Retrieve particle systems
  auto& bParticleSystems = sParticleSystem.getBlockParticleSystems();
  int noSys = bParticleSystems.size();
  auto& particleSystemA = *bParticleSystems[0];
  auto& particleSystemB = *bParticleSystems[1];
  auto dynamicsPtrA = particleSystemA.getDynamics();
  std::string nameA = dynamicsPtrA->getName();
  auto dynamicsPtrB = particleSystemB.getDynamics();
  std::string nameB = dynamicsPtrB->getName();
  EXPECT_EQ(noSys, noOfCuboids);
  EXPECT_EQ(nameA, "VerletParticleDynamics");
  EXPECT_EQ(nameB, "VerletParticleDynamics");
}

TEST(TestParticleSystem, VerletParticleDynamicsTranslationOnly)
{
  //Create Particlesystem
  ParticleSystem<T,PARTICLETYPESUB> particleSystem(3);
  //Add dynamics A
  std::shared_ptr<ParticleDynamics<T,PARTICLETYPESUB>> dynamicsASPtr =
    std::make_shared<VerletParticleDynamicsTranslationOnly<T,PARTICLETYPESUB>>();
  particleSystem.addDynamics(dynamicsASPtr);
  //Add dynamics B
  //Retrieve dynamics via particle interface
  auto particle = particleSystem.get(0);
  auto dynamicsPtrA = particle.getDynamics(0);
  std::string nameA = dynamicsPtrA->getName();
  EXPECT_EQ(nameA, "VerletParticleDynamicsTranslationOnly");
}

TEST(TestParticleSystem, VerletParticleDynamicsRotationOnly)
{
  //Create Particlesystem
  ParticleSystem<T,PARTICLETYPE> particleSystem(5);
  //Add dynamics A
  std::shared_ptr<ParticleDynamics<T,PARTICLETYPE>> dynamicsASPtr =
    std::make_shared<VerletParticleDynamicsRotationOnly<T,PARTICLETYPE>>();
  particleSystem.addDynamics(dynamicsASPtr);
  //Add dynamics B
  //Retrieve dynamics via particle interface
  auto particle = particleSystem.get(0);
  auto dynamicsPtrA = particle.getDynamics(0);
  std::string nameA = dynamicsPtrA->getName();
  EXPECT_EQ(nameA, "VerletParticleDynamicsRotationOnly");
}

TEST(TestParticleSystem, VerletParticleDynamicsWallCaptureDynamics)
{
  //Create Particlesystem
  ParticleSystem<T,PARTICLETYPESUB> particleSystem(2);
  //Add dynamics A
  Vector<T,3> domainMin(0.);
  Vector<T,3> domainExt(1.);
  IndicatorCuboid3D<T> cuboid( domainExt, domainMin );
  std::unique_ptr<IndicatorF<T,3>> wallInd =
      std::make_unique<IndicInverse<T,3>>(cuboid);
  SolidBoundary<T,3> wall(std::move(wallInd), 2, 0);
  std::shared_ptr<ParticleDynamics<T,PARTICLETYPESUB>> dynamicsASPtr =
    std::make_shared<VerletParticleDynamicsWallCapture<T,PARTICLETYPESUB>>(
      wall);
  particleSystem.addDynamics(dynamicsASPtr);
  //Add dynamics B
  //Retrieve dynamics via particle interface
  auto particle = particleSystem.get(0);
  auto dynamicsPtrA = particle.getDynamics(0);
  std::string nameA = dynamicsPtrA->getName();
  EXPECT_EQ(nameA, "VerletParticleDynamicsWallCapture");
}

//Factory functions via defineDynamics()

TEST(TestParticleSystem, AddAndRetrieveDynamicsFactory)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  particleSystem.defineDynamics<VerletParticleDynamics<T,PARTICLETYPE>>();
  particleSystem.defineDynamics<VerletParticleDynamicsTranslationOnly<T,PARTICLETYPE>>();
  auto& dynVec = particleSystem.getDynamicsVector();
  auto dynamicsPtr = particleSystem.getDynamics(0);
  auto dynamicsBPtr = particleSystem.getDynamics(1);
  std::string name = dynamicsPtr->getName();
  std::string nameB = dynamicsBPtr->getName();
  EXPECT_EQ(dynVec.size(), 2);
  EXPECT_EQ(name, "VerletParticleDynamics");
  EXPECT_EQ(nameB, "VerletParticleDynamicsTranslationOnly");
}

TEST(TestSuperParticleSystem, addAndRetrieveDynamicsFactory)
{
  Vector<T,3> suGEOo (0.);
  Vector<T,3> suGEOe (1.);
  IndicatorCuboid3D<T> cuboid( suGEOe, suGEOo );
  const int noOfCuboids = 2;
  CuboidDecomposition3D<T> cuboidDecomposition( cuboid, 0.1, noOfCuboids );
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,3> sGeometry( cuboidDecomposition, loadBalancer, 2 );
  SuperParticleSystem<T,PARTICLETYPEPARA> sParticleSystem(sGeometry);
  //Add dynamics
  sParticleSystem.defineDynamics<VerletParticleDynamics<T,PARTICLETYPEPARA>>();
  std::unique_ptr<IndicatorF<T,3>> wallInd =
      std::make_unique<IndicInverse<T,3>>(cuboid);
  SolidBoundary<T,3> wall(std::move(wallInd), 2, 0);
  sParticleSystem.defineDynamics<
    VerletParticleDynamicsVelocityWallReflection<T,PARTICLETYPEPARA>>(wall);
  //Retrieve particle systems
  auto& bParticleSystems = sParticleSystem.getBlockParticleSystems();
  int noSys = bParticleSystems.size();
  auto& particleSystemA = *bParticleSystems[0];
  auto& particleSystemB = *bParticleSystems[1];
  auto dynamicsPtrA = particleSystemA.getDynamics();
  std::string nameA = dynamicsPtrA->getName();
  auto dynamicsPtrB = particleSystemB.getDynamics(1);
  std::string nameB = dynamicsPtrB->getName();
  EXPECT_EQ(noSys, noOfCuboids);
  EXPECT_EQ(nameA, "VerletParticleDynamics");
  EXPECT_EQ(nameB, "VerletParticleDynamicsVelocityWallReflection");
}
