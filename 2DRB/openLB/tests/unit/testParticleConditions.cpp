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
#include <gtest/gtest.h>
#include "helper/unitTestHelper.h"
#include "helper/framework.h"

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::particles;
using namespace olb::particles::conditions;
using namespace olb::particles::creators;
using namespace olb::particles::communication;
using namespace olb::particles::access;

typedef double T;

::testing::Environment* const env = ::testing::AddGlobalTestEnvironment(new OlbInitTestEnvironment);

//#define VERBOSE_TESTS



/// ParticleSystem

TEST(TestParticleConditions, particle_matching_ID)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  ParticleSystem<T,PARTICLETYPE> particleSystem(10);

  using PCOND = particle_matching_ID<2>;
  forParticlesInParticleSystem<T,PARTICLETYPE,PCOND>( particleSystem,
    [](auto& particle){
#ifdef VERBOSE_TESTS
    particle.print();
#endif
  });
//  EXPECT_EQ(position, valInit);
//  EXPECT_EQ(position.getDim(), PARTICLETYPE::d);
}

TEST(TestParticleConditions, particle_mathing_ID_parallel)
{
  typedef ResolvedDecomposedParticle3D PARTICLETYPE;

  Vector<T,3> suGEOo (0.);
  Vector<T,3> suGEOe (1.);
  IndicatorCuboid3D<T> cuboid( suGEOe, suGEOo );
  const int noOfCuboids = 2;
  CuboidDecomposition3D<T> cuboidDecomposition( cuboid, 0.1, noOfCuboids );
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,3> sGeometry( cuboidDecomposition, loadBalancer, 2 );
  SuperParticleSystem<T,PARTICLETYPE> sParticleSystem(sGeometry);
  //TODO: add actual particles

  using PCOND = valid_particle_matching_ID<2>;
  communication::forParticlesInSuperParticleSystem<T,PARTICLETYPE,PCOND>( sParticleSystem,
    [](auto& particle,auto& particleSystem,int globiC){
#ifdef VERBOSE_TESTS
    particle.print();
#endif
  });
//  EXPECT_EQ(position, valInit);
//  EXPECT_EQ(position.getDim(), PARTICLETYPE::d);
}

//TODO: to be moved to particle utilities
TEST(TestParticleConditions, particle_matching_ID_util)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  ParticleSystem<T,PARTICLETYPE> particleSystem(5);
  //Do for particle mathing specific ID
  doForParticleMatchingID<T,PARTICLETYPE,3>( particleSystem,
    [](auto& particle){
#ifdef VERBOSE_TESTS
    particle.print();
#endif
  });

//  EXPECT_EQ(position, valInit);
//  EXPECT_EQ(position.getDim(), PARTICLETYPE::d);
}


/// SuperParticleSystem


TEST(TestParticleConditions, resolvedConditions)
{
  //Settings
  typedef ResolvedDecomposedParticle3D PARTICLETYPE;
  constexpr unsigned D = PARTICLETYPE::d;
  PhysR<T,D> origin (0.0);
  const T baseSize = 1.0;
  const T delX = 0.1;
  //Set up environment
  Vector<T,D> extend( baseSize );
  IndicatorCuboid<T,D> cuboid( extend, origin );
  const int noOfCuboids = 27;
  CuboidDecomposition<T,D> cuboidDecomposition( cuboid, delX, noOfCuboids );
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,D> sGeometry( cuboidDecomposition, loadBalancer, 2 );

  //Create super particle system
  SuperParticleSystem<T,PARTICLETYPE> sParticleSystem(sGeometry);

  //Define particle properties
  T radius = .2*baseSize;
  T epsilon = .5*delX;
  std::size_t noP = 9;
  for (std::size_t iP=0; iP<noP; ++iP){
    //Define position
    Vector<T,3> position(.1*iP*baseSize);
    //Add Sphere
    addResolvedSphere3D( sParticleSystem, position, radius, epsilon );
  }

  //Iterate over particles and randomly assign traits by globalId
  forParticlesInSuperParticleSystem(sParticleSystem,[](auto& particle){
      std::size_t globalId = getGlobalID( particle );
      //Randomly invalidate particle 3 and 7
      if (globalId==3 || globalId==7){
        particle.template setField<GENERAL,INVALID>( true );
      }
  });

  //Call size<>() with various conditions and evaluate
  std::size_t countAll = sParticleSystem.size<all_particles>();
  std::size_t countValid = sParticleSystem.size<valid_particles>();
  std::size_t countInvalid = sParticleSystem.size<invalid_particles>();
  std::size_t countValidCentres = sParticleSystem.size<valid_particle_centres>();
  std::size_t countValidSurfaces = sParticleSystem.size<valid_particle_surfaces>();

  //Test consistency
  EXPECT_EQ(countAll, countInvalid+countValid);
  EXPECT_EQ(countValid, countValidCentres+countValidSurfaces);

  //Test explicit
  EXPECT_EQ(countAll, 243);
  EXPECT_EQ(countValid, 189);
  EXPECT_EQ(countInvalid, 54);
  EXPECT_EQ(countValidCentres, 7);
  EXPECT_EQ(countValidSurfaces, 182);

}


TEST(TestParticleConditions, subgridConditions)
{
  //Settings
  typedef SubgridParticle3D PARTICLETYPE;
  constexpr unsigned D = PARTICLETYPE::d;
  PhysR<T,D> origin (0.0);
  const T baseSize = 1.0;
  const T delX = 0.1;
  //Set up environment
  Vector<T,D> extend( baseSize );
  IndicatorCuboid<T,D> cuboid( extend, origin );
  const int noOfCuboids = 27;
  CuboidDecomposition<T,D> cuboidDecomposition( cuboid, delX, noOfCuboids );
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,D> sGeometry( cuboidDecomposition, loadBalancer, 2 );

  //Create super particle system
  SuperParticleSystem<T,PARTICLETYPE> sParticleSystem(sGeometry);

  //Define particle properties
  T radius = .2*baseSize;
  std::size_t noP = 12;
  for (std::size_t iP=0; iP<noP; ++iP){
    //Define position
    Vector<T,3> position(.07*iP*baseSize);
    //Add Sphere
    addSubgridSphere3D( sParticleSystem, position, radius );
  }

  //Iterate over particles and randomly assign traits by globalId
  forParticlesInSuperParticleSystem(sParticleSystem,[](auto& particle){
      std::size_t globalId = getGlobalID( particle );
      //Randomly deactivate particle 1 and 8
      if (globalId==1 || globalId==5 || globalId==8){
        particle.template setField<DYNBEHAVIOUR,ACTIVE>( false );
      }
  });

  //Call size<>() with various conditions and evaluate
  std::size_t countAll = sParticleSystem.size<all_particles>();
  std::size_t countActive = sParticleSystem.size<active_particles>();

  //Test explicit
  EXPECT_EQ(countAll, 12);
  EXPECT_EQ(countActive, 9);
}
