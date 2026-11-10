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

#include <type_traits>

#include <gtest/gtest.h>
#include "helper/unitTestHelper.h"

#include <olb.h>


using namespace olb;
using namespace olb::descriptors;
using namespace olb::particles;
using namespace olb::particles::dynamics;
using SubgridParticle3DsingleCuboid  = PARTICLE_DESCRIPTOR<3, GENERAL_TMP<3>,
      MOBILITY_VERLET_NO_ANGLE<3>, FORCING_SUBGRID<3>,
      PHYSPROPERTIES_SUBGRID<3>,
      DYNBEHAVIOUR_BASIC>;
typedef double T;
typedef SubgridParticle3DsingleCuboid PARTICLETYPE;

//#define VERBOSE_TESTS




//Print order including activity and dummy dynamics
void printOrder( std::string heading,
                 ParticleSystem<T,PARTICLETYPE>& particleSystem)
{
  std::cout << heading << ":" << std::endl;
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
    auto particle = particleSystem.get(iP);
    bool active = particle.template getField<DYNBEHAVIOUR,ACTIVE>();
    auto dynamicsName = particleSystem.getDynamics(iP)->getName();
    std::cout << " iP=" << iP << ", active=" << active
              << ", dynamics=" << dynamicsName << std::endl;
  }
}




TEST(TestParticleSorting, testRandom)
{
  std::vector<bool> activity = {
    false,
    false,
    true,
    false,
    true,
    true,
    false,
    true,
    true,
    true,
    true
  };

  //Create particle system
  ParticleSystem<T,PARTICLETYPE> particleSystem(activity.size());

  //Create particle dynamics
  UnitConverter<T,D3Q19<>> converter( 1., 1., 1., 1., 1., 1.);
  particleSystem.defineDynamics<
    VerletParticleDynamics<T,PARTICLETYPE>>();

  //Set upt particle and count occurences of A and B for later reference
  int noApre = 0;
  int noBpre = 0;
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
    auto particle = particleSystem.get(iP);
    auto active = activity[iP];
    //Setting activity
    particle.setField<DYNBEHAVIOUR,ACTIVE>(active);
    if (active) {
      ++noBpre;
    }
    else {
      ++noApre;
    }
  }

#ifdef VERBOSE_TESTS
  printOrder("Original", particleSystem);
#endif

  //Sort particle system
  auto iPfirstB = sorting::partitionParticleSystem<
                  T,PARTICLETYPE,DYNBEHAVIOUR,ACTIVE>(particleSystem );

#ifdef VERBOSE_TESTS
  printOrder("Sorted", particleSystem);
  std::cout << "iPfirstB=" << iPfirstB << std::endl;
#endif


  int noApost = 0;
  int noBpost = 0;
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
    auto particle = particleSystem.get(iP);
    auto active = activity[iP];
    if (active) {
      ++noBpost;
    }
    else {
      ++noApost;
    }
  }

  //Test general
  EXPECT_EQ(noApre,noApost);
  EXPECT_EQ(noBpre,noBpost);

  //Test setting specifc
  EXPECT_EQ(iPfirstB,4);
}


TEST(TestParticleSorting, testFirstActive)
{
  std::vector<bool> activity = {
    true,
    false,
    true,
    false,
    true,
    false,
    true,
    true
  };

  //Create particle system
  ParticleSystem<T,PARTICLETYPE> particleSystem(activity.size());

  //Create particle dynamics
  UnitConverter<T,D3Q19<>> converter( 1., 1., 1., 1., 1., 1.);
  particleSystem.defineDynamics<
    VerletParticleDynamics<T,PARTICLETYPE>>();

  //Set upt particle and count occurences of A and B for later reference
  int noApre = 0;
  int noBpre = 0;
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
    auto particle = particleSystem.get(iP);
    auto active = activity[iP];
    //Setting activity
    particle.setField<DYNBEHAVIOUR,ACTIVE>(active);
    if (active) {
      ++noBpre;
    }
    else {
      ++noApre;
    }
  }

#ifdef VERBOSE_TESTS
  printOrder("Original", particleSystem);
#endif

  //Sort particle system
  auto iPfirstB = sorting::partitionParticleSystem<
                  T,PARTICLETYPE,DYNBEHAVIOUR,ACTIVE>(particleSystem );

#ifdef VERBOSE_TESTS
  printOrder("Sorted", particleSystem);
  std::cout << "iPfirstB=" << iPfirstB << std::endl;
#endif


  int noApost = 0;
  int noBpost = 0;
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
    auto particle = particleSystem.get(iP);
    auto active = activity[iP];
    if (active) {
      ++noBpost;
    }
    else {
      ++noApost;
    }
  }

  //Test general
  EXPECT_EQ(noApre,noApost);
  EXPECT_EQ(noBpre,noBpost);

  //Test setting specifc
  EXPECT_EQ(iPfirstB,3);
}


TEST(TestParticleSorting, testAllActive)
{
  std::vector<bool> activity = {
    true,
    true,
    true,
    true,
    true,
  };

  //Create particle system
  ParticleSystem<T,PARTICLETYPE> particleSystem(activity.size());

  //Create particle dynamics
  UnitConverter<T,D3Q19<>> converter( 1., 1., 1., 1., 1., 1.);
  particleSystem.defineDynamics<
    VerletParticleDynamics<T,PARTICLETYPE>>();

  //Set upt particle and count occurences of A and B for later reference
  int noApre = 0;
  int noBpre = 0;
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
    auto particle = particleSystem.get(iP);
    auto active = activity[iP];
    //Setting activity
    particle.setField<DYNBEHAVIOUR,ACTIVE>(active);
    if (active) {
      ++noBpre;
    }
    else {
      ++noApre;
    }
  }

#ifdef VERBOSE_TESTS
  printOrder("Original", particleSystem);
#endif

  //Sort particle system
  auto iPfirstB = sorting::partitionParticleSystem<
                  T,PARTICLETYPE,DYNBEHAVIOUR,ACTIVE>(particleSystem );

#ifdef VERBOSE_TESTS
  printOrder("Sorted", particleSystem);
  std::cout << "iPfirstB=" << iPfirstB << std::endl;
#endif


  int noApost = 0;
  int noBpost = 0;
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
    auto particle = particleSystem.get(iP);
    auto active = activity[iP];
    if (active) {
      ++noBpost;
    }
    else {
      ++noApost;
    }
  }

  //Test general
  EXPECT_EQ(noApre,noApost);
  EXPECT_EQ(noBpre,noBpost);

  //Test setting specifc
  EXPECT_EQ(iPfirstB,0);
}


TEST(TestParticleSorting, testAllInactive)
{
  std::vector<bool> activity = {
    false,
    false,
    false,
    false,
    false,
    false,
    false,
    false,
  };

  //Create particle system
  ParticleSystem<T,PARTICLETYPE> particleSystem(activity.size());

  //Create particle dynamics
  UnitConverter<T,D3Q19<>> converter( 1., 1., 1., 1., 1., 1.);
  particleSystem.defineDynamics<
    VerletParticleDynamics<T,PARTICLETYPE>>();

  //Set upt particle and count occurences of A and B for later reference
  int noApre = 0;
  int noBpre = 0;
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
    auto particle = particleSystem.get(iP);
    auto active = activity[iP];
    //Setting activity
    particle.setField<DYNBEHAVIOUR,ACTIVE>(active);
    if (active) {
      ++noBpre;
    }
    else {
      ++noApre;
    }
  }

#ifdef VERBOSE_TESTS
  printOrder("Original", particleSystem);
#endif

  //Sort particle system
  auto iPfirstB = sorting::partitionParticleSystem<
                  T,PARTICLETYPE,DYNBEHAVIOUR,ACTIVE>(particleSystem );

#ifdef VERBOSE_TESTS
  printOrder("Sorted", particleSystem);
  std::cout << "iPfirstB=" << iPfirstB << std::endl;
#endif


  int noApost = 0;
  int noBpost = 0;
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
    auto particle = particleSystem.get(iP);
    auto active = activity[iP];
    if (active) {
      ++noBpost;
    }
    else {
      ++noApost;
    }
  }

  //Test general
  EXPECT_EQ(noApre,noApost);
  EXPECT_EQ(noBpre,noBpost);

  //Test setting specifc
  EXPECT_EQ(iPfirstB,8);
}

