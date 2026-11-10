/*  This file is part of the OpenLB library
*
*  Copyright (C) 2022 Nicolas Hafen, Martin Sadric
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
#include "helper/framework.h"

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::particles;
using namespace olb::particles::dynamics;
typedef double T;

::testing::Environment* const env = ::testing::AddGlobalTestEnvironment(new OlbInitTestEnvironment);

//#define VERBOSE_TESTS

///2D

TEST(TestParticleUtilities, test2Dstd)
{
  const int D = 2;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>> DESCRIPTOR;

  DynState<T,DESCRIPTOR> state;

#ifdef VERBOSE_TESTS
  std::cout << "Position=        " << state.position << std::endl;
#endif
}

TEST(TestParticleUtilities, test2Dpos)
{
  const int D = 2;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>> DESCRIPTOR;

  Vector<T,D> position(3.0, 6.0);
  DynState<T,DESCRIPTOR> state( position );

#ifdef VERBOSE_TESTS
  std::cout << "Position=        " << state.position << std::endl;
#endif
}

TEST(TestParticleUtilities, test2DstdPos)
{
  const int D = 2;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>,SURFACE_RESOLVED<D>> DESCRIPTOR;

  Vector<T,D> position(3.0, 6.0);
  DynState<T,DESCRIPTOR> state( position );

#ifdef VERBOSE_TESTS
  std::cout << "Position=        " << state.position << std::endl;
  std::cout << "Angle=           " << state.angle << std::endl;
#endif
}

TEST(TestParticleUtilities, test2DvalAngle)
{
  const int D = 2;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>,SURFACE_RESOLVED<D>> DESCRIPTOR;

  Vector<T,D> position(3.0, 6.0);
  Vector<T,utilities::dimensions::convert<D>::rotation> angle(0.5);
  DynState<T,DESCRIPTOR> state( position, angle );

#ifdef VERBOSE_TESTS
  std::cout << "Position=        " << state.position << std::endl;
  std::cout << "Angle=           " << state.angle << std::endl;
#endif
}

///3D

TEST(TestParticleUtilities, test3Dstd)
{
  const int D = 3;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>> DESCRIPTOR;

  DynState<T,DESCRIPTOR> state;

#ifdef VERBOSE_TESTS
  std::cout << "Position=        " << state.position << std::endl;
#endif
}

TEST(TestParticleUtilities, test3Dval)
{
  const int D = 3;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>> DESCRIPTOR;

  Vector<T,D> position(2.2, 5.6, 4.2);
  DynState<T,DESCRIPTOR> state( position );

#ifdef VERBOSE_TESTS
  std::cout << "Position=        " << state.position << std::endl;
#endif
}

TEST(TestParticleUtilities, test3DstdPos)
{
  const int D = 3;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>,SURFACE_RESOLVED<D>> DESCRIPTOR;

  Vector<T,D> position(1.4, 4.6, 5.4);
  DynState<T,DESCRIPTOR> state( position );

#ifdef VERBOSE_TESTS
  std::cout << "Position=        " << state.position << std::endl;
  std::cout << "Angle=           " << state.angle << std::endl;
#endif
}

TEST(TestParticleUtilities, test3DvalAngle)
{
  const int D = 3;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>,SURFACE_RESOLVED<D>> DESCRIPTOR;

  Vector<T,D> position(1.4, 4.6, 5.4);
  Vector<T,utilities::dimensions::convert<D>::rotation> angle(0.2, 0.03, 1.9);
  DynState<T,DESCRIPTOR> state( position, angle );

#ifdef VERBOSE_TESTS
  std::cout << "Position=        " << state.position << std::endl;
  std::cout << "Angle=           " << state.angle << std::endl;
#endif
}



TEST(TestPurgeParticles, testPurgeNonParallel){

  typedef PARTICLE_DESCRIPTOR<3, GENERAL_EXTENDABLE<3>,
          MOBILITY_VERLET_NO_ANGLE<3>, FORCING_SUBGRID<3>,
          PHYSPROPERTIES_SUBGRID<3>,
          DYNBEHAVIOUR_BASIC > PARTICLETYPE;

  //Create particle system
  int noP = 11;
  ParticleSystem<T,PARTICLETYPE> particleSystem(noP);

  //Invalidate particles
  particleSystem.get(1).setField<GENERAL,INVALID>(true);
  particleSystem.get(4).setField<GENERAL,INVALID>(true);
  particleSystem.get(5).setField<GENERAL,INVALID>(true);
  particleSystem.get(6).setField<GENERAL,INVALID>(true);

  //Evaluate number of particles (before purge)
  auto sizePre = particleSystem.size();
  EXPECT_EQ(sizePre,noP);

  //Purge invalid particles
  purgeInvalidParticles<T,PARTICLETYPE>( particleSystem );

  //Evaluate number of particles (after purge)
  auto sizePost = particleSystem.size();
  EXPECT_EQ(sizePost,7);
}

TEST(TestPurgeParticles, testPurgeConsiderIndexShift){

  typedef PARTICLE_DESCRIPTOR<3, GENERAL_EXTENDABLE<3>,
          MOBILITY_VERLET_NO_ANGLE<3>, FORCING_SUBGRID<3>,
          PHYSPROPERTIES_SUBGRID<3>,
          DYNBEHAVIOUR_BASIC > PARTICLETYPE;

  //Create particle system
  int noP = 20;
  ParticleSystem<T,PARTICLETYPE> particleSystem(noP);

  //Invalidate particles
  particleSystem.get(1).setField<GENERAL,INVALID>(true);
  particleSystem.get(2).setField<GENERAL,INVALID>(true);
  particleSystem.get(3).setField<GENERAL,INVALID>(true);
  particleSystem.get(4).setField<GENERAL,INVALID>(true);
  particleSystem.get(5).setField<GENERAL,INVALID>(true);
  particleSystem.get(6).setField<GENERAL,INVALID>(true);
  particleSystem.get(15).setField<GENERAL,INVALID>(true);
  particleSystem.get(15).setField<GENERAL,INVALID>(true);
  particleSystem.get(19).setField<GENERAL,INVALID>(true);

  //Evaluate number of particles (before purge)
  auto sizePre = particleSystem.size();
  EXPECT_EQ(sizePre,noP);

  //Purge invalid particles
  purgeInvalidParticles<T,PARTICLETYPE>( particleSystem );

  //Evaluate number of particles (after purge)
  auto sizePost = particleSystem.size();
  EXPECT_EQ(sizePost,12);
}

TEST(TestPurgeParticles, testPurgeParallel){

  typedef SubgridParticle3D PARTICLETYPE;

  //Create super particle system
  Vector<T,3> suGEOo (0.);
  Vector<T,3> suGEOe (1.);
  IndicatorCuboid3D<T> cuboid( suGEOe, suGEOo );
  const int noOfCuboids = 2;
  CuboidDecomposition3D<T> cuboidDecomposition( cuboid, 0.1, noOfCuboids );
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,3> sGeometry( cuboidDecomposition, loadBalancer, 2 );
  SuperParticleSystem<T,PARTICLETYPE> sParticleSystem(sGeometry);

  //Add particles
  auto& particleSystemA = *sParticleSystem.getBlockParticleSystems()[0];
  auto& containerA = particleSystemA.get();
  containerA.push_back();
  containerA.push_back();
  containerA.push_back();
  auto& particleSystemB = *sParticleSystem.getBlockParticleSystems()[1];
  auto& containerB = particleSystemB.get();
  containerB.push_back();
  containerB.push_back();

  //Invalidate particles
  particleSystemA.get(1).setField<GENERAL,INVALID>(true);
  particleSystemA.get(2).setField<GENERAL,INVALID>(true);
  particleSystemB.get(0).setField<GENERAL,INVALID>(true);

  //Evaluate number of particles (before purge)
  auto sizePre = sParticleSystem.size();
  EXPECT_EQ(sizePre,5);

  //Purge invalid particles
  purgeInvalidParticles<T,PARTICLETYPE>( sParticleSystem );

  //Evaluate number of particles (after purge)
  auto sizePost = sParticleSystem.size();
  EXPECT_EQ(sizePost,2);
}



TEST(TestParticleUtilites, testDiscretePointsOnSphericalHull2D)
{
  Vector<T, 2> position(1,1);
  T radius(1);

  std::array<Vector<T, 2>, 8> pointsOnSphericalHull = {
    Vector<T, 2>(0.292893, 0.292893),
    Vector<T, 2>(1,0),  // right
    Vector<T, 2>(1.70711, 0.292893),
    Vector<T, 2>(0,1),  // back
    Vector<T, 2>(2,1),   // front
    Vector<T, 2>(0.292893, 1.70711),
    Vector<T, 2>(1, 2),   // left
    Vector<T, 2>(1.70711, 1.70711)
  };

  std::array<Vector<T, 2>, 8> actualPoints = discretePointsOnSphericalHull(position, radius);

  for (int i = 0; i < 8; ++i) {
    EXPECT_NEAR(pointsOnSphericalHull[i][0], actualPoints[i][0], 1e-5);
    EXPECT_NEAR(pointsOnSphericalHull[i][1], actualPoints[i][1], 1e-5);
  }
}

TEST(TestParticleUtilities, testDiscretePointsOnSphericalHull3D)
{
  Vector<T, 3> position(1, 1, 1);
  T radius(1);

  std::array<Vector<T, 3>, 26> pointsOnSphericalHull = {
    Vector<T, 3>(0.5, 0.5, 0.5),
    Vector<T, 3>(1, 0.292893, 0.292893),
    Vector<T, 3>(1.5, 0.5, 0.5),
    Vector<T, 3>(0.292893, 1, 0.292893),
    Vector<T, 3>(1, 1, 0),
    Vector<T, 3>(1.70711, 1, 0.292893),
    Vector<T, 3>(0.5, 1.5, 0.5),
    Vector<T, 3>(1, 1.70711, 0.292893),
    Vector<T, 3>(1.5, 1.5, 0.5),
    Vector<T, 3>(0.292893, 0.292893, 1),
    Vector<T, 3>(1, 0, 1),
    Vector<T, 3>(1.70711, 0.292893, 1),
    Vector<T, 3>(0, 1, 1),
    Vector<T, 3>(2, 1, 1),
    Vector<T, 3>(0.292893, 1.70711, 1),
    Vector<T, 3>(1, 2, 1),
    Vector<T, 3>(1.70711, 1.70711, 1),
    Vector<T, 3>(0.5, 0.5, 1.5),
    Vector<T, 3>(1, 0.292893, 1.70711),
    Vector<T, 3>(1.5, 0.5, 1.5),
    Vector<T, 3>(0.292893, 1, 1.70711),
    Vector<T, 3>(1, 1, 2),
    Vector<T, 3>(1.70711, 1, 1.70711),
    Vector<T, 3>(0.5, 1.5, 1.5),
    Vector<T, 3>(1, 1.70711, 1.70711),
    Vector<T, 3>(1.5, 1.5, 1.5)
  };

  std::array<Vector<T, 3>, 26> actualPoints = discretePointsOnSphericalHull(position, radius);


  for (int i = 0; i < 26; ++i) {
    EXPECT_NEAR(pointsOnSphericalHull[i][0], actualPoints[i][0], 1e-5);
    EXPECT_NEAR(pointsOnSphericalHull[i][1], actualPoints[i][1], 1e-5);
    EXPECT_NEAR(pointsOnSphericalHull[i][2], actualPoints[i][2], 1e-5);
  }
}

