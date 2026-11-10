/*  This file is part of the OpenLB library
*
*  Copyright (C) 2021 Jan E. Marquardt
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
#include <gtest/gtest.h>
#include "helper/unitTestHelper.h"

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::particles;

typedef double T;

// Test signed distance to particle
TEST(signedDistanceToParticle, testCircleCenter)
{
  using namespace descriptors;
  typedef ResolvedParticle2D PARTICLETYPE;
  T const radius = 1;
  PhysR<T,2> center(5.);
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  particles::creators::addResolvedCircle2D( particleSystem, center, radius, 0.5 );
  auto particle = particleSystem.get(0);
  particle.template setField<SURFACE,ANGLE>( 33. );

  T signedDistance;
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center);
  EXPECT_EQ(signedDistance, -1);
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center+PhysR<T,2>(0.5,0));
  EXPECT_EQ(signedDistance, -0.5);
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center+PhysR<T,2>(1.5,0));
  EXPECT_EQ(signedDistance, 0.5);
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center-PhysR<T,2>(0.5,0));
  EXPECT_EQ(signedDistance, -0.5);
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center-PhysR<T,2>(1.5,0));
  EXPECT_EQ(signedDistance, 0.5);
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center+PhysR<T,2>(0,0.5));
  EXPECT_EQ(signedDistance, -0.5);
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center+PhysR<T,2>(0,1.5));
  EXPECT_EQ(signedDistance, 0.5);
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center-PhysR<T,2>(0,0.5));
  EXPECT_EQ(signedDistance, -0.5);
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center-PhysR<T,2>(0,1.5));
  EXPECT_EQ(signedDistance, 0.5);
}

TEST(signedDistanceToParticle, testSphereCenter)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  T const radius = 1;
  PhysR<T,3> center(5.);
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  particles::creators::addResolvedSphere3D( particleSystem, center, radius, 0.5 );
  auto particle = particleSystem.get(0);
  particle.template setField<SURFACE,ANGLE>( Vector<T,3>(33.) );

  T signedDistance;
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center);
  EXPECT_EQ(signedDistance, -1);
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center+PhysR<T,3>(0.5,0,0));
  EXPECT_EQ(signedDistance, -0.5);
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center+PhysR<T,3>(1.5,0,0));
  EXPECT_EQ(signedDistance, 0.5);
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center-PhysR<T,3>(0.5,0,0));
  EXPECT_EQ(signedDistance, -0.5);
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center-PhysR<T,3>(1.5,0,0));
  EXPECT_EQ(signedDistance, 0.5);

  signedDistance = particles::resolved::signedDistanceToParticle(particle, center+PhysR<T,3>(0,0.5,0));
  EXPECT_EQ(signedDistance, -0.5);
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center+PhysR<T,3>(0,1.5,0));
  EXPECT_EQ(signedDistance, 0.5);
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center-PhysR<T,3>(0,0.5,0));
  EXPECT_EQ(signedDistance, -0.5);
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center-PhysR<T,3>(0,1.5,0));
  EXPECT_EQ(signedDistance, 0.5);

  signedDistance = particles::resolved::signedDistanceToParticle(particle, center+PhysR<T,3>(0,0,0.5));
  EXPECT_EQ(signedDistance, -0.5);
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center+PhysR<T,3>(0,0,1.5));
  EXPECT_EQ(signedDistance, 0.5);
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center-PhysR<T,3>(0,0,0.5));
  EXPECT_EQ(signedDistance, -0.5);
  signedDistance = particles::resolved::signedDistanceToParticle(particle, center-PhysR<T,3>(0,0,1.5));
  EXPECT_EQ(signedDistance, 0.5);
}

// Test distance to particle
TEST(distanceToParticle, testSphereCenter)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  T const radius = 1;
  PhysR<T,3> center(5.);
  Vector<short,3> directions[6] = {
    {-1, 0, 0}, {0,-1, 0},
    { 0, 0,-1}, {1, 0, 0},
    { 0, 1, 0}, {0, 0, 1}
  };
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  particles::creators::addResolvedSphere3D( particleSystem, center, radius, 0.5 );
  auto particle = particleSystem.get(0);

  for (unsigned i=0; i<6; ++i) {
    PhysR<T,3> direction(directions[i]);
    T distance = particles::resolved::distanceToParticle(particle, center, direction, 1e-16, 0.5);
    EXPECT_EQ(distance, radius);
  }
}

TEST(distanceToParticle, testSphereBoundary)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  T const radius = 1;
  PhysR<T,3> center(5.);
  PhysR<T,3> testPos(5.,6.,5.);
  Vector<short,3> directions[6] = {
    {-1, 0, 0}, {0,-1, 0},
    { 0, 0,-1}, {1, 0, 0},
    { 0, 1, 0}, {0, 0, 1}
  };
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  particles::creators::addResolvedSphere3D( particleSystem, center, radius, 0.5 );
  auto particle = particleSystem.get(0);

  for (unsigned i=0; i<6; ++i) {
    PhysR<T,3> direction(directions[i]);
    T distance = particles::resolved::distanceToParticle(particle, testPos, direction, 1e-16, 0.5);
    EXPECT_NEAR(distance, 0, 2e-8);
  }
}

TEST(distanceToParticle, testSphereOutside)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  T const radius = 1;
  PhysR<T,3> center(5.);
  PhysR<T,3> testPos(6.5,5.,5.);
  Vector<short,3> directions[6] = {
    {-1, 0, 0}, {0,-1, 0},
    { 0, 0,-1}, {1, 0, 0},
    { 0, 1, 0}, {0, 0, 1}
  };
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  particles::creators::addResolvedSphere3D( particleSystem, center, radius, 0.5 );
  auto particle = particleSystem.get(0);

  PhysR<T,3> direction(directions[0]);
  T distance = particles::resolved::distanceToParticle(particle, testPos, direction, 1e-16, 0.5);
  EXPECT_NEAR(distance, 0.5, 1e-16);
}


TEST(distanceToParticle, testCylinderCenter)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  PhysR<T,3> center(-5.);
  Vector<short,3> directions[6] = {
    {-1, 0, 0}, {0,-1, 0},
    { 0, 0,-1}, {1, 0, 0},
    { 0, 1, 0}, {0, 0, 1}
  };
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  T const radius = 2, height = 2*radius;
  particles::creators::addResolvedCylinder3D( particleSystem, center, Vector<T,3>(0.,1.,0.), height, radius, 0.);
  auto particle = particleSystem.get(0);

  for (unsigned i=0; i<6; ++i) {
    PhysR<T,3> direction(directions[i]);
    T distance = particles::resolved::distanceToParticle(particle, center, direction, 1e-16, 0.5);
    EXPECT_NEAR(distance, radius, 1e-15);
  }
}

TEST(distanceToParticle, testRotatedCylinderCenter)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  PhysR<T,3> center(-5.);
  Vector<T,3> angle(90., 0., 0.);
  Vector<short,3> directions[6] = {
    {-1, 0, 0}, {0,-1, 0},
    { 0, 0,-1}, {1, 0, 0},
    { 0, 1, 0}, {0, 0, 1}
  };
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  T const radius = 2, height = 2;
  particles::creators::addResolvedCylinder3D( particleSystem, center, Vector<T,3>(0.,1.,0.), height, radius, 0., 0., angle );
  auto particle = particleSystem.get(0);

  T distanceInd = 0;
  Vector<T,9> rotationMatrix = util::calculateRotationMatrix<T,3>( util::degreeToRadian(angle) );
  Vector<T,3> normal = util::executeRotation<T,3>(Vector<T,3>(0.,1.,0.), rotationMatrix, Vector<T,3>(0.));
  IndicatorCylinder3D<T> cyl(center, normal, radius, height);

  for (unsigned i=0; i<6; ++i) {
    PhysR<T,3> direction(directions[i]);
    T distance = particles::resolved::distanceToParticle(particle, center, direction, 1e-16, 0.5);
    cyl.distance(distanceInd, center, direction, 1e-16, 0.5);
    EXPECT_NEAR(distance, distanceInd, 1e-15);
  }
}
