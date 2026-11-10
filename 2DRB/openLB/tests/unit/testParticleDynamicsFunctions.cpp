/*  This file is part of the OpenLB library
*
*  Copyright (C) 2017 Markus Mohrhard, Davide Dapelo
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

#include "helper/framework.h"
#include "helper/unitTestHelper.h"
#include <olb.h>
#include <gtest/gtest.h>
#include <iostream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::particles;
using namespace olb::particles;

typedef double T;

template <typename T, typename PARTICLETYPE>
void setMobility(
    Particle<T, PARTICLETYPE>& particle, const Vector<T, PARTICLETYPE::d>& vel,
    const Vector<T, utilities::dimensions::convert<PARTICLETYPE::d>::rotation>&
        angVel)
{
  particle.template setField<MOBILITY, VELOCITY>(vel);
  particle.template setField<MOBILITY, ANG_VELOCITY>(
      utilities::dimensions::convert<PARTICLETYPE::d>::serialize_rotation(
          angVel));
}

TEST(LocalVelocity2D, restingCircle)
{
  typedef ResolvedCircle2D PARTICLETYPE;
  constexpr unsigned       D      = PARTICLETYPE::d;
  constexpr unsigned       Drot   = utilities::dimensions::convert<D>::rotation;
  T const                  radius = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  particles::creators::addResolvedCircle2D(particleSystem, Vector<T, D>(0.),
                                           radius, 0.);
  Particle<T, PARTICLETYPE> particle = particleSystem.get(0);
  setMobility(particle, Vector<T, D>(0.), Vector<T, Drot>(0.));

  PhysR<T, D>  pos(0.);
  Vector<T, D> expectedVel(0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedVel[iDim]);
  }

  pos         = PhysR<T, D>(radius);
  expectedVel = Vector<T, D>(0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedVel[iDim]);
  }

  pos         = PhysR<T, D>(-radius);
  expectedVel = Vector<T, D>(0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedVel[iDim]);
  }
}

TEST(LocalVelocity2D, movingCircle1)
{
  typedef ResolvedCircle2D PARTICLETYPE;
  constexpr unsigned       D      = PARTICLETYPE::d;
  constexpr unsigned       Drot   = utilities::dimensions::convert<D>::rotation;
  T const                  radius = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  particles::creators::addResolvedCircle2D(particleSystem, Vector<T, D>(0.),
                                           radius, 0.);
  Particle<T, PARTICLETYPE> particle = particleSystem.get(0);
  setMobility(particle, Vector<T, D>(0.1), Vector<T, Drot>(0.));

  PhysR<T, D>  pos(0.);
  Vector<T, D> expectedVel(0.1);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedVel[iDim]);
  }

  pos         = PhysR<T, D>(radius);
  expectedVel = Vector<T, D>(0.1);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedVel[iDim]);
  }

  pos         = PhysR<T, D>(-radius);
  expectedVel = Vector<T, D>(0.1);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedVel[iDim]);
  }
}

TEST(LocalVelocity2D, movingCircle2)
{
  typedef ResolvedCircle2D PARTICLETYPE;
  constexpr unsigned       D      = PARTICLETYPE::d;
  constexpr unsigned       Drot   = utilities::dimensions::convert<D>::rotation;
  T const                  radius = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  particles::creators::addResolvedCircle2D(particleSystem, Vector<T, D>(0.),
                                           radius, 0.);
  Particle<T, PARTICLETYPE> particle = particleSystem.get(0);
  setMobility(particle, Vector<T, D>(-0.137), Vector<T, Drot>(0.));

  PhysR<T, D>  pos(0.);
  Vector<T, D> expectedVel(-0.137);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedVel[iDim]);
  }

  pos         = PhysR<T, D>(radius);
  expectedVel = Vector<T, D>(-0.137);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedVel[iDim]);
  }

  pos         = PhysR<T, D>(-radius);
  expectedVel = Vector<T, D>(-0.137);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedVel[iDim]);
  }
}

TEST(LocalVelocity2D, movingCircle3)
{
  typedef ResolvedCircle2D PARTICLETYPE;
  constexpr unsigned       D      = PARTICLETYPE::d;
  constexpr unsigned       Drot   = utilities::dimensions::convert<D>::rotation;
  T const                  radius = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  particles::creators::addResolvedCircle2D(particleSystem, Vector<T, D>(0.),
                                           radius, 0.);
  Particle<T, PARTICLETYPE> particle = particleSystem.get(0);
  setMobility(particle, Vector<T, D>(-0.1, 0.2), Vector<T, Drot>(0.));

  PhysR<T, D>  pos(0.);
  Vector<T, D> expectedVel(-0.1, 0.2);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedVel[iDim]);
  }

  pos         = PhysR<T, D>(radius);
  expectedVel = Vector<T, D>(-0.1, 0.2);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedVel[iDim]);
  }

  pos         = PhysR<T, D>(-radius);
  expectedVel = Vector<T, D>(-0.1, 0.2);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedVel[iDim]);
  }
}

TEST(LocalVelocity2D, rotatingCircle1)
{
  typedef ResolvedCircle2D PARTICLETYPE;
  constexpr unsigned       D      = PARTICLETYPE::d;
  constexpr unsigned       Drot   = utilities::dimensions::convert<D>::rotation;
  T const                  radius = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  particles::creators::addResolvedCircle2D(particleSystem, Vector<T, D>(0.),
                                           radius, 0.);
  Particle<T, PARTICLETYPE> particle = particleSystem.get(0);
  setMobility(particle, Vector<T, D>(0.), Vector<T, Drot>(2.));

  PhysR<T, D>  pos(0.);
  Vector<T, D> expectedVel(0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedVel[iDim]);
  }

  pos         = PhysR<T, D>(radius, 0.);
  expectedVel = Vector<T, D>(0., 0.2);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedVel[iDim]);
  }

  pos         = PhysR<T, D>(0., radius);
  expectedVel = Vector<T, D>(-0.2, 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedVel[iDim]);
  }

  pos         = PhysR<T, D>(-radius, 0.);
  expectedVel = Vector<T, D>(0., -0.2);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedVel[iDim]);
  }

  pos         = PhysR<T, D>(0., -radius);
  expectedVel = Vector<T, D>(0.2, 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedVel[iDim]);
  }
}

TEST(LocalVelocity2D, rotatingCircle2)
{
  typedef ResolvedCircle2D PARTICLETYPE;
  constexpr unsigned       D      = PARTICLETYPE::d;
  constexpr unsigned       Drot   = utilities::dimensions::convert<D>::rotation;
  T const                  radius = 0.2;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  particles::creators::addResolvedCircle2D(particleSystem, Vector<T, D>(0.),
                                           radius, 0.);
  Particle<T, PARTICLETYPE> particle = particleSystem.get(0);
  setMobility(particle, Vector<T, D>(0.), Vector<T, Drot>(-1.5));

  PhysR<T, D>  pos(0.);
  Vector<T, D> expectedVel(0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedVel[iDim]);
  }

  pos         = PhysR<T, D>(radius, 0.);
  expectedVel = Vector<T, D>(0., -0.3);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_NEAR(
        particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
        expectedVel[iDim], 1e-16);
  }

  pos         = PhysR<T, D>(0., radius);
  expectedVel = Vector<T, D>(0.3, 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_NEAR(
        particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
        expectedVel[iDim], 1e-16);
  }

  pos         = PhysR<T, D>(-radius, 0.);
  expectedVel = Vector<T, D>(0., 0.3);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_NEAR(
        particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
        expectedVel[iDim], 1e-16);
  }

  pos         = PhysR<T, D>(0., -radius);
  expectedVel = Vector<T, D>(-0.3, 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_NEAR(
        particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
        expectedVel[iDim], 1e-16);
  }
}

TEST(LocalVelocity3D, restingSphere)
{
  typedef ResolvedSphere3D PARTICLETYPE;
  constexpr unsigned       D      = PARTICLETYPE::d;
  constexpr unsigned       Drot   = utilities::dimensions::convert<D>::rotation;
  T const                  radius = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  particles::creators::addResolvedSphere3D(particleSystem, Vector<T, D>(0.),
                                           radius, 0.);
  Particle<T, PARTICLETYPE> particle = particleSystem.get(0);
  setMobility(particle, Vector<T, D>(0.), Vector<T, Drot>(0.));

  PhysR<T, D>  pos(0.);
  Vector<T, D> expectedLocalVel(0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(radius);
  expectedLocalVel = Vector<T, D>(0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(-radius);
  expectedLocalVel = Vector<T, D>(0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }
}

TEST(LocalVelocity3D, movingSphere1)
{
  typedef ResolvedSphere3D PARTICLETYPE;
  constexpr unsigned       D      = PARTICLETYPE::d;
  constexpr unsigned       Drot   = utilities::dimensions::convert<D>::rotation;
  T const                  radius = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  particles::creators::addResolvedSphere3D(particleSystem, Vector<T, D>(0.),
                                           radius, 0.);
  Particle<T, PARTICLETYPE> particle = particleSystem.get(0);
  setMobility(particle, Vector<T, D>(0.1), Vector<T, Drot>(0.));

  PhysR<T, D>  pos(0.);
  Vector<T, D> expectedLocalVel(0.1);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(radius);
  expectedLocalVel = Vector<T, D>(0.1);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(-radius);
  expectedLocalVel = Vector<T, D>(0.1);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }
}

TEST(LocalVelocity3D, movingSphere2)
{
  typedef ResolvedSphere3D PARTICLETYPE;
  constexpr unsigned       D      = PARTICLETYPE::d;
  constexpr unsigned       Drot   = utilities::dimensions::convert<D>::rotation;
  T const                  radius = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  particles::creators::addResolvedSphere3D(particleSystem, Vector<T, D>(0.),
                                           radius, 0.);
  Particle<T, PARTICLETYPE> particle = particleSystem.get(0);
  setMobility(particle, Vector<T, D>(-0.11), Vector<T, Drot>(0.));

  PhysR<T, D>  pos(0.);
  Vector<T, D> expectedLocalVel(-0.11);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(radius);
  expectedLocalVel = Vector<T, D>(-0.11);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(-radius);
  expectedLocalVel = Vector<T, D>(-0.11);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }
}

TEST(LocalVelocity3D, movingSphere3)
{
  typedef ResolvedSphere3D PARTICLETYPE;
  constexpr unsigned       D      = PARTICLETYPE::d;
  constexpr unsigned       Drot   = utilities::dimensions::convert<D>::rotation;
  T const                  radius = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  particles::creators::addResolvedSphere3D(particleSystem, Vector<T, D>(0.),
                                           radius, 0.);
  Particle<T, PARTICLETYPE> particle = particleSystem.get(0);
  setMobility(particle, Vector<T, D>(-0.11, 0.22, -0.33), Vector<T, Drot>(0.));

  PhysR<T, D>  pos(0.);
  Vector<T, D> expectedLocalVel(-0.11, 0.22, -0.33);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(radius);
  expectedLocalVel = Vector<T, D>(-0.11, 0.22, -0.33);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(-radius);
  expectedLocalVel = Vector<T, D>(-0.11, 0.22, -0.33);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }
}

TEST(LocalVelocity3D, rotatingSphere1)
{
  typedef ResolvedSphere3D PARTICLETYPE;
  constexpr unsigned       D      = PARTICLETYPE::d;
  constexpr unsigned       Drot   = utilities::dimensions::convert<D>::rotation;
  T const                  radius = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  particles::creators::addResolvedSphere3D(particleSystem, Vector<T, D>(0.),
                                           radius, 0.);
  Particle<T, PARTICLETYPE> particle = particleSystem.get(0);
  setMobility(particle, Vector<T, D>(0.), Vector<T, Drot>(0., 0., 1.));

  PhysR<T, D>  pos(0.);
  Vector<T, D> expectedLocalVel(0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(radius, 0., 0.);
  expectedLocalVel = Vector<T, D>(0.0, 0.1, 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(-radius, 0., 0.);
  expectedLocalVel = Vector<T, D>(0., -0.1, 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(0., radius, 0.);
  expectedLocalVel = Vector<T, D>(-0.1, 0., 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(0., -radius, 0.);
  expectedLocalVel = Vector<T, D>(0.1, 0., 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }
}

TEST(LocalVelocity3D, rotatingSphere2)
{
  typedef ResolvedSphere3D PARTICLETYPE;
  constexpr unsigned       D      = PARTICLETYPE::d;
  constexpr unsigned       Drot   = utilities::dimensions::convert<D>::rotation;
  T const                  radius = 0.2;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  particles::creators::addResolvedSphere3D(particleSystem, Vector<T, D>(0.),
                                           radius, 0.);
  Particle<T, PARTICLETYPE> particle = particleSystem.get(0);
  setMobility(particle, Vector<T, D>(0.), Vector<T, Drot>(0., 0., -2.5));

  PhysR<T, D>  pos(0.);
  Vector<T, D> expectedLocalVel(0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(radius, 0., 0.);
  expectedLocalVel = Vector<T, D>(0.0, -0.5, 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(-radius, 0., 0.);
  expectedLocalVel = Vector<T, D>(0., 0.5, 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(0., radius, 0.);
  expectedLocalVel = Vector<T, D>(0.5, 0., 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(0., -radius, 0.);
  expectedLocalVel = Vector<T, D>(-0.5, 0., 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }
}

TEST(LocalVelocity3D, rotatingSphere3)
{
  typedef ResolvedSphere3D PARTICLETYPE;
  constexpr unsigned       D      = PARTICLETYPE::d;
  constexpr unsigned       Drot   = utilities::dimensions::convert<D>::rotation;
  T const                  radius = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  particles::creators::addResolvedSphere3D(particleSystem, Vector<T, D>(0.),
                                           radius, 0.);
  Particle<T, PARTICLETYPE> particle = particleSystem.get(0);
  setMobility(particle, Vector<T, D>(0.), Vector<T, Drot>(1., 0., 0.));

  PhysR<T, D>  pos(0.);
  Vector<T, D> expectedLocalVel(0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(0., radius, 0.);
  expectedLocalVel = Vector<T, D>(0., 0., 0.1);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(0., -radius, 0.);
  expectedLocalVel = Vector<T, D>(0., 0., -0.1);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(0., 0., radius);
  expectedLocalVel = Vector<T, D>(0., -0.1, 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(0., 0., -radius);
  expectedLocalVel = Vector<T, D>(0., 0.1, 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }
}

TEST(LocalVelocity3D, rotatingSphere4)
{
  typedef ResolvedSphere3D PARTICLETYPE;
  constexpr unsigned       D      = PARTICLETYPE::d;
  constexpr unsigned       Drot   = utilities::dimensions::convert<D>::rotation;
  T const                  radius = 0.2;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  particles::creators::addResolvedSphere3D(particleSystem, Vector<T, D>(0.),
                                           radius, 0.);
  Particle<T, PARTICLETYPE> particle = particleSystem.get(0);
  setMobility(particle, Vector<T, D>(0.), Vector<T, Drot>(-2.5, 0., 0.));

  PhysR<T, D>  pos(0.);
  Vector<T, D> expectedLocalVel(0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(0., radius, 0.);
  expectedLocalVel = Vector<T, D>(0., 0., -0.5);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(0., -radius, 0.);
  expectedLocalVel = Vector<T, D>(0., 0., 0.5);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(0., 0., radius);
  expectedLocalVel = Vector<T, D>(0., 0.5, 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(0., 0., -radius);
  expectedLocalVel = Vector<T, D>(0., -0.5, 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }
}

TEST(LocalVelocity3D, rotatingSphere5)
{
  typedef ResolvedSphere3D PARTICLETYPE;
  constexpr unsigned       D      = PARTICLETYPE::d;
  constexpr unsigned       Drot   = utilities::dimensions::convert<D>::rotation;
  T const                  radius = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  particles::creators::addResolvedSphere3D(particleSystem, Vector<T, D>(0.),
                                           radius, 0.);
  Particle<T, PARTICLETYPE> particle = particleSystem.get(0);
  setMobility(particle, Vector<T, D>(0.), Vector<T, Drot>(0., 1., 0.));

  PhysR<T, D>  pos(0.);
  Vector<T, D> expectedLocalVel(0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(radius, 0., 0.);
  expectedLocalVel = Vector<T, D>(0., 0., -0.1);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(-radius, 0., 0.);
  expectedLocalVel = Vector<T, D>(0., 0., 0.1);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(0., 0., radius);
  expectedLocalVel = Vector<T, D>(0.1, 0., 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(0., 0., -radius);
  expectedLocalVel = Vector<T, D>(-0.1, 0., 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }
}

TEST(LocalVelocity3D, rotatingSphere6)
{
  typedef ResolvedSphere3D PARTICLETYPE;
  constexpr unsigned       D      = PARTICLETYPE::d;
  constexpr unsigned       Drot   = utilities::dimensions::convert<D>::rotation;
  T const                  radius = 0.2;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  particles::creators::addResolvedSphere3D(particleSystem, Vector<T, D>(0.),
                                           radius, 0.);
  Particle<T, PARTICLETYPE> particle = particleSystem.get(0);
  setMobility(particle, Vector<T, D>(0.), Vector<T, Drot>(0., -2.5, 0.));

  PhysR<T, D>  pos(0.);
  Vector<T, D> expectedLocalVel(0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(radius, 0., 0.);
  expectedLocalVel = Vector<T, D>(0., 0., 0.5);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(-radius, 0., 0.);
  expectedLocalVel = Vector<T, D>(0., 0., -0.5);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(0., 0., radius);
  expectedLocalVel = Vector<T, D>(-0.5, 0., 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }

  pos              = PhysR<T, D>(0., 0., -radius);
  expectedLocalVel = Vector<T, D>(0.5, 0., 0.);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(particles::dynamics::calculateLocalVelocity(particle, pos)[iDim],
              expectedLocalVel[iDim]);
  }
}

// Torque calculation tests
TEST(TorqueFromForce2D, ZeroForce)
{
  constexpr unsigned D = 2;
  const Vector<T, D> force(T {0.});
  const Vector<T, D> lever(T {1.});

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], T {0}, 1e-16);
}

TEST(TorqueFromForce2D, ZeroLever)
{
  constexpr unsigned D = 2;
  const Vector<T, D> force(T {1.});
  const Vector<T, D> lever(T {0.});

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], T {0}, 1e-16);
}

TEST(TorqueFromForce2D, ForceAndLeverParallel1)
{
  constexpr unsigned D = 2;
  const Vector<T, D> force(T {1.});
  const Vector<T, D> lever(T {1.});

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], T {0}, 1e-16);
}

TEST(TorqueFromForce2D, ForceAndLeverParallel2)
{
  constexpr unsigned D = 2;
  const Vector<T, D> force(1., 2.);
  const Vector<T, D> lever(1., 2.);

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], T {0}, 1e-16);
}

TEST(TorqueFromForce2D, ForceAndLeverParallel3)
{
  constexpr unsigned D = 2;
  const Vector<T, D> force(-3.25, -2.5);
  const Vector<T, D> lever(-3.25, -2.5);

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], T {0}, 1e-16);
}

TEST(TorqueFromForce2D, ForceAndLeverNotParallel1)
{
  constexpr unsigned D = 2;
  const Vector<T, D> force(0., 2.5);
  const Vector<T, D> lever(1.4, 0.);

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], 3.5, 1e-16);
}

TEST(TorqueFromForce2D, ForceAndLeverNotParallel2)
{
  constexpr unsigned D = 2;
  const Vector<T, D> force(0., -2.5);
  const Vector<T, D> lever(1.4, 0.);

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], -3.5, 1e-16);
}

TEST(TorqueFromForce2D, ForceAndLeverNotParallel3)
{
  constexpr unsigned D = 2;
  const Vector<T, D> force(0., 2.5);
  const Vector<T, D> lever(-1.4, 0.);

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], -3.5, 1e-16);
}

TEST(TorqueFromForce2D, ForceAndLeverNotParallel4)
{
  constexpr unsigned D = 2;
  const Vector<T, D> force(0., -2.5);
  const Vector<T, D> lever(-1.4, 0.);

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], 3.5, 1e-16);
}

TEST(TorqueFromForce2D, ForceAndLeverNotParallel5)
{
  constexpr unsigned D = 2;
  const Vector<T, D> force(2.5, 0.);
  const Vector<T, D> lever(0., 1.4);

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], -3.5, 1e-16);
}

TEST(TorqueFromForce2D, ForceAndLeverNotParallel6)
{
  constexpr unsigned D = 2;
  const Vector<T, D> force(-2.5, 0.);
  const Vector<T, D> lever(0., 1.4);

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], 3.5, 1e-16);
}

TEST(TorqueFromForce2D, ForceAndLeverNotParallel7)
{
  constexpr unsigned D = 2;
  const Vector<T, D> force(2.5, 0.);
  const Vector<T, D> lever(0., -1.4);

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], 3.5, 1e-16);
}

TEST(TorqueFromForce2D, ForceAndLeverNotParallel8)
{
  constexpr unsigned D = 2;
  const Vector<T, D> force(-2.5, 0.);
  const Vector<T, D> lever(0., -1.4);

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], -3.5, 1e-16);
}

TEST(TorqueFromForce3D, ZeroForce)
{
  constexpr unsigned D = 3;
  const Vector<T, D> force(T {0.});
  const Vector<T, D> lever(T {1.});

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], T {0}, 1e-16);
  EXPECT_NEAR(torque[1], T {0}, 1e-16);
  EXPECT_NEAR(torque[2], T {0}, 1e-16);
}

TEST(TorqueFromForce3D, ZeroLever)
{
  constexpr unsigned D = 3;
  const Vector<T, D> force(T {1.});
  const Vector<T, D> lever(T {0.});

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], T {0}, 1e-16);
  EXPECT_NEAR(torque[1], T {0}, 1e-16);
  EXPECT_NEAR(torque[2], T {0}, 1e-16);
}

TEST(TorqueFromForce3D, ForceAndLeverParallel1)
{
  constexpr unsigned D = 3;
  const Vector<T, D> force(T {1.});
  const Vector<T, D> lever(T {1.});

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], T {0}, 1e-16);
  EXPECT_NEAR(torque[1], T {0}, 1e-16);
  EXPECT_NEAR(torque[2], T {0}, 1e-16);
}

TEST(TorqueFromForce3D, ForceAndLeverParallel2)
{
  constexpr unsigned D = 3;
  const Vector<T, D> force(1., 2., 3.);
  const Vector<T, D> lever(1., 2., 3.);

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], T {0}, 1e-16);
  EXPECT_NEAR(torque[1], T {0}, 1e-16);
  EXPECT_NEAR(torque[2], T {0}, 1e-16);
}

TEST(TorqueFromForce3D, ForceAndLeverParallel3)
{
  constexpr unsigned D = 3;
  const Vector<T, D> force(-3.25, -2.5, -5.65);
  const Vector<T, D> lever(-3.25, -2.5, -5.65);

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], T {0}, 1e-15);
  EXPECT_NEAR(torque[1], T {0}, 1e-15);
  EXPECT_NEAR(torque[2], T {0}, 1e-15);
}

TEST(TorqueFromForce3D, ForceAndLeverNotParallel1)
{
  constexpr unsigned D = 3;
  const Vector<T, D> force(1.4, 0.2, 8.5);
  const Vector<T, D> lever(4.8, 2.5, 3.1);

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], lever[1] * force[2] - lever[2] * force[1], 1e-16);
  EXPECT_NEAR(torque[1], lever[2] * force[0] - lever[0] * force[2], 1e-16);
  EXPECT_NEAR(torque[2], lever[0] * force[1] - lever[1] * force[0], 1e-16);
}

TEST(TorqueFromForce3D, ForceAndLeverNotParallel2)
{
  constexpr unsigned D = 3;
  const Vector<T, D> force(-4.2, -3.1, -7.5);
  const Vector<T, D> lever(-9.6, -2.3, -6.4);

  Vector<T, utilities::dimensions::convert<D>::rotation> torque =
      particles::dynamics::torque_from_force<D, T>::calculate(force, lever);

  EXPECT_NEAR(torque[0], lever[1] * force[2] - lever[2] * force[1], 1e-16);
  EXPECT_NEAR(torque[1], lever[2] * force[0] - lever[0] * force[2], 1e-16);
  EXPECT_NEAR(torque[2], lever[0] * force[1] - lever[1] * force[0], 1e-16);
}
