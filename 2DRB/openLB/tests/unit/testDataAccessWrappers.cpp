/*  This file is part of the OpenLB library
*
*  Copyright (C) 2022 Nicolas Hafen, Martin Sadric, Jan E. Marquardt
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

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::particles;
using namespace olb::particles::dynamics;
using namespace olb::particles::access;

typedef double T;

//#define VERBOSE_TESTS

// Field groups for testing
template<unsigned D>
struct FIELD_GROUP_WITH_DENSITY : public PARTICLE_DESCRIPTOR<D,DENSITY>, public PHYSPROPERTIES {};
template<unsigned D>
struct FIELD_GROUP_WITH_ELONGATION : public PARTICLE_DESCRIPTOR<D,ELONGATION>, public SURFACE {};
template<unsigned D>
struct FIELD_GROUP_WITH_RADIUS : public PARTICLE_DESCRIPTOR<D,RADIUS>, public PHYSPROPERTIES {};
template<unsigned D>
struct FIELD_GROUP_WITH_ADHESION : public PARTICLE_DESCRIPTOR<D,ADHESION>, public FORCING {};
template<unsigned D>
struct FIELD_GROUP_WITH_MATERIAL : public PARTICLE_DESCRIPTOR<D,MATERIAL>, public MECHPROPERTIES {};
template<unsigned D>
struct FIELD_GROUP_WITH_COR_OFFSET : public PARTICLE_DESCRIPTOR<D,COR_OFFSET>, public SURFACE {};
template<unsigned D>
struct FIELD_GROUP_WITH_ENLARGEMENT_FOR_CONTACT : public PARTICLE_DESCRIPTOR<D,ENLARGEMENT_FOR_CONTACT>, public NUMERICPROPERTIES {};

// TEST PROVIDES METHODS

TEST(TestDataAccessWrappers, providesID)
{
  typedef  PARTICLE_DESCRIPTOR<3, PARALLELIZATION_SUBGRID> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesID<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesID(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesID<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesID(particleB) );
}

TEST(TestDataAccessWrappers, providesInvalid)
{
  typedef  PARTICLE_DESCRIPTOR<3, GENERAL_EXTENDABLE<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesInvalid<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesInvalid(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesInvalid<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesInvalid(particleB) );
}

TEST(TestDataAccessWrappers, providesPosition)
{
  typedef  PARTICLE_DESCRIPTOR<3, GENERAL_TMP<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesPosition<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesPosition(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesPosition<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesPosition(particleB) );
}

TEST(TestDataAccessWrappers, providesRadius)
{
  typedef  PARTICLE_DESCRIPTOR<3, PHYSPROPERTIES_SUBGRID<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesRadius<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesRadius(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesRadius<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesRadius(particleB) );
}

TEST(TestDataAccessWrappers, providesMass)
{
  typedef  PARTICLE_DESCRIPTOR<3, PHYSPROPERTIES_RESOLVED<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesMass<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesMass(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesMass<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesMass(particleB) );
}

TEST(TestDataAccessWrappers, providesDensity)
{
  typedef  PARTICLE_DESCRIPTOR<3, FIELD_GROUP_WITH_DENSITY<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesDensity<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesDensity(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesDensity<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesDensity(particleB) );
}

TEST(TestDataAccessWrappers, providesMassOrDensity)
{
  typedef  PARTICLE_DESCRIPTOR<3, FIELD_GROUP_WITH_DENSITY<3>> PARTICLETYPE_WITH_DENSITY;
  ParticleSystem<T,PARTICLETYPE_WITH_DENSITY > particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesMassOrDensity<PARTICLETYPE_WITH_DENSITY>() );
  EXPECT_TRUE( providesMassOrDensity(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesMassOrDensity<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesMassOrDensity(particleB) );

  typedef  PARTICLE_DESCRIPTOR<3, PHYSPROPERTIES_RESOLVED<3>> PARTICLETYPE_WITH_MASS;
  ParticleSystem<T,PARTICLETYPE_WITH_MASS> particleSystemC(1);
  auto particleC = particleSystemC.get(0);
  EXPECT_TRUE( providesMass<PARTICLETYPE_WITH_MASS>() );
  EXPECT_TRUE( providesMass(particleC) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemD(1);
  auto particleD = particleSystemD.get(0);
  EXPECT_FALSE( providesMass<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesMass(particleD) );
}

TEST(TestDataAccessWrappers, providesElongation)
{
  typedef  PARTICLE_DESCRIPTOR<3, FIELD_GROUP_WITH_ELONGATION<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesElongation<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesElongation(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesElongation<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesElongation(particleB) );
}

TEST(TestDataAccessWrappers, providesAccelerationStrd)
{
  typedef  PARTICLE_DESCRIPTOR<3, MOBILITY_VERLET<3>> PARTICLETYPE_WITH;
  ParticleSystem<T, PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE(providesAccelerationStrd<PARTICLETYPE_WITH>());
  EXPECT_TRUE(providesAccelerationStrd(particleA));

  typedef PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T, PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE(providesAccelerationStrd<PARTICLETYPE_WITHOUT>());
  EXPECT_FALSE(providesAccelerationStrd(particleB));
}

TEST(TestDataAccessWrappers, providesAngle)
{
  typedef  PARTICLE_DESCRIPTOR<3, SURFACE_RESOLVED<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesAngle<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesAngle(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesAngle<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesAngle(particleB) );
}

TEST(TestDataAccessWrappers, providesVelocity)
{
  typedef  PARTICLE_DESCRIPTOR<3, MOBILITY_VERLET<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesVelocity<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesVelocity(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesVelocity<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesVelocity(particleB) );
}

TEST(TestDataAccessWrappers, providesAngVelocity)
{
  typedef  PARTICLE_DESCRIPTOR<3, MOBILITY_VERLET<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesVelocity<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesVelocity(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesVelocity<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesVelocity(particleB) );
}

TEST(TestDataAccessWrappers, providesAngAccelerationStrd)
{
  typedef  PARTICLE_DESCRIPTOR<3, MOBILITY_VERLET<3>> PARTICLETYPE_WITH;
  ParticleSystem<T, PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE(providesAngAccelerationStrd<PARTICLETYPE_WITH>());
  EXPECT_TRUE(providesAngAccelerationStrd(particleA));

  typedef PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T, PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE(providesAngAccelerationStrd<PARTICLETYPE_WITHOUT>());
  EXPECT_FALSE(providesAngAccelerationStrd(particleB));
}

TEST(TestDataAccessWrappers, providesForce)
{
  typedef  PARTICLE_DESCRIPTOR<3, FORCING_RESOLVED<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesForce<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesForce(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesForce<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesForce(particleB) );
}

TEST(TestDataAccessWrappers, providesTorque)
{
  typedef  PARTICLE_DESCRIPTOR<3, FORCING_RESOLVED<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesTorque<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesTorque(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesTorque<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesTorque(particleB) );
}

TEST(TestDataAccessWrappers, providesAdhesion)
{
  typedef  PARTICLE_DESCRIPTOR<3, FORCING_ADHESIVE<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesAdhesion<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesAdhesion(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesAdhesion<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesAdhesion(particleB) );
}

TEST(TestDataAccessWrappers, providesMomentOfInertia)
{
  typedef  PARTICLE_DESCRIPTOR<3, PHYSPROPERTIES_RESOLVED<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesMomentOfInertia<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesMomentOfInertia(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;


  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesMomentOfInertia<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesMomentOfInertia(particleB) );
}

TEST(TestDataAccessWrappers, providesValid)
{
  typedef  PARTICLE_DESCRIPTOR<3, GENERAL_EXTENDABLE <3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesValid<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesValid(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;


  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesValid<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesValid(particleB) );
}

TEST(TestDataAccessWrappers, providesDynamicsID)
{
  typedef  PARTICLE_DESCRIPTOR<3, DYNBEHAVIOUR_MULTI_DYN> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesDynamicsID<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesDynamicsID(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesDynamicsID<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesDynamicsID(particleB) );
}

TEST(TestDataAccessWrappers, providesComputeMotion)
{
  typedef  PARTICLE_DESCRIPTOR<3, DYNBEHAVIOUR_MULTI_DYN> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesComputeMotion<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesComputeMotion(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesComputeMotion<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesComputeMotion(particleB) );
}

TEST(TestDataAccessWrappers, providesComputeContact)
{
  typedef  PARTICLE_DESCRIPTOR<3, DYNBEHAVIOUR_DETACHABLE> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesComputeContact<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesComputeContact(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesComputeContact<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesComputeContact(particleB) );
}

TEST(TestDataAccessWrappers, providesSurface)
{
  typedef  PARTICLE_DESCRIPTOR<3, SURFACE_RESOLVED<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesSurface<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesSurface(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesSurface<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesSurface(particleB) );
}

TEST(TestDataAccessWrappers, providesParallelization)
{
  typedef  PARTICLE_DESCRIPTOR<3, PARALLELIZATION_SUBGRID> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesParallelization<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesParallelization(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesParallelization<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesParallelization(particleB) );
}

TEST(TestDataAccessWrappers, providesSpecies)
{
  typedef  PARTICLE_DESCRIPTOR<3, PHYSPROPERTIES_SUBGRID_REACTIVE<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesSpecies<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesSpecies(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesSpecies<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesSpecies(particleB) );
}

TEST(TestDataAccessWrappers, providesActive)
{
  typedef  PARTICLE_DESCRIPTOR<3, DYNBEHAVIOUR_BASIC> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesActive<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesActive(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesActive<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesActive(particleB) );
}

TEST(TestDataAccessWrappers, providesRotationMatrix)
{
  typedef  PARTICLE_DESCRIPTOR<3, SURFACE_RESOLVED<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesRotationMatrix<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesRotationMatrix(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesRotationMatrix<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesRotationMatrix(particleB) );
}

TEST(TestDataAccessWrappers, providesSmoothIndicator)
{
  typedef  PARTICLE_DESCRIPTOR<3, SURFACE_RESOLVED<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesSmoothIndicator<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesSmoothIndicator(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesSmoothIndicator<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesSmoothIndicator(particleB) );
}

TEST(TestDataAccessWrappers, providesDetaching)
{
  typedef  PARTICLE_DESCRIPTOR<3, DYNBEHAVIOUR_DETACHABLE> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesDetaching<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesDetaching(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesDetaching<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesDetaching(particleB) );
}

TEST(TestDataAccessWrappers, providesCORoffset)
{
  typedef  PARTICLE_DESCRIPTOR<3, SURFACE_RESOLVED_COR<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesCORoffset<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesCORoffset(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesCORoffset<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesCORoffset(particleB) );
}

TEST(TestDataAccessWrappers, providesEnlargementForContactTreatment)
{
  typedef  PARTICLE_DESCRIPTOR<3, NUMERICPROPERTIES_RESOLVED_CONTACT<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesEnlargementForContactTreatment<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesEnlargementForContactTreatment(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesEnlargementForContactTreatment<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesEnlargementForContactTreatment(particleB) );
}

TEST(TestDataAccessWrappers, providesContactMaterial)
{
  typedef  PARTICLE_DESCRIPTOR<3, MECHPROPERTIES_COLLISION<3>> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particleA = particleSystemA.get(0);
  EXPECT_TRUE( providesContactMaterial<PARTICLETYPE_WITH>() );
  EXPECT_TRUE( providesContactMaterial(particleA) );

  typedef  PARTICLE_DESCRIPTOR<3> PARTICLETYPE_WITHOUT;
  ParticleSystem<T,PARTICLETYPE_WITHOUT> particleSystemB(1);
  auto particleB = particleSystemB.get(0);
  EXPECT_FALSE( providesContactMaterial<PARTICLETYPE_WITHOUT>() );
  EXPECT_FALSE( providesContactMaterial(particleB) );
}

TEST(TestDataAccessWrappers, is2D)
{
  typedef ResolvedParticle2D PARTICLETYPE_1;
  ParticleSystem<T,PARTICLETYPE_1> particleSystem(2);
  auto particle_1 = particleSystem.get(1);

  typedef ResolvedParticle3D PARTICLETYPE_2;
  ParticleSystem<T,PARTICLETYPE_2> particleSystem_2(2);
  auto particle_2 = particleSystem_2.get(1);

  EXPECT_TRUE(is2D(particle_1));
  EXPECT_FALSE(is2D(particle_2));

}

TEST(TestDataAccessWrappers, is3D)
{
  typedef ResolvedParticle3D PARTICLETYPE_1;
  ParticleSystem<T,PARTICLETYPE_1> particleSystem(2);
  auto particle_1 = particleSystem.get(1);

  typedef ResolvedParticle2D PARTICLETYPE_2;
  ParticleSystem<T,PARTICLETYPE_2> particleSystem_2(2);
  auto particle_2 = particleSystem_2.get(1);

  EXPECT_TRUE(is3D(particle_1));
  EXPECT_FALSE(is3D(particle_2));

}

TEST(TestDataAccessWrappers, getPosition2D)
{
  typedef ResolvedParticle2D PARTICLETYPE;
  Vector<T,PARTICLETYPE::d> valInit(.5, 6.);

  ParticleSystem<T,PARTICLETYPE> particleSystem(2);
  auto particle = particleSystem.get(1);
  particle.setField<GENERAL,POSITION>( valInit );
  auto position = getPosition( particle );
  EXPECT_EQ(position, valInit);
  EXPECT_EQ(position.getDim(), PARTICLETYPE::d);
}

TEST(TestDataAccessWrappers, getPosition3D)
{
  typedef ResolvedSphere3D PARTICLETYPE;
  Vector<T,PARTICLETYPE::d> valInit(.2, 4., .9);

  ParticleSystem<T,PARTICLETYPE> particleSystem(3);
  auto particle = particleSystem.get(2);
  particle.setField<GENERAL,POSITION>( valInit );
  auto position = getPosition( particle );
  EXPECT_EQ(position, valInit);
  EXPECT_EQ(position.getDim(), PARTICLETYPE::d);
}

TEST(TestDataAccessWrappers, getVolume2D)
{
  typedef  PARTICLE_DESCRIPTOR<2, FIELD_GROUP_WITH_RADIUS<2>> PARTICLETYPE;
  ParticleSystem<T, PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  // Assuming the particle has a radius of 2.0
  setRadius(particle, 2.0);

  // Assuming the shape factor is 1.5
  T shapeFactor = 1.5;

  T expectedVolume = M_PI * 2.0 * 2.0 * shapeFactor;
  T volume = getVolume(particle, shapeFactor);

  EXPECT_EQ(volume, expectedVolume);
}

TEST(TestDataAccessWrappers, getVolume3D)
{
  typedef  PARTICLE_DESCRIPTOR<3, FIELD_GROUP_WITH_RADIUS<3>> PARTICLETYPE;
  ParticleSystem<T, PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  // Assuming the particle has a radius of 3.0
  setRadius(particle, 3.0);

  // Assuming the shape factor is 1.8
  T shapeFactor = 1.8;

  T expectedVolume = (4.0 / 3.0) * M_PI * 3.0 * 3.0 * 3.0 * shapeFactor;
  T volume = getVolume(particle, shapeFactor);

  EXPECT_NEAR(volume, expectedVolume,1e-9);
}

TEST(TestDataAccessWrappers, getAndSetAngle2D)
{
  typedef ResolvedParticle2D PARTICLETYPE;

  Vector<T, utilities::dimensions::convert<PARTICLETYPE::d>::rotation> valInit(0.5);
  ParticleSystem<T, PARTICLETYPE> particleSystem(2);
  auto particle = particleSystem.get(1);
  setAngle(particle, valInit);
  auto angle = getAngle(particle);

  EXPECT_EQ(angle, valInit);
  EXPECT_EQ(angle.getDim(), utilities::dimensions::convert<PARTICLETYPE::d>::rotation);
}

TEST(TestDataAccessWrappers, getAndSetAngle3D)
{
  typedef ResolvedParticle3D PARTICLETYPE;

  Vector<T, utilities::dimensions::convert<PARTICLETYPE::d>::rotation> valInit(0.5, 6.0, 0.9);
  ParticleSystem<T, PARTICLETYPE> particleSystem(4);
  auto particle = particleSystem.get(3);
  setAngle(particle, valInit);
  auto angle = getAngle(particle);

  EXPECT_EQ(angle, valInit);
  EXPECT_EQ(angle.getDim(), utilities::dimensions::convert<PARTICLETYPE::d>::rotation);
}

TEST(TestDataAccessWrappers, getAndSetAngleCORoffset2D)
{
  typedef  PARTICLE_DESCRIPTOR<2, FIELD_GROUP_WITH_COR_OFFSET<2>> PARTICLETYPE;

  Vector<T,PARTICLETYPE::d> valInit(0.5, 0.6);
  ParticleSystem<T, PARTICLETYPE> particleSystem(2);
  auto particle = particleSystem.get(1);
  setCORoffset(particle, valInit);

  EXPECT_EQ(getCORoffset(particle), valInit);
}

TEST(TestDataAccessWrappers, getAndSetAngleCORoffset3D)
{
  typedef  PARTICLE_DESCRIPTOR<3, FIELD_GROUP_WITH_COR_OFFSET<3>> PARTICLETYPE;

  Vector<T,PARTICLETYPE::d> valInit(0.5, 0.6, 0.7);
  ParticleSystem<T, PARTICLETYPE> particleSystem(2);
  auto particle = particleSystem.get(1);
  setCORoffset(particle, valInit);

  EXPECT_EQ(getCORoffset(particle), valInit);
}

TEST(TestDataAccessWrappers, getAndSetRotationMatrix2D)
{
  typedef ResolvedParticle2D PARTICLETYPE;
  Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::matrix> valInit(.5, 6., .9,.5);
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);
  auto particle = particleSystem.get(1);
  setRotationMatrix(particle, valInit);
  auto position = getRotationMatrix( particle );
  EXPECT_EQ(position, valInit);
  EXPECT_EQ(position.getDim(), utilities::dimensions::convert<PARTICLETYPE::d>::matrix);
}

TEST(TestDataAccessWrappers, getAndSetRotationMatrix3D)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::matrix> valInit(.5, 6.,.9,.5,0.,0.,0.,0.,8.);

  ParticleSystem<T,PARTICLETYPE> particleSystem(4);
  auto particle = particleSystem.get(3);
  setRotationMatrix(particle, valInit);
  auto position = getRotationMatrix( particle );
  EXPECT_EQ(position, valInit);
  EXPECT_EQ(position.getDim(), utilities::dimensions::convert<PARTICLETYPE::d>::matrix);
}

TEST(TestDataAccessWrappers, getAndSetContactMaterial2D)
{
  typedef  PARTICLE_DESCRIPTOR<2, FIELD_GROUP_WITH_MATERIAL<2>> PARTICLETYPE;
  int valInit(3);

  ParticleSystem<T,PARTICLETYPE> particleSystem(4);
  auto particle = particleSystem.get(3);
  setContactMaterial(particle, valInit);
  auto material = particle.template getField<MECHPROPERTIES, MATERIAL>();

  EXPECT_EQ(material, valInit);
}

TEST(TestDataAccessWrappers, getAndSetContactMaterial3D)
{
  typedef  PARTICLE_DESCRIPTOR<3, FIELD_GROUP_WITH_MATERIAL<3>> PARTICLETYPE;
  int valInit(3);

  ParticleSystem<T,PARTICLETYPE> particleSystem(4);
  auto particle = particleSystem.get(3);
  setContactMaterial(particle, valInit);
  auto material = particle.template getField<MECHPROPERTIES, MATERIAL>();

  EXPECT_EQ(material, valInit);
}

TEST(TestDataAccessWrappers, getSurfaceNormal3D)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::matrix> valInit(.5, 6.,.9,.5,0.,0.,0.,0.,8.);

  ParticleSystem<T,PARTICLETYPE> particleSystem(4);
  auto particle = particleSystem.get(3);
  particle.setField<SURFACE, ROT_MATRIX>( valInit );

  Vector<T,3> normal(.9,.0,8.);
  EXPECT_EQ(normal, getSurfaceNormal(particle));

}

TEST(TestDataAccessWrappers, getAndSetVelocity2D)
{
  typedef ResolvedParticle2D PARTICLETYPE;
  Vector<T,PARTICLETYPE::d> valInit(.5);
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);

  auto particle = particleSystem.get(1);
  setVelocity(particle, valInit );
  EXPECT_EQ(valInit, getVelocity(particle));

}

TEST(TestDataAccessWrappers, getAndSetVelocity3D)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  Vector<T,PARTICLETYPE::d> valInit(.5);
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);

  auto particle = particleSystem.get(1);
  particle.setField<MOBILITY,VELOCITY>( valInit );
  EXPECT_EQ(valInit, getVelocity(particle));

}

TEST(TestDataAccessWrappers, getAndSetAngularVelocity3D)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  Vector<T,PARTICLETYPE::d> valInit(.5);
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);

  auto particle = particleSystem.get(1);
  setAngularVelocity(particle,valInit);
  EXPECT_EQ(valInit, getAngularVelocity(particle));

}

TEST(TestDataAccessWrappers, getAndSetAngularVelocity2D)
{
  typedef ResolvedParticle2D PARTICLETYPE;
  Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> valInit(.5);
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);

  auto particle = particleSystem.get(1);
  setAngularVelocity(particle,valInit);
  EXPECT_EQ(valInit, getAngularVelocity(particle));

}

TEST(TestDataAccessWrappers, getForce2D)
{
  typedef ResolvedParticle2D PARTICLETYPE;
  Vector<T,PARTICLETYPE::d> valInit(.5);
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);

  auto particle = particleSystem.get(1);
  setForce(particle,valInit);
  EXPECT_EQ(valInit, getForce(particle));

}

TEST(TestDataAccessWrappers, getForce3D)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  Vector<T,PARTICLETYPE::d> valInit(.5);
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);

  auto particle = particleSystem.get(1);
  setForce(particle,valInit);
  EXPECT_EQ(valInit, getForce(particle));

}

TEST(TestDataAccessWrappers, getTorque2D)
{
  typedef ResolvedParticle2D PARTICLETYPE;
  Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> valInit(.2);
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);

  auto particle = particleSystem.get(1);
  setTorque(particle,valInit);
  EXPECT_EQ(valInit, getTorque(particle));

}

TEST(TestDataAccessWrappers, getTorque3D)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> valInit(.2e-14, 1.7e-13, 1.7e-13);
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);

  auto particle = particleSystem.get(1);
  setTorque(particle,valInit);
  EXPECT_EQ(valInit, getTorque(particle));

}

TEST(TestDataAccessWrappers, getAndSetEnlargementForContact2D)
{
  typedef  PARTICLE_DESCRIPTOR<2, FIELD_GROUP_WITH_ENLARGEMENT_FOR_CONTACT<2>> PARTICLETYPE;
  T valInit(.2);
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);

  auto particle = particleSystem.get(1);
  setEnlargementForContact(particle,valInit);
  EXPECT_EQ(valInit, getEnlargementForContact(particle));

}

TEST(TestDataAccessWrappers, getAndSetEnlargementForContact3D)
{
  typedef  PARTICLE_DESCRIPTOR<3, FIELD_GROUP_WITH_ENLARGEMENT_FOR_CONTACT<3>> PARTICLETYPE;
  T valInit(.2);
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);

  auto particle = particleSystem.get(1);
  setEnlargementForContact(particle,valInit);
  EXPECT_EQ(valInit, getEnlargementForContact(particle));

}

TEST(TestDataAccessWrappers, getAndSetMofi2D)
{
  typedef ResolvedParticle2D PARTICLETYPE;
  Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> valInit( 0.05 );
  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);
  setMomentOfInertia(particle, valInit );
  auto mofi = getMomentOfInertia( particle );
  EXPECT_EQ(mofi, valInit);
}

TEST(TestDataAccessWrappers, getAndSetMofi3D)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> valInit(.1, 0.15, .9);

  ParticleSystem<T,PARTICLETYPE> particleSystem(4);
  auto particle = particleSystem.get(2);
  setMomentOfInertia(particle,valInit);
  auto mofi = getMomentOfInertia( particle );
  EXPECT_EQ(mofi, valInit);
  EXPECT_EQ(mofi.getDim(), PARTICLETYPE::d);

}

TEST(TestDataAccessWrappers, getAndSetAdhesion2D)
{
  typedef  PARTICLE_DESCRIPTOR<2, FIELD_GROUP_WITH_ADHESION<2>> PARTICLETYPE;
  Vector<T,2> valInit(.5,.6);
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);

  auto particle = particleSystem.get(1);
  setAdhesion(particle,valInit);
  Vector<T,PARTICLETYPE::d> adhesion = getAdhesion(particle);
  EXPECT_EQ(valInit, adhesion);

}

TEST(TestDataAccessWrappers, getAndSetAccelerationStrd2D)
{
  typedef ResolvedParticle2D PARTICLETYPE;
  Vector<T,PARTICLETYPE::d> valInit(.5, 6.);
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);

  auto particle = particleSystem.get(1);
  setAccelerationStrd(particle,valInit);
  EXPECT_EQ(valInit, getAccelerationStrd(particle));
}

TEST(TestDataAccessWrappers, getAndSetAccelerationStrd3D)

{
  typedef ResolvedParticle3D PARTICLETYPE;
  Vector<T,PARTICLETYPE::d> valInit(.5, 6., .9);
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);

  auto particle = particleSystem.get(1);
  setAccelerationStrd(particle,valInit);
  EXPECT_EQ(valInit, getAccelerationStrd(particle));
}

TEST(TestDataAccessWrappers, getAcceleration2D)
{
  typedef ResolvedParticle2D PARTICLETYPE;
  Vector<T,PARTICLETYPE::d> forceInit(.2, 3.);
  T mass = 0.5;

  ParticleSystem<T,PARTICLETYPE> particleSystem(9);
  auto particle = particleSystem.get(5);
  setForce(particle,forceInit);
  particle.setField<PHYSPROPERTIES,MASS>( mass );
  auto acceleration = getAcceleration( particle );

  Vector<T,PARTICLETYPE::d> accelerationRef(
    forceInit[0] / mass,
    forceInit[1] / mass );
  EXPECT_EQ(acceleration, accelerationRef);
  EXPECT_EQ(acceleration.getDim(), PARTICLETYPE::d);
}

TEST(TestDataAccessWrappers, getAcceleration3D)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  Vector<T,PARTICLETYPE::d> forceInit(.5e-19, 2.6e-4,7.0e-12);
  T mass = 1.6;

  ParticleSystem<T,PARTICLETYPE> particleSystem(19);
  auto particle = particleSystem.get(12);
  particle.setField<FORCING,FORCE>( forceInit );
  particle.setField<PHYSPROPERTIES,MASS>( mass );
  auto acceleration = getAcceleration( particle );

  Vector<T,PARTICLETYPE::d> accelerationRef(
    forceInit[0] / mass,
    forceInit[1] / mass,
    forceInit[2] / mass );
  EXPECT_EQ(acceleration, accelerationRef);
  EXPECT_EQ(acceleration.getDim(), PARTICLETYPE::d);
}

TEST(TestDataAccessWrappers, getAngAcceleration2D)
{
  typedef ResolvedParticle2D PARTICLETYPE;
  T torqueInit = 4.2e-13;
  T mofi = 1.4e-15;

  ParticleSystem<T,PARTICLETYPE> particleSystem(4);
  auto particle = particleSystem.get(2);
  particle.setField<FORCING,TORQUE>( torqueInit );
  particle.setField<PHYSPROPERTIES,MOFI>( mofi );
  auto angAcceleration = getAngAcceleration( particle );

  T angAccelerationRef = torqueInit / mofi;
  EXPECT_EQ(angAcceleration[0], angAccelerationRef);
}

TEST(TestDataAccessWrappers, getAngAcceleration3D)

{
  typedef ResolvedParticle3D PARTICLETYPE;
  Vector<T,PARTICLETYPE::d> torqueInit(.4e-14, 1.5e-9,4.0e-18);
  Vector<T,PARTICLETYPE::d> mofi(1e-14,2e-15,3e-16);

  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);
  particle.setField<FORCING,TORQUE>( torqueInit );
  particle.setField<PHYSPROPERTIES,MOFI>( mofi );
  updateRotationMatrix( particle );
  auto angAcceleration = getAngAcceleration( particle );

  Vector<T,PARTICLETYPE::d> angAccelerationRef(
    torqueInit[0] / mofi[0],
    torqueInit[1] / mofi[1],
    torqueInit[2] / mofi[2] );

  T error = 1e-9;
  EXPECT_NEAR(angAcceleration[0], angAccelerationRef[0], error);
  EXPECT_NEAR(angAcceleration[1], angAccelerationRef[1], error);
  EXPECT_NEAR(angAcceleration[2], angAccelerationRef[2], error);
  EXPECT_EQ(angAcceleration.getDim(), PARTICLETYPE::d);

  T relErrorL2 = std::sqrt( //Represents change due to tensor adaption in #300
            ( std::pow(angAccelerationRef[0]-angAcceleration[0],2)
            + std::pow(angAccelerationRef[1]-angAcceleration[1],2)
            + std::pow(angAccelerationRef[2]-angAcceleration[2],2) )
          / ( std::pow(angAccelerationRef[0],2)
            + std::pow(angAccelerationRef[1],2)
            + std::pow(angAccelerationRef[2],2) ) );
  error = 1e-15;
  EXPECT_NEAR(relErrorL2, 0., error);
}

TEST(TestDataAccessWrappers, getAndSetAngAccelerationStrd2D)
{
  typedef ResolvedParticle2D PARTICLETYPE;
  Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> angAccelerationInit(.5);
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);

  auto particle = particleSystem.get(1);
  setAngAccelerationStrd(particle,angAccelerationInit);
  EXPECT_EQ(angAccelerationInit, getAngAccelerationStrd(particle));
}

TEST(TestDataAccessWrappers, getAndSetAngAccelerationStrd3D)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  Vector<T,utilities::dimensions::convert<PARTICLETYPE::d>::rotation> angAccelerationInit(.5);
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);

  auto particle = particleSystem.get(1);
  setAngAccelerationStrd(particle,angAccelerationInit);
  EXPECT_EQ(angAccelerationInit, getAngAccelerationStrd(particle));
}

TEST(TestDataAccessWrappers, isValid)
{
  typedef  PARTICLE_DESCRIPTOR<3, GENERAL_EXTENDABLE<3>> PARTICLETYPE;

  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setValid(particle, true);
  EXPECT_TRUE(isValid(particle));

  setValid(particle, false);
  EXPECT_FALSE(isValid(particle));
}

TEST(TestDataAccessWrappers, isActive)
{
  typedef  PARTICLE_DESCRIPTOR<3, DYNBEHAVIOUR_BASIC> PARTICLETYPE;

  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setActive(particle, true);
  EXPECT_TRUE(isActive(particle));

  setActive(particle, false);
  EXPECT_FALSE(isActive(particle));
}

TEST(TestDataAccessWrappers, getAndSetDensity2D)
{
  typedef  PARTICLE_DESCRIPTOR<2, FIELD_GROUP_WITH_DENSITY<2>> PARTICLETYPE;
  T valInit = 0.05;

  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setDensity(particle,valInit);
  EXPECT_EQ(getDensity(particle), valInit);
}

TEST(TestDataAccessWrappers, getAndSetDensity3D)
{
  typedef  PARTICLE_DESCRIPTOR<3, FIELD_GROUP_WITH_DENSITY<3>> PARTICLETYPE;
  T valInit = 0.05;

  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setDensity(particle,valInit);
  EXPECT_EQ(getDensity(particle), valInit);
}

TEST(TestDataAccessWrappers, getAndSetElongation2D)
{
  typedef  PARTICLE_DESCRIPTOR<2, FIELD_GROUP_WITH_ELONGATION<2>> PARTICLETYPE;
  Vector<T,PARTICLETYPE::d> valInit(0.5, 0.6);

  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setElongation(particle,valInit);
  EXPECT_EQ(getElongation(particle), valInit);
}

TEST(TestDataAccessWrappers, getAndSetElongation3D)
{
  typedef  PARTICLE_DESCRIPTOR<3, FIELD_GROUP_WITH_ELONGATION<3>> PARTICLETYPE;
  Vector<T,PARTICLETYPE::d> valInit(0.5, 0.6, 0.7);

  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setElongation(particle,valInit);
  EXPECT_EQ(getElongation(particle), valInit);
}

TEST(TestDataAccessWrappers, getAndSetGlobalID2D)
{
  std::size_t valInit = 42;
  typedef  PARTICLE_DESCRIPTOR<2, PARALLELIZATION_SUBGRID> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particle = particleSystemA.get(0);
  setGlobalID(particle, valInit);

  EXPECT_EQ(getGlobalID(particle),valInit);

}

TEST(TestDataAccessWrappers, getAndSetGlobalID3D)
{
  std::size_t valInit = 42;
  typedef  PARTICLE_DESCRIPTOR<3, PARALLELIZATION_SUBGRID> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particle = particleSystemA.get(0);
  setGlobalID(particle, valInit);

  EXPECT_EQ(getGlobalID(particle),valInit);

}

TEST(TestDataAccessWrappers, getAndSetGlobalIC2D)
{
  int valInit = 42;
  typedef  PARTICLE_DESCRIPTOR<2, PARALLELIZATION_RESOLVED> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particle = particleSystemA.get(0);
  setGlobalIC(particle, valInit);

  EXPECT_EQ(getGlobalIC(particle),valInit);

}

TEST(TestDataAccessWrappers, getAndSetGlobalIC3D)
{
  int valInit = 42;
  typedef  PARTICLE_DESCRIPTOR<3, PARALLELIZATION_RESOLVED> PARTICLETYPE_WITH;
  ParticleSystem<T,PARTICLETYPE_WITH> particleSystemA(1);
  auto particle = particleSystemA.get(0);
  setGlobalIC(particle, valInit);

  EXPECT_EQ(getGlobalIC(particle),valInit);

}

TEST(TestDataAccessWrappers, setAndGetMass2D)
{
  typedef ResolvedParticle2D PARTICLETYPE;
  T valInit = 0.05;

  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setMass(particle,valInit);
  EXPECT_EQ(getMass(particle), valInit);
}

TEST(TestDataAccessWrappers, setAndGetMass3D)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  T valInit = 0.05;

  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setMass(particle,valInit);
  EXPECT_EQ(getMass(particle), valInit);
}

TEST(TestDataAccessWrappers, setAndGetRadius2D)
{
  typedef  PARTICLE_DESCRIPTOR<2, FIELD_GROUP_WITH_RADIUS<2>> PARTICLETYPE;
  T valInit = 0.05;

  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setRadius(particle,valInit);
  EXPECT_EQ(getRadius(particle), valInit);
}

TEST(TestDataAccessWrappers, setAndGetRadius3D)
{
  typedef  PARTICLE_DESCRIPTOR<3, FIELD_GROUP_WITH_RADIUS<3>> PARTICLETYPE;
  T valInit = 0.05;

  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setRadius(particle,valInit);
  EXPECT_EQ(getRadius(particle), valInit);
}

TEST(TestDataAccessWrappers, setAndGetDensity2D)
{
  typedef  PARTICLE_DESCRIPTOR<3, FIELD_GROUP_WITH_DENSITY<3>> PARTICLETYPE;
  T valInit = 0.05;

  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setDensity(particle,valInit);
  EXPECT_EQ(getDensity(particle), valInit);
}

TEST(TestDataAccessWrappers, setAndGetDensity3D)
{
  typedef  PARTICLE_DESCRIPTOR<3, FIELD_GROUP_WITH_DENSITY<3>> PARTICLETYPE;
  T valInit = 0.05;

  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setDensity(particle,valInit);
  EXPECT_EQ(getDensity(particle), valInit);
}

TEST(TestDataAccessWrappers, setInvalid)
{
  typedef  PARTICLE_DESCRIPTOR<3, GENERAL_EXTENDABLE<3>> PARTICLETYPE;

  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setInvalid(particle, false);
  EXPECT_TRUE(isValid(particle));

  setInvalid(particle, true);
  EXPECT_FALSE(isValid(particle));

  setInvalid(particle);
  EXPECT_FALSE(isValid(particle));
}

TEST(TestDataAccessWrappers, setValid)
{
  typedef  PARTICLE_DESCRIPTOR<3, GENERAL_EXTENDABLE<3>> PARTICLETYPE;

  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setValid(particle, true);
  EXPECT_TRUE(isValid(particle));

  setValid(particle);
  EXPECT_TRUE(isValid(particle));

  setValid(particle, false);
  EXPECT_FALSE(isValid(particle));
}

TEST(TestDataAccessWrappers, setActive)
{
  typedef  PARTICLE_DESCRIPTOR<3, DYNBEHAVIOUR_BASIC> PARTICLETYPE;

  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setActive(particle, true);
  EXPECT_TRUE(isActive(particle));

  setActive(particle);
  EXPECT_TRUE(isActive(particle));

  setActive(particle, false);
  EXPECT_FALSE(isActive(particle));
}

TEST(TestDataAccessWrappers, setInactive)
{
  typedef  PARTICLE_DESCRIPTOR<3, DYNBEHAVIOUR_BASIC> PARTICLETYPE;

  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setInactive(particle, true);
  EXPECT_FALSE(isActive(particle));

  setInactive(particle);
  EXPECT_FALSE(isActive(particle));

  setInactive(particle, false);
  EXPECT_TRUE(isActive(particle));
}

TEST(TestDataAccessWrappers, setAndGetSmoothIndicatorPtr2D)
{
  typedef ResolvedParticle2D PARTICLETYPE;

  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  Vector<T,PARTICLETYPE::d> center(1.5, -2.3);
  T radius(3.75);
  std::unique_ptr<SmoothIndicatorF2D<T, T, true>> sIndicatorPtr(
    new SmoothIndicatorCircle2D<T, T, true>(center, radius, 0.));

  setSmoothIndicatorPtr(particle, sIndicatorPtr.get());

  EXPECT_EQ(getSmoothIndicatorPtr(particle), sIndicatorPtr.get());
}

TEST(TestDataAccessWrappers, setAndGetSmoothIndicatorPtr3D)
{
  typedef ResolvedParticle3D PARTICLETYPE;

  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  Vector<T,PARTICLETYPE::d> center(1.5, -2.3, 7.7);
  T radius(3.75);
  std::unique_ptr<SmoothIndicatorF3D<T, T, true>> sIndicatorPtr(
    new SmoothIndicatorSphere3D<T, T, true>(center, radius, 0.));

  setSmoothIndicatorPtr(particle, sIndicatorPtr.get());

  EXPECT_EQ(getSmoothIndicatorPtr(particle), sIndicatorPtr.get());
}

TEST(TestDataAccessWrappers, setAndGetDynamicsID2D)
{
  typedef  PARTICLE_DESCRIPTOR<2, DYNBEHAVIOUR_MULTI_DYN> PARTICLETYPE;
  unsigned short valInit = 42;
  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setDynamicsID(particle, valInit);
  EXPECT_EQ(getDynamicsID(particle), valInit);

}

TEST(TestDataAccessWrappers, setAndGetDynamicsID3D)
{
  typedef  PARTICLE_DESCRIPTOR<3, DYNBEHAVIOUR_MULTI_DYN> PARTICLETYPE;
  unsigned short valInit = 42;
  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setDynamicsID(particle, valInit);
  EXPECT_EQ(getDynamicsID(particle), valInit);

}

TEST(TestDataAccessWrappers, setAndIsDetaching2D)
{
  typedef  PARTICLE_DESCRIPTOR<2, DYNBEHAVIOUR_DETACHABLE> PARTICLETYPE;
  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setDetaching(particle, true);
  EXPECT_TRUE(isDetaching(particle));

  setDetaching(particle, false);
  EXPECT_FALSE(isDetaching(particle));

  setDetaching(particle, true);
  EXPECT_TRUE(isDetaching(particle));

}

TEST(TestDataAccessWrappers, setAndIsDetaching3D)
{
  typedef  PARTICLE_DESCRIPTOR<3, DYNBEHAVIOUR_DETACHABLE> PARTICLETYPE;
  ParticleSystem<T,PARTICLETYPE> particleSystem(1);
  auto particle = particleSystem.get(0);

  setDetaching(particle, true);
  EXPECT_TRUE(isDetaching(particle));

  setDetaching(particle, false);
  EXPECT_FALSE(isDetaching(particle));

  setDetaching(particle, true);
  EXPECT_TRUE(isDetaching(particle));

}

TEST(TestDataAccessWrappers, setPosition2D)
{
  typedef ResolvedParticle2D PARTICLETYPE;
  Vector<T,PARTICLETYPE::d> valInit(.5, 6.);

  ParticleSystem<T,PARTICLETYPE> particleSystem(2);
  auto particle = particleSystem.get(1);
  setPosition(particle, valInit);
  auto position = getPosition( particle );
  EXPECT_EQ(position, valInit);
  EXPECT_EQ(position.getDim(), PARTICLETYPE::d);
}

TEST(TestDataAccessWrappers, setPosition3D)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  Vector<T,PARTICLETYPE::d> valInit(.5, 6.,7.4);

  ParticleSystem<T,PARTICLETYPE> particleSystem(2);
  auto particle = particleSystem.get(1);
  setPosition(particle, valInit);
  auto position = getPosition( particle );
  EXPECT_EQ(position, valInit);
  EXPECT_EQ(position.getDim(), PARTICLETYPE::d);
}

TEST(TestDataAccessWrappers, EnableMotionComputation2D)
{
  typedef PARTICLE_DESCRIPTOR<2, DYNBEHAVIOUR_MULTI_DYN> PARTICLETYPE_2D;
  ParticleSystem<T, PARTICLETYPE_2D> particleSystem(1);
  auto particle = particleSystem.get(0);

  //enableMotionComputation(particle, true);
  //EXPECT_TRUE(isMotionComputationEnabled(particle));

  enableMotionComputation(particle, false);
  EXPECT_FALSE(isMotionComputationEnabled(particle));
}

TEST(TestDataAccessWrappers, EnableMotionComputation3D)
{
  typedef PARTICLE_DESCRIPTOR<3, DYNBEHAVIOUR_MULTI_DYN> PARTICLETYPE_3D;
  ParticleSystem<T, PARTICLETYPE_3D> particleSystem(1);
  auto particle = particleSystem.get(0);

  enableMotionComputation(particle, true);
  EXPECT_TRUE(isMotionComputationEnabled(particle));

  enableMotionComputation(particle, false);
  EXPECT_FALSE(isMotionComputationEnabled(particle));
}

TEST(TestDataAccessWrappers, EnableContactComputation2D)
{
  typedef PARTICLE_DESCRIPTOR<2, DYNBEHAVIOUR_MULTI_DYN> PARTICLETYPE_2D;
  ParticleSystem<T, PARTICLETYPE_2D> particleSystem(1);
  auto particle = particleSystem.get(0);

  enableContactComputation(particle, true);
  EXPECT_TRUE(isContactComputationEnabled(particle));

  enableContactComputation(particle, false);
  EXPECT_FALSE(isContactComputationEnabled(particle));
}

TEST(TestDataAccessWrappers, EnableContactComputation3D)
{
  typedef PARTICLE_DESCRIPTOR<3, DYNBEHAVIOUR_MULTI_DYN> PARTICLETYPE_3D;
  ParticleSystem<T, PARTICLETYPE_3D> particleSystem(1);
  auto particle = particleSystem.get(0);

  enableContactComputation(particle, true);
  EXPECT_TRUE(isContactComputationEnabled(particle));

  enableContactComputation(particle, false);
  EXPECT_FALSE(isContactComputationEnabled(particle));
}

TEST(TestDataAccessWrappers, DisableMotionComputation2D)
{
  typedef PARTICLE_DESCRIPTOR<2, DYNBEHAVIOUR_MULTI_DYN> PARTICLETYPE_2D;
  ParticleSystem<T, PARTICLETYPE_2D> particleSystem(1);
  auto particle = particleSystem.get(0);

  enableMotionComputation(particle, true);
  EXPECT_TRUE(isMotionComputationEnabled(particle));

  disableMotionComputation(particle);
  EXPECT_FALSE(isMotionComputationEnabled(particle));
}

TEST(TestDataAccessWrappers, DisableMotionComputation3D)
{
  typedef PARTICLE_DESCRIPTOR<3, DYNBEHAVIOUR_MULTI_DYN> PARTICLETYPE_3D;
  ParticleSystem<T, PARTICLETYPE_3D> particleSystem(1);
  auto particle = particleSystem.get(0);

  enableMotionComputation(particle, true);
  EXPECT_TRUE(isMotionComputationEnabled(particle));

  disableMotionComputation(particle);
  EXPECT_FALSE(isMotionComputationEnabled(particle));
}

TEST(TestDataAccessWrappers, DisableContactComputation2D)
{
  typedef PARTICLE_DESCRIPTOR<2, DYNBEHAVIOUR_MULTI_DYN> PARTICLETYPE_2D;
  ParticleSystem<T, PARTICLETYPE_2D> particleSystem(1);
  auto particle = particleSystem.get(0);

  enableContactComputation(particle, true);
  EXPECT_TRUE(isContactComputationEnabled(particle));

  disableContactComputation(particle);
  EXPECT_FALSE(isContactComputationEnabled(particle));
}

TEST(TestDataAccessWrappers, DisableContactComputation3D)
{
  typedef PARTICLE_DESCRIPTOR<3, DYNBEHAVIOUR_MULTI_DYN> PARTICLETYPE_3D;
  ParticleSystem<T, PARTICLETYPE_3D> particleSystem(1);
  auto particle = particleSystem.get(0);

  enableContactComputation(particle, true);
  EXPECT_TRUE(isContactComputationEnabled(particle));

  disableContactComputation(particle);
  EXPECT_FALSE(isContactComputationEnabled(particle));
}

//ADVANCED UNIT TESTS

TEST(TestDataAccessWrappers, getAngAcceleration3DofRotatedAroundX)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  Vector<T,PARTICLETYPE::d> torqueInit(.2e-14, 1.7e-13, 1.7e-13);
  Vector<T,PARTICLETYPE::d> mofi(1.3e-19,5e-17,4e-18);
  Vector<T,PARTICLETYPE::d> angleRad( .5*M_PI,0.,0.); //Rotate by 90 degree

  ParticleSystem<T,PARTICLETYPE> particleSystem(2);
  auto particle = particleSystem.get(0);
  particle.setField<FORCING,TORQUE>( torqueInit );
  particle.setField<PHYSPROPERTIES,MOFI>( mofi );
  updateRotationMatrix( particle );
  auto angAcceleration = getAngAcceleration( particle );

  //Now, rotate particle and again retrieve angular acceleration
  particle.setField<SURFACE,ANGLE>( angleRad );
  updateRotationMatrix( particle );
  auto angAcceleration90 = getAngAcceleration( particle );

  //Compare with Y,Z axis switched
  const T error = 1e-9;
  EXPECT_NEAR(angAcceleration90[0], angAcceleration[0], error);
  EXPECT_NEAR(angAcceleration90[1], angAcceleration[2], error);
  EXPECT_NEAR(angAcceleration90[2], angAcceleration[1], error);
  EXPECT_EQ(angAcceleration.getDim(), PARTICLETYPE::d);
}

TEST(TestDataAccessWrappers, setRestingParticle3DResolved)
{
  typedef ResolvedParticle3D PARTICLETYPE;
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);
  auto particle = particleSystem.get(0);

  setVelocity(particle, Vector<T,PARTICLETYPE::d>(1.));
  setAccelerationStrd(particle, Vector<T,PARTICLETYPE::d>(4.));
  setAngularVelocity(particle, Vector<T,PARTICLETYPE::d>(7.));
  setAngAccelerationStrd(particle, Vector<T,PARTICLETYPE::d>(10.));

  setRestingParticle(particle);

  Vector<T, PARTICLETYPE::d> zeroVec(0.);

  EXPECT_EQ(getVelocity(particle), zeroVec);
  EXPECT_EQ(getAccelerationStrd(particle), zeroVec);
  EXPECT_EQ(getAngularVelocity(particle), zeroVec);
  EXPECT_EQ(getAngAccelerationStrd(particle), zeroVec);

}

TEST(TestDataAccessWrappers, setRestingParticle3DSubgrid)
{
  typedef SubgridParticle3D PARTICLETYPE;
  ParticleSystem<T,PARTICLETYPE> particleSystem(2);
  auto particle = particleSystem.get(0);

  setVelocity(particle, Vector<T,PARTICLETYPE::d>(1.));
  setAccelerationStrd(particle, Vector<T,PARTICLETYPE::d>(4.));

  setRestingParticle(particle);

  Vector<T, PARTICLETYPE::d> zeroVec(0.);

  EXPECT_EQ(getVelocity(particle), zeroVec);
  EXPECT_EQ(getAccelerationStrd(particle), zeroVec);

}



