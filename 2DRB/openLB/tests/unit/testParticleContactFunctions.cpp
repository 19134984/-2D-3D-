/*  This file is part of the OpenLB library
*
*  Copyright (C) 2020 Jan E. Marquardt, Luca Heim
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
using namespace olb::particles::contact;

typedef double T;

TEST(effectiveYoungsModulus, testEffectiveYoungsModulus1)
{
  EXPECT_EQ(2. / 3.,
            particles::contact::evalEffectiveYoungModulus(1., 1., 0.5, 0.5));
}

TEST(effectiveYoungsModulus, testEffectiveYoungsModulus2)
{
  EXPECT_NEAR(4.129032258,
              particles::contact::evalEffectiveYoungModulus(8., 4., 0.75, 0.5),
              1e-10);
}

TEST(dampingFactor, testDampingFactor1)
{
  EXPECT_NEAR(1827. / 692., particles::contact::evalDampingFactor(0.85, 0.1),
              1e-15);
}

TEST(dampingFactor, testDampingFactor2)
{
  EXPECT_NEAR(1827. / 50516., particles::contact::evalDampingFactor(0.85, 7.3),
              1e-15);
}

TEST(dampingFactor, testDampingFactor3)
{
  EXPECT_NEAR(29.08434296, particles::contact::evalDampingFactor(0.05, 0.37),
              1e-8);
}

TEST(forceViaVolumeCoefficient, testCylinderWall)
{
  const T      radius = 1;
  const T      height = 1;
  Vector<T, 3> min(-radius, -radius, -0.5 * height);
  Vector<T, 3> max(radius, radius, 0.5 * height);
  const T      volume = M_PI * radius * radius * height;
  EXPECT_NEAR(
      2. / util::sqrt(M_PI),
      particles::contact::evalForceViaVolumeCoefficient(min, max, volume),
      1e-15);
}

TEST(forceViaVolumeCoefficient, testSphereWall)
{
  const T      radius      = 1.;
  const T      penetration = radius / 1000.;
  const T      a           = util::sqrt(radius * radius -
                                        (radius - penetration) * (radius - penetration));
  Vector<T, 3> min(-a, -a, -0.5 * penetration);
  Vector<T, 3> max(a, a, 0.5 * penetration);
  // Volume of spherical cap
  const T volume =
      M_PI / 3. * penetration * penetration * (3 * radius - penetration);
  // the k by Nassauer and Kuna (2013) is only valid for a small overlap, therefore, a difference makes sense
  EXPECT_NEAR(
      4. / (3 * util::sqrt(M_PI)),
      particles::contact::evalForceViaVolumeCoefficient(min, max, volume),
      1e-4);
}

TEST(boundingBoxCorrection, testWallCylinderZeroOverlap)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  T const                         radius = 2, height = 5;
  particles::creators::addResolvedCylinder3D(
      particleSystem, Vector<T, 3>(25., 0.5 * height, 25.),
      Vector<T, 3>(0., 1., 0.), height, radius, 0.);

  Vector<T, PARTICLETYPE::d>                      origin(0.);
  Vector<T, PARTICLETYPE::d>                      extend(50.);
  IndicatorCuboid3D<T>                            cuboid(extend, origin);
  std::unique_ptr<IndicatorF<T, PARTICLETYPE::d>> wallInd =
      std::make_unique<IndicInverse<T, PARTICLETYPE::d>>(cuboid, origin,
                                                         extend);
  std::vector<SolidBoundary<T, 3>> solidBoundaries;
  solidBoundaries.push_back(SolidBoundary<T, PARTICLETYPE::d>(std::move(wallInd), 2, 1));

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.wallContacts.push_back(WALLCONTACTTYPE(0, 0));
  PhysR<T, 3> min = Vector<T, 3>(25 - radius, -1, 25 - radius);
  contactContainer.wallContacts[0].updateMinMax(min);
  PhysR<T, 3> max = Vector<T, 3>(25 + radius, +1, 25 + radius);
  contactContainer.wallContacts[0].updateMinMax(max);

  particles::contact::particle_wall<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                    WALLCONTACTTYPE>::
      correctBoundingBoxNewContact(particleSystem, solidBoundaries,
                                   contactContainer.wallContacts[0], physDeltaX,
                                   10);

  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[0], 25 - radius, 5e-4);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[1], 0., 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[2], 25 - radius, 5e-4);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[0], 25 + radius, 5e-4);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[1], 0., 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[2], 25 + radius, 5e-4);
}

TEST(boundingBoxCorrection, testWallCylinderWithSmallOverlap)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  T const                         radius = 2, height = 5;
  T const                         depth = 1e-8 * height;
  particles::creators::addResolvedCylinder3D(
      particleSystem, Vector<T, 3>(25., 0.5 * height - depth, 25.),
      Vector<T, 3>(0., 1., 0.), height, radius, 0.);

  Vector<T, PARTICLETYPE::d>                      origin(0.);
  Vector<T, PARTICLETYPE::d>                      extend(50.);
  IndicatorCuboid3D<T>                            cuboid(extend, origin);
  std::unique_ptr<IndicatorF<T, PARTICLETYPE::d>> wallInd =
      std::make_unique<IndicInverse<T, PARTICLETYPE::d>>(
          cuboid, origin - 5 * physDeltaX, extend + 5 * physDeltaX);
  std::vector<SolidBoundary<T, 3>> solidBoundaries;
  solidBoundaries.push_back(SolidBoundary<T, PARTICLETYPE::d>(std::move(wallInd), 2, 1));

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.wallContacts.push_back(WALLCONTACTTYPE(0, 0));
  PhysR<T, 3> min = Vector<T, 3>(25 - radius, -1, 25 - radius);
  contactContainer.wallContacts[0].updateMinMax(min);
  PhysR<T, 3> max = Vector<T, 3>(25 + radius, +1, 25 + radius);
  contactContainer.wallContacts[0].updateMinMax(max);

  particles::contact::particle_wall<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                    WALLCONTACTTYPE>::
      correctBoundingBoxNewContact(particleSystem, solidBoundaries,
                                   contactContainer.wallContacts[0], physDeltaX,
                                   10);

  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[0], 25 - radius, 5e-4);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[1], -depth, 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[2], 25 - radius, 5e-4);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[0], 25 + radius, 5e-4);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[1], 0., 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[2], 25 + radius, 5e-4);
}

TEST(boundingBoxCorrection, testWallCylinderWithOverlap)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  T const                         radius = 2, height = 5;
  T const                         depth = 0.01 * height;
  particles::creators::addResolvedCylinder3D(
      particleSystem, Vector<T, 3>(25., 0.5 * height - depth, 25.),
      Vector<T, 3>(0., 1., 0.), height, radius, 0.);

  Vector<T, PARTICLETYPE::d>                      origin(0.);
  Vector<T, PARTICLETYPE::d>                      extend(50.);
  IndicatorCuboid3D<T>                            cuboid(extend, origin);
  std::unique_ptr<IndicatorF<T, PARTICLETYPE::d>> wallInd =
      std::make_unique<IndicInverse<T, PARTICLETYPE::d>>(
          cuboid, origin - 5 * physDeltaX, extend + 5 * physDeltaX);
  std::vector<SolidBoundary<T, 3>> solidBoundaries;
  solidBoundaries.push_back(SolidBoundary<T, PARTICLETYPE::d>(std::move(wallInd), 2, 1));

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.wallContacts.push_back(WALLCONTACTTYPE(0, 0));
  contactContainer.wallContacts.push_back(WALLCONTACTTYPE(0, 0));
  PhysR<T, 3> min = Vector<T, 3>(25 - radius, -1, 25 - radius);
  contactContainer.wallContacts[0].updateMinMax(min);
  PhysR<T, 3> max = Vector<T, 3>(25 + radius, +1, 25 + radius);
  contactContainer.wallContacts[0].updateMinMax(max);

  particles::contact::particle_wall<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                    WALLCONTACTTYPE>::
      correctBoundingBoxNewContact(particleSystem, solidBoundaries,
                                   contactContainer.wallContacts[0], physDeltaX,
                                   10);

  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[0], 25 - radius,
              7.5e-4);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[1], -depth, 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[2], 25 - radius,
              7.5e-4);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[0], 25 + radius,
              7.5e-4);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[1], 0., 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[2], 25 + radius,
              7.5e-4);
}

TEST(boundingBoxCorrection, testWallCylinderWithBigOverlap)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  T const                         radius = 2, height = 5;
  T const                         depth = 0.25 * height;
  particles::creators::addResolvedCylinder3D(
      particleSystem, Vector<T, 3>(25., 0.5 * height - depth, 25.),
      Vector<T, 3>(0., 1., 0.), height, radius, 0.);

  Vector<T, PARTICLETYPE::d>                      origin(0.);
  Vector<T, PARTICLETYPE::d>                      extend(50.);
  IndicatorCuboid3D<T>                            cuboid(extend, origin);
  std::unique_ptr<IndicatorF<T, PARTICLETYPE::d>> wallInd =
      std::make_unique<IndicInverse<T, PARTICLETYPE::d>>(
          cuboid, origin - 5 * physDeltaX, extend + 5 * physDeltaX);
  std::vector<SolidBoundary<T, 3>> solidBoundaries;
  solidBoundaries.push_back(SolidBoundary<T, PARTICLETYPE::d>(std::move(wallInd), 2, 1));

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.wallContacts.push_back(WALLCONTACTTYPE(0, 0));
  // In order to achieve higher accuracy, the rough bounding box must be bigger than the expected resulting bounding box
  // because otherwise the check might leave the bounding circumradius
  contactContainer.wallContacts.push_back(WALLCONTACTTYPE(0, 0));
  PhysR<T, 3> min = Vector<T, 3>(25 - radius, -depth, 25 - radius);
  contactContainer.wallContacts[0].updateMinMax(min);
  PhysR<T, 3> max = Vector<T, 3>(25 + radius, physDeltaX, 25 + radius);
  contactContainer.wallContacts[0].updateMinMax(max);

  particles::contact::particle_wall<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                    WALLCONTACTTYPE>::
      correctBoundingBoxNewContact(particleSystem, solidBoundaries,
                                   contactContainer.wallContacts[0], physDeltaX,
                                   10);

  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[0], 25 - radius, 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[1], -depth, 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[2], 25 - radius, 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[0], 25 + radius, 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[1], 0., 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[2], 25 + radius, 1e-3);
}

TEST(boundingBoxCorrection, testSphereSphereZeroOverlapBigBox)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  T const                         radius       = 2;
  T const                         displacement = 0;
  PhysR<T, 3>                     pos1(0., -radius + displacement, 0.);
  PhysR<T, 3>                     pos2(0., radius, 0.);

  particles::creators::addResolvedSphere3D(particleSystem, pos1, radius, 0.,
                                           0.);
  particles::creators::addResolvedSphere3D(particleSystem, pos2, radius, 0.,
                                           0.);

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.particleContacts.push_back(
      PARTICLECONTACTTYPE(std::array<size_t, 2>({0, 1})));
  PhysR<T, 3> min = Vector<T, 3>(-radius);
  contactContainer.particleContacts[0].updateMinMax(min);
  PhysR<T, 3> max = Vector<T, 3>(radius);
  contactContainer.particleContacts[0].updateMinMax(max);

  particles::contact::particle_particle<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                        WALLCONTACTTYPE>::
      correctBoundingBoxNewContact(
          particleSystem, contactContainer.particleContacts[0], physDeltaX, 10);

  T const a =
      std::sqrt(radius * radius - std::pow(radius - 0.5 * displacement, 2));

  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[0], -a, 6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[1], 0., 6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[2], -a, 6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[0], a, 6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[1], 0., 6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[2], a, 6e-4);
}

TEST(boundingBoxCorrection, testSphereSphereSmallOverlapBigBox)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  T const                         radius       = 2;
  T const                         displacement = 1e-8;
  PhysR<T, 3>                     pos1(0., -radius + displacement, 0.);
  PhysR<T, 3>                     pos2(0., radius, 0.);

  particles::creators::addResolvedSphere3D(particleSystem, pos1, radius, 0.,
                                           0.);
  particles::creators::addResolvedSphere3D(particleSystem, pos2, radius, 0.,
                                           0.);

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.particleContacts.push_back(
      PARTICLECONTACTTYPE(std::array<size_t, 2>({0, 1})));
  PhysR<T, 3> min = Vector<T, 3>(-radius);
  contactContainer.particleContacts[0].updateMinMax(min);
  PhysR<T, 3> max = Vector<T, 3>(+radius);
  contactContainer.particleContacts[0].updateMinMax(max);

  particles::contact::particle_particle<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                        WALLCONTACTTYPE>::
      correctBoundingBoxNewContact(
          particleSystem, contactContainer.particleContacts[0], physDeltaX, 10);

  T const a =
      std::sqrt(radius * radius - std::pow(radius - 0.5 * displacement, 2));

  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[0], -a, 5e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[1], 0., 6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[2], -a, 5e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[0], a, 5e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[1], displacement,
              6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[2], a, 5e-4);
}

TEST(boundingBoxCorrection, testSphereSphereOverlapBigBox)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  T const                         radius       = 2;
  T const                         displacement = 0.1;
  PhysR<T, 3>                     pos1(0., -radius + displacement, 0.);
  PhysR<T, 3>                     pos2(0., radius, 0.);

  particles::creators::addResolvedSphere3D(particleSystem, pos1, radius, 0.,
                                           0.);
  particles::creators::addResolvedSphere3D(particleSystem, pos2, radius, 0.,
                                           0.);

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.particleContacts.push_back(
      PARTICLECONTACTTYPE(std::array<size_t, 2>({0, 1})));
  PhysR<T, 3> min = Vector<T, 3>(-radius);
  contactContainer.particleContacts[0].updateMinMax(min);
  PhysR<T, 3> max = Vector<T, 3>(+radius);
  contactContainer.particleContacts[0].updateMinMax(max);

  particles::contact::particle_particle<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                        WALLCONTACTTYPE>::
      correctBoundingBoxNewContact(
          particleSystem, contactContainer.particleContacts[0], physDeltaX, 10);

  T const a =
      std::sqrt(radius * radius - std::pow(radius - 0.5 * displacement, 2));

  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[0], -a, 5e-2);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[1], 0., 6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[2], -a, 5e-2);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[0], a, 5e-2);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[1], displacement,
              6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[2], a, 5e-2);
}

TEST(boundingBoxCorrection, testSphereSphereBigOverlapBigBox)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  T const                         radius       = 2;
  T const                         displacement = 0.5;
  PhysR<T, 3>                     pos1(0., -radius + displacement, 0.);
  PhysR<T, 3>                     pos2(0., radius, 0.);

  particles::creators::addResolvedSphere3D(particleSystem, pos1, radius, 0.,
                                           0.);
  particles::creators::addResolvedSphere3D(particleSystem, pos2, radius, 0.,
                                           0.);

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.particleContacts.push_back(
      PARTICLECONTACTTYPE(std::array<size_t, 2>({0, 1})));
  PhysR<T, 3> min = Vector<T, 3>(-radius);
  contactContainer.particleContacts[0].updateMinMax(min);
  PhysR<T, 3> max = Vector<T, 3>(+radius);
  contactContainer.particleContacts[0].updateMinMax(max);

  particles::contact::particle_particle<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                        WALLCONTACTTYPE>::
      correctBoundingBoxNewContact(
          particleSystem, contactContainer.particleContacts[0], physDeltaX, 10);

  T const a =
      std::sqrt(radius * radius - std::pow(radius - 0.5 * displacement, 2));

  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[0], -a, 2e-2);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[1], 0., 6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[2], -a, 2e-2);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[0], a, 2e-2);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[1], displacement,
              6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[2], a, 2e-2);
}

TEST(boundingBoxCorrection, testWallCylinderZeroOverlapEnlargedParticle)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX                    = 0.1;
  T const particleEnlargementForContact = 3 * physDeltaX;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  T const                         radius = 2, height = 5;
  T const enlargedRadius = radius + particleEnlargementForContact;
  particles::creators::addResolvedCylinder3D(
      particleSystem,
      Vector<T, 3>(25., 0.5 * height + particleEnlargementForContact, 25.),
      Vector<T, 3>(0., 1., 0.), height, radius, 0.);
  auto particleA = particleSystem.get(0);
  particleA.setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(
      particleEnlargementForContact);

  Vector<T, PARTICLETYPE::d> origin(-particleEnlargementForContact);
  Vector<T, PARTICLETYPE::d> extend(50. + particleEnlargementForContact);
  IndicatorCuboid3D<T>       cuboid(extend, origin);
  std::unique_ptr<IndicatorF<T, PARTICLETYPE::d>> wallInd =
      std::make_unique<IndicInverse<T, PARTICLETYPE::d>>(cuboid, origin,
                                                         extend);
  std::vector<SolidBoundary<T, 3>> solidBoundaries;
  solidBoundaries.push_back(SolidBoundary<T, PARTICLETYPE::d>(
      std::move(wallInd), 2, 1, particleEnlargementForContact));

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.wallContacts.push_back(WALLCONTACTTYPE(0, 0));
  PhysR<T, 3> min =
      Vector<T, 3>(25 - enlargedRadius, -physDeltaX, 25 - enlargedRadius);
  contactContainer.wallContacts[0].updateMinMax(min);
  PhysR<T, 3> max =
      Vector<T, 3>(25 + enlargedRadius, +physDeltaX, 25 + enlargedRadius);
  contactContainer.wallContacts[0].updateMinMax(max);

  particles::contact::particle_wall<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                    WALLCONTACTTYPE>::
      correctBoundingBoxNewContact(particleSystem, solidBoundaries,
                                   contactContainer.wallContacts[0], physDeltaX,
                                   10);

  // std::rounding rounds edges, therefore we expect a smaller overlap in x- and z- direction
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[0], 22.982, 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[1], 0., 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[2], 22.982, 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[0], 27.018, 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[1], 0., 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[2], 27.018, 1e-3);
}

TEST(boundingBoxCorrection, testWallCylinderWithSmallOverlapEnlargedParticle)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX                    = 0.1;
  T const particleEnlargementForContact = 3 * physDeltaX;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  T const                         radius = 2, height = 5;
  T const enlargedRadius = radius + particleEnlargementForContact;
  T const depth          = 1e-8 * height;
  particles::creators::addResolvedCylinder3D(
      particleSystem,
      Vector<T, 3>(25., 0.5 * height - depth + particleEnlargementForContact,
                   25.),
      Vector<T, 3>(0., 1., 0.), height, radius, 0.);
  auto particleA = particleSystem.get(0);
  particleA.setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(
      particleEnlargementForContact);

  Vector<T, PARTICLETYPE::d> origin(-particleEnlargementForContact);
  Vector<T, PARTICLETYPE::d> extend(50. + particleEnlargementForContact);
  IndicatorCuboid3D<T>       cuboid(extend, origin);
  std::unique_ptr<IndicatorF<T, PARTICLETYPE::d>> wallInd =
      std::make_unique<IndicInverse<T, PARTICLETYPE::d>>(
          cuboid, origin - 5 * physDeltaX, extend + 5 * physDeltaX);
  std::vector<SolidBoundary<T, 3>> solidBoundaries;
  solidBoundaries.push_back(SolidBoundary<T, PARTICLETYPE::d>(
      std::move(wallInd), 2, 1, particleEnlargementForContact));

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.wallContacts.push_back(WALLCONTACTTYPE(0, 0));
  PhysR<T, 3> min =
      Vector<T, 3>(25 - enlargedRadius, -depth, 25 - enlargedRadius);
  contactContainer.wallContacts[0].updateMinMax(min);
  PhysR<T, 3> max = Vector<T, 3>(25 + enlargedRadius, 0., 25 + enlargedRadius);
  contactContainer.wallContacts[0].updateMinMax(max);

  particles::contact::particle_wall<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                    WALLCONTACTTYPE>::
      correctBoundingBoxNewContact(particleSystem, solidBoundaries,
                                   contactContainer.wallContacts[0], physDeltaX,
                                   10);

  // std::rounding rounds edges, therefore we expect a smaller overlap in x- and z- direction
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[0], 22.982, 5e-4);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[1], -depth, 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[2], 22.982, 5e-4);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[0], 27.018, 5e-4);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[1], 0., 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[2], 27.018, 5e-4);
}

TEST(boundingBoxCorrection, testWallCylinderWithOverlapEnlargedParticle)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX                    = 0.1;
  T const particleEnlargementForContact = 3 * physDeltaX;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  T const                         radius = 2, height = 5;
  T const enlargedRadius = radius + particleEnlargementForContact;
  T const depth          = 0.01 * height;
  particles::creators::addResolvedCylinder3D(
      particleSystem,
      Vector<T, 3>(25., 0.5 * height - depth + particleEnlargementForContact,
                   25.),
      Vector<T, 3>(0., 1., 0.), height, radius, 0.);
  auto particleA = particleSystem.get(0);
  particleA.setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(
      particleEnlargementForContact);

  Vector<T, PARTICLETYPE::d> origin(-particleEnlargementForContact);
  Vector<T, PARTICLETYPE::d> extend(50. + particleEnlargementForContact);
  IndicatorCuboid3D<T>       cuboid(extend, origin);
  std::unique_ptr<IndicatorF<T, PARTICLETYPE::d>> wallInd =
      std::make_unique<IndicInverse<T, PARTICLETYPE::d>>(
          cuboid, origin - 5 * physDeltaX, extend + 5 * physDeltaX);
  std::vector<SolidBoundary<T, 3>> solidBoundaries;
  solidBoundaries.push_back(SolidBoundary<T, PARTICLETYPE::d>(
      std::move(wallInd), 2, 1, particleEnlargementForContact));

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.wallContacts.push_back(WALLCONTACTTYPE(0, 0));
  PhysR<T, 3> min = Vector<T, 3>(25 - enlargedRadius, -1, 25 - enlargedRadius);
  contactContainer.wallContacts[0].updateMinMax(min);
  PhysR<T, 3> max = Vector<T, 3>(25 + enlargedRadius, 1, 25 + enlargedRadius);
  contactContainer.wallContacts[0].updateMinMax(max);

  particles::contact::particle_wall<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                    WALLCONTACTTYPE>::
      correctBoundingBoxNewContact(particleSystem, solidBoundaries,
                                   contactContainer.wallContacts[0], physDeltaX,
                                   10);

  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[0], 22.8326, 7.5e-4);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[1], -depth, 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[2], 22.8326, 7.5e-4);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[0], 27.1610, 7.5e-4);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[1], 0., 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[2], 27.1610, 7.5e-4);
}

TEST(boundingBoxCorrection, testWallCylinderWithBigOverlapEnlargedParticle)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX                    = 0.1;
  T const particleEnlargementForContact = 3 * physDeltaX;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  T const                         radius = 2, height = 5;
  T const enlargedRadius = radius + particleEnlargementForContact;
  T const depth          = 0.25 * height;
  particles::creators::addResolvedCylinder3D(
      particleSystem,
      Vector<T, 3>(25., 0.5 * height - depth + particleEnlargementForContact,
                   25.),
      Vector<T, 3>(0., 1., 0.), height, radius, 0.);
  auto particleA = particleSystem.get(0);
  particleA.setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(
      particleEnlargementForContact);

  Vector<T, PARTICLETYPE::d> origin(-particleEnlargementForContact);
  Vector<T, PARTICLETYPE::d> extend(50. + particleEnlargementForContact);
  IndicatorCuboid3D<T>       cuboid(extend, origin);
  std::unique_ptr<IndicatorF<T, PARTICLETYPE::d>> wallInd =
      std::make_unique<IndicInverse<T, PARTICLETYPE::d>>(
          cuboid, origin - 5 * physDeltaX, extend + 5 * physDeltaX);
  std::vector<SolidBoundary<T, 3>> solidBoundaries;
  solidBoundaries.push_back(SolidBoundary<T, PARTICLETYPE::d>(
      std::move(wallInd), 2, 1, particleEnlargementForContact));

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.wallContacts.push_back(WALLCONTACTTYPE(0, 0));
  // In order to achieve higher accuracy, the rough bounding box must be bigger than the expected resulting bounding box
  // because otherwise the check might leave the bounding circumradius
  PhysR<T, 3> min = Vector<T, 3>(25 - enlargedRadius, -depth - physDeltaX,
                                 25 - enlargedRadius);
  contactContainer.wallContacts[0].updateMinMax(min);
  PhysR<T, 3> max =
      Vector<T, 3>(25 + enlargedRadius, physDeltaX, 25 + enlargedRadius);
  contactContainer.wallContacts[0].updateMinMax(max);

  particles::contact::particle_wall<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                    WALLCONTACTTYPE>::
      correctBoundingBoxNewContact(

          particleSystem, solidBoundaries, contactContainer.wallContacts[0],
          physDeltaX, 10);

  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[0], 25 - enlargedRadius,
              1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[1], -depth, 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[2], 25 - enlargedRadius,
              1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[0], 25 + enlargedRadius,
              1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[1], 0., 1e-3);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[2], 25 + enlargedRadius,
              1e-3);
}

TEST(boundingBoxCorrection, testSphereSphereZeroOverlapBigBoxEnlargedParticle)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX                    = 0.1;
  T const particleEnlargementForContact = 3 * physDeltaX;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  T const                         radius = 2;
  T const     enlargedRadius = radius + particleEnlargementForContact;
  T const     displacement   = 0;
  PhysR<T, 3> pos1(
      0., -radius + displacement - 2 * particleEnlargementForContact, 0.);
  PhysR<T, 3> pos2(0., radius, 0.);

  particles::creators::addResolvedSphere3D(particleSystem, pos1, radius, 0.,
                                           0.);
  auto particleA = particleSystem.get(0);
  particleA.setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(
      particleEnlargementForContact);
  particles::creators::addResolvedSphere3D(particleSystem, pos2, radius, 0.,
                                           0.);
  auto particleB = particleSystem.get(1);
  particleB.setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(
      particleEnlargementForContact);

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.particleContacts.push_back(
      PARTICLECONTACTTYPE(std::array<size_t, 2>({0, 1})));
  PhysR<T, 3> min = Vector<T, 3>(-enlargedRadius);
  contactContainer.particleContacts[0].updateMinMax(min);
  PhysR<T, 3> max = Vector<T, 3>(+enlargedRadius);
  contactContainer.particleContacts[0].updateMinMax(max);

  particles::contact::particle_particle<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                        WALLCONTACTTYPE>::
      correctBoundingBoxNewContact(
          particleSystem, contactContainer.particleContacts[0], physDeltaX, 10);

  T const a = std::sqrt(enlargedRadius * enlargedRadius -
                        std::pow(enlargedRadius - 0.5 * displacement, 2));

  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[0], -a, 6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[1],
              -particleEnlargementForContact, 6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[2], -a, 6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[0], a, 6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[1],
              -particleEnlargementForContact, 6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[2], a, 6e-4);
}

TEST(boundingBoxCorrection, testSphereSphereSmallOverlapBigBoxEnlargedParticle)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX                    = 0.1;
  T const particleEnlargementForContact = 3 * physDeltaX;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  T const                         radius = 2;
  T const     enlargedRadius = radius + particleEnlargementForContact;
  T const     displacement   = 1e-8;
  PhysR<T, 3> pos1(
      0., -radius + displacement - 2 * particleEnlargementForContact, 0.);
  PhysR<T, 3> pos2(0., radius, 0.);

  particles::creators::addResolvedSphere3D(particleSystem, pos1, radius, 0.,
                                           0.);
  auto particleA = particleSystem.get(0);
  particleA.setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(
      particleEnlargementForContact);
  particles::creators::addResolvedSphere3D(particleSystem, pos2, radius, 0.,
                                           0.);
  auto particleB = particleSystem.get(1);
  particleB.setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(
      particleEnlargementForContact);

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.particleContacts.push_back(
      PARTICLECONTACTTYPE(std::array<size_t, 2>({0, 1})));
  PhysR<T, 3> min = Vector<T, 3>(-enlargedRadius);
  contactContainer.particleContacts[0].updateMinMax(min);
  PhysR<T, 3> max = Vector<T, 3>(+enlargedRadius);
  contactContainer.particleContacts[0].updateMinMax(max);

  particles::contact::particle_particle<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                        WALLCONTACTTYPE>::
      correctBoundingBoxNewContact(
          particleSystem, contactContainer.particleContacts[0], physDeltaX, 10);

  T const a = std::sqrt(enlargedRadius * enlargedRadius -
                        std::pow(enlargedRadius - 0.5 * displacement, 2));

  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[0], -a, 5e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[1],
              -particleEnlargementForContact, 6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[2], -a, 5e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[0], a, 5e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[1],
              displacement - particleEnlargementForContact, 6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[2], a, 5e-4);
}

TEST(boundingBoxCorrection, testSphereSphereOverlapBigBoxEnlargedParticle)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX                    = 0.1;
  T const particleEnlargementForContact = 3 * physDeltaX;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  T const                         radius = 2;
  T const     enlargedRadius = radius + particleEnlargementForContact;
  T const     displacement   = 0.1;
  PhysR<T, 3> pos1(
      0., -radius + displacement - 2 * particleEnlargementForContact, 0.);
  PhysR<T, 3> pos2(0., radius, 0.);

  particles::creators::addResolvedSphere3D(particleSystem, pos1, radius, 0.,
                                           0.);
  auto particleA = particleSystem.get(0);
  particleA.setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(
      particleEnlargementForContact);
  particles::creators::addResolvedSphere3D(particleSystem, pos2, radius, 0.,
                                           0.);
  auto particleB = particleSystem.get(1);
  particleB.setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(
      particleEnlargementForContact);

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.particleContacts.push_back(
      PARTICLECONTACTTYPE(std::array<size_t, 2>({0, 1})));
  PhysR<T, 3> min = Vector<T, 3>(-enlargedRadius);
  contactContainer.particleContacts[0].updateMinMax(min);
  PhysR<T, 3> max = Vector<T, 3>(+enlargedRadius);
  contactContainer.particleContacts[0].updateMinMax(max);

  particles::contact::particle_particle<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                        WALLCONTACTTYPE>::
      correctBoundingBoxNewContact(
          particleSystem, contactContainer.particleContacts[0], physDeltaX, 10);

  T const a = std::sqrt(enlargedRadius * enlargedRadius -
                        std::pow(enlargedRadius - 0.5 * displacement, 2));

  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[0], -a, 5e-2);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[1],
              -particleEnlargementForContact, 6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[2], -a, 5e-2);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[0], a, 5e-2);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[1],
              displacement - particleEnlargementForContact, 6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[2], a, 5e-2);
}

TEST(boundingBoxCorrection, testSphereSphereBigOverlapBigBoxEnlargedParticle)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX                    = 0.1;
  T const particleEnlargementForContact = 3 * physDeltaX;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  T const                         radius = 2;
  T const     enlargedRadius = radius + particleEnlargementForContact;
  T const     displacement   = 0.5;
  PhysR<T, 3> pos1(
      0., -radius + displacement - 2 * particleEnlargementForContact, 0.);
  PhysR<T, 3> pos2(0., radius, 0.);
  T const     a = std::sqrt(enlargedRadius * enlargedRadius -
                            std::pow(enlargedRadius - 0.5 * displacement, 2));

  particles::creators::addResolvedSphere3D(particleSystem, pos1, radius, 0.,
                                           0.);
  auto particleA = particleSystem.get(0);
  particleA.setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(
      particleEnlargementForContact);
  particles::creators::addResolvedSphere3D(particleSystem, pos2, radius, 0.,
                                           0.);
  auto particleB = particleSystem.get(1);
  particleB.setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(
      particleEnlargementForContact);

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.particleContacts.push_back(
      PARTICLECONTACTTYPE(std::array<size_t, 2>({0, 1})));
  PhysR<T, 3> min = Vector<T, 3>(-a, -particleEnlargementForContact, -a);
  contactContainer.particleContacts[0].updateMinMax(min);
  PhysR<T, 3> max =
      Vector<T, 3>(a, displacement - particleEnlargementForContact, a);
  contactContainer.particleContacts[0].updateMinMax(max);

  particles::contact::particle_particle<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                        WALLCONTACTTYPE>::
      correctBoundingBoxNewContact(
          particleSystem, contactContainer.particleContacts[0], physDeltaX, 10);

  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[0], -a, 2e-2);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[1],
              -particleEnlargementForContact, 6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[2], -a, 2e-2);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[0], a, 2e-2);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[1],
              displacement - particleEnlargementForContact, 6e-4);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[2], a, 2e-2);
}

// bounding box correction for existing contacts
TEST(boundingBoxCorrection, testCuboidCuboidZeroOverlap)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  PhysR<T, 3>                     sideLength(4);
  T const                         displacement = -physDeltaX;
  PhysR<T, 3> pos1(0., -0.5 * sideLength[1] + displacement, 0.);
  PhysR<T, 3> pos2(0., 0.5 * sideLength[1], 0.);

  particles::creators::addResolvedCuboid3D(particleSystem, pos1, sideLength, 0.,
                                           0.);
  particles::creators::addResolvedCuboid3D(particleSystem, pos2, sideLength, 0.,
                                           0.);

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.particleContacts.push_back(
      PARTICLECONTACTTYPE(std::array<size_t, 2>({0, 1})));
  PhysR<T, 3> min = Vector<T, 3>(-physDeltaX);
  contactContainer.particleContacts[0].updateMinMax(min);
  PhysR<T, 3> max = Vector<T, 3>(physDeltaX);
  contactContainer.particleContacts[0].updateMinMax(max);

  particles::contact::particle_particle<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                        WALLCONTACTTYPE>::
      correctBoundingBoxExistingContact(
          particleSystem, contactContainer.particleContacts[0], physDeltaX, 10);

  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[0],
              std::numeric_limits<olb::BaseType<T>>::max(), 1e-16);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[1],
              std::numeric_limits<olb::BaseType<T>>::max(), 1e-16);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[2],
              std::numeric_limits<olb::BaseType<T>>::max(), 1e-16);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[0],
              -std::numeric_limits<olb::BaseType<T>>::max(), 1e-16);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[1],
              -std::numeric_limits<olb::BaseType<T>>::max(), 1e-16);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[2],
              -std::numeric_limits<olb::BaseType<T>>::max(), 1e-16);
}

TEST(boundingBoxCorrection, testCuboidCuboidMinimalOverlap)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  PhysR<T, 3>                     sideLength(4);
  T const                         displacement = 0.;
  PhysR<T, 3> pos1(0., -0.5 * sideLength[1] + displacement, 0.);
  PhysR<T, 3> pos2(0., 0.5 * sideLength[1], 0.);

  particles::creators::addResolvedCuboid3D(particleSystem, pos1, sideLength, 0.,
                                           0.);
  particles::creators::addResolvedCuboid3D(particleSystem, pos2, sideLength, 0.,
                                           0.);

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.particleContacts.push_back(
      PARTICLECONTACTTYPE(std::array<size_t, 2>({0, 1})));
  PhysR<T, 3> min = Vector<T, 3>(-physDeltaX);
  contactContainer.particleContacts[0].updateMinMax(min);
  PhysR<T, 3> max = Vector<T, 3>(physDeltaX);
  contactContainer.particleContacts[0].updateMinMax(max);

  particles::contact::particle_particle<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                        WALLCONTACTTYPE>::
      correctBoundingBoxExistingContact(
          particleSystem, contactContainer.particleContacts[0], physDeltaX, 10);

  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[0], -2., 1e-10);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[1], 0., 1e-10);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[2], -2., 1e-10);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[0], 2., 1e-10);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[1], 0., 1e-10);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[2], 2., 1e-10);
}

TEST(boundingBoxCorrection, testCuboidCuboidSmallOverlap)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  PhysR<T, 3>                     sideLength(4);
  T const                         displacement = physDeltaX;
  PhysR<T, 3> pos1(0., -0.5 * sideLength[1] + displacement, 0.);
  PhysR<T, 3> pos2(0., 0.5 * sideLength[1], 0.);

  particles::creators::addResolvedCuboid3D(particleSystem, pos1, sideLength, 0.,
                                           0.);
  particles::creators::addResolvedCuboid3D(particleSystem, pos2, sideLength, 0.,
                                           0.);

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.particleContacts.push_back(
      PARTICLECONTACTTYPE(std::array<size_t, 2>({0, 1})));
  PhysR<T, 3> min = Vector<T, 3>(-physDeltaX);
  contactContainer.particleContacts[0].updateMinMax(min);
  PhysR<T, 3> max = Vector<T, 3>(physDeltaX);
  contactContainer.particleContacts[0].updateMinMax(max);

  particles::contact::particle_particle<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                        WALLCONTACTTYPE>::
      correctBoundingBoxExistingContact(
          particleSystem, contactContainer.particleContacts[0], physDeltaX, 10);

  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[0], -2., 1e-10);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[1], 0., 1e-10);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMin()[2], -2., 1e-10);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[0], 2., 1e-10);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[1], physDeltaX,
              1e-10);
  EXPECT_NEAR(contactContainer.particleContacts[0].getMax()[2], 2., 1e-10);
}

TEST(boundingBoxCorrection, testCuboidWallZeroOverlap)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  PhysR<T, 3>                     sideLength(4);
  T const                         displacement = -physDeltaX;
  PhysR<T, 3> pos(0., 0.5 * sideLength[1] - displacement, 0.);

  Vector<T, PARTICLETYPE::d>                      origin(-25., -50., -25.);
  Vector<T, PARTICLETYPE::d>                      extend(50.);
  std::unique_ptr<IndicatorF<T, PARTICLETYPE::d>> wallInd =
      std::make_unique<IndicatorCuboid3D<T>>(extend, origin);
  std::vector<SolidBoundary<T, PARTICLETYPE::d>> solidBoundaries;
  solidBoundaries.push_back(SolidBoundary<T, PARTICLETYPE::d>(std::move(wallInd), 2, 1));

  particles::creators::addResolvedCuboid3D(particleSystem, pos, sideLength, 0.,
                                           0.);

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.wallContacts.push_back(WALLCONTACTTYPE(0, 0));
  PhysR<T, 3> min = Vector<T, 3>(-physDeltaX);
  contactContainer.wallContacts[0].updateMinMax(min);
  PhysR<T, 3> max = Vector<T, 3>(physDeltaX);
  contactContainer.wallContacts[0].updateMinMax(max);

  particles::contact::particle_wall<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                    WALLCONTACTTYPE>::
      correctBoundingBoxExistingContact(particleSystem, solidBoundaries,
                                        contactContainer.wallContacts[0],
                                        physDeltaX, 10);

  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[0],
              std::numeric_limits<olb::BaseType<T>>::max(), 1e-16);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[1],
              std::numeric_limits<olb::BaseType<T>>::max(), 1e-16);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[2],
              std::numeric_limits<olb::BaseType<T>>::max(), 1e-16);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[0],
              -std::numeric_limits<olb::BaseType<T>>::max(), 1e-16);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[1],
              -std::numeric_limits<olb::BaseType<T>>::max(), 1e-16);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[2],
              -std::numeric_limits<olb::BaseType<T>>::max(), 1e-16);
}

TEST(boundingBoxCorrection, testCuboidWallMinimalOverlap)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  PhysR<T, 3>                     sideLength(4);
  T const                         displacement = 0.;
  PhysR<T, 3> pos(0., 0.5 * sideLength[1] - displacement, 0.);

  Vector<T, PARTICLETYPE::d>                      origin(-25., -50., -25.);
  Vector<T, PARTICLETYPE::d>                      extend(50.);
  std::unique_ptr<IndicatorF<T, PARTICLETYPE::d>> wallInd =
      std::make_unique<IndicatorCuboid3D<T>>(extend, origin);
  std::vector<SolidBoundary<T, PARTICLETYPE::d>> solidBoundaries;
  solidBoundaries.push_back(SolidBoundary<T, PARTICLETYPE::d>(std::move(wallInd), 2, 1));

  particles::creators::addResolvedCuboid3D(particleSystem, pos, sideLength, 0.,
                                           0.);

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.wallContacts.push_back(WALLCONTACTTYPE(0, 0));
  PhysR<T, 3> min = Vector<T, 3>(-physDeltaX);
  contactContainer.wallContacts[0].updateMinMax(min);
  PhysR<T, 3> max = Vector<T, 3>(physDeltaX);
  contactContainer.wallContacts[0].updateMinMax(max);

  particles::contact::particle_wall<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                    WALLCONTACTTYPE>::
      correctBoundingBoxExistingContact(particleSystem, solidBoundaries,
                                        contactContainer.wallContacts[0],
                                        physDeltaX, 10);

  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[0], -2., 1e-10);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[1], 0., 1e-10);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[2], -2., 1e-10);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[0], 2., 1e-10);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[1], 0., 1e-10);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[2], 2., 1e-10);
}

TEST(boundingBoxCorrection, testCuboidWallSmallOverlap)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
          WALLCONTACTTYPE;
  T const physDeltaX = 0.1;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  PhysR<T, 3>                     sideLength(4);
  T const                         displacement = physDeltaX;
  PhysR<T, 3> pos(0., 0.5 * sideLength[1] - displacement, 0.);

  Vector<T, PARTICLETYPE::d>                      origin(-25., -50., -25.);
  Vector<T, PARTICLETYPE::d>                      extend(50.);
  std::unique_ptr<IndicatorF<T, PARTICLETYPE::d>> wallInd =
      std::make_unique<IndicatorCuboid3D<T>>(extend, origin);
  std::vector<SolidBoundary<T, PARTICLETYPE::d>> solidBoundaries;
  solidBoundaries.push_back(SolidBoundary<T, PARTICLETYPE::d>(std::move(wallInd), 2, 1));

  particles::creators::addResolvedCuboid3D(particleSystem, pos, sideLength, 0.,
                                           0.);

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  contactContainer.wallContacts.push_back(WALLCONTACTTYPE(0, 0));
  PhysR<T, 3> min = Vector<T, 3>(-physDeltaX);
  contactContainer.wallContacts[0].updateMinMax(min);
  PhysR<T, 3> max = Vector<T, 3>(physDeltaX);
  contactContainer.wallContacts[0].updateMinMax(max);

  particles::contact::particle_wall<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                    WALLCONTACTTYPE>::
      correctBoundingBoxExistingContact(particleSystem, solidBoundaries,
                                        contactContainer.wallContacts[0],
                                        physDeltaX, 10);

  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[0], -2., 1e-10);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[1], -physDeltaX, 1e-10);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMin()[2], -2., 1e-10);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[0], 2., 1e-10);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[1], 0., 1e-10);
  EXPECT_NEAR(contactContainer.wallContacts[0].getMax()[2], 2., 1e-10);
}

// Normal force calculation
TEST(normalForce3D, testZeroOverlap)
{
  constexpr unsigned D = 3;

  const Vector<T, D> contactNormal(1);
  const Vector<T, D> relVel(1);
  const T            indentation   = 1;
  const T            overlapVolume = 0;
  const T            effectiveE    = 1;
  const T            dampingFactor = 1;
  const T            k             = 1;

  Vector<T, D> normalForce =
      particles::contact::getNormalForceFromOverlapVolume(
          contactNormal, indentation, overlapVolume, effectiveE, dampingFactor,
          k, relVel);

  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_EQ(0., normalForce[iD]);
  }
}

TEST(normalForce3D, testZeroIdentation)
{
  constexpr unsigned D = 3;

  const Vector<T, D> contactNormal(1);
  const Vector<T, D> relVel(1);
  const T            indentation   = 0;
  const T            overlapVolume = 1;
  const T            effectiveE    = 1;
  const T            dampingFactor = 1;
  const T            k             = 1;

  Vector<T, D> normalForce =
      particles::contact::getNormalForceFromOverlapVolume(
          contactNormal, indentation, overlapVolume, effectiveE, dampingFactor,
          k, relVel);

  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_EQ(0., normalForce[iD]);
  }
}

TEST(normalForce3D, testNoDamping1)
{
  constexpr unsigned D = 3;

  const Vector<T, D> contactNormal(normalize(Vector<T, D>(1)));
  const Vector<T, D> relVel(1);
  const T            indentation   = 1;
  const T            overlapVolume = 1;
  const T            effectiveE    = 1;
  const T            dampingFactor = 0;
  const T            k             = 1;

  Vector<T, D> normalForce =
      particles::contact::getNormalForceFromOverlapVolume(
          contactNormal, indentation, overlapVolume, effectiveE, dampingFactor,
          k, relVel);

  EXPECT_NEAR(norm(normalForce), 1, 1e-16);
  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_NEAR(normalize(Vector<T, D>(1))[iD], normalForce[iD], 1e-16);
  }
}

TEST(normalForce3D, testNoDamping2)
{
  constexpr unsigned D = 3;

  const Vector<T, D> contactNormal(normalize(Vector<T, D>(1)));
  const Vector<T, D> relVel(2);
  const T            indentation   = 2;
  const T            overlapVolume = 2;
  const T            effectiveE    = 2;
  const T            dampingFactor = 0;
  const T            k             = 2;

  Vector<T, D> normalForce =
      particles::contact::getNormalForceFromOverlapVolume(
          contactNormal, indentation, overlapVolume, effectiveE, dampingFactor,
          k, relVel);

  EXPECT_NEAR(norm(normalForce), 8, 1e-16);
  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_NEAR(8 * normalize(Vector<T, D>(1))[iD], normalForce[iD], 1e-16);
  }
}
TEST(normalForce3D, testNoDamping3)
{
  constexpr unsigned D = 3;

  const Vector<T, D> contactNormal(normalize(Vector<T, D>(-2)));
  const Vector<T, D> relVel(-2);
  const T            indentation   = -2;
  const T            overlapVolume = -2;
  const T            effectiveE    = -2;
  const T            dampingFactor = 0;
  const T            k             = -2;

  Vector<T, D> normalForce =
      particles::contact::getNormalForceFromOverlapVolume(
          contactNormal, indentation, overlapVolume, effectiveE, dampingFactor,
          k, relVel);

  EXPECT_NEAR(norm(normalForce), 8, 1e-16);
  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_NEAR(-8 * normalize(Vector<T, D>(1))[iD], normalForce[iD], 1e-16);
  }
}
TEST(normalForce3D, testWithDamping1)
{
  constexpr unsigned D = 3;

  const Vector<T, D> contactNormal(normalize(Vector<T, D>(1)));
  const Vector<T, D> relVel(-1);
  const T            indentation   = 1;
  const T            overlapVolume = 1;
  const T            effectiveE    = 1;
  const T            dampingFactor = 1;
  const T            k             = 1;

  Vector<T, D> normalForce =
      particles::contact::getNormalForceFromOverlapVolume(
          contactNormal, indentation, overlapVolume, effectiveE, dampingFactor,
          k, relVel);

  T const expectedMagnitude = 1 + util::sqrt(3);
  EXPECT_NEAR(norm(normalForce), expectedMagnitude, 1e-15);
  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_NEAR(expectedMagnitude * normalize(Vector<T, D>(1))[iD],
                normalForce[iD], 1e-15);
  }
}

TEST(normalForce3D, testWithDamping2)
{
  constexpr unsigned D = 3;

  const Vector<T, D> contactNormal(normalize(Vector<T, D>(-1)));
  const Vector<T, D> relVel(2);
  const T            indentation   = 2;
  const T            overlapVolume = 2;
  const T            effectiveE    = 2;
  const T            dampingFactor = 2;
  const T            k             = 2;

  Vector<T, D> normalForce =
      particles::contact::getNormalForceFromOverlapVolume(
          contactNormal, indentation, overlapVolume, effectiveE, dampingFactor,
          k, relVel);

  T const expectedMagnitude = 8 * (1 + 2 * util::sqrt(12));
  EXPECT_NEAR(norm(normalForce), expectedMagnitude, 1e-12);
  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_NEAR(-expectedMagnitude * normalize(Vector<T, D>(1))[iD],
                normalForce[iD], 1e-12);
  }
}

TEST(normalForce3D, testWithDamping3)
{
  constexpr unsigned D = 3;

  const Vector<T, D> contactNormal(normalize(Vector<T, D>(1)));
  const Vector<T, D> relVel(-2);
  const T            indentation   = -2;
  const T            overlapVolume = -2;
  const T            effectiveE    = -2;
  const T            dampingFactor = -2;
  const T            k             = -2;

  Vector<T, D> normalForce =
      particles::contact::getNormalForceFromOverlapVolume(
          contactNormal, indentation, overlapVolume, effectiveE, dampingFactor,
          k, relVel);

  // Negative magnitude is not allowed
  // T const expectedMagnitude = 8 * (1 - 2 * util::sqrt(12));
  // EXPECT_NEAR(-norm(normalForce), expectedMagnitude, 1e-12);
  EXPECT_NEAR(norm(normalForce), T {0}, 1e-12);
  for (unsigned iD = 0; iD < D; ++iD) {
    //EXPECT_NEAR(expectedMagnitude * normalize(Vector<T, D>(1))[iD],
    //            normalForce[iD], 1e-12);
    EXPECT_NEAR(T {0}, normalForce[iD], 1e-12);
  }
}
// Tangential force
TEST(tangentialForce3D, testWithoutFriction)
{
  constexpr unsigned D = 3;

  const T            coefficientStaticFriction       = 0;
  const T            coefficientKineticFriction      = 0;
  const T            staticKineticTransitionVelocity = 1;
  const T            k                               = 1;
  const Vector<T, D> relVel(1);
  const Vector<T, D> normalForce(1);
  const Vector<T, D> contactNormal(normalize(Vector<T, D>(1.)));

  const Vector<T, D> expectedTangentialForce(T {0});
  const Vector<T, D> tangentialForce =
      particles::contact::getTangentialForceFromOverlapVolume(
          contactNormal, coefficientStaticFriction, coefficientKineticFriction,
          staticKineticTransitionVelocity, k, relVel, normalForce);

  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_NEAR(tangentialForce[iD], expectedTangentialForce[iD], 1e-16);
  }
}

TEST(tangentialForce3D, testWithFriction1)
{
  constexpr unsigned D = 3;

  const T            coefficientStaticFriction       = 0.6;
  const T            coefficientKineticFriction      = 0.3;
  const T            staticKineticTransitionVelocity = 0.01;
  const T            k                               = 1;
  const Vector<T, D> relVel(0.1);
  const Vector<T, D> normalForce(1);
  const Vector<T, D> contactNormal(normalize(Vector<T, D>(0.23, 0.56, 0.18)));

  const T            expectedTangentialForce = 0.535757461650446;
  const Vector<T, D> expectedTangentialForceDirection(
      -0.550346973488248, 0.451735314416122, -0.702177623170728);
  const Vector<T, D> tangentialForce =
      particles::contact::getTangentialForceFromOverlapVolume(
          contactNormal, coefficientStaticFriction, coefficientKineticFriction,
          staticKineticTransitionVelocity, k, relVel, normalForce);

  EXPECT_NEAR(norm(tangentialForce), expectedTangentialForce, 1e-14);
  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_NEAR(tangentialForce[iD],
                expectedTangentialForce * expectedTangentialForceDirection[iD],
                1e-14);
  }
}

TEST(tangentialForce3D, testWithFriction2)
{
  constexpr unsigned D = 3;

  const T            coefficientStaticFriction       = 0.78;
  const T            coefficientKineticFriction      = 0.52;
  const T            staticKineticTransitionVelocity = 0.05;
  const T            k                               = 1;
  const Vector<T, D> relVel(0.25);
  const Vector<T, D> normalForce(1.3);
  const Vector<T, D> contactNormal(
      normalize(Vector<T, D>(-0.37, -0.96, -0.19)));

  const T            expectedTangentialForce = 1.22278358391258;
  const Vector<T, D> expectedTangentialForceDirection(
      -0.515586428764798, 0.353218361382272, -0.780645517284244);
  const Vector<T, D> tangentialForce =
      particles::contact::getTangentialForceFromOverlapVolume(
          contactNormal, coefficientStaticFriction, coefficientKineticFriction,
          staticKineticTransitionVelocity, k, relVel, normalForce);

  EXPECT_NEAR(norm(tangentialForce), expectedTangentialForce, 1e-14);
  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_NEAR(tangentialForce[iD],
                expectedTangentialForce * expectedTangentialForceDirection[iD],
                1e-14);
  }
}

// Force application

TEST(contactForceApplication3DsingleParticle, testWithoutForces)
{
  typedef ResolvedSphere3D PARTICLETYPE;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  const T                         radius = 2;
  const PhysR<T, 3>               pos(5.);
  const PhysR<T, 3>               contactPoint(1);
  const Vector<T, 3>              normalForce(0);
  const Vector<T, 3>              tangentialForce(0);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, 0., 0.);
  auto particle = particleSystem.get(0);

  particles::contact::applyForceOnParticle(dataMap, particle, particleSystem,
                                           contactPoint, normalForce,
                                           tangentialForce);

  Vector<T, 3> force  = particle.template getField<FORCING, FORCE>();
  Vector<T, 3> torque = particle.template getField<FORCING, TORQUE>();
  for (unsigned iD = 0; iD < 3; ++iD) {
    EXPECT_EQ(0., force[iD]);
    EXPECT_EQ(0., torque[iD]);
  }
}

TEST(contactForceApplication3DsingleParticle, testOnlyNormalForce1)
{
  typedef ResolvedSphere3D PARTICLETYPE;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  const T                         radius = 2;
  const PhysR<T, 3>               pos(5.);
  const PhysR<T, 3>               contactPoint(5.);
  const Vector<T, 3>              normalForce(5.);
  const Vector<T, 3>              tangentialForce(0);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, 0., 0.);
  auto particle = particleSystem.get(0);

  particles::contact::applyForceOnParticle(dataMap, particle, particleSystem,
                                           contactPoint, normalForce,
                                           tangentialForce);

  Vector<T, 3> force  = particle.template getField<FORCING, FORCE>();
  Vector<T, 3> torque = particle.template getField<FORCING, TORQUE>();
  for (unsigned iD = 0; iD < 3; ++iD) {
    EXPECT_EQ(5., force[iD]);
    EXPECT_EQ(0., torque[iD]);
  }
}

TEST(contactForceApplication3DsingleParticle, testOnlyNormalForce2)
{
  typedef ResolvedSphere3D PARTICLETYPE;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  const T                         radius = 2;
  const PhysR<T, 3>               pos(5.);
  const PhysR<T, 3>               contactPoint(.5);
  const Vector<T, 3>              normalForce(2.);
  const Vector<T, 3>              tangentialForce(0);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, 0., 0.);
  auto particle = particleSystem.get(0);

  particles::contact::applyForceOnParticle(dataMap, particle, particleSystem,
                                           contactPoint, normalForce,
                                           tangentialForce);

  Vector<T, 3> force  = particle.template getField<FORCING, FORCE>();
  Vector<T, 3> torque = particle.template getField<FORCING, TORQUE>();
  for (unsigned iD = 0; iD < 3; ++iD) {
    EXPECT_EQ(2., force[iD]);
    EXPECT_EQ(0., torque[iD]);
  }
}

TEST(contactForceApplication3DsingleParticle, testOnlyNormalForce3)
{
  typedef ResolvedSphere3D PARTICLETYPE;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  const T                         radius = 2;
  const PhysR<T, 3>               pos(5.);
  const PhysR<T, 3>               contactPoint(2.);
  const Vector<T, 3>              normalForce(-1, -2., -3);
  const Vector<T, 3>              tangentialForce(0);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, 0., 0.);
  auto particle = particleSystem.get(0);

  particles::contact::applyForceOnParticle(dataMap, particle, particleSystem,
                                           contactPoint, normalForce,
                                           tangentialForce);
  const Vector<T, 3> expectedTorque(3, -6, 3);

  Vector<T, 3> force  = particle.template getField<FORCING, FORCE>();
  Vector<T, 3> torque = particle.template getField<FORCING, TORQUE>();
  for (unsigned iD = 0; iD < 3; ++iD) {
    EXPECT_EQ(normalForce[iD], force[iD]);
    EXPECT_EQ(expectedTorque[iD], torque[iD]);
  }
}

TEST(contactForceApplication3DsingleParticle, testOnlyTangentialForce1)
{
  typedef ResolvedSphere3D PARTICLETYPE;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  const T                         radius = 2;
  const PhysR<T, 3>               pos(5.);
  const PhysR<T, 3>               contactPoint(5.);
  const Vector<T, 3>              normalForce(0);
  const Vector<T, 3>              tangentialForce(5);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, 0., 0.);
  auto particle = particleSystem.get(0);

  particles::contact::applyForceOnParticle(dataMap, particle, particleSystem,
                                           contactPoint, normalForce,
                                           tangentialForce);

  Vector<T, 3> force  = particle.template getField<FORCING, FORCE>();
  Vector<T, 3> torque = particle.template getField<FORCING, TORQUE>();
  for (unsigned iD = 0; iD < 3; ++iD) {
    EXPECT_EQ(5., force[iD]);
    EXPECT_EQ(0., torque[iD]);
  }
}

TEST(contactForceApplication3DsingleParticle, testOnlyTangentialForce2)
{
  typedef ResolvedSphere3D PARTICLETYPE;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  const T                         radius = 2;
  const PhysR<T, 3>               pos(5.);
  const PhysR<T, 3>               contactPoint(.5);
  const Vector<T, 3>              normalForce(0);
  const Vector<T, 3>              tangentialForce(2.);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, 0., 0.);
  auto particle = particleSystem.get(0);

  particles::contact::applyForceOnParticle(dataMap, particle, particleSystem,
                                           contactPoint, normalForce,
                                           tangentialForce);

  Vector<T, 3> force  = particle.template getField<FORCING, FORCE>();
  Vector<T, 3> torque = particle.template getField<FORCING, TORQUE>();
  for (unsigned iD = 0; iD < 3; ++iD) {
    EXPECT_EQ(2., force[iD]);
    EXPECT_EQ(0., torque[iD]);
  }
}

TEST(contactForceApplication3DsingleParticle, testOnlyTangentialForce3)
{
  typedef ResolvedSphere3D PARTICLETYPE;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  const T                         radius = 2;
  const PhysR<T, 3>               pos(5.);
  const PhysR<T, 3>               contactPoint(2.);
  const Vector<T, 3>              normalForce(0);
  const Vector<T, 3>              tangentialForce(-1, -2, -3);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, 0., 0.);
  auto particle = particleSystem.get(0);

  particles::contact::applyForceOnParticle(dataMap, particle, particleSystem,
                                           contactPoint, normalForce,
                                           tangentialForce);
  const Vector<T, 3> expectedTorque(3, -6, 3);

  Vector<T, 3> force  = particle.template getField<FORCING, FORCE>();
  Vector<T, 3> torque = particle.template getField<FORCING, TORQUE>();
  for (unsigned iD = 0; iD < 3; ++iD) {
    EXPECT_EQ(tangentialForce[iD], force[iD]);
    EXPECT_EQ(expectedTorque[iD], torque[iD]);
  }
}

TEST(contactForceApplication3DsingleParticle, testNormalAndTangentialForce)
{
  typedef ResolvedSphere3D PARTICLETYPE;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  const T                         radius = 2;
  const PhysR<T, 3>               pos(5.);
  const PhysR<T, 3>               contactPoint(8.5);
  const Vector<T, 3>              normalForce(3, 4, 1);
  const Vector<T, 3>              tangentialForce(-1, -2, -3);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, 0., 0.);
  auto particle = particleSystem.get(0);

  particles::contact::applyForceOnParticle(dataMap, particle, particleSystem,
                                           contactPoint, normalForce,
                                           tangentialForce);
  const Vector<T, 3> expectedTorque(-14, 14, 0);

  Vector<T, 3> force  = particle.template getField<FORCING, FORCE>();
  Vector<T, 3> torque = particle.template getField<FORCING, TORQUE>();
  for (unsigned iD = 0; iD < 3; ++iD) {
    EXPECT_EQ((normalForce + tangentialForce)[iD], force[iD]);
    EXPECT_EQ(expectedTorque[iD], torque[iD]);
  }
}

TEST(contactForceApplication3DtwoParticles, testWithZeroOverlap)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ContactProperties<T, 1>       CONTACTPROPERTIES;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  const T                         radius = 2;
  const PhysR<T, 3>               pos(5.);
  const PhysR<T, 3>               contactPoint(1);
  const Vector<T, 3>              contactNormal(normalize(Vector<T, 3>(1)));
  const T                         indentation                     = 1;
  const T                         overlapVolume                   = 0;
  const T                         dampingFactor                   = 0.5;
  const T                         coefficientKineticFriction      = 1;
  const T                         coefficientStaticFriction       = 1;
  const T                         staticKineticTransitionVelocity = 1;
  const T                         k                               = 1;

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T effectiveE =
      particles::contact::evalEffectiveYoungModulus(1., 1., 0.5, 0.5);
  CONTACTPROPERTIES contactProperties;
  contactProperties.set(0, 0, effectiveE, dampingFactor,
                        coefficientKineticFriction, coefficientStaticFriction,
                        staticKineticTransitionVelocity);

  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, 0., 0.);
  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, 0., 0.);
  auto particleA = particleSystem.get(0);
  particleA.template setField<MECHPROPERTIES, MATERIAL>(0);
  particleA.template setField<MOBILITY, VELOCITY>(Vector<T, 3>(0.25));
  auto particleB = particleSystem.get(1);
  particleB.template setField<MECHPROPERTIES, MATERIAL>(0);
  particleB.template setField<MOBILITY, VELOCITY>(Vector<T, 3>(0.25));

  particles::contact::applyForceFromOverlapVolume(
      dataMap, particleA, particleB, particleSystem, contactPoint,
      contactNormal, indentation, overlapVolume,
      particles::contact::evalRelativeVelocity(particleA, particleB,
                                               contactPoint),
      dampingFactor, contactProperties, k, 0, 0);

  Vector<T, 3> forceA  = particleA.template getField<FORCING, FORCE>();
  Vector<T, 3> torqueA = particleA.template getField<FORCING, TORQUE>();
  for (unsigned iD = 0; iD < 3; ++iD) {
    EXPECT_EQ(0., forceA[iD]);
    EXPECT_EQ(0., torqueA[iD]);
  }

  Vector<T, 3> forceB  = particleB.template getField<FORCING, FORCE>();
  Vector<T, 3> torqueB = particleB.template getField<FORCING, TORQUE>();
  for (unsigned iD = 0; iD < 3; ++iD) {
    EXPECT_EQ(0., forceB[iD]);
    EXPECT_EQ(0., torqueB[iD]);
  }
}

TEST(contactForceApplication3DtwoParticles, testWithZeroIndentation)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ContactProperties<T, 1>       CONTACTPROPERTIES;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  const T                         radius = 2;
  const PhysR<T, 3>               pos(5.);
  const PhysR<T, 3>               contactPoint(1);
  const Vector<T, 3>              contactNormal(normalize(Vector<T, 3>(1)));
  const T                         indentation                     = 0;
  const T                         overlapVolume                   = 1;
  const T                         dampingFactor                   = 0.5;
  const T                         coefficientKineticFriction      = 1;
  const T                         coefficientStaticFriction       = 1;
  const T                         staticKineticTransitionVelocity = 1;
  const T                         k                               = 1;

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T effectiveE =
      particles::contact::evalEffectiveYoungModulus(1., 1., 0.5, 0.5);
  CONTACTPROPERTIES contactProperties;
  contactProperties.set(0, 0, effectiveE, dampingFactor,
                        coefficientKineticFriction, coefficientStaticFriction,
                        staticKineticTransitionVelocity);

  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, 0., 0.);
  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, 0., 0.);
  auto particleA = particleSystem.get(0);
  particleA.template setField<MECHPROPERTIES, MATERIAL>(0);
  particleA.template setField<MOBILITY, VELOCITY>(Vector<T, 3>(0.25));
  auto particleB = particleSystem.get(1);
  particleB.template setField<MECHPROPERTIES, MATERIAL>(0);
  particleB.template setField<MOBILITY, VELOCITY>(Vector<T, 3>(0.25));

  particles::contact::applyForceFromOverlapVolume(
      dataMap, particleA, particleB, particleSystem, contactPoint,
      contactNormal, indentation, overlapVolume,
      particles::contact::evalRelativeVelocity(particleA, particleB,
                                               contactPoint),
      dampingFactor, contactProperties, k, 0, 0);

  Vector<T, 3> forceA  = particleA.template getField<FORCING, FORCE>();
  Vector<T, 3> torqueA = particleA.template getField<FORCING, TORQUE>();
  for (unsigned iD = 0; iD < 3; ++iD) {
    EXPECT_EQ(0., forceA[iD]);
    EXPECT_EQ(0., torqueA[iD]);
  }

  Vector<T, 3> forceB  = particleB.template getField<FORCING, FORCE>();
  Vector<T, 3> torqueB = particleB.template getField<FORCING, TORQUE>();
  for (unsigned iD = 0; iD < 3; ++iD) {
    EXPECT_EQ(0., forceB[iD]);
    EXPECT_EQ(0., torqueB[iD]);
  }
}

TEST(contactForceApplication3DtwoParticles, testTwoSpheres)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ContactProperties<T, 1>       CONTACTPROPERTIES;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  const T                         radius = 2;
  // The position doesn't influence the force, only the torque because of the cross product, all other values are fixed here and would be different if the following position was true
  const PhysR<T, 3>  pos(5.);
  const PhysR<T, 3>  contactPoint(1);
  const Vector<T, 3> contactNormal(normalize(Vector<T, 3>(-1)));
  const T            indentation                     = 0.5;
  const T            overlapVolume                   = 2.5;
  const T            dampingFactor                   = 0.5;
  const T            coefficientKineticFriction      = 0.3;
  const T            coefficientStaticFriction       = 0.6;
  const T            staticKineticTransitionVelocity = 0.01;
  const T            k                               = 0.9;

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T effectiveE =
      particles::contact::evalEffectiveYoungModulus(1., 1., 0.5, 0.5);
  CONTACTPROPERTIES contactProperties;
  contactProperties.set(0, 0, effectiveE, dampingFactor,
                        coefficientKineticFriction, coefficientStaticFriction,
                        staticKineticTransitionVelocity);

  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, 0., 0.);
  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, 0., 0.);
  auto particleA = particleSystem.get(0);
  particleA.template setField<MECHPROPERTIES, MATERIAL>(0);
  particleA.template setField<MOBILITY, VELOCITY>(Vector<T, 3>(0.25));
  auto particleB = particleSystem.get(1);
  particleB.template setField<MECHPROPERTIES, MATERIAL>(0);
  particleB.template setField<MOBILITY, VELOCITY>(Vector<T, 3>(0.25));
  particles::contact::applyForceFromOverlapVolume(
      dataMap, particleA, particleB, particleSystem, contactPoint,
      contactNormal, indentation, overlapVolume,
      particles::contact::evalRelativeVelocity(particleA, particleB,
                                               contactPoint),
      dampingFactor, contactProperties, k, 0, 0);

  Vector<T, 3> relVel =
      particles::dynamics::calculateLocalVelocity(particleA, contactPoint) -
      particles::dynamics::calculateLocalVelocity(particleB, contactPoint);
  const Vector<T, 3> normalForce =
      particles::contact::getNormalForceFromOverlapVolume(
          contactNormal, indentation, overlapVolume, effectiveE, dampingFactor,
          k, relVel);
  const Vector<T, 3> tangentialForce =
      particles::contact::getTangentialForceFromOverlapVolume(
          contactNormal, coefficientStaticFriction, coefficientKineticFriction,
          staticKineticTransitionVelocity, k, relVel, normalForce);
  const Vector<T, 3> expectedForce = normalForce + tangentialForce;
  const Vector<T, 3> expectedTorque =
      crossProduct3D(contactPoint - pos, expectedForce);

  Vector<T, 3> forceA  = particleA.template getField<FORCING, FORCE>();
  Vector<T, 3> torqueA = particleA.template getField<FORCING, TORQUE>();
  for (unsigned iD = 0; iD < 3; ++iD) {
    EXPECT_NEAR(expectedForce[iD], forceA[iD], 1e-16);
    EXPECT_NEAR(expectedTorque[iD], torqueA[iD], 1e-16);
  }

  Vector<T, 3> forceB  = particleB.template getField<FORCING, FORCE>();
  Vector<T, 3> torqueB = particleB.template getField<FORCING, TORQUE>();
  for (unsigned iD = 0; iD < 3; ++iD) {
    EXPECT_NEAR(-expectedForce[iD], forceB[iD], 1e-16);
    EXPECT_NEAR(-expectedTorque[iD], torqueB[iD], 1e-16);
  }
}

TEST(contactForceApplication3DParticleWall, testWithZeroOverlap)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ContactProperties<T, 1>       CONTACTPROPERTIES;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  const T                         radius = 2;
  const PhysR<T, 3>               pos(5.);
  const PhysR<T, 3>               contactPoint(1);
  const Vector<T, 3>              contactNormal(normalize(Vector<T, 3>(1)));
  const T                         indentation                     = 1;
  const T                         overlapVolume                   = 0;
  const T                         dampingFactor                   = 0.5;
  const T                         coefficientKineticFriction      = 1;
  const T                         coefficientStaticFriction       = 1;
  const T                         staticKineticTransitionVelocity = 1;
  const T                         k                               = 1;

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  // Pure dummy wall, only modulus of elasticity and poisson's ratio matter
  std::unique_ptr<IndicatorF<T, PARTICLETYPE::d>> wallInd =
      std::make_unique<IndicatorCuboid3D<T>>(Vector<T, PARTICLETYPE::d>(1),
                                             Vector<T, PARTICLETYPE::d>(0));
  SolidBoundary<T, PARTICLETYPE::d> wall(std::move(wallInd), 2, 0);
  const T                           effectiveE =
      particles::contact::evalEffectiveYoungModulus(1., 1., 0.5, 0.5);

  CONTACTPROPERTIES contactProperties;
  contactProperties.set(0, 0, effectiveE, dampingFactor,
                        coefficientKineticFriction, coefficientStaticFriction,
                        staticKineticTransitionVelocity);

  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, 0., 0.);
  auto particle = particleSystem.get(0);
  particle.template setField<MECHPROPERTIES, MATERIAL>(0);
  particle.template setField<MECHPROPERTIES, POISSON_RATIO>(0.5);
  particle.template setField<MOBILITY, VELOCITY>(Vector<T, 3>(0.25));

  particles::contact::applyForceFromOverlapVolume(
      dataMap, particle, particleSystem, wall, contactPoint, contactNormal,
      indentation, overlapVolume,
      particles::contact::evalRelativeVelocity(particle, contactPoint),
      dampingFactor, contactProperties, k, 0, 0);

  Vector<T, 3> force  = particle.template getField<FORCING, FORCE>();
  Vector<T, 3> torque = particle.template getField<FORCING, TORQUE>();
  for (unsigned iD = 0; iD < 3; ++iD) {
    EXPECT_EQ(0., force[iD]);
    EXPECT_EQ(0., torque[iD]);
  }
}

TEST(contactForceApplication3DParticleWall, testWithZeroIndentation)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ContactProperties<T, 1>       CONTACTPROPERTIES;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  const T                         radius = 2;
  const PhysR<T, 3>               pos(5.);
  const PhysR<T, 3>               contactPoint(1);
  const Vector<T, 3>              contactNormal(normalize(Vector<T, 3>(1)));
  const T                         indentation                     = 0;
  const T                         overlapVolume                   = 1;
  const T                         dampingFactor                   = 0.5;
  const T                         coefficientKineticFriction      = 1;
  const T                         coefficientStaticFriction       = 1;
  const T                         staticKineticTransitionVelocity = 1;
  const T                         k                               = 1;

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  // Pure dummy wall, only modulus of elasticity and poisson's ratio matter
  std::unique_ptr<IndicatorF<T, PARTICLETYPE::d>> wallInd =
      std::make_unique<IndicatorCuboid3D<T>>(Vector<T, PARTICLETYPE::d>(1),
                                             Vector<T, PARTICLETYPE::d>(0));
  SolidBoundary<T, PARTICLETYPE::d> wall(std::move(wallInd), 2, 0);
  const T                           effectiveE =
      particles::contact::evalEffectiveYoungModulus(1., 1., 0.5, 0.5);

  CONTACTPROPERTIES contactProperties;
  contactProperties.set(0, 0, effectiveE, dampingFactor,
                        coefficientKineticFriction, coefficientStaticFriction,
                        staticKineticTransitionVelocity);

  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, 0., 0.);
  auto particle = particleSystem.get(0);
  particle.template setField<MECHPROPERTIES, MATERIAL>(0);
  particle.template setField<MOBILITY, VELOCITY>(Vector<T, 3>(0.25));

  particles::contact::applyForceFromOverlapVolume(
      dataMap, particle, particleSystem, wall, contactPoint, contactNormal,
      indentation, overlapVolume,
      particles::contact::evalRelativeVelocity(particle, contactPoint),
      dampingFactor, contactProperties, k, 0, 0);

  Vector<T, 3> force  = particle.template getField<FORCING, FORCE>();
  Vector<T, 3> torque = particle.template getField<FORCING, TORQUE>();
  for (unsigned iD = 0; iD < 3; ++iD) {
    EXPECT_EQ(0., force[iD]);
    EXPECT_EQ(0., torque[iD]);
  }
}

TEST(contactForceApplication3DParticleWall, testMovingSphere)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ContactProperties<T, 1>       CONTACTPROPERTIES;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  const T                         radius = 2;
  // The position doesn't influence the force, only the torque because of the cross product, all other values are fixed here and would be different if the following position was true
  const PhysR<T, 3>  pos(5.);
  const PhysR<T, 3>  contactPoint(1);
  const Vector<T, 3> contactNormal(normalize(Vector<T, 3>(-1)));
  const T            indentation                     = 0.5;
  const T            overlapVolume                   = 2.5;
  const T            dampingFactor                   = 0.5;
  const T            coefficientKineticFriction      = 0.3;
  const T            coefficientStaticFriction       = 0.6;
  const T            staticKineticTransitionVelocity = 0.01;
  const T            k                               = 0.9;

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  // Pure dummy wall, only modulus of elasticity and poisson's ratio matter
  std::unique_ptr<IndicatorF<T, PARTICLETYPE::d>> wallInd =
      std::make_unique<IndicatorCuboid3D<T>>(Vector<T, PARTICLETYPE::d>(1),
                                             Vector<T, PARTICLETYPE::d>(0));
  SolidBoundary<T, PARTICLETYPE::d> wall(std::move(wallInd), 2, 0);
  const T                           effectiveE =
      particles::contact::evalEffectiveYoungModulus(1., 1., 0.5, 0.5);

  CONTACTPROPERTIES contactProperties;
  contactProperties.set(0, 0, effectiveE, dampingFactor,
                        coefficientKineticFriction, coefficientStaticFriction,
                        staticKineticTransitionVelocity);

  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, 0., 0.);
  auto particle = particleSystem.get(0);
  particle.template setField<MECHPROPERTIES, MATERIAL>(0);
  particle.template setField<MOBILITY, VELOCITY>(Vector<T, 3>(0.25));
  particles::contact::applyForceFromOverlapVolume(
      dataMap, particle, particleSystem, wall, contactPoint, contactNormal,
      indentation, overlapVolume,
      particles::contact::evalRelativeVelocity(particle, contactPoint),
      dampingFactor, contactProperties, k, 0, 0);

  Vector<T, 3> relVel =
      particles::dynamics::calculateLocalVelocity(particle, contactPoint);
  const Vector<T, 3> normalForce =
      particles::contact::getNormalForceFromOverlapVolume(
          contactNormal, indentation, overlapVolume, effectiveE, dampingFactor,
          k, relVel);
  const Vector<T, 3> tangentialForce =
      particles::contact::getTangentialForceFromOverlapVolume(
          contactNormal, coefficientStaticFriction, coefficientKineticFriction,
          staticKineticTransitionVelocity, k, relVel, normalForce);
  const Vector<T, 3> expectedForce = normalForce + tangentialForce;
  const Vector<T, 3> expectedTorque =
      crossProduct3D(contactPoint - pos, expectedForce);

  Vector<T, 3> force  = particle.template getField<FORCING, FORCE>();
  Vector<T, 3> torque = particle.template getField<FORCING, TORQUE>();
  for (unsigned iD = 0; iD < 3; ++iD) {
    EXPECT_NEAR(expectedForce[iD], force[iD], 1e-16);
    EXPECT_NEAR(expectedTorque[iD], torque[iD], 1e-16);
  }
}

TEST(contactForceApplication3DParticleWall, testRestingSphere)
{
  typedef ResolvedParticleWithContact3D PARTICLETYPE;
  typedef ContactProperties<T, 1>       CONTACTPROPERTIES;

  ParticleSystem<T, PARTICLETYPE> particleSystem;
  const T                         radius = 2;
  // The position doesn't influence the force, only the torque because of the cross product, all other values are fixed here and would be different if the following position was true
  const PhysR<T, 3>  pos(5.);
  const PhysR<T, 3>  contactPoint(1);
  const Vector<T, 3> contactNormal(normalize(Vector<T, 3>(-1)));
  const T            indentation                     = 0.5;
  const T            overlapVolume                   = 2.5;
  const T            dampingFactor                   = 0.5;
  const T            coefficientKineticFriction      = 0.3;
  const T            coefficientStaticFriction       = 0.6;
  const T            staticKineticTransitionVelocity = 0.01;
  const T            k                               = 0.9;

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  // Pure dummy wall, only modulus of elasticity and poisson's ratio matter
  std::unique_ptr<IndicatorF<T, PARTICLETYPE::d>> wallInd =
      std::make_unique<IndicatorCuboid3D<T>>(Vector<T, PARTICLETYPE::d>(1),
                                             Vector<T, PARTICLETYPE::d>(0));
  SolidBoundary<T, PARTICLETYPE::d> wall(std::move(wallInd), 2, 0);
  const T                           effectiveE =
      particles::contact::evalEffectiveYoungModulus(1., 1., 0.5, 0.5);

  CONTACTPROPERTIES contactProperties;
  contactProperties.set(0, 0, effectiveE, dampingFactor,
                        coefficientKineticFriction, coefficientStaticFriction,
                        staticKineticTransitionVelocity);

  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, 0., 0.);
  auto particle = particleSystem.get(0);
  particle.template setField<MECHPROPERTIES, MATERIAL>(0);
  particles::contact::applyForceFromOverlapVolume(
      dataMap, particle, particleSystem, wall, contactPoint, contactNormal,
      indentation, overlapVolume,
      particles::contact::evalRelativeVelocity(particle, contactPoint),
      dampingFactor, contactProperties, k, 0, 0);

  Vector<T, 3> relVel =
      particles::dynamics::calculateLocalVelocity(particle, contactPoint);
  const Vector<T, 3> normalForce =
      particles::contact::getNormalForceFromOverlapVolume(
          contactNormal, indentation, overlapVolume, effectiveE, dampingFactor,
          k, relVel);
  const Vector<T, 3> tangentialForce =
      particles::contact::getTangentialForceFromOverlapVolume(
          contactNormal, coefficientStaticFriction, coefficientKineticFriction,
          staticKineticTransitionVelocity, k, relVel, normalForce);
  const Vector<T, 3> expectedForce = normalForce + tangentialForce;
  const Vector<T, 3> expectedTorque =
      crossProduct3D(contactPoint - pos, expectedForce);

  Vector<T, 3> force  = particle.template getField<FORCING, FORCE>();
  Vector<T, 3> torque = particle.template getField<FORCING, TORQUE>();
  for (unsigned iD = 0; iD < 3; ++iD) {
    EXPECT_NEAR(expectedForce[iD], force[iD], 1e-16);
    EXPECT_NEAR(expectedTorque[iD], torque[iD], 1e-16);
  }
}

// "Process"-methods

TEST(cellProcessing, testInside)
{
  constexpr unsigned D = 3;

  const PhysR<T, D> contactPhysDeltaX(2.1);
  const PhysR<T, D> position(-3, -4, -5);
  PhysR<T, D>       max(-5);
  PhysR<T, D>       min(5);
  PhysR<T, D>       contactPoint(-1);
  Vector<T, D>      contactNormal(-1.5);
  unsigned          volumeCellCount(0);
  T                 surfaceCellCount(0);

  std::function<bool(PhysR<T, D>)> isInside = [](PhysR<T, D> x) {
    return true;
  };
  std::function<bool(PhysR<T, D>, PhysR<T, D>, bool&, bool&, Vector<T, D>,
                     Vector<T, D>)>
      isOnSurface = [](PhysR<T, D> x, PhysR<T, D> dx, bool& onSurfaceA,
                       bool& onSurfaceB, Vector<T, D> nA, Vector<T, D> nB) {
        return false;
      };
  std::function<Vector<T, D>(PhysR<T, D>, T)> calcNormal = [](PhysR<T, D> x,
                                                              T           m) {
    return Vector<T, D>(1.5);
  };
  std::function<void(const PhysR<T, D>&)> updateMinMax =
      [&min, &max](const PhysR<T, D>& pos) {
        particles::contact::updateMinMax(min, max, pos);
      };

  particles::contact::processCell(position, contactPhysDeltaX, contactPoint,
                                  contactNormal, volumeCellCount,
                                  surfaceCellCount, isInside, isOnSurface,
                                  calcNormal, calcNormal, updateMinMax);

  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_NEAR(min[iD], position[iD], 1e-16);
    EXPECT_NEAR(max[iD], position[iD], 1e-16);
    EXPECT_EQ(contactNormal[iD], -1.5);
    EXPECT_NEAR(contactPoint[iD], -1 + position[iD], 1e-16);
  }
  EXPECT_EQ(volumeCellCount, 1);
  EXPECT_EQ(surfaceCellCount, 0);
}

TEST(cellProcessing, testOutsideNotOnSurface)
{
  constexpr unsigned D = 3;

  const PhysR<T, D> contactPhysDeltaX(2.1);
  const PhysR<T, D> position(-3);
  PhysR<T, D>       max(-5);
  PhysR<T, D>       min(5);
  PhysR<T, D>       contactPoint(-1);
  Vector<T, D>      contactNormal(-1.5);
  unsigned          volumeCellCount(0);
  T                 surfaceCellCount(0);

  std::function<bool(PhysR<T, D>)> isInside = [](PhysR<T, D> x) {
    return false;
  };
  std::function<bool(PhysR<T, D>, PhysR<T, D>, bool&, bool&, Vector<T, D>,
                     Vector<T, D>)>
      isOnSurface = [](PhysR<T, D> x, PhysR<T, D> dx, bool& onSurfaceA,
                       bool& onSurfaceB, Vector<T, D> nA, Vector<T, D> nB) {
        return false;
      };
  std::function<Vector<T, D>(PhysR<T, D>, T)> calcNormal = [](PhysR<T, D> x,
                                                              T           m) {
    return Vector<T, D>(1.5);
  };
  std::function<void(const PhysR<T, D>&)> updateMinMax =
      [&min, &max](const PhysR<T, D>& pos) {
        particles::contact::updateMinMax(min, max, pos);
      };

  particles::contact::processCell(position, contactPhysDeltaX, contactPoint,
                                  contactNormal, volumeCellCount,
                                  surfaceCellCount, isInside, isOnSurface,
                                  calcNormal, calcNormal, updateMinMax);

  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_EQ(min[iD], 5);
    EXPECT_EQ(max[iD], -5);
    EXPECT_EQ(contactNormal[iD], -1.5);
    EXPECT_EQ(contactPoint[iD], -1);
  }
  EXPECT_EQ(volumeCellCount, 0);
  EXPECT_EQ(surfaceCellCount, 0);
}

TEST(cellProcessing, testOutsideOnSurfaceAzeroNorm)
{
  constexpr unsigned D = 3;

  const PhysR<T, D> contactPhysDeltaX(2.1);
  const PhysR<T, D> position(-3);
  PhysR<T, D>       max(-5);
  PhysR<T, D>       min(5);
  PhysR<T, D>       contactPoint(-1);
  Vector<T, D>      contactNormal(-1.5);
  unsigned          volumeCellCount(0);
  T                 surfaceCellCount(0);

  std::function<bool(PhysR<T, D>)> isInside = [](PhysR<T, D> x) {
    return false;
  };
  std::function<bool(PhysR<T, D>, PhysR<T, D>, bool&, bool&, Vector<T, D>,
                     Vector<T, D>)>
      isOnSurface = [](PhysR<T, D> x, PhysR<T, D> dx, bool& onSurfaceA,
                       bool& onSurfaceB, Vector<T, D> nA, Vector<T, D> nB) {
        onSurfaceA = true;
        return true;
      };
  std::function<Vector<T, D>(PhysR<T, D>, T)> calcNormal = [](PhysR<T, D> x,
                                                              T           m) {
    return Vector<T, D>(0.);
  };
  std::function<void(const PhysR<T, D>&)> updateMinMax =
      [&min, &max](const PhysR<T, D>& pos) {
        particles::contact::updateMinMax(min, max, pos);
      };

  particles::contact::processCell(position, contactPhysDeltaX, contactPoint,
                                  contactNormal, volumeCellCount,
                                  surfaceCellCount, isInside, isOnSurface,
                                  calcNormal, calcNormal, updateMinMax);

  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_EQ(min[iD], 5);
    EXPECT_EQ(max[iD], -5);
    EXPECT_EQ(contactNormal[iD], -1.5);
    EXPECT_EQ(contactPoint[iD], -1);
  }
  EXPECT_EQ(volumeCellCount, 0);
  EXPECT_EQ(surfaceCellCount, 0);
}

TEST(cellProcessing, testOutsideOnSurfaceAnonzeroNorm)
{
  constexpr unsigned D = 3;

  const PhysR<T, D> contactPhysDeltaX(2.1);
  const PhysR<T, D> position(-3);
  PhysR<T, D>       max(-5);
  PhysR<T, D>       min(5);
  PhysR<T, D>       contactPoint(-1);
  Vector<T, D>      contactNormal(-1.5);
  unsigned          volumeCellCount(0);
  T                 surfaceCellCount(0);

  const Vector<T, D>               normal(2.75);
  std::function<bool(PhysR<T, D>)> isInside = [](PhysR<T, D> x) {
    return false;
  };
  std::function<bool(PhysR<T, D>, PhysR<T, D>, bool&, bool&, Vector<T, D>,
                     Vector<T, D>)>
      isOnSurface = [](PhysR<T, D> x, PhysR<T, D> dx, bool& onSurfaceA,
                       bool& onSurfaceB, Vector<T, D> nA, Vector<T, D> nB) {
        onSurfaceA = true;
        return true;
      };
  std::function<Vector<T, D>(PhysR<T, D>, T)> calcNormal =
      [&normal](PhysR<T, D> x, T m) {
        return normal;
      };
  std::function<void(const PhysR<T, D>&)> updateMinMax =
      [&min, &max](const PhysR<T, D>& pos) {
        particles::contact::updateMinMax(min, max, pos);
      };

  particles::contact::processCell(position, contactPhysDeltaX, contactPoint,
                                  contactNormal, volumeCellCount,
                                  surfaceCellCount, isInside, isOnSurface,
                                  calcNormal, calcNormal, updateMinMax);
  const T surfaceApprox = particles::contact::evalApproxSurface(
      normalize(normal), contactPhysDeltaX);

  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_EQ(min[iD], 5);
    EXPECT_EQ(max[iD], -5);
    EXPECT_NEAR(contactNormal[iD], -1.5 + surfaceApprox * normalize(normal)[iD],
                2e-16);
    EXPECT_EQ(contactPoint[iD], -1);
  }
  EXPECT_EQ(volumeCellCount, 0);
  EXPECT_NEAR(surfaceCellCount, surfaceApprox, 1e-16);
}

TEST(cellProcessing, testOutsideOnSurfaceBzeroNorm)
{
  constexpr unsigned D = 3;

  const PhysR<T, D> contactPhysDeltaX(2.1);
  const PhysR<T, D> position(-3);
  PhysR<T, D>       max(-5);
  PhysR<T, D>       min(5);
  PhysR<T, D>       contactPoint(-1);
  Vector<T, D>      contactNormal(-1.5);
  unsigned          volumeCellCount(0);
  T                 surfaceCellCount(0);

  std::function<bool(PhysR<T, D>)> isInside = [](PhysR<T, D> x) {
    return false;
  };
  std::function<bool(PhysR<T, D>, PhysR<T, D>, bool&, bool&, Vector<T, D>,
                     Vector<T, D>)>
      isOnSurface = [](PhysR<T, D> x, PhysR<T, D> dx, bool& onSurfaceA,
                       bool& onSurfaceB, Vector<T, D> nA, Vector<T, D> nB) {
        onSurfaceB = true;
        return true;
      };
  std::function<Vector<T, D>(PhysR<T, D>, T)> calcNormal = [](PhysR<T, D> x,
                                                              T           m) {
    return Vector<T, D>(0.);
  };
  std::function<void(const PhysR<T, D>&)> updateMinMax =
      [&min, &max](const PhysR<T, D>& pos) {
        particles::contact::updateMinMax(min, max, pos);
      };

  particles::contact::processCell(position, contactPhysDeltaX, contactPoint,
                                  contactNormal, volumeCellCount,
                                  surfaceCellCount, isInside, isOnSurface,
                                  calcNormal, calcNormal, updateMinMax);

  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_EQ(min[iD], 5);
    EXPECT_EQ(max[iD], -5);
    EXPECT_EQ(contactNormal[iD], -1.5);
    EXPECT_EQ(contactPoint[iD], -1);
  }
  EXPECT_EQ(volumeCellCount, 0);
  EXPECT_EQ(surfaceCellCount, 0);
}

TEST(cellProcessing, testOutsideOnSurfaceBnonzeroNorm)
{
  constexpr unsigned D = 3;

  const PhysR<T, D> contactPhysDeltaX(2.1);
  const PhysR<T, D> position(-3);
  PhysR<T, D>       max(-5);
  PhysR<T, D>       min(5);
  PhysR<T, D>       contactPoint(-1);
  Vector<T, D>      contactNormal(-1.5);
  unsigned          volumeCellCount(0);
  T                 surfaceCellCount(0);

  const Vector<T, D>               normal(2.75);
  std::function<bool(PhysR<T, D>)> isInside = [](PhysR<T, D> x) {
    return false;
  };
  std::function<bool(PhysR<T, D>, PhysR<T, D>, bool&, bool&, Vector<T, D>,
                     Vector<T, D>)>
      isOnSurface = [](PhysR<T, D> x, PhysR<T, D> dx, bool& onSurfaceA,
                       bool& onSurfaceB, Vector<T, D> nA, Vector<T, D> nB) {
        onSurfaceB = true;
        return true;
      };
  std::function<Vector<T, D>(PhysR<T, D>, T)> calcNormal =
      [&normal](PhysR<T, D> x, T m) {
        return normal;
      };
  std::function<void(const PhysR<T, D>&)> updateMinMax =
      [&min, &max](const PhysR<T, D>& pos) {
        particles::contact::updateMinMax(min, max, pos);
      };

  particles::contact::processCell(position, contactPhysDeltaX, contactPoint,
                                  contactNormal, volumeCellCount,
                                  surfaceCellCount, isInside, isOnSurface,
                                  calcNormal, calcNormal, updateMinMax);
  const T surfaceApprox = particles::contact::evalApproxSurface(
      normalize(normal), contactPhysDeltaX);

  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_EQ(min[iD], 5);
    EXPECT_EQ(max[iD], -5);
    EXPECT_NEAR(contactNormal[iD], -1.5 - surfaceApprox * normalize(normal)[iD],
                1e-15);
    EXPECT_EQ(contactPoint[iD], -1);
  }
  EXPECT_EQ(volumeCellCount, 0);
  EXPECT_NEAR(surfaceCellCount, surfaceApprox, 1e-16);
}

TEST(cellProcessing, testOutsideOnSurfaceAOnSurfaceB)
{
  constexpr unsigned D = 3;

  const PhysR<T, D> contactPhysDeltaX(2.1);
  const PhysR<T, D> position(-3);
  PhysR<T, D>       max(-5);
  PhysR<T, D>       min(5);
  PhysR<T, D>       contactPoint(-1);
  Vector<T, D>      contactNormal(-1.5);
  unsigned          volumeCellCount(0);
  T                 surfaceCellCount(0);

  const Vector<T, D>               normal(2.75);
  std::function<bool(PhysR<T, D>)> isInside = [](PhysR<T, D> x) {
    return false;
  };
  std::function<bool(PhysR<T, D>, PhysR<T, D>, bool&, bool&, Vector<T, D>,
                     Vector<T, D>)>
      isOnSurface = [](PhysR<T, D> x, PhysR<T, D> dx, bool& onSurfaceA,
                       bool& onSurfaceB, Vector<T, D> nA, Vector<T, D> nB) {
        onSurfaceA = true;
        onSurfaceB = true;
        return true;
      };
  std::function<Vector<T, D>(PhysR<T, D>, T)> calcNormal =
      [&normal](PhysR<T, D> x, T m) {
        return normal;
      };
  std::function<void(const PhysR<T, D>&)> updateMinMax =
      [&min, &max](const PhysR<T, D>& pos) {
        particles::contact::updateMinMax(min, max, pos);
      };

  particles::contact::processCell(position, contactPhysDeltaX, contactPoint,
                                  contactNormal, volumeCellCount,
                                  surfaceCellCount, isInside, isOnSurface,
                                  calcNormal, calcNormal, updateMinMax);
  const T surfaceApprox = particles::contact::evalApproxSurface(
      normalize(normal), contactPhysDeltaX);

  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_EQ(min[iD], 5);
    EXPECT_EQ(max[iD], -5);
    EXPECT_EQ(contactNormal[iD], -1.5);
    EXPECT_EQ(contactPoint[iD], -1);
  }
  EXPECT_EQ(volumeCellCount, 0);
  EXPECT_NEAR(surfaceCellCount, 2 * surfaceApprox, 1e-16);
}
