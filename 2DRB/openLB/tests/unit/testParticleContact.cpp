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

typedef ResolvedParticleWithContact2D PARTICLE_2D;
typedef ResolvedParticleWithContact3D PARTICLE_3D;

using D2D =
    D2Q9<POROSITY, VELOCITY_NUMERATOR, VELOCITY_DENOMINATOR, CONTACT_DETECTION>;
using D3D = D3Q19<POROSITY, VELOCITY_NUMERATOR, VELOCITY_DENOMINATOR,
                  CONTACT_DETECTION>;

/// ensure olbInit is called for MPI
::testing::Environment* const foo_env =
    ::testing::AddGlobalTestEnvironment(new OlbInitTestEnvironment);

UnitConverterFromResolutionAndRelaxationTime<T, D2D> const converter2D(
    int {1}, // resolution: number of voxels per charPhysL
    (T)1, // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)1, // charPhysLength: reference length of simulation geometry
    (T)1, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)1, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)1  // physDensity: physical density in __kg / m^3__
);

UnitConverterFromResolutionAndRelaxationTime<T, D3D> const converter3D(
    int {1}, // resolution: number of voxels per charPhysL
    (T)1, // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)1, // charPhysLength: reference length of simulation geometry
    (T)1, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)1, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)1  // physDensity: physical density in __kg / m^3__
);

TEST(particleContactObject, testCopyConstructorAndAssignment)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                             PARTICLECONTACTTYPE;
  std::array<std::size_t, 2> particleIDs;
  particleIDs[0] = 79;
  particleIDs[1] = 43;
  PARTICLECONTACTTYPE contact(particleIDs);
  contact.setParticlePositions({PhysR<T, 3>(0.9), PhysR<T, 3>(1.3)});
  contact.setResponsibleRank(13);
  contact.setDampingFactor(2.2);
  contact.updateMinMax(PhysR<T, 3>(1.1));
  contact.increaseMinMax(Vector<T, 3>(0.5));
  contact.isNew(false);
  contact.setParticlePositionUpdated(true);

  PARTICLECONTACTTYPE constructed(contact);
  PARTICLECONTACTTYPE assigned;
  assigned = contact;

  EXPECT_EQ(contact.getIDs(), constructed.getIDs());
  EXPECT_EQ(contact.getResponsibleRank(), constructed.getResponsibleRank());
  EXPECT_EQ(contact.getParticlePositions(), constructed.getParticlePositions());
  EXPECT_EQ(contact.getDampingFactor(), constructed.getDampingFactor());
  EXPECT_EQ(contact.getMin(), constructed.getMin());
  EXPECT_EQ(contact.getMax(), constructed.getMax());
  EXPECT_EQ(contact.isNew(), constructed.isNew());
  EXPECT_EQ(contact.isParticlePositionUpdated(),
            constructed.isParticlePositionUpdated());

  EXPECT_EQ(contact.getIDs(), assigned.getIDs());
  EXPECT_EQ(contact.getResponsibleRank(), assigned.getResponsibleRank());
  EXPECT_EQ(contact.getParticlePositions(), assigned.getParticlePositions());
  EXPECT_EQ(contact.getDampingFactor(), assigned.getDampingFactor());
  EXPECT_EQ(contact.getMin(), assigned.getMin());
  EXPECT_EQ(contact.getMax(), assigned.getMax());
  EXPECT_EQ(contact.isNew(), assigned.isNew());
  EXPECT_EQ(contact.isParticlePositionUpdated(),
            assigned.isParticlePositionUpdated());
}

//Sphere tests

TEST(particleContactViaOverlapImproved3D, testTwoSpheresSmallOverlap)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 1>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T radius                   = 25;
  const int nCellX                   = 48. / converter3D.getPhysDeltaX();
  const T coefficientOfRestitution = 1;
  const T eps                      = 0.;
  const T expectedIndentationDepth = 2;
  const T density                  = 1;
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                               nCellX, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create sphere 1
  Vector<T, 3> pos1(25, 0, 25);
  particles::creators::addResolvedSphere3D(particleSystem, pos1, radius, eps,
                                           density);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);

  //Create sphere 2
  Vector<T, 3> pos2(25, 2 * radius - expectedIndentationDepth, 25);
  particles::creators::addResolvedSphere3D(particleSystem, pos2, radius, eps,
                                           density);
  auto particleB = particleSystem.get(1);
  particleB.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  contactProperties.set(particleMaterial, particleMaterial, effectiveE,
                        coefficientOfRestitution, 0, 0);

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D);

  // Force calculation
  particles::contact::particle_particle<
      T, PARTICLE_3D, PARTICLECONTACTTYPE,
      WALLCONTACTTYPE>::apply(particleSystem, contactContainer,
                              contactProperties, superGeometry, dataMap, 16);

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T alpha               = 3. / 2.;
  const T effectiveRadius     = 1. / ((1. / radius) + (1. / radius));
  const T expectedNormalForce = 4. / 3. * effectiveE *
                                std::pow(effectiveRadius, 2 - alpha) *
                                std::pow(expectedIndentationDepth, alpha);

  for (std::size_t iP = 0; iP < particleSystem.size(); ++iP) {
    auto                      particle = particleSystem.get(iP);
    Vector<T, PARTICLE_3D::d> force =
        particle.template getField<FORCING, FORCE>();
    const short sign = (iP == 1 ? 1 : -1);
    EXPECT_NEAR(sign * force[0], 0, expectedNormalForce * 1e-4);
    // deviation of up to 2% is allowed
    EXPECT_NEAR(sign * force[1], expectedNormalForce,
                expectedNormalForce * 2e-2);
    EXPECT_NEAR(sign * force[2], 0, expectedNormalForce * 1e-4);
  }

  EXPECT_EQ(contactContainer.particleContacts.size(), 1);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[0], 0);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[1], 1);
}

TEST(particleContactViaOverlapImproved3D, testTwoEnlargedSpheresSmallOverlap)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 1>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T radius                        = 25;
  const int nCellX                        = 48. / converter3D.getPhysDeltaX();
  const T particleEnlargementForContact = converter3D.getPhysDeltaX();
  const T coefficientOfRestitution      = 1;
  const T eps                           = 0.;
  const T expectedIndentationDepth      = 1;
  const T density                       = 1;
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     {nCellX, nCellX, nCellX}, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create sphere 1
  Vector<T, 3> pos1(25, 0, 25);
  particles::creators::addResolvedSphere3D(particleSystem, pos1, radius, eps,
                                           density);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  particleA.setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(
      particleEnlargementForContact);

  //Create sphere 2
  Vector<T, 3> pos2(25, 2 * radius - expectedIndentationDepth, 25);
  particles::creators::addResolvedSphere3D(particleSystem, pos2, radius, eps,
                                           density);
  auto particleB = particleSystem.get(1);
  particleB.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  particleB.setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(
      particleEnlargementForContact);
  contactProperties.set(particleMaterial, particleMaterial, effectiveE,
                        coefficientOfRestitution, 0, 0);
  auto solidBoundaries = std::vector<SolidBoundary<T, D3D::d>>();

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D,
      solidBoundaries, []() {
        return Vector<bool, 3>(false, false, false);
      });

  // Force calculation
  particles::contact::particle_particle<
      T, PARTICLE_3D, PARTICLECONTACTTYPE,
      WALLCONTACTTYPE>::apply(particleSystem, contactContainer,
                              contactProperties, superGeometry, dataMap, 16,
                              4. / (3. * util::sqrt(M_PI)));

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T alpha = 3. / 2.;
  const T effectiveRadius =
      1. / ((1. / (radius + particleEnlargementForContact)) +
            (1. / (radius + particleEnlargementForContact)));
  const T expectedNormalForce =
      4. / 3. * effectiveE * std::pow(effectiveRadius, 2 - alpha) *
      std::pow(expectedIndentationDepth + 2 * particleEnlargementForContact,
               alpha);

  for (std::size_t iP = 0; iP < particleSystem.size(); ++iP) {
    auto                      particle = particleSystem.get(iP);
    Vector<T, PARTICLE_3D::d> force =
        particle.template getField<FORCING, FORCE>();
    const short sign = (iP == 1 ? 1 : -1);
    EXPECT_NEAR(sign * force[0], 0, expectedNormalForce * 1e-3);
    // deviation of up to 2% is allowed
    EXPECT_NEAR(sign * force[1], expectedNormalForce,
                expectedNormalForce * 2e-2);
    EXPECT_NEAR(sign * force[2], 0, expectedNormalForce * 1e-3);
  }

  EXPECT_EQ(contactContainer.particleContacts.size(), 1);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[0], 0);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[1], 1);
}
TEST(particleContactViaOverlapImproved3D,
     testTwoSpheresSmallOverlapPeriodicBoundary)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 1>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T radius                   = 25;
  const int nCellX                   = 48. / converter3D.getPhysDeltaX();
  const T coefficientOfRestitution = 1;
  const T eps                      = 0.;
  const T expectedIndentationDepth = 2;
  const T density                  = 1;
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     {nCellX, nCellX, nCellX}, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create sphere 1
  Vector<T, 3> pos1(0, 0, 0);
  particles::creators::addResolvedSphere3D(particleSystem, pos1, radius, eps,
                                           density);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);

  //Create sphere 2
  Vector<T, 3> pos2(0, 2 * radius - expectedIndentationDepth, 0);
  particles::creators::addResolvedSphere3D(particleSystem, pos2, radius, eps,
                                           density);
  auto particleB = particleSystem.get(1);
  particleB.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  contactProperties.set(particleMaterial, particleMaterial, effectiveE,
                        coefficientOfRestitution, 0, 0);
  auto solidBoundaries = std::vector<SolidBoundary<T, D3D::d>>();

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D,
      solidBoundaries, []() {
        return Vector<bool, 3>(true, false, true);
      });

  // Force calculation
  particles::contact::particle_particle<
      T, PARTICLE_3D, PARTICLECONTACTTYPE,
      WALLCONTACTTYPE>::apply(particleSystem, contactContainer,
                              contactProperties, superGeometry, dataMap, 16);

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T alpha               = 3. / 2.;
  const T effectiveRadius     = 1. / ((1. / radius) + (1. / radius));
  const T expectedNormalForce = 4. / 3. * effectiveE *
                                std::pow(effectiveRadius, 2 - alpha) *
                                std::pow(expectedIndentationDepth, alpha);

  for (std::size_t iP = 0; iP < particleSystem.size(); ++iP) {
    auto                      particle = particleSystem.get(iP);
    Vector<T, PARTICLE_3D::d> force =
        particle.template getField<FORCING, FORCE>();
    const short sign = (iP == 1 ? 1 : -1);
    EXPECT_NEAR(sign * force[0], 0, 1e-16);
    // deviation of up to 2% is allowed
    EXPECT_NEAR(sign * force[1], expectedNormalForce,
                expectedNormalForce * 2e-2);
    EXPECT_NEAR(sign * force[2], 0, 1e-16);
  }

  EXPECT_EQ(contactContainer.particleContacts.size(), 1);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[0], 0);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[1], 1);
}

TEST(particleContactViaOverlapImproved3D, testWallSphereSmallOverlap)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 1>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T radius                   = 25;
  const int nCellX                   = 48. / converter3D.getPhysDeltaX();
  const T coefficientOfRestitution = 1;
  const T eps                      = 0.;
  const T expectedIndentationDepth = 2;
  const T density                  = 1;
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     {nCellX, nCellX, nCellX}, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create cuboid (wall)
  Vector<T, 3> pos1(25, 0, 25);
  Vector<T, 3> extend(50, 50, 50);
  particles::creators::addResolvedCuboid3D(particleSystem, pos1, extend, eps,
                                           density);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);

  //Create sphere
  Vector<T, 3> pos2(25, extend[1] / 2 + radius - expectedIndentationDepth, 25);
  particles::creators::addResolvedSphere3D(particleSystem, pos2, radius, eps,
                                           density);
  auto particleB = particleSystem.get(1);
  particleB.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  contactProperties.set(particleMaterial, particleMaterial, effectiveE,
                        coefficientOfRestitution, 0, 0);

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D);

  // Force calculation
  particles::contact::particle_particle<
      T, PARTICLE_3D, PARTICLECONTACTTYPE,
      WALLCONTACTTYPE>::apply(particleSystem, contactContainer,
                              contactProperties, superGeometry, dataMap, 16);

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T alpha               = 3. / 2.;
  const T expectedNormalForce = 4. / 3. * effectiveE *
                                std::pow(radius, 2 - alpha) *
                                std::pow(expectedIndentationDepth, alpha);

  for (std::size_t iP = 0; iP < particleSystem.size(); ++iP) {
    auto                      particle = particleSystem.get(iP);
    Vector<T, PARTICLE_3D::d> force =
        particle.template getField<FORCING, FORCE>();
    const short sign = (iP == 1 ? 1 : -1);
    EXPECT_NEAR(sign * force[0], 0, expectedNormalForce * 5e-3);
    // deviation of up to 2% is allowed
    EXPECT_NEAR(sign * force[1], expectedNormalForce,
                expectedNormalForce * 2e-2);
    EXPECT_NEAR(sign * force[2], 0, expectedNormalForce * 5e-3);
  }

  EXPECT_EQ(contactContainer.particleContacts.size(), 1);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[0], 0);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[1], 1);
}

TEST(particleContactViaOverlapImproved3D, testWallEnlargedSphereSmallOverlap)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 1>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T radius                        = 25;
  const int nCellX                        = 48. / converter3D.getPhysDeltaX();
  const T particleEnlargementForContact = 2 * converter3D.getPhysDeltaX();
  const T coefficientOfRestitution      = 1;
  const T eps                           = 0.;
  const T expectedIndentationDepth      = 1;
  const T density                       = 1;
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCellX, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create cuboid (wall)
  Vector<T, 3> pos1(25., 0., 25.);
  Vector<T, 3> extend(50, 50, 50);
  particles::creators::addResolvedCuboid3D(particleSystem, pos1, extend, eps,
                                           density);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  particleA.setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(
      particleEnlargementForContact);

  //Create sphere
  Vector<T, 3> pos2(25, 0.5 * extend[1] + radius - expectedIndentationDepth,
                    25);
  particles::creators::addResolvedSphere3D(particleSystem, pos2, radius, eps,
                                           density);
  auto particleB = particleSystem.get(1);
  particleB.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  particleB.setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(
      particleEnlargementForContact);
  contactProperties.set(particleMaterial, particleMaterial, effectiveE,
                        coefficientOfRestitution, 0, 0);
  auto solidBoundaries = std::vector<SolidBoundary<T, D3D::d>>();

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D,
      solidBoundaries, []() {
        return Vector<bool, 3>(false, false, false);
      });

  // Force calculation
  particles::contact::particle_particle<
      T, PARTICLE_3D, PARTICLECONTACTTYPE,
      WALLCONTACTTYPE>::apply(particleSystem, contactContainer,
                              contactProperties, superGeometry, dataMap, 16,
                              4. / (3 * util::sqrt(M_PI)));

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T alpha = 3. / 2.;
  const T expectedNormalForce =
      4. / 3. * effectiveE *
      std::pow(radius + particleEnlargementForContact, 2 - alpha) *
      std::pow(expectedIndentationDepth + 2 * particleEnlargementForContact,
               alpha);

  for (std::size_t iP = 0; iP < particleSystem.size(); ++iP) {
    auto                      particle = particleSystem.get(iP);
    Vector<T, PARTICLE_3D::d> force =
        particle.template getField<FORCING, FORCE>();
    const short sign = (iP == 1 ? 1 : -1);
    EXPECT_NEAR(sign * force[0], 0, expectedNormalForce * 5e-3);
    // deviation of up to 2% is allowed
    EXPECT_NEAR(sign * force[1], expectedNormalForce,
                expectedNormalForce * 2e-2);
    EXPECT_NEAR(sign * force[2], 0, expectedNormalForce * 5e-3);
  }

  EXPECT_EQ(contactContainer.particleContacts.size(), 1);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[0], 0);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[1], 1);
}

//different expectedIndentationDepth

TEST(particleContactViaOverlapImproved3D, testTwoSpheresSmallOverlap2)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 1>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T radius                   = 25;
  const int nCellX                   = 48. / converter3D.getPhysDeltaX();
  const T coefficientOfRestitution = 1;
  const T eps                      = 0.;
  const T expectedIndentationDepth = 1;
  const T density                  = 1;
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCellX, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create sphere 1
  Vector<T, 3> pos1(25, 0, 25);
  particles::creators::addResolvedSphere3D(particleSystem, pos1, radius, eps,
                                           density);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);

  //Create sphere 2
  Vector<T, 3> pos2(25, 2 * radius - expectedIndentationDepth, 25);
  particles::creators::addResolvedSphere3D(particleSystem, pos2, radius, eps,
                                           density);
  auto particleB = particleSystem.get(1);
  particleB.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  contactProperties.set(particleMaterial, particleMaterial, effectiveE,
                        coefficientOfRestitution, 0, 0);

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D);

  // Force calculation
  particles::contact::particle_particle<
      T, PARTICLE_3D, PARTICLECONTACTTYPE,
      WALLCONTACTTYPE>::apply(particleSystem, contactContainer,
                              contactProperties, superGeometry, dataMap, 16);

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T alpha               = 3. / 2.;
  const T effectiveRadius     = 1. / ((1. / radius) + (1. / radius));
  const T expectedNormalForce = 4. / 3. * effectiveE *
                                std::pow(effectiveRadius, 2 - alpha) *
                                std::pow(expectedIndentationDepth, alpha);
  // evaluates to 3.1426968052735447

  for (std::size_t iP = 0; iP < particleSystem.size(); ++iP) {
    auto                      particle = particleSystem.get(iP);
    Vector<T, PARTICLE_3D::d> force =
        particle.template getField<FORCING, FORCE>();
    const short sign = (iP == 1 ? 1 : -1);
    EXPECT_NEAR(sign * force[0], 0, expectedNormalForce * 1e-4);
    // deviation of up to 2% is allowed
    EXPECT_NEAR(sign * force[1], expectedNormalForce,
                expectedNormalForce * 2e-2);
    EXPECT_NEAR(sign * force[2], 0, expectedNormalForce * 1e-4);
  }

  EXPECT_EQ(contactContainer.particleContacts.size(), 1);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[0], 0);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[1], 1);
}

TEST(particleContactViaOverlapImproved3D, testWallSphereSmallOverlap2)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 1>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T radius                   = 25;
  const int nCellX                   = 48. / converter3D.getPhysDeltaX();
  const T coefficientOfRestitution = 1;
  const T eps                      = 0.;
  const T expectedIndentationDepth = 1;
  const T density                  = 1;
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCellX, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create cuboid (wall)
  Vector<T, 3> pos1(25, 0, 25);
  Vector<T, 3> extend(50, 50, 50);
  particles::creators::addResolvedCuboid3D(particleSystem, pos1, extend, eps,
                                           density);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);

  //Create sphere
  Vector<T, 3> pos2(25, extend[1] / 2 + radius - expectedIndentationDepth, 25);
  particles::creators::addResolvedSphere3D(particleSystem, pos2, radius, eps,
                                           density);
  auto particleB = particleSystem.get(1);
  particleB.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  contactProperties.set(particleMaterial, particleMaterial, effectiveE,
                        coefficientOfRestitution, 0, 0);

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D);

  // Force calculation
  particles::contact::particle_particle<
      T, PARTICLE_3D, PARTICLECONTACTTYPE,
      WALLCONTACTTYPE>::apply(particleSystem, contactContainer,
                              contactProperties, superGeometry, dataMap, 15);

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T alpha               = 3. / 2.;
  const T expectedNormalForce = 4. / 3. * effectiveE *
                                std::pow(radius, 2 - alpha) *
                                std::pow(expectedIndentationDepth, alpha);
  // evaluates to 4.4444444444444446

  for (std::size_t iP = 0; iP < particleSystem.size(); ++iP) {
    auto                      particle = particleSystem.get(iP);
    Vector<T, PARTICLE_3D::d> force =
        particle.template getField<FORCING, FORCE>();
    const short sign = (iP == 1 ? 1 : -1);
    EXPECT_NEAR(sign * force[0], 0, expectedNormalForce * 3e-3);
    // deviation of up to 4% is allowed
    EXPECT_NEAR(sign * force[1], expectedNormalForce,
                expectedNormalForce * 4e-2);
    EXPECT_NEAR(sign * force[2], 0, expectedNormalForce * 3e-3);
  }

  EXPECT_EQ(contactContainer.particleContacts.size(), 1);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[0], 0);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[1], 1);
}

//different radius

TEST(particleContactViaOverlapImproved3D, testTwoSpheresSmallOverlap3)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 1>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T radius1                  = 30;
  const T radius2                  = 20;
  const int nCellX                   = 48. / converter3D.getPhysDeltaX();
  const T coefficientOfRestitution = 1;
  const T eps                      = 0.;
  const T expectedIndentationDepth = 2;
  const T density                  = 1;
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCellX, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create sphere 1
  Vector<T, 3> pos1(25, 0, 25);
  particles::creators::addResolvedSphere3D(particleSystem, pos1, radius1, eps,
                                           density);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);

  //Create sphere 2
  Vector<T, 3> pos2(25, radius1 + radius2 - expectedIndentationDepth, 25);
  particles::creators::addResolvedSphere3D(particleSystem, pos2, radius2, eps,
                                           density);
  auto particleB = particleSystem.get(1);
  particleB.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  contactProperties.set(particleMaterial, particleMaterial, effectiveE,
                        coefficientOfRestitution, 0, 0);

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D);

  // Force calculation
  particles::contact::particle_particle<
      T, PARTICLE_3D, PARTICLECONTACTTYPE,
      WALLCONTACTTYPE>::apply(particleSystem, contactContainer,
                              contactProperties, superGeometry, dataMap, 16);

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T alpha               = 3. / 2.;
  const T effectiveRadius     = 1. / ((1. / radius1) + (1. / radius2));
  const T expectedNormalForce = 4. / 3. * effectiveE *
                                std::pow(effectiveRadius, 2 - alpha) *
                                std::pow(expectedIndentationDepth, alpha);
  // evaluates to 8.7092968632290777

  for (std::size_t iP = 0; iP < particleSystem.size(); ++iP) {
    auto                      particle = particleSystem.get(iP);
    Vector<T, PARTICLE_3D::d> force =
        particle.template getField<FORCING, FORCE>();
    const short sign = (iP == 1 ? 1 : -1);
    EXPECT_NEAR(sign * force[0], 0, expectedNormalForce * 5e-3);
    // deviation of up to 2% is allowed
    EXPECT_NEAR(sign * force[1], expectedNormalForce,
                expectedNormalForce * 2e-2);
    EXPECT_NEAR(sign * force[2], 0, expectedNormalForce * 5e-3);
  }

  EXPECT_EQ(contactContainer.particleContacts.size(), 1);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[0], 0);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[1], 1);
}

TEST(particleContactViaOverlapImproved3D, testWallSphereSmallOverlap3)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 1>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T radius                   = 30;
  const int nCellX                   = 48. / converter3D.getPhysDeltaX();
  const T coefficientOfRestitution = 1;
  const T eps                      = 0.;
  const T expectedIndentationDepth = 2;
  const T density                  = 1;
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCellX, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create cuboid (wall)
  Vector<T, 3> pos1(25, 0, 25);
  Vector<T, 3> extend(50, 50, 50);
  particles::creators::addResolvedCuboid3D(particleSystem, pos1, extend, eps,
                                           density);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);

  //Create sphere
  Vector<T, 3> pos2(25, extend[1] / 2 + radius - expectedIndentationDepth, 25);
  particles::creators::addResolvedSphere3D(particleSystem, pos2, radius, eps,
                                           density);
  auto particleB = particleSystem.get(1);
  particleB.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  contactProperties.set(particleMaterial, particleMaterial, effectiveE,
                        coefficientOfRestitution, 0, 0);

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D);

  // Force calculation
  particles::contact::particle_particle<
      T, PARTICLE_3D, PARTICLECONTACTTYPE,
      WALLCONTACTTYPE>::apply(particleSystem, contactContainer,
                              contactProperties, superGeometry, dataMap, 16);

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T alpha               = 3. / 2.;
  const T expectedNormalForce = 4. / 3. * effectiveE *
                                std::pow(radius, 2 - alpha) *
                                std::pow(expectedIndentationDepth, alpha);
  // evaluates to 13.770607453181926

  for (std::size_t iP = 0; iP < particleSystem.size(); ++iP) {
    auto                      particle = particleSystem.get(iP);
    Vector<T, PARTICLE_3D::d> force =
        particle.template getField<FORCING, FORCE>();
    const short sign = (iP == 1 ? 1 : -1);
    EXPECT_NEAR(sign * force[0], 0, expectedNormalForce * 5e-3);
    // deviation of up to 2% is allowed
    EXPECT_NEAR(sign * force[1], expectedNormalForce,
                expectedNormalForce * 2e-2);
    EXPECT_NEAR(sign * force[2], 0, expectedNormalForce * 5e-3);
  }

  EXPECT_EQ(contactContainer.particleContacts.size(), 1);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[0], 0);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[1], 1);
}

//Cylinder tests

TEST(particleContactViaOverlapImproved3D, testWallCylinderSmallOverlap)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 1>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T      radius                   = 25;
  const int    nCellX                   = 48. / converter3D.getPhysDeltaX();
  const T      coefficientOfRestitution = 1;
  const T      eps                      = 0.;
  const T      expectedIndentationDepth = 2;
  const T      density                  = 1;
  const T      height                   = 50;
  Vector<T, 3> angle(T(0));
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCellX, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create cuboid (wall)
  Vector<T, 3> pos1(25, 0, 25);
  Vector<T, 3> extend(50, 50, 50);
  particles::creators::addResolvedCuboid3D(particleSystem, pos1, extend, eps,
                                           density);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);

  //Create cylinder
  Vector<T, 3> pos2(25, extend[1] - expectedIndentationDepth, 25);
  Vector<T, 3> normal2(0, 1, 0);
  particles::creators::addResolvedCylinder3D(
      particleSystem, pos2, normal2, height, radius, eps, density, angle);
  auto particleB = particleSystem.get(1);
  particleB.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  contactProperties.set(particleMaterial, particleMaterial, effectiveE,
                        coefficientOfRestitution, 0, 0);

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D);

  // Force calculation
  particles::contact::particle_particle<
      T, PARTICLE_3D, PARTICLECONTACTTYPE,
      WALLCONTACTTYPE>::apply(particleSystem, contactContainer,
                              contactProperties, superGeometry, dataMap, 16,
                              2. / util::sqrt(M_PI));

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T expectedNormalForce =
      2 * radius * effectiveE * expectedIndentationDepth;
  // evaluates to 66.666666666666657

  for (std::size_t iP = 0; iP < particleSystem.size(); ++iP) {
    auto                      particle = particleSystem.get(iP);
    Vector<T, PARTICLE_3D::d> force =
        particle.template getField<FORCING, FORCE>();
    const short sign = (iP == 1 ? 1 : -1);
    EXPECT_NEAR(sign * force[0], 0, expectedNormalForce * 2e-3);
    // deviation of up to 2% is allowed
    EXPECT_NEAR(sign * force[1], expectedNormalForce,
                expectedNormalForce * 2e-2);
    EXPECT_NEAR(sign * force[2], 0, expectedNormalForce * 2e-3);
  }

  EXPECT_EQ(contactContainer.particleContacts.size(), 1);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[0], 0);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[1], 1);
}

TEST(particleContactViaOverlapImproved3D,
     testWallCylinderSmallOverlapPeriodicBoundary)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 1>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T      radius                   = 25;
  const int    nCellX                   = 48. / converter3D.getPhysDeltaX();
  const T      coefficientOfRestitution = 1;
  const T      eps                      = 0.;
  const T      expectedIndentationDepth = 2;
  const T      density                  = 1;
  const T      height                   = 50;
  Vector<T, 3> angle(T(0));
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCellX, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create cuboid (wall)
  Vector<T, 3> pos1(0, 0, 0);
  Vector<T, 3> extend(50, 50, 50);
  particles::creators::addResolvedCuboid3D(particleSystem, pos1, extend, eps,
                                           density);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);

  //Create cylinder
  Vector<T, 3> pos2(0, extend[1] - expectedIndentationDepth, 0);
  Vector<T, 3> normal2(0, 1, 0);
  particles::creators::addResolvedCylinder3D(
      particleSystem, pos2, normal2, height, radius, eps, density, angle);
  auto particleB = particleSystem.get(1);
  particleB.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  contactProperties.set(particleMaterial, particleMaterial, effectiveE,
                        coefficientOfRestitution, 0, 0);
  auto solidBoundaries = std::vector<SolidBoundary<T, D3D::d>>();

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D,
      solidBoundaries, []() {
        return Vector<bool, 3>(true, false, true);
      });

  // Force calculation
  particles::contact::particle_particle<
      T, PARTICLE_3D, PARTICLECONTACTTYPE,
      WALLCONTACTTYPE>::apply(particleSystem, contactContainer,
                              contactProperties, superGeometry, dataMap, 16,
                              2. / util::sqrt(M_PI));

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T expectedNormalForce =
      2 * radius * effectiveE * expectedIndentationDepth;
  // evaluates to 66.666666666666657

  // TODO: fix periodic boundary, which causes different forces here than in the test before
  for (std::size_t iP = 0; iP < particleSystem.size(); ++iP) {
    auto                      particle = particleSystem.get(iP);
    Vector<T, PARTICLE_3D::d> force =
        particle.template getField<FORCING, FORCE>();
    const short sign = (iP == 1 ? 1 : -1);
    EXPECT_NEAR(sign * force[0], 0, expectedNormalForce * 2e-3);
    // deviation of up to 2% is allowed
    EXPECT_NEAR(sign * force[1], expectedNormalForce,
                expectedNormalForce * 2.5e-2);
    EXPECT_NEAR(sign * force[2], 0, expectedNormalForce * 2e-3);
  }

  EXPECT_EQ(contactContainer.particleContacts.size(), 1);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[0], 0);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[1], 1);
}

//different expectedIndentationDepth

TEST(particleContactViaOverlapImproved3D, testWallCylinderSmallOverlap2)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 1>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T      radius                   = 25;
  const int    nCellX                   = 48. / converter3D.getPhysDeltaX();
  const T      coefficientOfRestitution = 1;
  const T      eps                      = 0.;
  const T      expectedIndentationDepth = 1;
  const T      density                  = 1;
  const T      height                   = 50;
  Vector<T, 3> angle(T(0));
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCellX, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create cuboid (wall)
  Vector<T, 3> pos1(25, 0, 25);
  Vector<T, 3> extend(50, 50, 50);
  particles::creators::addResolvedCuboid3D(particleSystem, pos1, extend, eps,
                                           density);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);

  //Create cylinder
  Vector<T, 3> pos2(25, extend[1] - expectedIndentationDepth, 25);
  Vector<T, 3> normal(0, 1, 0);
  particles::creators::addResolvedCylinder3D(
      particleSystem, pos2, normal, height, radius, eps, density, angle);
  auto particleB = particleSystem.get(1);
  particleB.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  contactProperties.set(particleMaterial, particleMaterial, effectiveE,
                        coefficientOfRestitution, 0, 0);

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D);

  // Force calculation
  particles::contact::particle_particle<
      T, PARTICLE_3D, PARTICLECONTACTTYPE,
      WALLCONTACTTYPE>::apply(particleSystem, contactContainer,
                              contactProperties, superGeometry, dataMap, 16,
                              2. / util::sqrt(M_PI));

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T expectedNormalForce =
      2 * radius * effectiveE * expectedIndentationDepth;

  for (std::size_t iP = 0; iP < particleSystem.size(); ++iP) {
    auto                      particle = particleSystem.get(iP);
    Vector<T, PARTICLE_3D::d> force =
        particle.template getField<FORCING, FORCE>();
    const short sign = (iP == 1 ? 1 : -1);
    EXPECT_NEAR(sign * force[0], 0, expectedNormalForce * 1e-3);
    // deviation of up to 4% is allowed
    EXPECT_NEAR(sign * force[1], expectedNormalForce,
                expectedNormalForce * 4e-2);
    EXPECT_NEAR(sign * force[2], 0, expectedNormalForce * 1e-3);
  }

  EXPECT_EQ(contactContainer.particleContacts.size(), 1);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[0], 0);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[1], 1);
}

//different radius

TEST(particleContactViaOverlapImproved3D, testWallCylinderSmallOverlap3)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 1>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T      radius                   = 10;
  const int    nCellX                   = 48. / converter3D.getPhysDeltaX();
  const T      coefficientOfRestitution = 1;
  const T      eps                      = 0.;
  const T      expectedIndentationDepth = 2;
  const T      density                  = 1;
  const T      height                   = 50;
  Vector<T, 3> angle(T(0));
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCellX, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create cuboid (wall)
  Vector<T, 3> pos1(25, 0, 25);
  Vector<T, 3> extend(50, 50, 50);
  particles::creators::addResolvedCuboid3D(particleSystem, pos1, extend, eps,
                                           density);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);

  //Create cylinder
  Vector<T, 3> pos2(25, extend[1] - expectedIndentationDepth, 25);
  Vector<T, 3> normal2(0, 1, 0);
  particles::creators::addResolvedCylinder3D(
      particleSystem, pos2, normal2, height, radius, eps, density, angle);
  auto particleB = particleSystem.get(1);
  particleB.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  contactProperties.set(particleMaterial, particleMaterial, effectiveE,
                        coefficientOfRestitution, 0, 0);

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D);

  // Force calculation
  particles::contact::particle_particle<
      T, PARTICLE_3D, PARTICLECONTACTTYPE,
      WALLCONTACTTYPE>::apply(particleSystem, contactContainer,
                              contactProperties, superGeometry, dataMap, 16,
                              2. / util::sqrt(M_PI));

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T expectedNormalForce =
      2 * radius * effectiveE * expectedIndentationDepth;
  // evaluates to 26.666666666666664

  for (std::size_t iP = 0; iP < particleSystem.size(); ++iP) {
    auto                      particle = particleSystem.get(iP);
    Vector<T, PARTICLE_3D::d> force =
        particle.template getField<FORCING, FORCE>();
    const short sign = (iP == 1 ? 1 : -1);
    EXPECT_NEAR(sign * force[0], 0, expectedNormalForce * 1e-8);
    // deviation of up to 2% is allowed
    EXPECT_NEAR(sign * force[1], expectedNormalForce,
                expectedNormalForce * 2e-2);
    EXPECT_NEAR(sign * force[2], 0, expectedNormalForce * 1e-8);
  }

  EXPECT_EQ(contactContainer.particleContacts.size(), 1);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[0], 0);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[1], 1);
}

//different height

TEST(particleContactViaOverlapImproved3D, testWallCylinderSmallOverlap4)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 1>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T      radius                   = 25;
  const int    nCellX                   = 48. / converter3D.getPhysDeltaX();
  const T      coefficientOfRestitution = 1;
  const T      eps                      = 0.;
  const T      expectedIndentationDepth = 2;
  const T      density                  = 1;
  const T      height                   = 40;
  Vector<T, 3> angle(T(0));
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCellX, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create cuboid (wall)
  Vector<T, 3> pos1(25, 0, 25);
  Vector<T, 3> extend(50, 50, 50);
  particles::creators::addResolvedCuboid3D(particleSystem, pos1, extend, eps,
                                           density);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);

  //Create cylinder
  Vector<T, 3> pos2(25, extend[1] - 5 - expectedIndentationDepth, 25);
  Vector<T, 3> normal2(0, 1, 0);
  particles::creators::addResolvedCylinder3D(
      particleSystem, pos2, normal2, height, radius, eps, density, angle);
  auto particleB = particleSystem.get(1);
  particleB.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  contactProperties.set(particleMaterial, particleMaterial, effectiveE,
                        coefficientOfRestitution, 0, 0);

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D);

  // Force calculation
  particles::contact::particle_particle<
      T, PARTICLE_3D, PARTICLECONTACTTYPE,
      WALLCONTACTTYPE>::apply(particleSystem, contactContainer,
                              contactProperties, superGeometry, dataMap, 16,
                              2. / util::sqrt(M_PI));

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T expectedNormalForce =
      2 * radius * effectiveE * expectedIndentationDepth;

  for (std::size_t iP = 0; iP < particleSystem.size(); ++iP) {
    auto                      particle = particleSystem.get(iP);
    Vector<T, PARTICLE_3D::d> force =
        particle.template getField<FORCING, FORCE>();
    const short sign = (iP == 1 ? 1 : -1);
    EXPECT_NEAR(sign * force[0], 0, expectedNormalForce * 2e-3);
    // deviation of up to 2% is allowed
    EXPECT_NEAR(sign * force[1], expectedNormalForce,
                expectedNormalForce * 2e-2);
    EXPECT_NEAR(sign * force[2], 0, expectedNormalForce * 2e-3);
  }

  EXPECT_EQ(contactContainer.particleContacts.size(), 1);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[0], 0);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[1], 1);
}

//Cone tests

TEST(particleContactViaOverlapImproved3D, testWallConeSmallOverlap)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 1>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T      radius1                  = 25;
  const T      radius2                  = 0;
  const int    nCellX                   = 48. / converter3D.getPhysDeltaX();
  const T      coefficientOfRestitution = 1;
  const T      eps                      = 0.;
  const T      expectedIndentationDepth = 2;
  const T      density                  = 1;
  const T      height                   = 50;
  Vector<T, 3> angle(T(0));
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCellX, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create cuboid (wall)
  Vector<T, 3> pos1(25, 0, 25);
  Vector<T, 3> extend(50, 50, 50);
  particles::creators::addResolvedCuboid3D(particleSystem, pos1, extend, eps,
                                           density);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);

  //Create cone
  Vector<T, 3> center1(25, extend[1] / 2 + height - expectedIndentationDepth,
                       25);
  Vector<T, 3> center2(25, extend[1] / 2 - expectedIndentationDepth, 25);
  particles::creators::addResolvedCone3D(particleSystem, center1, center2,
                                         radius1, radius2, eps, density, angle);
  auto particleB = particleSystem.get(1);
  particleB.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  contactProperties.set(particleMaterial, particleMaterial, effectiveE,
                        coefficientOfRestitution, 0, 0);

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D);

  // Force calculation
  particles::contact::particle_particle<
      T, PARTICLE_3D, PARTICLECONTACTTYPE,
      WALLCONTACTTYPE>::apply(particleSystem, contactContainer,
                              contactProperties, superGeometry, dataMap, 16,
                              0.65);

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T inclinationAngle   = 0.5 * M_PI - util::atan(radius1 / height);
  const T contactRegionDepth = expectedIndentationDepth * 2 / M_PI;
  const T contactRadius = contactRegionDepth / (util::tan(inclinationAngle));
  const T expectedNormalForce = 0.5 * M_PI * util::pow(contactRadius, 2) *
                                effectiveE * util::tan(inclinationAngle);
  // evaluates to 3.3953054526271012

  for (std::size_t iP = 0; iP < particleSystem.size(); ++iP) {
    auto                      particle = particleSystem.get(iP);
    Vector<T, PARTICLE_3D::d> force =
        particle.template getField<FORCING, FORCE>();
    const short sign = (iP == 1 ? 1 : -1);
    EXPECT_NEAR(sign * force[0], 0, expectedNormalForce * 2e-3);
    // deviation of up to 2% is allowed
    EXPECT_NEAR(sign * force[1], expectedNormalForce,
                expectedNormalForce * 2e-2);
    EXPECT_NEAR(sign * force[2], 0, expectedNormalForce * 2e-3);
  }

  EXPECT_EQ(contactContainer.particleContacts.size(), 1);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[0], 0);
  EXPECT_EQ(contactContainer.particleContacts[0].getIDs()[1], 1);
}
