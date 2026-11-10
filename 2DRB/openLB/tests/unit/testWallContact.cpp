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

TEST(wallContactObject, testCopyConstructorAndAssignment)
{
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  WALLCONTACTTYPE                                           contact(79, 43);
  contact.setParticlePosition(PhysR<T, 3>(0.9));
  contact.setResponsibleRank(13);
  contact.setDampingFactor(2.2);
  contact.updateMinMax(PhysR<T, 3>(1.1));
  contact.increaseMinMax(Vector<T, 3>(0.5));
  contact.isNew(false);
  contact.setParticlePositionUpdated(true);

  WALLCONTACTTYPE constructed(contact);
  WALLCONTACTTYPE assigned;
  assigned = contact;

  EXPECT_EQ(contact.getParticleID(), constructed.getParticleID());
  EXPECT_EQ(contact.getWallID(), constructed.getWallID());
  EXPECT_EQ(contact.getResponsibleRank(), constructed.getResponsibleRank());
  EXPECT_EQ(contact.getParticlePosition(), constructed.getParticlePosition());
  EXPECT_EQ(contact.getDampingFactor(), constructed.getDampingFactor());
  EXPECT_EQ(contact.getMin(), constructed.getMin());
  EXPECT_EQ(contact.getMax(), constructed.getMax());
  EXPECT_EQ(contact.isNew(), constructed.isNew());
  EXPECT_EQ(contact.isParticlePositionUpdated(),
            constructed.isParticlePositionUpdated());

  EXPECT_EQ(contact.getParticleID(), assigned.getParticleID());
  EXPECT_EQ(contact.getWallID(), assigned.getWallID());
  EXPECT_EQ(contact.getResponsibleRank(), assigned.getResponsibleRank());
  EXPECT_EQ(contact.getParticlePosition(), assigned.getParticlePosition());
  EXPECT_EQ(contact.getDampingFactor(), assigned.getDampingFactor());
  EXPECT_EQ(contact.getMin(), assigned.getMin());
  EXPECT_EQ(contact.getMax(), assigned.getMax());
  EXPECT_EQ(contact.isNew(), assigned.isNew());
  EXPECT_EQ(contact.isParticlePositionUpdated(),
            assigned.isParticlePositionUpdated());
}

TEST(wallContactViaOverlap3D, testSphere)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 2>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr unsigned wallMaterial        = 1;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T radius                   = 25;
  const T nCellX                   = 48. / converter3D.getPhysDeltaX();
  const T dampingConstant          = 0;
  const T eps                      = 0.;
  const T expectedIndentationDepth = 2;
  const T wallThickness            = 5 * expectedIndentationDepth;
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCellX, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  Vector<T, 3> wallExtend = {4 * radius, wallThickness, 4 * radius};
  Vector<T, 3> wallOrigin = {0., 0., 0.};
  std::unique_ptr<IndicatorF3D<T>> wallInd =
      std::make_unique<IndicatorCuboid3D<T>>(wallExtend, wallOrigin);
  superGeometry.rename(1, 3, wallInd);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  std::vector<SolidBoundary<T, 3>> solidBoundaries;
  solidBoundaries.push_back(SolidBoundary<T, 3>(std::move(wallInd), 3, 1));
  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create sphere 1
  Vector<T, 3> pos1(radius, radius + wallThickness - expectedIndentationDepth,
                    radius);
  particles::creators::addResolvedSphere3D(particleSystem, pos1, radius, eps);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  contactProperties.set(particleMaterial, wallMaterial, effectiveE,
                        dampingConstant, 0, 0);

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D,
      solidBoundaries);

  // Force calculation
  particles::contact::particle_wall<T, PARTICLE_3D, PARTICLECONTACTTYPE,
                                    WALLCONTACTTYPE>::apply(particleSystem,
                                                            solidBoundaries,
                                                            contactContainer,
                                                            contactProperties,
                                                            superGeometry,
                                                            dataMap, 16);

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T alpha               = 3. / 2.;
  const T effectiveRadius     = radius;
  const T expectedNormalForce = 4. / 3. * effectiveE *
                                std::pow(effectiveRadius, 2 - alpha) *
                                std::pow(expectedIndentationDepth, alpha);
  // evaluates to 12.57078722109418

  auto                      particle = particleSystem.get(0);
  Vector<T, PARTICLE_3D::d> force =
      particle.template getField<FORCING, FORCE>();
  EXPECT_NEAR(force[0], 0, expectedNormalForce * 5e-3);
  // deviation of up to 3% is allowed
  EXPECT_NEAR(force[1], expectedNormalForce, expectedNormalForce * 3e-2);
  EXPECT_NEAR(force[2], 0, expectedNormalForce * 4e-3);

  EXPECT_EQ(contactContainer.wallContacts.size(), 1);
  EXPECT_EQ(contactContainer.wallContacts[0].getParticlePosition(), pos1);
  EXPECT_EQ(contactContainer.wallContacts[0].getParticleID(), 0);
  EXPECT_EQ(contactContainer.wallContacts[0].getWallID(), 0);
}

TEST(wallContactViaOverlap3D, testEnlargedSphere)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 2>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr unsigned wallMaterial        = 1;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T radius                        = 25;
  const T nCellX                        = 48. / converter3D.getPhysDeltaX();
  const T particleEnlargementForContact = 1.5 * converter3D.getPhysDeltaX();
  const T dampingConstant               = 0;
  const T eps                           = 0.;
  const T expectedIndentationDepth      = 1.75;
  const T wallThickness                 = 5 * expectedIndentationDepth;
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCellX, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  Vector<T, 3> wallExtend = {4 * radius, wallThickness, 4 * radius};
  Vector<T, 3> wallOrigin = {0., 0., 0.};
  std::unique_ptr<IndicatorF3D<T>> wallInd =
      std::make_unique<IndicatorCuboid3D<T>>(wallExtend, wallOrigin);
  superGeometry.rename(1, 3, wallInd);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  std::vector<SolidBoundary<T, 3>> solidBoundaries;
  solidBoundaries.push_back(SolidBoundary<T, 3>(std::move(wallInd), 3, 1,
                                                particleEnlargementForContact));
  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create sphere 1
  Vector<T, 3> pos1(radius, radius + wallThickness - expectedIndentationDepth,
                    radius);
  particles::creators::addResolvedSphere3D(particleSystem, pos1, radius, eps);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  particleA.setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(
      particleEnlargementForContact);
  contactProperties.set(particleMaterial, wallMaterial, effectiveE,
                        dampingConstant, 0, 0);

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D,
      solidBoundaries, []() {
        return Vector<bool, 3>(false, false, false);
      });

  // Force calculation
  particles::contact::particle_wall<T, PARTICLE_3D, PARTICLECONTACTTYPE,
                                    WALLCONTACTTYPE>::apply(particleSystem,
                                                            solidBoundaries,
                                                            contactContainer,
                                                            contactProperties,
                                                            superGeometry,
                                                            dataMap, 16,
                                                            4. / (3 *
                                                                  util::sqrt(
                                                                      M_PI)));

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T alpha           = 3. / 2.;
  const T effectiveRadius = radius + particleEnlargementForContact;
  const T expectedNormalForce =
      4. / 3. * effectiveE * std::pow(effectiveRadius, 2 - alpha) *
      std::pow(expectedIndentationDepth + 2 * particleEnlargementForContact,
               alpha);

  auto                      particle = particleSystem.get(0);
  Vector<T, PARTICLE_3D::d> force =
      particle.template getField<FORCING, FORCE>();
  EXPECT_NEAR(force[0], 0, expectedNormalForce * 5e-3);
  // deviation of up to 2.5% is allowed
  EXPECT_NEAR(force[1], expectedNormalForce, expectedNormalForce * 2.5e-2);
  EXPECT_NEAR(force[2], 0, expectedNormalForce * 4e-3);

  EXPECT_EQ(contactContainer.wallContacts.size(), 1);
  EXPECT_EQ(contactContainer.wallContacts[0].getParticlePosition(), pos1);
  EXPECT_EQ(contactContainer.wallContacts[0].getParticleID(), 0);
  EXPECT_EQ(contactContainer.wallContacts[0].getWallID(), 0);
}
TEST(wallContactViaOverlap3D, testSpherePeriodicBoundary)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 2>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr unsigned wallMaterial        = 1;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T radius                   = 25;
  const T nCellX                   = 60. / converter3D.getPhysDeltaX();
  const T dampingConstant          = 0;
  const T eps                      = 0.;
  const T expectedIndentationDepth = 2;
  const T wallThickness            = 5 * expectedIndentationDepth;
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCellX, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  Vector<T, 3> wallExtend = {6 * radius, wallThickness, 6 * radius};
  Vector<T, 3> wallOrigin = {-.5 * wallExtend[0], 0., -.5 * wallExtend[0]};
  std::unique_ptr<IndicatorF3D<T>> wallInd =
      std::make_unique<IndicatorCuboid3D<T>>(wallExtend, wallOrigin);
  superGeometry.rename(1, 3, wallInd);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  std::vector<SolidBoundary<T, 3>> solidBoundaries;
  solidBoundaries.push_back(SolidBoundary<T, 3>(std::move(wallInd), 3, 1));
  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create sphere 1
  Vector<T, 3> pos1(0., radius + wallThickness - expectedIndentationDepth, 0.);
  particles::creators::addResolvedSphere3D(particleSystem, pos1, radius, eps);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  contactProperties.set(particleMaterial, wallMaterial, effectiveE,
                        dampingConstant, 0, 0);

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D,
      solidBoundaries, []() {
        return Vector<bool, 3>(true, false, true);
      });

  // Force calculation
  particles::contact::particle_wall<T, PARTICLE_3D, PARTICLECONTACTTYPE,
                                    WALLCONTACTTYPE>::apply(particleSystem,
                                                            solidBoundaries,
                                                            contactContainer,
                                                            contactProperties,
                                                            superGeometry,
                                                            dataMap, 16);

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T alpha               = 3. / 2.;
  const T effectiveRadius     = radius;
  const T expectedNormalForce = 4. / 3. * effectiveE *
                                std::pow(effectiveRadius, 2 - alpha) *
                                std::pow(expectedIndentationDepth, alpha);
  // evaluates to 12.57078722109418

  auto                      particle = particleSystem.get(0);
  Vector<T, PARTICLE_3D::d> force =
      particle.template getField<FORCING, FORCE>();
  EXPECT_NEAR(force[0], 0, expectedNormalForce * 5e-3);
  // deviation of up to 4% is allowed
  EXPECT_NEAR(force[1], expectedNormalForce, expectedNormalForce * 4e-2);
  EXPECT_NEAR(force[2], 0, expectedNormalForce * 4e-3);
  EXPECT_EQ(contactContainer.wallContacts.size(), 1);
  EXPECT_EQ(contactContainer.wallContacts[0].getParticlePosition(), pos1);
  EXPECT_EQ(contactContainer.wallContacts[0].getParticleID(), 0);
  EXPECT_EQ(contactContainer.wallContacts[0].getWallID(), 0);
}

TEST(wallContactViaOverlap3D, testCylinder)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 2>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr unsigned wallMaterial        = 1;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T      radius                   = 25;
  const T      nCellX                   = 48. / converter3D.getPhysDeltaX();
  const T      dampingConstant          = 0;
  const T      eps                      = 0.;
  const T      expectedIndentationDepth = 2;
  const T      height                   = 50;
  const T      density                  = 1.;
  Vector<T, 3> angle(0);
  const T      wallThickness = expectedIndentationDepth;
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCellX, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  Vector<T, 3> wallExtend = {2 * radius, 2 * wallThickness, 2 * radius};
  Vector<T, 3> wallOrigin = {0., 0., 0.};
  std::unique_ptr<IndicatorF3D<T>> wallInd =
      std::make_unique<IndicatorCuboid3D<T>>(wallExtend, wallOrigin);
  superGeometry.rename(1, 3, wallInd);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  std::vector<SolidBoundary<T, 3>> solidBoundaries;
  solidBoundaries.push_back(SolidBoundary<T, 3>(std::move(wallInd), 3, 1));
  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create cylinder 1
  Vector<T, 3> pos1(radius, wallThickness + height * 0.5, radius);
  Vector<T, 3> normal(0, 1, 0);
  particles::creators::addResolvedCylinder3D(
      particleSystem, pos1, normal, height, radius, eps, density, angle);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  contactProperties.set(particleMaterial, wallMaterial, effectiveE,
                        dampingConstant, 0, 0);

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D,
      solidBoundaries);

  // Force calculation
  particles::contact::particle_wall<T, PARTICLE_3D, PARTICLECONTACTTYPE,
                                    WALLCONTACTTYPE>::apply(particleSystem,
                                                            solidBoundaries,
                                                            contactContainer,
                                                            contactProperties,
                                                            superGeometry,
                                                            dataMap, 16,
                                                            2. / util::sqrt(
                                                                     M_PI));

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T expectedNormalForce =
      2 * radius * effectiveE * expectedIndentationDepth;
  // evaluates to 66.666666666666657

  auto                      particle = particleSystem.get(0);
  Vector<T, PARTICLE_3D::d> force =
      particle.template getField<FORCING, FORCE>();
  EXPECT_NEAR(force[0], 0, expectedNormalForce * 2e-3);
  // deviation of up to 2% is allowed
  EXPECT_NEAR(force[1], expectedNormalForce, expectedNormalForce * 2e-2);
  EXPECT_NEAR(force[2], 0, expectedNormalForce * 2e-3);

  EXPECT_EQ(contactContainer.wallContacts.size(), 1);
  EXPECT_EQ(contactContainer.wallContacts[0].getParticlePosition(), pos1);
  EXPECT_EQ(contactContainer.wallContacts[0].getParticleID(), 0);
  EXPECT_EQ(contactContainer.wallContacts[0].getWallID(), 0);
}

TEST(wallContactViaOverlap3D, testCylinderPeriodicBoundary)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 2>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr unsigned wallMaterial        = 1;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T      radius                   = 25;
  const T      nCellX                   = 60. / converter3D.getPhysDeltaX();
  const T      dampingConstant          = 0;
  const T      eps                      = 0.;
  const T      expectedIndentationDepth = 2;
  const T      height                   = 50;
  const T      density                  = 1.;
  Vector<T, 3> angle(0);
  const T      wallThickness = expectedIndentationDepth;
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCellX, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  Vector<T, 3> wallExtend = {4 * radius, 2 * wallThickness, 4 * radius};
  Vector<T, 3> wallOrigin = {-.5 * wallExtend[0], 0., -.5 * wallExtend[0]};
  std::unique_ptr<IndicatorF3D<T>> wallInd =
      std::make_unique<IndicatorCuboid3D<T>>(wallExtend, wallOrigin);
  superGeometry.rename(1, 3, wallInd);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  std::vector<SolidBoundary<T, 3>> solidBoundaries;
  solidBoundaries.push_back(SolidBoundary<T, 3>(std::move(wallInd), 3, 1));
  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create cylinder 1
  Vector<T, 3> pos1(converter3D.getPhysDeltaX(), wallThickness + height * 0.5,
                    converter3D.getPhysDeltaX());
  Vector<T, 3> normal(0, 1, 0);
  particles::creators::addResolvedCylinder3D(
      particleSystem, pos1, normal, height, radius, eps, density, angle);
  auto particleA = particleSystem.get(0);
  particleA.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  contactProperties.set(particleMaterial, wallMaterial, effectiveE,
                        dampingConstant, 0, 0);

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D,
      solidBoundaries, []() {
        return Vector<bool, 3>(true, false, true);
      });

  // Force calculation
  particles::contact::particle_wall<T, PARTICLE_3D, PARTICLECONTACTTYPE,
                                    WALLCONTACTTYPE>::apply(particleSystem,
                                                            solidBoundaries,
                                                            contactContainer,
                                                            contactProperties,
                                                            superGeometry,
                                                            dataMap, 16,
                                                            2. / util::sqrt(
                                                                     M_PI));

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T expectedNormalForce =
      2 * radius * effectiveE * expectedIndentationDepth;
  // evaluates to 66.666666666666657

  auto                      particle = particleSystem.get(0);
  Vector<T, PARTICLE_3D::d> force =
      particle.template getField<FORCING, FORCE>();
  EXPECT_NEAR(force[0], 0, expectedNormalForce * 2e-3);
  // deviation of up to 2% is allowed
  EXPECT_NEAR(force[1], expectedNormalForce, expectedNormalForce * 2e-2);
  EXPECT_NEAR(force[2], 0, expectedNormalForce * 2e-3);
  EXPECT_EQ(contactContainer.wallContacts.size(), 1);
  EXPECT_EQ(contactContainer.wallContacts[0].getParticlePosition(), pos1);
  EXPECT_EQ(contactContainer.wallContacts[0].getParticleID(), 0);
  EXPECT_EQ(contactContainer.wallContacts[0].getWallID(), 0);
}

//Cone tests

TEST(wallContactViaOverlap3D, testWallCone)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, 3, true>
                                                            PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, 3, true> WALLCONTACTTYPE;
  typedef ContactProperties<T, 2>                           CONTACTPROPERTIES;

  constexpr unsigned particleMaterial    = 0;
  constexpr unsigned wallMaterial        = 1;
  constexpr T        modulusOfElasticity = 1;
  constexpr T        poissonRatio        = 0.5;
  const T            effectiveE = particles::contact::evalEffectiveYoungModulus(
      modulusOfElasticity, modulusOfElasticity, poissonRatio, poissonRatio);

  std::multimap<int, std::unique_ptr<std::uint8_t[]>> dataMap;

  const T      radius1                  = 25;
  const T      radius2                  = 0;
  const T      nCellX                   = 48. / converter3D.getPhysDeltaX();
  const T      dampingConstant          = 0;
  const T      eps                      = 0.;
  const T      expectedIndentationDepth = 2;
  const T      density                  = 1;
  const T      height                   = 50;
  Vector<T, 3> angle(T(0));
  const T      wallThickness = 25;
  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCellX, 2);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);

  Vector<T, 3> wallExtend = {4 * radius1, wallThickness, 4 * radius1};
  Vector<T, 3> wallOrigin = {-.5 * wallExtend[0], 0., -.5 * wallExtend[0]};
  std::unique_ptr<IndicatorF3D<T>> wallInd =
      std::make_unique<IndicatorCuboid3D<T>>(wallExtend, wallOrigin);
  superGeometry.rename(1, 3, wallInd);
  std::vector<SolidBoundary<T, 3>> solidBoundaries;
  solidBoundaries.push_back(SolidBoundary<T, 3>(std::move(wallInd), 3, 1));

  SuperLattice<T, D3D>           sLattice(converter3D,superGeometry);
  ParticleSystem<T, PARTICLE_3D> particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  CONTACTPROPERTIES                                         contactProperties;

  //Create cone
  Vector<T, 3> center1(25, wallExtend[1] + height - expectedIndentationDepth,
                       25);
  Vector<T, 3> center2(25, wallExtend[1] - expectedIndentationDepth, 25);
  particles::creators::addResolvedCone3D(particleSystem, center1, center2,
                                         radius1, radius2, eps, density, angle);
  auto particle = particleSystem.get(0);
  particle.setField<MECHPROPERTIES, MATERIAL>(particleMaterial);
  contactProperties.set(particleMaterial, wallMaterial, effectiveE,
                        dampingConstant, 0, 0);

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D,
      solidBoundaries);

  // Force calculation
  particles::contact::particle_wall<T, PARTICLE_3D, PARTICLECONTACTTYPE,
                                    WALLCONTACTTYPE>::apply(particleSystem,
                                                            solidBoundaries,
                                                            contactContainer,
                                                            contactProperties,
                                                            superGeometry,
                                                            dataMap, 16, 0.65);

  // Calculation of expected result via Hertz theory (for E_i = E_j = 1 and sigma_i = sigma_j = 0.5)
  const T inclinationAngle   = 0.5 * M_PI - util::atan(radius1 / height);
  const T contactRegionDepth = expectedIndentationDepth * 2 / M_PI;
  const T contactRadius = contactRegionDepth / (util::tan(inclinationAngle));
  const T expectedNormalForce = 0.5 * M_PI * util::pow(contactRadius, 2) *
                                effectiveE * util::tan(inclinationAngle);
  // evaluates to 3.3953054526271012

  Vector<T, PARTICLE_3D::d> force =
      particle.template getField<FORCING, FORCE>();
  EXPECT_NEAR(force[0], 0, expectedNormalForce * 2e-3);
  // deviation of up to 2% is allowed
  EXPECT_NEAR(force[1], expectedNormalForce, expectedNormalForce * 2e-2);
  EXPECT_NEAR(force[2], 0, expectedNormalForce * 2e-3);
  EXPECT_EQ(contactContainer.wallContacts.size(), 1);
  EXPECT_EQ(contactContainer.wallContacts[0].getParticleID(), 0);
  EXPECT_EQ(contactContainer.wallContacts[0].getWallID(), 0);
}
