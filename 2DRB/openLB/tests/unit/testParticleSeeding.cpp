/*  This file is part of the OpenLB library
*
*  Copyright (C) 2022 Jan E. Marquardt
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

#include <cstdio>
#include <type_traits>

#include "helper/unitTestHelper.h"
#include <gtest/gtest.h>

#include <olb.h>

using namespace olb;
using namespace descriptors;
typedef double T;

// #define VERBOSE_TESTS

/// 3D Versions

TEST(TestResolvedParticleSeeding, CuboidDomainWithSpheres3D)
{
  using namespace particles;
  using namespace particles::creators;
  typedef ResolvedSphere3D        PARTICLETYPE;
  constexpr unsigned short        D = PARTICLETYPE::d;
  ParticleSystem<T, PARTICLETYPE> particleSystem;

  // Domain settings
  const Vector<T, D>   origin(-0.5);
  const Vector<T, D>   extend(1.0);
  IndicatorCuboid3D<T> cuboid(extend, origin);
  const T              domainVolume = extend[0] * extend[1] * extend[2];

  // Particle settings
  const T                       wantedParticleVolumeFraction = 0.2;
  const T                       sphereRadius                 = 0.02;
  const T                       epsilon = 0.0 * sphereRadius;
  SmoothIndicatorSphere3D<T, T> testSphere(PhysR<T, D>(), sphereRadius,
                                           epsilon);
  const std::function<T(const std::size_t& pID)> getCircumRadius =
      [&](const std::size_t& pID) {
        return testSphere.getCircumRadius();
      };

  const std::function<T(const std::size_t& pID)> getVolume =
      [&](const std::size_t& pID) {
        return 4. / 3. * M_PI * util::pow(sphereRadius, 3);
      };

  const std::function<void(const SpawnData<T, D>&, const std::size_t&)>
      createParticle =
          [&](const SpawnData<T, D>& spawnData, const std::size_t& pID) {
            creators::addResolvedSphere3D(particleSystem, spawnData.position,
                                          sphereRadius, epsilon);
          };

  std::vector<SpawnData<T, D>> spawnData = creators::setParticles<T, D>(
      wantedParticleVolumeFraction, cuboid, domainVolume, getCircumRadius,
      getVolume, createParticle);
  const T particleVolumeFraction =
      creators::calcParticleVolumeFraction(spawnData, domainVolume, getVolume);

#ifdef VERBOSE_TESTS
  std::cout << "particleVolumeFraction = " << particleVolumeFraction
            << std::endl;
#endif

  EXPECT_TRUE(particleVolumeFraction >= wantedParticleVolumeFraction);

  // Check that particles do not overlap with other particles
  std::size_t pID = 0;
  for (const SpawnData<T, D>& entry : spawnData) {
    for (unsigned i = ++pID; i < spawnData.size(); ++i) {
      T distanceCenterOfMass = norm(entry.position - spawnData[i].position);
      EXPECT_TRUE(distanceCenterOfMass >=
                  getCircumRadius(pID) + getCircumRadius(i));
    }
  }

  // Check that particles do not overlap with the wall
  pID = 0;
  for (const SpawnData<T, D>& entry : spawnData) {
    std::array<Vector<T, D>, 26> pointsOnHull =
        discretePointsOnSphericalHull(entry.position, getCircumRadius(pID));
    for (const Vector<T, D>& pointOnHull : pointsOnHull) {
      bool isInside;
      cuboid(&isInside, pointOnHull.data());
      EXPECT_TRUE(isInside);
    }
  }
}

TEST(TestResolvedParticleSeeding, CuboidDomainWithCuboids3D)
{
  using namespace particles;
  using namespace particles::creators;

  typedef ResolvedParticle3D      PARTICLETYPE;
  constexpr unsigned short        D = PARTICLETYPE::d;
  ParticleSystem<T, PARTICLETYPE> particleSystem;

  // Domain settings
  const Vector<T, D>   origin(-0.5);
  const Vector<T, D>   extend(1.0);
  IndicatorCuboid3D<T> cuboid(extend, origin);
  const T              domainVolume = extend[0] * extend[1] * extend[2];

  // Particle settings
  const T                                 wantedParticleVolumeFraction = 0.1;
  const T                                 cubeExtend                   = 0.02;
  const T                                 epsilon                      = 0.0;
  SmoothIndicatorCuboid3D<double, double> testCuboid(
      cubeExtend, cubeExtend, cubeExtend, PhysR<T, D>(), epsilon);
  const std::function<T(const std::size_t& pID)> getCircumRadius =
      [&](const std::size_t& pID) {
        return testCuboid.getCircumRadius();
      };

  const std::function<T(const std::size_t& pID)> getVolume =
      [&](const std::size_t& pID) {
        return util::pow(cubeExtend, D);
      };

  const std::function<void(const SpawnData<T, D>&, const std::size_t&)>
      createParticle =
          [&](const SpawnData<T, D>& spawnData, const std::size_t& pID) {
            creators::addResolvedCuboid3D(particleSystem, spawnData.position,
                                          Vector<T, D>(cubeExtend), epsilon);
          };

  std::vector<SpawnData<T, D>> spawnData = creators::setParticles<T, D>(
      wantedParticleVolumeFraction, cuboid, domainVolume, getCircumRadius,
      getVolume, createParticle);
  const T particleVolumeFraction =
      creators::calcParticleVolumeFraction(spawnData, domainVolume, getVolume);

#ifdef VERBOSE_TESTS
  std::cout << "particleVolumeFraction = " << particleVolumeFraction
            << std::endl;
#endif

  EXPECT_TRUE(particleVolumeFraction >= wantedParticleVolumeFraction);

  // Check that circum radius does not overlap
  std::size_t pID = 0;
  for (const SpawnData<T, D>& entry : spawnData) {
    for (unsigned i = ++pID; i < spawnData.size(); ++i) {
      T distanceCenterOfMass = norm(entry.position - spawnData[i].position);
      EXPECT_TRUE(distanceCenterOfMass >
                  getCircumRadius(pID) + getCircumRadius(i));
    }
  }

  // Check that particles do not overlap with the wall
  pID = 0;
  for (const SpawnData<T, D>& entry : spawnData) {
    std::array<Vector<T, D>, 26> pointsOnHull =
        discretePointsOnSphericalHull(entry.position, getCircumRadius(pID));
    for (const Vector<T, D>& pointOnHull : pointsOnHull) {
      bool isInside;
      cuboid(&isInside, pointOnHull.data());
      EXPECT_TRUE(isInside);
    }
  }
}

TEST(TestResolvedParticleSeeding, SaveAndReadParticlePositions3D)
{
  using namespace particles;
  using namespace particles::creators;

  constexpr unsigned           D        = 3;
  const std::string            filename = "test.csv";
  std::vector<SpawnData<T, D>> spawnData;
  const double epsilon = 1e-16;

  //Create randomizer
  util::Randomizer<T> randomizer;
  for (unsigned i = 0; i < 100; ++i) {
    const PhysR<T, D> randomPosition =
        randomizer.template generate<PhysR<T, D>>();
    spawnData.push_back(SpawnData<T, D>(randomPosition, Vector<T, D>(T {0.})));
  }

  const std::function<std::string(const std::size_t&)> evalIdentifier =
      [](const std::size_t& pID) {
        return "id";
      };
  saveParticlePositions(filename, spawnData, evalIdentifier);
  std::vector<std::vector<std::string>> particlePositionsRead =
      particles::creators::readParticlePositions(filename);
  std::remove(filename.c_str());

  std::size_t i = 0;
  for (const std::vector<std::string>& line : particlePositionsRead) {
    unsigned short j = 0;
    for (const std::string& word : line) {
      if (j < D) {
        EXPECT_NEAR(std::stod(word), spawnData[i].position[j], epsilon);
      }
      else if (j < D + D) {
        EXPECT_NEAR(std::stod(word), spawnData[i].angleInDegree[j - D], epsilon);
      }
      else {
        EXPECT_EQ(word, std::string("id"));
      }
      ++j;
    }
    ++i;
  }
}


/// 2D Versions

TEST(TestResolvedParticleSeeding, CuboidDomainWithSpheres2D)
{
  using namespace particles;
  using namespace particles::creators;

  typedef ResolvedCircle2D        PARTICLETYPE;
  constexpr unsigned short        D = PARTICLETYPE::d;
  ParticleSystem<T, PARTICLETYPE> particleSystem;

  // Domain settings
  const Vector<T, D>   origin(-0.5);
  const Vector<T, D>   extend(1.0);
  IndicatorCuboid2D<T> cuboid(extend, origin);
  const T              domainVolume = extend[0] * extend[1];

  // Particle settings
  const T                       wantedParticleVolumeFraction = 0.2;
  const T                       sphereRadius                 = 0.02;
  const T                       epsilon = 0.0 * sphereRadius;
  SmoothIndicatorCircle2D<T, T> testSphere(PhysR<T, D>(), sphereRadius,
                                           epsilon);
  const std::function<T(const std::size_t& pID)> getCircumRadius =
      [&](const std::size_t& pID) {
        return testSphere.getCircumRadius();
      };

  const std::function<T(const std::size_t& pID)> getVolume =
      [&](const std::size_t& pID) {
        return M_PI * util::pow(sphereRadius, 2);
      };

  const std::function<void(const SpawnData<T, D>&, const std::size_t&)>
      createParticle =
          [&](const SpawnData<T, D>& spawnData, const std::size_t& pID) {
            creators::addResolvedCircle2D(particleSystem, spawnData.position,
                                          sphereRadius, epsilon);
          };

  std::vector<SpawnData<T, D>> spawnData = creators::setParticles<T, D>(
      wantedParticleVolumeFraction, cuboid, domainVolume, getCircumRadius,
      getVolume, createParticle);
  const T particleVolumeFraction =
      creators::calcParticleVolumeFraction(spawnData, domainVolume, getVolume);

#ifdef VERBOSE_TESTS
  std::cout << "particleVolumeFraction = " << particleVolumeFraction
            << std::endl;
#endif

  EXPECT_TRUE(particleVolumeFraction >= wantedParticleVolumeFraction);

  // Check that particles do not overlap with other particles
  std::size_t pID = 0;
  for (const SpawnData<T, D>& entry : spawnData) {
    for (unsigned i = ++pID; i < spawnData.size(); ++i) {
      T distanceCenterOfMass = norm(entry.position - spawnData[i].position);
      EXPECT_TRUE(distanceCenterOfMass >=
                  getCircumRadius(pID) + getCircumRadius(i));
    }
  }

  // Check that particles do not overlap with the wall
  pID = 0;
  for (const SpawnData<T, D>& entry : spawnData) {
    std::array<Vector<T, D>, 8> pointsOnHull =
        discretePointsOnSphericalHull(entry.position, getCircumRadius(pID));
    for (const Vector<T, D>& pointOnHull : pointsOnHull) {
      bool isInside;
      cuboid(&isInside, pointOnHull.data());
      EXPECT_TRUE(isInside);
    }
  }
}

TEST(TestResolvedParticleSeeding, CuboidDomainWithCuboids2D)
{
  using namespace particles;
  using namespace particles::creators;

  typedef ResolvedParticle2D      PARTICLETYPE;
  constexpr unsigned short        D = PARTICLETYPE::d;
  ParticleSystem<T, PARTICLETYPE> particleSystem;

  // Domain settings
  const Vector<T, D>   origin(-0.5);
  const Vector<T, D>   extend(1.0);
  IndicatorCuboid2D<T> cuboid(extend, origin);
  const T              domainVolume = extend[0] * extend[1];

  // Particle settings
  const T                                 wantedParticleVolumeFraction = 0.2;
  const T                                 cubeExtend                   = 0.02;
  const T                                 epsilon                      = 0.0;
  SmoothIndicatorCuboid2D<double, double> testCuboid(PhysR<T, D>(), cubeExtend,
                                                     cubeExtend, epsilon);
  const std::function<T(const std::size_t& pID)> getCircumRadius =
      [&](const std::size_t& pID) {
        return testCuboid.getCircumRadius();
      };

  const std::function<T(const std::size_t& pID)> getVolume =
      [&](const std::size_t& pID) {
        return util::pow(cubeExtend, D);
      };

  const std::function<void(const SpawnData<T, D>&, const std::size_t&)>
      createParticle =
          [&](const SpawnData<T, D>& spawnData, const std::size_t& pID) {
            creators::addResolvedCuboid2D(particleSystem, spawnData.position,
                                          Vector<T, D>(cubeExtend), epsilon);
          };

  std::vector<SpawnData<T, D>> spawnData = creators::setParticles<T, D>(
      wantedParticleVolumeFraction, cuboid, domainVolume, getCircumRadius,
      getVolume, createParticle);
  const T particleVolumeFraction =
      creators::calcParticleVolumeFraction(spawnData, domainVolume, getVolume);

#ifdef VERBOSE_TESTS
  std::cout << "particleVolumeFraction = " << particleVolumeFraction
            << std::endl;
#endif

  EXPECT_TRUE(particleVolumeFraction >= wantedParticleVolumeFraction);

  // Check that circum radius does not overlap
  std::size_t pID = 0;
  for (const SpawnData<T, D>& entry : spawnData) {
    for (unsigned i = ++pID; i < spawnData.size(); ++i) {
      T distanceCenterOfMass = norm(entry.position - spawnData[i].position);
      EXPECT_TRUE(distanceCenterOfMass >=
                  getCircumRadius(pID) + getCircumRadius(i));
    }
  }

  // Check that particles do not overlap with the wall
  pID = 0;
  for (const SpawnData<T, D>& entry : spawnData) {
    std::array<Vector<T, D>, 8> pointsOnHull =
        discretePointsOnSphericalHull(entry.position, getCircumRadius(pID));
    for (const Vector<T, D>& pointOnHull : pointsOnHull) {
      bool isInside;
      cuboid(&isInside, pointOnHull.data());
      EXPECT_TRUE(isInside);
    }
  }
}

TEST(TestResolvedParticleSeeding, SaveAndReadParticlePositions2D)
{
  using namespace particles;
  using namespace particles::creators;

  constexpr unsigned           D        = 2;
  const std::string            filename = "test.csv";
  std::vector<SpawnData<T, D>> spawnData;
  const double epsilon = 1e-16;

  //Create randomizer
  util::Randomizer<T> randomizer;
  for (unsigned i = 0; i < 100; ++i) {
    const PhysR<T, D> randomPosition =
        randomizer.template generate<PhysR<T, D>>();
    spawnData.push_back(SpawnData<T, D>(randomPosition, Vector<T, 1>(T {0})));
  }

  const std::function<std::string(const std::size_t&)> evalIdentifier =
      [](const std::size_t& pID) {
        return "id";
      };
  saveParticlePositions(filename, spawnData, evalIdentifier);
  std::vector<std::vector<std::string>> particlePositionsRead =
      particles::creators::readParticlePositions(filename);
  std::remove(filename.c_str());

  std::size_t i = 0;
  for (const std::vector<std::string>& line : particlePositionsRead) {
    unsigned short j = 0;
    for (const std::string& word : line) {
      if (j < D) {
        EXPECT_NEAR(std::stod(word), spawnData[i].position[j], epsilon);
      }
      else if (j < D + 1) {
        EXPECT_NEAR(std::stod(word), spawnData[i].angleInDegree[j - D], epsilon);
      }
      else {
        EXPECT_EQ(word, std::string("id"));
      }
      ++j;
    }
    ++i;
  }
}

