#include "helper/unitTestHelper.h"
#include <gtest/gtest.h>
#include <iostream>

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::particles;

typedef ResolvedParticle3D PARTICLETYPE;
typedef double             T;

T getPorosity(Particle<T, PARTICLETYPE>& particle, T input[])
{
  T tmp[1];
  particles::resolved::evalSolidVolumeFraction(tmp, input, particle);
  return tmp[0];
}

//Sphere
TEST(TestSphere3D, testMomentOfInertiaAndMass)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  T                               radius(2.5), epsilon(0.1), density(0.1);
  particles::creators::addResolvedSphere3D(particleSystem, center, radius,
                                           epsilon, density);
  SmoothIndicatorSphere3D<T, T, false> sind(center, radius, epsilon);

  auto particle = particleSystem.get(0);
  for (unsigned iD = 0; iD < 3; ++iD) {
    const T mofi = particle.template getField<PHYSPROPERTIES, MOFI>()[iD];
    EXPECT_EQ(mofi, sind.calcMofiAndMass(density)[iD]);
  }
  const T mass = particle.template getField<PHYSPROPERTIES, MASS>();
  EXPECT_EQ(mass, sind.calcMofiAndMass(density)[3]);
}

TEST(TestSphere3D, testCenter)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  T                               radius(2.5), epsilon(0.1);
  particles::creators::addResolvedSphere3D(particleSystem, center, radius,
                                           epsilon);
  T    res[1];
  T    input[3] = {center[0], center[1], center[2]};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

TEST(TestSphere3D, testInnerBorder)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  T                               radius(2.5), epsilon(0.1);
  particles::creators::addResolvedSphere3D(particleSystem, center, radius,
                                           epsilon);
  T    res[1];
  T    input[3] = {center[0] + (radius - epsilon / 2.), center[1], center[2]};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

TEST(TestSphere3D, testAtRadius)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  T                               radius(2.5), epsilon(0.1);
  particles::creators::addResolvedSphere3D(particleSystem, center, radius,
                                           epsilon);
  T    res[1];
  T    input[3] = {center[0], center[1] - radius, center[2]};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_LT(0.0, res[0]);
  ASSERT_GT(1.0, res[0]);
}

TEST(TestSphere3D, testOuterBorder)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  T                               radius(2.5), epsilon(0.1);
  particles::creators::addResolvedSphere3D(particleSystem, center, radius,
                                           epsilon);
  T    res[1];
  T    input[3] = {center[0], center[1],
                   center[2] - (radius + epsilon / 2. +
                             2 * std::numeric_limits<T>::epsilon())};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_EQ(0.0, res[0]);
}

TEST(TestSphere3D, testFarOut)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  T                               radius(2.5), epsilon(0.1);
  particles::creators::addResolvedSphere3D(particleSystem, center, radius,
                                           epsilon);
  T    res[1];
  T    input[3] = {center[0] + 1000 * epsilon, center[1] - 500 * epsilon,
                   center[2] - 1000 * epsilon};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_EQ(0.0, res[0]);
}

TEST(TestSphere3D, testTranslatedCenter)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {20, 1, 4.3};
  T                               radius(2.5), epsilon(0.1);
  particles::creators::addResolvedSphere3D(particleSystem, center, radius,
                                           epsilon);
  T    res[1];
  T    input[3] = {center[0], center[1], center[2]};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

TEST(TestSphere3D, testTranslatedInnerBorder)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0.3, 80, -5};
  T                               radius(2.5), epsilon(0.1);
  particles::creators::addResolvedSphere3D(particleSystem, center, radius,
                                           epsilon);
  T    res[1];
  T    input[3] = {center[0] + (radius - epsilon / 2.), center[1], center[2]};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

TEST(TestSphere3D, testTranslatedAtRadius)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {-0.9, 1, 5};
  T                               radius(2.5), epsilon(0.1);
  particles::creators::addResolvedSphere3D(particleSystem, center, radius,
                                           epsilon);
  T    res[1];
  T    input[3] = {center[0], center[1] - radius, center[2]};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_LT(0.0, res[0]);
  ASSERT_GT(1.0, res[0]);
}

TEST(TestSphere3D, testTranslatedOuterBorder)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {-48.5, 22.5, 1};
  T                               radius(2.5), epsilon(0.1);
  particles::creators::addResolvedSphere3D(particleSystem, center, radius,
                                           epsilon);
  T    res[1];
  T    input[3] = {center[0], center[1],
                   center[2] - (radius + epsilon / 2. +
                             2 * std::numeric_limits<T>::epsilon())};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_EQ(0.0, res[0]);
}

TEST(TestSphere3D, testTranslatedFarOut)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {30.2, -400, 400};
  T                               radius(2.5), epsilon(0.1);
  particles::creators::addResolvedSphere3D(particleSystem, center, radius,
                                           epsilon);
  T    res[1];
  T    input[3] = {center[0] + 1000 * epsilon, center[1] - 500 * epsilon,
                   center[2] - 1000 * epsilon};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_EQ(0.0, res[0]);
}

TEST(TestSphere3D, testOuterBorderLargeRadius)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  T                               radius(20), epsilon(0.1);
  particles::creators::addResolvedSphere3D(particleSystem, center, radius,
                                           epsilon);
  T    res[1];
  T    input[3] = {center[0], center[1], center[2] - (radius + epsilon / 2.)};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_EQ(0.0, res[0]);
}

TEST(TestSphere3D, testCenterLargeEpsilon)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  T                               radius(2.5), epsilon(4);
  particles::creators::addResolvedSphere3D(particleSystem, center, radius,
                                           epsilon);
  T    res[1];
  T    input[3] = {center[0], center[1], center[2]};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

//Cuboid
TEST(TestCuboid3D, testMomentOfInertiaAndMass)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  Vector<T, 3>                    extend {20, 2, 2};
  T                               epsilon(0.5), density(1e2);
  particles::creators::addResolvedCuboid3D(particleSystem, center, extend,
                                           epsilon, density);
  SmoothIndicatorCuboid3D<T, T, false> sind(extend[0], extend[1], extend[2],
                                            center, epsilon);

  auto particle = particleSystem.get(0);
  for (unsigned iD = 0; iD < 3; ++iD) {
    const T mofi = particle.template getField<PHYSPROPERTIES, MOFI>()[iD];
    EXPECT_EQ(mofi, sind.calcMofiAndMass(density)[iD]);
  }
  const T mass = particle.template getField<PHYSPROPERTIES, MASS>();
  EXPECT_EQ(mass, sind.calcMofiAndMass(density)[3]);
}

TEST(TestCuboid3D, testCenter)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  Vector<T, 3>                    extend {20, 2, 2};
  T                               epsilon(0.5);
  particles::creators::addResolvedCuboid3D(particleSystem, center, extend,
                                           epsilon);
  T    res[1];
  T    input[3] = {center[0], center[1], center[2]};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

TEST(TestCuboid3D, testInnerBorder)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  Vector<T, 3>                    extend {20, 2, 2};
  T                               epsilon(0.5);
  particles::creators::addResolvedCuboid3D(particleSystem, center, extend,
                                           epsilon);
  T    res[1];
  T    input[3] = {center[0] + (extend[0] / 2 - epsilon / 2.), center[1],
                   center[2]};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

TEST(TestCuboid3D, testatLength)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  Vector<T, 3>                    extend {20, 2, 2};
  T                               epsilon(0.5);
  particles::creators::addResolvedCuboid3D(particleSystem, center, extend,
                                           epsilon);
  T    res[1];
  T    input[3] = {center[0], center[1] - extend[1] / 2, center[2]};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_LT(0.0, res[0]);
  ASSERT_GT(1.0, res[0]);
}

TEST(TestCuboid3D, testouterBorder)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  Vector<T, 3>                    extend {20, 2, 2};
  T                               epsilon(0.5);
  particles::creators::addResolvedCuboid3D(particleSystem, center, extend,
                                           epsilon);
  T    res[1];
  T    input[3] = {center[0], center[1],
                   center[2] - (extend[2] / 2 + epsilon / 2. +
                             2 * std::numeric_limits<T>::epsilon())};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_EQ(0, res[0]);
}

TEST(TestCuboid3D, testTranslateCenter)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {50, 40, 30};
  Vector<T, 3>                    extend {20, 2, 2};
  T                               epsilon(0.5);
  particles::creators::addResolvedCuboid3D(particleSystem, center, extend,
                                           epsilon);
  T    res[1];
  T    input[3] = {center[0], center[1], center[2]};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

TEST(TestCuboid3D, testRotateSmoothInd)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  Vector<T, 3>                    extend {20, 2, 2};
  Vector<T, 3>                    theta {0, 0, 45};
  T                               epsilon(0.5);
  particles::creators::addResolvedCuboid3D(particleSystem, center, extend,
                                           epsilon, 0., theta);
  T    res[1];
  T    input[3] = {5, 5, 0};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);
  ASSERT_EQ(1.0, res[0]);
}

TEST(TestCuboid3D, testRotateSmoothInd2)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  Vector<T, 3>                    extend {20, 2, 2};
  Vector<T, 3>                    theta {0, 0, 90};
  T                               epsilon(0.5);
  particles::creators::addResolvedCuboid3D(particleSystem, center, extend,
                                           epsilon, 0., theta);
  T    res[1];
  T    input[3] = {0, 9, 0};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);
  ASSERT_EQ(1.0, res[0]);
}

TEST(TestCuboid3D, testRotateandTranslateSmoothInd2)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {10, 10, 10};
  Vector<T, 3>                    extend {20, 2, 2};
  Vector<T, 3>                    theta {0, 0, 45};
  T                               epsilon(0.5);
  particles::creators::addResolvedCuboid3D(particleSystem, center, extend,
                                           epsilon, 0., theta);
  ;
  T    res[1];
  T    input[3] = {15, 15, 10};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);
  ASSERT_EQ(1.0, res[0]);
}

//Cylinder
TEST(TestSmoothIndicatorCylinder3D, testSmallCylinder)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  Vector<T, 3>                    normal {1, 0, 0};
  T                               length = 2;
  T                               radius = 1;
  Vector<T, 3>                    theta {0, 0, 0};
  T                               epsilon(0.5);
  particles::creators::addResolvedCylinder3D(
      particleSystem, center, normal, length, radius, epsilon, 0., theta);

  // Test of center
  T    res[1];
  T    input[3] = {0, 0, 0};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);
  ASSERT_EQ(1.0, res[0]);

  for (int i = -2; i < 3; ++i) {
    for (int j = -2; j < 3; ++j) {
      for (int k = -2; k < 3; ++k) {
        input[0] = i;
        input[1] = j;
        input[2] = k;
        if (input[0] == 0 && input[1] == 0 && input[2] == 0) {
          continue;
        }
        // Test of outside points
        else if ((util::abs(input[0]) > 1 || util::abs(input[1]) > 1 ||
                  util::abs(input[2]) > 1) ||
                 (util::abs(input[2]) == util::abs(input[1]) &&
                  util::abs(input[1]) > 0)) {
          res[0] = getPorosity(particle, input);
          ASSERT_EQ(0., res[0]);
        }
        // Test of points on smooth boundary
        else {
          res[0] = getPorosity(particle, input);
          ASSERT_NEAR(0.5, res[0], 1e-15);
        }
      }
    }
  }
}

//Cone
TEST(TestCone3D, testMomentOfInertiaAndMass)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center1 {0, 0, 0};
  Vector<T, 3>                    center2 {0, 5, 0};
  T            radius1(2.5), radius2(0), epsilon(0.5), density(1.);
  Vector<T, 3> angle {0, 0, 0};
  particles::creators::addResolvedCone3D(particleSystem, center1, center2,
                                         radius1, radius2, epsilon, density,
                                         angle);
  SmoothIndicatorCone3D<T, T, false> sind(center1, center2, radius1, radius2,
                                          epsilon, angle);

  auto particle = particleSystem.get(0);
  for (unsigned iD = 0; iD < 3; ++iD) {
    const T mofi = particle.template getField<PHYSPROPERTIES, MOFI>()[iD];
    EXPECT_EQ(mofi, sind.calcMofiAndMass(density)[iD]);
  }
  const T mass = particle.template getField<PHYSPROPERTIES, MASS>();
  EXPECT_EQ(mass, sind.calcMofiAndMass(density)[3]);
}

TEST(TestCone3D, testCenter1)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center1 {0, 0, 0};
  Vector<T, 3>                    center2 {0, 5, 0};
  T            radius1(2.5), radius2(0), epsilon(0.5), density(1.);
  Vector<T, 3> angle {0, 0, 0};
  particles::creators::addResolvedCone3D(particleSystem, center1, center2,
                                         radius1, radius2, epsilon, density,
                                         angle);
  T    res[1];
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, center1.data());

  EXPECT_NEAR(0.5, res[0], 1e-15);
}

TEST(TestCone3D, testCenter2)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center1 {0, 0, 0};
  Vector<T, 3>                    center2 {0, 5, 0};
  T            radius1(2.5), radius2(0), epsilon(0.5), density(1.);
  Vector<T, 3> angle {0, 0, 0};
  particles::creators::addResolvedCone3D(particleSystem, center1, center2,
                                         radius1, radius2, epsilon, density,
                                         angle);
  T    res[1];
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, center1.data());

  EXPECT_NEAR(0.5, res[0], 1e-15);
}

TEST(TestCone3D, testConeTanslate)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<double, 3>               center1 {0, 0, 0}, center2 {5, 0, 0};
  double       radius1(2.5), radius2(1), epsilon(.5), density(1.);
  Vector<T, 3> angle {0, 0, 0};
  particles::creators::addResolvedCone3D(particleSystem, center1, center2,
                                         radius1, radius2, epsilon, density,
                                         angle);
  T                 res[1];
  auto              particle = particleSystem.get(0);
  Vector<double, 3> pos {10, 0, 0};
  particle.template setField<GENERAL, POSITION>(pos);
  res[0] = getPorosity(particle, pos.data());

  ASSERT_EQ(1, res[0]);
}

TEST(TestCone3D, testConeTanslate2)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<double, 3>               center1 {0, 0, 0}, center2 {5, 0, 0};
  double       radius1(2.5), radius2(1), epsilon(.5), density(1.);
  Vector<T, 3> angle {0, 0, 0};
  particles::creators::addResolvedCone3D(particleSystem, center1, center2,
                                         radius1, radius2, epsilon, density,
                                         angle);
  T                 res[1];
  auto              particle = particleSystem.get(0);
  Vector<double, 3> pos {-10, 0, 0};
  particle.template setField<GENERAL, POSITION>(pos);
  res[0] = getPorosity(particle, pos.data());

  ASSERT_EQ(1, res[0]);
}

TEST(TestCone3D, testConeTanslate3)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<double, 3>               center1 {50, 0, 0}, center2 {55, 0, 0};
  double       radius1(2.5), radius2(1), epsilon(.5), density(1.);
  Vector<T, 3> angle {0, 0, 0};
  particles::creators::addResolvedCone3D(particleSystem, center1, center2,
                                         radius1, radius2, epsilon, density,
                                         angle);
  T                 res[1];
  auto              particle = particleSystem.get(0);
  Vector<double, 3> pos {10, 0, 0};
  particle.template setField<GENERAL, POSITION>(pos);
  res[0] = getPorosity(particle, pos.data());

  ASSERT_EQ(1, res[0]);
}

TEST(TestCone3D, testConeRotate)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<double, 3>               center1 {-10, 0, 0}, center2 {10, 0, 0};
  double       radius1(2), radius2(1), epsilon(.5), density(1.);
  Vector<T, 3> angle {0, 0, 0};
  particles::creators::addResolvedCone3D(particleSystem, center1, center2,
                                         radius1, radius2, epsilon, density,
                                         angle);
  auto particle = particleSystem.get(0);
  angle[2]      = M_PI / 2;
  particle.template setField<SURFACE, ROT_MATRIX>(
      util::calculateRotationMatrix<T, 3>(angle));

  T input[3] = {-2, 9, 0};
  T res[1];
  res[0] = getPorosity(particle, input);

  ASSERT_EQ(1, res[0]);
}

TEST(TestCone3D, testConeRotate2)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<double, 3>               center1 {-10, 0, 0}, center2 {10, 0, 0};
  double       radius1(2), radius2(1), epsilon(.5), density(1.);
  Vector<T, 3> angle {0, 0, 0};
  particles::creators::addResolvedCone3D(particleSystem, center1, center2,
                                         radius1, radius2, epsilon, density,
                                         angle);
  auto particle = particleSystem.get(0);
  angle[2]      = M_PI / 4;
  particle.template setField<SURFACE, ROT_MATRIX>(
      util::calculateRotationMatrix<T, 3>(angle));

  T input[3] = {2, 3, 0};
  T res[1];
  res[0] = getPorosity(particle, input);

  ASSERT_EQ(1, res[0]);
}

TEST(TestCone3D, testConeRotateAndTranslate)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<double, 3>               center1 {-10, 0, 0}, center2 {10, 0, 0};
  double       radius1(2), radius2(1), epsilon(.5), density(1.);
  Vector<T, 3> angle {0, 0, 0};
  particles::creators::addResolvedCone3D(particleSystem, center1, center2,
                                         radius1, radius2, epsilon, density,
                                         angle);
  auto particle = particleSystem.get(0);
  angle[2]      = M_PI / 4;
  particle.template setField<SURFACE, ROT_MATRIX>(
      util::calculateRotationMatrix<T, 3>(angle));
  Vector<T, 3> pos {10, 10, 10};
  particle.template setField<GENERAL, POSITION>(pos);

  T input[3] = {15, 15, 10};
  T res[1];
  res[0] = getPorosity(particle, input);

  ASSERT_EQ(1, res[0]);
}

//Ellipsoid
TEST(TestEllipsoid3D, testMomentOfInertiaAndMass)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  Vector<T, 3>                    radius(2.5, 5.1, 6.7);
  T                               epsilon(0.1), density(1.);
  Vector<T, 3>                    angle {0, 0, 0};
  particles::creators::addResolvedEllipsoid3D(particleSystem, center, radius,
                                              epsilon, density, angle);
  SmoothIndicatorEllipsoid3D<T, T, false> sind(center, radius, epsilon, angle);

  auto particle = particleSystem.get(0);
  for (unsigned iD = 0; iD < 3; ++iD) {
    const T mofi = particle.template getField<PHYSPROPERTIES, MOFI>()[iD];
    EXPECT_EQ(mofi, sind.calcMofiAndMass(density)[iD]);
  }
  const T mass = particle.template getField<PHYSPROPERTIES, MASS>();
  EXPECT_EQ(mass, sind.calcMofiAndMass(density)[3]);
}

TEST(TestEllipsoid3D, testCenter)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  Vector<T, 3>                    radius(2.5);
  T                               epsilon(0.1), density(1.);
  Vector<T, 3>                    angle {0, 0, 0};
  particles::creators::addResolvedEllipsoid3D(particleSystem, center, radius,
                                              epsilon, density, angle);
  T    res[1];
  T    input[3] = {center[0], center[1], center[2]};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

TEST(TestEllipsoid3D, testInnerBorder)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  Vector<T, 3>                    radius(2.5);
  T                               epsilon(0.1), density(1.);
  Vector<T, 3>                    angle {0, 0, 0};
  particles::creators::addResolvedEllipsoid3D(particleSystem, center, radius,
                                              epsilon, density, angle);
  T res[1];
  T input[3] = {center[0] + (radius[0] - epsilon / 2.), center[1], center[2]};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

TEST(TestEllipsoid3D, testAtRadius)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  Vector<T, 3>                    radius(2.5);
  T                               epsilon(0.1), density(1.);
  Vector<T, 3>                    angle {0, 0, 0};
  particles::creators::addResolvedEllipsoid3D(particleSystem, center, radius,
                                              epsilon, density, angle);
  T    res[1];
  T    input[3] = {center[0], center[1] - radius[1], center[2]};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_LT(0.0, res[0]);
  ASSERT_GT(1.0, res[0]);
}

TEST(TestEllipsoid3D, testOuterBorder)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  Vector<T, 3>                    radius(2.5);
  T                               epsilon(0.1), density(1.);
  Vector<T, 3>                    angle {0, 0, 0};
  particles::creators::addResolvedEllipsoid3D(particleSystem, center, radius,
                                              epsilon, density, angle);
  T    res[1];
  T    input[3] = {center[0], center[1],
                   center[2] - (radius[2] + epsilon / 2. +
                             2 * std::numeric_limits<T>::epsilon())};
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, input);

  ASSERT_EQ(0.0, res[0]);
}

TEST(TestEllipsoid3D, testMovedAndRotated)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {5, 5, 5};
  Vector<T, 3>                    radius {1, 2, 3};
  T                               epsilon(0.1), density(1.);
  Vector<T, 3>                    angle {0, 0, 90};
  particles::creators::addResolvedEllipsoid3D(particleSystem, center, radius,
                                              epsilon, density, angle);
  T    res[1];
  auto particle = particleSystem.get(0);

  T input[3] = {center[0], center[1], center[2]};
  res[0]     = getPorosity(particle, input);
  ASSERT_EQ(1.0, res[0]);

  input[1] = 6;
  res[0]   = getPorosity(particle, input);
  EXPECT_NEAR(0.5, res[0], 5e-16);

  input[1] = 4;
  res[0]   = getPorosity(particle, input);
  EXPECT_NEAR(0.5, res[0], 5e-16);

  input[0] = 3;
  input[1] = 5;
  res[0]   = getPorosity(particle, input);
  EXPECT_NEAR(0.5, res[0], 5e-16);

  input[0] = 7;
  res[0]   = getPorosity(particle, input);
  EXPECT_NEAR(0.5, res[0], 5e-16);

  input[0] = 5;
  input[2] = 2;
  res[0]   = getPorosity(particle, input);
  EXPECT_NEAR(0.5, res[0], 5e-16);

  input[2] = 8;
  res[0]   = getPorosity(particle, input);
  EXPECT_NEAR(0.5, res[0], 5e-16);
}

// Arbitrary shape
TEST(TestArbitraryShape, testMomentOfInertiaAndMass)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  T radius(2.5), epsilon(0.), latticeSpacing(radius / 10), density(1e3);
  std::shared_ptr<IndicatorF3D<T>> ind(
      new IndicatorSphere3D<T>(center, radius));
  SmoothIndicatorCustom3D<T, T, false> sind(latticeSpacing, ind, center,
                                            epsilon);
  particles::creators::addResolvedArbitraryShape3D(
      particleSystem, center, latticeSpacing, ind, epsilon, density);

  auto particle = particleSystem.get(0);
  for (unsigned iD = 0; iD < 3; ++iD) {
    const T mofi = particle.template getField<PHYSPROPERTIES, MOFI>()[iD];
    EXPECT_EQ(mofi, sind.calcMofiAndMass(density)[iD]);
  }
  const T mass = particle.template getField<PHYSPROPERTIES, MASS>();
  EXPECT_EQ(mass, sind.calcMofiAndMass(density)[3]);
}

TEST(TestArbitraryShape, testSphereInOrigin)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {0, 0, 0};
  T radius(2.5), epsilon(0.), latticeSpacing(radius / 10), density(1e3);
  std::shared_ptr<IndicatorF3D<T>> ind(
      new IndicatorSphere3D<T>(center, radius));
  particles::creators::addResolvedArbitraryShape3D(
      particleSystem, center, latticeSpacing, ind, epsilon, density);
  T    res[1];
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, center.data());
  EXPECT_EQ(1.0, res[0]);
}

TEST(TestArbitraryShape, testSphereMovedFromOrigin)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {10, 10, 10};
  T radius(2.5), epsilon(0.), latticeSpacing(radius / 10), density(1e3);
  std::shared_ptr<IndicatorF3D<T>> ind(
      new IndicatorSphere3D<T>(center, radius));
  particles::creators::addResolvedArbitraryShape3D(
      particleSystem, center, latticeSpacing, ind, epsilon, density);
  T    res[1];
  auto particle = particleSystem.get(0);
  res[0]        = getPorosity(particle, center.data());
  EXPECT_EQ(1.0, res[0]);
}

TEST(TestArbitraryShape, testCuboidMovedAndRotated1)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {5, 5, 5};
  Vector<T, 3>                    extend {6, 2, 2};
  Vector<T, 3>                    theta(0, 0, 90);
  T epsilon(0.), latticeSpacing(extend[1] / 32), density(1e3);
  std::shared_ptr<IndicatorF3D<T>> ind(
      new IndicatorCuboid3D<T>(extend, center));
  particles::creators::addResolvedArbitraryShape3D(
      particleSystem, center, latticeSpacing, ind, epsilon, density, theta);
  T    res[1];
  auto particle = particleSystem.get(0);

  Vector<T, 3> input = center + 0.5 * extend;
  input[2] = center[2];
  res[0]             = getPorosity(particle, input.data());
  EXPECT_EQ(0.0, res[0]);

  input[0] = center[0] - 0.5 * extend[0] + 0.5 * latticeSpacing;
  res[0]   = getPorosity(particle, input.data());
  EXPECT_EQ(0.0, res[0]);

  input[0] = center[0];
  input[1] = center[1] + 0.5 * extend[0] - 0.5 * latticeSpacing;
  res[0]   = getPorosity(particle, input.data());
  EXPECT_EQ(1.0, res[0]);
  input[1] += latticeSpacing;
  res[0] = getPorosity(particle, input.data());
  EXPECT_EQ(0.0, res[0]);

  input[1] = center[1] - 0.5 * extend[0] + 0.5 * latticeSpacing;
  res[0]   = getPorosity(particle, input.data());
  EXPECT_EQ(1.0, res[0]);
  input[1] -= latticeSpacing;
  res[0] = getPorosity(particle, input.data());
  EXPECT_EQ(0.0, res[0]);

  input[1] = center[0];
  input[2] = center[2] + 0.5 * extend[2] - 0.5 * latticeSpacing;
  res[0]   = getPorosity(particle, input.data());
  EXPECT_EQ(1.0, res[0]);
  input[2] += latticeSpacing;
  res[0] = getPorosity(particle, input.data());
  EXPECT_EQ(0.0, res[0]);

  input[2] = center[2] - 0.5 * extend[2] + 0.5 * latticeSpacing;
  res[0]   = getPorosity(particle, input.data());
  EXPECT_EQ(1.0, res[0]);
  input[2] -= latticeSpacing;
  res[0] = getPorosity(particle, input.data());
  EXPECT_EQ(0.0, res[0]);
}

TEST(TestArbitraryShape, testCuboidMovedAndRotated2)
{
  ParticleSystem<T, PARTICLETYPE> particleSystem;
  Vector<T, 3>                    center {5, 5, 5};
  Vector<T, 3>                    extend {0.7, 0.25, 0.25};
  Vector<T, 3>                    theta(0, 0, 90);
  T epsilon(0.), latticeSpacing(extend[1] / 32), density(1e3);
  std::shared_ptr<IndicatorF3D<T>> ind(
      new IndicatorCuboid3D<T>(extend, center));
  particles::creators::addResolvedArbitraryShape3D(
      particleSystem, center, latticeSpacing, ind, epsilon, density, theta);
  T    res[1];
  auto particle = particleSystem.get(0);

  Vector<T, 3> input = center + 0.5 * extend;
  input[2] = center[2];
  res[0]             = getPorosity(particle, input.data());
  EXPECT_EQ(0.0, res[0]);

  input[0] = center[0] - 0.5 * extend[0];
  res[0]   = getPorosity(particle, input.data());
  EXPECT_EQ(0.0, res[0]);

  input[0] = center[0];
  input[1] = center[1] + 0.5 * extend[0] - 0.5 * latticeSpacing;
  res[0]   = getPorosity(particle, input.data());
  EXPECT_EQ(1.0, res[0]);
  input[1] += latticeSpacing;
  res[0] = getPorosity(particle, input.data());
  EXPECT_EQ(0.0, res[0]);

  input[1] = center[1] - 0.5 * extend[0] + 0.5 * latticeSpacing;
  res[0]   = getPorosity(particle, input.data());
  EXPECT_EQ(1.0, res[0]);
  input[1] -= latticeSpacing;
  res[0] = getPorosity(particle, input.data());
  EXPECT_EQ(0.0, res[0]);

  input[1] = center[1];
  input[2] = center[2] + 0.5 * extend[2] - 0.5 * latticeSpacing;
  res[0]   = getPorosity(particle, input.data());
  EXPECT_EQ(1.0, res[0]);
  input[2] += latticeSpacing;
  res[0] = getPorosity(particle, input.data());
  EXPECT_EQ(0.0, res[0]);

  input[2] = center[2] - 0.5 * extend[2] + 0.5 * latticeSpacing;
  res[0]   = getPorosity(particle, input.data());
  EXPECT_EQ(1.0, res[0]);
  input[2] -= latticeSpacing;
  res[0] = getPorosity(particle, input.data());
  EXPECT_EQ(0.0, res[0]);
}
