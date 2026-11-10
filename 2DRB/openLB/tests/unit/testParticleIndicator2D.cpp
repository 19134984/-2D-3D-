#include <iostream>
#include <gtest/gtest.h>
#include "helper/unitTestHelper.h"

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::particles;

typedef ResolvedParticle2D PARTICLETYPE;
typedef double T;

T getPorosity(Particle<T,PARTICLETYPE>& particle, T input[])
{
  T tmp[1];
  particles::resolved::evalSolidVolumeFraction(tmp, input, particle);
  return tmp[0];
}

//Sphere
TEST(TestCircle2D, testMomentOfInertiaAndMass)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{0, 0};
  T radius(2.5), epsilon(0.1), density(1.);
  particles::creators::addResolvedCircle2D( particleSystem, center, radius, epsilon, density );
  SmoothIndicatorCircle2D<T,T,false> sind(center, radius, epsilon);

  auto particle = particleSystem.get(0);
  const T mofi = particle.template getField<PHYSPROPERTIES,MOFI>( );
  EXPECT_EQ(mofi, sind.calcMofiAndMass(density)[0]);
  const T mass = particle.template getField<PHYSPROPERTIES,MASS>( );
  EXPECT_EQ(mass, sind.calcMofiAndMass(density)[1]);
}

TEST(TestCircle2D, testCenter)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{0, 0};
  T radius(2.5), epsilon(0.1);
  particles::creators::addResolvedCircle2D( particleSystem, center, radius, epsilon );
  T res[1];
  T input[2] = {center[0], center[1]};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

TEST(TestCircle2D, testInnerBorder)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{0, 0};
  T radius(2.5), epsilon(0.1);
  particles::creators::addResolvedCircle2D( particleSystem, center, radius, epsilon );
  T res[1];
  T input[2] = {center[0]+(radius-epsilon/2.), center[1]};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

TEST(TestCircle2D, testAtRadius)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{0, 0};
  T radius(2.5), epsilon(0.1);
  particles::creators::addResolvedCircle2D( particleSystem, center, radius, epsilon );
  T res[1];
  T input[2] = {center[0], center[1]-radius};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_LT(0.0, res[0]);
  ASSERT_GT(1.0, res[0]);
}

TEST(TestCircle2D, testOuterBorder)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{0, 0};
  T radius(2.5), epsilon(0.1);
  particles::creators::addResolvedCircle2D( particleSystem, center, radius, epsilon );
  T res[1];
  T input[2] = {center[0], center[1]-(radius+epsilon/2.)};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_TRUE(util::nearZero(res[0]));
}

TEST(TestCircle2D, testFarOut)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{0, 0};
  T radius(2.5), epsilon(0.1);
  particles::creators::addResolvedCircle2D( particleSystem, center, radius, epsilon );
  T res[1];
  T input[2] = {center[0]+1000*epsilon, center[1]-500*epsilon};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_EQ(0.0, res[0]);
}

TEST(TestCircle2D, testTranslatedCenter)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{20, 4.3};
  T radius(2.5), epsilon(0.1);
  particles::creators::addResolvedCircle2D( particleSystem, center, radius, epsilon );
  T res[1];
  T input[2] = {center[0], center[1]};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

TEST(TestCircle2D, testTranslatedInnerBorder)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{0.3, -5};
  T radius(2.5), epsilon(0.1);
  particles::creators::addResolvedCircle2D( particleSystem, center, radius, epsilon );
  T res[1];
  T input[2] = {center[0]+(radius-epsilon/2.), center[1]};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

TEST(TestCircle2D, testTranslatedAtRadius)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{-0.9, 5};
  T radius(2.5), epsilon(0.1);
  particles::creators::addResolvedCircle2D( particleSystem, center, radius, epsilon );
  T res[1];
  T input[2] = {center[0], center[1]-radius};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_LT(0.0, res[0]);
  ASSERT_GT(1.0, res[0]);
}

TEST(TestCircle2D, testTranslatedOuterBorder)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{-48.5, 22.5};
  T radius(2.5), epsilon(0.1);
  particles::creators::addResolvedCircle2D( particleSystem, center, radius, epsilon );
  T res[1];
  T input[2] = {center[0], center[1]-(radius+epsilon/2.)};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);
  ASSERT_NEAR(0.0, res[0], 1e-26);
}

TEST(TestCircle2D, testTranslatedFarOut)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{30.2, -400};
  T radius(2.5), epsilon(0.1);
  particles::creators::addResolvedCircle2D( particleSystem, center, radius, epsilon );
  T res[1];
  T input[2] = {center[0]+1000*epsilon, center[1]-500*epsilon};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_EQ(0.0, res[0]);
}

TEST(TestCircle2D, testOuterBorderLargeRadius)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{0, 0};
  T radius(100), epsilon(0.1);
  particles::creators::addResolvedCircle2D( particleSystem, center, radius, epsilon );
  T res[1];
  T input[2] = {center[0],  center[1]-(radius+epsilon/2.)};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_TRUE(util::nearZero(res[0]));
}

TEST(TestCircle2D, testCenterLargeEpsilon)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{0, 0};
  T radius(2.5), epsilon(4);
  particles::creators::addResolvedCircle2D( particleSystem, center, radius, epsilon );
  T res[1];
  T input[2] = {center[0], center[1]};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

//Cuboid
TEST(TestCuboid2D, testMomentOfInertiaAndMass)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{.2, .1};
  Vector<T, 2> extend{.2, 0.4};
  T theta(0), epsilon(.05), density(0.5);
  particles::creators::addResolvedCuboid2D( particleSystem, center, extend, epsilon, density, theta );
  SmoothIndicatorCuboid2D<T,T,false> sind(center, extend[0], extend[1], epsilon, theta);

  auto particle = particleSystem.get(0);
  const T mofi = particle.template getField<PHYSPROPERTIES,MOFI>( );
  EXPECT_EQ(mofi, sind.calcMofiAndMass(density)[0]);
  const T mass = particle.template getField<PHYSPROPERTIES,MASS>( );
  EXPECT_EQ(mass, sind.calcMofiAndMass(density)[1]);
}

TEST(TestCuboid2D, testSmallCuboidCenter)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{.2, .1};
  Vector<T, 2> extend{.2, 0.4};
  T theta(0), epsilon(.05);
  particles::creators::addResolvedCuboid2D( particleSystem, center, extend, epsilon, 0., theta );
  T res[1];
  T input[2] = {center[0]+(extend[0]-epsilon)/2., center[1]+(extend[1]-epsilon)/2.};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

TEST(TestCuboid2D, testSmallCuboidUpperEdge)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{2, 11};
  Vector<T, 2> extend{4, 8};
  T theta(0), epsilon(.05);
  particles::creators::addResolvedCuboid2D( particleSystem, center, extend, epsilon, 0., theta );
  T res[1];
  T input[2] = {center[0], center[1]+(extend[1]-epsilon)/2.};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

TEST(TestCuboid2D, testSmallCuboidFarOut)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{2, 11};
  Vector<T, 2> extend{4, 8};
  T theta(0), epsilon(.05);
  particles::creators::addResolvedCuboid2D( particleSystem, center, extend, epsilon, 0., theta );
  T res[1];
  T input[2] = {center[0]+400, center[1]+16};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_EQ(0.0, res[0]);
}

TEST(TestCuboid2D, testLargeCuboidLeftEdge)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{2, 11};
  Vector<T, 2> extend{40, 800};
  T theta(0), epsilon(.05);
  particles::creators::addResolvedCuboid2D( particleSystem, center, extend, epsilon, 0., theta );
  T res[1];
  T input[2] = {center[0]-(extend[0]-epsilon)/2., center[1]};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

TEST(TestCuboid2D, testCuboidEdge)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{2, 11};
  Vector<T, 2> extend{4, 8};
  T theta(0), epsilon(.05);
  particles::creators::addResolvedCuboid2D( particleSystem, center, extend, epsilon, 0., theta );
  T res[1];
  T input[2] = {center[0]+(extend[0]-epsilon)/2., center[1]+(extend[1]-epsilon)/2.};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

TEST(TestCuboid2D, testCuboidEdgeHalfEpsilonOut)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{.2, -11};
  Vector<T, 2> extend{4, 8};
  T theta(0), epsilon(.05);
  particles::creators::addResolvedCuboid2D( particleSystem, center, extend, epsilon, 0., theta );
  T res[1];
  T input[2] = {center[0]+extend[0]/2., center[1]+extend[1]/2.};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_LT(0.0, res[0]);
  ASSERT_GT(1.0, res[0]);
}

TEST(TestCuboid2D, testCuboidEpsilonAboveEdge)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{2, 11};
  Vector<T, 2> extend{4, 8};
  T theta(0), epsilon(.05);
  particles::creators::addResolvedCuboid2D( particleSystem, center, extend, epsilon, 0., theta );
  T res[1];
  T input[2] = {center[0], center[1]+extend[1]/2.+epsilon};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_NEAR(0.0, res[0], 1e-15);
}

TEST(TestCuboid2D, testCuboidEpsilonOutCorner)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{2, 11};
  Vector<T, 2> extend{4, 8};
  T theta(0), epsilon(.05);
  particles::creators::addResolvedCuboid2D( particleSystem, center, extend, epsilon, 0., theta );
  T res[1];
  T input[2] = {center[0]-extend[0]/2.-epsilon, center[1]-extend[1]/2.-epsilon};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_NEAR(0.0, res[0], 1e-15);
}

TEST(TestCuboid2D, testCuboidLowerLeftCornerRotated)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{2, 11};
  Vector<T, 2> extend{4, 8};
  T theta(45), epsilon(.05);
  particles::creators::addResolvedCuboid2D( particleSystem, center, extend, epsilon, 0., theta );
  T res[1];
  T input[2] = {center[0]+util::sqrt(2), center[1]-3*util::sqrt(2)};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_NEAR(0.5, res[0], 1e-12);
}

//Triangle
TEST(TestTriangle2D, testMomentOfInertiaAndMass)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{12, -4};
  T radius(4), theta(180), epsilon(0.1), density(1e3);
  particles::creators::addResolvedTriangle2D( particleSystem, center, radius, epsilon, density, theta );
  SmoothIndicatorTriangle2D<T,T,false> sind(center, radius, epsilon);

  auto particle = particleSystem.get(0);
  const T mofi = particle.template getField<PHYSPROPERTIES,MOFI>( );
  EXPECT_EQ(mofi, sind.calcMofiAndMass(density)[0]);
  const T mass = particle.template getField<PHYSPROPERTIES,MASS>( );
  EXPECT_EQ(mass, sind.calcMofiAndMass(density)[1]);
}

TEST(TestTriangle2D, testTriangleCenter)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{12, -4};
  T radius(4), theta(180), epsilon(0.1);
  particles::creators::addResolvedTriangle2D( particleSystem, center, radius, epsilon, 0., theta );
  T res[1];
  T input[2] = {center[0], center[1]};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

TEST(TestTriangle2D, testTriangleRotatedCorner)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{12, -4};
  T radius(4), theta(180), epsilon(0.1);
  particles::creators::addResolvedTriangle2D( particleSystem, center, radius, epsilon, 0., theta );
  T res[1];
  T input[2] = {center[0], center[1] - radius};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_NEAR(0.5, res[0], 1e-14);
}


TEST(TestTriangle2D, testTriangleFarOut)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{12, -4};
  T radius(4), theta(180), epsilon(0.1);
  particles::creators::addResolvedTriangle2D( particleSystem, center, radius, epsilon, 0., theta );
  T res[1];
  T input[2] = {center[0]+4000, center[1]-200};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_EQ(0.0, res[0]);
}
TEST(TestTriangle2D, testTriangleEpsilonArea)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{12, -4};
  T radius(4), theta(180), epsilon(0.1);
  particles::creators::addResolvedTriangle2D( particleSystem, center, radius, epsilon, 0., theta );
  T res[1];
  T input[2] = {center[0], center[1]+(radius+epsilon)/2.};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_LT(0.0, res[0]);
  ASSERT_GT(1.0, res[0]);
}

TEST(TestTriangle2D, testTriangleEdge)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{12, -4};
  T radius(4), theta(180), epsilon(0.1);
  particles::creators::addResolvedTriangle2D( particleSystem, center, radius, epsilon, 0., theta );
  T res[1];
  T input[2] = {center[0], center[1] - (radius-epsilon)};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_EQ(1.0, res[0]);
}

TEST(TestTriangle2D, testTriangleLargeRadiusFarOut)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{12, -4};
  T radius(400), theta(180), epsilon(0.1);
  particles::creators::addResolvedTriangle2D( particleSystem, center, radius, epsilon, 0., theta );
  T res[1];
  T input[2] = {center[0]+800, center[1] - 300};
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, input);

  ASSERT_EQ(0.0, res[0]);
}


// Arbitrary shape
TEST(TestArbitraryShape, testMomentOfInertiaAndMass)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{0, 0};
  T radius(2.5), theta(0.), epsilon(0.), latticeSpacing(1e-2), density(1e3);
  std::shared_ptr<IndicatorF2D<T>> ind (new IndicatorCircle2D<T>(center, radius));
  SmoothIndicatorCustom2D<T, T, false> sind(latticeSpacing, ind, center, epsilon, theta);
  particles::creators::addResolvedArbitraryShape2D( particleSystem, center, latticeSpacing, ind, epsilon, density, theta );

  auto particle = particleSystem.get(0);
  T mofi = particle.template getField<PHYSPROPERTIES,MOFI>( );
  EXPECT_EQ(mofi, sind.calcMofiAndMass(density)[0]);
  const T mass = particle.template getField<PHYSPROPERTIES,MASS>( );
  EXPECT_EQ(mass, sind.calcMofiAndMass(density)[1]);
}

TEST(TestArbitraryShape, testCircleInOrigin)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{0, 0};
  T radius(2.5), theta(0.), epsilon(0.), latticeSpacing(1e-2), density(1e3);
  std::shared_ptr<IndicatorF2D<T>> ind (new IndicatorCircle2D<T>(center, radius));
  particles::creators::addResolvedArbitraryShape2D( particleSystem, center, latticeSpacing, ind, epsilon, density, theta );
  T res[1];
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, center.data());
  EXPECT_EQ(1.0, res[0]);
}

TEST(TestArbitraryShape, testCircleMovedFromOrigin)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{10, 10};
  T radius(2.5), theta(0.), epsilon(0.), latticeSpacing(1e-2), density(1e3);
  std::shared_ptr<IndicatorF2D<double>> ind (new IndicatorCircle2D<double>(center, radius));
  particles::creators::addResolvedArbitraryShape2D( particleSystem, center, latticeSpacing, ind, epsilon, density, theta );
  T res[1];
  auto particle = particleSystem.get(0);
  res[0] = getPorosity(particle, center.data());
  EXPECT_EQ(1.0, res[0]);
}

TEST(TestArbitraryShape, testCuboidMovedAndRotated)
{
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  Vector<T, 2> center{5, 5};
  Vector<T, 2> extend{6, 2};
  T theta(90), epsilon(0.), latticeSpacing(1e-2), density(1e3);
  std::shared_ptr<IndicatorF2D<T>> ind (new IndicatorCuboid2D<T>(extend, center, 0.));
  particles::creators::addResolvedArbitraryShape2D( particleSystem, center, latticeSpacing, ind, epsilon, density, theta );
  T res[1];
  auto particle = particleSystem.get(0);

  T input[2] = {8,5.};
  res[0] = getPorosity(particle, input);
  EXPECT_EQ(0.0, res[0]);

  input[0] = 2;
  res[0] = getPorosity(particle, input);
  EXPECT_EQ(0.0, res[0]);

  input[0] = 5;
  input[1] = 8.-0.5*latticeSpacing;
  res[0] = getPorosity(particle, input);
  EXPECT_EQ(1.0, res[0]);

  input[1] = 2+0.5*latticeSpacing;
  res[0] = getPorosity(particle, input);
  EXPECT_EQ(1.0, res[0]);
}
