/*  This file is part of the OpenLB library
*
*  Copyright (C) 2025 Christoph Gaul, Mathias J. Krause
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

#include <olb.h>
#include <olb.h>

#include "../unit/helper/framework.h"
#include <iostream>
#include <gtest/gtest.h>
#include "../unit/helper/unitTestHelper.h"

/// ensure olbInit is called for MPI
::testing::Environment* const foo_env = ::testing::AddGlobalTestEnvironment(new OlbInitTestEnvironment);

using namespace olb;
using namespace olb::descriptors;

TEST(FloatingPointRelationalOperatorsTest, testEpsilonMin)
{
  EXPECT_TRUE(util::numericLimits::epsilon<float>() > 0);
  EXPECT_TRUE(util::numericLimits::epsilon<double>() > 0);
  EXPECT_TRUE(util::numericLimits::minLimit<float>() > 0);
  EXPECT_TRUE(util::numericLimits::minLimit<double>() > 0);
}

TEST(FloatingPointRelationalOperatorsTest, testEqual)
{
  double a       = 1.0;
  double epsilon = std::numeric_limits<double>::epsilon();
  double b       = a - 100 * epsilon;

  for (int i = -100; i < 100; ++i) {
    b += epsilon;
    if (i <= -7) {
      EXPECT_FALSE(util::equal(b, a));
      EXPECT_TRUE(util::less(b, a));
      EXPECT_TRUE(util::lessOrEqual(b, a));
      EXPECT_FALSE(util::greaterOrEqual(b, a));
      EXPECT_FALSE(util::greater(b, a));
    }
    if (i > -7 && i < 6) {
      EXPECT_TRUE(util::equal(a, b));
      EXPECT_TRUE(util::lessOrEqual(b, a));
      EXPECT_TRUE(util::greaterOrEqual(b, a));
      EXPECT_FALSE(util::greater(b, a));
      EXPECT_FALSE(util::less(b, a));
    }
    else if (i >= 6) {
      EXPECT_FALSE(util::equal(b, a));
      EXPECT_TRUE(util::greater(b, a));
      EXPECT_TRUE(util::greaterOrEqual(b, a));
      EXPECT_FALSE(util::lessOrEqual(b, a));
      EXPECT_FALSE(util::less(b, a));
    }
  }
}

TEST(FloatingPointRelationalOperatorsTest, testOperatorsFloat)
{
  float a       = 1.0;
  float epsilon = std::numeric_limits<float>::epsilon();
  float b       = a - 100 * epsilon;

  for (int i = -100; i < 100; ++i) {
    b += epsilon;
    if (i <= -7) {
      EXPECT_FALSE(util::equal(b, a));
      EXPECT_TRUE(util::less(b, a));
      EXPECT_TRUE(util::lessOrEqual(b, a));
      EXPECT_FALSE(util::greaterOrEqual(b, a));
      EXPECT_FALSE(util::greater(b, a));
    }
    if (i > -7 && i < 6) {
      EXPECT_TRUE(util::equal(a, b));
      EXPECT_TRUE(util::lessOrEqual(b, a));
      EXPECT_TRUE(util::greaterOrEqual(b, a));
      EXPECT_FALSE(util::greater(b, a));
      EXPECT_FALSE(util::less(b, a));
    }
    else if (i >= 6) {
      EXPECT_FALSE(util::equal(b, a));
      EXPECT_TRUE(util::greater(b, a));
      EXPECT_TRUE(util::greaterOrEqual(b, a));
      EXPECT_FALSE(util::lessOrEqual(b, a));
      EXPECT_FALSE(util::less(b, a));
    }
  }
}

TEST(FloatingPointRelationalOperatorsTest, testInsideOutsideFunctions)
{
  double a       = 1.;
  double epsilon = std::numeric_limits<double>::epsilon() * a;
  double b       = a - 100 * epsilon;
  double c       = a + 100 * epsilon;
  for (int i = -100; i < 100; i++) {
    b += epsilon;
    c -= epsilon;
    if (i < -7) {
      EXPECT_TRUE(util::isInside(a, b, c));
      EXPECT_FALSE(util::isOutside(a, b, c));
    }
    else if (i >= -6 && i < 5) {
      EXPECT_FALSE(util::isInside(a, b, c));
      EXPECT_FALSE(util::isOutside(a, b, c));
    }
    else if (i >= 6) {
      EXPECT_FALSE(util::isInside(a, b, c));
      EXPECT_TRUE(util::isOutside(a, b, c));
    }
  }
}

TEST(FloatingPointRelationalOperatorsTest, testInsideOutsideFunctionsFloat)
{
  float a       = 1.;
  float epsilon = std::numeric_limits<float>::epsilon() * a;
  float b       = a - 100 * epsilon;
  float c       = a + 100 * epsilon;
  for (int i = -100; i < 100; i++) {
    b += epsilon;
    c -= epsilon;
    if (i < -7) {
      EXPECT_TRUE(util::isInside(a, b, c));
      EXPECT_FALSE(util::isOutside(a, b, c));
    }
    else if (i >= -6 && i < 5) {
      EXPECT_FALSE(util::isInside(a, b, c));
      EXPECT_FALSE(util::isOutside(a, b, c));
    }
    else if (i >= 6) {
      EXPECT_FALSE(util::isInside(a, b, c));
      EXPECT_TRUE(util::isOutside(a, b, c));
    }
  }
}

TEST(FloatingPointRelationalOperatorsTest, testCuboidSizes)
{
  std::vector<std::string> strategies     = {"standard"};
  std::vector<std::string> indicatorNames = {"Cuboid"};
  double                   size           = 1.0;

  std::vector<int> nCuboidsList    = {2, 10, 25, 50, 100};
  std::vector<int> nResolutionList = {10, 20, 40, 80};
  double           deltaX;
  int              sizeExponentMax = 15;
  int              sizeExponentMin = -15;
  for (int sizeExponent = 0; sizeExponent < sizeExponentMax; ++sizeExponent) {
    size *= 10.0;
    IndicatorCuboid3D<double> indicator(Vector<double, 3> {size, size, size}, Vector<double, 3> {0, 0, 0});
    for (std::size_t jCuboids : nCuboidsList) {
      for (std::size_t iResolution : nResolutionList) {
        deltaX = size / double(iResolution);
        //std::cout << "Testing size" << size << " with " << jCuboids
        //          << " cuboids and resolution " << iResolution << std::endl;
        CuboidDecomposition3D<double> cuboid(indicator, deltaX, jCuboids);
        EXPECT_DOUBLE_EQ(indicator.getMin()[0], cuboid.getMinPhysR()[0]);
        EXPECT_DOUBLE_EQ(indicator.getMin()[1], cuboid.getMinPhysR()[1]);
        EXPECT_DOUBLE_EQ(indicator.getMin()[2], cuboid.getMinPhysR()[2]);
        EXPECT_DOUBLE_EQ(indicator.getMax()[0] + deltaX, cuboid.getMaxPhysR()[0]);
        EXPECT_DOUBLE_EQ(indicator.getMax()[1] + deltaX, cuboid.getMaxPhysR()[1]);
        EXPECT_DOUBLE_EQ(indicator.getMax()[2] + deltaX, cuboid.getMaxPhysR()[2]);
      }
    }
  }

  size = 1;
  for (int sizeExponent = 0; sizeExponent < -sizeExponentMin; ++sizeExponent) {
    size *= 0.1;
    IndicatorCuboid3D<double> indicator(Vector<double, 3> {size, size, size}, Vector<double, 3> {0, 0, 0});
    for (std::size_t jCuboids : nCuboidsList) {
      for (std::size_t iResolution : nResolutionList) {
        deltaX = size / double(iResolution);
        //std::cout << "Testing size" << size << " with " << jCuboids
        //          << " cuboids and resolution " << iResolution << std::endl;
        CuboidDecomposition3D<double> cuboid(indicator, deltaX, jCuboids);
        EXPECT_DOUBLE_EQ(indicator.getMin()[0], cuboid.getMinPhysR()[0]);
        EXPECT_DOUBLE_EQ(indicator.getMin()[1], cuboid.getMinPhysR()[1]);
        EXPECT_DOUBLE_EQ(indicator.getMin()[2], cuboid.getMinPhysR()[2]);
        EXPECT_DOUBLE_EQ(indicator.getMax()[0] + deltaX, cuboid.getMaxPhysR()[0]);
        EXPECT_DOUBLE_EQ(indicator.getMax()[1] + deltaX, cuboid.getMaxPhysR()[1]);
        EXPECT_DOUBLE_EQ(indicator.getMax()[2] + deltaX, cuboid.getMaxPhysR()[2]);
      }
    }
  }
}
