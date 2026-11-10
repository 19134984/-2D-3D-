/*  This file is part of the OpenLB library
*
*  Copyright (C) 2017 Markus Mohrhard, 2025 Christoph Gaul
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
#include <string>

/// ensure olbInit is called for MPI
::testing::Environment* const foo_env = ::testing::AddGlobalTestEnvironment(new OlbInitTestEnvironment);

using namespace olb;
using namespace olb::descriptors;

std::unique_ptr<IndicatorF3D<double>> createIndicator(const std::string& name, double size)
{
  if (name == "Sphere") {
    return std::make_unique<IndicatorSphere3D<double>>(Vector<double, 3> {0, 0, 0}, size);
  }
  else if (name == "Cuboid") {
    return std::make_unique<IndicatorCuboid3D<double>>(Vector<double, 3> {size, size, size},
                                                       Vector<double, 3> {0, 0, 0});
  }
  else if (name == "Cylinder") {
    return std::make_unique<IndicatorCylinder3D<double>>(Vector<double, 3> {0, 0, 0}, Vector<double, 3> {1, 0, 0},
                                                         size);
  }
  else if (name == "Ellipsoid") {
    return std::make_unique<IndicatorEllipsoid3D<double>>(Vector<double, 3> {0, 0, 0}, Vector<double, 3> {1, 1, 1});
  }
  else if (name == "SuperEllipsoid") {
    return std::make_unique<IndicatorSuperEllipsoid3D<double>>(Vector<double, 3> {0, 0, 0}, size, size, size, size,
                                                               size);
  }
  else if (name == "Polygon") {
    return std::make_unique<IndicatorPolygon3D<double>>(std::vector<Vector<double, 4>> {
        {0, 0, 0, 0},
        {1, 0, 0, 0},
        {1, 1, 0, 0},
        {0, 1, 0, 0}
    });
  }
  else if (name == "Cone") {
    return std::make_unique<IndicatorCone3D<double>>(Vector<double, 3> {0, 0, 0}, Vector<double, 3> {1, 0, 0}, size,
                                                     size);
  }
  return nullptr;
}

void indicatorDecompositionTests(CuboidDecomposition3D<double>& cuboid,
                                 IndicatorF3D<double>& indicator, double deltaX)
{
  EXPECT_DOUBLE_EQ(indicator.getMin()[0], cuboid.getMinPhysR()[0]);
  EXPECT_DOUBLE_EQ(indicator.getMin()[1], cuboid.getMinPhysR()[1]);
  EXPECT_DOUBLE_EQ(indicator.getMin()[2], cuboid.getMinPhysR()[2]);
  EXPECT_DOUBLE_EQ(indicator.getMax()[0] + deltaX, cuboid.getMaxPhysR()[0]);
  EXPECT_DOUBLE_EQ(indicator.getMax()[1] + deltaX, cuboid.getMaxPhysR()[1]);
  EXPECT_DOUBLE_EQ(indicator.getMax()[2] + deltaX, cuboid.getMaxPhysR()[2]);
}

TEST(CuboidDecomposition3DTest, testShrinkIndicatorLoop)
{
  std::vector<std::string> strategies     = {"standard", "weight", "volume"};
  std::vector<std::string> indicatorNames = {"Sphere",  "Cuboid",    "Cylinder", "Ellipsoid"};
  double                   size           = 1.0;

  std::vector<int> nCuboidsList    = {2, 10, 25, 50, 100};
  std::vector<int> nResolutionList = {10, 20, 40, 80, 160};
  double           deltaX;
  for (std::string& indicatorName : indicatorNames) {
    auto indicator = createIndicator(indicatorName, size);
    for (std::string& strategy : strategies) {
      std::cout << std::endl;
      std::cout << "Testing indicator: " << indicatorName << " with decomposition strategy: " << strategy << std::endl;
      std::cout << std::endl;
      for (std::size_t jCuboids : nCuboidsList) {
        for (std::size_t iResolution : nResolutionList) {
          deltaX = 1.0 / double(iResolution);
          if (strategy == "standard") {
            CuboidDecomposition3D<double> cuboid(*indicator, deltaX, jCuboids);
            indicatorDecompositionTests(cuboid, *indicator, deltaX);
          }
          else {
            CuboidDecomposition3D<double> cuboid(*indicator, deltaX, jCuboids, strategy);
            indicatorDecompositionTests(cuboid, *indicator, deltaX);
          }
        }
      }
    }
  }
}

TEST(CuboidDecomposition3DTest, testShrinkBifurcationGeometry)
{
  STLreader<double>             stlReader("../../examples/particles/bifurcation3d/bifurcation3d.stl", 0.1, 100.0);
  CuboidDecomposition3D<double> cuboid(stlReader, 0.1, 100);

  EXPECT_EQ(std::size_t(70), cuboid.size());
  EXPECT_NEAR(double(-2.974935), cuboid.getMinPhysR()[0], 1e-6);
  EXPECT_NEAR(double(3.025065), cuboid.getMaxPhysR()[0], 1e-6);
  EXPECT_NEAR(double(0.075), cuboid.getMinPhysVolume(), 1e-6);
  EXPECT_NEAR(double(1.2), cuboid.getMaxPhysVolume(), 1e-6);
}

TEST(CuboidDecomposition3DTest, testShrinkAortaGeometry)
{
  STLreader<double>             stlReader("../../examples/turbulence/aorta3d/aorta3d.stl", 0.1, 0.1);
  CuboidDecomposition3D<double> cuboid(stlReader, 0.1, 100);

  EXPECT_EQ(std::size_t(62), cuboid.size());
  EXPECT_NEAR(double(20.09925), cuboid.getMinPhysR()[0], 1e-6);
  EXPECT_NEAR(double(24.39925), cuboid.getMaxPhysR()[0], 1e-6);
  EXPECT_NEAR(double(0.004), cuboid.getMinPhysVolume(), 1e-6);
  EXPECT_NEAR(double(1.848), cuboid.getMaxPhysVolume(), 1e-6);
}

TEST(CuboidDecomposition3DTest, testShrink2Dloop)
{
  Vector<double, 2>         extend(1, 1);
  Vector<double, 2>         origin(0, 0);
  IndicatorCuboid2D<double> indicator(extend, origin);
  double                    deltaX;
  int                       numCuboids = 2;
  for (int i = 5; i < 200; i += 1) {
    deltaX = 1.0 / double(i);
    CuboidDecomposition2D<double> cuboid(indicator, deltaX, numCuboids);

    EXPECT_NEAR(double(0.0), cuboid.getMinPhysR()[0], 1e-6);
    EXPECT_NEAR(double(0.0), cuboid.getMinPhysR()[1], 1e-6);
    EXPECT_NEAR(double(extend[0]) + deltaX, cuboid.getMaxPhysR()[0], 1e-6);
    EXPECT_NEAR(double(extend[1]) + deltaX, cuboid.getMaxPhysR()[1], 1e-6);
  }
}
