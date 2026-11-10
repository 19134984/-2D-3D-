/*  This file is part of the OpenLB library
*
*  Copyright (C) 2018 Markus Mohrhard
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

#include <iostream>
#include <gtest/gtest.h>
#include "helper/unitTestHelper.h"

#include <olb.h>

#include "core/vector.h"
#include "geometry/cuboid.h"

using namespace olb;
using namespace olb::descriptors;

//Tests for Hyperplane2D

//tests for method originAt
TEST(Hyperplane2D, testHyperplane2DOrigin1)
{
  Vector<double, 2> origin = {1, 1};
  Hyperplane2D<double> hyperPlane;
  hyperPlane.originAt(origin);

  EXPECT_NEAR(hyperPlane.origin[0], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[1], 1, 1e-10);
}
TEST(Hyperplane2D, testHyperplane2DOrigin2)
{
  Vector<double, 2> origin = {0, -1};
  Hyperplane2D<double> hyperPlane;
  hyperPlane.originAt(origin);

  EXPECT_NEAR(hyperPlane.origin[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[1], -1, 1e-10);
}
TEST(Hyperplane2D, testHyperplane2DOrigin3)
{
  Vector<double, 2> origin = {0, 0};
  Hyperplane2D<double> hyperPlane;
  hyperPlane.originAt(origin);

  EXPECT_NEAR(hyperPlane.origin[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[1], 0, 1e-10);
}
TEST(Hyperplane2D, testHyperplane2DOrigin4)
{
  Vector<double, 2> origin = {-34.2, 78};
  Hyperplane2D<double> hyperPlane;
  hyperPlane.originAt(origin);

  EXPECT_NEAR(hyperPlane.origin[0], -34.2, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[1], 78, 1e-10);
}
//tests for method centeredIn
TEST(Hyperplane2D, testHyperplane2DCenter1)
{
  Cuboid2D<double> cuboid({0, 0}, 1,{ 4, 2});
  Hyperplane2D<double> hyperPlane;
  hyperPlane.centeredIn(cuboid);

  EXPECT_NEAR(hyperPlane.origin[0], 2, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[1], 1, 1e-10);
}
TEST(Hyperplane2D, testHyperplane2DCenter2)
{
  Cuboid2D<double> cuboid({0, 0}, 0.5,{ 4, 2});
  Hyperplane2D<double> hyperPlane;
  hyperPlane.centeredIn(cuboid);

  EXPECT_NEAR(hyperPlane.origin[0], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[1], 0.5, 1e-10);
}
TEST(Hyperplane2D, testHyperplane2DCenter3)
{
  Cuboid2D<double> cuboid({1, 1}, 1,{ 9, 6});
  Hyperplane2D<double> hyperPlane;
  hyperPlane.centeredIn(cuboid);

  EXPECT_NEAR(hyperPlane.origin[0], 5.5, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[1], 4, 1e-10);
}
TEST(Hyperplane2D, testHyperplane2DCenter4)
{
  Cuboid2D<double> cuboid({-6.1, 4.8}, 2.5,{ 352, 76});
  Hyperplane2D<double> hyperPlane;
  hyperPlane.centeredIn(cuboid);

  EXPECT_NEAR(hyperPlane.origin[0], 433.9, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[1], 99.8, 1e-10);
}
TEST(Hyperplane2D, testHyperplane2DCenter5)
{
  Cuboid2D<double> cuboid({-6, -4}, 1,{ 12, 8});
  Hyperplane2D<double> hyperPlane;
  hyperPlane.centeredIn(cuboid);

  EXPECT_NEAR(hyperPlane.origin[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[1], 0, 1e-10);
}
//tests for method parallelTo
TEST(Hyperplane2D, testHyperplane2DDirection1)
{
  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> direction(1, 2);

  hyperPlane.parallelTo(direction);

  EXPECT_NEAR(hyperPlane.u[0], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 2, 1e-10);

  EXPECT_NEAR(hyperPlane.normal[0], 0.894427191, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], -0.4472135955, 1e-10);
}
TEST(Hyperplane2D, testHyperplane2DDirection2)
{
  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> direction(4, -0.5);

  hyperPlane.parallelTo(direction);

  EXPECT_NEAR(hyperPlane.u[0], 4, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], -0.5, 1e-10);

  EXPECT_NEAR(hyperPlane.normal[0], -0.1240347346, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], -0.99227787671366774, 1e-10);
}
TEST(Hyperplane2D, testHyperplane2DDirection3)
{
  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> direction(12, 0);

  hyperPlane.parallelTo(direction);

  EXPECT_NEAR(hyperPlane.u[0], 12, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 0, 1e-10);

  EXPECT_NEAR(hyperPlane.normal[0], 0., 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], -1, 1e-10);

  bool isParallelToX = hyperPlane.isParallelToX();
  EXPECT_TRUE(isParallelToX);
}
TEST(Hyperplane2D, testHyperplane2DDirection4)
{
  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> direction(0, 5);

  hyperPlane.parallelTo(direction);

  EXPECT_NEAR(hyperPlane.u[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 5, 1e-10);

  EXPECT_NEAR(hyperPlane.normal[0], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], 0, 1e-10);

  bool isParallelToY = hyperPlane.isParallelToY();
  EXPECT_TRUE(isParallelToY);
}
//tests for method normalTo
TEST(Hyperplane2D, testHyperplane2DNormal1)
{
  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> normalVec(2, 0);

  hyperPlane.normalTo(normalVec);

  EXPECT_NEAR(hyperPlane.u[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[0], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], 0, 1e-10);
  bool isParallelToY = hyperPlane.isParallelToY();
  EXPECT_TRUE(isParallelToY);
}
TEST(Hyperplane2D, testHyperplane2DNormal2)
{
  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> normalVec(0, -3);

  hyperPlane.normalTo(normalVec);

  EXPECT_NEAR(hyperPlane.u[0], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], -1, 1e-10);
  bool isParallelToX = hyperPlane.isParallelToX();
  EXPECT_TRUE(isParallelToX);
}
TEST(Hyperplane2D, testHyperplane2DNormal3)
{
  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> normalVec(2, 2);

  hyperPlane.normalTo(normalVec);

  EXPECT_NEAR(hyperPlane.u[0], 0.7071067812, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], -0.7071067812, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[0], 0.70710678118654746, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], 0.70710678118654746, 1e-10);
}

//tests for Hyperplane3D

//tests for method originAt
TEST(Hyperplane3D, testHyperplan3DOrigin1)
{
  Vector<double, 3> originVec(1, 1, 1);
  Hyperplane3D<double> hyperPlane;
  hyperPlane.originAt(originVec);

  EXPECT_NEAR(hyperPlane.origin[0], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[1], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[2], 1, 1e-10);
}
TEST(Hyperplane3D, testHyperplan3DOrigin2)
{
  Vector<double, 3> originVec(0, 0, 34);
  Hyperplane3D<double> hyperPlane;
  hyperPlane.originAt(originVec);

  EXPECT_NEAR(hyperPlane.origin[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[2], 34, 1e-10);
}
TEST(Hyperplane3D, testHyperplan3DOrigin3)
{
  Vector<double, 3> originVec(0.4, 45, -1);
  Hyperplane3D<double> hyperPlane;
  hyperPlane.originAt(originVec);

  EXPECT_NEAR(hyperPlane.origin[0], 0.4, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[1], 45, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[2], -1, 1e-10);
}
TEST(Hyperplane3D, testHyperplan3DOrigin4)
{
  Vector<double, 3> originVec(-8.5, 78, -21.8);
  Hyperplane3D<double> hyperPlane;
  hyperPlane.originAt(originVec);

  EXPECT_NEAR(hyperPlane.origin[0], -8.5, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[1], 78, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[2], -21.8, 1e-10);
}

//tests for method centeredIn
TEST(Hyperplane3D, testHyperplane3DCenter1)
{
  Cuboid3D<double> cuboid({0, 0, 0}, 1,{ 2, 4, 6});
  Hyperplane3D<double> hyperPlane;
  hyperPlane.centeredIn(cuboid);

  EXPECT_NEAR(hyperPlane.origin[0], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[1], 2, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[2], 3, 1e-10);
}
TEST(Hyperplane3D, testHyperplane3DCenter2)
{
  Cuboid3D<double> cuboid({0, 0, 0}, 0.5,{ 2, 4, 6});
  Hyperplane3D<double> hyperPlane;
  hyperPlane.centeredIn(cuboid);

  EXPECT_NEAR(hyperPlane.origin[0], 0.5, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[1], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[2], 1.5, 1e-10);
}
TEST(Hyperplane3D, testHyperplane3DCenter3)
{
  Cuboid3D<double> cuboid({0, 0, -3}, 2,{ 3, 1, 8});
  Hyperplane3D<double> hyperPlane;
  hyperPlane.centeredIn(cuboid);

  EXPECT_NEAR(hyperPlane.origin[0], 3, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[1], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[2], 5, 1e-10);
}
TEST(Hyperplane3D, testHyperplane3DCenter4)
{
  Cuboid3D<double> cuboid({34.2, -2, -53.2}, 2,{ 63, 2, 7});
  Hyperplane3D<double> hyperPlane;
  hyperPlane.centeredIn(cuboid);

  EXPECT_NEAR(hyperPlane.origin[0], 97.2, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[2], -46.2, 1e-10);
}
TEST(Hyperplane3D, testHyperplane3DCenter5)
{
  Cuboid3D<double> cuboid({-0.1, -9, 36.5}, 1.8,{ 2, 4, 6});
  Hyperplane3D<double> hyperPlane;
  hyperPlane.centeredIn(cuboid);

  EXPECT_NEAR(hyperPlane.origin[0], 1.7, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[1], -5.4, 1e-10);
  EXPECT_NEAR(hyperPlane.origin[2], 41.9, 1e-10);
}

//tests for method spannedBy
TEST(Hyperplane3D, testHyperplane3DSpan1)
{
  Vector<double, 3> uVec(1,0,0);
  Vector<double, 3> vVec(0,1,0);
  Hyperplane3D<double> hyperPlane;
  hyperPlane.spannedBy(uVec, vVec);

  EXPECT_NEAR(hyperPlane.u[0], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[2], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[1], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.v[2], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[2], 1, 1e-10);

  EXPECT_TRUE(hyperPlane.isXYPlane());
}
TEST(Hyperplane3D, testHyperplane3DSpan2)
{
  Vector<double, 3> uVec(0,3,0);
  Vector<double, 3> vVec(0,0,4);
  Hyperplane3D<double> hyperPlane;
  hyperPlane.spannedBy(uVec, vVec);

  EXPECT_NEAR(hyperPlane.u[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 3, 1e-10);
  EXPECT_NEAR(hyperPlane.u[2], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[2], 4, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[0], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[2], 0, 1e-10);

  EXPECT_TRUE(hyperPlane.isYZPlane());
}
TEST(Hyperplane3D, testHyperplane3DSpan3)
{
  Vector<double, 3> uVec(0,0,8.6);
  Vector<double, 3> vVec(1.3,0,0);
  Hyperplane3D<double> hyperPlane;
  hyperPlane.spannedBy(uVec, vVec);

  EXPECT_NEAR(hyperPlane.u[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[2], 8.6, 1e-10);
  EXPECT_NEAR(hyperPlane.v[0], 1.3, 1e-10);
  EXPECT_NEAR(hyperPlane.v[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[2], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[2], 0, 1e-10);

  EXPECT_TRUE(hyperPlane.isXZPlane());
}

//if you choose a more complex vector, the POSTCONDITION is not fulfilled, because of annihilation
/*
TEST(Hyperplane3D, testHyperplane3DSpan4)
{
  Vector<double, 3> uVec(-8, 22, 4);
  Vector<double, 3> vVec(0.6, 7.9, -42.25);
  Hyperplane3D<double> hyperPlane;
  hyperPlane.spannedBy(uVec, vVec);

  EXPECT_NEAR(hyperPlane.u[0], -8, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 22, 1e-10);
  EXPECT_NEAR(hyperPlane.u[2], 4, 1e-10);
  EXPECT_NEAR(hyperPlane.v[0], 0.6, 1e-10);
  EXPECT_NEAR(hyperPlane.v[1], 7.9, 1e-10);
  EXPECT_NEAR(hyperPlane.v[2], -42.25, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[0], -961.1, 1e-10); //normal[0] with/without normalize = -0.941451/-961.1
  EXPECT_NEAR(hyperPlane.normal[1], -335.6, 1e-10); //normal[1] with/without normalize = -0.328739/-335.6
  EXPECT_NEAR(hyperPlane.normal[2], -76.4, 1e-10); //normal[2] with/without normalize = -0.0748381/-76.4

  //u*normal without normalize = -5.68434e-13 ->violates the POSTCONDITION
  //v*normal without normalize = 0

  //u*normal with normalize = 0
  //v*normal with normalize = -4.44089e-16 ->violates the POSTCONDITION
}*/

//tests for method normalTo
TEST(Hyperplane3D, testHyperplane3DNormal1)
{
  Vector<double, 3> normalVec(1,0,0);
  Hyperplane3D<double> hyperPlane;
  hyperPlane.normalTo(normalVec);

  EXPECT_NEAR(hyperPlane.u[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.u[2], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[2], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[0], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[2], 0, 1e-10);

  EXPECT_TRUE(hyperPlane.isYZPlane());
}
TEST(Hyperplane3D, testHyperplane3DNormal2)
{
  Vector<double, 3> normalVec(0,2.3,0);
  Hyperplane3D<double> hyperPlane;
  hyperPlane.normalTo(normalVec);

  EXPECT_NEAR(hyperPlane.u[0], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[2], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[2], -1, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[2], 0, 1e-10);

  EXPECT_TRUE(hyperPlane.isXZPlane());
}
TEST(Hyperplane3D, testHyperplane3DNormal3)
{
  Vector<double, 3> normalVec(0,0,-87);
  Hyperplane3D<double> hyperPlane;
  hyperPlane.normalTo(normalVec);

  EXPECT_NEAR(hyperPlane.u[0], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[2], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[1], -1, 1e-10);
  EXPECT_NEAR(hyperPlane.v[2], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[2], -1, 1e-10);

  EXPECT_TRUE(hyperPlane.isXYPlane());
}
TEST(Hyperplane3D, testHyperplane3DNormal4)
{
  Vector<double, 3> normalVec(1,1,1);
  Hyperplane3D<double> hyperPlane;
  hyperPlane.normalTo(normalVec);

  EXPECT_NEAR(hyperPlane.u[0], 0.70710678118654746, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[2], -0.70710678118654746, 1e-10);
  EXPECT_NEAR(hyperPlane.v[0], -0.40824829046386307, 1e-10);
  EXPECT_NEAR(hyperPlane.v[1], 0.81649658092772615, 1e-10);
  EXPECT_NEAR(hyperPlane.v[2], -0.40824829046386307, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[0], 0.57735026918962584, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], 0.57735026918962584, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[2], 0.57735026918962584, 1e-10);
}
TEST(Hyperplane3D, testHyperplane3DNormal5)
{
  Vector<double, 3> normalVec(1.8,6.5,-2);
  Hyperplane3D<double> hyperPlane;
  hyperPlane.normalTo(normalVec);

  EXPECT_NEAR(hyperPlane.u[0], -0.74329414624716628, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[2], -0.66896473162244963, 1e-10);
  EXPECT_NEAR(hyperPlane.v[0], -0.61809873069581167, 1e-10);
  EXPECT_NEAR(hyperPlane.v[1], 0.38248160771262191, 1e-10);
  EXPECT_NEAR(hyperPlane.v[2], 0.68677636743979076, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[0], 0.25586670605399719, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], 0.92396310519498981, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[2], -0.28429634005999688, 1e-10);
}

//tests for method applyMatrixToSpan
TEST(Hyperplane3D, testHyperplane3DMatrixToSpan1)
{
  Hyperplane3D<double> hyperPlane;

  Vector<double, 3> normalVec(0,6,0);
  hyperPlane.normalTo(normalVec);

  hyperPlane.applyMatrixToSpan(
  {1, 0, 0},
  {0, 1, 0},
  {0, 0, 1});

  EXPECT_NEAR(hyperPlane.u[0], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[2], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[2], -1, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[2], 0, 1e-10);
}
TEST(Hyperplane3D, testHyperplane3DMatrixToSpan2)
{
  Hyperplane3D<double> hyperPlane;

  Vector<double, 3> normalVec(2, 0, 0);
  hyperPlane.normalTo(normalVec);
  hyperPlane.applyMatrixToSpan(
  {2, 0, 0},
  {0, 4, 0},
  {0, 0, 6});

  EXPECT_NEAR(hyperPlane.u[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 4, 1e-10);
  EXPECT_NEAR(hyperPlane.u[2], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[2], 6, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[0], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[2], 0, 1e-10);

  EXPECT_TRUE(hyperPlane.isYZPlane());
}
TEST(Hyperplane3D, testHyperplane3DMatrixToSpan3)
{
  Hyperplane3D<double> hyperPlane;

  Vector<double, 3> uVec(1,0,1);
  Vector<double, 3> vVec(0,1,0);
  hyperPlane.spannedBy(uVec, vVec);

  hyperPlane.applyMatrixToSpan(
  {24.2, 35.1, -152},
  {96,    3.3,   53},
  {0,    24.1, -8.8});

  EXPECT_NEAR(hyperPlane.u[0], -127.8, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 149, 1e-10);
  EXPECT_NEAR(hyperPlane.u[2], -8.8, 1e-10);
  EXPECT_NEAR(hyperPlane.v[0], 35.1, 1e-10);
  EXPECT_NEAR(hyperPlane.v[1], 3.3, 1e-10);
  EXPECT_NEAR(hyperPlane.v[2], 24.1, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[0], -0.70710678118654746, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[2], 0.70710678118654746, 1e-10);
}

//rotate Span around X
TEST(Hyperplane3D, testHyperplane3DRotateAroundX1)
{
  Hyperplane3D<double> hyperPlane;

  Vector<double, 3> normalVec(0,0,1);
  hyperPlane.normalTo(normalVec);

  hyperPlane.rotateSpanAroundX(M_PI/2);

  EXPECT_NEAR(hyperPlane.u[0], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[2], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[2], 1, 1e-10);

  EXPECT_TRUE(hyperPlane.isXYPlane());
}
TEST(Hyperplane3D, testHyperplane3DRotateAroundX2)
{
  Hyperplane3D<double> hyperPlane;

  Vector<double, 3> normalVec(0,1,0);
  hyperPlane.normalTo(normalVec);

  EXPECT_TRUE(hyperPlane.isXZPlane());

  hyperPlane.rotateSpanAroundX(M_PI/2);

  EXPECT_NEAR(hyperPlane.u[0], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[2], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[1], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.v[2], 0, 1e-10);

  //Attention: rotateSpanAround doesn't create a new normal, so you have to do it manually
  hyperPlane.spannedBy(hyperPlane.u, hyperPlane.v);
  EXPECT_TRUE(hyperPlane.isXYPlane());
}
TEST(Hyperplane3D, testHyperplane3DRotateAroundX3)
{
  Hyperplane3D<double> hyperPlane;

  Vector<double, 3> normalVec(0,0,1);
  hyperPlane.normalTo(normalVec);

  hyperPlane.rotateSpanAroundX(2*M_PI);

  EXPECT_NEAR(hyperPlane.u[0], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[2], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[1], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.v[2], 0, 1e-10);

  EXPECT_TRUE(hyperPlane.isXYPlane());
}

//rotate Span around Y
//if you rotate a XZSpan around Y the normal should stay the same
TEST(Hyperplane3D, testHyperplane3DRotateAroundY1)
{
  Hyperplane3D<double> hyperPlane;

  Vector<double, 3> uVec(1,0,0);
  Vector<double, 3> vVec(0,0,1);
  hyperPlane.spannedBy(uVec, vVec);

  hyperPlane.rotateSpanAroundY(M_PI);

  EXPECT_NEAR(hyperPlane.u[0], -1, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[2], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[2], -1, 1e-10);

  hyperPlane.spannedBy(hyperPlane.u, hyperPlane.v);
  EXPECT_NEAR(hyperPlane.normal[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], -1, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[2], 0, 1e-10);

  EXPECT_TRUE(hyperPlane.isXZPlane());
}
TEST(Hyperplane3D, testHyperplane3DRotateAroundY2)
{
  Hyperplane3D<double> hyperPlane;

  Vector<double, 3> uVec(0,1,0);
  Vector<double, 3> vVec(0,0,1);
  hyperPlane.spannedBy(uVec, vVec);

  hyperPlane.rotateSpanAroundY(M_PI);

  EXPECT_NEAR(hyperPlane.u[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.u[2], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[2], -1, 1e-10);

  EXPECT_TRUE(hyperPlane.isYZPlane());
}

//rotate Span around Z
//if you rotate a XYSpan around Z the normal should stay the same
TEST(Hyperplane3D, testHyperplane3DRotateAroundZ1)
{
  Hyperplane3D<double> hyperPlane;

  Vector<double, 3> normalVec(0,0,1);
  hyperPlane.normalTo(normalVec);

  hyperPlane.rotateSpanAroundZ(-8.12);

  hyperPlane.spannedBy(hyperPlane.u, hyperPlane.v);
  EXPECT_NEAR(hyperPlane.normal[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[2], 1, 1e-10);

  EXPECT_TRUE(hyperPlane.isXYPlane());
}
TEST(Hyperplane3D, testHyperplane3DRotateAroundZ2)
{
  Hyperplane3D<double> hyperPlane;

  Vector<double, 3> normalVec(1, 0, 0);
  hyperPlane.normalTo(normalVec);

  hyperPlane.rotateSpanAroundZ(M_PI/2);

  EXPECT_NEAR(hyperPlane.u[0], -1, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[2], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[2], 1, 1e-10);

  EXPECT_TRUE(hyperPlane.isYZPlane());
}
TEST(Hyperplane3D, testHyperplane3DRotateAroundZ3)
{
  Hyperplane3D<double> hyperPlane;

  Vector<double, 3> normalVec(0, 1, 0);
  hyperPlane.normalTo(normalVec);

  EXPECT_TRUE(hyperPlane.isXZPlane());

  hyperPlane.rotateSpanAroundZ(3*M_PI/2);

  EXPECT_NEAR(hyperPlane.u[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.u[1], -1, 1e-10);
  EXPECT_NEAR(hyperPlane.u[2], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[0], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.v[2], -1, 1e-10);

  hyperPlane.spannedBy(hyperPlane.u, hyperPlane.v);
  EXPECT_NEAR(hyperPlane.normal[0], 1, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[1], 0, 1e-10);
  EXPECT_NEAR(hyperPlane.normal[2], 0, 1e-10);

  EXPECT_TRUE(hyperPlane.isYZPlane());
}



//tests for HyperplaneLattice2D

//tests for method getVectorU

TEST(HyperplaneLattice2D, testHyperplaneLattice2DDirection1)
{
  CuboidDecomposition2D<double> geometry(0, 1, 2);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(-1, 0);
  Vector<double, 2> direction(1, 0);
  hyperPlane.parallelTo(direction);
  hyperPlane.originAt(origin);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane, 0.362624215);

  Vector<double, 2> vecU = hyperPlaneLattice.getVectorU();

  EXPECT_NEAR(vecU[0], 0.362624215, 1e-10);
  EXPECT_NEAR(vecU[1], 0, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DDirection2)
{
  CuboidDecomposition2D<double> geometry(0, 1, 2);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> direction(1, 1);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane, 1.0);

  Vector<double, 2> vecU = hyperPlaneLattice.getVectorU();

  EXPECT_NEAR(vecU[0], 0.70710678118654746, 1e-10);
  EXPECT_NEAR(vecU[1], 0.70710678118654746, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DDirection3)
{
  CuboidDecomposition2D<double> geometry(0, 1, 2);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> direction(3, 6);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 2> vecU = hyperPlaneLattice.getVectorU();

  EXPECT_NEAR(vecU[0], 0.4472135955, 1e-10);
  EXPECT_NEAR(vecU[1], 0.894427191, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DDirection4)
{
  CuboidDecomposition2D<double> geometry(0, 2*0.70710678118654746, 2);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> direction(-34, -34);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 2> vecU = hyperPlaneLattice.getVectorU();

  EXPECT_NEAR(vecU[0], -1, 1e-10);
  EXPECT_NEAR(vecU[1], -1, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DDirection5)
{
  CuboidDecomposition2D<double> geometry(0, 0.5, 2);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> direction(-5, 36);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 2> vecU = hyperPlaneLattice.getVectorU();

  EXPECT_NEAR(vecU[0], -0.06878418563734387, 1e-10);
  EXPECT_NEAR(vecU[1], 0.49524613658887584, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperPlaneLattice2DDirection6)
{
  CuboidDecomposition2D<double> geometry(0, 1, 2);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> direction(1, 1);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane, 2.0);

  Vector<double, 2> vecU = hyperPlaneLattice.getVectorU();

  EXPECT_NEAR(vecU[0], 1.4142135623730949, 1e-10);
  EXPECT_NEAR(vecU[1], 1.4142135623730949, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperPlaneLattice2DDirection7)
{
  CuboidDecomposition2D<double> geometry(0, 1, 2);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(0, 0);
  Vector<double, 2> direction(1, 1);
  hyperPlane.parallelTo(direction);
  hyperPlane.originAt(origin);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane, 3);

  Vector<double, 2> vecU = hyperPlaneLattice.getVectorU();

  EXPECT_NEAR(vecU[0], 0.70710678118654746, 1e-10);
  EXPECT_NEAR(vecU[1], 0.70710678118654746, 1e-10);
}

//tests for method getPhysOrigin

TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysOrigin1)
{
  CuboidDecomposition2D<double> geometry(0, 1, 2);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(0, 0);
  Vector<double, 2> direction(1, 0);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 2> physOrigin = hyperPlaneLattice.getPhysOrigin();
  Vector<double, 2> vecU = hyperPlaneLattice.getVectorU();

  EXPECT_NEAR(physOrigin[0], 0, 1e-10);
  EXPECT_NEAR(physOrigin[1], 0, 1e-10);
  EXPECT_NEAR(vecU[0], 1, 1e-10);
  EXPECT_NEAR(vecU[1], 0, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysOrigin2)
{
  CuboidDecomposition2D<double> geometry(8, 1, 4);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(8, 8);
  Vector<double, 2> direction(1, 1);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 2> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physOrigin[0], 8, 1e-10);
  EXPECT_NEAR(physOrigin[1], 8, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysOrigin3)
{
  CuboidDecomposition2D<double> geometry(8, 1, 4);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(0, 0);
  Vector<double, 2> direction(1, 1);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 2> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physOrigin[0], 7.77817459305202163, 1e-10);
  EXPECT_NEAR(physOrigin[1], 7.77817459305202163, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysOrigin4)
{
  CuboidDecomposition2D<double> geometry(8, 1, 4);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(1, 1);
  Vector<double, 2> direction(1, 1);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 2> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physOrigin[0], 8.0710678118654737, 1e-10);
  EXPECT_NEAR(physOrigin[1], 8.0710678118654737, 1e-10);
}
//be aware that your direction is fitted for deltaR if the spacing is default
//in that case vecU=2*0.5*util::sqrt(2)
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysOrigin5)
{
  CuboidDecomposition2D<double> geometry(8, 2, 4);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(0, 0);
  Vector<double, 2> direction(1, 1);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 2> physOrigin = hyperPlaneLattice.getPhysOrigin();
  Vector<double, 2> vecU = hyperPlaneLattice.getVectorU();

  EXPECT_NEAR(physOrigin[0], 7.0710678118654746, 1e-10);
  EXPECT_NEAR(physOrigin[1], 7.0710678118654746, 1e-10);
  EXPECT_NEAR(vecU[0], util::sqrt(2), 1e-10);
  EXPECT_NEAR(vecU[1], util::sqrt(2), 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysOrigin6)
{
  CuboidDecomposition2D<double> geometry(8, 1, 2);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(8, 8);
  Vector<double, 2> direction(0, 1);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 2> physOrigin = hyperPlaneLattice.getPhysOrigin();
  Vector<double, 2> vecU = hyperPlaneLattice.getVectorU();

  EXPECT_NEAR(physOrigin[0], 8, 1e-10);
  EXPECT_NEAR(physOrigin[1], 8, 1e-10);
  EXPECT_NEAR(vecU[0], 0, 1e-10);
  EXPECT_NEAR(vecU[1], 1, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysOrigin7)
{
  CuboidDecomposition2D<double> geometry({9, 8}, 1, {346, 3});

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(12, -345.);
  Vector<double, 2> direction(0, 1);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 2> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physOrigin[0], 12, 1e-10);
  EXPECT_NEAR(physOrigin[1], 8, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysOrigin8)
{
  CuboidDecomposition2D<double> geometry({12, 8}, 1, {5, 20});

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(0, 16);
  Vector<double, 2> direction(-1, 0);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 2> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physOrigin[0], 16, 1e-10);
  EXPECT_NEAR(physOrigin[1], 16, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysOrigin9)
{
  CuboidDecomposition2D<double> geometry({12, 8}, 1.5, {5, 20});

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(-1, 20);
  Vector<double, 2> direction(1, 0);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 2> physOrigin = hyperPlaneLattice.getPhysOrigin();
  Vector<double, 2> vecU = hyperPlaneLattice.getVectorU();

  EXPECT_NEAR(physOrigin[0], 12.5, 1e-10);
  EXPECT_NEAR(physOrigin[1], 20, 1e-10);
  EXPECT_NEAR(vecU[0], 1.5,1e-10);
  EXPECT_NEAR(vecU[1], 0,1e-10);
}
//resolution has no influence on the origin
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysOrigin10)
{
  CuboidDecomposition2D<double> geometry({12, 8}, 1, {5, 20});

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(0, 16);
  Vector<double, 2> direction(-1, 0);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane, 256);
  Vector<double, 2> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physOrigin[0], 16, 1e-10);
  EXPECT_NEAR(physOrigin[1], 16, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysOrigin11)
{
  CuboidDecomposition2D<double> geometry({12, 8}, 1, {5, 20});

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(0, 16);
  Vector<double, 2> direction(-1, 0);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane, 2.);
  Vector<double, 2> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physOrigin[0], 16, 1e-10);
  EXPECT_NEAR(physOrigin[1], 16, 1e-10);
}

//tests for method getPhysSpacing

//in default the Spacing is fitted for deltaR
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysSpacing1)
{
  CuboidDecomposition2D<double> geometry(0, 1, 2);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(0, 0);
  Vector<double, 2> direction(1, 0);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);
  double physSpacing = hyperPlaneLattice.getPhysSpacing();

  EXPECT_NEAR(physSpacing, 1, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysSpacing2)
{
  CuboidDecomposition2D<double> geometry(0, 4.945, 2);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(0, 0);
  Vector<double, 2> direction(35, 0);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);
  double physSpacing = hyperPlaneLattice.getPhysSpacing();
  Vector<double, 2> vecU = hyperPlaneLattice.getVectorU();

  EXPECT_NEAR(physSpacing, 4.945, 1e-10);
  EXPECT_NEAR(vecU[0], 4.945, 1e-10);
  EXPECT_NEAR(vecU[1], 0, 1e-10);
}
//be careful with a physSpacing too large, you could "jump" over the cuboid
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysSpacing3)
{
  CuboidDecomposition2D<double> geometry(0, 1, 2);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(-24, 0);
  Vector<double, 2> direction(-5, 0);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane, 25.);
  double physSpacing = hyperPlaneLattice.getPhysSpacing();
  Vector<double, 2> vecU = hyperPlaneLattice.getVectorU();

  EXPECT_NEAR(physSpacing, 25., 1e-10);
  EXPECT_NEAR(vecU[0], -25, 1e-10);
  EXPECT_NEAR(vecU[1], 0, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysSpacing4)
{
  CuboidDecomposition2D<double> geometry(0, 1, 2);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(0, 0);
  Vector<double, 2> direction(0, 9.57);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane, 10);
  double physSpacing = hyperPlaneLattice.getPhysSpacing();
  Vector<double, 2> vecU = hyperPlaneLattice.getVectorU();

  EXPECT_NEAR(physSpacing, 0.2, 1e-10);
  EXPECT_NEAR(vecU[0], 0, 1e-10);
  EXPECT_NEAR(vecU[1], 0.2, 1e-10);
}
//physSpacing is calculated after _origin
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysSpacing5)
{
  CuboidDecomposition2D<double> geometry(0, 1, 2);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(1, -1);
  Vector<double, 2> direction(0, 1);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane, 1);
  double physSpacing = hyperPlaneLattice.getPhysSpacing();
  Vector<double, 2> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physSpacing, 2, 1e-10);
  EXPECT_NEAR(physOrigin[0], 1, 1e-10);
  EXPECT_NEAR(physOrigin[1], 0, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysSpacing6)
{
  CuboidDecomposition2D<double> geometry(0, 1, 2);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(1, -0.5);
  Vector<double, 2> direction(0, 1);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane, 1);
  double physSpacing = hyperPlaneLattice.getPhysSpacing();
  Vector<double, 2> vecU = hyperPlaneLattice.getVectorU();

  EXPECT_NEAR(physSpacing, 2, 1e-10);
  EXPECT_NEAR(vecU[0], 0, 1e-10);
  EXPECT_NEAR(vecU[1], 2, 1e-10);
}

//tests for method getPhysR

TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysR1)
{
  CuboidDecomposition2D<double> geometry(0, 1, 2);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(0, 0);
  Vector<double, 2> direction(1, 0);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);
  Vector<double, 2> physR = hyperPlaneLattice.getPhysR(0);

  EXPECT_NEAR(physR[0], 0, 1e-10);
  EXPECT_NEAR(physR[1], 0, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysR2)
{
  CuboidDecomposition2D<double> geometry(0, 1, 2);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(41, 1);
  Vector<double, 2> direction(-1, 0);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);
  Vector<double, 2> physR = hyperPlaneLattice.getPhysR(1);

  EXPECT_NEAR(physR[0], 0, 1e-10);
  EXPECT_NEAR(physR[1], 1, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysR3)
{
  CuboidDecomposition2D<double> geometry(0, 3, 2);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(1414, 3);
  Vector<double, 2> direction(-4, 0);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);
  Vector<double, 2> physR = hyperPlaneLattice.getPhysR(0);

  EXPECT_NEAR(physR[0], 4, 1e-10);
  EXPECT_NEAR(physR[1], 3, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysR4)
{
  CuboidDecomposition2D<double> geometry(0, 3, 2);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(1414, 3);
  Vector<double, 2> direction(-4, 0);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane, 10);
  Vector<double, 2> physR = hyperPlaneLattice.getPhysR(7);

  EXPECT_NEAR(physR[0], -0.2, 1e-10);
  EXPECT_NEAR(physR[1], 3, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysR5)
{
  CuboidDecomposition2D<double> geometry(0, 1, 4);

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(2, 3);
  Vector<double, 2> direction(-1, 0);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane, 2.);
  Vector<double, 2> physR = hyperPlaneLattice.getPhysR(0);

  EXPECT_NEAR(physR[0], 2, 1e-10);
  EXPECT_NEAR(physR[1], 3, 1e-10);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DphysR6)
{
  CuboidDecomposition2D<double> geometry({8, 12}, 1, {6, 20});

  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(2, 3);
  Vector<double, 2> direction(1, 1);
  hyperPlane.originAt(origin);
  hyperPlane.parallelTo(direction);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);
  Vector<double, 2> physR = hyperPlaneLattice.getPhysR(0);

  EXPECT_NEAR(physR[0], 2+13*0.5*util::sqrt(2), 1e-10);
  EXPECT_NEAR(physR[1], 3+13*0.5*util::sqrt(2), 1e-10);

  physR =hyperPlaneLattice.getPhysR(1);
  EXPECT_NEAR(physR[0], 2+14*0.5*util::sqrt(2), 1e-10);
  EXPECT_NEAR(physR[1], 3+14*0.5*util::sqrt(2), 1e-10);
}

//tests for method getN

TEST(HyperplaneLattice2D, testHyperplaneLattice2DgetN1)
{
  CuboidDecomposition2D<double> geometry(8, 1, 4);
  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(0, 0);
  Vector<double, 2> direction(1, 1);
  hyperPlane.parallelTo(direction);
  hyperPlane.originAt(origin);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane, 1);

  int n = hyperPlaneLattice.getN();
  EXPECT_EQ(n, 1);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DgetN2)
{
  CuboidDecomposition2D<double> geometry({8, 8}, 1,{ 4, 4});
  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(0, 0);
  Vector<double, 2> direction(1, 1);
  hyperPlane.parallelTo(direction);
  hyperPlane.originAt(origin);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane, 462);

  int n = hyperPlaneLattice.getN();
  EXPECT_EQ(n, 462);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DgetN3)
{
  CuboidDecomposition2D<double> geometry({8, 8}, 1,{ 4, 4});
  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(8, 0);
  Vector<double, 2> direction(0, 1);
  hyperPlane.parallelTo(direction);
  hyperPlane.originAt(origin);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);

  int n = hyperPlaneLattice.getN();
  EXPECT_EQ(n, 4);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DgetN4)
{
  CuboidDecomposition2D<double> geometry({8, 8}, 3,{ 4, 4});
  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(8, 0);
  Vector<double, 2> direction(0, 1);
  hyperPlane.parallelTo(direction);
  hyperPlane.originAt(origin);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);

  int n = hyperPlaneLattice.getN();
  EXPECT_EQ(n, 4);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DgetN5)
{
  CuboidDecomposition2D<double> geometry({8, 8}, 1,{ 4, 4});
  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(8, 0);
  Vector<double, 2> direction(0, 1);
  hyperPlane.parallelTo(direction);
  hyperPlane.originAt(origin);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane, 2.);

  int n = hyperPlaneLattice.getN();
  EXPECT_EQ(n, 2);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DgetN6)
{
  CuboidDecomposition2D<double> geometry({8, 8}, 3,{ 4, 4});
  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(8, -1);
  Vector<double, 2> direction(0, 52);
  hyperPlane.parallelTo(direction);
  hyperPlane.originAt(origin);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane, 2.);

  int n = hyperPlaneLattice.getN();
  EXPECT_EQ(n, 6);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DgetN7)
{
  CuboidDecomposition2D<double> geometry({8, 8}, 1,{ 4, 4});
  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(0, 0);
  Vector<double, 2> direction(1, 1);
  hyperPlane.parallelTo(direction);
  hyperPlane.originAt(origin);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);

  int n = hyperPlaneLattice.getN();
  EXPECT_EQ(n, 6);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DgetN8)
{
  CuboidDecomposition2D<double> geometry({8, 8}, 2,{ 4, 4});
  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(0, 0);
  Vector<double, 2> direction(1, 1);
  hyperPlane.parallelTo(direction);
  hyperPlane.originAt(origin);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);

  int n = hyperPlaneLattice.getN();
  EXPECT_EQ(n, 6);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DgetN9)
{
  CuboidDecomposition2D<double> geometry({8, 8}, 3,{ 4, 4});
  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(0, 0);
  Vector<double, 2> direction(1, 1);
  hyperPlane.parallelTo(direction);
  hyperPlane.originAt(origin);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);

  int n = hyperPlaneLattice.getN();
  EXPECT_EQ(n, 5);
}
TEST(HyperplaneLattice2D, testHyperplaneLattice2DgetN10)
{
  CuboidDecomposition2D<double> geometry({8, 8}, 2,{ 100, 80});
  Hyperplane2D<double> hyperPlane;
  Vector<double, 2> origin(0, 0);
  Vector<double, 2> direction(208, 168);
  hyperPlane.parallelTo(direction);
  hyperPlane.originAt(origin);

  HyperplaneLattice2D<double> hyperPlaneLattice(geometry, hyperPlane);

  int n = hyperPlaneLattice.getN();
  EXPECT_EQ(n, 127);
}



//tests for HyperplaneLattice3D

//tests for method getVectorU and getVectorV

TEST(HyperplaneLattice3D, testHyperplaneLattice3Ddirection1)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 2, 2, 2});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> direction1(0, 1, 0);
  Vector<double, 3> direction2(0, 0, 1);
  hyperPlane.spannedBy(direction1, direction2);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 3> vecU = hyperPlaneLattice.getVectorU();
  Vector<double, 3> vecV = hyperPlaneLattice.getVectorV();

  EXPECT_NEAR(vecU[0], 0, 1e-10);
  EXPECT_NEAR(vecU[1], 1, 1e-10);
  EXPECT_NEAR(vecU[2], 0, 1e-10);
  EXPECT_NEAR(vecV[0], 0, 1e-10);
  EXPECT_NEAR(vecV[1], 0, 1e-10);
  EXPECT_NEAR(vecV[2], 1, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3Ddirection2)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 2, 2, 2});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> normal(1, 0, 0);
  hyperPlane.normalTo(normal);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 3> vecU = hyperPlaneLattice.getVectorU();
  Vector<double, 3> vecV = hyperPlaneLattice.getVectorV();

  EXPECT_NEAR(vecU[0], 0, 1e-10);
  EXPECT_NEAR(vecU[1], 1, 1e-10);
  EXPECT_NEAR(vecU[2], 0, 1e-10);
  EXPECT_NEAR(vecV[0], 0, 1e-10);
  EXPECT_NEAR(vecV[1], 0, 1e-10);
  EXPECT_NEAR(vecV[2], 1, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3Ddirection3)
{
  CuboidDecomposition3D<double> geometry({41, 4, -4}, 1,{ 53, 4, 34});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> direction1(1, 1, 0);
  Vector<double, 3> direction2(0, 0, 4);
  hyperPlane.spannedBy(direction1, direction2);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 3> vecU = hyperPlaneLattice.getVectorU();
  Vector<double, 3> vecV = hyperPlaneLattice.getVectorV();

  EXPECT_NEAR(vecU[0], 0.70710678118654746, 1e-10);
  EXPECT_NEAR(vecU[1], 0.70710678118654746, 1e-10);
  EXPECT_NEAR(vecU[2], 0, 1e-10);
  EXPECT_NEAR(vecV[0], 0, 1e-10);
  EXPECT_NEAR(vecV[1], 0, 1e-10);
  EXPECT_NEAR(vecV[2], 1, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3Ddirection4)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 0.5,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> direction1(2, 0, 0);
  Vector<double, 3> direction2(0, -1, 0);
  hyperPlane.spannedBy(direction1, direction2);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 3> vecU = hyperPlaneLattice.getVectorU();
  Vector<double, 3> vecV = hyperPlaneLattice.getVectorV();

  EXPECT_NEAR(vecU[0], 0.5, 1e-10);
  EXPECT_NEAR(vecU[1], 0, 1e-10);
  EXPECT_NEAR(vecU[2], 0, 1e-10);
  EXPECT_NEAR(vecV[0], 0, 1e-10);
  EXPECT_NEAR(vecV[1], -0.5, 1e-10);
  EXPECT_NEAR(vecV[2], 0, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3Ddirection5)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 0.342545,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> direction1(2, 0, 0);
  Vector<double, 3> direction2(0, 1, 0);
  hyperPlane.spannedBy(direction1, direction2);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 3> vecU = hyperPlaneLattice.getVectorU();
  Vector<double, 3> vecV = hyperPlaneLattice.getVectorV();

  EXPECT_NEAR(vecU[0], 0.342545, 1e-10);
  EXPECT_NEAR(vecU[1], 0, 1e-10);
  EXPECT_NEAR(vecU[2], 0, 1e-10);
  EXPECT_NEAR(vecV[0], 0, 1e-10);
  EXPECT_NEAR(vecV[1], 0.342545, 1e-10);
  EXPECT_NEAR(vecV[2], 0, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3Ddirection6)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> direction1(2, 0, 0);
  Vector<double, 3> direction2(0, 1, 0);
  hyperPlane.spannedBy(direction1, direction2);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane, 0.342545);

  Vector<double, 3> vecU = hyperPlaneLattice.getVectorU();
  Vector<double, 3> vecV = hyperPlaneLattice.getVectorV();

  EXPECT_NEAR(vecU[0], 0.342545, 1e-10);
  EXPECT_NEAR(vecU[1], 0, 1e-10);
  EXPECT_NEAR(vecU[2], 0, 1e-10);
  EXPECT_NEAR(vecV[0], 0, 1e-10);
  EXPECT_NEAR(vecV[1], 0.342545, 1e-10);
  EXPECT_NEAR(vecV[2], 0, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3Ddirection7)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> direction1(2, 0, 0);
  Vector<double, 3> direction2(0, 1, 0);
  hyperPlane.spannedBy(direction1, direction2);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane, 4);

  Vector<double, 3> vecU = hyperPlaneLattice.getVectorU();
  Vector<double, 3> vecV = hyperPlaneLattice.getVectorV();

  EXPECT_NEAR(vecU[0], 1, 1e-10);
  EXPECT_NEAR(vecU[1], 0, 1e-10);
  EXPECT_NEAR(vecU[2], 0, 1e-10);
  EXPECT_NEAR(vecV[0], 0, 1e-10);
  EXPECT_NEAR(vecV[1], 1, 1e-10);
  EXPECT_NEAR(vecV[2], 0, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3Ddirection8)
{
  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(1, -2, 5);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(3, 0, 0);
  Vector<double, 3> directionV(0, 0, -8);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(hyperPlane, 2., 2, 2);

  Vector<double, 3> vecU = hyperPlaneLattice.getVectorU();
  Vector<double, 3> vecV = hyperPlaneLattice.getVectorV();

  EXPECT_NEAR(vecU[0], 2, 1e-10);
  EXPECT_NEAR(vecU[1], 0, 1e-10);
  EXPECT_NEAR(vecU[2], 0, 1e-10);
  EXPECT_NEAR(vecV[0], 0, 1e-10);
  EXPECT_NEAR(vecV[1], 0, 1e-10);
  EXPECT_NEAR(vecV[2], -2, 1e-10);
}

//tests for method getPhysOrigin

TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysOrigin1)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 3> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physOrigin[0], 0, 1e-10);
  EXPECT_NEAR(physOrigin[1], 0, 1e-10);
  EXPECT_NEAR(physOrigin[2], 0, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysOrigin2)
{
  CuboidDecomposition3D<double> geometry({3, 0, 0}, 1,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(0, 0, 0);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 3> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physOrigin[0], 3, 1e-10);
  EXPECT_NEAR(physOrigin[1], 0, 1e-10);
  EXPECT_NEAR(physOrigin[2], 0, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysOrigin3)
{
  CuboidDecomposition3D<double> geometry({4, 5, 0}, 1,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(6, 6, 0);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(4, 0, 0);
  Vector<double, 3> directionV(0, 3, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 3> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physOrigin[0], 4, 1e-10);
  EXPECT_NEAR(physOrigin[1], 5, 1e-10);
  EXPECT_NEAR(physOrigin[2], 0, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysOrigin4)
{
  CuboidDecomposition3D<double> geometry({4, 5, 2}, 1,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(-2, -4, 3);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(4, 0, 0);
  Vector<double, 3> directionV(0, 3, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 3> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physOrigin[0], 4, 1e-10);
  EXPECT_NEAR(physOrigin[1], 5, 1e-10);
  EXPECT_NEAR(physOrigin[2], 3, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysOrigin5)
{
  CuboidDecomposition3D<double> geometry({4, 5, 2}, 1,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(5, 6, 3);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(0, -2, 0);
  Vector<double, 3> directionV(0, 0, 1);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 3> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physOrigin[0], 5, 1e-10);
  EXPECT_NEAR(physOrigin[1], 8, 1e-10);
  EXPECT_NEAR(physOrigin[2], 2, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysOrigin6)
{
  CuboidDecomposition3D<double> geometry({4, 5, 2}, 3,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(5, 6, 3);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(0, -2, 0);
  Vector<double, 3> directionV(0, 0, 1);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 3> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physOrigin[0], 5, 1e-10);
  EXPECT_NEAR(physOrigin[1], 15, 1e-10);
  EXPECT_NEAR(physOrigin[2], 3, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysOrigin7)
{
  CuboidDecomposition3D<double> geometry({4, 5, 2}, 1,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(0, 6, 0);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 0, 1);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane, 465);

  Vector<double, 3> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physOrigin[0], 4, 1e-10);
  EXPECT_NEAR(physOrigin[1], 6, 1e-10);
  EXPECT_NEAR(physOrigin[2], 2, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysOrigin8)
{
  CuboidDecomposition3D<double> geometry({4, 5, 2}, 1,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(0, 6, 0);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 0, 1);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane, 2.);

  Vector<double, 3> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physOrigin[0], 4, 1e-10);
  EXPECT_NEAR(physOrigin[1], 6, 1e-10);
  EXPECT_NEAR(physOrigin[2], 2, 1e-10);
}
//be aware of rounding, it can have a big impact on your origin
//as you can see in the two following tests
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysOrigin9)
{
  CuboidDecomposition3D<double> geometry({4, 5, 2}, 1,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(0.5, 6, 0);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 0, 1);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane, 3.);

  Vector<double, 3> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physOrigin[0], 3.5, 1e-10);
  EXPECT_NEAR(physOrigin[1], 6, 1e-10);
  EXPECT_NEAR(physOrigin[2], 3, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysOrigin10)
{
  CuboidDecomposition3D<double> geometry({4, 5, 2}, 1,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(0.4, 6, 0);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 0, 1);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane, 3.);

  Vector<double, 3> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physOrigin[0], 6.4, 1e-10);
  EXPECT_NEAR(physOrigin[1], 6, 1e-10);
  EXPECT_NEAR(physOrigin[2], 3, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysOrigin11)
{
  CuboidDecomposition3D<double> geometry({3, 5, 2}, 1,{ 4, 6, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(0, 0, 0);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 1, 0);
  Vector<double, 3> directionV(0, 0, 1);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 3> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physOrigin[0], 3.5*util::sqrt(2), 1e-10);
  EXPECT_NEAR(physOrigin[1], 3.5*util::sqrt(2), 1e-10);
  EXPECT_NEAR(physOrigin[2], 2, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysOrigin12)
{
  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(0.5, 781, -32.3);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(hyperPlane, 1., 2, 2);

  Vector<double, 3> physOrigin = hyperPlaneLattice.getPhysOrigin();

  EXPECT_NEAR(physOrigin[0], 0.5, 1e-10);
  EXPECT_NEAR(physOrigin[1], 781, 1e-10);
  EXPECT_NEAR(physOrigin[2], -32.3, 1e-10);
}

//tests for method getPhysSpacing

TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysSpacing1)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 3, 3, 3});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(0, 0, 0);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  double physSpacing = hyperPlaneLattice.getPhysSpacing();

  EXPECT_NEAR(physSpacing, 1, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysSpacing2)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1.87,{ 3, 3, 3});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(0, 0, 0);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  double physSpacing = hyperPlaneLattice.getPhysSpacing();
  Vector<double, 3> vecU = hyperPlaneLattice.getVectorU();
  Vector<double, 3> vecV = hyperPlaneLattice.getVectorV();

  EXPECT_NEAR(physSpacing, 1.87, 1e-10);
  EXPECT_NEAR(vecU[0], 1.87, 1e-10);
  EXPECT_NEAR(vecU[1], 0, 1e-10);
  EXPECT_NEAR(vecU[2], 0, 1e-10);
  EXPECT_NEAR(vecV[0], 0, 1e-10);
  EXPECT_NEAR(vecV[1], 1.87, 1e-10);
  EXPECT_NEAR(vecV[2], 0, 1e-10);
}
//be careful with a physSpacing too large, you could "jump" over the cuboid
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysSpacing3)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 3, 3, 3});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(22, -21, 1);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(-1, 0, 0);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane, 21.);

  double physSpacing = hyperPlaneLattice.getPhysSpacing();
  Vector<double, 3> vecU = hyperPlaneLattice.getVectorU();
  Vector<double, 3> vecV = hyperPlaneLattice.getVectorV();

  EXPECT_NEAR(physSpacing, 21, 1e-10);
  EXPECT_NEAR(vecU[0], -21, 1e-10);
  EXPECT_NEAR(vecU[1], 0, 1e-10);
  EXPECT_NEAR(vecU[2], 0, 1e-10);
  EXPECT_NEAR(vecV[0], 0, 1e-10);
  EXPECT_NEAR(vecV[1], 21, 1e-10);
  EXPECT_NEAR(vecV[2], 0, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysSpacing4)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 3, 3, 3});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(0, 0, 1);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(-1, 0, 0);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane, 9);

  double physSpacing = hyperPlaneLattice.getPhysSpacing();
  Vector<double, 3> vecU = hyperPlaneLattice.getVectorU();
  Vector<double, 3> vecV = hyperPlaneLattice.getVectorV();

  EXPECT_NEAR(physSpacing, 0.333333333333333333, 1e-10);
  EXPECT_NEAR(vecU[0], -0.333333333333333333, 1e-10);
  EXPECT_NEAR(vecU[1], 0, 1e-10);
  EXPECT_NEAR(vecU[2], 0, 1e-10);
  EXPECT_NEAR(vecV[0], 0, 1e-10);
  EXPECT_NEAR(vecV[1], 0.333333333333333333, 1e-10);
  EXPECT_NEAR(vecV[2], 0, 1e-10);
}
//the resolution is set for the direction with more latticepoints
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysSpacing5)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 1, 4, 2});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(0, 0, 1);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(-1, 0, 0);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane, 10);

  double physSpacing = hyperPlaneLattice.getPhysSpacing();
  Vector<double, 3> vecU = hyperPlaneLattice.getVectorU();
  Vector<double, 3> vecV = hyperPlaneLattice.getVectorV();

  EXPECT_NEAR(physSpacing, 0.4, 1e-10);
  EXPECT_NEAR(vecU[0], -0.4, 1e-10);
  EXPECT_NEAR(vecU[1], 0, 1e-10);
  EXPECT_NEAR(vecU[2], 0, 1e-10);
  EXPECT_NEAR(vecV[0], 0, 1e-10);
  EXPECT_NEAR(vecV[1], 0.4, 1e-10);
  EXPECT_NEAR(vecV[2], 0, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysSpacing6)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 1, 4, 2});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(0, 0, 1);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 1);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane, 10);

  double physSpacing = hyperPlaneLattice.getPhysSpacing();
  Vector<double, 3> vecU = hyperPlaneLattice.getVectorU();
  Vector<double, 3> vecV = hyperPlaneLattice.getVectorV();

  EXPECT_NEAR(physSpacing, 0.4, 1e-10);
  EXPECT_NEAR(vecU[0], 0.28284271247461906, 1e-10);
  EXPECT_NEAR(vecU[1], 0, 1e-10);
  EXPECT_NEAR(vecU[2], 0.28284271247461906, 1e-10);
  EXPECT_NEAR(vecV[0], 0, 1e-10);
  EXPECT_NEAR(vecV[1], 0.4, 1e-10);
  EXPECT_NEAR(vecV[2], 0, 1e-10);
}
//when physSpacing is given, _origin is calculated with it
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysSpacing7)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 0.5,{ 4, 8, 6});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(0, -2.3, -1);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(0, 0, 1);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane, 2.);

  double physSpacing = hyperPlaneLattice.getPhysSpacing();
  Vector<double, 3> vecU = hyperPlaneLattice.getVectorU();
  Vector<double, 3> vecV = hyperPlaneLattice.getVectorV();

  EXPECT_NEAR(physSpacing, 2, 1e-10);
  EXPECT_NEAR(vecU[0], 0, 1e-10);
  EXPECT_NEAR(vecU[1], 0, 1e-10);
  EXPECT_NEAR(vecU[2], 2, 1e-10);
  EXPECT_NEAR(vecV[0], 0, 1e-10);
  EXPECT_NEAR(vecV[1], 2, 1e-10);
  EXPECT_NEAR(vecV[2], 0, 1e-10);
}
//when physSpacing is not given, it is calculated after _origin
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysSpacing8)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(0, -2, -1);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(0, 0, 1);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane, 1);

  double physSpacing = hyperPlaneLattice.getPhysSpacing();
  Vector<double, 3> vecU = hyperPlaneLattice.getVectorU();
  Vector<double, 3> vecV = hyperPlaneLattice.getVectorV();

  EXPECT_NEAR(physSpacing, 4, 1e-10);
  EXPECT_NEAR(vecU[0], 0, 1e-10);
  EXPECT_NEAR(vecU[1], 0, 1e-10);
  EXPECT_NEAR(vecU[2], 4, 1e-10);
  EXPECT_NEAR(vecV[0], 0, 1e-10);
  EXPECT_NEAR(vecV[1], 4, 1e-10);
  EXPECT_NEAR(vecV[2], 0, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysSpacing9)
{
  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(6, -3, 3);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(hyperPlane, 3.8, 2, 2);

  double physSpacing = hyperPlaneLattice.getPhysSpacing();

  EXPECT_NEAR(physSpacing, 3.8, 1e-10);
}

//tests for method getPhysR

TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysR1)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(0, 0, 0);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 3> physR = hyperPlaneLattice.getPhysR(0, 0);

  EXPECT_NEAR(physR[0], 0, 1e-10);
  EXPECT_NEAR(physR[1], 0, 1e-10);
  EXPECT_NEAR(physR[2], 0, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysR2)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(-324, -2, 1);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 3> physR = hyperPlaneLattice.getPhysR(0, 1);

  EXPECT_NEAR(physR[0], 0, 1e-10);
  EXPECT_NEAR(physR[1], 1, 1e-10);
  EXPECT_NEAR(physR[2], 1, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysR3)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(-324, -2, 1);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(-1, 0, 0);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 3> physR = hyperPlaneLattice.getPhysR(1, 2);

  EXPECT_NEAR(physR[0], 2, 1e-10);
  EXPECT_NEAR(physR[1], 2, 1e-10);
  EXPECT_NEAR(physR[2], 1, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysR4)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 4, 25, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(3, 0, 41);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(0, 0, 1);
  Vector<double, 3> directionV(0, -1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 3> physR = hyperPlaneLattice.getPhysR(2, 2);

  EXPECT_NEAR(physR[0], 3, 1e-10);
  EXPECT_NEAR(physR[1], 22, 1e-10);
  EXPECT_NEAR(physR[2], 2, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysR5)
{
  CuboidDecomposition3D<double> geometry({1, -2, 0}, 1,{ 4, 25, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(3, 0, 4);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 1);
  Vector<double, 3> directionV(0, -1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 3> physR = hyperPlaneLattice.getPhysR(2, 2);

  EXPECT_NEAR(physR[0], 2.2928932188134512, 1e-10);
  EXPECT_NEAR(physR[1], 20, 1e-10);
  EXPECT_NEAR(physR[2], 3.2928932188134508, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysR6)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 2,{ 4, 25, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(3, 0, 42);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(0, 0, 1);
  Vector<double, 3> directionV(0, -1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  Vector<double, 3> physR = hyperPlaneLattice.getPhysR(3, 24);

  EXPECT_NEAR(physR[0], 3, 1e-10);
  EXPECT_NEAR(physR[1], 0, 1e-10);
  EXPECT_NEAR(physR[2], 6, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysR7)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 4, 25, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(3, 0, 41);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(0, 0, 1);
  Vector<double, 3> directionV(0, -1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane, 5);

  Vector<double, 3> physR = hyperPlaneLattice.getPhysR(2, 2);

  EXPECT_NEAR(physR[0], 3, 1e-10);
  EXPECT_NEAR(physR[1], 14, 1e-10);
  EXPECT_NEAR(physR[2], 10, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysR8)
{
  CuboidDecomposition3D<double> geometry({8, 8, 8}, 1,{ 4, 4, 2});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(3, 10, 2);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(10, 0, 0);
  Vector<double, 3> directionV(0, 0, 1);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane, 3.);

  Vector<double, 3> physR = hyperPlaneLattice.getPhysR(0, 1);

  EXPECT_NEAR(physR[0], 9, 1e-10);
  EXPECT_NEAR(physR[1], 10, 1e-10);
  EXPECT_NEAR(physR[2], 11, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysR9)
{
  CuboidDecomposition3D<double> geometry({8, 8, 8}, 1,{ 4, 4, 2});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(3, 10, 2);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(0, 10, 0);
  Vector<double, 3> directionV(1, 0, 1);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane, util::sqrt(2));

  Vector<double, 3> physR = hyperPlaneLattice.getPhysR(0, 2);

  EXPECT_NEAR(physR[0], 11, 1e-10);
  EXPECT_NEAR(physR[1], 8.5857864376269006, 1e-10);
  EXPECT_NEAR(physR[2], 10, 1e-10);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DphysR10)
{
  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(-7, 21., 2.4);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(hyperPlane, 1.5, 35, 74);

  Vector<double, 3> physR = hyperPlaneLattice.getPhysR(8, 5);

  EXPECT_NEAR(physR[0], 5, 1e-10);
  EXPECT_NEAR(physR[1], 28.5, 1e-10);
  EXPECT_NEAR(physR[2], 2.4, 1e-10);
}

//tests for method getN

TEST(HyperplaneLattice3D, testHyperplaneLattice3DgetN1)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(0, 0, 0);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  double nX = hyperPlaneLattice.getNx();
  double nY = hyperPlaneLattice.getNy();

  EXPECT_EQ(nX, 4);
  EXPECT_EQ(nY, 4);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DgetN2)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 4, 2, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(0, 0, 0);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  double nX = hyperPlaneLattice.getNx();
  double nY = hyperPlaneLattice.getNy();

  EXPECT_EQ(nX, 4);
  EXPECT_EQ(nY, 2);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DgetN3)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 4, 2, 8});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(0, 0, 0);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 0, 8);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  double nX = hyperPlaneLattice.getNx();
  double nY = hyperPlaneLattice.getNy();

  EXPECT_EQ(nX, 4);
  EXPECT_EQ(nY, 8);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DgetN4)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 2,{ 4, 4, 4});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(3, 2, 0);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(0, 3, 0);
  Vector<double, 3> directionV(0, 0, -2);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  double nX = hyperPlaneLattice.getNx();
  double nY = hyperPlaneLattice.getNy();

  EXPECT_EQ(nX, 4);
  EXPECT_EQ(nY, 4);
}
//be aware how the origin can affect the latticepoints, as you can see in the two following tests
TEST(HyperplaneLattice3D, testHyperplaneLattice3DgetN5)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 3,{ 4, 6, 8});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(-1.4, 2, 0);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  double nX = hyperPlaneLattice.getNx();
  double nY = hyperPlaneLattice.getNy();

  EXPECT_EQ(nX, 4);
  EXPECT_EQ(nY, 6);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DgetN6)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 3,{ 4, 6, 8});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(-1.5, 2, 0);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  double nX = hyperPlaneLattice.getNx();
  double nY = hyperPlaneLattice.getNy();

  EXPECT_EQ(nX, 3);
  EXPECT_EQ(nY, 6);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DgetN7)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 4, 6, 10});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(-1.5, 2, 9);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane, 1);

  double nX = hyperPlaneLattice.getNx();
  double nY = hyperPlaneLattice.getNy();

  EXPECT_EQ(nX, 1);
  EXPECT_EQ(nY, 1);
}
//the resolution is fitted for the direction with more latticepoints
TEST(HyperplaneLattice3D, testHyperplaneLattice3DgetN8)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 4, 6, 8});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(-10, 2, -3);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 0, 1);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane, 356);

  double nX = hyperPlaneLattice.getNx();
  double nY = hyperPlaneLattice.getNy();

  EXPECT_EQ(nX, 179);
  EXPECT_EQ(nY, 356);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DgetN9)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 4, 6, 8});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(-10, 2, -3);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 0, 1);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane, 2.);

  double nX = hyperPlaneLattice.getNx();
  double nY = hyperPlaneLattice.getNy();

  EXPECT_EQ(nX, 2);
  EXPECT_EQ(nY, 4);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DgetN10)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 4, 6, 8});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(-10, 2, -3);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 0, 0);
  Vector<double, 3> directionV(0, 0, 1);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane, 4.);

  double nX = hyperPlaneLattice.getNx();
  double nY = hyperPlaneLattice.getNy();

  EXPECT_EQ(nX, 1);
  EXPECT_EQ(nY, 2);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DgetN11)
{
  CuboidDecomposition3D<double> geometry({0, 0, 0}, 1,{ 200, 200, 200});

  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(0, 0, 0);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(1, 1, 0);
  Vector<double, 3> directionV(0, 0, 1);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(geometry, hyperPlane);

  double nX = hyperPlaneLattice.getNx();
  double nY = hyperPlaneLattice.getNy();

  EXPECT_EQ(nX, 283);
  EXPECT_EQ(nY, 200);
}
TEST(HyperplaneLattice3D, testHyperplaneLattice3DgetN12)
{
  Hyperplane3D<double> hyperPlane;
  Vector<double, 3> origin(4, -35, 81);
  hyperPlane.originAt(origin);
  Vector<double, 3> directionU(-25, 0, 0);
  Vector<double, 3> directionV(0, 1, 0);
  hyperPlane.spannedBy(directionU, directionV);

  HyperplaneLattice3D<double> hyperPlaneLattice(hyperPlane, 1., 213, 8);

  double nX = hyperPlaneLattice.getNx();
  double nY = hyperPlaneLattice.getNy();

  EXPECT_EQ(nX, 213);
  EXPECT_EQ(nY, 8);
}
