/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Nicolas Hafen
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
using namespace util;

typedef double T;

//maxElement

TEST(TestVectorHelpers, maxElementDoubleLate)
{
  Vector<T,3> vector(0.2, 0.4, 0.9);
  T res = max_element(vector);
  EXPECT_EQ(res,0.9);
}

TEST(TestVectorHelpers, maxElementDoubleMid)
{
  Vector<T,3> vector(0.2, 1.1, 0.9);
  T res = max_element(vector);
  EXPECT_EQ(res,1.1);
}

TEST(TestVectorHelpers, maxElementDoubleEarly)
{
  Vector<T,4> vector(2.003,0.2, -10.1, 0.9);
  T res = max_element(vector);
  EXPECT_EQ(res,2.003);
}

TEST(TestVectorHelpers, maxElementIntLate)
{
  Vector<int,3> vector(2, 4, 9);
  T res = max_element(vector);
  EXPECT_EQ(res,9);
}

TEST(TestVectorHelpers, maxElementIntMid)
{
  Vector<int,3> vector(2, 11, 9);
  T res = max_element(vector);
  EXPECT_EQ(res,11);
}

TEST(TestVectorHelpers, maxElementIntEarly)
{
  Vector<int,4> vector(12,2,1,9);
  T res = max_element(vector);
  EXPECT_EQ(res,12);
}

//minElement

TEST(TestVectorHelpers, minElementDoubleLate)
{
  Vector<T,3> vector(0.2, 0.4, 0.05);
  T res = min_element(vector);
  EXPECT_EQ(res,0.05);
}

TEST(TestVectorHelpers, minElementDoubleMid)
{
  Vector<T,3> vector(0.2, 0.01, 0.9);
  T res = min_element(vector);
  EXPECT_EQ(res,0.01);
}

TEST(TestVectorHelpers, minElementDoubleEarly)
{
  Vector<T,5> vector(0.002,2.003,0.2, 1.1, 0.9);
  T res = min_element(vector);
  EXPECT_EQ(res,0.002);
}

TEST(TestVectorHelpers, minElementIntLate)
{
  Vector<int,3> vector(2, 4, 1);
  T res = min_element(vector);
  EXPECT_EQ(res,1);
}

TEST(TestVectorHelpers, minElementIntMid)
{
  Vector<int,3> vector(2, 1, 9);
  T res = min_element(vector);
  EXPECT_EQ(res,1);
}

TEST(TestVectorHelpers, minElementIntEarly)
{
  Vector<int,6> vector(0,1,1,2,1,9);
  T res = min_element(vector);
  EXPECT_EQ(res,0);
}

//maxElementPos

TEST(TestVectorHelpers, maxElementPosDoubleLate)
{
  Vector<T,3> vector(0.2, 0.4, 0.9);
  unsigned res = maxElementPos(vector);
  EXPECT_EQ(res,2);
}

TEST(TestVectorHelpers, maxElementPosDoubleMid)
{
  Vector<T,3> vector(-3.2, 1.1, 0.9);
  unsigned res = maxElementPos(vector);
  EXPECT_EQ(res,1);
}

TEST(TestVectorHelpers, maxElementPosDoubleEarly)
{
  Vector<T,4> vector(2.003, 0.2, 1.1, 0.9);
  unsigned res = maxElementPos(vector);
  EXPECT_EQ(res,0);
}

TEST(TestVectorHelpers, maxElementPosIntLate)
{
  Vector<int,3> vector(2, -4, 9);
  unsigned res = maxElementPos(vector);
  EXPECT_EQ(res,2);
}

TEST(TestVectorHelpers, maxElementPosIntMid)
{
  Vector<int,3> vector(2, 11, 9);
  unsigned res = maxElementPos(vector);
  EXPECT_EQ(res,1);
}

TEST(TestVectorHelpers, maxElementPosIntEarly)
{
  Vector<int,4> vector(12,-2,1,-9);
  unsigned res = maxElementPos(vector);
  EXPECT_EQ(res,0);
}

//minElementPos

TEST(TestVectorHelpers, minElementPosDoubleLate)
{
  Vector<T,3> vector(0.2, 0.4, 0.1);
  unsigned res = minElementPos(vector);
  EXPECT_EQ(res,2);
}

TEST(TestVectorHelpers, minElementPosDoubleMid)
{
  Vector<T,3> vector(0.2, 0.01, 0.9);
  unsigned res = minElementPos(vector);
  EXPECT_EQ(res,1);
}

TEST(TestVectorHelpers, minElementPosDoubleEarly)
{
  Vector<T,5> vector(0.002,2.003,0.2, 1.1, 0.9);
  unsigned res = minElementPos(vector);
  EXPECT_EQ(res,0);
}

TEST(TestVectorHelpers, minElementPosIntLate)
{
  Vector<int,3> vector(2, 4, 1);
  unsigned res = minElementPos(vector);
  EXPECT_EQ(res,2);
}

TEST(TestVectorHelpers, minElementPosIntMid)
{
  Vector<int,3> vector(2, 1, 9);
  unsigned res = minElementPos(vector);
  EXPECT_EQ(res,1);
}

TEST(TestVectorHelpers, minElementPosIntEarly)
{
  Vector<int,6> vector(0,1,1,2,1,9);
  T res = minElementPos(vector);
  EXPECT_EQ(res,0);
}

//maxElementAbs

TEST(TestVectorHelpers, maxElementAbsDoubleLate)
{
  Vector<T,3> vector(0.2, -0.4, 0.9);
  T res = maxElementAbs(vector);
  EXPECT_EQ(res,0.9);
}

TEST(TestVectorHelpers, maxElementAbsDoubleMid)
{
  Vector<T,3> vector(0.2, -1.1, 0.9);
  T res = maxElementAbs(vector);
  EXPECT_EQ(res,-1.1);
}

TEST(TestVectorHelpers, maxElementAbsDoubleEarly)
{
  Vector<T,4> vector(2.003,0.2, -0.1, -0.9);
  T res = maxElementAbs(vector);
  EXPECT_EQ(res,2.003);
}

TEST(TestVectorHelpers, maxElementAbsIntLate)
{
  Vector<int,3> vector(-2, -4, -9);
  T res = maxElementAbs(vector);
  EXPECT_EQ(res,-9);
}

TEST(TestVectorHelpers, maxElementAbsIntMid)
{
  Vector<int,3> vector(2, 11, 9);
  T res = maxElementAbs(vector);
  EXPECT_EQ(res,11);
}

TEST(TestVectorHelpers, maxElementAbsIntEarly)
{
  Vector<int,4> vector(-12,2,-1,9);
  T res = maxElementAbs(vector);
  EXPECT_EQ(res,-12);
}

//maxElementAbsPos

TEST(TestVectorHelpers, maxElementAbsPosDoubleLate)
{
  Vector<T,3> vector(0.2, -0.4, 0.9);
  unsigned res = maxElementAbsPos(vector);
  EXPECT_EQ(res,2);
}

TEST(TestVectorHelpers, maxElementAbsPosDoubleMid)
{
  Vector<T,3> vector(0.2, -1.1, 0.9);
  unsigned res = maxElementAbsPos(vector);
  EXPECT_EQ(res,1);
}

TEST(TestVectorHelpers, maxElementAbsPosDoubleEarly)
{
  Vector<T,4> vector(2.003,0.2, -0.1, -0.9);
  unsigned res = maxElementAbsPos(vector);
  EXPECT_EQ(res,0);
}

TEST(TestVectorHelpers, maxElementAbsPosIntLate)
{
  Vector<int,3> vector(-2, -4, -9);
  unsigned res = maxElementAbsPos(vector);
  EXPECT_EQ(res,2);
}

TEST(TestVectorHelpers, maxElementAbsPosIntMid)
{
  Vector<int,3> vector(2, 11, 9);
  unsigned res = maxElementAbsPos(vector);
  EXPECT_EQ(res,1);
}

TEST(TestVectorHelpers, maxElementAbsPosIntEarly)
{
  Vector<int,4> vector(-12,2,-1,9);
  unsigned res = maxElementAbsPos(vector);
  EXPECT_EQ(res,0);
}


TEST(TestVectorHelpers, Determinant2d)
{
  constexpr Vector<double, 2> V1(1.,2.);
  constexpr Vector<double, 2> V2(-1.,6.);

  constexpr auto det = util::determinant(V1, V2);
  EXPECT_EQ(det, 8);
}

TEST(TestVectorHelpers, Determinant3d)
{
  constexpr Vector<double, 3> V1(1.,2.,3.);
  constexpr Vector<double, 3> V2(-1.,6.,2.);
  constexpr Vector<double, 3> V3(-1.,2.,1.);

  constexpr auto det = util::determinant(V1, V2, V3);
  EXPECT_EQ(det, 12);
}

TEST(TestVectorHelpers, SolveLinear2d)
{
  constexpr Vector<double, 2> V1(1.,2.);
  constexpr Vector<double, 2> V2(-1.,6.);
  constexpr Vector<double, 2> rhs(-1.,30.);

  const auto x = util::solveLinearSystem(V1, V2, rhs);
  EXPECT_EQ(x[0], 3.);
  EXPECT_EQ(x[1], 4.);
}

TEST(TestVectorHelpers, SolveLinear3d)
{
  constexpr Vector<double, 3> V1(1.,2.,3.);
  constexpr Vector<double, 3> V2(-1.,6.,2.);
  constexpr Vector<double, 3> V3(3.,1.,1.);
  constexpr Vector<double, 3> rhs(14.,35.,22.);

  const auto x = util::solveLinearSystem(V1, V2, V3, rhs);
  EXPECT_EQ(x[0], 3.);
  EXPECT_EQ(x[1], 4.);
  EXPECT_EQ(x[2], 5.);
}

TEST(TestVectorHelpers, serialSymmetricTensorIndex2d)
{
  EXPECT_EQ(serialSymmetricTensorIndex<2>(0,0), 0);
  EXPECT_EQ(serialSymmetricTensorIndex<2>(0,1), 1);
  EXPECT_EQ(serialSymmetricTensorIndex<2>(1,0), 1);
  EXPECT_EQ(serialSymmetricTensorIndex<2>(1,1), 2);
}

TEST(TestVectorHelpers, serialSymmetricTensorIndex3d)
{
  EXPECT_EQ(serialSymmetricTensorIndex<3>(0,0), 0);
  EXPECT_EQ(serialSymmetricTensorIndex<3>(0,1), 1);
  EXPECT_EQ(serialSymmetricTensorIndex<3>(0,2), 2);
  EXPECT_EQ(serialSymmetricTensorIndex<3>(1,0), 1);
  EXPECT_EQ(serialSymmetricTensorIndex<3>(1,1), 3);
  EXPECT_EQ(serialSymmetricTensorIndex<3>(1,2), 4);
  EXPECT_EQ(serialSymmetricTensorIndex<3>(2,0), 2);
  EXPECT_EQ(serialSymmetricTensorIndex<3>(2,1), 4);
  EXPECT_EQ(serialSymmetricTensorIndex<3>(2,2), 5);
}
