/*  This file is part of the OpenLB library
*
*  Copyright (C) 2022 Jan E. Marquardt, 2023 Julius Jessberger
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
using namespace olb::particles;

typedef double T;

TEST(TestMatrix, ConstructFromVector)
{
  Vector<T,6> vec(1,2,3,4,5,6);
  Matrix<T, 2, 3> matrix(vec);

  EXPECT_EQ(matrix[0][0], 1);
  EXPECT_EQ(matrix[0][1], 2);
  EXPECT_EQ(matrix[0][2], 3);
  EXPECT_EQ(matrix[1][0], 4);
  EXPECT_EQ(matrix[1][1], 5);
  EXPECT_EQ(matrix[1][2], 6);
}

TEST(TestMatrix, Transpose1)
{
  T data[2][2] = {
      {1, 2},
      {3, 4}
  };
  Matrix<T, 2, 2> matrix(data);
  Matrix<T, 2, 2> transpose = matrix.transpose();

  EXPECT_EQ(matrix[0][0], data[0][0]);
  EXPECT_EQ(matrix[1][0], data[1][0]);
  EXPECT_EQ(matrix[0][1], data[0][1]);
  EXPECT_EQ(matrix[1][1], data[1][1]);

  EXPECT_EQ(matrix[0][0], transpose[0][0]);
  EXPECT_EQ(matrix[1][0], transpose[0][1]);
  EXPECT_EQ(matrix[0][1], transpose[1][0]);
  EXPECT_EQ(matrix[1][1], transpose[1][1]);
}

TEST(TestMatrix, Transpose2)
{
  T data[3][2] = {
      {1, 2},
      {3, 4},
      {5, 6}
  };
  Matrix<T, 3, 2> matrix(data);
  Matrix<T, 2, 3> transpose = matrix.transpose();

  for(unsigned i=0; i<3; ++i) {
    for(unsigned j=0; j<2; ++j) {
      EXPECT_EQ(matrix[i][j], data[i][j]);
    }
  }

  EXPECT_EQ(matrix[0][0], transpose[0][0]);
  EXPECT_EQ(matrix[1][0], transpose[0][1]);
  EXPECT_EQ(matrix[2][0], transpose[0][2]);
  EXPECT_EQ(matrix[0][1], transpose[1][0]);
  EXPECT_EQ(matrix[1][1], transpose[1][1]);
  EXPECT_EQ(matrix[2][1], transpose[1][2]);
}

TEST(TestMatrix, MatrixVectorMultiplication)
{
  T data[3][3] = {
    {2, 3, -4},
    {11,8,7},
    {2,5,3}
  };
  const Matrix<T,3,3> matrix(data);
  const Vector<T,3> vec(3,7,5);

  const Vector<T,3> result = matrix * vec;
  EXPECT_EQ(result[0], 7);
  EXPECT_EQ(result[1], 124);
  EXPECT_EQ(result[2], 56);
}

TEST(TestMatrix, Determinant2d)
{
  constexpr T data[2][2] = {
    {1, -1},
    {2, 6}
  };
  constexpr Matrix<T,2,2> matrix(data);

  constexpr auto det = util::determinant(matrix);
  EXPECT_EQ(det, 8);
}

TEST(TestMatrix, Determinant3d)
{
  constexpr T data[3][3] = {
    {1, -1, -1},
    {2, 6, 2},
    {3, 2, 1}
  };
  constexpr Matrix<T,3,3> matrix(data);

  constexpr auto det = util::determinant(matrix);
  EXPECT_EQ(det, 12);
}

TEST(TestMatrix, SolveLinear2d)
{
  constexpr T data[2][2] = {
    {1,  2},
    {-1, 6}
  };
  constexpr Matrix<T,2,2> matrix(data);
  constexpr Vector<double, 2> rhs(-1.,30.);

  const auto x = util::solveLinearSystem(matrix, rhs);
  EXPECT_EQ(x[0], 3.);
  EXPECT_EQ(x[1], 4.);
}

TEST(TestMatrix, SolveLinear2dIll1)
{
  constexpr T data[2][2] = {
    {1, 2},
    {2, 4}
  };
  constexpr Matrix<T,2,2> matrix(data);
  constexpr Vector<double, 2> rhs(-1.,30.);

  const auto x = util::solveLinearSystem(matrix, rhs);
  EXPECT_TRUE(isinf(x[0]));
  EXPECT_TRUE(isinf(x[1]));
}

TEST(TestMatrix, SolveLinear2dIll2)
{
  constexpr T data[2][2] = {
    {1, 2},
    {2, 4}
  };
  constexpr Matrix<T,2,2> matrix(data);
  constexpr Vector<double, 2> rhs(4.,8.);

  const auto x = util::solveLinearSystem(matrix, rhs);
  EXPECT_TRUE(isnan(x[0]));
  EXPECT_TRUE(isnan(x[1]));
}

TEST(TestMatrix, SolveLinear3d)
{
  constexpr T data[3][3] = {
    {1, 2, 3},
    {-1, 6, 2},
    {3, 1, 1}
  };
  constexpr Matrix<T,3,3> matrix(data);
  constexpr Vector<T,3> rhs(14,35,22);

  const auto x = util::solveLinearSystem(matrix, rhs);
  EXPECT_EQ(x[0], 3.);
  EXPECT_EQ(x[1], 4.);
  EXPECT_EQ(x[2], 5.);
}

TEST(TestMatrix, SolveLinear3dIll1)
{
  constexpr T data[3][3] = {
    {1, 2, 3},
    {1, 0, 0},
    {0, 4, 6}
  };
  constexpr Matrix<T,3,3> matrix(data);
  constexpr Vector<T,3> rhs(14,35,22);

  const auto x = util::solveLinearSystem(matrix, rhs);
  EXPECT_TRUE(isinf(x[0]));
  EXPECT_TRUE(isinf(x[1]));
  EXPECT_TRUE(isinf(x[2]));
}

TEST(TestMatrix, SolveLinear3dIll2)
{
  constexpr T data[3][3] = {
    {1, 2, 3},
    {1, 0, 0},
    {0, 4, 6}
  };
  constexpr Matrix<T,3,3> matrix(data);
  constexpr Vector<T,3> rhs(2,-8,-12);

  const auto x = util::solveLinearSystem(matrix, rhs);
  EXPECT_TRUE(isnan(x[0]));
  EXPECT_TRUE(isnan(x[1]));
  EXPECT_TRUE(isnan(x[2]));
}
