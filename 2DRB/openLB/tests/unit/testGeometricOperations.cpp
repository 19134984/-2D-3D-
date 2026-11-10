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

Vector<T, 4> matrixMultiplication(Vector<T, 4> rotMat, Vector<T, 4> invRotMat)
{
  return Vector<T, 4>(
      // Row 1
      rotMat[0] * invRotMat[0] + rotMat[1] * invRotMat[2],
      rotMat[0] * invRotMat[1] + rotMat[1] * invRotMat[3],
      // Row 2
      rotMat[2] * invRotMat[0] + rotMat[3] * invRotMat[2],
      rotMat[2] * invRotMat[1] + rotMat[3] * invRotMat[3]);
}

Vector<T, 9> matrixMultiplication(Vector<T, 9> rotMat, Vector<T, 9> invRotMat)
{
  return Vector<T, 9>(
      // Row 1
      rotMat[0] * invRotMat[0] + rotMat[1] * invRotMat[3] +
          rotMat[2] * invRotMat[6],
      rotMat[0] * invRotMat[1] + rotMat[1] * invRotMat[4] +
          rotMat[2] * invRotMat[7],
      rotMat[0] * invRotMat[2] + rotMat[1] * invRotMat[5] +
          rotMat[2] * invRotMat[8],
      // Row 2
      rotMat[3] * invRotMat[0] + rotMat[4] * invRotMat[3] +
          rotMat[5] * invRotMat[6],
      rotMat[3] * invRotMat[1] + rotMat[4] * invRotMat[4] +
          rotMat[5] * invRotMat[7],
      rotMat[3] * invRotMat[2] + rotMat[4] * invRotMat[5] +
          rotMat[5] * invRotMat[8],
      // Row 2
      rotMat[6] * invRotMat[0] + rotMat[7] * invRotMat[3] +
          rotMat[8] * invRotMat[6],
      rotMat[6] * invRotMat[1] + rotMat[7] * invRotMat[4] +
          rotMat[8] * invRotMat[7],
      rotMat[6] * invRotMat[2] + rotMat[7] * invRotMat[5] +
          rotMat[8] * invRotMat[8]);
}

// Angle conversion
TEST(DegreeToRadian, Scalar)
{
  for (int i = -24; i <= 24; ++i) {
    T angle = i * 15;
    EXPECT_NEAR(M_PI * angle / 180., util::degreeToRadian(angle), 1e-11);
  }
}

TEST(DegreeToRadian, Vector2d)
{
  for (int i = -24; i <= 24; ++i) {
    Vector<T, 1> angle(i * 15);
    EXPECT_NEAR(M_PI * angle[0] / 180., util::degreeToRadian(angle)[0], 1e-11);
  }
}

TEST(DegreeToRadian, Vector3d)
{
  for (int i = -24; i <= 24; ++i) {
    Vector<T, 3> angle(i * 15);
    for (unsigned iDim = 0; iDim < 3; ++iDim) {
      EXPECT_NEAR(M_PI * angle[iDim] / 180., util::degreeToRadian(angle)[iDim],
                  1e-11);
    }
  }
}

TEST(RadianToDegree, Scalar)
{
  for (int i = -24; i <= 24; ++i) {
    T angle = i * 15;
    EXPECT_NEAR(180 * angle / M_PI, util::radianToDegree(angle), 1e-11);
  }
}

TEST(RadianToDegree, Vector2d)
{
  for (int i = -24; i <= 24; ++i) {
    Vector<T, 1> angle(i * 15);
    EXPECT_NEAR(180 * angle[0] / M_PI, util::radianToDegree(angle)[0], 1e-11);
  }
}

TEST(RadianToDegree, Vector3d)
{
  for (int i = -24; i <= 24; ++i) {
    Vector<T, 3> angle(i * 15);
    for (unsigned iDim = 0; iDim < 3; ++iDim) {
      EXPECT_NEAR(180 * angle[iDim] / M_PI, util::radianToDegree(angle)[iDim],
                  1e-11);
    }
  }
}

// Rotation Matrix
// 2D
TEST(RotationMatrix2D, NoRotation)
{
  Vector<T, 1> angle(0.);
  Vector<T, 4> rotMat = util::calculateRotationMatrix<T, 2>(angle);

  // row 1
  EXPECT_EQ(rotMat[0], 1);
  EXPECT_EQ(rotMat[1], 0);
  // row 2
  EXPECT_EQ(rotMat[2], 0);
  EXPECT_EQ(rotMat[3], 1);
}

TEST(RotationMatrix2D, Rotation1)
{
  Vector<T, 1> angle(0.5 * M_PI);
  Vector<T, 4> rotMat = util::calculateRotationMatrix<T, 2>(angle);

  // row 1
  EXPECT_NEAR(rotMat[0], 0, 1e-16);
  EXPECT_NEAR(rotMat[1], -1, 1e-16);
  // row 2
  EXPECT_NEAR(rotMat[2], 1, 1e-16);
  EXPECT_NEAR(rotMat[3], 0, 1e-16);
}

TEST(RotationMatrix2D, Rotation2)
{
  Vector<T, 1> angle(-0.5 * M_PI);
  Vector<T, 4> rotMat = util::calculateRotationMatrix<T, 2>(angle);

  // row 1
  EXPECT_NEAR(rotMat[0], 0, 1e-16);
  EXPECT_NEAR(rotMat[1], 1, 1e-16);
  // row 2
  EXPECT_NEAR(rotMat[2], -1, 1e-16);
  EXPECT_NEAR(rotMat[3], 0, 1e-16);
}

TEST(RotationMatrix2D, Rotation3)
{
  Vector<T, 1> angle(0.25 * M_PI);
  Vector<T, 4> rotMat = util::calculateRotationMatrix<T, 2>(angle);

  T const tmp = 0.5 * util::sqrt(2);
  // row 1
  EXPECT_NEAR(rotMat[0], tmp, 1e-15);
  EXPECT_NEAR(rotMat[1], -tmp, 1e-15);
  // row 2
  EXPECT_NEAR(rotMat[2], tmp, 1e-15);
  EXPECT_NEAR(rotMat[3], tmp, 1e-15);
}

TEST(RotationMatrix2D, Rotation4)
{
  Vector<T, 1> angle(-0.25 * M_PI);
  Vector<T, 4> rotMat = util::calculateRotationMatrix<T, 2>(angle);

  T const tmp = 0.5 * util::sqrt(2);
  // row 1
  EXPECT_NEAR(rotMat[0], tmp, 1e-15);
  EXPECT_NEAR(rotMat[1], tmp, 1e-15);
  // row 2
  EXPECT_NEAR(rotMat[2], -tmp, 1e-15);
  EXPECT_NEAR(rotMat[3], tmp, 1e-15);
}

TEST(InverseRotationMatrix2D, NoRotation)
{
  Vector<T, 1> angle(0.);
  Vector<T, 4> rotMat    = util::calculateRotationMatrix<T, 2>(angle);
  Vector<T, 4> invRotMat = util::invertRotationMatrix<T, 2>(rotMat);

  Vector<T, 4> product = matrixMultiplication(rotMat, invRotMat);

  // row 1
  EXPECT_EQ(product[0], 1);
  EXPECT_EQ(product[1], 0);
  // row 2
  EXPECT_EQ(product[2], 0);
  EXPECT_EQ(product[3], 1);
}

TEST(InverseRotationMatrix2D, Rotation1)
{
  Vector<T, 1> angle(13.87);
  Vector<T, 4> rotMat    = util::calculateRotationMatrix<T, 2>(angle);
  Vector<T, 4> invRotMat = util::invertRotationMatrix<T, 2>(rotMat);

  Vector<T, 4> product = matrixMultiplication(rotMat, invRotMat);

  // row 1
  EXPECT_NEAR(product[0], 1, 1e-15);
  EXPECT_NEAR(product[1], 0, 1e-15);
  // row 2
  EXPECT_NEAR(product[2], 0, 1e-15);
  EXPECT_NEAR(product[3], 1, 1e-15);
}

TEST(InverseRotationMatrix2D, Rotation2)
{
  Vector<T, 1> angle(-77.8);
  Vector<T, 4> rotMat    = util::calculateRotationMatrix<T, 2>(angle);
  Vector<T, 4> invRotMat = util::invertRotationMatrix<T, 2>(rotMat);

  Vector<T, 4> product = matrixMultiplication(rotMat, invRotMat);

  // row 1
  EXPECT_NEAR(product[0], 1, 1e-15);
  EXPECT_NEAR(product[1], 0, 1e-15);
  // row 2
  EXPECT_NEAR(product[2], 0, 1e-15);
  EXPECT_NEAR(product[3], 1, 1e-15);
}

TEST(CalculateInverseRotationMatrix2D, Rotation1)
{
  Vector<T, 1> angle(13.87);
  Vector<T, 4> rotMat    = util::calculateRotationMatrix<T, 2>(angle);
  Vector<T, 4> invRotMat = util::calculateInverseRotationMatrix<T, 2>(angle);

  Vector<T, 4> product = matrixMultiplication(rotMat, invRotMat);

  // row 1
  EXPECT_NEAR(product[0], 1, 1e-15);
  EXPECT_NEAR(product[1], 0, 1e-15);
  // row 2
  EXPECT_NEAR(product[2], 0, 1e-15);
  EXPECT_NEAR(product[3], 1, 1e-15);
}

TEST(CalculateInverseRotationMatrix2D, Rotation2)
{
  Vector<T, 1> angle(-77.8);
  Vector<T, 4> rotMat    = util::calculateRotationMatrix<T, 2>(angle);
  Vector<T, 4> invRotMat = util::calculateInverseRotationMatrix<T, 2>(angle);

  Vector<T, 4> product = matrixMultiplication(rotMat, invRotMat);

  // row 1
  EXPECT_NEAR(product[0], 1, 1e-15);
  EXPECT_NEAR(product[1], 0, 1e-15);
  // row 2
  EXPECT_NEAR(product[2], 0, 1e-15);
  EXPECT_NEAR(product[3], 1, 1e-15);
}

// 3D
TEST(RotationMatrix3D, NoRotation)
{
  Vector<T, 3> angle(0., 0., 0.);
  Vector<T, 9> rotMat = util::calculateRotationMatrix<T, 3>(angle);

  // row 1
  EXPECT_EQ(rotMat[0], 1);
  EXPECT_EQ(rotMat[1], 0);
  EXPECT_EQ(rotMat[2], 0);
  // row 2
  EXPECT_EQ(rotMat[3], 0);
  EXPECT_EQ(rotMat[4], 1);
  EXPECT_EQ(rotMat[5], 0);
  // row 3
  EXPECT_EQ(rotMat[6], 0);
  EXPECT_EQ(rotMat[7], 0);
  EXPECT_EQ(rotMat[8], 1);
}

TEST(RotationMatrix3D, Rotation1)
{
  Vector<T, 3> angle(0.5 * M_PI, 0.5 * M_PI, 0.5 * M_PI);
  Vector<T, 9> rotMat = util::calculateRotationMatrix<T, 3>(angle);

  // row 1
  EXPECT_NEAR(rotMat[0], 0, 1e-16);
  EXPECT_NEAR(rotMat[1], 0, 1e-16);
  EXPECT_NEAR(rotMat[2], 1, 1e-16);
  // row 2
  EXPECT_NEAR(rotMat[3], 0, 1e-16);
  EXPECT_NEAR(rotMat[4], 1, 1e-16);
  EXPECT_NEAR(rotMat[5], 0, 1e-16);
  // row 3
  EXPECT_NEAR(rotMat[6], -1, 1e-16);
  EXPECT_NEAR(rotMat[7], 0, 1e-16);
  EXPECT_NEAR(rotMat[8], 0, 1e-16);
}

TEST(RotationMatrix3D, Rotation2)
{
  Vector<T, 3> angle(-0.5 * M_PI, -0.5 * M_PI, -0.5 * M_PI);
  Vector<T, 9> rotMat = util::calculateRotationMatrix<T, 3>(angle);

  // row 1
  EXPECT_NEAR(rotMat[0], 0, 1e-15);
  EXPECT_NEAR(rotMat[1], 0, 1e-15);
  EXPECT_NEAR(rotMat[2], 1, 1e-15);
  // row 2
  EXPECT_NEAR(rotMat[3], 0, 1e-15);
  EXPECT_NEAR(rotMat[4], -1, 1e-15);
  EXPECT_NEAR(rotMat[5], 0, 1e-15);
  // row 3
  EXPECT_NEAR(rotMat[6], 1, 1e-15);
  EXPECT_NEAR(rotMat[7], 0, 1e-15);
  EXPECT_NEAR(rotMat[8], 0, 1e-15);
}

TEST(RotationMatrix3D, Rotation3)
{
  Vector<T, 3> angle(0.25 * M_PI, 0.25 * M_PI, 0.25 * M_PI);
  Vector<T, 9> rotMat = util::calculateRotationMatrix<T, 3>(angle);

  T const tmp = 0.5 * util::sqrt(2);
  // row 1
  EXPECT_NEAR(rotMat[0], tmp * tmp, 1e-15);
  EXPECT_NEAR(rotMat[1], tmp * tmp * tmp - tmp * tmp, 1e-15);
  EXPECT_NEAR(rotMat[2], tmp * tmp * tmp + tmp * tmp, 1e-15);
  // row 2
  EXPECT_NEAR(rotMat[3], tmp * tmp, 1e-15);
  EXPECT_NEAR(rotMat[4], tmp * tmp * tmp + tmp * tmp, 1e-15);
  EXPECT_NEAR(rotMat[5], tmp * tmp * tmp - tmp * tmp, 1e-15);
  // row 3
  EXPECT_NEAR(rotMat[6], -tmp, 1e-15);
  EXPECT_NEAR(rotMat[7], tmp * tmp, 1e-15);
  EXPECT_NEAR(rotMat[8], tmp * tmp, 1e-15);
}

TEST(RotationMatrix3D, Rotation4)
{
  Vector<T, 3> angle(-0.25 * M_PI, -0.25 * M_PI, -0.25 * M_PI);
  Vector<T, 9> rotMat = util::calculateRotationMatrix<T, 3>(angle);

  T const tmp = 0.5 * util::sqrt(2);
  // row 1
  EXPECT_NEAR(rotMat[0], tmp * tmp, 1e-15);
  EXPECT_NEAR(rotMat[1], -tmp * -tmp * tmp - tmp * -tmp, 1e-15);
  EXPECT_NEAR(rotMat[2], tmp * -tmp * tmp + (-tmp * -tmp), 1e-15);
  // row 2
  EXPECT_NEAR(rotMat[3], tmp * -tmp, 1e-15);
  EXPECT_NEAR(rotMat[4], -tmp * -tmp * -tmp + tmp * tmp, 1e-15);
  EXPECT_NEAR(rotMat[5], tmp * -tmp * -tmp - (-tmp * tmp), 1e-15);
  // row 3
  EXPECT_NEAR(rotMat[6], -(-tmp), 1e-15);
  EXPECT_NEAR(rotMat[7], -tmp * tmp, 1e-15);
  EXPECT_NEAR(rotMat[8], tmp * tmp, 1e-15);
}

TEST(InverseRotationMatrix3D, NoRotation)
{
  Vector<T, 3> angle(0., 0., 0.);
  Vector<T, 9> rotMat    = util::calculateRotationMatrix<T, 3>(angle);
  Vector<T, 9> invRotMat = util::invertRotationMatrix<T, 3>(rotMat);

  Vector<T, 9> product = matrixMultiplication(rotMat, invRotMat);

  // row 1
  EXPECT_EQ(product[0], 1);
  EXPECT_EQ(product[1], 0);
  EXPECT_EQ(product[2], 0);
  // row 2
  EXPECT_EQ(product[3], 0);
  EXPECT_EQ(product[4], 1);
  EXPECT_EQ(product[5], 0);
  // row 3
  EXPECT_EQ(product[6], 0);
  EXPECT_EQ(product[7], 0);
  EXPECT_EQ(product[8], 1);
}

TEST(InverseRotationMatrix3D, Rotation1)
{
  Vector<T, 3> angle(-10., 20., 1.);
  Vector<T, 9> rotMat    = util::calculateRotationMatrix<T, 3>(angle);
  Vector<T, 9> invRotMat = util::invertRotationMatrix<T, 3>(rotMat);

  Vector<T, 9> product = matrixMultiplication(rotMat, invRotMat);

  // row 1
  EXPECT_NEAR(product[0], 1, 1e-15);
  EXPECT_NEAR(product[1], 0, 1e-15);
  EXPECT_NEAR(product[2], 0, 1e-15);
  // row 2
  EXPECT_NEAR(product[3], 0, 1e-15);
  EXPECT_NEAR(product[4], 1, 1e-15);
  EXPECT_NEAR(product[5], 0, 1e-15);
  // row 3
  EXPECT_NEAR(product[6], 0, 1e-15);
  EXPECT_NEAR(product[7], 0, 1e-15);
  EXPECT_NEAR(product[8], 1, 1e-15);
}

TEST(InverseRotationMatrix3D, Rotation2)
{
  Vector<T, 3> angle(23.3, -11.1, -39.27);
  Vector<T, 9> rotMat    = util::calculateRotationMatrix<T, 3>(angle);
  Vector<T, 9> invRotMat = util::invertRotationMatrix<T, 3>(rotMat);

  Vector<T, 9> product = matrixMultiplication(rotMat, invRotMat);

  // row 1
  EXPECT_NEAR(product[0], 1, 1e-15);
  EXPECT_NEAR(product[1], 0, 1e-15);
  EXPECT_NEAR(product[2], 0, 1e-15);
  // row 2
  EXPECT_NEAR(product[3], 0, 1e-15);
  EXPECT_NEAR(product[4], 1, 1e-15);
  EXPECT_NEAR(product[5], 0, 1e-15);
  // row 3
  EXPECT_NEAR(product[6], 0, 1e-15);
  EXPECT_NEAR(product[7], 0, 1e-15);
  EXPECT_NEAR(product[8], 1, 1e-15);
}

TEST(CalculateInverseRotationMatrix3D, Rotation1)
{
  Vector<T, 3> angle(-10., 20., 1.);
  Vector<T, 9> rotMat    = util::calculateRotationMatrix<T, 3>(angle);
  Vector<T, 9> invRotMat = util::calculateInverseRotationMatrix<T, 3>(angle);

  Vector<T, 9> product = matrixMultiplication(rotMat, invRotMat);

  // row 1
  EXPECT_NEAR(product[0], 1, 1e-15);
  EXPECT_NEAR(product[1], 0, 1e-15);
  EXPECT_NEAR(product[2], 0, 1e-15);
  // row 2
  EXPECT_NEAR(product[3], 0, 1e-15);
  EXPECT_NEAR(product[4], 1, 1e-15);
  EXPECT_NEAR(product[5], 0, 1e-15);
  // row 3
  EXPECT_NEAR(product[6], 0, 1e-15);
  EXPECT_NEAR(product[7], 0, 1e-15);
  EXPECT_NEAR(product[8], 1, 1e-15);
}

TEST(CalculateInverseRotationMatrix3D, Rotation2)
{
  Vector<T, 3> angle(23.3, -11.1, -39.27);
  Vector<T, 9> rotMat    = util::calculateRotationMatrix<T, 3>(angle);
  Vector<T, 9> invRotMat = util::calculateInverseRotationMatrix<T, 3>(angle);

  Vector<T, 9> product = matrixMultiplication(rotMat, invRotMat);

  // row 1
  EXPECT_NEAR(product[0], 1, 1e-15);
  EXPECT_NEAR(product[1], 0, 1e-15);
  EXPECT_NEAR(product[2], 0, 1e-15);
  // row 2
  EXPECT_NEAR(product[3], 0, 1e-15);
  EXPECT_NEAR(product[4], 1, 1e-15);
  EXPECT_NEAR(product[5], 0, 1e-15);
  // row 3
  EXPECT_NEAR(product[6], 0, 1e-15);
  EXPECT_NEAR(product[7], 0, 1e-15);
  EXPECT_NEAR(product[8], 1, 1e-15);
}

// Rotation Execution
// 2D
TEST(RotationExecution2D, NoRotation)
{
  Vector<T, 1> angle(0.);
  Vector<T, 4> rotMat = util::calculateRotationMatrix<T, 2>(angle);
  Vector<T, 2> rotCenter(5., 5.);

  PhysR<T, 2> input(1.);
  PhysR<T, 2> output;

  output = util::executeRotation<T, 2, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_EQ(input[iDim], output[iDim] + rotCenter[iDim]);
  }

  output = util::executeRotation<T, 2, false>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_EQ(input[iDim], output[iDim]);
  }
}

TEST(RotationExecution2D, Rotation1)
{
  Vector<T, 1> angle(45.);
  Vector<T, 4> rotMat =
      util::calculateRotationMatrix<T, 2>(util::degreeToRadian(angle));
  Vector<T, 2> rotCenter(5., 5.);

  PhysR<T, 2> input(1., 0.);
  PhysR<T, 2> output;
  PhysR<T, 2> expectation(util::sqrt(0.5), util::sqrt(0.5));

  output =
      util::executeRotation<T, 2, true>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_NEAR(expectation[iDim], output[iDim], 1e-14);
  }

  output =
      util::executeRotation<T, 2, false>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_NEAR(expectation[iDim] + rotCenter[iDim], output[iDim], 1e-14);
  }
}

TEST(RotationExecution2D, Rotation2)
{
  Vector<T, 1> angle(360.);
  Vector<T, 4> rotMat =
      util::calculateRotationMatrix<T, 2>(util::degreeToRadian(angle));
  Vector<T, 2> rotCenter(5., 5.);

  PhysR<T, 2> input(23.1, 9.37);
  PhysR<T, 2> output;
  PhysR<T, 2> expectation(input);

  output =
      util::executeRotation<T, 2, true>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_NEAR(expectation[iDim], output[iDim], 1e-14);
  }

  output =
      util::executeRotation<T, 2, false>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_NEAR(expectation[iDim] + rotCenter[iDim], output[iDim], 1e-14);
  }
}

TEST(RotationExecution2D, Rotation3)
{
  Vector<T, 1> angle(-360.);
  Vector<T, 4> rotMat =
      util::calculateRotationMatrix<T, 2>(util::degreeToRadian(angle));
  Vector<T, 2> rotCenter(5., 5.);

  PhysR<T, 2> input(23.1, 9.37);
  PhysR<T, 2> output;
  PhysR<T, 2> expectation(input);

  output =
      util::executeRotation<T, 2, true>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_NEAR(expectation[iDim], output[iDim], 1e-14);
  }

  output =
      util::executeRotation<T, 2, false>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_NEAR(expectation[iDim] + rotCenter[iDim], output[iDim], 1e-14);
  }
}

TEST(RotationExecution2D, Rotation4)
{
  Vector<T, 1> angle(-45.);
  Vector<T, 4> rotMat =
      util::calculateRotationMatrix<T, 2>(util::degreeToRadian(angle));
  Vector<T, 2> rotCenter(7., 7.);

  PhysR<T, 2> input(2., 0.);
  PhysR<T, 2> output;
  PhysR<T, 2> expectation(util::sqrt(2), -util::sqrt(2));

  output =
      util::executeRotation<T, 2, true>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_NEAR(expectation[iDim], output[iDim], 1e-14);
  }

  output =
      util::executeRotation<T, 2, false>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_NEAR(expectation[iDim] + rotCenter[iDim], output[iDim], 1e-14);
  }
}

TEST(RotationExecution2D, Rotation5)
{
  Vector<T, 1> angle(-180);
  Vector<T, 4> rotMat =
      util::calculateRotationMatrix<T, 2>(util::degreeToRadian(angle));
  Vector<T, 2> rotCenter(7., 7.);

  PhysR<T, 2> input(0., 5.);
  PhysR<T, 2> output;
  PhysR<T, 2> expectation(0., -5.);

  output =
      util::executeRotation<T, 2, true>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_NEAR(expectation[iDim], output[iDim], 1e-14);
  }

  output =
      util::executeRotation<T, 2, false>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_NEAR(expectation[iDim] + rotCenter[iDim], output[iDim], 1e-14);
  }
}

TEST(RotationExecution2D, Rotation6)
{
  Vector<T, 1> angle(225);
  Vector<T, 4> rotMat =
      util::calculateRotationMatrix<T, 2>(util::degreeToRadian(angle));
  Vector<T, 2> rotCenter(7., 7.);

  PhysR<T, 2> input(0., 2.);
  PhysR<T, 2> output;
  PhysR<T, 2> expectation(util::sqrt(2), -util::sqrt(2));

  output =
      util::executeRotation<T, 2, true>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_NEAR(expectation[iDim], output[iDim], 1e-14);
  }

  output =
      util::executeRotation<T, 2, false>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_NEAR(expectation[iDim] + rotCenter[iDim], output[iDim], 1e-14);
  }
}

TEST(RotationExecution2D, Rotation7)
{
  Vector<T, 1> angle(135);
  Vector<T, 4> rotMat =
      util::calculateRotationMatrix<T, 2>(util::degreeToRadian(angle));
  Vector<T, 2> rotCenter(7., 7.);

  PhysR<T, 2> input(0., 2.);
  PhysR<T, 2> output;
  PhysR<T, 2> expectation(-util::sqrt(2), -util::sqrt(2));

  output =
      util::executeRotation<T, 2, true>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_NEAR(expectation[iDim], output[iDim], 1e-14);
  }

  output =
      util::executeRotation<T, 2, false>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_NEAR(expectation[iDim] + rotCenter[iDim], output[iDim], 1e-14);
  }
}

TEST(InvertRotationExecution2D, NoRotation)
{
  Vector<T, 1> angle(0.);
  Vector<T, 4> rotMat = util::calculateRotationMatrix<T, 2>(angle);
  Vector<T, 2> rotCenter(5., 5.);

  PhysR<T, 2> input(1.);
  PhysR<T, 2> output;

  output = util::executeRotation<T, 2, true>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 2, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }

  output = util::executeRotation<T, 2, false>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 2, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }
}

TEST(InvertRotationExecution2D, Rotation1)
{
  Vector<T, 1> angle(45.);
  Vector<T, 4> rotMat =
      util::calculateRotationMatrix<T, 2>(util::degreeToRadian(angle));
  Vector<T, 2> rotCenter(5., 5.);

  PhysR<T, 2> input(1., 0.);
  PhysR<T, 2> output;

  output = util::executeRotation<T, 2, true>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 2, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }

  output = util::executeRotation<T, 2, false>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 2, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }
}

TEST(InvertRotationExecution2D, Rotation2)
{
  Vector<T, 1> angle(360.);
  Vector<T, 4> rotMat =
      util::calculateRotationMatrix<T, 2>(util::degreeToRadian(angle));
  Vector<T, 2> rotCenter(5., 5.);

  PhysR<T, 2> input(23.1, 9.37);
  PhysR<T, 2> output;

  output = util::executeRotation<T, 2, true>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 2, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }

  output = util::executeRotation<T, 2, false>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 2, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }
}

TEST(InvertRotationExecution2D, Rotation3)
{
  Vector<T, 1> angle(-360.);
  Vector<T, 4> rotMat =
      util::calculateRotationMatrix<T, 2>(util::degreeToRadian(angle));
  Vector<T, 2> rotCenter(5., 5.);

  PhysR<T, 2> input(23.1, 9.37);
  PhysR<T, 2> output;

  output = util::executeRotation<T, 2, true>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 2, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }

  output = util::executeRotation<T, 2, false>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 2, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }
}

TEST(InvertRotationExecution2D, Rotation4)
{
  Vector<T, 1> angle(-45.);
  Vector<T, 4> rotMat =
      util::calculateRotationMatrix<T, 2>(util::degreeToRadian(angle));
  Vector<T, 2> rotCenter(7., 7.);

  PhysR<T, 2> input(2., 0.);
  PhysR<T, 2> output;

  output = util::executeRotation<T, 2, true>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 2, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }

  output = util::executeRotation<T, 2, false>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 2, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }
}

TEST(InvertRotationExecution2D, Rotation5)
{
  Vector<T, 1> angle(-180);
  Vector<T, 4> rotMat =
      util::calculateRotationMatrix<T, 2>(util::degreeToRadian(angle));
  Vector<T, 2> rotCenter(7., 7.);

  PhysR<T, 2> input(0., 5.);
  PhysR<T, 2> output;

  output = util::executeRotation<T, 2, true>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 2, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }

  output = util::executeRotation<T, 2, false>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 2, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }
}

TEST(InvertRotationExecution2D, Rotation6)
{
  Vector<T, 1> angle(225);
  Vector<T, 4> rotMat =
      util::calculateRotationMatrix<T, 2>(util::degreeToRadian(angle));
  Vector<T, 2> rotCenter(7., 7.);

  PhysR<T, 2> input(0., 2.);
  PhysR<T, 2> output;

  output = util::executeRotation<T, 2, true>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 2, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }

  output = util::executeRotation<T, 2, false>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 2, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }
}

TEST(InvertRotationExecution2D, Rotation7)
{
  Vector<T, 1> angle(135);
  Vector<T, 4> rotMat =
      util::calculateRotationMatrix<T, 2>(util::degreeToRadian(angle));
  Vector<T, 2> rotCenter(7., 7.);

  PhysR<T, 2> input(0., 2.);
  PhysR<T, 2> output;

  output = util::executeRotation<T, 2, true>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 2, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }

  output = util::executeRotation<T, 2, false>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 2, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 2; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }
}

// 3D
TEST(RotationExecution3D, NoRotation)
{
  Vector<T, 3> angle(0.);
  Vector<T, 9> rotMat = util::calculateRotationMatrix<T, 3>(angle);
  Vector<T, 3> rotCenter(5., 5., 5.);

  PhysR<T, 3> input(1.);
  PhysR<T, 3> output;

  output = util::executeRotation<T, 3, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_EQ(input[iDim], output[iDim] + rotCenter[iDim]);
  }

  output = util::executeRotation<T, 3, false>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_EQ(input[iDim], output[iDim]);
  }
}

TEST(RotationExecution3D, Rotation1)
{
  Vector<T, 3> angle(0., 0., 135);
  Vector<T, 9> rotMat =
      util::calculateRotationMatrix<T, 3>(util::degreeToRadian(angle));
  Vector<T, 3> rotCenter(7., 7., 7.);

  PhysR<T, 3> input(0., 2., 0.);
  PhysR<T, 3> output;
  PhysR<T, 3> expectation(-util::sqrt(2), -util::sqrt(2), 0.);

  output =
      util::executeRotation<T, 3, true>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_NEAR(expectation[iDim], output[iDim], 1e-14);
  }

  output =
      util::executeRotation<T, 3, false>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_NEAR(expectation[iDim] + rotCenter[iDim], output[iDim], 1e-14);
  }
}

TEST(RotationExecution3D, Rotation2)
{
  Vector<T, 3> angle(0., 0., 225);
  Vector<T, 9> rotMat =
      util::calculateRotationMatrix<T, 3>(util::degreeToRadian(angle));
  Vector<T, 3> rotCenter(7., 7., 7.);

  PhysR<T, 3> input(2., 0., 0.);
  PhysR<T, 3> output;
  PhysR<T, 3> expectation(-util::sqrt(2), -util::sqrt(2), 0.);

  output =
      util::executeRotation<T, 3, true>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_NEAR(expectation[iDim], output[iDim], 1e-14);
  }

  output =
      util::executeRotation<T, 3, false>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_NEAR(expectation[iDim] + rotCenter[iDim], output[iDim], 1e-14);
  }
}

TEST(RotationExecution3D, Rotation3)
{
  Vector<T, 3> angle(0., 45., 0.);
  Vector<T, 9> rotMat =
      util::calculateRotationMatrix<T, 3>(util::degreeToRadian(angle));
  Vector<T, 3> rotCenter(7., 7., 7.);

  PhysR<T, 3> input(2., 0., 0.);
  PhysR<T, 3> output;
  PhysR<T, 3> expectation(util::sqrt(2), 0., -util::sqrt(2));

  output =
      util::executeRotation<T, 3, true>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_NEAR(expectation[iDim], output[iDim], 1e-14);
  }

  output =
      util::executeRotation<T, 3, false>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_NEAR(expectation[iDim] + rotCenter[iDim], output[iDim], 1e-14);
  }
}

TEST(RotationExecution3D, Rotation4)
{
  Vector<T, 3> angle(0., -45., 0.);
  Vector<T, 9> rotMat =
      util::calculateRotationMatrix<T, 3>(util::degreeToRadian(angle));
  Vector<T, 3> rotCenter(7., 7., 7.);

  PhysR<T, 3> input(0., 0., 2.);
  PhysR<T, 3> output;
  PhysR<T, 3> expectation(-util::sqrt(2), 0., util::sqrt(2));

  output =
      util::executeRotation<T, 3, true>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_NEAR(expectation[iDim], output[iDim], 1e-14);
  }

  output =
      util::executeRotation<T, 3, false>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_NEAR(expectation[iDim] + rotCenter[iDim], output[iDim], 1e-14);
  }
}

TEST(RotationExecution3D, Rotation5)
{
  Vector<T, 3> angle(135., 0., 0.);
  Vector<T, 9> rotMat =
      util::calculateRotationMatrix<T, 3>(util::degreeToRadian(angle));
  Vector<T, 3> rotCenter(7., 7., 7.);

  PhysR<T, 3> input(0., 0., 2.);
  PhysR<T, 3> output;
  PhysR<T, 3> expectation(0., -util::sqrt(2), -util::sqrt(2));

  output =
      util::executeRotation<T, 3, true>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_NEAR(expectation[iDim], output[iDim], 1e-14);
  }

  output =
      util::executeRotation<T, 3, false>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_NEAR(expectation[iDim] + rotCenter[iDim], output[iDim], 1e-14);
  }
}

TEST(RotationExecution3D, Rotation6)
{
  Vector<T, 3> angle(225., 0., 0.);
  Vector<T, 9> rotMat =
      util::calculateRotationMatrix<T, 3>(util::degreeToRadian(angle));
  Vector<T, 3> rotCenter(7., 7., 7.);

  PhysR<T, 3> input(0., -2., 0.);
  PhysR<T, 3> output;
  PhysR<T, 3> expectation(0., util::sqrt(2), util::sqrt(2));

  output =
      util::executeRotation<T, 3, true>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_NEAR(expectation[iDim], output[iDim], 1e-14);
  }

  output =
      util::executeRotation<T, 3, false>(input + rotCenter, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_NEAR(expectation[iDim] + rotCenter[iDim], output[iDim], 1e-14);
  }
}

TEST(InvertRotationExecution3D, NoRotation)
{
  Vector<T, 3> angle(0.);
  Vector<T, 9> rotMat = util::calculateRotationMatrix<T, 3>(angle);
  Vector<T, 3> rotCenter(5., 5., 5.);

  PhysR<T, 3> input(1.);
  PhysR<T, 3> output;

  output = util::executeRotation<T, 3, true>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 3, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }

  output = util::executeRotation<T, 3, false>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 3, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }
}

TEST(InvertRotationExecution3D, Rotation1)
{
  Vector<T, 3> angle(0., 0., 135);
  Vector<T, 9> rotMat =
      util::calculateRotationMatrix<T, 3>(util::degreeToRadian(angle));
  Vector<T, 3> rotCenter(7., 7., 7.);

  PhysR<T, 3> input(0., 2., 0.);
  PhysR<T, 3> output;

  output = util::executeRotation<T, 3, true>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 3, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }

  output = util::executeRotation<T, 3, false>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 3, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }
}

TEST(InvertRotationExecution3D, Rotation2)
{
  Vector<T, 3> angle(0., 0., 225);
  Vector<T, 9> rotMat =
      util::calculateRotationMatrix<T, 3>(util::degreeToRadian(angle));
  Vector<T, 3> rotCenter(7., 7., 7.);

  PhysR<T, 3> input(2., 0., 0.);
  PhysR<T, 3> output;

  output = util::executeRotation<T, 3, true>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 3, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }

  output = util::executeRotation<T, 3, false>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 3, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }
}

TEST(InvertRotationExecution3D, Rotation3)
{
  Vector<T, 3> angle(0., 45., 0.);
  Vector<T, 9> rotMat =
      util::calculateRotationMatrix<T, 3>(util::degreeToRadian(angle));
  Vector<T, 3> rotCenter(7., 7., 7.);

  PhysR<T, 3> input(2., 0., 0.);
  PhysR<T, 3> output;

  output = util::executeRotation<T, 3, true>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 3, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }

  output = util::executeRotation<T, 3, false>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 3, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }
}

TEST(InvertRotationExecution3D, Rotation4)
{
  Vector<T, 3> angle(0., -45., 0.);
  Vector<T, 9> rotMat =
      util::calculateRotationMatrix<T, 3>(util::degreeToRadian(angle));
  Vector<T, 3> rotCenter(7., 7., 7.);

  PhysR<T, 3> input(0., 0., 2.);
  PhysR<T, 3> output;

  output = util::executeRotation<T, 3, true>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 3, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }

  output = util::executeRotation<T, 3, false>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 3, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }
}

TEST(InvertRotationExecution3D, Rotation5)
{
  Vector<T, 3> angle(135., 0., 0.);
  Vector<T, 9> rotMat =
      util::calculateRotationMatrix<T, 3>(util::degreeToRadian(angle));
  Vector<T, 3> rotCenter(7., 7., 7.);

  PhysR<T, 3> input(0., 0., 2.);
  PhysR<T, 3> output;

  output = util::executeRotation<T, 3, true>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 3, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }

  output = util::executeRotation<T, 3, false>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 3, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }
}

TEST(InvertRotationExecution3D, Rotation6)
{
  Vector<T, 3> angle(225., 0., 0.);
  Vector<T, 9> rotMat =
      util::calculateRotationMatrix<T, 3>(util::degreeToRadian(angle));
  Vector<T, 3> rotCenter(7., 7., 7.);

  PhysR<T, 3> input(0., -2., 0.);
  PhysR<T, 3> output;

  output = util::executeRotation<T, 3, true>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 3, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }

  output = util::executeRotation<T, 3, false>(input, rotMat, rotCenter);
  output = util::invertRotation<T, 3, true>(input, rotMat, rotCenter);
  for (unsigned iDim = 0; iDim < 3; ++iDim) {
    EXPECT_EQ(input[iDim], input[iDim]);
  }
}

