/*  This file is part of the OpenLB library
*
*  Copyright (C) 2020 Clara Schragmann, 2023 Julius Jessberger
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
#include "helper/unitTestHelper.h"

#include <olb.h>

using namespace olb;
typedef int T;


TEST(testVector, VectorInitialization)
{
  Vector<T, 3> V1(1,2,3);
  EXPECT_EQ(1,V1[0]);
  EXPECT_EQ(2,V1[1]);
  EXPECT_EQ(3,V1[2]);

  EXPECT_EQ(1,V1.data()[0]);
  EXPECT_EQ(2,V1.data()[1]);
  EXPECT_EQ(3,V1.data()[2]);

  Vector<T, 2> V2(1,2);
  EXPECT_EQ(1,V2[0]);
  EXPECT_EQ(2,V2[1]);

  Vector<T, 3> V3(V1);
  EXPECT_EQ(1,V3[0]);
  EXPECT_EQ(2,V3[1]);
  EXPECT_EQ(3,V3[2]);

  T array[3] {1,2,3};
  Vector<T, 3> V4(array);
  EXPECT_EQ(1,V4[0]);
  EXPECT_EQ(2,V4[1]);
  EXPECT_EQ(3,V4[2]);

  std::vector<T> stdV {1,2,3};
  Vector<T, 3> V5(stdV);
  EXPECT_EQ(1,V5[0]);
  EXPECT_EQ(2,V5[1]);
  EXPECT_EQ(3,V5[2]);

  Vector<T, 3> V6(1);
  EXPECT_EQ(1,V6[0]);
  EXPECT_EQ(1,V6[1]);
  EXPECT_EQ(1,V6[2]);

  Vector<T, 3> V7 = V1;
  EXPECT_EQ(1,V7[0]);
  EXPECT_EQ(2,V7[1]);
  EXPECT_EQ(3,V7[2]);

  Vector<T, 3> V8 = {1,2,3};
  EXPECT_EQ(1,V8[0]);
  EXPECT_EQ(2,V8[1]);
  EXPECT_EQ(3,V8[2]);
}


TEST(testVector, VectorComparison)
{
  Vector<T, 3> V1(1,2,3);
  Vector<T, 3> V2(1,2,3);
  Vector<T, 3> V3(2,3,4);

  EXPECT_TRUE(V2==V1);
  EXPECT_FALSE(V3==V1);

  EXPECT_FALSE(V2!=V1);
  EXPECT_TRUE(V3!=V1);

  EXPECT_FALSE(closeToZero(V1));

  Vector<double, 3> V5(2,3,4);
  V5 = normalize(V5, 1.);
  EXPECT_TRUE(norm(V5)==1);
  Vector<double, 3> V6(2,3,4);
  V6 = normalize(V6);
  EXPECT_TRUE(norm_squared(V6)==1);

  Vector<T, 3> V7(crossProduct3D(V2,V3));
  EXPECT_EQ(-1,V7[0]);
  EXPECT_EQ(2,V7[1]);
  EXPECT_EQ(-1,V7[2]);

  Vector<T, 3> V8(-1,2,-3);
  V8 = olb::abs(V8);
  EXPECT_EQ(1,V8[0]);
  EXPECT_EQ(2,V8[1]);
  EXPECT_EQ(3,V8[2]);

  Vector<T, 3> V9(2,3,4);
  std::vector<T> V10 = toStdVector(V9);
  EXPECT_EQ(2,V10[0]);
  EXPECT_EQ(3,V10[1]);
  EXPECT_EQ(4,V10[2]);

  Vector<double, 3> V11(1,2,3);
  EXPECT_EQ(norm(V11),std::sqrt(14));
  EXPECT_EQ(norm_squared(V11),14);
}

TEST(testVector, VectorAddition)
{
  Vector<T, 3> V1(1,2,3);

  Vector<T, 3> V2(2,3,4);
  V2+=V1;
  EXPECT_EQ(3,V2[0]);
  EXPECT_EQ(5,V2[1]);
  EXPECT_EQ(7,V2[2]);

  Vector<T, 3> V3(2,3,4);
  V3+=2;
  EXPECT_EQ(4,V3[0]);
  EXPECT_EQ(5,V3[1]);
  EXPECT_EQ(6,V3[2]);

  Vector<T, 3> V4(2,3,4);
  Vector<T, 3> V5(3,5,7);

  EXPECT_TRUE((V1+1)==V4);
  EXPECT_TRUE((1+V1)==V4);
  EXPECT_TRUE((V1+V4)==V5);

  Vector<double, 3> V6(0.1,0.2,0.3);
  EXPECT_EQ((V1+V6)[0],1.1);
  EXPECT_EQ((V1+V6)[1],2.2);
  EXPECT_EQ((V1+V6)[2],3.3);

  Vector<double, 3> V7(2.5,3.5,4.5);
  Vector<decltype(T{}+double{}), 3> V8 = V1 + V7;
  EXPECT_EQ(3.5,V8[0]);
  EXPECT_EQ(5.5,V8[1]);
  EXPECT_EQ(7.5,V8[2]);
}

TEST(testVector, VectorSubtraction)
{
  Vector<T, 3> V1(1,2,3);

  Vector<T, 3> V4(2,3,4);
  V4-=V1;
  EXPECT_EQ(1,V4[0]);
  EXPECT_EQ(1,V4[1]);
  EXPECT_EQ(1,V4[2]);

  Vector<T, 3> V5(2,3,4);
  V5-=2;
  EXPECT_EQ(0,V5[0]);
  EXPECT_EQ(1,V5[1]);
  EXPECT_EQ(2,V5[2]);

  Vector<T, 3> V6(0,1,2);
  Vector<T, 3> V7(0,-1,-2);
  Vector<T, 3> V8(1,1,1);

  EXPECT_TRUE((V1-1)==V6);
  EXPECT_TRUE((1-V1)==V7);
  EXPECT_TRUE((V1-V6)==V8);

  Vector<double, 3> V9(0.1,0.2,0.3);
  EXPECT_EQ((V1-V9)[0],0.9);
  EXPECT_EQ((V1-V9)[1],1.8);
  EXPECT_EQ((V1-V9)[2],2.7);

  Vector<double, 3> V10(0.5,1.5,2.5);
  Vector<decltype(T{}-double{}), 3> V11 = V1 - V10;
  EXPECT_EQ(0.5,V11[0]);
  EXPECT_EQ(0.5,V11[1]);
  EXPECT_EQ(0.5,V11[2]);
}

TEST(testVector, MultiplyAndDivide)
{
  Vector<T, 3> V1(1,2,3);
  Vector<T, 3> V2(2,4,6);

  EXPECT_TRUE((V1*2)==V2);
  EXPECT_TRUE((2*V1)==V2);
  EXPECT_EQ((V1*V2),28);

  EXPECT_TRUE((V2/2)==V1);

  Vector<T, 3> V3(2,3,4);
  V3*=2;
  EXPECT_EQ(4,V3[0]);
  EXPECT_EQ(6,V3[1]);
  EXPECT_EQ(8,V3[2]);

  Vector<T, 3> V4(4,6,4);
  V4/=2;
  EXPECT_EQ(2,V4[0]);
  EXPECT_EQ(3,V4[1]);
  EXPECT_EQ(2,V4[2]);
}

TEST(testVector, ScalarVectorMultiply)
{
  {
    Vector<std::int32_t,3> lhs{5,5,5};
    std::int32_t rhs(2);
    auto result = lhs * rhs;
    EXPECT_EQ(result[0], 10);
    EXPECT_EQ(result[1], 10);
    EXPECT_EQ(result[2], 10);
  }

  {
    std::int32_t lhs(2);
    Vector<std::int32_t,3> rhs{5,5,5};
    auto result = lhs * rhs;
    EXPECT_EQ(result[0], 10);
    EXPECT_EQ(result[1], 10);
    EXPECT_EQ(result[2], 10);
  }
}

TEST(testVector, MixedTypeScalarVectorMultiply)
{
  {
    Vector<std::int32_t,3> lhs{5,5,5};
    double rhs(0.9);
    auto result = lhs * rhs;
    EXPECT_EQ(result[0], 4.5);
    EXPECT_EQ(result[1], 4.5);
    EXPECT_EQ(result[2], 4.5);
  }

  {
    double lhs(0.9);
    Vector<std::int32_t,3> rhs{5,5,5};
    auto result = lhs * rhs;
    EXPECT_EQ(result[0], 4.5);
    EXPECT_EQ(result[1], 4.5);
    EXPECT_EQ(result[2], 4.5);
  }
}

TEST(testVector, SizeComparison)
{
  Vector<T, 3> V1(1,2,3);
  Vector<T, 3> V2(2,4,6);
  Vector<T, 3> V3(2,5,7);

  EXPECT_TRUE(V1<V2);
  EXPECT_TRUE(V2>V1);
  EXPECT_TRUE(V1<=V2);
  EXPECT_TRUE(V2>=V1);
  EXPECT_FALSE(V3>V2);
  EXPECT_TRUE(V3>=V2);
  EXPECT_FALSE(V2<V3);
  EXPECT_TRUE(V2<=V3);

  EXPECT_TRUE((maxv(V1,V2))==V2);
}

TEST(testVector, ComplexArithmetics)
{
  Vector<double, 3> V1(1.,2.,3.);
  Vector<double, 3> V2(-1.,6.,9.);

  EXPECT_EQ(0.,V1*(crossProduct3D(V1, V2)));
  EXPECT_EQ(0.,V2*(crossProduct3D(V1, V2)));

  Vector<double, 3> V3 = ((V1+V2)/2+5)*3;
  EXPECT_EQ(V3[0],15);
  EXPECT_EQ(V3[1],27);
  EXPECT_EQ(V3[2],33);

  Vector<double, 3> V4 = ((olb::abs(V1-V2)-5.)*2);
  EXPECT_EQ(V4[0],-6);
  EXPECT_EQ(V4[1],-2);
  EXPECT_EQ(V4[2],2);

  Vector<double, 3> V5 = norm(V1+V2-7.)*V1;
  EXPECT_EQ(V5[0],std::sqrt(75));
  EXPECT_EQ(V5[1],2*std::sqrt(75));
  EXPECT_EQ(V5[2],3*std::sqrt(75));
}

TEST(testVector, Angles)
{
  Vector<double, 3> V1(1.,2.,3.);
  Vector<double, 3> V2(-1.,6.,9.);

  Vector<double, 3> angles = util::angleBetweenVectors(V1, V2);
  EXPECT_EQ(angles[0]*180/M_PI, 0);
  EXPECT_NEAR(angles[1]*180/M_PI, 24.775140568832, 1e-12);
  EXPECT_NEAR(angles[2]*180/M_PI, 36.027373385104, 1e-12);
}
