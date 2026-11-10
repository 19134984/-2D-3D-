/*  This file is part of the OpenLB library
*
*  Copyright (C) 2020 Clara Schragmann
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

using namespace olb;
using namespace olb::descriptors;
using namespace olb::util;

// AD Tests utilizing type double

TEST(ADVectorTest, testVectorInitialization)
{
  typedef ADf<double,1> T;

  T f1;
  f1 = 1.;
  f1.setDiffVariable(0);

  T f2;
  f2 = 2.;
  f2.setDiffVariable(0);

  T f3;
  f3 = 3.;
  f3.setDiffVariable(0);

  Vector<T, 3> V1(f1,f2,f3);
  ASSERT_EQ(1,V1[0].v());
  ASSERT_EQ(2,V1[1].v());
  ASSERT_EQ(3,V1[2].v());

  ASSERT_EQ(1.0, V1[0].d(0));
  ASSERT_EQ(1.0, V1[1].d(0));
  ASSERT_EQ(1.0, V1[2].d(0));
}

TEST(ADVectorTest, testVectorAddition)
{
  typedef ADf<double,1> T;

  T f1;
  f1 = 1.;
  f1.setDiffVariable(0);

  T f2;
  f2 = 2.;
  f2.setDiffVariable(0);

  T f3;
  f3 = 3.;
  f3.setDiffVariable(0);

  //addition of a vector
  Vector<T, 3> V1(f1,f2,f3);
  V1 += V1;
  ASSERT_EQ(2,V1[0].v());
  ASSERT_EQ(4,V1[1].v());
  ASSERT_EQ(6,V1[2].v());

  ASSERT_EQ(2.0, V1[0].d(0));
  ASSERT_EQ(2.0, V1[1].d(0));
  ASSERT_EQ(2.0, V1[2].d(0));

  //addition of a scalar
  Vector<T, 3> V2(f1,f2,f3);
  V2 += 3;
  ASSERT_EQ(4,V2[0].v());
  ASSERT_EQ(5,V2[1].v());
  ASSERT_EQ(6,V2[2].v());

  ASSERT_EQ(1.0, V2[0].d(0));
  ASSERT_EQ(1.0, V2[1].d(0));
  ASSERT_EQ(1.0, V2[2].d(0));
}

TEST(ADVectorTest, testVectorSubtraction)
{
  typedef ADf<double,1> T;

  T f1;
  f1 = 1.;
  f1.setDiffVariable(0);

  T f2;
  f2 = 2.;
  f2.setDiffVariable(0);

  T f3;
  f3 = 3.;
  f3.setDiffVariable(0);

  //subtraction of a vector
  Vector<T, 3> V1(f1,f2,f3);
  Vector<T, 3> V2(f3,f2,f1);

  V1 -= V2;
  ASSERT_EQ(-2,V1[0].v());
  ASSERT_EQ(0,V1[1].v());
  ASSERT_EQ(2,V1[2].v());

  ASSERT_EQ(0.0, V1[0].d(0));
  ASSERT_EQ(0.0, V1[1].d(0));
  ASSERT_EQ(0.0, V1[2].d(0));

  //subtraction of a scalar
  Vector<T, 3> V3(f1,f2,f3);
  V3 -= 1;
  ASSERT_EQ(0,V3[0].v());
  ASSERT_EQ(1,V3[1].v());
  ASSERT_EQ(2,V3[2].v());

  ASSERT_EQ(1.0, V3[0].d(0));
  ASSERT_EQ(1.0, V3[1].d(0));
  ASSERT_EQ(1.0, V3[2].d(0));
}

TEST(ADVectorTest, testVectorMultiplication)
{
  typedef ADf<double,1> T;

  T f1;
  f1 = 1.;
  f1.setDiffVariable(0);

  T f2;
  f2 = 2.;
  f2.setDiffVariable(0);

  T f3;
  f3 = 3.;
  f3.setDiffVariable(0);

  //multiplication of a vector
  Vector<T, 3> V1(f1,f2,f3);
  Vector<T, 3> V2(f3,f2,f1);

  V1 *= V2;
  ASSERT_EQ(3,V1[0].v());
  ASSERT_EQ(4,V1[1].v());
  ASSERT_EQ(3,V1[2].v());

  ASSERT_EQ(4.0, V1[0].d(0));
  ASSERT_EQ(4.0, V1[1].d(0));
  ASSERT_EQ(4.0, V1[2].d(0));

  //multiplication of a scalar
  Vector<T, 3> V3(f1,f2,f3);
  V3 *= 3;
  ASSERT_EQ(3,V3[0].v());
  ASSERT_EQ(6,V3[1].v());
  ASSERT_EQ(9,V3[2].v());

  ASSERT_EQ(3.0, V3[0].d(0));
  ASSERT_EQ(3.0, V3[1].d(0));
  ASSERT_EQ(3.0, V3[2].d(0));
}

TEST(ADVectorTest, testVectorDivision)
{
  typedef ADf<double,1> T;

  T f1;
  f1 = 1.;
  f1.setDiffVariable(0);

  T f2;
  f2 = 2.;
  f2.setDiffVariable(0);

  T f3;
  f3 = 3.;
  f3.setDiffVariable(0);

  //division of a scalar
  Vector<T, 3> V3(f1,f2,f3);
  V3 /= 2;

  ASSERT_EQ(0.5,V3[0].v());
  ASSERT_EQ(1,V3[1].v());
  ASSERT_EQ(1.5,V3[2].v());

  ASSERT_EQ(0.5, V3[0].d(0));
  ASSERT_EQ(0.5, V3[1].d(0));
  ASSERT_EQ(0.5, V3[2].d(0));

}

TEST(ADVectorTest, testVectorArithmetic)
{
  typedef ADf<double,1> T;

  T f1;
  f1 = 1.;
  f1.setDiffVariable(0);

  T f2;
  f2 = 2.;
  f2.setDiffVariable(0);

  T f3;
  f3 = 3.;
  f3.setDiffVariable(0);

  Vector<T, 3> V1(f1,f2,f3);
  Vector<T, 3> V2(f1,f2,f3);
  Vector<T, 3> V3(f3,f1,f2);
  Vector<T, 3> V4(-f1,f2,-f3);

  ASSERT_TRUE(V1==V2);
  ASSERT_FALSE(V1==V3);
  ASSERT_FALSE(closeToZero(V1));

  Vector<T, 3> V5(crossProduct3D(V1,V3));
  ASSERT_EQ(1,V5[0]);
  ASSERT_EQ(7,V5[1]);
  ASSERT_EQ(-5,V5[2]);

  Vector<T, 3> V6(maxv(V1, V3));
  ASSERT_EQ(3,V6[0]);
  ASSERT_EQ(2,V6[1]);
  ASSERT_EQ(3,V6[2]);

  //ASSERT_EQ(norm(V1),sqrt(14));
  //Vector<T, 3> V7(normalize(V1));
  //ASSERT_TRUE(norm(V7)==1.);

  /*V4 = olb::abs(V4);
  ASSERT_EQ(1,V4[0]);
  ASSERT_EQ(2,V4[1]);
  ASSERT_EQ(3,V4[2]);*/
}

TEST(ADVectorTest, testVectorComparison)
{
  typedef ADf<double,1> T;

  T f1;
  f1 = 1.;
  f1.setDiffVariable(0);

  T f2;
  f2 = 2.;
  f2.setDiffVariable(0);

  T f3;
  f3 = 3.;
  f3.setDiffVariable(0);

  Vector<T, 3> V1(f1,f2,f3);
  Vector<T, 3> V2(f1,f2,f3);
  Vector<T, 3> V3(f3,f1,f2);

  ASSERT_TRUE(V1==V2);
  ASSERT_FALSE(V1!=V2);
  ASSERT_FALSE(V1==V3);
  ASSERT_TRUE(V1!=V3);

  Vector<T, 3> V4(f3,f3,f3);

  ASSERT_TRUE(V4>=V1);
  ASSERT_FALSE(V4>V1);
  ASSERT_TRUE(V1<=V4);
  ASSERT_FALSE(V1<V4);
}
