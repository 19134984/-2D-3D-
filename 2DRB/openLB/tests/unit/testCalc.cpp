/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Julius Jessberger
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
using namespace olb::util;

TEST(TestCalc, testFmodPos)
{
  double a = 1.5;
  long double b = 0.4;
  ADf<double,1> c = 1.5;
  ADf<double,1> d = 0.4;
  float e1 = -1.4;
  double e2 = -1.4;
  ADf<double,1> f = -1.7;

  EXPECT_NEAR(util::fmod_pos(a,b), 0.3, 1.e-15);
  EXPECT_NEAR(util::fmod_pos(a,c), 0.0, 1.e-15);
  EXPECT_NEAR(util::fmod_pos(b,c), 0.4, 1.e-15);
  EXPECT_NEAR(util::fmod_pos(c,b).v(), 0.3, 1.e-15);
  EXPECT_NEAR(util::fmod_pos(c,d).v(), 0.3, 1.e-15);

  EXPECT_NEAR(util::fmod_pos(e1,b), 0.2, 1.e-7);
  EXPECT_NEAR(util::fmod_pos(e2,b), 0.2, 1.e-15);
  EXPECT_NEAR(util::fmod_pos(f,a), 1.3, 1.e-15);
}

TEST(TestKahan, testFloat)
{
  using T = float;
  const T s {static_cast<T>(util::sqrt(2))};
  const int N {1000000};

  KahanSummator<T> k;
  for (int i = 0; i < N; ++i) {
    k.add(s);
  }
  EXPECT_EQ(k.getSum(), T(N) * s);
}

TEST(TestKahan, testDouble)
{
  using T = double;
  const T s {util::sqrt(2)};
  const int N {1000000};

  KahanSummator<T> k;
  for (int i = 0; i < N; ++i) {
    k.add(s);
  }
  EXPECT_EQ(k.getSum(), T(N) * s);
}

TEST(TestKahan, testLongDouble)
{
  using T = long double;
  const T s {util::sqrt(2)};
  const int N {1000000};

  KahanSummator<T> k;
  for (int i = 0; i < N; ++i) {
    k.add(s);
  }
  EXPECT_EQ(k.getSum(), T(N) * s);
}

TEST(TestCrossSum, testPositive)
{
  using K = unsigned;
  const K k {12345};
  auto r = crossSum(k);
  EXPECT_EQ(r, 1+2+3+4+5);
}

TEST(TestCrossSum, testZero)
{
  using K = unsigned;
  const K k {0};
  auto r = crossSum(k);
  EXPECT_EQ(r, 0);
}

TEST(TestCrossSum, testNegative)
{
  using K = int;
  const K k {-123};
  auto r = crossSum(k);
  EXPECT_EQ(r, 1+2+3);
}
