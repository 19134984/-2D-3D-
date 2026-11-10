/*  This file is part of the OpenLB library
*
*  Copyright (C) 2021 Julius Jeßberger
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


#define AdWarnings

#include <iostream>
#include <gtest/gtest.h>
#include "helper/unitTestHelper.h"

#include <olb.h>

using namespace olb;
using namespace olb::util;

using S = double;

TEST(testIniDiagonal, testSingleEntry)
{
  constexpr int dim = 1;
  S x = 10.5;
  util::ADf<S,dim> f(x);
  EXPECT_EQ(f.d(0), 0);

  iniDiagonal(&f, dim);
  EXPECT_EQ(f.v(), x);
  EXPECT_EQ(f.d(0), 1);
}

TEST(testIniDiagonal, testTwoEntries)
{
  constexpr int dim = 3;
  ADf<S,dim> f[dim];
  for (int i = 0; i < dim; ++i) {
    f[i] = 2.0;
    EXPECT_EQ(f[i].v(), 2.0);
    for (int j = 0; j < dim; ++j) {
      EXPECT_EQ(f[i].d(j), 0);
    }
  }

  iniDiagonal(f, dim);
  for (int i = 0; i < dim; ++i) {
    EXPECT_EQ(f[i].v(), 2.0);
    for (int j = 0; j < dim; ++j) {
      bool r = (i == j);
      EXPECT_EQ(f[i].d(j), r);
    }
  }
}

TEST(testCopyDerivatives, testSixEntries)
{
  constexpr int dim = 3;
  constexpr int length = 2;
  ADf<S,dim> f[length];

  f[0] = 0.0;
  f[0].d(0) = 0.0;
  f[0].d(1) = 0.1;
  f[0].d(2) = 0.2;
  f[1] = 1.0;
  f[1].d(0) = 1.0;
  f[1].d(1) = 1.1;
  f[1].d(2) = 1.2;

  S target[dim*length];
  copyDerivatives(target, f, length);

  EXPECT_EQ(target[0], 0.0);
  EXPECT_EQ(target[1], 0.1);
  EXPECT_EQ(target[2], 0.2);
  EXPECT_EQ(target[3], 1.0);
  EXPECT_EQ(target[4], 1.1);
  EXPECT_EQ(target[5], 1.2);
}

TEST(testCopyAs, testDoubleFloat)
{
  std::vector<double> source ({0,1});
  auto res = util::copyAs <float,double> (source);

  // EXPECT_TRUE((std::is_same_v<decltype(res),float>));
  EXPECT_EQ(res.size(), 2);
  EXPECT_EQ(res[0], 0);
  EXPECT_EQ(res[1], 1);
}

TEST(testCopyAs, testDoubleAD)
{
  using S = double;
  using T = ADf<S,3>;
  std::vector<S> source ({0,1});
  auto res = util::copyAs <T,S> (source);

  // EXPECT_TRUE((std::is_same_v<decltype(res),T>));
  EXPECT_EQ(res.size(), 2);
  EXPECT_EQ(res[0], 0);
  EXPECT_EQ(res[1], 1);
  EXPECT_TRUE(res[0].hasZeroDerivative());
  EXPECT_TRUE(res[1].hasZeroDerivative());
}

TEST(testCopyAs, testADdouble)
{
  using S = ADf<double,3>;
  using T = double;
  std::vector<S> source ({0,1});
  auto res = util::copyAs <T,S> (source);

  // EXPECT_TRUE((std::is_same_v<decltype(res),T>));
  EXPECT_EQ(res.size(), 2);
  EXPECT_EQ(res[0], 0);
  EXPECT_EQ(res[1], 1);
}

TEST(testIniAD, testSingle)
{
  double x {1.1};
  auto y = iniAD(x);
  EXPECT_EQ(y.v(), 1.1);
  EXPECT_EQ(y.d(0), 1.0);
}

TEST(testIniAD, testArray)
{
  using S = double;
  constexpr unsigned d = 2;

  S x[2] {0,1};
  auto res = iniAD<d>(x);

  EXPECT_EQ(res[0].v(), 0);
  EXPECT_EQ(res[1].v(), 1);
  EXPECT_EQ(res[0].d(0), 1);
  EXPECT_EQ(res[0].d(1), 0);
  EXPECT_EQ(res[1].d(0), 0);
  EXPECT_EQ(res[1].d(1), 1);
}

TEST(testIniAD, testVector)
{
  using S = double;
  constexpr unsigned d = 2;

  std::vector<S> x {0,1};
  auto res = iniAD<d,S,util::StdVector>(x);

  EXPECT_EQ(res.size(), 2);
  EXPECT_EQ(res[0].v(), 0);
  EXPECT_EQ(res[1].v(), 1);
  EXPECT_EQ(res[0].d(0), 1);
  EXPECT_EQ(res[0].d(1), 0);
  EXPECT_EQ(res[1].d(0), 0);
  EXPECT_EQ(res[1].d(1), 1);
}

TEST(testDerivativeFAD_single, testSquare)
{
  S input {2.5};
  auto result = derivativeFAD ([](auto x){
    return x*x;
  }, input);
  EXPECT_EQ(result, 2.0*2.5);
}

TEST(testDerivativeFAD_return, testSquare)
{
  S input[1] {2.5};
  auto result = derivativeFAD<1>([](auto x){
    return x[0]*x[0];
  }, input);
  EXPECT_EQ(result[0], 2.0*2.5);
}

TEST(testDerivativeFAD_return, testSumSquare)
{
  constexpr unsigned dim = 2;
  S input[dim] {2.5, 6.6};
  auto result = derivativeFAD<dim>([](auto x){
    return x[0]*x[0] + x[0]*x[1];
  }, input);
  EXPECT_EQ(result[0], 2.0*2.5+6.6);
  EXPECT_EQ(result[1], 2.5);
}

TEST(testDerivativeFAD_inplace, testSquare)
{
  S input[1] {2.5};
  S output[1];
  derivativeFAD<1>([](auto y, auto x){
    y[0] = x[0]*x[0];
  }, output, input, 1);

  EXPECT_EQ(output[0], 2.0*2.5);
}

TEST(testDerivativeFAD_inplace, testSumSquare)
{
  constexpr unsigned dim = 2;
  S input[dim] {2.5, 6.6};
  S output[dim];
  derivativeFAD<dim>([](auto y, auto x){
    y[0] = x[0]*x[0] + x[0]*x[1];
  }, output, input, 1);

  EXPECT_EQ(output[0], 2.0*2.5+6.6);
  EXPECT_EQ(output[1], 2.5);
}

TEST(testDerivativeFAD_inplace, testSixEntries)
{
  constexpr unsigned dim = 2;
  constexpr unsigned comp = 3;
  S input[dim] {2.5, 6.6};
  S output[dim*comp];
  derivativeFAD<dim>([](auto y, auto x){
    y[0] = x[0]+x[1];
    y[1] = 3.0 * x[0] * x[1];
    y[2] = x[1] * x[1];
  }, output, input, comp);

  EXPECT_EQ(output[0], 1.0);
  EXPECT_EQ(output[1], 1.0);
  EXPECT_EQ(output[2], 3.0*6.6);
  EXPECT_EQ(output[3], 3.0*2.5);
  EXPECT_EQ(output[4], 0.0);
  EXPECT_EQ(output[5], 2.0*6.6);
}
