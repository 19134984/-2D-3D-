/*  This file is part of the OpenLB library
*
*  Copyright (C) 2020 Jan E. Marquardt, Luca Heim
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
using namespace olb::descriptors;

constexpr unsigned DIM2 = 2;

TEST(ADTest, testTypeTraits)
{
  struct S { };

  EXPECT_FALSE( (util::is_adf<double>::value) );
  EXPECT_FALSE( (util::is_adf_v<double>) );
  EXPECT_FALSE( (util::is_adf_v<int>) );
  EXPECT_FALSE( (util::is_adf_v<S>) );

  EXPECT_TRUE( (util::is_adf_v<util::ADf<double,1>>) );
  EXPECT_TRUE( (util::is_adf_v<util::ADf<util::ADf<float,0>,1>>) );
  EXPECT_TRUE( (util::is_adf_v<util::ADf<S,1>>) );
}

// AD Tests utilizing type double

TEST(doubleADTest, testConstructFromNoArgument)
{
  typedef double S;
  S x = 10.5;
  util::ADf<S,DIM2> f;
  f = x;
  f.setDiffVariable(0);
  EXPECT_EQ(f.v(), 10.5);
  EXPECT_EQ(f.d(0), 1.0);
}

TEST(doubleADTest, testConstructFromADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f;
  f = x;
  f.setDiffVariable(0);
  T g(f);
  EXPECT_EQ(g.v(), 10.5);
  EXPECT_EQ(g.d(0), 1.0);
}

TEST(doubleADTest, testConstructFromDouble)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f(x);
  f.setDiffVariable(0);
  EXPECT_EQ(f.v(), 10.5);
  EXPECT_EQ(f.d(0), 1.0);
}

TEST(doubleADTest, testConstructFromDoubleWithDerivative)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  Vector<S,2> d = {5,0};
  T f(x,d);
  EXPECT_EQ(f.v(), 10.5);
  EXPECT_EQ(f.d(0), 5);
}

TEST(doubleADTest, testConstructFromDoubleWithMultiDimDerivative)
{
  typedef double S;
  typedef util::ADf<S,5> T;
  S x = 10.5;
  Vector<S,5> d = {5,2,8,1.1,-100};
  T f(x,d);
  EXPECT_EQ(f.v(), 10.5);
  EXPECT_EQ(f.d(), d);
}

TEST(doubleADTest, testSingleDimConst)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[1];
  f[0] = x;
  f[0].setDiffVariable(0);
  EXPECT_EQ(f[0].v(), 10.5);
  EXPECT_EQ(f[0].d(0), 1.0);
}

TEST(doubleADTest, testSingleDimAdditionADfdouble)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = f[0] + 9.25;
  EXPECT_EQ(f[1].v(), 19.75);
  EXPECT_EQ(f[1].d(0), 1.0);
}

TEST(doubleADTest, testSingleDimAdditionADfADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = f[0]+f[0];
  EXPECT_EQ(f[1].v(), 21.0);
  EXPECT_EQ(f[1].d(0), 2.0);
}

TEST(doubleADTest, testSingleDimAdditionDoubleADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] =  9.25 + f[0];
  EXPECT_EQ(f[1].v(), 19.75);
  EXPECT_EQ(f[1].d(0), 1.0);
}

TEST(doubleADTest, testSingleDimAdditionDoubleDouble)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T f = 2.5 + 1.2;
  EXPECT_EQ(f.v(), 3.7);
  EXPECT_EQ(f.d(0), 0.0);
}

TEST(doubleADTest, testSingleDimSubtraktionADfdouble)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = f[0] - 19.25;
  EXPECT_EQ(f[1].v(), -8.75);
  EXPECT_EQ(f[1].d(0), 1.0);
}

TEST(doubleADTest, testSingleDimSubtraktionADfADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = f[0] - f[0];
  EXPECT_EQ(f[1].v(), 0.0);
  EXPECT_EQ(f[1].d(0), 0.0);
}

TEST(doubleADTest, testSingleDimSubtraktionDoubleADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] =  9.25 - f[0];
  EXPECT_EQ(f[1].v(), -1.25);
  EXPECT_EQ(f[1].d(0), -1.0);
}

TEST(doubleADTest, testSingleDimSubtraktionDoubleDouble)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T f = 2.5 - 1.2;
  EXPECT_EQ(f.v(), 1.3);
  EXPECT_EQ(f.d(0), 0.0);
}

TEST(doubleADTest, testSingleDimMultiplicationADfdouble)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = f[0] * 9.25;
  EXPECT_EQ(f[1].v(), 97.125);
  EXPECT_EQ(f[1].d(0), 9.25);
}

TEST(doubleADTest, testSingleDimMultiplicationADfADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = f[0]*f[0];
  EXPECT_EQ(f[1].v(), 110.25);
  EXPECT_EQ(f[1].d(0), 21.0);
}

TEST(doubleADTest, testSingleDimMultiplicationDoubleADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] =  9.25 * +f[0];
  EXPECT_EQ(f[1].v(), +97.125);
  EXPECT_EQ(f[1].d(0), +9.25);
}

TEST(doubleADTest, testSingleDimMultiplicationDoubleDouble)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T f = 2.5 * 1.2;
  EXPECT_EQ(f.v(), 3.0);
  EXPECT_EQ(f.d(0), 0.0);
}

TEST(doubleADTest, testSingleDimDivisionDoubleADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 2.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = 1.0/f[0];
  EXPECT_NEAR(f[1].v(), 0.4, 1e-12);
  EXPECT_NEAR(f[1].d(0), -0.16, 1e-12);
}

TEST(doubleADTest, testSingleDimDivisionADfADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = f[0]/f[0];
  EXPECT_EQ(f[1].v(), 1.0);
  EXPECT_EQ(f[1].d(0), 0.0);
}

TEST(doubleADTest, testSingleDimDivisionADfDouble)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 2.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = f[0]/-5.0;
  EXPECT_NEAR(f[1].v(), -0.5, 1e-12);
  EXPECT_NEAR(f[1].d(0), -0.2, 1e-12);
}

TEST(doubleADTest, testSingleDimDivisionDoubleDouble)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T f = 3.4 / -1.6;
  EXPECT_EQ(f.v(), -2.125);
  EXPECT_EQ(f.d(0), 0.0);
}

TEST(doubleADTest, testSingleDimPowADfint)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 2.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::pow(f[0],-1);
  EXPECT_NEAR(f[1].v(), 0.4, 1e-15);
  EXPECT_EQ(f[1].d(0), -0.16);
}

TEST(doubleADTest, testSingleDimPowADfdouble)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 2.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::pow(f[0],-1.0);
  EXPECT_NEAR(f[1].v(), 0.4, 1e-15);
  EXPECT_EQ(f[1].d(0), -0.16);
}

TEST(doubleADTest, testSingleDimPowdoubleADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = -1;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::pow(2.5, f[0]);
  EXPECT_NEAR(f[1].v(), 0.4, 1e-15);
  EXPECT_NEAR(f[1].d(0), 0.366516292749662026073410884707204428580040487963304987116, 1e-15);
}

TEST(doubleADTest, testSingleDimPowdoubleZeroBaseADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 0;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::pow(f[0], 5.5);
  EXPECT_NEAR(f[1].v(), 0, 1e-15);
  EXPECT_NEAR(f[1].d(0), 0, 1e-15);
}

TEST(doubleADTest, testSingleDimPowdoubleNegativeBaseADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = -2.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::pow(f[0], 3.0);
  EXPECT_NEAR(f[1].v(), -15.625, 1e-15);
  EXPECT_NEAR(f[1].d(0), 18.75, 1e-15);
}

TEST(doubleADTest, testSingleDimPowAdfAdf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 3.0;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::pow(f[0], f[0]);
  EXPECT_EQ(f[1].v(), 27.0);
  EXPECT_NEAR(f[1].d(0), 56.6625317940389617, 1e-12);
}

TEST(doubleADTest, testMultiDimPowAdfAdf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 3.0;
  S y = 1.77;
  T f[3];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = y;
  f[1].setDiffVariable(1);
  f[2] = util::pow(f[0], f[1]);
  EXPECT_NEAR(f[2].v(), 6.990441744888698651441358430347227855484708956073113669125, 1e-12);
  EXPECT_NEAR(f[2].d(0), 4.124360629484332204350401473904864434735978284083137064784, 1e-12);
  EXPECT_NEAR(f[2].d(1), 7.679785204153267407696533372050107003013422483333443401365, 1e-12);
}

TEST(doubleADTest, testSingleDimExp)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::exp(f[0]);
  EXPECT_NEAR(f[1].v(), 36315.502674246636, 1e-12);
  EXPECT_NEAR(f[1].d(0), 36315.502674246636, 1e-12);
}

TEST(doubleADTest, testSingleDimSqr)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::sqr(f[0]);
  EXPECT_EQ(f[1].v(), 110.25);
  EXPECT_EQ(f[1].d(0), 21.0);
}

TEST(doubleADTest, testSingleDimLog)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 5.25;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::log(f[0]);
  EXPECT_NEAR(f[1].v(), 1.65822807660353, 1e-12);
  EXPECT_NEAR(f[1].d(0), 0.1904761904761905, 1e-12);
}

TEST(doubleADTest, testSingleDimLog10)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 53.7;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::log10(f[0]);
  EXPECT_NEAR(f[1].v(), 1.7299742856995556068654105950281739286312009548692244, 1e-12);
  EXPECT_NEAR(f[1].d(0), 0.0080874205196136280754400171120410629849980820447610, 1e-12);
}

TEST(doubleADTest, testSingleDimLog2)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 94.6;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::log2(f[0]);
  EXPECT_NEAR(f[1].v(), 6.5637682784520328472131852577044702148066953193247361, 1e-12);
  EXPECT_NEAR(f[1].d(0), 0.0152504761193336512405911699894491769284000629403064, 1e-12);
}

TEST(doubleADTest, testSingleDimLog1P)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 4.9;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::log1p(f[0]);
  EXPECT_NEAR(f[1].v(), 1.7749523509116737665980589190353334164622453007016815, 1e-12);
  EXPECT_NEAR(f[1].d(0), 0.169491525423728813559322033898305084745762711864406779661, 1e-12);
}

TEST(doubleADTest, testSingleDimSqrt)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 25;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::sqrt(f[0]);
  EXPECT_EQ(f[1].v(), 5);
  EXPECT_EQ(f[1].d(0), 0.1);
}

TEST(doubleADTest, testSingleDimSinus)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::sin(f[0]);
  EXPECT_NEAR(f[1].v(), -0.87969575997167, 1e-12);
  EXPECT_NEAR(f[1].d(0), -0.475536927995993, 1e-12);
}

TEST(doubleADTest, testSingleDimCosinus)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::cos(f[0]);
  EXPECT_NEAR(f[1].v(), -0.475536927995993, 1e-12);
  EXPECT_NEAR(f[1].d(0), 0.87969575997167, 1e-12);
}

TEST(doubleADTest, testSingleDimTangens)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::tan(f[0]);
  EXPECT_NEAR(f[1].v(), 1.84989999342193, 1e-12);
  EXPECT_NEAR(f[1].d(0), 4.42212998566244, 1e-12);
}

TEST(doubleADTest, testSingleDimArcsine)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 0.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::asin(f[0]);
  EXPECT_NEAR(f[1].v(), 0.523598775598299, 1e-12);
  EXPECT_NEAR(f[1].d(0), 1.15470053837925, 1e-12);
}

TEST(doubleADTest, testSingleDimArccosine)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 0.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::acos(f[0]);
  EXPECT_NEAR(f[1].v(), 1.0471975511966, 1e-12);
  EXPECT_NEAR(f[1].d(0), -1.15470053837925, 1e-12);
}

TEST(doubleADTest, testSingleDimArctangent)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 0.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::atan(f[0]);
  EXPECT_NEAR(f[1].v(), 0.463647609000806, 1e-12);
  EXPECT_NEAR(f[1].d(0), 0.8, 1e-12);
}

TEST(doubleADTest, testSingleDimArctangent2)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 0.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::atan2(f[0], 1.0);
  EXPECT_NEAR(f[1].v(), 0.463647609000806, 1e-12);
  EXPECT_NEAR(f[1].d(0), 0.8, 1e-12);

  f[1] = util::atan2(1.0, f[0]);
  EXPECT_NEAR(f[1].v(), 1.107148717794, 1e-12);
  EXPECT_NEAR(f[1].d(0), -0.8, 1e-12);

  f[1] = util::atan2(f[0], f[0]);
  EXPECT_NEAR(f[1].v(), 0.785398163397, 1e-12);
  EXPECT_NEAR(f[1].d(0), 0, 1e-12);
}

TEST(doubleADTest, testSingleDimHyperbolicsine)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 1.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::sinh(f[0]);
  EXPECT_NEAR(f[1].v(), 2.1292794550948,1e-6);
  EXPECT_NEAR(f[1].d(0), 2.3524096152432,1e-6);
}

TEST(doubleADTest, testSingleDimHyperboliccosine)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 1.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::cosh(f[0]);
  EXPECT_NEAR(f[1].v(), 2.3524096152432,1e-6);
  EXPECT_NEAR(f[1].d(0), 2.1292794550948,1e-6);
}

TEST(doubleADTest, testSingleDimHyperbolictangent)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 1.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::tanh(f[0]);
  EXPECT_NEAR(f[1].v(), 0.90514825364487,1e-6);
  EXPECT_NEAR(f[1].d(0), 0.18070663892,1e-6);
}

TEST(doubleADTest, testSingleDimArcyperbolicsine)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 1.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::asinh(f[0]);
  EXPECT_NEAR(f[1].v(), 1.194763217,1e-6);
  EXPECT_NEAR(f[1].d(0), 0.5547001962,1e-6);
}

TEST(doubleADTest, testSingleDimArcyperboliccosine)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 1.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::acosh(f[0]);
  EXPECT_NEAR(f[1].v(), 0.9624236501,1e-6);
  EXPECT_NEAR(f[1].d(0), 0.894427191,1e-6);
}

TEST(doubleADTest, testSingleDimArcyperbolictangent)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 0.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::atanh(f[0]);
  EXPECT_NEAR(f[1].v(), 0.5493061443,1e-6);
  EXPECT_NEAR(f[1].d(0), 1.33333333333,1e-6);
}

TEST(doubleADTest, testMultiDimFmod)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  S y = 3.0;
  T f[3];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = y;
  f[1].setDiffVariable(1);
  f[2] = util::fmod(f[0], f[1]);
  EXPECT_NEAR(f[2].v(), 1.5, 1e-12);
  EXPECT_NEAR(f[2].d(0), 1.0, 1e-12);
  EXPECT_NEAR(f[2].d(1), -3.0, 1e-12);
}

TEST(doubleADTest, testMultiDimMultiplication)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = f[0]*f[0];
  EXPECT_EQ(f[1].v(), 110.25);
  EXPECT_EQ(f[1].d(0), 21.0);
  EXPECT_EQ(f[1].d(1), 0.0);
}

// AD Tests utilizing type double with multiple variables

TEST(doubleADTest, testSingleDimAdditionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 10.5;
  S x2 = 11.25;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = f[0]+f[1];
  EXPECT_EQ(f[2].v(), 21.75);
  EXPECT_EQ(f[2].d(0), 1.0);
  EXPECT_EQ(f[2].d(1), 1.0);
}

TEST(doubleADTest, testSingleDimSubtractionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 10.5;
  S x2 = 11.25;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = f[0]-f[1];
  EXPECT_EQ(f[2].v(), -0.75);
  EXPECT_EQ(f[2].d(0), 1.0);
  EXPECT_EQ(f[2].d(1), -1.0);
}

TEST(doubleADTest, testSingleDimMultiplicationMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 10.5;
  S x2 = 11.25;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = f[0]*f[1];
  EXPECT_EQ(f[2].v(), 118.125);
  EXPECT_EQ(f[2].d(0), 11.25);
  EXPECT_EQ(f[2].d(1), 10.5);
}

TEST(doubleADTest, testSingleDimDivisionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 2.5;
  S x2 = 2;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = f[1]/f[0];
  EXPECT_NEAR(f[2].v(), 0.8, 1e-12);
  EXPECT_NEAR(f[2].d(0), -0.32, 1e-12);
  EXPECT_NEAR(f[2].d(1), 0.4, 1e-12);
}

TEST(doubleADTest, testSingleDimPowMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 2.5;
  S x2 = 2.0;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::pow(f[0],f[1]);
  EXPECT_NEAR(f[2].v(), 6.25, 1e-15);
  EXPECT_NEAR(f[2].d(0), 5, 1e-15);
}

TEST(doubleADTest, testSingleDimExpAdditionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::exp(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), 6.685894442, 1e-9);
  EXPECT_NEAR(f[2].d(0), 6.685894442, 1e-9);
}

TEST(doubleADTest, testSingleDimExpMultiplicationMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::exp(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), 1.8221188, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.7288475202, 1e-9);
}

TEST(doubleADTest, testSingleDimExpDivisionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::exp(f[0] / f[1]);
  EXPECT_NEAR(f[2].v(), 42.521082, 1e-9);
  EXPECT_NEAR(f[2].d(0), 106.302705, 1e-9);
}

TEST(doubleADTest, testSingleDimExpSubtractionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::exp(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 3.004166024, 1e-9);
  EXPECT_NEAR(f[2].d(0), 3.004166024, 1e-9);
}

TEST(doubleADTest, testSingleDimSqrAdditionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sqr(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), 3.61, 1e-10);
  EXPECT_NEAR(f[2].d(0), 3.8, 1e-10);
}

TEST(doubleADTest, testSingleDimSqrMultiplicationMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sqr(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), 0.36, 1e-10);
  EXPECT_NEAR(f[2].d(0), 0.48, 1e-10);
}

TEST(doubleADTest, testSingleDimSqrDivisionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sqr(f[0] / f[1]);
  EXPECT_NEAR(f[2].v(), 14.0625, 1e-10);
  EXPECT_NEAR(f[2].d(0), 18.75, 1e-10);
}

TEST(doubleADTest, testSingleDimSqrSubtractionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sqr(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 1.21, 1e-10);
  EXPECT_NEAR(f[2].d(0), 2.2, 1e-10);
}

TEST(doubleADTest, testSingleDimLogAdditionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::log(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), 0.6418538862, 1e-10);
  EXPECT_NEAR(f[2].d(0),0.5263157895, 1e-10);
}

TEST(doubleADTest, testSingleDimLogMultiplicationMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::log(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), -0.5108256238, 1e-10);
  EXPECT_NEAR(f[2].d(0), 0.66666666666, 1e-10);
}

TEST(doubleADTest, testSingleDimLogDivisionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::log(f[0] / f[1]);
  EXPECT_NEAR(f[2].v(), 1.32175584, 1e-10);
  EXPECT_NEAR(f[2].d(0), 0.66666666666, 1e-10);
}

TEST(doubleADTest, testSingleDimLogSubtractionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::log(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 0.0953101798, 1e-10);
  EXPECT_NEAR(f[2].d(0), 0.90909090909090, 1e-10);
}

TEST(doubleADTest, testSingleDimSqrtAdditionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sqrt(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), 1.378404875, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.3627381251, 1e-9);
}

TEST(doubleADTest, testSingleDimSqrtMultiplicationMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sqrt(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), 0.7745966692, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.2581988897, 1e-9);
}

TEST(doubleADTest, testSingleDimSqrtDivisionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sqrt(f[0] / f[1]);
  EXPECT_NEAR(f[2].v(), 1.936491673, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.6454972244, 1e-9);
}

TEST(doubleADTest, testSingleDimSqrtSubtractionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sqrt(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 1.048808848, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.4767312946, 1e-9);
}

TEST(doubleADTest, testSingleDimSinusAdditionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sin(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), 0.9463000877, 1e-9);
  EXPECT_NEAR(f[2].d(0), -0.3232895669, 1e-9);
}

TEST(doubleADTest, testSingleDimSinusMultiplicationMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sin(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), 0.5646424734, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.330134246, 1e-9);
}

TEST(doubleADTest, testSingleDimSinusDivisionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sin(f[0] / f[1]);
  EXPECT_NEAR(f[2].v(), -0.5715613187, 1e-9);
  EXPECT_NEAR(f[2].d(0), -2.051398393, 1e-9);
}

TEST(doubleADTest, testSingleDimSinusSubtractionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sin(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 0.8912073601, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.4535961214, 1e-9);
}

TEST(doubleADTest, testSingleDimCosinusAdditionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::cos(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), -0.3232895669, 1e-9);
  EXPECT_NEAR(f[2].d(0), -0.9463000877, 1e-9);
}

TEST(doubleADTest, testSingleDimCosinusMultiplicationMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::cos(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), 0.8253356149, 1e-9);
  EXPECT_NEAR(f[2].d(0), -0.2258569894, 1e-9);
}

TEST(doubleADTest, testSingleDimCosinusDivisionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::cos(f[0] / f[1]);
  EXPECT_NEAR(f[2].v(), -0.8205593573, 1e-9);
  EXPECT_NEAR(f[2].d(0), 1.428903297, 1e-9);
}

TEST(doubleADTest, testSingleDimCosinusSubtractionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::cos(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 0.4535961214, 1e-9);
  EXPECT_NEAR(f[2].d(0), -0.8912073601, 1e-9);
}

TEST(doubleADTest, testSingleDimTangensAdditionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::tan(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), -2.927097515, 1e-9);
  EXPECT_NEAR(f[2].d(0), 9.56789986, 1e-9);
}

TEST(doubleADTest, testSingleDimTangensMultiplicationMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::tan(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), 0.6841368083, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.587217269, 1e-9);
}

TEST(doubleADTest, testSingleDimTangensDivisionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::tan(f[0] / f[1]);
  EXPECT_NEAR(f[2].v(), 0.6965508511, 1e-9);
  EXPECT_NEAR(f[2].d(0), 3.71295772, 1e-9);
}

TEST(doubleADTest, testSingleDimTangensSubtractionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::tan(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 1.964759657, 1e-9);
  EXPECT_NEAR(f[2].d(0), 4.860280511, 1e-9);
}

TEST(doubleADTest, testSingleDimArcsineAdditionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::asin(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), 1.119769515, 1e-9);
  EXPECT_NEAR(f[2].d(0), 2.294157339, 1e-9);
}

TEST(doubleADTest, testSingleDimArcsineMultiplicationMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::asin(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), 0.2013579208, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.4082482905, 1e-9);
}

TEST(doubleADTest, testSingleDimArcsineDivisionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::asin(f[1] / f[0]);
  EXPECT_NEAR(f[2].v(), 0.927295218, 1e-9);
  EXPECT_NEAR(f[2].d(0), -2.66666666666, 1e-9);
}

TEST(doubleADTest, testSingleDimArcsineSubtractionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::asin(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 0.1001674212, 1e-9);
  EXPECT_NEAR(f[2].d(0), 1.005037815, 1e-9);
}

TEST(doubleADTest, testSingleDimArccosineAdditionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::acos(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), 0.4510268118, 1e-9);
  EXPECT_NEAR(f[2].d(0), -2.294157339, 1e-9);
}

TEST(doubleADTest, testSingleDimArccosineMultplicationMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::acos(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), 1.369438406, 1e-9);
  EXPECT_NEAR(f[2].d(0), -0.4082482905, 1e-9);
}

TEST(doubleADTest, testSingleDimArccosineDivisionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::acos(f[1] / f[0]);
  EXPECT_NEAR(f[2].v(), 0.6435011088, 1e-9);
  EXPECT_NEAR(f[2].d(0), 2.66666666666, 1e-9);
}

TEST(doubleADTest, testSingleDimArccosineSubtractionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::acos(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 1.470628906, 1e-9);
  EXPECT_NEAR(f[2].d(0), -1.005037815, 1e-9);
}

TEST(doubleADTest, testSingleDimArctangentAdditionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::atan(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), 0.7328151018, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.5524861878, 1e-9);
}

TEST(doubleADTest, testSingleDimArctangentMultiplicationMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::atan(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), 0.1973955598, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.3846153846, 1e-9);
}

TEST(doubleADTest, testSingleDimArctangentDivisionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::atan(f[1] / f[0]);
  EXPECT_NEAR(f[2].v(), 0.6747409422, 1e-9);
  EXPECT_NEAR(f[2].d(0), -0.9756097561, 1e-9);
}

TEST(doubleADTest, testSingleDimArctangentSubtractionMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::atan(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 0.09966865249, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.99009900990099, 1e-9);
}


TEST(doubleADTest, testMultiDimMultiplicationMultipleVariables)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 10.5;
  S x2 = 11.25;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = f[0]*f[1];
  EXPECT_EQ(f[2].v(), 118.125);
  EXPECT_EQ(f[2].d(0), 11.25);
  EXPECT_EQ(f[2].d(1), 10.5);
}

TEST(doubleADTest, testhasZeroDerivative)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = 5;
  f[1] = x;
  f[1].setDiffVariable(0);
  EXPECT_TRUE(f[0].hasZeroDerivative());
  EXPECT_FALSE(f[1].hasZeroDerivative());
}

TEST(doubleADTest, testEQADfADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  T y = 2.5;
  EXPECT_TRUE(x == y);
}

TEST(doubleADTest, testEQdoubleADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(2.5 == x);
}

TEST(doubleADTest, testEQADfdouble)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(x == 2.5);
}

TEST(doubleADTest, testNotEQADfADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x1 = 2.5;
  T x2 = 2.4;
  EXPECT_TRUE(x1 != x2);
}

TEST(doubleADTest, testNotEQdoubleADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(2.6 != x);
}

TEST(doubleADTest, testNotEQADfdouble)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(x != 2.6);
}

TEST(doubleADTest, testSmallerADfADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x1 = 2.5;
  T x2 = 2.4;
  EXPECT_TRUE(x1 > x2);
}

TEST(doubleADTest, testSmallerdoubleADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(2.6 > x);
}

TEST(doubleADTest, testSmallerADfdouble)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(x > 2.4);
}

TEST(doubleADTest, testSmallerEQADfADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x1 = 2.5;
  T x2 = 2.4;
  EXPECT_TRUE(x1 >= x2);
}

TEST(doubleADTest, testSmallerEQdoubleADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(2.5 >= x);
}

TEST(doubleADTest, testSmallerEQADfdouble)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(x >= 2.4);
}

TEST(doubleADTest, testBiggerADfADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x1 = 2.5;
  T x2 = 2.6;
  EXPECT_TRUE(x1 < x2);
}

TEST(doubleADTest, testBiggerdoubleADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(2.4 < x);
}

TEST(doubleADTest, testBiggerADfdouble)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(x < 2.6);
}

TEST(doubleADTest, testBiggerEQADfADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x1 = 2.5;
  T x2 = 2.6;
  EXPECT_TRUE(x1 <= x2);
}

TEST(doubleADTest, testBiggerEQdoubleADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(2.5 <= x);
}

TEST(doubleADTest, testBiggerEQADfdouble)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(x <= 2.6);
}

TEST(doubleADTest, testFloor)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  x.setDiffVariable(0);
  EXPECT_TRUE(util::floor(x) == 2);
  EXPECT_EQ(util::floor(x).d(0), 0.0);
}

TEST(doubleADTest, testFabs1)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = -2.5;
  x.setDiffVariable(0);
  EXPECT_TRUE(util::fabs(x) == 2.5);
  EXPECT_EQ(x.d(0), 1.0);
}

TEST(doubleADTest, testFabs2)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 5.25;
  x.setDiffVariable(0);
  EXPECT_TRUE(util::fabs(x) == 5.25);
  EXPECT_EQ(x.d(0), 1.0);
}

TEST(doubleADTest, testFabs3)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 0;
  x.setDiffVariable(0);
  EXPECT_TRUE(util::fabs(x) == 0);
}

TEST(doubleADTest, testMaxADfdouble)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;

  EXPECT_TRUE(util::max(x, 5.25) == 5.25);
}

TEST(doubleADTest, testMaxADfADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x1 = 2.5;
  T x2 = 2.0;
  EXPECT_TRUE(util::max(x1, x2) == 2.5);
}

TEST(doubleADTest, testMaxdoubleADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(util::max(7.25, x) == 7.25);
}

TEST(doubleADTest, testMinADfdouble)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(util::min(x, 2.125) == 2.125);
}

TEST(doubleADTest, testMinADfADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x1 = 2.5;
  T x2 = 2.0;
  EXPECT_TRUE(util::min(x1, x2) == 2.0);
}

TEST(doubleADTest, testMindoubleADf)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(util::min(2.125, x) == 2.125);
}

TEST(doubleADTest, testNearZero)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  T x = 1e-18;
  EXPECT_TRUE(nearZero(x));
}

// AD Tests utilizing type long double

TEST(longdoubleADTest, testConstructFromNoArgument)
{
  typedef double S;
  S x = 10.5;
  util::ADf<S,DIM2> f;
  f = x;
  f.setDiffVariable(0);
  EXPECT_EQ(f.v(), 10.5);
  EXPECT_EQ(f.d(0), 1.0);
}

TEST(longdoubleADTest, testConstructFromADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f;
  f = x;
  f.setDiffVariable(0);
  T g(f);
  EXPECT_EQ(g.v(), 10.5);
  EXPECT_EQ(g.d(0), 1.0);
}

TEST(longdoubleADTest, testConstructFromLongDouble)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f(x);
  f.setDiffVariable(0);
  EXPECT_EQ(f.v(), 10.5);
  EXPECT_EQ(f.d(0), 1.0);
}

TEST(longdoubleADTest, testConstructFromLongDoubleWithDerivative)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  Vector<S,2> d = {5,0};
  T f(x,d);
  EXPECT_EQ(f.v(), 10.5);
  EXPECT_EQ(f.d(0), 5);
}

TEST(longdoubleADTest, testConstructFromLongDoubleWithMultiDimDerivative)
{
  typedef long double S;
  typedef util::ADf<S,5> T;
  S x = 10.5;
  Vector<S,5> d = {5,2,8,1.1,-100};
  T f(x,d);
  EXPECT_EQ(f.v(), 10.5);
  EXPECT_EQ(f.d(), d);
}

TEST(longdoubleADTest, testSingleDimConst)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[1];
  f[0] = x;
  f[0].setDiffVariable(0);
  EXPECT_EQ(f[0].v(), 10.5);
  EXPECT_EQ(f[0].d(0), 1.0);
}

TEST(longdoubleADTest, testSingleDimAdditionADfdouble)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = -f[0] + 9.25;
  EXPECT_EQ(f[1].v(), -1.25);
  EXPECT_EQ(f[1].d(0), -1.0);
}

TEST(longdoubleADTest, testSingleDimAdditionADfADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = f[0]+f[0];
  EXPECT_EQ(f[1].v(), 21.0);
  EXPECT_EQ(f[1].d(0), 2.0);
}

TEST(longdoubleADTest, testSingleDimAdditionDoubleADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 0.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] =  9.25 + -f[0];
  EXPECT_EQ(f[1].v(), 8.75);
  EXPECT_EQ(f[1].d(0), -1.0);
}

TEST(longdoubleADTest, testSingleDimAdditionDoubleDouble)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T f = 2.5 + -1.2;
  EXPECT_EQ(f.v(), 1.3);
  EXPECT_EQ(f.d(0), 0.0);
}

TEST(longdoubleADTest, testSingleDimSubtraktionADfdouble)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = f[0] - 19.25;
  EXPECT_EQ(f[1].v(), -8.75);
  EXPECT_EQ(f[1].d(0), 1.0);
}

TEST(longdoubleADTest, testSingleDimSubtraktionADfADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = f[0] - f[0];
  EXPECT_EQ(f[1].v(), 0.0);
  EXPECT_EQ(f[1].d(0), 0.0);
}

TEST(longdoubleADTest, testSingleDimSubtraktionDoubleADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = -10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] =  9.25 - f[0];
  EXPECT_EQ(f[1].v(), 19.75);
  EXPECT_EQ(f[1].d(0), -1.0);
}

TEST(longdoubleADTest, testSingleDimSubtraktionDoubleDouble)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T f = 2.5 - 1.2;
  EXPECT_EQ(f.v(), 1.3);
  EXPECT_EQ(f.d(0), 0.0);
}

TEST(longdoubleADTest, testSingleDimMultiplicationADfdouble)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = f[0] * 9.25;
  EXPECT_EQ(f[1].v(), 97.125);
  EXPECT_EQ(f[1].d(0), 9.25);
}

TEST(longdoubleADTest, testSingleDimMultiplicationADfADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = f[0]*f[0];
  EXPECT_EQ(f[1].v(), 110.25);
  EXPECT_EQ(f[1].d(0), 21.0);
}

TEST(longdoubleADTest, testSingleDimMultiplicationDoubleADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = -10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] =  9.25 * -f[0];
  EXPECT_EQ(f[1].v(), 97.125);
  EXPECT_EQ(f[1].d(0), -9.25);
}

TEST(longdoubleADTest, testSingleDimMultiplicationDoubleDouble)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T f = 2.5 * 1.2;
  EXPECT_EQ(f.v(), 3.0);
  EXPECT_EQ(f.d(0), 0.0);
}

TEST(longdoubleADTest, testSingleDimDivisionDoubleADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 2.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = 1.0/f[0];
  EXPECT_NEAR(f[1].v(), 0.4, 1e-12);
  EXPECT_NEAR(f[1].d(0), -0.16, 1e-12);
}

TEST(longdoubleADTest, testSingleDimDivisionADfADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = f[0]/f[0];
  EXPECT_EQ(f[1].v(), 1.0);
  EXPECT_EQ(f[1].d(0), 0.0);
}

TEST(longdoubleADTest, testSingleDimDivisionADfDouble)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 2.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = f[0]/-5.0;
  EXPECT_NEAR(f[1].v(), -0.5, 1e-12);
  EXPECT_NEAR(f[1].d(0), -0.2, 1e-12);
}

TEST(longdoubleADTest, testSingleDimDivisionDoubleDouble)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T f = 3.4 / -1.6;
  EXPECT_EQ(f.v(), -2.125);
  EXPECT_EQ(f.d(0), 0.0);
}

TEST(longdoubleADTest, testSingleDimPow)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 2.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::pow(f[0],-1);
  EXPECT_NEAR(f[1].v(), 0.4, 1e-12);
  EXPECT_NEAR(f[1].d(0), -0.16, 1e-12);
}

TEST(longdoubleADTest, testSingleDimPowZeroBaseADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 0;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::pow(f[0], 5.5);
  EXPECT_NEAR(f[1].v(), 0, 1e-15);
  EXPECT_NEAR(f[1].d(0), 0, 1e-15);
}

TEST(longdoubleADTest, testSingleDimPowNegativeBaseADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = -2.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::pow(f[0], 3.0);
  EXPECT_NEAR(f[1].v(), -15.625, 1e-15);
  EXPECT_NEAR(f[1].d(0), 18.75, 1e-15);
}


TEST(longdoubleADTest, testSingleDimExpl)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = expl(f[0]);
  EXPECT_NEAR(f[1].v(), 36315.502674246636, 1e-12);
  EXPECT_NEAR(f[1].d(0), 36315.502674246636, 1e-12);
}

TEST(longdoubleADTest, testSingleDimExp)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::exp(f[0]);
  EXPECT_NEAR(f[1].v(), 36315.502674246636, 1e-12);
  EXPECT_NEAR(f[1].d(0), 36315.502674246636, 1e-12);
}

TEST(longdoubleADTest, testSingleDimSqr)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::sqr(f[0]);
  EXPECT_EQ(f[1].v(), 110.25);
  EXPECT_EQ(f[1].d(0), 21.0);
}

TEST(longdoubleADTest, testSingleDimLogl)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 4.7;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = logl(f[0]);
  EXPECT_NEAR(f[1].v(), 1.54756250871601, 1e-12);
  EXPECT_NEAR(f[1].d(0), 0.212765957446809, 1e-12);
}

TEST(longdoubleADTest, testSingleDimLog)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 4.7;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::log(f[0]);
  EXPECT_NEAR(f[1].v(), 1.54756250871601, 1e-12);
  EXPECT_NEAR(f[1].d(0), 0.212765957446809, 1e-12);
}

TEST(longdoubleADTest, testSingleDimLog10l)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 53.7;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = log10l(f[0]);
  EXPECT_NEAR(f[1].v(), 1.7299742856995556068654105950281739286312009548692244, 1e-12);
  EXPECT_NEAR(f[1].d(0), 0.0080874205196136280754400171120410629849980820447610, 1e-12);
}

TEST(longdoubleADTest, testSingleDimLog10)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 53.7;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::log10(f[0]);
  EXPECT_NEAR(f[1].v(), 1.7299742856995556068654105950281739286312009548692244, 1e-12);
  EXPECT_NEAR(f[1].d(0), 0.0080874205196136280754400171120410629849980820447610, 1e-12);
}

TEST(longdoubleADTest, testSingleDimLog2l)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 94.6;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = log2l(f[0]);
  EXPECT_NEAR(f[1].v(), 6.5637682784520328472131852577044702148066953193247361, 1e-12);
  EXPECT_NEAR(f[1].d(0), 0.0152504761193336512405911699894491769284000629403064, 1e-12);
}

TEST(longdoubleADTest, testSingleDimLog2)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 94.6;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::log2(f[0]);
  EXPECT_NEAR(f[1].v(), 6.5637682784520328472131852577044702148066953193247361, 1e-12);
  EXPECT_NEAR(f[1].d(0), 0.0152504761193336512405911699894491769284000629403064, 1e-12);
}

TEST(longdoubleADTest, testSingleDimLog1Pl)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 4.9;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = log1pl(f[0]);
  EXPECT_NEAR(f[1].v(), 1.7749523509116737665980589190353334164622453007016815, 1e-12);
  EXPECT_NEAR(f[1].d(0), 0.169491525423728813559322033898305084745762711864406779661, 1e-12);
}

TEST(longdoubleADTest, testSingleDimLog1P)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 4.9;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::log1p(f[0]);
  EXPECT_NEAR(f[1].v(), 1.7749523509116737665980589190353334164622453007016815, 1e-12);
  EXPECT_NEAR(f[1].d(0), 0.169491525423728813559322033898305084745762711864406779661, 1e-12);
}

TEST(longdoubleADTest, testSingleDimSqrt)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 25.25;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::sqrt(f[0]);
  EXPECT_NEAR(f[1].v(), 5.02493781056044, 1e-12);
  EXPECT_NEAR(f[1].d(0), 0.09950371902099892, 1e-12);
}

TEST(longdoubleADTest, testSingleDimSinus)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::sin(f[0]);
  EXPECT_NEAR(f[1].v(), -0.87969575997167, 1e-12);
  EXPECT_NEAR(f[1].d(0), -0.475536927995993, 1e-12);
}

TEST(longdoubleADTest, testSingleDimCosinus)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::cos(f[0]);
  EXPECT_NEAR(f[1].v(), -0.475536927995993, 1e-12);
  EXPECT_NEAR(f[1].d(0), 0.87969575997167, 1e-12);
}

TEST(longdoubleADTest, testSingleDimTangens)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::tan(f[0]);
  EXPECT_NEAR(f[1].v(), 1.84989999342193, 1e-12);
  EXPECT_NEAR(f[1].d(0), 4.42212998566244, 1e-12);
}

TEST(longdoubleADTest, testSingleDimArcsine)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 0.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::asin(f[0]);
  EXPECT_NEAR(f[1].v(), 0.523598775598299, 1e-12);
  EXPECT_NEAR(f[1].d(0), 1.15470053837925, 1e-12);
}

TEST(longdoubleADTest, testSingleDimArccosine)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 0.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::acos(f[0]);
  EXPECT_NEAR(f[1].v(), 1.0471975511966, 1e-12);
  EXPECT_NEAR(f[1].d(0), -1.15470053837925, 1e-12);
}

TEST(longdoubleADTest, testSingleDimArctangent)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 0.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::atan(f[0]);
  EXPECT_NEAR(f[1].v(), 0.463647609000806, 1e-12);
  EXPECT_NEAR(f[1].d(0), 0.8, 1e-12);
}

TEST(longdoubleADTest, testSingleDimHyperbolicsine)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 1.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::sinh(f[0]);
  EXPECT_NEAR(f[1].v(), 2.1292794550948,1e-6);
  EXPECT_NEAR(f[1].d(0), 2.3524096152432,1e-6);
}

TEST(longdoubleADTest, testSingleDimHyperboliccosine)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 1.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::cosh(f[0]);
  EXPECT_NEAR(f[1].v(), 2.3524096152432,1e-6);
  EXPECT_NEAR(f[1].d(0), 2.1292794550948,1e-6);
}

TEST(longdoubleADTest, testSingleDimHyperbolictangent)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 1.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::tanh(f[0]);
  EXPECT_NEAR(f[1].v(), 0.90514825364487,1e-6);
  EXPECT_NEAR(f[1].d(0), 0.18070663892,1e-6);
}

TEST(longdoubleADTest, testSingleDimArcyperbolicsine)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 1.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::asinh(f[0]);
  EXPECT_NEAR(f[1].v(), 1.194763217,1e-6);
  EXPECT_NEAR(f[1].d(0), 0.5547001962,1e-6);
}

TEST(longdoubleADTest, testSingleDimArcyperboliccosine)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 1.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::acosh(f[0]);
  EXPECT_NEAR(f[1].v(), 0.9624236501,1e-6);
  EXPECT_NEAR(f[1].d(0), 0.894427191,1e-6);
}

TEST(longdoubleADTest, testSingleDimArcyperbolictangent)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 0.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = util::atanh(f[0]);
  EXPECT_NEAR(f[1].v(), 0.5493061443,1e-6);
  EXPECT_NEAR(f[1].d(0), 1.33333333333,1e-6);
}

TEST(longdoubleADTest, testMultiDimFmod)
{
  typedef double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  S y = 3.0;
  T f[3];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = y;
  f[1].setDiffVariable(1);
  f[2] = util::fmod(f[0], f[1]);
  EXPECT_NEAR(f[2].v(), 1.5, 1e-12);
  EXPECT_NEAR(f[2].d(0), 1.0, 1e-12);
  EXPECT_NEAR(f[2].d(1), -3.0, 1e-12);
}

TEST(longdoubleADTest, testMultiDimMultiplication)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = f[0]*f[0];
  EXPECT_EQ(f[1].v(), 110.25);
  EXPECT_EQ(f[1].d(0), 21.0);
  EXPECT_EQ(f[1].d(1), 0.0);
}

// AD Tests utilizing type long double with multiple variables

TEST(longdoubleADTest, testSingleDimAdditionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 10.5;
  S x2 = 11.25;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = f[0]+f[1];
  EXPECT_EQ(f[2].v(), 21.75);
  EXPECT_EQ(f[2].d(0), 1.0);
  EXPECT_EQ(f[2].d(1), 1.0);
}

TEST(longdoubleADTest, testSingleDimSubtractionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 10.5;
  S x2 = 11.25;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = f[0]-f[1];
  EXPECT_EQ(f[2].v(), -0.75);
  EXPECT_EQ(f[2].d(0), 1.0);
  EXPECT_EQ(f[2].d(1), -1.0);
}

TEST(longdoubleADTest, testSingleDimMultiplicationMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 10.5;
  S x2 = 11.25;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = f[0]*f[1];
  EXPECT_EQ(f[2].v(), 118.125);
  EXPECT_EQ(f[2].d(0), 11.25);
  EXPECT_EQ(f[2].d(1), 10.5);
}

TEST(longdoubleADTest, testSingleDimDivisionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 2.5;
  S x2 = 2;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = f[1]/f[0];
  EXPECT_NEAR(f[2].v(), 0.8, 1e-12);
  EXPECT_NEAR(f[2].d(0), -0.32, 1e-12);
  EXPECT_NEAR(f[2].d(1), 0.4, 1e-12);
}

TEST(longdoubleADTest, testSingleDimPowMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 2.5;
  S x2 = 2.0;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::pow(f[0],f[1]);
  EXPECT_EQ(f[2].v(), 6.25);
  EXPECT_EQ(f[2].d(0), 5);
}

TEST(longdoubleADTest, testSingleDimExpAdditionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::exp(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), 6.685894442, 1e-9);
  EXPECT_NEAR(f[2].d(0), 6.685894442, 1e-9);
}

TEST(longdoubleADTest, testSingleDimExpMultiplicationMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::exp(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), 1.8221188, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.7288475202, 1e-9);
}

TEST(longdoubleADTest, testSingleDimExpDivisionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::exp(f[0] / f[1]);
  EXPECT_NEAR(f[2].v(), 42.521082, 1e-9);
  EXPECT_NEAR(f[2].d(0), 106.302705, 1e-9);
}

TEST(longdoubleADTest, testSingleDimExpSubtractionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::exp(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 3.004166024, 1e-9);
  EXPECT_NEAR(f[2].d(0), 3.004166024, 1e-9);
}

TEST(longdoubleADTest, testSingleDimSqrAdditionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sqr(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), 3.61, 1e-10);
  EXPECT_NEAR(f[2].d(0), 3.8, 1e-10);
}

TEST(longdoubleADTest, testSingleDimSqrMultiplicationMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sqr(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), 0.36, 1e-10);
  EXPECT_NEAR(f[2].d(0), 0.48, 1e-10);
}

TEST(longdoubleADTest, testSingleDimSqrDivisionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sqr(f[0] / f[1]);
  EXPECT_NEAR(f[2].v(), 14.0625, 1e-10);
  EXPECT_NEAR(f[2].d(0), 18.75, 1e-10);
}

TEST(longdoubleADTest, testSingleDimSqrSubtractionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sqr(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 1.21, 1e-10);
  EXPECT_NEAR(f[2].d(0), 2.2, 1e-10);
}

TEST(longdoubleADTest, testSingleDimLogAdditionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::log(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), 0.6418538862, 1e-10);
  EXPECT_NEAR(f[2].d(0),0.5263157895, 1e-10);
}

TEST(longdoubleADTest, testSingleDimLogMultiplicationMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::log(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), -0.5108256238, 1e-10);
  EXPECT_NEAR(f[2].d(0), 0.66666666666, 1e-10);
}

TEST(longdoubleADTest, testSingleDimLogDivisionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::log(f[0] / f[1]);
  EXPECT_NEAR(f[2].v(), 1.32175584, 1e-10);
  EXPECT_NEAR(f[2].d(0), 0.66666666666, 1e-10);
}

TEST(longdoubleADTest, testSingleDimLogSubtractionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::log(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 0.0953101798, 1e-10);
  EXPECT_NEAR(f[2].d(0), 0.90909090909090, 1e-10);
}

TEST(longdoubleADTest, testSingleDimSqrtAdditionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sqrt(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), 1.378404875, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.3627381251, 1e-9);
}

TEST(longdoubleADTest, testSingleDimSqrtMultiplicationMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sqrt(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), 0.7745966692, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.2581988897, 1e-9);
}

TEST(longdoubleADTest, testSingleDimSqrtDivisionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sqrt(f[0] / f[1]);
  EXPECT_NEAR(f[2].v(), 1.936491673, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.6454972244, 1e-9);
}

TEST(longdoubleADTest, testSingleDimSqrtSubtractionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sqrt(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 1.048808848, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.4767312946, 1e-9);
}

TEST(longdoubleADTest, testSingleDimSinusAdditionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sin(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), 0.9463000877, 1e-9);
  EXPECT_NEAR(f[2].d(0), -0.3232895669, 1e-9);
}

TEST(longdoubleADTest, testSingleDimSinusMultiplicationMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sin(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), 0.5646424734, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.330134246, 1e-9);
}

TEST(longdoubleADTest, testSingleDimSinusDivisionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sin(f[0] / f[1]);
  EXPECT_NEAR(f[2].v(), -0.5715613187, 1e-9);
  EXPECT_NEAR(f[2].d(0), -2.051398393, 1e-9);
}

TEST(longdoubleADTest, testSingleDimSinusSubtractionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::sin(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 0.8912073601, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.4535961214, 1e-9);
}

TEST(longdoubleADTest, testSingleDimCosinusAdditionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::cos(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), -0.3232895669, 1e-9);
  EXPECT_NEAR(f[2].d(0), -0.9463000877, 1e-9);
}

TEST(longdoubleADTest, testSingleDimCosinusMultiplicationMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::cos(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), 0.8253356149, 1e-9);
  EXPECT_NEAR(f[2].d(0), -0.2258569894, 1e-9);
}

TEST(longdoubleADTest, testSingleDimCosinusDivisionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::cos(f[0] / f[1]);
  EXPECT_NEAR(f[2].v(), -0.8205593573, 1e-9);
  EXPECT_NEAR(f[2].d(0), 1.428903297, 1e-9);
}

TEST(longdoubleADTest, testSingleDimCosinusSubtractionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::cos(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 0.4535961214, 1e-9);
  EXPECT_NEAR(f[2].d(0), -0.8912073601, 1e-9);
}

TEST(longdoubleADTest, testSingleDimTangensAdditionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::tan(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), -2.927097515, 1e-9);
  EXPECT_NEAR(f[2].d(0), 9.56789986, 1e-9);
}

TEST(longdoubleADTest, testSingleDimTangensMultiplicationMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::tan(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), 0.6841368083, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.587217269, 1e-9);
}

TEST(longdoubleADTest, testSingleDimTangensDivisionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::tan(f[0] / f[1]);
  EXPECT_NEAR(f[2].v(), 0.6965508511, 1e-9);
  EXPECT_NEAR(f[2].d(0), 3.71295772, 1e-9);
}

TEST(longdoubleADTest, testSingleDimTangensSubtractionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 1.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::tan(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 1.964759657, 1e-9);
  EXPECT_NEAR(f[2].d(0), 4.860280511, 1e-9);
}

TEST(longdoubleADTest, testSingleDimArcsineAdditionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::asin(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), 1.119769515, 1e-9);
  EXPECT_NEAR(f[2].d(0), 2.294157339, 1e-9);
}

TEST(longdoubleADTest, testSingleDimArcsineMultiplicationMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::asin(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), 0.2013579208, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.4082482905, 1e-9);
}

TEST(longdoubleADTest, testSingleDimArcsineDivisionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::asin(f[1] / f[0]);
  EXPECT_NEAR(f[2].v(), 0.927295218, 1e-9);
  EXPECT_NEAR(f[2].d(0), -2.66666666666, 1e-9);
}

TEST(longdoubleADTest, testSingleDimArcsineSubtractionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::asin(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 0.1001674212, 1e-9);
  EXPECT_NEAR(f[2].d(0), 1.005037815, 1e-9);
}

TEST(longdoubleADTest, testSingleDimArccosineAdditionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::acos(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), 0.4510268118, 1e-9);
  EXPECT_NEAR(f[2].d(0), -2.294157339, 1e-9);
}

TEST(longdoubleADTest, testSingleDimArccosineMultplicationMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::acos(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), 1.369438406, 1e-9);
  EXPECT_NEAR(f[2].d(0), -0.4082482905, 1e-9);
}

TEST(longdoubleADTest, testSingleDimArccosineDivisionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::acos(f[1] / f[0]);
  EXPECT_NEAR(f[2].v(), 0.6435011088, 1e-9);
  EXPECT_NEAR(f[2].d(0), 2.66666666666, 1e-9);
}

TEST(longdoubleADTest, testSingleDimArccosineSubtractionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::acos(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 1.470628906, 1e-9);
  EXPECT_NEAR(f[2].d(0), -1.005037815, 1e-9);
}

TEST(longdoubleADTest, testSingleDimArctangentAdditionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::atan(f[0] + f[1]);
  EXPECT_NEAR(f[2].v(), 0.7328151018, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.5524861878, 1e-9);
}

TEST(longdoubleADTest, testSingleDimArctangentMultiplicationMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::atan(f[0] * f[1]);
  EXPECT_NEAR(f[2].v(), 0.1973955598, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.3846153846, 1e-9);
}

TEST(longdoubleADTest, testSingleDimArctangentDivisionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::atan(f[1] / f[0]);
  EXPECT_NEAR(f[2].v(), 0.6747409422, 1e-9);
  EXPECT_NEAR(f[2].d(0), -0.9756097561, 1e-9);
}

TEST(longdoubleADTest, testSingleDimArctangentSubtractionMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 0.5;
  S x2 = 0.4;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = util::atan(f[0] - f[1]);
  EXPECT_NEAR(f[2].v(), 0.09966865249, 1e-9);
  EXPECT_NEAR(f[2].d(0), 0.99009900990099, 1e-9);
}

TEST(longdoubleADTest, testMultiDimMultiplicationMultipleVariables)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x1 = 10.5;
  S x2 = 11.25;
  T f[3];
  f[0] = x1;
  f[0].setDiffVariable(0);
  f[1] = x2;
  f[1].setDiffVariable(1);
  f[2] = f[0]*f[1];
  EXPECT_EQ(f[2].v(), 118.125);
  EXPECT_EQ(f[2].d(0), 11.25);
  EXPECT_EQ(f[2].d(1), 10.5);
}

//Test for certain functions

TEST(longdoubleADTest, testhasZeroDerivative)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  S x = 10.5;
  T f[2];
  f[0] = 5;
  f[1] = x;
  f[1].setDiffVariable(0);
  EXPECT_TRUE(f[0].hasZeroDerivative());
  EXPECT_FALSE(f[1].hasZeroDerivative());
}

TEST(longdoubleADTest, testEQADfADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  T y = 2.5;
  EXPECT_TRUE(x == y);
}

TEST(longdoubleADTest, testEQdoubleADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(2.5 == x);
}

TEST(longdoubleADTest, testEQADfdouble)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(x == 2.5);
}

TEST(longdoubleADTest, testNotEQADfADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x1 = 2.5;
  T x2 = 2.4;
  EXPECT_TRUE(x1 != x2);
}

TEST(longdoubleADTest, testNotEQdoubleADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(2.6 != x);
}

TEST(longdoubleADTest, testNotEQADfdouble)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(x != 2.6);
}

TEST(longdoubleADTest, testSmallerADfADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x1 = 2.5;
  T x2 = 2.4;
  EXPECT_TRUE(x1 > x2);
}

TEST(longdoubleADTest, testSmallerdoubleADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(2.6 > x);
}

TEST(longdoubleADTest, testSmallerADfdouble)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(x > 2.4);
}

TEST(longdoubleADTest, testSmallerEQADfADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x1 = 2.5;
  T x2 = 2.4;
  EXPECT_TRUE(x1 >= x2);
}

TEST(longdoubleADTest, testSmallerEQdoubleADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(2.5 >= x);
}

TEST(longdoubleADTest, testSmallerEQADfdouble)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(x >= 2.4);
}

TEST(longdoubleADTest, testBiggerADfADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x1 = 2.5;
  T x2 = 2.6;
  EXPECT_TRUE(x1 < x2);
}

TEST(longdoubleADTest, testBiggerdoubleADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(2.4 < x);
}

TEST(longdoubleADTest, testBiggerADfdouble)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(x < 2.6);
}

TEST(longdoubleADTest, testBiggerEQADfADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x1 = 2.5;
  T x2 = 2.6;
  EXPECT_TRUE(x1 <= x2);
}

TEST(longdoubleADTest, testBiggerEQdoubleADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(2.5 <= x);
}

TEST(longdoubleADTest, testBiggerEQADfdouble)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(x <= 2.6);
}

TEST(longdoubleADTest, testFloor)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(util::floor(x) == 2);
}

TEST(longdoubleADTest, testFabs1)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = -2.5;
  EXPECT_TRUE(util::fabs(x) == 2.5);
}

TEST(longdoubleADTest, testFabs2)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 5.25;
  EXPECT_TRUE(util::fabs(x) == 5.25);
}

TEST(longdoubleADTest, testFabs3)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 0;
  x.setDiffVariable(0);
  EXPECT_TRUE(util::fabs(x) == 0);
}

TEST(longdoubleADTest, testMaxADfdouble)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;

  EXPECT_TRUE(util::max(x, 5.25) == 5.25);
}

TEST(longdoubleADTest, testMaxADfADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x1 = 2.5;
  T x2 = 2.0;
  EXPECT_TRUE(util::max(x1, x2) == 2.5);
}

TEST(longdoubleADTest, testMaxdoubleADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(util::max<T>(7.25, x) == 7.25);
}

TEST(longdoubleADTest, testMinADfdouble)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(util::min<T>(x, 2.125) == 2.125);
}

TEST(longdoubleADTest, testMinADfADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x1 = 2.5;
  T x2 = 2.0;
  EXPECT_TRUE(util::min(x1, x2) == 2.0);
}

TEST(longdoubleADTest, testMindoubleADf)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 2.5;
  EXPECT_TRUE(util::min<T>(2.125, x) == 2.125);
}

TEST(longdoubleADTest, testNearZero)
{
  typedef long double S;
  typedef util::ADf<S,DIM2> T;
  T x = 1e-28;
  EXPECT_TRUE(nearZero(x));
}

// AD Tests utilizing type float

TEST(floatADTest, testSingleDimExpf)
{
  typedef float S;
  typedef util::ADf<S,DIM2> T;
  S x = 2.5;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = expf(f[0]);
  EXPECT_NEAR(f[1].v(), 12.1824939607035, 1e-6);
  EXPECT_NEAR(f[1].d(0), 12.1824939607035, 1e-6);
}

TEST(floatADTest, testSingleDimLogf)
{
  typedef float S;
  typedef util::ADf<S,DIM2> T;
  S x = 4.7;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = logf(f[0]);
  EXPECT_NEAR(f[1].v(), 1.54756250871601, 1e-6);
  EXPECT_NEAR(f[1].d(0), 0.212765957446809, 1e-6);
}

TEST(floatADTest, testSingleDimLog10f)
{
  typedef float S;
  typedef util::ADf<S,DIM2> T;
  S x = 53.7;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = log10f(f[0]);
  EXPECT_NEAR(f[1].v(), 1.7299742856995556068654105950281739286312009548692244, 1e-6);
  EXPECT_NEAR(f[1].d(0), 0.0080874205196136280754400171120410629849980820447610, 1e-6);
}

TEST(floatADTest, testSingleDimLog2f)
{
  typedef float S;
  typedef util::ADf<S,DIM2> T;
  S x = 94.6;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = log2f(f[0]);
  EXPECT_NEAR(f[1].v(), 6.5637682784520328472131852577044702148066953193247361, 1e-6);
  EXPECT_NEAR(f[1].d(0), 0.0152504761193336512405911699894491769284000629403064, 1e-6);
}

TEST(floatADTest, testSingleDimLog1Pf)
{
  typedef float S;
  typedef util::ADf<S,DIM2> T;
  S x = 4.9;
  T f[2];
  f[0] = x;
  f[0].setDiffVariable(0);
  f[1] = log1pf(f[0]);
  EXPECT_NEAR(f[1].v(), 1.7749523509116737665980589190353334164622453007016815, 1e-6);
  EXPECT_NEAR(f[1].d(0), 0.169491525423728813559322033898305084745762711864406779661, 1e-6);
}
