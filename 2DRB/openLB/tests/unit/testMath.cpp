/*  This file is part of the OpenLB library
*
*  Copyright (C) 2020 Luca Heim, Jan E. Marquardt
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

TEST(omath, testPow)
{
  int x = 3;
  float y1 = 2.5;
  float y2 = 3.25;
  double z1 = 2.5;
  double z2 = 3.25;
  long double k1 = 2.5;
  long double k2 = 3.25;

  EXPECT_EQ(util::pow(y1, y2), std::pow(y1, y2));
  EXPECT_EQ(typeid(util::pow(y1, y2)), typeid(std::pow(y1, y2)));
  EXPECT_EQ(util::powf(y1, y2), std::pow(y1, y2));
  EXPECT_EQ(typeid(util::powf(y1, y2)), typeid(std::pow(y1, y2)));

  EXPECT_EQ(util::pow(k1, k2), std::pow(k1, k2));
  EXPECT_EQ(typeid(util::pow(k1, k2)), typeid(std::pow(k1, k2)));
  EXPECT_EQ(util::pow(z1, z2), std::pow(z1, z2));
  EXPECT_EQ(typeid(util::pow(z1, z2)), typeid(std::pow(z1, z2)));
  EXPECT_EQ(util::pow(z1, x), std::pow(z1, x));
  EXPECT_EQ(typeid(util::pow(z1, x)), typeid(std::pow(z1, x)));
  EXPECT_EQ(util::pow(y1, x), std::pow(y1, x));
  EXPECT_EQ(typeid(util::pow(y1, x)), typeid(std::pow(y1, x)));
  EXPECT_EQ(util::pow(k1, x), std::pow(k1, x));
  EXPECT_EQ(typeid(util::pow(k1, x)), typeid(std::pow(k1, x)));

  std::cout << typeid(util::pow(y1, x)).name() << std::endl;

  EXPECT_EQ(util::pow(x, x), std::pow(x, x));
  EXPECT_EQ(typeid(util::pow(x, x)), typeid(std::pow(x, x)));
}

TEST(omath, testExp)
{
  int x = 3;
  float y = 2.5;
  double z = 2.5;
  long double k = 2.5;

  EXPECT_EQ(util::exp(y), std::exp(y));
  EXPECT_EQ(typeid(util::exp(y)), typeid(std::exp(y)));
  EXPECT_EQ(util::expf(y), std::exp(y));
  EXPECT_EQ(typeid(util::expf(y)), typeid(std::exp(y)));

  EXPECT_EQ(util::exp(z), std::exp(z));
  EXPECT_EQ(typeid(util::exp(z)), typeid(std::exp(z)));

  EXPECT_EQ(util::exp(k), std::exp(k));
  EXPECT_EQ(typeid(util::exp(k)), typeid(std::exp(k)));

  EXPECT_EQ(util::exp(x), std::exp(x));
  EXPECT_EQ(typeid(util::exp(x)), typeid(std::exp(x)));
}

TEST(omath, testLog)
{
  int x = 3;
  float y = 2.5;
  double z = 2.5;
  long double k = 2.5;

  EXPECT_EQ(util::log(y), std::log(y));
  EXPECT_EQ(typeid(util::log(y)), typeid(std::log(y)));
  EXPECT_EQ(util::logf(y), std::log(y));
  EXPECT_EQ(typeid(util::logf(y)), typeid(std::log(y)));

  EXPECT_EQ(util::log(z), std::log(z));
  EXPECT_EQ(typeid(util::log(z)), typeid(std::log(z)));

  EXPECT_EQ(util::log(k), std::log(k));
  EXPECT_EQ(typeid(util::log(k)), typeid(std::log(k)));
  EXPECT_EQ(util::log(x), std::log(x));
  EXPECT_EQ(typeid(util::log(x)), typeid(std::log(x)));
}

TEST(omath, testLog10)
{
  int x = 3;
  float y = 2.5;
  double z = 2.5;
  long double k = 2.5;

  EXPECT_EQ(util::log10(y), std::log10(y));
  EXPECT_EQ(typeid(util::log10(y)), typeid(std::log10(y)));
  EXPECT_EQ(util::log10f(y), std::log10(y));
  EXPECT_EQ(typeid(util::log10f(y)), typeid(std::log10(y)));

  EXPECT_EQ(util::log10(z), std::log10(z));
  EXPECT_EQ(typeid(util::log10(z)), typeid(std::log10(z)));

  EXPECT_EQ(util::log10(k), std::log10(k));
  EXPECT_EQ(typeid(util::log10(k)), typeid(std::log10(k)));

  EXPECT_EQ(util::log10(x), std::log10(x));
  EXPECT_EQ(typeid(util::log10(x)), typeid(std::log10(x)));
}

TEST(omath, testLog2)
{
  int x = 3;
  float y = 2.5;
  double z = 2.5;
  long double k = 2.5;

  EXPECT_EQ(util::log2(y), std::log2(y));
  EXPECT_EQ(typeid(util::log2(y)), typeid(std::log2(y)));
  EXPECT_EQ(util::log2f(y), std::log2(y));
  EXPECT_EQ(typeid(util::log2f(y)), typeid(std::log2(y)));

  EXPECT_EQ(util::log2(z), std::log2(z));
  EXPECT_EQ(typeid(util::log2(z)), typeid(std::log2(z)));

  EXPECT_EQ(util::log2(k), std::log2(k));
  EXPECT_EQ(typeid(util::log2(k)), typeid(std::log2(k)));

  EXPECT_EQ(util::log2(x), std::log2(x));
  EXPECT_EQ(typeid(util::log2(x)), typeid(std::log2(x)));
}

TEST(omath, testLog1p)
{
  int x = 3;
  float y = 2.5;
  double z = 2.5;
  long double k = 2.5;

  EXPECT_EQ(util::log1p(y), std::log1p(y));
  EXPECT_EQ(typeid(util::log1p(y)), typeid(std::log1p(y)));
  EXPECT_EQ(util::log1pf(y), std::log1p(y));
  EXPECT_EQ(typeid(util::log1pf(y)), typeid(std::log1p(y)));

  EXPECT_EQ(util::log1p(z), std::log1p(z));
  EXPECT_EQ(typeid(util::log1p(z)), typeid(std::log1p(z)));

  EXPECT_EQ(util::log1p(k), std::log1p(k));
  EXPECT_EQ(typeid(util::log1p(k)), typeid(std::log1p(k)));
  EXPECT_EQ(util::log1pl(k), std::log1p(k));
  EXPECT_EQ(typeid(util::log1pl(k)), typeid(std::log1p(k)));

  EXPECT_EQ(util::log1p(x), std::log1p(x));
  EXPECT_EQ(typeid(util::log1p(x)), typeid(std::log1p(x)));
}

TEST(omath, testSqrt)
{
  int x = 3;
  float y = 2.5;
  double z = 2.5;
  long double k = 2.5;

  EXPECT_EQ(util::sqrt(y), std::sqrt(y));
  EXPECT_EQ(typeid(util::sqrt(y)), typeid(std::sqrt(y)));
  EXPECT_EQ(util::sqrtf(y), std::sqrt(y));
  EXPECT_EQ(typeid(util::sqrtf(y)), typeid(std::sqrt(y)));

  EXPECT_EQ(util::sqrt(z), std::sqrt(z));
  EXPECT_EQ(typeid(util::sqrt(z)), typeid(std::sqrt(z)));

  EXPECT_EQ(util::sqrt(k), std::sqrt(k));
  EXPECT_EQ(typeid(util::sqrt(k)), typeid(std::sqrt(k)));

  EXPECT_EQ(util::sqrt(x), std::sqrt(x));
  EXPECT_EQ(typeid(util::sqrt(x)), typeid(std::sqrt(x)));
}

TEST(omath, testSin)
{
  int x = 3;
  float y = 2.5;
  double z = 2.5;
  long double k = 2.5;

  EXPECT_EQ(util::sin(y), std::sin(y));
  EXPECT_EQ(typeid(util::sin(y)), typeid(std::sin(y)));
  EXPECT_EQ(util::sinf(y), std::sin(y));
  EXPECT_EQ(typeid(util::sinf(y)), typeid(std::sin(y)));

  EXPECT_EQ(util::sin(z), std::sin(z));
  EXPECT_EQ(typeid(util::sin(z)), typeid(std::sin(z)));

  EXPECT_EQ(util::sin(k), std::sin(k));
  EXPECT_EQ(typeid(util::sin(k)), typeid(std::sin(k)));
  EXPECT_EQ(util::sinl(k), std::sin(k));
  EXPECT_EQ(typeid(util::sinl(k)), typeid(std::sin(k)));

  EXPECT_EQ(util::sin(x), std::sin(x));
  EXPECT_EQ(typeid(util::sin(x)), typeid(std::sin(x)));
}

TEST(omath, testSinh)
{
  int x = 3;
  float y = 2.5;
  double z = 2.5;
  long double k = 2.5;

  EXPECT_EQ(util::sinh(y), std::sinh(y));
  EXPECT_EQ(typeid(util::sinh(y)), typeid(std::sinh(y)));
  EXPECT_EQ(util::sinhf(y), std::sinh(y));
  EXPECT_EQ(typeid(util::sinhf(y)), typeid(std::sinh(y)));

  EXPECT_EQ(util::sinh(z), std::sinh(z));
  EXPECT_EQ(typeid(util::sinh(z)), typeid(std::sinh(z)));

  EXPECT_EQ(util::sinh(k), std::sinh(k));
  EXPECT_EQ(typeid(util::sinh(k)), typeid(std::sinh(k)));
  EXPECT_EQ(util::sinhl(k), std::sinh(k));
  EXPECT_EQ(typeid(util::sinhl(k)), typeid(std::sinh(k)));

  EXPECT_EQ(util::sinh(x), std::sinh(x));
  EXPECT_EQ(typeid(util::sinh(x)), typeid(std::sinh(x)));
}

TEST(omath, testCos)
{
  int x = 3;
  float y = 2.5;
  double z = 2.5;
  long double k = 2.5;

  EXPECT_EQ(util::cos(y), std::cos(y));
  EXPECT_EQ(typeid(util::cos(y)), typeid(std::cos(y)));
  EXPECT_EQ(util::cosf(y), std::cos(y));
  EXPECT_EQ(typeid(util::cosf(y)), typeid(std::cos(y)));

  EXPECT_EQ(util::cos(z), std::cos(z));
  EXPECT_EQ(typeid(util::cos(z)), typeid(std::cos(z)));

  EXPECT_EQ(util::cos(k), std::cos(k));
  EXPECT_EQ(typeid(util::cos(k)), typeid(std::cos(k)));
  EXPECT_EQ(util::cosl(k), std::cos(k));
  EXPECT_EQ(typeid(util::cosl(k)), typeid(std::cos(k)));

  EXPECT_EQ(util::cos(x), std::cos(x));
  EXPECT_EQ(typeid(util::cos(x)), typeid(std::cos(x)));
}

TEST(omath, testCosh)
{
  int x = 3;
  float y = 2.5;
  double z = 2.5;
  long double k = 2.5;

  EXPECT_EQ(util::cosh(y), std::cosh(y));
  EXPECT_EQ(typeid(util::cosh(y)), typeid(std::cosh(y)));
  EXPECT_EQ(util::coshf(y), std::cosh(y));
  EXPECT_EQ(typeid(util::coshf(y)), typeid(std::cosh(y)));

  EXPECT_EQ(util::cosh(z), std::cosh(z));
  EXPECT_EQ(typeid(util::cosh(z)), typeid(std::cosh(z)));

  EXPECT_EQ(util::cosh(k), std::cosh(k));
  EXPECT_EQ(typeid(util::cosh(k)), typeid(std::cosh(k)));
  EXPECT_EQ(util::coshl(k), std::cosh(k));
  EXPECT_EQ(typeid(util::coshl(k)), typeid(std::cosh(k)));

  EXPECT_EQ(util::cosh(x), std::cosh(x));
  EXPECT_EQ(typeid(util::cosh(x)), typeid(std::cosh(x)));
}

TEST(omath, testTan)
{
  int x = 3;
  float y = 2.5;
  double z = 2.5;
  long double k = 2.5;

  EXPECT_EQ(util::tan(y), std::tan(y));
  EXPECT_EQ(typeid(util::tan(y)), typeid(std::tan(y)));
  EXPECT_EQ(util::tanf(y), std::tan(y));
  EXPECT_EQ(typeid(util::tanf(y)), typeid(std::tan(y)));

  EXPECT_EQ(util::tan(z), std::tan(z));
  EXPECT_EQ(typeid(util::tan(z)), typeid(std::tan(z)));

  EXPECT_EQ(util::tan(k), std::tan(k));
  EXPECT_EQ(typeid(util::tan(k)), typeid(std::tan(k)));
  EXPECT_EQ(util::tanl(k), std::tan(k));
  EXPECT_EQ(typeid(util::tanl(k)), typeid(std::tan(k)));

  EXPECT_EQ(util::tan(x), std::tan(x));
  EXPECT_EQ(typeid(util::tan(x)), typeid(std::tan(x)));
}

TEST(omath, testTanh)
{
  int x = 3;
  float y = 2.5;
  double z = 2.5;
  long double k = 2.5;

  EXPECT_EQ(util::tanh(y), std::tanh(y));
  EXPECT_EQ(typeid(util::tanh(y)), typeid(std::tanh(y)));
  EXPECT_EQ(util::tanhf(y), std::tanh(y));
  EXPECT_EQ(typeid(util::tanhf(y)), typeid(std::tanh(y)));

  EXPECT_EQ(util::tanh(z), std::tanh(z));
  EXPECT_EQ(typeid(util::tanh(z)), typeid(std::tanh(z)));

  EXPECT_EQ(util::tanh(k), std::tanh(k));
  EXPECT_EQ(typeid(util::tanh(k)), typeid(std::tanh(k)));
  EXPECT_EQ(util::tanhl(k), std::tanh(k));
  EXPECT_EQ(typeid(util::tanhl(k)), typeid(std::tanh(k)));

  EXPECT_EQ(util::tanh(x), std::tanh(x));
  EXPECT_EQ(typeid(util::tanh(x)), typeid(std::tanh(x)));
}

TEST(omath, testAsin)
{
  int x = 1;
  float y = 0.5;
  double z = 0.5;
  long double k = 0.5;

  EXPECT_EQ(util::asin(y), std::asin(y));
  EXPECT_EQ(typeid(util::asin(y)), typeid(std::asin(y)));
  EXPECT_EQ(util::asinf(y), std::asin(y));
  EXPECT_EQ(typeid(util::asinf(y)), typeid(std::asin(y)));

  EXPECT_EQ(util::asin(z), std::asin(z));
  EXPECT_EQ(typeid(util::asin(z)), typeid(std::asin(z)));

  EXPECT_EQ(util::asin(k), std::asin(k));
  EXPECT_EQ(typeid(util::asin(k)), typeid(std::asin(k)));
  EXPECT_EQ(util::asinl(k), std::asin(k));
  EXPECT_EQ(typeid(util::asinl(k)), typeid(std::asin(k)));

  EXPECT_EQ(util::asin(x), std::asin(x));
  EXPECT_EQ(typeid(util::asin(x)), typeid(std::asin(x)));
}

TEST(omath, testAsinh)
{
  int x = 1;
  float y = 0.5;
  double z = 0.5;
  long double k = 0.5;

  EXPECT_EQ(util::asinh(y), std::asinh(y));
  EXPECT_EQ(typeid(util::asinh(y)), typeid(std::asinh(y)));
  EXPECT_EQ(util::asinhf(y), std::asinh(y));
  EXPECT_EQ(typeid(util::asinhf(y)), typeid(std::asinh(y)));

  EXPECT_EQ(util::asinh(z), std::asinh(z));
  EXPECT_EQ(typeid(util::asinh(z)), typeid(std::asinh(z)));

  EXPECT_EQ(util::asinh(k), std::asinh(k));
  EXPECT_EQ(typeid(util::asinh(k)), typeid(std::asinh(k)));
  EXPECT_EQ(util::asinhl(k), std::asinh(k));
  EXPECT_EQ(typeid(util::asinhl(k)), typeid(std::asinh(k)));

  EXPECT_EQ(util::asinh(x), std::asinh(x));
  EXPECT_EQ(typeid(util::asinh(x)), typeid(std::asinh(x)));
}

TEST(omath, testAcos)
{
  int x = 1;
  float y = 0.5;
  double z = 0.5;
  long double k = 0.5;

  EXPECT_EQ(util::acos(y), std::acos(y));
  EXPECT_EQ(typeid(util::acos(y)), typeid(std::acos(y)));
  EXPECT_EQ(util::acosf(y), std::acos(y));
  EXPECT_EQ(typeid(util::acosf(y)), typeid(std::acos(y)));

  EXPECT_EQ(util::acos(z), std::acos(z));
  EXPECT_EQ(typeid(util::acos(z)), typeid(std::acos(z)));

  EXPECT_EQ(util::acos(k), std::acos(k));
  EXPECT_EQ(typeid(util::acos(k)), typeid(std::acos(k)));
  EXPECT_EQ(util::acosl(k), std::acos(k));
  EXPECT_EQ(typeid(util::acosl(k)), typeid(std::acos(k)));

  EXPECT_EQ(util::acos(x), std::acos(x));
  EXPECT_EQ(typeid(util::acos(x)), typeid(std::acos(x)));
}

TEST(omath, testAcosh)
{
  int x = 1;
  float y = 1.5;
  double z = 1.5;
  long double k = 1.5;

  EXPECT_EQ(util::acosh(y), std::acosh(y));
  EXPECT_EQ(typeid(util::acosh(y)), typeid(std::acosh(y)));
  EXPECT_EQ(util::acoshf(y), std::acosh(y));
  EXPECT_EQ(typeid(util::acoshf(y)), typeid(std::acosh(y)));

  EXPECT_EQ(util::acosh(z), std::acosh(z));
  EXPECT_EQ(typeid(util::acosh(z)), typeid(std::acosh(z)));

  EXPECT_EQ(util::acosh(k), std::acosh(k));
  EXPECT_EQ(typeid(util::acosh(k)), typeid(std::acosh(k)));
  EXPECT_EQ(util::acoshl(k), std::acosh(k));
  EXPECT_EQ(typeid(util::acoshl(k)), typeid(std::acosh(k)));

  EXPECT_EQ(util::acosh(x), std::acosh(x));
  EXPECT_EQ(typeid(util::acosh(x)), typeid(std::acosh(x)));
}

TEST(omath, testAtan)
{
  int x = 1;
  float y = 0.5;
  double z = 0.5;
  long double k = 0.5;

  EXPECT_EQ(util::atan(y), std::atan(y));
  EXPECT_EQ(typeid(util::atan(y)), typeid(std::atan(y)));
  EXPECT_EQ(util::atanf(y), std::atan(y));
  EXPECT_EQ(typeid(util::atanf(y)), typeid(std::atan(y)));

  EXPECT_EQ(util::atan(z), std::atan(z));
  EXPECT_EQ(typeid(util::atan(z)), typeid(std::atan(z)));

  EXPECT_EQ(util::atan(k), std::atan(k));
  EXPECT_EQ(typeid(util::atan(k)), typeid(std::atan(k)));
  EXPECT_EQ(util::atanl(k), std::atan(k));
  EXPECT_EQ(typeid(util::atanl(k)), typeid(std::atan(k)));

  EXPECT_EQ(util::atan(x), std::atan(x));
  EXPECT_EQ(typeid(util::atan(x)), typeid(std::atan(x)));
}

TEST(omath, testAtanh)
{
  int x = 1;
  float y = 0.5;
  double z = 0.5;
  long double k = 0.5;

  EXPECT_EQ(util::atanh(y), std::atanh(y));
  EXPECT_EQ(typeid(util::atanh(y)), typeid(std::atanh(y)));
  EXPECT_EQ(util::atanhf(y), std::atanh(y));
  EXPECT_EQ(typeid(util::atanhf(y)), typeid(std::atanh(y)));

  EXPECT_EQ(util::atanh(z), std::atanh(z));
  EXPECT_EQ(typeid(util::atanh(z)), typeid(std::atanh(z)));

  EXPECT_EQ(util::atanh(k), std::atanh(k));
  EXPECT_EQ(typeid(util::atanh(k)), typeid(std::atanh(k)));
  EXPECT_EQ(util::atanhl(k), std::atanh(k));
  EXPECT_EQ(typeid(util::atanhl(k)), typeid(std::atanh(k)));

  EXPECT_EQ(util::atanh(x), std::atanh(x));
  EXPECT_EQ(typeid(util::atanh(x)), typeid(std::atanh(x)));
}

TEST(omath, testFloor)
{
  int x = 2;
  float y = 2.5;
  double z = 2.5;
  long double k = 2;

  EXPECT_EQ(util::floor(y), std::floor(y));
  EXPECT_EQ(typeid(util::floor(y)), typeid(std::floor(y)));
  EXPECT_EQ(util::floorf(y), std::floor(y));
  EXPECT_EQ(typeid(util::floorf(y)), typeid(std::floor(y)));

  EXPECT_EQ(util::floor(z), std::floor(z));
  EXPECT_EQ(typeid(util::floor(z)), typeid(std::floor(z)));

  EXPECT_EQ(util::floor(k), std::floor(k));
  EXPECT_EQ(typeid(util::floor(k)), typeid(std::floor(k)));
  EXPECT_EQ(util::floorl(k), std::floor(k));
  EXPECT_EQ(typeid(util::floorl(k)), typeid(std::floor(k)));

  EXPECT_EQ(util::floor(x), std::floor(x));
  EXPECT_EQ(typeid(util::floor(x)), typeid(std::floor(x)));
}

TEST(omath, testCeil)
{
  int x = 2;
  float y = 2.5;
  double z = 2.5;
  long double k = 2.5;

  EXPECT_EQ(util::ceil(y), std::ceil(y));
  EXPECT_EQ(typeid(util::ceil(y)), typeid(std::ceil(y)));
  EXPECT_EQ(util::ceilf(y), std::ceil(y));
  EXPECT_EQ(typeid(util::ceilf(y)), typeid(std::ceil(y)));
  EXPECT_EQ(util::ceil(z), std::ceil(z));
  EXPECT_EQ(typeid(util::ceil(z)), typeid(std::ceil(z)));
  EXPECT_EQ(util::ceil(k), std::ceil(k));
  EXPECT_EQ(typeid(util::ceil(k)), typeid(std::ceil(k)));
  EXPECT_EQ(util::ceil(x), std::ceil(x));
  EXPECT_EQ(typeid(util::ceil(x)), typeid(std::ceil(x)));
}

TEST(omath, testTrunc)
{
  int x = 2;
  float y = 2.5;
  double z = 2.5;
  long double k = 2;

  EXPECT_EQ(util::trunc(y), std::trunc(y));
  EXPECT_EQ(typeid(util::trunc(y)), typeid(std::trunc(y)));
  EXPECT_EQ(util::truncf(y), std::trunc(y));
  EXPECT_EQ(typeid(util::truncf(y)), typeid(std::trunc(y)));
  EXPECT_EQ(util::trunc(z), std::trunc(z));
  EXPECT_EQ(typeid(util::trunc(z)), typeid(std::trunc(z)));
  EXPECT_EQ(util::trunc(k), std::trunc(k));
  EXPECT_EQ(typeid(util::trunc(k)), typeid(std::trunc(k)));
  EXPECT_EQ(util::trunc(x), std::trunc(x));
  EXPECT_EQ(typeid(util::trunc(x)), typeid(std::trunc(x)));
}

TEST(omath, testRound)
{
  int x = 2;
  float y = 2.3;
  long double k = 2.7;

  EXPECT_EQ(util::round(x), std::round(x));
  EXPECT_EQ(typeid(util::round(x)), typeid(std::round(x)));
  EXPECT_EQ(util::roundf(y), std::round(y));
  EXPECT_EQ(typeid(util::roundf(y)), typeid(std::round(y)));
}

TEST(omath, testAbsolute)
{
  int x = -2;
  long int m = -4;
  long long int l = -3;
  float y = -2.3;
  double z = -2.3;
  long double k = -2.7;

  EXPECT_EQ(util::abs(x), std::abs(x));
  EXPECT_EQ(typeid(util::abs(x)), typeid(std::abs(x)));
  EXPECT_EQ(util::fabs(y), std::fabs(y));
  EXPECT_EQ(typeid(util::fabs(y)), typeid(std::fabs(y)));
  EXPECT_EQ(util::fabsf(y), std::fabs(y));
  EXPECT_EQ(typeid(util::fabsf(y)), typeid(std::fabs(y)));
  EXPECT_EQ(util::fabs(z), std::fabs(z));
  EXPECT_EQ(typeid(util::fabs(z)), typeid(std::fabs(z)));
  EXPECT_EQ(util::fabs(k), std::fabs(k));
  EXPECT_EQ(typeid(util::fabs(k)), typeid(std::fabs(k)));
  EXPECT_EQ(util::fabs(x), std::fabs(x));
  EXPECT_EQ(typeid(util::fabs(x)), typeid(std::fabs(x)));
  EXPECT_EQ(util::fabs(m), std::fabs(m));
  EXPECT_EQ(typeid(util::fabs(m)), typeid(std::fabs(m)));
  EXPECT_EQ(util::fabs(l), std::fabs(l));
  EXPECT_EQ(typeid(util::fabs(l)), typeid(std::fabs(l)));
}
