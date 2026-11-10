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


#include <iostream>
#include <gtest/gtest.h>
#include "helper/unitTestHelper.h"

#include <olb.h>

using namespace olb;
using namespace olb::util;


TEST(doubleL1Test, testConst)
{
  using T = double;
  TimeIntegrator<T> integrator(0.2, 3.2, 0.5);

  // Result is initialized to zero
  EXPECT_DOUBLE_EQ(integrator.getResult(), T(0));

  // Integrator ignores values outside the interval
  integrator.takeValue(0, 10.0);
  EXPECT_DOUBLE_EQ(integrator.getResult(), T(0));
  integrator.takeValue(11, 10.0);
  EXPECT_DOUBLE_EQ(integrator.getResult(), T(0));

  // Correct handling of the first time step
  integrator.takeValue(1, 1.0);
  EXPECT_DOUBLE_EQ(integrator.getResult(), 0.55);  // 0.25 + (2.0 + 0.3 / 0.5) * 0.3 / 2.0); // result for linear extrapolation
  // Correct handling of the second time step
  integrator.takeValue(2, 1.0);
  EXPECT_DOUBLE_EQ(integrator.getResult(), 1.25 - 0.2);

  // Correct handling of intermediate time steps
  integrator.takeValue(3, 1.0);
  EXPECT_DOUBLE_EQ(integrator.getResult(), 1.75 - 0.2);
  integrator.takeValue(4, 1.0);
  EXPECT_DOUBLE_EQ(integrator.getResult(), 2.25 - 0.2);

  // Correct handling of the prelast time step
  integrator.takeValue(5, 1.0);
  EXPECT_DOUBLE_EQ(integrator.getResult(), 2.55);  // 2.75 - 0.2 - 0.2 * 0.2);  // see above
  // Correct handling of the last time step
  integrator.takeValue(6, 1.0);
  EXPECT_DOUBLE_EQ(integrator.getResult(), 3.2 - 0.2);
}

TEST(doubleL1Test, testLinear)
{
  using T = double;
  TimeIntegrator<T> integrator(0.2, 3.2, 0.5);

  // Result is initialized to zero
  EXPECT_DOUBLE_EQ(integrator.getResult(), T(0));

  // Integrator ignores values outside the interval
  integrator.takeValue(0, 10.0);
  EXPECT_DOUBLE_EQ(integrator.getResult(), T(0));
  integrator.takeValue(11, 10.0);
  EXPECT_DOUBLE_EQ(integrator.getResult(), T(0));

  // Correct handling of the first time step
  integrator.takeValue(1, 1.5);
  EXPECT_DOUBLE_EQ(integrator.getResult(), 1.5 * 0.55);  // 1.5 * (0.25 + (2.0 + 0.3 / 0.5) * 0.3 / 2.0));
  // Correct handling of the second time step
  integrator.takeValue(2, 2.0);
  EXPECT_DOUBLE_EQ(integrator.getResult(), 1.5 * 0.55 + 2.0 * 0.5);  // 1.8 - 0.02);

  // Correct handling of intermediate time steps
  integrator.takeValue(3, 2.5);
  EXPECT_DOUBLE_EQ(integrator.getResult(), 1.5*0.55 + 0.5*(2.0+2.5));  // 0.5*1.5*1.5 + 1.5 - 0.02 - 0.2 + 2.5 * 0.25);
  integrator.takeValue(4, 3.0);
  EXPECT_DOUBLE_EQ(integrator.getResult(), 1.5*0.55 + 0.5*(2.0+2.5+3.0));  // 0.5*2.0*2.0 + 2.0 - 0.02 - 0.2 + 3.0 * 0.25);

  // Correct handling of the prelast time step
  integrator.takeValue(5, 3.5);
  EXPECT_DOUBLE_EQ(integrator.getResult(), 1.5*0.55 + 0.5*(2.0+2.5+3.0+3.5));  // 0.5*2.5*2.5 + 2.5 - 0.02 - 0.2 + 3.5 * 0.25 - 3.5 * 0.2 * 0.2);
  // Correct handling of the last time step
  integrator.takeValue(6, 4.0);
  EXPECT_DOUBLE_EQ(integrator.getResult(), 1.5*0.55 + 0.5*(2.0+2.5+3.0+3.5) + 4.0*0.45);  // 0.5*3.2*3.2 + 3.2 - 0.02 - 0.2);
}

TEST(doubleL1Test, testSinus)
{
  using T = double;
  const T dt = 1.0e-5;
  TimeIntegrator<T> integrator(-0.2, 2.5, dt);

  // Result is initialized to zero
  EXPECT_DOUBLE_EQ(integrator.getResult(), T(0));

  // Integrator ignores values outside the interval
  integrator.takeValue(int(-1e+8), 10.0);
  EXPECT_DOUBLE_EQ(integrator.getResult(), T(0));
  integrator.takeValue(int(12e+8), 10.0);
  EXPECT_DOUBLE_EQ(integrator.getResult(), T(0));

  for (int iT(-1.0 / dt); iT <= int(6.0 / dt); iT++) {
    integrator.takeValue(iT, util::sin(T(iT) * dt));
  }
  EXPECT_NEAR(integrator.getResult(), - util::cos(2.5) + util::cos(-0.2), 1.0e-10);
}

TEST(doubleL2Test, testLinear)
{
  using T = double;
  const T dt = 1.0e-5;
  TimeIntegrator<T,2> integrator(-0.2, 2.5, dt);

  // Result is initialized to zero
  EXPECT_DOUBLE_EQ(integrator.getResult(), T(0));

  for (int iT(-1.0 / dt); iT <= int(6.0 / dt); iT++) {
    integrator.takeValue(iT, T(iT) * dt);
  }
  EXPECT_NEAR(
    integrator.getResult(),
    util::sqrt( (util::pow(2.5, 3) - util::pow(-0.2, 3)) / 3),
    1.0e-10);
}

TEST(doubleL2Test, testSquare)
{
  using T = double;
  const T dt = 1.0e-5;
  TimeIntegrator<T,2> integrator(-0.4, 2.0, dt);

  // Result is initialized to zero
  EXPECT_DOUBLE_EQ(integrator.getResult(), T(0));

  for (int iT(-1.0 / dt); iT <= int(3.0 / dt); iT++) {
    integrator.takeValue(iT, T(iT) * dt * T(iT) * dt);
  }
  EXPECT_NEAR(
    integrator.getResult(),
    util::sqrt( (util::pow(2.0, 5) - util::pow(-0.4, 5)) / 5),
    1.0e-9);
}

TEST(doubleLinfTest, testSquare)
{
  using T = double;
  const T dt = 1.0e-6;
  TimeIntegrator<T,0> integrator(-0.4, 2.0, dt);

  // Result is initialized to zero
  EXPECT_DOUBLE_EQ(integrator.getResult(), T(0));

  for (int iT(-1.0 / dt); iT <= int(3.0 / dt); iT++) {
    integrator.takeValue(iT, T(iT) * dt * T(iT) * dt);
  }
  EXPECT_NEAR(integrator.getResult(), 2.0 * 2.0, 1.0e-5);
}

TEST(doubleLinfTest, testSinus)
{
  using T = double;
  const T dt = 1.0e-6;
  TimeIntegrator<T,0> integrator(-0.1, 2.0, dt);

  // Result is initialized to zero
  EXPECT_DOUBLE_EQ(integrator.getResult(), T(0));

  for (int iT(-1.0 / dt); iT <= int(3.0 / dt); iT++) {
    integrator.takeValue(iT, util::sin(T(iT) * dt));
  }
  EXPECT_NEAR(integrator.getResult(), 1.0, 1.0e-5);
}

TEST(adfL1Test, testConst1)
{
  using S = double;
  using T = ADf<S,1>;
  const S dt = 1.0e-1;

  TimeIntegrator<T> integrator(0.2, 3.2, dt);

  // Result is initialized to zero
  EXPECT_DOUBLE_EQ(integrator.getResult().v(), S(0));
  EXPECT_NEAR(integrator.getResult().d(0), S(0), 1.0e-16);

  for (int iT(-1.0 / dt); iT <= int(4.0 / dt); iT++) {
    integrator.takeValue(iT, T(1.0, {{ 1.0 }}));
  }
  EXPECT_NEAR(integrator.getResult().v(), 3.2 - 0.2, 1.0e-10);
  EXPECT_NEAR(integrator.getResult().d(0), 3.2 - 0.2, 1.0e-10);
}

TEST(adfL1Test, testConst2)
{
  using S = double;
  using T = ADf<S,1>;
  const S dt = 1.0e-1;

  T upperBound(3.2, {{1.0}});
  TimeIntegrator<T> integrator(0.2, upperBound, dt);

  // Result is initialized to zero
  EXPECT_DOUBLE_EQ(integrator.getResult().v(), S(0));
  EXPECT_NEAR(integrator.getResult().d(0), S(0), 1.0e-16);

  for (int iT(-1.0 / dt); iT <= upperBound / dt + 1; iT++) {
    integrator.takeValue(iT, T(1));
  }
  EXPECT_NEAR(integrator.getResult().v(), 3.2 - 0.2, 1.0e-10);
  EXPECT_NEAR(integrator.getResult().d(0), 1.0, 1.0e-10);
}

TEST(adfL1Test, testAverage)
{
  using S = double;
  using T = ADf<S,2>;
  const S dt = 1.0e-1;

  T lowerBound(0.2, {{1.0, 0.0}});
  T upperBound(3.2, {{0.0, 1.0}});
  TimeIntegrator<T> integrator(lowerBound, upperBound, dt);

  // Result is initialized to zero
  EXPECT_DOUBLE_EQ(integrator.getResult().v(), S(0));
  EXPECT_NEAR(integrator.getResult().d(0), S(0), 1.0e-16);
  EXPECT_NEAR(integrator.getResult().d(1), S(0), 1.0e-16);

  for (int iT(-1.0 / dt); iT <= upperBound / dt + 2; iT++) {
    integrator.takeValue(iT, T(1) / (upperBound - lowerBound));
  }
  EXPECT_NEAR(integrator.getResult().v(), 1.0, 1.0e-10);
  EXPECT_NEAR(integrator.getResult().d(0), 0.0, 1.0e-10);
  EXPECT_NEAR(integrator.getResult().d(1), 0.0, 1.0e-10);
}

TEST(adfL1Test, testSinus1)
{
  // Derive an integral by its limits via fund. thm. of calculus
  using S = double;
  using T = ADf<S,2>;
  const S dt = 1.0e-5;

  T lowerBound(0.2, {{1.0, 0.0}});
  T upperBound(3.2, {{0.0, 1.0}});
  TimeIntegrator<T> integrator(lowerBound, upperBound, dt);

  // Result is initialized to zero
  EXPECT_DOUBLE_EQ(integrator.getResult().v(), S(0));
  EXPECT_NEAR(integrator.getResult().d(0), S(0), 1.0e-16);
  EXPECT_NEAR(integrator.getResult().d(1), S(0), 1.0e-16);

  for (int iT(-1.0 / dt); iT <= upperBound / dt + 2; iT++) {
    integrator.takeValue(iT, util::sin(iT * dt));
  }
  EXPECT_NEAR(integrator.getResult().v(), - util::cos(upperBound) + util::cos(lowerBound), 1.0e-10);
  EXPECT_NEAR(integrator.getResult().d(0), - util::sin(lowerBound), 1.0e-5);
  EXPECT_NEAR(integrator.getResult().d(1), util::sin(upperBound), 1.0e-10);
}

TEST(adfL1Test, testSinus2)
{
  // Derive a parameter dependent integral
  using S = double;
  using T = ADf<S,2>;
  const S dt = 1.0e-5;

  T p1(0.5, {{1.0, 0.0}});
  T p2(-2.9, {{0.0, 1.0}});
  T lowerBound(0.2);
  T upperBound(3.2);
  TimeIntegrator<T> integrator(lowerBound, upperBound, dt);

  // Result is initialized to zero
  EXPECT_DOUBLE_EQ(integrator.getResult().v(), S(0));
  EXPECT_NEAR(integrator.getResult().d(0), S(0), 1.0e-16);
  EXPECT_NEAR(integrator.getResult().d(1), S(0), 1.0e-16);

  for (int iT(-1.0 / dt); iT <= upperBound / dt + 2; iT++) {
    integrator.takeValue(iT, util::sin(iT * dt * p1 + p1 * p2));
  }
  EXPECT_NEAR(
    integrator.getResult().v(),
    - util::cos(p1*upperBound + p1*p2) / p1
    + util::cos(p1*lowerBound + p1*p2) / p1,
    1.0e-10);
  EXPECT_NEAR(
    integrator.getResult().d(0),
    util::sin(p1*upperBound + p1*p2) * (upperBound+p2) / p1
    - util::sin(p1*lowerBound + p1*p2) * (lowerBound+p2) / p1
    + util::cos(p1*upperBound + p1*p2) / (p1*p1)
    - util::cos(p1*lowerBound + p1*p2) / (p1*p1),
    1.0e-10);
  EXPECT_NEAR(
    integrator.getResult().d(1),
    util::sin(p1*upperBound + p1*p2)
    - util::sin(p1*lowerBound + p1*p2),
    1.0e-10);
}

TEST(adfL1Test, testSinus3)
{
  // Both integrand and limits depend on the derivation variable
  using S = double;
  using T = ADf<S,1>;
  const S dt = 1.0e-5;

  T p(0.5, {{1.0}});
  S lowerBound(0.2);
  T upperBound = T(2.7) + p;
  TimeIntegrator<T> integrator(lowerBound, upperBound, dt);

  // Result is initialized to zero
  EXPECT_DOUBLE_EQ(integrator.getResult().v(), S(0));
  EXPECT_NEAR(integrator.getResult().d(0), S(0), 1.0e-16);

  for (int iT(-1.0 / dt); iT <= upperBound / dt + 2; iT++) {
    integrator.takeValue(iT, util::sin(p*iT*dt + p*p));
  }
  EXPECT_NEAR(
    integrator.getResult().v(),
    - util::cos(p*upperBound + p*p) / p + util::cos(p*lowerBound + p*p) / p,
    1.0e-10);
  EXPECT_NEAR(
    integrator.getResult().d(0),
    util::sin(p*upperBound + p*p) * (upperBound + 2*p) / p
    - util::sin(p*lowerBound + p*p) * (lowerBound + 2*p) / p
    + util::cos(p*upperBound + p*p) / (p*p)
    - util::cos(p*lowerBound + p*p) / (p*p)
    + util::sin(p*upperBound + p*p),
    1.0e-10);
}

TEST(arrayL1test, testConstLinear)
{
  // integrate a constant and a linear function
  using T = double;
  const T dt = 1.0e-5;
  TimeIntegratorsArray<T,2,1> integrator(0.5, 2.0, dt);

  auto result0 = integrator.getResult();
  EXPECT_DOUBLE_EQ(result0[0], 0.0);
  EXPECT_DOUBLE_EQ(result0[1], 0.0);

  for (int iT = 0; iT < 2.0 / dt + 2; ++iT) {
    integrator.takeValues(iT, {1.0, T(2*iT*dt)});
  }

  auto result1 = integrator.getResult();
  EXPECT_DOUBLE_EQ(result1[0], 1.5);
  EXPECT_NEAR(result1[1], 3.75, 1.e-8);
}
