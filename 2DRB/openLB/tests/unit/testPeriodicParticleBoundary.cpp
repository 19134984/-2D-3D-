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

#include "helper/unitTestHelper.h"
#include <gtest/gtest.h>
#include <iostream>

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::particles;

typedef double T;

TEST(PeriodicResolvedParticleBoundary, checkIfParticleIsOutsideDomainOutideX)
{
  constexpr unsigned D = 3;

  // Input
  const PhysR<T, D>     cellMin(-5.25);
  const PhysR<T, D>     cellMax(5.25);
  const PhysR<T, D>     position(4, 2, 2);
  const T               circumRadius(2);
  const Vector<bool, D> periodic(true, false, false);

  // Resulting data
  bool        surfaceOutOfGeometry = false;
  PhysR<T, D> ghostPos(0.);

  checkSmoothIndicatorOutOfGeometry(surfaceOutOfGeometry, ghostPos, cellMin,
                                    cellMax, position, circumRadius, periodic);

  EXPECT_TRUE(surfaceOutOfGeometry);
  EXPECT_EQ(ghostPos[0], -6 - 0.5);
  EXPECT_EQ(ghostPos[1], 2);
  EXPECT_EQ(ghostPos[2], 2);
}

TEST(PeriodicResolvedParticleBoundary, checkIfParticleIsOutsideDomainOutideY)
{
  constexpr unsigned D = 3;

  // Input
  const PhysR<T, D>     cellMin(-5.25);
  const PhysR<T, D>     cellMax(5.25);
  const PhysR<T, D>     position(2, 4, 2);
  const T               circumRadius(2);
  const Vector<bool, D> periodic(false, true, false);

  // Resulting data
  bool        surfaceOutOfGeometry = false;
  PhysR<T, D> ghostPos(0.);

  checkSmoothIndicatorOutOfGeometry(surfaceOutOfGeometry, ghostPos, cellMin,
                                    cellMax, position, circumRadius, periodic);

  EXPECT_TRUE(surfaceOutOfGeometry);
  EXPECT_EQ(ghostPos[0], 2);
  EXPECT_EQ(ghostPos[1], -6 - 0.5);
  EXPECT_EQ(ghostPos[2], 2);
}

TEST(PeriodicResolvedParticleBoundary, checkIfParticleIsOutsideDomainOutideZ)
{
  constexpr unsigned D = 3;

  // Input
  const PhysR<T, D>     cellMin(-5.25);
  const PhysR<T, D>     cellMax(5.25);
  const PhysR<T, D>     position(2, 2, 4);
  const T               circumRadius(2);
  const Vector<bool, D> periodic(false, false, true);

  // Resulting data
  bool        surfaceOutOfGeometry = false;
  PhysR<T, D> ghostPos(0.);

  checkSmoothIndicatorOutOfGeometry(surfaceOutOfGeometry, ghostPos, cellMin,
                                    cellMax, position, circumRadius, periodic);

  EXPECT_TRUE(surfaceOutOfGeometry);
  EXPECT_EQ(ghostPos[0], 2);
  EXPECT_EQ(ghostPos[1], 2);
  EXPECT_EQ(ghostPos[2], -6 - 0.5);
}

TEST(PeriodicResolvedParticleBoundary, checkIfParticleIsOutsideDomainOutide)
{
  constexpr unsigned D = 3;

  // Input
  const PhysR<T, D>     cellMin(-5.25);
  const PhysR<T, D>     cellMax(5.25);
  const PhysR<T, D>     position(-4, -4, -4);
  const T               circumRadius(2);
  const Vector<bool, D> periodic(true);

  // Resulting data
  bool        surfaceOutOfGeometry = false;
  PhysR<T, D> ghostPos(0.);

  checkSmoothIndicatorOutOfGeometry(surfaceOutOfGeometry, ghostPos, cellMin,
                                    cellMax, position, circumRadius, periodic);

  EXPECT_TRUE(surfaceOutOfGeometry);
  EXPECT_EQ(ghostPos[0], 6 + 0.5);
  EXPECT_EQ(ghostPos[1], 6 + 0.5);
  EXPECT_EQ(ghostPos[2], 6 + 0.5);
}

TEST(PeriodicResolvedParticleBoundary,
     checkIfParticleIsOutsideDomainOutideNonPeriodic)
{
  constexpr unsigned D = 3;

  // Input
  const PhysR<T, D>     cellMin(-5.25);
  const PhysR<T, D>     cellMax(5.25);
  const PhysR<T, D>     position(-4, -4, -4);
  const T               circumRadius(2);
  const Vector<bool, D> periodic(false);

  // Resulting data
  bool        surfaceOutOfGeometry = true;
  PhysR<T, D> ghostPos(0.);

  checkSmoothIndicatorOutOfGeometry(surfaceOutOfGeometry, ghostPos, cellMin,
                                    cellMax, position, circumRadius, periodic);

  EXPECT_FALSE(surfaceOutOfGeometry);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(ghostPos[iDim], position[iDim]);
  }
}

TEST(PeriodicResolvedParticleBoundary, checkIfParticleIsOutsideDomainInside)
{
  constexpr unsigned D = 3;

  // Input
  const PhysR<T, D>     cellMin(-5.25);
  const PhysR<T, D>     cellMax(5.25);
  const PhysR<T, D>     position(-3, -3, -3);
  const T               circumRadius(2);
  const Vector<bool, D> periodic(true);

  // Resulting data
  bool        surfaceOutOfGeometry = true;
  PhysR<T, D> ghostPos(0.);

  checkSmoothIndicatorOutOfGeometry(surfaceOutOfGeometry, ghostPos, cellMin,
                                    cellMax, position, circumRadius, periodic);

  EXPECT_FALSE(surfaceOutOfGeometry);
  for (unsigned iDim = 0; iDim < D; ++iDim) {
    EXPECT_EQ(ghostPos[iDim], position[iDim]);
  }
}

