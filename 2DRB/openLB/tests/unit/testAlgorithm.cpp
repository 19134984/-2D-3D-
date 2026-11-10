/*  This file is part of the OpenLB library
*
*  Copyright (C) 2020 Julius Jessberger, Jan E. Marquardt
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
using T = int;

TEST(oalgorithm, testMax)
{
  T z1 = 2;
  T z2 = 3;
  T z3 = -z1;
  T z4 = -z2;

  EXPECT_EQ(util::max(z1, z2), std::max(z1, z2));
  EXPECT_EQ(util::max(z1, z4), std::max(z1, z4));
  EXPECT_EQ(util::max(z3, z4), std::max(z3, z4));

  EXPECT_EQ(util::max(std::initializer_list<T>{1, -4, 3}),
            std::max(std::initializer_list<T>{1, -4, 3}));
}


TEST(oalgorithm, testMin)
{
  T z1 = 2;
  T z2 = 3;
  T z3 = -z1;
  T z4 = -z2;

  EXPECT_EQ(util::min(z1, z2), std::min(z1, z2));
  EXPECT_EQ(util::min(z1, z4), std::min(z1, z4));
  EXPECT_EQ(util::min(z3, z4), std::min(z3, z4));

  EXPECT_EQ(util::min(std::initializer_list<T> {1, -4, 3}),
            std::min(std::initializer_list<T> {1, -4, 3}));
}

TEST(populationsContributing, toDirection1)
{
  using DESCRIPTOR = descriptors::D2Q9<>;
  auto pops = util::populationsContributingToDirection<DESCRIPTOR,1,0>();
  EXPECT_EQ(pops.size(), 3);
  EXPECT_EQ(pops[0], 5);
  EXPECT_EQ(pops[1], 6);
  EXPECT_EQ(pops[2], 7);
}

TEST(populationsContributing, toDirection2)
{
  using DESCRIPTOR = descriptors::D2Q9<>;
  auto pops = util::populationsContributingToDirection<DESCRIPTOR,-1,0>();
  EXPECT_EQ(pops.size(), 3);
  EXPECT_EQ(pops[0], 1);
  EXPECT_EQ(pops[1], 2);
  EXPECT_EQ(pops[2], 3);
}

TEST(populationsContributing, toDirection3)
{
  using DESCRIPTOR = descriptors::D2Q9<>;
  auto pops = util::populationsContributingToDirection<DESCRIPTOR,0,1>();
  EXPECT_EQ(pops.size(), 3);
  EXPECT_EQ(pops[0], 1);
  EXPECT_EQ(pops[1], 7);
  EXPECT_EQ(pops[2], 8);
}
