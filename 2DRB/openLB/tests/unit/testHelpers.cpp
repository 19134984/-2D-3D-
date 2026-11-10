/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Julius Jessberger
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

TEST(TestIsContained, testDouble)
{
  using T = double;
  std::vector<T> c{1,2,3};
  EXPECT_TRUE(util::isContained(c, T{1}));
  EXPECT_FALSE(util::isContained(c, T{4}));
}

TEST(TestIsContained, testInt)
{
  using T = int;
  std::vector<T> c{1,2,3};
  EXPECT_TRUE(util::isContained(c, T{1}));
  EXPECT_FALSE(util::isContained(c, T{4}));
}

TEST(TestContainerCreator, testStdVector)
{
  using T = double;
  constexpr unsigned n = 2;
  using C = std::vector<T>;

  const auto c = util::ContainerCreator<C>::create(n);

  EXPECT_EQ(c.size(), 2);
  EXPECT_EQ(c[0], T());
  EXPECT_EQ(c[1], T());
}

TEST(TestContainerCreator, testStdArray)
{
  using T = double;
  constexpr unsigned n = 3;
  using C = std::array<T,n>;

  const auto c = util::ContainerCreator<C>::create(n);

  EXPECT_EQ(c.size(), 3);
  EXPECT_EQ(c[0], T());
  EXPECT_EQ(c[1], T());
  EXPECT_EQ(c[2], T());
}

TEST(TestContainerCreator, testOlbVector)
{
  using T = double;
  constexpr unsigned n = 2;
  using C = Vector<T,n>;

  const auto c = util::ContainerCreator<C>::create(n);

  EXPECT_EQ(c.size(), 2);
  EXPECT_EQ(c[0], T());
  EXPECT_EQ(c[1], T());
}
