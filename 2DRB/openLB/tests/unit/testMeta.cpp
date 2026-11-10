/*  This file is part of the OpenLB library
*
*  Copyright (C) 2021 Adrian Kummerlaender
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

#include "core/meta.h"

using namespace olb;

template <unsigned I>
using item = std::integral_constant<unsigned, I>;

TEST(meta_list, get)
{
  using list = meta::list<item<0>,item<1>,item<2>,item<3>,item<4>>;
  EXPECT_EQ(list::get<0>::value, 0);
  EXPECT_EQ(list::get<1>::value, 1);
  EXPECT_EQ(list::get<2>::value, 2);
  EXPECT_EQ(list::get<3>::value, 3);
  EXPECT_EQ(list::get<4>::value, 4);
}

TEST(meta_list, for_each_single)
{
  unsigned i=0;
  meta::list<item<0>>::for_each([&i](auto id) {
    EXPECT_EQ(id.get(), i);
    i += 1;
  });
  EXPECT_EQ(i, 1);
}

TEST(meta_list, for_each_multiple)
{
  unsigned i=0;
  meta::list<item<0>,item<1>,item<2>,item<3>,item<4>>::for_each([&i](auto id) {
    EXPECT_EQ(id.get(), i);
    i += 1;
  });
  EXPECT_EQ(i, 5);
}

TEST(meta_list, push)
{
  EXPECT_EQ(meta::list<item<0>>::get<0>::value, 0);
  EXPECT_EQ(meta::list<item<0>>::push<item<1>>::get<0>::value, 1);
  EXPECT_EQ(meta::list<item<0>>::push<item<1>>::get<1>::value, 0);
}

TEST(meta_list, append)
{
  EXPECT_EQ(meta::list<item<0>>::get<0>::value, 0);
  EXPECT_EQ(meta::list<item<0>>::append<item<1>>::get<0>::value, 0);
  EXPECT_EQ(meta::list<item<0>>::append<item<1>>::get<1>::value, 1);
}

TEST(meta_list, include_disjoint)
{
  using odd  = meta::list<item<1>,item<3>>;
  using even = meta::list<item<0>,item<2>,item<4>>;

  using all = even::decompose_into<odd::include>;

  EXPECT_EQ(all::get<0>(), 1);
  EXPECT_EQ(all::get<1>(), 3);
  EXPECT_EQ(all::get<2>(), 0);
  EXPECT_EQ(all::get<3>(), 2);
  EXPECT_EQ(all::get<4>(), 4);
}

TEST(meta_list, merge_disjoint)
{
  using odd  = meta::list<item<1>,item<3>>;
  using even = meta::list<item<0>,item<2>,item<4>>;

  using all = meta::merge<even,odd>;

  EXPECT_EQ(all::get<0>(), 1);
  EXPECT_EQ(all::get<1>(), 3);
  EXPECT_EQ(all::get<2>(), 0);
  EXPECT_EQ(all::get<3>(), 2);
  EXPECT_EQ(all::get<4>(), 4);
}

TEST(meta_list, include_non_disjoint)
{
  using a = meta::list<item<1>,item<2>,item<3>,item<5>>;
  using b = meta::list<item<1>,item<3>,item<4>>;

  using all = a::decompose_into<b::include>;

  EXPECT_EQ(all::get<0>(), 1);
  EXPECT_EQ(all::get<1>(), 3);
  EXPECT_EQ(all::get<2>(), 4);
  EXPECT_EQ(all::get<3>(), 2);
  EXPECT_EQ(all::get<4>(), 5);
}

TEST(meta_list, merge_non_disjoint)
{
  using a = meta::list<item<1>,item<2>,item<3>,item<5>>;
  using b = meta::list<item<1>,item<3>,item<4>>;

  using all = meta::merge<a,b>;

  EXPECT_EQ(all::get<0>(), 1);
  EXPECT_EQ(all::get<1>(), 3);
  EXPECT_EQ(all::get<2>(), 4);
  EXPECT_EQ(all::get<3>(), 2);
  EXPECT_EQ(all::get<4>(), 5);
}
