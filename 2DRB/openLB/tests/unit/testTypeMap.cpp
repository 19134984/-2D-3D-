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

#include "utilities/typeMap.h"

using namespace olb::meta;

struct key1 { };
struct key2 { };
struct key3 { };

struct value1 { };
struct value2 { };
struct value3 { };

TEST(TypeMap, keys)
{
  EXPECT_TRUE((std::is_same_v<map<key1, value1>::keys_t,                             list<key1>>));
  EXPECT_TRUE((std::is_same_v<map<key1, value1, key2, value2>::keys_t,               list<key1,key2>>));
  EXPECT_TRUE((std::is_same_v<map<key1, value1, key2, value2, key3, value3>::keys_t, list<key1,key2,key3>>));
}

TEST(TypeMap, values)
{
  EXPECT_TRUE((std::is_same_v<map<key1, value1>::values_t,                             list<value1>>));
  EXPECT_TRUE((std::is_same_v<map<key1, value1, key2, value2>::values_t,               list<value1,value2>>));
  EXPECT_TRUE((std::is_same_v<map<key1, value1, key2, value2, key3, value3>::values_t, list<value1,value2,value3>>));
}

TEST(TypeMap, size)
{
  EXPECT_EQ((map<key1, value1>::size),                             1);
  EXPECT_EQ((map<key1, value1, key2, value2>::size),               2);
  EXPECT_EQ((map<key1, value1, key2, value2, key3, value3>::size), 3);
}

TEST(TypeMap, index)
{
  EXPECT_EQ((map<key1, value1, key2, value2, key3, value3>::keys_t::index<key1>()), 0);
  EXPECT_EQ((map<key1, value1, key2, value2, key3, value3>::keys_t::index<key2>()), 1);
  EXPECT_EQ((map<key1, value1, key2, value2, key3, value3>::keys_t::index<key3>()), 2);
}

TEST(TypeMap, invalid_index)
{
  EXPECT_EQ((map<key1, value1, key2, value2, key3, value3>::keys_t::index<value1>()), 3);
}

TEST(TypeMap, element)
{
  EXPECT_TRUE((std::is_same_v<map<key1, value1, key2, value2, key3, value3>::keys_t::get<0>, key1>));
  EXPECT_TRUE((std::is_same_v<map<key1, value1, key2, value2, key3, value3>::keys_t::get<1>, key2>));
  EXPECT_TRUE((std::is_same_v<map<key1, value1, key2, value2, key3, value3>::keys_t::get<2>, key3>));
}
