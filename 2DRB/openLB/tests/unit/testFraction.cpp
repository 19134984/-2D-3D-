/*  This file is part of the OpenLB library
*
*  Copyright (C) 2019 Adrian Kummerlaender
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
using namespace olb::utilities;

// Explicitly check for equality as `Fraction` should not introduce any additional
// errors compared to direct floating point usage. i.e. do not use `EXPECT_NEAR`

TEST(Fraction, testFractionAsDouble)
{
  EXPECT_EQ((Fraction{0  }.as<double>()), double{0} / double{1});
  EXPECT_EQ((Fraction{1  }.as<double>()), double{1} / double{1});
  EXPECT_EQ((Fraction{1,2}.as<double>()), double{1} / double{2});
  EXPECT_EQ((Fraction{2,1}.as<double>()), double{2} / double{1});
  EXPECT_EQ((Fraction{1,3}.as<double>()), double{1} / double{3});
  EXPECT_EQ((Fraction{3,1}.as<double>()), double{3} / double{1});
}

TEST(Fraction, testFractionAsFloat)
{
  EXPECT_EQ((Fraction{0  }.as<float>()), float{0} / float{1});
  EXPECT_EQ((Fraction{1  }.as<float>()), float{1} / float{1});
  EXPECT_EQ((Fraction{1,2}.as<float>()), float{1} / float{2});
  EXPECT_EQ((Fraction{2,1}.as<float>()), float{2} / float{1});
  EXPECT_EQ((Fraction{1,3}.as<float>()), float{1} / float{3});
  EXPECT_EQ((Fraction{3,1}.as<float>()), float{3} / float{1});
}

TEST(Fraction, testFractionInverseAsDouble)
{
  EXPECT_EQ((Fraction{1  }.inverseAs<double>()), double{1} / double{1});
  EXPECT_EQ((Fraction{1,2}.inverseAs<double>()), double{2} / double{1});
  EXPECT_EQ((Fraction{2,1}.inverseAs<double>()), double{1} / double{2});
  EXPECT_EQ((Fraction{1,3}.inverseAs<double>()), double{3} / double{1});
  EXPECT_EQ((Fraction{3,1}.inverseAs<double>()), double{1} / double{3});
}

TEST(Fraction, testFractionInverseAsFloat)
{
  EXPECT_EQ((Fraction{1  }.inverseAs<float>()), float{1} / float{1});
  EXPECT_EQ((Fraction{1,2}.inverseAs<float>()), float{2} / float{1});
  EXPECT_EQ((Fraction{2,1}.inverseAs<float>()), float{1} / float{2});
  EXPECT_EQ((Fraction{1,3}.inverseAs<float>()), float{3} / float{1});
  EXPECT_EQ((Fraction{3,1}.inverseAs<float>()), float{1} / float{3});
}
