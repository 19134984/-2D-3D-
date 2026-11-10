/*  This file is part of the OpenLB library
*
*  Copyright (C) 2021 Nicolas Hafen
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

#include <type_traits>

#include <gtest/gtest.h>
#include "helper/unitTestHelper.h"

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
typedef double T;

//#define VERBOSE_TESTS



/// CORRECT

TEST(TestFieldDescriptors, testDescriptor)
{
  const int D = 3;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>> DESCRIPTOR;

  using type = typename meta::derived_type_in_nested<DESCRIPTOR>::type;

#ifdef VERBOSE_TESTS
  cout << "derivedType: " << typeid(type).name() << std::endl;
#endif

  //Test setting specifc
  bool isExpectedType = std::is_same<type,DESCRIPTOR>::value;
  EXPECT_EQ(isExpectedType,true);
}


TEST(TestFieldDescriptors, testDescriptorGroup)
{
  const int D = 3;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>> DESCRIPTOR;

  using type = typename meta::derived_type_in_nested<DESCRIPTOR,GENERAL>::type;

#ifdef VERBOSE_TESTS
  cout << "derivedType: " << typeid(type).name() << std::endl;
#endif

  //Test setting specifc
  bool isExpectedType = std::is_same<type,GENERAL_TMP<D>>::value;
  EXPECT_EQ(isExpectedType,true);
}


TEST(TestFieldDescriptors, testDescriptorGroupField)
{
  const int D = 3;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>> DESCRIPTOR;

  using type = typename meta::derived_type_in_nested<DESCRIPTOR,GENERAL,POSITION>::type;

#ifdef VERBOSE_TESTS
  cout << "derivedType: " << typeid(type).name() << std::endl;
#endif

  //Test setting specifc
  bool isExpectedType = std::is_same<type,POSITION>::value;
  EXPECT_EQ(isExpectedType,true);
}


///INCORRECT


TEST(TestFieldDescriptors, testDescriptorWrong)
{
  using type = typename meta::derived_type_in_nested<D2Q9<>>::type;

#ifdef VERBOSE_TESTS
  cout << "derivedType: " << typeid(type).name() << std::endl;
#endif

  //Test setting specifc
  bool isExpectedType = std::is_same<type,D2Q9<>>::value;
  EXPECT_EQ(isExpectedType,true);
}


TEST(TestFieldDescriptors, testDescriptorGroupWrong)
{
  const int D = 3;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>> DESCRIPTOR;

  using type = typename meta::derived_type_in_nested<DESCRIPTOR,MOBILITY>::type;

#ifdef VERBOSE_TESTS
  cout << "derivedType: " << typeid(type).name() << std::endl;
#endif

  //Test setting specifc
  bool isExpectedType = std::is_same<type,void>::value;
  EXPECT_EQ(isExpectedType,true);
}


TEST(TestFieldDescriptors, testDescriptorGroupFieldWrong)
{
  const int D = 3;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>> DESCRIPTOR;

  using type = typename meta::derived_type_in_nested<DESCRIPTOR,GENERAL,ADHESION>::type;

#ifdef VERBOSE_TESTS
  cout << "derivedType: " << typeid(type).name() << std::endl;
#endif

  //Test setting specifc
  bool isExpectedType = std::is_same<type,void>::value;
  EXPECT_EQ(isExpectedType,true);
}

TEST(TestFieldDescriptors, testDescriptorGroupWrongField)
{
  const int D = 3;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>> DESCRIPTOR;

  using type = typename meta::derived_type_in_nested<DESCRIPTOR,MOBILITY,POSITION>::type;

#ifdef VERBOSE_TESTS
  cout << "derivedType: " << typeid(type).name() << std::endl;
#endif

  //Test setting specifc
  bool isExpectedType = std::is_same<type,void>::value;
  EXPECT_EQ(isExpectedType,true);
}



///CHECK PROVIDED



template<typename DESCRIPTOR, typename GROUP, typename FIELD>
bool evalProvided(){
  bool provides = DESCRIPTOR::template providesNested<GROUP,FIELD>();
#ifdef VERBOSE_TESTS
  cout << "Evaluate provided fields:" << std::endl;
  cout << " DESCRIPTOR: " << typeid(DESCRIPTOR).name() << std::endl;
  cout << " GROUP:      " << typeid(GROUP).name() << std::endl;
  cout << " FIELD:      " << typeid(FIELD).name() << std::endl;
  cout << " ok:         " << provides << std::endl;
#endif
  return provides;
}



TEST(TestFieldDescriptors, testPosition)
{
  const int D = 3;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>> DESCRIPTOR;
  bool provides = evalProvided<DESCRIPTOR,GENERAL,POSITION>();
  //Test setting specifc
  EXPECT_EQ(provides,true);
}


TEST(TestFieldDescriptors, testAdhesion)
{
  const int D = 3;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>> DESCRIPTOR;
  bool provides = evalProvided<DESCRIPTOR,GENERAL,ADHESION>();
  //Test setting specifc
  EXPECT_EQ(provides,false);
}

TEST(TestFieldDescriptors, testForcing)
{
  const int D = 3;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>> DESCRIPTOR;
  bool provides = evalProvided<DESCRIPTOR,MOBILITY,ADHESION>();
  //Test setting specifc
  EXPECT_EQ(provides,false);
}


//ONE SPECIAL TEST TO CHECK ORIGINAL FUNCTIONALITY
TEST(TestFieldDescriptors, testVelocityNumerator)
{
  typedef D3Q19<VELOCITY_NUMERATOR> DESCRIPTOR;
  bool provides = DESCRIPTOR::template provides<VELOCITY_NUMERATOR>();
  //Test setting specifc
  EXPECT_EQ(provides,true);
}
