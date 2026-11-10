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
using namespace olb::particles;
using namespace olb::particles::dynamics;
typedef double T;

//#define VERBOSE_TESTS


TEST(TestParticleFunctor, test2DdynamicFieldGroupContainer)
{
  const int D = 2;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>> DESCRIPTOR;
  std::size_t initialSize = 11;
  typedef DynamicFieldGroupsD<T,typename DESCRIPTOR::fields_t> FIELD_ARRAY_TYPE;

  Container<T,DESCRIPTOR,FIELD_ARRAY_TYPE> fieldContainer(initialSize);

  auto sizeA = fieldContainer.size();
  auto capacityA = fieldContainer.capacity();
  fieldContainer.push_back();
  auto sizeB = fieldContainer.size();
  auto capacityB = fieldContainer.capacity();
  fieldContainer.push_back();
  auto sizeC = fieldContainer.size();
  auto capacityC = fieldContainer.capacity();

#ifdef VERBOSE_TESTS
  std::cout << "SizeA=        " << sizeA << std::endl;
  std::cout << "CapacityA=    " << capacityA << std::endl;
  std::cout << "->push_back() " << std::endl;
  std::cout << "SizeB=        " << sizeB << std::endl;
  std::cout << "CapacityB=    " << capacityB << std::endl;
  std::cout << "->push_back() " << std::endl;
  std::cout << "SizeC=        " << sizeC << std::endl;
  std::cout << "CapacityC=    " << capacityC << std::endl;
#endif

  //Test general
  EXPECT_EQ(sizeA,initialSize);

  //Test setting specifc
  EXPECT_EQ(sizeC,13);
  EXPECT_EQ(capacityC,22);
}

TEST(TestParticleFunctor, test3DdynamicFieldGroupContainer)
{
  const int D = 3;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>> DESCRIPTOR;
  std::size_t initialSize = 20;
  typedef DynamicFieldGroupsD<T,typename DESCRIPTOR::fields_t> FIELD_ARRAY_TYPE;

  Container<T,DESCRIPTOR,FIELD_ARRAY_TYPE> fieldContainer(initialSize);

  auto sizeA = fieldContainer.size();
  auto capacityA = fieldContainer.capacity();
  fieldContainer.push_back();
  auto sizeB = fieldContainer.size();
  auto capacityB = fieldContainer.capacity();
  fieldContainer.push_back();
  auto sizeC = fieldContainer.size();
  auto capacityC = fieldContainer.capacity();

#ifdef VERBOSE_TESTS
  std::cout << "SizeA=        " << sizeA << std::endl;
  std::cout << "CapacityA=    " << capacityA << std::endl;
  std::cout << "->push_back() " << std::endl;
  std::cout << "SizeB=        " << sizeB << std::endl;
  std::cout << "CapacityB=    " << capacityB << std::endl;
  std::cout << "->push_back() " << std::endl;
  std::cout << "SizeC=        " << sizeC << std::endl;
  std::cout << "CapacityC=    " << capacityC << std::endl;
#endif

  //Test general
  EXPECT_EQ(sizeA,initialSize);

  //Test setting specifc
  EXPECT_EQ(sizeC,22);
  EXPECT_EQ(capacityC,40);
}


