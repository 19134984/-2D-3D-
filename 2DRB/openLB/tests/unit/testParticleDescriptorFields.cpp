/*  This file is part of the OpenLB library
*
*  Copyright (C) 2023 Nicolas Hafen
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
#include "helper/framework.h"

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::particles;
//using namespace olb::particles::dynamics;
typedef double T;

::testing::Environment* const env = ::testing::AddGlobalTestEnvironment(new OlbInitTestEnvironment);

//#define VERBOSE_TESTS




TEST(TestParticleDescriptorFields, testActive)
{
  //Settings
  const int D = 2;
  typedef PARTICLE_DESCRIPTOR<
    D,
    GENERAL_TMP<D>,
    DYNBEHAVIOUR_BASIC
  > PARTICLETYPE;
  //Create particle system
  ParticleSystem<T,PARTICLETYPE> particleSystem(4);
  //Retrieve particle
  auto particle = particleSystem.get(2);
  //Check activity
  bool active = access::isActive( particle );
  //Evaluate
  EXPECT_EQ(active,true);
}

