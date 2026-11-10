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
#include "helper/framework.h"
#include "helper/unitTestHelper.h"
#include <gtest/gtest.h>
#include <iostream>
#include <type_traits>

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::particles;
using namespace olb::particles::dynamics;
//using namespace olb::particles::access;

using SubgridParticle3DsingleCuboid  = PARTICLE_DESCRIPTOR<3, GENERAL_TMP<3>,
      MOBILITY_VERLET_NO_ANGLE<3>, FORCING_SUBGRID<3>,
      PHYSPROPERTIES_SUBGRID<3>,
      DYNBEHAVIOUR_BASIC>;

typedef double                         T;
typedef ResolvedParticle3D             PARTICLETYPE;
typedef ResolvedDecomposedParticle3D     PARTICLETYPEPARA;
typedef SubgridParticle3DsingleCuboid  PARTICLETYPESUB;

::testing::Environment* const env =
    ::testing::AddGlobalTestEnvironment(new OlbInitTestEnvironment);

//#define VERBOSE_TESTS

TEST(TestSuperParticleSystem, cuboidNeighborhood)
{
  //Settings
  typedef ResolvedDecomposedParticle3D PARTICLETYPE;
  constexpr unsigned                 D = PARTICLETYPE::d;
  PhysR<T, D>                        origin(0.0);
  const T                            baseSize = 1.0;
  const T                            delX     = 0.1;
  //Set up environment
  Vector<T, D>             extend(baseSize);
  IndicatorCuboid<T, D>    cuboid(extend, origin);
  const int                noOfCuboids = 27;
  CuboidDecomposition<T, D>     cuboidDecomposition(cuboid, delX, noOfCuboids);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, D>      sGeometry(cuboidDecomposition, loadBalancer, 2);

  //Create super particle system
  SuperParticleSystem<T, PARTICLETYPE> sParticleSystem(sGeometry);

  // Test with default circumRadius
  {
    std::vector<std::set<int> > expectedCuboidNeighborhood
    = { { 13, 12, 10, 9, 4, 3, 1 },
        { 14, 12, 11, 10, 9, 5, 4, 3, 2, 13, 0 },
        { 13, 11, 10, 5, 4, 14, 1 },
        { 16, 15, 12, 10, 9, 7, 6, 4, 1, 13, 0 },
        { 17, 16, 15, 14, 0, 13, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12 },
        { 16, 13, 11, 10, 8, 7, 17, 4, 2, 14, 1 },
        { 15, 13, 12, 7, 4, 16, 3 },
        { 15, 14, 13, 12, 8, 6, 5, 17, 4, 16, 3 },
        { 16, 14, 13, 7, 5, 17, 4 },
        { 22, 21, 19, 18, 12, 10, 4, 3, 1, 13, 0 },
        { 23, 22, 21, 20, 0, 13, 1, 14, 2, 3, 4, 5, 18, 9, 11, 12, 19 },
        { 22, 20, 19, 13, 23, 10, 5, 4, 2, 14, 1 },
        { 25, 24, 22, 21, 0, 13, 1, 3, 16, 4, 6, 19, 7, 9, 10, 15, 18 },
        { 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
        { 26, 25, 23, 22, 1, 2, 4, 17, 5, 7, 20, 8, 10, 11, 13, 16, 19 },
        { 24, 22, 21, 13, 25, 12, 7, 6, 4, 16, 3 },
        { 26, 25, 24, 23, 3, 4, 17, 5, 6, 7, 8, 21, 12, 13, 14, 15, 22 },
        { 25, 23, 22, 16, 14, 26, 13, 8, 7, 5, 4 },
        { 21, 19, 13, 12, 10, 22, 9 },
        { 21, 20, 18, 14, 13, 12, 11, 23, 10, 22, 9 },
        { 22, 19, 14, 13, 11, 23, 10 },
        { 24, 19, 18, 16, 15, 13, 25, 12, 10, 22, 9 },
        { 26, 25, 24, 23, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21 },
        { 25, 22, 20, 19, 17, 16, 14, 26, 13, 11, 10 },
        { 22, 21, 16, 15, 13, 25, 12 },
        { 24, 23, 22, 21, 17, 16, 15, 14, 26, 13, 12 },
        { 25, 23, 22, 17, 16, 14, 13 } };

    std::vector<std::set<int>> cuboidNeighborhood =
        sParticleSystem.getCuboidNeighborhood();

    EXPECT_EQ(noOfCuboids, cuboidNeighborhood.size());
    ASSERT_EQ(expectedCuboidNeighborhood.size(), cuboidNeighborhood.size());

    for (unsigned i = 0; i < cuboidNeighborhood.size(); ++i) {
      EXPECT_EQ(expectedCuboidNeighborhood[i], cuboidNeighborhood[i]);
      EXPECT_EQ(particles::communication::getSurfaceTouchingICs(
            sParticleSystem, Vector<T, PARTICLETYPE::d>() , T() ,
            Vector<bool, PARTICLETYPE::d>(false) , i ),
          cuboidNeighborhood[i]);
    }
  }

  // Test with a bigger circumRadius
  sParticleSystem.updateOffsetFromCircumRadius(baseSize * 0.21);
  {
    std::vector<std::set<int>> expectedCuboidNeighborhood = {
      { 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 },
      { 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 0, 13, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
      { 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 0, 13, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
      { 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 0, 13, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
      { 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 0, 13, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12 },
      { 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 0, 13, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12 },
      { 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 0, 13, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12 },
      { 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 0, 13, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12 },
      { 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 0, 13, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12 },
      { 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 0, 13, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12 },
      { 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 0, 13, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12 },
      { 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 0, 13, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12 },
      { 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 0, 13, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 },
      { 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
      { 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
      { 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 14, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
      { 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 15, 14, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
      { 26, 25, 24, 23, 22, 21, 20, 19, 18, 16, 15, 14, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
      { 26, 25, 24, 23, 22, 21, 20, 19, 17, 16, 15, 14, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
      { 26, 25, 24, 23, 22, 21, 20, 18, 17, 16, 15, 14, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
      { 26, 25, 24, 23, 22, 21, 19, 18, 17, 16, 15, 14, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
      { 26, 25, 24, 23, 22, 20, 19, 18, 17, 16, 15, 14, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
      { 26, 25, 24, 23, 21, 20, 19, 18, 17, 16, 15, 14, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
      { 26, 25, 24, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
      { 26, 25, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
      { 26, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 },
      { 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 }
    };

    std::vector<std::set<int>> cuboidNeighborhood =
        sParticleSystem.getCuboidNeighborhood();

    EXPECT_EQ(noOfCuboids, cuboidNeighborhood.size());
    ASSERT_EQ(expectedCuboidNeighborhood.size(), cuboidNeighborhood.size());

    for (unsigned i = 0; i < cuboidNeighborhood.size(); ++i) {
      EXPECT_EQ(expectedCuboidNeighborhood[i], cuboidNeighborhood[i]);
      EXPECT_EQ(particles::communication::getSurfaceTouchingICs(
            sParticleSystem, Vector<T, PARTICLETYPE::d>() , T() ,
            Vector<bool, PARTICLETYPE::d>(false) , i ),
          cuboidNeighborhood[i]);
    }
  }

  // Test with a gigantic circumRadius
  sParticleSystem.updateOffsetFromCircumRadius(baseSize);
  {
    std::set<int> allCuboids;
    for(int i=0; i<noOfCuboids; ++i) {
      allCuboids.insert(i);
    }

    std::vector<std::set<int>> expectedCuboidNeighborhood;
    expectedCuboidNeighborhood.resize(noOfCuboids);
    for(int i=0; i<expectedCuboidNeighborhood.size(); ++i) {
      expectedCuboidNeighborhood[i] = allCuboids;
      expectedCuboidNeighborhood[i].erase(expectedCuboidNeighborhood[i].find(i));
    }

    std::vector<std::set<int>> cuboidNeighborhood =
        sParticleSystem.getCuboidNeighborhood();

    EXPECT_EQ(noOfCuboids, cuboidNeighborhood.size());
    ASSERT_EQ(expectedCuboidNeighborhood.size(), cuboidNeighborhood.size());

    for (unsigned i = 0; i < cuboidNeighborhood.size(); ++i) {
      EXPECT_EQ(expectedCuboidNeighborhood[i], cuboidNeighborhood[i]);
      EXPECT_EQ(particles::communication::getSurfaceTouchingICs(
            sParticleSystem, Vector<T, PARTICLETYPE::d>() , T() ,
            Vector<bool, PARTICLETYPE::d>(false) , i ),
          cuboidNeighborhood[i]);
    }
  }

}
