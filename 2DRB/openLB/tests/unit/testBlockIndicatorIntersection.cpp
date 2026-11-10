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

#include "helper/framework.h"
#include "helper/unitTestHelper.h"
#include <gtest/gtest.h>

#include <olb.h>

using namespace olb;
using namespace olb::particles;
typedef double T;

//#define VERBOSE_TESTS

/// ensure olbInit is called for MPI
::testing::Environment* const foo_env =
    ::testing::AddGlobalTestEnvironment(new OlbInitTestEnvironment);

/// Create BlockGeometry environment
template <typename T, unsigned D>
BlockGeometry<T, D> createBlockGeometry(PhysR<T, D> origin, PhysR<T, D> extent,
                                        T delX, int padding = 2,
                                        int noCuboids = 1)
{
  int                   iCglob = 0;
  IndicatorCuboid<T, D> cuboidIndicator(extent, origin);
  CuboidDecomposition<T, D>  cuboidDecomposition(cuboidIndicator, delX, noCuboids);
  auto&                 cuboid = cuboidDecomposition.get(iCglob);
  BlockGeometry<T, D>   blockGeometry(cuboid, padding, iCglob);
  return blockGeometry;
}

/// Evaluate intersection
template <typename T, typename DESCRIPTOR>
void evaluateIntersection(BlockGeometry<T, DESCRIPTOR::d>& blockGeometry,
                          BlockLattice<T, DESCRIPTOR>&     blockLattice,
                          PhysR<T, DESCRIPTOR::d>&         particlePosition,
                          T                                particleCircumRadius,
                          LatticeR<DESCRIPTOR::d>&         start,
                          LatticeR<DESCRIPTOR::d>& end, bool& intersection,
                          int& numOfCells)
{
  constexpr unsigned D       = DESCRIPTOR::d;
  T                  padding = blockGeometry.getPadding();
  T                  delX    = blockGeometry.getDeltaR();

  //Check intersection
  intersection = getBlockParticleIntersection<T, D>(
      blockGeometry, 1. / T(delX), start, end, particlePosition,
      particleCircumRadius);

  //Count cells
  blockLattice.forSpatialLocations(start, end, [&](LatticeR<D> latticeR) {
    ++numOfCells;
  });

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);


#ifdef VERBOSE_TESTS

  //Write field to console
  ConsoleWriter<T,D,decltype(blockLattice)> consoleWriter( std::cout, blockLattice );
  consoleWriter.print();
  if constexpr(D==2){
    IndicatorCircle2D<T> indicator(particlePosition, particleCircumRadius);
    consoleWriter.writeIndicator(blockGeometry, indicator);
  } else {
    unsigned sliceAxis=2;
    int slicePosOnAxis=1;
    IndicatorSphere3D<T> indicator(particlePosition, particleCircumRadius);
    consoleWriter.writeIndicator(blockGeometry, indicator, sliceAxis, slicePosOnAxis);
  }

  //Summary
  std::cout << "BLOCKGEOMETRY:     " << std::endl;
  std::cout << " PhysR Origin=     " << blockGeometry.getOrigin() << std::endl;
  std::cout << " LatticeR Extent=  " << blockGeometry.getExtent() << std::endl;
  std::cout << " PhysR End=        "
            << blockGeometry.getOrigin() + (blockGeometry.getExtent() * delX)
            << std::endl;
  std::cout << "BLOCKLATTICE:      " << std::endl;
  std::cout << " LatticeR Extent=  " << blockLattice.getExtent() << std::endl;
  std::cout << "PARTICLE:          " << std::endl;
  std::cout << " PhysR Min=        "
            << particlePosition - PhysR<T, D>(particleCircumRadius)
            << std::endl;
  std::cout << " PhysR Max=        "
            << particlePosition + PhysR<T, D>(particleCircumRadius)
            << std::endl;
  std::cout << "RANGE:             " << std::endl;
  std::cout << " Intersection=     " << intersection << std::endl;
  std::cout << " LatticeR Start=   " << start << std::endl;
  std::cout << " LatticeR End=     " << end << std::endl;
  std::cout << " LatticeR Extent=  " << rangeExtentR << std::endl;
  std::cout << " PhysR Extent=     " << rangeExtent << std::endl;
  std::cout << " Number of cells=  " << numOfCells << std::endl;
#endif
}

/// 3D Versions

TEST(TestBlockParticleIntersection, 3DsmallRadiusPassLowerBounds)
{
  //Settings
  typedef descriptors::D3Q19<> DESCRIPTOR;
  constexpr unsigned           D = DESCRIPTOR::d;
  PhysR<T, D>                  origin(0.5);
  T                            baseSize = 1.;
  T                            delX     = 0.1;
  int                          padding  = 2;
  PhysR<T, D>                  particlePosition(origin);
  T                            particleCircumRadius = 0.2 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(0, 0, 0));
  EXPECT_EQ(end, LatticeR<D>(2, 2, 2));
  EXPECT_EQ(numOfCells, 27);
}

TEST(TestBlockParticleIntersection, 3DmediumRadiusPassLowerBounds)
{
  //Settings
  typedef descriptors::D3Q19<> DESCRIPTOR;
  constexpr unsigned           D = DESCRIPTOR::d;
  PhysR<T, D>                  origin(0.5);
  T                            baseSize = 1.;
  T                            delX     = 0.1;
  int                          padding  = 2;
  PhysR<T, D>                  particlePosition(origin);
  T                            particleCircumRadius = 0.5 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(0, 0, 0));
  EXPECT_EQ(end, LatticeR<D>(5, 5, 5));
  EXPECT_EQ(numOfCells, 216);
}

TEST(TestBlockParticleIntersection, 3DbigRadiusPassLowerBounds)
{
  //Settings
  typedef descriptors::D3Q19<> DESCRIPTOR;
  constexpr unsigned           D = DESCRIPTOR::d;
  PhysR<T, D>                  origin(0.5);
  T                            baseSize = 1.;
  T                            delX     = 0.1;
  int                          padding  = 2;
  PhysR<T, D>                  particlePosition(origin);
  T                            particleCircumRadius = 1.2 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(0, 0, 0));
  EXPECT_EQ(end, LatticeR<D>(10, 10, 10));
  EXPECT_EQ(numOfCells, 1331);
}

TEST(TestBlockParticleIntersection, 3DsmallRadiusPassUpperBounds)
{
  //Settings
  typedef descriptors::D3Q19<> DESCRIPTOR;
  constexpr unsigned           D = DESCRIPTOR::d;
  PhysR<T, D>                  origin(0.5);
  T                            baseSize = 1.;
  T                            delX     = 0.1;
  int                          padding  = 2;
  PhysR<T, D>                  particlePosition(origin + baseSize);
  T                            particleCircumRadius = 0.1 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(9, 9, 9));
  EXPECT_EQ(end, LatticeR<D>(10, 10, 10));
  EXPECT_EQ(numOfCells, 8);
}

TEST(TestBlockParticleIntersection, 3DmediumRadiusPassUpperBounds)
{
  //Settings
  typedef descriptors::D3Q19<> DESCRIPTOR;
  constexpr unsigned           D = DESCRIPTOR::d;
  PhysR<T, D>                  origin(0.5);
  T                            baseSize = 1.;
  T                            delX     = 0.1;
  int                          padding  = 2;
  PhysR<T, D>                  particlePosition(origin + baseSize);
  T                            particleCircumRadius = 0.65 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(3, 3, 3));
  EXPECT_EQ(end, LatticeR<D>(10, 10, 10));
  EXPECT_EQ(numOfCells, 512);
}

TEST(TestBlockParticleIntersection, 3DbigRadiusPassUpperBounds)
{
  //Settings
  typedef descriptors::D3Q19<> DESCRIPTOR;
  constexpr unsigned           D = DESCRIPTOR::d;
  PhysR<T, D>                  origin(0.5);
  T                            baseSize = 1.;
  T                            delX     = 0.1;
  int                          padding  = 2;
  PhysR<T, D>                  particlePosition(origin + baseSize);
  T                            particleCircumRadius = 0.9 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(0, 0, 0));
  EXPECT_EQ(end, LatticeR<D>(10, 10, 10));
  EXPECT_EQ(numOfCells, 1331);
}

TEST(TestBlockParticleIntersection, 3DmediumRadiusPassUpperBoundX)
{
  //Settings
  typedef descriptors::D3Q19<> DESCRIPTOR;
  constexpr unsigned           D = DESCRIPTOR::d;
  PhysR<T, D>                  origin(0.5);
  T                            baseSize = 1.;
  T                            delX     = 0.1;
  int                          padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(0.94, 0.5, 0.5) * baseSize);
  T           particleCircumRadius = 0.4 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(5, 0, 0));
  EXPECT_EQ(end, LatticeR<D>(10, 9, 9));
  EXPECT_EQ(numOfCells, 600);
}

TEST(TestBlockParticleIntersection, 3DmediumRadiusPassLowerBoundX)
{
  //Settings
  typedef descriptors::D3Q19<> DESCRIPTOR;
  constexpr unsigned           D = DESCRIPTOR::d;
  PhysR<T, D>                  origin(0.5);
  T                            baseSize = 1.;
  T                            delX     = 0.1;
  int                          padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(-0.1, 0.5, 0.5) * baseSize);
  T           particleCircumRadius = 0.45 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(0, 0, 0));
  EXPECT_EQ(end, LatticeR<D>(4, 10, 10));
  EXPECT_EQ(numOfCells, 605);
}

TEST(TestBlockParticleIntersection, 3DmediumRadiusPassUpperBoundY)
{
  //Settings
  typedef descriptors::D3Q19<> DESCRIPTOR;
  constexpr unsigned           D = DESCRIPTOR::d;
  PhysR<T, D>                  origin(0.5);
  T                            baseSize = 1.;
  T                            delX     = 0.1;
  int                          padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(0.5, 1.1, 0.5) * baseSize);
  T           particleCircumRadius = 0.43 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(0, 6, 0));
  EXPECT_EQ(end, LatticeR<D>(10, 10, 10));
  EXPECT_EQ(numOfCells, 605);
}

TEST(TestBlockParticleIntersection, 3DmediumRadiusPassLowerBoundY)
{
  //Settings
  typedef descriptors::D3Q19<> DESCRIPTOR;
  constexpr unsigned           D = DESCRIPTOR::d;
  PhysR<T, D>                  origin(0.5);
  T                            baseSize = 1.;
  T                            delX     = 0.1;
  int                          padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(0.5, 0.2, 0.5) * baseSize);
  T           particleCircumRadius = 0.35 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(1, 0, 1));
  EXPECT_EQ(end, LatticeR<D>(9, 6, 9));
  EXPECT_EQ(numOfCells, 567);
}

TEST(TestBlockParticleIntersection, 3DmediumRadiusPassUpperBoundZ)
{
  //Settings
  typedef descriptors::D3Q19<> DESCRIPTOR;
  constexpr unsigned           D = DESCRIPTOR::d;
  PhysR<T, D>                  origin(0.5);
  T                            baseSize = 1.;
  T                            delX     = 0.1;
  int                          padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(0.5, 0.5, 0.9) * baseSize);
  T           particleCircumRadius = 0.4 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(0, 0, 4));
  EXPECT_EQ(end, LatticeR<D>(9, 9, 10));
  EXPECT_EQ(numOfCells, 700);
}

TEST(TestBlockParticleIntersection, 3DmediumRadiusPassLowerBoundZ)
{
  //Settings
  typedef descriptors::D3Q19<> DESCRIPTOR;
  constexpr unsigned           D = DESCRIPTOR::d;
  PhysR<T, D>                  origin(0.5);
  T                            baseSize = 1.;
  T                            delX     = 0.1;
  int                          padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(0.5, 0.5, 0.1) * baseSize);
  T           particleCircumRadius = 0.4 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(0, 0, 0));
  EXPECT_EQ(end, LatticeR<D>(9, 9, 5));
  EXPECT_EQ(numOfCells, 600);
}

TEST(TestBlockParticleIntersection, 3DmediumRadiusOutOfUpperBoundX)
{
  //Settings
  typedef descriptors::D3Q19<> DESCRIPTOR;
  constexpr unsigned           D = DESCRIPTOR::d;
  PhysR<T, D>                  origin(0.5);
  T                            baseSize = 1.;
  T                            delX     = 0.1;
  int                          padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(1.8, 0.5, 0.5) * baseSize);
  T           particleCircumRadius = 0.45 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, false);
  EXPECT_EQ(start, LatticeR<D>(13, 0, 0));
  EXPECT_EQ(end, LatticeR<D>(10, 10, 10));
  EXPECT_EQ(numOfCells, 0);
}

TEST(TestBlockParticleIntersection, 3DmediumRadiusOutOfLowerBoundX)
{
  //Settings
  typedef descriptors::D3Q19<> DESCRIPTOR;
  constexpr unsigned           D = DESCRIPTOR::d;
  PhysR<T, D>                  origin(0.5);
  T                            baseSize = 1.;
  T                            delX     = 0.1;
  int                          padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(-0.8, 0.5, 0.5) * baseSize);
  T           particleCircumRadius = 0.45 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, false);
  EXPECT_EQ(start, LatticeR<D>(0, 0, 0));
  EXPECT_EQ(end, LatticeR<D>(-3, 10, 10));
  EXPECT_EQ(numOfCells, 0);
}

TEST(TestBlockParticleIntersection, 3DmediumRadiusOutOfUpperBoundY)
{
  //Settings
  typedef descriptors::D3Q19<> DESCRIPTOR;
  constexpr unsigned           D = DESCRIPTOR::d;
  PhysR<T, D>                  origin(0.5);
  T                            baseSize = 1.;
  T                            delX     = 0.1;
  int                          padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(0.5, 1.9, 0.5) * baseSize);
  T           particleCircumRadius = 0.43 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, false);
  EXPECT_EQ(start, LatticeR<D>(0, 14, 0));
  EXPECT_EQ(end, LatticeR<D>(10, 10, 10));
  EXPECT_EQ(numOfCells, 0);
}

TEST(TestBlockParticleIntersection, 3DmediumRadiusOutOfLowerBoundY)
{
  //Settings
  typedef descriptors::D3Q19<> DESCRIPTOR;
  constexpr unsigned           D = DESCRIPTOR::d;
  PhysR<T, D>                  origin(0.5);
  T                            baseSize = 1.;
  T                            delX     = 0.1;
  int                          padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(0.5, -1.0, 0.5) * baseSize);
  T           particleCircumRadius = 0.5 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, false);
  EXPECT_EQ(start, LatticeR<D>(0, 0, 0));
  EXPECT_EQ(end, LatticeR<D>(10, -5, 10));
  EXPECT_EQ(numOfCells, 0);
}

TEST(TestBlockParticleIntersection, 3DmediumRadiusOutOfUpperBoundZ)
{
  //Settings
  typedef descriptors::D3Q19<> DESCRIPTOR;
  constexpr unsigned           D = DESCRIPTOR::d;
  PhysR<T, D>                  origin(0.5);
  T                            baseSize = 1.;
  T                            delX     = 0.1;
  int                          padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(0.5, 0.5, 1.7) * baseSize);
  T           particleCircumRadius = 0.4 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, false);
  EXPECT_EQ(start, LatticeR<D>(0, 0, 13));
  EXPECT_EQ(end, LatticeR<D>(9, 9, 10));
  EXPECT_EQ(numOfCells, 0);
}

TEST(TestBlockParticleIntersection, 3DmediumRadiusOutOfLowerBoundZ)
{
  //Settings
  typedef descriptors::D3Q19<> DESCRIPTOR;
  constexpr unsigned           D = DESCRIPTOR::d;
  PhysR<T, D>                  origin(0.5);
  T                            baseSize = 1.;
  T                            delX     = 0.1;
  int                          padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(0.5, 0.5, -1.6) * baseSize);
  T           particleCircumRadius = 0.4 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, false);
  EXPECT_EQ(start, LatticeR<D>(0, 0, 0));
  EXPECT_EQ(end, LatticeR<D>(9, 9, -12));
  EXPECT_EQ(numOfCells, 0);
}

/// 2D Versions

TEST(TestBlockParticleIntersection, 2DsmallRadiusPassLowerBounds)
{
  //Settings
  typedef descriptors::D2Q9<> DESCRIPTOR;
  constexpr unsigned          D = DESCRIPTOR::d;
  PhysR<T, D>                 origin(0.5);
  T                           baseSize = 1.;
  T                           delX     = 0.1;
  int                         padding  = 2;
  PhysR<T, D>                 particlePosition(origin);
  T                           particleCircumRadius = 0.2 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(0, 0));
  EXPECT_EQ(end, LatticeR<D>(2, 2));
  EXPECT_EQ(numOfCells, 9);
}

TEST(TestBlockParticleIntersection, 2DmediumRadiusPassLowerBounds)
{
  //Settings
  typedef descriptors::D2Q9<> DESCRIPTOR;
  constexpr unsigned          D = DESCRIPTOR::d;
  PhysR<T, D>                 origin(0.5);
  T                           baseSize = 1.;
  T                           delX     = 0.1;
  int                         padding  = 2;
  PhysR<T, D>                 particlePosition(origin);
  T                           particleCircumRadius = 0.5 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(0, 0));
  EXPECT_EQ(end, LatticeR<D>(5, 5));
  EXPECT_EQ(numOfCells, 36);
}

TEST(TestBlockParticleIntersection, 2DbigRadiusPassLowerBounds)
{
  //Settings
  typedef descriptors::D2Q9<> DESCRIPTOR;
  constexpr unsigned          D = DESCRIPTOR::d;
  PhysR<T, D>                 origin(0.5);
  T                           baseSize = 1.;
  T                           delX     = 0.1;
  int                         padding  = 2;
  PhysR<T, D>                 particlePosition(origin);
  T                           particleCircumRadius = 1.2 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(0, 0));
  EXPECT_EQ(end, LatticeR<D>(10, 10));
  EXPECT_EQ(numOfCells, 121);
}

TEST(TestBlockParticleIntersection, 2DsmallRadiusPassUpperBounds)
{
  //Settings
  typedef descriptors::D2Q9<> DESCRIPTOR;
  constexpr unsigned          D = DESCRIPTOR::d;
  PhysR<T, D>                 origin(0.5);
  T                           baseSize = 1.;
  T                           delX     = 0.1;
  int                         padding  = 2;
  PhysR<T, D>                 particlePosition(origin + baseSize);
  T                           particleCircumRadius = 0.1 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(9, 9));
  EXPECT_EQ(end, LatticeR<D>(10, 10));
  EXPECT_EQ(numOfCells, 4);
}

TEST(TestBlockParticleIntersection, 2DmediumRadiusPassUpperBounds)
{
  //Settings
  typedef descriptors::D2Q9<> DESCRIPTOR;
  constexpr unsigned          D = DESCRIPTOR::d;
  PhysR<T, D>                 origin(0.5);
  T                           baseSize = 1.;
  T                           delX     = 0.1;
  int                         padding  = 2;
  PhysR<T, D>                 particlePosition(origin + baseSize);
  T                           particleCircumRadius = 0.65 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(3, 3));
  EXPECT_EQ(end, LatticeR<D>(10, 10));
  EXPECT_EQ(numOfCells, 64);
}

TEST(TestBlockParticleIntersection, 2DbigRadiusPassUpperBounds)
{
  //Settings
  typedef descriptors::D2Q9<> DESCRIPTOR;
  constexpr unsigned          D = DESCRIPTOR::d;
  PhysR<T, D>                 origin(0.5);
  T                           baseSize = 1.;
  T                           delX     = 0.1;
  int                         padding  = 2;
  PhysR<T, D>                 particlePosition(origin + baseSize);
  T                           particleCircumRadius = 0.9 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(0, 0));
  EXPECT_EQ(end, LatticeR<D>(10, 10));
  EXPECT_EQ(numOfCells, 121);
}

TEST(TestBlockParticleIntersection, 2DmediumRadiusPassUpperBoundX)
{
  //Settings
  typedef descriptors::D2Q9<> DESCRIPTOR;
  constexpr unsigned          D = DESCRIPTOR::d;
  PhysR<T, D>                 origin(0.5);
  T                           baseSize = 1.;
  T                           delX     = 0.1;
  int                         padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(0.94, 0.5) * baseSize);
  T           particleCircumRadius = 0.4 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(5, 0));
  EXPECT_EQ(end, LatticeR<D>(10, 9));
  EXPECT_EQ(numOfCells, 60);
}

TEST(TestBlockParticleIntersection, 2DmediumRadiusPassLowerBoundX)
{
  //Settings
  typedef descriptors::D2Q9<> DESCRIPTOR;
  constexpr unsigned          D = DESCRIPTOR::d;
  PhysR<T, D>                 origin(0.5);
  T                           baseSize = 1.;
  T                           delX     = 0.1;
  int                         padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(-0.1, 0.5) * baseSize);
  T           particleCircumRadius = 0.45 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(0, 0));
  EXPECT_EQ(end, LatticeR<D>(4, 10));
  EXPECT_EQ(numOfCells, 55);
}

TEST(TestBlockParticleIntersection, 2DmediumRadiusPassUpperBoundY)
{
  //Settings
  typedef descriptors::D2Q9<> DESCRIPTOR;
  constexpr unsigned          D = DESCRIPTOR::d;
  PhysR<T, D>                 origin(0.5);
  T                           baseSize = 1.;
  T                           delX     = 0.1;
  int                         padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(0.5, 1.1) * baseSize);
  T           particleCircumRadius = 0.43 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(0, 6));
  EXPECT_EQ(end, LatticeR<D>(10, 10));
  EXPECT_EQ(numOfCells, 55);
}

TEST(TestBlockParticleIntersection, 2DmediumRadiusPassLowerBoundY)
{
  //Settings
  typedef descriptors::D2Q9<> DESCRIPTOR;
  constexpr unsigned          D = DESCRIPTOR::d;
  PhysR<T, D>                 origin(0.5);
  T                           baseSize = 1.;
  T                           delX     = 0.1;
  int                         padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(0.5, 0.2) * baseSize);
  T           particleCircumRadius = 0.35 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, true);
  EXPECT_EQ(start, LatticeR<D>(1, 0));
  EXPECT_EQ(end, LatticeR<D>(9, 6));
  EXPECT_EQ(numOfCells, 63);
}

TEST(TestBlockParticleIntersection, 2DmediumRadiusOutOfUpperBoundX)
{
  //Settings
  typedef descriptors::D2Q9<> DESCRIPTOR;
  constexpr unsigned          D = 2;
  PhysR<T, D>                 origin(0.5);
  T                           baseSize = 1.;
  T                           delX     = 0.1;
  int                         padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(1.8, 0.5) * baseSize);
  T           particleCircumRadius = 0.45 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, false);
  EXPECT_EQ(start, LatticeR<D>(13, 0));
  EXPECT_EQ(end, LatticeR<D>(10, 10));
  EXPECT_EQ(numOfCells, 0);
}

TEST(TestBlockParticleIntersection, 2DmediumRadiusOutOfLowerBoundX)
{
  //Settings
  typedef descriptors::D2Q9<> DESCRIPTOR;
  constexpr unsigned          D = 2;
  PhysR<T, D>                 origin(0.5);
  T                           baseSize = 1.;
  T                           delX     = 0.1;
  int                         padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(-.8, 0.5) * baseSize);
  T           particleCircumRadius = 0.45 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, false);
  EXPECT_EQ(start, LatticeR<D>(0, 0));
  EXPECT_EQ(end, LatticeR<D>(-3, 10));
  EXPECT_EQ(numOfCells, 0);
}

TEST(TestBlockParticleIntersection, 2DmediumRadiusOutOfUpperBoundY)
{
  //Settings
  typedef descriptors::D2Q9<> DESCRIPTOR;
  constexpr unsigned          D = 2;
  PhysR<T, D>                 origin(0.5);
  T                           baseSize = 1.;
  T                           delX     = 0.1;
  int                         padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(0.5, 1.9) * baseSize);
  T           particleCircumRadius = 0.43 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, 0);
  EXPECT_EQ(start, LatticeR<D>(0, 14));
  EXPECT_EQ(end, LatticeR<D>(10, 10));
  EXPECT_EQ(numOfCells, 0);
}

TEST(TestBlockParticleIntersection, 2DmediumRadiusOutOfLowerBoundY)
{
  //Settings
  typedef descriptors::D2Q9<> DESCRIPTOR;
  constexpr unsigned          D = DESCRIPTOR::d;
  PhysR<T, D>                 origin(0.5);
  T                           baseSize = 1.;
  T                           delX     = 0.1;
  int                         padding  = 2;
  PhysR<T, D> particlePosition(origin + PhysR<T, D>(0.5, -1.0) * baseSize);
  T           particleCircumRadius = 0.5 * baseSize;

  //Create BlockGeometry and BlockLattice
  PhysR<T, D>         extent(baseSize);
  BlockGeometry<T, D> blockGeometry =
      createBlockGeometry<T, D>(origin, extent, delX, padding);
  ConcreteBlockLattice<T, DESCRIPTOR> blockLattice(blockGeometry.getExtent(),
                                                   padding);

  //Evaluate intersection
  LatticeR<D> start, end;
  bool        intersection;
  int         numOfCells = 0;
  evaluateIntersection<T, DESCRIPTOR>(blockGeometry, blockLattice,
                                      particlePosition, particleCircumRadius,
                                      start, end, intersection, numOfCells);

  //Test general
  LatticeR<D> rangeExtentR(end - start);
  PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
  PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
  EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
  EXPECT_GE(start, LatticeR<D>(-padding));

  //Test setting specifc
  EXPECT_EQ(intersection, false);
  EXPECT_EQ(start, LatticeR<D>(0, 0));
  EXPECT_EQ(end, LatticeR<D>(10, -5));
  EXPECT_EQ(numOfCells, 0);
}

TEST(TestBlockParticleIntersection, multipleBlocks2D)
{
  //Settings
  typedef descriptors::D2Q9<> DESCRIPTOR;
  constexpr unsigned          D         = DESCRIPTOR::d;
  constexpr unsigned          noCuboids = 4;
  T                           baseSize  = 2.;
  T                           delX      = 0.1;
  int                         padding   = 0;
  PhysR<T, D> particlePosition(PhysR<T, D>(0.5, 0.5) * baseSize);
  T           particleCircumRadius = 0.25 * baseSize;

  //Create BlockGeometry and BlockLattice
  int            nCell = baseSize / delX;
  CuboidDecomposition2D<T> cuboidDecomposition(0., delX, nCell, noCuboids);
  HeuristicLoadBalancer<T>    loadBalancer(cuboidDecomposition);
  SuperGeometry<T, D>         superGeometry(cuboidDecomposition, loadBalancer);
  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter2D(
    int {nCell}, // resolution: number of voxels per charPhysL
    (T)1, // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)baseSize, // charPhysLength: reference length of simulation geometry
    (T)1, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)1, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)1  // physDensity: physical density in __kg / m^3__
);
  SuperLattice<T, DESCRIPTOR> sLattice(converter2D, superGeometry);

  std::array<Vector<T, D>, noCuboids> startPerBlock;
  std::array<Vector<T, D>, noCuboids> endPerBlock;
  unsigned                            numCells = 0;

  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    //Evaluate intersection
    LatticeR<D> start, end;
    bool        intersection;
    int         numOfCells = 0;
    evaluateIntersection<T, DESCRIPTOR>(superGeometry.getBlockGeometry(iC),
                                        sLattice.getBlock(iC), particlePosition,
                                        particleCircumRadius, start, end,
                                        intersection, numOfCells);

    //Test general
    LatticeR<D> rangeExtentR(end - start);
    PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
    PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
    EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
    EXPECT_GE(start, LatticeR<D>(-padding));

    //Test setting specifc
    EXPECT_EQ(intersection, true);

    numCells += numOfCells;
    startPerBlock[iC] =
        superGeometry.getBlockGeometry(iC).getOrigin() / delX + start;
    endPerBlock[iC] =
        superGeometry.getBlockGeometry(iC).getOrigin() / delX + end;
  }

  Vector<T, D> startComb = startPerBlock[0];
  Vector<T, D> endComb   = endPerBlock[0];

  EXPECT_EQ(numCells, std::pow(2*std::round(particleCircumRadius / delX) + 1, D));
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    for (int jC = iC + 1; jC < iC; ++jC) {
      for (unsigned iD = 0; iD < D; ++iD) {
        EXPECT_NE(startPerBlock[iC][iD], startPerBlock[jC][iD]);
        EXPECT_NE(endPerBlock[iC][iD], endPerBlock[jC][iD]);
      }
    }
    for (unsigned iD = 0; iD < D; ++iD) {
      startComb[iD] = std::min(startComb[iD], startPerBlock[iC][iD]);
      endComb[iD]   = std::max(endComb[iD], endPerBlock[iC][iD]);
    }
  }

  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_EQ(startComb[iD],
              (particlePosition - PhysR<T, D>(particleCircumRadius))[iD] /
                  delX);
    EXPECT_EQ(endComb[iD],
              (particlePosition + PhysR<T, D>(particleCircumRadius))[iD] /
                  delX);
  }
}

TEST(TestBlockParticleIntersection, multipleBlocks3D)
{
  //Settings
  typedef descriptors::D3Q19<> DESCRIPTOR;
  constexpr unsigned           D         = DESCRIPTOR::d;
  constexpr unsigned           noCuboids = 8;
  T                            baseSize  = 2.;
  T                            delX      = 0.1;
  int                          padding   = 0;
  PhysR<T, D> particlePosition(PhysR<T, D>(0.5, 0.5, 0.5) * baseSize);
  T           particleCircumRadius = 0.25 * baseSize;

  //Create BlockGeometry and BlockLattice
  int                 nCell = baseSize / delX;
  CuboidDecomposition3D<T>      cuboidDecomposition(0., delX, nCell,
                                          noCuboids);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, D>      superGeometry(cuboidDecomposition, loadBalancer);
  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter3D(
    int {nCell}, // resolution: number of voxels per charPhysL
    (T)1, // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)baseSize, // charPhysLength: reference length of simulation geometry
    (T)1, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)1, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)1  // physDensity: physical density in __kg / m^3__
);
  SuperLattice<T, DESCRIPTOR> sLattice(converter3D, superGeometry);

  std::array<Vector<T, D>, noCuboids> startPerBlock;
  std::array<Vector<T, D>, noCuboids> endPerBlock;
  unsigned                            numCells = 0;

  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    //Evaluate intersection
    LatticeR<D> start, end;
    bool        intersection;
    int         numOfCells = 0;
    evaluateIntersection<T, DESCRIPTOR>(superGeometry.getBlockGeometry(iC),
                                        sLattice.getBlock(iC), particlePosition,
                                        particleCircumRadius, start, end,
                                        intersection, numOfCells);

    //Test general
    LatticeR<D> rangeExtentR(end - start);
    PhysR<T, D> rangeExtent(PhysR<T, D>(rangeExtentR) * delX);
    PhysR<T, D> rangeExtent_maxExpected(2. * particleCircumRadius + 2. * delX);
    EXPECT_LE(rangeExtent, rangeExtent_maxExpected);
    EXPECT_GE(start, LatticeR<D>(-padding));

    //Test setting specifc
    EXPECT_EQ(intersection, true);

    numCells += numOfCells;
    startPerBlock[iC] =
        superGeometry.getBlockGeometry(iC).getOrigin() / delX + start;
    endPerBlock[iC] =
        superGeometry.getBlockGeometry(iC).getOrigin() / delX + end;
  }

  Vector<T, D> startComb = startPerBlock[0];
  Vector<T, D> endComb   = endPerBlock[0];

  EXPECT_EQ(numCells, std::pow(2*std::round(particleCircumRadius / delX) + 1, D));
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    for (int jC = iC + 1; jC < iC; ++jC) {
      for (unsigned iD = 0; iD < D; ++iD) {
        EXPECT_NE(startPerBlock[iC][iD], startPerBlock[jC][iD]);
        EXPECT_NE(endPerBlock[iC][iD], endPerBlock[jC][iD]);
      }
    }
    for (unsigned iD = 0; iD < D; ++iD) {
      startComb[iD] = std::min(startComb[iD], startPerBlock[iC][iD]);
      endComb[iD]   = std::max(endComb[iD], endPerBlock[iC][iD]);
    }
  }

  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_EQ(startComb[iD],
              (particlePosition - PhysR<T, D>(particleCircumRadius))[iD] /
                  delX);
    EXPECT_EQ(endComb[iD],
              (particlePosition + PhysR<T, D>(particleCircumRadius))[iD] /
                  delX);
  }
}
