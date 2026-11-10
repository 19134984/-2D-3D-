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
#include <olb.h>
#include <gtest/gtest.h>
#include <iostream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::particles;
using namespace olb::particles::contact;

typedef double                        T;
typedef ResolvedParticleWithContact2D PARTICLE_2D;
typedef ResolvedParticleWithContact3D PARTICLE_3D;

using D2D =
    D2Q9<POROSITY, VELOCITY_NUMERATOR, VELOCITY_DENOMINATOR, CONTACT_DETECTION>;
using D3D = D3Q19<POROSITY, VELOCITY_NUMERATOR, VELOCITY_DENOMINATOR,
                  CONTACT_DETECTION>;

/// ensure olbInit is called for MPI
::testing::Environment* const foo_env =
    ::testing::AddGlobalTestEnvironment(new OlbInitTestEnvironment);

UnitConverterFromResolutionAndRelaxationTime<T, D2D> const converter2D(
    int {1}, // resolution: number of voxels per charPhysL
    (T)1, // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)1, // charPhysLength: reference length of simulation geometry
    (T)1, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)1, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)1  // physDensity: physical density in __kg / m^3__
);

UnitConverterFromResolutionAndRelaxationTime<T, D3D> const converter3D(
    int {1}, // resolution: number of voxels per charPhysL
    (T)1, // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)1, // charPhysLength: reference length of simulation geometry
    (T)1, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)1, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)1  // physDensity: physical density in __kg / m^3__
);

//#define VERBOSE_TESTS

TEST(onLatticeContactDetection, resetOfHelperField)
{
  typedef ParticleContactArbitraryFromOverlapVolume<T, D3D::d, true>
      PARTICLECONTACTTYPE;
  typedef WallContactArbitraryFromOverlapVolume<T, D3D::d, true>
      WALLCONTACTTYPE;

  const T radius  = 5;
  const int nCellX  = 6. / converter3D.getPhysDeltaX();
  const T eps     = 0.;
  const T density = 1;
  // Creation of lattice
  CuboidDecomposition<T, D3D::d> cuboidDecomposition(
      0., converter3D.getPhysDeltaX(), {nCellX, nCellX, nCellX}, 4);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, D3D::d> superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
#ifdef VERBOSE_TESTS
  cuboidDecomposition.printExtended();
  superGeometry.print();
  converter3D.print();
#endif
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);
  sLattice.defineDynamics<PorousParticleBGKdynamics<T, D3D>>(superGeometry, 1);
  sLattice.setParameter<descriptors::OMEGA>(
      converter3D.getLatticeRelaxationFrequency());
  // Set initial condition
  AnalyticalConst3D<T, T> one(1.);
  sLattice.defineField<descriptors::POROSITY>(
      superGeometry.getMaterialIndicator(1), one);
  AnalyticalConst3D<T, T>    zero(0.);
  AnalyticalConst3D<T, T>    rho(1.);
  AnalyticalComposed3D<T, T> u(zero, zero, zero);
  //Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU(superGeometry, 1, rho, u);
  sLattice.iniEquilibrium(superGeometry, 1, rho, u);
  // Make the lattice ready for simulation
  sLattice.initialize();

  ParticleSystem<T, PARTICLE_3D>                            particleSystem;
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;

  //Create sphere 1
  Vector<T, D3D::d> pos1(T {0});
  particles::creators::addResolvedSphere3D(particleSystem, pos1, radius, eps,
                                           density);
  auto particleA       = particleSystem.get(0);
  auto solidBoundaries = std::vector<SolidBoundary<T, D3D::d>>();

  particles::dynamics::coupleResolvedParticlesToLattice<
      T, D3D, PARTICLE_3D, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
      particleSystem, contactContainer, superGeometry, sLattice, converter3D,
      solidBoundaries, []() {
        return Vector<bool, D3D::d>(false, false, false);
      });

  // Expect that entries are written to the CONTACT_DETECTION field
  unsigned numFieldEntries = 0;
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    BlockLattice<T, D3D>&           blockLattice = sLattice.getBlock(iC);
    const BlockGeometry<T, D3D::d>& blockGeometry =
        superGeometry.getBlockGeometry(iC);
    const int padding = blockGeometry.getPadding();
    forSpatialLocationsInBlockParticleIntersection(
        blockGeometry, blockLattice, padding,
        particles::access::getPosition(particleA),
        particles::access::getRadius(particleA),
        [&](const LatticeR<D3D::d>& latticeR) {
          auto contactHelper =
              blockLattice.get(latticeR)
                  .template getFieldPointer<descriptors::CONTACT_DETECTION>();
          if (contactHelper[0] > 0) {
            ++numFieldEntries;
          }
        });
  }
  EXPECT_GT(numFieldEntries, 0);

  // Now expect that the field resets during collide
  numFieldEntries = 0;
  sLattice.collideAndStream();
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    BlockLattice<T, D3D>&           blockLattice = sLattice.getBlock(iC);
    const BlockGeometry<T, D3D::d>& blockGeometry =
        superGeometry.getBlockGeometry(iC);
    const int padding = blockGeometry.getPadding();
    forSpatialLocationsInBlockParticleIntersection(
        blockGeometry, blockLattice, padding,
        particles::access::getPosition(particleA),
        particles::access::getRadius(particleA),
        [&](const LatticeR<D3D::d>& latticeR) {
          auto contactHelper =
              blockLattice.get(latticeR)
                  .template getFieldPointer<descriptors::CONTACT_DETECTION>();
          if (contactHelper[0] > 0) {
            ++numFieldEntries;
          }
        });
  }
  EXPECT_EQ(numFieldEntries, 0);
}
