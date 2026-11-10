/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2021 Jan E. Marquardt
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

// see https://gitlab.com/openlb/olb/-/issues/114
// In [1] on p. 41, the analytical solution for a rotating sphere in a linear shear flow is given. We should use it as a first (simple) test.
// [1] H. Nirschl, “Partikelbewegungen in Scherströmungen ohne und mit Berücksichtigung des Einflusses angrenzender Wände,” VDI Verlag, no. 7, 1997.

#include <olb.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "../unit/helper/framework.h"

#include <gtest/gtest.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace olb::particles;
using namespace olb::particles::dynamics;

typedef double T;
//Define lattice type
typedef PorousParticleD3Q19Descriptor DESCRIPTOR;
//Define particleType
typedef ResolvedSphere3D PARTICLETYPE;

/// ensure olbInit is called for MPI
::testing::Environment* const foo_env =
    ::testing::AddGlobalTestEnvironment(new OlbInitTestEnvironment);

// Geometry
const Vector<T, 3>    channelExtend(0.3, 0.1, 0.25);
const Vector<bool, 3> periodicity(true, false, true);

struct TestParam {
  T uMax;
  T xPosition;
};

TestParam params[] = {
    {0.005, 0.5 * channelExtend[0]},
    {0.005, 0.0 * channelExtend[0]},
    { 0.01, 0.5 * channelExtend[0]},
    { 0.01, 0.0 * channelExtend[0]}
};

class ParticleMigrationTest : public ::testing::TestWithParam<TestParam> {};

// Parameters for the simulation setup
const int N        = 10;
const T   tau      = 1.0;
const T   maxPhysT = 0.4; // max. simulation time in s, SI unit
const T   iTwrite  = 0.2 * maxPhysT;

// Fluid
const T rhoF      = 1200;
const T viscosity = 1e4 / (rhoF * 1e3);

// Particle
const T            rhoP                 = rhoF;
const T            sphereRadius         = 0.015;
const Vector<T, 3> externalAcceleration = {0., 0., 0.};

void prepareGeometry(UnitConverter<T, DESCRIPTOR> const& converter,
                     SuperGeometry<T, 3>& superGeometry, int N)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0, 2);
  superGeometry.rename(2, 1, {1, 1, 1});

  const T deltaX = converter.getPhysDeltaX();

  Vector<T, 3> boundaryOrigin(-0.5 * deltaX);
  Vector<T, 3> boundaryExtend(channelExtend[0] + deltaX, deltaX,
                              channelExtend[2] + deltaX);

  IndicatorCuboid3D inflow(boundaryExtend, boundaryOrigin);
  superGeometry.rename(2, 3, 1, inflow);

  boundaryOrigin[1] = channelExtend[1] - 0.5 * deltaX;
  IndicatorCuboid3D outflow(boundaryExtend, boundaryOrigin);
  superGeometry.rename(2, 4, 1, outflow);

  superGeometry.clean();
  superGeometry.innerClean();

  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
  return;
}

void prepareLattice(SuperLattice<T, DESCRIPTOR>&        sLattice,
                    UnitConverter<T, DESCRIPTOR> const& converter,
                    SuperGeometry<T, 3>& superGeometry, const T uMax)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  sLattice.defineDynamics<NoDynamics>(superGeometry, 0);
  sLattice.defineDynamics<PorousParticleBGKdynamics>(superGeometry, 1);
  sLattice.defineDynamics<PorousParticleBGKdynamics>(superGeometry, 3);
  sLattice.defineDynamics<PorousParticleBGKdynamics>(superGeometry, 4);

  {
    auto& communicator = sLattice.getCommunicator(stage::PostPostProcess());
    communicator
        .requestFields<POROSITY, VELOCITY_NUMERATOR, VELOCITY_DENOMINATOR>();
    communicator.requestOverlap(sLattice.getOverlap());
    communicator.exchangeRequests();
  }

  boundary::set<boundary::LocalVelocity>(sLattice, superGeometry, 3);
  boundary::set<boundary::LocalVelocity>(sLattice, superGeometry, 4);

  sLattice.setParameter<descriptors::OMEGA>(
      converter.getLatticeRelaxationFrequency());

  AnalyticalConst3D<T, T> one(1.);
  AnalyticalConst3D<T, T> zero(0.);
  AnalyticalConst3D<T, T> maxU(uMax / converter.getConversionFactorVelocity());
  AnalyticalConst3D<T, T> minU(-uMax / converter.getConversionFactorVelocity());
  sLattice.defineField<POROSITY>(
      superGeometry.getMaterialIndicator({0, 1, 2, 3, 4}), one);

  AnalyticalConst3D<T, T> rho(1.);

  AnalyticalComposed3D<T, T> u(zero, zero, zero);
  AnalyticalComposed3D<T, T> uUpper(maxU, zero, zero);
  AnalyticalComposed3D<T, T> uLower(minU, zero, zero);

  //Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU(superGeometry, 1, rho, u);
  sLattice.iniEquilibrium(superGeometry, 1, rho, u);
  sLattice.defineRhoU(superGeometry, 3, rho, uUpper);
  sLattice.iniEquilibrium(superGeometry, 3, rho, uUpper);
  sLattice.defineRhoU(superGeometry, 4, rho, uLower);
  sLattice.iniEquilibrium(superGeometry, 4, rho, uLower);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

TEST_P(ParticleMigrationTest, testCouetteFlow3d)
{

  // === 1st Step: Initialization ===
  TestParam param                   = GetParam();
  const T   uMax                    = param.uMax;
  const T   expectedAngularVelocity = 0.5 * 2 * uMax / channelExtend[1];

  const Vector<T, 3> spherePos(param.xPosition, 0.5 * channelExtend[1],
                               0.5 * channelExtend[2]);

  OstreamManager clout(std::cout, "main");

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> converter(
      (int)N,              //resolution
      (T)tau,              //charLatticeVelocity
      (T)2 * sphereRadius, //charPhysLength
      (T)uMax,             //charPhysVelocity
      (T)viscosity,        //physViscosity
      (T)rhoF              //physDensity
  );
  converter.print();

  /// === 2rd Step: Prepare Geometry ===
  Vector<T, 3>         channelOrigin(0., 0., 0.);
  IndicatorCuboid3D<T> channel(channelExtend, channelOrigin);

#ifdef PARALLEL_MODE_MPI
  CuboidDecomposition3D<T> cuboidDecomposition(channel,
                                     converter.getConversionFactorLength(),
                                     singleton::mpi().getSize());
#else
  CuboidDecomposition3D<T> cuboidDecomposition(channel,
                                     converter.getConversionFactorLength(), 1);
#endif

  cuboidDecomposition.setPeriodicity(periodicity);

  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, 3>      superGeometry(cuboidDecomposition, loadBalancer, 2);
  prepareGeometry(converter, superGeometry, N);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice(converter, superGeometry);

  prepareLattice(sLattice, converter, superGeometry, uMax);

  /// === 4th Step: Main Loop with Timer ===
  Timer<double> timer(converter.getLatticeTime(maxPhysT),
                      superGeometry.getStatistics().getNvoxel());
  timer.start();
  // Create ParticleSystem
  ParticleSystem<T, PARTICLETYPE> particleSystem;

  //Create particle manager handling coupling, gravity and particle dynamics
  ParticleManager<T, DESCRIPTOR, PARTICLETYPE> particleManager(
      particleSystem, superGeometry, sLattice, converter, externalAcceleration,
      periodicity);

  // Calculate particle quantities
  const T epsilon = 0.5 * converter.getConversionFactorLength();

  // Create Particle 1
  creators::addResolvedSphere3D(particleSystem, spherePos, sphereRadius,
                                epsilon, rhoP);

  // Create and assign resolved particle dynamics
  particleSystem.defineDynamics<
    VerletParticleDynamics<T, PARTICLETYPE>>();

  clout << "MaxIT: " << converter.getLatticeTime(maxPhysT) << std::endl;
  for (size_t iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {
    // Execute particle manager
    particleManager.execute<
        couple_lattice_to_particles_single_cuboid<T, DESCRIPTOR, PARTICLETYPE>,
        apply_gravity<T, PARTICLETYPE>, process_dynamics_single_cuboid<T, PARTICLETYPE>,
        couple_particles_to_lattice_single_cuboid<T, DESCRIPTOR, PARTICLETYPE>>();

    /// Writes output on the console
    if (iT % converter.getLatticeTime(iTwrite) == 0) {
      timer.update(iT);
      timer.printStep();
      sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
      for (std::size_t iP = 0; iP < particleSystem.size(); ++iP) {
        auto particle = particleSystem.get(iP);
        particles::io::printResolvedParticleInfo(particle);
      }
    }

    // Collide and stream
    sLattice.collideAndStream();
  }

  auto               particle = particleSystem.get(0);
  const Vector<T, 3> angVel =
      particle.template getField<MOBILITY, ANG_VELOCITY>();
  // we allow a relative error of 3% or less
  EXPECT_NEAR(expectedAngularVelocity, angVel[2],
              0.03 * expectedAngularVelocity);

  timer.stop();
}

INSTANTIATE_TEST_CASE_P(ParticleMigration, ParticleMigrationTest,
                        testing::ValuesIn(params));
