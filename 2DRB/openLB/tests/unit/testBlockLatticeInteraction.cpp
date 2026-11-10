#include "helper/framework.h"
#include "helper/unitTestHelper.h"
#include <olb.h>
#include <gtest/gtest.h>
#include <iostream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::particles;

using D2D = D2Q9<POROSITY, VELOCITY_NUMERATOR, VELOCITY_DENOMINATOR>;
using D3D = D3Q19<POROSITY, VELOCITY_NUMERATOR, VELOCITY_DENOMINATOR>;

typedef double T;

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

template <typename PARTICLETYPE>
T getPorosity(Particle<T, PARTICLETYPE>& particle, T input[])
{
  T tmp[1];
  particles::resolved::evalSolidVolumeFraction(tmp, input, particle);
  return tmp[0];
}

// 2D Tests
TEST(testBlockLatticeInteraction, testParticleFarOut2D)
{
  constexpr unsigned D            = 2;
  int                nCell[D]     = {4, 4};
  T const            radius       = 2;
  T const            eps          = 1;
  unsigned const     numberBlocks = 2;

  // Creation of lattice
  CuboidDecomposition2D<T>      cuboidDecomposition(0., converter2D.getPhysDeltaX(),
                                          nCell, numberBlocks);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, D>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D2D> sLattice(converter2D,superGeometry);

  ParticleSystem<T, ResolvedParticle2D> particleSystem;
  Vector<T, D>                          pos(-500);
  particles::creators::addResolvedCircle2D(particleSystem, pos, radius, eps);
  auto particle = particleSystem.get(0);

  //Write particle field
  setSuperParticleField(superGeometry, sLattice, converter2D, particle,
                        Vector<bool, D>(false));

  for (unsigned iC = 0; iC < numberBlocks; ++iC) {
    sLattice.getBlock(iC).forCoreSpatialLocations([&](LatticeR<D> latticeR) {
      T porosity =
          sLattice.getBlock(iC).get(latticeR).template getField<POROSITY>();
      Vector<T, D> velocityDenominator =
          sLattice.getBlock(iC)
              .get(latticeR)
              .template getField<VELOCITY_DENOMINATOR>();
      Vector<T, D> velocityNumerator =
          sLattice.getBlock(iC)
              .get(latticeR)
              .template getField<VELOCITY_NUMERATOR>();

      EXPECT_EQ(1., porosity);
      for (unsigned iDim = 0; iDim < D; ++iDim) {
        EXPECT_EQ(0., velocityDenominator[iDim]);
        EXPECT_EQ(0., velocityNumerator[iDim]);
      }
    });
  }
}

TEST(testBlockLatticeInteraction, testParticleInsideWithZeroVelocity2D)
{
  constexpr unsigned D            = 2;
  int                nCell[D]     = {4, 4};
  T const            radius       = 2;
  T const            eps          = 1;
  unsigned const     numberBlocks = 2;

  // Creation of lattice
  CuboidDecomposition2D<T>      cuboidDecomposition(0., converter2D.getPhysDeltaX(),
                                          nCell, numberBlocks);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, D>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D2D> sLattice(converter2D, superGeometry);

  ParticleSystem<T, ResolvedParticle2D> particleSystem;
  Vector<T, D>                          pos(0);
  particles::creators::addResolvedCircle2D(particleSystem, pos, radius, eps);
  auto particle = particleSystem.get(0);

  //Write particle field
  setSuperParticleField(superGeometry, sLattice, converter2D, particle,
                        Vector<bool, D>(false));

  for (unsigned iC = 0; iC < numberBlocks; ++iC) {
    sLattice.getBlock(iC).forCoreSpatialLocations([&](LatticeR<D> latticeR) {
      PhysR<T, D> physR = superGeometry.getBlockGeometry(iC).getPhysR(latticeR);
      T           porosity =
          sLattice.getBlock(iC).get(latticeR).template getField<POROSITY>();
      Vector<T, D> velocityDenominator =
          sLattice.getBlock(iC)
              .get(latticeR)
              .template getField<VELOCITY_DENOMINATOR>();
      Vector<T, D> velocityNumerator =
          sLattice.getBlock(iC)
              .get(latticeR)
              .template getField<VELOCITY_NUMERATOR>();
      T const expectedVelocityDenominator = getPorosity(particle, physR.data());
      T const expectedPorosity            = 1 - expectedVelocityDenominator;

      EXPECT_EQ(expectedPorosity, porosity);
      for (unsigned iDim = 0; iDim < D; ++iDim) {
        EXPECT_EQ(expectedVelocityDenominator, velocityDenominator[iDim]);
        EXPECT_EQ(0., velocityNumerator[iDim]);
      }
    });
  }
}

TEST(testBlockLatticeInteraction, testParticleInsideWithVelocity2D)
{
  constexpr unsigned D            = 2;
  int                nCell[D]     = {4, 4};
  T const            radius       = 2;
  T const            eps          = 1;
  unsigned const     numberBlocks = 2;

  // Creation of lattice
  CuboidDecomposition2D<T>      cuboidDecomposition(0., converter2D.getPhysDeltaX(),
                                          nCell, numberBlocks);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, D>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D2D> sLattice(converter2D, superGeometry);

  ParticleSystem<T, ResolvedParticle2D> particleSystem;
  Vector<T, D>                          pos(0);
  Vector<T, D>                          vel(0.1, 0.2);
  particles::creators::addResolvedCircle2D(particleSystem, pos, radius, eps, 0.,
                                           vel);
  auto particle = particleSystem.get(0);

  //Create eccentric velocity field
  EccentricLatticeVelocityField<T, T, D2D> eccentricVelocityField(
      particle.getField<GENERAL, POSITION>(),
      particle.getField<MOBILITY, VELOCITY>(),
      particle.getField<MOBILITY, ANG_VELOCITY>(), converter2D);

  //Write particle field
  setSuperParticleField(superGeometry, sLattice, converter2D, particle,
                        Vector<bool, D>(false));

  for (unsigned iC = 0; iC < numberBlocks; ++iC) {
    sLattice.getBlock(iC).forCoreSpatialLocations([&](LatticeR<D> latticeR) {
      PhysR<T, D> physR = superGeometry.getBlockGeometry(iC).getPhysR(latticeR);
      T           porosity =
          sLattice.getBlock(iC).get(latticeR).template getField<POROSITY>();
      Vector<T, D> velocityDenominator =
          sLattice.getBlock(iC)
              .get(latticeR)
              .template getField<VELOCITY_DENOMINATOR>();
      Vector<T, D> velocityNumerator =
          sLattice.getBlock(iC)
              .get(latticeR)
              .template getField<VELOCITY_NUMERATOR>();
      T const expectedVelocityDenominator = getPorosity(particle, physR.data());
      T const expectedPorosity            = 1 - expectedVelocityDenominator;
      T       expectedVelocityNumerator[3] = {0., 0., 0.};
      if (!util::nearZero(expectedVelocityDenominator)) {
        eccentricVelocityField(expectedVelocityNumerator, physR.data());
        //Apply weighting by porosity
        for (unsigned iDim = 0; iDim != D; ++iDim) {
          expectedVelocityNumerator[iDim] *= 1 - porosity;
        }
      }

      EXPECT_EQ(expectedPorosity, porosity);
      for (unsigned iDim = 0; iDim < D; ++iDim) {
        EXPECT_EQ(expectedVelocityDenominator, velocityDenominator[iDim]);
        EXPECT_NEAR(expectedVelocityNumerator[iDim], velocityNumerator[iDim],
                    1e-16);
      }
    });
  }
}

TEST(testBlockLatticeInteraction,
     testParticleInsideWithVelocityAndAngularVelocity2D)
{
  constexpr unsigned D            = 2;
  int                nCell[D]     = {4, 4};
  T const            radius       = 2;
  T const            eps          = 1;
  unsigned const     numberBlocks = 2;

  // Creation of lattice
  CuboidDecomposition2D<T>      cuboidDecomposition(0., converter2D.getPhysDeltaX(),
                                          nCell, numberBlocks);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, D>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D2D> sLattice(converter2D, superGeometry);

  ParticleSystem<T, ResolvedParticle2D> particleSystem;
  Vector<T, D>                          pos(0);
  Vector<T, D>                          vel(0.2, 0.3);
  particles::creators::addResolvedCircle2D(particleSystem, pos, radius, eps, 0.,
                                           vel);
  auto particle = particleSystem.get(0);
  particle.template setField<MOBILITY, ANG_VELOCITY>(0.1);

  //Create eccentric velocity field
  EccentricLatticeVelocityField<T, T, D2D> eccentricVelocityField(
      particle.getField<GENERAL, POSITION>(),
      particle.getField<MOBILITY, VELOCITY>(),
      particle.getField<MOBILITY, ANG_VELOCITY>(), converter2D);

  //Write particle field
  setSuperParticleField(superGeometry, sLattice, converter2D, particle,
                        Vector<bool, D>(false));

  for (unsigned iC = 0; iC < numberBlocks; ++iC) {
    sLattice.getBlock(iC).forCoreSpatialLocations([&](LatticeR<D> latticeR) {
      PhysR<T, D> physR = superGeometry.getBlockGeometry(iC).getPhysR(latticeR);
      T           porosity =
          sLattice.getBlock(iC).get(latticeR).template getField<POROSITY>();
      Vector<T, D> velocityDenominator =
          sLattice.getBlock(iC)
              .get(latticeR)
              .template getField<VELOCITY_DENOMINATOR>();
      Vector<T, D> velocityNumerator =
          sLattice.getBlock(iC)
              .get(latticeR)
              .template getField<VELOCITY_NUMERATOR>();
      T const expectedVelocityDenominator = getPorosity(particle, physR.data());
      T const expectedPorosity            = 1 - expectedVelocityDenominator;
      T       expectedVelocityNumerator[3] = {0., 0., 0.};
      if (!util::nearZero(expectedVelocityDenominator)) {
        eccentricVelocityField(expectedVelocityNumerator, physR.data());
        //Apply weighting by porosity
        for (unsigned iDim = 0; iDim != D; ++iDim) {
          expectedVelocityNumerator[iDim] *= 1 - porosity;
        }
      }

      EXPECT_EQ(expectedPorosity, porosity);
      for (unsigned iDim = 0; iDim < D; ++iDim) {
        EXPECT_EQ(expectedVelocityDenominator, velocityDenominator[iDim]);
        EXPECT_NEAR(expectedVelocityNumerator[iDim], velocityNumerator[iDim],
                    1e-16);
      }
    });
  }
}

TEST(testBlockLatticeInteraction,
     testTwoParticleInsideWithVelocityAndAngularVelocity2D)
{
  constexpr unsigned D            = 2;
  int                nCell[D]     = {4, 4};
  T const            radius       = 2;
  T const            eps          = 1;
  unsigned const     numberBlocks = 2;

  // Creation of lattice
  CuboidDecomposition2D<T>      cuboidDecomposition(0., converter2D.getPhysDeltaX(),
                                          nCell, numberBlocks);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, D>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D2D> sLattice(converter2D, superGeometry);

  ParticleSystem<T, ResolvedParticle2D> particleSystem;
  Vector<T, D>                          pos(0);
  Vector<T, D>                          vel(0.2, 0.3);
  particles::creators::addResolvedCircle2D(particleSystem, pos, radius, eps, 0.,
                                           vel);
  particles::creators::addResolvedCircle2D(particleSystem, pos, radius, eps, 0.,
                                           vel);
  auto particle = particleSystem.get(0);
  particle.template setField<MOBILITY, ANG_VELOCITY>(0.1);
  auto particle2 = particleSystem.get(1);

  //Create eccentric velocity field
  EccentricLatticeVelocityField<T, T, D2D> eccentricVelocityField(
      particle.getField<GENERAL, POSITION>(),
      particle.getField<MOBILITY, VELOCITY>(),
      particle.getField<MOBILITY, ANG_VELOCITY>(), converter2D);
  EccentricLatticeVelocityField<T, T, D2D> eccentricVelocityField2(
      particle2.getField<GENERAL, POSITION>(),
      particle2.getField<MOBILITY, VELOCITY>(),
      particle2.getField<MOBILITY, ANG_VELOCITY>(), converter2D);

  //Write particle field
  setSuperParticleField(superGeometry, sLattice, converter2D, particle,
                        Vector<bool, D>(false));
  setSuperParticleField(superGeometry, sLattice, converter2D, particle2,
                        Vector<bool, D>(false));

  for (unsigned iC = 0; iC < numberBlocks; ++iC) {
    sLattice.getBlock(iC).forCoreSpatialLocations([&](LatticeR<D> latticeR) {
      PhysR<T, D> physR = superGeometry.getBlockGeometry(iC).getPhysR(latticeR);
      T           porosity =
          sLattice.getBlock(iC).get(latticeR).template getField<POROSITY>();
      Vector<T, D> velocityDenominator =
          sLattice.getBlock(iC)
              .get(latticeR)
              .template getField<VELOCITY_DENOMINATOR>();
      Vector<T, D> velocityNumerator =
          sLattice.getBlock(iC)
              .get(latticeR)
              .template getField<VELOCITY_NUMERATOR>();
      T expectedVelocityDenominator = getPorosity(particle, physR.data()) +
                                      getPorosity(particle2, physR.data());
      T expectedPorosity = (1 - getPorosity(particle, physR.data())) *
                           (1 - getPorosity(particle2, physR.data()));
      T expectedVelocityNumerator[3] = {0., 0., 0.};
      if (!util::nearZero(expectedVelocityDenominator)) {
        T vel2[D];
        eccentricVelocityField(expectedVelocityNumerator, physR.data());
        eccentricVelocityField2(vel2, physR.data());
        //Apply weighting by porosity
        for (unsigned iDim = 0; iDim != D; ++iDim) {
          expectedVelocityNumerator[iDim] *=
              getPorosity(particle, physR.data());
          vel2[iDim] *= getPorosity(particle, physR.data());
          expectedVelocityNumerator[iDim] += vel2[iDim];
        }
      }

      EXPECT_EQ(expectedPorosity, porosity);
      for (unsigned iDim = 0; iDim < D; ++iDim) {
        EXPECT_EQ(expectedVelocityDenominator, velocityDenominator[iDim]);
        EXPECT_NEAR(expectedVelocityNumerator[iDim], velocityNumerator[iDim],
                    1e-16);
      }
    });
  }
}

// 3D Tests
TEST(testBlockLatticeInteraction, testParticleFarOut3D)
{
  constexpr unsigned D            = 3;
  int                nCell[D]     = {4, 4, 4};
  T const            radius       = 2;
  T const            eps          = 1;
  unsigned const     numberBlocks = 2;

  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCell,
                                     numberBlocks);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, D>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D,superGeometry);

  ParticleSystem<T, ResolvedParticle3D> particleSystem;
  Vector<T, D>                          pos(-500);
  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, eps);
  auto particle = particleSystem.get(0);

  //Create eccentric velocity field
  EccentricLatticeVelocityField<T, T, D3D> eccentricVelocityField(
      particle.getField<GENERAL, POSITION>(),
      particle.getField<MOBILITY, VELOCITY>(),
      particle.getField<MOBILITY, ANG_VELOCITY>(), converter3D);

  //Write particle field
  setSuperParticleField(superGeometry, sLattice, converter3D, particle,
                        Vector<bool, D>(false));

  for (unsigned iC = 0; iC < numberBlocks; ++iC) {
    sLattice.getBlock(iC).forCoreSpatialLocations([&](LatticeR<D> latticeR) {
      T porosity =
          sLattice.getBlock(iC).get(latticeR).template getField<POROSITY>();
      Vector<T, D> velocityDenominator =
          sLattice.getBlock(iC)
              .get(latticeR)
              .template getField<VELOCITY_DENOMINATOR>();
      Vector<T, D> velocityNumerator =
          sLattice.getBlock(iC)
              .get(latticeR)
              .template getField<VELOCITY_NUMERATOR>();

      EXPECT_EQ(1., porosity);
      for (unsigned iDim = 0; iDim < D; ++iDim) {
        EXPECT_EQ(0., velocityDenominator[iDim]);
        EXPECT_EQ(0., velocityNumerator[iDim]);
      }
    });
  }
}

TEST(testBlockLatticeInteraction, testParticleInsideWithZeroVelocity3D)
{
  constexpr unsigned D            = 3;
  int                nCell[D]     = {4, 4, 4};
  T const            radius       = 2;
  T const            eps          = 1;
  unsigned const     numberBlocks = 2;

  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCell,
                                     numberBlocks);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, D>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D, superGeometry);

  ParticleSystem<T, ResolvedParticle3D> particleSystem;
  Vector<T, D>                          pos(0);
  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, eps);
  auto particle = particleSystem.get(0);

  //Create eccentric velocity field
  EccentricLatticeVelocityField<T, T, D3D> eccentricVelocityField(
      particle.getField<GENERAL, POSITION>(),
      particle.getField<MOBILITY, VELOCITY>(),
      particle.getField<MOBILITY, ANG_VELOCITY>(), converter3D);

  //Write particle field
  setSuperParticleField(superGeometry, sLattice, converter3D, particle,
                        Vector<bool, D>(false));

  for (unsigned iC = 0; iC < numberBlocks; ++iC) {
    sLattice.getBlock(iC).forCoreSpatialLocations([&](LatticeR<D> latticeR) {
      PhysR<T, D> physR = superGeometry.getBlockGeometry(iC).getPhysR(latticeR);
      T           porosity =
          sLattice.getBlock(iC).get(latticeR).template getField<POROSITY>();
      Vector<T, D> velocityDenominator =
          sLattice.getBlock(iC)
              .get(latticeR)
              .template getField<VELOCITY_DENOMINATOR>();
      Vector<T, D> velocityNumerator =
          sLattice.getBlock(iC)
              .get(latticeR)
              .template getField<VELOCITY_NUMERATOR>();
      T const expectedVelocityDenominator = getPorosity(particle, physR.data());
      T const expectedPorosity            = 1 - expectedVelocityDenominator;

      EXPECT_EQ(expectedPorosity, porosity);
      for (unsigned iDim = 0; iDim < D; ++iDim) {
        EXPECT_EQ(expectedVelocityDenominator, velocityDenominator[iDim]);
        EXPECT_EQ(0., velocityNumerator[iDim]);
      }
    });
  }
}

TEST(testBlockLatticeInteraction, testParticleInsideWithVelocity3D)
{
  constexpr unsigned D            = 3;
  int                nCell[D]     = {4, 4, 4};
  T const            radius       = 2;
  T const            eps          = 1;
  unsigned const     numberBlocks = 2;

  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCell,
                                     numberBlocks);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, D>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D, superGeometry);

  ParticleSystem<T, ResolvedParticle3D> particleSystem;
  Vector<T, D>                          pos(0);
  Vector<T, D>                          vel(0.1, 0.2, 0.3);
  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, eps, 0.,
                                           vel);
  auto particle = particleSystem.get(0);

  //Create eccentric velocity field
  EccentricLatticeVelocityField<T, T, D3D> eccentricVelocityField(
      particle.getField<GENERAL, POSITION>(),
      particle.getField<MOBILITY, VELOCITY>(),
      particle.getField<MOBILITY, ANG_VELOCITY>(), converter3D);

  //Write particle field
  setSuperParticleField(superGeometry, sLattice, converter3D, particle,
                        Vector<bool, D>(false));

  for (unsigned iC = 0; iC < numberBlocks; ++iC) {
    sLattice.getBlock(iC).forCoreSpatialLocations([&](LatticeR<D> latticeR) {
      PhysR<T, D> physR = superGeometry.getBlockGeometry(iC).getPhysR(latticeR);
      T           porosity =
          sLattice.getBlock(iC).get(latticeR).template getField<POROSITY>();
      Vector<T, D> velocityDenominator =
          sLattice.getBlock(iC)
              .get(latticeR)
              .template getField<VELOCITY_DENOMINATOR>();
      Vector<T, D> velocityNumerator =
          sLattice.getBlock(iC)
              .get(latticeR)
              .template getField<VELOCITY_NUMERATOR>();
      T const expectedVelocityDenominator = getPorosity(particle, physR.data());
      T const expectedPorosity            = 1 - expectedVelocityDenominator;
      T       expectedVelocityNumerator[3] = {0., 0., 0.};
      if (!util::nearZero(expectedVelocityDenominator)) {
        eccentricVelocityField(expectedVelocityNumerator, physR.data());
        //Apply weighting by porosity
        for (unsigned iDim = 0; iDim != D; ++iDim) {
          expectedVelocityNumerator[iDim] *= 1 - porosity;
        }
      }

      EXPECT_EQ(expectedPorosity, porosity);
      for (unsigned iDim = 0; iDim < D; ++iDim) {
        EXPECT_EQ(expectedVelocityDenominator, velocityDenominator[iDim]);
        EXPECT_NEAR(expectedVelocityNumerator[iDim], velocityNumerator[iDim],
                    1e-16);
      }
    });
  }
}

TEST(testBlockLatticeInteraction,
     testParticleInsideWithVelocityAndAngularVelocity3D)
{
  constexpr unsigned D            = 3;
  int                nCell[D]     = {4, 4, 4};
  T const            radius       = 2;
  T const            eps          = 1;
  unsigned const     numberBlocks = 2;

  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCell,
                                     numberBlocks);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, D>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D, superGeometry);

  ParticleSystem<T, ResolvedParticle3D> particleSystem;
  Vector<T, D>                          pos(0);
  Vector<T, D>                          vel(0.1, 0.2, 0.3);
  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, eps, 0.,
                                           vel);
  auto particle = particleSystem.get(0);
  particle.template setField<MOBILITY, ANG_VELOCITY>(vel);

  //Create eccentric velocity field
  EccentricLatticeVelocityField<T, T, D3D> eccentricVelocityField(
      particle.getField<GENERAL, POSITION>(),
      particle.getField<MOBILITY, VELOCITY>(),
      particle.getField<MOBILITY, ANG_VELOCITY>(), converter3D);

  //Write particle field
  setSuperParticleField(superGeometry, sLattice, converter3D, particle,
                        Vector<bool, D>(false));

  for (unsigned iC = 0; iC < numberBlocks; ++iC) {
    sLattice.getBlock(iC).forCoreSpatialLocations([&](LatticeR<D> latticeR) {
      PhysR<T, D> physR = superGeometry.getBlockGeometry(iC).getPhysR(latticeR);
      T           porosity =
          sLattice.getBlock(iC).get(latticeR).template getField<POROSITY>();
      Vector<T, D> velocityDenominator =
          sLattice.getBlock(iC)
              .get(latticeR)
              .template getField<VELOCITY_DENOMINATOR>();
      Vector<T, D> velocityNumerator =
          sLattice.getBlock(iC)
              .get(latticeR)
              .template getField<VELOCITY_NUMERATOR>();
      T const expectedVelocityDenominator = getPorosity(particle, physR.data());
      T const expectedPorosity            = 1 - expectedVelocityDenominator;
      T       expectedVelocityNumerator[3] = {0., 0., 0.};
      if (!util::nearZero(expectedVelocityDenominator)) {
        eccentricVelocityField(expectedVelocityNumerator, physR.data());
        //Apply weighting by porosity
        for (unsigned iDim = 0; iDim != D; ++iDim) {
          expectedVelocityNumerator[iDim] *= 1 - porosity;
        }
      }

      EXPECT_EQ(expectedPorosity, porosity);
      for (unsigned iDim = 0; iDim < D; ++iDim) {
        EXPECT_EQ(expectedVelocityDenominator, velocityDenominator[iDim]);
        EXPECT_NEAR(expectedVelocityNumerator[iDim], velocityNumerator[iDim],
                    1e-16);
      }
    });
  }
}

TEST(testBlockLatticeInteraction,
     testTwoParticleInsideWithVelocityAndAngularVelocity3D)
{
  constexpr unsigned D            = 3;
  int                nCell[D]     = {4, 4, 4};
  T const            radius       = 2;
  T const            eps          = 1;
  unsigned const     numberBlocks = 2;

  // Creation of lattice
  CuboidDecomposition3D<T> cuboidDecomposition(0., converter3D.getPhysDeltaX(),
                                     nCell,
                                     numberBlocks);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  SuperGeometry<T, D>      superGeometry(cuboidDecomposition, loadBalancer);
  superGeometry.rename(0, 1);
  SuperLattice<T, D3D> sLattice(converter3D, superGeometry);

  ParticleSystem<T, ResolvedParticle3D> particleSystem;
  Vector<T, D>                          pos(0);
  Vector<T, D>                          vel(0.1, 0.2, 0.3);
  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, eps, 0.,
                                           vel);
  particles::creators::addResolvedSphere3D(particleSystem, pos, radius, eps, 0.,
                                           vel);
  auto particle = particleSystem.get(0);
  particle.template setField<MOBILITY, ANG_VELOCITY>(vel);
  auto particle2 = particleSystem.get(1);

  //Create eccentric velocity field
  EccentricLatticeVelocityField<T, T, D3D> eccentricVelocityField(
      particle.getField<GENERAL, POSITION>(),
      particle.getField<MOBILITY, VELOCITY>(),
      particle.getField<MOBILITY, ANG_VELOCITY>(), converter3D);
  EccentricLatticeVelocityField<T, T, D3D> eccentricVelocityField2(
      particle2.getField<GENERAL, POSITION>(),
      particle2.getField<MOBILITY, VELOCITY>(),
      particle2.getField<MOBILITY, ANG_VELOCITY>(), converter3D);

  //Write particle field
  setSuperParticleField(superGeometry, sLattice, converter3D, particle,
                        Vector<bool, D>(false));
  setSuperParticleField(superGeometry, sLattice, converter3D, particle2,
                        Vector<bool, D>(false));

  for (unsigned iC = 0; iC < numberBlocks; ++iC) {
    sLattice.getBlock(iC).forCoreSpatialLocations([&](LatticeR<D> latticeR) {
      PhysR<T, D> physR = superGeometry.getBlockGeometry(iC).getPhysR(latticeR);
      T           porosity =
          sLattice.getBlock(iC).get(latticeR).template getField<POROSITY>();
      T velocityDenominator = sLattice.getBlock(iC)
                                  .get(latticeR)
                                  .template getField<VELOCITY_DENOMINATOR>();
      Vector<T, D> velocityNumerator =
          sLattice.getBlock(iC)
              .get(latticeR)
              .template getField<VELOCITY_NUMERATOR>();
      T expectedVelocityDenominator = getPorosity(particle, physR.data()) +
                                      getPorosity(particle2, physR.data());
      T expectedPorosity = (1 - getPorosity(particle, physR.data())) *
                           (1 - getPorosity(particle2, physR.data()));
      T expectedVelocityNumerator[3] = {0., 0., 0.};
      if (!util::nearZero(expectedVelocityDenominator)) {
        T vel2[D];
        eccentricVelocityField(expectedVelocityNumerator, physR.data());
        eccentricVelocityField2(vel2, physR.data());
        //Apply weighting by porosity
        for (unsigned iDim = 0; iDim != D; ++iDim) {
          expectedVelocityNumerator[iDim] *=
              getPorosity(particle, physR.data());
          vel2[iDim] *= getPorosity(particle, physR.data());
          expectedVelocityNumerator[iDim] += vel2[iDim];
        }
      }

      EXPECT_EQ(expectedPorosity, porosity);
      EXPECT_EQ(expectedVelocityDenominator, velocityDenominator);
      for (unsigned iDim = 0; iDim < D; ++iDim) {
        EXPECT_NEAR(expectedVelocityNumerator[iDim], velocityNumerator[iDim],
                    1e-16);
      }
    });
  }
}
