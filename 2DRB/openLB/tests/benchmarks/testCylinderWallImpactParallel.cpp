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
#ifdef PARALLEL_MODE_MPI

#include <olb.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "../unit/helper/framework.h"

#include <gtest/gtest.h>

namespace olb {
namespace descriptors {
struct OLD_VELOCITY : public FIELD_BASE<0, 1, 0> {};
struct OLD_OLD_VELOCITY : public FIELD_BASE<0, 1, 0> {};
struct OLD_ANG_VELOCITY : public FIELD_BASE<0, 1, 0> {};
struct OLD_OLD_ANG_VELOCITY : public FIELD_BASE<0, 1, 0> {};
struct CUSTOM {};
template <unsigned D>
struct CUSTOMFIELDS
    : public PARTICLE_DESCRIPTOR<D, OLD_VELOCITY, OLD_OLD_VELOCITY,
                                 OLD_ANG_VELOCITY, OLD_OLD_ANG_VELOCITY>
    , public CUSTOM {};
} // namespace descriptors
} // namespace olb

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::particles;
using namespace olb::particles::dynamics;
using namespace olb::particles::contact;
using namespace olb::particles::communication;
using namespace olb::particles::access;
using namespace olb::util;

typedef double T;

//Define lattice type
typedef D3Q19<POROSITY, VELOCITY_NUMERATOR, VELOCITY_DENOMINATOR,
              CONTACT_DETECTION>
    DESCRIPTOR;
//Define particle type
typedef PARTICLE_DESCRIPTOR<
    DESCRIPTOR::d, GENERAL_EXTENDABLE<DESCRIPTOR::d>,
    MOBILITY_VERLET<DESCRIPTOR::d>, SURFACE_RESOLVED_PARALLEL<DESCRIPTOR::d>,
    FORCING_RESOLVED<DESCRIPTOR::d>, PHYSPROPERTIES_RESOLVED<DESCRIPTOR::d>,
    MECHPROPERTIES_COLLISION<DESCRIPTOR::d>,
    NUMERICPROPERTIES_RESOLVED_CONTACT<DESCRIPTOR::d>,
    CUSTOMFIELDS<DESCRIPTOR::d>, PARALLELIZATION_RESOLVED>
    PARTICLETYPE;
// Define particle-particle contact type
typedef ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
    PARTICLECONTACTTYPE;
// Define particle-wall contact type
typedef WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>
    WALLCONTACTTYPE;

/// ensure olbInit is called for MPI
::testing::Environment* const foo_env =
    ::testing::AddGlobalTestEnvironment(new OlbInitTestEnvironment);

struct TestParam {
  T angle;
  T allowedErrorVel;
  T allowedErrorAngularVel;
};

/*
angle: 0;   velocity: 0.08500000000000000   angular velocity: 0
angle: 10;  velocity: -0.03639435700970071  angular velocity: -34.89282749498069
angle: 20;  velocity: -0.01896392232442171  angular velocity: -36.44765774368297
angle: 30;  velocity: 0.007463935608702979  angular velocity: -36.24702250823359
angle: 40;  velocity: 0.04319049868509287   angular velocity: -30.72441939691441
angle: 50;  velocity: 0.07683504100490955   angular velocity: -15.08864559659039
angle: 60;  velocity: 0.08249631274130925   angular velocity: 8.488015751869153
angle: 70;  velocity: 0.05435463777033014   angular velocity: 27.31060438322266
angle: 80;  velocity: 0.01721179789900002   angular velocity: 35.39580645875379
angle: 90;  velocity: 0.08500000000000000   angular velocity: 0)
*/

TestParam params[] = {
    { 0,  0.0035, 0.60},
    {10,  0.0015, 0.55},
    {20,  0.0012,  0.2},
    {30,  0.0007, 0.25},
    {40,  0.0004, 0.32},
    {50, 0.00085, 0.75},
    {60,  0.0004,  1.8},
    {70,   0.002,  1.5},
    {80,   0.001, 0.65},
    {90, 0.00225, 0.45}
};

class CylinderWallImpactTest : public ::testing::TestWithParam<TestParam> {};

// Discretization Settings
int const res                 = 8;
T const   charLatticeVelocity = 0.0019;

// Time Settings
T const        maxPhysT     = 0.015;  // max. simulation time in s
T const        iTwrite      = 0.0004; // write out intervall in s
unsigned const subTimesteps = 10;

// Domain Settings
T const lengthX = 0.02;
T const lengthY = 0.012;
T const lengthZ = 0.02;

// Fluid Settings
T const physDensity        = 1.2922;
T const dynamicViscosity   = 18.2e-6;
T const kinematicViscosity = dynamicViscosity / physDensity;

//Particle Settings
T const radius    = 4e-3;
T const cylHeight = 5.3e-3;
T const distanceImpactAndCenter =
    util::sqrt(radius * radius + 0.25 * cylHeight * cylHeight);
T const density = 3.1e-4 / (util::pow(radius, 2) * M_PI * cylHeight);

Vector<T, DESCRIPTOR::d> cylAngle(0., 0., 0.);
Vector<T, DESCRIPTOR::d> cylVelocity(0., -0.1, 0.);
Vector<T, DESCRIPTOR::d> cylNormal(0., 1., 0.);
Vector<T, DESCRIPTOR::d> cylCenter(lengthX * .5,
                                   util::sin(M_PI* cylAngle[2] / 180) * radius +
                                       util::cos(M_PI * cylAngle[2] / 180) *
                                           0.5 * cylHeight +
                                       1.1 * 2 * radius / res,
                                   lengthZ * .5);
Vector<T, DESCRIPTOR::d> externalAcceleration(0.);

// Contact Settings
constexpr T const signedCoeffRestitution = 0.85;
constexpr T const YoungsModulus          = 1.0e7;
constexpr T const PoissonRation          = 0.4;

// Characteristic Quantities
T const charPhysLength = 2 * radius;
T charPhysVelocity     = util::abs(cylVelocity[1]); // assumed maximal velocity;

unsigned const contactResolutionPerDirection = 9;

Vector<T, 3> oldAngularVelocity(0.);
Vector<T, 3> oldVelocity(0.);

T angleBetween(T x)
{
  return (x <= 0 ? 90
          : x >= 90
              ? 0
              : util::acos(radius / distanceImpactAndCenter) * 180. / M_PI);
}

T sumAngles(T x) { return x + angleBetween(x); }

T cosSumAngles(T x) { return util::cos(util::degreeToRadian(sumAngles(x))); }

T angularVel(T x)
{
  T const mass = 0.00031;
  T const mofi = 1.966e-9;
  return mass * cylVelocity[1] * (1 + signedCoeffRestitution) *
         distanceImpactAndCenter * cosSumAngles(x) /
         (mofi + mass * distanceImpactAndCenter * distanceImpactAndCenter *
                     cosSumAngles(x) * cosSumAngles(x));
}

T vel(T x)
{
  return angularVel(x) * distanceImpactAndCenter * cosSumAngles(x) -
         signedCoeffRestitution * cylVelocity[1];
}

bool isFinished(XParticleSystem<T, PARTICLETYPE>& particleSystem)
{
  const auto evalParticle = [&](Particle<T, PARTICLETYPE>& particle) {
    const Vector<T, 3> oldVelocity(
        particle.template getField<CUSTOM, OLD_VELOCITY>());
    const Vector<T, 3> oldoldVelocity(
        particle.template getField<CUSTOM, OLD_OLD_VELOCITY>());
    const Vector<T, 3> oldAngularVelocity(
        particle.template getField<CUSTOM, OLD_ANG_VELOCITY>());
    const Vector<T, 3> oldoldAngularVelocity(
        particle.template getField<CUSTOM, OLD_OLD_ANG_VELOCITY>());

    const Vector<T, 3> velocity(
        particle.template getField<MOBILITY, VELOCITY>());
    const Vector<T, 3> force(particle.template getField<FORCING, FORCE>());
    const Vector<T, 3> angularVelocity(
        particle.template getField<MOBILITY, ANG_VELOCITY>());

    bool const angularVelDone = angularVelocity == oldAngularVelocity &&
                                oldoldAngularVelocity == angularVelocity;
    bool const velDone = velocity != cylVelocity && velocity == oldVelocity &&
                         velocity == oldoldVelocity;
    bool const forceZero = force == Vector<T, 3>(0.);

    particle.template setField<CUSTOM, OLD_OLD_ANG_VELOCITY>(
        oldAngularVelocity);
    particle.template setField<CUSTOM, OLD_ANG_VELOCITY>(angularVelocity);
    particle.template setField<CUSTOM, OLD_OLD_VELOCITY>(oldVelocity);
    particle.template setField<CUSTOM, OLD_VELOCITY>(velocity);

    return angularVelDone && velDone && forceZero;
  };

#ifdef FEATURE_PARALLEL
  int numberOnTrue = 555;
  int number       = 0;

  communication::forParticlesInSuperParticleSystem(
      particleSystem,
      [&](Particle<T, PARTICLETYPE>&       particle,
          ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
        const std::size_t particleIC =
            particle.template getField<PARALLELIZATION, IC>();
        if (particleIC == globiC) {
          if (evalParticle(particle)) {
            number = numberOnTrue;
          }
        }
      });
  singleton::mpi().reduceAndBcast(number, MPI_MAX);
  return number == numberOnTrue;
#else
  auto particle = particleSystem.get(0);
  return evalParticle(particle);
#endif
  __builtin_unreachable();
}

void resetForceOnParticle(SuperParticleSystem<T, PARTICLETYPE>& particleSystem)
{
  particleSystem.get(0).template setField<FORCING, FORCE>(
      Vector<T, DESCRIPTOR::d>(0.));
  particleSystem.get(0).template setField<FORCING, TORQUE>(
      Vector<T, DESCRIPTOR::d>(0.));
}

// Prepare geometry
void prepareGeometry(UnitConverter<T, DESCRIPTOR> const& converter,
                     SuperGeometry<T, DESCRIPTOR::d>&    superGeometry,
                     SolidBoundary<T, DESCRIPTOR::d>&    wall)
{
  superGeometry.rename(0, 1);
  superGeometry.rename(1, wall.getLatticeMaterial(), *wall.getIndicator());
}

// Set up the geometry of the simulation
void prepareLattice(SuperLattice<T, DESCRIPTOR>&        sLattice,
                    UnitConverter<T, DESCRIPTOR> const& converter,
                    SuperGeometry<T, DESCRIPTOR::d>&    superGeometry)
{
  /// Material=0 -->do nothing
  sLattice.defineDynamics<NoDynamics>(superGeometry, 0);
  sLattice.defineDynamics<PorousParticleBGKdynamics>(superGeometry, 1);
  sLattice.defineDynamics<BounceBack>(superGeometry, 2);

  sLattice.setParameter<descriptors::OMEGA>(
      converter.getLatticeRelaxationFrequency());

  {
    auto& communicator = sLattice.getCommunicator(stage::PostPostProcess());
    communicator
        .requestFields<POROSITY, VELOCITY_NUMERATOR, VELOCITY_DENOMINATOR>();
    communicator.requestOverlap(sLattice.getOverlap());
    communicator.exchangeRequests();
  }
}

//Set Boundary Values
void setBoundaryValues(SuperLattice<T, DESCRIPTOR>&        sLattice,
                       UnitConverter<T, DESCRIPTOR> const& converter, int iT,
                       SuperGeometry<T, DESCRIPTOR::d>& superGeometry)
{
  if (iT == 0) {
    AnalyticalConst<DESCRIPTOR::d, T, T> zero(0.);
    AnalyticalConst<DESCRIPTOR::d, T, T> one(1.);
    sLattice.defineField<descriptors::POROSITY>(
        superGeometry.getMaterialIndicator({0, 1, 2}), one);
    // Set initial condition
    AnalyticalConst<DESCRIPTOR::d, T, T> ux(0.);
    AnalyticalConst<DESCRIPTOR::d, T, T> uy(0.);
    AnalyticalConst3D<T, T>              uz(0.);
    AnalyticalComposed3D<T, T>           u(ux, uy, uz);
    AnalyticalConst<DESCRIPTOR::d, T, T> rho(1.);

    //Initialize all values of distribution functions to their local equilibrium
    sLattice.defineRhoU(superGeometry, 1, rho, u);
    sLattice.iniEquilibrium(superGeometry, 1, rho, u);

    // Make the lattice ready for simulation
    sLattice.initialize();
  }
}

TEST_P(CylinderWallImpactTest, testCylWallImpactParallel)
{
  // === 1st Step: Initialization ===
  TestParam param = GetParam();

  cylAngle  = Vector<T, DESCRIPTOR::d>(0., 0., param.angle);
  cylCenter = Vector<T, DESCRIPTOR::d>(
      lengthX * .5,
      util::sin(M_PI * cylAngle[2] / 180) * radius +
          util::cos(M_PI * cylAngle[2] / 180) * 0.5 * cylHeight +
          1.1 * 2 * radius / res,
      lengthZ * .5);

  UnitConverterFromResolutionAndLatticeVelocity<T, DESCRIPTOR> converter(
      (int)res,               //resolution
      (T)charLatticeVelocity, //charLatticeVelocity
      (T)charPhysLength,      //charPhysLength
      (T)charPhysVelocity,    //charPhysVelocity
      (T)kinematicViscosity,  //physViscosity
      (T)physDensity          //physDensity
  );

  /// === 2rd Step: Prepare Geometry ===
  /// Instantiation of a cuboidDecomposition with weights
  Vector<T, DESCRIPTOR::d> origin(0.);
  Vector<T, DESCRIPTOR::d> extend;
  extend[0] = lengthX;
  extend[1] = lengthY;
  extend[2] = lengthZ;
  IndicatorCuboid3D<T> cuboid(extend, origin);
  // Setting (big) bounding box is necessary for distance calculation
  std::unique_ptr<IndicatorF<T, DESCRIPTOR::d>> wallInd =
      std::make_unique<IndicInverse<T, DESCRIPTOR::d>>(
          cuboid, origin - 5 * charPhysLength / res,
          extend + 5 * charPhysLength / res);
  // Create vector containing all solidBoundaries
  std::vector<SolidBoundary<T, 3>> solidBoundaries;
  solidBoundaries.push_back(
      SolidBoundary<T, DESCRIPTOR::d>(std::move(wallInd), 2, 0));

  CuboidDecomposition3D<T> cuboidDecomposition(cuboid,
                                     converter.getConversionFactorLength(), 1);

  HeuristicLoadBalancer<T>        loadBalancer(cuboidDecomposition);
  SuperGeometry<T, DESCRIPTOR::d> superGeometry(cuboidDecomposition, loadBalancer,
                                                2);
  prepareGeometry(converter, superGeometry, solidBoundaries[0]);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice(converter, superGeometry);

  // Prepare lattice
  prepareLattice(sLattice, converter, superGeometry);

  // Create ParticleSystem
  SuperParticleSystem<T, PARTICLETYPE> sParticleSystem(superGeometry);
  // Create ParticleContact container
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  // Container containg contact properties
  constexpr ContactProperties<T, 1> contactProperties =
      createContactProperties<T, 1>(
          [&](ContactProperties<T, 1>& contactProperties) {
            contactProperties.set(
                0, 0,
                evalEffectiveYoungModulus(YoungsModulus, YoungsModulus,
                                          PoissonRation, PoissonRation),
                signedCoeffRestitution, 0, 0);
          });

  // Calculate particle quantities
  T epsilon = 0.5 * converter.getConversionFactorLength();

  // Create Particle 1
  creators::addResolvedCylinder3D(sParticleSystem, cylCenter, cylNormal,
                                  cylHeight, radius, epsilon, density, cylAngle,
                                  cylVelocity);
  forParticlesInSuperParticleSystem(
      sParticleSystem,
      [](Particle<T, PARTICLETYPE>&       particle,
         ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
        particle.setField<MECHPROPERTIES, MATERIAL>(0);
      });

  // Create and assign resolved particle dynamics
  sParticleSystem.defineDynamics<VerletParticleDynamics<T, PARTICLETYPE>>();

  // Create communicators
  communication::ParticleCommunicator particleCommunicator;

  /// === 4th Step: Main Loop with Timer ===
  Timer<T> timer(converter.getLatticeTime(maxPhysT),
                 superGeometry.getStatistics().getNvoxel());
  timer.start();

  /// === 5th Step: Definition of Initial and Boundary Conditions ===
  setBoundaryValues(sLattice, converter, 0, superGeometry);

  // Setup memory for saving pre contact forces
  std::vector<Vector<T, DESCRIPTOR::d>> forceWithoutContact;
  forceWithoutContact.resize(1);
  std::vector<Vector<T, DESCRIPTOR::d>> torqueWithoutContact;
  torqueWithoutContact.resize(1);

  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT) + 10; ++iT) {
    // Reset force on particle
    resetForceOnParticle(sParticleSystem);

    unsigned currSteps = 1;
    if (contactContainer.wallContacts.size() != 0 ||
        contactContainer.particleContacts.size() != 0) {
      currSteps = subTimesteps;
    }
    T currDeltaT = converter.getPhysDeltaT() / currSteps;

    // Processing with sub-timesteps
    forParticlesInSuperParticleSystem(
        sParticleSystem,
        [&forceWithoutContact, &torqueWithoutContact](
            Particle<T, PARTICLETYPE>&       particle,
            ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
          forceWithoutContact[0] = particle.template getField<FORCING, FORCE>();
          torqueWithoutContact[0] =
              particle.template getField<FORCING, TORQUE>();
        });

    for (unsigned i = 0; i < currSteps; ++i) {
      // Remove old contact force
      forParticlesInSuperParticleSystem(
          sParticleSystem,
          [&forceWithoutContact, &torqueWithoutContact](
              Particle<T, PARTICLETYPE>&       particle,
              ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
            particle.template setField<FORCING, FORCE>(forceWithoutContact[0]);
            particle.template setField<FORCING, TORQUE>(
                torqueWithoutContact[0]);
          });

      particles::contact::processContacts<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                          WALLCONTACTTYPE,
                                          ContactProperties<T, 1>>(
          sParticleSystem, solidBoundaries, contactContainer, contactProperties,
          superGeometry, particleCommunicator.contactTreatmentComm,
          contactResolutionPerDirection, 2. / util::sqrt(M_PI));

      // Process particles (Verlet algorithm)
      communication::forParticlesInSuperParticleSystem<
          T, PARTICLETYPE,
          conditions::valid_particle_centres //only consider center for resolved
          >(sParticleSystem,
            [&](Particle<T, PARTICLETYPE>&       particle,
                ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
              particle.process(currDeltaT);
            });

      particles::communication::communicateEquationsOfMotionResults(
          sParticleSystem, contactContainer,
          particleCommunicator.equationsOfMotionComm);
    }
    communicatePostContactTreatmentContacts(
        contactContainer, sParticleSystem, converter.getPhysDeltaX(),
        particleCommunicator.particleContactDetectionComm,
        particleCommunicator.wallContactDetectionComm, Vector<bool, 3>(false));

    // Couple particles to lattice
    particles::dynamics::coupleResolvedParticlesToLattice<
        T, DESCRIPTOR, PARTICLETYPE, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
        sParticleSystem, contactContainer, superGeometry, sLattice, converter,
        solidBoundaries);

    // Communicate found contacts
    communicateContactsParallel(
        contactContainer, sParticleSystem, converter.getPhysDeltaX(),
        particleCommunicator.particleContactDetectionComm,
        particleCommunicator.wallContactDetectionComm, Vector<bool, 3>(false));

    if (iT % converter.getLatticeTime(iTwrite) == 0 &&
        isFinished(sParticleSystem)) {
      break;
    }
  }

  timer.stop();

  auto         particle = sParticleSystem.get(0);
  Vector<T, 3> velocity(particle.template getField<MOBILITY, VELOCITY>());
  Vector<T, 3> angularVelocity(
      particle.template getField<MOBILITY, ANG_VELOCITY>());

  // using max in expected accuracy since some results are close to 0
  std::cout << "angle: " << param.angle << std::endl;
  EXPECT_NEAR(vel(param.angle), velocity[1], param.allowedErrorVel);
  EXPECT_NEAR(angularVel(param.angle), angularVelocity[2],
              param.allowedErrorAngularVel);
}

INSTANTIATE_TEST_CASE_P(ParticleWallInteraction, CylinderWallImpactTest,
                        testing::ValuesIn(params));

#endif
