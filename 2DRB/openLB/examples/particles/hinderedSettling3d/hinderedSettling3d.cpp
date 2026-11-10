/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Jan Marquardt
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
/*
 * This example is based on 10.1016/j.cpc.2024.109321.
 * The limestone shapes are from the online particle database PARROT
 * (https://parrot.tu-freiberg.de/).
 * The number of triangles has been reduced.
 */

#include "olb.h"

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <memory>
#include <vector>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace olb::names;


#define WithContact

// === Step 1: Declarations ===

using MyCase = Case<
                NavierStokes, Lattice<double, D3Q19<POROSITY, VELOCITY_NUMERATOR, VELOCITY_DENOMINATOR,
                                                          CONTACT_DETECTION, FORCE>>
>;

namespace olb::parameters {

struct COEFFICIENT_OF_RESTITUTION : public descriptors::FIELD_BASE<1> { };
struct COEFFICIENT_STATIC_FRICTION : public descriptors::FIELD_BASE<1> { };
struct COEFFICIENT_KINETIC_FRICTION : public descriptors::FIELD_BASE<1> { };
struct EPS : public descriptors::FIELD_BASE<1> { };
struct PART_POISSON_RATIO : public descriptors::FIELD_BASE<1> { };
struct PARTICLE_TYPE : public descriptors::TYPED_FIELD_BASE<int,1> { };
struct WANTED_PARTICLE_VOLUME_FRACTION : public descriptors::FIELD_BASE<1> { };
struct DENSITY_RATIO : public descriptors::FIELD_BASE<1> { };
struct STATIC_KINETIC_TRANSITION_VELOCITY : public descriptors::FIELD_BASE<1> { };
struct CONTACT_BOX_RESOLUTION_PER_DIRECTION : public descriptors::FIELD_BASE<1> { };
struct PARTICLE_CONTACT_MATERIAL : public descriptors::FIELD_BASE<1> { };
struct WALL_CONTACT_MATERIAL : public descriptors::FIELD_BASE<1> { };
struct EXTENT_TO_DIAMETER : public descriptors::FIELD_BASE<1> { };
struct PHYS_PURGE_ITER_T : public descriptors::FIELD_BASE<1> { };
struct PARTICLE_ENLARGEMENT_FOR_CONTACT : public descriptors::FIELD_BASE<1> { };
struct PARTICLE_NUMBER : public descriptors::TYPED_FIELD_BASE<std::size_t,1> { };
struct PARTICLE_VOLUME_FRACTION : public descriptors::FIELD_BASE<1> { };

}


Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters)
{
  using T           = MyCase::value_t;
  const T      size = parameters.get<parameters::PHYS_CHAR_LENGTH>() * parameters.get<parameters::EXTENT_TO_DIAMETER>();
  const Vector origin {0, 0, 0};
  IndicatorCuboid3D<T> cuboid({size, size, size}, origin);

  const T            physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();
  Mesh<T, MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}



#ifdef PARALLEL_MODE_MPI
MPI_Comm averageParticleVelocityComm; /// Communicator for calculation of average particle velocity
MPI_Comm numberParticleCellsComm;     /// Communicator for calculation of current number of particle cells


std::string getParticleIdentifier(const std::size_t& pID) { return std::to_string(pID); }

template <typename PARTICLESYSTEM>
MyCase::value_t evalAverageSettlingVelocity(MyCase& myCase, PARTICLESYSTEM& xParticleSystem)
{
  using T = MyCase::value_t;
  using namespace olb::particles;
  using namespace olb::particles::access;
#ifdef WithContact
  using PARTICLETYPE = ResolvedDecomposedParticleWithContact3D;
#else
  using PARTICLETYPE = ResolvedDecomposedParticle3D;
#endif
  T averageSettlingVelocity {0};
  auto& parameters = myCase.getParameters();

  communication::forParticlesInSuperParticleSystem<T, PARTICLETYPE, conditions::valid_particle_centres>(
      xParticleSystem,
      [&](Particle<T, PARTICLETYPE>& particle, ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
        averageSettlingVelocity += getVelocity(particle)[2];
      });

  singleton::mpi().reduceAndBcast(averageSettlingVelocity, MPI_SUM, singleton::mpi().bossId(),
                                  averageParticleVelocityComm);

  return averageSettlingVelocity / parameters.get<parameters::PARTICLE_NUMBER>();
}

template<typename PARTICLESYSTEM>
void updateBodyForce(MyCase& myCase, PARTICLESYSTEM& xParticleSystem)
{
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes {});
  auto& geometry   = myCase.getGeometry();
  auto& converter  = lattice.getUnitConverter();
  SuperLatticePhysExternalPorosity3D<T, DESCRIPTOR> porosity(lattice, converter);
  SuperAverage3D<T>                                 avPorosityF(porosity, geometry, 1);
  int                                               input[1] {};
  T                                                 fluidVolumeFraction[1];
  avPorosityF(fluidVolumeFraction, input);
  const T volumeRatio = parameters.get<parameters::PARTICLE_VOLUME_FRACTION>() / fluidVolumeFraction[0];

  // Apply equal to submerged weight of the particles to the fluid
  std::vector<T> balancingAcceleration(3, T(0.));
  balancingAcceleration[2] = parameters.get<parameters::GRAVITATIONAL_ACC>() * volumeRatio *
                            (parameters.get<parameters::PART_RHO>() / parameters.get<parameters::PHYS_DENSITY>() - 1);

  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(lattice, converter);
  SuperAverage3D<T>                         avgVel(velocity, geometry, 1);
  T                                         vel[3] {};
  avgVel(vel, input);
  // Avoid numerical drift
  balancingAcceleration[2] -= vel[2] * converter.getCharPhysVelocity() / converter.getCharPhysLength();

  const T conversionFactor =
      converter.getConversionFactorTime() * converter.getConversionFactorTime() / converter.getConversionFactorLength();
  balancingAcceleration[2] *= conversionFactor;

  AnalyticalConst3D<T, T> acc(balancingAcceleration);
  fields::set<descriptors::FORCE>(lattice, geometry.getMaterialIndicator({1}), acc);

}

template<typename T>
T calculateCubeEdgeLengthFromSphereRadius(MyCase& myCase)
{
  auto& parameters = myCase.getParameters();
  // Calculate the volume of the sphere
  const T sphereVolume = (4.0 / 3.0) * M_PI * util::pow( parameters.get<parameters::PART_RADIUS>(), 3);

  // Calculate the edge length of the cube
  const T cubeEdgeLength = util::pow(sphereVolume, 1.0 / 3.0);

  return cubeEdgeLength;
}

void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  auto& geometry = myCase.getGeometry();
  auto& cuboidDecomposition = geometry.getCuboidDecomposition();
  clout << "Prepare Geometry ..." << std::endl;

  cuboidDecomposition.setPeriodicity({true, true, true});

  geometry.rename(0, 1);

  geometry.clean();
  geometry.innerClean();

  geometry.checkForErrors();
  geometry.getStatistics().print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes {});
  auto& geometry   = myCase.getGeometry();
  lattice.setUnitConverter<
      UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR>>(
        parameters.get<parameters::RESOLUTION>(),
        parameters.get<parameters::LATTICE_RELAXATION_TIME>(),
        parameters.get<parameters::PHYS_CHAR_LENGTH>(),
        parameters.get<parameters::PHYS_CHAR_VELOCITY>(),
        parameters.get<parameters::PHYS_CHAR_VISCOSITY>(),
        parameters.get<parameters::PHYS_DENSITY>()
  );
  auto& converter  = lattice.getUnitConverter();
  converter.print();

  dynamics::set<PorousParticleKupershtokhForcedBGKdynamics<T, DESCRIPTOR>>(lattice, geometry.getMaterialIndicator({1}));

  lattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  {
    auto& communicator = lattice.getCommunicator(stage::PostPostProcess());
    communicator.requestFields<POROSITY, VELOCITY_NUMERATOR, VELOCITY_DENOMINATOR>();
    communicator.requestOverlap(lattice.getOverlap());
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(MyCase& myCase)
{
  OstreamManager clout(std::cout, "setBoundaryValues");
  auto&                   lattice    = myCase.getLattice(NavierStokes {});

  lattice.initialize();
}

template<typename PARTICLESYSTEM>
void getResults(MyCase& myCase, std::size_t iT, Timer<double>& timer, PARTICLESYSTEM& xParticleSystem)
{
  OstreamManager clout(std::cout, "getResults");
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes {});
  auto& converter  = lattice.getUnitConverter();

  std::size_t iTwrite = converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T>());
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(lattice, converter);
  if (parameters.get<parameters::VTK_ENABLED>()) {
    SuperVTMwriter3D<T>                               vtkWriter("sedimentation");
    SuperLatticePhysPressure3D<T, DESCRIPTOR>         pressure(lattice, converter);
    SuperLatticePhysExternalPorosity3D<T, DESCRIPTOR> externalPor(lattice, converter);

    vtkWriter.addFunctor(pressure);
    vtkWriter.addFunctor(velocity);
    vtkWriter.addFunctor(externalPor);

    if (iT == 0) {
      /// Writes the converter log file
      SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(lattice);
      SuperLatticeRank3D<T, DESCRIPTOR>   rank(lattice);
      vtkWriter.write(cuboid);
      vtkWriter.write(rank);
      vtkWriter.createMasterFile();
    }

    if (iT % iTwrite == 0) {
      vtkWriter.write(iT);
    }
  }

  if (iT % iTwrite == 0) {
    clout << "Average settling velocity: " << evalAverageSettlingVelocity(myCase, xParticleSystem) << " in m/s" << std::endl;

    timer.update(iT);
    timer.printStep();
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
  }
}

void simulate(MyCase& myCase)
{
  using T                             = MyCase::value_t;
  using DESCRIPTOR                    = MyCase::descriptor_t_of<NavierStokes>;
  using namespace olb::particles;
  using namespace olb::particles::dynamics;
  using namespace olb::particles::contact;
  using namespace olb::particles::access;
#ifdef WithContact
  using PARTICLETYPE = ResolvedDecomposedParticleWithContact3D;
#else
  using PARTICLETYPE = ResolvedDecomposedParticle3D;
#endif
  using PARTICLECONTACTTYPE = ParticleContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>;
  using WALLCONTACTTYPE     = WallContactArbitraryFromOverlapVolume<T, PARTICLETYPE::d, true>;
  auto&             parameters        = myCase.getParameters();
  auto&             lattice           = myCase.getLattice(NavierStokes {});
  auto&             geometry          = myCase.getGeometry();
  auto&             converter         = lattice.getUnitConverter();
  OstreamManager clout(std::cout, "simulate");

  std::vector<std::pair<std::string, MyCase::value_t>> limestoneStlFiles = {
    std::make_pair("./limestone/312_reduced.stl", 2.518e-4), std::make_pair("./limestone/529_reduced.stl", 2.625e-4),
    std::make_pair("./limestone/1076_reduced.stl", 1.012e-4), std::make_pair("./limestone/1270_reduced.stl", 1.214e-4),
    std::make_pair("./limestone/1810_reduced.stl", 0.862e-4)};


  std::string positionsfilename = [&parameters] {
    if (parameters.get<parameters::PARTICLE_TYPE>() == 2) {
      return std::string("limestone/particlepositions_limestone_15_0.150021");
    }
    else
      return std::string("particlepositions_sphere_1.5mm_15_0.30");
  }();

  if (MPI_Comm_dup(MPI_COMM_WORLD, &averageParticleVelocityComm) != MPI_SUCCESS) {
    throw std::runtime_error("Unable to duplicate MPI communicator");
  }
  if (MPI_Comm_dup(MPI_COMM_WORLD, &numberParticleCellsComm) != MPI_SUCCESS) {
    throw std::runtime_error("Unable to duplicate MPI communicator");
  }

  std::vector<SolidBoundary<T, DESCRIPTOR::d>> solidBoundaries;
  const T                                      epsilon     = T {0.5} * converter.getConversionFactorLength();
  const T                                      halfEpsilon = T {0.5} * epsilon;
  const T                                      size        = parameters.get<parameters::PHYS_CHAR_LENGTH>() * parameters.get<parameters::EXTENT_TO_DIAMETER>();
  IndicatorCuboid3D<T> cuboid({size, size, size}, {0.,0.,0.});
  constexpr auto       getPeriodicity = []() {
    return Vector<bool, 3>(true, true, true);
  };


  T maxCircumRadius = T {0};

  // Create smooth indicators for limestones
  std::vector<std::shared_ptr<STLreader<T>>>                             limestoneSTLreaders;
  std::vector<std::unique_ptr<olb::SmoothIndicatorCustom3D<T, T, true>>> limestoneIndicators;
  const T latticeSpacingDiscreteParticle = T {0.2} * converter.getConversionFactorLength();

  if (parameters.get<parameters::PARTICLE_TYPE>() == 2) {
    clout << "Initializing limestone ..." << std::endl;

    {
      unsigned i = 0;
      for (auto& STLfile : limestoneStlFiles) {
        limestoneSTLreaders.push_back(
            std::make_shared<STLreader<T>>(STLfile.first, converter.getConversionFactorLength(), STLfile.second));
        limestoneIndicators.push_back(std::make_unique<olb::SmoothIndicatorCustom3D<T, T, true>>(
            latticeSpacingDiscreteParticle, limestoneSTLreaders[i], olb::Vector<T, 3>(T {}), epsilon,
            olb::Vector<T, 3>(T {})));
        limestoneIndicators.back()->calcMofiAndMass(parameters.get<parameters::PART_RHO>());
        ++i;
      }
    }

    clout << "Initializing limestone ... OK" << std::endl;

    for (auto& STLsurface : limestoneIndicators) {
      maxCircumRadius = util::max(maxCircumRadius, STLsurface->getCircumRadius());
    }
  }
  if  (parameters.get<parameters::PARTICLE_TYPE>() == 0) {
    maxCircumRadius = parameters.get<parameters::PART_RADIUS>() + halfEpsilon;
  }

  if constexpr (access::providesContactMaterial<PARTICLETYPE>()) {
    const T detectionDistance = T {0.5} * util::sqrt(PARTICLETYPE::d) * converter.getPhysDeltaX();
    maxCircumRadius           = maxCircumRadius - halfEpsilon + util::max(halfEpsilon, detectionDistance);
  }

  SuperParticleSystem<T, PARTICLETYPE> xParticleSystem(geometry, maxCircumRadius);
  particles::communication::checkCuboidSizes(xParticleSystem);

  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;

  ContactProperties<T, 1> contactProperties;
  contactProperties.set(0, 0,
                        evalEffectiveYoungModulus(parameters.get<parameters::YOUNGS_MODULUS>(), parameters.get<parameters::YOUNGS_MODULUS>(),
                                                  parameters.get<parameters::PART_POISSON_RATIO>(), parameters.get<parameters::PART_POISSON_RATIO>()),
                                                  parameters.get<parameters::COEFFICIENT_OF_RESTITUTION>(),
                                                  parameters.get<parameters::COEFFICIENT_KINETIC_FRICTION>(),
                                                  parameters.get<parameters::COEFFICIENT_STATIC_FRICTION>(),
                                                  parameters.get<parameters::STATIC_KINETIC_TRANSITION_VELOCITY>());

  // Create and assign resolved particle dynamics
  xParticleSystem.defineDynamics<VerletParticleDynamics<T, PARTICLETYPE>>();

  //Create particle manager handling coupling, gravity and particle dynamics
  ParticleManager<T, DESCRIPTOR, PARTICLETYPE> particleManager(xParticleSystem, geometry, lattice, converter,
                                                               parameters.get<parameters::GRAVITY>(), getPeriodicity());

  // Create Communicators

  const communication::ParticleCommunicator& particleCommunicator = particleManager.getParticleCommunicator();

  const std::function<T(const std::size_t&)> getCircumRadius = [&](const std::size_t& pID) {
    if (parameters.get<parameters::PARTICLE_TYPE>() == 2) {
      return limestoneIndicators[pID % limestoneStlFiles.size()]->getCircumRadius();
    }
    if (parameters.get<parameters::PARTICLE_TYPE>() == 0) {
      return parameters.get<parameters::PART_RADIUS>() + T {0.5} * epsilon;
    }
    if (parameters.get<parameters::PARTICLE_TYPE>() == 1) {
      return T {0.5} * (calculateCubeEdgeLengthFromSphereRadius<T>(myCase) * util::sqrt(3) + epsilon);
    }
    throw std::runtime_error("Invalid PARTICLE_TYPE set");
  };
  const std::function<T(const std::size_t&)> getParticleVolume = [&](const std::size_t& pID) {
    if (parameters.get<parameters::PARTICLE_TYPE>() == 2) {
      return limestoneIndicators[pID % limestoneStlFiles.size()]->getVolume();
    }
    if (parameters.get<parameters::PARTICLE_TYPE>() == 0) {
      return T {4} / T {3} * M_PI * parameters.get<parameters::PART_RADIUS>() * parameters.get<parameters::PART_RADIUS>() * parameters.get<parameters::PART_RADIUS>();
    }
    if (parameters.get<parameters::PARTICLE_TYPE>() == 1) {
      return util::pow(calculateCubeEdgeLengthFromSphereRadius<T>(myCase), 3);
    }
    throw std::runtime_error("Invalid PARTICLE_TYPE set");
  };

  const std::function<void(const particles::creators::SpawnData<T, DESCRIPTOR::d>&, const std::string&)>
      createParticleFromString =
          [&](const particles::creators::SpawnData<T, DESCRIPTOR::d>& data, const std::string& pID) {
            const PhysR<T, 3>  physPosition   = data.position;
            const Vector<T, 3> angleInDegrees = data.angleInDegree;

            if (parameters.get<parameters::PARTICLE_TYPE>() == 2) {
              if (parameters.get<parameters::PARTICLE_NUMBER>() < limestoneStlFiles.size()) {
                std::shared_ptr<STLreader<T>> limestoneIndicator = std::make_shared<STLreader<T>>(
                    limestoneStlFiles[parameters.get<parameters::PARTICLE_NUMBER>()].first, converter.getConversionFactorLength(),
                    limestoneStlFiles[parameters.get<parameters::PARTICLE_NUMBER>()].second);
                creators::addResolvedArbitraryShape3D<T, PARTICLETYPE>(xParticleSystem, physPosition,
                                                                       latticeSpacingDiscreteParticle,
                                                                       limestoneIndicator, epsilon,
                                                                       parameters.get<parameters::PART_RHO>());
              }
              else {
                creators::addResolvedObject<T, PARTICLETYPE>(xParticleSystem, parameters.get<parameters::PARTICLE_NUMBER>() % limestoneStlFiles.size(),
                                                             physPosition, parameters.get<parameters::PART_RHO>(), angleInDegrees);
              }
            }
            if (parameters.get<parameters::PARTICLE_TYPE>() == 0) {
              creators::addResolvedSphere3D(xParticleSystem, physPosition, parameters.get<parameters::PART_RADIUS>(), epsilon, parameters.get<parameters::PART_RHO>());
            }
            if (parameters.get<parameters::PARTICLE_TYPE>() == 1) {
              creators::addResolvedCuboid3D(xParticleSystem, physPosition,
                                            Vector<T, 3>(calculateCubeEdgeLengthFromSphereRadius<T>(myCase)), epsilon,
                                            parameters.get<parameters::PART_RHO>());
            }
            parameters.set<parameters::PARTICLE_NUMBER>(parameters.get<parameters::PARTICLE_NUMBER>() + 1);
          };

  if (positionsfilename != "") {
    clout << "Spawning particles from " << positionsfilename << " ..." << std::endl;
    if (std::filesystem::exists(positionsfilename)) {
      std::vector<particles::creators::SpawnData<T, DESCRIPTOR::d>> tmpSpawnData;
      tmpSpawnData = particles::creators::setParticles<T, 3>(positionsfilename, parameters.get<parameters::WANTED_PARTICLE_VOLUME_FRACTION>(),
                                                            cuboid, size*size*size, getParticleVolume, createParticleFromString);

      T tmpVolume = T {0};
      for (unsigned pID = 0; pID < tmpSpawnData.size(); ++pID) {
        tmpVolume += getParticleVolume(pID);
      }
      parameters.set<parameters::PARTICLE_VOLUME_FRACTION>(tmpVolume / (size*size*size));
    }
    else {
      OLB_ASSERT(false, positionsfilename + " does not exist.");
    }
  }
  else {
    OLB_ASSERT(false, "No particle positions file given.");
  }
  if (parameters.get<parameters::PARTICLE_TYPE>() == 2) {
    limestoneIndicators.clear();
  }

  clout << "Spawning particles from " << positionsfilename << " ... OK" << std::endl;

  maxCircumRadius = 0.;
  communication::forParticlesInSuperParticleSystem(xParticleSystem, [&](Particle<T, PARTICLETYPE>&       particle,
                                                         ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
#ifdef WithContact
    particle.setField<MECHPROPERTIES, MATERIAL>(0);
    particle.setField<NUMERICPROPERTIES, ENLARGEMENT_FOR_CONTACT>(
        parameters.get<parameters::PARTICLE_ENLARGEMENT_FOR_CONTACT>());
#endif // WithContact

    const T currCircumRadius = access::getRadius(particle);
    maxCircumRadius          = util::max(currCircumRadius, maxCircumRadius);
  });

  singleton::mpi().reduceAndBcast(maxCircumRadius, MPI_MAX, singleton::mpi().bossId(),
                                  particleCommunicator.contactTreatmentComm);

  xParticleSystem.updateOffsetFromCircumRadius(maxCircumRadius);

  setBoundaryValues(myCase);

  clout << "MaxIT: " << converter.getLatticeTime(parameters.get<parameters::PHYS_MAX_T>()) << std::endl;
  std::size_t iTpurge = converter.getLatticeTime(parameters.get<parameters::PHYS_PURGE_ITER_T>());

  Timer<T> timer(converter.getLatticeTime(parameters.get<parameters::PHYS_MAX_T>()), geometry.getStatistics().getNvoxel());
  timer.start();
  for (std::size_t iT = 0; iT < converter.getLatticeTime(parameters.get<parameters::PHYS_MAX_T>()); ++iT) {
#ifndef WithContact
    particleManager.execute<
        couple_lattice_to_parallel_particles<T, DESCRIPTOR, PARTICLETYPE>,
        communicate_parallel_surface_force<T, PARTICLETYPE>,
        apply_gravity<T, PARTICLETYPE>,
        process_dynamics_parallel<T, PARTICLETYPE>,
        update_particle_core_distribution<T, PARTICLETYPE>>();

    particles::dynamics::coupleResolvedParticlesToLattice<
        T, DESCRIPTOR, PARTICLETYPE, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
        xParticleSystem, contactContainer, geometry, lattice, converter,
        solidBoundaries, getPeriodicity);

#else  // WithContact
    // Couple lattice to particles
    couple_lattice_to_parallel_particles<T, DESCRIPTOR, PARTICLETYPE>::execute(
        xParticleSystem, geometry, lattice, converter, getPeriodicity());

    communicate_parallel_surface_force<T, PARTICLETYPE>::execute(
        xParticleSystem, particleCommunicator);

    // Apply contacts
    particles::contact::processContacts<T, PARTICLETYPE, PARTICLECONTACTTYPE,
                                        WALLCONTACTTYPE,
                                        ContactProperties<T, 1>>(
        xParticleSystem, solidBoundaries, contactContainer, contactProperties,
        geometry,
        particleCommunicator.contactTreatmentComm,
        parameters.get<parameters::CONTACT_BOX_RESOLUTION_PER_DIRECTION>(), T {4. / (3. * util::sqrt(M_PI))},
        getPeriodicity);

    // Apply gravity
    communication::forParticlesInSuperParticleSystem<
        T, PARTICLETYPE,
        conditions::valid_particle_centres //only consider center for resolved
        >(xParticleSystem,
          [&](Particle<T, PARTICLETYPE>&       particle,
              ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
            apply_gravity<T, PARTICLETYPE>::execute(xParticleSystem, particle,
                                                    parameters.get<parameters::GRAVITY>(),
                                                    converter.getPhysDeltaT());
          });


    // Process particles (Verlet algorithm)
    communication::forParticlesInSuperParticleSystem<
        T, PARTICLETYPE,
        conditions::valid_particle_centres //only consider center for resolved
        >(xParticleSystem,
          [&](Particle<T, PARTICLETYPE>&       particle,
              ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
            particle.process(converter.getPhysDeltaT());
          });

    communicatePostContactTreatmentContacts(
        contactContainer, xParticleSystem, converter.getPhysDeltaX(),
        particleCommunicator.particleContactDetectionComm,
        particleCommunicator.wallContactDetectionComm,
        getPeriodicity());

    update_particle_core_distribution<T, PARTICLETYPE>::execute(
        xParticleSystem, converter.getPhysDeltaX(), particleCommunicator,
        getPeriodicity());

    // Couple particles to lattice
    particles::dynamics::coupleResolvedParticlesToLattice<
        T, DESCRIPTOR, PARTICLETYPE, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
        xParticleSystem, contactContainer, geometry, lattice, converter,
        solidBoundaries, getPeriodicity);

    // Communicate found contacts
    communicateContactsParallel(
        contactContainer, xParticleSystem, converter.getPhysDeltaX(),
        particleCommunicator.particleContactDetectionComm,
        particleCommunicator.wallContactDetectionComm,
        getPeriodicity());

    if constexpr (isPeriodic(getPeriodicity())) {
      accountForPeriodicParticleBoundary(xParticleSystem, contactContainer,
                                         geometry, getPeriodicity);
    }


#endif // WithContact

    if (iT % iTpurge == 0) {
      purgeInvalidParticles<T, PARTICLETYPE>(xParticleSystem);
    }

    getResults(myCase, iT, timer, xParticleSystem);

    updateBodyForce(myCase, xParticleSystem);

    lattice.collideAndStream();

    timer.update(iT);
  }
  timer.stop();
  timer.printSummary();
}
#endif                                // PARALLEL_MODE_MPI

int main(int argc, char* argv[])
{
  /// === 1st Step: Initialization ===
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");
#ifdef PARALLEL_MODE_MPI //Check if MPI is available, otherwise throw error
  std::vector<std::string> cmdInput;
  MyCase::ParametersD      myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<PARTICLE_NUMBER>(0);
    myCaseParameters.set<PARTICLE_TYPE>(0); // 0: spheres, 1: cubes, 2: limestones
    myCaseParameters.set<WANTED_PARTICLE_VOLUME_FRACTION>(0.15);
    myCaseParameters.set<RESOLUTION>(9);
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.55);
    myCaseParameters.set<ARCHIMEDES>(1000.);
    myCaseParameters.set<GRAVITATIONAL_ACC>(9.81);
    myCaseParameters.set<DENSITY_RATIO>(3.3);
    myCaseParameters.set<PART_RADIUS>(0.0015);
    myCaseParameters.set<COEFFICIENT_OF_RESTITUTION>(0.926);
    myCaseParameters.set<COEFFICIENT_KINETIC_FRICTION>(0.16);
    myCaseParameters.set<COEFFICIENT_STATIC_FRICTION>(0.32);
    myCaseParameters.set<STATIC_KINETIC_TRANSITION_VELOCITY>(0.001);
    myCaseParameters.set<YOUNGS_MODULUS>(5.0e3);
    myCaseParameters.set<PART_POISSON_RATIO>(0.245);
    myCaseParameters.set<CONTACT_BOX_RESOLUTION_PER_DIRECTION>(8);
    myCaseParameters.set<EXTENT_TO_DIAMETER>(15.0);
    myCaseParameters.set<PHYS_DENSITY>(1000);
    myCaseParameters.set<VTK_ENABLED>(true);
    myCaseParameters.set<PHYS_CHAR_LENGTH>([] (MyCase::ParametersD& p) {
      return p.get<PART_RADIUS>() * 2.0;
    });
    myCaseParameters.set<PART_RHO>([] (MyCase::ParametersD& p) {
      return p.get<DENSITY_RATIO>() * p.get<PHYS_DENSITY>();
    });
    myCaseParameters.set<PHYS_DELTA_X>([] (MyCase::ParametersD& p) {
      return p.get<PART_RADIUS>() * 2.0 / p.get<RESOLUTION>();
    });
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>([] (MyCase::ParametersD& p) {
      return util::sqrt(p.get<GRAVITATIONAL_ACC>() * (p.get<DENSITY_RATIO>() - 1) *
                        util::pow(p.get<PART_RADIUS>() * 2.0, 3) /
                        p.get<ARCHIMEDES>());
    });
    myCaseParameters.set<GRAVITY>([] (MyCase::ParametersD& p) {
      return Vector<double, 3> {0., 0.,
                                -p.get<GRAVITATIONAL_ACC>() *
                                    (1. - p.get<PHYS_DENSITY>() / p.get<PART_RHO>())};
    });
    myCaseParameters.set<PHYS_DYNAMIC_VISCOSITY>([] (MyCase::ParametersD& p) {
      return p.get<PHYS_CHAR_VISCOSITY>() * p.get<PHYS_DENSITY>();
    });
    myCaseParameters.set<PHYS_CHAR_VELOCITY>([] (MyCase::ParametersD& p) {
      return (2. / 9.) * (p.get<PART_RHO>() - p.get<PHYS_DENSITY>()) *
             p.get<GRAVITATIONAL_ACC>() * p.get<PART_RADIUS>() *
             p.get<PART_RADIUS>() / p.get<PHYS_DYNAMIC_VISCOSITY>();
    });
    myCaseParameters.set<PHYS_MAX_T>([] (MyCase::ParametersD& p) {
      return 200.0 * 2. * p.get<PART_RADIUS>() / p.get<PHYS_CHAR_VELOCITY>();
    });
    myCaseParameters.set<PHYS_VTK_ITER_T>([] (MyCase::ParametersD& p) {
      return 0.02 * p.get<PHYS_MAX_T>();
    });
    myCaseParameters.set<PHYS_PURGE_ITER_T>([] (MyCase::ParametersD& p) {
      return 0.06 * p.get<PHYS_MAX_T>();
    });
    myCaseParameters.set<PARTICLE_ENLARGEMENT_FOR_CONTACT>([] (MyCase::ParametersD& p) {
      return p.get<PHYS_DELTA_X>() / 5.0;
    });
  }

  Mesh   mesh = createMesh(myCaseParameters);
  MyCase myCase(myCaseParameters, mesh);

  prepareGeometry(myCase);

  prepareLattice(myCase);

  simulate(myCase);

#else
  throw std::runtime_error("This example is not designed to run in serial mode.");
#endif // PARALLEL_MODE_MPI
}
