/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006-2021 Nicolas Hafen, Mathias J. Krause
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

/* settlingCube3d.cpp:
 * The case examines the settling of a cubical silica particle
 * under the influence of gravity.
 * The object is surrounded by water in a rectangular domain
 * limited by no-slip boundary conditions.
 * For the calculation of forces an DNS approach is chosen
 * which also leads to a back-coupling of the particle on the fluid,
 * inducing a flow.
 *
 * The simulation is based on the homogenised lattice Boltzmann approach
 * (HLBM) introduced in "Particle flow simulations with homogenised
 * lattice Boltzmann methods" by Krause et al.
 * and extended in "Towards the simulation of arbitrarily shaped 3D particles
 * using a homogenised lattice Boltzmann method" by Trunk et al.
 * for the simulation of 3D particles.
 *
 * This example demonstrates the usage of HLBM in the OpenLB framework.
 * To improve parallel performance, the particle decomposition scheme
 * described in 10.48550/arXiv.2312.14172 is used.
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::graphics;
using namespace olb::util;

using MyCase = Case<
NavierStokes, Lattice<double, descriptors::PorousParticleD3Q19Descriptor>
>;

namespace olb::parameters {
struct CUBE_DENSITY : public descriptors::FIELD_BASE<1> {};
struct CUBE_1_EDGE_LENGTH : public descriptors::FIELD_BASE<1> {};
struct CUBE_2_EDGE_LENGTH : public descriptors::FIELD_BASE<1> {};
struct CUBE_1_ORIENTATION : public descriptors::FIELD_BASE<0,1,0> {};
struct CUBE_1_POSITION : public descriptors::FIELD_BASE<0,1,0> {};
struct CUBE_2_ORIENTATION : public descriptors::FIELD_BASE<0, 1, 0> {};
struct CUBE_2_POSITION : public descriptors::FIELD_BASE<0,1,0> {};
struct EXTERNAL_ACCELERATION : public descriptors::FIELD_BASE<0, 1, 0> {};
struct CUBE_1_VELOCITY : public descriptors::FIELD_BASE<0, 1, 0> {};
struct CUBE_2_VELOCITY : public descriptors::FIELD_BASE<0, 1, 0> {};
}

Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters)
{
  using T = MyCase::value_t;

  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector origin {0., 0., 0.};

  IndicatorCuboid3D<T> cuboid(extent, origin);

  const T physDeltaX = parameters.get<parameters::DOMAIN_EXTENT>()[0] / parameters.get<parameters::RESOLUTION>();

  #ifdef PARALLEL_MODE_MPI
    Mesh<T, MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  #else
    Mesh<T, MyCase::d> mesh(cuboid, physDeltaX, 7);
  #endif

  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  auto& geometry = myCase.getGeometry();

  geometry.rename(0, 2);
  geometry.rename(2, 1, {1, 1, 1});

  geometry.clean();
  geometry.innerClean();

  geometry.checkForErrors();
  geometry.getStatistics().print();

  clout << "Prepare Geometry ... OK" << std::endl;
  return;
}

void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;
  clout << "setting Velocity Boundaries ..." << std::endl;

  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

  auto& params   = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice  = myCase.getLattice(NavierStokes {});

  const int res                 = params.get<parameters::RESOLUTION>();
  const T   charLatticeVelocity = params.get<parameters::LATTICE_CHAR_VELOCITY>();
  const T   charPhysLength      = params.get<parameters::PHYS_CHAR_LENGTH>();
  const T   charPhysVelocity    = params.get<parameters::PHYS_CHAR_VELOCITY>();
  const T   physViscosity       = params.get<parameters::PHYS_KINEMATIC_VISCOSITY>();
  const T   physDensity         = params.get<parameters::PHYS_DENSITY>();

  lattice.setUnitConverter<UnitConverterFromResolutionAndLatticeVelocity<T, DESCRIPTOR>>(
      (int)res,
      (T)charLatticeVelocity,
      (T)charPhysLength,
      (T)charPhysVelocity,
      (T)physViscosity,
      (T)physDensity
  );

  auto& converter = lattice.getUnitConverter();

  dynamics::set<PorousParticleBGKdynamics>(lattice, geometry.getMaterialIndicator({1}));
  boundary::set<boundary::BounceBack>(lattice, geometry, 2);

  lattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  {
    auto& communicator = lattice.getCommunicator(stage::PostPostProcess());
    communicator
        .requestFields<descriptors::POROSITY, descriptors::VELOCITY_NUMERATOR, descriptors::VELOCITY_DENOMINATOR>();
    communicator.requestOverlap(lattice.getOverlap());
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase)
{
  OstreamManager clout(std::cout, "setBoundaryValues");

  auto& geometry = myCase.getGeometry();
  auto& lattice  = myCase.getLattice(NavierStokes {});

  fields::set<descriptors::POROSITY>(lattice, geometry.getMaterialIndicator({0,1,2}), 1.);

  lattice.initialize();
}

template <typename PARTICLESYSTEM>
void getResults(MyCase& myCase, int iT, Timer<MyCase::value_t>& timer, PARTICLESYSTEM& xParticleSystem)
{
  #ifdef PARALLEL_MODE_MPI
    using PARTICLETYPE = descriptors::ResolvedDecomposedParticle3D;
  #endif
  using namespace olb::particles;


  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

  auto& params    = myCase.getParameters();
  auto& lattice   = myCase.getLattice(NavierStokes {});
  auto& converter = lattice.getUnitConverter();

  const T iTwrite = params.get<parameters::PHYS_VTK_ITER_T>();

  if(params.get<parameters::VTK_ENABLED>()){
    SuperVTMwriter3D<T>                               vtkWriter("sedimentation");
    SuperLatticePhysVelocity3D<T, DESCRIPTOR>         velocity(lattice, converter);
    SuperLatticePhysPressure3D<T, DESCRIPTOR>         pressure(lattice, converter);
    SuperLatticePhysExternalPorosity3D<T, DESCRIPTOR> externalPor(lattice, converter);
    vtkWriter.addFunctor(velocity);
    vtkWriter.addFunctor(pressure);
    vtkWriter.addFunctor(externalPor);

    if (iT == 0) {
      SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(lattice);
      SuperLatticeRank3D<T, DESCRIPTOR>   rank(lattice);
      vtkWriter.write(cuboid);
      vtkWriter.write(rank);
      vtkWriter.createMasterFile();
    }

    if (iT % converter.getLatticeTime(iTwrite) == 0) {
      vtkWriter.write(iT);
    }
  }

  if (iT % converter.getLatticeTime(iTwrite) == 0) {
    timer.update(iT);
    timer.printStep();
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
    #ifdef PARALLEL_MODE_MPI
      communication::forParticlesInSuperParticleSystem<T, PARTICLETYPE, conditions::valid_particle_centres>(
          xParticleSystem,
          [&](Particle<T, PARTICLETYPE>& particle, ParticleSystem<T, PARTICLETYPE>& particleSystem, int globiC) {
            particles::io::printResolvedParticleInfo(particle);
          });
    #else
    for (std::size_t iP = 0; iP < xParticleSystem.size(); ++iP) {
      auto particle = xParticleSystem.get(iP);
      particles::io::printResolvedParticleInfo(particle);
    }
    #endif
  }
}

void simulate(MyCase& myCase)
{
  OstreamManager clout(std::cout, "getResults");

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using namespace olb::particles;
  using namespace olb::particles::dynamics;

  auto& params    = myCase.getParameters();
  auto& geometry  = myCase.getGeometry();
  auto& lattice   = myCase.getLattice(NavierStokes {});
  auto& converter = lattice.getUnitConverter();

  const Vector externalAcceleration  = params.get<parameters::EXTERNAL_ACCELERATION>();
  const T maxPhysT              = params.get<parameters::PHYS_MAX_T>();
  const Vector cube1Center      = params.get<parameters::CUBE_1_POSITION>();
  const Vector cube2Center      = params.get<parameters::CUBE_2_POSITION>();
  const T      cube1EdgeLength  = params.get<parameters::CUBE_1_EDGE_LENGTH>();
  const T      cube2EdgeLength  = params.get<parameters::CUBE_2_EDGE_LENGTH>();
  const T      cubeDensity      = params.get<parameters::CUBE_DENSITY>();
  const Vector cube1Orientation = params.get<parameters::CUBE_1_ORIENTATION>();
  const Vector cube2Orientation = params.get<parameters::CUBE_2_ORIENTATION>();

  #ifdef PARALLEL_MODE_MPI
    using PARTICLETYPE = descriptors::ResolvedDecomposedParticle3D;
  #else
    using PARTICLETYPE = descriptors::ResolvedParticle3D;
  #endif

  #ifdef PARALLEL_MODE_MPI
    SuperParticleSystem<T, PARTICLETYPE> particleSystem(geometry);
  #else
    ParticleSystem<T, PARTICLETYPE> particleSystem;
  #endif

  ParticleManager<T, DESCRIPTOR, PARTICLETYPE> particleManager(particleSystem, geometry, lattice, converter,
                                                               externalAcceleration);

  particleSystem.defineDynamics<VerletParticleDynamics<T, PARTICLETYPE>>();

  T            epsilon = 0.5 * converter.getPhysDeltaX();
  Vector<T, 3> cube1Extend(cube1EdgeLength);
  Vector<T, 3> cube2Extend(cube2EdgeLength);

  creators::addResolvedCuboid3D<T, PARTICLETYPE>(particleSystem, cube1Center, cube1Extend, epsilon, cubeDensity, cube1Orientation);

  creators::addResolvedCuboid3D<T, PARTICLETYPE>(particleSystem, cube2Center, cube2Extend, epsilon, cubeDensity, cube2Orientation);

  particleSystem.checkForErrors();

  Timer<T> timer(converter.getLatticeTime(maxPhysT), geometry.getStatistics().getNvoxel());
  timer.start();

  clout << "MaxIT: " << converter.getLatticeTime(maxPhysT) << std::endl;

  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT) + 10; ++iT) {

    particleManager.execute<couple_lattice_to_particles<T, DESCRIPTOR, PARTICLETYPE>,

    #ifdef PARALLEL_MODE_MPI
      communicate_surface_force<T, PARTICLETYPE>,
    #endif
      apply_gravity<T, PARTICLETYPE>, process_dynamics<T, PARTICLETYPE>,
    #ifdef PARALLEL_MODE_MPI
      update_particle_core_distribution<T, PARTICLETYPE>,
    #endif
    couple_particles_to_lattice<T, DESCRIPTOR, PARTICLETYPE>>();

    getResults(myCase, iT, timer, particleSystem);

    lattice.collideAndStream();
  }

  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");

  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    using T = MyCase::value_t;
    myCaseParameters.set<RESOLUTION>(30);
    myCaseParameters.set<LATTICE_CHAR_VELOCITY>(0.01);
    myCaseParameters.set<PHYS_MAX_T>(0.5);
    myCaseParameters.set<PHYS_VTK_ITER_T>(0.02);
    myCaseParameters.set<VTK_ENABLED>(true);
    myCaseParameters.set<DOMAIN_EXTENT>({0.01, 0.01, 0.05});
    myCaseParameters.set<PHYS_DENSITY>(1000);
    myCaseParameters.set<PHYS_KINEMATIC_VISCOSITY>(1.e-5);
    myCaseParameters.set<CUBE_DENSITY>(2500);
    myCaseParameters.set<CUBE_1_EDGE_LENGTH>(0.0025);
    myCaseParameters.set<CUBE_2_EDGE_LENGTH>(0.0025);
    myCaseParameters.set<CUBE_1_ORIENTATION>({0., 15., 0.});
    myCaseParameters.set<CUBE_2_ORIENTATION>({0., 0., 15.});
    myCaseParameters.set<CUBE_1_VELOCITY>({0., 0., 0.});
    myCaseParameters.set<CUBE_2_VELOCITY>({0., 0., 0.});
    myCaseParameters.set<PHYS_CHAR_LENGTH>(myCaseParameters.get<DOMAIN_EXTENT>()[0]);
    myCaseParameters.set<PHYS_CHAR_VELOCITY>(0.15);
    myCaseParameters.set<EXTERNAL_ACCELERATION>([] (MyCase::ParametersD& p) -> Vector<T,3> {
      return {.0,.0, -T(9.81) * (T(1) - p.get<PHYS_DENSITY>() / p.get<CUBE_DENSITY>())};
    });
    myCaseParameters.set<CUBE_1_POSITION>([] (MyCase::ParametersD& p) -> Vector<T,3> {
      return {p.get<DOMAIN_EXTENT>()[0] * (T).5, p.get<DOMAIN_EXTENT>()[1] * (T).5,
              p.get<DOMAIN_EXTENT>()[2] * (T).9};
    });
    myCaseParameters.set<CUBE_2_POSITION>([] (MyCase::ParametersD& p) -> Vector<T,3> {
      return {p.get<DOMAIN_EXTENT>()[0] * (T).5, p.get<DOMAIN_EXTENT>()[1] * (T).51,
              p.get<DOMAIN_EXTENT>()[2] * (T).7};
    });
  }

  Mesh mesh = createMesh(myCaseParameters);

  MyCase myCase(myCaseParameters, mesh);

  prepareGeometry(myCase);

  prepareLattice(myCase);

  setInitialValues(myCase);

  simulate(myCase);
}
