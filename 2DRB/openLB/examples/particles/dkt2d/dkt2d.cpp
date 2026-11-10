/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006-2023 Nicolas Hafen, Robin Trunk, Jan E. Marquardt, Mathias J. Krause
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



/* dkt2d.cpp:
 * The case examines the settling of two circles under gravity
 * in a surrounding fluid. The rectangular domain is limited
 * by no-slip boundary conditions.
 * For the calculation of forces a DNS approach is chosen
 * which also leads to a back-coupling of the particle on the fluid,
 * inducing a flow.
 * The simulation is based on the homogenised lattice Boltzmann approach
 * (HLBM) introduced by Krause et al. in "Particle flow simulations
 * with homogenised lattice Boltzmann methods".
 * The contact treatment follows 10.1016/j.partic.2022.12.005.
 * The drafting-kissing-tumbling benchmark case is e.g. described
 * in "Drafting, kissing and tumbling process of two particles
 * with different sizes" by Wang et al.
 * or "The immersed boundary-lattice Boltzmann method
 * for solving fluid-particles interaction problems" by Feng and Michaelides.
 * The example demonstrates the usage of HLBM in the OpenLB framework
 * as well as the utilisation of the Gnuplot-writer
 * to print simulation results.
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace olb::names;


// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::PorousParticleWithContactD2Q9Descriptor>
>;

namespace olb::parameters {

struct COEFFICIENT_OF_RESTITUTION : public descriptors::FIELD_BASE<1> { };
struct COEFFICIENT_STATIC_FRICTION : public descriptors::FIELD_BASE<1> { };
struct COEFFICIENT_KINETIC_FRICTION : public descriptors::FIELD_BASE<1> { };
struct CENTER_1 : public descriptors::FIELD_BASE<0,1> { };
struct CENTER_2 : public descriptors::FIELD_BASE<0,1> { };
struct CONTACT_BOX_RESOLUTION_PER_DIRECTION : public descriptors::FIELD_BASE<1> { };
struct PARTICLE_CONTACT_MATERIAL : public descriptors::FIELD_BASE<1> { };
struct WALL_CONTACT_MATERIAL : public descriptors::FIELD_BASE<1> { };
struct EPS : public descriptors::FIELD_BASE<1> { };
struct PART_POISSON_RATIO : public descriptors::FIELD_BASE<1> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector origin{0, 0};
  IndicatorCuboid2D<T> cuboid(extent, origin);

  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;
  auto& geometry = myCase.getGeometry();

  geometry.rename(0, 2);
  geometry.rename(2, 1, {1, 1});

  geometry.clean();
  geometry.innerClean();

  geometry.checkForErrors();
  geometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  lattice.setUnitConverter<UnitConverter<T,DESCRIPTOR>>(
    ( T )   parameters.get<parameters::PHYS_DELTA_X>(),
    ( T )   parameters.get<parameters::PHYS_DELTA_T>(),
    ( T )   parameters.get<parameters::PHYS_CHAR_LENGTH>(),
    ( T )   parameters.get<parameters::PHYS_CHAR_VELOCITY>(),
    ( T )   parameters.get<parameters::PHYS_CHAR_VISCOSITY>(),
    ( T )   parameters.get<parameters::PHYS_CHAR_DENSITY>()
  );
  auto& converter = lattice.getUnitConverter();
  converter.print();
  clout << "Prepare Lattice ..." << std::endl;

  dynamics::set<PorousParticleBGKdynamics>(lattice, geometry.getMaterialIndicator({1}));
  boundary::set<boundary::BounceBack>(lattice, geometry, 2);

  lattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(MyCase& myCase)
{
  OstreamManager clout(std::cout, "setBoundaryValues");
  clout << "Set Boundary Values ..." << std::endl;

  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();

  fields::set<POROSITY>(lattice, geometry.getMaterialIndicator({1,2}), 1);

  lattice.initialize();
  clout << "Set Boundary Values ... OK" << std::endl;
}

template<typename PARTICLESYSTEM>
void getResults(MyCase& myCase, int iT,
                Timer<MyCase::value_t>& timer,
                PARTICLESYSTEM& particleSystem )
{
  OstreamManager clout(std::cout, "getResults");
  using T = MyCase::value_t;
  using namespace olb::particles;
  using namespace olb::particles::access;

  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& converter = lattice.getUnitConverter();
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using PARTICLETYPE = ResolvedCircleWithContact2D;

  if(parameters.get<parameters::VTK_ENABLED>()){
    SuperVTMwriter2D<T> vtkWriter("sedimentation");
    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(lattice, converter);
    SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure(lattice, converter);
    SuperLatticePhysExternalPorosity2D<T, DESCRIPTOR> externalPor(lattice, converter);
    SuperLatticeMomentumExchangeForceLocal<T, DESCRIPTOR, PARTICLETYPE> momentumExchange(
      lattice, converter, geometry, particleSystem);

    vtkWriter.addFunctor(velocity);
    vtkWriter.addFunctor(pressure);
    vtkWriter.addFunctor(externalPor);
    vtkWriter.addFunctor(momentumExchange);

    if (iT == 0) {
      converter.write("dkt");
      SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid(lattice);
      SuperLatticeRank2D<T, DESCRIPTOR> rank(lattice);
      vtkWriter.write(cuboid);
      vtkWriter.write(rank);
      vtkWriter.createMasterFile();
    }

    if (iT % converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T>()) == 0) {
      vtkWriter.write(iT);
    }
  }


  auto particleA = particleSystem.get( 0 );
  auto particleB = particleSystem.get( 1 );

  if (parameters.get<parameters::GNUPLOT_ENABLED>()) {
    std::string gnuplotFilename = "gnuplot.dat";
    if (iT % converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T>()) == 0) {
      if (singleton::mpi().getRank() == 0) {

        std::ofstream myfile;
        myfile.open (gnuplotFilename.c_str(), std::ios::app);
        T p2PosY = particleB.template getField<GENERAL,POSITION>()[1];
        T p1PosY = particleA.template getField<GENERAL,POSITION>()[1];
        T p2PosX = particleB.template getField<GENERAL,POSITION>()[0];
        T p1PosX = particleA.template getField<GENERAL,POSITION>()[0];
        myfile
            << converter.getPhysTime(iT) << " "
            << std::setprecision(9)
            << p2PosY << " "
            << p1PosY << " "
            << p2PosX << " "
            << p1PosX << std::endl;
        myfile.close();
      }
    }
  }

  if (iT % converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T>()) == 0) {
    timer.update(iT);
    timer.printStep();
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
    for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
      auto particle = particleSystem.get(iP);
      particles::io::printResolvedParticleInfo(particle);
    }
  }

  return;
}

void simulate(MyCase& myCase )
{
  using T = MyCase::value_t;
  using namespace olb::particles;
  using namespace olb::particles::access;
  using namespace olb::particles::contact;
  using namespace olb::particles::dynamics;
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& converter = lattice.getUnitConverter();
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using PARTICLETYPE = ResolvedCircleWithContact2D;
  using PARTICLECONTACTTYPE = ParticleContactArbitraryFromOverlapVolume<T, DESCRIPTOR::d, true>;
  using WALLCONTACTTYPE = WallContactArbitraryFromOverlapVolume<T, DESCRIPTOR::d, true>;
  Timer<T> timer(converter.getLatticeTime(parameters.get<parameters::PHYS_MAX_T>()), geometry.getStatistics().getNvoxel());
  timer.start();
  OstreamManager clout(std::cout, "simulate");

  ParticleSystem<T,PARTICLETYPE> particleSystem;

  ParticleManager<T,DESCRIPTOR,PARTICLETYPE> particleManager(
    particleSystem, geometry, lattice, converter, parameters.get<parameters::GRAVITY>());

  particleSystem.defineDynamics<VerletParticleDynamics<T,PARTICLETYPE>>();

  std::vector<SolidBoundary<T, DESCRIPTOR::d>> solidBoundaries;
  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  IndicatorCuboid2D<T> cuboid(extent, {0, 0});

  solidBoundaries.push_back(  SolidBoundary<T, DESCRIPTOR::d>(
        std::make_unique<IndicInverse<T, DESCRIPTOR::d>>(
          cuboid, cuboid.getMin() - 5 * converter.getPhysDeltaX(),
          cuboid.getMax() + 5 * converter.getPhysDeltaX()), 2, parameters.get<parameters::WALL_CONTACT_MATERIAL>()));

  // Create objects for contact treatment
  ContactContainer<T, PARTICLECONTACTTYPE, WALLCONTACTTYPE> contactContainer;
  // Generate lookup table for contact properties
  ContactProperties<T, 1> contactProperties;
  contactProperties.set(parameters.get<parameters::PARTICLE_CONTACT_MATERIAL>(), parameters.get<parameters::WALL_CONTACT_MATERIAL>(),
                        evalEffectiveYoungModulus(parameters.get<parameters::YOUNGS_MODULUS>(), parameters.get<parameters::YOUNGS_MODULUS>(),
                                                  parameters.get<parameters::POISSON_RATIO>(), parameters.get<parameters::POISSON_RATIO>()),
                                                  parameters.get<parameters::COEFFICIENT_OF_RESTITUTION>(), parameters.get<parameters::COEFFICIENT_KINETIC_FRICTION>(),
                                                  parameters.get<parameters::COEFFICIENT_STATIC_FRICTION>());


  T epsilon = parameters.get<parameters::EPS>() * converter.getPhysDeltaX();

  creators::addResolvedCircle2D( particleSystem, parameters.get<parameters::CENTER_1>(),
                                 parameters.get<parameters::PART_RADIUS>(), epsilon, parameters.get<parameters::PART_RHO>() );

  creators::addResolvedCircle2D( particleSystem, parameters.get<parameters::CENTER_2>(),
                                 parameters.get<parameters::PART_RADIUS>(), epsilon, parameters.get<parameters::PART_RHO>() );

  particleSystem.checkForErrors();

  for (std::size_t iP = 0; iP < particleSystem.size(); ++iP) {
    auto particle = particleSystem.get(iP);
    setContactMaterial(particle, parameters.get<parameters::PARTICLE_CONTACT_MATERIAL>());
  }

  setBoundaryValues(myCase);

  {
    auto& communicator = lattice.getCommunicator(stage::PostPostProcess());
    communicator.requestOverlap(lattice.getOverlap());
    communicator.requestFields<POROSITY,VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR>();
    communicator.exchangeRequests();
  }

  clout << "MaxIT: " << converter.getLatticeTime( parameters.get<parameters::PHYS_MAX_T>() ) << std::endl;
  for (std::size_t iT = 0; iT < converter.getLatticeTime( parameters.get<parameters::PHYS_MAX_T>() ) + 10; ++iT) {

    particleManager.execute<
      couple_lattice_to_particles<T,DESCRIPTOR,PARTICLETYPE>,
      apply_gravity<T,PARTICLETYPE>
    >();

    processContacts<T, PARTICLETYPE, PARTICLECONTACTTYPE, WALLCONTACTTYPE, ContactProperties<T, 1>>(
        particleSystem, solidBoundaries, contactContainer, contactProperties,
        geometry, parameters.get<parameters::CONTACT_BOX_RESOLUTION_PER_DIRECTION>());

    particleManager.execute<process_dynamics<T,PARTICLETYPE>>();

    coupleResolvedParticlesToLattice<T, DESCRIPTOR, PARTICLETYPE, PARTICLECONTACTTYPE, WALLCONTACTTYPE>(
        particleSystem, contactContainer, geometry, lattice, converter, solidBoundaries);

    getResults(myCase, iT, timer, particleSystem);

    lattice.collideAndStream();
  }

  // Run Gnuplot
  if (singleton::mpi().getRank() == 0) {
    if (!system(NULL)) {
      exit (EXIT_FAILURE);
    }
    int ret = system("gnuplot dkt.p");
    if (ret == -1) {
      clout << "Writing Gnuplot failed!" << std::endl;
    }
  }

  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[])
{
  /// === 1st Step: Initialization ===
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");

   MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<EPS                                 >(0.5);
    myCaseParameters.set<PHYS_MAX_T                          >(6.);
    myCaseParameters.set<DOMAIN_EXTENT                       >({0.02, 0.08});
    myCaseParameters.set<CENTER_1                            >({0.01, 0.068});
    myCaseParameters.set<CENTER_2                            >({0.00999, 0.072});
    myCaseParameters.set<PART_RHO                            >(1010.);
    myCaseParameters.set<PART_RADIUS                         >(0.001);
    myCaseParameters.set<CONTACT_BOX_RESOLUTION_PER_DIRECTION>(8);
    myCaseParameters.set<PARTICLE_CONTACT_MATERIAL           >(0);
    myCaseParameters.set<WALL_CONTACT_MATERIAL               >(0);
    myCaseParameters.set<YOUNGS_MODULUS                      >(1e6);
    myCaseParameters.set<PART_POISSON_RATIO                  >(0.3);
    myCaseParameters.set<COEFFICIENT_OF_RESTITUTION          >(0.9);
    myCaseParameters.set<COEFFICIENT_STATIC_FRICTION         >(0.6);
    myCaseParameters.set<COEFFICIENT_KINETIC_FRICTION        >(0.3);
    myCaseParameters.set<PHYS_CHAR_DENSITY                   >(1000);
    myCaseParameters.set<PHYS_CHAR_VELOCITY                  >(0.2);
    myCaseParameters.set<PHYS_CHAR_LENGTH                    >(0.002);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY                 >(1E-6);
    myCaseParameters.set<PHYS_DELTA_X                        >(0.0001);
    myCaseParameters.set<PHYS_DELTA_T                        >(5.e-4);
    myCaseParameters.set<PHYS_VTK_ITER_T                     >(0.125);
    myCaseParameters.set<VTK_ENABLED                         >(true);
    myCaseParameters.set<GNUPLOT_ENABLED                     >(true);
    myCaseParameters.set<GRAVITATIONAL_ACC                   >(-9.81);
    myCaseParameters.set<GRAVITY                             >([](MyCase::ParametersD& p) {
      return Vector<double, 2>{0., p.get<GRAVITATIONAL_ACC>()*
       (1.-p.get<PHYS_CHAR_DENSITY>() / p.get<PART_RHO>())};
    });
  }

  Mesh mesh = createMesh(myCaseParameters);
  MyCase myCase(myCaseParameters, mesh);

  prepareGeometry(myCase);

  prepareLattice(myCase);

  simulate(myCase);
}
