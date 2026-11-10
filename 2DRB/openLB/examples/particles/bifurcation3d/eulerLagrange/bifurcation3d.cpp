/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C)
 *  2023      Nicolas Hafen, Frantisek Prinz, Mathias J. Krause
 *  2011-2016 Thomas Henn, Mathias J. Krause, Marie-Luise Maier
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

/* bifurcation3d.cpp:
 * This example examines a steady particulate flow past a bifurcation. At the inlet,
 * an inflow condition with grad_n u = 0 and rho = 1 is implemented.
 * At both outlets, a Poiseuille profile is imposed on the velocity.
 * After a start time, particles are put into the bifurcation at the
 * inlet and experience a stokes drag force.
 *
 * A publication using the same geometry can be found here:
 * http://link.springer.com/chapter/10.1007/978-3-642-36961-2_5
 *  *
 */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace olb::names;

using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q19<>>
>;

// === Step 1: Declarations ===
namespace olb::parameters {

struct INLET_RADIUS : public descriptors::FIELD_BASE<1> {};
struct INLET_CENTER : public descriptors::FIELD_BASE<0, 1> {};
struct INLET_NORMAL : public descriptors::FIELD_BASE<0, 1> {};
struct OUTLET_RADIUS0 : public descriptors::FIELD_BASE<1> {};
struct OUTLET_CENTER0 : public descriptors::FIELD_BASE<0, 1> {};
struct OUTLET_NORMAL0 : public descriptors::FIELD_BASE<0, 1> {};
struct OUTLET_RADIUS1 : public descriptors::FIELD_BASE<1> {};
struct OUTLET_CENTER1 : public descriptors::FIELD_BASE<0, 1> {};
struct OUTLET_NORMAL1 : public descriptors::FIELD_BASE<0, 1> {};
struct FLUID_PHYS_MAX_T : public descriptors::FIELD_BASE<1> {};
struct FLUID_MAX_PHYS_T : public descriptors::FIELD_BASE<1> {};
struct PARTICLE_MAX_PHYS_T : public descriptors::FIELD_BASE<1> {};
struct NO_OF_PARTICLES : public descriptors::FIELD_BASE<1> {};

struct PARTICLE_DYNAMICS_SETUP : public descriptors::TYPED_FIELD_BASE<bool, 1> {};
struct FLUID_EXISTS : public descriptors::TYPED_FIELD_BASE<bool, 1> {};
struct BOUZIDI : public descriptors::TYPED_FIELD_BASE<bool, 1> {};

}

Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters)
{
  using T                 = MyCase::value_t;
  const T      physDeltaX = 2 * parameters.get<parameters::INLET_RADIUS>() / parameters.get<parameters::RESOLUTION>();
  STLreader<T> stlReader("../bifurcation3d.stl", physDeltaX);
  IndicatorLayer3D<T> extendedDomain(stlReader, physDeltaX);
  std::size_t         noOfCuboids = util::max(16, 4 * singleton::mpi().getSize());
  Mesh<T, MyCase::d>  mesh(extendedDomain, physDeltaX, noOfCuboids);
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

//Ensure that parallel mode is used
#ifdef PARALLEL_MODE_MPI

void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T          = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();

  const T      physDeltaX = 2 * parameters.get<parameters::INLET_RADIUS>() / parameters.get<parameters::RESOLUTION>();
  STLreader<T> stlReader("../bifurcation3d.stl", physDeltaX);
  IndicatorLayer3D<T> extendedDomain(stlReader, physDeltaX);
  geometry.rename(0, 2, extendedDomain);
  geometry.rename(2, 1, stlReader);
  geometry.clean();

  // rename the material at the inlet
  IndicatorCircle3D<T>   inletCircle(parameters.get<parameters::INLET_CENTER>(),
                                     parameters.get<parameters::INLET_NORMAL>(),
                                     parameters.get<parameters::INLET_RADIUS>());
  IndicatorCylinder3D<T> inlet(inletCircle, 2 * physDeltaX);

  geometry.rename(2, 3, 1, inlet);

  // rename the material at the outlet0
  IndicatorCircle3D<T>   outletCircle0(parameters.get<parameters::OUTLET_CENTER0>(),
                                       parameters.get<parameters::OUTLET_NORMAL0>(),
                                       0.95 * parameters.get<parameters::OUTLET_RADIUS0>());
  IndicatorCylinder3D<T> outlet0(outletCircle0, 4 * physDeltaX);
  geometry.rename(2, 4, outlet0);

  // rename the material at the outlet1
  IndicatorCircle3D<T>   outletCircle1(parameters.get<parameters::OUTLET_CENTER1>(),
                                       parameters.get<parameters::OUTLET_NORMAL1>(),
                                       0.95 * parameters.get<parameters::OUTLET_RADIUS1>());
  IndicatorCylinder3D<T> outlet1(outletCircle1, 4 * physDeltaX);
  geometry.rename(2, 5, outlet1);

  IndicatorCircle3D<T>   inletCircleExtended(parameters.get<parameters::INLET_CENTER>(),
                                             parameters.get<parameters::INLET_NORMAL>(),
                                             parameters.get<parameters::INLET_RADIUS>() + 2 * physDeltaX);
  IndicatorCylinder3D<T> inletExtended(inletCircleExtended, 2 * physDeltaX);
  geometry.rename(2, 6, inletExtended);

  geometry.clean();
  geometry.innerClean(true);
  geometry.checkForErrors();

  clout << "Prepare Geometry ... OK" << std::endl;
  return;
}

void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(NavierStokes {});


  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR>>(
      parameters.get<parameters::RESOLUTION>(),
      (T)0.557646, // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
      (T)parameters.get<parameters::INLET_RADIUS>() * 2., // charPhysLength: reference length of simulation geometry
      (T)parameters.get<parameters::REYNOLDS>() * 1.5e-5 / (parameters.get<parameters::INLET_RADIUS>() * 2.),
      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
      (T)parameters.get<parameters::PHYS_CHAR_VISCOSITY>() , // physViscosity: physical kinematic viscosity in __m^2 / s__
      (T)parameters.get<parameters::PHYS_CHAR_DENSITY>()   // physDensity: physical density in __kg / m^3__
  );
  auto& converter   = lattice.getUnitConverter();
  converter.print();
  const T omega     = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  dynamics::set<BGKdynamics>(lattice, geometry, 1);


  if (parameters.get<parameters::BOUZIDI>()) {
    STLreader<T> stlReader("../bifurcation3d.stl", converter.getPhysDeltaX());
    setBouzidiBoundary(lattice, geometry, 2, stlReader);
  }
  else {
    boundary::set<boundary::BounceBack>(lattice, geometry, 2);
  }

  // Material=3 -->bulk dynamics (inflow)
  dynamics::set<BGKdynamics>(lattice, geometry, 3);

  // Material=4 -->bulk dynamics (outflow)
  dynamics::set<BGKdynamics>(lattice, geometry, 4);
  dynamics::set<BGKdynamics>(lattice, geometry, 5);

  // Setting of the boundary conditions
  boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 3);
  boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 4);
  boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 5);

  lattice.setParameter<descriptors::OMEGA>(omega);

  clout << "Prepare Lattice ... OK" << std::endl;
  return;
}

void setInitialValues(MyCase& myCase)
{
  OstreamManager clout(std::cout, "setInitialValues");
  clout << "Set Initial Values ..." << std::endl;
  auto& lattice  = myCase.getLattice(NavierStokes {});

  lattice.initialize();

  clout << "Set Initial Values ... OK" << std::endl;
}

void setTemporalValues(MyCase& myCase, std::size_t iT)
{
  using T          = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(NavierStokes {});
  auto& converter  = lattice.getUnitConverter();
  // No of time steps for smooth start-up
  std::size_t             iTmaxStart = converter.getLatticeTime(0.8 * parameters.get<parameters::FLUID_MAX_PHYS_T>());
  SinusStartScale<T, int> startScale(iTmaxStart, T(1));
  int                     iTvec[1] = {int(iT)};
  T                       frac[1]  = {T(0)};
  startScale(frac, iTvec);
  T maxVelocity = frac[0] * converter.getCharPhysVelocity() * 3. / 4. *
                  util::pow(parameters.get<parameters::INLET_RADIUS>(), 2) /
                  util::pow(parameters.get<parameters::OUTLET_RADIUS0>(), 2);

  Vector<T, 3>          outletNormal0 = parameters.get<parameters::OUTLET_NORMAL0>();
  Vector<T, 3>          outletNormal1 = parameters.get<parameters::OUTLET_NORMAL1>();
  Vector<T, 3>          outletCenter0 = parameters.get<parameters::OUTLET_CENTER0>();
  Vector<T, 3>          outletCenter1 = parameters.get<parameters::OUTLET_CENTER1>();
  CirclePoiseuille3D<T> poiseuilleU4(outletCenter0[0], outletCenter0[1], outletCenter0[2], outletNormal0[0],
                                     outletNormal0[1], outletNormal0[2], parameters.get<parameters::OUTLET_RADIUS0>() * 0.95, -maxVelocity);

  CirclePoiseuille3D<T> poiseuilleU5(outletCenter1[0], outletCenter1[1], outletCenter1[2], outletNormal1[0],
                                     outletNormal1[1], outletNormal1[2], parameters.get<parameters::OUTLET_RADIUS1>() * 0.95, -maxVelocity);

  momenta::setVelocity(lattice, geometry.getMaterialIndicator(4), poiseuilleU4);
  momenta::setVelocity(lattice, geometry.getMaterialIndicator(5), poiseuilleU5);
}

template<typename PARTICLESYSTEM>
void prepareParticles(MyCase& myCase,
                      PARTICLESYSTEM& superParticleSystem,
                      SolidBoundary<MyCase::value_t, 3>& wall,
                      STLreader<MyCase::value_t>& stlReader)
{
  OstreamManager clout(std::cout, "prepareParticles");
  clout << "Prepare Particles ..." << std::endl;
  using namespace particles;
  using namespace particles::subgrid;
  using namespace particles::dynamics;
  using namespace particles::creators;
  using namespace particles::io;
  using T          = MyCase::value_t;
  using PARTICLETYPE = SubgridParticle3D;
  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(NavierStokes {});
  auto& converter  = lattice.getUnitConverter();

  //Create material indicators for particle boundaries
  std::vector<int>                              wallMaterials {2};
  std::shared_ptr<SuperIndicatorMaterial<T, 3>> wallMaterialIndicator =
      std::make_shared<SuperIndicatorMaterial<T, 3>>(geometry, wallMaterials);
  std::vector<int>                              outletMaterials {4, 5};
  std::shared_ptr<SuperIndicatorMaterial<T, 3>> outletMaterialIndicator =
      std::make_shared<SuperIndicatorMaterial<T, 3>>(geometry, outletMaterials);

  //Add selected particle dynamics
  if (parameters.get<parameters::PARTICLE_DYNAMICS_SETUP>() == 1) {
    superParticleSystem.template defineDynamics<
        VerletParticleDynamicsMaterialAwareWallCaptureAndEscape<T, PARTICLETYPE>
        >( wall, wallMaterialIndicator, outletMaterialIndicator );
  }
  else {
    superParticleSystem.template defineDynamics<
        VerletParticleDynamicsMaterialCaptureAndEscape<T, PARTICLETYPE>
        >( wallMaterialIndicator, outletMaterialIndicator );
  }

  // particles generation at inlet
  Vector<T, 3> center(parameters.get<parameters::INLET_CENTER>());
  center[2] = 0.074;
  IndicatorCircle3D<T>   inflowCircle(center, parameters.get<parameters::INLET_NORMAL>(),
                                      parameters.get<parameters::INLET_RADIUS>() - converter.getPhysDeltaX() * 2.5);
  IndicatorCylinder3D<T> inletCylinder(inflowCircle, 0.01 * converter.getPhysDeltaX());

  Randomizer<T> randomizer;
  addParticles(superParticleSystem, inletCylinder, parameters.get<parameters::PART_RHO>(),
               parameters.get<parameters::PART_RADIUS>(), parameters.get<parameters::NO_OF_PARTICLES>(), randomizer);

  //Print super particle system summary
  superParticleSystem.print();
  VTUwriter<T, PARTICLETYPE, true> superParticleWriter("particles_master", false, false);

  SuperParticleGroupedFieldF<T, PARTICLETYPE, GENERAL, POSITION>    particlePosition(superParticleSystem);
  SuperParticleGroupedFieldF<T, PARTICLETYPE, PHYSPROPERTIES, MASS> particleMass(superParticleSystem);
  clout << "Prepare Particles ... OK" << std::endl;
}

template<typename PARTICLESYSTEM>
bool getResults(MyCase& myCase, std::size_t iT, Timer<MyCase::value_t>& fluidTimer,
                STLreader<MyCase::value_t>& stlReader,
                PARTICLESYSTEM& superParticleSystem,
                Timer<MyCase::value_t>& particleTimer)
{
  OstreamManager clout(std::cout, "getResults");
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using PARTICLETYPE = SubgridParticle3D;
  using namespace particles;
  using namespace particles::subgrid;
  using namespace particles::dynamics;
  using namespace particles::io;
  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(NavierStokes {});
  auto& converter  = lattice.getUnitConverter();

  SuperVTMwriter3D<T> vtmWriter("bifurcation3d");
  SuperVTMwriter3D<T> vtmWriterStartTime("startingTimeBifurcation3d");

  //Create material indicators for particle boundaries
  std::vector<int>                              wallMaterials {2};
  std::shared_ptr<SuperIndicatorMaterial<T, 3>> wallMaterialIndicator =
      std::make_shared<SuperIndicatorMaterial<T, 3>>(geometry, wallMaterials);
  std::vector<int>                              outletMaterials {4, 5};
  std::shared_ptr<SuperIndicatorMaterial<T, 3>> outletMaterialIndicator =
      std::make_shared<SuperIndicatorMaterial<T, 3>>(geometry, outletMaterials);

  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(lattice, converter);
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(lattice, converter);
  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(pressure);

  vtmWriterStartTime.addFunctor(velocity);
  vtmWriterStartTime.addFunctor(pressure);

  VTUwriter<T, PARTICLETYPE, true> superParticleWriter("particles_master", false, false);

  SuperParticleGroupedFieldF<T, PARTICLETYPE, GENERAL, POSITION>      particlePosition(superParticleSystem);
  SuperParticleGroupedFieldF<T, PARTICLETYPE, PHYSPROPERTIES, MASS>   particleMass(superParticleSystem);
  SuperParticleGroupedFieldF<T, PARTICLETYPE, PHYSPROPERTIES, RADIUS> particleRadius(superParticleSystem);
  SuperParticleGroupedFieldF<T, PARTICLETYPE, MOBILITY, VELOCITY>     particleVelocity(superParticleSystem);
  SuperParticleGroupedFieldF<T, PARTICLETYPE, FORCING, FORCE>         particleForce(superParticleSystem);
  SuperParticleGroupedFieldF<T, PARTICLETYPE, DYNBEHAVIOUR, ACTIVE>   particleActivity(superParticleSystem);
  superParticleWriter.addFunctor(particlePosition, "Position");
  superParticleWriter.addFunctor(particleMass, "Mass");
  superParticleWriter.addFunctor(particleRadius, "Radius");
  superParticleWriter.addFunctor(particleVelocity, "Velocity");
  superParticleWriter.addFunctor(particleForce, "Force");
  superParticleWriter.addFunctor(particleActivity, "Active");

  std::size_t fluidMaxT = converter.getLatticeTime(parameters.get<parameters::FLUID_MAX_PHYS_T>());

  if (iT == 0) {
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(lattice);
    SuperLatticeRank3D<T, DESCRIPTOR>   rank(lattice);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();
    superParticleWriter.createMasterFile();
    vtmWriterStartTime.createMasterFile();

    clout << "N=" << parameters.get<parameters::RESOLUTION>()
          << "; maxTimeSteps(fluid)=" << converter.getLatticeTime(parameters.get<parameters::FLUID_MAX_PHYS_T>())
          << "; noOfCuboid=" << geometry.getCuboidDecomposition().size()
          << "; Re=" << parameters.get<parameters::REYNOLDS>()
          << "; noOfParticles=" << parameters.get<parameters::NO_OF_PARTICLES>()
          << "; maxTimeSteps(particle)=" << converter.getLatticeTime(parameters.get<parameters::PARTICLE_MAX_PHYS_T>())
          << "; St="
          << parameters.get<parameters::STOKES>()
          << std::endl;
  }

  // Writes the .vtk and .gif files
  if (iT < fluidMaxT) {
    if (!parameters.get<parameters::FLUID_EXISTS>()) {
      vtmWriterStartTime.write(iT);
      SuperEuklidNorm3D<T>  normVel(velocity);
      BlockReduction3D2D<T> planeReduction(normVel, {0, -1, 0}, 600, BlockDataSyncMode::ReduceOnly);
      // write output as JPEG
      heatmap::write(planeReduction, iT);
    }

    // Timer statics
    fluidTimer.update(iT);
    fluidTimer.printStep();

    vtmWriter.write(iT);
    // Lattice statistics
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));

    // Flux at the inlet and outlet regions
    const std::vector<int> materials = {1, 3, 4, 5};

    IndicatorCircle3D<T>                inlet(parameters.get<parameters::INLET_CENTER>() +
                                                  2. * converter.getPhysDeltaX() * parameters.get<parameters::INLET_NORMAL>(),
                                              parameters.get<parameters::INLET_NORMAL>(),
                                              parameters.get<parameters::INLET_RADIUS>() + 2. * converter.getPhysDeltaX());
    SuperPlaneIntegralFluxVelocity3D<T> vFluxInflow(lattice, converter, geometry, inlet, materials);
    vFluxInflow.print("inflow", "ml/s");
    SuperPlaneIntegralFluxPressure3D<T> pFluxInflow(lattice, converter, geometry, inlet, materials);
    pFluxInflow.print("inflow", "N", "Pa");

    IndicatorCircle3D<T>                outlet0(parameters.get<parameters::OUTLET_CENTER0>() +
                                                    2. * converter.getPhysDeltaX() * parameters.get<parameters::OUTLET_NORMAL0>(),
                                                parameters.get<parameters::OUTLET_NORMAL0>(),
                                                parameters.get<parameters::OUTLET_RADIUS0>() + 2. * converter.getPhysDeltaX());
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOutflow0(lattice, converter, geometry, outlet0, materials);
    vFluxOutflow0.print("outflow0", "ml/s");
    SuperPlaneIntegralFluxPressure3D<T> pFluxOutflow0(lattice, converter, geometry, outlet0, materials);
    pFluxOutflow0.print("outflow0", "N", "Pa");

    IndicatorCircle3D<T>                outlet1(parameters.get<parameters::OUTLET_CENTER1>() +
                                                    2. * converter.getPhysDeltaX() * parameters.get<parameters::OUTLET_NORMAL1>(),
                                                parameters.get<parameters::OUTLET_NORMAL1>(),
                                                parameters.get<parameters::OUTLET_RADIUS1>() + 2. * converter.getPhysDeltaX());
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOutflow1(lattice, converter, geometry, outlet1, materials);
    vFluxOutflow1.print("outflow1", "ml/s");
    SuperPlaneIntegralFluxPressure3D<T> pFluxOutflow1(lattice, converter, geometry, outlet1, materials);
    pFluxOutflow1.print("outflow1", "N", "Pa");
  }

  // Particle output
  if (iT >= fluidMaxT) {
    // advance particle timer
    particleTimer.print(iT - fluidMaxT);

    //delete invalidated particles
    purgeInvalidParticles<T, PARTICLETYPE>(superParticleSystem);

    //Define materials for capture rate
    std::vector<int>             materialsOutput {4, 5};
    SuperIndicatorMaterial<T, 3> materialIndicatorOutput(geometry, materialsOutput);
    std::size_t                  noActive;
    captureStatistics(superParticleSystem, materialIndicatorOutput, noActive);

    superParticleWriter.write(iT);

    // true as long as certain amount of active particles
    if ((noActive < 0.001 * parameters.get<parameters::NO_OF_PARTICLES>() &&
        iT > 0.9 * converter.getLatticeTime(parameters.get<parameters::FLUID_MAX_PHYS_T>() +
                                            parameters.get<parameters::PARTICLE_MAX_PHYS_T>()))
        || noActive == 0) {
      return false;
    }
  }
  return true;
}

void simulate( MyCase& myCase )
{
  using T                   = MyCase::value_t;
  using DESCRIPTOR          = MyCase::descriptor_t_of<NavierStokes>;
  using PARTICLETYPE        = SubgridParticle3D;
  using namespace particles;
  using namespace particles::subgrid;
  using namespace particles::communication;
  using namespace particles::dynamics;

  auto& parameters          = myCase.getParameters();
  auto& geometry            = myCase.getGeometry();
  auto& lattice             = myCase.getLattice(NavierStokes {});
  auto& converter           = lattice.getUnitConverter();

  OstreamManager clout(std::cout, "simulate");

    /// === Step 9: Create SuperParticleSystem and ParticleManager ===
  SuperParticleSystem<T, PARTICLETYPE> superParticleSystem(geometry);

  ParticleManager<T, DESCRIPTOR, PARTICLETYPE> particleManager(
      superParticleSystem, geometry,  lattice, converter);

  STLreader<T>   stlReader("../bifurcation3d.stl", converter.getPhysDeltaX());

  const unsigned latticeMaterial = 2; //Material number of wall
  const unsigned contactMaterial = 0; //Material identifier (only relevant for contact model)
  SolidBoundary<T,3> wall( std::make_unique<IndicInverse<T,3>>(stlReader),
                           latticeMaterial, contactMaterial );

  Timer<T>       fluidTimer(converter.getLatticeTime(parameters.get<parameters::FLUID_MAX_PHYS_T>()),
                            geometry.getStatistics().getNvoxel());

  Timer<T> particleTimer(converter.getLatticeTime(parameters.get<parameters::PARTICLE_MAX_PHYS_T>()),
                         parameters.get<parameters::NO_OF_PARTICLES>());
  fluidTimer.start();

  std::size_t iT         = 0;
  std::size_t iTperiod   = converter.getLatticeTime(0.5 * parameters.get<parameters::PHYS_STAT_ITER_T>());
  std::size_t iTmaxStart = converter.getLatticeTime(0.8 * parameters.get<parameters::FLUID_MAX_PHYS_T>());


  /// === Step 10.1 loading fluid solution or simulating ===
  if (!(lattice.load("fluidSolution_N" + std::to_string(parameters.get<parameters::RESOLUTION>())))) {
    clout << "No fluid data available, generating new data" << std::endl;
    parameters.set<parameters::FLUID_EXISTS>(false);

    for (; iT <= converter.getLatticeTime(parameters.get<parameters::FLUID_MAX_PHYS_T>()); ++iT) {
      if (iT <= iTmaxStart && iT % iTperiod == 0) {
        setTemporalValues(myCase, iT);
      }
      lattice.collideAndStream();
      if (iT % iTperiod == 0) {
        getResults(myCase, iT, fluidTimer, stlReader, superParticleSystem, particleTimer);
      }
    }

    fluidTimer.stop();
    fluidTimer.printSummary();
    lattice.communicate();

    // calculated results are written in a file
    lattice.save("fluidSolution_N" + std::to_string(parameters.get<parameters::RESOLUTION>()));
  }

  // if there exists already data of the fluid from an earlier calculation, this is used
  else {
    clout << "using existing fluid data" << std::endl;
    iT = converter.getLatticeTime(parameters.get<parameters::FLUID_MAX_PHYS_T>());
    getResults(myCase, iT, fluidTimer, stlReader, superParticleSystem, particleTimer);
    fluidTimer.stop();
  }

  /// === Step 10.2: Prepare Particle dynamics and add particles ===
  prepareParticles(myCase, superParticleSystem, wall, stlReader);

  /// === Step 10.3 Particle simulation ===
  initializeParticleVelocity(lattice, geometry, converter, superParticleSystem);
  particleTimer.start();
  clout << " starting particle simulation " << std::endl;

  for (; iT <= converter.getLatticeTime(parameters.get<parameters::FLUID_MAX_PHYS_T>() +
                                        parameters.get<parameters::PARTICLE_MAX_PHYS_T>());
       ++iT) {
    // particles simulation starts after run up time is over
    particleManager.execute<
                    couple_lattice_to_particles<T, DESCRIPTOR, PARTICLETYPE>,
                    process_dynamics<T, PARTICLETYPE>,
                    update_particle_core_distribution<T, PARTICLETYPE>
    >();
    if (iT % iTperiod == 0) {
      if(getResults(myCase, iT, fluidTimer, stlReader, superParticleSystem, particleTimer)) {
        continue;
      } else {
        break;
      }
    }
  }
  particleTimer.stop();
  particleTimer.printSummary();
}

#endif //PARALLEL_MODE_MPI

int main(int argc, char* argv[])
{

  initialize(&argc, &argv);

  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<REYNOLDS>(50);
    myCaseParameters.set<RESOLUTION>(19);
    myCaseParameters.set<PHYS_STAT_ITER_T>(0.2); // time interval for writeout and updating boundary values in s
    myCaseParameters.set<PART_RADIUS>(1.5e-4);
    myCaseParameters.set<PART_RHO>(998.2);
    myCaseParameters.set<FLUID_MAX_PHYS_T>(5);     // max. simulation time in s
    myCaseParameters.set<PARTICLE_MAX_PHYS_T>(20); // max. particle simulation time in s
    myCaseParameters.set<NO_OF_PARTICLES>(1000);
    myCaseParameters.set<PARTICLE_DYNAMICS_SETUP>(1);// 0: based on material number,  1: based on more accurate stl description
    myCaseParameters.set<INLET_CENTER>({0., 0., 0.0786395});
    myCaseParameters.set<OUTLET_CENTER0>({-0.0235929682287551, -0.000052820468762797, -0.021445708949909});
    myCaseParameters.set<OUTLET_CENTER1>({0.0233643529416147, 0.00000212439067050152, -0.0211994104877918});
    myCaseParameters.set<INLET_RADIUS>(0.00999839);
    myCaseParameters.set<OUTLET_RADIUS0>(0.007927);
    myCaseParameters.set<OUTLET_RADIUS1>(0.00787134);
    myCaseParameters.set<INLET_NORMAL>({0., 0., -1.});
    myCaseParameters.set<OUTLET_NORMAL0>({0.505126, -0.04177, 0.862034});
    myCaseParameters.set<OUTLET_NORMAL1>({-0.483331, -0.0102764, 0.875377});
    myCaseParameters.set<BOUZIDI>(1);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(1.5e-5);
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1.225);
    myCaseParameters.set<STOKES>([] (MyCase::ParametersD& p) {
              return  p.get<PART_RHO>() * p.get<PART_RADIUS>()
                    * p.get<PART_RADIUS>() * p.get<REYNOLDS>()
                    / (18. *p.get<INLET_RADIUS>() * p.get<PHYS_CHAR_DENSITY>()
                    * p.get<INLET_RADIUS>() );
                });
    myCaseParameters.set<FLUID_EXISTS> (1);

  }
  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase = MyCase(myCaseParameters, mesh);

#ifdef PARALLEL_MODE_MPI

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase);

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);

#else
  std::cerr << std::endl << "ERROR: Subgrid particles can only be used with MPI!" << std::endl << std::endl;
#endif //PARALLEL_MODE_MPI
}
