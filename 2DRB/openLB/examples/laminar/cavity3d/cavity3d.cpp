/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2014, 2025 Mathias J. Krause, Yuji (Sam) Shimojima
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

/* cavity3d.cpp:
 * This example illustrates a flow in a cuboid, lid-driven cavity.
 * This version is for parallel use. A version for sequential use
 * is also available.
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;
// === Step 1: Declarations ===
using MyCase = Case<NavierStokes, Lattice<double, descriptors::D3Q19<>>>;


/// @brief Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters)
{
  using T                 = MyCase::value_t;
  const T      physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();
  const T      length = parameters.get<parameters::PHYS_CHAR_LENGTH>(); // length of the cavity in x- and z-direction
  Vector<T, 3> origin(T(0));
  Vector<T, 3> extend(length + 0.5 * physDeltaX);
  IndicatorCuboid3D<T> cuboid(extend, origin);

  Mesh<T, MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

/// @brief Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
/// @note The material numbers will be used to assign physics to lattice nodes
void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;
  using T          = MyCase::value_t;
  auto& sGeometry  = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();

  const T length = parameters.get<parameters::PHYS_CHAR_LENGTH>(); // length of the cavity in x- and z-direction

  Vector<T, 3>         origin((T)0);
  Vector<T, 3>         extend(length + 0.5 * physDeltaX);
  IndicatorCuboid3D<T> cuboid(extend, origin);

  // Sets material number for fluid and boundary
  sGeometry.rename(0, 2, cuboid);
  sGeometry.rename(2, 1, {1, 1, 1});

  T                    eps = physDeltaX;
  Vector<T, 3>         libOrigin(-eps, length - eps, -eps);
  Vector<T, 3>         libExtend(length + 2 * eps, 2 * eps, length + 2 * eps);
  IndicatorCuboid3D<T> lid(libExtend, libOrigin);

  sGeometry.rename(2, 3, 1, lid);

  // Removes all not needed boundary voxels outside the surface
  sGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  sGeometry.innerClean();
  sGeometry.checkForErrors();

  sGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}
/// @brief Set lattice dynamics
/// @param myCase The Case instance which keeps the simulation data
void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  using T          = MyCase::value_t;
  auto& sGeometry  = myCase.getGeometry();
  auto& sLattice   = myCase.getLattice(NavierStokes {});
  auto& parameters = myCase.getParameters();

  clout << "Prepare Lattice ..." << std::endl;
  {
    using namespace olb::parameters;
    sLattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, MyCase::descriptor_t_of<NavierStokes>>>(
        parameters.get<RESOLUTION>(), // resolution: number of voxels per charPhysL
        parameters
            .get<LATTICE_RELAXATION_TIME>(), // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
        parameters.get<PHYS_CHAR_LENGTH>(),  // charPhysLength: reference length of simulation geometry
        parameters.get<
            PHYS_CHAR_VELOCITY>(), // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
        parameters.get<PHYS_CHAR_VISCOSITY>(), // physViscosity: physical kinematic viscosity in __m^2 / s__
        parameters.get<PHYS_CHAR_DENSITY>()    // physDensity: physical density in __kg / m^3__
    );
  }
  sLattice.getUnitConverter().print();

  // Material=1 -->bulk dynamics
  dynamics::set<ConstRhoBGKdynamics>(sLattice, sGeometry.getMaterialIndicator(1));

  // Material=2,3 -->bulk dynamics, velocity boundary
  boundary::set<boundary::InterpolatedVelocity>(sLattice, sGeometry, 2);
  boundary::set<boundary::InterpolatedVelocity>(sLattice, sGeometry, 3);
  clout << "Prepare Lattice ... OK" << std::endl;
}

/// Set initial condition for primal variables (velocity and density)
/// @param myCase The Case instance which keeps the simulation data
/// @note Be careful: initial values have to be set using lattice units
void setInitialValues(MyCase& myCase)
{
  OstreamManager clout(std::cout, "Initialization");
  clout << "lattice initialization ..." << std::endl;

  using T = MyCase::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes {});
  auto& sGeometry = myCase.getGeometry();

  Vector<T,3> uTop{sLattice.getUnitConverter().getCharPhysVelocity(), (T)0, (T)0};
  momenta::setVelocity(sLattice, sGeometry.getMaterialIndicator(3), uTop);
  const T omega = sLattice.getUnitConverter().getLatticeRelaxationFrequency();

  sLattice.setParameter<descriptors::OMEGA>(omega);

  sLattice.initialize();
  clout << "Initialization ... OK" << std::endl;
}

void setTemporalValues(MyCase& myCase, std::size_t iT) {}

void getResults(MyCase& myCase, std::size_t iT, util::Timer<MyCase::value_t>& timer, bool converged)
{
  OstreamManager clout(std::cout, "getResults");
  using T                       = MyCase::value_t;
  auto&               sLattice  = myCase.getLattice(NavierStokes {});
  const auto&         converter = sLattice.getUnitConverter();
  auto&               sGeometry = myCase.getGeometry();
  SuperVTMwriter3D<T> vtmWriter("cavity3d");

  const T logT  = (T)1.0;
  const T saveT = (T)1.0;

  if (iT == 0) {
    SuperLatticeCuboid3D                                                   cuboid(sLattice);
    SuperLatticeRank3D                                                     rank(sLattice);
    SuperLatticeDiscreteNormal3D<T, MyCase::descriptor_t_of<NavierStokes>> discreteNormal(
        sLattice, sGeometry, sGeometry.getMaterialIndicator({2, 3}));
    SuperLatticeDiscreteNormalType3D<T, MyCase::descriptor_t_of<NavierStokes>> discreteNormalType(
        sLattice, sGeometry, sGeometry.getMaterialIndicator({2, 3}));
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.write(discreteNormal);
    vtmWriter.write(discreteNormalType);
    vtmWriter.createMasterFile();
  }

  if ((iT % converter.getLatticeTime(logT) == 0 && iT > 0) || converged) {
    timer.update(iT);
    timer.printStep(2);
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
  }

  if ((iT % converter.getLatticeTime(saveT) == 0 && iT > 0) || converged) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticePhysVelocity3D velocity(sLattice, converter);
    SuperLatticePhysPressure3D pressure(sLattice, converter);

    vtmWriter.addFunctor(velocity);
    vtmWriter.addFunctor(pressure);

    vtmWriter.write(iT);

    // define vector which span the plane
    Vector<T, 3> u(1, 0, 0);
    Vector<T, 3> v(0, 1, 0);
    T            tmp       = (T)converter.getCharPhysLength() / 2.0;
    T            origin[3] = {tmp, tmp, tmp};

    SuperEuklidNorm3D<T>  normVel(velocity);
    BlockReduction3D2D<T> planeReduction(normVel, origin, u, v, 600, BlockDataSyncMode::ReduceOnly);

    heatmap::plotParam<T> plotParam;
    plotParam.maxValue = (T)1.0;
    plotParam.name     = "velocity";
    heatmap::write(planeReduction, iT, plotParam);
  }
}

/// @brief Execute simulation: set initial values and run time loop
/// @param myCase The Case instance which keeps the simulation data
void simulate(MyCase& myCase)
{
  OstreamManager clout(std::cout, "Time marching");
  using T           = MyCase::value_t;
  auto& sLattice    = myCase.getLattice(NavierStokes {});
  auto& parameters  = myCase.getParameters();
  auto  maxLatticeT = sLattice.getUnitConverter().getLatticeTime(parameters.get<parameters::PHYS_MAX_T>());

  util::ValueTracer<T> converge(sLattice.getUnitConverter().getLatticeTime(parameters.get<parameters::PHYS_INTERVAL_T>()),
                                parameters.get<parameters::RESIDUUM>());

  util::Timer<T> timer(maxLatticeT, util::pow<int>(sLattice.getUnitConverter().getResolution(), 3));

  clout << "starting simulation..." << std::endl;
  timer.start();

  for (std::size_t iT = 0; iT <= maxLatticeT; ++iT) {
    if (converge.hasConverged()) {
      clout << "Simulation converged." << std::endl;
      getResults(myCase, iT, timer, converge.hasConverged());
      break;
    }

    setTemporalValues(myCase, iT);

    sLattice.collideAndStream();

    getResults(myCase, iT, timer, converge.hasConverged());
    converge.takeValue(sLattice.getStatistics().getAverageEnergy(), true);
  }
  sLattice.setProcessingContext(ProcessingContext::Evaluation);
  timer.stop();
  timer.printSummary();
  return;
}

// Setup and run a simulation
int main(int argc, char** argv)
{

  initialize(&argc, &argv);
  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION>(30);
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.509);
    myCaseParameters.set<PHYS_CHAR_LENGTH>(1.0);
    myCaseParameters.set<PHYS_CHAR_VELOCITY>(1.0);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(0.001);
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1.0);
    myCaseParameters.set<PHYS_MAX_T>(100);
    myCaseParameters.set<PHYS_INTERVAL_T>(1.0);
    myCaseParameters.set<RESIDUUM>(1.0e-3);
    myCaseParameters.set<PHYS_DELTA_X>([] (MyCase::ParametersD& p) {
      return p.get<PHYS_CHAR_LENGTH>() / p.get<RESOLUTION>();
    });
  }

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);
  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);
  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase);

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);

  return 0;
}
