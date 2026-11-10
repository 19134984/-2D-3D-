/*  Lattice Boltzmann sample, written in C++, using the OpenLB library
 *
 *  Copyright (C) 2025 Fedor Bukreev, Adrian Kummerländer,
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

/* nozzle2d.cpp:
 * This is learning example for dynamics and BCs changes in a turbulent nozzle
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<float, descriptors::D2Q9<>>
>;

namespace olb::parameters {

struct INLET_HEIGHT : public descriptors::FIELD_BASE<1> {};
struct START_FACTOR : public descriptors::FIELD_BASE<1> {};

}

// Create mesh for simulation case
Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& params)
{
  using T                     = MyCase::value_t;
  Vector               extent = params.get<parameters::DOMAIN_EXTENT>();
  std::vector<T>       origin(2, T());
  IndicatorCuboid2D<T> cuboid(extent, origin);

  const T            physDeltaX = params.get<parameters::INLET_HEIGHT>() / params.get<parameters::RESOLUTION>();
  Mesh<T, MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry(MyCase& myCase)
{
  using T        = MyCase::value_t;
  auto& params   = myCase.getParameters();
  auto& geometry = myCase.getGeometry();

  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  // Get relevant parameters for geometry creation
  T lengthX = params.get<parameters::DOMAIN_EXTENT>()[0];
  T lengthY = params.get<parameters::DOMAIN_EXTENT>()[1];
  T inletHeight  = params.get<parameters::INLET_HEIGHT>();
  T deltaX  = inletHeight / params.get<parameters::RESOLUTION>();

  Vector<T, 2> extend(lengthX, lengthY);
  Vector<T, 2> origin;

  geometry.rename(0, 2);
  geometry.rename(2, 1, {1, 1});

  // Set material number for inflow
  IndicatorCuboid2D<T> inlet({T(2)*deltaX, inletHeight}, {-deltaX, (lengthY-inletHeight)/T(2)});
  geometry.rename(2, 3, 1, inlet);

  // Set material number for outflow
  origin[0] = lengthX - deltaX;
  IndicatorCuboid2D<T> outflow(extend, origin);
  geometry.rename(2, 4, outflow);

  // Removes all not needed boundary voxels outside the surface
  geometry.clean(false, {1, 5});
  geometry.checkForErrors();

  clout << "Prepare Geometry ... OK" << std::endl;

  geometry.getStatistics().print();
}

void prepareLattice(MyCase& myCase)
{
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& params     = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes {});
  auto& geometry   = myCase.getGeometry();

  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  // Get parameters for Unit converter
  int resolution      = params.get<parameters::RESOLUTION>();
  int Re              = params.get<parameters::REYNOLDS>();
  T   latticeVelocity = params.get<parameters::LATTICE_CHAR_VELOCITY>();
  T   physViscosity   = params.get<parameters::PHYS_CHAR_VISCOSITY>();
  T   physDensity     = params.get<parameters::PHYS_CHAR_DENSITY>();
  T   inletHeight     = params.get<parameters::INLET_HEIGHT>();

  // Set UnitConverter
  lattice.setUnitConverter<UnitConverterFromResolutionAndLatticeVelocity<T, DESCRIPTOR>>(
      (T)resolution,  // resolution: number of voxels per charPhysL
      (T)latticeVelocity,        // Max CFL
      (T)inletHeight, // charPhysLength: reference length of simulation geometry
      (T)Re * physViscosity /
          inletHeight,  // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
      (T)physViscosity, // physViscosity: physical kinematic viscosity in __m^2 / s__
      (T)physDensity
  );
  const auto& converter = lattice.getUnitConverter();
  converter.print();

  // Material=1 --> Bulk dynamics
  dynamics::set<BGKdynamics>(lattice, geometry, 1);

  // Material=2 -->No Slip
  boundary::set<boundary::FullSlip>(lattice, geometry, 2);

  // Material=3 -->Fixed Velocity
  boundary::set<boundary::LocalVelocity>(lattice, geometry, 3);

  boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 4);

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase)
{
  using T          = MyCase::value_t;
  auto& lattice    = myCase.getLattice(NavierStokes {});
  auto& converter  = lattice.getUnitConverter();

  const T omega = converter.getLatticeRelaxationFrequency();
  lattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  lattice.initialize();
}

void setTemporalValues(MyCase& myCase, std::size_t iT)
{
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});
  const auto& converter = lattice.getUnitConverter();

  const T startUpFactor = params.get<parameters::START_FACTOR>();
  const T maxStartTPhys = params.get<parameters::PHYS_MAX_T>() * startUpFactor;
  const auto maxStartT  = converter.getLatticeTime(maxStartTPhys);
  const size_t iTupdate   = (size_t) params.get<parameters::LATTICE_INTERVAL_T>();

  if (iT%iTupdate == 0 && iT <= maxStartT) {
    PolynomialStartScale<T,std::size_t> scale(maxStartT, 1);
    T frac{};
    scale(&frac, &iT);
    PowerLawTurbulent2D<T> velocity(geometry, 3, converter.getCharPhysVelocity()*frac, T(), 2, 0.03);
    momenta::setVelocity(lattice, geometry.getMaterialIndicator(3), velocity);
    lattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

void getResults(MyCase& myCase, std::size_t iT)
{

  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& params     = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes {});

  OstreamManager clout(std::cout, "getResults");

  const int vtkIter  = lattice.getUnitConverter().getLatticeTime(params.get<parameters::PHYS_VTK_ITER_T>());

  if (iT == 0) {
    SuperVTMwriter2D<T> vtmWriter("nozzle2d");
    vtmWriter.createMasterFile();
  }

  if (iT % vtkIter == 0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperVTMwriter2D<T>                         vtmWriter("nozzle2d");
    SuperLatticePhysVelocity2D<T, DESCRIPTOR>   velocity(lattice, lattice.getUnitConverter());
    SuperLatticePhysPressure2D<T, DESCRIPTOR>   pressure(lattice, lattice.getUnitConverter());
    vtmWriter.addFunctor(velocity);
    vtmWriter.addFunctor(pressure);
    vtmWriter.write(iT);
  }

  myCase.report(iT);
}

void simulate(MyCase& myCase)
{
  using T           = MyCase::value_t;
  auto&   params    = myCase.getParameters();
  auto&   lattice   = myCase.getLattice(NavierStokes {});
  const T physMaxT  = params.get<parameters::PHYS_MAX_T>();

  const std::size_t iTmax = myCase.getLattice(NavierStokes {}).getUnitConverter().getLatticeTime(physMaxT);

  for (std::size_t iT = 0; iT < iTmax; ++iT) {
    // === Update Inlet Values ===
    setTemporalValues(myCase, iT);

    // === Collide and Stream Execution ===
    lattice.collideAndStream();

    // === Computation and Output of the Results ===
    getResults(myCase, iT);
  }

  myCase.finalize(iTmax);
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    using T          = MyCase::value_t;
    myCaseParameters.set<RESOLUTION>(30);
    myCaseParameters.set<REYNOLDS>(500);
    myCaseParameters.set<LATTICE_CHAR_VELOCITY>(0.05);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>((T)1.e-5);
    myCaseParameters.set<PHYS_CHAR_DENSITY>((T)1.);
    myCaseParameters.set<PHYS_MAX_T>((T)500.);
    myCaseParameters.set<DOMAIN_EXTENT>(Vector<T,2> {6, 2});
    myCaseParameters.set<INLET_HEIGHT>((T)0.2);
    myCaseParameters.set<START_FACTOR>((T)0.1);
    myCaseParameters.set<LATTICE_INTERVAL_T>(1);
    myCaseParameters["reporter"]["statistics"].set<PHYS_ITER_T>(1.0);
    myCaseParameters["reporter"]["timer"].set<PHYS_ITER_T>(1.0);
    myCaseParameters.set<PHYS_VTK_ITER_T>(1.0);
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
