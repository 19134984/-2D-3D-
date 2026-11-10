/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2019 2025 Mathias J. Krause, Adrian Kummerlaender, Yuji (Sam) Shimojima
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net
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

#include "olb.h"

// Undefine to test a minimal bounce back cavity
#define LID_DRIVEN
//// Use bounce back (velocity) boundaries instead of local velocity
#define LID_DRIVEN_BOUNCE_BACK

using namespace olb;
using namespace olb::names;

// === Step 1: Declarations ===
// Build with FEATURES=FP16 to store POPULATION in f16; compute stays float.
#ifdef FEATURE_FP16
using MyCase = Case<NavierStokes,
                    Lattice<float, descriptors::D3Q19<descriptors::VELOCITY,
                              descriptors::tag::store_precision<descriptors::POPULATION,
                                                                FieldPrecision::F16>>>>;
#else
using MyCase = Case<NavierStokes,
                    Lattice<float, descriptors::D3Q19<descriptors::VELOCITY>>>;
#endif

namespace olb::parameters {

struct CUBOIDS_PER_PROCESS : public descriptors::TYPED_FIELD_BASE<std::size_t, 1> {};

}

Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T                 = MyCase::value_t;
  const T      physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();
  const T      length = parameters.get<parameters::PHYS_CHAR_LENGTH>(); // length of the cavity in x- and z-direction
  Vector<T, 3> origin((T)0);
  Vector<T, 3> extend(length + 0.5 * physDeltaX);
  IndicatorCuboid3D<T> cuboid(extend, origin);

  OstreamManager clout(std::cout, "debug");
  clout << "Meshing…" << std::endl;;

  Mesh<T,MyCase::d> mesh(cuboid,
                         physDeltaX,
                           parameters.get<parameters::CUBOIDS_PER_PROCESS>()
                         * singleton::mpi().getSize(),
                         "splitOnly");
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().buildNeighborCache(parameters.get<parameters::OVERLAP>());

  clout << "Done." << std::endl;;

  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  OstreamManager clout(std::cout, "debug");
  clout << "Prepare geometry…" << std::endl;

  using T          = MyCase::value_t;
  auto& sGeometry  = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();
  const T length     = parameters.get<parameters::PHYS_CHAR_LENGTH>(); // length of the cavity in x- and z-direction

  Vector<T, 3>         origin((T)0);
  Vector<T, 3>         extend(length + 0.25 * physDeltaX);
  IndicatorCuboid3D<T> cuboid(extend, origin);

  sGeometry.rename(0, 2, cuboid);
  sGeometry.rename(2, 1, {1, 1, 1});

  T                    eps = physDeltaX;
  Vector<T, 3>         libOrigin(-eps, length - eps, -eps);
  Vector<T, 3>         libExtend(length + 2 * eps, 2 * eps, length + 2 * eps);
  IndicatorCuboid3D<T> lid(libExtend, libOrigin);

  sGeometry.rename(2, 3, 1, lid);

  //sGeometry.clean(true);
  //sGeometry.innerClean(true);
  //sGeometry.checkForErrors(true);

  clout << "Done." << std::endl;;
}

void prepareLattice(MyCase& myCase) {
  OstreamManager clout(std::cout, "debug");
  clout << "Prepare lattice…" << std::endl;;

  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& sGeometry  = myCase.getGeometry();
  auto& sLattice   = myCase.getLattice(NavierStokes {});
  auto& parameters = myCase.getParameters();

  sLattice.setStatisticsOff();

  {
    using namespace olb::parameters;
    sLattice.setUnitConverter<UnitConverterFromResolutionAndLatticeVelocity<T, DESCRIPTOR>>(
        parameters.get<RESOLUTION>(),       // resolution: number of voxels per charPhysL
        0.01,                               //lattice_char_velocity : characteristic non-dimensional(lattice) valocity
        parameters.get<PHYS_CHAR_LENGTH>(), // charPhysLength: reference length of simulation geometry
        1.0,        // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
        1.0 / 1000.0, // physViscosity: physical kinematic viscosity in __m^2 / s__
        1.0         // physDensity: physical density in __kg / m^3__
    );
  }

  /// Material=1 -->bulk dynamics
  dynamics::set<BGKdynamics>(sLattice, sGeometry, 1);

#ifdef LID_DRIVEN
#ifdef LID_DRIVEN_BOUNCE_BACK
  boundary::set<boundary::BounceBack>(sLattice, sGeometry, 2);
  dynamics::set<BounceBackVelocity>(sLattice, sGeometry, 3);
#else // Local velocity boundaries
  boundary::set<boundary::LocalVelocity<T, DESCRIPTOR, BGKdynamics<T, DESCRIPTOR>>>(sLattice, sGeometry, 2);
  boundary::set<boundary::LocalVelocity<T, DESCRIPTOR, BGKdynamics<T, DESCRIPTOR>>>(sLattice, sGeometry, 3);
#endif
#else
  boundary::set<boundary::BounceBack>(sLattice, sGeometry, 2);
  boundary::set<boundary::BounceBack>(sLattice, sGeometry, 3);
#endif

  clout << "Done." << std::endl;;
}

void setInitialValues(MyCase& myCase) {
  OstreamManager clout(std::cout, "debug");
  clout << "Set initial values…" << std::endl;

  auto& sLattice = myCase.getLattice(NavierStokes {});
  auto& sGeometry = myCase.getGeometry();

#ifdef LID_DRIVEN
  Vector uTop{sLattice.getUnitConverter().getCharPhysVelocity(), 0, 0};
  momenta::setVelocity(sLattice, sGeometry.getMaterialIndicator(3), uTop);
#endif

  sLattice.setParameter<descriptors::OMEGA>(
    sLattice.getUnitConverter().getLatticeRelaxationFrequency());

  sLattice.initialize();

  clout << "Done." << std::endl;;
}

void getResults(MyCase& myCase) {
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

  auto& sLattice = myCase.getLattice(NavierStokes{});

  SuperVTMwriter3D<T> vtmWriter("cavity3d");
  const auto&         converter = sLattice.getUnitConverter();

  SuperLatticeCuboid3D       cuboidF(sLattice);
  SuperLatticeRank3D         rankF(sLattice);
  SuperLatticePhysVelocity3D velocityF(sLattice, converter);
  SuperLatticePhysField3D<T,DESCRIPTOR,fields::moments::U> velocityFieldF(sLattice, converter.getConversionFactorVelocity());
  velocityFieldF.getName() = "moments::U";
  SuperLatticePhysPressure3D pressureF(sLattice, converter);

  vtmWriter.write(cuboidF);
  vtmWriter.write(rankF);
  vtmWriter.write(velocityF);
  vtmWriter.write(velocityFieldF);
  vtmWriter.write(pressureF);
}

void simulate(MyCase& myCase) {
  OstreamManager clout(std::cout, "debug");

  clout << "Simulating…" << std::endl;

  auto& sLattice   = myCase.getLattice(NavierStokes {});
  auto& parameters = myCase.getParameters();

  if (parameters.get<parameters::VTK_ENABLED>()) {
    sLattice.writeSummary();
    getResults(myCase);
  }


  myCase.setReporter<reporter::Timer>("timer", /*writeCsv=*/true);

  const std::size_t maxT = parameters.get<parameters::TIME_STEPS>();
  for (std::size_t iT = 0; iT < maxT; ++iT) {
    sLattice.collideAndStream();
    myCase.report(iT);
  }
  myCase.finalize(maxT);

  if (parameters.get<parameters::VTK_ENABLED>()) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    getResults(myCase);
  }

  clout << "Done." << std::endl;

  singleton::mpi().barrier();
  singleton::exit(0);
}

int main(int argc, char** argv) {
  initialize(&argc, &argv, false, false);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION>(100);
    myCaseParameters.set<TIME_STEPS>(1000);
    myCaseParameters["reporter"]["timer"].set<LATTICE_ITER_T>(100);
    myCaseParameters.set<PHYS_CHAR_LENGTH>(1.0);
    myCaseParameters.set<VTK_ENABLED>(false);
    myCaseParameters.set<CUBOIDS_PER_PROCESS>(1);
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
