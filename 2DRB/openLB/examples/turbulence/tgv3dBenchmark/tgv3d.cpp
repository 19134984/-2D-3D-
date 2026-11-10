/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Adrian Kummerlaender
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

// Taylor-Green vortex throughput benchmark: the cavity3dBenchmark harness
// (same --RESOLUTION/--TIME_STEPS CLI, same [result] MLUPs output) applied to a
// triply-periodic Smagorinsky-LES case. Exercises the periodic halo-exchange path
// (every domain face wraps around, crossing rank boundaries under decomposition).

#include <olb.h>

#include <numbers>

using namespace olb;
using namespace olb::names;
using namespace olb::descriptors;

// === Step 1: Declarations ===
// Build with FEATURES=FP16 to store POPULATION in f16; compute stays float.
#ifdef FEATURE_FP16
using MyCase = Case<NavierStokes,
                    Lattice<FLOATING_POINT_TYPE,
                            D3Q19<tag::store_precision<POPULATION, FieldPrecision::F16>>>>;
#else
using MyCase = Case<NavierStokes,
                    Lattice<FLOATING_POINT_TYPE, D3Q19<>>>;
#endif

// TGV_PLAIN_BGK swaps in plain BGK (same D3Q19<> layout, same populations traffic,
// no strain-rate / effective-omega compute) to isolate the LES collision cost.
#ifdef TGV_PLAIN_BGK
using BulkDynamics = BGKdynamics<MyCase::value_t, MyCase::descriptor_t>;
#else
using BulkDynamics = SmagorinskyBGKdynamics<MyCase::value_t, MyCase::descriptor_t>;
#endif

namespace olb::parameters {

struct CUBOIDS_PER_PROCESS : public descriptors::TYPED_FIELD_BASE<std::size_t, 1> {};

}


Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const T N          = parameters.get<parameters::RESOLUTION>();
  const T physDeltaX = 1. / int(std::nearbyint(N / (2 * std::numbers::pi_v<T>)));

  Mesh<T,MyCase::d> mesh(0, physDeltaX, N,
                           parameters.get<parameters::CUBOIDS_PER_PROCESS>()
                         * singleton::mpi().getSize());
  mesh.getCuboidDecomposition().setPeriodicity({true, true, true});
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().buildNeighborCache(parameters.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  auto& geometry = myCase.getGeometry();
  geometry.rename(0, 1);                 // the whole periodic box is bulk fluid
  geometry.communicate();
}

void prepareLattice(MyCase& myCase) {
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes{});

  const T N   = parameters.get<parameters::RESOLUTION>();
  const T Re  = parameters.get<parameters::REYNOLDS>();
  [[maybe_unused]] const T smagoConst = parameters.get<parameters::SMAGORINSKY>();

  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    int(std::nearbyint(N / (2 * std::numbers::pi_v<T>))),
    0.507639,   // latticeRelaxationTime
    1,          // charPhysLength
    1,          // charPhysVelocity
    1. / Re,    // physViscosity
    1.0);       // physDensity

  lattice.setStatisticsOff();

  dynamics::set<BulkDynamics>(lattice, myCase.getGeometry(), 1);
#ifndef TGV_PLAIN_BGK
  lattice.setParameter<collision::LES::SMAGORINSKY>(smagoConst);
#endif
  lattice.setParameter<descriptors::OMEGA>(lattice.getUnitConverter().getLatticeRelaxationFrequency());
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});

  AnalyticalConst3D<T,T> rho(1);
  const T u0 = lattice.getUnitConverter().getCharPhysVelocity();
  AnalyticalFfromCallableF<MyCase::d,T,T> u([u0](Vector<T,3> p) -> Vector<T,3> {
    return { u0*util::sin(p[0])*util::cos(p[1])*util::cos(p[2]),
            -u0*util::cos(p[0])*util::sin(p[1])*util::cos(p[2]),
             0};
  });

  auto domain = myCase.getGeometry().getMaterialIndicator(1);
  momenta::setDensity(lattice, domain, rho);
  momenta::setVelocity(lattice, domain, u);
  lattice.initialize();
}

void getResults(MyCase& myCase) {
  using T        = MyCase::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  const auto& converter = sLattice.getUnitConverter();

  SuperVTMwriter3D<T> vtmWriter("tgv3d");
  SuperLatticeCuboid3D        cuboidF(sLattice);
  SuperLatticeRank3D          rankF(sLattice);
  SuperLatticePhysVelocity3D  velocityF(sLattice, converter);
  vtmWriter.write(cuboidF);
  vtmWriter.write(rankF);
  vtmWriter.write(velocityF);
}

void simulate(MyCase& myCase) {
  OstreamManager clout(std::cout, "debug");
  clout << "Simulating…" << std::endl;

  auto& sLattice   = myCase.getLattice(NavierStokes{});
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

  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION>(100);
    myCaseParameters.set<TIME_STEPS>(1000);
    myCaseParameters["reporter"]["timer"].set<LATTICE_ITER_T>(100);
    myCaseParameters.set<REYNOLDS>(800);
    myCaseParameters.set<SMAGORINSKY>(0.1);
    myCaseParameters.set<VTK_ENABLED>(false);
    myCaseParameters.set<CUBOIDS_PER_PROCESS>(1);
  }

  Mesh mesh = createMesh(myCaseParameters);
  MyCase myCase(myCaseParameters, mesh);

  prepareGeometry(myCase);
  prepareLattice(myCase);
  setInitialValues(myCase);
  simulate(myCase);

  return 0;
}
