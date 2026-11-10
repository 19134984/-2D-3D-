/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2019 Mathias J. Krause, Adrian Kummerlaender
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

#include <olb.h>

using namespace olb;

using DESCRIPTOR = descriptors::D2Q9<>;
using T = double;
using BulkDynamics = BGKdynamics<T,DESCRIPTOR>;

#define LID_DRIVEN

void prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter,
                     IndicatorF2D<T>& indicator,
                     SuperGeometry<T,2>& superGeometry)
{
  // Sets material number for fluid and boundary
  superGeometry.rename(0,2,indicator);
  superGeometry.rename(2,1,{1,1});

  T eps = converter.getConversionFactorLength();
  Vector<T,2> origin(-eps, converter.getCharPhysLength() - eps);
  Vector<T,2> extend(converter.getCharPhysLength() + 2*eps, 2*eps);
  IndicatorCuboid2D<T> lid(extend,origin);

  superGeometry.rename(2,3,1,lid);

  const bool verbose = false;
  superGeometry.clean(verbose);
  superGeometry.innerClean(verbose);
  superGeometry.checkForErrors(verbose);
}

void prepareLattice(SuperLattice<T,DESCRIPTOR>& superLattice,
                    SuperGeometry<T,2>& superGeometry,
                    UnitConverter<T,DESCRIPTOR> const& converter)
{
  const T omega = converter.getLatticeRelaxationFrequency();

  /// Material=0 -->do nothing
  superLattice.defineDynamics<NoDynamics>(superGeometry, 0);

  /// Material=1 -->bulk dynamics
  superLattice.defineDynamics<BulkDynamics>(superGeometry, 1);

#ifdef LID_DRIVEN
  boundary::set<boundary::LocalVelocity<T,DESCRIPTOR,BulkDynamics>>(superLattice, superGeometry, 2);
  boundary::set<boundary::LocalVelocity<T,DESCRIPTOR,BulkDynamics>>(superLattice, superGeometry, 3);
#else
  superLattice.defineDynamics<BounceBack>(superGeometry, 2);
  superLattice.defineDynamics<BounceBack>(superGeometry, 3);
#endif

  superLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
}

void setBoundaryValues(SuperLattice<T,DESCRIPTOR>& superLattice,
                       SuperGeometry<T,2>& superGeometry,
                       UnitConverter<T,DESCRIPTOR> const& converter,
                       int iT)
{
  if (iT==0) {
    AnalyticalConst2D<T,T> rhoF(1);
    Vector<T,2> velocity{};
    AnalyticalConst2D<T,T> uF(velocity);

    superLattice.iniEquilibrium(superGeometry, 1, rhoF, uF);
    superLattice.iniEquilibrium(superGeometry, 2, rhoF, uF);
    superLattice.iniEquilibrium(superGeometry, 3, rhoF, uF);

    superLattice.defineRhoU(superGeometry, 1, rhoF, uF);
    superLattice.defineRhoU(superGeometry, 2, rhoF, uF);
    superLattice.defineRhoU(superGeometry, 3, rhoF, uF);

#ifdef LID_DRIVEN
    velocity[0] = converter.getCharLatticeVelocity();
    AnalyticalConst2D<T,T> u(velocity);
    superLattice.defineU(superGeometry,3,u);
#endif

    superLattice.initialize();
  }
}

void getResults(SuperLattice<T,DESCRIPTOR>& superLattice,
                SuperGeometry<T,2>& superGeometry,
                UnitConverter<T,DESCRIPTOR> const& converter)
{
  SuperVTMwriter2D<T> vtmWriter("cavity2d");

  SuperGeometryF<T,DESCRIPTOR> geometryF(superGeometry);
  SuperLatticeCuboid2D<T,DESCRIPTOR> cuboidF(superLattice);
  SuperLatticeRank2D<T,DESCRIPTOR> rankF(superLattice);
  SuperLatticePhysVelocity2D<T,DESCRIPTOR> velocityF(superLattice, converter);
  SuperLatticePhysPressure2D<T,DESCRIPTOR> pressureF(superLattice, converter);

  vtmWriter.write(geometryF);
  vtmWriter.write(cuboidF);
  vtmWriter.write(rankF);
  vtmWriter.write(velocityF);
  vtmWriter.write(pressureF);
}

int main(int argc, char **argv)
{
  initialize(&argc, &argv, false, false);

  singleton::directories().setOlbDir("../../../");
  singleton::directories().setOutputDir("./tmp/");

  std::size_t size  = argc >= 2 ? atoi(argv[1]) : 100;
  std::size_t steps = argc >= 3 ? atoi(argv[2]) : 10;
  std::size_t cuboids_per_process = argc == 4 ? atoi(argv[3]) : 1;

  UnitConverterFromResolutionAndLatticeVelocity<T,DESCRIPTOR> const converter(
    size,      // resolution: number of voxels per charPhysL
    0.01,      // lattice velocity
    1.0,       // charPhysLength: reference length of simulation geometry
    1.0,       // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    1./1000.,  // physViscosity: physical kinematic viscosity in __m^2 / s__
    1.0        // physDensity: physical density in __kg / m^3__
  );

  Vector<T,2> origin{};
  Vector<T,2> extend(converter.getCharPhysLength());
  IndicatorCuboid2D<T> cube(extend, origin);

#ifdef PARALLEL_MODE_MPI
  CuboidDecomposition2D<T> cuboidDecomposition(cube, converter.getConversionFactorLength(), cuboids_per_process*singleton::mpi().getSize());
#else
  CuboidDecomposition2D<T> cuboidDecomposition(cube, converter.getConversionFactorLength(), cuboids_per_process);
#endif

  BlockLoadBalancer<T> loadBalancer(singleton::mpi().getRank(), singleton::mpi().getSize(), cuboidDecomposition.size(), 0);

  SuperGeometry<T,2> superGeometry(cuboidDecomposition, loadBalancer, 2);

  prepareGeometry(converter, cube, superGeometry);

  SuperLattice<T,DESCRIPTOR> superLattice(converter, superGeometry);
  superLattice.setStatisticsOff();

  prepareLattice(superLattice, superGeometry, converter);
  setBoundaryValues(superLattice, superGeometry, converter, 0);

  getResults(superLattice, superGeometry, converter);

  util::Timer<T> timer(steps, superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < steps; ++iT) {
    superLattice.collideAndStream();
  }

  timer.stop();
  timer.update(steps);

  if (singleton::mpi().isMainProcessor()) {
    std::cout << size << ", "
              << steps << ", "
              << cuboids_per_process << ", "
              << singleton::mpi().getSize() << ", "
#ifdef PARALLEL_MODE_OMP
              << singleton::omp().getSize() << ", "
#else
              << 1 << ", "
#endif
              << timer.getTotalMLUPs() << std::endl;
  }

  getResults(superLattice, superGeometry, converter);

  return 0;
}
