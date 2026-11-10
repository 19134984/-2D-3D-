/* Lattice Boltzmann sample, written in C++, using the OpenLB
 * library
 *
 * Copyright (C) 2019-2025 Mathias J. Krause, Julius Jeßberger,
 *                         Shota Ito
 * E-mail contact: info@openlb.net
 * The most recent release of OpenLB can be downloaded at
 * <http:  //www.openlb.net/>
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#include <olb.h>

using namespace olb;

using T = FLOATING_POINT_TYPE;
using DESCRIPTOR = descriptors::D3Q19<>;
using ADDESCRIPTOR = descriptors::D3Q7<>;

const T latticeU = 0.1;
const T charU = 0.0996;   // charU
const T charNu = 1.e-6;   // charNu
const T charL = 0.00133;  // hydraulic diameter = 4 * surface / perimeter
const T charRho = 1000;   // charRhoFluid
const int N = 20;         // resolution of the hydraulic diameter
const T physDeltaX = charL/N;
const T physDeltaT = 400*6.67670682730923725550428726727e-5/(N*N);

// TODO this value has to be adapted to real values
// todo realistic value 1.e-9
const T diffusion = 30e-9;
const T physTotalTime = 8;   // time for fluid simulation
const T physStartTime = .6;    // time to start fluid pulsation
const int latticeOverlap = 2;
const T stlSize = 1.0;

void prepareGeometry(UnitConverter<T, DESCRIPTOR> const& converterNS,
                     IndicatorF3D<T>& indicator, SuperGeometry<T,3>& superGeometry)
{

  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0, 2, indicator);
  superGeometry.rename(2, 1, {1, 1, 1});

  // Returns the minimum phys position in each direction for material 2
  Vector<T,3> minR = superGeometry.getStatistics().getMinPhysR(2);
  Vector<T,3> maxR = superGeometry.getStatistics().getMaxPhysR(2);
  Vector<T,3> centerR = superGeometry.getStatistics().getCenterPhysR(2);

  // sets circle of both inflows and the outflow, with direction and radius
  IndicatorCircle3D<T> inflow1(minR[0], minR[1], centerR[2], 0., -1., 0.,
                               (maxR[0] - minR[0]) / 3.);
  IndicatorCircle3D<T> inflow2(maxR[0], minR[1], centerR[2], 0., -1., 0.,
                               (maxR[0] - minR[0]) / 3.);
  IndicatorCircle3D<T> outflow((maxR[0] + minR[0]) / 2., maxR[1], centerR[2],
                               0., 1., 0., (maxR[0] - minR[0]) / 3.);

  // sets cylinder on that in-/out-flow circles with length
  IndicatorCylinder3D<T> layerInflow1(inflow1, converterNS.getConversionFactorLength());
  IndicatorCylinder3D<T> layerInflow2(inflow2, converterNS.getConversionFactorLength());
  IndicatorCylinder3D<T> layerOutflow(outflow, converterNS.getConversionFactorLength());
  // renames all boundary voxels of material fromBcMat to toBcMat if two neighbour voxel
  // in the direction of the discrete normal are fluid voxel with material fluidM in the region
  // where the indicator function is fulfilled
  superGeometry.rename(2, 3, 1, layerInflow1); // layer of inflow1 gets mat = 3
  superGeometry.rename(2, 4, 1, layerInflow2); // layer of inflow2 gets mat = 4
  superGeometry.rename(2, 5, 1, layerOutflow); // layer of outflow gets mat = 5

  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(
  SuperGeometry<T,3>& superGeometry,
  SuperLattice<T, DESCRIPTOR>& sLattice,
  const T omega,
  SuperLattice<T, ADDESCRIPTOR>& sLatticeAD,
  T omegaAD)
{

  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  auto bulkIndicator = superGeometry.getMaterialIndicator({1,3,4,5});
  // dynamics for fluid
  sLattice.defineDynamics<BGKdynamics>(bulkIndicator);

  // boundary conditions for fluid
  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry, 4);
  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry, 5);


  // dynamics for adsorptive
  auto bulkIndicatorAD = superGeometry.getMaterialIndicator({1,3,5});
  dynamics::set<ParticleAdvectionDiffusionBGKdynamics>(sLatticeAD, bulkIndicatorAD);

  // boundary for adsorptive
  boundary::set<boundary::BounceBack>(sLatticeAD, superGeometry, 2);
  boundary::set<boundary::BounceBack>(sLatticeAD, superGeometry, 4);
  boundary::set<boundary::AdvectionDiffusionDirichlet>(sLatticeAD, superGeometry, 3);
  setZeroGradientBoundary<T,ADDESCRIPTOR>(sLatticeAD, superGeometry.getMaterialIndicator({5}));

  // initialisation for adsorptive
  auto initIndicatorAD = superGeometry.getMaterialIndicator({1,2,4,5});
  AnalyticalConst3D<T, T> u0(0.);
  AnalyticalConst3D<T, T> rhoSmall(1.e-8);
  sLatticeAD.iniEquilibrium(initIndicatorAD, rhoSmall, u0);

  sLattice.setParameter<descriptors::OMEGA>(omega);
  sLatticeAD.setParameter<descriptors::OMEGA>(omegaAD);

  // Make the lattice ready for simulation
  sLattice.initialize();
  sLatticeAD.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
  return;
}

void setBoundaryValues(
  SuperGeometry<T,3>& superGeometry,
  SuperLattice<T, DESCRIPTOR>& sLattice,
  UnitConverter<T, DESCRIPTOR> const& converterNS,
  int iT)
{
  OstreamManager clout(std::cout, "setBoundaryValues");

  std::vector < T > maxVelocity(3, T());
  const T distanceToBoundary = converterNS.getConversionFactorLength() / 2.;
  const T velNS = converterNS.getCharPhysVelocity();
  const int itStartTime = converterNS.getLatticeTime(physStartTime);

  if (iT <= itStartTime && iT % 50 == 0) {

    SinusStartScale<T, int> startScale(itStartTime, T(1)); //1+-amplitude
    int help[1] = { iT };
    T frac[3] = { T() };
    startScale(frac, help);

    // set lattice velocity on boundary
    maxVelocity[1] = velNS * frac[0];

    RectanglePoiseuille3D<T> u5(superGeometry, 5, maxVelocity,
                                distanceToBoundary, distanceToBoundary, distanceToBoundary);
    momenta::setVelocity(sLattice, superGeometry.getMaterialIndicator({5}), u5);
  }
}

void getResults(
  SuperLattice<T, DESCRIPTOR>& sLattice,
  UnitConverter<T, DESCRIPTOR>& converterNS,
  SuperLattice<T, ADDESCRIPTOR>& sLatticeAD,
  size_t iT, SuperGeometry<T,3>& superGeometry,
  util::Timer<double>& timer)
{
  OstreamManager clout(std::cout, "getResults");

  const int itTotalTime = converterNS.getLatticeTime(physTotalTime);
  const int iTperiodConsole = itTotalTime / 100;   // Console output
  const int iTperiodVTK = itTotalTime / 200;   // Writes the vtk files
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocityNS(sLattice, converterNS);
  SuperLatticeDensity3D<T, ADDESCRIPTOR> adsorptive(sLatticeAD);

  SuperVTMwriter3D<T> vtmWriter("microMixer3d");
  vtmWriter.addFunctor(velocityNS);
  vtmWriter.addFunctor(adsorptive);

  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  // console output
  if (iT % iTperiodConsole == 0) {
    timer.update(iT);
    timer.printStep();

    // output for latticeStatistics
    sLattice.getStatistics().print(iT, converterNS.getPhysTime(iT));
    sLatticeAD.getStatistics().print(iT, converterNS.getPhysTime(iT));
  }

  // vtk and gif output
  if (iT % iTperiodVTK == 0) {
    vtmWriter.write(iT);
  }
}


int main(int argc, char* argv[])
{
  /// === 1st Step: Initialization ===
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");

  UnitConverter<T,DESCRIPTOR> converterNS(
    ( T )   physDeltaX,     //physDeltaX
    ( T )   physDeltaT,     //physDeltaT
    ( T )   charL,          //charPhysLength
    ( T )   charU,          //charPhysVelocity
    ( T )   charNu,         //physViscosity
    ( T )   charRho         //physDensity
  );
  clout << "navier-stokes converter " << std::endl;
  converterNS.print();

  UnitConverter<T,ADDESCRIPTOR> converterAD(
    ( T )   physDeltaX,     //physDeltaX
    ( T )   physDeltaT,     //physDeltaT
    ( T )   charL,          //charPhysLength
    ( T )   charU,          //charPhysVelocity
    ( T )   diffusion,      //physViscosity
    ( T )   charRho         //physDensity
  );
  clout << "advection-diffusion converter " << std::endl;
  converterAD.print();

  /// compute relaxation parameter to solve the advection-diffusion
  /// equation in the lattice Boltzmann context
  T omegaAD = converterAD.getLatticeRelaxationFrequency();
  clout << "diffusion = " << diffusion << std::endl;
  clout << "omegaAD = " << omegaAD << std::endl;
  clout << "tauAD = " << 1./omegaAD << std::endl;
  /// === 2nd Step: Prepare Geometry ===

  // Instantiation of an empty cuboidDecomposition
  const int noOfCuboids = util::max(16, 4 * singleton::mpi().getSize());

  STLreader<T> stlReader("microMixer3d.stl",
                         converterNS.getConversionFactorLength(),
                         stlSize);

  IndicatorLayer3D<T> extendedDomain( stlReader,
                                      converterNS.getConversionFactorLength() );

  CuboidDecomposition3D<T> cuboidDecomposition( extendedDomain,
                                                converterNS.getConversionFactorLength(),
                                                noOfCuboids );

  // Instantiation of an empty loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  // Default instantiation of superGeometry
  SuperGeometry<T,3> superGeometry(cuboidDecomposition, loadBalancer, latticeOverlap);
  prepareGeometry(converterNS, stlReader, superGeometry);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice(converterNS, superGeometry);
  SuperLattice<T, ADDESCRIPTOR> sLatticeAD(converterAD, superGeometry);
  prepareLattice(superGeometry, sLattice, converterNS.getLatticeRelaxationFrequency(),
                 sLatticeAD, omegaAD);

  SuperLatticeCoupling coupling(
    NavierStokesAdvectionDiffusionCoupling{},
    names::NavierStokes{}, sLattice,
    names::Temperature{}, sLatticeAD
  );
  coupling.restrictTo(superGeometry.getMaterialIndicator({1}));

  /// === 4th Step: Fluid Main Loop with Timer ===
  util::Timer<double> timer(converterNS.getLatticeTime(physTotalTime),
                      superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (size_t iT = 0; iT < converterNS.getLatticeTime(physTotalTime); ++iT) {

    setBoundaryValues(superGeometry, sLattice, converterNS, iT);

    sLattice.collideAndStream();
    coupling.apply();
    sLatticeAD.collideAndStream();

    getResults(sLattice, converterNS, sLatticeAD, iT, superGeometry, timer);
  }

  timer.stop();
  timer.printSummary();
}
