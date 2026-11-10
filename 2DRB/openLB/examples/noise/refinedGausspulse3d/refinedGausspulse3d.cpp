/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Philipp Spelten
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

#include "olb.h"
#include "../gausspulse3d/gausspulse3d.h"
#include "refinedGausspulse3d.h"

using namespace olb;

void simulateRefined( MyCase& myCase ) {
  // Here, lattice and geometry are recreated, independend of the mesh in myCase
  OstreamManager clout(std::cout, "simulateRefined");
  clout << "simulateRefined ..." << std::endl;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto&         parameters    = myCase.getParameters();
  auto&         converterL0   = myCase.getLattice(names::NavierStokes{}).getUnitConverter();
  const size_t  iTgraph       = parameters.get<parameters::LATTICE_VTK_ITER_T>();
  const size_t  iTlog         = std::min(iTgraph, (std::size_t) parameters.get<parameters::LATTICE_STAT_ITER_T>());
  const bool    doGraphicalOutput = parameters.get<parameters::DO_GRAPHICAL_OUTPUT>();
  const bool    doL2Plot      = parameters.get<parameters::DO_L2_PLOT>();
  const Vector  extentDomain  = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector  originDomain  = parameters.get<parameters::DOMAIN_ORIGIN>();
  const Vector  extentFluid   = parameters.get<parameters::CORE_EXTENT>();
  const Vector  originFluid   = -0.5*extentFluid;
  const T       physDeltaX    = parameters.get<parameters::PHYS_DELTA_X>();
  const T       spongeDepthPU = physDeltaX * parameters.get<parameters::SPONGE_DEPTH_LU>();
  std::string   suffix        = parameters.get<parameters::OUTDIR_SUFFIX>();

  IndicatorCuboid3D<T> domainI( extentDomain, originDomain );
  IndicatorCuboid3D<T> fluidI(  extentFluid, originFluid );

  /// === Coarse Cuboid ===
  CuboidDecomposition3D<T> cuboidDecompositionL0( domainI, physDeltaX );
  // CAREFUL: Periodicity is not set by default. This causes issues after refinement
  cuboidDecompositionL0.setPeriodicity( {true, true, true} );
  std::vector<T> fractions(3);
  fractions[0] = spongeDepthPU/extentDomain[0]*.9;
  fractions[1] = (extentDomain[0]-2*spongeDepthPU)*1.2/extentDomain[0];
  fractions[2] = spongeDepthPU/extentDomain[0]*.9;
  clout << "Splitting level 0 into fractions " << fractions << std::endl;
  cuboidDecompositionL0.splitFractional(0, 0, fractions);
  cuboidDecompositionL0.splitFractional(1, 1, fractions);
  cuboidDecompositionL0.splitFractional(3, 2, fractions);
  cuboidDecompositionL0.print();
  /// === Coarse Cuboid ===

  /// === Refined Cuboid ===
  auto cuboidDecompositionL1 = cuboidDecompositionL0;
  cuboidDecompositionL1.remove( fluidI );
  cuboidDecompositionL1.refine(2);
  cuboidDecompositionL1.print();
  /// === Refined Cuboid ===

  // Adjust weights for balancing
  for (int iC=0; iC < cuboidDecompositionL0.size(); ++iC) {
    auto& cuboid = cuboidDecompositionL0.get(iC);
    cuboid.setWeight(cuboid.getLatticeVolume());
    auto origin = cuboidDecompositionL0.get(iC).getOrigin();
    for (int jC=0; jC < cuboidDecompositionL1.size(); ++jC) {
      if (cuboidDecompositionL1.get(jC).getOrigin() == origin) {
        cuboid.setWeight(cuboid.getWeight() + 2*cuboidDecompositionL1.get(jC).getLatticeVolume());
      }
    }
  }

  // === Coarse geometry and lattice ===
  HeuristicLoadBalancer<T> loadBalancerL0( cuboidDecompositionL0 );
  SuperGeometry<T,3> geometryL0( cuboidDecompositionL0, loadBalancerL0 );
  geometryL0.setWriteIncrementalVTK( false );
  geometryL0.rename(0,1);
  geometryL0.clean();
  geometryL0.checkForErrors();
  geometryL0.print();
  SuperLattice<T,DESCRIPTOR> latticeL0( converterL0, cuboidDecompositionL0, loadBalancerL0, 3 );
  latticeL0.defineDynamics<BGKdynamics>( geometryL0, 1 );
  latticeL0.setParameter<descriptors::OMEGA>( converterL0.getLatticeRelaxationFrequency() );
  latticeL0.initialize();
  // === Coarse geometry and lattice ===

  // === Fine geometry and lattice ===
  RefinedLoadBalancer<T,3> loadBalancerL1( cuboidDecompositionL0, loadBalancerL0, cuboidDecompositionL1);
  auto converterL1 = convectivelyRefineUnitConverter( converterL0, 2 );
  converterL1.print();
  SuperGeometry<T,3> geometryL1( cuboidDecompositionL1, loadBalancerL1 );
  geometryL1.setWriteIncrementalVTK( false );
  geometryL1.rename(0,1);
  geometryL1.checkForErrors();
  geometryL1.print();
  SuperLattice<T,DESCRIPTOR> latticeL1( converterL1, cuboidDecompositionL1, loadBalancerL1, 3 );
  latticeL1.defineDynamics<BGKdynamics>( geometryL1, 1 );
  latticeL1.setParameter<descriptors::OMEGA>( converterL1.getLatticeRelaxationFrequency() );
  latticeL1.initialize();
  // === Fine geometry and lattice ===

  // === Coarse initial values ===
  const Vector uInfty = parameters.get<parameters::PULSE_PHYS_VELOCITY>();
  AnalyticalConst3D<T,T> u( converterL0.getLatticeVelocity(uInfty[0]),
                            converterL0.getLatticeVelocity(uInfty[1]),
                            converterL0.getLatticeVelocity(uInfty[2]) );
  AnalyticalConst3D<T,T> rho0( 1. );
  latticeL0.defineRhoU( geometryL0.getMaterialIndicator( 1 ), rho0, u);
  latticeL0.iniEquilibrium( geometryL0.getMaterialIndicator( 1 ), rho0, u);
  latticeL0.initialize();
  // === Coarse initial values ===

  // === Fine initial values ===
  setFineInitialValues( parameters, latticeL1, geometryL1 );
  // === Fine initial values ===

  getInitialSetup( latticeL0, geometryL0, suffix+"_farField" );
  getInitialSetup( latticeL1, geometryL1, suffix+"_core" );

  // === Rohde coupling ===
  clout << "Rohde coupling ..." << std::endl;
  auto coupler = refinement::rohde::makeCoupler(
    latticeL0, geometryL0,
    latticeL1, geometryL1,
    geometryL1.getMaterialIndicator( 1 )
  );
  clout << "Rohde coupling ... OK" << std::endl;
  // === Rohde coupling ===

  // === Initialize pressure L2 norm plot
  Gnuplot<T> gplot_L0("l2_absolute");
  gplot_L0.setLabel("time []", "absolute L2 norm []");
  SuperIndicatorFfromIndicatorF3D<T> domainFluidL0( new IndicatorCuboid3D<T>( extentFluid, originFluid ), geometryL0 );
  T Lp0_L0 = L2Norm<3, T, DESCRIPTOR>( latticeL0, converterL0, domainFluidL0 );
  Gnuplot<T> gplot_L1("l2_absolute_fine");
  gplot_L1.setLabel("time []", "absolute L2 norm []");
  SuperIndicatorFfromIndicatorF3D<T> domainFluidL1( new IndicatorCuboid3D<T>( extentFluid, originFluid ), geometryL1 );
  T Lp0_L1 = L2Norm<3,T,DESCRIPTOR>( latticeL1, latticeL1.getUnitConverter(), domainFluidL1 );
  // === Initialize pressure L2 norm plot

  // === Timer setup ===
  size_t iTmaxT = converterL0.getLatticeTime(parameters.get<parameters::PHYS_MAX_T>());
  size_t iTmax  = parameters.get<parameters::LATTICE_MAX_T>();
  iTmax = ( iTmax > 0 ) ? std::min(iTmax+1, iTmaxT) : iTmaxT;
  util::Timer<T> timer(iTmax, geometryL0.getStatistics().getNvoxel()
                              + 2*geometryL1.getStatistics().getNvoxel());
  timer.start();
  // === Timer setup ===

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    /// === Computation and Output of the Results ===
    if ( doL2Plot ) {
      setPlotData( myCase, iT, gplot_L0, Lp0_L0 );
      setFinePlotData( myCase, latticeL1, geometryL1, iT*2, gplot_L1, Lp0_L1 );
    }
    if ( doGraphicalOutput ) {
      getGraphicalResults( myCase, iT );
      getFineGraphicalResults( myCase, iT*2, latticeL1 );
    }

    /// === Collide and Stream Execution ===
    latticeL0.setProcessingContext(ProcessingContext::Simulation);
    latticeL1.setProcessingContext(ProcessingContext::Simulation);
    latticeL0.collide();
      coupler->apply(meta::id<refinement::rohde::CoarseToFineO>{});
      latticeL1.collideAndStream();
      latticeL1.collideAndStream();
    latticeL0.stream();
    coupler->apply(meta::id<refinement::rohde::FineToCoarseO>{});

    /// === Print some (numerical and computational) statistics ===
    if ( iT%iTlog == 0 ) {
      latticeL1.getStatistics().print(iT*2, converterL0.getPhysTime(iT));
      timer.print(iT);
    }

    singleton::pool().wait();

    gplot_L0.setYrange(1e-5, 1);
    gplot_L0.setLogScale(2);
    gplot_L0.writePNG(-1, -1, "gplot_l2_abs");
    gplot_L1.setYrange(1e-5, 1);
    gplot_L1.setLogScale(2);
    gplot_L1.writePNG(-1, -1, "gplot_l2_abs");
  }

  timer.stop();
  timer.printSummary();
  clout << "simulateRefined ... OK" << std::endl;
}

int main(int argc, char* argv[])
{
  // === Step 1: Initialization ===
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");
  clout << "Starting gausspulse3d ..." << std::endl;

  /// === Step 2a: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  setGetParameters(myCaseParameters, argc, argv);

  // === Step 2b: Set coarse-specific parameters ===
  myCaseParameters.set<parameters::FAR_FIELD_TYPE>( COARSE );
  myCaseParameters.set<parameters::RESOLUTION>( int(myCaseParameters.get<parameters::RESOLUTION>() / 2. ) );
  myCaseParameters.set<parameters::LATTICE_VTK_ITER_T>( int(myCaseParameters.get<parameters::LATTICE_VTK_ITER_T>() / 2. ) + 1 );
  myCaseParameters.set<parameters::LATTICE_MAX_T>( myCaseParameters.get<parameters::LATTICE_MAX_T>() / 2 );

  // === Step 2b: set output directory depending on input values ===
  setOutDir(myCaseParameters);

  /// === Step 3: Create Coarse Mesh ===
  Mesh<MyCase::value_t,MyCase::d> mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Coarse Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Coarse Lattice ===
  prepareLattice(myCase);

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulateRefined(myCase);
}
