/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Fedor Bukreev, Adrian Kummerländer, Mathias J. Krause
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

 /* laminarReactiveTmixer.cpp:
 * Benchmark Simulation of Laminar Reactive Micromixing Using Lattice Boltzmann Methods
 * Fedor Bukreev, Adrian Kummerländer, Julius Jeßberger, Dennis Teutscher, Stephan Simonis,
 * Dieter Bothe, and Mathias J. Krause, AIAA Journal, https://doi.org/10.2514/1.J064234
 *
 * This example is a new benchmark for laminar reactive micromixing in a T-shaped mixer
 * with A+B->C reaction and special Smagorinsky-like stabilization approach.
 */


#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::descriptors;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q19<>>,
  Concentration<0>,  Lattice<double, descriptors::D3Q7<VELOCITY, SOURCE, OMEGA>>,
  Concentration<1>,  Lattice<double, descriptors::D3Q7<VELOCITY, SOURCE, OMEGA>>,
  Concentration<2>,  Lattice<double, descriptors::D3Q7<VELOCITY, SOURCE, OMEGA>>
>;

namespace olb::parameters {

struct LX_INLET  : public descriptors::FIELD_BASE<1> { };
struct LY_OUTLET : public descriptors::FIELD_BASE<1> { };
struct PLOT_OVER_LINE : public TYPED_FIELD_BASE<bool,1> { };
struct CONCENTRATION : public descriptors::FIELD_BASE<6> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T lx_inlet = params.get<parameters::LX_INLET>();
  const T ly_outlet = params.get<parameters::LY_OUTLET>();

  const T hydrInlet  = T(2)*ly_outlet*extent[2]/(ly_outlet + extent[2]);
  const T physDeltaX = hydrInlet / params.get<parameters::RESOLUTION>();

  Vector<T,3> extendReactor( extent[0] + T(2)*physDeltaX, ly_outlet + T(2)*physDeltaX, extent[2] + T(2)*physDeltaX );
  Vector<T,3> originReactor( -physDeltaX, -(ly_outlet + T(2)*physDeltaX)/T(2), -(extent[2] + T(2)*physDeltaX)/T(2) );
  std::shared_ptr<IndicatorCuboid3D<T> > reactor( new IndicatorCuboid3D<T> ( extendReactor, originReactor ));

  Vector<T,3> extendInflow( lx_inlet + T(2)*physDeltaX, extent[1] + T(2)*physDeltaX, extent[2] + T(2)*physDeltaX );
  Vector<T,3> originInflow( -physDeltaX, -(extent[1] + T(2)*physDeltaX)/T(2), -(extent[2] + T(2)*physDeltaX)/T(2) );
  std::shared_ptr<IndicatorCuboid3D<T> > inflow( new IndicatorCuboid3D<T> ( extendInflow, originInflow ));

  IndicatorIdentity3D<T> mixer(reactor + inflow);

  Mesh<T,MyCase::d> mesh(mixer, physDeltaX, singleton::mpi().getSize(), "volume");
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  return mesh;
}

// Stores geometry information in form of material numbers
void prepareGeometry(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T lx_inlet = params.get<parameters::LX_INLET>();
  const T ly_outlet = params.get<parameters::LY_OUTLET>();

  const T hydrInlet  = T(2)*ly_outlet*extent[2]/(ly_outlet + extent[2]);
  const T physDeltaX = hydrInlet / params.get<parameters::RESOLUTION>();

  geometry.rename( 0,2 );

  Vector<T,3> extendReactor( extent[0], ly_outlet, extent[2] );
  Vector<T,3> originReactor( 0., -ly_outlet/T(2), -extent[2]/T(2) );
  IndicatorCuboid3D<T> reactor( extendReactor, originReactor );

  Vector<T,3> extendInflow( lx_inlet, extent[1], extent[2] );
  Vector<T,3> originInflow( 0., -extent[1]/T(2), -extent[2]/T(2) );
  IndicatorCuboid3D<T> inflow( extendInflow, originInflow );

  geometry.rename( 2, 1, reactor );
  geometry.rename( 2, 1, inflow );

  Vector<T,3> extendInlet1( lx_inlet, physDeltaX, extent[2]);
  Vector<T,3> originInlet1( 0., -extent[1]/T(2), -extent[2]/T(2) );
  IndicatorCuboid3D<T> inlet1( extendInlet1, originInlet1 );
  geometry.rename( 1, 3, inlet1 );

  Vector<T,3> extendInlet2( lx_inlet, physDeltaX, extent[2]);
  Vector<T,3> originInlet2( 0., extent[1]/T(2)-physDeltaX, -extent[2]/T(2) );
  IndicatorCuboid3D<T> inlet2( extendInlet2, originInlet2 );
  geometry.rename( 1, 4, inlet2 );

  Vector<T,3> extendOutlet( physDeltaX, ly_outlet, extent[2] );
  Vector<T,3> originOutlet( extent[0]-physDeltaX, -ly_outlet/T(2), -extent[2]/T(2) );
  IndicatorCuboid3D<T> outlet( extendOutlet, originOutlet );
  geometry.rename( 1, 5, outlet );

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  // Removes all not needed boundary voxels inside the surface
  geometry.innerClean();
  geometry.checkForErrors();
  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the NS lattice of the simulation
void prepareLatticeNS(MyCase& myCase) {
  OstreamManager clout(std::cout,"prepareLattice");

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& lattice = myCase.getLattice(NavierStokes{});

  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using RADDESCRIPTOR = MyCase::descriptor_t_of<Concentration<0>>;

  const int N = params.get<parameters::RESOLUTION>();
  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  const T ly_outlet = params.get<parameters::LY_OUTLET>();
  const T hydrInlet  = T(2)*ly_outlet*extent[2]/(ly_outlet + extent[2]);
  const T physCharVelocity = params.get<parameters::PHYS_CHAR_VELOCITY>();
  const T latticeCharVelocity = params.get<parameters::LATTICE_CHAR_VELOCITY>();
  const T physCharViscosity = params.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T physCharDensity = params.get<parameters::PHYS_CHAR_DENSITY>();

  lattice.setUnitConverter<UnitConverterFromResolutionAndLatticeVelocity<T,DESCRIPTOR>>(
    (T)   N,                   // resolution
    (T)   latticeCharVelocity, // maximal lattice velocity
    (T)   hydrInlet,           // charPhysLength: reference length of simulation geometry
    (T)   physCharVelocity,    // charPhysVelocity
    (T)   physCharViscosity,   // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   physCharDensity      // physDensity: physical density in __kg / m^3__
  );
  const auto& converter = lattice.getUnitConverter();
  converter.print();

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  // Material=3 -->bulk dynamics (inlet1 (-Y))
  // Material=4 -->bulk dynamics (inlet2 (+Y))
  // Material=5 -->bulk dynamics (outlet)
  auto bulkIndicator = geometry.getMaterialIndicator({1, 3, 4, 5});
  dynamics::set<BGKdynamics>(lattice, bulkIndicator);

  // Material=2 -->bounce back
  dynamics::set<BounceBack>(lattice, geometry.getMaterialIndicator(2));

  //if interpolated boundary conditions are chosen
  boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 3);
  boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 4);
  boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 5);

  lattice.setParameter<descriptors::OMEGA>(omega);

  auto DIFFUSION = params.get<parameters::DIFFUSION_COMPONENTS>();

  T omegaA = converter.getLatticeRelaxationFrequencyFromDiffusivity<RADDESCRIPTOR>( DIFFUSION[0] );
  T omegaB = converter.getLatticeRelaxationFrequencyFromDiffusivity<RADDESCRIPTOR>( DIFFUSION[1] );
  T omegaC = converter.getLatticeRelaxationFrequencyFromDiffusivity<RADDESCRIPTOR>( DIFFUSION[2] );

  auto& coupling = myCase.setCouplingOperator(
    "LESReaction",
    LESReactionCoupling<T,3>{},
    names::NavierStokes{}, lattice,
    names::Concentration0{}, myCase.getLattice(Concentration<0>{}),
    names::Concentration1{}, myCase.getLattice(Concentration<1>{}),
    names::Concentration2{}, myCase.getLattice(Concentration<2>{})
  );
  coupling.setParameter<LESReactionCoupling<T,3>::LATTICE_REACTION_COEFF>(converter.getConversionFactorTime()*1.e5);
  coupling.setParameter<LESReactionCoupling<T,3>::STOCH_COEFF>({-1, -1., 1.});
  coupling.setParameter<LESReactionCoupling<T,3>::REACTION_ORDER>({1., 1., 0.});
  coupling.setParameter<LESReactionCoupling<T,3>::SMAGORINSKY_PREFACTOR>(0.2);
  coupling.setParameter<LESReactionCoupling<T,3>::SCHMIDT>({0.03, 0.03, 0.03});
  coupling.setParameter<LESReactionCoupling<T,3>::OMEGA_NSE>(1./converter.getLatticeRelaxationFrequency());
  coupling.setParameter<LESReactionCoupling<T,3>::OMEGAS_ADE>({omegaA, omegaB, omegaC});

  clout << "Prepare Lattice for NSE... OK" << std::endl;
}

// Set up the CRAD lattice of the simulation
template<int ID>
void prepareLatticeCRAD(MyCase& myCase) {
  OstreamManager clout(std::cout,"prepareLatticeCRAD");
  clout << "Prepare Lattice for Coupled Reaction simualtion..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& NSlattice = myCase.getLattice(NavierStokes{});
  auto& CRADlattice = myCase.getLattice(Concentration<ID>{});
  using RADDESCRIPTOR = MyCase::descriptor_t_of<Concentration<ID>>;

  const auto& NSconverter = NSlattice.getUnitConverter();
  CRADlattice.setUnitConverter(NSconverter);

  // bulk
  dynamics::set<SourcedLimitedAdvectionDiffusionBGKdynamics>(CRADlattice, geometry.getMaterialIndicator({1, 3, 4}));

  // Material=2 -->bounce back
  dynamics::set<BounceBack>(CRADlattice, geometry.getMaterialIndicator(2));

  // Setting of the boundary conditions, inflow and outflow with Dirichlet and Neumann BCs
  boundary::set<boundary::AdvectionDiffusionDirichlet>(CRADlattice, geometry.getMaterialIndicator({3}));
  boundary::set<boundary::AdvectionDiffusionDirichlet>(CRADlattice, geometry.getMaterialIndicator({4}));
  setZeroGradientBoundary<T, RADDESCRIPTOR>(CRADlattice, geometry.getMaterialIndicator({5}));

  auto DIFFUSION = params.get<parameters::DIFFUSION_COMPONENTS>();
  T omega = NSconverter.getLatticeRelaxationFrequencyFromDiffusivity<RADDESCRIPTOR>( DIFFUSION[ID] );
  AnalyticalConst3D<T,T> omegaD( omega );
  fields::set<descriptors::OMEGA>(CRADlattice, geometry.getMaterialIndicator({1, 2, 3, 4, 5}), omegaD);
  CRADlattice.template setParameter<descriptors::OMEGA>(omega);

  {
    auto& communicatorCR = CRADlattice.getCommunicator(stage::PreCollide());
    communicatorCR.requestOverlap(params.get<parameters::OVERLAP>());
    communicatorCR.exchangeRequests();
  }

  clout << "Prepare Lattice for Coupled Reaction simualtion... OK" << std::endl;
}

void setInitialValuesNS(MyCase& myCase) {
  using T = MyCase::value_t;

  auto& lattice = myCase.getLattice(NavierStokes{});

  // Make the lattice ready for simulation
  lattice.initialize();
}

template<int ID>
void setInitialValuesCRAD(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();
  auto& CRADlattice = myCase.getLattice(Concentration<ID>{});
  const auto& converter = myCase.getLattice(NavierStokes{}).getUnitConverter();

  auto rhos = params.get<parameters::CONCENTRATION>();
  const T physMaxVelocity = params.get<parameters::PHYS_MAX_VELOCITY>();

  AnalyticalConst3D<T,T> rho0(1.e-8);
  AnalyticalConst3D<T,T> rhoA(rhos[2*ID+0]);
  AnalyticalConst3D<T,T> rhoB(rhos[2*ID+1]);

  T distance2Wall = converter.getConversionFactorLength()/T(2);
  std::vector<T> maxVelocity1( 3,0 );
  std::vector<T> maxVelocity2( 3,0 );
  maxVelocity1[1] = converter.getLatticeVelocity(physMaxVelocity);
  maxVelocity2[1] = converter.getLatticeVelocity(physMaxVelocity);
  RectanglePoiseuille3D<T> poiseuilleU1( geometry, 3, maxVelocity1, distance2Wall, distance2Wall, distance2Wall );
  RectanglePoiseuille3D<T> poiseuilleU2( geometry, 4, maxVelocity2, distance2Wall, distance2Wall, distance2Wall );

  momenta::setDensity(CRADlattice, geometry.getMaterialIndicator({1, 2, 5}), rho0);
  momenta::setDensity(CRADlattice, geometry.getMaterialIndicator(3), rhos[2*ID+0]);
  momenta::setDensity(CRADlattice, geometry.getMaterialIndicator(4), rhos[2*ID+1]);
  momenta::setVelocity(CRADlattice, geometry.getMaterialIndicator(3), poiseuilleU1);
  momenta::setVelocity(CRADlattice, geometry.getMaterialIndicator(4), poiseuilleU2);

  CRADlattice.iniEquilibrium( geometry, 3, rhoA, poiseuilleU1 );
  CRADlattice.iniEquilibrium( geometry, 4, rhoB, poiseuilleU2 );

  fields::set<descriptors::SOURCE>(CRADlattice, geometry.getMaterialIndicator({1, 2, 3, 4, 5 }), rho0);

  CRADlattice.initialize();
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{
  OstreamManager clout( std::cout,"setTemporalValues" );
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});
  const auto& converter = lattice.getUnitConverter();

  // No of time steps for smooth start-up
  const T startTime = params.get<parameters::PHYS_START_T>();
  const int iTmaxStart = converter.getLatticeTime( startTime );
  const T physMaxVelocity = params.get<parameters::PHYS_MAX_VELOCITY>();

  if ( int(iT)%100 == 0 && int(iT) <= iTmaxStart ) {

    // Smooth start curve, sinus
    SinusStartScale<T,int> startScale(iTmaxStart, T(1));

    // Creates and sets the Poiseuille inflow profile using functors
    int iTvec[1]= {int(iT)};
    T frac[1]= {};
    startScale( frac,iTvec );
    std::vector<T> maxVelocity1( 3,0 );
    std::vector<T> maxVelocity2( 3,0 );
    maxVelocity1[1] = frac[0]*converter.getLatticeVelocity(physMaxVelocity);
    maxVelocity2[1] = -frac[0]*converter.getLatticeVelocity(physMaxVelocity);

    T distance2Wall = converter.getConversionFactorLength()/T(2);
    RectanglePoiseuille3D<T> poiseuilleU1( geometry, 3, maxVelocity1, distance2Wall, distance2Wall, distance2Wall );
    RectanglePoiseuille3D<T> poiseuilleU2( geometry, 4, maxVelocity2, distance2Wall, distance2Wall, distance2Wall );

    momenta::setVelocity(lattice, geometry.getMaterialIndicator(3), poiseuilleU1);
    momenta::setVelocity(lattice, geometry.getMaterialIndicator(4), poiseuilleU2);

    lattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

// Output to console and files
void getResultsNS(MyCase& myCase,
                  util::Timer<MyCase::value_t>& timer,
                  std::size_t iT)
{
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& params = myCase.getParameters();
  const auto& converter = lattice.getUnitConverter();
  OstreamManager clout( std::cout,"getResultsNS" );

  SuperVTMwriter3D<T> vtmWriterNS( "mixer_fluid" );
  SuperGeometryF<T, DESCRIPTOR::d> materials(geometry);
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( lattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( lattice, converter );
  vtmWriterNS.addFunctor( materials );
  vtmWriterNS.addFunctor( velocity );
  vtmWriterNS.addFunctor( pressure );

  const T maxPhysT = params.get<parameters::PHYS_MAX_T>();
  const int  vtkIter  = converter.getLatticeTime( maxPhysT/10 );
  const int  statIter = converter.getLatticeTime( maxPhysT/100 );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( lattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( lattice );
    vtmWriterNS.write( cuboid );
    vtmWriterNS.write( rank );
    vtmWriterNS.createMasterFile();
  }

  // Writes the ppm files
  if ( iT%vtkIter==0 ) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriterNS.write( iT );
  }

  // Writes output on the console
  if ( iT%statIter==0 && iT>=0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();
    lattice.getStatistics().print(iT,converter.getPhysTime(iT));
  }

  const T lx_inlet = params.get<parameters::LX_INLET>();
  const bool plotOverLine = params.get<parameters::PLOT_OVER_LINE>();
  if(plotOverLine) {
    if (iT == (converter.getLatticeTime( maxPhysT ) - 1)) {
      // Preparation
      lattice.setProcessingContext(ProcessingContext::Evaluation);
      Gnuplot<T> csvWriter("flowField");
      int dummy[4];
      T result[4];
      const int nX = geometry.getCuboidDecomposition().getMotherCuboid().getNx();  // number of voxels in x-direction
      const int minX = converter.getLatticeLength(lx_inlet);                       // beginning of the mixing channel
      auto mat = new SuperIndicatorIdentity3D<T> (geometry.getMaterialIndicator({1,5}));
      const Vector<T,3> normal( 1, 0, 0 );

      for (int iX = minX; iX < nX - 3; ++iX) {
        const T x = converter.getPhysLength(iX);
        const Vector<T,3> center( x, 0, 0 );
        auto circle = new SuperIndicatorFfromIndicatorF3D<T>(
          std::shared_ptr<IndicatorF3D<T>>(new IndicatorCircle3D<T> ( center, normal, 0.0003)), geometry);
        SuperIndicatorMultiplication3D<T> plane(circle, mat);

        // Average concentration: use average functor
        SuperAverage3D<T> (velocity, plane).operator()(result, dummy);
        const T avVel = result[0];
        SuperAverage3D<T> (pressure, plane).operator()(result, dummy);
        const T avPres = result[0];

        csvWriter.setData(x, {avVel, avPres});
      }
      csvWriter.writePNG();
    }
  }
}

void getResultsCRAD(MyCase& myCase,
                    util::Timer<MyCase::value_t>& timer,
                    std::size_t iT)
{
  using T = MyCase::value_t;
  using RADDESCRIPTOR = MyCase::descriptor_t_of<Concentration<0>>;
  auto& geometry = myCase.getGeometry();
  const auto& converter = myCase.getLattice(NavierStokes{}).getUnitConverter();
  const auto& converterAD = myCase.getLattice(Concentration<0>{}).getUnitConverter();
  auto& params = myCase.getParameters();
  OstreamManager clout( std::cout,"getResultsCRAD" );

  SuperVTMwriter3D<T> vtmWriterCRAD( "mixer_wholeReaction" );

  // insert concentrations and velocoties for each component into vector and vtk
  std::vector<SuperLatticeDensity3D<T, RADDESCRIPTOR>*> densities;
  std::vector<SuperLatticePhysViscosity3D<T, RADDESCRIPTOR>*> diffs;
  std::vector<std::string> names = {"A","B","P"};

  densities.emplace_back( new SuperLatticeDensity3D<T, RADDESCRIPTOR> (myCase.getLattice(Concentration<0>{})));
  densities.emplace_back( new SuperLatticeDensity3D<T, RADDESCRIPTOR> (myCase.getLattice(Concentration<1>{})));
  densities.emplace_back( new SuperLatticeDensity3D<T, RADDESCRIPTOR> (myCase.getLattice(Concentration<2>{})));
  for (int i = 0; i<3;i++){
    densities[i]->getName() = "Concentration " + names[i];
    vtmWriterCRAD.addFunctor(*densities[i]);
  }
  diffs.emplace_back( new SuperLatticePhysViscosity3D<T, RADDESCRIPTOR> (myCase.getLattice(Concentration<0>{}), converterAD));
  diffs.emplace_back( new SuperLatticePhysViscosity3D<T, RADDESCRIPTOR> (myCase.getLattice(Concentration<1>{}), converterAD));
  diffs.emplace_back( new SuperLatticePhysViscosity3D<T, RADDESCRIPTOR> (myCase.getLattice(Concentration<2>{}), converterAD));
  for (int i = 0; i<3;i++){
    diffs[i]->getName() = "Diffusisvity " + names[i];
    vtmWriterCRAD.addFunctor(*diffs[i]);
  }

  const T maxPhysT = params.get<parameters::PHYS_MAX_T>();
  const T maxPhysTAD = params.get<parameters::PHYS_MAX_T_AD>();
  const int  vtkIter  = converter.getLatticeTime( maxPhysTAD/10 );
  const int  statIter = converter.getLatticeTime( maxPhysTAD/100 );

  if ( iT==0 ) {
    vtmWriterCRAD.createMasterFile();
  }

  // Writes the ppm files
  if ( iT%vtkIter==0 ) {
    myCase.getLattice(Concentration<0>{}).setProcessingContext(ProcessingContext::Evaluation);
    myCase.getLattice(Concentration<1>{}).setProcessingContext(ProcessingContext::Evaluation);
    myCase.getLattice(Concentration<2>{}).setProcessingContext(ProcessingContext::Evaluation);

    SuperGeometryF<T, RADDESCRIPTOR::d> materials(geometry);
    vtmWriterCRAD.addFunctor( materials );
    vtmWriterCRAD.write( iT );
  }

  // Writes output on the console
  if ( iT%statIter==0 && iT>=0 ) {
    // Lattice statistics console output
    myCase.getLattice(Concentration<1>{}).getStatistics().print(iT,converter.getPhysTime(iT));
  }

  const bool plotOverLine = params.get<parameters::PLOT_OVER_LINE>();
  if(plotOverLine) {
    // Compute mixing quantities along central axis
    if (iT == (converter.getLatticeTime( maxPhysT ) - 1)) {
      myCase.getLattice(Concentration<0>{}).setProcessingContext(ProcessingContext::Evaluation);
      myCase.getLattice(Concentration<1>{}).setProcessingContext(ProcessingContext::Evaluation);
      myCase.getLattice(Concentration<2>{}).setProcessingContext(ProcessingContext::Evaluation);
      // Preparation
      Gnuplot<T> csvWriter("mixing");
      int dummy[4];
      T result[4];
      const int nX = geometry.getCuboidDecomposition().getMotherCuboid().getNx();  // number of voxels in x-direction
      const int minX = 3;                      // beginning of the mixing channel
      auto mat = new SuperIndicatorIdentity3D<T> (geometry.getMaterialIndicator({1/*,5*/}));
      const Vector<T,3> normal( 1, 0, 0 );

      for (int iX = minX; iX < nX - 3; ++iX) {
        const T x = converter.getPhysLength(iX);
        const Vector<T,3> center( x, 0, 0 );
        auto circle = new SuperIndicatorFfromIndicatorF3D<T>(
          std::shared_ptr<IndicatorF3D<T>>(new IndicatorCircle3D<T> ( center, normal, 0.0003)), geometry);
        SuperIndicatorMultiplication3D<T> plane(circle, mat);

        // Average concentration: use average functor
        SuperAverage3D<T> (densities[2], plane).operator()(result, dummy);
        const T avConc = result[0];
        if (iX == nX - 8) {  // at outlet
          clout << "av. conc. at outlet = " << avConc << std::endl;
        }
        csvWriter.setData(x, avConc);
      }
      csvWriter.writePNG();
    }
  }
}

void simulate(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& params = myCase.getParameters();
  const auto& converter = myCase.getLattice(NavierStokes{}).getUnitConverter();
  const T maxPhysT = params.get<parameters::PHYS_MAX_T>();
  const T maxPhysTAD = params.get<parameters::PHYS_MAX_T_AD>();
  const std::size_t iTmax = converter.getLatticeTime(maxPhysT);
  const int  statIterNS = converter.getLatticeTime( maxPhysT/100 );
  const int  statIterAD = converter.getLatticeTime( maxPhysTAD/100 );

  myCase.getLattice(NavierStokes{}).setStatisticsOff();
  myCase.getLattice(Concentration<0>{}).setStatisticsOff();
  myCase.getLattice(Concentration<1>{}).setStatisticsOff();
  myCase.getLattice(Concentration<2>{}).setStatisticsOff();

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    if(iT%statIterNS == 0){
      myCase.getLattice(NavierStokes{}).setStatisticsOn();
    }
    if(iT == 0){
      myCase.getLattice(Concentration<1>{}).setStatisticsOn();
    }

    setTemporalValues(myCase, iT);

    myCase.getLattice(NavierStokes{}).collideAndStream();

    getResultsNS(myCase, timer, iT);

    if(iT == 0){getResultsCRAD(myCase, timer, iT);}
    if(iT >= converter.getLatticeTime( maxPhysT - maxPhysTAD )){
      if(iT%statIterAD == 0){
        myCase.getLattice(Concentration<1>{}).setStatisticsOn();
      }
      if(iT%100 == 0){
        myCase.getOperator("LESReaction").apply();
      }
      myCase.getLattice(Concentration<0>{}).collideAndStream();
      myCase.getLattice(Concentration<1>{}).collideAndStream();
      myCase.getLattice(Concentration<2>{}).collideAndStream();
      getResultsCRAD(myCase, timer, iT);
    }
    if(iT%statIterNS == 0){
      myCase.getLattice(NavierStokes{}).setStatisticsOff();
    }
    if(iT%statIterAD == 0){
      myCase.getLattice(Concentration<1>{}).setStatisticsOff();
    }
  }

  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[]) {
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<DOMAIN_EXTENT        >({1600e-6,800e-6,100e-6});
    myCaseParameters.set<RESOLUTION           >(40);
    myCaseParameters.set<PHYS_MAX_T           >(25*(1600e-6 + 800e-6/2)/1.4);
    myCaseParameters.set<PHYS_MAX_T_AD        >(0.01);
    myCaseParameters.set<PHYS_CHAR_VELOCITY   >(1.4);
    myCaseParameters.set<PHYS_MAX_VELOCITY    >(1.4*2.1);
    myCaseParameters.set<LATTICE_CHAR_VELOCITY>(0.05);
    myCaseParameters.set<PHYS_CHAR_DENSITY    >(1000.);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY  >(1.002e-6);
    myCaseParameters.set<LX_INLET             >(100.e-6);
    myCaseParameters.set<LY_OUTLET            >(200.e-6);
    myCaseParameters.set<DIFFUSION_COMPONENTS >({1.6e-9,2.e-10,2.e-10});
    myCaseParameters.set<CONCENTRATION        >({1., 1.e-8,1.e-8,0.001,1.e-8,1.e-8});
    myCaseParameters.set<PLOT_OVER_LINE       >(true);
    myCaseParameters.set<PHYS_START_T         >(0.02*25*(1600e-6 + 800e-6/2)/1.4);
  }

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Lattice ===
  prepareLatticeNS(myCase);
  prepareLatticeCRAD<0>(myCase);
  prepareLatticeCRAD<1>(myCase);
  prepareLatticeCRAD<2>(myCase);

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValuesNS(myCase);
  setInitialValuesCRAD<0>(myCase);
  setInitialValuesCRAD<1>(myCase);
  setInitialValuesCRAD<2>(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);
}
