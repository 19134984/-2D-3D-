/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Fedor Bukreev, Jan E. Marquardt, Mathias J. Krause
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

// Experimental: Setup hasn't been exntensively validated.

/* How to convert TIFF to VTI:
 * 1. Convert the TIFF series to VTI using ParaView
 *    (Point data are mandatory, cell data do not work)
 *     - Open TIFF series in ParaView
 *     - Apply "cell data to point data" filter (if it's grayed out, applying the filter isn't necessary)
 *     - Save the result as VTI
 *     - Use ASCII formatting
 *     - Solely export one array and use the name as an input when running the simulation (see <arrayname>)
 * 2. Change "vtiName" and "arrayName" in the main() function
 * 3. Start the simulation
 */

 // The used VTI µCT-scan is taken from: https://doi.org/10.1016/j.energy.2021.122151

#include <olb.h>

using namespace olb;
using namespace olb::names;

namespace olb::parameters {

struct CONVERGENCE_CHECK_T  : public descriptors::FIELD_BASE<1> { };
struct SCALING_FACTOR       : public descriptors::FIELD_BASE<1> { };
struct TIME_SCALING_FACTOR  : public descriptors::FIELD_BASE<1> { };

}

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q19<>>
>;

std::shared_ptr<BlockVTIreader3D<MyCase::value_t, MyCase::value_t>> vtiReader; // to ensure that the data persists

std::shared_ptr<IndicatorBlockData3D<MyCase::value_t>>
generateIndicatorFromVTI(const std::string vtiFile, const std::string arrayName, MyCase::ParametersD& parameters)
{
  const MyCase::value_t sourceScale = parameters.get<parameters::SCALING_FACTOR>();

  vtiReader = std::make_shared<BlockVTIreader3D<MyCase::value_t, MyCase::value_t>>(vtiFile, arrayName);

  auto           cuboidSample     = vtiReader->getCuboid();
  MyCase::value_t              deltaRsample     = cuboidSample.getDeltaR() * sourceScale;
  Vector<int, 3> extentSample     = cuboidSample.getExtent();
  Vector<MyCase::value_t, 3>   originSamplePhys = cuboidSample.getOrigin() * sourceScale;
  Vector<MyCase::value_t, 3>   extentSamplePhys = {deltaRsample * MyCase::value_t(extentSample[0] + 0.5),
                                     deltaRsample * MyCase::value_t(extentSample[1] + 0.5),
                                     deltaRsample * MyCase::value_t(extentSample[2] + 0.5)};
  parameters.set<parameters::DOMAIN_EXTENT>(extentSamplePhys);

  std::shared_ptr<IndicatorBlockData3D<MyCase::value_t>> ind(
      new IndicatorBlockData3D<MyCase::value_t>(vtiReader->getBlockData(), extentSamplePhys,
                                  originSamplePhys, deltaRsample, true));

  return ind;
}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters, std::string vtiFile, std::string arrayName) {
  using T = MyCase::value_t;
  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();
  std::shared_ptr<IndicatorBlockData3D<T>> rock =
      generateIndicatorFromVTI(vtiFile, arrayName, parameters);
  IndicatorLayer3D<T> layer(*rock, physDeltaX);

  Mesh<T,MyCase::d> mesh(layer, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({false, true, true});
  return mesh;
}

void prepareGeometry( MyCase& myCase, std::string vtiFile,  std::string arrayName ) {
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();

  // Set the cells inside of the layer indicator to MN 2
  std::shared_ptr<IndicatorBlockData3D<T>> rock =
      generateIndicatorFromVTI(vtiFile, arrayName, parameters);
  IndicatorLayer3D<T> layer(*rock, physDeltaX);
  geometry.rename(0, 2, layer);

  // Set material number 1 for cells in the rock cavity with 1 cell offset in each direction
  geometry.rename(2, 1, {1,1,1});

  Vector<T, 3> origin = geometry.getStatistics().getMinPhysR(2);
  Vector<T, 3> extend = geometry.getStatistics().getMaxPhysR(2);

  // Set material number for inflow
  origin[0] = geometry.getStatistics().getMinPhysR(2)[0] -
              T {0.5} * physDeltaX;
  extend[0] = T {1} * physDeltaX;
  IndicatorCuboid3D<T> inflow(extend, origin);
  // MN 3 is set only if the neighbor cell has MN 1
  geometry.rename(2, 3, 1, inflow);

  // Set material number for outflow
  origin[0] = geometry.getStatistics().getMaxPhysR(2)[0] -
              T {0.5} * physDeltaX;
  extend[0] = T {1} * physDeltaX;
  IndicatorCuboid3D<T> outflow(extend, origin);
  // MN 4 is set only if the neighbor cell has MN 1
  geometry.rename(2, 4, 1, outflow);

  geometry.clean();

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  geometry.checkForErrors();

  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( MyCase& myCase ) {
  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});
  Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  T       tau = parameters.get<parameters::LATTICE_RELAXATION_TIME>();
  T       fluidDensity = parameters.get<parameters::PHYS_CHAR_DENSITY>();
  T       kinematicViscosity = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  T       pressureDrop = parameters.get<parameters::PRESSURE_DROP>();
  T       permeabilityForVelocityEstimation = T {1e-7};
  T       charPhysVelocityEstimation = permeabilityForVelocityEstimation * pressureDrop /
                                       (kinematicViscosity * fluidDensity * extent[0]);

  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR>>(
      parameters.get<parameters::RESOLUTION>(), // resolution: number of voxels per charPhysL
      (T)tau, // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
      (T)extent[1], // charPhysLength: reference length of simulation geometry
      (T)charPhysVelocityEstimation, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
      (T)kinematicViscosity, // physViscosity: physical kinematic viscosity in __m^2 / s__
      (T)fluidDensity // physDensity: physical density in __kg / m^3__
  );
  auto& converter = lattice.getUnitConverter();
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("rock");

  // Material=1 -->bulk dynamics
  auto bulkIndicator = geometry.getMaterialIndicator({1, 3, 4});
  dynamics::set<BGKdynamics>(lattice, bulkIndicator);

  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(lattice, geometry, 2);

  // Setting of the boundary conditions
  // local boundary conditions
  boundary::set<boundary::LocalPressure>(lattice, geometry, 3);
  boundary::set<boundary::LocalPressure>(lattice, geometry, 4);

  lattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  clout << "Prepare Lattice ... OK" << std::endl;
}
void setInitialValues( MyCase& myCase ) {
  auto& lattice   = myCase.getLattice(NavierStokes{});

  // Make the lattice ready for simulation
  lattice.initialize();
}

void setTemporalValues( MyCase& myCase,
                        std::size_t iT )
{
  OstreamManager clout(std::cout, "setTemporalValues");
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();
  const int iTmaxStart = converter.getLatticeTime(parameters.get<parameters::PHYS_START_T>());
  const T pressureDrop = parameters.get<parameters::PRESSURE_DROP>();

  // No of time steps for smooth start-up
  const int iTupdate = std::max(iTmaxStart / 200, 1);

  if (iT % iTupdate == 0 && int(iT) <= iTmaxStart) {
    // Smooth start curve, sinus
    // SinusStartScale<T,int> StartScale(iTmaxStart, T(1));

    // Smooth start curve, polynomial
    PolynomialStartScale<T, std::size_t> StartScale(iTmaxStart, T(1));

    // Creates and sets the Poiseuille inflow profile using functors
    T   frac[1]  = {};
    StartScale(frac, &iT);
    const T                 currentPressure = frac[0] * pressureDrop;

    momenta::setPressure(lattice, geometry.getMaterialIndicator({3}), currentPressure);

    clout << "step=" << iT << "; currentPressureDifference=" << currentPressure
          << std::endl;

    lattice.setProcessingContext<Array<momenta::FixedDensity::RHO>>(
        ProcessingContext::Simulation);
  }
}

void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  OstreamManager clout( std::cout,"getResults" );
  using T = MyCase::value_t;
  using DESCRIPTOR  = MyCase::descriptor_t_of<NavierStokes>;
  auto& lattice     = myCase.getLattice(NavierStokes{});
  auto& geometry    = myCase.getGeometry();
  auto& converter   = lattice.getUnitConverter();
  auto& parameters  = myCase.getParameters();
  const T maxPhysT            = parameters.get<parameters::PHYS_MAX_T>();
  const T lx                  = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T kinematicViscosity  = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T fluidDensity        = parameters.get<parameters::PHYS_CHAR_DENSITY>();
  const T pressureDrop        = parameters.get<parameters::PRESSURE_DROP>();
  bool converged              = parameters.get<parameters::CONVERGED>();
  const int statIter          = converter.getLatticeTime(parameters.get<parameters::PHYS_STAT_ITER_T>());
  const int vtkIter           = converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T>());
  T currentPermeability;

  SuperVTMwriter3D<T>                       vtmWriter("resolvedRock3d");
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(lattice, converter);
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(lattice, converter);

  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(pressure);

  if (iT == 0) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, DESCRIPTOR>   cuboid(lattice);
    SuperLatticeRank3D<T, DESCRIPTOR>     rank(lattice);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);

    vtmWriter.createMasterFile();
  }

  // Writes output on the console
  if ( iT % statIter == 0 || converged ) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    // Timer console output
    timer.update(iT);
    timer.printStep();

    // Lattice statistics console output
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));

    SuperAverage3D<T> avgVel(velocity, geometry, 1);
    int               input[1] {};
    Vector<T,3> currentAverageVelocity;
    avgVel(currentAverageVelocity.data(), input);
    clout << "Average x-velocity: " << currentAverageVelocity[0] << " m/s"
          << std::endl;

    currentPermeability = currentAverageVelocity[0] * kinematicViscosity *
                          fluidDensity * lx / pressureDrop;

    clout << "Permeability: " << currentPermeability << " m^2" << std::endl;
  }

  // Writes the vtk files
  if ( iT % vtkIter == 0 || maxPhysT ) {
    vtmWriter.write(iT);

    {
      SuperEuklidNorm3D<T>  normVel(velocity);
      BlockReduction3D2D<T> planeReduction(normVel, Vector<T, 3>({0, 0, 1}));
      // write output as JPEG
      heatmap::write(planeReduction, iT);
    }
  }
}

void simulate( MyCase& myCase ) {
  OstreamManager clout( std::cout, "simulation" );
  using T = MyCase::value_t;
  auto& lattice     = myCase.getLattice(NavierStokes{});
  auto& geometry    = myCase.getGeometry();
  auto& converter   = lattice.getUnitConverter();
  auto& parameters  = myCase.getParameters();
  const T   epsilon             = parameters.get<parameters::CONVERGENCE_PRECISION>();
  const T   pressureDrop        = parameters.get<parameters::PRESSURE_DROP>();
  const T   kinematicViscosity  = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T   fluidDensity        = parameters.get<parameters::PHYS_CHAR_DENSITY>();
  const T   timeScalingFactor   = parameters.get<parameters::TIME_SCALING_FACTOR>();
  const T   lx                  = parameters.get<parameters::DOMAIN_EXTENT>()[0];

  constexpr T permeabilityForVelocityEstimation = T {1e-7};
  // Estimating velocity using Darcy's law
  const T charPhysVelocityEstimation =
      permeabilityForVelocityEstimation * pressureDrop /
      (kinematicViscosity * fluidDensity * lx);
  // Assumption. Likely needs an update after some experience is gained
  const T maxPhysT = timeScalingFactor * lx / charPhysVelocityEstimation;

  parameters.set<parameters::PHYS_MAX_T             >( maxPhysT );
  parameters.set<parameters::PHYS_START_T           >( maxPhysT * 0.2 );
  parameters.set<parameters::PHYS_STAT_ITER_T       >( maxPhysT * .01 );
  parameters.set<parameters::PHYS_VTK_ITER_T        >( maxPhysT * .05 );
  parameters.set<parameters::CONVERGENCE_CHECK_T    >( maxPhysT / 100. );

  const T   iTmax               = converter.getLatticeTime(parameters.get<parameters::PHYS_MAX_T>());
  const int iTmaxStart          = converter.getLatticeTime(parameters.get<parameters::PHYS_START_T>());
  const int interval            = converter.getLatticeTime(parameters.get<parameters::CONVERGENCE_CHECK_T>());
  util::ValueTracer<T> converge(interval, epsilon);

  clout << "Timeframe to be simulated: " << maxPhysT << " s" << std::endl;

  clout << "starting simulation..." << std::endl;
  clout << "MaxIT: " << iTmax << std::endl;
  util::Timer<T> timer(iTmax,
                       geometry.getStatistics().getNvoxel());
  timer.start();

  parameters.set<parameters::CONVERGED>(false);
  for (std::size_t iT = 0; iT < iTmax; ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setTemporalValues( myCase, iT );

    // === 6th Step: Collide and Stream Execution ===
    lattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( myCase, timer, iT );

    if ( int(iT) > iTmaxStart ) {
      converge.takeValue(lattice.getStatistics().getAverageEnergy(), false);
    }

    if (converge.hasConverged()) {
      parameters.set<parameters::CONVERGED>( true );
      getResults( myCase, timer, iT );
      break;
    }
  }

  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<PHYS_DELTA_X               >(     5e-6 );
    myCaseParameters.set<LATTICE_RELAXATION_TIME    >(      .75 );
    myCaseParameters.set<PHYS_CHAR_DENSITY          >(     1000 );
    myCaseParameters.set<PHYS_CHAR_VISCOSITY        >(     1e-6 );
    myCaseParameters.set<CONVERGENCE_PRECISION      >(     1e-6 );
    myCaseParameters.set<RESOLUTION                 >(       50 );
    myCaseParameters.set<REYNOLDS                   >(       1. );
    myCaseParameters.set<CONVERGED                  >(    false );
    myCaseParameters.set<PRESSURE_DROP              >(       1. );
    myCaseParameters.set<SCALING_FACTOR             >(       1. );
    myCaseParameters.set<TIME_SCALING_FACTOR        >(       4. );
  }

  std::string vtiFile = "rock.vti";
  std::string arrayName = "scalars";

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters, vtiFile, arrayName);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase, vtiFile, arrayName);

  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase);

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);
}
