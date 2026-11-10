/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2016 Vojtech Cvrcek, Mathias J. Krause
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

/* powerLaw2d.cpp:
 * This example examines a steady flow of a non-newtonian fluid in a channel.
 * At the inlet, a profile for non-newtonian fluid is imposed on the velocity,
 * where as the outlet implements an outflow condition grad_x p = 0.
 * The power law model is for n=1 and m=charNu in fact the classical poiseuille flow.
 * One can validate the error with using functors in void error.
 *
 *
 */

#define _SMAGORINSKY // Smagorinsky power-law instead of power-law

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<FLOATING_POINT_TYPE, descriptors::D2Q9<OMEGA>>
>;

namespace olb::parameters {

  struct PERIODIC_BC        : public descriptors::TYPED_FIELD_BASE<bool,1> { };
  struct POWER_LAW_EXPONENT : public descriptors::FIELD_BASE<1> { };
  struct NU_MIN             : public descriptors::FIELD_BASE<1> { };
  struct NU_MAX             : public descriptors::FIELD_BASE<1> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const T physLengthX   = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T physDeltaX    = parameters.get<parameters::DOMAIN_EXTENT>()[1]/parameters.get<parameters::RESOLUTION>();
  const T physLengthY   = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const bool bcPeriodic = parameters.get<parameters::PERIODIC_BC>();

  const Vector extent{physLengthX, physLengthY};
  const Vector origin{0, 0};
  IndicatorCuboid2D<T> cuboid(extent, origin);

  #ifdef PARALLEL_MODE_MPI
    Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  #else
    Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, 1);
  #endif

  if (bcPeriodic) parameters.set<parameters::OVERLAP>(2);
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());

  // Periodic boundaries in x-direction
  if (bcPeriodic) mesh.getCuboidDecomposition().setPeriodicity({true, false});
  return mesh;
}

void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();

  const T physLengthX   = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T physLengthY   = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const T physDeltaX    = parameters.get<parameters::DOMAIN_EXTENT>()[1] / parameters.get<parameters::RESOLUTION>();
  const bool bcPeriodic = parameters.get<parameters::PERIODIC_BC>();

  Vector<T,2> extent(physLengthX, physLengthY);
  Vector<T,2> origin;

  geometry.rename(0, 2);
  geometry.rename(2, 1, {1, 1});

  // Set material number for inflow
  extent[0] = 1.2*physDeltaX;
  origin[0] = -physDeltaX;
  IndicatorCuboid2D<T> inflow(extent, origin);

  if (bcPeriodic) geometry.rename(1, 3, inflow);
  else geometry.rename(2, 3, 1, inflow);

  // Set material number for outflow
  origin[0] = physLengthX - .5*physDeltaX;
  IndicatorCuboid2D<T> outflow(extent, origin);

  if (bcPeriodic) geometry.rename(1, 4, outflow);
  else geometry.rename(2, 4, 1, outflow);

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  // Removes all not needed boundary voxels inside the surface
  geometry.innerClean();
  geometry.checkForErrors();
  geometry.getStatistics().print();

  clout << "Prepare Geometry ... OK" << std::endl;
  return;
}

// Set up the geometry of the simulation
void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes{});
  auto& geometry   = myCase.getGeometry();

  const T physDeltaX    = parameters.get<parameters::DOMAIN_EXTENT>()[1] / parameters.get<parameters::RESOLUTION>();
  const T physLengthX   = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T physLengthY   = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const bool bcPeriodic = parameters.get<parameters::PERIODIC_BC>();
  const int resolution  = parameters.get<parameters::RESOLUTION>();
  const T latticeRelaxationTime = parameters.get<parameters::LATTICE_RELAXATION_TIME>();
  const T maxVelocity   = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T Re            = parameters.get<parameters::REYNOLDS>();
  const T n             = parameters.get<parameters::POWER_LAW_EXPONENT>();
  const T physDensity   = parameters.get<parameters::PHYS_CHAR_DENSITY>();


  lattice.setUnitConverter<PowerLawUnitConverterFrom_Resolution_RelaxationTime_Reynolds_PLindex<T, DESCRIPTOR>>(
    resolution,                           // resolution: number of voxels per charPhysL
    latticeRelaxationTime,                // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    physLengthY,                          // charPhysLength: reference length of simulation geometry
    maxVelocity,                          // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    Re,                                   // Reynolds number
    n,                                    // power-law index
    physDensity                           // physDensity: physical density in __kg / m^3__
  );

  const auto& converter = lattice.getUnitConverter();
  lattice.getUnitConverter().print();

#ifdef _SMAGORINSKY
  dynamics::set<SmagorinskyPowerLawBGKdynamics>(lattice, geometry.getMaterialIndicator(1));
#else
  dynamics::set<PowerLawBGKdynamics>(lattice, geometry.getMaterialIndicator(1));
#endif

  boundary::set<boundary::BounceBack>(lattice, geometry.getMaterialIndicator(2));

  T distance2Wall = physDeltaX/2.;
  T p0 = converter.getPhysConsistencyCoeff()*util::pow(converter.getCharPhysVelocity(), n)*util::pow((n + 1.)/n, n)*util::pow(2./(physLengthY-distance2Wall*2), n + 1.);

  if (bcPeriodic) {
    dynamics::set(lattice,
                  geometry.getMaterialIndicator(3),
                  meta::id<typename PowerLawBGKdynamics<T,DESCRIPTOR>::template exchange_combination_rule<
                    powerlaw::PeriodicPressureOffset<-1,0>
                  >>{});
    lattice.setParameter<powerlaw::PRESSURE_OFFSET<-1,0>>(
      -converter.getLatticeDensityFromPhysPressure(p0*(physLengthX + distance2Wall*2.))+1);
  }
  else {
    dynamics::set<PowerLawBGKdynamics>(lattice, geometry.getMaterialIndicator(3));
    boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 3);
  }

  if (bcPeriodic) {
    dynamics::set(lattice,
                  geometry.getMaterialIndicator(4),
                  meta::id<typename PowerLawBGKdynamics<T,DESCRIPTOR>::template exchange_combination_rule<
                    powerlaw::PeriodicPressureOffset<1,0>
                  >>{});
    lattice.setParameter<powerlaw::PRESSURE_OFFSET<1,0>>(
      converter.getLatticeDensityFromPhysPressure(p0*(physLengthX + distance2Wall*2.))-1);
  }
  else {
    dynamics::set<PowerLawBGKdynamics>(lattice, geometry.getMaterialIndicator(4));
    boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 4);
  }

  lattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  lattice.setParameter<powerlaw::M>(converter.getLatticeConsistencyCoeff());
  lattice.setParameter<powerlaw::N>(n);

  const T nuMin = parameters.get<parameters::NU_MIN>();
  const T nuMax = parameters.get<parameters::NU_MAX>();
  lattice.setParameter<powerlaw::OMEGA_MIN>(1./(nuMax*descriptors::invCs2<T,DESCRIPTOR>() + 0.5));
  lattice.setParameter<powerlaw::OMEGA_MAX>(1./(nuMin*descriptors::invCs2<T,DESCRIPTOR>() + 0.5));
#ifdef _SMAGORINSKY
  lattice.setParameter<collision::LES::SMAGORINSKY>(parameters.get<parameters::SMAGORINSKY>());
#endif

  fields::set<descriptors::OMEGA>(lattice,
                                  geometry.getMaterialIndicator({1,3,4}),
                                  converter.getLatticeRelaxationFrequency());

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  OstreamManager clout(std::cout, "setInitialValues");
  using T = MyCase::value_t;
  auto& lattice    = myCase.getLattice(NavierStokes{});
  auto& geometry   = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  const auto& converter = lattice.getUnitConverter();

  const T n             = parameters.get<parameters::POWER_LAW_EXPONENT>();
  const T physLengthX   = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T physLengthY   = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const T physDeltaX    = parameters.get<parameters::DOMAIN_EXTENT>()[1] / parameters.get<parameters::RESOLUTION>();

  T distance2Wall = physDeltaX / 2;
  T p0 = converter.getPhysConsistencyCoeff()
       * util::pow(converter.getCharPhysVelocity(), n)
       * util::pow((n + 1.)/n, n)
       * util::pow(2./(physLengthY-distance2Wall*2), n + 1.);

  AnalyticalLinear2D<T,T> pressureF(-p0,
                                    0,
                                    p0*(physLengthX + distance2Wall*2)/2);
  PowerLaw2D<T> velocityF(geometry, 3, converter.getCharPhysVelocity(), distance2Wall, (n + 1.)/n);

  momenta::setPressureAndVelocity(lattice,
                                  geometry.getMaterialIndicator({1,2,3,4}),
                                  pressureF,
                                  velocityF);
  equilibria::setPressureAndVelocity(lattice,
                                     geometry.getMaterialIndicator({1,2,3,4}),
                                     pressureF,
                                     velocityF);

  lattice.initialize();

  clout << "Set Initial Values ... OK" << std::endl;
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{ }

// Compute error norms
void errorNorms(MyCase& myCase)
{
  OstreamManager clout(std::cout, "error");
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes{});
  auto& geometry   = myCase.getGeometry();

  const UnitConverter<T,DESCRIPTOR>& converter = lattice.getUnitConverter();

  const T physDeltaX    = parameters.get<parameters::DOMAIN_EXTENT>()[1] / parameters.get<parameters::RESOLUTION>();
  const T physLengthX   = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T physLengthY   = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const T n             = parameters.get<parameters::POWER_LAW_EXPONENT>();

  int input[1] = { };
  T result[1]  = { };

  T distance2Wall = physDeltaX/2.;

  PowerLaw2D<T> uSol(geometry, 3, converter.getCharPhysVelocity(), distance2Wall, (n + 1.)/n);
  SuperLatticePhysVelocity2D<T,DESCRIPTOR> u(lattice, converter);

  T p0 = converter.getPhysConsistencyCoeff()*util::pow(converter.getCharPhysVelocity(), n)*util::pow((n + 1.)/n, n)*util::pow(2./(physLengthY-distance2Wall*2), n + 1.);
  AnalyticalLinear2D<T,T> pressureSol(-p0, 0, p0*(physLengthX + distance2Wall*2.)/2.);
  SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure(lattice, converter);

  auto indicatorF = geometry.getMaterialIndicator(1);

  // velocity error
  SuperAbsoluteErrorL1Norm2D<T> absVelocityErrorNormL1(u, uSol, indicatorF);
  absVelocityErrorNormL1(result, input);
  clout << "velocity-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm2D<T> relVelocityErrorNormL1(u, uSol, indicatorF);
  relVelocityErrorNormL1(result, input);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm2D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, input);
  clout << "velocity-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm2D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, input);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm2D<T> absVelocityErrorNormLinf(u, uSol, indicatorF);
  absVelocityErrorNormLinf(result, input);
  clout << "velocity-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm2D<T> relVelocityErrorNormLinf(u, uSol, indicatorF);
  relVelocityErrorNormLinf(result, input);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;

  // pressure error
  SuperAbsoluteErrorL1Norm2D<T> absPressureErrorNormL1(pressure, pressureSol, indicatorF);
  absPressureErrorNormL1(result, input);
  clout << "pressure-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm2D<T> relPressureErrorNormL1(pressure, pressureSol, indicatorF);
  relPressureErrorNormL1(result, input);
  clout << "; pressure-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm2D<T> absPressureErrorNormL2(pressure, pressureSol, indicatorF);
  absPressureErrorNormL2(result, input);
  clout << "pressure-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm2D<T> relPressureErrorNormL2(pressure, pressureSol, indicatorF);
  relPressureErrorNormL2(result, input);
  clout << "; pressure-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm2D<T> absPressureErrorNormLinf(pressure, pressureSol, indicatorF);
  absPressureErrorNormLinf(result, input);
  clout << "pressure-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm2D<T> relPressureErrorNormLinf(pressure, pressureSol, indicatorF);
  relPressureErrorNormLinf(result, input);
  clout << "; pressure-Linf-error(rel)=" << result[0] << std::endl;
}

// Output to console and files
void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT
)
{
  OstreamManager clout(std::cout, "getResults");
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes{});

  const UnitConverter<T,DESCRIPTOR>& converter = lattice.getUnitConverter();

  const int vtkIter  = converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T>());

  SuperVTMwriter2D<T> vtmWriter("powerLaw2d");
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(lattice, converter);
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure(lattice, converter);
  SuperLatticePhysViscosity2D<T, DESCRIPTOR> viscosity(lattice, converter);

  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(pressure);
  vtmWriter.addFunctor(viscosity);

  if (iT==0) {
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid(lattice);
    SuperLatticeRank2D<T, DESCRIPTOR> rank(lattice);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();
  }

  if (iT%vtkIter==0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    vtmWriter.write(iT);

    SuperEuklidNorm2D<T, DESCRIPTOR> normVel(velocity);
    BlockReduction2D2D<T> planeReduction(normVel, 600, BlockDataSyncMode::ReduceOnly);
    // write output of velocity as JPEG
    heatmap::write(planeReduction, iT);

    timer.update(iT);
    timer.printStep();
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
    errorNorms(myCase);
  }
}

void simulate(MyCase& myCase){
  OstreamManager clout(std::cout,"simulate");
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;

  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(NavierStokes{});

  const UnitConverter<T,DESCRIPTOR>& converter = lattice.getUnitConverter();

  const int resolution  = parameters.get<parameters::RESOLUTION>();
  const T physDeltaX    = parameters.get<parameters::DOMAIN_EXTENT>()[1] / parameters.get<parameters::RESOLUTION>();
  const T physLengthX   = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T physLengthY   = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const T n             = parameters.get<parameters::POWER_LAW_EXPONENT>();
  const T maxVelocity   = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const bool bcPeriodic = parameters.get<parameters::PERIODIC_BC>();

  const std::size_t iTmax = converter.getLatticeTime(
    parameters.get<parameters::PHYS_MAX_T>());

  clout << "Starting simulation.. " << std::endl;
  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  util::ValueTracer<T> converge(converter.getLatticeTime(0.1 * parameters.get<parameters::PHYS_VTK_ITER_T>()),
                                parameters.get<parameters::CONVERGENCE_PRECISION>());

  for (std::size_t iT=0; iT<iTmax; ++iT) {

    if (converge.hasConverged()) {
      clout << "Convergence reached." << std::endl;
      parameters.set<parameters::CONVERGED>(true);
      getResults(myCase, timer, iT);
      clout << "Converged after " << iT << " iterations." << std::endl;
      break;
    }

    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2 & Step 8.3: Computation and Output of the Results, Collide and Stream Execution ===
    getResults(myCase, timer, iT);
    converge.takeValue(lattice.getStatistics().getAverageEnergy(), true);
    lattice.collideAndStream();

    if (bcPeriodic) {
      lattice.stripeOffDensityOffset(lattice.getStatistics().getAverageRho()-(T)1);
    }
  }
  timer.stop();
  timer.printSummary();

  // === 7th Step: Gnuplot ===
  Gnuplot<T> gplot("centerVelocity");
  T charLengthY = converter.getLatticeLength(physLengthY);
  for (int iY=0; iY<=charLengthY; ++iY) {
    T dx = 1. / resolution;
    T point[2]= {T(),T()};
    point[0] = .9*physLengthX;
    point[1] = (T)iY/charLengthY;
    std::vector<T> axisPoint(2, T());
    axisPoint[0] = physLengthX/2.;
    axisPoint[1] = physLengthY/2.;
    std::vector<T> axisDirection(2, T());
    axisDirection[0] = 1;
    axisDirection[1] = 0;
    T distance2Wall = physDeltaX/2.;
    PowerLaw2D<T> uSol(geometry, 3, maxVelocity, distance2Wall, (n + 1.)/n);
    T analytical[2] = {T(), T()};
    uSol(analytical, point);
    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(lattice, converter);
    AnalyticalFfromSuperF2D<T> intpolateVelocity(velocity, true);
    T numerical[2] = {T(),T()};
    intpolateVelocity(numerical,point);
    gplot.setData(iY*dx, {analytical[0], numerical[0]}, {"analytical", "numerical"});
  }
  // Create PNG file
  gplot.writePNG();
}


int main( int argc, char* argv[] )
{
  initialize( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  XMLreader config("input.xml");

  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    using T = MyCase::value_t;
    myCaseParameters.set<RESOLUTION>(config["geometry"]["N"].get<int>());
    myCaseParameters.set<DOMAIN_EXTENT>({config["geometry"]["lx"].get<T>(), config["geometry"]["ly"].get<T>()});
    myCaseParameters.set<PHYS_CHAR_VELOCITY>(config["dynamics"]["maxU"].get<T>());
    myCaseParameters.set<REYNOLDS>(config["dynamics"]["Re"].get<T>());
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(config["dynamics"]["tau"].get<T>());
    myCaseParameters.set<POWER_LAW_EXPONENT>(config["dynamics"]["n"].get<T>());
    myCaseParameters.set<PHYS_MAX_T>(config["time"]["Tmax"].get<T>());
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1.0);
    myCaseParameters.set<NU_MIN>(2.9686e-3);
    myCaseParameters.set<NU_MAX>(3.1667);
    myCaseParameters.set<SMAGORINSKY>(0.15);
    myCaseParameters.set<PERIODIC_BC>(true);

    myCaseParameters.set<PHYS_VTK_ITER_T>(config["time"]["Tprint"].get<T>());
    myCaseParameters.set<CONVERGENCE_PRECISION>(1e-6);
    myCaseParameters.set<CONVERGED>(false);

  }

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase);

  /// === Step 7: Set Initial Conditions ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);
}
