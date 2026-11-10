/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Florian Kaiser
 *                2008 Orestis Malaspinas
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

/* Header related to squareCavity2dLaminar.cpp
 * The reference is the paper in "Gaedtke, M., Wachter, S., Raedle, M., Nirschl, H., & Krause, M. J. (2018).
 * Application of a lattice Boltzmann method combined with a Smagorinsky turbulence model to spatially resolved heat flux inside a refrigerated vehicle.
 * Computers & Mathematics with Applications, 76(10), 2315-2329."
 */

// natural convection of air in a square cavity in 2D

#include <olb.h>

using namespace olb;
using namespace olb::graphics;
using namespace olb::names;

namespace olb::parameters {
  struct OUTPUT_VALUE_WRAPPER  : public descriptors::FIELD_BASE<5> { };
}

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<float, descriptors::D2Q9<descriptors::FORCE>>,
  Temperature,  Lattice<float, descriptors::D2Q5<descriptors::VELOCITY>>
>;

/// @brief Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters)
{
  using T = MyCase::value_t_of<NavierStokes>;

  const T dx = parameters.get<parameters::PHYS_DELTA_X>();

  const Vector domainExtend = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector extend = {domainExtend[0] + dx, domainExtend[1] + dx};

  const Vector domainOrigin = parameters.get<parameters::DOMAIN_ORIGIN>();
  const Vector origin = {domainOrigin[0] - (T) dx / 2, domainOrigin[1] - (T) dx / 2};
  IndicatorCuboid2D<T> cuboid(extend, origin);

  Mesh<T,MyCase::d> mesh(cuboid, dx, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({ false, false });

  return mesh;
}

void prepareGeometry(MyCase& myCase){
  OstreamManager clout(std::cout, "preprareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t_of<NavierStokes>;

  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  const T dx                = parameters.get<parameters::PHYS_DELTA_X>();
  const Vector extend       = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector origin       = parameters.get<parameters::DOMAIN_ORIGIN>();

  geometry.rename(0, 4, {1, 1});

  IndicatorCuboid2D<T> bulkCuboid(extend, origin);
  geometry.rename(4, 1, bulkCuboid);

  Vector extendwallleft{dx, extend[1] + dx};
  Vector originwallleft{origin[0] - (T) dx / 2, origin[1] - (T) dx / 2};
  IndicatorCuboid2D<T> wallleft(extendwallleft, originwallleft);
  geometry.rename( 4, 2, 1, wallleft  );

  Vector extendwallright{dx, extend[1] + dx};
  Vector originwallright{origin[0] + extend[0], origin[1] - (T) dx / 2};
  IndicatorCuboid2D<T> wallright(extendwallright, originwallright);

  geometry.rename( 4, 3, 1, wallright );

  geometry.clean();
  geometry.checkForErrors();

  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase){
  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T             = MyCase::value_t_of<NavierStokes>;
  using NSEDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ADEDESCRIPTOR = MyCase::descriptor_t_of<Temperature>;

  auto& geometry    = myCase.getGeometry();
  auto& NSElattice  = myCase.getLattice(NavierStokes{});
  auto& ADElattice  = myCase.getLattice(Temperature{});
  auto& parameters  = myCase.getParameters();

  const T physCharLength          = parameters.get<parameters::PHYS_CHAR_LENGTH>();
  const T tau                     = parameters.get<parameters::LATTICE_RELAXATION_TIME>();
  const T physViscosity           = parameters.get<parameters::PHYS_KINEMATIC_VISCOSITY>();
  const T physDeltaX              = parameters.get<parameters::PHYS_DELTA_X>();
  const T physDeltaT              = (tau - (T) 1 / 2) / descriptors::invCs2<T,NSEDESCRIPTOR>() * physDeltaX * physDeltaX / physViscosity;
  const T physCharVelocity        = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T physDensity             = parameters.get<parameters::PHYS_CHAR_DENSITY>();
  const T physThermalExpansion    = parameters.get<parameters::PHYS_THERMAL_EXPANSION>();
  const T physThermalConductivity = parameters.get<parameters::PHYS_THERMAL_CONDUCTIVITY>();
  const T physHeatCapacity        = parameters.get<parameters::PHYS_HEAT_CAPACITY>();
  const T g                       = parameters.get<parameters::GRAVITATIONAL_ACC>();
  const T Tcold                   = parameters.get<parameters::T_COLD>();
  const T Thot                    = parameters.get<parameters::T_HOT>();

  NSElattice.setUnitConverter<ThermalUnitConverter<T,NSEDESCRIPTOR,ADEDESCRIPTOR>>(
    (T) physDeltaX,
    (T) physDeltaT,
    (T) physCharLength,
    (T) physCharVelocity,
    (T) physViscosity,
    (T) physDensity,
    (T) physThermalConductivity,
    (T) physHeatCapacity,
    (T) physThermalExpansion,
    (T) Tcold,
    (T) Thot
  );
  const auto& converter = NSElattice.getUnitConverter();
  converter.print();

  ADElattice.setUnitConverter(converter);

  dynamics::set<ForcedBGKdynamics>(NSElattice, geometry.getMaterialIndicator({ 1, 2, 3 }));
  dynamics::set<AdvectionDiffusionBGKdynamics>(ADElattice, geometry.getMaterialIndicator({ 1, 2, 3 }));

  boundary::set<boundary::BounceBack>(ADElattice, geometry.getMaterialIndicator(4));
  boundary::set<boundary::BounceBack>(NSElattice, geometry.getMaterialIndicator(4));
  boundary::set<boundary::AdvectionDiffusionDirichlet>(ADElattice, geometry.getMaterialIndicator({ 2, 3 }));
  boundary::set<boundary::LocalVelocity>(NSElattice, geometry.getMaterialIndicator({ 2, 3 }));

  T boussinesqForcePrefactor = g * converter.getConversionFactorTime()
                                 * converter.getCharPhysTemperatureDifference()
                                 * converter.getPhysThermalExpansionCoefficient()
                                 / converter.getConversionFactorVelocity();

  auto& coupling = myCase.setCouplingOperator(
    "Boussinesq",
    NavierStokesAdvectionDiffusionCoupling{},
    names::NavierStokes{}, NSElattice,
    names::Temperature{},  ADElattice
  );
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::T0>(
    converter.getLatticeTemperature(Tcold)
  );
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::FORCE_PREFACTOR>(
    boussinesqForcePrefactor * Vector<T,2>{0.0,1.0}
  );

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase){
  OstreamManager clout(std::cout,"setInitialValues");
  clout << "Set initial values ..." << std::endl;

  using T               = MyCase::value_t_of<NavierStokes>;

  auto& geometry        = myCase.getGeometry();
  auto& NSElattice      = myCase.getLattice(NavierStokes{});
  auto& ADElattice      = myCase.getLattice(Temperature{});
  const auto& converter = NSElattice.getUnitConverter();

  const T NSEomega      = converter.getLatticeRelaxationFrequency();
  const T ADEomega      = converter.getLatticeThermalRelaxationFrequency();
  const T Tcold         = converter.getCharPhysLowTemperature();
  const T Thot          = converter.getCharPhysHighTemperature();
  const T Tmean         = (Thot + Tcold) / 2.;

  /// define initial conditions

  momenta::setTemperature(ADElattice, geometry.getMaterialIndicator(1), Tmean);
  momenta::setTemperature(ADElattice, geometry.getMaterialIndicator(2), Thot);
  momenta::setTemperature(ADElattice, geometry.getMaterialIndicator(3), Tcold);

  NSElattice.setParameter<descriptors::OMEGA>(NSEomega);
  ADElattice.setParameter<descriptors::OMEGA>(ADEomega);

  /// Make the lattice ready for simulation
  NSElattice.initialize();
  ADElattice.initialize();

  clout << "Set initial values ... OK" << std::endl;
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{
  // Nothing to do
}

void computeNusselt(MyCase& myCase){
  OstreamManager clout(std::cout, "computeNusselt");

  using T = MyCase::value_t_of<NavierStokes>;

  auto& geometry      = myCase.getGeometry();
  auto& NSElattice    = myCase.getLattice(NavierStokes{});
  auto& ADElattice    = myCase.getLattice(Temperature{});
  auto& converter     = NSElattice.getUnitConverter();
  auto& parameters    = myCase.getParameters();

  const int N         = converter.getResolution();

  ADElattice.setProcessingContext(ProcessingContext::Evaluation);
  NSElattice.setProcessingContext(ProcessingContext::Evaluation);

  int material = 0, voxel = 0;
  T T_x = 0, T_xplus1 = 0, T_xplus2 = 0, q = 0;

  for (int iC = 0; iC < NSElattice.getLoadBalancer().size(); iC++)
  {
    int ny = NSElattice.getBlock(iC).getNy();
    int iX = 0;
    for (int iY = 0; iY < ny; ++iY)
    {
      material = geometry.getBlockGeometry(iC).getMaterial(iX,iY);

      T_x = ADElattice.getBlock(iC).get(iX,iY).computeRho();
      T_xplus1 = ADElattice.getBlock(iC).get(iX+1,iY).computeRho();
      T_xplus2 = ADElattice.getBlock(iC).get(iX+2,iY).computeRho();

      if ( material == 2 )
      {
        q += (3.0*T_x - 4.0*T_xplus1 + 1.0*T_xplus2)/2.0*N;
        voxel++;
      }
    }
  }

  #ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(q, MPI_SUM);
    singleton::mpi().reduceAndBcast(voxel, MPI_SUM);
  #endif

  parameters.set<parameters::NUSSELT>(q / (T) voxel);

  ADElattice.setProcessingContext(ProcessingContext::Simulation);
  NSElattice.setProcessingContext(ProcessingContext::Simulation);
}

void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                int iT)
{
  using T = MyCase::value_t_of<NavierStokes>;
  using NSEDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ADEDESCRIPTOR = MyCase::descriptor_t_of<Temperature>;

  auto& NSElattice        = myCase.getLattice(NavierStokes{});
  auto& ADElattice        = myCase.getLattice(Temperature{});
  const auto& converter   = NSElattice.getUnitConverter();
  auto& parameters        = myCase.getParameters();

  const int statIter      = converter.getLatticeTime(parameters.get<parameters::PHYS_STAT_ITER_T>());
  const int vtkIter       = converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T>());
  const bool converged    = parameters.get<parameters::CONVERGED>();

  if (iT == 0)
  {
    NSElattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperVTMwriter2D<T> vtkWriter("squareCavity2dLaminar");
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid2D  cuboid(NSElattice);
    SuperLatticeRank2D    rank(NSElattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.createMasterFile();

    NSElattice.setProcessingContext(ProcessingContext::Simulation);
  }

  if ((iT % vtkIter == 0 && iT > 0) || converged)
  {
    NSElattice.setProcessingContext(ProcessingContext::Evaluation);
    ADElattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperVTMwriter2D<T> vtkWriter("squareCavity2dLaminar");
    SuperLatticePhysVelocity2D<T, NSEDESCRIPTOR>                    velocity(NSElattice, converter);
    SuperLatticePhysPressure2D<T, NSEDESCRIPTOR>                    pressure(NSElattice, converter);
    SuperLatticePhysTemperature2D<T, NSEDESCRIPTOR, ADEDESCRIPTOR>  temperature(ADElattice, converter);

    vtkWriter.addFunctor(velocity);
    vtkWriter.addFunctor(pressure);
    vtkWriter.addFunctor(temperature);

    vtkWriter.write(iT);

    NSElattice.setProcessingContext(ProcessingContext::Simulation);
    ADElattice.setProcessingContext(ProcessingContext::Simulation);
  }

  if (iT % statIter == 0 || converged)
  {
    NSElattice.setProcessingContext(ProcessingContext::Evaluation);
    ADElattice.setProcessingContext(ProcessingContext::Evaluation);

    const T Thot            = parameters.get<parameters::T_HOT>();
    const T Tcold           = parameters.get<parameters::T_COLD>();

    SuperLatticePhysVelocity2D<T, NSEDESCRIPTOR>                    velocity(NSElattice, converter);
    SuperLatticePhysTemperature2D<T, NSEDESCRIPTOR, ADEDESCRIPTOR>  temperature(ADElattice, converter);

    timer.printStep();
    /// NSELattice statistics console output
    NSElattice.getStatistics().print(iT,converter.getPhysTime(iT));
    /// ADElattice statistics console output
    ADElattice.getStatistics().print(iT,converter.getPhysTime(iT));

    BlockReduction2D2D<T> planeReduction(temperature, 600, BlockDataSyncMode::ReduceOnly);
    BlockGifWriter<T> gifWriter;
    gifWriter.write(planeReduction, Tcold-1, Thot+1, iT, "temperature");

    SuperEuklidNorm2D<T, NSEDESCRIPTOR> normVel( velocity );
    BlockReduction2D2D<T> planeReduction2(normVel, 600, BlockDataSyncMode::ReduceOnly);
    BlockGifWriter<T> gifWriter2;
    gifWriter2.write( planeReduction2, iT, "velocity" );

    NSElattice.setProcessingContext(ProcessingContext::Simulation);
    ADElattice.setProcessingContext(ProcessingContext::Simulation);
  }

  if ( converged )
  {
    NSElattice.setProcessingContext(ProcessingContext::Evaluation);
    ADElattice.setProcessingContext(ProcessingContext::Evaluation);

    computeNusselt(myCase);

    SuperLatticePhysVelocity2D<T, NSEDESCRIPTOR>  velocity(NSElattice, converter);
    AnalyticalFfromSuperF2D<T>                    interpolation(velocity, true);

    const T lx              = parameters.get<parameters::PHYS_CHAR_LENGTH>();
    /// Initialize vectors for data output
    T xVelocity[2] = { T() };
    T outputVelX[2] = { T() };
    T yVelocity[2] = { T() };
    T outputVelY[2] = { T() };

    const int outputSize = 512;
    Vector<T, outputSize> velX;
    Vector<T, outputSize> posX;
    Vector<T, outputSize> velY;
    Vector<T, outputSize> posY;

    /// loop for the resolution of the cavity at x = lx/2 in yDirection and vice versa
    for (int n = 0; n < outputSize; ++n)
    {
      T yPosition[2] = { lx / 2, lx * n / (T) outputSize };
      T xPosition[2] = { lx * n / (T) outputSize, lx / 2 };

      /// Interpolate xVelocity at x = lx/2 for each yPosition
      interpolation(xVelocity, yPosition);
      interpolation(yVelocity, xPosition);
      /// Store the interpolated values to compare them among each other in order to detect the maximum
      velX[n] = xVelocity[0];
      posY[n] = yPosition[1];
      velY[n] = yVelocity[1];
      posX[n] = xPosition[0];

      /// Initialize output with the corresponding velocities and positions at the origin
      if (n == 0)
      {
        outputVelX[0] = velX[0];
        outputVelX[1] = posY[0];
        outputVelY[0] = velY[0];
        outputVelY[1] = posX[0];
      }
      /// look for the maximum velocity in xDirection and the corresponding position in yDirection
      if (n > 0 && velX[n] > outputVelX[0])
      {
        outputVelX[0] = velX[n];
        outputVelX[1] = posY[n];
      }
      /// look for the maximum velocity in yDirection and the corresponding position in xDirection
      if (n > 0 && velY[n] > outputVelY[0])
      {
        outputVelY[0] = velY[n];
        outputVelY[1] = posX[n];
      }
    }

    parameters.set<parameters::OUTPUT_VALUE_WRAPPER>(
      {
        outputVelX[0],
        outputVelY[0],
        outputVelX[1],
        outputVelY[1],
        parameters.get<parameters::NUSSELT>()
      }
    );
  }
}

void simulate(MyCase& myCase){
  OstreamManager clout(std::cout,"Simulation");
  clout << "Starting Simulation ..." << std::endl;

  using T = MyCase::value_t_of<NavierStokes>;

  auto& NSElattice = myCase.getLattice(NavierStokes{});
  auto& ADElattice = myCase.getLattice(Temperature{});
  const auto& converter = NSElattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  const std::size_t iTmax = converter.getLatticeTime(parameters.get<parameters::PHYS_MAX_T>());

  const std::size_t convIter = parameters.get<parameters::CONV_ITER>();
  util::ValueTracer<T> converge(6, parameters.get<parameters::CONVERGENCE_PRECISION>());

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();
  for (std::size_t iT = 0; iT < iTmax; ++iT)
  {
    setTemporalValues(myCase, iT);

    if (converge.hasConverged() && !parameters.get<parameters::CONVERGED>())
    {
      parameters.set<parameters::CONVERGED>(true);
      clout << "Simulation converged." << std::endl;
      clout << "Time " << iT << "." << std::endl;

      getResults(myCase, timer, iT);
      break;
    }

    NSElattice.collideAndStream();
    ADElattice.collideAndStream();
    myCase.getOperator("Boussinesq").apply();

    getResults(myCase, timer, iT);

    if(iT % convIter == 0 && !parameters.get<parameters::CONVERGED>())
    {
      computeNusselt(myCase);
      converge.takeValue(parameters.get<parameters::NUSSELT>(), true);
    }

    timer.update(iT);
  }

  timer.stop();
  timer.printSummary();
}
