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

/* Header related to squareCavity2dTurbulent.cpp
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
  struct OUTPUT_VALUE_WRAPPER             : public descriptors::FIELD_BASE<5            > { };
  struct STATISTICS_ENSEMBLES   : public descriptors::FIELD_BASE<1            > { };
  struct STOP_SIM               : public descriptors::TYPED_FIELD_BASE<bool, 1> { };
}

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D2Q9<descriptors::FORCE, descriptors::TAU_EFF>>,
  Temperature,  Lattice<double, descriptors::D2Q5<descriptors::VELOCITY, descriptors::TAU_EFF>>
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

void prepareGeometry(MyCase& myCase)
{
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

void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T = MyCase::value_t_of<NavierStokes>;
  using NSEDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ADEDESCRIPTOR = MyCase::descriptor_t_of<Temperature>;

  auto& geometry = myCase.getGeometry();
  auto& NSElattice = myCase.getLattice(NavierStokes{});
  auto& ADElattice = myCase.getLattice(Temperature{});
  auto& parameters = myCase.getParameters();

  const T physCharLength          = parameters.get<parameters::PHYS_CHAR_LENGTH>();
  const int N                     = parameters.get<parameters::RESOLUTION>();
  const T physViscosity           = parameters.get<parameters::PHYS_KINEMATIC_VISCOSITY>();
  const T physDeltaX              = parameters.get<parameters::PHYS_DELTA_X>();
  const T physCharVelocity        = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T physDeltaT              = 2. * 0.056 / physCharVelocity * physCharLength / N;
  const T physDensity             = parameters.get<parameters::PHYS_CHAR_DENSITY>();
  const T physThermalExpansion    = parameters.get<parameters::PHYS_THERMAL_EXPANSION>();
  const T physThermalConductivity = parameters.get<parameters::PHYS_THERMAL_CONDUCTIVITY>();
  const T physHeatCapacity        = parameters.get<parameters::PHYS_HEAT_CAPACITY>();
  const T g                       = parameters.get<parameters::GRAVITATIONAL_ACC>();
  const T smagoConst              = parameters.get<parameters::SMAGORINSKY>();
  const T prTurb                  = parameters.get<parameters::PRANDTL_TURB>();
  const T Tcold                   = parameters.get<parameters::T_COLD>();
  const T Thot                    = parameters.get<parameters::T_HOT>();
  const T Tmean                   = parameters.get<parameters::T_MEAN>();

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

  dynamics::set<ExternalTauEffLESForcedBGKdynamics>(NSElattice, geometry.getMaterialIndicator({ 1, 2, 3 }));
  dynamics::set<ExternalTauEffLESBGKadvectionDiffusionDynamics>(ADElattice, geometry.getMaterialIndicator({ 1, 2, 3 }));

  boundary::set<boundary::BounceBack>(NSElattice, geometry, 4);
  boundary::set<boundary::BounceBack>(ADElattice, geometry, 4);

  /// sets boundary
  boundary::set<boundary::LocalVelocity>(NSElattice, geometry.getMaterialIndicator({ 2, 3 }));
  boundary::set<boundary::AdvectionDiffusionDirichlet>(ADElattice, geometry.getMaterialIndicator({ 2, 3 }));

  const T omegaNSE  =  converter.getLatticeRelaxationFrequency();
  const T omegaADE  =  converter.getLatticeThermalRelaxationFrequency();

  fields::set<descriptors::TAU_EFF>(NSElattice, geometry.getMaterialIndicator({ 1, 2, 3 }), 1. / omegaNSE);
  fields::set<descriptors::TAU_EFF>(ADElattice, geometry.getMaterialIndicator({ 1, 2, 3 }), 1. / omegaADE);

  NSElattice.setParameter<descriptors::OMEGA>( omegaNSE );
  ADElattice.setParameter<descriptors::OMEGA>( omegaADE );

  const T boussinesqForcePrefactor = g
                                   * converter.getConversionFactorTime()
                                   * converter.getCharPhysTemperatureDifference()
                                   * converter.getPhysThermalExpansionCoefficient()
                                   / converter.getConversionFactorVelocity();

  const T preFactor = smagoConst * smagoConst
                    * descriptors::invCs2<T,NSEDESCRIPTOR>() * descriptors::invCs2<T,NSEDESCRIPTOR>()
                    * 2 * util::sqrt(2);

  auto& coupling = myCase.setCouplingOperator(
    "Boussinesq",
    SmagorinskyBoussinesqCoupling{},
    names::NavierStokes{}, NSElattice,
    names::Temperature{},  ADElattice
  );
  coupling.setParameter<SmagorinskyBoussinesqCoupling::T0>(
    converter.getLatticeTemperature(Tmean)
  );
  coupling.setParameter<SmagorinskyBoussinesqCoupling::FORCE_PREFACTOR>(
    boussinesqForcePrefactor * Vector<T, MyCase::d>{ 0.0, 1.0 }
  );
  coupling.setParameter<SmagorinskyBoussinesqCoupling::SMAGORINSKY_PREFACTOR>(
    preFactor
  );
  coupling.setParameter<SmagorinskyBoussinesqCoupling::PR_TURB>(
    prTurb
  );
  coupling.setParameter<SmagorinskyBoussinesqCoupling::OMEGA_NSE>(
    omegaNSE
  );
  coupling.setParameter<SmagorinskyBoussinesqCoupling::OMEGA_ADE>(
    omegaADE
  );

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase)
{
    OstreamManager clout(std::cout,"setInitialValues");
    clout << "Set initial values ..." << std::endl;

    using T = MyCase::value_t_of<NavierStokes>;

    auto& geometry = myCase.getGeometry();
    auto& NSElattice = myCase.getLattice(NavierStokes{});
    auto& ADElattice = myCase.getLattice(Temperature{});
    const auto& converter = NSElattice.getUnitConverter();

    const T NSEomega = converter.getLatticeRelaxationFrequency();
    const T ADEomega = converter.getLatticeThermalRelaxationFrequency();
    const T Tcold = converter.getCharPhysLowTemperature();
    const T Thot  = converter.getCharPhysHighTemperature();
    const T Tmean = (Thot + Tcold) / 2.;

    /// define initial conditions
    momenta::setTemperature(ADElattice, geometry.getMaterialIndicator(1), Tmean);
    momenta::setTemperature(ADElattice, geometry.getMaterialIndicator(2), Thot);
    momenta::setTemperature(ADElattice, geometry.getMaterialIndicator(3), Tcold);

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

void computeNusselt(MyCase& myCase)
{
  OstreamManager clout(std::cout, "computeNusselt");

  using T = MyCase::value_t_of<NavierStokes>;

  auto& geometry      = myCase.getGeometry();
  auto& NSElattice    = myCase.getLattice(NavierStokes{});
  auto& ADElattice    = myCase.getLattice(Temperature{});
  auto& converter     = NSElattice.getUnitConverter();
  auto& parameters    = myCase.getParameters();

  const int N = converter.getResolution();
  int material = 0, voxel = 0;
  T T_x = 0, T_xplus1 = 0, T_xplus2 = 0, q = 0;

  for (int iC = 0; iC < NSElattice.getLoadBalancer().size(); iC++) {
    int ny = NSElattice.getBlock(iC).getNy();
    int iX = 0;
    for (int iY = 0; iY < ny; ++iY) {
      material = geometry.getBlockGeometry(iC).getMaterial(iX,iY);

      T_x = ADElattice.getBlock(iC).get(iX,iY).computeRho();
      T_xplus1 = ADElattice.getBlock(iC).get(iX+1,iY).computeRho();
      T_xplus2 = ADElattice.getBlock(iC).get(iX+2,iY).computeRho();

      if ( material == 2 ) {
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
}

void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                int iT)
{
  OstreamManager clout(std::cout, "getResults");

  using NSEDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ADEDESCRIPTOR = MyCase::descriptor_t_of<Temperature>;
  using T = MyCase::value_t_of<NavierStokes>;

  auto& NSElattice        = myCase.getLattice(NavierStokes{});
  auto& ADElattice        = myCase.getLattice(Temperature{});
  const auto& converter   = NSElattice.getUnitConverter();
  auto& parameters        = myCase.getParameters();

  const int statIter      = converter.getLatticeTime(parameters.get<parameters::PHYS_STAT_ITER_T>());
  const int vtkIter       = converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T>());
  const bool converged    = parameters.get<parameters::CONVERGED>();
  const bool stopSim      = parameters.get<parameters::STOP_SIM>();

  if (iT == 0)
  {
    NSElattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperVTMwriter2D<T>                     vtkWriter("squareCavity2dTurbulent");

    SuperLatticeCuboid2D<T, NSEDESCRIPTOR>  cuboid(NSElattice);
    SuperLatticeRank2D<T, NSEDESCRIPTOR>    rank(NSElattice);

    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.createMasterFile();

    NSElattice.setProcessingContext(ProcessingContext::Simulation);
  }

  if ((iT % vtkIter == 0 && iT > 0) || stopSim)
  {
    NSElattice.setProcessingContext(ProcessingContext::Evaluation);
    ADElattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperVTMwriter2D<T> vtkWriter("squareCavity2dTurbulent");

    SuperLatticePhysVelocity2D<T, NSEDESCRIPTOR> velocity(NSElattice, converter);
    SuperLatticePhysPressure2D<T, NSEDESCRIPTOR> pressure(NSElattice, converter);
    SuperLatticePhysTemperature2D<T, NSEDESCRIPTOR, ADEDESCRIPTOR> temperature(ADElattice, converter);

    vtkWriter.addFunctor(pressure);
    vtkWriter.addFunctor(velocity);
    vtkWriter.addFunctor(temperature);

    vtkWriter.write(iT);

    NSElattice.setProcessingContext(ProcessingContext::Simulation);
    ADElattice.setProcessingContext(ProcessingContext::Simulation);
  }

  if (iT % statIter == 0 || stopSim)
  {
    NSElattice.setProcessingContext(ProcessingContext::Evaluation);
    ADElattice.setProcessingContext(ProcessingContext::Evaluation);

    const T Thot            = parameters.get<parameters::T_HOT>();
    const T Tcold           = parameters.get<parameters::T_COLD>();

    SuperLatticePhysVelocity2D<T, NSEDESCRIPTOR> velocity(NSElattice, converter);
    SuperLatticePhysTemperature2D<T, NSEDESCRIPTOR, ADEDESCRIPTOR> temperature(ADElattice, converter);

    BlockReduction2D2D<T> planeReduction(temperature, 600, BlockDataSyncMode::ReduceOnly);
    BlockGifWriter<T> gifWriter;
    gifWriter.write(planeReduction, Tcold-1, Thot+1, iT, "temperature");

    SuperEuklidNorm2D<T, NSEDESCRIPTOR> normVel( velocity );
    BlockReduction2D2D<T> planeReduction2(normVel, 600, BlockDataSyncMode::ReduceOnly);
    BlockGifWriter<T> gifWriter2;
    gifWriter2.write( planeReduction2, iT, "velocity" );

    timer.printStep();
    /// NSELattice statistics console output
    NSElattice.getStatistics().print(iT,converter.getPhysTime(iT));
    /// ADElattice statistics console output
    ADElattice.getStatistics().print(iT,converter.getPhysTime(iT));

    NSElattice.setProcessingContext(ProcessingContext::Simulation);
    ADElattice.setProcessingContext(ProcessingContext::Simulation);
  }

  if ( converged )
  {
    NSElattice.setProcessingContext(ProcessingContext::Evaluation);
    ADElattice.setProcessingContext(ProcessingContext::Evaluation);

    NSElattice.communicate();
    ADElattice.communicate();

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
      if (n == 0) {
        outputVelX[0] = velX[0];
        outputVelX[1] = posY[0];
        outputVelY[0] = velY[0];
        outputVelY[1] = posX[0];
      }
      /// look for the maximum velocity in xDirection and the corresponding position in yDirection
      if (n > 0 && velX[n] > outputVelX[0]) {
        outputVelX[0] = velX[n];
        outputVelX[1] = posY[n];
      }
      /// look for the maximum velocity in yDirection and the corresponding position in xDirection
      if (n > 0 && velY[n] > outputVelY[0]) {
        outputVelY[0] = velY[n];
        outputVelY[1] = posX[n];
      }
    }

    parameters.set<parameters::OUTPUT_VALUE_WRAPPER>(
      {outputVelX[0],
       outputVelY[0],
       outputVelX[1],
       outputVelY[1],
       parameters.get<parameters::NUSSELT>()}
    );

    timer.printStep();
    /// NSELattice statistics console output
    NSElattice.getStatistics().print(iT,converter.getPhysTime(iT));
    /// ADElattice statistics console output
    ADElattice.getStatistics().print(iT,converter.getPhysTime(iT));
  }
}

void simulate(MyCase& myCase)
{
  OstreamManager clout(std::cout,"Simulation");
  clout << "Starting Simulation ..." << std::endl;

  using NSEDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ADEDESCRIPTOR = MyCase::descriptor_t_of<Temperature>;
  using T             = MyCase::value_t_of<NavierStokes>;

  auto& NSElattice              = myCase.getLattice(NavierStokes{});
  auto& ADElattice              = myCase.getLattice(Temperature{});
  const auto& converter         = NSElattice.getUnitConverter();
  auto& parameters              = myCase.getParameters();

  const std::size_t iTmax       = converter.getLatticeTime(parameters.get<parameters::PHYS_MAX_T>());

  const int convIter = parameters.get<parameters::CONV_ITER>();
  util::ValueTracer<T> converge(6, parameters.get<parameters::CONVERGENCE_PRECISION>());

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  SuperLatticePhysTemperature2D<T,NSEDESCRIPTOR, ADEDESCRIPTOR> sTemp(ADElattice, converter);
  SuperLatticePhysVelocity2D<T,NSEDESCRIPTOR>                   sVel(NSElattice, converter);

  SuperLatticeTimeAveragedF2D<T>                                sAveragedTemp(sTemp);
  SuperLatticeTimeAveragedF2D<T>                                sAveragedVel(sVel);
  SuperLatticeTimeAveragedCrossCorrelationF2D<T>                sAveragedTempVelCross(sTemp, sVel);
  SuperLatticeTimeAveragedCrossCorrelationF2D<T>                sAveragedVelVelCross(sVel, sVel);

  for (std::size_t iT = 0; iT < iTmax; ++iT)
  {
    setTemporalValues(myCase, iT);

    NSElattice.collideAndStream();
    ADElattice.collideAndStream();
    myCase.getOperator("Boussinesq").apply();

    if (iT % convIter == 0 && !parameters.get<parameters::CONVERGED>())
    {
      computeNusselt(myCase);
      converge.takeValue(parameters.get<parameters::NUSSELT>(), true);
    }

    getResults(myCase, timer, iT);

    if (converge.hasConverged() && !parameters.get<parameters::STOP_SIM>())
    {
      parameters.set<parameters::CONVERGED>(true);

      SuperVTMwriter2D<T> vtkWriter("squareCavity2dTurbulent");

      sAveragedTemp.addEnsemble();
      sAveragedVel.addEnsemble();
      sAveragedTempVelCross.addEnsemble();
      sAveragedVelVelCross.addEnsemble();

      vtkWriter.write(sAveragedTemp);
      vtkWriter.write(sAveragedVel);
      vtkWriter.write(sTemp);
      vtkWriter.write(sVel);

      if (sAveragedTemp.getEnsembles() >= parameters.get<parameters::STATISTICS_ENSEMBLES>())
      {
        parameters.set<parameters::STOP_SIM>(true);
      }
    }

    if (parameters.get<parameters::STOP_SIM>())
    {
      clout << "Stopping simulation." << std::endl;
      break;
    }

    timer.update(iT);
  }

  timer.stop();
  timer.printSummary();
}
