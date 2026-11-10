/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006-2025 Mathias J. Krause, Fabian Klemens,
 *  Julius Jessberger, Shota Ito
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


/** @file testFlow3d.cpp
 * @brief In this example, a 3d flow in a cube is simulated.
 *
 * It shows the basic structure of an OpenLB simulation and should be suitable for beginners.
 *
 * As every OpenLB simulation, it consists of eight steps:
 * @li Step 1: Declarations
 * @li Step 2: Initialization
 * @li Step 3: Create mesh
 * @li Step 4: Create case
 * @li Step 5: Prepare geometry
 * @li Step 6: Prepare lattice
 * @li Step 7: Set initial values
 * @li Step 8: Simulate
 *
 * For more information, we refer to the OpenLB user guide.
 */

/// Include OpenLB library and load analytical solution

#include <olb.h>
#include "analyticalSolutionTestFlow3D.h"

using namespace olb;
using namespace olb::names;

enum class GeometryType : int {
  cube = 0,
  sphere = 1
};

enum class BoundaryType : int {
  interpolated = 0,
  local = 1,
  bouzidi = 2
};

namespace olb::parameters {

struct ERROR_VELOCITY_L1 : public descriptors::FIELD_BASE<1> { };
struct ERROR_VELOCITY_L2 : public descriptors::FIELD_BASE<1> { };
struct ERROR_VELOCITY_LINF : public descriptors::FIELD_BASE<1> { };
struct ERROR_PRESSURE_L1 : public descriptors::FIELD_BASE<1> { };
struct ERROR_PRESSURE_L2 : public descriptors::FIELD_BASE<1> { };
struct ERROR_PRESSURE_LINF : public descriptors::FIELD_BASE<1> { };
struct ERROR_STRAIN_RATE_L1 : public descriptors::FIELD_BASE<1> { };
struct ERROR_STRAIN_RATE_L2 : public descriptors::FIELD_BASE<1> { };
struct ERROR_STRAIN_RATE_LINF : public descriptors::FIELD_BASE<1> { };
struct ERROR_DISSIPATION_L1 : public descriptors::FIELD_BASE<1> { };
struct ERROR_DISSIPATION_L2 : public descriptors::FIELD_BASE<1> { };
struct ERROR_DISSIPATION_LINF : public descriptors::FIELD_BASE<1> { };
struct ERROR_STRESS_L1 : public descriptors::FIELD_BASE<1> { };
struct ERROR_STRESS_L2 : public descriptors::FIELD_BASE<1> { };
struct ERROR_STRESS_LINF : public descriptors::FIELD_BASE<1> { };

struct GEOMETRY_TYPE : public descriptors::TYPED_FIELD_BASE<GeometryType,1> { };
struct BOUNDARY_TYPE : public descriptors::TYPED_FIELD_BASE<BoundaryType,1> { };

}

/// @brief Step 1: Declare simulation structure.
/// Model name and lattice type are collected in a Case class
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q19<descriptors::FORCE>>
>;

/// @brief Step 3: Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  std::shared_ptr<IndicatorF3D<T>> domain, sphere;
  Vector<T,3> extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const T physDeltaX = extent[0] /
                       parameters.get<parameters::RESOLUTION>();
  switch (parameters.get<parameters::GEOMETRY_TYPE>()) {
    case GeometryType::sphere: {
      Vector<T,3> center{extent[0]/2., extent[1]/2., extent[2]/2.};
      sphere = std::make_shared<IndicatorSphere3D<T>>(center, extent[0] / 2.);
      domain = std::make_shared<IndicatorLayer3D<T>>(*sphere, physDeltaX);
      break;
    }
    default: {
      Vector<T,3> origin{0, 0, 0};
      domain = std::make_shared<IndicatorCuboid3D<T>>(extent, origin);
    }
  }
  Mesh<T,MyCase::d> mesh(*domain, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

/// @brief Step 5: Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
/// @note The material numbers are used to assign physics to lattice nodes
void prepareGeometry(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  Vector<T,3> extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const T physDeltaX = extent[0] /
                       parameters.get<parameters::RESOLUTION>();

  switch (parameters.get<parameters::GEOMETRY_TYPE>()) {
    case GeometryType::sphere: {
      Vector<T,3> center{extent[0]/2., extent[1]/2., extent[2]/2.};
      IndicatorSphere3D<T> sphere(center, extent[0] / 2.);
      IndicatorLayer3D<T> sphereLayer(sphere, physDeltaX);
      geometry.rename(0,2,sphereLayer);
      geometry.rename(2,1,sphere);

      geometry.clean();
      geometry.checkForErrors();
      break;
    }
    default: {
      geometry.rename(0, 2);
      geometry.rename(2, 1, {1,1,1});
    }
  }
  geometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

/// @brief Step 6: Set lattice dynamics and boundary conditions
/// @param myCase The Case instance which keeps the simulation data
void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& parameters = myCase.getParameters();

  /// @li Set up a unit converter with the characteristic physical units
  lattice.setUnitConverter<UnitConverterFromResolutionAndLatticeVelocity<MyCase::value_t,
                                                                         MyCase::descriptor_t_of<NavierStokes>>>(
    parameters.get<parameters::RESOLUTION>(),               // resolution
    parameters.get<parameters::LATTICE_CHAR_VELOCITY>(),    // charLatticeVelocity
    parameters.get<parameters::DOMAIN_EXTENT>()[0],         // charPhysLength: reference length of simulation geometry in [m]
    parameters.get<parameters::PHYS_CHAR_VELOCITY>(),       // charPhysVelocity: highest expected velocity during simulation in [m/s]
    parameters.get<parameters::PHYS_CHAR_VISCOSITY>(),      // physViscosity: physical kinematic viscosity in [m^2/s]
    parameters.get<parameters::PHYS_CHAR_DENSITY>()         // physDensity: physical density [kg/m^3]
  );
  auto& converter = lattice.getUnitConverter();
  converter.print();

  /// @li Material=1 --> bulk dynamics
  dynamics::set<ForcedBGKdynamics>(lattice, geometry.getMaterialIndicator({1}));
  /// @li Material=2,3 --> velocity boundary
  switch (parameters.get<parameters::BOUNDARY_TYPE>()) {
    case BoundaryType::interpolated: {
      boundary::set<boundary::InterpolatedVelocity>(lattice, geometry.getMaterialIndicator({2}));
      break;
    }
    case BoundaryType::local: {
      boundary::set<boundary::LocalVelocity>(lattice, geometry.getMaterialIndicator({2}));
      break;
    }
    case BoundaryType::bouzidi: {
      Vector<T,3> extent = parameters.get<parameters::DOMAIN_EXTENT>();
      Vector<T,3> center{extent[0]/2., extent[1]/2., extent[2]/2.};
      IndicatorSphere3D<T> sphere(center, extent[0] / 2.);
      setBouzidiBoundary<T,MyCase::descriptor_t,BouzidiVelocityPostProcessor>(lattice, geometry, 2, sphere);
      break;
    }
    default: {
      std::runtime_error("Invalid Boundary Type");
    }
  }
  /// @li Set lattice relaxation frequency
  lattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
}

// Step 7: Set initial values for primal variables (e.g. velocity, density) and additional fields
/// @param myCase The Case instance which keeps the simulation data
/// @note Initial values have to be set using lattice units
void setInitialValues(MyCase& myCase)
{
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();

  /// @li Initialize density to one everywhere
  /// @li Initialize velocity to be tangential at the lid and zero in the bulk and at the walls
  AnalyticalConst3D<T,T> rhoF(1);
  const Vector<T,3> u{0.,0.,0.};
  AnalyticalConst3D<T,T> uF(u);

  /// @li Set force field
  ForceTestFlow3D<T,T,MyCase::descriptor_t> forceF(converter);
  const T latticeScaling(converter.getConversionFactorMass() / converter.getConversionFactorForce());
  AnalyticalScaled3D<T,T> scaledForceF(forceF, latticeScaling);  // conversion to lattice units
  fields::set<descriptors::FORCE>(lattice, geometry.getMaterialIndicator({1}), scaledForceF);

  /// @li Initialize lattice
  lattice.initialize();
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT) {
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& converter = lattice.getUnitConverter();
  auto& parameters = myCase.getParameters();
  const T physStartT = parameters.get<parameters::PHYS_MAX_T>() * (1.0 / 3.0);
  const T physBoundaryValueUpdateT = parameters.get<parameters::PHYS_BOUNDARY_VALUE_UPDATE_T>();

  const std::size_t itStart = converter.getLatticeTime(physStartT);
  const std::size_t itLatticeUpdate = converter.getLatticeTime(physBoundaryValueUpdateT);
  const std::size_t itUpdate = (itLatticeUpdate == 0) ? 1 : itLatticeUpdate;

  if (iT <= itStart && iT % itUpdate == 0) {
    /// @li Compute scaling factor
    PolynomialStartScale<T,T> startScaleF(itStart, T(1));
    const T iTvec[1] = {T(iT)};
    T frac[1] = {};
    startScaleF(frac, iTvec);

    /// @li Take analytical velocity solution, scale it to lattice units, set the boundary data
    VelocityTestFlow3D<T,T,MyCase::descriptor_t> velocityF(converter);
    AnalyticalScaled3D<T,T> uBoundaryStartF(velocityF, frac[0] / converter.getConversionFactorVelocity());
    AnalyticalScaled3D<T,T> uBoundaryPhysStartF(velocityF, frac[0]);
    switch (parameters.get<parameters::BOUNDARY_TYPE>()) {
      case BoundaryType::bouzidi: {
        setBouzidiVelocity(lattice, geometry, 2, uBoundaryStartF);
        lattice.template setProcessingContext<Array<descriptors::BOUZIDI_VELOCITY>>(ProcessingContext::Simulation);
    break;
      }
      default: {
        momenta::setVelocity(lattice, geometry.getMaterialIndicator({2}), uBoundaryPhysStartF);
        lattice.template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
      }
    }
  }
}

void computeErrorNorms(MyCase& myCase) {
  using T = MyCase::value_t;
  OstreamManager clout(std::cout, "ErrorNorm");
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& converter = sLattice.getUnitConverter();
  sLattice.setProcessingContext(ProcessingContext::Evaluation);

  // Compute error to analytic solution
  // Analytic solutions
  VelocityTestFlow3D<T,T,MyCase::descriptor_t> analytic_vel(converter);
  PressureTestFlow3D<T,T,MyCase::descriptor_t> analytic_pre(converter);
  StrainRateTestFlow3D<T,T,MyCase::descriptor_t> analytic_str(converter);
  DissipationTestFlow3D<T,T,MyCase::descriptor_t> analytic_dis(converter);
  StressTensorTestFlow3D<T,T,MyCase::descriptor_t> analytic_stress(converter);

  // Interpolated analytic solutions
  SuperLatticeFfromAnalyticalF3D<T,MyCase::descriptor_t> lattice_analytical_vel(analytic_vel, sLattice);
  SuperLatticeFfromAnalyticalF3D<T,MyCase::descriptor_t> lattice_analytical_pre(analytic_pre, sLattice);
  SuperLatticeFfromAnalyticalF3D<T,MyCase::descriptor_t> lattice_analytical_str(analytic_str, sLattice);
  SuperLatticeFfromAnalyticalF3D<T,MyCase::descriptor_t> lattice_analytical_dis(analytic_dis, sLattice);
  SuperLatticeFfromAnalyticalF3D<T,MyCase::descriptor_t> lattice_analytical_stress(analytic_stress, sLattice);

  // Simulated solutions
  SuperLatticePhysVelocity3D<T,MyCase::descriptor_t> vel(sLattice, converter);
  SuperLatticePhysPressure3D<T,MyCase::descriptor_t> pre(sLattice, converter);
  SuperLatticePhysStrainRate3D<T,MyCase::descriptor_t> str(sLattice, converter);
  SuperLatticePhysDissipation3D<T,MyCase::descriptor_t> dis(sLattice, converter);
  SuperLatticePhysStressTensor3D<T,MyCase::descriptor_t> stress(sLattice, converter);

  int tmp[4]{0}; T norm[1]{0}; int tmp2[4]{0}; T norm2[1]{0};
  SuperL1Norm3D<T> uL1Norm(lattice_analytical_vel - vel, myCase.getGeometry(), 1);
  SuperL1Norm3D<T> uAnaL1Norm(lattice_analytical_vel, myCase.getGeometry(), 1);
  uL1Norm(norm, tmp); uAnaL1Norm(norm2, tmp2);
  clout << "velocity-L1-error(abs)=" << norm[0] << "; velocity-L1-error(rel)=" << norm[0] / norm2[0] << std::endl;
  myCase.getParameters().set<parameters::ERROR_VELOCITY_L1>(norm[0]);
  SuperL2Norm3D<T> uL2Norm(lattice_analytical_vel - vel, myCase.getGeometry(), 1);
  SuperL2Norm3D<T> uAnaL2Norm(lattice_analytical_vel, myCase.getGeometry(), 1);
  uL2Norm(norm, tmp); uAnaL2Norm(norm2, tmp2);
  clout << "velocity-L2-error(abs)=" << norm[0] << "; velocity-L2-error(rel)=" << norm[0] / norm2[0] << std::endl;
  myCase.getParameters().set<parameters::ERROR_VELOCITY_L2>(norm[0]);
  SuperLinfNorm3D<T> uLinfNorm(lattice_analytical_vel - vel, myCase.getGeometry(), 1);
  SuperLinfNorm3D<T> uAnaLinfNorm(lattice_analytical_vel, myCase.getGeometry(), 1);
  uLinfNorm(norm, tmp); uAnaLinfNorm(norm2, tmp2);
  clout << "velocity-Linf-error(abs)=" << norm[0] << "; velocity-Linf-error(rel)=" << norm[0] / norm2[0] << std::endl;
  myCase.getParameters().set<parameters::ERROR_VELOCITY_LINF>(norm[0]);
  SuperL1Norm3D<T> pL1Norm(lattice_analytical_pre - pre, myCase.getGeometry(), 1);
  SuperL1Norm3D<T> pAnaL1Norm(lattice_analytical_pre, myCase.getGeometry(), 1);
  pL1Norm(norm, tmp); pAnaL1Norm(norm2, tmp2);
  clout << "pressure-L1-error(abs)=" << norm[0] << "; pressure-L1-error(rel)=" << norm[0] / norm2[0] << std::endl;
  myCase.getParameters().set<parameters::ERROR_PRESSURE_L1>(norm[0]);
  SuperL2Norm3D<T> pL2Norm(lattice_analytical_pre - pre, myCase.getGeometry(), 1);
  SuperL2Norm3D<T> pAnaL2Norm(lattice_analytical_pre, myCase.getGeometry(), 1);
  pL2Norm(norm, tmp); pAnaL2Norm(norm2, tmp2);
  clout << "pressure-L2-error(abs)=" << norm[0] << "; pressure-L2-error(rel)=" << norm[0] / norm2[0] << std::endl;
  myCase.getParameters().set<parameters::ERROR_PRESSURE_L2>(norm[0]);
  SuperLinfNorm3D<T> pLinfNorm(lattice_analytical_pre - pre, myCase.getGeometry(), 1);
  SuperLinfNorm3D<T> pAnaLinfNorm(lattice_analytical_pre, myCase.getGeometry(), 1);
  pLinfNorm(norm, tmp); pAnaLinfNorm(norm2, tmp2);
  clout << "pressure-Linf-error(abs)=" << norm[0] << "; pressure-Linf-error(rel)=" << norm[0] / norm2[0] << std::endl;
  myCase.getParameters().set<parameters::ERROR_PRESSURE_LINF>(norm[0]);
  SuperL1Norm3D<T> sL1Norm(lattice_analytical_str - str, myCase.getGeometry(), 1);
  SuperL1Norm3D<T> sAnaL1Norm(lattice_analytical_str, myCase.getGeometry(), 1);
  sL1Norm(norm, tmp); sAnaL1Norm(norm2, tmp2);
  clout << "strainRate-L1-error(abs)=" << norm[0] << "; strainRate-L1-error(rel)=" << norm[0] / norm2[0] << std::endl;
  myCase.getParameters().set<parameters::ERROR_STRAIN_RATE_L1>(norm[0]);
  SuperL2Norm3D<T> sL2Norm(lattice_analytical_str - str, myCase.getGeometry(), 1);
  SuperL2Norm3D<T> sAnaL2Norm(lattice_analytical_str, myCase.getGeometry(), 1);
  sL2Norm(norm, tmp); sAnaL2Norm(norm2, tmp2);
  clout << "strainRate-L2-error(abs)=" << norm[0] << "; strainRate-L2-error(rel)=" << norm[0] / norm2[0] << std::endl;
  myCase.getParameters().set<parameters::ERROR_STRAIN_RATE_L2>(norm[0]);
  SuperLinfNorm3D<T> sLinfNorm(lattice_analytical_str - str, myCase.getGeometry(), 1);
  SuperLinfNorm3D<T> sAnaLinfNorm(lattice_analytical_str, myCase.getGeometry(), 1);
  sLinfNorm(norm, tmp); sAnaLinfNorm(norm2, tmp2);
  clout << "strainRate-Linf-error(abs)=" << norm[0] << "; strainRate-Linf-error(rel)" << norm[0] / norm2[0] << std::endl;
  myCase.getParameters().set<parameters::ERROR_STRAIN_RATE_LINF>(norm[0]);
  SuperL1Norm3D<T> dL1Norm(lattice_analytical_dis - dis, myCase.getGeometry(), 1);
  SuperL1Norm3D<T> dAnaL1Norm(lattice_analytical_dis, myCase.getGeometry(), 1);
  dL1Norm(norm, tmp); dAnaL1Norm(norm2, tmp2);
  clout << "dissipation-L1-error(abs)=" << norm[0] << "; dissipation-L1-error(rel)=" << norm[0] / norm2[0] << std::endl;
  myCase.getParameters().set<parameters::ERROR_DISSIPATION_L1>(norm[0]);
  SuperL2Norm3D<T> dL2Norm(lattice_analytical_dis - dis, myCase.getGeometry(), 1);
  SuperL2Norm3D<T> dAnaL2Norm(lattice_analytical_dis, myCase.getGeometry(), 1);
  dL2Norm(norm, tmp); dAnaL2Norm(norm2, tmp2);
  clout << "dissipation-L2-error(abs)=" << norm[0] << "; dissipation-L2-error(rel)=" << norm[0] / norm2[0] << std::endl;
  myCase.getParameters().set<parameters::ERROR_DISSIPATION_L2>(norm[0]);
  SuperLinfNorm3D<T> dLinfNorm(lattice_analytical_dis - dis, myCase.getGeometry(), 1);
  SuperLinfNorm3D<T> dAnaLinfNorm(lattice_analytical_dis, myCase.getGeometry(), 1);
  dLinfNorm(norm, tmp); dAnaLinfNorm(norm2, tmp2);
  clout << "dissipation-Linf-error(abs)=" << norm[0] << "; dissipation-Linf-error(rel)=" << norm[0] / norm2[0] << std::endl;
  myCase.getParameters().set<parameters::ERROR_DISSIPATION_LINF>(norm[0]);
  SuperL1Norm3D<T> stL1Norm(lattice_analytical_stress - stress, myCase.getGeometry(), 1);
  SuperL1Norm3D<T> stAnaL1Norm(lattice_analytical_stress, myCase.getGeometry(), 1);
  stL1Norm(norm, tmp); stAnaL1Norm(norm2, tmp2);
  clout << "stress-L1-error(abs)=" << norm[0] << "; stress-L1-error(rel)=" << norm[0] / norm2[0] << std::endl;
  myCase.getParameters().set<parameters::ERROR_STRESS_L1>(norm[0]);
  SuperL2Norm3D<T> stL2Norm(lattice_analytical_stress - stress, myCase.getGeometry(), 1);
  SuperL2Norm3D<T> stAnaL2Norm(lattice_analytical_stress, myCase.getGeometry(), 1);
  stL2Norm(norm, tmp); stAnaL2Norm(norm2, tmp2);
  clout << "stress-L2-error(abs)=" << norm[0] << "; stress-L2-error(rel)=" << norm[0] / norm2[0] << std::endl;
  myCase.getParameters().set<parameters::ERROR_STRESS_L2>(norm[0]);
  SuperLinfNorm3D<T> stLinfNorm(lattice_analytical_stress - stress, myCase.getGeometry(), 1);
  SuperLinfNorm3D<T> stAnaLinfNorm(lattice_analytical_stress, myCase.getGeometry(), 1);
  stLinfNorm(norm, tmp); stAnaLinfNorm(norm2, tmp2);
  clout << "stress-Linf-error(abs)=" << norm[0] << "; stress-Linf-error(rel)=" << norm[0] / norm2[0] << std::endl;
  myCase.getParameters().set<parameters::ERROR_STRESS_LINF>(norm[0]);
}

/// Step 8.3: Compute simulation results at times
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();
  const T physMaxT = myCase.getParameters().get<parameters::PHYS_MAX_T>();
  const std::size_t iTvtk = converter.getLatticeTime(physMaxT/1);
  const std::size_t iTlog = converter.getLatticeTime(physMaxT/10);

  SuperVTMwriter3D<T> vtmWriter("testFlow3d");
  SuperLatticePhysVelocity3D velocityF(lattice, converter);
  SuperLatticePhysPressure3D pressureF(lattice, converter);
  vtmWriter.addFunctor(velocityF);
  vtmWriter.addFunctor(pressureF);

  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  // Writes the VTK files
  if (iT%iTvtk == 0 && iT > 0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }

  /// @li Print some (numerical and computational) statistics
  if (iT%iTlog == 0) {
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
    computeErrorNorms(myCase);
    timer.print(iT);
  }
}

/// Step 8: Execute simulation
/// @param myCase The Case instance which keeps the simulation data
/// Run time loop
void simulate(MyCase& myCase)
{
  using T = MyCase::value_t;
  const T physMaxT = myCase.getParameters().get<parameters::PHYS_MAX_T>();
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& converter = sLattice.getUnitConverter();
  const std::size_t iTmax = converter.getLatticeTime(physMaxT);
  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
#ifndef PLATFORM_GPU_CUDA
  util::ValueTracer<T> converge(converter.getLatticeTime(1.0), 1e-8);
#endif
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT)
  {
#ifndef PLATFORM_GPU_CUDA
    if(converge.hasConverged()) {
      getResults(myCase, timer, iT);
      computeErrorNorms(myCase);
      break;
    }
#endif
    setTemporalValues(myCase, iT);

    /// @li Step 8.2: Collide and Stream Execution
    sLattice.collideAndStream();

    /// @li Stripe off density offset due to Dirichlet boundary conditions
    sLattice.stripeOffDensityOffset(
      sLattice.getStatistics().getAverageRho() - T{1});

    /// @li Step 8.3: Computation and Output of the Results
    getResults(myCase, timer, iT);
#ifndef PLATFORM_GPU_CUDA
    converge.takeValue(sLattice.getStatistics().getAverageEnergy(), true);
#endif
  }

  /// @li Evaluate timer
  timer.stop();
  timer.printSummary();
}
