/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2018 Marc Haußmann, Mathias J. Krause, Jonathan Jeppener-Haltenhoff
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

/* poiseuille3d.cpp:
 * This example examines a 3D Poseuille flow
 * It illustrates the usage of different flow setups, boundary conditions
 * and computation of error norms.
 * The following simulation options can be combined freely:
 * - if the compiler flag ENABLE_MRT is set, mrt collision operators are used
 * - forced/ nonForced flow
 * - different boundary conditions
 * - simulation only or eoc convergence analysis
 *
 * The main code of the simulation is in case.h as it is also used by the
 * example ../../pdeSolverEoc/poiseuille3dEoc
 *
 * Set flag in order to use mrt collision operators instead of bgk
 * #define ENABLE_MRT
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;

enum class FlowType : int {
  FORCED = 0,
  NON_FORCED = 1
};

enum class BoundaryType : int {
  BOUNCE_BACK = 0,
  LOCAL = 1,
  INTERPOLATED = 2,
  BOUZIDI= 3 ,
  FREE_SLIP = 4,
  PARTIAL_SLIP = 5
};

namespace olb::parameters {

struct FLOW_TYPE      : public descriptors::TYPED_FIELD_BASE<FlowType,1> { };
struct BOUNDARY_TYPE  : public descriptors::TYPED_FIELD_BASE<BoundaryType,1> { };
struct PARTIAL_SLIP_TUNER   : public descriptors::TYPED_FIELD_BASE<int,1> { };
struct EOC            : public descriptors::TYPED_FIELD_BASE<bool,1> { };

struct VELOCITY_L1_ABS_ERROR    : public descriptors::FIELD_BASE<1> { };
struct VELOCITY_L2_ABS_ERROR    : public descriptors::FIELD_BASE<1> { };
struct VELOCITY_LINF_ABS_ERROR  : public descriptors::FIELD_BASE<1> { };
struct STRAIN_RATE_L1_ABS_ERROR    : public descriptors::FIELD_BASE<1> { };
struct STRAIN_RATE_L2_ABS_ERROR    : public descriptors::FIELD_BASE<1> { };
struct STRAIN_RATE_LINF_ABS_ERROR  : public descriptors::FIELD_BASE<1> { };
struct WSS_L1_ABS_ERROR    : public descriptors::FIELD_BASE<1> { };
struct WSS_L2_ABS_ERROR    : public descriptors::FIELD_BASE<1> { };
struct WSS_LINF_ABS_ERROR  : public descriptors::FIELD_BASE<1> { };
struct PRESSURE_L1_ABS_ERROR    : public descriptors::FIELD_BASE<1> { };
struct PRESSURE_L2_ABS_ERROR    : public descriptors::FIELD_BASE<1> { };
struct PRESSURE_LINF_ABS_ERROR  : public descriptors::FIELD_BASE<1> { };

struct DIAMETER       : public descriptors::FIELD_BASE<1> { };  // diameter of the pipe
struct LENGTH         : public descriptors::FIELD_BASE<1> { };  // length of the pipe

}

using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q19<descriptors::FORCE>>
>;

Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  OstreamManager clout(std::cout, "prepareGeometry");
  const T radius = parameters.get<parameters::DIAMETER>()/2;
  const T length = parameters.get<parameters::LENGTH>();
  Vector<T, 3>            center0(0, radius, radius);
  Vector<T, 3>            center1(length + 0.5 * parameters.get<parameters::PHYS_DELTA_X>(), radius, radius);
  IndicatorCylinder3D<T>  pipe(center0, center1, radius);
  IndicatorLayer3D<T>     extendedDomain(pipe, parameters.get<parameters::PHYS_DELTA_X>());

  Mesh<T, MyCase::d> mesh(extendedDomain, parameters.get<parameters::PHYS_DELTA_X>(), singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  if (parameters.get<parameters::FLOW_TYPE>() == FlowType::FORCED) {
    mesh.getCuboidDecomposition().setPeriodicity({true, false, false});
  }
  return mesh;
}

// Stores geometry information in form of material numbers
void prepareGeometry(MyCase& myCase)
{
  using T = MyCase::value_t;
  OstreamManager clout(std::cout, "prepareGeometry");

  clout << "Prepare Geometry ..." << std::endl;

  auto&parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();

  const T radius = parameters.get<parameters::DIAMETER>()/2.0;
  const T length = parameters.get<parameters::LENGTH>();
  FlowType flowType = parameters.get<parameters::FLOW_TYPE>();

  const T dx = parameters.get<parameters::PHYS_DELTA_X>();

  Vector<T, 3> center0(-dx * 0.2, radius, radius);
  Vector<T, 3> center1(length, radius, radius);
  if (flowType == FlowType::FORCED) {
    center0[0] -= 3.*dx;
    center1[0] += 3.*dx;
  }
  IndicatorCylinder3D<T> pipe(center0, center1, radius);

  geometry.rename(0, 2);

  geometry.rename(2, 1, pipe);

  if (flowType == FlowType::NON_FORCED) {
    geometry.clean();
    Vector<T, 3> origin(0, radius, radius);
    Vector<T, 3> extend = origin;

    // Set material number for inflow
    origin[0] = -dx * 2;
    extend[0] = dx * 2;
    IndicatorCylinder3D<T> inflow(origin, extend, radius);
    geometry.rename(2, 3, 1, inflow);

    // Set material number for outflow
    origin[0] = length - 2 * dx;
    extend[0] = length + 2 * dx;
    IndicatorCylinder3D<T> outflow(extend, origin, radius);
    geometry.rename(2, 4, 1, outflow);
  }

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  // Removes all not needed boundary voxels inside the surface
  geometry.innerClean();
  geometry.checkForErrors();

  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(MyCase& myCase)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});

  const T radius = parameters.get<parameters::DIAMETER>()/2.0;
  const T length = parameters.get<parameters::LENGTH>();

  const T diameter = parameters.get<parameters::DIAMETER>();
  const T physU = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();

  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    int {parameters.get<parameters::RESOLUTION>()},                  // resolution: number of voxels per charPhysL
    (T)   parameters.get<parameters::LATTICE_RELAXATION_TIME>(), // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   diameter,           // charPhysLength: reference length of simulation geometry
    (T)   physU,// charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   diameter*physU/parameters.get<parameters::REYNOLDS>(),  // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   parameters.get<parameters::PHYS_CHAR_DENSITY>() // physDensity: physical density in __kg / m^3__
  );
  auto& converter = lattice.getUnitConverter();
  converter.print();

  FlowType      flowType      = parameters.get<parameters::FLOW_TYPE>();
  BoundaryType  boundaryType  = parameters.get<parameters::BOUNDARY_TYPE>();

  if (flowType == FlowType::NON_FORCED) {
    dynamics::set<BGKdynamics>(lattice, myCase.getGeometry(), 1);
  } else {
    dynamics::set<ForcedBGKdynamics>(lattice, myCase.getGeometry(), 1);
  }

  Vector<T, 3> center0(0, radius, radius);
  Vector<T, 3> center1(length, radius, radius);
  if (boundaryType == BoundaryType::BOUZIDI) {
    center0[0] -= 0.5*physDeltaX;
    center1[0] += 0.5*physDeltaX;
    if (flowType == FlowType::FORCED) {
      center0[0] -= 3.*physDeltaX;
      center1[0] += 3.*physDeltaX;
    }
  }
  IndicatorCylinder3D<T> pipe(center0, center1, radius);

  switch (boundaryType) {
    case BoundaryType::BOUNCE_BACK:
      boundary::set<boundary::BounceBack>(lattice, myCase.getGeometry(), 2);
      break;
    case BoundaryType::FREE_SLIP:
      boundary::set<boundary::FullSlip>(lattice, myCase.getGeometry(), 2);
      break;
    case BoundaryType::PARTIAL_SLIP:
      boundary::set<boundary::PartialSlip>(lattice, myCase.getGeometry(), 2);
      lattice.template setParameter<descriptors::TUNER>(parameters.get<parameters::PARTIAL_SLIP_TUNER>());
      break;
    case BoundaryType::BOUZIDI:
      setBouzidiBoundary<T, DESCRIPTOR, BouzidiPostProcessor>(lattice, myCase.getGeometry(), 2, pipe);
      break;
    case BoundaryType::LOCAL:
      boundary::set<boundary::LocalVelocity>(lattice, myCase.getGeometry(), 2);
      break;
    case BoundaryType::INTERPOLATED:
    default:
      boundary::set<boundary::InterpolatedVelocity>(lattice, myCase.getGeometry(), 2);
      break;
  }

  if (flowType == FlowType::NON_FORCED) {
    switch (boundaryType) {
      case BoundaryType::BOUZIDI:
        setBouzidiBoundary<T, DESCRIPTOR, BouzidiVelocityPostProcessor>(lattice, myCase.getGeometry(), 3, pipe);
        break;
      case BoundaryType::LOCAL:
        boundary::set<boundary::LocalVelocity>(lattice, myCase.getGeometry(), 3);
        break;
      case BoundaryType::BOUNCE_BACK:
      case BoundaryType::INTERPOLATED:
      case BoundaryType::FREE_SLIP:
      case BoundaryType::PARTIAL_SLIP:
      default:
        boundary::set<boundary::InterpolatedVelocity>(lattice, myCase.getGeometry(), 3);
        break;
    }
    if (boundaryType == BoundaryType::LOCAL) {
      boundary::set<boundary::LocalPressure>(lattice, myCase.getGeometry(), 4);
    }
    else {
      boundary::set<boundary::InterpolatedPressure>(lattice, myCase.getGeometry(), 4);
    }
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  OstreamManager clout( std::cout,"setInitialValues" );
  clout << "setInitialValues ..." << std::endl;
  using T             = MyCase::value_t;
  auto& parameters    = myCase.getParameters();
  auto& lattice       = myCase.getLattice(NavierStokes{});
  auto& geometry      = myCase.getGeometry();
  auto& converter     = lattice.getUnitConverter();

  const T length      = parameters.get<parameters::LENGTH>();
  const T diameter    = parameters.get<parameters::DIAMETER>();
  const T physU       = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T radius      = diameter / 2.;

  BoundaryType boundaryType = parameters.get<parameters::BOUNDARY_TYPE>();
  FlowType flowType   = parameters.get<parameters::FLOW_TYPE>();

  std::vector<T> origin = { length, radius, radius};
  std::vector<T> axis = { 1, 0, 0 };
  CirclePoiseuille3D<T> poiseuilleU(origin, axis, lattice.getUnitConverter().getCharPhysVelocity(), radius);

  if (flowType == FlowType::FORCED) {
    const T D = lattice.getUnitConverter().getLatticeLength(diameter);
    Vector<T,3> poiseuilleF = 0;
    poiseuilleF[0] = 4*lattice.getUnitConverter().getLatticeViscosity()
                   * lattice.getUnitConverter().getCharLatticeVelocity()
                   / (D * D / 4);
    fields::set<descriptors::FORCE>(lattice, myCase.getGeometry().getMaterialIndicator({1,2}), poiseuilleF);

    momenta::setVelocity(lattice, myCase.getGeometry().getMaterialIndicator({1,2}), poiseuilleU);
  } else {
    T p0 = 4. * lattice.getUnitConverter().getPhysViscosity() * physU * length / (radius * radius);

    AnalyticalLinear3D<T,T> pressureF(-p0 / length, 0, 0, p0);

    momenta::setPressureAndVelocity(lattice,
                                    geometry.getMaterialIndicator({0,1,2,3}),
                                    pressureF,
                                    poiseuilleU);


    if (boundaryType == BoundaryType::BOUZIDI) {
      setBouzidiVelocity(lattice, myCase.getGeometry(), 3, poiseuilleU);
    }
  }

  lattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  lattice.initialize();

  clout << "setInitialValues ... OK" << std::endl;
}

SuperLatticePhysWallShearStress3D<MyCase::value_t,MyCase::descriptor_t_of<NavierStokes>> createWallShearStressF(MyCase& myCase)
{
  using T           = MyCase::value_t;
  using DESCRIPTOR  = MyCase::descriptor_t_of<NavierStokes>;
  auto& lattice     = myCase.getLattice(NavierStokes{});
  auto& geometry    = myCase.getGeometry();
  auto& parameters  = myCase.getParameters();
  auto& converter   = lattice.getUnitConverter();

  const T radius    = parameters.get<parameters::DIAMETER>()/2.0;
  const T length    = parameters.get<parameters::LENGTH>();

  // set up size-increased indicator and instantiate wall shear stress functor (wss)
  const T center0x = (parameters.get<parameters::FLOW_TYPE>() == FlowType::FORCED) ? -converter.getPhysDeltaX() * 4.2 : -converter.getPhysDeltaX() * 0.2;
  const Vector<T, 3> center0Extended(center0x, radius, radius);
  const T center1x = (parameters.get<parameters::FLOW_TYPE>() == FlowType::FORCED) ? length + 4.*converter.getPhysDeltaX() : length;
  const Vector<T, 3> center1Extended(center1x, radius, radius);
  IndicatorCylinder3D<T> pipeExtended(center0Extended, center1Extended, radius);
  IndicatorLayer3D<T> indicatorExtended (pipeExtended, 0.9*converter.getPhysDeltaX()*parameters.get<parameters::RESOLUTION>()/11.);
  return SuperLatticePhysWallShearStress3D<T,DESCRIPTOR>(lattice, geometry, 2, converter, indicatorExtended);
}

void calculateError(MyCase& myCase)
{
  OstreamManager clout( std::cout,"error" );
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();

  const FlowType  flowType  = parameters.get<parameters::FLOW_TYPE>();
  const T         radius    = parameters.get<parameters::DIAMETER>()/2.0;
  const T         length    = parameters.get<parameters::LENGTH>();

  lattice.setProcessingContext(ProcessingContext::Evaluation);

  // velocity error
  const T maxVelocity = converter.getCharPhysVelocity();
  const std::vector<T> axisPoint = {length, radius, radius};
  const std::vector<T> axisDirection = { 1, 0, 0 };
  CirclePoiseuille3D<T> uSol(axisPoint, axisDirection, maxVelocity, radius);
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> u( lattice,converter );
  auto indicatorF = geometry.getMaterialIndicator(1);

  int tmp[]= { };
  T result[2]= { };

  SuperAbsoluteErrorL1Norm3D<T> absVelocityErrorNormL1(u, uSol, indicatorF);
  absVelocityErrorNormL1(result, tmp);
  clout << "velocity-L1-error(abs)=" << result[0];
  parameters.set<parameters::VELOCITY_L1_ABS_ERROR>(result[0]);
  SuperRelativeErrorL1Norm3D<T> relVelocityErrorNormL1(u, uSol, indicatorF);
  relVelocityErrorNormL1(result, tmp);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm3D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, tmp);
  clout << "velocity-L2-error(abs)=" << result[0];
  parameters.set<parameters::VELOCITY_L2_ABS_ERROR>(result[0]);
  SuperRelativeErrorL2Norm3D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, tmp);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm3D<T> absVelocityErrorNormLinf(u, uSol, indicatorF);
  absVelocityErrorNormLinf(result, tmp);
  clout << "velocity-Linf-error(abs)=" << result[0];
  parameters.set<parameters::VELOCITY_LINF_ABS_ERROR>(result[0]);
  SuperRelativeErrorLinfNorm3D<T> relVelocityErrorNormLinf(u, uSol, indicatorF);
  relVelocityErrorNormLinf(result, tmp);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;

  // strainRate error
  CirclePoiseuilleStrainRate3D<T, DESCRIPTOR> sSol( converter, radius );
  SuperLatticePhysStrainRate3D<T,DESCRIPTOR> s( lattice,converter );

  SuperAbsoluteErrorL1Norm3D<T> absStrainRateErrorNormL1(s, sSol, indicatorF);
  absStrainRateErrorNormL1(result, tmp);
  clout << "strainRate-L1-error(abs)=" << result[0];
  parameters.set<parameters::STRAIN_RATE_L1_ABS_ERROR>(result[0]);
  SuperRelativeErrorL1Norm3D<T> relStrainRateErrorNormL1(s, sSol, indicatorF);
  relStrainRateErrorNormL1(result, tmp);
  clout << "; strainRate-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm3D<T> absStrainRateErrorNormL2(s, sSol, indicatorF);
  absStrainRateErrorNormL2(result, tmp);
  clout << "strainRate-L2-error(abs)=" << result[0];
  parameters.set<parameters::STRAIN_RATE_L2_ABS_ERROR>(result[0]);
  SuperRelativeErrorL2Norm3D<T> relStrainRateErrorNormL2(s, sSol, indicatorF);
  relStrainRateErrorNormL2(result, tmp);
  clout << "; strainRate-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm3D<T> absStrainRateErrorNormLinf(s, sSol, indicatorF);
  absStrainRateErrorNormLinf(result, tmp);
  clout << "strainRate-Linf-error(abs)=" << result[0];
  parameters.set<parameters::STRAIN_RATE_LINF_ABS_ERROR>(result[0]);
  SuperRelativeErrorLinfNorm3D<T> relStrainRateErrorNormLinf(s, sSol, indicatorF);
  relStrainRateErrorNormLinf(result, tmp);
  clout << "; strainRate-Linf-error(rel)=" << result[0] << std::endl;

  // wallShearStress error
  AnalyticalConst3D<T,T> wssSol(4. * converter.getPhysViscosity() * converter.getPhysDensity() * maxVelocity / parameters.get<parameters::DIAMETER>());
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> wssSolLattice (wssSol, lattice);
  SuperLatticePhysWallShearStress3D<T,DESCRIPTOR> wss = createWallShearStressF(myCase);

  auto indicatorB = geometry.getMaterialIndicator(2);

  SuperAbsoluteErrorL1Norm3D<T> absWallShearStressErrorNormL1(wss, wssSol, indicatorB);
  absWallShearStressErrorNormL1(result, tmp);
  clout << "wss-L1-error(abs)=" << result[0];
  parameters.set<parameters::WSS_L1_ABS_ERROR>(result[0]);
  SuperRelativeErrorL1Norm3D<T> relWallShearStressErrorNormL1(wss, wssSol, indicatorB);
  relWallShearStressErrorNormL1(result, tmp);
  clout << "; wss-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm3D<T> absWallShearStressErrorNormL2(wss, wssSol, indicatorB);
  absWallShearStressErrorNormL2(result, tmp);
  clout << "wss-L2-error(abs)=" << result[0];
  parameters.set<parameters::WSS_L2_ABS_ERROR>(result[0]);
  SuperRelativeErrorL2Norm3D<T> relWallShearStressErrorNormL2(wss, wssSol, indicatorB);
  relWallShearStressErrorNormL2(result, tmp);
  clout << "; wss-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm3D<T> absWallShearStressErrorNormLinf(wss, wssSol, indicatorB);
  absWallShearStressErrorNormLinf(result, tmp);
  clout << "wss-Linf-error(abs)=" << result[0];
  parameters.set<parameters::WSS_LINF_ABS_ERROR>(result[0]);
  SuperRelativeErrorLinfNorm3D<T> relWallShearStressErrorNormLinf(wss, wssSol, indicatorB);
  relWallShearStressErrorNormLinf(result, tmp);
  clout << "; wss-Linf-error(rel)=" << result[0] << std::endl;

  if (flowType == FlowType::NON_FORCED) {
    // pressure error
    T p0 = 4. * converter.getPhysViscosity() * maxVelocity * length / (radius * radius);
    AnalyticalLinear3D<T, T> pressureSol(-p0 / length, 0, 0, p0);
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(lattice, converter);

    SuperAbsoluteErrorL1Norm3D<T> absPressureErrorNormL1(pressure, pressureSol, indicatorF);
    absPressureErrorNormL1(result, tmp);
    clout << "pressure-L1-error(abs)=" << result[0];
    parameters.set<parameters::PRESSURE_L1_ABS_ERROR>(result[0]);
    SuperRelativeErrorL1Norm3D<T> relPressureErrorNormL1(pressure, pressureSol, indicatorF);
    relPressureErrorNormL1(result, tmp);
    clout << "; pressure-L1-error(rel)=" << result[0] << std::endl;

    SuperAbsoluteErrorL2Norm3D<T> absPressureErrorNormL2(pressure, pressureSol, indicatorF);
    absPressureErrorNormL2(result, tmp);
    clout << "pressure-L2-error(abs)=" << result[0];
    parameters.set<parameters::PRESSURE_L2_ABS_ERROR>(result[0]);
    SuperRelativeErrorL2Norm3D<T> relPressureErrorNormL2(pressure, pressureSol, indicatorF);
    relPressureErrorNormL2(result, tmp);
    clout << "; pressure-L2-error(rel)=" << result[0] << std::endl;

    SuperAbsoluteErrorLinfNorm3D<T> absPressureErrorNormLinf(pressure, pressureSol, indicatorF);
    absPressureErrorNormLinf(result, tmp);
    clout << "pressure-Linf-error(abs)=" << result[0];
    parameters.set<parameters::PRESSURE_LINF_ABS_ERROR>(result[0]);
    SuperRelativeErrorLinfNorm3D<T> relPressureErrorNormLinf(pressure, pressureSol, indicatorF);
    relPressureErrorNormLinf(result, tmp);
    clout << "; pressure-Linf-error(rel)=" << result[0] << std::endl;
  }
}

void getResults(MyCase& myCase, std::size_t iT, util::Timer<MyCase::value_t>& timer, bool converged=false)
{
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

  OstreamManager clout( std::cout,"getResults" );
  static Gnuplot<T> gplot("centerVelocity");

  auto& parameters  = myCase.getParameters();
  auto& geometry    = myCase.getGeometry();
  auto& lattice     = myCase.getLattice(NavierStokes{});
  auto& converter   = lattice.getUnitConverter();

  const bool    eoc             = parameters.get<parameters::EOC>();
  const BoundaryType boundaryType = parameters.get<parameters::BOUNDARY_TYPE>();
  const bool    noslipBoundary  = (   (boundaryType != BoundaryType::FREE_SLIP)
                                   && (boundaryType != BoundaryType::PARTIAL_SLIP));

  const size_t  statIter        = converter.getLatticeTime( parameters.get<parameters::PHYS_STAT_ITER_T>() );
  const int     vtmIter         = converter.getLatticeTime( parameters.get<parameters::PHYS_VTK_ITER_T>() );
  bool          lastTimeStep    = (converged || (iT + 1 == converter.getLatticeTime( parameters.get<parameters::PHYS_MAX_T>() )));

  const T       radius          = parameters.get<parameters::DIAMETER>()/2.0;
  const T       length          = parameters.get<parameters::LENGTH>();

  // VTK and image output only if no EOC analysis
  if (!eoc) {
    // set up size-increased indicator and instantiate wall shear stress functor (wss)
    SuperLatticePhysWallShearStress3D<T,DESCRIPTOR> wss = createWallShearStressF(myCase);

    SuperVTMwriter3D<T> vtmWriter( "poiseuille3d" );
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( lattice, converter );
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( lattice, converter );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    vtmWriter.addFunctor( wss );

    const T maxVelocity = converter.getCharPhysVelocity();
    std::vector<T> axisPoint = {length, radius, radius};
    std::vector<T> axisDirection = { 1, 0, 0 };
    CirclePoiseuille3D<T> uSol(axisPoint, axisDirection, maxVelocity, radius);
    SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> analyticalVelocityLattice(uSol, lattice);
    analyticalVelocityLattice.getName() = "analytical solution";
    vtmWriter.addFunctor(analyticalVelocityLattice);

    if ( iT==0 ) {
      // Writes the geometry, cuboid no. and rank no. as vti file for visualization
      SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( lattice );
      SuperLatticeRank3D<T, DESCRIPTOR> rank( lattice );
      SuperLatticeDiscreteNormal3D<T, DESCRIPTOR> discreteNormal( lattice, geometry, geometry.getMaterialIndicator({2, 3}) );
      SuperLatticeDiscreteNormalType3D<T, DESCRIPTOR> discreteNormalType( lattice, geometry, geometry.getMaterialIndicator({2, 3, 4, 5}) );
      vtmWriter.write( cuboid );
      vtmWriter.write( rank );
      vtmWriter.write( discreteNormal );
      vtmWriter.write( discreteNormalType );

      vtmWriter.createMasterFile();
    }

    // Writes the vtm files and profile text file
    if ( iT%vtmIter==0 || lastTimeStep ) {
      lattice.setProcessingContext(ProcessingContext::Evaluation);

      vtmWriter.write( iT );

      SuperEuklidNorm3D<T> normVel( velocity );
      BlockReduction3D2D<T> planeReduction( normVel, Vector<T,3>{0,0,1}, 600, BlockDataSyncMode::ReduceOnly );
      // write output as JPEG
      heatmap::write(planeReduction, iT);
    }
  }

  // Writes output on the console
  if ( iT%statIter==0 || lastTimeStep ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    lattice.getStatistics().printWithoutPlot( iT,converter.getPhysTime( iT ) );

    // Error norms
    if (noslipBoundary) {
      if ( (!eoc) || lastTimeStep ) {
        calculateError(myCase);
      }
    }
  }

  // Gnuplot output
  if ((noslipBoundary) && (lastTimeStep) && (!eoc)) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    // plot velocity magnitude over line through the center of the simulation domain
    const T maxVelocity = converter.getCharPhysVelocity();
    T D = converter.getLatticeLength( parameters.get<parameters::DIAMETER>() );
    T dx = 1. / T(converter.getResolution());
    T point[3] { };
    point[0] = length/2.;
    point[2] = ( T )radius;
    std::vector<T> axisPoint {length, radius, radius};
    std::vector<T> axisDirection { 1, 0, 0 };
    CirclePoiseuille3D<T> uSol(axisPoint, axisDirection, maxVelocity, radius);
    T analytical[3] { };
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( lattice, converter );
    AnalyticalFfromSuperF3D<T> intpolateVelocity( velocity, true, 1 );
    T numerical[3] { };
    for ( int iY=0; iY<=D; ++iY ) {
      point[1] = ( T )converter.getPhysLength(iY);
      uSol( analytical,point );
      intpolateVelocity( numerical,point );
      gplot.setData( iY*dx, {analytical[0],numerical[0]}, {"analytical","numerical"} );
    }
    // Create PNG file
    gplot.writePNG();
  }
}

void simulate(MyCase& myCase) {
  OstreamManager clout( std::cout,"simulate" );
  using T           = MyCase::value_t;
  auto& parameters  = myCase.getParameters();
  auto& lattice     = myCase.getLattice(NavierStokes{});
  auto& geometry    = myCase.getGeometry();
  auto& converter   = lattice.getUnitConverter();

  const size_t  iTmax   = converter.getLatticeTime( parameters.get<parameters::PHYS_MAX_T>() );
  const size_t  iTcheck = converter.getLatticeTime( parameters.get<parameters::PHYS_CONVERGENCE_CHECK_T>() );
  const T       epsilon = parameters.get<parameters::RESIDUUM>();

  clout << "starting simulation..." << std::endl;
  util::Timer<T>        timer(iTmax, geometry.getStatistics().getNvoxel());
  util::ValueTracer<T>  converge( iTcheck, epsilon );
  timer.start();

  for ( std::size_t iT = 0; iT < iTmax; ++iT ) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << std::endl;
      getResults( myCase, iT, timer, true );
      break;
    }

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    // in this application no boundary conditions have to be adjusted

    // === 6th Step: Collide and Stream Execution ===
    lattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( myCase, iT, timer );
    converge.takeValue( lattice.getStatistics().getAverageEnergy(), false );
  }

  timer.stop();
  timer.printSummary();
}

void setDefaultParameters(MyCase::ParametersD& params) {
  using namespace olb::parameters;
  params.set<FLOW_TYPE               >(FlowType::NON_FORCED);
  params.set<BOUNDARY_TYPE           >(BoundaryType::INTERPOLATED);
  params.set<RESOLUTION              >(           21 );
  params.set<PHYS_CHAR_VELOCITY      >(           1. );
  params.set<PHYS_CHAR_DENSITY       >(           1. );
  params.set<PHYS_MAX_T              >(          20. );  // max. simulation time in s, SI unit
  params.set<LENGTH                  >(           2. );  // length of the pie
  params.set<DIAMETER                >(           1. );  // diameter of the pipe
  params.set<REYNOLDS                >(          10. );
  params.set<LATTICE_RELAXATION_TIME >(          0.8 );
  params.set<RESIDUUM                >(         1e-5 );
  params.set<PARTIAL_SLIP_TUNER      >(           0. );  // for partialSlip only: 0->bounceBack, 1->freeSlip
  params.set<EOC                     >(        false );
  params.set<OVERLAP                 >(            3 );
  params.set<PHYS_CONVERGENCE_CHECK_T>([] (MyCase::ParametersD& p) {  // interval for the convergence check in s
    return p.get<PHYS_MAX_T>()*0.0125;
  });
  params.set<PHYS_STAT_ITER_T>([] (MyCase::ParametersD& p) {  // interval for the statistics output in s
    return p.get<PHYS_MAX_T>()/20.;
  });
  params.set<PHYS_VTK_ITER_T>([] (MyCase::ParametersD& p) {  // interval for the vtk output in s
    return p.get<PHYS_MAX_T>()/20.;
  });
  params.set<PHYS_DELTA_X>([] (MyCase::ParametersD& p) {
    return p.get<DIAMETER>()/p.get<RESOLUTION>();
  });
}
