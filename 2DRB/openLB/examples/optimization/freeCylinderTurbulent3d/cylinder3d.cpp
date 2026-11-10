/*  Lattice Boltzmann sample, written in C++, using the OpenLB library
 *
 *  Copyright (C) 2026 Shota Ito
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

#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::opti;

template <typename T, typename DESCRIPTOR>
using SmagorinskyHLBMdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Porous<momenta::BulkTuple>,
  equilibria::ThirdOrder,
  collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB>
>;

template <typename T, typename DESCRIPTOR>
using FrozenPrimalDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Porous<momenta::BulkTuple>,
  equilibria::ThirdOrder,
  collision::FrozenOmegaFromCell<collision::ThirdOrderRLB>
>;

template <typename T, typename DESCRIPTOR>
using FringeDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Porous<momenta::BulkTuple>,
  equilibria::ThirdOrder,
  collision::ParameterFromCell<collision::LES::SMAGORINSKY,
                               collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB>>
>;

enum class AdjointMode : int {
  FULL_ADJOINT_LES = 0,
  FROZEN_TURBULENCE = 1
};

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<float, descriptors::D3Q19<>>
>;

using ObjectiveF = functors::FarFieldDragF<0>;
using ObjectiveParameters = functors::FunctorParameters<ObjectiveF, MyCase::value_t, MyCase::descriptor_t>;
using PrimalDynamics = SmagorinskyHLBMdynamics<MyCase::value_t, MyCase::descriptor_t>;
using FrozenPrimalDynamicsT = FrozenPrimalDynamics<MyCase::value_t, MyCase::descriptor_t>;
using ControlledField = descriptors::POROSITY;

using MyOptiCase = OptiCaseAdjoint<
  Controlled, MyCase,
  Adjoint, MyCase,
  Objective, ObjectiveParameters
>;

namespace olb::parameters {

struct RADIUS_CYLINDER                    : public descriptors::FIELD_BASE<1> { };
struct DOMAIN_EXTENT_X                    : public descriptors::FIELD_BASE<1> { };
struct DOMAIN_EXTENT_Y                    : public descriptors::FIELD_BASE<1> { };
struct DOMAIN_EXTENT_Z                    : public descriptors::FIELD_BASE<1> { };
struct CENTER_CYLINDER_X                  : public descriptors::FIELD_BASE<1> { };
struct CENTER_CYLINDER_Y                  : public descriptors::FIELD_BASE<1> { };
struct AVERAGING_START_FRACTION           : public descriptors::FIELD_BASE<1> { };
struct RAMP_UP_END_FRACTION               : public descriptors::FIELD_BASE<1> { };
struct ADJOINT_RAMP_UP_UPDATE             : public descriptors::FIELD_BASE<1> { };
struct ADJOINT_RAMP_UP_END_FRACTION       : public descriptors::FIELD_BASE<1> { };
struct CONTROL_SHIFT                      : public descriptors::FIELD_BASE<1> { };
struct SAMPLE_POS                         : public descriptors::FIELD_BASE<0,1> { };
struct DISCRETE_ADJOINT_GRADIENTS         : public descriptors::FIELD_BASE<1> { };
struct ADF_GRADIENTS                      : public descriptors::FIELD_BASE<1> { };
struct FDQ_GRADIENTS                      : public descriptors::FIELD_BASE<1> { };
struct NUMBER_SAMPLES                     : public descriptors::FIELD_BASE<1> { };
struct OBJECTIVE_DOMAIN_X                 : public descriptors::FIELD_BASE<1> { };
struct OBJECTIVE_DOMAIN_Y                 : public descriptors::FIELD_BASE<1> { };
struct OBJECTIVE_DOMAIN_Z                 : public descriptors::FIELD_BASE<1> { };
struct FRINGE_THICKNESS                   : public descriptors::FIELD_BASE<1> { };
struct MAX_FRINGE_SMAGORINSKY_CONSTANT    : public descriptors::FIELD_BASE<1> { };
struct ADJOINT_PHYS_VTK_ITER_T            : public descriptors::FIELD_BASE<1> { };
struct ADJOINT_LATTICE_RELAXATION_TIME    : public descriptors::FIELD_BASE<1> { };
struct ADJOINT_SGS_TREATMENT              : public descriptors::TYPED_FIELD_BASE<AdjointMode,1> { };

}

template <typename CASE>
Mesh<typename CASE::value_t, CASE::d> createMesh(typename CASE::ParametersD& parameters) {
  using T = typename CASE::value_t;
  const T physLengthX = parameters.template get<parameters::DOMAIN_EXTENT_X>();
  const T physLengthY = parameters.template get<parameters::DOMAIN_EXTENT_Y>();
  const T physLengthZ = parameters.template get<parameters::DOMAIN_EXTENT_Z>();
  const T physDeltaX  = 2*parameters.template get<parameters::RADIUS_CYLINDER>()/parameters.template get<parameters::RESOLUTION>();
  const Vector extent{physLengthX, physLengthY, physLengthZ};
  const Vector origin{0, 0, 0};
  IndicatorCuboid3D<T> cuboid(extent, origin);
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.template get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({false, true, true});
  return mesh;
}

// Stores geometry information in form of material numbers
template <typename CASE>
void prepareGeometry(CASE& myCase)
{
  using T = CASE::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry  = myCase.getGeometry();
  const T physLengthX = parameters.template get<parameters::DOMAIN_EXTENT_X>();
  const T physDeltaX  = 2*parameters.template get<parameters::RADIUS_CYLINDER>()
                      /   parameters.template get<parameters::RESOLUTION>();
  const T physLengthY = parameters.template get<parameters::DOMAIN_EXTENT_Y>();
  const T physLengthZ = parameters.template get<parameters::DOMAIN_EXTENT_Z>();
  const T radiusCylinder  = parameters.template get<parameters::RADIUS_CYLINDER>();

  geometry.rename(0, 2);
  geometry.rename(2, 1, {1, 1, 1});
  Vector<T, 3> extent(physLengthX, physLengthY, physLengthZ);
  Vector<T, 3> origin;
  // Set material number for inflow
  extent[0] = 2.*physDeltaX;
  origin[0] = -physDeltaX;
  IndicatorCuboid3D<T> inflow(extent, origin);
  geometry.rename(2, 3, 1, inflow);
  // Set material number for outflow
  origin[0] = physLengthX-physDeltaX;
  IndicatorCuboid3D<T> outflow(extent, origin);
  geometry.rename(2, 4, 1, outflow);
  // Set material number for cylinder
  Vector<T,3> center1(parameters.template get<parameters::CENTER_CYLINDER_X>(), parameters.template get<parameters::CENTER_CYLINDER_Y>(), -physDeltaX);
  Vector<T,3> center2(parameters.template get<parameters::CENTER_CYLINDER_X>(), parameters.template get<parameters::CENTER_CYLINDER_Y>(), physLengthZ+physDeltaX);
  std::shared_ptr<IndicatorF3D<T>> cylinder = std::make_shared<IndicatorCylinder3D<T>>(center1, center2, radiusCylinder);
  geometry.rename(1, 5, cylinder);

  Vector<T,3> center(parameters.template get<parameters::CENTER_CYLINDER_X>() + 10. * radiusCylinder, parameters.template get<parameters::CENTER_CYLINDER_Y>(), physLengthZ / T(2.0));
  IndicatorCuboid3D<T> objectiveDomain(parameters.template get<parameters::OBJECTIVE_DOMAIN_X>(),
                                       parameters.template get<parameters::OBJECTIVE_DOMAIN_Y>(),
                                       parameters.template get<parameters::OBJECTIVE_DOMAIN_Z>(),
                                       center);
  geometry.rename(1, 6, objectiveDomain);

  // Add fringe zone
  origin[0] = physLengthX - parameters.template get<parameters::FRINGE_THICKNESS>() * physDeltaX;
  extent[0] = parameters.template get<parameters::FRINGE_THICKNESS>() * physDeltaX;
  origin[1] = -physDeltaX;
  extent[1] = physLengthY + 2.0*physDeltaX;
  origin[2] = -physDeltaX;
  extent[2] = physLengthZ + 2.0*physDeltaX;
  IndicatorCuboid3D<T> fringeI(extent, origin);
  geometry.rename(1, 7, fringeI);

  // Removes all not needed boundary voxels outside the surface
  geometry.print();
}

void setObjectiveParameters(auto& converter, auto& objectiveParameters) {
  objectiveParameters.template set<parameters::CONVERSION_FORCE>(converter.getConversionFactorForce());
  objectiveParameters.template set<parameters::LATTICE_REFERENCE_VELOCITY>(converter.getCharLatticeVelocity());
  objectiveParameters.template set<parameters::LATTICE_REFERENCE_PRESSURE>(0.);
}

// Set up the geometry of the simulation
template <typename CASE>
void prepareLattice(CASE& myCase, auto& objectiveParameters)
{
  using T = typename CASE::value_t;
  using DESCRIPTOR = typename CASE::descriptor_t;
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes{});
  auto& geometry   = myCase.getGeometry();
  lattice.template setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    parameters.template get<parameters::RESOLUTION>(),              // resolution: number of voxels per charPhysL
    parameters.template get<parameters::LATTICE_RELAXATION_TIME>(), // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    2.*parameters.template get<parameters::RADIUS_CYLINDER>(),      // charPhysLength: reference length of simulation geometry
    parameters.template get<parameters::PHYS_CHAR_VELOCITY>(),      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    parameters.template get<parameters::PHYS_CHAR_VELOCITY>() *     // physViscosity: physical kinematic viscosity in __m^2 / s__
    2.*parameters.template get<parameters::RADIUS_CYLINDER>() /
    parameters.template get<parameters::REYNOLDS>(),
    parameters.template get<parameters::PHYS_CHAR_DENSITY>()        // physDensity: physical density in __kg / m^3__
  );
  lattice.getUnitConverter().print();
  auto& converter = lattice.getUnitConverter();

  // Set Dynamics and BC
  auto bulkIndicator = geometry.getMaterialIndicator({1,5,6});
  dynamics::set<typename SmagorinskyHLBMdynamics<T,DESCRIPTOR>::template wrap_collision<collision::TrackAveragePopulation>>(lattice, bulkIndicator);
  boundary::set<boundary::LocalVelocity>(lattice, geometry.getMaterialIndicator({3}));
  boundary::set<T,DESCRIPTOR,boundary::InterpolatedPressure<T,DESCRIPTOR,FringeDynamics<T,DESCRIPTOR>>>(lattice,
                                                                                                        geometry.getMaterialIndicator({4}),
                                                                                                        geometry.getMaterialIndicator({7}),
                                                                                                        geometry.getMaterialIndicator({0}));
  dynamics::set<typename FringeDynamics<T,DESCRIPTOR>::template wrap_collision<collision::TrackAveragePopulation>>(lattice, geometry.getMaterialIndicator({7}));
  lattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  lattice.template setParameter<collision::LES::SMAGORINSKY>(parameters.template get<parameters::SMAGORINSKY_CONSTANT>());

  // Porosity field
  const T physDeltaX  = converter.getPhysDeltaX();
  fields::set<descriptors::POROSITY>(lattice, geometry.getMaterialIndicator({1,3,4,6,7}), 1.);
  fields::set<descriptors::POROSITY>(lattice, geometry.getMaterialIndicator({5}), 0);

  // Fringezone
  const T thickness = parameters.template get<parameters::FRINGE_THICKNESS>();
  const T bulkSmagoConst = parameters.template get<parameters::SMAGORINSKY_CONSTANT>();
  const T maxSmagoConst = parameters.template get<parameters::MAX_FRINGE_SMAGORINSKY_CONSTANT>();
  const T physLengthX = parameters.template get<parameters::DOMAIN_EXTENT_X>();
  const T physLengthY = parameters.template get<parameters::DOMAIN_EXTENT_Y>();
  const T physLengthZ = parameters.template get<parameters::DOMAIN_EXTENT_Z>();
  Vector<T,3> origin(physLengthX - thickness * physDeltaX, -physDeltaX, -physDeltaX);
  Vector<T,3> extent(thickness * physDeltaX, physLengthY + 2*physDeltaX, physLengthZ + 2*physDeltaX);
  IndicatorCuboid3D<T> fringe(extent, origin);
  SuperIndicatorFfromIndicatorF3D<T> fringeI(fringe, geometry);
  AnalyticalFfromCallableF<3,T,T> fringeF([&](Vector<T,3> physR) -> Vector<T,1> {
    return bulkSmagoConst + ((physR[0] - origin[0]) / (thickness * physDeltaX)) * (maxSmagoConst - bulkSmagoConst);
  });
  fields::set<collision::LES::SMAGORINSKY>(lattice, fringeI, fringeF);

  // Tracking type post-processors for averaged primal fields and their gradients
  auto objectiveDomain = geometry.getMaterialIndicator({6});
  lattice.template addPostProcessor<stage::Evaluation>(objectiveDomain, meta::id<operators::TrackAverageFunctorO<ObjectiveF, opti::AVERAGE_J>>());

  switch (parameters.template get<parameters::ADJOINT_SGS_TREATMENT>()) {
    case AdjointMode::FULL_ADJOINT_LES: {
      lattice.template addPostProcessor<stage::Evaluation>(objectiveDomain,
        meta::id<operators::TrackAverageFunctorO<functors::DerivativeF<ObjectiveF, descriptors::POPULATION, PrimalDynamics>, opti::AVERAGE_DJDF>>());
      lattice.template addPostProcessor<stage::Evaluation>(objectiveDomain,
        meta::id<operators::TrackAverageFunctorO<functors::DerivativeF<ObjectiveF, ControlledField, PrimalDynamics>, opti::AVERAGE_DJDALPHA<ControlledField>>>());
      break;
    }
    case AdjointMode::FROZEN_TURBULENCE: {
      lattice.template addPostProcessor<stage::Evaluation>(objectiveDomain,
        meta::id<operators::TrackAverageFunctorO<functors::DerivativeF<ObjectiveF, descriptors::POPULATION, FrozenPrimalDynamicsT>, opti::AVERAGE_DJDF>>());
      lattice.template addPostProcessor<stage::Evaluation>(objectiveDomain,
        meta::id<operators::TrackAverageFunctorO<functors::DerivativeF<ObjectiveF, ControlledField, FrozenPrimalDynamicsT>, opti::AVERAGE_DJDALPHA<ControlledField>>>());
      auto fluidDomain = geometry.getMaterialIndicator({1,5,6,7});
      lattice.template addPostProcessor<stage::Evaluation>(fluidDomain,
        meta::id<operators::TrackAverageFunctorO<functors::SmagorinskyTurbulentTauF<PrimalDynamics>, descriptors::AVERAGE_TURBULENT_TAU>>());
      lattice.template setParameter<parameters::OMEGA>(converter.getLatticeRelaxationFrequency());
      lattice.template setParameter<parameters::SMAGORINSKY_CONSTANT>(parameters.template get<parameters::SMAGORINSKY_CONSTANT>());
      break;
      }
  }

  // Populate the parameters required for the objective functional here once
  setObjectiveParameters(converter, objectiveParameters);
  functors::setAllParametersIn(lattice, objectiveParameters);
  lattice.template setParameter<parameters::FACTOR>(util::pow(converter.getPhysDeltaX(),3));

  // Operator for ramping-up inflow velocity on the GPU
  auto& smoothInflowO = myCase.setCouplingOperator("SmoothInflowO",
                                                   SingleLatticeO<operators::SmoothConstInflowO>{},
                                                   names::Lattice1{}, lattice);
  smoothInflowO.restrictTo(geometry.getMaterialIndicator({3}));
  smoothInflowO.template setParameter<operators::SmoothConstInflowO::INFLOW_VELOCITY>({0.0, 0.0, 0.0});
}

template <typename CASE>
void setInitialValues(CASE& myCase) {
  auto& lattice = myCase.getLattice(NavierStokes{});

  lattice.initialize();
  lattice.writeSummary();
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
template <typename CASE>
void setTemporalValues(CASE& myCase,
                      std::size_t iT)
{
  using T = typename CASE::value_t;
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes{});
  auto& geometry   = myCase.getGeometry();
  auto& converter  = lattice.getUnitConverter();

  const T maxPhysT   = parameters.template get<parameters::PHYS_MAX_T>();

  if (myCase.hasName(Controlled{})) {
    std::size_t iTmaxStart = converter.getLatticeTime(
      maxPhysT*parameters.template get<parameters::RAMP_UP_END_FRACTION>());
    std::size_t iTupdate   = converter.getLatticeTime(
      parameters.template get<parameters::RAMP_UP_T>());
    if (iT % iTupdate == 0 && iT <= iTmaxStart) {
      using Coupling_t = SuperLatticeCoupling<SingleLatticeO<operators::SmoothConstInflowO>,
                                              meta::map<names::Lattice1,
                                                        descriptors::VALUED_DESCRIPTOR<T,typename CASE::descriptor_t>>>;
      Coupling_t& smoothInflowO = dynamic_cast<Coupling_t&>(myCase.getOperator("SmoothInflowO"));
      const T frac = 1. - util::cos(0.5 * std::numbers::pi_v<T> * static_cast<T>(iT) / static_cast<T>(iTmaxStart));
      const T targetU = frac * converter.getCharLatticeVelocity();
      smoothInflowO.template setParameter<operators::SmoothConstInflowO::INFLOW_VELOCITY>({targetU, 0.0, 0.0});
      smoothInflowO.apply();
    }
  } else {
    std::size_t iTmaxStart = lattice.getUnitConverter().getLatticeTime(
      maxPhysT*parameters.template get<parameters::ADJOINT_RAMP_UP_END_FRACTION>());
    std::size_t iTupdate   = lattice.getUnitConverter().getLatticeTime(
      parameters.template get<parameters::ADJOINT_RAMP_UP_UPDATE>());
    if (iT % iTupdate == 0 && iT <= iTmaxStart) {
      T frac{};
      PolynomialStartScale<T,std::size_t>(iTmaxStart, T(1))(&frac, &iT);
      auto scaleF = makeWriteFunctorO<functors::ScaleFieldF<opti::DJDF_TMP>, opti::DJDF>(lattice);
      scaleF->template setParameter<parameters::FACTOR>(frac);
      scaleF->restrictTo(geometry.getMaterialIndicator({5,6}));
      scaleF->apply();
    }
  }
}

// Computes the pressure drop between the voxels before and after the cylinder
template <typename CASE>
void getResults(CASE& myCase,
                util::Timer<typename CASE::value_t>& timer,
                std::size_t iT,
                auto& objectiveParameters)
{
  OstreamManager clout( std::cout,"getResults" );
  using T          = typename CASE::value_t;
  using DESCRIPTOR = typename CASE::descriptor_t;
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes{});
  auto& geometry   = myCase.getGeometry();

  const UnitConverter<T,DESCRIPTOR>& converter = lattice.getUnitConverter();
  const std::size_t statIter = converter.getLatticeTime(parameters.template get<parameters::PHYS_STAT_ITER_T>());
  const std::size_t iTstartAvg = converter.getLatticeTime(parameters.template get<parameters::PHYS_MAX_T>() *
                                                          parameters.template get<parameters::AVERAGING_START_FRACTION>());
  if (iT == iTstartAvg) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    SuperLatticeField3D<T,DESCRIPTOR,descriptors::POPULATION> pop(lattice);
    fields::set<descriptors::AVERAGE_POPULATION>(lattice, geometry.getMaterialIndicator({1,2,3,4,5,6,7}), pop);
    lattice.setProcessingContext(ProcessingContext::Simulation);
    auto objectiveDomain = geometry.getMaterialIndicator({6});
    auto objectiveF = makeWriteFunctorO<ObjectiveF, opti::AVERAGE_J>(lattice);
    functors::setAllParametersIn(*objectiveF, objectiveParameters);
    objectiveF->restrictTo(objectiveDomain);
    objectiveF->apply();

    switch (parameters.template get<parameters::ADJOINT_SGS_TREATMENT>()) {
      case AdjointMode::FULL_ADJOINT_LES: {
        break;
      }
      case AdjointMode::FROZEN_TURBULENCE: {
        auto turbTauF = makeWriteFunctorO<functors::SmagorinskyTurbulentTauF<PrimalDynamics>, descriptors::AVERAGE_TURBULENT_TAU>(lattice);
        turbTauF->template setParameter<parameters::OMEGA>(converter.getLatticeRelaxationFrequency());
        turbTauF->template setParameter<parameters::SMAGORINSKY_CONSTANT>(parameters.template get<parameters::SMAGORINSKY_CONSTANT>());
        turbTauF->restrictTo(geometry.getMaterialIndicator({1,5,6,7}));
        turbTauF->apply();
        break;
      }
    }

    if constexpr(std::is_same_v<CASE,MyCase>) {
      switch (parameters.template get<parameters::ADJOINT_SGS_TREATMENT>()) {
        case AdjointMode::FULL_ADJOINT_LES: {
          auto objectiveDF = makeWriteFunctorO<functors::DerivativeF<ObjectiveF, descriptors::POPULATION, PrimalDynamics>,
                                               opti::AVERAGE_DJDF>(lattice);
          functors::setAllParametersIn(*objectiveDF, objectiveParameters);
          objectiveDF->template setParameter<parameters::FACTOR>(util::pow(converter.getPhysDeltaX(),3));
          objectiveDF->restrictTo(objectiveDomain);
          objectiveDF->apply();

          auto objectiveDA = makeWriteFunctorO<functors::DerivativeF<ObjectiveF, ControlledField, PrimalDynamics>,
                                               opti::AVERAGE_DJDALPHA<ControlledField>>(lattice);
          functors::setAllParametersIn(*objectiveDA, objectiveParameters);
          objectiveDA->template setParameter<parameters::FACTOR>(util::pow(converter.getPhysDeltaX(),3));
          objectiveDA->restrictTo(objectiveDomain);
          objectiveDA->apply();
          break;
        }
        case AdjointMode::FROZEN_TURBULENCE: {
          auto objectiveDF = makeWriteFunctorO<functors::DerivativeF<ObjectiveF, descriptors::POPULATION, FrozenPrimalDynamicsT>,
                                               opti::AVERAGE_DJDF>(lattice);
          functors::setAllParametersIn(*objectiveDF, objectiveParameters);
          objectiveDF->template setParameter<parameters::FACTOR>(util::pow(converter.getPhysDeltaX(),3));
          objectiveDF->restrictTo(objectiveDomain);
          objectiveDF->apply();

          auto objectiveDA = makeWriteFunctorO<functors::DerivativeF<ObjectiveF, ControlledField, FrozenPrimalDynamicsT>,
                                               opti::AVERAGE_DJDALPHA<ControlledField>>(lattice);
          functors::setAllParametersIn(*objectiveDA, objectiveParameters);
          objectiveDA->template setParameter<parameters::FACTOR>(util::pow(converter.getPhysDeltaX(),3));
          objectiveDA->restrictTo(objectiveDomain);
          objectiveDA->apply();
          break;
        }
      }
    }
  } else if (iT < iTstartAvg) {
    lattice.template setParameter<descriptors::LATTICE_TIME_T>(T(2));
  } else {
    lattice.template setParameter<descriptors::LATTICE_TIME_T>(T(iT - iTstartAvg + 1));
  }

  if (iT%statIter == 0) {
    // Timer console output
    timer.update(iT);
    timer.printStep();

    // Lattice statistics console output
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
  }

  const std::size_t vtkIter  = converter.getLatticeTime(parameters.template get<parameters::PHYS_VTK_ITER_T>());
  const std::size_t adjointVtkIter  = converter.getLatticeTime(parameters.template get<parameters::ADJOINT_PHYS_VTK_ITER_T>());
  if (myCase.hasName(Controlled{})) {
    if (iT == 0) {
      SuperVTMwriter3D<T> vtmWriter("cylinder3d");
      vtmWriter.createMasterFile();
    }

    // Writes the vtk files
    if ((iT % vtkIter == 0 && iT >= 0) || iT == 0) {
      lattice.setProcessingContext(ProcessingContext::Evaluation);
      lattice.scheduleBackgroundOutputVTK([&, iT](auto task) {
        SuperVTMwriter3D<T> vtmWriter("cylinder3d");
        SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(lattice, converter);
        SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(lattice, converter);
        SuperLatticeField3D<T,DESCRIPTOR,descriptors::AVERAGE_POPULATION> popAvg(lattice);
        SuperLatticeField3D<T,DESCRIPTOR,descriptors::POROSITY> porosity(lattice);
        SuperLatticeField3D<T,DESCRIPTOR,opti::AVERAGE_J> jFAvg(lattice);
        SuperLatticeField3D<T,DESCRIPTOR,opti::AVERAGE_DJDF> djdfFAvg(lattice);
        SuperLatticeField3D<T,DESCRIPTOR,opti::AVERAGE_DJDALPHA<ControlledField>> djdaFAvg(lattice);
        djdaFAvg.getName() = "djdaFAvg";

        vtmWriter.addFunctor(porosity);
        vtmWriter.addFunctor(velocity);
        vtmWriter.addFunctor(pressure);
        vtmWriter.addFunctor(popAvg);
        vtmWriter.addFunctor(jFAvg);
        vtmWriter.addFunctor(djdfFAvg);
        vtmWriter.addFunctor(djdaFAvg);
        task(vtmWriter, iT);
      });
    }
  } else {
    if (iT == 0) {
      SuperVTMwriter3D<T> vtmWriter("adjointCylinder3d");
      vtmWriter.createMasterFile();
    }

    // Writes the vtk files
    if ((iT % adjointVtkIter == 0 && iT >= 0) || iT == 0) {
      lattice.setProcessingContext(ProcessingContext::Evaluation);
      lattice.scheduleBackgroundOutputVTK([&, iT](auto task) {
        SuperVTMwriter3D<T> vtmWriter("adjointCylinder3d");
        SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(lattice, converter);
        SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(lattice, converter);
        SuperLatticeField3D<T,DESCRIPTOR,descriptors::POPULATION> adjointsF(lattice);
        SuperLatticeField3D<T,DESCRIPTOR,opti::DJDF> djdfF(lattice);
        SuperLatticeField3D<T,DESCRIPTOR,opti::F> primalF(lattice);

        vtmWriter.addFunctor(velocity);
        vtmWriter.addFunctor(pressure);
        vtmWriter.addFunctor(adjointsF);
        vtmWriter.addFunctor(djdfF);
        vtmWriter.addFunctor(primalF);
        task(vtmWriter, iT);
      });
    }
  }
}

template <typename CASE>
void simulate(CASE& myCase, auto& objectiveParameters){
  OstreamManager clout(std::cout, "simulate");
  using T = typename CASE::value_t;
  auto& parameters = myCase.getParameters();
  auto& converter = myCase.getLattice(NavierStokes{}).getUnitConverter();
  const T physMaxT = parameters.template get<parameters::PHYS_MAX_T>();
  const std::size_t iTmax = converter.getLatticeTime(physMaxT);
  const std::size_t iTstartAvg = converter.getLatticeTime(parameters.template get<parameters::PHYS_MAX_T>() *
                                                          parameters.template get<parameters::AVERAGING_START_FRACTION>());

  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    myCase.getLattice(NavierStokes{}).collideAndStream();

    if(myCase.hasName(Controlled{}) && iT >= iTstartAvg) {
      myCase.getLattice(NavierStokes{}).executePostProcessors(stage::Evaluation());
    }

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, timer, iT, objectiveParameters);
  }

  timer.stop();
  timer.printSummary();
}

void setInitialControl(MyOptiCase& optiCase) {
  OstreamManager clout(std::cout, "setInitialControl");
  using T = MyCase::value_t;
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& lattice = controlledCase.getLattice(NavierStokes{});
  auto& geometry = controlledCase.getGeometry();
  auto& control = optiCase.getController();
  auto& converter = lattice.getUnitConverter();
  auto& parameters = controlledCase.getParameters();

  // Intialize controlled field in superLattice
  T porosity = projection::permeabilityToPorosity(parameters.get<parameters::INITIAL_CONTROL_SCALAR>(), converter);
  clout << "GridTerm: " << projection::gridTerm(converter) << std::endl;
  clout << "Initial porosity: " << porosity << std::endl;
  fields::set<ControlledField>(lattice, geometry.getMaterialIndicator({5}), porosity);
  control.template setProjection<projection::PermeabilityToPorosity<T>>(converter);
  control.template set<ControlledField>(geometry, 5, lattice);
}

void applyControl(MyOptiCase& optiCase) {
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& lattice = controlledCase.getLattice(NavierStokes{});
  optiCase.getController().template setUpdatedControlsOnField<ControlledField>(lattice);
  lattice.template setProcessingContext<Array<ControlledField>>(ProcessingContext::Simulation);
}

MyCase::value_t computeObjective(MyOptiCase& optiCase) {
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& converter = controlledLattice.getUnitConverter();
  auto objectiveDomain = optiCase.getCase(Controlled{}).getGeometry().getMaterialIndicator({6});
  return integrateField<opti::AVERAGE_J>(controlledLattice, objectiveDomain, converter.getPhysDeltaX())[0];
}

MyCase::value_t solvePrimalAndComputeObjective(MyOptiCase& optiCase) {
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& objectiveParameters = optiCase.get(Objective{});

  // Reset prior simulation lattice
  controlledCase.resetLattices();

  // Prepare new lattices
  prepareLattice(controlledCase, objectiveParameters);
  setInitialValues(controlledCase);

  // Set updated controls in the simulation
  applyControl(optiCase);

  // Execute simulation
  simulate(controlledCase, objectiveParameters);

  // Compute Objective value from simulation results
  return computeObjective(optiCase);
}

void prepareAdjointLattice(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& adjointLattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& geometry = optiCase.getCase(Adjoint{}).getGeometry();
  auto& parameters = optiCase.getCase(Adjoint{}).getParameters();

  adjointLattice.setUnitConverter(controlledLattice.getUnitConverter());

  // Define dual dynamics
  switch (parameters.template get<parameters::ADJOINT_SGS_TREATMENT>()) {
    case AdjointMode::FULL_ADJOINT_LES: {
      dynamics::set<olb::Dual<SmagorinskyHLBMdynamics<T,DESCRIPTOR>>,T,DESCRIPTOR>(adjointLattice, geometry.getMaterialIndicator({1,5,6,7}));
      adjointLattice.template setParameter<descriptors::OMEGA>(1.0 / parameters.template get<parameters::ADJOINT_LATTICE_RELAXATION_TIME>());
      adjointLattice.template setParameter<collision::LES::SMAGORINSKY>(parameters.template get<parameters::SMAGORINSKY_CONSTANT>());
      break;
    }
    case AdjointMode::FROZEN_TURBULENCE: {
      dynamics::set<olb::Dual<FrozenPrimalDynamics<T,DESCRIPTOR>>,T,DESCRIPTOR>(adjointLattice, geometry.getMaterialIndicator({1,5,6,7}));
      // This should be only a dummy parameter
      adjointLattice.template setParameter<descriptors::OMEGA>(adjointLattice.getUnitConverter().getLatticeRelaxationFrequency());
      break;
    }
  }
  boundary::set<boundary::ZeroPopulation>(adjointLattice, geometry.getMaterialIndicator({3,4}));
  adjointLattice.revertStreaming();
}

void setAdjointInitialValues(MyOptiCase& optiCase) {
  auto& adjointCase = optiCase.getCase(Adjoint{});
  auto& adjointLattice = adjointCase.getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& parameters = adjointCase.getParameters();
  auto& geometry = adjointCase.getGeometry();

  setInitialValues(adjointCase);

  copyFields<opti::AVERAGE_DJDF, opti::DJDF>(controlledLattice, adjointLattice);
  copyFields<ControlledField, ControlledField>(controlledLattice, adjointLattice);
  // Copy adjoint source term into tmp field to apply ramp-up
  copyFields<opti::DJDF, opti::DJDF_TMP>(adjointLattice, adjointLattice);
  copyFields<descriptors::AVERAGE_POPULATION, opti::F>(controlledLattice, adjointLattice);
  switch (parameters.template get<parameters::ADJOINT_SGS_TREATMENT>()) {
    case AdjointMode::FULL_ADJOINT_LES: {
      break;
    }
    case AdjointMode::FROZEN_TURBULENCE: {
      copyFields<descriptors::AVERAGE_TURBULENT_TAU, descriptors::AVERAGE_TURBULENT_TAU>(controlledLattice, adjointLattice);
      fields::set<descriptors::OMEGA>(adjointLattice, geometry.getMaterialIndicator({1,5,6,7}), parameters.template get<parameters::ADJOINT_LATTICE_RELAXATION_TIME>());
      adjointLattice.template setProcessingContext<Array<descriptors::OMEGA>>(ProcessingContext::Simulation);
      break;
    }
  }
}

std::vector<MyCase::value_t> computeDerivative(MyOptiCase& optiCase) {
  auto& adjointLattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& converter = controlledLattice.getUnitConverter();
  auto objectiveDomain = optiCase.getCase(Controlled{}).getGeometry().getMaterialIndicator({6});
  auto& parameters = optiCase.getCase(Adjoint{}).getParameters();
  auto& geometry = optiCase.getCase(Adjoint{}).getGeometry();

  auto optimalityO = makeWriteFunctorO<functors::OptimalityF<PrimalDynamics, ControlledField>, opti::SENSITIVITY<ControlledField>>(adjointLattice);
  optiCase.getController().template setUpdatedProjectionDerivativesOnField<opti::DPROJECTIONDALPHA<ControlledField>>(adjointLattice);
  adjointLattice.template setProcessingContext<Array<opti::DPROJECTIONDALPHA<ControlledField>>>(ProcessingContext::Simulation);

  switch (parameters.template get<parameters::ADJOINT_SGS_TREATMENT>()) {
    case AdjointMode::FULL_ADJOINT_LES: {
      auto dCdalphaO = makeWriteFunctorO<functors::DerivativeF<functors::CollisionUsingAveragePopulationF<PrimalDynamics>, ControlledField, PrimalDynamics>, opti::DCDALPHA<ControlledField>>(controlledLattice);
      dCdalphaO->template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
      dCdalphaO->template setParameter<collision::LES::SMAGORINSKY>(parameters.template get<parameters::SMAGORINSKY_CONSTANT>());
      dCdalphaO->template setParameter<parameters::FACTOR>(1.0);
      dCdalphaO->apply();
      break;
    }
    case AdjointMode::FROZEN_TURBULENCE: {
      using T = MyCase::value_t;
      using DESCRIPTOR = MyCase::descriptor_t;
      fields::set<descriptors::OMEGA>(controlledLattice, geometry.getMaterialIndicator({1,5,6,7}), parameters.template get<parameters::LATTICE_RELAXATION_TIME>());
      controlledLattice.template setProcessingContext<Array<descriptors::OMEGA>>(ProcessingContext::Simulation);
      auto dCdalphaO = makeWriteFunctorO<functors::DerivativeF<functors::CollisionUsingAveragePopulationF<FrozenPrimalDynamics<T,DESCRIPTOR>>, ControlledField, FrozenPrimalDynamicsT>, opti::DCDALPHA<ControlledField>>(controlledLattice);
      dCdalphaO->template setParameter<parameters::FACTOR>(1.0);
      dCdalphaO->apply();
      break;
    }
  }

  copyFields<opti::DCDALPHA<ControlledField>,opti::DCDALPHA<ControlledField>>(controlledLattice, adjointLattice);
  copyFields<opti::AVERAGE_DJDALPHA<ControlledField>,opti::DJDALPHA<ControlledField>>(controlledLattice, adjointLattice);
  optimalityO->apply();

  // Return serial vector containing total derivatives of objective regarding controls
  return getSerializedFromField<opti::SENSITIVITY<ControlledField>>(adjointLattice, objectiveDomain);
}

auto getShiftedControl(MyOptiCase& optiCase, MyCase::value_t shiftedLocalPorosity) {
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& lattice = controlledCase.getLattice(NavierStokes{});
  auto& geometry = controlledCase.getGeometry();
  auto& control = optiCase.getController();
  auto& parameters = controlledCase.getParameters();

  setInitialControl(optiCase);
  Vector<MyCase::value_t,3> samplePhysPos(parameters.get<parameters::SAMPLE_POS>());

  if ( auto samplePos = controlledCase.getMesh().getCuboidDecomposition().getLatticeR(samplePhysPos)) {
    if (controlledCase.getMesh().getLoadBalancer().isLocal((*samplePos)[0])) {
      auto cell = lattice.get(*samplePos);
      cell.template setField<ControlledField>(shiftedLocalPorosity);
    }
    control.template set<ControlledField>(geometry, 5, lattice);
  } else {
    throw std::runtime_error("Invalid sample position.");
  }
  return control.get();
}

void sampleDerivativesWithAdjoints(MyCase& adjointCase) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& sLattice = adjointCase.getLattice(NavierStokes{});
  auto& parameters = adjointCase.getParameters();

  sLattice.template setProcessingContext<Array<opti::SENSITIVITY<ControlledField>>>(ProcessingContext::Evaluation);
  SuperLatticeField3D<T,DESCRIPTOR,opti::SENSITIVITY<ControlledField>> gradF(sLattice);
  AnalyticalFfromSuperF3D<T,T> grad(gradF);

  T output[1], input[3];
  input[0] = parameters.get<parameters::SAMPLE_POS>()[0];
  input[1] = parameters.get<parameters::SAMPLE_POS>()[1];
  input[2] = parameters.get<parameters::SAMPLE_POS>()[2];
  grad(output, input);
  std::cout << "Adjoint Gradient at " << parameters.get<parameters::SAMPLE_POS>()
                                      << ": (" << output[0]
                                      << ")" << std::endl;
  parameters.set<parameters::DISCRETE_ADJOINT_GRADIENTS>(output[0]);
}

std::vector<MyCase::value_t> solveAdjointAndComputeDerivative(MyOptiCase& optiCase) {
  auto& adjointCase = optiCase.getCase(Adjoint{});
  auto& objectiveParameters = optiCase.get(Objective{});
  // Reset prior simulation lattice
  adjointCase.resetLattices();

  // Prepare new lattices
  prepareAdjointLattice(optiCase);
  setAdjointInitialValues(optiCase);

  // Execute simulation
  simulate(adjointCase, objectiveParameters);

  // Compute Derivative value from simulation results
  return computeDerivative(optiCase);
}

void plotDerivatives(MyCase::value_t adjoint, unsigned index) {
  using T = MyCase::value_t;

  static Gnuplot<T> gplot("derivatives");
  gplot.setData(T(index), {adjoint}, {"adjoints"}, "bottom right", {'l', 'l'});
  gplot.writePDF();
}

auto sampleCircle(MyCase::value_t cx, MyCase::value_t cy, MyCase::value_t r, MyCase::value_t dx, MyCase::value_t z, int nSamples) {
  using T = MyCase::value_t;
  std::vector<std::vector<T>> points;
  points.reserve(nSamples);

  const T twoPi = 2.0 * M_PI;
  const T  rShifted = r - dx;
  for (int i = 0; i < nSamples; ++i) {
    const T theta = 0.5*M_PI - twoPi * i / nSamples;
    std::vector<T> point({0.0, 0.0, 0.0});
    point[0] = cx + rShifted*util::cos(theta);
    point[1] = cy + rShifted*util::sin(theta);
    point[2] = z;
    points.push_back(point);
  }
  return points;
}

void printSampledGradients(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  auto& myCaseParameters = optiCase.getCase(Controlled{}).getParameters();
  auto& adjointCase = optiCase.getCase(Adjoint{});

  // Evaluating gradients on the cylinder surface
  const T dx  = myCaseParameters.get<parameters::RADIUS_CYLINDER>()/myCaseParameters.get<parameters::RESOLUTION>();
  const T radius = myCaseParameters.get<parameters::RADIUS_CYLINDER>();
  const T sampleDepth = myCaseParameters.get<parameters::DOMAIN_EXTENT_Z>() / 2.0;
  auto samplePoints = sampleCircle(myCaseParameters.get<parameters::CENTER_CYLINDER_X>(),
                                   myCaseParameters.get<parameters::CENTER_CYLINDER_Y>(),
                                   radius, dx, sampleDepth, myCaseParameters.get<parameters::NUMBER_SAMPLES>());

  unsigned index = 0;
  for (auto point : samplePoints) {
    myCaseParameters.set<parameters::SAMPLE_POS>(point);
    sampleDerivativesWithAdjoints(adjointCase);

    std::cout << "Summary:" << std::endl;
    std::cout << "Averaged Dis.   Gradient at " << optiCase.getCase(Controlled{}).getParameters().get<parameters::SAMPLE_POS>()
                                                << ": " << optiCase.getCase(Controlled{}).getParameters().get<parameters::DISCRETE_ADJOINT_GRADIENTS>()
                                                << std::endl;
    plotDerivatives(optiCase.getCase(Controlled{}).getParameters().get<parameters::DISCRETE_ADJOINT_GRADIENTS>(), index);
    ++index;
  }
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION                     >(21);
    myCaseParameters.set<REYNOLDS                       >(3900.);
    myCaseParameters.set<PHYS_MAX_T                     >(800.);
    myCaseParameters.set<DOMAIN_EXTENT_X                >(4.0);
    myCaseParameters.set<DOMAIN_EXTENT_Y                >(1.01);
    myCaseParameters.set<DOMAIN_EXTENT_Z                >(.51);
    myCaseParameters.set<PHYS_CHAR_VELOCITY             >(0.2);
    myCaseParameters.set<LATTICE_RELAXATION_TIME        >(0.501);
    myCaseParameters.set<ADJOINT_LATTICE_RELAXATION_TIME>(0.65);
    myCaseParameters.set<PHYS_CHAR_DENSITY              >(1.0);
    myCaseParameters.set<RADIUS_CYLINDER                >(0.05);
    myCaseParameters.set<CENTER_CYLINDER_X              >(0.7);
    myCaseParameters.set<CENTER_CYLINDER_Y              >(0.5);
    myCaseParameters.set<PHYS_VTK_ITER_T                >(200.);
    myCaseParameters.set<ADJOINT_PHYS_VTK_ITER_T        >(100.);
    myCaseParameters.set<PHYS_STAT_ITER_T               >(1.5);
    myCaseParameters.set<AVERAGING_START_FRACTION       >(0.1);
    myCaseParameters.set<RAMP_UP_T                      >(0.01);
    myCaseParameters.set<RAMP_UP_END_FRACTION           >(0.075);
    myCaseParameters.set<ADJOINT_RAMP_UP_UPDATE         >(0.01);
    myCaseParameters.set<ADJOINT_RAMP_UP_END_FRACTION   >(0.8);
    myCaseParameters.set<SAMPLE_POS                     >({0.7, 0.2, 0.205});
    myCaseParameters.set<INITIAL_CONTROL_SCALAR         >(3.85488e-8);
    myCaseParameters.set<CONTROL_SHIFT                  >(1.e-10);
    myCaseParameters.set<NUMBER_SAMPLES                 >(50);
    myCaseParameters.set<OBJECTIVE_DOMAIN_X             >(0.004);
    myCaseParameters.set<OBJECTIVE_DOMAIN_Y             >(0.62);
    myCaseParameters.set<OBJECTIVE_DOMAIN_Z             >(0.62);
    myCaseParameters.set<SMAGORINSKY_CONSTANT           >(0.1);
    myCaseParameters.set<FRINGE_THICKNESS               >(200);
    myCaseParameters.set<MAX_FRINGE_SMAGORINSKY_CONSTANT>(1.);
    myCaseParameters.set<ADJOINT_SGS_TREATMENT          >(AdjointMode::FULL_ADJOINT_LES);
  }

  MyOptiCase optiCase;
  ObjectiveParameters objectiveF;
  optiCase.set<Objective>(objectiveF);

  Mesh mesh = createMesh<MyCase>(myCaseParameters);
  MyCase myCase(myCaseParameters, mesh);
  prepareGeometry(myCase);
  prepareLattice(myCase, objectiveF);

  MyCase adjointCase(myCaseParameters, mesh);
  prepareGeometry(adjointCase);

  optiCase.setCase<Controlled>(myCase);
  optiCase.setCase<Adjoint>(adjointCase);

  setInitialControl(optiCase);
  optiCase.setObjective(solvePrimalAndComputeObjective);
  optiCase.setDerivative(solveAdjointAndComputeDerivative);

  auto& controls = optiCase.getController().get();
  // Evaluate objective
  optiCase.computeObjective(controls);
  // Compute gradients with adjoints
  optiCase.computeDerivatives(controls);

  printSampledGradients(optiCase);

  // Prevent segfaults due to asyncrone vtk outputs
  optiCase.getCase(Controlled{}).getLattice(NavierStokes{}).setProcessingContext(ProcessingContext::Evaluation);
  optiCase.getCase(Controlled{}).getLattice(NavierStokes{}).setProcessingContext(ProcessingContext::Simulation);
  optiCase.getCase(Adjoint{}).getLattice(NavierStokes{}).setProcessingContext(ProcessingContext::Evaluation);
  optiCase.getCase(Adjoint{}).getLattice(NavierStokes{}).setProcessingContext(ProcessingContext::Simulation);
}
