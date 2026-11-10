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

// This app is based on the work in https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7211144.

// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::opti;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D2Q9<descriptors::POROSITY>>
>;

using MyADfCase = Case<
  NavierStokes, Lattice<util::ADf<MyCase::value_t,1>, MyCase::descriptor_t>
>;

using DissipationF = functors::AddF<functors::DissipationF, functors::PorousDissipationF>;
using FarFieldDragF = functors::FarFieldDragF<0>;

using PrimalDynamics = PorousBGKdynamics<MyCase::value_t, MyCase::descriptor_t>;
using ControlledField = descriptors::POROSITY;

using MyOptiCase = OptiCaseAdjoint<
  Controlled, MyCase,
  Adjoint, MyCase,
  Derivatives, MyADfCase
>;

enum class ObjectiveType : int {
  DISSIPATION = 0,
  FAR_FIELD_DRAG = 1
};

namespace olb::parameters {
struct RADIUS_CYLINDER               : public descriptors::FIELD_BASE<1> { };
struct DOMAIN_EXTENT_X               : public descriptors::FIELD_BASE<1> { };
struct DOMAIN_EXTENT_Y               : public descriptors::FIELD_BASE<1> { };
struct CENTER_CYLINDER_X             : public descriptors::FIELD_BASE<1> { };
struct CENTER_CYLINDER_Y             : public descriptors::FIELD_BASE<1> { };
struct AVERAGING_START_FRACTION      : public descriptors::FIELD_BASE<1> { };
struct RAMP_UP_END_FRACTION          : public descriptors::FIELD_BASE<1> { };
struct ADJOINT_RAMP_UP_UPDATE        : public descriptors::FIELD_BASE<1> { };
struct ADJOINT_RAMP_UP_END_FRACTION  : public descriptors::FIELD_BASE<1> { };
struct CONTROL_SHIFT                 : public descriptors::FIELD_BASE<1> { };
struct SAMPLE_POS                    : public descriptors::FIELD_BASE<0,1> { };
struct DISCRETE_ADJOINT_GRADIENTS    : public descriptors::FIELD_BASE<1> { };
struct ADF_GRADIENTS                 : public descriptors::FIELD_BASE<1> { };
struct FDQ_GRADIENTS                 : public descriptors::FIELD_BASE<1> { };
struct NUMBER_SAMPLES                : public descriptors::FIELD_BASE<1> { };
struct OBJECTIVE_DOMAIN_X            : public descriptors::FIELD_BASE<1> { };
struct OBJECTIVE_DOMAIN_Y            : public descriptors::FIELD_BASE<1> { };
struct VELOCITY_SAMPLING_POS_X       : public descriptors::FIELD_BASE<1> { };
struct VELOCITY_SAMPLING_POS_Y       : public descriptors::FIELD_BASE<1> { };
struct OBJECTIVE_TYPE                : public descriptors::TYPED_FIELD_BASE<ObjectiveType,1> { };
}

template <typename CASE>
Mesh<typename CASE::value_t, CASE::d> createMesh(typename CASE::ParametersD& parameters) {
  using T = typename CASE::value_t;
  const T physLengthX = parameters.template get<parameters::DOMAIN_EXTENT_X>();
  const T physLengthY = parameters.template get<parameters::DOMAIN_EXTENT_Y>();
  const T physDeltaX  = 2*parameters.template get<parameters::RADIUS_CYLINDER>()/parameters.template get<parameters::RESOLUTION>();
  const Vector extent{physLengthX, physLengthY};
  const Vector origin{0, 0};
  IndicatorCuboid2D<T> cuboid(extent, origin);
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.template get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({false, true});
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
  const T radiusCylinder  = parameters.template get<parameters::RADIUS_CYLINDER>();

  geometry.rename(0, 2);
  geometry.rename(2, 1, {1, 1});
  Vector<T, 2> extent(physLengthX, physLengthY);
  Vector<T, 2> origin;
  // Set material number for inflow
  extent[0] = 2.*physDeltaX;
  origin[0] = -physDeltaX;
  IndicatorCuboid2D<T> inflow(extent, origin );
  geometry.rename(2, 3, 1, inflow);
  // Set material number for outflow
  origin[0] = physLengthX-physDeltaX;
  IndicatorCuboid2D<T> outflow(extent, origin);
  geometry.rename(2, 4, 1, outflow);
  // Set material number for cylinder
  Vector<T,2> center(parameters.template get<parameters::CENTER_CYLINDER_X>(), parameters.template get<parameters::CENTER_CYLINDER_Y>());
  std::shared_ptr<IndicatorF2D<T>> circle = std::make_shared<IndicatorCircle2D<T>>(center, radiusCylinder);
  geometry.rename(1, 5, circle);

  switch (parameters.template get<parameters::OBJECTIVE_TYPE>()) {
    case ObjectiveType::DISSIPATION: {
      // Dissipation objective
      IndicatorCuboid2D<T> objectiveDomain(parameters.template get<parameters::OBJECTIVE_DOMAIN_X>(),
                                           parameters.template get<parameters::OBJECTIVE_DOMAIN_Y>(),
                                           center);
      geometry.rename(1, 6, objectiveDomain);
      break;
    }
    case ObjectiveType::FAR_FIELD_DRAG: {
      // FarFieldDrag objective
      center[0]+= radiusCylinder * 10;
      IndicatorCuboid2D<T> objectiveDomain(physDeltaX * 2.,
                                           parameters.template get<parameters::OBJECTIVE_DOMAIN_Y>(),
                                           center);
      geometry.rename(1, 6, objectiveDomain);
      break;
    }
  }

  // Remove the walls again
  geometry.rename(2,1);
  // Removes all not needed boundary voxels outside the surface
  geometry.print();
}

// Set up the geometry of the simulation
template <typename CASE>
void prepareLattice(CASE& myCase)
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

  // Material=1 -->bulk dynamics
  auto bulkIndicator = geometry.getMaterialIndicator({1,5,6});
  dynamics::set<typename PorousBGKdynamics<T,DESCRIPTOR>::template wrap_collision<collision::TrackAveragePopulation>>(lattice, bulkIndicator);
  boundary::set<boundary::LocalVelocity>(lattice, geometry.getMaterialIndicator({3}));
  boundary::set<boundary::LocalPressure>(lattice, geometry.getMaterialIndicator({4}));

  fields::set<descriptors::POROSITY>(lattice, geometry.getMaterialIndicator({1,3,4,6}), 1.);
  fields::set<descriptors::POROSITY>(lattice, geometry.getMaterialIndicator({5}), 0);

  auto& converter = lattice.getUnitConverter();
  lattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  switch (parameters.template get<parameters::OBJECTIVE_TYPE>()) {
    case ObjectiveType::DISSIPATION: {
      auto objectiveDomain = geometry.getMaterialIndicator({5,6});
      lattice.template addPostProcessor<stage::Evaluation>(objectiveDomain, meta::id<operators::TrackAverageFunctorO<DissipationF, opti::AVERAGE_J>>());
      if constexpr(std::is_same_v<CASE,MyCase>) {
        lattice.template addPostProcessor<stage::Evaluation>(objectiveDomain,
          meta::id<operators::TrackAverageFunctorO<functors::DerivativeF<DissipationF, descriptors::POPULATION, PrimalDynamics>, opti::AVERAGE_DJDF>>());
        lattice.template addPostProcessor<stage::Evaluation>(objectiveDomain,
          meta::id<operators::TrackAverageFunctorO<functors::DerivativeF<DissipationF, ControlledField, PrimalDynamics>, opti::AVERAGE_DJDALPHA<ControlledField>>>());
        lattice.template setParameter<parameters::FACTOR>(util::pow(converter.getPhysDeltaX(),2));
      }
      lattice.template setParameter<parameters::OMEGA>(converter.getLatticeRelaxationFrequency());
      lattice.template setParameter<parameters::DT>(converter.getPhysDeltaT());
      lattice.template setParameter<parameters::PHYS_CHAR_VISCOSITY>(converter.getPhysViscosity());
      lattice.template setParameter<parameters::DX>(converter.getPhysDeltaX());
      lattice.template setParameter<parameters::LATTICE_VISCOSITY>(converter.getLatticeViscosity());
      lattice.template setParameter<parameters::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
      break;
    }
    case ObjectiveType::FAR_FIELD_DRAG: {
      auto objectiveDomain = geometry.getMaterialIndicator({6});
      lattice.template addPostProcessor<stage::Evaluation>(objectiveDomain, meta::id<operators::TrackAverageFunctorO<FarFieldDragF, opti::AVERAGE_J>>());
      if constexpr(std::is_same_v<CASE,MyCase>) {
        lattice.template addPostProcessor<stage::Evaluation>(objectiveDomain,
          meta::id<operators::TrackAverageFunctorO<functors::DerivativeF<FarFieldDragF, descriptors::POPULATION, PrimalDynamics>, opti::AVERAGE_DJDF>>());
        lattice.template addPostProcessor<stage::Evaluation>(objectiveDomain,
          meta::id<operators::TrackAverageFunctorO<functors::DerivativeF<FarFieldDragF, ControlledField, PrimalDynamics>, opti::AVERAGE_DJDALPHA<ControlledField>>>());
        lattice.template setParameter<parameters::FACTOR>(util::pow(converter.getPhysDeltaX(),2));
      }
      lattice.template setParameter<parameters::CONVERSION_FORCE>(converter.getConversionFactorForce());
      lattice.template setParameter<parameters::LATTICE_REFERENCE_VELOCITY>(converter.getCharLatticeVelocity());
      lattice.template setParameter<parameters::LATTICE_REFERENCE_PRESSURE>(T(0.));
      break;
    }
  }
}

template <typename CASE>
void setInitialValues(CASE& myCase) {
  using T = typename CASE::value_t;
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

  const T maxPhysT   = parameters.template get<parameters::PHYS_MAX_T>();

  if (myCase.hasName(Controlled{}) || myCase.hasName(Derivatives{})) {
    std::size_t iTmaxStart = lattice.getUnitConverter().getLatticeTime(
      maxPhysT*parameters.template get<parameters::RAMP_UP_END_FRACTION>());
    std::size_t iTupdate   = lattice.getUnitConverter().getLatticeTime(
      parameters.template get<parameters::RAMP_UP_T>());
    if (iT % iTupdate == 0 && iT <= iTmaxStart) {
      T frac{};
      PolynomialStartScale<T,std::size_t>(iTmaxStart, T(1))(&frac, &iT);
      T maxVelocity = lattice.getUnitConverter().getCharPhysVelocity()*frac;
      momenta::setVelocity(lattice, geometry.getMaterialIndicator({3}), maxVelocity);
      lattice.template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
        ProcessingContext::Simulation);
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

void plotObjectiveDerivative(MyADfCase::value_t objective, MyADfCase::value_t trackAvgObj, MyADfCase::value_t physT) {
  using T = MyADfCase::value_t;
  static Gnuplot<T> gplot("j");
  gplot.setData(physT, {objective.v()}, {"j"}, "bottom right", {'l'});
  gplot.writePDF();
  static Gnuplot<T> gplot2("djdalpha");
  gplot2.setData(physT, {util::abs(objective.d(0)), util::abs(trackAvgObj.d(0))}, {"djdalpha", "djdalpha tr. avg."}, "bottom right", {'l', 'l'});
  gplot2.writePDF();
}

void plotVerticalWakeVelocity(MyCase::value_t vel_y, MyCase::value_t physT) {
  using T = MyCase::value_t;
  static Gnuplot<T> gplot("wake_vel_y");
  gplot.setData(physT, {vel_y}, {"vel y"}, "bottom right", {'l'});
  gplot.writePDF();
}

// Computes the pressure drop between the voxels before and after the cylinder
template <typename CASE>
void getResults(CASE& myCase,
                util::Timer<typename CASE::value_t>& timer,
                std::size_t iT)
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
    SuperLatticeField2D<T,DESCRIPTOR,descriptors::POPULATION> pop(lattice);
    fields::set<descriptors::AVERAGE_POPULATION>(lattice, geometry.getMaterialIndicator({1,2,3,4,5,6}), pop);
    lattice.setProcessingContext(ProcessingContext::Simulation);
    switch (parameters.template get<parameters::OBJECTIVE_TYPE>()) {
      case ObjectiveType::DISSIPATION: {
        auto objectiveDomain = geometry.getMaterialIndicator({5,6});
        auto objectiveF = makeWriteFunctorO<DissipationF, opti::AVERAGE_J>(lattice);
        objectiveF->template setParameter<parameters::OMEGA>(converter.getLatticeRelaxationFrequency());
        objectiveF->template setParameter<parameters::DT>(converter.getPhysDeltaT());
        objectiveF->template setParameter<parameters::PHYS_CHAR_VISCOSITY>(converter.getPhysViscosity());
        objectiveF->template setParameter<parameters::DX>(converter.getPhysDeltaX());
        objectiveF->template setParameter<parameters::LATTICE_VISCOSITY>(converter.getLatticeViscosity());
        objectiveF->template setParameter<parameters::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
        objectiveF->restrictTo(objectiveDomain);
        objectiveF->apply();
        if constexpr(std::is_same_v<CASE,MyCase>) {
          auto objectiveDF = makeWriteFunctorO<functors::DerivativeF<DissipationF, descriptors::POPULATION, PrimalDynamics>,
                                               opti::AVERAGE_DJDF>(lattice);
          objectiveDF->template setParameter<parameters::OMEGA>(converter.getLatticeRelaxationFrequency());
          objectiveDF->template setParameter<parameters::DT>(converter.getPhysDeltaT());
          objectiveDF->template setParameter<parameters::PHYS_CHAR_VISCOSITY>(converter.getPhysViscosity());
          objectiveDF->template setParameter<parameters::DX>(converter.getPhysDeltaX());
          objectiveDF->template setParameter<parameters::LATTICE_VISCOSITY>(converter.getLatticeViscosity());
          objectiveDF->template setParameter<parameters::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
          objectiveDF->template setParameter<parameters::FACTOR>(util::pow(converter.getPhysDeltaX(),2));
          objectiveDF->restrictTo(objectiveDomain);
          objectiveDF->apply();

          auto objectiveDA = makeWriteFunctorO<functors::DerivativeF<DissipationF, ControlledField, PrimalDynamics>,
                                               opti::AVERAGE_DJDALPHA<ControlledField>>(lattice);
          objectiveDA->template setParameter<parameters::OMEGA>(converter.getLatticeRelaxationFrequency());
          objectiveDA->template setParameter<parameters::DT>(converter.getPhysDeltaT());
          objectiveDA->template setParameter<parameters::PHYS_CHAR_VISCOSITY>(converter.getPhysViscosity());
          objectiveDA->template setParameter<parameters::DX>(converter.getPhysDeltaX());
          objectiveDA->template setParameter<parameters::LATTICE_VISCOSITY>(converter.getLatticeViscosity());
          objectiveDA->template setParameter<parameters::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
          objectiveDA->template setParameter<parameters::FACTOR>(util::pow(converter.getPhysDeltaX(),2));
          objectiveDA->restrictTo(objectiveDomain);
          objectiveDA->apply();
        }
        break;
      }
      case ObjectiveType::FAR_FIELD_DRAG: {
        auto objectiveDomain = geometry.getMaterialIndicator({6});
        auto objectiveF = makeWriteFunctorO<FarFieldDragF, opti::AVERAGE_J>(lattice);
        objectiveF->template setParameter<parameters::CONVERSION_FORCE>(converter.getConversionFactorForce());
        objectiveF->template setParameter<parameters::LATTICE_REFERENCE_VELOCITY>(converter.getCharLatticeVelocity());
        objectiveF->template setParameter<parameters::LATTICE_REFERENCE_PRESSURE>(T(0.));
        objectiveF->restrictTo(objectiveDomain);
        objectiveF->apply();
        if constexpr(std::is_same_v<CASE,MyCase>) {
          auto objectiveDF = makeWriteFunctorO<functors::DerivativeF<FarFieldDragF, descriptors::POPULATION, PrimalDynamics>,
                                               opti::AVERAGE_DJDF>(lattice);
          objectiveDF->template setParameter<parameters::CONVERSION_FORCE>(converter.getConversionFactorForce());
          objectiveDF->template setParameter<parameters::LATTICE_REFERENCE_VELOCITY>(converter.getCharLatticeVelocity());
          objectiveDF->template setParameter<parameters::LATTICE_REFERENCE_PRESSURE>(T(0.));
          objectiveDF->template setParameter<parameters::FACTOR>(util::pow(converter.getPhysDeltaX(),2));
          objectiveDF->restrictTo(objectiveDomain);
          objectiveDF->apply();
          auto objectiveDA = makeWriteFunctorO<functors::DerivativeF<FarFieldDragF, ControlledField, PrimalDynamics>,
                                               opti::AVERAGE_DJDALPHA<ControlledField>>(lattice);
          objectiveDA->template setParameter<parameters::CONVERSION_FORCE>(converter.getConversionFactorForce());
          objectiveDA->template setParameter<parameters::LATTICE_REFERENCE_VELOCITY>(converter.getCharLatticeVelocity());
          objectiveDA->template setParameter<parameters::LATTICE_REFERENCE_PRESSURE>(T(0.));
          objectiveDA->template setParameter<parameters::FACTOR>(util::pow(converter.getPhysDeltaX(),2));
          objectiveDA->restrictTo(objectiveDomain);
          objectiveDA->apply();
        }
        break;
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

    if constexpr (std::is_same_v<CASE, MyADfCase>) {
      switch (parameters.template get<parameters::OBJECTIVE_TYPE>()) {
        case ObjectiveType::DISSIPATION: {
          auto objectiveDomain = geometry.getMaterialIndicator({5,6});
          auto objectiveF = makeWriteFunctorO<DissipationF, opti::J>(lattice);
          objectiveF->template setParameter<parameters::OMEGA>(converter.getLatticeRelaxationFrequency());
          objectiveF->template setParameter<parameters::DT>(converter.getPhysDeltaT());
          objectiveF->template setParameter<parameters::PHYS_CHAR_VISCOSITY>(converter.getPhysViscosity());
          objectiveF->template setParameter<parameters::DX>(converter.getPhysDeltaX());
          objectiveF->template setParameter<parameters::LATTICE_VISCOSITY>(converter.getLatticeViscosity());
          objectiveF->template setParameter<parameters::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
          objectiveF->restrictTo(objectiveDomain);
          objectiveF->apply();
          T objective = integrateField<opti::J>(lattice, objectiveDomain, converter.getPhysDeltaX())[0];
          T objectiveAvg = integrateField<opti::AVERAGE_J>(lattice, objectiveDomain, converter.getPhysDeltaX())[0];
          plotObjectiveDerivative(objective, objectiveAvg, converter.getPhysTime(iT));
          break;
        }
        case ObjectiveType::FAR_FIELD_DRAG: {
          auto objectiveDomain = geometry.getMaterialIndicator({6});
          auto objectiveF = makeWriteFunctorO<FarFieldDragF, opti::J>(lattice);
          objectiveF->template setParameter<parameters::CONVERSION_FORCE>(converter.getConversionFactorForce());
          objectiveF->template setParameter<parameters::LATTICE_REFERENCE_VELOCITY>(converter.getCharLatticeVelocity());
          objectiveF->template setParameter<parameters::LATTICE_REFERENCE_PRESSURE>(T(0.));
          objectiveF->restrictTo(objectiveDomain);
          objectiveF->apply();
          T objective = integrateField<opti::J>(lattice, objectiveDomain, converter.getPhysDeltaX())[0];
          T objectiveAvg = integrateField<opti::AVERAGE_J>(lattice, objectiveDomain, converter.getPhysDeltaX())[0];
          plotObjectiveDerivative(objective, objectiveAvg, converter.getPhysTime(iT));
          break;
        }
      }
    } else {
      if (myCase.hasName(Controlled{})) {
        lattice.setProcessingContext(ProcessingContext::Evaluation);
        SuperLatticePhysVelocity2D<T,DESCRIPTOR> velocity(lattice, converter);
        AnalyticalFfromSuperF2D<T,T> vel(velocity);
        T output[2], input[2];
        input[0] = parameters.template get<parameters::VELOCITY_SAMPLING_POS_X>();
        input[1] = parameters.template get<parameters::VELOCITY_SAMPLING_POS_Y>();
        vel(output, input);
        std::cout << "Wake velocity: " << output[0] << " | " << output[1] << std::endl;
        plotVerticalWakeVelocity(output[1], converter.getPhysTime(iT));
      }
    }
  }

  const std::size_t vtkIter  = converter.getLatticeTime(parameters.template get<parameters::PHYS_VTK_ITER_T>());
  if (myCase.hasName(Controlled{})) {
    SuperVTMwriter2D<T> vtmWriter("cylinder2d");
    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(lattice, converter);
    SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure(lattice, converter);
    SuperLatticeField2D<T,DESCRIPTOR,descriptors::POROSITY> porosityF(lattice);
    SuperLatticeField2D<T,DESCRIPTOR,descriptors::POPULATION> pop(lattice);
    SuperLatticeField2D<T,DESCRIPTOR,descriptors::AVERAGE_POPULATION> popAvg(lattice);
    SuperLatticeField2D<T,DESCRIPTOR,opti::AVERAGE_J> jFAvg(lattice);
    SuperLatticeField2D<T,DESCRIPTOR,opti::AVERAGE_DJDF> djdfFAvg(lattice);
    SuperLatticeField2D<T,DESCRIPTOR,opti::AVERAGE_DJDALPHA<ControlledField>> djdaFAvg(lattice);
    djdaFAvg.getName() = "djdaFAvg";

    vtmWriter.addFunctor(velocity);
    vtmWriter.addFunctor(pressure);
    vtmWriter.addFunctor(porosityF);
    vtmWriter.addFunctor(pop);
    vtmWriter.addFunctor(popAvg);
    vtmWriter.addFunctor(jFAvg);
    vtmWriter.addFunctor(djdfFAvg);
    vtmWriter.addFunctor(djdaFAvg);

    if (iT == 0) {
      vtmWriter.createMasterFile();
    }

    // Writes the vtk files
    if (iT % vtkIter == 0 && iT > 0) {
      lattice.setProcessingContext(ProcessingContext::Evaluation);
      vtmWriter.write(iT);
    }
  } else {
    SuperVTMwriter2D<T> vtmWriter("adjointCylinder2d");
    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(lattice, converter);
    SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure(lattice, converter);
    SuperLatticeField2D<T,DESCRIPTOR,descriptors::POROSITY> porosityF(lattice);
    SuperLatticeField2D<T,DESCRIPTOR,opti::DJDF> djdfF(lattice);
    SuperLatticeField2D<T,DESCRIPTOR,opti::F> fF(lattice);
    SuperLatticeField2D<T,DESCRIPTOR,descriptors::POPULATION> adjointsF(lattice);
    SuperLatticeField2D<T,DESCRIPTOR,momenta::FixedVelocityMomentumGeneric::VELOCITY> boundaryVel(lattice);

    vtmWriter.addFunctor(velocity);
    vtmWriter.addFunctor(pressure);
    vtmWriter.addFunctor(porosityF);
    vtmWriter.addFunctor(djdfF);
    vtmWriter.addFunctor(fF);
    vtmWriter.addFunctor(adjointsF);
    vtmWriter.addFunctor(boundaryVel);

    if (iT == 0) {
      vtmWriter.createMasterFile();
    }

    // Writes the vtk files
    if (iT % vtkIter == 0 && iT > 0) {
      lattice.setProcessingContext(ProcessingContext::Evaluation);
      vtmWriter.write(iT);
    }
  }
}

template <typename CASE>
void simulate(CASE& myCase){
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

    if((myCase.hasName(Controlled{}) || myCase.hasName(Derivatives{})) && iT >= iTstartAvg) {
      myCase.getLattice(NavierStokes{}).executePostProcessors(stage::Evaluation());
    }

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, timer, iT);
  }

  timer.stop();
  timer.printSummary();
}

void setInitialControl(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& lattice = controlledCase.getLattice(NavierStokes{});
  auto& geometry = controlledCase.getGeometry();
  auto& control = optiCase.getController();
  auto& converter = lattice.getUnitConverter();
  auto& parameters = controlledCase.getParameters();

  // Intialize controlled field in superLattice
  T porosity = projection::permeabilityToPorosity(parameters.get<parameters::INITIAL_CONTROL_SCALAR>(), converter);
  std::cout << "GridTerm: " << projection::gridTerm(converter) << std::endl;
  std::cout << "Initial porosity: " << porosity << std::endl;
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
  auto& geometry = optiCase.getCase(Controlled{}).getGeometry();
  auto& parameters = optiCase.getCase(Controlled{}).getParameters();
  switch (parameters.template get<parameters::OBJECTIVE_TYPE>()) {
    case ObjectiveType::DISSIPATION: {
      auto objectiveDomain = geometry.getMaterialIndicator({5,6});
      return integrateField<opti::AVERAGE_J>(controlledLattice, objectiveDomain, converter.getPhysDeltaX())[0];
      break;
    }
    case ObjectiveType::FAR_FIELD_DRAG: {
      auto objectiveDomain = geometry.getMaterialIndicator({6});
      return integrateField<opti::AVERAGE_J>(controlledLattice, objectiveDomain, converter.getPhysDeltaX())[0];
      break;
    }
  }
}

MyCase::value_t solvePrimalAndComputeObjective(MyOptiCase& optiCase) {
  auto& controlledCase = optiCase.getCase(Controlled{});
  // Reset prior simulation lattice
  controlledCase.resetLattices();

  // Prepare new lattices
  prepareLattice(controlledCase);
  setInitialValues(controlledCase);

  // Set updated controls in the simulation
  applyControl(optiCase);

  // Execute simulation
  simulate(controlledCase);

  // Compute Objective value from simulation results
  return computeObjective(optiCase);
}

void prepareAdjointLattice(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& adjointLattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& geometry = optiCase.getCase(Adjoint{}).getGeometry();

  adjointLattice.setUnitConverter(controlledLattice.getUnitConverter());

  // Define dual dynamics
  dynamics::set<olb::Dual<PorousBGKdynamics<T,DESCRIPTOR>>,T,DESCRIPTOR>(adjointLattice, geometry.getMaterialIndicator({1,5,6}));
  boundary::set<boundary::ZeroPopulation>(adjointLattice, geometry.getMaterialIndicator({2,3,4}));
  adjointLattice.template setParameter<descriptors::OMEGA>(adjointLattice.getUnitConverter().getLatticeRelaxationFrequency());
  adjointLattice.revertStreaming();
}

void setAdjointInitialValues(MyOptiCase& optiCase) {
  auto& adjointCase = optiCase.getCase(Adjoint{});
  auto& adjointLattice = adjointCase.getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& converter = controlledLattice.getUnitConverter();

  setInitialValues(adjointCase);

  copyFields<opti::AVERAGE_DJDF, opti::DJDF>(controlledLattice, adjointLattice);
  copyFields<ControlledField, ControlledField>(controlledLattice, adjointLattice);
  // Copy adjoint source term into tmp field to apply ramp-up
  copyFields<opti::DJDF, opti::DJDF_TMP>(adjointLattice, adjointLattice);
  copyFields<momenta::FixedVelocityMomentumGeneric::VELOCITY, momenta::FixedVelocityMomentumGeneric::VELOCITY>(controlledLattice, adjointLattice);
  copyFields<descriptors::AVERAGE_POPULATION, opti::F>(controlledLattice, adjointLattice);
}

std::vector<MyCase::value_t> computeDerivative(MyOptiCase& optiCase) {
  auto& adjointLattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& converter = controlledLattice.getUnitConverter();
  auto& geometry = optiCase.getCase(Controlled{}).getGeometry();
  auto& parameters = optiCase.getCase(Adjoint{}).getParameters();
  auto optimalityO = makeWriteFunctorO<functors::OptimalityF<PrimalDynamics, ControlledField>, opti::SENSITIVITY<ControlledField>>(adjointLattice);
  optiCase.getController().template setUpdatedProjectionDerivativesOnField<opti::DPROJECTIONDALPHA<ControlledField>>(adjointLattice);
  adjointLattice.template setProcessingContext<Array<opti::DPROJECTIONDALPHA<ControlledField>>>(ProcessingContext::Simulation);
  auto dCdalphaO = makeWriteFunctorO<functors::DerivativeF<functors::CollisionUsingAveragePopulationF<PrimalDynamics>, ControlledField, PrimalDynamics>, opti::DCDALPHA<ControlledField>>(controlledLattice);
  dCdalphaO->template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  dCdalphaO->template setParameter<parameters::FACTOR>(1.0);
  dCdalphaO->apply();

  copyFields<opti::DCDALPHA<ControlledField>,opti::DCDALPHA<ControlledField>>(controlledLattice, adjointLattice);
  copyFields<opti::AVERAGE_DJDALPHA<ControlledField>,opti::DJDALPHA<ControlledField>>(controlledLattice, adjointLattice);
  optimalityO->apply();

  // Return serial vector containing total derivatives of objective regarding controls
  switch (parameters.template get<parameters::OBJECTIVE_TYPE>()) {
    case ObjectiveType::DISSIPATION: {
      auto objectiveDomain = geometry.getMaterialIndicator({5,6});
      return getSerializedFromField<opti::SENSITIVITY<ControlledField>>(adjointLattice, objectiveDomain);
      break;
    }
    case ObjectiveType::FAR_FIELD_DRAG: {
      auto objectiveDomain = geometry.getMaterialIndicator({6});
      return getSerializedFromField<opti::SENSITIVITY<ControlledField>>(adjointLattice, objectiveDomain);
      break;
    }
  }
}

void getOptiResults(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& adjointLattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& converter = controlledLattice.getUnitConverter();
  std::size_t iT = optiCase.getOptimizationStep();

  SuperVTMwriter2D<T> vtmWriter("cylinderTopology2d");
  SuperLatticePhysVelocity2D velocityF(adjointLattice, converter);
  SuperLatticePhysPressure2D pressureF(adjointLattice, converter);
  SuperLatticeField2D<T,DESCRIPTOR,descriptors::POROSITY> porosityF(adjointLattice);
  SuperLatticeField2D<T,DESCRIPTOR,opti::J> jF(controlledLattice);
  SuperLatticeField2D<T,DESCRIPTOR,opti::DJDF> djdfF(adjointLattice);
  SuperLatticeField2D<T,DESCRIPTOR,opti::SENSITIVITY<ControlledField>> djdalpha(adjointLattice);
  SuperLatticeField2D<T,DESCRIPTOR,opti::DPROJECTIONDALPHA<ControlledField>> dproj(adjointLattice);
  vtmWriter.addFunctor(velocityF);
  vtmWriter.addFunctor(pressureF);
  vtmWriter.addFunctor(porosityF);
  vtmWriter.addFunctor(jF);
  vtmWriter.addFunctor(djdfF);
  vtmWriter.addFunctor(djdalpha);
  djdalpha.getName() = "djdalpha";
  vtmWriter.addFunctor(dproj);
  dproj.getName() = "dproj";

  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  // Writes the VTK files
  if (iT >= 0) {
    controlledLattice.setProcessingContext(ProcessingContext::Evaluation);
    adjointLattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }
}

auto getShiftedControl(MyOptiCase& optiCase, MyCase::value_t shiftedLocalPorosity) {
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& lattice = controlledCase.getLattice(NavierStokes{});
  auto& geometry = controlledCase.getGeometry();
  auto& control = optiCase.getController();
  auto& parameters = controlledCase.getParameters();

  setInitialControl(optiCase);
  Vector<MyCase::value_t,2> samplePhysPos(parameters.get<parameters::SAMPLE_POS>());

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
  auto& converter = sLattice.getUnitConverter();

  SuperLatticeField2D<T,DESCRIPTOR,opti::SENSITIVITY<ControlledField>> gradF(sLattice);
  AnalyticalFfromSuperF2D<T,T> grad(gradF);

  T output[1], input[2];
  input[0] = parameters.get<parameters::SAMPLE_POS>()[0];
  input[1] = parameters.get<parameters::SAMPLE_POS>()[1];
  grad(output, input);
  std::cout << "Adjoint Gradient at " << parameters.get<parameters::SAMPLE_POS>()
                                      << ": (" << output[0]
                                      << ")" << std::endl;
  parameters.set<parameters::DISCRETE_ADJOINT_GRADIENTS>(output[0]);
}

std::vector<MyCase::value_t> solveAdjointAndComputeDerivative(MyOptiCase& optiCase) {
  auto& adjointCase = optiCase.getCase(Adjoint{});
  // Reset prior simulation lattice
  adjointCase.resetLattices();

  // Prepare new lattices
  prepareAdjointLattice(optiCase);
  setAdjointInitialValues(optiCase);

  // Execute simulation
  simulate(adjointCase);

  // Compute Derivative value from simulation results
  auto derivatives = computeDerivative(optiCase);

  // Optimization related IO
  getOptiResults(optiCase);

  return derivatives;
}

void computeDerivativesWithFDQ(MyOptiCase& optiCase, MyCase::value_t objective) {
  using T = MyCase::value_t;
  auto& lattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();
  auto& parameters = optiCase.getCase(Adjoint{}).getParameters();

  T gradient(0.0);
  const T permeability = parameters.get<parameters::INITIAL_CONTROL_SCALAR>();
  const T shiftedPermeability = parameters.get<parameters::CONTROL_SHIFT>() +
                                parameters.get<parameters::INITIAL_CONTROL_SCALAR>();
  const T porosity = projection::permeabilityToPorosity(permeability, converter);
  const T shiftedPorosity = projection::permeabilityToPorosity(shiftedPermeability, converter);
  auto shiftedControl = getShiftedControl(optiCase, shiftedPorosity);
  auto shiftedObjective = optiCase.computeObjective(shiftedControl, 0);
  gradient = (shiftedObjective - objective) / (shiftedPermeability - permeability);
  optiCase.getCase(Controlled{}).getParameters().set<parameters::FDQ_GRADIENTS>(gradient);
}

void plotDerivatives(MyCase::value_t adjoint, MyCase::value_t FDQ, MyCase::value_t ADf, unsigned index) {
  using T = MyCase::value_t;

  static Gnuplot<T> gplot("derivatives");
  gplot.setData(T(index), {adjoint, FDQ, ADf}, {"adjoints", "FDQ", "ADf"}, "bottom right", {'l', 'l', 'l'});
  gplot.writePDF();
}

void plotAveragedDerivatives(MyCase::value_t adjoint, MyCase::value_t ADf, unsigned index) {
  using T = MyCase::value_t;

  static Gnuplot<T> gplot("averaged_derivatives");
  gplot.setData(T(index), {adjoint, ADf}, {"adjoints", "ADf"}, "bottom right", {'l', 'l'});
  gplot.writePDF();
}

void sampleDerivativesWithADf(MyADfCase& myCase) {
  myCase.resetLattices();
  prepareLattice(myCase);
  setInitialValues(myCase);

  using T = MyADfCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  // reset the controlled field
  T permeability = parameters.get<parameters::INITIAL_CONTROL_SCALAR>();
  T porosity = projection::permeabilityToPorosity(permeability, converter);
  fields::set<ControlledField>(lattice, geometry.getMaterialIndicator({5}), porosity);
  // seeding the permeability for that single cell
  permeability.setDiffVariable(0);
  Vector<T,1> control = projection::permeabilityToPorosity(permeability, converter);
  // set the seeded control in the controlled field
  Vector<MyADfCase::value_t,2> samplePhysPos(parameters.get<parameters::SAMPLE_POS>());
  if (auto samplePos = myCase.getMesh().getCuboidDecomposition().getLatticeR(samplePhysPos)) {
    if (myCase.getMesh().getLoadBalancer().isLocal((*samplePos)[0])) {
      auto cell = lattice.get(*samplePos);
      cell.template setField<ControlledField>(control);
    }
  }
  lattice.template setProcessingContext<Array<ControlledField>>(ProcessingContext::Simulation);

  simulate(myCase);
  // Compute objective
  switch (parameters.template get<parameters::OBJECTIVE_TYPE>()) {
    case ObjectiveType::DISSIPATION: {
      auto objectiveF = makeWriteFunctorO<DissipationF, opti::J>(lattice);
      auto objectiveDomain = geometry.getMaterialIndicator({5,6});
      objectiveF->template setParameter<parameters::OMEGA>(converter.getLatticeRelaxationFrequency());
      objectiveF->template setParameter<parameters::DT>(converter.getPhysDeltaT());
      objectiveF->template setParameter<parameters::PHYS_CHAR_VISCOSITY>(converter.getPhysViscosity());
      objectiveF->template setParameter<parameters::DX>(converter.getPhysDeltaX());
      objectiveF->template setParameter<parameters::LATTICE_VISCOSITY>(converter.getLatticeViscosity());
      objectiveF->template setParameter<parameters::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
      objectiveF->restrictTo(objectiveDomain);
      objectiveF->apply();
      T objective = integrateField<opti::J>(lattice, objectiveDomain, converter.getPhysDeltaX())[0];
      parameters.set<parameters::ADF_GRADIENTS>(T(objective.d(0)));
      break;
    }
    case ObjectiveType::FAR_FIELD_DRAG: {
      auto objectiveF = makeWriteFunctorO<FarFieldDragF, opti::J>(lattice);
      auto objectiveDomain = geometry.getMaterialIndicator({6});
      objectiveF->template setParameter<parameters::CONVERSION_FORCE>(converter.getConversionFactorForce());
      objectiveF->template setParameter<parameters::LATTICE_REFERENCE_VELOCITY>(converter.getCharLatticeVelocity());
      objectiveF->template setParameter<parameters::LATTICE_REFERENCE_PRESSURE>(T(0.));
      objectiveF->restrictTo(objectiveDomain);
      objectiveF->apply();
      T objective = integrateField<opti::J>(lattice, objectiveDomain, converter.getPhysDeltaX())[0];
      parameters.set<parameters::ADF_GRADIENTS>(T(objective.d(0)));
      break;
    }
  }
}

void sampleAveragedDerivativesWithADf(MyADfCase& myCase) {
  myCase.resetLattices();
  prepareLattice(myCase);
  setInitialValues(myCase);

  OstreamManager clout(std::cout, "simulate");
  using T = MyADfCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  // reset the controlled field
  T permeability = parameters.get<parameters::INITIAL_CONTROL_SCALAR>();
  T porosity = projection::permeabilityToPorosity(permeability, converter);
  fields::set<ControlledField>(lattice, geometry.getMaterialIndicator({5}), porosity);

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

    if (iT == iTstartAvg) {
      // seeding the permeability for that single cell
      permeability.setDiffVariable(0);
      Vector<T,1> control = projection::permeabilityToPorosity(permeability, converter);
      // set the seeded control in the controlled field
      Vector<MyADfCase::value_t,2> samplePhysPos(parameters.get<parameters::SAMPLE_POS>());
      if (auto samplePos = myCase.getMesh().getCuboidDecomposition().getLatticeR(samplePhysPos)) {
        if (myCase.getMesh().getLoadBalancer().isLocal((*samplePos)[0])) {
          auto cell = lattice.get(*samplePos);
          cell.template setField<ControlledField>(control);
        }
      }
      lattice.template setProcessingContext<Array<ControlledField>>(ProcessingContext::Simulation);
    }

    if((myCase.hasName(Controlled{}) || myCase.hasName(Derivatives{})) && iT >= iTstartAvg) {
      myCase.getLattice(NavierStokes{}).executePostProcessors(stage::Evaluation());
    }

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, timer, iT);
  }

  timer.stop();
  timer.printSummary();

  // set in parameters
  switch (parameters.template get<parameters::OBJECTIVE_TYPE>()) {
    case ObjectiveType::DISSIPATION: {
      auto objectiveDomain = geometry.getMaterialIndicator({5,6});
      T objective = integrateField<opti::AVERAGE_J>(lattice, objectiveDomain, converter.getPhysDeltaX())[0];
      parameters.set<parameters::ADF_GRADIENTS>(T(objective.d(0)));
      break;
    }
    case ObjectiveType::FAR_FIELD_DRAG: {
      auto objectiveDomain = geometry.getMaterialIndicator({6});
      T objective = integrateField<opti::AVERAGE_J>(lattice, objectiveDomain, converter.getPhysDeltaX())[0];
      parameters.set<parameters::ADF_GRADIENTS>(T(objective.d(0)));
      break;
    }
  }
}

auto sampleCircle(MyCase::value_t cx, MyCase::value_t cy, MyCase::value_t r, MyCase::value_t dx, int nSamples) {
  using T = MyCase::value_t;
  std::vector<std::vector<T>> points;
  points.reserve(nSamples);

  const T twoPi = 2.0 * M_PI;
  const T  rShifted = r - dx;
  for (int i = 0; i < nSamples; ++i) {
    const T theta = 0.5*M_PI - twoPi * i / nSamples;
    std::vector<T> point({0.0, 0.0});
    point[0] = cx + rShifted*util::cos(theta);
    point[1] = cy + rShifted*util::sin(theta);
    points.push_back(point);
  }
  return points;
}

void compareFinalGradients(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  using U = util::ADf<T, 1>;
  auto& myCaseParameters = optiCase.getCase(Controlled{}).getParameters();
  auto& adjointCase = optiCase.getCase(Adjoint{});
  auto& adfCase = optiCase.getCase(Derivatives{});
  auto& myADfCaseParameters = adfCase.getParameters();
  auto& controls = optiCase.getController().get();

  // Evaluate objective
  T objective = optiCase.computeObjective(controls);
  // Compute gradients with adjoints
  optiCase.computeDerivatives(controls);

  // Evaluating gradients on the cylinder surface
  const T dx  = myCaseParameters.get<parameters::RADIUS_CYLINDER>()/myCaseParameters.get<parameters::RESOLUTION>();
  const T radius = myCaseParameters.get<parameters::RADIUS_CYLINDER>();
  auto samplePoints = sampleCircle(1.2, 0.4, radius, dx, myCaseParameters.get<parameters::NUMBER_SAMPLES>());

  unsigned index = 0;
  for (auto point : samplePoints) {
    myCaseParameters.set<parameters::SAMPLE_POS>(point);
    myADfCaseParameters.set<parameters::SAMPLE_POS>({U(point[0]), U(point[1])});
    sampleDerivativesWithAdjoints(adjointCase);
    computeDerivativesWithFDQ(optiCase, objective);
    sampleDerivativesWithADf(adfCase);

    std::cout << "Summary:" << std::endl;
    std::cout << "Dis.   Gradient at " << optiCase.getCase(Controlled{}).getParameters().get<parameters::SAMPLE_POS>()
                                       << ": " << optiCase.getCase(Controlled{}).getParameters().get<parameters::DISCRETE_ADJOINT_GRADIENTS>()
                                       << std::endl;
    std::cout << "FDQ.   Gradient at " << optiCase.getCase(Controlled{}).getParameters().get<parameters::SAMPLE_POS>()
                                       << ": " << optiCase.getCase(Controlled{}).getParameters().get<parameters::FDQ_GRADIENTS>()
                                       << std::endl;
    std::cout << "ADF.   Gradient at " << adfCase.getParameters().get<parameters::SAMPLE_POS>()
                                       << ": " << adfCase.getParameters().get<parameters::ADF_GRADIENTS>()
                                       << std::endl;
    plotDerivatives(optiCase.getCase(Controlled{}).getParameters().get<parameters::DISCRETE_ADJOINT_GRADIENTS>(),
                    optiCase.getCase(Controlled{}).getParameters().get<parameters::FDQ_GRADIENTS>(),
                    adfCase.getParameters().get<parameters::ADF_GRADIENTS>(), index);
    ++index;
  }
}

void compareAveragedGradients(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  using U = util::ADf<T, 1>;
  auto& myCaseParameters = optiCase.getCase(Controlled{}).getParameters();
  auto& adjointCase = optiCase.getCase(Adjoint{});
  auto& adfCase = optiCase.getCase(Derivatives{});
  auto& myADfCaseParameters = adfCase.getParameters();
  auto& converter = optiCase.getCase(Controlled{}).getLattice(NavierStokes{}).getUnitConverter();
  auto& controls = optiCase.getController().get();

  // Evaluate objective
  T objective = optiCase.computeObjective(controls);
  // Compute gradients with adjoints
  optiCase.computeDerivatives(controls);

  std::cout << "Compare averaged gradients..." << std::endl;
  T frac = myCaseParameters.template get<parameters::AVERAGING_START_FRACTION>();

  std::cout << "Averaging start fraction: " << std::setprecision(16) << frac << std::endl;
  const T physMaxT = myCaseParameters.template get<parameters::PHYS_MAX_T>();
  const std::size_t iTmax = converter.getLatticeTime(physMaxT);
  const std::size_t iTstartAvg = converter.getLatticeTime(myCaseParameters.template get<parameters::PHYS_MAX_T>() *
                                                          myCaseParameters.template get<parameters::AVERAGING_START_FRACTION>());
  std::cout << "Discrete averaging time steps: " << iTmax - iTstartAvg << std::endl;

  // Evaluating gradients on the cylinder surface
  const T dx  = myCaseParameters.get<parameters::RADIUS_CYLINDER>()/myCaseParameters.get<parameters::RESOLUTION>();
  const T radius = myCaseParameters.get<parameters::RADIUS_CYLINDER>();
  auto samplePoints = sampleCircle(1.2, 0.4, radius, dx, myCaseParameters.get<parameters::NUMBER_SAMPLES>());

  unsigned index = 0;
  for (auto point : samplePoints) {
    myCaseParameters.set<parameters::SAMPLE_POS>(point);
    myADfCaseParameters.set<parameters::SAMPLE_POS>({U(point[0]), U(point[1])});
    sampleDerivativesWithAdjoints(adjointCase);
    sampleAveragedDerivativesWithADf(adfCase);

    std::cout << "Summary:" << std::endl;
    std::cout << "Averaged Dis.   Gradient at " << optiCase.getCase(Controlled{}).getParameters().get<parameters::SAMPLE_POS>()
                                                << ": " << optiCase.getCase(Controlled{}).getParameters().get<parameters::DISCRETE_ADJOINT_GRADIENTS>()
                                                << std::endl;
    std::cout << "Averaged ADF.   Gradient at " << adfCase.getParameters().get<parameters::SAMPLE_POS>()
                                                << ": " << adfCase.getParameters().get<parameters::ADF_GRADIENTS>()
                                                << std::endl;
    plotAveragedDerivatives(optiCase.getCase(Controlled{}).getParameters().get<parameters::DISCRETE_ADJOINT_GRADIENTS>(),
                            adfCase.getParameters().get<parameters::ADF_GRADIENTS>(), index);
    ++index;
  }
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);

  using T = MyCase::value_t;
  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION                  >(21);
    myCaseParameters.set<REYNOLDS                    >(100.);
    myCaseParameters.set<PHYS_MAX_T                  >(300.);
    myCaseParameters.set<DOMAIN_EXTENT_X             >(3.6);
    myCaseParameters.set<DOMAIN_EXTENT_Y             >(.81);
    myCaseParameters.set<PHYS_CHAR_VELOCITY          >(0.2);
    myCaseParameters.set<LATTICE_RELAXATION_TIME     >(0.55);
    myCaseParameters.set<PHYS_CHAR_DENSITY           >(1.0);
    myCaseParameters.set<RADIUS_CYLINDER             >(0.05);
    myCaseParameters.set<CENTER_CYLINDER_X           >(1.2);
    myCaseParameters.set<CENTER_CYLINDER_Y           >(0.4);
    myCaseParameters.set<PHYS_VTK_ITER_T             >(10);
    myCaseParameters.set<PHYS_STAT_ITER_T            >(2.5);
    myCaseParameters.set<AVERAGING_START_FRACTION    >(0.8);
    myCaseParameters.set<RAMP_UP_T                   >(0.01);
    myCaseParameters.set<RAMP_UP_END_FRACTION        >(0.4);
    myCaseParameters.set<ADJOINT_RAMP_UP_UPDATE      >(0.01);
    myCaseParameters.set<ADJOINT_RAMP_UP_END_FRACTION>(0.5);
    myCaseParameters.set<SAMPLE_POS                  >({1.2, 0.45});
    myCaseParameters.set<INITIAL_CONTROL_SCALAR      >(1.e-6);
    myCaseParameters.set<CONTROL_SHIFT               >(1.e-10);
    myCaseParameters.set<NUMBER_SAMPLES              >(16);
    myCaseParameters.set<OBJECTIVE_DOMAIN_X          >(2.2);
    myCaseParameters.set<OBJECTIVE_DOMAIN_Y          >(0.8);
    myCaseParameters.set<VELOCITY_SAMPLING_POS_X     >(2.4);
    myCaseParameters.set<VELOCITY_SAMPLING_POS_Y     >(0.4);
    myCaseParameters.set<OBJECTIVE_TYPE              >(ObjectiveType::FAR_FIELD_DRAG);
  }
  auto myADfCaseParameters = myCaseParameters.copyAs<MyADfCase::value_t>();

  Mesh mesh = createMesh<MyCase>(myCaseParameters);
  MyCase myCase(myCaseParameters, mesh);
  prepareGeometry(myCase);
  prepareLattice(myCase);

  MyCase adjointCase(myCaseParameters, mesh);
  prepareGeometry(adjointCase);

  Mesh adfmesh = createMesh<MyADfCase>(myADfCaseParameters);
  MyADfCase adfCase(myADfCaseParameters, adfmesh);
  adfCase.setName(Controlled{}.name);
  prepareGeometry(adfCase);
  prepareLattice(adfCase);

  MyOptiCase optiCase;
  optiCase.setCase<Controlled>(myCase);
  optiCase.setCase<Adjoint>(adjointCase);
  optiCase.setCase<Derivatives>(adfCase);

  setInitialControl(optiCase);
  optiCase.setObjective(solvePrimalAndComputeObjective);
  optiCase.setDerivative(solveAdjointAndComputeDerivative);

  //compareFinalGradients(optiCase);
  compareAveragedGradients(optiCase);
}
