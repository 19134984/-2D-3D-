/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
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

// This app is based on the work in doi: 10.1016/j.jcp.2026.115151.
// The datasets can be downloaded from https://bwsyncandshare.kit.edu/f/6013224286

// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::opti;

using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::DualPorousD3Q19Descriptor>
>;
using MyOptiCase = OptiCaseAdjoint<
  Controlled, MyCase,
  Adjoint, MyCase,
  Reference, MyCase
>;

using BulkDynamics = PorousBGKdynamics<MyCase::value_t,
                                       MyCase::descriptor_t_of<NavierStokes>>;
using ObjectiveF = functors::L2DistanceF<functors::VelocityF>;
using ControlledField = descriptors::POROSITY;

namespace olb::parameters{

struct DOMAIN_CENTER1           : public descriptors::FIELD_BASE<0,1>{ };
struct DOMAIN_CENTER2           : public descriptors::FIELD_BASE<0,1>{ };
struct DOMAIN_RADIUS            : public descriptors::FIELD_BASE<1>{ };
struct DESIGN_DOMAIN_CENTER1    : public descriptors::FIELD_BASE<0,1>{ };
struct DESIGN_DOMAIN_CENTER2    : public descriptors::FIELD_BASE<0,1>{ };
struct DESIGN_DOMAIN_RADIUS     : public descriptors::FIELD_BASE<1>{ };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const Vector center1 = parameters.get<parameters::DOMAIN_CENTER1>();
  const Vector center2 = parameters.get<parameters::DOMAIN_CENTER2>();
  const T radius = parameters.get<parameters::DOMAIN_RADIUS>();
  IndicatorCylinder3D<T> cylinder(center1, center2, radius);

  const T physDeltaX = radius*2. / parameters.get<parameters::RESOLUTION>();
  Mesh<T,MyCase::d> mesh(cylinder, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& superGeometry = myCase.getGeometry();
  superGeometry.setWriteIncrementalVTK(true);
  const Vector center1 = parameters.get<parameters::DOMAIN_CENTER1>();
  const Vector center2 = parameters.get<parameters::DOMAIN_CENTER2>();
  const T radius = parameters.get<parameters::DOMAIN_RADIUS>();
  IndicatorCylinder3D<T> cylinder(center1, center2, radius);

  superGeometry.rename(0, 2, cylinder);
  superGeometry.rename(2, 1, {1,1,1});

  const T dz = radius / parameters.get<parameters::RESOLUTION>();
  Vector<T,3> center1_minus_dz = cylinder.getCenter1();
  center1_minus_dz[2] -= dz;
  Vector<T,3> center1_plus_dz = center1_minus_dz;
  center1_plus_dz[2] += 2.*dz;
  Vector<T,3> center2_minus_dz = cylinder.getCenter2();
  center2_minus_dz[2] -= dz;
  Vector<T,3> center2_plus_dz = center2_minus_dz;
  center2_plus_dz[2] += 2.*dz;

  IndicatorCylinder3D<T> inflow(center2_minus_dz, center2_plus_dz, radius);
  superGeometry.rename(2, 3, 1, inflow);
  superGeometry.clean();

  IndicatorCylinder3D<T> outflow(center1_minus_dz, center1_plus_dz, radius);
  superGeometry.rename(2, 4, 1, outflow);

  // Indicators for material 6 (designDomain)
  const Vector designCenter1 = parameters.get<parameters::DESIGN_DOMAIN_CENTER1>();
  const Vector designCenter2 = parameters.get<parameters::DESIGN_DOMAIN_CENTER2>();
  const T designRadius = parameters.get<parameters::DESIGN_DOMAIN_RADIUS>();
  IndicatorCylinder3D<T> designDomain(designCenter1, designCenter2, designRadius);
  superGeometry.rename(1, 6, designDomain);

  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
}

void prepareLattice(MyCase& myCase) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  sLattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>>(
    parameters.get<parameters::RESOLUTION>(),               // resolution
    parameters.get<parameters::LATTICE_RELAXATION_TIME>(),  // relaxationTime
    parameters.get<parameters::DOMAIN_RADIUS>() * 2.,       // charPhysLength: reference length of simulation geometry in [m]
    parameters.get<parameters::PHYS_CHAR_VELOCITY>(),       // charPhysVelocity: highest expected velocity during simulation in [m/s]
    parameters.get<parameters::PHYS_CHAR_VISCOSITY>(),      // physViscosity: physical kinematic viscosity in [m^2/s]
    parameters.get<parameters::PHYS_CHAR_DENSITY>()         // physDensity: physical density [kg/m^3]
  );
  auto& converter = sLattice.getUnitConverter();
  converter.print();

  auto bulkIndicator = superGeometry.getMaterialIndicator({1,3,4,6});
  dynamics::set<PorousBGKdynamics>(sLattice, bulkIndicator);

  boundary::set<boundary::BounceBack>(sLattice, superGeometry, 2);
  boundary::set<boundary::InterpolatedVelocity>(sLattice, superGeometry.getMaterialIndicator({3}));
  boundary::set<boundary::InterpolatedPressure>(sLattice, superGeometry.getMaterialIndicator({4}));

  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  auto& converter = sLattice.getUnitConverter();
  auto bulkIndicator = superGeometry.getMaterialIndicator({1,2,3,4,6});

  STLreader<T> stlReader("foam_transformed.stl", converter.getConversionFactorLength(), 1.0);
  SuperIndicatorFfromIndicatorF3D<T> stlF(stlReader, superGeometry);
  fields::set<descriptors::POROSITY>(sLattice, bulkIndicator, 1.);
  fields::set<descriptors::POROSITY>(sLattice, stlF, 0.);

  sLattice.initialize();
}

void setTemporalValues(MyCase& myCase, std::size_t iT) {
  using T = MyCase::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& converter = sLattice.getUnitConverter();
  auto& superGeometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  const T physStartT = parameters.get<parameters::PHYS_START_T>();
  const T physBoundaryT = parameters.get<parameters::PHYS_BOUNDARY_VALUE_UPDATE_T>();

  const std::size_t itStart = converter.getLatticeTime(physStartT);
  const std::size_t itUpdate = converter.getLatticeTime(physBoundaryT);
  if (iT <= itStart && iT % itUpdate == 0) {
    PolynomialStartScale<T,std::size_t> StartScale( itStart, T( 1 ) );
    std::size_t iTvec[1] = {iT};
    T frac[1] = {};
    StartScale(frac, iTvec);

    const T uF = 2.0 * frac[0] * converter.getCharPhysVelocity();
    CirclePoiseuille3D<T> poiseuilleU(superGeometry, 3, uF, T());
    momenta::setVelocity(sLattice, superGeometry.getMaterialIndicator({3}), poiseuilleU);
    sLattice.template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }
}

void getResults(MyCase& myCase, std::size_t iT, util::Timer<MyCase::value_t> timer) {
  using T = MyCase::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& parameters = myCase.getParameters();
  const T physMaxT = parameters.get<parameters::PHYS_MAX_T>();
  auto& converter = sLattice.getUnitConverter();
  const std::size_t iTlog = converter.getLatticeTime(physMaxT/5);

  // Get statistics
  if (iT%iTlog == 0) {
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    timer.print(iT);
  }
}

void simulate(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& sLattice = myCase.getLattice(NavierStokes{});
  auto& superGeometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  auto& converter = sLattice.getUnitConverter();
  const T physMaxT = parameters.get<parameters::PHYS_MAX_T>();

  // === 4th Step: Main Loop with Timer ===
  const std::size_t iTmax = converter.getLatticeTime(physMaxT);
  util::Timer<T> timer(iTmax, superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    if (!myCase.hasName(Adjoint{})) {
      setTemporalValues(myCase, iT);
    }
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    if (myCase.hasName(Adjoint{})) {
      myCase.getLattice(NavierStokes{}).stripeOffDensityOffset(
        myCase.getLattice(NavierStokes{}).getStatistics().getAverageRho() - T{1});
    }
    // === 7th Step: Computation and Output of the Results ===
    getResults(myCase, iT, timer);
  }
  timer.stop();
  timer.printSummary();
  myCase.getLattice(NavierStokes{}).setProcessingContext(ProcessingContext::Evaluation);
}

void prepareAdjointLattice(MyOptiCase& optiCase) {
  auto& adjointLattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& geometry = optiCase.getCase(Adjoint{}).getGeometry();

  // Set adjoint unit converter as the controlled case
  adjointLattice.setUnitConverter(controlledLattice.getUnitConverter());

  // Define dual physics
  dynamics::set<DualPorousBGKDynamics>(adjointLattice, geometry.getMaterialIndicator({1,6}));
  boundary::set<boundary::BounceBack>(adjointLattice, geometry.getMaterialIndicator({2,3,4}));
  adjointLattice.template setParameter<descriptors::OMEGA>(adjointLattice.getUnitConverter().getLatticeRelaxationFrequency());

  adjointLattice.revertStreaming();
}

void setAdjointInitialValues(MyOptiCase& optiCase) {
  auto& adjointLattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& converter = adjointLattice.getUnitConverter();
  auto& geometry = optiCase.getCase(Adjoint{}).getGeometry();
  auto bulkIndicator = geometry.getMaterialIndicator({1,2,3,4,6});
  auto objectiveDomain = geometry.getMaterialIndicator({6});

  adjointLattice.initialize();

  // Compute source term for the dual simulation
  auto objectiveDerivativeO = makeWriteFunctorO<functors::DerivativeF<ObjectiveF,descriptors::POPULATION,BulkDynamics>,
                                                opti::DJDF>(controlledLattice);
  objectiveDerivativeO->restrictTo(objectiveDomain);
  objectiveDerivativeO->template setParameter<parameters::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
  objectiveDerivativeO->template setParameter<parameters::NORMALIZE>(computeL2Norm<ObjectiveF::Reference>(controlledLattice, objectiveDomain, converter.getPhysDeltaX()));
  objectiveDerivativeO->template setParameter<parameters::FACTOR>(util::pow(converter.getPhysDeltaX(),3));
  objectiveDerivativeO->apply();

  // Write primal population in dedicated fields
  writeFunctorTo<functors::PopulationF,opti::F>(controlledLattice, bulkIndicator);

  // Provide fields required by the dual collision operator
  copyFields<ControlledField,ControlledField>(controlledLattice, adjointLattice);
  copyFields<opti::F,opti::F>(controlledLattice, adjointLattice);
  copyFields<opti::DJDF,opti::DJDF>(controlledLattice, adjointLattice);
}

void setInitialControl(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& lattice = controlledCase.getLattice(NavierStokes{});
  auto& geometry = controlledCase.getGeometry();
  auto& control = optiCase.getController();
  auto& converter = lattice.getUnitConverter();
  auto& parameters = controlledCase.getParameters();
  auto designDomain = geometry.getMaterialIndicator({6});

  // Intialize controlled field in superLattice
  T porosity = projection::permeabilityToPorosity(parameters.get<parameters::INITIAL_CONTROL_SCALAR>(), converter);
  AnalyticalConst3D<T,T> controls(porosity);
  fields::set<ControlledField>(lattice, geometry.getMaterialIndicator({1,2,3,4,6}), 1.0);
  fields::set<ControlledField>(lattice, designDomain, controls);
  control.template setProjection<projection::Sigmoid<T>>();
  control.template set<ControlledField>(geometry, 6, lattice);
}

void applyControl(MyOptiCase& optiCase) {
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& lattice = controlledCase.getLattice(NavierStokes{});
  auto& geometry = controlledCase.getGeometry();

  fields::set<ControlledField>(lattice, geometry.getMaterialIndicator({1,2,3,4,6}), 1.0);
  optiCase.getController().template setUpdatedControlsOnField<ControlledField>(lattice);
  lattice.initialize();
}

MyCase::value_t objectiveF(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& converter = controlledLattice.getUnitConverter();
  auto& geometry = optiCase.getCase(Controlled{}).getGeometry();
  auto objectiveDomain = geometry.getMaterialIndicator({6});

  // Evaluate functor for objective computation
  auto objectiveO = makeWriteFunctorO<ObjectiveF,opti::J>(controlledLattice);
  objectiveO->restrictTo(objectiveDomain);

  // Write velocity field in the reference field
  AnalyticalDataVolumeF<T> velocityMRI("GUI_data_231113_7_medium_quality_vel.vti", 3);
  fields::set<ObjectiveF::Reference>(controlledLattice, objectiveDomain, velocityMRI);
  controlledLattice.template setProcessingContext<Array<ObjectiveF::Reference>>(ProcessingContext::Simulation);

  // Get solution from the reference simulation for the inverse problem
  objectiveO->template setParameter<parameters::CONVERSION_VELOCITY>(converter.getConversionFactorVelocity());
  objectiveO->template setParameter<parameters::NORMALIZE>(computeL2Norm<ObjectiveF::Reference>(controlledLattice, objectiveDomain, converter.getPhysDeltaX()));
  objectiveO->apply();
  return integrateField<opti::J>(controlledLattice, objectiveDomain, converter.getPhysDeltaX())[0];
}

std::vector<MyCase::value_t> derivativeF(MyOptiCase& optiCase) {
  auto& adjointLattice = optiCase.getCase(Adjoint{}).getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& converter = controlledLattice.getUnitConverter();

  // Evaluate optimality condition
  auto optimalityO = makeWriteFunctorO<functors::OptimalityF<BulkDynamics,ControlledField>,
                                       opti::SENSITIVITY<ControlledField>>(adjointLattice);
  optiCase.getController().template setUpdatedProjectionDerivativesOnField<opti::DPROJECTIONDALPHA<ControlledField>>(adjointLattice);
  adjointLattice.template setProcessingContext<Array<opti::DPROJECTIONDALPHA<ControlledField>>>(ProcessingContext::Simulation);
  // Compute jacobian of collision operator regarding control variable
  auto dCDalphaO = makeWriteFunctorO<functors::DerivativeF<functors::CollisionF<BulkDynamics>,ControlledField,BulkDynamics>,
                                     opti::DCDALPHA<ControlledField>>(controlledLattice);
  dCDalphaO->template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  dCDalphaO->template setParameter<parameters::FACTOR>(1.0);
  dCDalphaO->apply();
  // Jacobian is computed on primal lattice as jacobian is evaluated for primal populations
  copyFields<opti::DCDALPHA<ControlledField>,opti::DCDALPHA<ControlledField>>(controlledLattice, adjointLattice);
  optimalityO->apply();

  // Return serial vector containing total derivatives of objective regarding controls
  adjointLattice.setProcessingContext(ProcessingContext::Evaluation);
  return getSerializedFromField<opti::SENSITIVITY<ControlledField>>(adjointLattice, optiCase.getController().getDesignDomain<MyCase::descriptor_t>());
}

void getOptiResults(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  auto& referenceLattice = optiCase.getCase(Reference{}).getLattice(NavierStokes{});
  auto& controlledLattice = optiCase.getCase(Controlled{}).getLattice(NavierStokes{});
  auto& converter = controlledLattice.getUnitConverter();
  std::size_t iT = optiCase.getOptimizationStep();

  SuperVTMwriter3D<T> vtmWriter("foamId3d");
  SuperLatticePhysPressure3D pressureF(controlledLattice, converter);
  SuperLatticePhysVelocity3D velocityF(controlledLattice, converter);
  SuperLatticeField3D<T,MyCase::descriptor_t,ControlledField> porosityF(controlledLattice);
  AnalyticalDataVolumeF<T> velocityMRI("GUI_data_231113_7_medium_quality_vel.vti", 3);
  SuperLatticeFfromAnalyticalF3D<T,MyCase::descriptor_t> velocityMRIF(velocityMRI, controlledLattice);
  AnalyticalDataVolumeF<T> porosityMRI("GUI_data_231113_7_medium_quality_geom.vti", 1);
  SuperLatticeFfromAnalyticalF3D<T,MyCase::descriptor_t> porosityMRIF(porosityMRI, controlledLattice);
  SuperLatticePhysVelocity3D velocityRefF(referenceLattice, converter);
  SuperLatticeField3D<T,MyCase::descriptor_t,ControlledField> porosityRefF(referenceLattice);

  pressureF.getName() = "physPressure (CFD-MRI)";
  vtmWriter.addFunctor(pressureF);
  velocityF.getName() = "physVelocity (CFD-MRI)";
  vtmWriter.addFunctor(velocityF);
  porosityF.getName() = "porosity (CFD-MRI)";
  vtmWriter.addFunctor(porosityF);
  velocityMRIF.getName() = "physVelocity (MRV)";
  vtmWriter.addFunctor(velocityMRIF);
  porosityMRIF.getName() = "porosity (MRI)";
  vtmWriter.addFunctor(porosityMRIF);
  velocityRefF.getName() = "physVelocity (CT-CFD)";
  vtmWriter.addFunctor(velocityRefF);
  porosityRefF.getName() = "porosity (CT)";
  vtmWriter.addFunctor(porosityRefF);

  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  // Writes the VTK files
  if (iT > 0) {
    referenceLattice.setProcessingContext(ProcessingContext::Evaluation);
    controlledLattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }
}

// Evaluate how much a distribution has two distinct peaks (doi: 10.22237/jmasm/1177992120)
MyCase::value_t computeBimodalityIndex(auto& porosityF, auto&& indicatorF) {
  using T = MyCase::value_t;
  auto& loadBalancer = porosityF.getSuperStructure().getLoadBalancer();
  auto& decomposition = porosityF.getSuperStructure().getCuboidDecomposition();

  // compute first the average porosity
  int numCells = 0;
  T sumPorosity = 0.0;
  for (int iC=0; iC < decomposition.size(); ++iC) {
    if (loadBalancer.isLocal(iC)) {
      const int loc = loadBalancer.loc(iC);
      auto& blockStructure = indicatorF->getBlockIndicatorF(loc).getBlockGeometry();
      blockStructure.forCoreSpatialLocations([&](auto iX, auto iY, auto iZ) {
        const int input[3] = {iX, iY, iZ};
        T curPor[4];
        if (indicatorF->getBlockIndicatorF(loc)(iX, iY, iZ)) {
          porosityF.getBlockF(loc)(curPor, input);
          sumPorosity += curPor[0];
          ++numCells;
        }
      });
    }
  }

#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(numCells, MPI_SUM);
  singleton::mpi().reduceAndBcast(sumPorosity, MPI_SUM);
#endif
  const T averagePorosity = sumPorosity / (numCells * 1.);

  // compute the moments
  T secondMoment = 0.0;
  T thirdMoment = 0.0;
  T fourthMoment = 0.0;
  for (int iC=0; iC < decomposition.size(); ++iC) {
    if (loadBalancer.isLocal(iC)) {
      const int loc = loadBalancer.loc(iC);
      auto& blockStructure = indicatorF->getBlockIndicatorF(loc).getBlockGeometry();
      blockStructure.forCoreSpatialLocations([&](auto iX, auto iY, auto iZ) {
        const int input[3] = {iX, iY, iZ};
        T curPor[4];
        if (indicatorF->getBlockIndicatorF(loc)(iX, iY, iZ)) {
          porosityF.getBlockF(loc)(curPor, input);
          const T arg = curPor[0] - averagePorosity;
          secondMoment += util::pow(arg, 2.0);
          thirdMoment += util::pow(arg, 3.0);
          fourthMoment += util::pow(arg, 4.0);
        }
      });
    }
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(secondMoment, MPI_SUM);
  singleton::mpi().reduceAndBcast(thirdMoment, MPI_SUM);
  singleton::mpi().reduceAndBcast(fourthMoment, MPI_SUM);
#endif
  // compute standard deviation
  const T sd = util::sqrt((1.0 / (numCells - 1.0)) * secondMoment);
  // compute skewness
  const T tmp = numCells / ((numCells - 1.0) * (numCells - 2.0));
  const T m3 = tmp * (1.0 / util::pow(sd, 3.0)) * thirdMoment;
  // compute kurtosis minus summand 3(n-1)
  const T m4 = tmp * ((numCells + 1.0) / (numCells - 3.0)) * (1.0 / util::pow(sd, 4.0)) * fourthMoment;
  // return bimodality index
  return (util::pow(m3, 2.0) + 1.0) / m4;
}

// Correlate two distributions (https://www.jstor.org/stable/115794)
std::vector<MyCase::value_t> computePearsonCorrelationCoefficient(auto& pairs, auto&& indicatorF) {
  using T = MyCase::value_t;
  auto& referenceF = *(pairs[0].first);
  auto& loadBalancer = referenceF.getSuperStructure().getLoadBalancer();
  auto& decomposition = referenceF.getSuperStructure().getCuboidDecomposition();

  // Compute average porosity first
  int numCells = 0;
  std::vector<std::vector<T>> sumPorosities(pairs.size(), {2, T()});
  for (int iC=0; iC < decomposition.size(); ++iC) {
    if (loadBalancer.isLocal(iC)) {
      const int loc = loadBalancer.loc(iC);
      auto& blockStructure = indicatorF->getBlockIndicatorF(loc).getBlockGeometry();
      blockStructure.forCoreSpatialLocations([&](auto iX, auto iY, auto iZ) {
        const int input[3] = {iX, iY, iZ};
        T curPor1[4]; T curPor2[4];
        if (indicatorF->getBlockIndicatorF(loc)(iX, iY, iZ)) {
          int i = 0;
          for (auto [F1, F2] : pairs) {
            F1->getBlockF(loc)(curPor1, input);
            sumPorosities[i][0] += curPor1[0];
            F2->getBlockF(loc)(curPor2, input);
            sumPorosities[i][1] += curPor2[0];
            ++i;
          }
          ++numCells;
        }
      });
    }
  }

#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(numCells, MPI_SUM);
  for (auto& e : sumPorosities) {
    singleton::mpi().reduceAndBcast(e[0], MPI_SUM);
    singleton::mpi().reduceAndBcast(e[1], MPI_SUM);
  }
#endif
  std::vector<std::vector<T>> averagePorosities(pairs.size(), {2, T()});
  for (std::size_t i=0; i<sumPorosities.size(); ++i) {
    averagePorosities[i][0] = sumPorosities[i][0] / (numCells * 1.);
    averagePorosities[i][1] = sumPorosities[i][1] / (numCells * 1.);
  }

  // Compute the moments and the inner product
  std::vector<std::vector<T>> secondMoments(pairs.size(), {2, T()});
  std::vector<T> covariances(pairs.size(), T());
  for (int iC=0; iC < decomposition.size(); ++iC) {
    if (loadBalancer.isLocal(iC)) {
      const int loc = loadBalancer.loc(iC);
      auto& blockStructure = indicatorF->getBlockIndicatorF(loc).getBlockGeometry();
      blockStructure.forCoreSpatialLocations([&](auto iX, auto iY, auto iZ) {
        const int input[3] = {iX, iY, iZ};
        T curPor1[4]; T curPor2[4];
        if (indicatorF->getBlockIndicatorF(loc)(iX, iY, iZ)) {
          std::vector<std::vector<T>> arg(pairs.size(), {2, T()}); int i=0;
          for (auto [F1, F2] : pairs) {
            F1->getBlockF(loc)(curPor1, input);
            F2->getBlockF(loc)(curPor2, input);
            arg[i] = {curPor1[0] - averagePorosities[i][0], curPor2[0] - averagePorosities[i][1]};
            covariances[i] += arg[i][0]*arg[i][1];
            secondMoments[i][0] += util::pow(arg[i][0], 2.0);
            secondMoments[i][1] += util::pow(arg[i][1], 2.0);
            ++i;
          }
        }
      });
    }
  }
#ifdef PARALLEL_MODE_MPI
  for (std::size_t i=0; i<pairs.size(); ++i) {
    singleton::mpi().reduceAndBcast(secondMoments[i][0], MPI_SUM);
    singleton::mpi().reduceAndBcast(secondMoments[i][1], MPI_SUM);
    singleton::mpi().reduceAndBcast(covariances[i], MPI_SUM);
  }
#endif
  // compute standard deviation
  const T factor = 1. / (numCells - 1.);
  std::vector<std::vector<T>> standardDeviations(pairs.size(), {2, T()});
  {int i=0;
  for (auto& e : standardDeviations) {
    e = {util::sqrt(factor * secondMoments[i][0]), util::sqrt(factor * secondMoments[i][1])};
    ++i;
  }}

  // Compute prarson index
  std::vector<T> PCC(pairs.size(), T());
  for (std::size_t i=0; i<PCC.size(); ++i) {
    PCC[i] = (factor * covariances[i]) / (standardDeviations[i][0] * standardDeviations[i][1]);
  }

  return PCC;
}

// Computes the similarity index according to binary classification (doi: 10.1016/j.patrec.2005.10.010)
// thresholds are used to define all porosity values below that value to consider being a solid node
std::vector<MyCase::value_t> computeJaccardIndex(auto& pairs, auto&& indicatorF, auto thresholds) {
  using T = MyCase::value_t;
  auto& referenceF = *(pairs[0].first);
  auto& loadBalancer = referenceF.getSuperStructure().getLoadBalancer();
  auto& decomposition = referenceF.getSuperStructure().getCuboidDecomposition();

  // All possible cases are:
  // tP (truePositive): Identified as solid is actual solid node
  // tN (trueNegative): Identified as fluid is actual fluid node
  // fP (falsePositive): Identified as solid is actual fluid node
  // fN (falseNegative): Identified as fluid is actual solid node
  std::vector<int> tP(pairs.size(), 0);
  std::vector<int> tN(pairs.size(), 0);
  std::vector<int> fP(pairs.size(), 0);
  std::vector<int> fN(pairs.size(), 0);
  int numCells = 0;
  for (int iC=0; iC < decomposition.size(); ++iC) {
    if (loadBalancer.isLocal(iC)) {
      const int loc = loadBalancer.loc(iC);
      auto& blockStructure = indicatorF->getBlockIndicatorF(loc).getBlockGeometry();
      blockStructure.forCoreSpatialLocations([&](auto iX, auto iY, auto iZ) {
        const int input[3] = {iX, iY, iZ};
        T curPor1[4]; T curPor2[4];
        if (indicatorF->getBlockIndicatorF(loc)(iX, iY, iZ)) {
          int i = 0;
          for (auto [F1, F2] : pairs) {
            F1->getBlockF(loc)(curPor1, input);
            F2->getBlockF(loc)(curPor2, input);
            const int cP = (curPor1[0] < thresholds[i][0]) ? 0 : 1;
            const int rP = (curPor2[0] < thresholds[i][1]) ? 0 : 1;
            if ((cP == 0) && (rP == 0)) {
              ++tP[i];
            }
            else if ((cP == 1) && (rP == 1)) {
              ++tN[i];
            }
            else if ((cP == 0) && (rP == 1)) {
              ++fP[i];
            }
            else if ((cP == 1) && (rP == 0)) {
              ++fN[i];
            }
            ++i;
          }
          ++numCells;
        }
      });
    }
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(numCells, MPI_SUM);
  for (std::size_t i=0; i<pairs.size(); ++i) {
    singleton::mpi().reduceAndBcast(tP[i], MPI_SUM);
    singleton::mpi().reduceAndBcast(tN[i], MPI_SUM);
    singleton::mpi().reduceAndBcast(fP[i], MPI_SUM);
    singleton::mpi().reduceAndBcast(fN[i], MPI_SUM);
  }
#endif

  std::vector<T> J(pairs.size(), T());
  for (std::size_t i=0; i<pairs.size(); ++i) {
    J[i] = (tP[i] * 1.) / ((tP[i] + fP[i] + fN[i]) * 1.);
  }
  return J;
}

void computeQualityMeasures(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& referenceCase = optiCase.get(Reference{});
  auto& controlledCase = optiCase.get(Controlled{});
  auto& referenceLattice = referenceCase.getLattice(NavierStokes{});
  auto& controlledLattice = controlledCase.getLattice(NavierStokes{});
  auto& geometry = controlledCase.getGeometry();

  const std::size_t iT = optiCase.getOptimizationStep();
  auto porosity_cfd_mri = std::make_shared<SuperLatticePorosity3D<T,DESCRIPTOR>>(controlledLattice);
  auto porosity_ct = std::make_shared<SuperLatticePorosity3D<T,DESCRIPTOR>>(referenceLattice);
  auto porosityMRI = std::make_shared<AnalyticalDataVolumeF<T>>("GUI_data_231113_7_medium_quality_geom.vti", 1);
  auto porosity_mri = std::make_shared<SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR>>(*porosityMRI, controlledLattice);

  std::vector<std::pair<std::shared_ptr<SuperLatticeF3D<T,DESCRIPTOR>>,
                        std::shared_ptr<SuperLatticeF3D<T,DESCRIPTOR>>>> pairs;
  pairs.push_back({porosity_cfd_mri, porosity_ct});
  pairs.push_back({porosity_mri, porosity_cfd_mri});
  pairs.push_back({porosity_ct, porosity_mri});

  std::vector<std::vector<T>> porosityThresholdsForJaccard{{0.97, 0.50},
                                                           {0.50, 0.97},
                                                           {0.50, 0.50}};

  const T B = computeBimodalityIndex(*porosity_cfd_mri, geometry.getMaterialIndicator({6}));
  const std::vector<T> PCC = computePearsonCorrelationCoefficient(pairs, geometry.getMaterialIndicator({6}));
  const std::vector<T> J = computeJaccardIndex(pairs, geometry.getMaterialIndicator({6}), porosityThresholdsForJaccard);

  OstreamManager clout(std::cout, "QualityMeasures");
  clout << "========== QUALITY MEASURES ==========" << std::endl;
  clout << "Optimization step: " << iT << std::endl;
  clout << "Bimodality Index                (CFD-MRI          ): " << B << std::endl;
  clout << "Pearson Correlation Coefficient (CFD-MRI / CT     ): " << PCC[0] << std::endl;
  clout << "Pearson Correlation Coefficient (MRI     / CFD-MRI): " << PCC[1] << std::endl;
  clout << "Pearson Correlation Coefficient (CT      / MRI    ): " << PCC[2] << std::endl;
  clout << "Jaccard Index                   (CFD-MRI / CT     ): " << J[0] << std::endl;
  clout << "Jaccard Index                   (MRI     / CFD-MRI): " << J[1] << std::endl;
  clout << "Jaccard Index                   (CT      / MRI    ): " << J[2] << std::endl;
  clout << "======================================" << std::endl;
}

MyCase::value_t computeObjective(MyOptiCase& optiCase) {
  auto& controlledCase = optiCase.getCase(Controlled{});
  // Reset prior simulation lattice
  controlledCase.resetLattices();

  // Prepare new lattices
  prepareLattice(controlledCase);

  // Set updated controls in the simulation
  applyControl(optiCase);

  // Execute simulation
  simulate(controlledCase);

  // Optimization related IO
  getOptiResults(optiCase);

  // Compute quality measures
  computeQualityMeasures(optiCase);

  // Compute Objective value from simulation results
  return objectiveF(optiCase);
}

std::vector<MyCase::value_t> computeDerivative(MyOptiCase& optiCase) {
  auto& adjointCase = optiCase.getCase(Adjoint{});
  // Reset prior simulation lattice
  adjointCase.resetLattices();

  // Prepare new lattices
  prepareAdjointLattice(optiCase);
  setAdjointInitialValues(optiCase);

  // Execute simulation
  simulate(adjointCase);

  // Compute Derivative value from simulation results
  return derivativeF(optiCase);
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);

  MyCase::ParametersD myCaseParameters;
  {
    using namespace parameters;
    myCaseParameters.set<DOMAIN_CENTER1              >(     {0.017, 0.0155, -0.02});
    myCaseParameters.set<DOMAIN_CENTER2              >(     {0.017, 0.0155, 0.130});
    myCaseParameters.set<DOMAIN_RADIUS               >(                     0.0125);
    myCaseParameters.set<PHYS_CHAR_VELOCITY          >(                  0.0071641);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY         >(                       1e-6);
    myCaseParameters.set<PHYS_CHAR_DENSITY           >(                      0.997);
    myCaseParameters.set<PHYS_MAX_T                  >(                       50.0);
    myCaseParameters.set<RESOLUTION                  >(                         51);
    myCaseParameters.set<LATTICE_RELAXATION_TIME     >(                       0.51);
    myCaseParameters.set<PHYS_START_T                >(                        4.0);
    myCaseParameters.set<PHYS_BOUNDARY_VALUE_UPDATE_T>(                       0.05);
    myCaseParameters.set<INITIAL_CONTROL_SCALAR      >(                      1.e-7);
    myCaseParameters.set<DESIGN_DOMAIN_CENTER1       >({0.017, 0.0155, 0.00451379});
    myCaseParameters.set<DESIGN_DOMAIN_CENTER2       >( {0.017, 0.0155, 0.0285546});
    myCaseParameters.set<DESIGN_DOMAIN_RADIUS        >(                     0.0125);
  }

  /// Step 3: Create Mesh
  Mesh mesh = createMesh(myCaseParameters);

  /// Step B: Create Cases and prepare
  // Prepare controlled case
  MyCase myCase(myCaseParameters, mesh);
  prepareGeometry(myCase);
  prepareLattice(myCase);

  // Create adjoint case for gradient computation
  MyCase adjointCase(myCaseParameters, mesh);
  prepareGeometry(adjointCase);

  // Compute solution for the objective functional
  MyCase referenceCase(myCaseParameters, mesh);
  prepareGeometry(referenceCase);
  prepareLattice(referenceCase);
  setInitialValues(referenceCase);
  simulate(referenceCase);

  /// Step C: Create OptiCase and set cases
  MyOptiCase optiCase;
  optiCase.setCase<Controlled>(myCase);
  optiCase.setCase<Adjoint>(adjointCase);
  optiCase.setCase<Reference>(referenceCase);

  /// Step D: Set initial control
  setInitialControl(optiCase);

  /// Step E: Define objective routine
  optiCase.setObjective(computeObjective);

  /// Step F: Define derivative routine
  optiCase.setDerivative(computeDerivative);

  OptimizerParametersD<MyCase::value_t> optiParameters;
  {
    using namespace parameters;
    optiParameters.set<optimizer::EPS>(1.e-75);
    optiParameters.set<linesearch::INITIAL_STEP_SIZE>(500.);
    optiParameters.set<linesearch::STEP_CONDITION>("Wolfe");
  }

  OptimizerLBFGS<MyCase::value_t,std::vector<MyCase::value_t>> optimizer(
    optiCase.getController().size(), optiParameters);
  optimizer.optimize(optiCase);
}
