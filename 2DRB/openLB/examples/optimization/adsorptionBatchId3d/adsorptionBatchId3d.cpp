/*  Lattice Boltzmann sample, written in C++, using the OpenLBlibrary
 *
 *  Copyright (C) 2026 Simon Großmann, Shota Ito
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

// This app is based on the work in doi: 10.1016/j.ces.2025.121467.

// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::opti;

namespace olb::parameters {
  struct BATCH_REACTOR_LENGTH            : public descriptors::FIELD_BASE<1> {};
  struct PARTICLE_RADIUS                 : public descriptors::FIELD_BASE<1> {};
  struct PARTICLE_DENSITY                : public descriptors::FIELD_BASE<1> {};
  struct INITIAL_CONCENTRATION           : public descriptors::FIELD_BASE<1> {};
  struct HENRY_ISOTHERM_CONSTANT         : public descriptors::FIELD_BASE<1> {};
  struct SURFACE_DIFFUSION_COEFFICIENT   : public descriptors::FIELD_BASE<1> {};
  struct FILM_MASS_TRANSFER_COEFFICIENT  : public descriptors::FIELD_BASE<1> {};

  struct SURFACE_TRANSFER_CONSTANT       : public descriptors::FIELD_BASE<1> {};
  struct DISTRIBUTION_FACTOR             : public descriptors::FIELD_BASE<1> {};

  struct L2Distance                      : public descriptors::FIELD_BASE<1> {};
}

using MyCase = Case<
  Concentration0, Lattice<double, descriptors::D3Q7<descriptors::SOURCE>>,
  Concentration1, Lattice<double, descriptors::D3Q7<descriptors::SOURCE>>,
  Concentration2, Lattice<double, descriptors::D3Q7<descriptors::SOURCE>>
>;
using MyADfCase = MyCase::exchange_value_t<util::ADf<double, 1>>;
using MyOptiCase = OptiCaseADf<
  Controlled, MyCase,
  Derivatives, MyADfCase
>;

template <typename T>
struct AnalyticalSolution : public AnalyticalF3D<T,T> {
  T t;
  T D_b;

  AnalyticalSolution(T time, T D_b, T ks) : AnalyticalF3D<T,T>(1),
                                            t(time * ks),
                                            D_b(D_b) {}
  bool operator() (T output[], const T input[]) override {
    output[0] = 1. / (D_b+1.) + D_b / (D_b+1.) * util::exp(-(D_b+1)*t);
    return true;
  }
};

template <typename PARAMETERS>
auto createMesh(PARAMETERS& parameters) {
  using T = PARAMETERS::value_t;
  const Vector<T,3> origin = parameters.template get<parameters::DOMAIN_EXTENT>() * (-0.5);
  const Vector<T,3> extend = parameters.template get<parameters::DOMAIN_EXTENT>();
  IndicatorCuboid3D<T> cuboid(extend, origin);

  Mesh<T,3> mesh(cuboid, parameters.template get<parameters::PHYS_DELTA_X>(), singleton::mpi().getSize());

  mesh.setOverlap(parameters.template get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({true, true, true});

  return mesh;
 }

// prepare numerical domain
template <typename CASE>
void prepareGeometry(CASE& myCase){
  auto& superGeometry   = myCase.getGeometry();
  superGeometry.rename(0, 1);
  superGeometry.print();
}

template <typename CASE>
void prepareLattice(CASE& myCase){
  using T = CASE::value_t;
  auto& superGeometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  auto& ParticleLattice = myCase.getLattice(Concentration0{});
  auto& LoadingLattice = myCase.getLattice(Concentration1{});
  auto& SoluteLattice = myCase.getLattice(Concentration2{});

  ParticleLattice.template setUnitConverter<AdsorptionConverterFromSchmidtNumberAndRelaxation<T,MyCase::descriptor_t_of<Concentration0>>>(
    (T)   parameters.template get<parameters::SCHMIDT>(),
    (T)   parameters.template get<parameters::REYNOLDS>(),
    (T)   parameters.template get<parameters::LATTICE_RELAXATION_TIME>(),
    (T)   parameters.template get<parameters::BATCH_REACTOR_LENGTH>(),
    (T)   parameters.template get<parameters::BATCH_REACTOR_LENGTH>(),
    (T)   parameters.template get<parameters::RESOLUTION>(),
    (T)   0.02, // physCharVelocity (dummy)
    (T)   parameters.template get<parameters::PHYS_CHAR_DENSITY>(),
    (T)   parameters.template get<parameters::PARTICLE_DENSITY>()
  );
  ParticleLattice.getUnitConverter().print();
  const auto& converter = ParticleLattice.getUnitConverter();
  LoadingLattice.setUnitConverter(converter);
  SoluteLattice.setUnitConverter(converter);

  // define dynamics and boundary conditions
  const auto bulkIndicator = superGeometry.getMaterialIndicator({1});
  dynamics::set<SourcedAdvectionDiffusionBGKdynamics>(ParticleLattice, bulkIndicator);
  dynamics::set<SourcedAdvectionDiffusionBGKdynamics>(LoadingLattice, bulkIndicator);
  dynamics::set<SourcedAdvectionDiffusionBGKdynamics>(SoluteLattice, bulkIndicator);

  auto& coupling = myCase.setCouplingOperator(
    "AdsorptionReaction",
    AdsorptionBatchCoupling3D<AdsorptionReaction<Isotherm::LinearIsotherm>>{},
    Concentration0{}, ParticleLattice,
    Concentration1{}, SoluteLattice,
    Concentration2{}, LoadingLattice);
  coupling.template setParameter<AdsorptionReaction<Isotherm::LinearIsotherm>::K_F>(
    parameters.template get<parameters::FILM_MASS_TRANSFER_COEFFICIENT>());
  coupling.template setParameter<AdsorptionReaction<Isotherm::LinearIsotherm>::D_S>(
    parameters.template get<parameters::SURFACE_DIFFUSION_COEFFICIENT>());
  coupling.template setParameter<AdsorptionReaction<Isotherm::LinearIsotherm>::C_0>(
    parameters.template get<parameters::INITIAL_CONCENTRATION>());
  coupling.template setParameter<AdsorptionReaction<Isotherm::LinearIsotherm>::R_P>(
    parameters.template get<parameters::PARTICLE_RADIUS>());

  Isotherm::LinearIsotherm::setParameters<T>(
    parameters.template get<parameters::HENRY_ISOTHERM_CONSTANT>(), 0, coupling);
  AdsorptionReaction<Isotherm::LinearIsotherm>::computeParameters<T>(coupling, converter);

  OstreamManager clout(std::cout, "AdsorptionReaction");
  AdsorptionReaction<Isotherm::LinearIsotherm>::print<T>(clout, coupling, converter);
}

template <typename CASE>
void setInitialValues(CASE& myCase) {
  auto& superGeometry = myCase.getGeometry();
  auto& ParticleLattice = myCase.getLattice(Concentration0{});
  auto& LoadingLattice = myCase.getLattice(Concentration1{});
  auto& SoluteLattice = myCase.getLattice(Concentration2{});
  const auto& converter = ParticleLattice.getUnitConverter();
  auto& parameters = myCase.getParameters();
  auto factor = converter.getConversionFactorDensity();

  momenta::setDensity(ParticleLattice, superGeometry.getMaterialIndicator({1}), 1.* factor);
  momenta::setDensity(LoadingLattice, superGeometry.getMaterialIndicator({1}), 0.);
  momenta::setDensity(SoluteLattice, superGeometry.getMaterialIndicator({1}),
    parameters.template get<parameters::INITIAL_CONCENTRATION>() * factor);

  ParticleLattice.template setParameter<descriptors::OMEGA>(
    converter.getLatticeRelaxationFrequency());
  LoadingLattice.template setParameter<descriptors::OMEGA>(
    converter.getLatticeRelaxationFrequency());
  SoluteLattice.template setParameter<descriptors::OMEGA>(
    converter.getLatticeRelaxationFrequency());

  ParticleLattice.initialize();
  LoadingLattice.initialize();
  SoluteLattice.initialize();
}

template <typename CASE>
void objectiveContribution(CASE& myCase, auto& solute, auto& analytical) {
  OstreamManager clout(std::cout, "Objective");
  using T = CASE::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& SoluteLattice = myCase.getLattice(Concentration2{});

  SoluteLattice.setProcessingContext(ProcessingContext::Evaluation);
  SuperAbsoluteErrorL2Norm3D<T> norm(solute, analytical, geometry.getMaterialIndicator({1}));
  int input [3]; T output[1];
  norm(output, input);
  parameters.template set<parameters::L2Distance>(
    util::pow(parameters.template get<parameters::L2Distance>() + output[0], 2));
}

// compute results and generate output
template <typename CASE, typename TIMER>
void getResults(CASE& myCase, TIMER& timer, std::size_t iT)
{
  using T = CASE::value_t;
  auto& ParticleLattice = myCase.getLattice(Concentration0{});
  auto& LoadingLattice = myCase.getLattice(Concentration1{});
  auto& SoluteLattice = myCase.getLattice(Concentration2{});
  auto& converter = ParticleLattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  // generate vtm files for visualisation
  SuperVTMwriter3D<T> vtmWriter("adsorptionBatchId3d");
  SuperLatticeDensity3D<T,MyCase::descriptor_t_of<Concentration0>> particle(ParticleLattice);
  SuperLatticeDensity3D<T,MyCase::descriptor_t_of<Concentration1>> solute(SoluteLattice);
  SuperLatticeDensity3D<T,MyCase::descriptor_t_of<Concentration2>> loading(LoadingLattice);
  AnalyticalSolution<T> solution(converter.getPhysTime(iT),
                                 parameters.template get<parameters::DISTRIBUTION_FACTOR>(),
                                 parameters.template get<parameters::SURFACE_TRANSFER_CONSTANT>());
  SuperLatticeFfromAnalyticalF3D<T,MyCase::descriptor_t_of<Concentration0>> analytical(solution, ParticleLattice);

  particle.getName() = "Particle Concentration";
  loading.getName() = "Loading Concentration";
  solute.getName() = "Solute Concentration";
  analytical.getName() = "Analytical Solution";
  vtmWriter.addFunctor(particle);
  vtmWriter.addFunctor(loading);
  vtmWriter.addFunctor(solute);
  vtmWriter.addFunctor(analytical);

  const std::size_t vtkIter = converter.getLatticeTime(parameters.template get<parameters::PHYS_MAX_T>() / 100);
  const std::size_t statIter = converter.getLatticeTime(parameters.template get<parameters::PHYS_MAX_T>() / 5);

  if (iT == 0) {
    // write the geometry
    vtmWriter.createMasterFile();

    // write inital field
    vtmWriter.write(iT);
  }

  if (iT%statIter == 0) {
    // Timer console output
    timer.update(iT);
    timer.printStep();

    // Lattice statistics console output
    SoluteLattice.getStatistics().print(iT, converter.getPhysTime(iT));
  }
  // write the vtk files
  if (iT % vtkIter == 0 && iT > 0) {
    SoluteLattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iT);
  }

  objectiveContribution(myCase, solute, solution);
}

template <typename CASE>
void simulate(CASE& myCase){
  OstreamManager clout(std::cout,"simulate");
  using T = CASE::value_t;
  auto& parameters = myCase.getParameters();
  auto& ParticleLattice = myCase.getLattice(Concentration0{});
  auto& LoadingLattice = myCase.getLattice(Concentration1{});
  auto& SoluteLattice = myCase.getLattice(Concentration2{});
  auto& converter = ParticleLattice.getUnitConverter();
  auto& superGeometry = myCase.getGeometry();

  util::Timer<T> timer(converter.getLatticeTime(parameters.template get<parameters::PHYS_MAX_T>()), superGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT < converter.getLatticeTime(parameters.template get<parameters::PHYS_MAX_T>()); ++iT) {
    myCase.getOperator("AdsorptionReaction").apply();
    ParticleLattice.collideAndStream();
    LoadingLattice.collideAndStream();
    SoluteLattice.collideAndStream();
    getResults(myCase, timer, iT);
  }

  timer.stop();
  timer.printSummary();
}

void setInitialControl(MyOptiCase& optiCase) {
  // control[0] : Henry isotherm constant
  optiCase.getController().set({15});
}

template <typename CASE, typename CONTROLLER>
void applyControl(CASE& myCase, CONTROLLER& controller) {
  auto& parameters = myCase.getParameters();
  auto controls = controller.get();
  parameters.template set<parameters::HENRY_ISOTHERM_CONSTANT>(controls[0]);
}

template <typename NAME>
auto computeObjective(MyOptiCase& optiCase, NAME) {
  using CASE = MyOptiCase::case_t<NAME>;
  using T = CASE::value_t;
  auto& myCase = optiCase.getCase(NAME{});
  auto& parameters = myCase.getParameters();

  myCase.resetLattices();
  parameters.template set<parameters::L2Distance>(0.);

  applyControl(myCase, optiCase.getController<T>());
  prepareLattice(myCase);
  setInitialValues(myCase);
  simulate(myCase);
  return parameters.template get<parameters::L2Distance>();
}

int main(int argc, char* argv[]) {
    /// === 1st Step: Initialization ===
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION                    >(       20);
    myCaseParameters.set<BATCH_REACTOR_LENGTH          >(      0.1);
    myCaseParameters.set<REYNOLDS                      >(      17.);
    myCaseParameters.set<SCHMIDT                       >(      1.5);
    myCaseParameters.set<PARTICLE_RADIUS               >(   1.5e-4);
    myCaseParameters.set<PARTICLE_DENSITY              >(     0.94);
    myCaseParameters.set<INITIAL_CONCENTRATION         >(       1.);
    myCaseParameters.set<LATTICE_RELAXATION_TIME       >(    0.695);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY           >(  1.18e-4);
    myCaseParameters.set<PHYS_CHAR_DENSITY             >(       1.);
    myCaseParameters.set<HENRY_ISOTHERM_CONSTANT       >(      45.);
    myCaseParameters.set<SURFACE_DIFFUSION_COEFFICIENT >(    5e-11);
    myCaseParameters.set<FILM_MASS_TRANSFER_COEFFICIENT>(       0.);
    myCaseParameters.set<PHYS_DELTA_X                  >([](MyCase::ParametersD& p) {
      return p.get<BATCH_REACTOR_LENGTH>() /
             p.get<RESOLUTION>();
    });
    myCaseParameters.set<DOMAIN_EXTENT>([](MyCase::ParametersD& p) {
      return Vector{p.get<BATCH_REACTOR_LENGTH>(),
                    p.get<BATCH_REACTOR_LENGTH>(),
                    p.get<BATCH_REACTOR_LENGTH>()};
    });
    myCaseParameters.set<PHYS_MAX_T                    >(      15.);
    // Used for the analytical solution, computed once here based on
    // reference parameters for D_s and r_p.
    myCaseParameters.set<SURFACE_TRANSFER_CONSTANT>([](MyCase::ParametersD& p) {
      return 15.*
             p.get<SURFACE_DIFFUSION_COEFFICIENT>() /
             util::pow(p.get<PARTICLE_RADIUS>(), 2);
    });
    myCaseParameters.set<DISTRIBUTION_FACTOR>([](MyCase::ParametersD& p) {
      return p.get<PARTICLE_DENSITY>() *
             p.get<HENRY_ISOTHERM_CONSTANT>();
    });
  }
  auto myADfParameters = myCaseParameters.copyAs<MyADfCase::value_t>();

  OptimizerParametersD<MyCase::value_t> optiParameters;
  {
    using namespace parameters;
    optiParameters.set<optimizer::EPS>(1e-10);
    optiParameters.set<optimizer::MAX_OPTIMIZATION_STEPS>(100);
    optiParameters.set<linesearch::MAX_STEP_ATTEMPTS>(40);
    optiParameters.set<linesearch::STEP_CONDITION>("StrongWolfe");
    optiParameters.set<optimizer::VERBOSE>(true);
    optiParameters.set<optimizer::WITH_LOWER_BOUND>(true);
    optiParameters.set<optimizer::LOWER_BOUND>(1e-3);
    optiParameters.set<optimizer::CONTROL_EPS>(0.);
  }

  auto mesh = createMesh(myCaseParameters);
  auto aDfMesh = createMesh(myADfParameters);

  MyCase myCase(myCaseParameters, mesh);
  prepareGeometry(myCase);

  MyADfCase myADfCase(myADfParameters, aDfMesh);
  prepareGeometry(myADfCase);

  MyOptiCase optiCase;
  optiCase.setCase<Controlled>(myCase);
  optiCase.setCase<Derivatives>(myADfCase);

  setInitialControl(optiCase);

  optiCase.setObjective([](auto& optiCase, auto name) {
    return computeObjective(optiCase, name);
  });

  OptimizerLBFGS<double,std::vector<double>> optimizer(
    optiCase.getController().size(), optiParameters,
    {OptimizerLogType::objective, OptimizerLogType::control, OptimizerLogType::derivative}
  );

  optimizer.optimize(optiCase);
}
