/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Florian Raichle, Luiz Czelusniak, Fedor Bukreev
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

/**
 * adsorption3D
 * Simulating adsorption in a batch reactor using an Euler-Euler approach.
 * The model is based on the linear driving force model and uses advection diffusion reaction lattices for particles,
 * solute and particle loading.
 *
 * Different isotherms and mass transfer models can be used.
 * An analytical solution is implemented when using the linear isotherm and surface diffusion.
 */

#include <olb.h>

//#include "../isotherms.h"

using namespace olb;
using namespace olb::names;
using namespace olb::descriptors;

template <typename T>
class BatchSolution: public AnalyticalF3D<T, T> {
protected:
  T t;
  T D_b_;

public:
  BatchSolution(T time, T D_b, T ks, T k_f, T particleRadius) : AnalyticalF3D<T, T>(1), t(time*3.*k_f/particleRadius/D_b), D_b_(D_b) {}
  bool operator () (T output[], const T input[]) override {
    output[0] = T(1)/(D_b_+T(1)) + D_b_ / (D_b_+T(1)) * util::exp(-(D_b_+T(1))*t);
    return true;
  }
};

using MyCase = Case<
  NavierStokes,       Lattice<double, descriptors::D3Q19<>>,
  Concentration<0>,     Lattice<double, descriptors::D3Q7<VELOCITY,SOURCE>>,
  Concentration<1>,     Lattice<double, descriptors::D3Q7<VELOCITY,SOURCE>>,
  Concentration<2>,     Lattice<double, descriptors::D3Q7<VELOCITY,SOURCE>>
>;

namespace olb::parameters {

struct PARTICLE_RADIUS : public descriptors::FIELD_BASE<1> { };
struct PARTICLE_DENSITY : public descriptors::FIELD_BASE<1> { };
struct ISO_CONST_A : public descriptors::FIELD_BASE<1> { };
struct ISO_CONST_B : public descriptors::FIELD_BASE<1> { };
struct K_F : public descriptors::FIELD_BASE<1> { };
struct C_0 : public descriptors::FIELD_BASE<1> { };
struct D_S : public descriptors::FIELD_BASE<1> { };
struct D_B : public descriptors::FIELD_BASE<1> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const T physDeltaX = extend[0]/parameters.get<parameters::RESOLUTION>();

  Vector<T,3> origin(-extend[0] / 2, -extend[0] / 2, -extend[0] / 2);
  IndicatorCuboid3D<T> cuboid(extend, origin);

  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({true,true,true});
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  auto& geometry = myCase.getGeometry();

  geometry.rename(0, 1);
  geometry.communicate();

  geometry.clean();
  geometry.checkForErrors();

  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

template<size_t ID>
void prepareLatticeAD(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareLatticeAD");
  clout << "Prepare ADE Lattice ..." << std::endl;
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<Concentration<ID>>;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(Concentration<ID>{});

  Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const T Sc = parameters.get<parameters::SCHMIDT>();
  const T Re = parameters.get<parameters::REYNOLDS>();
  const T tau_ads = parameters.get<parameters::LATTICE_RELAXATION_TIME>();
  const T partRho = parameters.get<parameters::PARTICLE_DENSITY>();
  const int N = parameters.get<parameters::RESOLUTION>();

  // Set up a unit converter with the characteristic physical units
  lattice.template setUnitConverter<AdsorptionConverterFromSchmidtNumberAndRelaxation<T, DESCRIPTOR>>(
   (T) Sc,
   (T) Re,
   (T) tau_ads,
   (T) extend[0],
   (T) extend[0],
   (T) N,
   (T) 0.02,
   (T) 1,
   (T) partRho
 );
  lattice.getUnitConverter().print();

  const T omega = lattice.getUnitConverter().getLatticeRelaxationFrequency();

  // Material=1 --> bulk dynamics
  dynamics::set<SourcedAdvectionDiffusionBGKdynamics>(lattice, geometry.getMaterialIndicator(1));

  // Lattice initialize
  lattice.template setParameter<descriptors::OMEGA>(omega);

  {
    auto& communicator = lattice.getCommunicator(stage::Full());
    communicator.template requestField<descriptors::VELOCITY>();
    communicator.template requestField<descriptors::SOURCE>();
    communicator.requestOverlap(parameters.get<parameters::OVERLAP>());
    communicator.exchangeRequests();
  }

  clout << "Prepare ADE Lattice ... OK" << std::endl;
}

void prepareLatticeNS(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareLatticeNS");
  clout << "Prepare NSE Lattice ..." << std::endl;
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& latticePAD = myCase.getLattice(Concentration<0>{});
  auto& latticeCAD = myCase.getLattice(Concentration<1>{});
  auto& latticeQAD = myCase.getLattice(Concentration<2>{});

  Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const T physDeltaX = extend[0]/parameters.get<parameters::RESOLUTION>();
  const T fluidViscosity = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T Re = parameters.get<parameters::REYNOLDS>();

  // Set up a unit converter with the characteristic physical units
  lattice.setUnitConverter<UnitConverter<T,DESCRIPTOR>>(
    physDeltaX,
    latticePAD.getUnitConverter().getPhysDeltaT(),
    extend[0],
    Re * fluidViscosity / extend[0],
    fluidViscosity,
    1.225
  );
  lattice.getUnitConverter().print();

  const T omega = lattice.getUnitConverter().getLatticeRelaxationFrequency();

  // Material=1 --> bulk dynamics
  dynamics::set<BGKdynamics>(lattice, geometry.getMaterialIndicator(1));

  // Lattice initialize
  lattice.setParameter<descriptors::OMEGA>(omega);

  {
    auto &communicator = lattice.getCommunicator(stage::Full());
    communicator.requestField<descriptors::VELOCITY>();
    communicator.requestOverlap(lattice.getOverlap());
    communicator.exchangeRequests();
  }

  auto& coupling = myCase.setCouplingOperator(
    "ParticleTransport",
    AdsorptionFullCoupling3D<AdsorptionReaction<Isotherm::LinearIsotherm>,
                             ade_forces::AdvDiffDragForce3D>{},
    names::NavierStokes{}, lattice,
    names::Concentration0{}, latticePAD,
    names::Concentration1{}, latticeCAD,
    names::Concentration2{}, latticeQAD );

  // Setting Isotherm Parameters
  const T isoConstA = parameters.get<parameters::ISO_CONST_A>();
  const T isoConstB = parameters.get<parameters::ISO_CONST_B>();
  Isotherm::LinearIsotherm::setParameters<T>(isoConstA, isoConstB, coupling);

  // Setting Adsorption Reaction Parameters
  const T k_f = parameters.get<parameters::K_F>();
  const T D_s = parameters.get<parameters::D_S>();
  const T c_0 = parameters.get<parameters::C_0>();
  const T particleRadius = parameters.get<parameters::PARTICLE_RADIUS>();
  coupling.template setParameter<AdsorptionReaction<Isotherm::LinearIsotherm>::K_F>(k_f);
  coupling.template setParameter<AdsorptionReaction<Isotherm::LinearIsotherm>::D_S>(D_s);
  coupling.template setParameter<AdsorptionReaction<Isotherm::LinearIsotherm>::C_0>(c_0);
  coupling.template setParameter<AdsorptionReaction<Isotherm::LinearIsotherm>::R_P>(particleRadius);

  // Compute the interaction parameters
  AdsorptionReaction<Isotherm::LinearIsotherm>::computeParameters<T>(coupling, latticePAD.getUnitConverter());

  AdsorptionReaction<Isotherm::LinearIsotherm>::print<T>(clout, coupling, latticePAD.getUnitConverter());

  // Compute the drag force parameters
  const T partRho = parameters.get<parameters::PARTICLE_DENSITY>();
  ade_forces::AdvDiffDragForce3D::computeParametersFromRhoAndRadius<T>(partRho, particleRadius, coupling, lattice.getUnitConverter());

  T V_l = extend[0] * extend[1] * extend[2];
  T D_b = V_l*partRho * Isotherm::LinearIsotherm::getLoadingFromCoupling<T>(c_0, coupling)/(V_l*c_0);
  parameters.set<parameters::D_B>(D_b);

  clout << "Prepare NSE Lattice ... OK" << std::endl;
}

void setInitialValuesNSE(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});

  momenta::setDensity( lattice, geometry.getMaterialIndicator({0, 1}), T(1.0));

  // Make the lattice ready for simulation
  lattice.initialize();
}

template<size_t ID>
void setInitialValuesADE(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(Concentration<ID>{});
  auto& parameters = myCase.getParameters();

  // Initial conditions
  const T rho_ = parameters.get<parameters::PHYS_DENSITY_COMPONENTS>()[ID];
  momenta::setDensity( lattice, geometry.getMaterialIndicator({0, 1}), rho_);

  // Make the lattice ready for simulation
  lattice.initialize();
}

void getResults(MyCase& myCase,
                util::Timer<double>& timer,
                size_t iT)
{
  OstreamManager clout(std::cout, "getResults");
  using T = MyCase::value_t;
  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ADEDESCRIPTOR = MyCase::descriptor_t_of<Concentration<0>>;
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();
  auto& latticeNS = myCase.getLattice(NavierStokes{});
  auto& latticePAD = myCase.getLattice(Concentration<0>{});
  auto& latticeCAD = myCase.getLattice(Concentration<1>{});
  auto& latticeQAD = myCase.getLattice(Concentration<2>{});

  SuperVTMwriter3D<T> vtmWriter("adsorption3D");

  const int N = parameters.get<parameters::RESOLUTION>();
  const T maxPhysT = parameters.get<parameters::PHYS_MAX_T>();
  const int vtkIter = latticeCAD.getUnitConverter().getLatticeTime(maxPhysT/T(100))+1;

  const T particleRadius = parameters.get<parameters::PARTICLE_RADIUS>();
  const T Re = parameters.get<parameters::REYNOLDS>();
  const T partRho = parameters.get<parameters::PARTICLE_DENSITY>();

  if (iT == 0) {
    latticeNS.setProcessingContext(ProcessingContext::Evaluation);
    SuperLatticeCuboid3D<T, NSDESCRIPTOR> cuboid(latticeNS);
    SuperLatticeRank3D<T, NSDESCRIPTOR> rank(latticeNS);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();

    // Print some output of the chosen simulation setup
    clout << "N=" << N << "; maxTimeSteps=" << latticeNS.getUnitConverter().getLatticeTime(maxPhysT)
          << "; noOfCuboid=" << geometry.getCuboidDecomposition().size() << "; Re=" << Re
          << "; St="
          << (T(2) * partRho * particleRadius * particleRadius * latticeNS.getUnitConverter().getCharPhysVelocity())
              / (T(9) * latticeNS.getUnitConverter().getPhysViscosity() * latticeNS.getUnitConverter().getPhysDensity() * latticeNS.getUnitConverter().getCharPhysLength())
          << std::endl;
  }


  if (iT % vtkIter == 0) {
    // Writes the vtk files
    latticeNS.setProcessingContext(ProcessingContext::Evaluation);
    latticePAD.setProcessingContext(ProcessingContext::Evaluation);
    latticeQAD.setProcessingContext(ProcessingContext::Evaluation);
    latticeCAD.setProcessingContext(ProcessingContext::Evaluation);
    SuperGeometryF<T,3> materials(geometry);
    SuperLatticeDensity3D<T, ADEDESCRIPTOR> particles(latticePAD);
    SuperLatticeDensity3D<T, ADEDESCRIPTOR> loading(latticeQAD);
    SuperLatticeDensity3D<T, ADEDESCRIPTOR> phosphateConcentration(latticeCAD);
    SuperLatticePhysField3D<T, ADEDESCRIPTOR, descriptors::VELOCITY> extFieldParticles(latticePAD, latticeNS.getUnitConverter().getConversionFactorVelocity());
    SuperLatticePhysField3D<T, ADEDESCRIPTOR, descriptors::VELOCITY> extFieldPhosphate(latticeCAD, latticeNS.getUnitConverter().getConversionFactorVelocity());
    SuperLatticePhysField3D<T, ADEDESCRIPTOR, descriptors::SOURCE> sourcePhosphate(latticeCAD, 1);
    AnalyticalFfromSuperF3D<T> concentrationInterpolation(phosphateConcentration, true, true);
    AnalyticalFfromSuperF3D<T> loadingInterpolation(loading, true, true);

    const T D_s = parameters.get<parameters::D_S>();
    const T D_b = parameters.get<parameters::D_B>();
    const T k_f = parameters.get<parameters::K_F>();
    T k_s= T(15.) * D_s / ( particleRadius * particleRadius );
    BatchSolution<T> concentrationSol(latticeCAD.getUnitConverter().getPhysTime(iT), D_b, k_s, k_f, particleRadius);
    SuperLatticeFfromAnalyticalF3D<T, ADEDESCRIPTOR> solution(concentrationSol, latticeQAD);

    vtmWriter.addFunctor(materials);
    vtmWriter.addFunctor(particles, "particle density");
    vtmWriter.addFunctor(extFieldParticles, "particle velocity");
    vtmWriter.addFunctor(loading, "particle loading");
    vtmWriter.addFunctor(phosphateConcentration, "solute concentration");
    vtmWriter.addFunctor(sourcePhosphate, "solute source");
    vtmWriter.addFunctor(solution, "solution");
    vtmWriter.write(iT);

    timer.update(iT);
    timer.printStep();

    auto indicatorF = geometry.getMaterialIndicator(1);
    int tmp{};
    T result[2] { };
    SuperRelativeErrorL2Norm3D<T> relConcentrationError(phosphateConcentration, concentrationSol, indicatorF);
    relConcentrationError(result, &tmp);
    clout << "concentration-L2-error(rel)=" << result[0] << std::endl;

    latticeNS.getStatistics().print(iT, latticeNS.getUnitConverter().getPhysTime(iT));
  }
}

void simulate(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();

  const std::size_t iTmax = myCase.getLattice(Concentration<0>{}).getUnitConverter().getLatticeTime(
    parameters.get<parameters::PHYS_MAX_T>());

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {

    myCase.getOperator("ParticleTransport").apply();

    /// === Step 7.2.2: Computation and Output of the Results ===
    getResults(myCase, timer, iT);

    /// === Step 7.2.3: Collide and Stream Execution ===
    myCase.getLattice(Concentration<0>{}).collideAndStream();
    myCase.getLattice(Concentration<1>{}).collideAndStream();
    myCase.getLattice(Concentration<2>{}).collideAndStream();
    myCase.getLattice(NavierStokes{}).collideAndStream();
  }

  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[]) {
  /// === Step 2: Initialization ===
  initialize(&argc, &argv);

  /// === Step 2.1: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION             >(21);       // resolution of the hydraulic diameter // 20
    myCaseParameters.set<REYNOLDS               >(17);
    myCaseParameters.set<SCHMIDT                >(1.5);
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.6952380952380952);
    myCaseParameters.set<PARTICLE_RADIUS        >(1.5e-04);
    myCaseParameters.set<PARTICLE_DENSITY       >(0.940);
    myCaseParameters.set<DOMAIN_EXTENT          >({0.1,0.1,0.1});
    myCaseParameters.set<ISO_CONST_A            >(45);
    myCaseParameters.set<ISO_CONST_B            >(0.5);
    myCaseParameters.set<K_F                    >(5.37E-5);
    myCaseParameters.set<C_0                    >(1);
    myCaseParameters.set<D_S                    >(5.e-11);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY    >(0.00011764705882352942);
    myCaseParameters.set<parameters::PHYS_DENSITY_COMPONENTS  >({1,1,0});
    myCaseParameters.set<PHYS_MAX_T   >(15);    // time for fluid simulation
  }

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Lattice ===
  prepareLatticeAD<0>(myCase);
  prepareLatticeAD<1>(myCase);
  prepareLatticeAD<2>(myCase);
  prepareLatticeNS(myCase);

  setInitialValuesADE<0>(myCase);
  setInitialValuesADE<1>(myCase);
  setInitialValuesADE<2>(myCase);
  setInitialValuesNSE(myCase);

  /// === Step 7: Simulate ===
  simulate(myCase);
}
