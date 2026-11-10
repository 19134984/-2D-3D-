/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006, 2007, 2008 Jonas Latt, Orestis Malaspina, Andrea Parmigiani
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

/* rayleighBenard3d.cpp:
 * Rayleigh-Benard convection rolls in 3D, simulated with
 * the thermal LB model by Z. Guo e.a., between a hot plate at
 * the bottom and a cold plate at the top.
 */


#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<FLOATING_POINT_TYPE, descriptors::D3Q19<descriptors::FORCE>>,
  Temperature,  Lattice<FLOATING_POINT_TYPE, descriptors::D3Q7<descriptors::VELOCITY>>
>;

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  std::vector<T> origin(3, T());
  IndicatorCuboid3D<T> cuboid(extent, origin);

  const T physDeltaX = parameters.get<parameters::PHYS_CHAR_LENGTH>() / parameters.get<parameters::RESOLUTION>();
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({true, false, true});
  return mesh;
}

void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();

  const T physDeltaX  = parameters.get<parameters::PHYS_CHAR_LENGTH>() / parameters.get<parameters::RESOLUTION>();
  const T physLengthX = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T physLengthY = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const T physLengthZ = parameters.get<parameters::DOMAIN_EXTENT>()[2];

  // Sets material number for fluid and boundary
  geometry.rename(0, 2);
  geometry.rename(2, 1, {0, 1, 0});

  std::vector<T> extent(3, T(0));
  extent[0] = physLengthX;
  extent[1] = physDeltaX;
  extent[2] = physLengthZ;
  std::vector<T> origin(3, T(0));
  IndicatorCuboid3D<T> bottom(extent, origin);

  origin[1] = physLengthY-physDeltaX;
  IndicatorCuboid3D<T> top(extent, origin);

  origin[0] = physLengthX/2.;
  origin[1] = physDeltaX;
  origin[2] = physLengthZ/2.;
  extent = {physDeltaX, physDeltaX, physDeltaX};
  IndicatorCuboid3D<T> perturbation(extent, origin);

  /// Set material numbers for bottom, top and pertubation
  geometry.rename(2, 2, 1, bottom);
  geometry.rename(2, 3, 1, top);
  geometry.rename(1, 4, perturbation);

  /// Removes all not needed boundary voxels outside the surface
  geometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  geometry.innerClean();
  geometry.checkForErrors();
  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout,"prepareLattice");

  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();

  auto& NSlattice = myCase.getLattice(NavierStokes{});
  auto& ADlattice = myCase.getLattice(Temperature{});

  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using TDESCRIPTOR  = MyCase::descriptor_t_of<Temperature>;

  const T charPhysLength          = parameters.get<parameters::PHYS_CHAR_LENGTH>();
  const T charPhysDensity         = parameters.get<parameters::PHYS_CHAR_DENSITY>();
  const T physViscosity           = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T physThermalConductivity = parameters.get<parameters::PHYS_THERMAL_CONDUCTIVITY>();
  const T physCp                  = parameters.get<parameters::PHYS_CP>();
  const T physDeltaX              = charPhysLength / parameters.get<parameters::RESOLUTION>();
  const T Ra                      = parameters.get<parameters::RAYLEIGH>();
  const T Pr                      = parameters.get<parameters::PRANDTL>();
  const T Tcold                   = parameters.get<parameters::T_COLD>();
  const T Thot                    = parameters.get<parameters::T_HOT>();

  NSlattice.setUnitConverter<ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR>>(
    physDeltaX,
    charPhysLength / (physViscosity / charPhysLength* util::sqrt(Ra/Pr)) * physDeltaX, // physDeltaT = charLatticeVelocity / charPhysVelocity * physDeltaX
    charPhysLength,
    physViscosity / charPhysLength * util::sqrt(Ra/Pr), // charPhysVelocity
    physViscosity,
    charPhysDensity,
    physThermalConductivity,
    physCp,
    Ra * physViscosity * physViscosity / Pr / 9.81 / (Thot - Tcold) / util::pow(charPhysLength, 3),  // physThermalExpansionCoefficient
    Tcold,
    Thot
  );

  const auto& converter = NSlattice.getUnitConverter();
  converter.print();

  ADlattice.setUnitConverter(converter);

  dynamics::set<ForcedBGKdynamics>(NSlattice, geometry, 1);
  boundary::set<boundary::BounceBack>(NSlattice, geometry, 2);
  boundary::set<boundary::BounceBack>(NSlattice, geometry, 3);
  dynamics::set<ForcedBGKdynamics>(NSlattice, geometry, 4);

  dynamics::set<AdvectionDiffusionBGKdynamics>(ADlattice,
                                               geometry.getMaterialIndicator({1, 2, 3, 4}));
  boundary::set<boundary::AdvectionDiffusionDirichlet>(ADlattice, geometry, 2);
  boundary::set<boundary::AdvectionDiffusionDirichlet>(ADlattice, geometry, 3);

  NSlattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  ADlattice.setParameter<descriptors::OMEGA>(converter.getLatticeThermalRelaxationFrequency());

  clout << "Setup NSE-ADE coupling" << std::endl;

  T boussinesqForcePrefactor = 9.81 / converter.getConversionFactorVelocity() * converter.getConversionFactorTime() *
                               converter.getCharPhysTemperatureDifference() * converter.getPhysThermalExpansionCoefficient();

  auto& coupling = myCase.setCouplingOperator(
    "Boussinesq",
    NavierStokesAdvectionDiffusionCoupling{},
    names::NavierStokes{}, NSlattice,
    names::Temperature{}, ADlattice);

  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::T0>(
    converter.getLatticeTemperature(Tcold));
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::FORCE_PREFACTOR>(
    boussinesqForcePrefactor * Vector<T,3>{0.0, 1.0, 0.0});

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry   = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  auto& NSlattice = myCase.getLattice(NavierStokes{});
  auto& ADlattice = myCase.getLattice(Temperature{});

  const T Tcold    = parameters.get<parameters::T_COLD>();
  const T Thot     = parameters.get<parameters::T_HOT>();
  const T Tperturb = parameters.get<parameters::T_PERTURBATION>();

  momenta::setTemperature(ADlattice, geometry.getMaterialIndicator(1), Tcold);
  momenta::setTemperature(ADlattice, geometry.getMaterialIndicator(2), Thot);
  momenta::setTemperature(ADlattice, geometry.getMaterialIndicator(3), Tcold);
  momenta::setTemperature(ADlattice, geometry.getMaterialIndicator(4), Tperturb);

  NSlattice.initialize();
  ADlattice.initialize();
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{ }

void getResults(MyCase& myCase,
                std::size_t iT)

{
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();

  auto& NSlattice = myCase.getLattice(NavierStokes{});
  auto& ADlattice = myCase.getLattice(Temperature{});
  const auto& converter = NSlattice.getUnitConverter();

  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using TDESCRIPTOR  = MyCase::descriptor_t_of<Temperature>;

  const bool converged = parameters.get<parameters::CONVERGED>();
  const T Tcold        = parameters.get<parameters::T_COLD>();
  const T Thot         = parameters.get<parameters::T_HOT>();
  const T physLengthZ  = parameters.get<parameters::DOMAIN_EXTENT>()[2];

  if (iT == 0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperVTMwriter3D<T> vtkWriter("rayleighBenard3d");
    SuperLatticeCuboid3D<T, NSDESCRIPTOR> cuboid(NSlattice);
    SuperLatticeRank3D<T, NSDESCRIPTOR> rank(NSlattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.createMasterFile();
  }

  const int vtkIter  = converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T>());

  /// Writes the VTK files and prints statistics
  if (iT%vtkIter == 0 || converged) {
    ADlattice.setProcessingContext(ProcessingContext::Evaluation);
    NSlattice.setProcessingContext(ProcessingContext::Evaluation);

    using V = MyCase::value_t_of<NavierStokes>;

    NSlattice.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtkWriter("rayleighBenard3d");
      SuperLatticePhysVelocity3D velocity(NSlattice, converter);
      SuperLatticePhysPressure3D pressure(NSlattice, converter);
      SuperLatticePhysTemperature3D<V,NSDESCRIPTOR,TDESCRIPTOR> temperature(ADlattice, static_cast<const ThermalUnitConverter<V,NSDESCRIPTOR,TDESCRIPTOR>&>(converter));
      vtkWriter.addFunctor(pressure);
      vtkWriter.addFunctor(velocity);
      vtkWriter.addFunctor(temperature);
      task(vtkWriter, iT);
    });

    SuperLatticePhysTemperature3D<T, NSDESCRIPTOR, TDESCRIPTOR> temperature(ADlattice, static_cast<const ThermalUnitConverter<V,NSDESCRIPTOR,TDESCRIPTOR>&>(converter));
    BlockReduction3D2D<T> planeReduction(temperature, {0, 0, physLengthZ/T(2)}, {0, 0, 1});
    BlockGifWriter<T> gifWriter;
    gifWriter.write(planeReduction, Tcold-T(0.1), Thot+T(0.1), iT, "temperature");
  }

  myCase.report(iT);
}

void simulate(MyCase& myCase) {
  OstreamManager clout(std::cout,"simulate");
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();

  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(
    parameters.get<parameters::PHYS_MAX_T>());

  util::ValueTracer<T> converge(myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(50.),
                                parameters.get<parameters::CONVERGENCE_PRECISION>());

  for (std::size_t iT=0; iT < iTmax; ++iT) {

    if (converge.hasConverged()){
      parameters.set<parameters::CONVERGED>(true);
      getResults(myCase, iT);
      clout << "Converged after " << iT << " iterations." << std::endl;
      break;
    }

    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    myCase.getLattice(NavierStokes{}).collideAndStream();
    myCase.getLattice(Temperature{}).collideAndStream();

    myCase.getOperator("Boussinesq").apply();

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, iT);
    converge.takeValue(myCase.getLattice(Temperature{}).getStatistics().getAverageEnergy(), true);
  }

  myCase.getLattice(NavierStokes{}).setProcessingContext(ProcessingContext::Evaluation);
  myCase.getLattice(Temperature{}).setProcessingContext(ProcessingContext::Evaluation);

  myCase.finalize(iTmax);
}

int main(int argc, char *argv[])
{
  /// === 1st Step: Initialization ===
  OstreamManager clout(std::cout,"main");
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<DOMAIN_EXTENT>({0.2, 0.1, 0.1});
    myCaseParameters.set<RESOLUTION>(40);
    myCaseParameters.set<RAYLEIGH>(1e6);
    myCaseParameters.set<PRANDTL>(0.71);
    myCaseParameters.set<PHYS_MAX_T>(20.);
    myCaseParameters.set<T_HOT>(274.15);
    myCaseParameters.set<T_COLD>(273.15);
    myCaseParameters.set<PHYS_CHAR_LENGTH>(0.1);
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1.0);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(1e-5);
    myCaseParameters.set<PHYS_THERMAL_CONDUCTIVITY>(0.03);

    myCaseParameters.set<T_PERTURBATION>([] (MyCase::ParametersD& p) {
      return p.get<T_COLD>()/5. + 4./5.*p.get<T_HOT>();
    });

    myCaseParameters.set<PHYS_CP>([] (MyCase::ParametersD& p) {
      return p.get<PRANDTL>() * p.get<PHYS_THERMAL_CONDUCTIVITY>()
              / (p.get<PHYS_CHAR_VISCOSITY>()/p.get<PHYS_CHAR_DENSITY>());
    });

    myCaseParameters.set<CONVERGENCE_PRECISION>(1e-5);
    myCaseParameters.set<CONVERGED>(false);

    myCaseParameters.set<PHYS_VTK_ITER_T>(1.0);
    myCaseParameters["reporter"]["statistics"].set<PHYS_ITER_T>(0.1);
    myCaseParameters["reporter"]["timer"].set<PHYS_ITER_T>(0.1);
  }

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase);

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);
}
