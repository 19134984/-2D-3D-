/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2008 Orestis Malaspinas
 *                2025 Adrian Kummerlaender
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

namespace olb::parameters {
  struct T_PERTURB : public descriptors::FIELD_BASE<1> { };
}

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D2Q9<descriptors::FORCE>>,
  Temperature,  Lattice<double, descriptors::D2Q5<descriptors::VELOCITY>>
>;

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;
  Vector extent = params.get<parameters::DOMAIN_EXTENT>();
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extent, origin);

  const T physDeltaX = 0.1 / params.get<parameters::RESOLUTION>();
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({true,false});
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  geometry.rename(0,2);
  geometry.rename(2,1,{0,1});

  const T physDeltaX = params.get<parameters::PHYS_DELTA_X>();

  std::vector<T> extend( 2, T(0) );
  extend[0] = params.get<parameters::DOMAIN_EXTENT>()[0];
  extend[1] = physDeltaX;
  std::vector<T> origin( 2, T(0) );
  origin[1] = -physDeltaX/2.;
  IndicatorCuboid2D<T> bottom(extend, origin);

  origin[1] = params.get<parameters::DOMAIN_EXTENT>()[1]-physDeltaX/2.;
  IndicatorCuboid2D<T> top(extend, origin);

  origin[0] = params.get<parameters::DOMAIN_EXTENT>()[0]/2.;
  origin[1] = physDeltaX;
  extend[0] = physDeltaX;
  extend[1] = physDeltaX;
  IndicatorCuboid2D<T> perturbation(extend, origin);

  /// Set material numbers for bottom, top and pertubation
  geometry.rename(2,2,1,bottom);
  geometry.rename(2,3,1,top);
  geometry.rename(1,4,perturbation);

  geometry.clean();
  geometry.innerClean();
  geometry.checkForErrors();

  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase) {
  OstreamManager clout(std::cout,"prepareLattice");

  using T = MyCase::value_t;
  auto& geometry  = myCase.getGeometry();
  auto& params    = myCase.getParameters();

  auto& NSlattice = myCase.getLattice(NavierStokes{});
  auto& ADlattice = myCase.getLattice(Temperature{});

  using NSDESCRIPTOR  = MyCase::descriptor_t_of<NavierStokes>;
  using TDESCRIPTOR   = MyCase::descriptor_t_of<Temperature>;

  const T Tcold = params.get<parameters::T_COLD>();

  NSlattice.setUnitConverter(params);
  ADlattice.setUnitConverter(params);

  const auto& NSconverter = NSlattice.getUnitConverter();
  NSconverter.print();

  const auto& ADconverter = ADlattice.getUnitConverter();
  ADconverter.print();

  dynamics::set<ForcedBGKdynamics>(NSlattice, geometry.getMaterialIndicator({1,4}));
  boundary::set<boundary::BounceBack>(NSlattice, geometry.getMaterialIndicator({2,3}));

  dynamics::set<AdvectionDiffusionBGKdynamics>(ADlattice, geometry.getMaterialIndicator({1, 2, 3, 4}));

  boundary::set<boundary::AdvectionDiffusionDirichlet>(ADlattice, geometry, 2);
  boundary::set<boundary::AdvectionDiffusionDirichlet>(ADlattice, geometry, 3);

  NSlattice.setParameter<descriptors::OMEGA>(NSlattice.getUnitConverter().getLatticeRelaxationFrequency());
  ADlattice.setParameter<descriptors::OMEGA>(ADlattice.getUnitConverter().getLatticeThermalRelaxationFrequency());

  clout << "Setup NSE-ADE coupling" << std::endl;

  const T boussinesqForcePrefactor = 9.81 / ADconverter.getConversionFactorVelocity()
                                   * ADconverter.getConversionFactorTime()
                                   * ADconverter.getCharPhysTemperatureDifference()
                                   * ADconverter.getPhysThermalExpansionCoefficient();

  auto& coupling = myCase.setCouplingOperator(
    "Boussinesq",
    NavierStokesAdvectionDiffusionCoupling{},
    names::NavierStokes{}, NSlattice,
    names::Temperature{},  ADlattice);
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::T0>(
    ADconverter.getLatticeTemperature(Tcold));
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::FORCE_PREFACTOR>(
    boussinesqForcePrefactor * Vector<T,2>{0.0,1.0});

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& NSlattice = myCase.getLattice(NavierStokes{});
  auto& ADlattice = myCase.getLattice(Temperature{});

  const T Tcold     = params.get<parameters::T_COLD>();
  const T Thot      = params.get<parameters::T_HOT>();
  const T Tperturb  = params.get<parameters::T_PERTURB>();

  auto indicator = geometry.getMaterialIndicator({1,2,3,4});
  const auto& converter = ADlattice.getUnitConverter();

  /// define initial conditions
  AnalyticalConst2D<T,T> u0(0.0, 0.0);
  AnalyticalConst2D<T,T> PHYS_CHAR_LOW_TEMPERATURE(converter.getLatticeTemperature(Tcold));
  AnalyticalConst2D<T,T> PHYS_CHAR_HIGH_TEMPERATURE(converter.getLatticeTemperature(Thot));
  AnalyticalConst2D<T,T> T_perturb(converter.getLatticeTemperature(Tperturb));

  momenta::setTemperature(ADlattice, geometry.getMaterialIndicator({1, 3}), Tcold);
  momenta::setTemperature(ADlattice, geometry.getMaterialIndicator({2}), Thot);
  momenta::setTemperature(ADlattice, geometry.getMaterialIndicator({4}), Tperturb);

  ADlattice.iniEquilibrium(geometry, 1, PHYS_CHAR_LOW_TEMPERATURE, u0);
  ADlattice.iniEquilibrium(geometry, 2, PHYS_CHAR_HIGH_TEMPERATURE, u0);
  ADlattice.iniEquilibrium(geometry, 3, PHYS_CHAR_LOW_TEMPERATURE, u0);
  ADlattice.iniEquilibrium(geometry, 4, T_perturb, u0);

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
  auto& NSlattice = myCase.getLattice(NavierStokes{});
  auto& ADlattice = myCase.getLattice(Temperature{});
  const auto& NSconverter = NSlattice.getUnitConverter();
  const auto& ADconverter = ADlattice.getUnitConverter();
  auto& params    = myCase.getParameters();

  SuperVTMwriter2D<T> vtkWriter("rayleighBenard2d");

  const int saveIter = NSconverter.getLatticeTime(params.get<parameters::PHYS_VTK_ITER_T>());

  if (iT == 0) {
    SuperLatticeCuboid2D cuboid(NSlattice);
    SuperLatticeRank2D rank(NSlattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.createMasterFile();
  }

  if (iT % saveIter == 0) {
    NSlattice.setProcessingContext(ProcessingContext::Evaluation);
    ADlattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticePhysVelocity2D velocity(NSlattice, NSconverter);
    vtkWriter.addFunctor(velocity);
    SuperLatticePhysPressure2D pressure(NSlattice, NSconverter);
    vtkWriter.addFunctor(pressure);

    using V = MyCase::value_t_of<NavierStokes>;
    using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
    using TDESCRIPTOR = MyCase::descriptor_t_of<Temperature>;
    SuperLatticePhysTemperature2D temperature(
      ADlattice,
      ADconverter);
    vtkWriter.addFunctor(temperature);

    vtkWriter.write(iT);
  }

  myCase.report(iT);
}

void simulate(MyCase& myCase) {
  auto& parameters = myCase.getParameters();

  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(
    parameters.get<parameters::PHYS_MAX_T>());

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    myCase.getLattice(NavierStokes{}).collideAndStream();
    myCase.getLattice(Temperature{}).collideAndStream();

    myCase.getOperator("Boussinesq").apply();

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, iT);
  }

  myCase.finalize(iTmax);
}

int main(int argc, char* argv[]) {
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<DOMAIN_EXTENT            >({2, 1});
    myCaseParameters.set<RESOLUTION               >(10);
    myCaseParameters.set<RAYLEIGH                 >(1e4);
    myCaseParameters.set<PRANDTL                  >(0.71);
    myCaseParameters.set<PHYS_MAX_T               >(1000);
    myCaseParameters.set<PHYS_CHAR_LENGTH         >(0.1);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY      >(1e-5);
    myCaseParameters.set<PHYS_DENSITY             >(1.);
    myCaseParameters.set<PHYS_THERMAL_CONDUCTIVITY>(0.03);
    myCaseParameters.set<LATTICE_CHAR_VELOCITY    >(0.1);
    myCaseParameters.set<T_HOT                    >(274.15);
    myCaseParameters.set<T_COLD                   >(273.15);
    myCaseParameters.set<T_PERTURB                >(
      1. / 5. * myCaseParameters.get<T_COLD>() + 4./5. * myCaseParameters.get<T_HOT>()
    );
    myCaseParameters.set<PHYS_DELTA_X             >(
      myCaseParameters.get<PHYS_CHAR_LENGTH>() / myCaseParameters.get<RESOLUTION>()
    );
    myCaseParameters.set<PHYS_CHAR_VELOCITY       >(
      myCaseParameters.get<PHYS_CHAR_VISCOSITY>() / myCaseParameters.get<PHYS_CHAR_LENGTH>() * util::sqrt(myCaseParameters.get<RAYLEIGH>() / myCaseParameters.get<PRANDTL>())
    );
    myCaseParameters.set<PHYS_DELTA_T             >(
      myCaseParameters.get<LATTICE_CHAR_VELOCITY>() / myCaseParameters.get<PHYS_CHAR_VELOCITY>() * myCaseParameters.get<PHYS_DELTA_X>()
    );
    myCaseParameters.set<PHYS_HEAT_CAPACITY       >(
      myCaseParameters.get<PRANDTL>() * myCaseParameters.get<PHYS_THERMAL_CONDUCTIVITY>() /
      (myCaseParameters.get<PHYS_CHAR_VISCOSITY>() * myCaseParameters.get<PHYS_DENSITY>())
    );
    myCaseParameters.set<PHYS_THERMAL_EXPANSION   >(
      myCaseParameters.get<RAYLEIGH>() * util::pow(myCaseParameters.get<PHYS_CHAR_VISCOSITY>(), 2) /
      (myCaseParameters.get<PRANDTL>() * 9.81 * (myCaseParameters.get<T_HOT>() - myCaseParameters.get<T_COLD>()) * util::pow(myCaseParameters.get<PHYS_CHAR_LENGTH>(), 3))
    );
    myCaseParameters.set<PHYS_VTK_ITER_T          >(10.0);
    myCaseParameters["reporter"]["statistics"].set<PHYS_ITER_T>(10.0);
    myCaseParameters["reporter"]["timer"].set<PHYS_ITER_T>(10.0);
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
