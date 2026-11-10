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

namespace olb::parameters {
  struct T_PERTURB  : public descriptors::FIELD_BASE<1> { };
}

// === Step 1: Declarations ===
using MyCase = Case<
  names::NavierStokes, Lattice<double, descriptors::D2Q9<descriptors::FORCE>>,
  names::Temperature,  Lattice<double, descriptors::D2Q5<descriptors::VELOCITY>>
>;

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extent, origin);

  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({true,false});
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;

  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  geometry.rename(0,2);
  geometry.rename(2,1,{0,1});

  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();

  std::vector<T> extend( 2, T(0) );
  extend[0] = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  extend[1] = physDeltaX;
  std::vector<T> origin( 2, T(0) );
  IndicatorCuboid2D<T> bottom(extend, origin);

  origin[1] = parameters.get<parameters::DOMAIN_EXTENT>()[1]-physDeltaX;
  IndicatorCuboid2D<T> top(extend, origin);

  origin[0] = parameters.get<parameters::DOMAIN_EXTENT>()[0]/2.;
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

  using T             = MyCase::value_t;
  using NSEDESCRIPTOR  = MyCase::descriptor_t_of<names::NavierStokes>;
  using ADEDESCRIPTOR   = MyCase::descriptor_t_of<names::Temperature>;

  auto& geometry      = myCase.getGeometry();
  auto& NSElattice    = myCase.getLattice(names::NavierStokes{});
  auto& ADElattice    = myCase.getLattice(names::Temperature{});
  auto& parameters    = myCase.getParameters();

  const T Ra                              = parameters.get<parameters::RAYLEIGH>();
  const T Pr                              = parameters.get<parameters::PRANDTL>();
  const T g                               = parameters.get<parameters::GRAVITATIONAL_ACC>();

  const T physCharLength                  = parameters.get<parameters::PHYS_CHAR_LENGTH>();
  const T physDensity                     = parameters.get<parameters::PHYS_CHAR_DENSITY>();
  const T physThermalConductivity         = parameters.get<parameters::PHYS_THERMAL_CONDUCTIVITY>();
  const T physViscosity                   = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T Tcold                           = parameters.get<parameters::T_COLD>();
  const T Thot                            = parameters.get<parameters::T_HOT>();

  const T physDeltaX                      = parameters.get<parameters::PHYS_DELTA_X>();
  const T latticeCharVelocity             = parameters.get<parameters::LATTICE_CHAR_VELOCITY>();

  const T physCharVelocity                = physViscosity / physCharLength * util::sqrt( Ra / Pr );
  const T physDeltaT                      = latticeCharVelocity / physCharVelocity * physDeltaX;
  const T physThermalExpansionCoefficient = Ra * physViscosity * physViscosity / (Pr * g * (Thot - Tcold) * physCharLength * physCharLength * physCharLength);
  const T physSpecificHeatCapacity        = Pr * physThermalConductivity / (physViscosity * physDensity);

  NSElattice.setUnitConverter<ThermalUnitConverter<T,NSEDESCRIPTOR,ADEDESCRIPTOR>>(
    physDeltaX,
    physDeltaT,
    physCharLength,
    physCharVelocity,
    physViscosity,
    physDensity,
    physThermalConductivity,
    physSpecificHeatCapacity,
    physThermalExpansionCoefficient,
    Tcold,
    Thot
  );
  const auto& converter = NSElattice.getUnitConverter();
  converter.print();

  ADElattice.setUnitConverter(converter);

  dynamics::set<ForcedBGKdynamics>(NSElattice, geometry.getMaterialIndicator({1,4}));
  boundary::set<boundary::BounceBack>(NSElattice, geometry.getMaterialIndicator({2,3}));

  dynamics::set<AdvectionDiffusionBGKdynamics>(ADElattice, geometry.getMaterialIndicator({1, 2, 3, 4}));

  boundary::set<boundary::AdvectionDiffusionDirichlet>(ADElattice, geometry, 2);
  boundary::set<boundary::AdvectionDiffusionDirichlet>(ADElattice, geometry, 3);

  NSElattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  ADElattice.setParameter<descriptors::OMEGA>(converter.getLatticeThermalRelaxationFrequency());

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  auto& NSElattice = myCase.getLattice(names::NavierStokes{});
  auto& ADElattice = myCase.getLattice(names::Temperature{});

  const T Tcold     = parameters.get<parameters::T_COLD>();
  const T Thot      = parameters.get<parameters::T_HOT>();
  const T Tperturb  = parameters.get<parameters::T_PERTURB>();

  auto indicator = geometry.getMaterialIndicator({1,2,3,4});
  const auto& converter = NSElattice.getUnitConverter();

  AnalyticalConst2D<T,T> u0(0.0, 0.0);
  AnalyticalConst2D<T,T> T_cold(converter.getLatticeTemperature(Tcold));
  AnalyticalConst2D<T,T> T_hot(converter.getLatticeTemperature(Thot));
  AnalyticalConst2D<T,T> T_perturb(converter.getLatticeTemperature(Tperturb));

  momenta::setTemperature(ADElattice, geometry.getMaterialIndicator({1, 3}), Tcold);
  momenta::setTemperature(ADElattice, geometry.getMaterialIndicator({2}), Thot);
  momenta::setTemperature(ADElattice, geometry.getMaterialIndicator({4}), Tperturb);

  ADElattice.iniEquilibrium(geometry, 1, T_cold, u0);
  ADElattice.iniEquilibrium(geometry, 2, T_hot, u0);
  ADElattice.iniEquilibrium(geometry, 3, T_cold, u0);
  ADElattice.iniEquilibrium(geometry, 4, T_perturb, u0);

  NSElattice.initialize();
  ADElattice.initialize();
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{ }

void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  using T = MyCase::value_t;
  auto& NSElattice = myCase.getLattice(names::NavierStokes{});
  auto& ADElattice = myCase.getLattice(names::Temperature{});
  const auto& converter = NSElattice.getUnitConverter();

  SuperVTMwriter2D<T> vtkWriter("rayleighBenard2d");

  const int statIter = converter.getLatticeTime(10.0);
  const int saveIter = converter.getLatticeTime(10.0);

  if (iT == 0) {
    SuperLatticeCuboid2D cuboid(NSElattice);
    SuperLatticeRank2D rank(NSElattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.createMasterFile();
  }

  if (iT % statIter == 0) {
    timer.update(iT);
    timer.printStep();
    NSElattice.getStatistics().print(iT, converter.getPhysTime(iT));
  }

  if (iT % saveIter == 0) {
    NSElattice.setProcessingContext(ProcessingContext::Evaluation);
    ADElattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticePhysVelocity2D velocity(NSElattice, converter);
    vtkWriter.addFunctor(velocity);
    SuperLatticePhysPressure2D pressure(NSElattice, converter);
    vtkWriter.addFunctor(pressure);

    using V = MyCase::value_t_of<names::NavierStokes>;
    using NSEDESCRIPTOR = MyCase::descriptor_t_of<names::NavierStokes>;
    using ADEDESCRIPTOR = MyCase::descriptor_t_of<names::Temperature>;
    SuperLatticePhysTemperature2D temperature(
      ADElattice,
      static_cast<const ThermalUnitConverter<V,NSEDESCRIPTOR,ADEDESCRIPTOR>&>(converter));
    vtkWriter.addFunctor(temperature);

    vtkWriter.write(iT);
  }
}

void simulate(MyCase& myCase) {
  using T = MyCase::value_t;

  auto& geometry   = myCase.getGeometry();
  auto& NSElattice = myCase.getLattice(names::NavierStokes{});
  auto& ADElattice = myCase.getLattice(names::Temperature{});
  const auto& converter = NSElattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  const std::size_t iTmax = converter.getLatticeTime(
    parameters.get<parameters::PHYS_MAX_T>()
  );

  util::Timer<T> timer(iTmax, geometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT <= iTmax; ++iT) {

    setTemporalValues(myCase, iT);

    NSElattice.collideAndStream();
    ADElattice.collideAndStream();

    getResults(myCase, timer, iT);
  }

  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[]) {
  initialize(&argc, &argv);

  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION               >(10    );
    myCaseParameters.set<PHYS_MAX_T               >(1000  );
    myCaseParameters.set<DOMAIN_EXTENT            >({2, 1});
    myCaseParameters.set<PHYS_VTK_ITER_T          >(1.0   );
    myCaseParameters.set<PHYS_STAT_ITER_T         >(0.1);

    myCaseParameters.set<RAYLEIGH                 >(1e4   );
    myCaseParameters.set<PRANDTL                  >(0.71  );
    myCaseParameters.set<GRAVITATIONAL_ACC        >(9.81);
    myCaseParameters.set<PHYS_CHAR_LENGTH         >(    0.1  );
    myCaseParameters.set<PHYS_CHAR_VELOCITY       >(    1.0  );
    myCaseParameters.set<PHYS_CHAR_DENSITY        >(    1.0  );
    myCaseParameters.set<PHYS_CHAR_VISCOSITY      >(    1e-5 );
    myCaseParameters.set<PHYS_THERMAL_CONDUCTIVITY>(    0.03  );
    myCaseParameters.set<T_HOT                    >(274.15);
    myCaseParameters.set<T_COLD                   >(273.15);
    myCaseParameters.set<LATTICE_CHAR_VELOCITY    >(  0.1);

    myCaseParameters.set<PHYS_DELTA_X             >(
      myCaseParameters.get<PHYS_CHAR_LENGTH>() / myCaseParameters.get<RESOLUTION>()
    );
    myCaseParameters.set<T_PERTURB        >(
      1./5. * myCaseParameters.get<T_COLD>() + 4./5. * myCaseParameters.get<T_HOT>()
    );
  }

  Mesh mesh = createMesh(myCaseParameters);

  MyCase myCase(myCaseParameters, mesh);

  prepareGeometry(myCase);

  prepareLattice(myCase);

  setInitialValues(myCase);

  simulate(myCase);
}
