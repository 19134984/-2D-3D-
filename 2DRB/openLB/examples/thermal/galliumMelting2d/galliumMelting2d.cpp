/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2020 Maximilian Gaedtke, Larissa Dietz
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

/* The solution for the melting problem (solid-liquid phase change)
coupled with natural convection is found using the lattice Boltzmann
method after Rongzong Huang and Huiying Wu (2015)[1]. The equilibrium
distribution function for the temperature is modified in order to deal
with the latent-heat source term. That way, iteration steps or solving
a group of linear equations is avoided, which results in enhanced efficiency.
The phase interface is located by the current total enthalpy, and
its movement is considered by the immersed moving boundary scheme after
Noble and Torczynski (1998)[2]. This method was validated by comparison
with experimental values (e.g. Gau and Viskanta (1986) [3]).

[1] Rongzong Huang, Huiying Wu, Phase interface effects in the total enthalpy-based lattice
Boltzmann model for solid–liquid phase change, Journal of Computational Physics 294 (2015) 346–362.

[2] D. Noble, J. Torczynski, A lattice-Boltzmann method for partially saturated
computational cells, Int. J. Modern Phys. C 9 (8) (1998) 1189–1202.

[3] C. Gau, R. Viskanta, Melting and Solidification of a Pure Metal on a
Vertikal Wall, Journal of Heat Transfer  108(1) (1986): 174–181.
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;

using MyCase = Case<NavierStokes,
Lattice<double, descriptors::D2Q9<descriptors::POROSITY, descriptors::VELOCITY_SOLID, descriptors::FORCE, descriptors::OMEGA>>,
Temperature, Lattice<double, descriptors::D2Q5<descriptors::VELOCITY, descriptors::TEMPERATURE>>
>;

namespace olb::parameters {
struct CP_S : public descriptors::FIELD_BASE<1> {};
struct CP_L : public descriptors::FIELD_BASE<1> {};
struct LAMBDA_S : public descriptors::FIELD_BASE<1> {};
struct LAMBDA_L : public descriptors::FIELD_BASE<1> {};
} // namespace olb::parameters

Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& params)
{
  using T            = MyCase::value_t;
  Vector  extent     = params.get<parameters::DOMAIN_EXTENT>();
  const T physDeltaX = params.get<parameters::DOMAIN_EXTENT>()[0] / params.get<parameters::RESOLUTION>();
  extent[0] += 2 * physDeltaX;
  extent[1] += physDeltaX;
  std::vector<T>       origin(2, T());
  IndicatorCuboid2D<T> cuboid(extent, origin);

  Mesh<T, MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;

  auto& params   = myCase.getParameters();
  auto& geometry = myCase.getGeometry();

  const T N  = params.get<parameters::RESOLUTION>();
  const T lx = params.get<parameters::DOMAIN_EXTENT>()[0];
  const T ly = params.get<parameters::DOMAIN_EXTENT>()[1];
  const T dx = lx / N;

  geometry.rename(0, 4);

  std::vector<T> extend(2, T());
  extend[0] = lx;
  extend[1] = ly;
  std::vector<T> origin(2, T());
  origin[0] = dx;
  origin[1] = 0.5 * dx;
  IndicatorCuboid2D<T> cuboid2(extend, origin);

  geometry.rename(4, 1, cuboid2);

  std::vector<T> extendwallleft(2, T(0));
  extendwallleft[0] = dx;
  extendwallleft[1] = ly;
  std::vector<T> originwallleft(2, T(0));
  originwallleft[0] = 0.0;
  originwallleft[1] = 0.0;
  IndicatorCuboid2D<T> wallleft(extendwallleft, originwallleft);

  std::vector<T> extendwallright(2, T(0));
  extendwallright[0] = dx;
  extendwallright[1] = ly;
  std::vector<T> originwallright(2, T(0));
  originwallright[0] = lx + dx;
  originwallright[1] = 0.0;
  IndicatorCuboid2D<T> wallright(extendwallright, originwallright);

  geometry.rename(4, 2, 1, wallleft);
  geometry.rename(4, 3, 1, wallright);

  geometry.clean();
  geometry.innerClean();
  geometry.checkForErrors();

  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T                                       = MyCase::value_t;
  using NSDESCRIPTOR                            = MyCase::descriptor_t_of<NavierStokes>;
  using TDESCRIPTOR                             = MyCase::descriptor_t_of<Temperature>;
  using TotalEnthalpyAdvectionDiffusionDynamics = TotalEnthalpyAdvectionDiffusionTRTdynamics<T, TDESCRIPTOR>;

  auto& params    = myCase.getParameters();
  auto& geometry  = myCase.getGeometry();
  auto& NSlattice = myCase.getLattice(NavierStokes {});
  auto& ADlattice = myCase.getLattice(Temperature {});

  const T N          = params.get<parameters::RESOLUTION>();
  const T tau        = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T lx         = params.get<parameters::DOMAIN_EXTENT>()[0];
  const T physDeltaX = lx / N;
  const T cp_s       = params.get<parameters::CP_S>();
  const T cp_l       = params.get<parameters::CP_L>();
  const T cp_ref     = 2.0 * cp_s * cp_l / (cp_s + cp_l);
  const T Thot       = params.get<parameters::T_HOT>();
  const T Tcold      = params.get<parameters::T_COLD>();
  const T Tmelt      = (302.8 - 301.3) / (311.0 - 301.3) + 0.5;
  const T Ra         = params.get<parameters::RAYLEIGH>();
  const T Pr         = params.get<parameters::PRANDTL>();
  const T Ste        = params.get<parameters::STEFAN>();
  const T lambda_l   = params.get<parameters::LAMBDA_L>();
  const T lambda_s   = params.get<parameters::LAMBDA_S>();
  const T R_lambda   = lambda_s / lambda_l;
  const T L          = cp_l * (Thot - Tmelt) / Ste; // J / kg
  const T physDeltaT = 6093. / 1.81e-3 / descriptors::invCs2<T, NSDESCRIPTOR>() * (tau - 0.5) * physDeltaX * physDeltaX;
  const T char_u     = util::sqrt(9.81 * 1.2e-4 * (311. - 302.8) * 6093.);

  NSlattice.setUnitConverter<ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR>>(
      (T)physDeltaX,      // physDeltaX
      (T)physDeltaT,      // physDeltaT
      (T)lx,              // charPhysLength
      (T)char_u,          // charPhysVelocity
      (T)1.81e-3 / 6093., // physViscosity
      (T)6093.,           // physDensity
      (T)32.,             // physThermalConductivity
      (T)381.,            // physSpecificHeatCapacity
      (T)1.2e-4,          // physThermalExpansionCoefficient
      (T)Tcold,           // charPhysLowTemperature
      (T)Thot             // charPhysHighTemperature
  );

  auto& converter = NSlattice.getUnitConverter();
  converter.print();

  ADlattice.setUnitConverter(converter);

  auto& coupling = myCase.setCouplingOperator("Boussinesq", TotalEnthalpyPhaseChangeCoupling {}, names::NavierStokes {},
                                              NSlattice, names::Temperature {}, ADlattice);
  coupling.restrictTo(geometry.getMaterialIndicator({1}));

  dynamics::set<ForcedPSMBGKdynamics>(NSlattice, geometry.getMaterialIndicator({1, 2, 3, 4}));
  dynamics::set<TotalEnthalpyAdvectionDiffusionDynamics>(ADlattice, geometry.getMaterialIndicator({1, 2, 3}));
  boundary::set<boundary::BounceBack>(ADlattice, geometry, 4);

  boundary::set<boundary::RegularizedTemperature>(ADlattice, geometry.getMaterialIndicator({2, 3}));
  boundary::set<boundary::InterpolatedVelocity>(NSlattice, geometry.getMaterialIndicator({2, 3, 4}));

  T omega  = converter.getLatticeRelaxationFrequency();
  T Tomega = converter.getLatticeThermalRelaxationFrequency();

  std::vector<T> dir {0.0, 1.0};
  std::vector<T> forcePrefactor {0, 0};
  T              boussinesqForcePrefactor =
      Ra / util::pow(T(N), 3) * Pr *
      util::pow(cp_ref / descriptors::invCs2<T, TDESCRIPTOR>() * (converter.getLatticeThermalRelaxationTime() - 0.5),
                2);
  clout << "boussinesq " << Ra / util::pow(T(N), 3) * Pr * lambda_l * lambda_l << std::endl;

  T normDir = T();
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    normDir += dir[iD] * dir[iD];
  }
  normDir = util::sqrt(normDir);
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    dir[iD] /= normDir;
  }

  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    forcePrefactor[iD] = boussinesqForcePrefactor * dir[iD];
  }

  NSlattice.setParameter<descriptors::OMEGA>(omega);
  ADlattice.setParameter<descriptors::OMEGA>(Tomega);

  ADlattice.setParameter<collision::TRT::MAGIC>(T(0.25));

  ADlattice.setParameter<TotalEnthalpy::T_S>(Tmelt);
  ADlattice.setParameter<TotalEnthalpy::T_L>(Tmelt);
  ADlattice.setParameter<TotalEnthalpy::CP_S>(cp_s);
  ADlattice.setParameter<TotalEnthalpy::CP_L>(cp_l);
  ADlattice.setParameter<TotalEnthalpy::LAMBDA_S>(cp_ref / descriptors::invCs2<T, TDESCRIPTOR>() *
                                                  (converter.getLatticeThermalRelaxationTime() - 0.5) * R_lambda);
  ADlattice.setParameter<TotalEnthalpy::LAMBDA_L>(cp_ref / descriptors::invCs2<T, TDESCRIPTOR>() *
                                                  (converter.getLatticeThermalRelaxationTime() - 0.5));
  ADlattice.setParameter<TotalEnthalpy::L>(L);

  coupling.template setParameter<TotalEnthalpyPhaseChangeCoupling::T_S>(Tmelt);
  coupling.template setParameter<TotalEnthalpyPhaseChangeCoupling::T_L>(Tmelt);
  coupling.template setParameter<TotalEnthalpyPhaseChangeCoupling::CP_S>(cp_s);
  coupling.template setParameter<TotalEnthalpyPhaseChangeCoupling::CP_L>(cp_l);
  coupling.template setParameter<TotalEnthalpyPhaseChangeCoupling::L>(L);
  coupling.template setParameter<TotalEnthalpyPhaseChangeCoupling::FORCE_PREFACTOR>(forcePrefactor);
  coupling.template setParameter<TotalEnthalpyPhaseChangeCoupling::T_COLD>(Tcold);
  coupling.template setParameter<TotalEnthalpyPhaseChangeCoupling::DELTA_T>(T(1.));
}

void setInitialValues(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");

  using T         = MyCase::value_t;
  auto& params    = myCase.getParameters();
  auto& geometry  = myCase.getGeometry();
  auto& NSlattice = myCase.getLattice(NavierStokes {});
  auto& ADlattice = myCase.getLattice(Temperature {});
  auto& converter = NSlattice.getUnitConverter();

  const T cp_s  = params.get<parameters::CP_S>();
  const T cp_l  = params.get<parameters::CP_L>();
  const T Thot  = params.get<parameters::T_HOT>();
  const T Tcold = params.get<parameters::T_COLD>();

  const T lattice_Hcold = cp_s * Tcold;
  const T lattice_Hhot  = cp_l * Thot;

  T omega = converter.getLatticeRelaxationFrequency();

  AnalyticalConst2D<T, T> rho(1.);
  AnalyticalConst2D<T, T> u0(0.0, 0.0);
  AnalyticalConst2D<T, T> T_cold(lattice_Hcold);
  AnalyticalConst2D<T, T> T_hot(lattice_Hhot);
  AnalyticalConst2D<T, T> omegaField(omega);
  fields::set<descriptors::OMEGA>(NSlattice, geometry.getMaterialIndicator({1, 2, 3, 4}), omegaField);

  fields::set<descriptors::VELOCITY>(ADlattice, geometry.getMaterialIndicator({1, 2, 3}), u0);
  momenta::setTemperature(ADlattice, geometry.getMaterialIndicator({1, 3}), T_cold);
  ADlattice.iniEquilibrium(geometry.getMaterialIndicator({1, 3}), T_cold, u0);
  momenta::setTemperature(ADlattice, geometry.getMaterialIndicator({2}), T_hot);
  ADlattice.iniEquilibrium(geometry, 2, T_hot, u0);

  NSlattice.initialize();
  ADlattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void getResults(MyCase& myCase, util::Timer<MyCase::value_t>& timer, std::size_t iT)
{
  OstreamManager clout(std::cout, "getResults");

  using T            = MyCase::value_t;
  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using TDESCRIPTOR  = MyCase::descriptor_t_of<Temperature>;

  auto&       NSlattice = myCase.getLattice(NavierStokes {});
  auto&       ADlattice = myCase.getLattice(Temperature {});
  const auto& converter = NSlattice.getUnitConverter();

  SuperVTMwriter2D<T>                                        vtkWriter("galliumMelting2d");
  SuperLatticeField2D<T, TDESCRIPTOR, descriptors::VELOCITY> velocity(ADlattice);
  SuperLatticePhysPressure2D<T, NSDESCRIPTOR>                pressure(NSlattice, converter);

  SuperLatticeDensity2D<T, TDESCRIPTOR> enthalpy(ADlattice);
  enthalpy.getName() = "enthalpy";
  SuperLatticeField2D<T, NSDESCRIPTOR, descriptors::POROSITY> liquid_frac(NSlattice);
  liquid_frac.getName() = "liquid fraction";
  SuperLatticeField2D<T, TDESCRIPTOR, descriptors::TEMPERATURE> temperature(ADlattice);
  temperature.getName() = "temperature";
  SuperLatticeField2D<T, NSDESCRIPTOR, descriptors::FORCE> force(NSlattice);
  force.getName() = "force";
  vtkWriter.addFunctor(pressure);
  vtkWriter.addFunctor(velocity);
  vtkWriter.addFunctor(enthalpy);
  vtkWriter.addFunctor(liquid_frac);
  vtkWriter.addFunctor(temperature);
  vtkWriter.addFunctor(force);

  const int vtkIter = converter.getLatticeTime(0.5);

  if (iT == 0) {
    SuperLatticeCuboid2D<T, NSDESCRIPTOR> cuboid(NSlattice);
    SuperLatticeRank2D<T, NSDESCRIPTOR>   rank(NSlattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
  }

  if (iT % vtkIter == 0) {
    NSlattice.setProcessingContext(ProcessingContext::Evaluation);
    ADlattice.setProcessingContext(ProcessingContext::Evaluation);

    timer.update(iT);
    timer.printStep();

    NSlattice.getStatistics().print(iT, converter.getPhysTime(iT));

    ADlattice.getStatistics().print(iT, converter.getPhysTime(iT));

    vtkWriter.write(iT);
  }
}
void simulate(MyCase& myCase)
{
  using T = MyCase::value_t;

  auto& params    = myCase.getParameters();
  auto& geometry  = myCase.getGeometry();
  auto& NSlattice = myCase.getLattice(NavierStokes {});
  auto& ADlattice = myCase.getLattice(Temperature {});
  auto& converter = NSlattice.getUnitConverter();

  const T maxPhysT = params.get<parameters::PHYS_MAX_T>();

  util::Timer<T> timer(converter.getLatticeTime(maxPhysT), geometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT) + 1; ++iT) {
    myCase.getOperator("Boussinesq").apply();
    NSlattice.collideAndStream();
    ADlattice.collideAndStream();

    getResults(myCase, timer, iT);
  }

  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[])
{
  OstreamManager clout(std::cout, "main");
  initialize(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    using T = MyCase::value_t;
    myCaseParameters.set<RESOLUTION>(128);
    myCaseParameters.set<DOMAIN_EXTENT>({(T)88.9e-3, (T)63.5e-3});
    myCaseParameters.set<LATTICE_RELAXATION_TIME>((T)0.51);
    myCaseParameters.set<RAYLEIGH>((T)2e6);
    myCaseParameters.set<PRANDTL>((T)0.0216);
    myCaseParameters.set<STEFAN>((T)0.039);
    myCaseParameters.set<PHYS_MAX_T>((T)1140.);
    myCaseParameters.set<T_COLD>((T)0.5);
    myCaseParameters.set<T_HOT>((T)1.5);
    myCaseParameters.set<LAMBDA_S>((T)33.5);
    myCaseParameters.set<LAMBDA_L>((T)32.0);
    myCaseParameters.set<CP_S>((T)1.0);
    myCaseParameters.set<CP_L>((T)1.0);
  }

  Mesh mesh = createMesh(myCaseParameters);

  MyCase myCase(myCaseParameters, mesh);

  prepareGeometry(myCase);

  prepareLattice(myCase);

  setInitialValues(myCase);

  simulate(myCase);
}
