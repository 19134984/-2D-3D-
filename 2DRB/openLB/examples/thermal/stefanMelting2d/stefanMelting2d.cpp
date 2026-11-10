/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2020 Maximilian Gaedke, Larissa Dietz
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
is found using the lattice Boltzmann method after Rongzong Huang
and Huiying Wu (2015)[1]. The equilibrium distribution
function for the temperature is modified in order to deal with the
latent-heat source term. That way, iteration steps or solving a group
of linear equations is avoided, which results in enhanced efficiency.
The phase interface is located by the current total enthalpy, and
its movement is considered by the immersed moving boundary scheme after
Noble and Torczynski (1998)[2]. This method was validated by the
problem of conduction-induced melting in a semi-infinite space,
comparing its results to analytical solutions.

[1] Rongzong Huang, Huiying Wu, Phase interface effects in the total enthalpy-based lattice
Boltzmann model for solid–liquid phase change, Journal of Computational Physics 294 (2015) 346–362.

[2] D. Noble, J. Torczynski, A lattice-Boltzmann method for partially saturated
computational cells, Int. J. Modern Phys. C 9 (8) (1998) 1189–1202.
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;

using MyCase = Case<
NavierStokes, Lattice<double, descriptors::D2Q9<descriptors::POROSITY, descriptors::VELOCITY_SOLID, descriptors::FORCE, descriptors::OMEGA>>,
Temperature, Lattice<double, descriptors::D2Q5<descriptors::VELOCITY, descriptors::TEMPERATURE>>
>;

namespace olb::parameters {
struct CP_S : public descriptors::FIELD_BASE<1> {};
struct CP_L : public descriptors::FIELD_BASE<1> {};
struct NUM_ERROR_SAMPLES : public descriptors::FIELD_BASE<1> {};
struct K : public descriptors::FIELD_BASE<1> {};
struct ITERATIONS : public descriptors::FIELD_BASE<1> {};
struct MAX_ERROR : public descriptors::FIELD_BASE<1> {};
struct STARTING_VALUE : public descriptors::FIELD_BASE<1> {};
struct SUM_ERROR_MELT : public descriptors::FIELD_BASE<1> {};
struct SUM_ERROR_TEMP : public descriptors::FIELD_BASE<1> {};
struct SUM_ERROR_MELT_INF : public descriptors::FIELD_BASE<1> {};
struct SUM_ERROR_TEMP_INF : public descriptors::FIELD_BASE<1> {};
struct NUM_ERRORS : public descriptors::FIELD_BASE<1> {};
} // namespace olb::parameters

Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& params)
{
  using T            = MyCase::value_t;
  Vector  extent     = params.get<parameters::DOMAIN_EXTENT>();
  const T physDeltaX = params.get<parameters::DOMAIN_EXTENT>()[0] / params.get<parameters::RESOLUTION>();

  std::vector<T>       origin(2, T());
  IndicatorCuboid2D<T> cuboid(extent, origin);

  Mesh<T, MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({false, true});
  return mesh;
}

void getError(MyCase& myCase, int iT)
{
  using T            = MyCase::value_t;
  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using TDESCRIPTOR  = MyCase::descriptor_t_of<Temperature>;

  auto&   params     = myCase.getParameters();
  auto&   geometry   = myCase.getGeometry();
  auto&   NSlattice  = myCase.getLattice(NavierStokes {});
  auto&   ADlattice  = myCase.getLattice(Temperature {});
  const T N          = params.get<parameters::RESOLUTION>();
  const T tau        = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T lx         = params.get<parameters::DOMAIN_EXTENT>()[0];
  const T cp_s       = params.get<parameters::CP_S>();
  const T cp_l       = params.get<parameters::CP_L>();
  const T Thot       = params.get<parameters::T_HOT>();
  const T Tcold      = params.get<parameters::T_COLD>();
  const T lambda     = params.get<parameters::PHYS_LAMBDA>();
  const T density    = params.get<parameters::PHYS_DENSITY>();
  const T k          = params.get<parameters::K>();
  const T cp_ref     = 2.0 * cp_s * cp_l / (cp_s + cp_l);
  const T physDeltaX = lx / N;
  const T physDeltaT =
      density * cp_ref / lambda / descriptors::invCs2<T, NSDESCRIPTOR>() * (tau - 0.5) * physDeltaX * physDeltaX;
  OstreamManager clout(std::cout, "error");

  int input[1]  = {};
  T   result[1] = {};

  auto indicatorF = geometry.getMaterialIndicator({1});

  SuperLatticeField2D<T, NSDESCRIPTOR, descriptors::POROSITY> liquid_frac(NSlattice);
  liquid_frac.getName() = "liquid fraction";

  AnalyticalFfromCallableF<MyCase::d, T, T> melt_solution([&](Vector<T, 2> x) -> Vector<T, 1> {
    T X = 2.0 * k * util::sqrt(lambda / density / cp_l * iT * physDeltaT);
    if (util::nearZero(x[0] - X)) {
      return 0.5;
    }
    else if (x[0] < X) {
      return 1.0;
    }
    else if (x[0] > X) {
      return 0.0;
    }
    throw std::invalid_argument("ERROR!");
  });

  SuperAbsoluteErrorL2Norm2D<T> absMeltFractionErrorNormL2(liquid_frac, melt_solution, indicatorF);
  absMeltFractionErrorNormL2(result, input);
  params.set<parameters::SUM_ERROR_MELT>(params.get<parameters::SUM_ERROR_MELT>() + result[0]);
  clout << "melt-fraction-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm2D<T> relMeltFractionErrorNormL2(liquid_frac, melt_solution, indicatorF);
  relMeltFractionErrorNormL2(result, input);
  clout << "; melt-fraction-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm2D<T> absMeltFractionErrorNormLinf(liquid_frac, melt_solution, indicatorF);
  absMeltFractionErrorNormLinf(result, input);
  params.set<parameters::SUM_ERROR_MELT_INF>(params.get<parameters::SUM_ERROR_MELT_INF>() + result[0]);
  clout << "melt-fraction-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm2D<T> relMeltFractionErrorNormLinf(liquid_frac, melt_solution, indicatorF);
  relMeltFractionErrorNormLinf(result, input);
  clout << "; melt-fraction-Linf-error(rel)=" << result[0] << std::endl;

  SuperLatticeField2D<T, TDESCRIPTOR, descriptors::TEMPERATURE> temperature(ADlattice);
  temperature.getName() = "temperature";

  AnalyticalFfromCallableF<MyCase::d, T, T> temp_solution([&](Vector<T, 2> x) -> Vector<T,1>{
    T X = 2.0 * k * util::sqrt(lambda / density / cp_l * iT * physDeltaT);
    if (x[0] < X) {
      return Thot - (Thot - Tcold) / erf(k) * erf(x[0] * k / X);
    }
    if (x[0] >= X) {
      return Tcold;
    }
    throw std::invalid_argument("ERROR!");
  });

  SuperAbsoluteErrorL2Norm2D<T> absTemperatureErrorNormL2(temperature, temp_solution, indicatorF);
  absTemperatureErrorNormL2(result, input);
  params.set<parameters::SUM_ERROR_TEMP>(params.get<parameters::SUM_ERROR_TEMP>() + result[0]);
  clout << "temperature-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm2D<T> relTemperatureErrorNormL2(temperature, temp_solution, indicatorF);
  relTemperatureErrorNormL2(result, input);
  clout << "; temperature-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm2D<T> absTemperatureErrorNormLinf(temperature, temp_solution, indicatorF);
  absTemperatureErrorNormLinf(result, input);
  params.set<parameters::SUM_ERROR_TEMP_INF>(params.get<parameters::SUM_ERROR_TEMP_INF>() + result[0]);
  clout << "temperature-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm2D<T> relTemperatureErrorNormLinf(temperature, temp_solution, indicatorF);
  relTemperatureErrorNormLinf(result, input);
  clout << "; temperature-Linf-error(rel)=" << result[0] << std::endl;

  params.set<parameters::NUM_ERRORS>(params.get<parameters::NUM_ERRORS>() + 1);
}

void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;

  auto& params   = myCase.getParameters();
  auto& geometry = myCase.getGeometry();

  const T ly         = params.get<parameters::DOMAIN_EXTENT>()[1];
  const T physDeltaX = params.get<parameters::DOMAIN_EXTENT>()[0] / params.get<parameters::RESOLUTION>();

  geometry.rename(0, 2);
  geometry.rename(2, 1, {1, 0});
  std::vector<T> extend(2, T(0));
  extend[0] = physDeltaX;
  extend[1] = ly + 2. * physDeltaX;
  std::vector<T> origin(2, T(0));
  origin[0] = -physDeltaX;
  IndicatorCuboid2D<T> left(extend, origin);
  geometry.rename(2, 3, 1, left);

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
  using TotalEnthalpyAdvectionDiffusionDynamics = TotalEnthalpyAdvectionDiffusionBGKdynamics<T, TDESCRIPTOR>;

  auto& params    = myCase.getParameters();
  auto& geometry  = myCase.getGeometry();
  auto& NSlattice = myCase.getLattice(NavierStokes {});
  auto& ADlattice = myCase.getLattice(Temperature {});

  const T N          = params.get<parameters::RESOLUTION>();
  const T tau        = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T lx         = params.get<parameters::DOMAIN_EXTENT>()[0];
  const T cp_s       = params.get<parameters::CP_S>();
  const T cp_l       = params.get<parameters::CP_L>();
  const T Thot       = params.get<parameters::T_HOT>();
  const T Tcold      = params.get<parameters::T_COLD>();
  const T Ste        = params.get<parameters::STEFAN>();
  const T lambda     = params.get<parameters::PHYS_LAMBDA>();
  const T density    = params.get<parameters::PHYS_DENSITY>();
  const T k          = params.get<parameters::K>();
  const T cp_ref     = 2.0 * cp_s * cp_l / (cp_s + cp_l);
  const T L          = cp_l * (Thot - Tcold) / Ste;
  const T physDeltaX = lx / N;
  const T physDeltaT =
      density * cp_ref / lambda / descriptors::invCs2<T, NSDESCRIPTOR>() * (tau - 0.5) * physDeltaX * physDeltaX;

  NSlattice.setUnitConverter<ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR>>(
      (T)physDeltaX,              // physDeltaX
      (T)physDeltaT,              // physDeltaT
      (T)lx,                      // charPhysLength
      (T)1.0,                     // charPhysVelocity
      (T)lambda / density / cp_l, // physViscosity
      (T)density,                 // physDensity
      (T)lambda,                  // physThermalConductivity
      (T)cp_l,                    // physSpecificHeatCapacity
      (T)1.0,                     // physThermalExpansionCoefficient
      (T)Tcold,                   // charPhysLowTemperature
      (T)Thot                     // charPhysHighTemperature
  );
  auto& converter = NSlattice.getUnitConverter();
  ADlattice.setUnitConverter(converter);
  converter.print();

  const T lattice_Hcold = cp_s * Tcold;
  const T lattice_Hhot  = cp_l * Thot;

  clout << "H_cold: " << lattice_Hcold << std::endl;
  clout << "H_hot: " << lattice_Hhot << std::endl;
  clout << "k: " << std::setprecision(17) << k << std::endl;
  clout << "Ste: " << Ste << std::endl;
  clout << "lattice cp: " << converter.getLatticeSpecificHeatCapacity(cp_l) << std::endl;

  auto& coupling = myCase.setCouplingOperator("Boussinesq", TotalEnthalpyPhaseChangeCoupling {}, names::NavierStokes {},
                                              NSlattice, names::Temperature {}, ADlattice);
  coupling.restrictTo(geometry.getMaterialIndicator({1}));

  T Tomega  = converter.getLatticeThermalRelaxationFrequency();
  T NSomega = converter.getLatticeRelaxationFrequency();

  dynamics::set<ForcedPSMBGKdynamics>(NSlattice, geometry.getMaterialIndicator({1, 2, 3}));
  dynamics::set<TotalEnthalpyAdvectionDiffusionDynamics>(ADlattice, geometry.getMaterialIndicator({1, 3}));
  boundary::set<boundary::BounceBack>(ADlattice, geometry, 2);

  boundary::set<boundary::LocalVelocity>(NSlattice, geometry, 2);
  boundary::set<boundary::LocalVelocity>(NSlattice, geometry, 3);

  boundary::set<boundary::RegularizedTemperature>(ADlattice, geometry.getMaterialIndicator(3));

  NSlattice.setParameter<descriptors::OMEGA>(NSomega);

  ADlattice.setParameter<descriptors::OMEGA>(Tomega);
  ADlattice.setParameter<TotalEnthalpy::T_S>(Tcold);
  ADlattice.setParameter<TotalEnthalpy::T_L>(Tcold);
  ADlattice.setParameter<TotalEnthalpy::CP_S>(cp_s);
  ADlattice.setParameter<TotalEnthalpy::CP_L>(cp_l);
  ADlattice.setParameter<TotalEnthalpy::LAMBDA_S>(cp_ref / descriptors::invCs2<T, TDESCRIPTOR>() * (tau - 0.5));
  ADlattice.setParameter<TotalEnthalpy::LAMBDA_L>(cp_ref / descriptors::invCs2<T, TDESCRIPTOR>() * (tau - 0.5));
  ADlattice.setParameter<TotalEnthalpy::L>(L);

  std::vector<T> dir {0.0, 1.0};
  std::vector<T> forcePrefactor {0, 0};
  T boussinesqForcePrefactor = 9.81 / converter.getConversionFactorVelocity() * converter.getConversionFactorTime() *
                               converter.getCharPhysTemperatureDifference() *
                               converter.getPhysThermalExpansionCoefficient();

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

  coupling.template setParameter<TotalEnthalpyPhaseChangeCoupling::T_S>(Tcold);
  coupling.template setParameter<TotalEnthalpyPhaseChangeCoupling::T_L>(Tcold);
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

  using T = MyCase::value_t;

  auto& params    = myCase.getParameters();
  auto& geometry  = myCase.getGeometry();
  auto& NSlattice = myCase.getLattice(NavierStokes {});
  auto& ADlattice = myCase.getLattice(Temperature {});
  auto& converter = NSlattice.getUnitConverter();

  const T Re    = params.get<parameters::REYNOLDS>();
  const T Thot  = params.get<parameters::T_HOT>();
  const T Tcold = params.get<parameters::T_COLD>();
  const T cp_s  = params.get<parameters::CP_S>();
  const T cp_l  = params.get<parameters::CP_L>();

  const T lattice_Hcold = cp_s * Tcold;
  const T lattice_Hhot  = cp_l * Thot;

  std::vector<T>          zero(2, T());
  AnalyticalConst2D<T, T> u(zero);
  AnalyticalConst2D<T, T> rho(1.);
  AnalyticalConst2D<T, T> force(zero);

  T u_Re = converter.getLatticeVelocity(Re * converter.getPhysViscosity() / converter.getCharPhysLength());

  AnalyticalConst2D<T, T> u_top(converter.getCharLatticeVelocity(), u_Re);
  AnalyticalConst2D<T, T> u_bot(0.0, u_Re);
  T                       omega = converter.getLatticeRelaxationFrequency();
  AnalyticalConst2D<T, T> omegaField(omega);
  fields::set<descriptors::OMEGA>(NSlattice, geometry.getMaterialIndicator({1,2,3,4}), omegaField);

  AnalyticalConst2D<T, T> Cold(lattice_Hcold);
  AnalyticalConst2D<T, T> Hot(lattice_Hhot);

  fields::set<descriptors::VELOCITY>(ADlattice, geometry.getMaterialIndicator({1,2,3}), u);

  momenta::setTemperature(ADlattice, geometry.getMaterialIndicator({1,2}), Tcold);
  momenta::setTemperature(ADlattice, geometry.getMaterialIndicator({3}), Thot);
  ADlattice.iniEquilibrium(geometry.getMaterialIndicator({1, 2}), Cold, u);
  ADlattice.iniEquilibrium(geometry, 3, Hot, u);

  NSlattice.initialize();
  ADlattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void getResults(MyCase& myCase, util::Timer<MyCase::value_t>& timer, int iT)
{
  OstreamManager clout(std::cout, "getResults");

  using T            = MyCase::value_t;
  using NSDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using TDESCRIPTOR  = MyCase::descriptor_t_of<Temperature>;

  auto& params    = myCase.getParameters();
  auto& NSlattice = myCase.getLattice(NavierStokes {});
  auto& ADlattice = myCase.getLattice(Temperature {});
  auto& converter = NSlattice.getUnitConverter();

  const T   Thot              = params.get<parameters::T_HOT>();
  const T   Tcold             = params.get<parameters::T_COLD>();
  const T   lambda            = params.get<parameters::PHYS_LAMBDA>();
  const T   density           = params.get<parameters::PHYS_DENSITY>();
  const T   cp_l              = params.get<parameters::CP_L>();
  const T   cp_s              = params.get<parameters::CP_S>();
  const T   tau               = params.get<parameters::LATTICE_RELAXATION_TIME>();
  const T   physDeltaX        = params.get<parameters::DOMAIN_EXTENT>()[0] / params.get<parameters::RESOLUTION>();
  const int num_error_samples = params.get<parameters::NUM_ERROR_SAMPLES>();
  const T   k                 = params.get<parameters::K>();
  const T   cp_ref            = 2.0 * cp_s * cp_l / (cp_s + cp_l); // J /kg K
  const T   physDeltaT =
      density * cp_ref / lambda / descriptors::invCs2<T, NSDESCRIPTOR>() * (tau - 0.5) * physDeltaX * physDeltaX;

  const T maxPhysT = util::pow(0.25 / 2.0 / k, 2) * density * cp_l / lambda;

  SuperVTMwriter2D<T>                                      vtkWriter("stefanMelting2d");
  SuperLatticePhysVelocity2D<T, NSDESCRIPTOR>              velocity(NSlattice, converter);
  SuperLatticePhysPressure2D<T, NSDESCRIPTOR>              pressure(NSlattice, converter);
  SuperLatticePhysHeatFlux2D<T, NSDESCRIPTOR, TDESCRIPTOR> heatflux(ADlattice, converter);
  SuperLatticeDensity2D<T, TDESCRIPTOR>                    enthalpy(ADlattice);
  enthalpy.getName() = "enthalpy";
  SuperLatticeField2D<T, NSDESCRIPTOR, descriptors::POROSITY> liquid_frac(NSlattice);
  liquid_frac.getName() = "liquid fraction";
  SuperLatticeField2D<T, TDESCRIPTOR, descriptors::TEMPERATURE> temperature(ADlattice);
  temperature.getName() = "temperature";

  AnalyticalFfromCallableF<MyCase::d, T, T> melt_solution([&](Vector<T, 2> x) -> Vector<T, 1> {
    T X = 2.0 * k * util::sqrt(lambda / density / cp_l * iT * physDeltaT);
    if (util::nearZero(x[0] - X)) {
      return 0.5;
    }
    else if (x[0] < X) {
      return 1.0;
    }
    else if (x[0] > X) {
      return 0.0;
    }
    throw std::invalid_argument("ERROR!");
  });
  melt_solution.getName() = "Analytical liquid fraction";

  SuperLatticeFfromAnalyticalF2D<T, TDESCRIPTOR> melt_solution_lattice(melt_solution, ADlattice);

  AnalyticalFfromCallableF<MyCase::d, T, T> temp_solution([&](Vector<T, 2> x) -> Vector<T, 1> {
    T X = 2.0 * k * util::sqrt(lambda / density / cp_l * iT * physDeltaT);
    if (x[0] < X) {
      return Thot - (Thot - Tcold) / std::erf(k) * std::erf(x[0] * k / X);
    }
    if (x[0] >= X) {
      return Tcold;
    }
    throw std::invalid_argument("ERROR!");
  });
  temp_solution.getName() = "Analytical temperature";

  SuperLatticeFfromAnalyticalF2D<T, TDESCRIPTOR> temp_solution_lattice(temp_solution, ADlattice);

  vtkWriter.addFunctor(pressure);
  vtkWriter.addFunctor(velocity);
  vtkWriter.addFunctor(enthalpy);
  vtkWriter.addFunctor(temperature);
  vtkWriter.addFunctor(heatflux);
  vtkWriter.addFunctor(liquid_frac);
  vtkWriter.addFunctor(melt_solution_lattice);
  vtkWriter.addFunctor(temp_solution_lattice);

  const int vtkIter = converter.getLatticeTime(maxPhysT / T(num_error_samples));

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
    ADlattice.getStatistics().print(iT, converter.getPhysTime(iT));
    timer.print(iT);
    getError(myCase, iT);

    vtkWriter.write(iT);
  }
}

void simulate(MyCase& myCase)
{
  OstreamManager clout(std::cout, "getResults");
  using T = MyCase::value_t;

  auto& params    = myCase.getParameters();
  auto& geometry  = myCase.getGeometry();
  auto& NSlattice = myCase.getLattice(NavierStokes {});
  auto& ADlattice = myCase.getLattice(Temperature {});
  auto& converter = NSlattice.getUnitConverter();

  const T lambda   = params.get<parameters::PHYS_LAMBDA>();
  const T density  = params.get<parameters::PHYS_DENSITY>();
  const T cp_l     = params.get<parameters::CP_L>();
  const T k        = params.get<parameters::K>();
  const T maxPhysT = util::pow(0.25 / 2.0 / k, 2) * density * cp_l / lambda;

  util::Timer<T> timer(converter.getLatticeTime(maxPhysT), geometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT) + 1; ++iT) {
    myCase.getOperator("Boussinesq").apply();
    std::vector<T>          zero(2, T());
    AnalyticalConst2D<T, T> u(zero);
    ADlattice.defineField<descriptors::VELOCITY>(geometry, 1, u);
    getResults(myCase, timer, iT);

    // NSlattice.collideAndStream();
    ADlattice.collideAndStream();
  }
  const int num_errors         = params.get<parameters::NUM_ERRORS>();
  const T   sum_error_melt     = params.get<parameters::SUM_ERROR_MELT>();
  const T   sum_error_temp     = params.get<parameters::SUM_ERROR_TEMP>();
  const T   sum_error_melt_inf = params.get<parameters::SUM_ERROR_MELT_INF>();
  const T   sum_error_temp_inf = params.get<parameters::SUM_ERROR_TEMP_INF>();

  clout << std::setprecision(17);
  clout << "liquid_fraction-L2-error(abs,mean)=" << sum_error_melt / num_errors << std::endl;
  clout << "temperature-L2-error(abs,mean)=" << sum_error_temp / num_errors << std::endl;
  clout << "liquid_fraction-Linf-error(abs,mean)=" << sum_error_melt_inf / num_errors << std::endl;
  clout << "temperature-Linf-error(abs,mean)=" << sum_error_temp_inf / num_errors << std::endl;

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
    myCaseParameters.set<RESOLUTION>((int)128);
    myCaseParameters.set<DOMAIN_EXTENT>({1.0, 1 / 8.});
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(1.);
    myCaseParameters.set<REYNOLDS>(20.);
    myCaseParameters.set<STEFAN>(0.01);
    myCaseParameters.set<T_COLD>(0.5);
    myCaseParameters.set<T_HOT>(1.5);
    myCaseParameters.set<PHYS_LAMBDA>(1. / 6.);
    myCaseParameters.set<CP_S>(1.0);
    myCaseParameters.set<CP_L>(1.0);
    myCaseParameters.set<PHYS_DENSITY>(1.0);
    myCaseParameters.set<NUM_ERROR_SAMPLES>((int)25);
    myCaseParameters.set<STARTING_VALUE>(1.e-5);
    myCaseParameters.set<ITERATIONS>(100);
    myCaseParameters.set<MAX_ERROR>(2e-7);
    myCaseParameters.set<K>([] (MyCase::ParametersD& p) -> T  {
      int       i             = 0;
      T       x             = p.get<STARTING_VALUE>();
      int     maxIterations = p.get<ITERATIONS>();
      const T maxError      = p.get<MAX_ERROR>();
      const T Ste           = p.get<STEFAN>();
      T       frac;
      do {
        frac = (Ste / util::exp(x * x) / std::erf(x) - x * util::sqrt(M_PI)) /
               (-(2 * Ste * util::exp(-x * x) * x / erf(x)) -
                (2 * Ste * util::exp(-2 * x * x) / (util::sqrt(M_PI) * (std::erf(x) * std::erf(x)))) - util::sqrt(M_PI));
        x = x - frac;
        i++;
      } while (i < maxIterations && util::abs(frac) > maxError);
      return x;
    });
    myCaseParameters.set<PHYS_MAX_T>([] (MyCase::ParametersD& p) {
      return util::pow(0.25 / 2.0 / p.get<K>(), 2) * p.get<PHYS_DENSITY>() * p.get<CP_L>() / p.get<PHYS_LAMBDA>();
    });
  }

  myCaseParameters.set<parameters::SUM_ERROR_MELT>(0);
  myCaseParameters.set<parameters::SUM_ERROR_TEMP>(0);
  myCaseParameters.set<parameters::SUM_ERROR_MELT_INF>(0);
  myCaseParameters.set<parameters::SUM_ERROR_TEMP_INF>(0);

  Mesh mesh = createMesh(myCaseParameters);

  MyCase myCase(myCaseParameters, mesh);

  prepareGeometry(myCase);

  prepareLattice(myCase);

  setInitialValues(myCase);

  simulate(myCase);
}
