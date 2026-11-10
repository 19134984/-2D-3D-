/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2014-2018 Albert Mink
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

/* sphere3d.cpp:
 * A 3D implementation of radiation being emitted from a sphere and
 * spreading in a bigger sphere through highly scattering media. The macroscopic
 * method is based on the diffusion approximation.
 * The theoretical background and validation of said method are detailed in
 * [A. Mink, G. Thäter, H. Nirschl, and M. J. Krause. “A 3D Lattice Boltzmann method
 * for light simulation in participating media”. In: Journal of Computational Science
 * 17.Part 2 (2016). Discrete Simulation of Fluid Dynamics 2015, pp. 431–437. DOI:
 * 10.1016/j.jocs.2016.03.014.]
 */

#include <olb.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <format>

using namespace olb;
using namespace olb::names;

using MyCase = Case<Radiation, Lattice<double, descriptors::D3Q7<descriptors::VELOCITY>>>;

namespace olb::parameters {

struct WRITE_SEC : public descriptors::FIELD_BASE<1> {};

} // namespace olb::parameters

// PostProcessor for updating the source term in each cell
struct ComputeSourceTerm {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct C_T : public descriptors::FIELD_BASE<1> {};
  struct MU_EFF : public descriptors::FIELD_BASE<1> {};

  using parameters = meta::list<C_T, MU_EFF>;

  int getPriority() const { return 0; }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform
  {
    using V  = typename CELL::value_t;
    V phi    = cell.computeRho();
    V mu_eff = parameters.template get<MU_EFF>();
    V source = -mu_eff * mu_eff * phi * parameters.template get<C_T>();
    cell.template setField<descriptors::SOURCE>(source);
  }
};

// Compute the analytical solution and error
void calculate_error(MyCase& myCase)
{
  OstreamManager clout(std::cout, "error");

  using T        = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params   = myCase.getParameters();

  auto& Rlattice = myCase.getLattice(Radiation {});

  using RDESCRIPTOR = MyCase::descriptor_t_of<Radiation>;

  const T absorption = params.get<parameters::ABSORPTION>();
  const T scattering = params.get<parameters::SCATTERING>();

  int input[1];
  T   normAnaSol[1], absErr[1], relErr[1];
  // analytical solution
  PLSsolution3D<T, T, RDESCRIPTOR>               dSol(absorption, scattering);
  SuperLatticeFfromAnalyticalF3D<T, RDESCRIPTOR> dSolLattice(dSol, Rlattice);
  // approximated solution
  SuperLatticeDensity3D<T, RDESCRIPTOR> d(Rlattice);

  SuperL2Norm3D<T> dL2Norm(dSolLattice - d, geometry, 1);
  SuperL2Norm3D<T> dSolL2Norm(dSolLattice, geometry, 1);
  dL2Norm(absErr, input);
  dSolL2Norm(normAnaSol, input);
  relErr[0] = absErr[0] / normAnaSol[0];
  clout << "density-L2-error(util::abs)=" << absErr[0] << "; density-L2-error(rel)=" << relErr[0] << std::endl;

  SuperLinfNorm3D<T> dLinfNorm(dSolLattice - d, geometry, 1);
  dLinfNorm(absErr, input);
  clout << "density-Linf-error(util::abs)=" << absErr[0] << std::endl;
}

Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& params)
{
  using T = MyCase::value_t;

  const auto lx = params.get<parameters::DOMAIN_L>();
  const auto dx = params.get<parameters::PHYS_DELTA_X>();

  Vector<T, 3>         origin;
  IndicatorSphere3D<T> sphereAll(origin, 1 * lx);

  Mesh<T, MyCase::d> mesh(sphereAll, dx, 2 * singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());

  return mesh;
}

void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepraringGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T        = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params   = myCase.getParameters();

  const T lx = params.get<parameters::DOMAIN_L>();

  Vector<T, 3> origin;

  //Domain, Sphere of Size lx
  IndicatorSphere3D<T> sphereAll(origin, 1 * lx);
  geometry.rename(0, 2, sphereAll);
  geometry.rename(2, 1, {1, 1, 1});

  //Emitter Spehere of 10% domain:l
  IndicatorSphere3D<T> sphereEm(origin, 0.1 * lx);
  geometry.rename(1, 3, sphereEm);

  geometry.clean();
  geometry.innerClean();
  geometry.checkForErrors();

  geometry.print();
}

void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T        = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params   = myCase.getParameters();

  auto& Rlattice = myCase.getLattice(Radiation {});

  using RDESCRIPTOR = MyCase::descriptor_t_of<Radiation>;

  Rlattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, RDESCRIPTOR>>(
      (int)params.get<parameters::RESOLUTION>(), // resolution: number of voxels per charPhysL
      (T)params.get<
          parameters::LATTICE_RELAXATION_TIME>(), // latticeRelaxationTime: relaxation time, has to be greater than 0.5!
      (T)1.0,                                     // charPhysLength: reference length of simulation geometry
      (T)0.0, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__  || none for radiation
      (T)1.0, // physViscosity: physical kinematic viscosity in __m^2 / s__  || radiative diffusion coefficient D
      (T)1.0  // physDensity: physical density in __kg / m^3__
  );

  const auto& converter = Rlattice.getUnitConverter();
  converter.print();

  const T omega = converter.getLatticeRelaxationFrequency();

  //define dynamics
  // sourced advection diffusion dynamics OR poisson dynamics with sink term
  // when using poisson dynamics, setting the sink parameter is necessary
  // when using sourced advection diffusion dynamics, the sink term is handeled by the post processor
  auto bulkIndicator = geometry.getMaterialIndicator({1});
  //dynamics::set<PoissonDynamics<T,DESCRIPTOR>>(Rlattice, bulkIndicator);
  dynamics::set<SourcedAdvectionDiffusionBGKdynamics<T, RDESCRIPTOR>>(Rlattice, bulkIndicator);

  // Material=2,3 -->first order equilibrium
  dynamics::set<EquilibriumBoundaryFirstOrder>(Rlattice, geometry, 2);
  dynamics::set<EquilibriumBoundaryFirstOrder>(Rlattice, geometry, 3);

  // set the parameters for the dynamics
  Rlattice.setParameter<descriptors::OMEGA>(omega);
  // Rlattice.setParameter<collision::Poisson::SINK>(latticeSink);

  // Add PostProcessor and set parameters
  Rlattice.addPostProcessor<stage::PreCollide>(meta::id<ComputeSourceTerm> {});
  Rlattice.setParameter<ComputeSourceTerm::C_T>(converter.getConversionFactorTime());
  Rlattice.setParameter<ComputeSourceTerm::MU_EFF>(params.get<parameters::MU_EFF>());

  clout << "Prepare Lattice ... OK" << std::endl;
  return;
}

void setInitialValues(MyCase& myCase)
{
  OstreamManager clout(std::cout, "setInitialValues");
  clout << "Setting Initial Values and BC ..." << std::endl;

  using T        = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params   = myCase.getParameters();

  auto& Rlattice = myCase.getLattice(Radiation {});

  using RDESCRIPTOR = MyCase::descriptor_t_of<Radiation>;

  const T absorption = params.get<parameters::ABSORPTION>();
  const T scattering = params.get<parameters::SCATTERING>();

  // since adv-diffusion model is used, the velocity it set to 0
  AnalyticalConst3D<T, T>          u0(0.0, 0.0, 0.0);            // 3D -> 3D
  AnalyticalConst3D<T, T>          rho0(0.0);                    // 3D -> 3D
  PLSsolution3D<T, T, RDESCRIPTOR> dSol(absorption, scattering); // analytical Solution

  // initialize media with density from analytical solution
  // at iT=0 the error is given by the maschinen genauigkeit
  Rlattice.iniEquilibrium(geometry, 1, rho0, u0);
  Rlattice.iniEquilibrium(geometry, 2, dSol, u0);
  Rlattice.iniEquilibrium(geometry, 3, dSol, u0);

  momenta::setDensity(Rlattice, geometry.getMaterialIndicator({2}), dSol);
  momenta::setDensity(Rlattice, geometry.getMaterialIndicator({3}), dSol);

  Rlattice.initialize();

  clout << "Setting Initial Values and BC ... OK" << std::endl;
}

void getResults(MyCase& myCase, util::Timer<MyCase::value_t>& timer, std::size_t iT)
{
  OstreamManager clout(std::cout, "getResults");

  using T      = MyCase::value_t;
  auto& params = myCase.getParameters();

  auto&       Rlattice  = myCase.getLattice(Radiation {});
  const auto& converter = Rlattice.getUnitConverter();

  using RDESCRIPTOR = MyCase::descriptor_t_of<Radiation>;

  const T absorption = params.get<parameters::ABSORPTION>();
  const T scattering = params.get<parameters::SCATTERING>();

  const auto saveIter = converter.getLatticeTime(params.get<parameters::WRITE_SEC>() / 32.);

  SuperVTMwriter3D<T>                   vtmWriter("sphere3d");
  SuperLatticeDensity3D<T, RDESCRIPTOR> density(Rlattice);
  vtmWriter.addFunctor(density);

  if (iT == 0) {
    SuperLatticeCuboid3D cuboid(Rlattice);
    SuperLatticeRank3D   rank(Rlattice);

    vtmWriter.write(cuboid);
    vtmWriter.write(rank);

    vtmWriter.createMasterFile();
  }

  if (iT % saveIter == 0) {
    std::cout << "Write added functor at time " << iT << "." << std::endl;
    vtmWriter.write(iT);

    timer.update(iT);
    timer.printStep();
    Rlattice.getStatistics().print(iT, converter.getPhysTime(iT));
  }

  if (iT % 500 == 0) {
    SuperLatticeDensity3D<T, RDESCRIPTOR> density(Rlattice);
    AnalyticalFfromSuperF3D<T>            analyDensity(density, true, 1);
    PLSsolution3D<T, T, RDESCRIPTOR>      dSol(absorption, scattering);

    for (int nY = -25; nY < 25; ++nY) {
      double position[3] = {0, ((double)nY) / 25.0, 0.0};

      double densityEval[1] = {0};
      analyDensity(densityEval, position);
      double densitySol[1] = {0};
      dSol(densitySol, position);

      // printf("%1.3f\t%1.3f\t%1.3f\n", position[1], densityEval[0], densitySol[0] );
      // fout << postition[1] << ", " << densityEval[0] << ", " << densitySol[0]  << std::endl;
    }
  }
}

void writeToFile(MyCase& myCase, std::size_t iT)
{
  OstreamManager clout(std::cout, "writeToFile");
  using T      = MyCase::value_t;
  auto& params = myCase.getParameters();

  auto&       Rlattice  = myCase.getLattice(Radiation {});
  const auto& converter = Rlattice.getUnitConverter();

  using RDESCRIPTOR = MyCase::descriptor_t_of<Radiation>;

  const T absorption = params.get<parameters::ABSORPTION>();
  const T scattering = params.get<parameters::SCATTERING>();

  std::string   name = std::format("Results_along_line_at_time_{}.txt", std::to_string(iT));
  std::ofstream fout(name, std::ios::trunc);

  if (!fout) {
    clout << "Error: could not open " << name << std::endl;
  }

  SuperLatticeDensity3D<T, RDESCRIPTOR> density(Rlattice);
  AnalyticalFfromSuperF3D<T>            analytDensity(density, true, 1);
  PLSsolution3D<T, T, RDESCRIPTOR>      dSol(absorption, scattering);

  int n = converter.getResolution() / 2;

  for (int nY = 0; nY < n; ++nY) {
    double position[3]    = {0, double(nY) / n, 0};
    double densityEval[1] = {0};
    analytDensity(densityEval, position);
    double densitySol[1] = {0};
    dSol(densitySol, position);
    if (singleton::mpi().getRank() == 0) {
      fout << position[1] << ", " << densityEval[0] << ", " << densitySol[0] << ", "
           << util::abs(densityEval[0] - densitySol[0]) / densitySol[0] << std::endl;
      printf("%1.3f\t%1.3f\t%1.3f\t%1.4f\n", position[1], densityEval[0], densitySol[0],
             util::abs(densityEval[0] - densitySol[0]) / densitySol[0]);
    }
  }

  fout.close();
}

void simulate(MyCase& myCase)
{
  OstreamManager clout(std::cout, "simulate");
  clout << "Starting Simulation ..." << std::endl;

  using T        = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& params   = myCase.getParameters();

  auto&       Rlattice  = myCase.getLattice(Radiation {});
  const auto& converter = Rlattice.getUnitConverter();

  const auto maxPhysT = params.get<parameters::PHYS_MAX_T>();
  const auto saveIter = converter.getLatticeTime(params.get<parameters::WRITE_SEC>() / 32.);

  util::Timer<T> timer(converter.getLatticeTime(maxPhysT), geometry.getStatistics().getNvoxel());
  timer.start();

  util::ValueTracer<T> converge(params.get<parameters::RESOLUTION>(), 1e-6);

  std::size_t iT;
  for (iT = 0; iT <= converter.getLatticeTime(maxPhysT); ++iT) {
    converge.takeValue(Rlattice.getStatistics().getAverageRho(), true);

    if (converge.hasConverged()) {
      clout << "Simulation converged. -- " << iT << std::endl;
      break;
    }

    if (iT % saveIter == 0) {
      Rlattice.getStatistics().print(iT, converter.getPhysTime(iT));
      timer.print(iT);
      calculate_error(myCase);
    }
    Rlattice.collideAndStream();
    getResults(myCase, timer, iT);
  }

  writeToFile(myCase, converter.getPhysTime(iT));

  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<WRITE_SEC>(1.0);
    myCaseParameters.set<RESOLUTION>(40);
    myCaseParameters.set<PHYS_MAX_T>(100);
    myCaseParameters.set<DOMAIN_L>(1.0);
    myCaseParameters.set<ABSORPTION>(0.5);
    myCaseParameters.set<SCATTERING>(1.5);

    myCaseParameters.set<LATTICE_RELAXATION_TIME>(1.0);
    myCaseParameters.set<MU_EFF>([] (MyCase::ParametersD& p) {
      return util::sqrt(3 * p.get<ABSORPTION>() *
                        (p.get<ABSORPTION>() + p.get<SCATTERING>()));
    });

    myCaseParameters.set<PHYS_DELTA_X>([] (MyCase::ParametersD& p) {
      return p.get<DOMAIN_L>() / p.get<RESOLUTION>();
    });
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
