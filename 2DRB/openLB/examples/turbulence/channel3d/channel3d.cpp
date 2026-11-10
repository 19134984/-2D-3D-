/*  Lattice Boltzmann sample, written in C++, using the OpenLB library
 *
 *  Copyright (C) 2025 Fedor Bukreev, Adrian Kummerlaender, Yuji (Sam) Shimojima,
 *                     Jonathan Jeppener-Haltenhoff, Marc Haußmann,
 *                     Mathias J. Krause
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

/* channel3d.cpp:
 * This example examines a wall-bounded flow at high Reynolds numbers.
 * The dynamics follow the RLB collision operator of 3rd order and Smagorinsky-Lilly turbulence model.
 * The near wall region is modelled by a wall function.
 *
 * The example shows the usage of turbulent wall model with different parameters.
 * All WM paramertes are integrated into HLBM interface.
 * As showcase for stability improvement HRR collision is possible.
 *
 * The reference paper with setup and results description:
 * Myoungkyu Lee and Robert D. Moser,
 * Direct numerical simulation of turbulent channel flow up to Re_tau = 5200,
 * 2015, Journal of Fluid Mechanics, vol. 774, pp. 395-415
 *
 * The paper can be found at http://journals.cambridge.org/article_S0022112015002682
 * The data can be found in: https://turbulence.oden.utexas.edu/channel2015/data/
 *
 * The comparative plots can be generated using python script plotResults.py
 * The referenceData of Re_tau 1000 case are stored in referenceData.csv
 */

#include <olb.h>

#include "superForceTermApplyInChannel3D.h"

using namespace olb;
using namespace names;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<float, descriptors::D3Q19<descriptors::FORCE, descriptors::Y1, descriptors::U_TAU>>
>;

namespace olb::parameters {

struct WITH_WALL_MODEL : public descriptors::TYPED_FIELD_BASE<bool, 1> {};
struct HRR_COLLISION : public descriptors::TYPED_FIELD_BASE<bool, 1> {};

// physical simulated length adapted for lattice distance to boundary in meters
struct PHYS_ADAPTED_CHAR_LENGTH : public descriptors::FIELD_BASE<1> {};

// (lattice-)distance from one point away from the wall to the wall
struct LATTICE_DISTANCE_TO_WALL : public descriptors::FIELD_BASE<1> {};

//friction reynolds number ReTau
struct REYNOLDS_TAU : public descriptors::FIELD_BASE<1> {};

//domain size base; this vector filed will be multiplied by std::numbers::pi_v<T> * physRefL in the directions of the wall parallel and by physRefL in the wall normal direction.
struct DOMAIN_EXTENT_BASE : public descriptors::FIELD_BASE<0, 1> {};

/// Wallfunction parameters

//  Used method for density reconstruction
//  0: use local density
//  1: extrapolation (Guo)
//  2: constant (rho = 1.)
struct METHOD_OF_RHO : public descriptors::TYPED_FIELD_BASE<int, 1> {};

//  Used method for non-equilibrium population reconstruction
//  0: extrapolation NEQ (Guo Zhaoli)
//  1: first or second order finite differnce (Malaspinas)
//  2: equilibrium scheme (no fNeq)
struct METHOD_OF_FNEQ : public descriptors::TYPED_FIELD_BASE<int, 1> {};

//  Used wall profile
//  0: power law profile
//  1: Spalding profile
struct METHOD_WALL_FUNCTION_PROFILE : public descriptors::TYPED_FIELD_BASE<int, 1> {};

// check if descriptor with body force is used
struct ENABLE_BODY_FORCE : public descriptors::TYPED_FIELD_BASE<bool, 1> {};

// interpolate sampling velocity along given normal between lattice voxels
struct ENABLE_INTERPOLATE_SAMPLE_VELOCITY : public descriptors::TYPED_FIELD_BASE<bool, 1> {};

// use van Driest damping function for turbulent viscosity in boundary cell
struct ENABLE_VANDRIEST : public descriptors::TYPED_FIELD_BASE<bool, 1> {};

struct ENABLE_MOVING_WALL : public descriptors::TYPED_FIELD_BASE<bool, 1> {};

struct ENABLE_TURBULENCE_STATISTICS : public descriptors::TYPED_FIELD_BASE<bool, 1> {};

//  distance from cell to velocity sampling point in lattice units
struct LATTICE_DISTANCE_SAMPLING : public descriptors::FIELD_BASE<1> {};

// time until until statistics sampling in seconds. This value will be multiple by characteristic physical time in main function.
struct PHYS_CONVERGE_TIME : public descriptors::FIELD_BASE<1> {};

// statistics sampling time in seconds.
struct PHYS_STATISTICS_TIME : public descriptors::FIELD_BASE<1> {};

// This value will be multiple by characteristic physical time and stored in PHYS_CONVERGE_TIME in main function.
struct PHYS_CONVERGE_TIME_BASE : public descriptors::FIELD_BASE<1> {};

// This value will be multiple by characteristic physical time and stored in PHYS_STATISTICS_TIME in main function.
struct PHYS_STATISTICS_TIME_BASE : public descriptors::FIELD_BASE<1> {};

//HRR Collision paramertes
// very strong influence of numerical diffusion. 0.99 meand 1% FDM strain rate tensor, which is already enough!
struct CONST_HRR_HYBRID : public descriptors::FIELD_BASE<1> {};

} // namespace olb::parameters

// seed the rng with time if SEED_WITH_TIME is set, otherwise just use a fixed seed.
#if defined SEED_WITH_TIME
#include <chrono>
auto                       seed = std::chrono::system_clock().now().time_since_epoch().count();
std::default_random_engine generator(seed);
#else
std::default_random_engine generator(0x1337533DAAAAAAAA);
#endif

template <typename T, typename S>
class Channel3D : public AnalyticalF3D<T, S> {
protected:
  T turbulenceIntensity;
  T maxVelocity;
  T distanceToWall;
  T obst_z;
  T obst_r;
  T a;
  T b;

public:
  Channel3D(MyCase& myCase, T frac)
      : AnalyticalF3D<T, S>(3)
  {
    const auto converter = myCase.getLattice(NavierStokes {}).getUnitConverter();
    const auto physRefL = myCase.getParameters().get<parameters::PHYS_CHAR_LENGTH>();

    turbulenceIntensity  = 0.05;
    maxVelocity          = converter.getLatticeVelocity(converter.getCharPhysVelocity() * (8.0 / 7.0)); // Centerline Velocity
    obst_r               = physRefL;
    a                    = -1.;
    b                    = 1.;
  };

  bool operator()(T output[], const S input[])
  {
    std::uniform_real_distribution<BaseType<T>> distribution(a, b);
    T                                           nRandom1 = distribution(generator);
    T                                           nRandom2 = distribution(generator);
    T                                           nRandom3 = distribution(generator);

    if ((util::abs(input[2] - obst_r)) < obst_r) {
      T u_calc = maxVelocity * util::pow(((obst_r - util::abs(input[2] - obst_r)) / obst_r), 1. / 7.);

      output[0] = turbulenceIntensity * nRandom1 * maxVelocity + u_calc;
      output[1] = turbulenceIntensity * nRandom2 * maxVelocity;
      output[2] = turbulenceIntensity * nRandom3 * maxVelocity;
    }
    return true;
  };
};

/// @brief Create a simulation mesh, based on user-specific geometry
/// @return An instance of Mesh, which keeps the relevant information
Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters)
{
  using T = MyCase::value_t;
  using namespace parameters;
  const T              physDeltaX                 = parameters.get<parameters::PHYS_DELTA_X>();
  const auto           lx                         = parameters.get<DOMAIN_EXTENT>()[0];
  const T              ly                         = parameters.get<DOMAIN_EXTENT>()[1];
  const auto           adaptedPhysSimulatedLength = parameters.get<PHYS_ADAPTED_CHAR_LENGTH>();
  const auto           latticeWallDistance        = parameters.get<LATTICE_DISTANCE_TO_WALL>();
  Vector<T, 3>         extend(lx, ly, adaptedPhysSimulatedLength);
  Vector<T, 3>         origin((T)0.0, (T)0.0, -((T)1.0 - latticeWallDistance) * physDeltaX);
  IndicatorCuboid3D<T> cuboid(extend, origin);

  Mesh<T, MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.getCuboidDecomposition().setPeriodicity({true, true, false});
  mesh.setOverlap(parameters.get<OVERLAP>());
  return mesh;
}

/// @brief Set material numbers for different parts of the domain
/// @param myCase The Case instance which keeps the simulation data
/// @note The material numbers will be used to assign physics to lattice nodes
void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;
  using namespace parameters;
  auto& parameters                      = myCase.getParameters();
  using T                               = MyCase::value_t;
  auto&      sGeometry                  = myCase.getGeometry();
  const auto physDeltaX                 = parameters.get<PHYS_DELTA_X>(); // lattice spacing
  const auto lx                         = parameters.get<DOMAIN_EXTENT>()[0];
  const T    ly                         = parameters.get<DOMAIN_EXTENT>()[1];
  const auto adaptedPhysSimulatedLength = parameters.get<PHYS_ADAPTED_CHAR_LENGTH>();
  const auto latticeWallDistance        = parameters.get<LATTICE_DISTANCE_TO_WALL>();

  Vector<T, 3>         extend(lx, ly, adaptedPhysSimulatedLength);
  Vector<T, 3>         origin((T)0.0, (T)0.0, -((T)1.0 - latticeWallDistance) * physDeltaX);
  IndicatorCuboid3D<T> cuboid(extend, origin);

  sGeometry.rename(0, 2, cuboid);
  sGeometry.rename(2, 1, {0, 0, 1});
  sGeometry.clean();
  sGeometry.innerClean();
  sGeometry.checkForErrors();

  sGeometry.print();

  olb::Vector<T, 3> PhyMax = sGeometry.getStatistics().getMaxPhysR(2);
  olb::Vector<T, 3> PhyMin = sGeometry.getStatistics().getMinPhysR(2);
  clout << "Dimension of the channel in meters: x = " << PhyMax[0] - PhyMin[0];
  clout << " ; y = " << PhyMax[1] - PhyMin[1];
  clout << " ; z = " << PhyMax[2] - PhyMin[2] << std::endl;

  clout << "Prepare Geometry ... OK" << std::endl;
}

/// Set initial condition for primal variables (velocity and density)
/// @param myCase The Case instance which keeps the simulation data
/// @note Be careful: initial values have to be set using lattice units
void setInitialConditions(MyCase& myCase)
{
  OstreamManager clout(std::cout, "setInitialConditions");
  clout << "Set initial conditions ..." << std::endl;
  using T         = MyCase::value_t;
  auto& sGeometry = myCase.getGeometry();
  auto& sLattice  = myCase.getLattice(NavierStokes {});
  auto& parameters = myCase.getParameters();

  AnalyticalConst3D<T, T> rho(1.0);
  AnalyticalConst3D<T, T> rho0(0.0);
  AnalyticalConst3D<T, T> u0(0.0, 0.0, 0.0);
  Channel3D<T, T>         uSol(myCase, 1.0);

  momenta::setVelocity(sLattice, sGeometry.getMaterialIndicator(1), uSol);
  sLattice.iniEquilibrium(sGeometry, 1, rho, uSol);

  momenta::setVelocity(sLattice, sGeometry.getMaterialIndicator(2), uSol);
  sLattice.iniEquilibrium(sGeometry, 2, rho, uSol);

  sLattice.defineField<descriptors::VELOCITY2>(sGeometry.getMaterialIndicator({0, 1, 2}), uSol);
  sLattice.defineField<descriptors::AVERAGE_VELOCITY>(sGeometry.getMaterialIndicator({0, 1, 2}), u0);
  sLattice.defineField<descriptors::AVERAGE_PRESSURE>(sGeometry.getMaterialIndicator({0, 1, 2}), rho0);
  sLattice.defineField<descriptors::AVERAGE_SQUARE_PRESSURE>(sGeometry.getMaterialIndicator({0, 1, 2}), rho0);
  sLattice.defineField<descriptors::FORCE>(sGeometry.getMaterialIndicator({0, 1, 2}), u0);
  if(parameters.get<parameters::WITH_WALL_MODEL>() || parameters.get<parameters::HRR_COLLISION>()) {
    sLattice.defineField<descriptors::POROSITY>(sGeometry.getMaterialIndicator({0, 2}), rho0);
    sLattice.defineField<descriptors::POROSITY>(sGeometry, 1, rho);
  }
  sLattice.setParameter<descriptors::OMEGA>(sLattice.getUnitConverter().getLatticeRelaxationFrequency());
  sLattice.setParameter<collision::LES::SMAGORINSKY>(T(0.2));

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Set initial conditions ... OK" << std::endl;
}

/// @brief Set lattice dynamics
/// @param myCase The Case instance which keeps the simulation data
void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;
  using T          = MyCase::value_t;
  auto& sGeometry  = myCase.getGeometry();
  auto& sLattice   = myCase.getLattice(NavierStokes {});
  auto& parameters = myCase.getParameters();
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

  sLattice.setUnitConverter<UnitConverter<T, DESCRIPTOR>>(
      // cell size
      parameters.get<parameters::PHYS_DELTA_X>(),
      // time step size
      parameters.get<parameters::PHYS_DELTA_T>(),
      // charPhysLength: reference length of simulation geometry
      T(2) * parameters.get<parameters::PHYS_CHAR_LENGTH>(),
      // charPhysVelocity: characteristic mean bulk velocity in __m / s__
      parameters.get<parameters::PHYS_CHAR_VELOCITY>(),
      // charphysViscosity:characteristic physical kinematic viscosity in __m^2 / s__
      parameters.get<parameters::PHYS_CHAR_VISCOSITY>(),
      // physDensity:characteristic physical density in __kg / m^3__
      parameters.get<parameters::PHYS_CHAR_DENSITY>());
  sLattice.getUnitConverter().print();

  const auto physConvergeTime           = parameters.get<parameters::PHYS_CONVERGE_TIME>();
  const auto maxPhysT                   = parameters.get<parameters::PHYS_MAX_T>();
  const auto adaptedPhysSimulatedLength = parameters.get<parameters::PHYS_ADAPTED_CHAR_LENGTH>();
  const auto ReTau                      = parameters.get<parameters::REYNOLDS_TAU>();
  const auto physRefL                   = parameters.get<parameters::PHYS_CHAR_LENGTH>();
  const auto N                          = parameters.get<parameters::RESOLUTION>();
  const auto latticeWallDistance        = parameters.get<parameters::LATTICE_DISTANCE_TO_WALL>();

  clout << "----------------------------------------------------------------------" << std::endl;
  clout << "Converge time(s): " << physConvergeTime << std::endl;
  clout << "Lattice converge time: " << sLattice.getUnitConverter().getLatticeTime(physConvergeTime) << std::endl;
  clout << "Max. Phys. simulation time(s): " << maxPhysT << std::endl;
  clout << "Max. Lattice simulation time: " << sLattice.getUnitConverter().getLatticeTime(maxPhysT) << std::endl;
  clout << "----------------------------------------------------------------------" << std::endl;
  clout << "Channel height(m): " << adaptedPhysSimulatedLength << std::endl;
  clout << "y+ value: "
        << (ReTau * sLattice.getUnitConverter().getPhysViscosity() / (physRefL)) *
               ((2.0 / (T)(N + 2 * latticeWallDistance)) * latticeWallDistance) /
               sLattice.getUnitConverter().getPhysViscosity()
        << std::endl;
  clout << "y+ value spacing: "
        << (ReTau * sLattice.getUnitConverter().getPhysViscosity() / (physRefL)) *
               (sLattice.getUnitConverter().getPhysDeltaX()) / sLattice.getUnitConverter().getPhysViscosity()
        << std::endl;
  clout << "----------------------------------------------------------------------" << std::endl;

  //sLattice.defineDynamics<BounceBack>(sGeometry, 2);
  const T              physDeltaX                 = parameters.get<parameters::PHYS_DELTA_X>();
  const T              lx                         = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T              ly                         = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  Vector<T, 3>         extend(lx + T(6) * physDeltaX, ly + T(6) * physDeltaX, T(2) * physRefL);
  Vector<T, 3>         origin(-T(3) * physDeltaX, -T(3) * physDeltaX, T(0));
  IndicatorCuboid3D<T> channel(extend, origin);
  setBouzidiBoundary<T,DESCRIPTOR>(sLattice, sGeometry, 2, channel);
  if(parameters.get<parameters::WITH_WALL_MODEL>()) {
    WallModelParameters<T> wallModelParameters;
    {
      using namespace olb::parameters;
      wallModelParameters.bodyForce                 = parameters.get<ENABLE_BODY_FORCE>();
      wallModelParameters.rhoMethod                 = parameters.get<METHOD_OF_RHO>();
      wallModelParameters.fNeqMethod                = parameters.get<METHOD_OF_FNEQ>();
      wallModelParameters.samplingCellDistance      = parameters.get<LATTICE_DISTANCE_SAMPLING>();
      wallModelParameters.interpolateSampleVelocity = parameters.get<ENABLE_INTERPOLATE_SAMPLE_VELOCITY>();
      wallModelParameters.useVanDriest              = parameters.get<ENABLE_VANDRIEST>();
      wallModelParameters.wallFunctionProfile       = parameters.get<METHOD_WALL_FUNCTION_PROFILE>();
      wallModelParameters.latticeWallDistance       = parameters.get<LATTICE_DISTANCE_TO_WALL>();
      wallModelParameters.movingWall                = parameters.get<ENABLE_MOVING_WALL>();
      wallModelParameters.turbulenceStatistics      = parameters.get<ENABLE_TURBULENCE_STATISTICS>();
    }
    /// Material=1 -->bulk dynamics
    setTurbulentWallModelDynamics(sLattice, sGeometry, 1, wallModelParameters);
    setTurbulentWallModel(sLattice, sGeometry, 2, wallModelParameters, &channel);
  }else{
    sLattice.defineDynamics<SmagorinskyForcedBGKdynamics<T, DESCRIPTOR>::template wrap_collision<collision::TrackTurbulenceStatistics>>(sGeometry, 1);
  }
  if(parameters.get<parameters::HRR_COLLISION>()) {
    const T hybridConst = parameters.get<parameters::CONST_HRR_HYBRID>();
    sLattice.addPostProcessor<stage::PostStream>(sGeometry.getMaterialIndicator({1}),
                                                 meta::id<FDMstrainRateTensorPostProcessor> {});
    AnalyticalConst3D<T, T> hybrid(hybridConst);
    sLattice.defineField<collision::HYBRID>(sGeometry, 1, hybrid);
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

void getResults(MyCase& myCase, std::size_t iT, util::Timer<MyCase::value_t>& timer)
{
  OstreamManager clout(std::cout, "getResults");
  using T                     = MyCase::value_t;
  using DESCRIPTOR            = MyCase::descriptor_t_of<NavierStokes>;
  auto&      sGeometry        = myCase.getGeometry();
  auto&      sLattice         = myCase.getLattice(NavierStokes {});
  auto&      parameters       = myCase.getParameters();
  const auto maxPhysT         = parameters.get<parameters::PHYS_MAX_T>();
  const auto ReTau            = parameters.get<parameters::REYNOLDS_TAU>();
  const auto charPhysNu       = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const auto physConvergeTime = parameters.get<parameters::PHYS_CONVERGE_TIME>();
  const T    checkstatistics  = (T)maxPhysT / 200.0;

  if (iT == 0) {
    // Writes the geometry for visualization
    SuperVTMwriter3D<T> vtmWriter("channel3d");
    vtmWriter.createMasterFile();
  }

  // Writes output on the console
  if (iT % sLattice.getUnitConverter().getLatticeTime(checkstatistics) == 0) {
    // Timer console output
    timer.update(iT);
    timer.printStep(2);
    // Lattice statistics console output
    sLattice.getStatistics().print(iT, iT * sLattice.getUnitConverter().getPhysDeltaT());
    clout << "Max. physical velocity(m/s): "
          << sLattice.getUnitConverter().getPhysVelocity(sLattice.getStatistics().getMaxU()) << std::endl;
    clout << "Max u+:"
          << sLattice.getUnitConverter().getPhysVelocity(sLattice.getStatistics().getMaxU()) /
                 (ReTau * charPhysNu / (sLattice.getUnitConverter().getCharPhysLength() / 2.))
          << std::endl;
  }

  std::size_t iTstartAvg = sLattice.getUnitConverter().getLatticeTime(physConvergeTime);
  if (iT == iTstartAvg) {
    SuperLatticeVelocity3D<T, DESCRIPTOR> latticeVelocity(sLattice);
    SuperLatticePressure3D<T, DESCRIPTOR> latticePressure(sLattice);
    SuperLatticeSquarePressure3D<T, DESCRIPTOR> latticeSquarePressure(sLattice);
    SuperLatticeSquareVelocity3D<T, DESCRIPTOR> latticeSquareVelocity(sLattice);
    SuperLatticeField3D<T, DESCRIPTOR, descriptors::U_TAU> uTauL(sLattice);
    sLattice.defineField<descriptors::AVERAGE_VELOCITY>(sGeometry.getMaterialIndicator({1}), latticeVelocity);
    sLattice.defineField<descriptors::AVERAGE_PRESSURE>(sGeometry.getMaterialIndicator({1}), latticePressure);
    sLattice.defineField<descriptors::AVERAGE_SQUARE_PRESSURE>(sGeometry.getMaterialIndicator({1}), latticeSquarePressure);
    sLattice.defineField<descriptors::AVERAGE_VELOCITY_X_VELOCITY>(sGeometry.getMaterialIndicator({1}), latticeSquareVelocity);
    sLattice.defineField<descriptors::AVERAGE_U_TAU>(sGeometry.getMaterialIndicator({1}), uTauL);
  }
  if (iT < iTstartAvg) {
    sLattice.setParameter<descriptors::LATTICE_TIME>(2);
  }
  else {
    sLattice.setParameter<descriptors::LATTICE_TIME>(iT - iTstartAvg + 1);
  }

  if (/*iT % sLattice.getUnitConverter().getLatticeTime(checkstatistics) == 0 ||*/
      iT == 0 ||
      iT == sLattice.getUnitConverter().getLatticeTime(maxPhysT) - 1) {
    // Writes the vtk files
    sLattice.executePostProcessors(stage::Evaluation {});
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    if(parameters.get<parameters::WITH_WALL_MODEL>()) {
      sLattice.scheduleBackgroundOutputVTK([&, iT](auto task) {
        SuperVTMwriter3D<T>        vtmWriter("channel3d");
        SuperLatticePhysVelocity3D velocity(sLattice, sLattice.getUnitConverter());
        SuperLatticeRMSPhysPressure3D rmsPressure(sLattice, sLattice.getUnitConverter());
        SuperLatticePhysPressure3D pressure(sLattice, sLattice.getUnitConverter());
        SuperLatticePhysField3D<T, DESCRIPTOR, descriptors::AVERAGE_VELOCITY> sAveragedVel(
            sLattice, sLattice.getUnitConverter().getConversionFactorVelocity());
        SuperLatticePhysField3D<T, DESCRIPTOR, descriptors::WMVELOCITY> wmvelocity(
              sLattice, sLattice.getUnitConverter().getConversionFactorVelocity());
        wmvelocity.getName() = "wmvelocity";
        SuperLatticeField3D<T, DESCRIPTOR, descriptors::Y1> y1(sLattice);
        y1.getName() = "y1";
        SuperGeometryF<T, 3> geom(sGeometry);
        vtmWriter.addFunctor(velocity);
        vtmWriter.addFunctor(pressure);
        vtmWriter.addFunctor(sAveragedVel);
        vtmWriter.addFunctor(geom);
        vtmWriter.addFunctor(wmvelocity);
        vtmWriter.addFunctor(rmsPressure);
        vtmWriter.addFunctor(y1);

        task(vtmWriter, iT);
      });
    } else {
      sLattice.scheduleBackgroundOutputVTK([&, iT](auto task) {
        SuperVTMwriter3D<T>        vtmWriter("channel3d");
        SuperLatticePhysVelocity3D velocity(sLattice, sLattice.getUnitConverter());
        SuperLatticeRMSPhysPressure3D rmsPressure(sLattice, sLattice.getUnitConverter());
        SuperLatticePhysPressure3D pressure(sLattice, sLattice.getUnitConverter());
        SuperLatticePhysField3D<T, DESCRIPTOR, descriptors::AVERAGE_VELOCITY> sAveragedVel(
            sLattice, sLattice.getUnitConverter().getConversionFactorVelocity());
        SuperLatticeField3D<T, DESCRIPTOR, descriptors::Y1> y1(sLattice);
        y1.getName() = "y1";
        SuperGeometryF<T, 3> geom(sGeometry);
        vtmWriter.addFunctor(velocity);
        vtmWriter.addFunctor(pressure);
        vtmWriter.addFunctor(sAveragedVel);
        vtmWriter.addFunctor(geom);
        vtmWriter.addFunctor(rmsPressure);
        vtmWriter.addFunctor(y1);

        task(vtmWriter, iT);
      });
    }
  }

  if (iT == sLattice.getUnitConverter().getLatticeTime(maxPhysT) - 1) {
    // Preparation
    const T lx                  = parameters.get<parameters::DOMAIN_EXTENT>()[0];
    const T ly                  = parameters.get<parameters::DOMAIN_EXTENT>()[1];
    const T physDeltaX          = parameters.get<parameters::PHYS_DELTA_X>();
    const T latticeWallDistance = parameters.get<parameters::LATTICE_DISTANCE_TO_WALL>();
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    SuperLatticeRMSPhysPressure3D<T, DESCRIPTOR>    rmsPressure(sLattice, sLattice.getUnitConverter());
    SuperLatticePhysField3D<T, DESCRIPTOR, descriptors::AVERAGE_VELOCITY> sAveragedVel(
        sLattice, sLattice.getUnitConverter().getConversionFactorVelocity());
    sAveragedVel.getName() = "sAveragedVel";
    SuperLatticePhysField3D<T, DESCRIPTOR, descriptors::AVERAGE_U_TAU> sAveragedUTau(
        sLattice, sLattice.getUnitConverter().getConversionFactorVelocity());
    sAveragedUTau.getName() = "sAveragedUTau";
    CSV<T> csvWriter("flowField", ';', {"y+", "u_tau", "uAv+", "uu++", "uv++", "uw++", "vv++", "vw++", "ww++", "pRMS"},
                     ".csv");
    int    dummy[1] = {};
    const int nZ    = sGeometry.getCuboidDecomposition().getMotherCuboid().getNz(); // number of voxels in x-direction
    auto      mat   = new SuperIndicatorIdentity3D<T>(sGeometry.getMaterialIndicator({1}));
    T         avUTau[sAveragedUTau.getTargetDim() + 1];
    SuperNonZeroAverage3D<T>(sAveragedUTau, sGeometry.getMaterialIndicator({1})).operator()(avUTau, dummy);

    if(avUTau[0] != T(0)){
      for (int iZ = 1; iZ < int(nZ/2); ++iZ) {
        const T z = sLattice.getUnitConverter().getPhysLength(iZ);
        auto cube = new SuperIndicatorFfromIndicatorF3D<T>(
          std::shared_ptr<IndicatorF3D<T>>(new IndicatorCuboid3D<T> ({T(0.8)*lx, T(0.8)*ly, T(0.25)*physDeltaX},{T(0.1)*lx, T(0.1)*ly, z-latticeWallDistance*physDeltaX})), sGeometry);
        SuperIndicatorMultiplication3D<T> plane(cube, mat);

        // Average concentration: use average functor
        T avRmsPres[rmsPressure.getTargetDim() + 1];
        SuperAverage3D<T>(rmsPressure, plane).operator()(avRmsPres, dummy);
        T avVel[sAveragedVel.getTargetDim() + 1];
        SuperAverage3D<T>(sAveragedVel, plane).operator()(avVel, dummy);

        csvWriter.writeDataFile(
            (z - latticeWallDistance * physDeltaX) * avUTau[0] / sLattice.getUnitConverter().getPhysViscosity(),
            {avUTau[0], avVel[0] / avUTau[0], avRmsPres[0]});
      }
    }
  }
}

/// @brief Execute simulation: set initial values and run time loop
/// @param myCase The Case instance which keeps the simulation data
void simulate(MyCase& myCase)
{

  OstreamManager clout(std::cout, "Time marching");
  using T               = MyCase::value_t;
  auto&      sGeometry  = myCase.getGeometry();
  auto&      sLattice   = myCase.getLattice(NavierStokes {});
  auto&      parameters = myCase.getParameters();
  auto&      converter  = sLattice.getUnitConverter();
  const auto maxPhysT   = parameters.get<parameters::PHYS_MAX_T>();
  const auto physRefL = parameters.get<parameters::PHYS_CHAR_LENGTH>();
  const auto ReTau                      = parameters.get<parameters::REYNOLDS_TAU>();
  const auto lx                         = sGeometry.getStatistics().getMaxPhysR(1)[0] -  sGeometry.getStatistics().getMinPhysR(1)[0];
  const auto ly                         = sGeometry.getStatistics().getMaxPhysR(1)[1] -  sGeometry.getStatistics().getMinPhysR(1)[1];
  const auto adaptedPhysSimulatedLength = parameters.get<parameters::PHYS_ADAPTED_CHAR_LENGTH>();


  //forcing of the channel
  SuperForceTermApplyInChannel3d forcingProcessor(sLattice, sGeometry, converter, ReTau, physRefL ,
               parameters.get<parameters::LATTICE_DISTANCE_TO_WALL>(), parameters.get<parameters::WITH_WALL_MODEL>());

  forcingProcessor.applyForcing();

  util::Timer<T> timer(converter.getLatticeTime(maxPhysT), sGeometry.getStatistics().getNvoxel());
  clout << "starting simulation..." << std::endl;
  timer.start();

  for (size_t iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {
    sLattice.setParameter<descriptors::LATTICE_TIME>(iT);
    getResults(myCase, iT, timer);
    forcingProcessor.applyForcing();
    sLattice.collideAndStream();
    sLattice.stripeOffDensityOffset(sLattice.getStatistics().getAverageRho()-(T)1);
  }
  sLattice.setProcessingContext(ProcessingContext::Evaluation);
  timer.stop();
  timer.printSummary();
}
int main(int argc, char* argv[])
{

  initialize(&argc, &argv);
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    using T = MyCase::value_t;
    using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
    myCaseParameters.set<WITH_WALL_MODEL>(false);
    myCaseParameters.set<HRR_COLLISION>(false);
    myCaseParameters.set<RESOLUTION>(40);
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1.0);
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.50025);
    myCaseParameters.set<REYNOLDS_TAU>(1000.512);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(5.0 / 100000.0);
    myCaseParameters.set<PHYS_CHAR_VELOCITY>(1.0);
    myCaseParameters.set<PHYS_CHAR_LENGTH>(1.0);
    myCaseParameters.set<LATTICE_DISTANCE_TO_WALL>(0.5);
    myCaseParameters.set<PHYS_ADAPTED_CHAR_LENGTH>([] (MyCase::ParametersD& p) {
      T   latticeWallDistance = p.get<LATTICE_DISTANCE_TO_WALL>();
      int N                   = p.get<RESOLUTION>();
      T   physRefL            = p.get<PHYS_CHAR_LENGTH>();
      return 2.0 * physRefL / (1.0 - 2.0 / N * (1.0 - latticeWallDistance));
    });

    myCaseParameters.set<PHYS_DELTA_X>([] (MyCase::ParametersD& p) {
      return p.get<PHYS_ADAPTED_CHAR_LENGTH>() / p.get<RESOLUTION>();
    });
    myCaseParameters.set<PHYS_DELTA_T>([] (MyCase::ParametersD& p) {
      return (p.get<LATTICE_RELAXATION_TIME>() - T(0.5)) / p.get<PHYS_CHAR_VISCOSITY>()
              / descriptors::invCs2<T,DESCRIPTOR>()
              * p.get<PHYS_DELTA_X>() * p.get<PHYS_DELTA_X>();
    });
    myCaseParameters.set<PHYS_CHAR_TIME>([] (MyCase::ParametersD& p) {
      T ReTau      = p.get<REYNOLDS_TAU>();
      T physRefL   = p.get<PHYS_CHAR_LENGTH>();
      T charPhysNu = p.get<PHYS_CHAR_VISCOSITY>();
      return physRefL / (ReTau * charPhysNu / physRefL);
    });

    myCaseParameters.set<PHYS_CONVERGE_TIME_BASE>(40.0);

    myCaseParameters.set<PHYS_CONVERGE_TIME>([] (MyCase::ParametersD& p) {
      T multiple  = p.get<PHYS_CONVERGE_TIME_BASE>();
      T charPhysT = p.get<PHYS_CHAR_TIME>();
      return charPhysT * multiple;
    });

    myCaseParameters.set<PHYS_STATISTICS_TIME_BASE>(150.0);
    myCaseParameters.set<PHYS_STATISTICS_TIME>([] (MyCase::ParametersD& p) {
      T multiple  = p.get<PHYS_STATISTICS_TIME_BASE>();
      T charPhysT = p.get<PHYS_CHAR_TIME>();
      return charPhysT * multiple;
    });

    myCaseParameters.set<PHYS_MAX_T>([] (MyCase::ParametersD& p) {
      return p.get<PHYS_CONVERGE_TIME>() + p.get<PHYS_STATISTICS_TIME>();
    });

    myCaseParameters.set<DOMAIN_EXTENT_BASE>({2.0, 2.0, 2.0});

    myCaseParameters.set<DOMAIN_EXTENT>([] (MyCase::ParametersD& p) -> Vector<T, 3>  {
      T physRefL = p.get<PHYS_CHAR_LENGTH>();
      return {p.get<DOMAIN_EXTENT_BASE>()[0] * std::numbers::pi_v<T> * physRefL,
              p.get<DOMAIN_EXTENT_BASE>()[1] * std::numbers::pi_v<T> * physRefL,
              p.get<DOMAIN_EXTENT_BASE>()[2] * physRefL};
    });

    /// Wallfunction parameters
    myCaseParameters.set<METHOD_OF_RHO>(0);
    myCaseParameters.set<METHOD_OF_FNEQ>(2);
    myCaseParameters.set<METHOD_WALL_FUNCTION_PROFILE>(0);
    myCaseParameters.set<ENABLE_BODY_FORCE>(true);
    myCaseParameters.set<ENABLE_INTERPOLATE_SAMPLE_VELOCITY>(false);
    myCaseParameters.set<ENABLE_VANDRIEST>(false);
    myCaseParameters.set<LATTICE_DISTANCE_SAMPLING>(3.5);
    myCaseParameters.set<OVERLAP>([] (MyCase::ParametersD& p) {
      return int(p.get<WITH_WALL_MODEL>()*p.get<LATTICE_DISTANCE_SAMPLING>() + T(1));
    });
    myCaseParameters.set<ENABLE_MOVING_WALL>(false);
    myCaseParameters.set<ENABLE_TURBULENCE_STATISTICS>(true);
  }

  {
    using namespace olb::parameters;
    if (!(myCaseParameters.get<METHOD_OF_RHO>() >= 0 && myCaseParameters.get<METHOD_OF_RHO>() <= 2)) {
      std::cerr << "Error: density reconstruction methods are 0 to 2" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (!(myCaseParameters.get<METHOD_OF_FNEQ>() >= 0 && myCaseParameters.get<METHOD_OF_FNEQ>() <= 2)) {
      std::cerr << "Error: non-equilibrium population reconstruction methods are 0 to 2" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (!(myCaseParameters.get<METHOD_WALL_FUNCTION_PROFILE>() >= 0 &&
          myCaseParameters.get<METHOD_WALL_FUNCTION_PROFILE>() <= 1)) {
      std::cerr << "Error:  wall function profile methods are 0 to 1" << std::endl;
      std::exit(EXIT_FAILURE);
    }
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
  setInitialConditions(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);

  return 0;
}
