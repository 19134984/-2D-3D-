/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2026 Dominik Ebert, Arsh Kumbhat, Stephan Simonis,
 *                     Fedor Bukreev, Adrian Kummerländer, Shota Ito
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

/* airfoilOpti2d.cpp:
 * This optimization example examines laminar flow past a NACA0012 airfoil placed in a
 * large domain. The drag (CD) and lift (CL) coefficient are accurate upto
 * Re = 1000 for angles of attack upto 10 degrees. Cp and Cf have also
 * been validated for Re = 1000 at 10 degrees AOA.
 * The objective function (J) contains the mean of force coefficients in a
 * specified averaging intervall [t_0, t_max]. The averaged Drag coefficient is minimized
 * or the Lift-to-Drag ratio is maximized according to the angle of attack.
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::opti;

using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D2Q9<descriptors::VELOCITY>>
>;

using MyOptiCase = OptiCaseFDQ<
  Controlled, MyCase
>;

enum class ControlType : int {
  AOA = 0,
  AOA_WITH_CHAMBER = 1
};

enum class ObjectiveType : int {
  LIFT_DRAG_RATIO = 0,
  DRAG = 1
};

namespace olb::parameters {
  struct ANGLE_OF_ATTACK              : public descriptors::FIELD_BASE<1> {};
  struct ANGLE_OF_ATTACK_IN_DEGREES   : public descriptors::FIELD_BASE<1> {};
  struct CHORD_LENGTH                 : public descriptors::FIELD_BASE<1> {};
  struct MAX_CAMBER                   : public descriptors::FIELD_BASE<1> {};
  struct CENTER_AIRFOIL               : public descriptors::FIELD_BASE<0,1> {};
  struct CST_COEFFICIENTS_TOP         : public descriptors::FIELD_BASE<8> {};
  struct CST_COEFFICIENTS_BOTTOM      : public descriptors::FIELD_BASE<8> {};
  struct PHYS_START_AVERAGE_T         : public descriptors::FIELD_BASE<1> {};

  struct CONTROL_TYPE                 : public descriptors::TYPED_FIELD_BASE<ControlType,1> {};
  struct OBJECTIVE_TYPE               : public descriptors::TYPED_FIELD_BASE<ObjectiveType,1> {};

  struct AVERAGE_CD                   : public descriptors::FIELD_BASE<1> {};
  struct AVERAGE_CL                   : public descriptors::FIELD_BASE<1> {};
  struct AVERAGE_COUNT                : public descriptors::TYPED_FIELD_BASE<unsigned, 1> {};

  struct FDQ_STENCIL_WIDTH_IN_DEGREES : public descriptors::FIELD_BASE<1> {};
}

struct SmoothInflowUpdateO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct VELOCITY : public descriptors::FIELD_BASE<0,1> { };

  using parameters = meta::list<VELOCITY>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    auto& cell = cells.template get<names::NavierStokes>();
    auto u = parameters.template get<VELOCITY>();
    cell.defineU(u.data());
  }
};

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const Vector origin{0, 0};
  IndicatorCuboid2D<T> cuboid(extent, origin);

  const T physDeltaX = parameters.get<parameters::PHYS_DELTA_X>();
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  geometry.setWriteIncrementalVTK(true);

  Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  Vector origin{0., 0.};
  const T L = parameters.get<parameters::PHYS_DELTA_X>();

  geometry.rename(0, 2);
  geometry.rename(2, 1, {1, 1});

  // Set material number for inflow
  extent[0] = 2. * L;
  origin[0] = -L;
  IndicatorCuboid2D<T> inflow(extent, origin);
  geometry.rename(2, 3, 1, inflow);

  // Set material number for outflow
  origin[0] = parameters.get<parameters::DOMAIN_EXTENT>()[0] - L;
  IndicatorCuboid2D<T> outflow(extent, origin);
  geometry.rename(2, 4, 1, outflow);

  // Set material number for airfoil
  const Vector center = parameters.get<parameters::CENTER_AIRFOIL>();
  const T angleOfAttack = parameters.get<parameters::ANGLE_OF_ATTACK>();
  IndicatorParametricAirfoil2D<T> airfoilI(center,
                                           parameters.get<parameters::CHORD_LENGTH>(),
                                           util::toStdVector(parameters.get<parameters::CST_COEFFICIENTS_TOP>()),
                                           util::toStdVector(parameters.get<parameters::CST_COEFFICIENTS_BOTTOM>()),
                                           angleOfAttack,
                                           parameters.get<parameters::MAX_CAMBER>(),
                                           0.4);

  IndicatorRotate<T,2> airfoil(airfoilI.getAerodynamicCenter(),
                               angleOfAttack,
                               airfoilI);
  geometry.rename(1, 5, airfoil);

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  geometry.checkForErrors();
  geometry.print();
}

void prepareLattice(MyCase& myCase)
{
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto bulkIndicator = geometry.getMaterialIndicator({1});

  /// Set up a unit converter with the characteristic physical units
  lattice.setUnitConverter<UnitConverterFromResolutionAndLatticeVelocity<T,DESCRIPTOR>>(
    parameters.get<parameters::RESOLUTION>(),            // resolution: number of voxels
    parameters.get<parameters::MACH>() / sqrt(3),        // latticeU
    parameters.get<parameters::CHORD_LENGTH>(),          // charPhysLength
    parameters.get<parameters::PHYS_CHAR_VELOCITY>(),    // charPhysVelocity
    parameters.get<parameters::PHYS_CHAR_VELOCITY>() *   // charPhysViscosity
    parameters.get<parameters::CHORD_LENGTH>() /
    parameters.get<parameters::REYNOLDS>(),
    parameters.get<parameters::PHYS_CHAR_DENSITY>()      // physDensity
  );
  lattice.getUnitConverter().print();

  // Material=1 --> Bulk dynamics
  lattice.defineDynamics<BGKdynamics>(bulkIndicator);

  // Material=2 --> Full Slip
  boundary::set<boundary::FullSlip>(lattice, geometry, 2);

  // Setting of the boundary conditions
  boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 3);
  boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 4);

  const Vector center = parameters.get<parameters::CENTER_AIRFOIL>();
  const T angleOfAttack = parameters.get<parameters::ANGLE_OF_ATTACK>();
  IndicatorParametricAirfoil2D<T> airfoilI(center,
                                           parameters.get<parameters::CHORD_LENGTH>(),
                                           util::toStdVector(parameters.get<parameters::CST_COEFFICIENTS_TOP>()),
                                           util::toStdVector(parameters.get<parameters::CST_COEFFICIENTS_BOTTOM>()),
                                           angleOfAttack,
                                           parameters.get<parameters::MAX_CAMBER>(),
                                           0.4);

  IndicatorRotate<T,2> airfoil(airfoilI.getAerodynamicCenter(),
                               angleOfAttack,
                               airfoilI);

  setBouzidiBoundary(lattice, geometry, 5, airfoil);

  // Initial conditions
  AnalyticalConst2D<T, T> rhoF(1);
  AnalyticalConst2D<T, T> rho0(0);
  std::vector<T>          velocity(2, T(0));
  AnalyticalConst2D<T, T> uF(velocity);

  // Initialize all values of distribution functions to their local equilibrium
  lattice.defineRhoU(geometry.getMaterialIndicator({1,2,3,4,5}), rhoF, uF);
  lattice.iniEquilibrium(geometry.getMaterialIndicator({1,2,3,4,5}), rhoF, uF);
  lattice.defineField<descriptors::VELOCITY>(geometry.getMaterialIndicator({1,2,3,4,5}), uF);
  lattice.defineField<descriptors::MATERIAL>(geometry.getMaterialIndicator({1,2,3,4}), rhoF);
  lattice.defineField<descriptors::MATERIAL>(geometry, 5, rho0);

  lattice.setParameter<descriptors::OMEGA>(lattice.getUnitConverter().getLatticeRelaxationFrequency());
  lattice.initialize();

  // Inflow velocity ramp-up
  auto& smoothInflowUpdateO = myCase.setCouplingOperator("SmoothInflow",
                                                         SmoothInflowUpdateO{},
                                                         names::NavierStokes{}, lattice);
  smoothInflowUpdateO.restrictTo(geometry.getMaterialIndicator(3));
}

void setTemporalValues(MyCase& myCase, std::size_t iT) {
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();
  auto& parameters = myCase.getParameters();
  auto& smoothInflowO = dynamic_cast<SuperLatticeCoupling<SmoothInflowUpdateO,
                                                          meta::map<NavierStokes, descriptors::VALUED_DESCRIPTOR<T,MyCase::descriptor_t>>
                                                         >&>(myCase.getOperator("SmoothInflow"));

  const std::size_t iTmaxStart =
    converter.getLatticeTime(parameters.get<parameters::PHYS_START_PERIOD>());
  const std::size_t iTupdate = parameters.get<parameters::PHYS_BOUNDARY_VALUE_UPDATE_T>();
  if (iT % iTupdate == 0 && iT <= iTmaxStart) {
    const T tau = (T)iT / (T)iTmaxStart;
    T smooth_factor = tau * tau * tau * (10.0 - 15.0 * tau + 6.0 * tau * tau);
    smoothInflowO.setParameter<SmoothInflowUpdateO::VELOCITY>({smooth_factor * converter.getCharLatticeVelocity(), 0});
    smoothInflowO.apply();
  }
}

void trackAveragedCoefficients(MyCase& myCase, std::size_t iT) {
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& converter = lattice.getUnitConverter();

  const T startAveragePhysT = parameters.get<parameters::PHYS_START_AVERAGE_T>();
  if (iT >= (converter.getLatticeTime(startAveragePhysT) - 1) && iT > 0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    // Drag and Lift via Momentum Exchange
    SuperLatticePhysDrag2D drag(lattice, geometry, 5, converter, parameters.get<parameters::CHORD_LENGTH>());

    int input[3] = {};
    T   _drag[drag.getTargetDim()];
    drag(_drag, input);

    const unsigned count = parameters.get<parameters::AVERAGE_COUNT>() + 1;
    parameters.set<parameters::AVERAGE_COUNT>(count);
    const T CD_avg = parameters.get<parameters::AVERAGE_CD>();
    const T CL_avg = parameters.get<parameters::AVERAGE_CL>();
    parameters.set<parameters::AVERAGE_CD>(CD_avg + (_drag[0] - CD_avg) / (count * 1.0));
    parameters.set<parameters::AVERAGE_CL>(CL_avg + (_drag[1] - CL_avg) / (count * 1.0));
  }
}

void getResults(MyCase& myCase, util::Timer<MyCase::value_t>& timer, std::size_t iT) {
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();
  const std::size_t iTlog = converter.getLatticeTime(parameters.get<parameters::PHYS_STAT_ITER_T>());

  /// Print some (numerical and computational) statistics
  if (iT%iTlog == 0) {
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
    timer.print(iT);

    // Compute CD / CL
    trackAveragedCoefficients(myCase, iT);
  }
}

void simulate(MyCase& myCase)
{
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();

  parameters.set<parameters::AVERAGE_CD>(0.);
  parameters.set<parameters::AVERAGE_CL>(0.);
  parameters.set<parameters::AVERAGE_COUNT>(0.);

  const T physMaxT = parameters.get<parameters::PHYS_MAX_T>();

  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(physMaxT);
  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    setTemporalValues(myCase, iT);
    myCase.getLattice(NavierStokes{}).collideAndStream();
    getResults(myCase, timer, iT);
  }

  timer.stop();
  timer.printSummary();
}

void setInitialControl(MyOptiCase& optiCase) {
  using T = MyOptiCase::value_t;
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& parameters = controlledCase.getParameters();
  OstreamManager clout(std::cout, "Control");

  switch (parameters.get<parameters::CONTROL_TYPE>()) {
    case ControlType::AOA: {
      const T aoa = parameters.get<parameters::ANGLE_OF_ATTACK>();
      clout << "Control angle of attack." << std::endl;
      optiCase.getController().set({aoa});
      break;
    }
    case ControlType::AOA_WITH_CHAMBER:
    default: {
      const T aoa = parameters.get<parameters::ANGLE_OF_ATTACK>();
      const T chamber = parameters.get<parameters::MAX_CAMBER>();
      clout << "Control angle of attack & max chamber." << std::endl;
      optiCase.getController().set({aoa, chamber});
      break;
    }
  }
}

void applyControl(MyOptiCase& optiCase) {
  auto& control = optiCase.getController();
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& parameters = controlledCase.getParameters();

  switch (parameters.get<parameters::CONTROL_TYPE>()) {
    case ControlType::AOA: {
      parameters.set<parameters::ANGLE_OF_ATTACK>(control[0]);
      break;
    }
    case ControlType::AOA_WITH_CHAMBER:
    default: {
      parameters.set<parameters::ANGLE_OF_ATTACK>(control[0]);
      parameters.set<parameters::MAX_CAMBER>(control[1]);
      break;
    }
  }
}

MyOptiCase::value_t objectiveF(MyOptiCase& optiCase) {
  using T = MyOptiCase::value_t;
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& parameters = controlledCase.getParameters();
  OstreamManager clout(std::cout, "Objective");
  T objective = 0.0;

  switch (parameters.get<parameters::OBJECTIVE_TYPE>()) {
    case ObjectiveType::LIFT_DRAG_RATIO:
    default: {
      objective = -parameters.get<parameters::AVERAGE_CL>() /
                   parameters.get<parameters::AVERAGE_CD>();
      clout << "Lift-to-drag ratio: " << objective << std::endl;
      break;
    }
    case ObjectiveType::DRAG: {
      objective = parameters.get<parameters::AVERAGE_CD>();
      clout << "Drag: " << objective << std::endl;
      break;
    }
  }
  return objective;
}

void getOptiResults(MyOptiCase& optiCase) {
  using T = MyCase::value_t;
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& lattice = controlledCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();
  std::size_t iM = optiCase.getObjectiveCount();

  SuperVTMwriter2D<T> vtmWriter("airfoilOpti2d");
  SuperLatticePhysVelocity2D velocityF(lattice, converter);
  SuperLatticePhysPressure2D pressureF(lattice, converter);
  SuperLatticeField2D<T,MyCase::descriptor_t,descriptors::MATERIAL> materialF(lattice);
  vtmWriter.addFunctor(velocityF);
  vtmWriter.addFunctor(pressureF);
  vtmWriter.addFunctor(materialF);

  if (iM == 0) {
    vtmWriter.createMasterFile();
  }

  if (iM > 0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    vtmWriter.write(iM);
  }
}

MyOptiCase::value_t computeObjective(MyOptiCase& optiCase) {
  auto& controlledCase = optiCase.getCase(Controlled{});
  controlledCase.resetGeometry();
  controlledCase.resetLattices();

  applyControl(optiCase);
  prepareGeometry(controlledCase);
  prepareLattice(controlledCase);
  simulate(controlledCase);
  getOptiResults(optiCase);
  return objectiveF(optiCase);
}

/// Steps B-G: Set up and run the optimization routine
int main(int argc, char* argv[])
{
  using T = MyCase::value_t;
  // Step 2: Initialization
  initialize(&argc, &argv);

  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION                  >(        60); // 120
    myCaseParameters.set<REYNOLDS                    >(      1000);
    myCaseParameters.set<PHYS_MAX_T                  >(      200.);
    myCaseParameters.set<PHYS_START_PERIOD           >([](MyCase::ParametersD& p){
      return p.get<PHYS_MAX_T>() * 0.4;
    });
    myCaseParameters.set<PHYS_START_AVERAGE_T        >([](MyCase::ParametersD& p){
      return p.get<PHYS_MAX_T>() * 0.5;
    });
    myCaseParameters.set<PHYS_STAT_ITER_T            >(       0.5);
    myCaseParameters.set<PHYS_BOUNDARY_VALUE_UPDATE_T>(         5);
    myCaseParameters.set<ANGLE_OF_ATTACK_IN_DEGREES  >(       25.);
    myCaseParameters.set<ANGLE_OF_ATTACK             >([](MyCase::ParametersD& p){
      return p.get<ANGLE_OF_ATTACK_IN_DEGREES>() *
             M_PI / 180.;
    });
    myCaseParameters.set<MAX_CAMBER                  >(        0.);
    myCaseParameters.set<MACH                        >(     0.069);
    myCaseParameters.set<PHYS_CHAR_VELOCITY          >(       0.2);
    myCaseParameters.set<CHORD_LENGTH                >(        1.);
    myCaseParameters.set<PHYS_DELTA_X                >([](MyCase::ParametersD& p){
      return p.get<CHORD_LENGTH>() /
             p.get<RESOLUTION>();
    });
    myCaseParameters.set<DOMAIN_EXTENT               >([](MyCase::ParametersD& p){
      return Vector{16., 4. + p.get<PHYS_DELTA_X>()};
    });
    myCaseParameters.set<CENTER_AIRFOIL              >([](MyCase::ParametersD& p){
      return Vector{6., p.get<DOMAIN_EXTENT>()[1] / 2.};
    });
    myCaseParameters.set<PHYS_CHAR_DENSITY           >(       1.0);
    myCaseParameters.set<CST_COEFFICIENTS_TOP        >(
      {0.17461638362047,0.144580464119799,0.187466427348086,0.11360028706947,
       0.170833407038172,0.12860181567567,0.141970597786503,0.141375015320354});
    myCaseParameters.set<CST_COEFFICIENTS_BOTTOM     >([](MyCase::ParametersD& p){
      return p.get<CST_COEFFICIENTS_TOP>() * -1.;
    });
    myCaseParameters.set<CONTROL_TYPE                >(ControlType::AOA_WITH_CHAMBER);
    myCaseParameters.set<OBJECTIVE_TYPE              >(ObjectiveType::LIFT_DRAG_RATIO);
    myCaseParameters.set<FDQ_STENCIL_WIDTH_IN_DEGREES>(0.4);
  }

  OptimizerParametersD<MyCase::value_t> optiParameters;
  {
    using namespace parameters;
    optiParameters.set<optimizer::EPS>(1e-4);
    optiParameters.set<optimizer::MAX_OPTIMIZATION_STEPS>(20);
    optiParameters.set<linesearch::INITIAL_STEP_SIZE>(0.01);
    optiParameters.set<optimizer::CONTROL_EPS>(0.);
    optiParameters.set<linesearch::MAX_STEP_ATTEMPTS>(10);
    optiParameters.set<linesearch::STEP_CONDITION>("Wolfe");
    optiParameters.set<optimizer::VERBOSE>(true);
    optiParameters.set<optimizer::LOG_FILE_NAME>("log");
  }

  // Step 3: Create Mesh
  Mesh mesh = createMesh(myCaseParameters);

  // Step 4: Create Case
  MyCase myCase(myCaseParameters, mesh);

  MyOptiCase optiCase;
  optiCase.setCase<Controlled>(myCase);
  // Set the stencil width for the FDQ
  const T stencilWidth = myCaseParameters.get<parameters::FDQ_STENCIL_WIDTH_IN_DEGREES>() *
                         M_PI / 180.0;
  optiCase.setParameter<parameters::FDQ_STENCIL_WIDTH>(stencilWidth);

  setInitialControl(optiCase);

  optiCase.setObjective(computeObjective);

  OptimizerSteepestDescent<T,std::vector<T>> optimizer(
    optiCase.getController().size(),
    optiParameters,
    {OptimizerLogType::objective, OptimizerLogType::control, OptimizerLogType::derivative}
  );

  optimizer.optimize(optiCase);
}
