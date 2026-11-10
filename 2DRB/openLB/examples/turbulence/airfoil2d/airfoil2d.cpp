/*  Lattice Boltzmann sample, written in C++, using the OpenLBlibrary
 *
 *  Copyright (C) 2025 Mathias J. Krause, David Heidenthal, Adrian Kummerländer, Michael Grinschewski, Fedor Bukreev
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

/* airfoil2d.cpp:
 * This example examines a steady flow past a NACA airfoil placed in a channel.
 * At the inlet, a random Turbulent profile is imposed on the velocity,
 * whereas the outlet implements a Dirichlet pressure condition set by p = 0.
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<float, descriptors::D2Q9<descriptors::POROSITY, descriptors::VELOCITY>>
>;

namespace olb::parameters {
struct CHORD_LENGTH : public descriptors::FIELD_BASE<1> {};
struct AIRFOIL_PARAMETERS : public descriptors::TYPED_FIELD_BASE<int,3> {};
struct ANGLE_OF_ATTACK : public descriptors::FIELD_BASE<1> {};
}

// Create mesh for simulation case
Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& params)
{
  using T                     = MyCase::value_t;
  Vector               extent = params.get<parameters::DOMAIN_EXTENT>();
  std::vector<T>       origin(2, T());
  IndicatorCuboid2D<T> cuboid(extent, origin);

  const T            physDeltaX = params.get<parameters::CHORD_LENGTH>() / params.get<parameters::RESOLUTION>();
  Mesh<T, MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  return mesh;
}

struct SmoothInflowUpdateO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct VELOCITY : public descriptors::FIELD_BASE<0,1> { };

  using parameters = meta::list<VELOCITY>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    auto& cell = cells.template get<names::NavierStokes>();
    auto  u    = parameters.template get<VELOCITY>();
    cell.defineU(u.data());
  }
};

void prepareGeometry(MyCase& myCase)
{
  using T        = MyCase::value_t;
  auto& params   = myCase.getParameters();
  auto& geometry = myCase.getGeometry();

  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  // Get relevant parameters for geometry creation
  T lengthX = params.get<parameters::DOMAIN_EXTENT>()[0];
  T lengthY = params.get<parameters::DOMAIN_EXTENT>()[1];

  T deltaX = params.get<parameters::CHORD_LENGTH>() / params.get<parameters::RESOLUTION>();

  Vector<T, 2> extend(lengthX, lengthY);
  Vector<T, 2> origin;

  // Get airfoil parameters
  Vector<T, 2>   airfoilCenter = params.get<parameters::CENTER>();
  Vector<int,3> airfoilParams = {params.get<parameters::AIRFOIL_PARAMETERS>()[0],
                                 params.get<parameters::AIRFOIL_PARAMETERS>()[1],
                                 params.get<parameters::AIRFOIL_PARAMETERS>()[2]};
  T              angleOfAttack = params.get<parameters::ANGLE_OF_ATTACK>();
  T              chordLength   = params.get<parameters::CHORD_LENGTH>();

  // Create airfoil indicator
  IndicatorAirfoil2D<T> airfoilI(airfoilCenter, chordLength, airfoilParams[0]/100., airfoilParams[1]/10., airfoilParams[2]/100.);
  IndicatorRotate<T, 2> airfoil(airfoilCenter, angleOfAttack * std::numbers::pi / 90., airfoilI);

  geometry.rename(0, 2);
  geometry.rename(2, 1, {1, 1});

  // Set material number for inflow
  extend[0] = 2. * deltaX;
  origin[0] = -deltaX;
  IndicatorCuboid2D<T> inflow(extend, origin);
  geometry.rename(2, 3, 1, inflow);

  {
    IndicatorCuboid2D<T> cornerI(deltaX, Vector {-deltaX / 2, lengthY});
    geometry.rename(2, 6, cornerI);
  }
  {
    IndicatorCuboid2D<T> cornerI(deltaX, Vector {-deltaX / 2, -deltaX / 2});
    geometry.rename(2, 6, cornerI);
  }

  {
    IndicatorCuboid2D<T> fringeI(Vector {20 * deltaX, lengthY}, Vector {lengthX - 20 * deltaX, 0});
    geometry.rename(1, 7, fringeI);
  }

  // Set material number for outflow
  origin[0] = lengthX - deltaX;
  IndicatorCuboid2D<T> outflow(extend, origin);
  geometry.rename(2, 4, outflow);
  // Set material number for airfoil
  geometry.rename(1, 5, airfoil);

  // Removes all not needed boundary voxels outside the surface
  geometry.clean(false, {1, 7});
  geometry.checkForErrors();

  clout << "Prepare Geometry ... OK" << std::endl;

  geometry.getStatistics().print();
}

void prepareLattice(MyCase& myCase)
{
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& params     = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes {});
  auto& geometry   = myCase.getGeometry();

  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  // Wallfunction parameters
  //  Used method for density reconstruction
  //  0: use local density
  //  1: extrapolation (Guo)
  //  2: constant (rho = 1.)
  const int rhoMethod = 0;
  //  Used wall profile
  //  0: power law profile
  //  1: Spalding profile
  const int wallFunctionProfile = 1;
  // check if descriptor with body force is used
  const bool bodyForce = false;
  // interpolate sampling velocity along given normal between lattice voxels
  const bool interpolateSampleVelocity = true;
  // use van Driest damping function for turbulent viscosity in boundary cell
  const bool useVanDriest = false;
  //  distance from cell to real wall in lattice units if no geometry indicator is given as input
  const T latticeWallDistance = 0.5;
  //  distance from cell to velocity sampling point in lattice units
  const T    samplingCellDistance = 3.5;
  const bool movingWall           = false;
  const bool averageVelocity      = false;

  WallModelParameters<T> wallModelParameters;
  wallModelParameters.bodyForce                 = bodyForce;
  wallModelParameters.rhoMethod                 = rhoMethod;
  wallModelParameters.samplingCellDistance      = samplingCellDistance;
  wallModelParameters.interpolateSampleVelocity = interpolateSampleVelocity;
  wallModelParameters.useVanDriest              = useVanDriest;
  wallModelParameters.wallFunctionProfile       = wallFunctionProfile;
  wallModelParameters.latticeWallDistance       = latticeWallDistance;
  wallModelParameters.movingWall                = movingWall;
  wallModelParameters.averageVelocity           = averageVelocity;

  // Get parameters for Unit converter
  int resolution    = params.get<parameters::RESOLUTION>();
  int Re            = params.get<parameters::REYNOLDS>();
  T   physViscosity = params.get<parameters::PHYS_CHAR_VISCOSITY>();
  T   chordLength   = params.get<parameters::CHORD_LENGTH>();

  // Set UnitConverter
  lattice.setUnitConverter<UnitConverterFromResolutionAndLatticeVelocity<T, DESCRIPTOR>>(
      (T)resolution,  // resolution: number of voxels per charPhysL
      (T)0.03,        // Max cell speed?
      (T)chordLength, // charPhysLength: reference length of simulation geometry
      (T)Re * physViscosity /
          chordLength,  // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
      (T)physViscosity, // physViscosity: physical kinematic viscosity in __m^2 / s__
      (T)1.0
  );

  // Get airfoil parameters
  Vector<T,2> airfoilCenter = params.get<parameters::CENTER>();
  Vector<int,3> airfoilParams = {params.get<parameters::AIRFOIL_PARAMETERS>()[0],
                                 params.get<parameters::AIRFOIL_PARAMETERS>()[1],
                                 params.get<parameters::AIRFOIL_PARAMETERS>()[2]};
  T              angleOfAttack = params.get<parameters::ANGLE_OF_ATTACK>();

  // Create airfoil indicator
  IndicatorAirfoil2D<T> airfoilI(airfoilCenter, chordLength, airfoilParams[0]/100., airfoilParams[1]/10., airfoilParams[2]/100.);
  IndicatorRotate<T, 2> airfoil(airfoilCenter, angleOfAttack * std::numbers::pi / 90., airfoilI);

  // Bulk material --> Wall Model Dynamics
  setTurbulentWallModelDynamics(lattice, geometry, 1, wallModelParameters);

  // Material=2 -->Full Slip
  boundary::set<boundary::FullSlip>(lattice, geometry, 2); // <-

  // Material=6 -->Bounce Back
  boundary::set<boundary::BounceBack>(lattice, geometry, 6);

  boundary::set<boundary::LocalVelocity>(lattice, geometry, 3);

  using FringeDynamics = dynamics::Tuple<
    MyCase::value_t,
    MyCase::descriptor_t_of<NavierStokes>,
    momenta::BulkTuple,
    equilibria::ThirdOrder,
    collision::ParameterFromCell<collision::LES::SMAGORINSKY,
                                 collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB>>
  >;

  boundary::set<T, DESCRIPTOR, boundary::InterpolatedPressure<T,DESCRIPTOR,FringeDynamics>>(
      lattice, geometry.getMaterialIndicator(4), geometry.getMaterialIndicator(7), geometry.getMaterialIndicator(0));

  // Turbulent Wall Model on airfoil indicator
  setBouzidiBoundary(lattice, geometry, 5, airfoil);
  setTurbulentWallModel(lattice, geometry, 5, wallModelParameters);

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase)
{
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& params     = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes {});
  auto& geometry   = myCase.getGeometry();
  auto& converter  = lattice.getUnitConverter();

  T lengthX = params.get<parameters::DOMAIN_EXTENT>()[0];
  T lengthY = params.get<parameters::DOMAIN_EXTENT>()[1];

  // Initialize all values of distribution functions to their local equilibrium
  auto bulkMaterialIndicator = geometry.getMaterialIndicator({1, 2, 3, 4, 5, 6, 7});

  fields::set<descriptors::VELOCITY>(lattice, bulkMaterialIndicator, std::vector<T>(2, T(0)));
  fields::set<descriptors::POROSITY>(lattice, geometry.getMaterialIndicator({1, 2, 3, 4, 6, 7}), 1);
  fields::set<descriptors::POROSITY>(lattice, geometry.getMaterialIndicator({5}), 0);

  using FringeDynamics = dynamics::Tuple<
    MyCase::value_t,
    MyCase::descriptor_t_of<NavierStokes>,
    momenta::BulkTuple,
    equilibria::ThirdOrder,
    collision::ParameterFromCell<collision::LES::SMAGORINSKY,
                                 collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB>>
  >;

  // Define fringe zone for end of simulation domain
  lattice.defineDynamics<FringeDynamics>(geometry, 7);
  {
    Vector<T, 2>                       origin(lengthX - 20 * converter.getPhysDeltaX(), -converter.getPhysDeltaX());
    Vector<T, 2>                       extend(20 * converter.getPhysDeltaX(), lengthY + 2 * converter.getPhysDeltaX());
    IndicatorCuboid2D<T>               out(extend, origin);
    SuperIndicatorFfromIndicatorF2D<T> outI(out, geometry);

    AnalyticalFfromCallableF<2, T, T> smagorinskyF([&](Vector<T, 2> physR) -> Vector<T, 1> {
      return 0.3 + ((physR[0] - origin[0]) / (20 * converter.getPhysDeltaX())) * (2 - 0.3);
    });

    fields::set<collision::LES::SMAGORINSKY>(lattice, outI, smagorinskyF);
  }

  const T omega = converter.getLatticeRelaxationFrequency();
  lattice.setParameter<descriptors::OMEGA>(omega);
  lattice.setParameter<collision::LES::SMAGORINSKY>(0.3);

  // Make the lattice ready for simulation
  lattice.initialize();
}

void getResults(MyCase& myCase, util::Timer<MyCase::value_t>& timer, std::size_t iT)
{

  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& params     = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes {});
  auto& geometry   = myCase.getGeometry();
  auto& converter  = lattice.getUnitConverter();

  OstreamManager clout(std::cout, "getResults");

  const int vtkIter  = lattice.getUnitConverter().getLatticeTime(.3);
  const int statIter = lattice.getUnitConverter().getLatticeTime(.1);

  // Create AirfoilIndicator
  Vector<T,2> airfoilCenter = params.get<parameters::CENTER>();
  Vector<int,3> airfoilParams = {params.get<parameters::AIRFOIL_PARAMETERS>()[0],
                                 params.get<parameters::AIRFOIL_PARAMETERS>()[1],
                                 params.get<parameters::AIRFOIL_PARAMETERS>()[2]};
  T   angleOfAttack = params.get<parameters::ANGLE_OF_ATTACK>();
  T   chordLength   = params.get<parameters::CHORD_LENGTH>();

  // Create airfoil indicator
  IndicatorAirfoil2D<T> airfoilI(airfoilCenter, chordLength, airfoilParams[0]/100., airfoilParams[1]/10., airfoilParams[2]/100.);
  IndicatorRotate<T, 2> airfoil(airfoilCenter, angleOfAttack * std::numbers::pi / 90., airfoilI);

  if (iT == 0) {
    SuperVTMwriter2D<T> vtmWriter("airfoil2d");
    vtmWriter.createMasterFile();
  }

  if (iT % statIter == 0) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    T point[2] = {};
    point[0]   = airfoilCenter[0] + chordLength;
    point[1]   = airfoilCenter[1];
    SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure(lattice, lattice.getUnitConverter());
    AnalyticalFfromSuperF2D<T>                intpolateP(pressure, true);
    T                                         p;
    intpolateP(&p, point);

    // Timer console output
    timer.update(iT);
    timer.printStep();

    // Lattice statistics console output
    lattice.getStatistics().print(iT, lattice.getUnitConverter().getPhysTime(iT));

    // Drag, lift, pressure drop
    AnalyticalFfromSuperF2D<T>            intpolatePressure(pressure, true);
    SuperLatticePhysDrag2D<T, DESCRIPTOR> drag(lattice, geometry, 5, lattice.getUnitConverter());

    T point1[2] = {};
    T point2[2] = {};

    point1[0] = airfoilCenter[0] - chordLength / 2;
    point1[1] = airfoilCenter[1];

    point2[0] = airfoilCenter[0] + chordLength / 2;
    point2[1] = airfoilCenter[1];

    T p1, p2;
    intpolatePressure(&p1, point1);
    intpolatePressure(&p2, point2);

    clout << "pressure1=" << p1;
    clout << "; pressure2=" << p2;

    T pressureDrop = p1 - p2;
    clout << "; pressureDrop=" << pressureDrop;

    int input[3] = {};
    T   _drag[drag.getTargetDim()];
    drag(_drag, input);
    clout << "; drag=" << _drag[0] << "; lift=" << _drag[1] << std::endl;
  }

  if (iT % vtkIter == 0) {
    SuperVTMwriter2D<T>                                             vtmWriter("airfoil2d");
    SuperLatticePhysVelocity2D<T, DESCRIPTOR>                       velocity(lattice, lattice.getUnitConverter());
    SuperLatticePhysPressure2D<T, DESCRIPTOR>                       pressure(lattice, lattice.getUnitConverter());
    SuperLatticeField2D<T, DESCRIPTOR, descriptors::Y1>             y1(lattice);
    SuperLatticePhysField2D<T, DESCRIPTOR, descriptors::WMVELOCITY> wmvelocity(
        lattice, lattice.getUnitConverter().getConversionFactorVelocity());
    wmvelocity.getName() = "wmvelocity";
    SuperLatticePhysWallShearStress2D<T, DESCRIPTOR> wss(lattice, geometry, 5, lattice.getUnitConverter(), airfoil);

    vtmWriter.addFunctor(velocity);
    vtmWriter.addFunctor(pressure);
    vtmWriter.addFunctor(y1);
    vtmWriter.addFunctor(wmvelocity);
    vtmWriter.addFunctor(wss);

    vtmWriter.write(iT);
  }
}

void simulate(MyCase& myCase)
{
  using T           = MyCase::value_t;
  auto&   params    = myCase.getParameters();
  auto&   lattice   = myCase.getLattice(NavierStokes {});
  auto&   geometry  = myCase.getGeometry();
  auto&   converter = lattice.getUnitConverter();
  const T physMaxT  = params.get<parameters::PHYS_MAX_T>();

  const std::size_t iTmax = myCase.getLattice(NavierStokes {}).getUnitConverter().getLatticeTime(physMaxT);
  util::Timer<T>    timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  SuperLatticeCoupling smoothInflowUpdateO(SmoothInflowUpdateO {}, names::NavierStokes {}, lattice);
  smoothInflowUpdateO.restrictTo(geometry.getMaterialIndicator(3));

  for (std::size_t iT = 0; iT < iTmax; ++iT) {
    {

      const std::size_t iTmaxStart = converter.getLatticeTime(0.1 * physMaxT);
      const std::size_t iTupdate   = converter.getLatticeTime(0.001);

      if (iT % iTupdate == 0 && iT <= iTmaxStart) {
        T                          frac[1] = {};
        PolynomialStartScale<T, T> scale(iTmaxStart, T(1));
        T                          iTvec[1] = {T(iT)};
        scale(frac, iTvec);
        smoothInflowUpdateO.setParameter<SmoothInflowUpdateO::VELOCITY>(
            {frac[0] * converter.getCharLatticeVelocity(), 0});
        smoothInflowUpdateO.apply();
      }
    }

    // === Collide and Stream Execution ===
    lattice.collideAndStream();

    // === Computation and Output of the Results ===
    getResults(myCase, timer, iT);
  }

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
    using T          = MyCase::value_t;
    myCaseParameters.set<RESOLUTION>(100);
    myCaseParameters.set<REYNOLDS>(1400000);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>((T)1.e-5);
    myCaseParameters.set<PHYS_MAX_T>((T)16.);
    myCaseParameters.set<CHORD_LENGTH>((T)1.);
    myCaseParameters.set<AIRFOIL_PARAMETERS>({2, 4, 12});
    myCaseParameters.set<ANGLE_OF_ATTACK>(3.);
    myCaseParameters.set<DOMAIN_EXTENT>([] (MyCase::ParametersD& p) {
      return Vector{6,2} * p.get<CHORD_LENGTH>();
    });
    myCaseParameters.set<CENTER>([] (MyCase::ParametersD& p) {
      auto pos = p.get<DOMAIN_EXTENT>();
      pos[0] /= 4;
      pos[1] /= 2;
      return pos;
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
  return 0;
}
