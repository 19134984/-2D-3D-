/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024-2025 Danial Khazaeipoul
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

/*
 * This example presents dynamics of a single bubble rising in a column of stagnant liquid.
 * The bubble is initialized at the bottom-center of the domain and starts to rise due to the
 * buoyancy force. The free surface Lattice Boltzmann Method (FSLBM) is used to simulate the
 * bubble dynamics. It is important to note that this is a basic version and is only used for
 * demonstration purpose. The actual case must implement a Bubble Model to track and update
 * bubble density for a more realistic simulation.
 */

#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::descriptors;
using namespace olb::graphics;


// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<double, descriptors::D3Q27<  descriptors::FORCE,
  FreeSurface::MASS,
  FreeSurface::EPSILON,
  FreeSurface::CELL_TYPE,
  FreeSurface::CELL_FLAGS,
  FreeSurface::TEMP_MASS_EXCHANGE,
  FreeSurface::PREVIOUS_VELOCITY,
  FreeSurface::HAS_INTERFACE_NBRS>>
>;


/// @brief Create a simulation mesh, based on user-specified geometry
/// @return An instance of a mesh with the relevant information
Mesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;

  const T length = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T height = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const T width  = parameters.get<parameters::DOMAIN_EXTENT>()[2];

  const std::vector<T> origin(3, T(0));
  const std::vector<T> extend = {length, height, width};
  IndicatorCuboid3D<T> cuboid(extend, origin);

  const T physDeltaX = parameters.get<parameters::PHYS_CHAR_LENGTH>() / parameters.get<parameters::RESOLUTION>();

  Mesh<T, MyCase::d> mesh = Mesh<T, MyCase::d>(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(3);

  mesh.getCuboidDecomposition().setPeriodicity({true, true, false});

  return mesh;
}


//- Set the fields based on the initial volume fraction: a spherical gas bubble in a liquid
template<typename T>
class bubble final : public AnalyticalF3D<T, T>
{
protected:
  MyCase& myCase_;

  //- Lattice grid spacing in SI units [m]
  T deltaX_;

  //- Bubble radius in SI units [m]
  T radius_;

  //- Bubble center in SI units [m]
  //std::vector<T> center_;
  std::array<T,3> center_;

  //- Sets field values, i.e., independent of the descriptor dimension; always uses 3 components.
  std::array<T, 3> values_;

public:

  //- Constructor from the input file, Navier-Stokes unit converter, and field values
  bubble
  (
    MyCase& myCase,
    UnitConverter<MyCase::value_t, MyCase::descriptor_t_of<NavierStokes>> const& NavierStokesConverter,
    const std::array<MyCase::value_t,3>& values
  ) : AnalyticalF3D<MyCase::value_t, MyCase::value_t>(1), myCase_(myCase), deltaX_(NavierStokesConverter.getPhysDeltaX()), radius_(MyCase::value_t(0)), values_(values)
  {
    auto& parameters = myCase_.getParameters();
    radius_ =   static_cast<T>(parameters.get<parameters::RADIUS>());
    center_ = { static_cast<T>(parameters.get<parameters::CENTER>()[0]),
                static_cast<T>(parameters.get<parameters::CENTER>()[1]),
                static_cast<T>(parameters.get<parameters::CENTER>()[2]) };
  }

  //- Set the initial volume fraction field based on the bubble geometry, i.e., gas (0), interface (1), fluid (2)
  bool operator()(MyCase::value_t output[], const MyCase::value_t x[]) override {
    //- Initially, set the current node as fluid
    output[0] = values_[2];

    //- Calculate the distance from the bubble center and verify if the node is within the bubble.
    //- Additionally, check if the node lies within the specified limit to be considered as interface.
    Vector<T, 3> fromCtr = {x[0] - center_[0], x[1] - center_[1], x[2] - center_[2]};

    if (norm(fromCtr) <= radius_)
    {
      output[0] = values_[0];
    }
    else
    {
      const T cutoff = deltaX_ * T{1.1};
      bool found = false;

      util::forEachOffset<3>([&](int i, int j, int k) {
        if (!found) {
          Vector<T, 3> offset = {fromCtr[0] + i * cutoff, fromCtr[1] + j * cutoff, fromCtr[2] + k * cutoff};

          if (norm(offset) <= radius_) {
            output[0] = values_[1];
            found = true;
          }
        }
      });
    }

    return true;
  }
};

//- Initializing hydrostatic pressure within the domain
template<typename T>
class HydroStaticPressure3D final : public AnalyticalF3D<T, T>
{
protected:
  MyCase& myCase_;

  //- Navier-Stokes unit converter
  UnitConverter<T, MyCase::descriptor_t_of<NavierStokes>> const& converter_;

  //- Liquid column width in z-direction [m]
  T width_;

  //- The gravity vector [m/s^2]
  std::array<T,3> gravity_;

public:

  //- Constructor from the Advection-Diffusion unit converter
  HydroStaticPressure3D
  (
    MyCase& myCase,
    UnitConverter<T, MyCase::descriptor_t_of<NavierStokes>> const& NavierStokesConverter
  ) : AnalyticalF3D<T, T>(1), myCase_(myCase), converter_(NavierStokesConverter)
  {
    auto& parameters = myCase_.getParameters();
    //- Reading the liquid column width [m] and gravity [m/s^2] from the input file
    width_ = parameters.get<parameters::DOMAIN_EXTENT>()[2];
    //configuration["Application"]["PhysParameters"]["PhysGravity"].read(gravity_);
    gravity_ = { static_cast<T>(parameters.get<parameters::GRAVITY>()[0]),
                 static_cast<T>(parameters.get<parameters::GRAVITY>()[1]),
                 static_cast<T>(parameters.get<parameters::GRAVITY>()[2]) };
  }

  //- Density computed using hydrostatic pressure (p = p0 + rho * g * h),
  //  with the reference density established at the midpoint in the z-direction
  bool operator()(T output[], const T x[]) override {
    output[0] = converter_.getLatticeDensityFromPhysPressure
    (
      converter_.getCharPhysPressure() + converter_.getPhysDensity() * util::abs(gravity_[2]) * (T(0.5) * width_ - x[2])
    );

    return true;
  }
};


//- Calculate the center of mass of the bubble, i.e., works only for a single bubble
Vector<MyCase::value_t, 3> bubbleCenterOfMass(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();

  //- Initialize the center of mass of the bubble
  Vector<T, 3> centerOfMass(T(0));
  T nCells = T(0);

  for (int iC = 0; iC < lattice.getLoadBalancer().size(); ++iC)
  {
    const auto& block = lattice.getBlock(iC);
    block.forCoreSpatialLocations([&](auto iX, auto iY, auto iZ) {
      auto cell = block.get({iX, iY, iZ});

      //- Calculate the center of mass of the bubble, i.e., is Gas Type (0)
      if (isCellType(cell, FreeSurface::Type::Gas))
      {
        centerOfMass += geometry.getBlockGeometry(iC).getPhysR({iX, iY, iZ});
        nCells += T(1);
      }
    });
  }

#ifdef PARALLEL_MODE_MPI
  for (size_t i = 0; i < centerOfMass.size(); ++i)
  {
    singleton::mpi().reduceAndBcast(centerOfMass[i], MPI_SUM);
  }
  singleton::mpi().reduceAndBcast(nCells, MPI_SUM);
#endif

  centerOfMass[0] /= nCells;
  centerOfMass[1] /= nCells;
  centerOfMass[2] /= nCells;

  return centerOfMass;
}

//- Prepare the geometry using the indicator method
void prepareGeometry(MyCase& myCase)
{
  //using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  OstreamManager clout(std::cout, "prepareGeometry");

  //- Setting the material ID to (1) for all voxels, i.e. domain-wide
  //- Note that the side walls in x- and y-directions are periodic
  geometry.rename(0, 1);

  //- Reading cuboid length, height, and width [m] from the input file
  const T length = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T height = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const T width  = parameters.get<parameters::DOMAIN_EXTENT>()[2];
  const T epsilon = T(0.5) * parameters.get<parameters::PHYS_CHAR_LENGTH>() / parameters.get<parameters::RESOLUTION>();
  std::vector<T> origin = {-epsilon, -epsilon, -epsilon};
  std::vector<T> extend = {length + T(2) * epsilon, height + T(2) * epsilon, T(2) * epsilon};

  //- Setting material ID 2 for the bottom wall boundary voxels
  IndicatorCuboid3D<T> bottomWall(extend, origin);
  geometry.rename(1, 2, bottomWall);

  //- Setting material ID 3 for the top wall boundary voxels
  origin[2] += width;
  IndicatorCuboid3D<T> topWall(extend, origin);
  geometry.rename(1, 3, topWall);

  //- Removes all not needed boundary voxels outside the surface
  //- Removes all not needed boundary voxels inside the surface
  //- Checks the geometry for possible errors
  geometry.clean();
  geometry.innerClean();
  geometry.checkForErrors();
  geometry.print();

  clout << "Geometry created successfully." << std::endl;
}

//- Prepare lattice for the Navier-Stokes equations with the free-surface model
void prepareLattice (MyCase& myCase) {
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& parameters = myCase.getParameters();

  OstreamManager clout(std::cout, "prepareLattice");

  const int resolution = parameters.get<parameters::RESOLUTION>();
  const T lattice_relaxation_time = parameters.get<parameters::LATTICE_RELAXATION_TIME>();
  const T charPhysLength = parameters.get<parameters::PHYS_CHAR_LENGTH>();
  const T charPhysVelocity = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T physViscosity = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T physDensity = parameters.get<parameters::PHYS_CHAR_DENSITY>();

  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR>>(
    int {resolution},              // resolution: number of voxels per charPhysL
    (T) lattice_relaxation_time,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T) charPhysLength,            // charPhysLength: reference length of simulation geometry
    (T) charPhysVelocity,          // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T) physViscosity,             // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T) physDensity                // physDensity: physical density in __kg / m^3__
  );

  auto& converter = lattice.getUnitConverter();

  //converter.print();
  //converter.write("NavierStokesConverter");

  //- Lattice relaxation frequency for the Navier-Stokes equations
  T NavierStokesOmega = converter.getLatticeRelaxationFrequency();

  //- rho0 is the initial density profile calculated from the hydro-static pressure
  HydroStaticPressure3D<T> rho0(myCase, converter);

  //- velocity0 is the uniform initial velocity profile
  AnalyticalConst3D<T, T> velocity0(T(0), T(0), T(0));

  //- Setting NSEs bulk dynamics for the defined material IDs
  //- Bulk materials: (1) fluid, (2) bottomWall, and (3) topWall
  auto bulkIndicator = geometry.getMaterialIndicator({1, 2, 3});
  lattice.defineDynamics<NoDynamics>(geometry, 0);
  lattice.defineDynamics<SmagorinskyForcedBGKdynamics<T, DESCRIPTOR>>(bulkIndicator);

  //- Setting NSEs velocity boundary condition for the material (2), i.e. bottomWall
  // boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 2);
  boundary::set<boundary::BounceBack>(lattice, geometry, 2);

  //- Setting NSEs velocity boundary condition for the material (3), i.e. topWall
  // boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 3);
  boundary::set<boundary::BounceBack>(lattice, geometry, 3);

  //- Setting the lattice relaxation frequency and Smagorinsky constant
  lattice.setParameter<descriptors::OMEGA>(NavierStokesOmega);
  lattice.setParameter<collision::LES::SMAGORINSKY>(T(0.2));

  //- Define user-defined constants for field initialization
  AnalyticalConst3D<T, T> zero(T(0));
  AnalyticalConst3D<T, T> one(T(1));
  AnalyticalConst3D<T, T> two(T(2));
  AnalyticalConst3D<T, T> four(T(4));
  AnalyticalConst3D<T, T> zeros(T(0), T(0), T(0));

  //- Create a single bubble with the specified center and radius
  bubble<T> types( myCase, converter, {T(0), T(1), T(2)});
  bubble<T> mass( myCase, converter, {T(0), T(0.5), T(1)});

  //- Initialize the free surface fields on the lattice
  for (auto i : {0, 1, 2, 3}) {
    lattice.defineField<FreeSurface::MASS>(geometry, i, zero);
    lattice.defineField<FreeSurface::EPSILON>(geometry, i, zero);
    lattice.defineField<FreeSurface::CELL_TYPE>(geometry, i, zero);
    lattice.defineField<FreeSurface::CELL_FLAGS>(geometry, i, zero);
    lattice.defineField<FreeSurface::PREVIOUS_VELOCITY>(geometry, i, zeros);
    lattice.defineField<FORCE>(geometry, i, zeros);
  }

  //- Setting the initial types, volume fraction, and mass based on the bubble geometry for the material IDs
  //- Material ID 1 denotes the internal fluid nodes
  lattice.defineField<FreeSurface::CELL_TYPE>(geometry, 1, types);
  lattice.defineField<FreeSurface::MASS>(geometry, 1, mass);
  lattice.defineField<FreeSurface::EPSILON>(geometry, 1, mass);

  //- Materials ID 2 and 3 denote the bottom and top wall boundary voxels
  //- Note that epsilon is set to 1.0 for the wall boundary voxels, i.e., a wetting wall condition
  for (auto i : {2, 3}) {
    lattice.defineField<FreeSurface::CELL_TYPE>(geometry, i, four);
    lattice.defineField<FreeSurface::EPSILON>(geometry, i, one);
  }

  //- Setting the gravity force in the z-direction, i.e., in lattice units
  std::array<T,3> gravity;
  gravity[0] = parameters.get<parameters::GRAVITY>()[0];
  gravity[1] = parameters.get<parameters::GRAVITY>()[1];
  gravity[2] = parameters.get<parameters::GRAVITY>()[2];
  T gScale = T(1.0) / converter.getConversionFactorForce() * converter.getConversionFactorMass();
  AnalyticalConst3D<T,T> gForce{gravity[0] * gScale, gravity[1] * gScale, gravity[2] * gScale};
  lattice.defineField<FORCE>(geometry, 1, gForce);

  //- Setting the initial conditions for the Navier-Stokes equations
  lattice.defineRhoU(bulkIndicator, rho0, velocity0);
  lattice.iniEquilibrium(bulkIndicator, rho0, velocity0);

  FreeSurface::initialize(lattice);
  lattice.initialize();

  static FreeSurface3DSetup<T, DESCRIPTOR> freeSurfaceHandler(lattice);
  freeSurfaceHandler.addPostProcessor();

  T tScale = util::pow(converter.getConversionFactorTime(), 2)
    / (converter.getPhysDensity() * util::pow(converter.getPhysDeltaX(), 3));

  bool drop_isolated_cells = parameters.get<FreeSurface::DROP_ISOLATED_CELLS>();
  bool hasSurfaceTension = parameters.get<FreeSurface::HAS_SURFACE_TENSION>();
  T surfaceTensionCoeff = parameters.get<FreeSurface::SURFACE_TENSION_PARAMETER>();
  T transitionThreshold = parameters.get<FreeSurface::TRANSITION>();
  T lonelyThreshold = parameters.get<FreeSurface::LONELY_THRESHOLD>();

  //- Setting the free-surface parameters
  lattice.setParameter<FreeSurface::DROP_ISOLATED_CELLS>(drop_isolated_cells);
  lattice.setParameter<FreeSurface::TRANSITION>(transitionThreshold);
  lattice.setParameter<FreeSurface::LONELY_THRESHOLD>(lonelyThreshold);
  lattice.setParameter<FreeSurface::HAS_SURFACE_TENSION>(hasSurfaceTension);
  lattice.setParameter<FreeSurface::SURFACE_TENSION_PARAMETER>(tScale * surfaceTensionCoeff);
  lattice.setParameter<FreeSurface::FORCE_DENSITY>({gravity[0] * gScale, gravity[1] * gScale, gravity[2] * gScale});

  clout << "Lattice created successfully." << std::endl;
}


//- Set the output data for post-processing and visualization
void postprocess(MyCase& myCase, const std::size_t& iT, util::Timer<MyCase::value_t>& timer, std::size_t& iTPrev, Vector<MyCase::value_t, 3>& CenterOfMassPrev)
{
  OstreamManager clout(std::cout, "postprocess");
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& converter = lattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  SuperVTMwriter3D<T> vtmWriter("basicRisingBubble3d");
  static CSV<T> csvWriter("basicRisingBubble3d", ',', {"Time", "X", "Y", "Z", "RiseVelocity"});
  T maxPhysTime = parameters.get<parameters::PHYS_MAX_T>();
  const bool lastTimeStep = (iT + 1 == converter.getLatticeTime(maxPhysTime));

  const T vtkSaveTime = parameters.get<parameters::PHYS_VTK_ITER_T>();
  const T logTime     = parameters.get<parameters::PHYS_STAT_ITER_T>();

  if (iT == 0) {
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(lattice);
    SuperLatticeRank3D<T, DESCRIPTOR> rank(lattice);

    vtmWriter.write(cuboid);
    vtmWriter.write(rank);
    vtmWriter.createMasterFile();

    T surfaceTensionCoeff = parameters.get<FreeSurface::SURFACE_TENSION_PARAMETER>();

    //- Calculate the Bond number for the rising bubble
    const T Bond = util::abs(T(4) * converter.getPhysDensity()
                                  * parameters.get<parameters::GRAVITY>()[2]
                                  * util::pow(parameters.get<parameters::RADIUS>(), 2)
                                  / surfaceTensionCoeff);
    clout << "Bond number = " << Bond << std::endl;

    //- Calculate the Morton number for the rising bubble
    const T Morton = util::abs(parameters.get<parameters::GRAVITY>()[2] * util::pow(converter.getPhysViscosity(), 4)
                               * util::pow(parameters.get<parameters::PHYS_CHAR_DENSITY>() / surfaceTensionCoeff, 3));
    clout << "Morton number = " << Morton << std::endl;
  }

  //- Write the output data for visualization in VTK format
  if (iT % converter.getLatticeTime(vtkSaveTime) == 0 || lastTimeStep) {
    //- Triggers data transfer between the host (CPU) and the device (GPU)
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(lattice, converter);
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(lattice, converter);
    SuperLatticeDensity3D<T, DESCRIPTOR> density(lattice);
    SuperLatticeField3D<T, DESCRIPTOR, FreeSurface::CELL_TYPE> types(lattice);
    SuperLatticeField3D<T, DESCRIPTOR, FreeSurface::EPSILON> epsilon(lattice);
    SuperLatticeField3D<T, DESCRIPTOR, FreeSurface::MASS> mass(lattice);

    vtmWriter.addFunctor(velocity, "Velocity");
    vtmWriter.addFunctor(pressure, "Pressure");
    vtmWriter.addFunctor(density, "Density");
    vtmWriter.addFunctor(types, "Types");
    vtmWriter.addFunctor(epsilon, "Epsilon");
    vtmWriter.addFunctor(mass, "Mass");
    vtmWriter.write(iT);

    //- Write bubble statistics in SI units: id, center of mass [m], and rise velocity [m/s]
    Vector<T, 3> centerOfMass = bubbleCenterOfMass(myCase);
    const T time = converter.getPhysTime(iT);
    const T timePrev = converter.getPhysTime(iTPrev);
    const T riseVelocity = (centerOfMass[2] - CenterOfMassPrev[2]) / (time - timePrev + util::numericLimits::epsilon<T>());
    csvWriter.writeDataFile(time, {centerOfMass[0], centerOfMass[1], centerOfMass[2], riseVelocity});
    iTPrev = iT;
    CenterOfMassPrev = centerOfMass;
    const T Reynolds = riseVelocity * converter.getCharPhysLength() / converter.getPhysViscosity();
    clout << "Calculated Re = " << Reynolds << ", Experimental Re = 94" << std::endl;
    clout << "Bubble's center of mass = (" << centerOfMass[0] << ", " << centerOfMass[1] << ", " << centerOfMass[2] << ")" << std::endl;
  }

  //- Write the output data for logging in the console
  if (iT % converter.getLatticeTime(logTime) == 0 || lastTimeStep) {
    timer.update(iT);
    timer.printStep();
    lattice.getStatistics().print(iT, converter.getPhysTime(iT));
  }
}


/// Set initial condition for primal variables (velocity and density)
/// @param myCase The Case instance which keeps the simulation data
/// @note Be careful: initial values have to be set using lattice units
void setInitialValues( MyCase& myCase ){
  // Nothing to do here
}

/// Update boundary values at times (and external fields, if they exist)
/// @param myCase The Case instance which keeps the simulation data
/// @param iT The time step
/// @note Be careful: boundary values have to be set using lattice units
void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{
  // Nothing to do here, because simulation does not depend on time
}


/// @brief Execute simulation: set initial values and run time loop
/// @param myCase The Case instance which keeps the simulation data
void simulate( MyCase& myCase ){

  using T = MyCase::value_t;
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& geometry = myCase.getGeometry();
  auto& converter = lattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  // Prepare the main time loop using Timer and convergence check
  T maxPhysTime = parameters.get<parameters::PHYS_MAX_T>();
  util::Timer<T> timer(converter.getLatticeTime(maxPhysTime), geometry.getStatistics().getNvoxel());
  timer.start();

  //- Enable statistics by tracking bubble's center of mass
  std::size_t iTPrev = 0;
  Vector<T, 3> CenterOfMassPrev = bubbleCenterOfMass(myCase);

  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysTime); ++iT) {
    //-- 5th Step: write the output data for post-processing and visualization
    postprocess( myCase, iT, timer, iTPrev, CenterOfMassPrev);

    //-- 6th Step: evolve the Navier-Stokes equations using FSLBM approximation
    lattice.collideAndStream();
  }

  timer.stop();
  timer.printSummary();

}


/// Setup and run a simulation
int main(int argc, char* argv[])
{
  //-- 1st Step: Reading input files and solver initialization
  initialize(&argc, &argv, false, false);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION>(16);
    myCaseParameters.set<CENTER>({0.1044, 0.1044, 0.0261});
    myCaseParameters.set<olb::parameters::RADIUS>(0.01305);
    myCaseParameters.set<DOMAIN_EXTENT>({0.2088, 0.2088, 0.522});

    myCaseParameters.set<PHYS_CHAR_LENGTH>(0.0261);
    myCaseParameters.set<PHYS_CHAR_PRESSURE>(0.0);
    myCaseParameters.set<PHYS_CHAR_VELOCITY>(0.3533088193);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(0.0000980996);
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1095.0);
    myCaseParameters.set<PHYS_MAX_T>(1.0);
    myCaseParameters.set<GRAVITY>({0., 0., -9.81});
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(.501);

    myCaseParameters.set<PHYS_VTK_ITER_T>(0.05);
    myCaseParameters.set<PHYS_STAT_ITER_T>(0.05);

    myCaseParameters.set<FreeSurface::DROP_ISOLATED_CELLS>(true);
    myCaseParameters.set<FreeSurface::HAS_SURFACE_TENSION>(true);
    myCaseParameters.set<FreeSurface::SURFACE_TENSION_PARAMETER>(0.0636306414);
    myCaseParameters.set<FreeSurface::TRANSITION>(0.01);
    myCaseParameters.set<FreeSurface::LONELY_THRESHOLD>(0.1);
  }

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry( myCase );

  /// === Step 6: Prepare Lattice ===
  prepareLattice( myCase );

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);

}
