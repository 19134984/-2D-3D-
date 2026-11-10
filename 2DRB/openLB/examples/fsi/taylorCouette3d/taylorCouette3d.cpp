/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Fedor Bukreev, Adrian Kummerlaende
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
 *  Boston, MA  02110-1201, USA.
 */

// Taylor-Couette flow with rotating inner and static outer cylinders
// Literature with parameters and results for this example can be found here:
// https://doi.org/10.1146/annurev-fluid-122414-034353
// https://doi.org/10.1017/jfm.2014.134
// https://doi.org/10.1063/1.4863312

#include <olb.h>

using namespace olb;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<float, descriptors::D3Q19<
    descriptors::POROSITY,
    descriptors::VELOCITY,
    fields::fsi::ELEMENT_TAG,
    descriptors::AVERAGE_VELOCITY
  >>
>;

namespace olb::parameters {

struct ETA : public descriptors::FIELD_BASE<1> { };
struct INNER_R : public descriptors::FIELD_BASE<1> { };
struct OUTER_R : public descriptors::FIELD_BASE<1> { };
struct SIGMA : public  descriptors::FIELD_BASE<1> { };
struct ROTATE_OMEGA : public  descriptors::FIELD_BASE<1> { };
struct PHYS_AVERAGING_START_T : public  descriptors::FIELD_BASE<1> { };

}

template <typename T, typename S>
class TaylorCouette3D : public AnalyticalF3D<T,S> {
private:
  T turbulenceIntensity;
  T maxVelocity;
  T a;
  T b;
  T innerR;
  T outerR;

  std::default_random_engine generator;

public:
  TaylorCouette3D(MyCase& myCase)
    : AnalyticalF3D<T,S>(3)
    , generator(0x1337533DAAAAAAAA)
  {
    auto& params = myCase.getParameters();
    const auto& converter = myCase.getLattice(NavierStokes()).getUnitConverter();
    innerR = params.get<parameters::INNER_R>();
    outerR = params.get<parameters::OUTER_R>();
    turbulenceIntensity = params.get<parameters::TURBULENCE_INTENSITY>();
    maxVelocity = converter.getLatticeVelocity(params.get<parameters::PHYS_CHAR_VELOCITY>());
    a = -1.;
    b =  1.;
  };

  bool operator()(T output[], const S input[])
  {
    std::uniform_real_distribution<BaseType<T>> distribution(a, b);
    T nRandom1 = distribution(generator);
    T nRandom2 = distribution(generator);
    T nRandom3 = distribution(generator);

    output[0] = 0.;
    output[1] = 0.;
    output[2] = 0.;

    T radius = util::sqrt(util::pow(input[0] - outerR,2) + util::pow(input[1] - outerR,2));

    if( radius >= innerR && radius <= outerR ){

    T u_calc_tang = maxVelocity*(outerR - radius)/(outerR-innerR);
    T u_calc_x = u_calc_tang*util::sin(std::acos((input[0] - outerR)/radius));
    T u_calc_y = -u_calc_tang*util::cos(std::acos((input[0] - outerR)/radius));

    if(input[1] - outerR < T(0)) u_calc_x *= T(-1);

    output[0] = turbulenceIntensity*nRandom1*u_calc_x + u_calc_x;
    output[1] = turbulenceIntensity*nRandom2*u_calc_y + u_calc_y;
    output[2] = turbulenceIntensity*nRandom3*u_calc_tang;
    }else{
    output[0] = 0.;
    output[1] = 0.;
    output[2] = 0.;
    }

    return true;
  };
};

std::shared_ptr<IndicatorF3D<MyCase::value_t>> createDomainI(MyCase::ParametersD& params) {
  using T = MyCase::value_t;

  const T outerR = params.get<parameters::OUTER_R>();
  const T length = params.get<parameters::DOMAIN_L>();

  std::shared_ptr<IndicatorF3D<T>> boxI(new IndicatorCuboid3D<T>({2*outerR,2*outerR,length},{0,0,0}));
  return boxI;
}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;

  const T gapW = params.get<parameters::OUTER_R>() - params.get<parameters::INNER_R>();
  const T deltaX = gapW / params.get<parameters::RESOLUTION>();

  auto domainI = createDomainI(params);

  Mesh<T,MyCase::d> mesh(*domainI, deltaX, singleton::mpi().getSize());
  mesh.getCuboidDecomposition().setPeriodicity({false, false, true});
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  auto& params = myCase.getParameters();
  auto domainI = createDomainI(params);
  auto& sGeometry = myCase.getGeometry();

  sGeometry.rename(0, 3);
  sGeometry.rename(3, 2, {1,1,0});

  sGeometry.checkForErrors();
  sGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase) {
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

  auto& sGeometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& sLattice = myCase.getLattice(NavierStokes());

  const T gapW = params.get<parameters::OUTER_R>() - params.get<parameters::INNER_R>();

  sLattice.setUnitConverter<UnitConverterFromResolutionAndLatticeVelocity<T,DESCRIPTOR>>(
    (int) params.get<parameters::RESOLUTION>(), // resolution of gap
    (T)   0.1,              // CFL
    (T)   gapW,             // charPhysLength: reference length of simulation geometry
    (T)   params.get<parameters::PHYS_CHAR_VELOCITY>(),   // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   params.get<parameters::PHYS_CHAR_VISCOSITY>(), // physViscosity: physical kinematic viscosity in __m^3 / s__
    (T)   1.                // physDensity: physical density in __kg / m^3__
  );
  const auto& converter = sLattice.getUnitConverter();
  converter.print();

  WallModelParameters<T> wallModelParameters;
  wallModelParameters.bodyForce = false;
  wallModelParameters.rhoMethod = 0;
  wallModelParameters.fNeqMethod = 0;
  wallModelParameters.samplingCellDistance = 2.5;
  wallModelParameters.interpolateSampleVelocity = true;
  wallModelParameters.useVanDriest = false;
  wallModelParameters.wallFunctionProfile = 1;
  wallModelParameters.movingWall = true;
  wallModelParameters.averageVelocity = true;
  setTurbulentWallModelDynamics(sLattice,
                                sGeometry.getMaterialIndicator({1,2}),
                                wallModelParameters);

  sLattice.addPostProcessor(sGeometry.getMaterialIndicator(2),
                            meta::id<TurbulentWallModelPostProcessor<false,true,false,1,true>>{});
  sLattice.addPostProcessor(sGeometry.getMaterialIndicator({1,2}),
                            meta::id<TurbulentWallModelPorousFneqFDMPostProcessor<0>>{});

  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());
  sLattice.setParameter<collision::LES::SMAGORINSKY>(0.2);
  sLattice.setParameter<descriptors::SAMPLING_DISTANCE>(wallModelParameters.samplingCellDistance);

  boundary::set<boundary::BounceBack>(sLattice, sGeometry, 3);
  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;

  auto& params  = myCase.getParameters();
  auto& sGeometry = myCase.getGeometry();
  auto& sLattice  = myCase.getLattice(NavierStokes());
  const auto& converter  = sLattice.getUnitConverter();

  const T outerR = params.get<parameters::OUTER_R>();
  const T length = params.get<parameters::DOMAIN_L>();

  AnalyticalConst3D<T,T> rhoF(1);
  TaylorCouette3D<T,T> uCalc(myCase);
  auto bulkIndicator = sGeometry.getMaterialIndicator({1,2});
  sLattice.iniEquilibrium(bulkIndicator, rhoF, uCalc);
  sLattice.defineRhoU(bulkIndicator, rhoF, uCalc);
  fields::set<descriptors::AVERAGE_VELOCITY>(sLattice, sGeometry.getMaterialIndicator({0,1,2,3}), 0);

  // set everywhere 0 porosity
  fields::set<descriptors::POROSITY>(sLattice, sGeometry.getMaterialIndicator({0,1,2,3}), 0);

  {
    // set porosity & wall-model-porosity 1 in the bulk and rotating inner cylinder
    IndicatorCylinder3D<T> cylinderI({outerR - T{0.25}*converter.getPhysDeltaX(), outerR - T{0.25}*converter.getPhysDeltaX(), -T{5.}*converter.getPhysDeltaX()},{outerR - T{0.25}*converter.getPhysDeltaX(), outerR - T{0.25}*converter.getPhysDeltaX(), length+T{10.}*converter.getPhysDeltaX()}, outerR);
    SuperIndicatorFfromIndicatorF3D<T> cylinderIndicatorF(cylinderI, sGeometry);
    fields::set<descriptors::POROSITY>(sLattice, cylinderIndicatorF, 1);
    fields::set<descriptors::WMPOROSITY>(sLattice, bulkIndicator, 1);
  }
  {
    // set Y1 vectors on the outer cylinder wall
    IndicatorCuboid3D<T> cuboid({2*outerR,2*outerR,length},{0,0,0});
    std::shared_ptr<IndicatorF3D<T>> cuboidLayer(new IndicatorLayer3D<T>(cuboid, T{5.}*converter.getPhysDeltaX()));
    std::shared_ptr<IndicatorF3D<T>> cylinderOutI(new IndicatorCylinder3D<T>({outerR - T{0.25}*converter.getPhysDeltaX(), outerR - T{0.25}*converter.getPhysDeltaX(), -T{5.}*converter.getPhysDeltaX()},{outerR - T{0.25}*converter.getPhysDeltaX(), outerR - T{0.25}*converter.getPhysDeltaX(), length+T{10.}*converter.getPhysDeltaX()}, outerR));
    SuperIndicatorFfromIndicatorF3D<T> cylinderIndicatorOutF(cylinderOutI, sGeometry);
    SuperIndicatorFfromIndicatorF3D<T> outerF((cuboidLayer - cylinderOutI), sGeometry);

    setWallDistance(sLattice, outerF, cylinderIndicatorOutF, cylinderOutI.get());
  }

  sLattice.initialize();
}

void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  OstreamManager clout(std::cout, "getResults");

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

  auto& sGeometry = myCase.getGeometry();
  auto& sLattice  = myCase.getLattice(NavierStokes());
  const auto& converter  = sLattice.getUnitConverter();
  auto& params  = myCase.getParameters();

  const int vtkIter  = converter.getLatticeTime(1);
  const int statIter = converter.getLatticeTime(1./30);

  if (iT == 0) {
    {
      SuperVTMwriter3D<T> vtmWriter("flow");
      vtmWriter.createMasterFile();
    }
  }

  if (iT % statIter == 0) {
    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    if (std::isnan(sLattice.getStatistics().getAverageRho())) {
      std::exit(-1);
    }
  }

  const T gapW = params.get<parameters::OUTER_R>() - params.get<parameters::INNER_R>();
  const T charPhysU = params.get<parameters::PHYS_CHAR_VELOCITY>();
  std::size_t iTstartAvg = converter.getLatticeTime(200.*gapW/charPhysU);

  if (iT == iTstartAvg) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    SuperLatticeVelocity3D<T,DESCRIPTOR> latticeVelocityF(sLattice);
    fields::set<descriptors::AVERAGE_VELOCITY>(sLattice,
                                               sGeometry.getMaterialIndicator({1,2,3,4,5,6,7}),
                                               latticeVelocityF);
  }
  if (iT < iTstartAvg) {
    sLattice.setParameter<descriptors::LATTICE_TIME>(2);
  }
  else {
    sLattice.setParameter<descriptors::LATTICE_TIME>(iT - iTstartAvg + 1);
  }

  if (iT % vtkIter == 0 && iT > converter.getLatticeTime(1.)) {
    sLattice.executePostProcessors(stage::Evaluation{});
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    sLattice.scheduleBackgroundOutputVTK([&,iT,iTstartAvg](auto task) {
      SuperVTMwriter3D<T> vtkWriter("flow");
      SuperLatticePhysVelocity3D velocity(sLattice, converter);
      SuperLatticePhysPressure3D pressure(sLattice, converter);
      vtkWriter.addFunctor(pressure);
      vtkWriter.addFunctor(velocity);

      SuperLatticePhysField3D<T,DESCRIPTOR,descriptors::AVERAGE_VELOCITY> velocityAv(
        sLattice, converter.getConversionFactorVelocity());
      velocityAv.getName() = "avgU";

      if (iT > iTstartAvg) {
        vtkWriter.addFunctor(velocityAv);
      }

      task(vtkWriter, iT);
    });
  }
}

void simulate(MyCase& myCase) {
  using T = MyCase::value_t;

  auto& params = myCase.getParameters();
  auto& sGeometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes());
  const auto& converter = sLattice.getUnitConverter();

  SuperPorousElementEmbeddingO fsiEmbeddingO(sLattice);
  fsiEmbeddingO.registerElementType<CylinderWithWallModelPorosityF>(sGeometry.getMaterialIndicator(2));

  const T innerR = params.get<parameters::INNER_R>();
  const T outerR = params.get<parameters::OUTER_R>();
  const T length = params.get<parameters::DOMAIN_L>();
  const T omega = params.get<parameters::ROTATE_OMEGA>();
  Vector center{outerR - T{0.25}*converter.getPhysDeltaX(),
                outerR - T{0.25}*converter.getPhysDeltaX(),
                       - T{   5}*converter.getPhysDeltaX()};
  Vector<T,3> axis{0,0,1};
  auto valveParams = makeParametersD<T,descriptors::SPATIAL_DESCRIPTOR<3>>(
    fields::fsi::ELEMENT_TAG{},      1,
    fields::fsi::ELEMENT_PIVOT{},    center,
    fields::fsi::ELEMENT_ROTATION{}, 0,
    fields::fsi::ELEMENT_U_ROTATION{}, Vector<T,3>{0, 0, -omega},
    CylinderPorosityF::RADIUS{},     innerR,
    CylinderPorosityF::CENTER{},     center,
    CylinderPorosityF::LENGTH{},     (length+T{10}*converter.getPhysDeltaX()),
    CylinderPorosityF::AXIS{},       axis
  );
  fsiEmbeddingO.add(valveParams);

  sLattice.setParameter<fields::converter::PHYS_LENGTH>(
    1/converter.getConversionFactorLength());
  sLattice.setParameter<fields::converter::PHYS_VELOCITY>(
    1/converter.getConversionFactorVelocity());
  sLattice.setParameter<fields::converter::PHYS_DELTA_X>(
    converter.getPhysDeltaX());

  sLattice.setProcessingContext(ProcessingContext::Simulation);
  fsiEmbeddingO.initialize();

  const std::size_t iTmax = converter.getLatticeTime(
    params.get<parameters::PHYS_MAX_T>());

  util::Timer<T> timer(iTmax, sGeometry.getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    /// === Step 8.1: Collide and Stream Execution ===
    myCase.getLattice(NavierStokes{}).collideAndStream();

    sLattice.stripeOffDensityOffset(sLattice.getStatistics().getAverageRho()-(T)1);

    /// === Step 8.2: Computation and Output of the Results ===
    getResults(myCase, timer, iT);
  }

  sLattice.setProcessingContext(ProcessingContext::Evaluation);

  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[]) {
  initialize( &argc, &argv );

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION>(20);
    myCaseParameters.set<ETA>(0.5);
    myCaseParameters.set<INNER_R>(0.3);
    myCaseParameters.set<TAYLOR>(1e10);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(1e-5);
    myCaseParameters.set<TURBULENCE_INTENSITY>(0.05);

    myCaseParameters.set<OUTER_R>([] (MyCase::ParametersD& p) {
      return p.get<INNER_R>() / p.get<ETA>();
    });
    myCaseParameters.set<DOMAIN_L>([] (MyCase::ParametersD& p) {
      const auto gapW = p.get<OUTER_R>() - p.get<INNER_R>();
      return 2 * std::numbers::pi * gapW / 3;
    });
    myCaseParameters.set<ROTATE_OMEGA>([] (MyCase::ParametersD& p) {
      const auto innerR = p.get<INNER_R>();
      const auto outerR = p.get<OUTER_R>();
      const auto sigma = util::pow((innerR + outerR)/2./util::sqrt(innerR*outerR), 4);
      return util::sqrt(p.get<TAYLOR>()*util::pow(p.get<PHYS_CHAR_VISCOSITY>(), 2)/(0.25*sigma*util::pow(outerR-innerR, 3)));
    });
    myCaseParameters.set<PHYS_CHAR_VELOCITY>([] (MyCase::ParametersD& p) {
      return p.get<ROTATE_OMEGA>() * p.get<INNER_R>();
    });
    myCaseParameters.set<PHYS_MAX_T>([] (MyCase::ParametersD& p) {
      const auto gapW = p.get<OUTER_R>() - p.get<INNER_R>();
      return 1000 * gapW / p.get<PHYS_CHAR_VELOCITY>();
    });
    myCaseParameters.set<PHYS_AVERAGING_START_T>([] (MyCase::ParametersD& p) {
      const auto gapW = p.get<OUTER_R>() - p.get<INNER_R>();
      return 200 * gapW / p.get<PHYS_CHAR_VELOCITY>();
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
