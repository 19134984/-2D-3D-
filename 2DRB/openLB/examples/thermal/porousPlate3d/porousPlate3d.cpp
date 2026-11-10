/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Florian Kaiser
 *                2008 Orestis Malaspinas
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

// #define AdvectionDiffusionDirichletBC
// #define RegularizedTemperatureBoundary
#define RegularizedHeatFluxBC

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<float, descriptors::D3Q19<descriptors::FORCE>>,
  Temperature,  Lattice<float, descriptors::D3Q7<descriptors::VELOCITY>>
>;

template <typename T, typename S>
class AnalyticalVelocity : public AnalyticalF3D<T, S>
{
  private:
    T _Re;
    T _u0;
    T _v0;
    T _ly;
  public:
    AnalyticalVelocity(T Re, T u0, T v0, T ly) : AnalyticalF3D<T, S>(3),
      _Re(Re), _u0(u0), _v0(v0), _ly(ly)
    {
      this->getName() = "AnalyticalVelocity";
    };

    bool operator()(T output[3], const S x[3]) override
    {
      output[0] = _u0*((util::exp(_Re* x[1] / _ly) - 1) / (util::exp(_Re) - 1));
      output[1] = _v0;
      output[2] = 0.0;
      return true;
    };
};

template <typename T, typename S>
class AnalyticalTemperature : public AnalyticalF3D<T, S>
{
  private:
    T _Re;
    T _Pr;
    T _ly;
    T _T0;
    T _deltaT;
  public:
  AnalyticalTemperature(T Re, T Pr, T ly, T T0, T deltaT) : AnalyticalF3D<T, S>(1),
    _Re(Re), _Pr(Pr), _ly(ly), _T0(T0), _deltaT(deltaT)
  {
    this->getName() = "AnalyticalTemperature";
  };

  bool operator()(T output[1], const S x[3]) override
  {
    output[0] = _T0 + _deltaT*((util::exp(_Pr*_Re*x[1] / _ly) - 1) / (util::exp(_Pr*_Re) - 1));
    return true;
  };
};

template <typename T, typename S>
class AnalyticalHeatFlux : public AnalyticalF3D<T, S>
{
  private:
    T _Re;
    T _Pr;
    T _deltaT;
    T _ly;
    T _lambda;
  public:
  AnalyticalHeatFlux(T Re, T Pr, T deltaT, T ly, T lambda) : AnalyticalF3D<T, S>(3),
    _Re(Re), _Pr(Pr), _deltaT(deltaT), _ly(ly), _lambda(lambda)
  {
    this->getName() = "AnalyticalHeatFlux";
  };

  bool operator()(T output[3], const S x[3]) override
  {
    output[0] = 0;
    output[1] = - _lambda * _Re * _Pr * _deltaT / _ly * (util::exp(_Pr * _Re * x[1] / _ly))/(util::exp(_Pr * _Re) - 1);
    output[2] = 0;
    return true;
  };
};

void computeError(MyCase& myCase)
{
  OstreamManager clout( std::cout, "computeError" );

  using T = MyCase::value_t_of<NavierStokes>;
  using NSEDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ADEDESCRIPTOR = MyCase::descriptor_t_of<Temperature>;

  auto& geometry = myCase.getGeometry();
  auto& NSElattice = myCase.getLattice(NavierStokes{});
  auto& ADElattice = myCase.getLattice(Temperature{});
  const auto& converter = NSElattice.getUnitConverter();

  const T Re = converter.getReynoldsNumber();
  const T Pr = converter.getPrandtlNumber();

  int input[1] = { };
  T result[1]  = { };

  auto indicatorF = geometry.getMaterialIndicator({1, 2, 3});

  T u_Re = Re * converter.getPhysViscosity() / converter.getCharPhysLength();
  AnalyticalVelocity<T,T> uSol(Re, converter.getCharPhysVelocity(), u_Re, converter.getCharPhysLength());
  SuperLatticePhysVelocity3D<T,NSEDESCRIPTOR> u(NSElattice,converter);
  SuperAbsoluteErrorL2Norm3D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, input);
  clout << "velocity-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm3D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, input);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;

  AnalyticalTemperature<T,T> tSol(Re, Pr, converter.getCharPhysLength(), converter.getCharPhysLowTemperature(), converter.getCharPhysTemperatureDifference());
  SuperLatticePhysTemperature3D<T,NSEDESCRIPTOR,ADEDESCRIPTOR> t(ADElattice,converter);
  SuperAbsoluteErrorL2Norm3D<T> absTemperatureErrorNormL2(t, tSol, indicatorF);
  absTemperatureErrorNormL2(result, input);
  clout << "temperature-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm3D<T> relTemperatureErrorNormL2(t, tSol, indicatorF);
  relTemperatureErrorNormL2(result, input);
  clout << "; temperature-L2-error(rel)=" << result[0] << std::endl;

  AnalyticalHeatFlux<T,T> HeatFluxSol(Re, Pr, converter.getCharPhysTemperatureDifference(), converter.getCharPhysLength(), converter.getThermalConductivity());
  SuperLatticePhysHeatFlux3D<T,NSEDESCRIPTOR,ADEDESCRIPTOR> HeatFlux(ADElattice,converter);
  SuperAbsoluteErrorL2Norm3D<T> absHeatFluxErrorNormL2(HeatFlux, HeatFluxSol, indicatorF);
  absHeatFluxErrorNormL2(result, input);
  clout << "heatFlux-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm3D<T> relHeatFluxErrorNormL2(HeatFlux, HeatFluxSol, indicatorF);
  relHeatFluxErrorNormL2(result, input);
  clout << "; heatFlux-L2-error(rel)=" << result[0] << std::endl;
}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters)
{
  using T             = MyCase::value_t_of<NavierStokes>;

  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const T physDeltaX  = parameters.get<parameters::PHYS_DELTA_X>();
  const Vector origin = parameters.get<parameters::DOMAIN_ORIGIN>();
  IndicatorCuboid3D<T> cuboid(extent, origin);

  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({true, false, true});
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t_of<NavierStokes>;
  auto& geometry  = myCase.getGeometry();
  auto& parameters    = myCase.getParameters();

  geometry.rename(0, 2);
  geometry.rename(2, 1, {0, 1, 0});

  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const T physDeltaX  = parameters.get<parameters::PHYS_DELTA_X>();
  const Vector origin = parameters.get<parameters::DOMAIN_ORIGIN>();
  const Vector bottomExtent{extent[0], physDeltaX, extent[2]};

  IndicatorCuboid3D<T> bottom(bottomExtent, origin);

  geometry.rename(2, 3, 1, bottom);

  geometry.clean();
  geometry.innerClean();
  geometry.checkForErrors();

  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase) {
  OstreamManager clout(std::cout,"prepareLattice");

  using T = MyCase::value_t_of<NavierStokes>;
  using NSEDESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  using ADEDESCRIPTOR = MyCase::descriptor_t_of<Temperature>;

  auto& geometry      = myCase.getGeometry();
  auto& NSElattice    = myCase.getLattice(NavierStokes{});
  auto& ADElattice    = myCase.getLattice(Temperature{});
  auto& parameters    = myCase.getParameters();

  const T Ra                              = parameters.get<parameters::RAYLEIGH>();
  const T Re                              = parameters.get<parameters::REYNOLDS>();
  const T Pr                              = parameters.get<parameters::PRANDTL>();
  const T physDeltaX                      = parameters.get<parameters::PHYS_DELTA_X>();
  const T g                               = parameters.get<parameters::GRAVITATIONAL_ACC>();
  const T physCharLength                  = parameters.get<parameters::PHYS_CHAR_LENGTH>();
  const T physCharVelocity                = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T physDensity                     = parameters.get<parameters::PHYS_CHAR_DENSITY>();
  const T Tcold                           = parameters.get<parameters::T_COLD>();
  const T Thot                            = parameters.get<parameters::T_HOT>();
  const T physThermalConductivity         = parameters.get<parameters::PHYS_THERMAL_CONDUCTIVITY>();
  const T tau                             = parameters.get<parameters::LATTICE_RELAXATION_TIME>();

  const T dT                              = Thot - Tcold;
  const T physViscosity                   = physCharLength * physCharVelocity / Re;
  const T physThermalDiffusivity          = physViscosity / Pr;
  const T physThermalExpansionCoefficient = physViscosity * physThermalDiffusivity * Ra / (g * dT * physCharLength * physCharLength * physCharLength);
  const T physSpecificHeatCapacity        = physThermalConductivity / (physDensity * physThermalDiffusivity);
  const T physDeltaT                      = (tau - (T) 1 / 2) / descriptors::invCs2<T,NSEDESCRIPTOR>() * physDeltaX * physDeltaX / physViscosity;

  NSElattice.setUnitConverter<ThermalUnitConverter<T,NSEDESCRIPTOR,ADEDESCRIPTOR>>(
    (T) physDeltaX,
    (T) physDeltaT,
    (T) physCharLength,
    (T) physCharVelocity,
    (T) physViscosity,
    (T) physDensity,
    (T) physThermalConductivity,
    (T) physSpecificHeatCapacity,
    (T) physThermalExpansionCoefficient,
    (T) Tcold,
    (T) Thot
  );
  const auto& converter = NSElattice.getUnitConverter();
  converter.print();

  ADElattice.setUnitConverter(converter);

  ADElattice.defineDynamics<AdvectionDiffusionBGKdynamics>(geometry.getMaterialIndicator({ 1, 2, 3 }));
  NSElattice.defineDynamics<ForcedBGKdynamics>(geometry.getMaterialIndicator({ 1, 2, 3 }));


  boundary::set<boundary::LocalVelocity>(NSElattice, geometry.getMaterialIndicator({ 2, 3 }));

  #ifdef AdvectionDiffusionDirichletBC
  boundary::set<boundary::AdvectionDiffusionDirichlet>(ADElattice, geometry.getMaterialIndicator({ 2, 3 }));
  #endif
  #ifdef RegularizedTemperatureBC
  boundary::set<boundary::RegularizedTemperature>(ADElattice, geometry.getMaterialIndicator({ 2, 3 }));
  #endif
  #ifdef RegularizedHeatFluxBC
    boundary::set<boundary::RegularizedHeatFlux>(ADElattice, geometry, 2);
    boundary::set<boundary::RegularizedTemperature>(ADElattice, geometry, 3);
  #endif

  ADElattice.setParameter<descriptors::OMEGA>(converter.getLatticeThermalRelaxationFrequency());
  NSElattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  const T boussinesqForcePrefactor = g
                                   / converter.getConversionFactorVelocity()
                                   * converter.getConversionFactorTime()
                                   * converter.getCharPhysTemperatureDifference()
                                   * converter.getPhysThermalExpansionCoefficient();

  auto& coupling = myCase.setCouplingOperator(
    "Boussinesq",
    NavierStokesAdvectionDiffusionCoupling{},
    names::NavierStokes{}, NSElattice,
    names::Temperature{},  ADElattice
  );
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::T0>(
    converter.getLatticeTemperature(Tcold)
  );
  coupling.setParameter<NavierStokesAdvectionDiffusionCoupling::FORCE_PREFACTOR>(
    boussinesqForcePrefactor * Vector<T,3>{ 0.0, 1.0, 0.0}
  );

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase)
{
  OstreamManager clout(std::cout, "setInitialValues");
  clout << "Set Initial Values ..." << std::endl;

  using T = MyCase::value_t_of<NavierStokes>;

  auto& geometry    = myCase.getGeometry();
  auto& NSElattice  = myCase.getLattice(NavierStokes{});
  auto& ADElattice  = myCase.getLattice(Temperature{});
  auto& converter   = NSElattice.getUnitConverter();

  const T Tcold     = converter.getCharPhysLowTemperature();

  Vector u_top(converter.getCharPhysVelocity(), converter.getCharPhysVelocity(), 0.0);
  Vector u_bot(0.0, converter.getCharPhysVelocity(), 0.0);

  momenta::setVelocity(NSElattice, geometry.getMaterialIndicator(3), u_bot);
  momenta::setVelocity(NSElattice, geometry.getMaterialIndicator(2), u_top);

  momenta::setTemperature(ADElattice, geometry.getMaterialIndicator( { 1, 3 }), Tcold);

  #ifdef RegularizedHeatFluxBC
    Vector<T,3> heatFlux;
    T input[3] = {0., 1., 0.};

    const T Re = converter.getReynoldsNumber();
    const T Pr = converter.getPrandtlNumber();

    AnalyticalHeatFlux<T,T> HeatFluxSol(Re, Pr, converter.getCharPhysTemperatureDifference(), converter.getCharPhysLength(), converter.getThermalConductivity());
    HeatFluxSol(heatFlux.data(), input);
    momenta::setHeatFlux(ADElattice, geometry.getMaterialIndicator(2), heatFlux);
  #else
    const T Thot      = converter.getCharPhysHighTemperature();
    momenta::setTemperature(ADElattice, geometry.getMaterialIndicator(2), Thot);
  #endif

  NSElattice.initialize();
  ADElattice.initialize();

  clout << "Set Initial Values ... OK" << std::endl;
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{
  // Nothing to do here
}

void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT,
                bool converged)
{
  OstreamManager clout(std::cout,"getResults");

  using T               = MyCase::value_t_of<NavierStokes>;
  using NSEDESCRIPTOR   = MyCase::descriptor_t_of<NavierStokes>;
  using ADEDESCRIPTOR   = MyCase::descriptor_t_of<Temperature>;

  auto& NSElattice      = myCase.getLattice(NavierStokes{});
  auto& ADElattice      = myCase.getLattice(Temperature{});
  const auto& converter = NSElattice.getUnitConverter();
  auto& parameters      = myCase.getParameters();

  const int statIter    = converter.getLatticeTime(parameters.get<parameters::PHYS_STAT_ITER_T>());
  const int vtkIter     = converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T>());

  if (iT == 0)
  {
    SuperVTMwriter3D<T> vtkWriter("porousPlate3d");

    const T Re            = converter.getReynoldsNumber();
    const T Pr            = converter.getPrandtlNumber();
    const T charL         = converter.getCharPhysLength();

    AnalyticalHeatFlux<T,T> HeatFluxSol(Re, Pr, converter.getCharPhysTemperatureDifference(), charL, converter.getThermalConductivity());

    T tmpIn[3] = {0., 1., 0.};
    T tmpOut[3];
    HeatFluxSol(tmpOut,tmpIn);
    clout << converter.getLatticeHeatFlux(tmpOut[0]) << " " << converter.getLatticeHeatFlux(tmpOut[1]) << std::endl;
    clout << tmpOut[0] << " " << tmpOut[1] << std::endl;

    SuperLatticeCuboid3D<T, NSEDESCRIPTOR> cuboid(NSElattice);
    SuperLatticeRank3D<T, NSEDESCRIPTOR> rank(NSElattice);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
  }

  if (iT % vtkIter == 0 || converged)
  {
    ADElattice.setProcessingContext(ProcessingContext::Evaluation);
    NSElattice.setProcessingContext(ProcessingContext::Evaluation);

    computeError(myCase);

    SuperVTMwriter3D<T> vtkWriter("porousPlate3d");

    SuperLatticePhysVelocity3D<T, NSEDESCRIPTOR>                    velocity(NSElattice, converter);
    SuperLatticePhysPressure3D<T, NSEDESCRIPTOR>                    pressure(NSElattice, converter);
    SuperLatticePhysTemperature3D<T, NSEDESCRIPTOR, ADEDESCRIPTOR>  temperature(ADElattice, converter);
    SuperLatticePhysHeatFlux3D<T, NSEDESCRIPTOR, ADEDESCRIPTOR>     heatflux(ADElattice, converter);

    const T Re            = converter.getReynoldsNumber();
    const T Pr            = converter.getPrandtlNumber();
    const T charL         = converter.getCharPhysLength();
    const T charU         = converter.getCharPhysVelocity();
    const T Tcold         = converter.getCharPhysLowTemperature();
    const T Thot          = converter.getCharPhysHighTemperature();

    AnalyticalVelocity<T,T>                         uSol(Re, charU, charU, charL);
    SuperLatticeFfromAnalyticalF3D<T,NSEDESCRIPTOR> uSolLattice(uSol, NSElattice);
    AnalyticalTemperature<T,T>                      TSol(Re, Pr, charL, Tcold, converter.getCharPhysTemperatureDifference());
    SuperLatticeFfromAnalyticalF3D<T,ADEDESCRIPTOR> TSolLattice(TSol,ADElattice);
    AnalyticalHeatFlux<T,T>                         HeatFluxSol(Re, Pr, converter.getCharPhysTemperatureDifference(), charL, converter.getThermalConductivity());
    SuperLatticeFfromAnalyticalF3D<T,ADEDESCRIPTOR> HeatFluxSolLattice(HeatFluxSol, ADElattice);

    vtkWriter.addFunctor( pressure );
    vtkWriter.addFunctor( velocity );
    vtkWriter.addFunctor( temperature );
    vtkWriter.addFunctor( heatflux );
    vtkWriter.addFunctor( uSolLattice );
    vtkWriter.addFunctor( TSolLattice );
    vtkWriter.addFunctor( HeatFluxSolLattice );

    vtkWriter.write(iT);

    BlockReduction3D2D<T> planeReduction( temperature, {0,0,1}, 600, BlockDataSyncMode::ReduceOnly );
    heatmap::plotParam<T> jpeg_Param;
    jpeg_Param.maxValue = Thot;
    jpeg_Param.minValue = Tcold;
    heatmap::write(planeReduction, iT, jpeg_Param);
  }

  if (iT % statIter == 0 || converged)
  {
    NSElattice.getStatistics().print(iT,converter.getPhysTime(iT));
    timer.print(iT);
  }
}

void simulate(MyCase& myCase)
{
  OstreamManager clout(std::cout,"simulate");
  clout << "Starting simulation" << std::endl;

  using T = MyCase::value_t;

  auto& NSElattice = myCase.getLattice(NavierStokes{});
  auto& ADElattice = myCase.getLattice(Temperature{});
  auto& coupling   = myCase.getOperator("Boussinesq");
  auto& converter  = NSElattice.getUnitConverter();
  auto& parameters = myCase.getParameters();

  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(
    parameters.get<parameters::PHYS_MAX_T>()
  );

  const T epsilon       = parameters.get<parameters::EPSILON>();
  const T convCheckTime = parameters.get<parameters::CONV_ITER>();
  util::ValueTracer<T> converge(converter.getLatticeTime(convCheckTime), epsilon);

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();
  for (std::size_t iT = 0; iT < iTmax; ++iT)
  {
    if (converge.hasConverged())
    {
      clout << "Simulation converged." << std::endl;
      getResults(myCase, timer, iT, converge.hasConverged());
      break;
    }

    setTemporalValues(myCase, iT);

    NSElattice.collideAndStream();
    coupling.apply();
    ADElattice.collideAndStream();

    getResults(myCase, timer, iT, converge.hasConverged());

    converge.takeValue(NSElattice.getStatistics().getAverageEnergy());
  }

  timer.stop();
  timer.printSummary();
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);

  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION               >(    20     );
    myCaseParameters.set<REYNOLDS                 >(     5     );
    myCaseParameters.set<RAYLEIGH                 >(   100     );
    myCaseParameters.set<PRANDTL                  >(     0.71  );
    myCaseParameters.set<DOMAIN_ORIGIN            >( {0, 0, 0} );
    myCaseParameters.set<PHYS_CHAR_LENGTH         >(     1.0   );
    myCaseParameters.set<PHYS_CHAR_VELOCITY       >(     1.0   );
    myCaseParameters.set<PHYS_CHAR_DENSITY        >(     1.0   );
    myCaseParameters.set<PHYS_THERMAL_CONDUCTIVITY>(     1.0   );
    myCaseParameters.set<LATTICE_RELAXATION_TIME  >(     1.0   );
    myCaseParameters.set<LATTICE_CHAR_VELOCITY    >(     0.1   );
    myCaseParameters.set<GRAVITATIONAL_ACC        >(     9.81  );
    myCaseParameters.set<PHYS_MAX_T               >( 10000     );
    myCaseParameters.set<T_HOT                    >(   274.15  );
    myCaseParameters.set<T_COLD                   >(   273.15  );
    myCaseParameters.set<EPSILON                  >(     1e-5  );
    myCaseParameters.set<CONV_ITER                >(     1.0   );
    myCaseParameters.set<PHYS_STAT_ITER_T         >(    10.0   );
    myCaseParameters.set<PHYS_VTK_ITER_T          >(   100.0   );
  }
  myCaseParameters.set<parameters::PHYS_DELTA_X>([] (MyCase::ParametersD& p) -> MyCase::value_t  {
    return {p.get<parameters::PHYS_CHAR_LENGTH>() / p.get<parameters::RESOLUTION>()};
  });
  myCaseParameters.set<parameters::DOMAIN_EXTENT>([] (MyCase::ParametersD& p) -> Vector<MyCase::value_t, 3>  {
    return {p.get<parameters::PHYS_CHAR_LENGTH>(),
            p.get<parameters::PHYS_CHAR_LENGTH>(),
            p.get<parameters::PHYS_DELTA_X>() * 3,
            };
  });


  Mesh mesh = createMesh(myCaseParameters);

  MyCase myCase(myCaseParameters, mesh);

  prepareGeometry(myCase);

  prepareLattice(myCase);

  setInitialValues(myCase);

  simulate(myCase);
}
