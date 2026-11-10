/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
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

/**
 * This example implements the model for predicting the dynamic
 * response of a mechanical heart valve described by Stijnen et al.
 * DOI: 10.1016/j.jfluidstructs.2004.04.007. 2004.
 **/

#define DISABLE_CSE

#include <olb.h>

using namespace olb;
using namespace olb::names;

// === Step 1: Declarations ===
using MyCase = Case<
  NavierStokes, Lattice<float, descriptors::D2Q9<
    descriptors::POROSITY,
    descriptors::VELOCITY,
    fields::fsi::ELEMENT_TAG,
    fields::fsi::ELEMENT_TORQUE,

    fields::moments::RHO,
    fields::moments::U
  >>
>;

namespace olb::parameters {

struct T_P : public descriptors::FIELD_BASE<1> { };
struct U_AMP : public descriptors::FIELD_BASE<1> { };
struct U_REF : public descriptors::FIELD_BASE<1> { };

}

std::shared_ptr<IndicatorF2D<MyCase::value_t>> createDomainI(MyCase::ParametersD& params) {
  using T = MyCase::value_t;

  const T channelH = params.get<parameters::DOMAIN_EXTENT>()[0];
  const T channelL = params.get<parameters::DOMAIN_EXTENT>()[1];

  const T cavityR  = channelH;
  const T cavityX  = channelH+cavityR;
  const T cavityY  = channelH;

  std::shared_ptr<IndicatorF2D<T>> channelI(new IndicatorCuboid2D<T>({channelL,channelH},{0,0}));
  std::shared_ptr<IndicatorF2D<T>> cavityI(new IndicatorCircle2D<T>({cavityX, cavityY}, cavityR));

  std::shared_ptr<IndicatorF2D<T>> cornerBoxI(new IndicatorCuboid2D<T>({0.00721568,0.00545155},{0.0592418,channelH}));
  std::shared_ptr<IndicatorF2D<T>> cornerCircleI(new IndicatorCircle2D<T>({0.0664575,0.0275},0.0075));
  auto cornerI = cornerBoxI - (cornerBoxI * cornerCircleI);

  return channelI+cavityI+cornerI;
}

std::shared_ptr<IndicatorF2D<MyCase::value_t>> createFsiRegionI(MyCase::ParametersD& params) {
  using T = MyCase::value_t;

  const T channelH = params.get<parameters::DOMAIN_EXTENT>()[0];

  const T cavityR  = channelH;

  std::shared_ptr<IndicatorF2D<T>> regionI(new IndicatorCuboid2D<T>({T{1.2}*cavityR,2*channelH},{channelH-T{0.1}*cavityR,0}));
  return regionI;
}


Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& params) {
  using T = MyCase::value_t;

  const T deltaX = params.get<parameters::DOMAIN_EXTENT>()[0] / params.get<parameters::RESOLUTION>();

  auto domainI = createDomainI(params);
  std::shared_ptr<IndicatorF2D<T>> paddedDomainI(
    new IndicatorLayer2D<T>(*domainI, deltaX));

  Mesh<T,MyCase::d> mesh(*paddedDomainI, deltaX, singleton::mpi().getSize());
  mesh.setOverlap(params.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().buildNeighborCache(params.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  using T = MyCase::value_t;

  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  auto& params = myCase.getParameters();

  auto domainI = createDomainI(params);
  auto fsiRegionI = createFsiRegionI(params);

  auto& sGeometry = myCase.getGeometry();

  sGeometry.rename(0, 3);
  sGeometry.rename(3, 1, domainI);
  sGeometry.rename(1, 2, fsiRegionI);

  const T deltaX = params.get<parameters::DOMAIN_EXTENT>()[0] / params.get<parameters::RESOLUTION>();
  const T channelH = params.get<parameters::DOMAIN_EXTENT>()[0];
  const T channelL = params.get<parameters::DOMAIN_EXTENT>()[1];

  {
    Vector<T,2> extend(channelL, channelH);
    Vector<T,2> origin;
    extend[0] = 2*deltaX;
    origin[0] = -deltaX;
    IndicatorCuboid2D<T> inflow(extend, origin);
    sGeometry.rename(3,4,1, inflow);
  }

  // Set material number for outflow
  {
    Vector<T,2> extend(channelL, channelH);
    Vector<T,2> origin;
    origin[0] = channelL-deltaX;
    IndicatorCuboid2D<T> outflow(extend, origin);
    sGeometry.rename(3,5,1, outflow);
  }

  sGeometry.clean(true, {1,2});
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

  const T channelH = params.get<parameters::DOMAIN_EXTENT>()[0];
  const T deltaX = channelH / params.get<parameters::RESOLUTION>();
  const T deltaT = 9.6e-5 * (T{40}/params.get<parameters::RESOLUTION>());

  const UnitConverter<T,DESCRIPTOR> converter(
    (T)   deltaX,           // deltaX
    (T)   deltaT,           // deltaT
    (T)   channelH,         // charPhysLength: reference length of simulation geometry
    (T)   0.04,             // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   4e-6,             // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1090              // physDensity: physical density in __kg / m^3__
  );
  converter.print();
  sLattice.setUnitConverter(converter);

  using BulkDynamics = dynamics::Tuple<
    T, DESCRIPTOR,
    momenta::BulkTuple,
    equilibria::SecondOrder,
    collision::BGK,
    forcing::fsi::HLBM
  >;

  sLattice.defineDynamics<NoDynamics>(sGeometry, 0);
  sLattice.defineDynamics<BulkDynamics>(sGeometry, 1);
  sLattice.defineDynamics<BulkDynamics>(sGeometry, 2);

  boundary::set<boundary::BounceBack>(sLattice, sGeometry, 3);
  boundary::set<boundary::InterpolatedVelocity>(sLattice, sGeometry, 4);
  boundary::set<boundary::InterpolatedPressure>(sLattice, sGeometry, 5);

  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  sLattice.writeSummary();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;

  auto& sGeometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes());

  fields::set<descriptors::POROSITY>(sLattice, sGeometry.getMaterialIndicator({1,2}), 1);
  fields::set<descriptors::POROSITY>(sLattice, sGeometry.getMaterialIndicator({0,3}), 0);

  // Track min/max/mean of the HLBM fields at statistics time
  sLattice.getStatistics().template registerField<descriptors::POROSITY>("POROSITY");
  sLattice.getStatistics().template registerField<descriptors::VELOCITY>("VELOCITY");

  sLattice.initialize();
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{
  using T = MyCase::value_t;

  auto& params = myCase.getParameters();
  const T T_p = params.get<parameters::T_P>();
  const T u_amp = params.get<parameters::U_AMP>();
  const T u_ref = params.get<parameters::U_REF>();

  auto& sGeometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes());
  const auto& converter = sLattice.getUnitConverter();

  const std::size_t iTupdate = converter.getLatticeTime(T_p/1000);
  const std::size_t iTstart = converter.getLatticeTime(T_p);

  if (iT % iTupdate == 0 && iT < iTstart) {
    T frac{};
    PolynomialStartScale<T,std::size_t>(iTstart, T{1})(&frac, &iT);
    momenta::setVelocity(sLattice, sGeometry.getMaterialIndicator(4), Vector{frac*u_ref, 0});

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
  else if (iT % iTupdate == 0 && iT >= iTstart) {
    T u { };
    const T t = std::fmod(converter.getPhysTime(iT), T_p);
    if (t / T_p <= 0.37) {
      u = u_ref + u_amp * util::sin(((2*std::numbers::pi_v<T>*t)/T_p)/0.74);
    } else if (t / T_p > 0.37 && t / T_p <= 1) {
      u = u_ref + 0.5 * u_amp * util::sin(2*std::numbers::pi_v<T>*(t/T_p+0.26)/1.26);
    }

    momenta::setVelocity(sLattice, sGeometry.getMaterialIndicator(4), Vector{u, 0});
    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

void getResults(MyCase& myCase,
                std::size_t iT)
{
  OstreamManager clout(std::cout, "getResults");

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

  auto& params = myCase.getParameters();

  auto& sLattice = myCase.getLattice(NavierStokes());
  const auto& converter = sLattice.getUnitConverter();

  const int vtkIter  = converter.getLatticeTime(0.10);

  const bool vtkEnabled = params.get<parameters::VTK_ENABLED>();

  if (iT == 0 && vtkEnabled) {
    {
      SuperVTMwriter2D<T> vtmWriter("flow");
      vtmWriter.createMasterFile();
    }
    {
      SuperVTMwriter2D<T> vtmWriter("structure");
      vtmWriter.createMasterFile();
    }
  }

  if (iT % vtkIter == 0 && vtkEnabled) {
    sLattice.executePostProcessors(stage::Evaluation{});
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    {
      SuperVTMwriter2D<T> vtkWriter("flow");
      SuperLatticePhysVelocity2D velocity(sLattice, converter);
      SuperLatticePhysPressure2D pressure(sLattice, converter);
      vtkWriter.addFunctor(pressure);
      vtkWriter.addFunctor(velocity);
      vtkWriter.write(iT);
    }
    {
      SuperVTMwriter2D<T> vtkWriter("structure");
      SuperLatticeField2D<T,DESCRIPTOR,descriptors::POROSITY> porosity(sLattice);
      SuperLatticeField2D<T,DESCRIPTOR,fields::fsi::ELEMENT_TAG> tag(sLattice);
      SuperLatticeField2D<T,DESCRIPTOR,fields::fsi::ELEMENT_FORCE> force(sLattice);
      SuperLatticeField2D<T,DESCRIPTOR,fields::fsi::ELEMENT_TORQUE> torque(sLattice);
      vtkWriter.addFunctor(porosity);
      vtkWriter.addFunctor(tag);
      vtkWriter.addFunctor(force);
      vtkWriter.addFunctor(torque);
      vtkWriter.write(iT);
    }
  }
}

void simulate(MyCase& myCase) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

  auto& params = myCase.getParameters();
  auto& sGeometry = myCase.getGeometry();
  auto& sLattice = myCase.getLattice(NavierStokes());
  const auto& converter = sLattice.getUnitConverter();

  SuperPorousElementEmbeddingO fsiEmbeddingO(sLattice);
  fsiEmbeddingO.registerElementType<CuboidPorosityF>(sGeometry.getMaterialIndicator(2));

  const T channelH = params.get<parameters::DOMAIN_EXTENT>()[0];
  const T valveL = 1.07*channelH;
  const T valveMomentaOfInertia = 1100*util::pow(valveL,3)*(0.001)/3;
  const T T_p = params.get<parameters::T_P>();

  T valveW = 0.5*converter.getPhysDeltaX();

  T angle = 0;
  T angularU = 0.0;
  Vector lower{channelH, channelH+T{0.5}*valveW};

  auto valveParams = makeParametersD<T,descriptors::SPATIAL_DESCRIPTOR<2>>(
    fields::fsi::ELEMENT_TAG{},      1,
    fields::fsi::ELEMENT_LOWER{},    lower,
    fields::fsi::ELEMENT_UPPER{},    lower + Vector{valveL,valveW},
    fields::fsi::ELEMENT_PIVOT{},    lower + Vector{valveW/2,valveW/2},
    fields::fsi::ELEMENT_ROTATION{}, angle
  );
  fsiEmbeddingO.add(valveParams);

  SuperPorousElementReductionO<T,DESCRIPTOR,fields::fsi::ELEMENT_TORQUE> fsiReductionO(
    sLattice,
    sGeometry.getMaterialIndicator(2));
  fsiReductionO.resize(1);
  fsiReductionO.addCollectionO(meta::id<CollectPorousBoundaryTorqueO>{});

  sLattice.setParameter<fields::converter::PHYS_LENGTH>(
    1/converter.getConversionFactorLength());
  sLattice.setParameter<fields::converter::PHYS_VELOCITY>(
    1/converter.getConversionFactorVelocity());
  sLattice.setParameter<fields::converter::PHYS_DELTA_X>(
    converter.getPhysDeltaX());

  sLattice.setProcessingContext(ProcessingContext::Simulation);
  fsiEmbeddingO.initialize();

  CSV<T> valveWriter("valve_n" + std::to_string(params.get<parameters::RESOLUTION>()),
                     ';', {"time", "angle"}, ".csv");
  OstreamManager clresult(std::cout, "result");

  const std::size_t iTmax = converter.getLatticeTime(
    params.get<parameters::PHYS_MAX_T>());

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    myCase.getLattice(NavierStokes{}).collideAndStream();

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, iT);
    myCase.report(iT);

    fsiReductionO.apply();

    if (fsiReductionO.rankDoesFSI()) {
      for (unsigned iElement=1; iElement <= fsiReductionO.getElementCount(); ++iElement) {
        // Divide by deltaX due to 3D mass conversion factor choice in UnitConverter
        auto t = converter.getPhysTorque(fsiReductionO.getField<fields::fsi::ELEMENT_TORQUE>(iElement))
               / converter.getPhysDeltaX();

        auto angularAccel = t / valveMomentaOfInertia;
        if (iT >= converter.getLatticeTime(T_p)) {
          angle += angularU * converter.getPhysDeltaT();
          angularU += angularAccel * converter.getPhysDeltaT();
        }

        fsiEmbeddingO.setField<fields::fsi::ELEMENT_ROTATION>(iElement-1, angle);
        fsiEmbeddingO.setField<fields::fsi::ELEMENT_U_ROTATION>(iElement-1, angularU);

        if (iT % converter.getLatticeTime(0.01) == 0) {
          clresult << "t=" << converter.getPhysTime(iT) << "; torque=" << t << "; angle=" << angle << "; angularU=" << angularU << "; angularA=" << angularAccel << std::endl;
          valveWriter.writeDataFile(converter.getPhysTime(iT), T{angle});
        }
      }

      fsiEmbeddingO.setProcessingContext<Array<fields::fsi::ELEMENT_ROTATION>>(ProcessingContext::Simulation);
      fsiEmbeddingO.setProcessingContext<Array<fields::fsi::ELEMENT_U_ROTATION>>(ProcessingContext::Simulation);
    }

    fsiEmbeddingO.apply();
  }

  sLattice.setProcessingContext(ProcessingContext::Evaluation);

  // Final full statistics report
  myCase.finalize(iTmax);
}

int main(int argc, char* argv[]) {
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    using T = MyCase::value_t;
    myCaseParameters.set<RESOLUTION   >(40);
    myCaseParameters.set<DOMAIN_L     >(0.02);
    myCaseParameters.set<DOMAIN_EXTENT>([](MyCase::ParametersD& p) {
      return Vector { p.get<DOMAIN_L>(), p.get<DOMAIN_L>()*T{6.0} };
    });
    myCaseParameters.set<T_P        >(2.45);
    myCaseParameters.set<U_AMP      >(0.11);
    myCaseParameters["reporter"]["statistics"].set<PHYS_ITER_T>(0.01);
    myCaseParameters["reporter"]["timer"].set<PHYS_ITER_T>(0.01);
    myCaseParameters.set<U_REF      >(0.04);
    myCaseParameters.set<PHYS_MAX_T >(myCaseParameters.get<T_P>() * 10);
    myCaseParameters.set<VTK_ENABLED>(true);
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
