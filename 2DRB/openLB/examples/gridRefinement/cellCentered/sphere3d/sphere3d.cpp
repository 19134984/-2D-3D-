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

#include <olb.h>

using namespace olb;
using namespace olb::names;

using MyCase = Case<
  NavierStokes, RefinedLattice<float, descriptors::D3Q19<>>
>;

namespace olb::parameters {

struct RADIUS_SPHERE : public descriptors::FIELD_BASE<1> {};

}

RefinedMesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& params)
{
  if (auto meshFile = singleton::input().getCaseFile()) {
    return *meshFile;
  } else {
    throw std::runtime_error("No case file given, use \"--case-file\".");
  }
}

auto makeConverter(MyCase::ParametersD& params)
{
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  const int resolution = params.get<parameters::RESOLUTION>();
  const T   Re         = params.get<parameters::REYNOLDS>();
  const T   charL      = params.get<parameters::PHYS_CHAR_LENGTH>();
  return UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>(
    resolution,      // resolution: number of voxels per charPhysL
    (T) 0.505,       // latticeRelaxationTime
    charL,           // charPhysLength
    charL,           // charPhysVelocity
    charL*charL/Re,  // physViscosity
    (T) 1.0          // physDensity
  );
}

void prepareGeometry(MyCase& myCase)
{
  using T        = MyCase::value_t;
  auto& params   = myCase.getParameters();

  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  T physDeltaX = params.get<parameters::PHYS_CHAR_LENGTH>() / params.get<parameters::RESOLUTION>();
  Vector<T, 3> extent = params.get<parameters::DOMAIN_EXTENT>();
  Vector<T, 3> origin = params.get<parameters::ORIGIN>();

  auto& refinedGeometry = myCase.getGeometry();
  const std::size_t levelCount = refinedGeometry.getLevelCount();

  for (std::size_t iLevel = 0; iLevel < levelCount; ++iLevel) {
    auto& sGeometry = refinedGeometry.getLevel(iLevel);

    IndicatorCuboid3D<T> domainI(extent, origin);
    IndicatorLayer3D<T>  domainLayerI(domainI, physDeltaX);

    sGeometry.rename(0, 2, domainLayerI);
    sGeometry.rename(2, 1, domainI);

    // Set material number for inflow
    Vector<T,3> inflowOrigin = origin;
    inflowOrigin[0] -= physDeltaX;
    Vector<T,3> inflowExtent = extent;
    inflowExtent[0] = 2*physDeltaX;
    IndicatorCuboid3D<T> inflow(inflowExtent, inflowOrigin);
    sGeometry.rename(2, 3, 1, inflow, {1, 0, 0});

    // Set material number for outflow
    Vector<T,3> outflowOrigin = origin;
    outflowOrigin[0] = origin[0] + extent[0] - physDeltaX;
    Vector<T,3> outflowExtent = extent;
    outflowExtent[0] = 2*physDeltaX;
    IndicatorCuboid3D<T> outflow(outflowExtent, outflowOrigin);
    sGeometry.rename(2, 4, 1, outflow, {-1, 0, 0});

    // Set material number for sphere
    IndicatorSphere3D<T> sphereI({0.55, 0.25, 0.25}, params.get<parameters::RADIUS_SPHERE>());
    sGeometry.rename(1,5, sphereI);
  }

  refinedGeometry.markCoarseBorders(11); // inactive coarse-border blocks

  for (std::size_t iLevel = 0; iLevel < levelCount; ++iLevel) {
    auto& sGeometry = refinedGeometry.getLevel(iLevel);
    sGeometry.getStatistics().print();
    sGeometry.disableStatisticsAndCommunication();
  }

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase)
{
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  auto& refinedLattice = myCase.getLattice(names::NavierStokes{});
  auto& refinedGeometry = myCase.getGeometry();
  const std::size_t levelCount = refinedGeometry.getLevelCount();

  for (std::size_t iLevel = 0; iLevel < levelCount; ++iLevel) {
    auto& sGeometry = refinedGeometry.getLevel(iLevel);
    auto& sLattice = refinedLattice.getLevel(iLevel);

    const T omega = sLattice.getUnitConverter().getLatticeRelaxationFrequency();

    // Material=1 -->bulk dynamics
    auto bulkIndicator = sGeometry.getMaterialIndicator({1});

    using BulkDynamics = dynamics::Tuple<
      T, DESCRIPTOR,
      momenta::BulkTuple,
      equilibria::ThirdOrder,
      collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB>
    >;
    dynamics::set<BulkDynamics>(sLattice, bulkIndicator);

    // Material=2 -->bounce back (domain boundaries if any)
    boundary::set<boundary::BounceBack>(sLattice, sGeometry, 2);

    // Material=3 -->Inflow
    boundary::set<boundary::InterpolatedVelocity>(sLattice, sGeometry, 3);
    // Material=4 -->Outflow
    boundary::set<boundary::InterpolatedPressure>(sLattice, sGeometry, 4);

    // Material=5 -->BounceBack for sphere
    boundary::set<boundary::BounceBack>(sLattice, sGeometry, 5);

    sLattice.setParameter<descriptors::OMEGA>(omega);
    sLattice.setParameter<collision::LES::SMAGORINSKY>(0.18);
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase)
{
  auto& refinedLattice = myCase.getLattice(names::NavierStokes{});
  const std::size_t levelCount = refinedLattice.getLevelCount();

  for (std::size_t iLevel = 0; iLevel < levelCount; ++iLevel) {
    auto& sLattice = refinedLattice.getLevel(iLevel);
    sLattice.initialize();
  }
}

void simulate(MyCase& myCase)
{
  using T           = MyCase::value_t;
  using DESCRIPTOR  = MyCase::descriptor_t_of<NavierStokes>;
  auto&   params    = myCase.getParameters();
  auto&   refinedLattice = myCase.getLattice(names::NavierStokes{});
  auto&   refinedGeometry = myCase.getGeometry();
  auto&   mesh      = myCase.getMesh();

  auto&   converter0 = refinedLattice.getLevel(0).getUnitConverter();
  const T physMaxT  = params.get<parameters::PHYS_MAX_T>();

  const std::size_t iTmax      = converter0.getLatticeTime(physMaxT);
  const std::size_t iTmaxStart = converter0.getLatticeTime(1.0);
  const std::size_t iTupdate   = 5;

  using SliceF = functors::pointwise::Composed<functors::pointwise::VorticityMagnitudeF,
                                               functors::pointwise::PhysVelocityF>;
  const T dxFine = refinedLattice.getLevel(refinedLattice.getLevelCount()-1).getUnitConverter().getPhysDeltaX();

  // Cross-level mid-plane as adaptive VTU; single-level uniform VTIs at their own resolution: a
  // finest window at the sphere and the coarse inlet plane
  Hyperplane3D<T> midPlane;
  midPlane.originAt({0, 0, 0.25}).normalTo({0, 0, 1});
  refinement::PointExtraction<T,DESCRIPTOR,SliceF> midSlice(refinedLattice, mesh,
    refinement::sampling::slice(refinedLattice, mesh, midPlane));
  auto fineLattice  = sampling::planeLattice<T>({0.46, 0.16, 0.25}, {0.18, 0.18, 0}, dxFine);
  auto inletLattice = sampling::planeLattice<T>({0.05, 0, 0}, {0, 0.5, 0.5}, converter0.getPhysDeltaX());
  refinement::PointExtraction<T,DESCRIPTOR,SliceF> fineSlice (refinedLattice, mesh, fineLattice);
  refinement::PointExtraction<T,DESCRIPTOR,SliceF> inletSlice(refinedLattice, mesh, inletLattice);

  for (auto* s : {&midSlice, &fineSlice, &inletSlice}) {
    s->template setParameter<fields::converter::PHYS_VELOCITY>(converter0.getConversionFactorVelocity());
  }
  midSlice.template setParameter<fields::converter::PHYS_DELTA_X>(converter0.getPhysDeltaX());
  fineSlice.template setParameter<fields::converter::PHYS_DELTA_X>(dxFine);
  inletSlice.template setParameter<fields::converter::PHYS_DELTA_X>(converter0.getPhysDeltaX());

  auto midFields = io::pointFieldsFor(midSlice);
  midFields.add("Level", 1, [&](std::size_t i, unsigned){ return T(midSlice.getPointLevel(i)); });
  auto fineFields  = io::pointFieldsFor(fineSlice);
  auto inletFields = io::pointFieldsFor(inletSlice);
  const std::size_t sliceIter = util::max(converter0.getLatticeTime(0.1), std::size_t{1});

  auto preStepHook = [&](unsigned iLevel, std::size_t iT) {
    if (iLevel == 0 && iT % iTupdate == 0 && iT <= iTmaxStart) {
      auto& sLattice = refinedLattice.getLevel(0);
      auto& sGeometry = refinedGeometry.getLevel(0);

      T frac = 0;
      PolynomialStartScale<T,std::size_t>(iTmaxStart, 1)(&frac, &iT);
      std::vector<T> maxVelocity(3,0);
      maxVelocity[0] = 2.25*frac*converter0.getCharPhysVelocity();
      const T distance2Wall = converter0.getPhysDeltaX()/2.;
      RectanglePoiseuille3D<T> poiseuilleU(sGeometry, 3, maxVelocity,
                                           distance2Wall, distance2Wall, distance2Wall);
      momenta::setVelocity(sLattice, sGeometry.getMaterialIndicator(3), poiseuilleU);

      sLattice.template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
        ProcessingContext::Simulation);
    }
  };

  for (std::size_t iT = 0; iT < iTmax; ++iT) {
    refinement::rohde::executeStep(refinedLattice, 0, iT, preStepHook);
    myCase.report(iT);

    if (iT % sliceIter == 0) {
      const std::string dir = singleton::directories().getVtkOutDir();
      midSlice.update();
      io::writeVtu(midSlice, midFields, dir + "slice", iT);
      fineSlice.update();
      io::writeVti(fineSlice, fineFields, io::regularGrid(fineLattice), dir + "fineSlice", iT);
      inletSlice.update();
      io::writeVti(inletSlice, inletFields, io::regularGrid(inletLattice), dir + "inletSlice", iT);
    }
  }

  myCase.finalize(iTmax);
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");

  using T           = MyCase::value_t;

  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<PHYS_CHAR_LENGTH>(0.2);
    myCaseParameters.set<RESOLUTION>(16);
    myCaseParameters.set<REYNOLDS>(100);
    myCaseParameters.set<PHYS_MAX_T>(16);
    myCaseParameters.set<RADIUS_SPHERE>(0.1);
    myCaseParameters.set<ORIGIN>(Vector<T, 3>{0.0, 0.0, 0.0});
    myCaseParameters.set<DOMAIN_EXTENT>(Vector<T, 3>{2.5, 0.5, 0.5});
    myCaseParameters["reporter"]["statistics"].set<PHYS_ITER_T>(0.5);
    myCaseParameters["reporter"]["timer"].set<PHYS_ITER_T>(0.01);
    myCaseParameters["reporter"]["vtk"].set<PHYS_ITER_T>(0.1);
  }

  RefinedMesh<MyCase::value_t, MyCase::d> mesh = createMesh(myCaseParameters);

  MyCase myCase(myCaseParameters, mesh);

  auto& refinedLattice = myCase.getLattice(names::NavierStokes{});
  auto& refinedGeometry = myCase.getGeometry();

  auto converter = makeConverter(myCaseParameters);
  refinedLattice.setUnitConverter(converter);

  prepareGeometry(myCase);
  prepareLattice(myCase);
  setInitialValues(myCase);

  refinedLattice.initialize(refinedGeometry);

  auto vtkFunctors = [&myCase](auto& writer, std::size_t iLevel, auto& sLattice, std::size_t iT) {
    auto& converter = sLattice.getUnitConverter();
    SuperLatticePhysVelocity3D velocityF(sLattice, converter);
    SuperLatticePhysPressure3D pressureF(sLattice, converter);
    SuperLatticeRefinementMetricKnudsen3D qualityF(sLattice, converter);
    SuperGeometryF<MyCase::value_t, 3> geometryF(myCase.getGeometry().getLevel(iLevel));
    if (iT == 0) {
      writer.addFunctor(geometryF);
    }
    writer.addFunctor(velocityF);
    writer.addFunctor(pressureF);
    writer.addFunctor(qualityF);
    writer.write();
  };
  refinement::setVTKreporter(myCase, "vtk", names::NavierStokes{}, "sphere3d", vtkFunctors);

  simulate(myCase);
}
