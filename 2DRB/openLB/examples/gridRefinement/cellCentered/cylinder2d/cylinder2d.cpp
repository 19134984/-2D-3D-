/*  Lattice Boltzmann sample, written in C++, using the OpenLB library
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
#include <refinement/octree.h>

using namespace olb;
using namespace olb::names;

using MyCase = Case<
  NavierStokes, RefinedLattice<FLOATING_POINT_TYPE, descriptors::D2Q9<>>
>;

namespace olb::parameters {

struct RADIUS_CYLINDER : public descriptors::FIELD_BASE<1> {};

struct MAX_REFINEMENT : public descriptors::TYPED_FIELD_BASE<int,1> {};
struct REFINEMENT_BLOCK_SIZE : public descriptors::TYPED_FIELD_BASE<int,1> {};

}

RefinedMesh<MyCase::value_t, MyCase::d> createMesh(MyCase::ParametersD& params)
{
  using T                     = MyCase::value_t;
  T L = 0.1 / params.get<parameters::RESOLUTION>();
  Vector<T, 2>         extent(2.2, 0.41 + L);
  std::vector<T>       origin(2, T());
  IndicatorCuboid2D<T> cuboid(extent, origin);

  const T            physDeltaX = L;

  const int maxRef = params.get<parameters::MAX_REFINEMENT>();
  const int refBlockSize = params.get<parameters::REFINEMENT_BLOCK_SIZE>();

  Vector<T, 2>   center(0.2, 0.2 + L / 2.0);
  T              radius = params.get<parameters::RADIUS_CYLINDER>();

  IndicatorCircle2D<T> circle(center, radius);
  IndicatorLayer2D<T> paddedCircle(circle, physDeltaX);
  IndicatorCuboid2D<T> refineCuboidI({1.0, 0.25}, {0.08,0.08});

  OctreeRefinedMeshBlueprint<T,2> blueprint(cuboid, physDeltaX, refBlockSize, maxRef);
  blueprint.refine(refineCuboidI, 1, 0);
  blueprint.refine(paddedCircle, 2, 2);

  blueprint.balance();
  blueprint.markBorders();
  blueprint.reduceWeights(circle);
  blueprint.loadBalance(singleton::mpi().getSize());
  blueprint.coarse();

  return RefinedMesh(blueprint);
}

auto makeConverter(MyCase::ParametersD& params)
{
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  const int resolution = params.get<parameters::RESOLUTION>();
  const T   Re         = params.get<parameters::REYNOLDS>();
  const T   radius     = params.get<parameters::RADIUS_CYLINDER>();
  return UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR>(
    resolution,             // resolution: number of voxels per charPhysL
    (T) 0.51,               // latticeRelaxationTime
    (T) 2.0*radius,         // charPhysLength
    (T) 0.2,                // charPhysVelocity
    (T) 0.2*2.*radius/Re,   // physViscosity
    (T) 1.0                 // physDensity
  );
}

void prepareGeometry(MyCase& myCase, SuperGeometry<MyCase::value_t, 2>& sGeometry)
{
  using T        = MyCase::value_t;
  auto& params   = myCase.getParameters();

  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  T L = 0.1 / params.get<parameters::RESOLUTION>();
  Vector<T, 2> extend(2.2, 0.41 + L);
  Vector<T, 2> origin;

  IndicatorCuboid2D<T> domain(extend, origin);
  sGeometry.rename(0,2,domain);
  sGeometry.rename(2,1,{1,1});

  Vector<T,2> innerExtent(extend);
  Vector<T,2> innerOrigin(origin);
  innerOrigin += L*2;
  innerExtent -= L*4;
  IndicatorCuboid2D<T> innerCuboid(innerExtent, innerOrigin);
  sGeometry.rename(2,1, innerCuboid);

  // Set material number for inflow
  extend[0] = 2.*L;
  origin[0] = -L;
  IndicatorCuboid2D<T> inflow( extend, origin );
  sGeometry.rename( 2,3,1,inflow );

  // Set material number for outflow
  origin[0] = 2.2 - L;
  IndicatorCuboid2D<T> outflow( extend, origin );
  sGeometry.rename( 2,4,1,outflow );

  // Set material number for cylinder
  Vector<T, 2> center(0.2, 0.2 + L / 2.0);
  T radius = params.get<parameters::RADIUS_CYLINDER>();
  IndicatorCircle2D<T> circle(center, radius);

  sGeometry.rename( 1,5, circle );

  sGeometry.clean();
  sGeometry.checkForErrors();

  sGeometry.getStatistics().print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase,
                    SuperGeometry<MyCase::value_t, 2>& sGeometry,
                    SuperLattice<MyCase::value_t,MyCase::descriptor_t_of<names::NavierStokes>>& sLattice)
{
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<names::NavierStokes>;
  auto& params     = myCase.getParameters();

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = sLattice.getUnitConverter().getLatticeRelaxationFrequency();

  // Material=1 -->bulk dynamics
  auto bulkIndicator = sGeometry.getMaterialIndicator({1});
  dynamics::set<BGKdynamics>(sLattice, bulkIndicator);

  // Material=2 -->bounce back
  boundary::set<boundary::BounceBack>(sLattice, sGeometry, 2);

  boundary::set<boundary::InterpolatedVelocity>(sLattice, sGeometry, 3);
  boundary::set<boundary::InterpolatedPressure>(sLattice, sGeometry, 4);

  T L = 0.1 / params.get<parameters::RESOLUTION>();
  Vector<T, 2> center(0.2, 0.2 + L / 2.0);
  T radius = params.get<parameters::RADIUS_CYLINDER>();
  IndicatorCircle2D<T> circle(center, radius);

  // Material=5 -->bouzidi / bounce back
  #ifdef BOUZIDI
  setBouzidiBoundary(sLattice, sGeometry, 5, circle);
  #else
  boundary::set<boundary::BounceBack>(sLattice, sGeometry, 5);
  #endif

  sLattice.setParameter<descriptors::OMEGA>(omega);

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase, SuperGeometry<MyCase::value_t, 2>& sGeometry, SuperLattice<MyCase::value_t,MyCase::descriptor_t_of<names::NavierStokes>>& sLattice)
{
  sLattice.initialize();
}

void getResults(MyCase& myCase, std::size_t iT)
{
  using T          = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<names::NavierStokes>;
  auto& refinedLattice = myCase.getLattice(names::NavierStokes{});

  const T physVtkInterval = myCase.getParameters()["reporter"]["vtk"].template get<parameters::PHYS_ITER_T>();
  const int vtkIter = refinedLattice.getLevel(0).getUnitConverter().getLatticeTime(physVtkInterval);
  if (vtkIter <= 0) {
    return;
  }

  if (iT == 0) {
    for (std::size_t iLevel = 0; iLevel < refinedLattice.getLevelCount(); iLevel++) {
      auto name = "level" + std::to_string(iLevel);
      SuperVTMwriter2D<T> vtmWriter(name, 1);
      if (iLevel == 0) {
        SuperLatticeRank2D rank( refinedLattice.getLevel(0) );
        vtmWriter.write(rank);
      }
      SuperGeometryF<T,DESCRIPTOR::d> geometryF(myCase.getGeometry().getLevel(iLevel));
      geometryF.getName() = name + "_geometry";
      vtmWriter.write(geometryF);
      SuperLatticeCuboid2D<T,DESCRIPTOR> cuboidF(refinedLattice.getLevel(iLevel));
      cuboidF.getName() = name + "_cuboid";
      vtmWriter.write(cuboidF);
      vtmWriter.createMasterFile();
    }
  }

  if (iT % vtkIter == 0) {
    for (std::size_t iLevel = 0; iLevel < refinedLattice.getLevelCount(); iLevel++) {
      auto& sLattice = refinedLattice.getLevel(iLevel);
      auto& converter = sLattice.getUnitConverter();
      sLattice.setProcessingContext(ProcessingContext::Evaluation);

      auto name = "level" + std::to_string(iLevel);
      SuperVTMwriter2D<T> vtmWriter(name);
      SuperLatticePhysVelocity2D velocityF(sLattice, converter);
      SuperLatticePhysPressure2D pressureF(sLattice, converter);
      SuperLatticeRefinementMetricKnudsen2D qualityF(sLattice, converter);
      vtmWriter.addFunctor(qualityF);
      vtmWriter.addFunctor(velocityF);
      vtmWriter.addFunctor(pressureF);
      vtmWriter.write(iT);
    }
  }
}

void simulate(MyCase& myCase)
{
  using T           = MyCase::value_t;
  auto&   params    = myCase.getParameters();
  auto&   refinedLattice = myCase.getLattice(names::NavierStokes{});
  auto&   refinedGeometry = myCase.getGeometry();

  auto&   converter0 = refinedLattice.getLevel(0).getUnitConverter();
  const T physMaxT  = params.get<parameters::PHYS_MAX_T>();

  const std::size_t iTmax      = converter0.getLatticeTime(physMaxT);
  const std::size_t iTmaxStart = converter0.getLatticeTime(6.4);
  const std::size_t iTupdate   = 5;

  auto preStepHook = [&](unsigned iLevel, std::size_t iT) {
    if (iLevel == 0 && iT % iTupdate == 0 && iT <= iTmaxStart) {
      auto& sLattice = refinedLattice.getLevel(0);
      auto& sGeometry = refinedGeometry.getLevel(0);

      PolynomialStartScale<T,T> StartScale( iTmaxStart, T( 1 ) );
      T iTvec[1] = {T( iT )};
      T frac[1] = {};
      StartScale( frac,iTvec );
      T maxVelocity = converter0.getCharPhysVelocity()*frac[0]*3./2.;
      T distance2Wall = converter0.getPhysDeltaX()/2.;
      Poiseuille2D<T> poiseuilleU( sGeometry, 3, maxVelocity, distance2Wall );

      momenta::setVelocity(sLattice, sGeometry.getMaterialIndicator(3), poiseuilleU);

      sLattice.template setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
        ProcessingContext::Simulation);
    }
  };

  for (std::size_t iT = 0; iT < iTmax; ++iT) {
    refinement::rohde::executeStep(refinedLattice, 0, iT, preStepHook);
    getResults(myCase, iT);
    myCase.report(iT);
  }

  myCase.finalize(iTmax);
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");

  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION>(11);
    myCaseParameters.set<REYNOLDS>(100);
    myCaseParameters.set<PHYS_MAX_T>(16);
    myCaseParameters.set<RADIUS_CYLINDER>(0.05);
    myCaseParameters.set<MAX_REFINEMENT>(2);
    myCaseParameters.set<REFINEMENT_BLOCK_SIZE>(8);
    myCaseParameters["reporter"]["statistics"].set<PHYS_ITER_T>(0.1);
    myCaseParameters["reporter"]["timer"].set<PHYS_ITER_T>(0.1);
    myCaseParameters["reporter"]["vtk"].set<PHYS_ITER_T>(1.0);
  }

  RefinedMesh<MyCase::value_t, MyCase::d> mesh = createMesh(myCaseParameters);

  MyCase myCase(myCaseParameters, mesh);

  auto& refinedLattice = myCase.getLattice(names::NavierStokes{});
  auto& refinedGeometry = myCase.getGeometry();

  auto converter = makeConverter(myCaseParameters);
  refinedLattice.setUnitConverter(converter);

  for (std::size_t iLevel = 0; iLevel < refinedLattice.getLevelCount(); iLevel++) {
    prepareGeometry(myCase, refinedGeometry.getLevel(iLevel));
  }

  for (std::size_t iLevel = 0; iLevel < refinedLattice.getLevelCount(); iLevel++) {
    prepareLattice(myCase, refinedGeometry.getLevel(iLevel), refinedLattice.getLevel(iLevel));
  }

  for (std::size_t iLevel = 0; iLevel < refinedLattice.getLevelCount(); iLevel++) {
    setInitialValues(myCase, refinedGeometry.getLevel(iLevel), refinedLattice.getLevel(iLevel));
  }

  refinedLattice.initialize(refinedGeometry);

  simulate(myCase);
}
