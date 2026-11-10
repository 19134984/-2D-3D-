/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2007, 2012 Jonas Latt, Mathias J. Krause,
 *  Vojtech Cvrcek, Peter Weisbrod,
 *  2025 Fedor Bukreev
 *  2026 Johannes Grafen
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

/* poiseuille2d.cpp:
 * This example examines a 2D Poiseuille flow
 * It illustrates the usage of different flow setups, boundary conditions
 * and computation of error norms.
 * The following simulation options can be combined freely:
 * - if the compiler flag ENABLE_MRT is set, mrt collision operators are used
 * - forced/ nonForced flow
 * - different boundary conditions
 * - simulation only or eoc convergence analysis
 */

// the main code of the simulation is in poiseuille2d.h as it is also used by the
// example poiseuille2dEoc

// set flag in order to use mrt collision operators instead of bgk
//#define ENABLE_MRT

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::names;

#ifdef ENABLE_MRT
  using MyCase = Case<
    NavierStokes, Lattice<FLOATING_POINT_TYPE, descriptors::D2Q9<tag::MRT, FORCE>>
  >;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;

  using BulkDynamics       = MRTdynamics<T,DESCRIPTOR>;
  using ForcedBulkDynamics = ForcedMRTdynamics<T,DESCRIPTOR>;
#else
  using MyCase = Case<
    NavierStokes, Lattice<FLOATING_POINT_TYPE, descriptors::D2Q9<FORCE>>
  >;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  using BulkDynamics       = BGKdynamics<T,DESCRIPTOR>;
  using ForcedBulkDynamics = ForcedBGKdynamics<T,DESCRIPTOR>;
#endif

enum class FlowType : int {
  forced = 0,
  nonForced = 1,
};

enum class BoundaryType : int {
  bounceBack   = 0,
  local        = 1,
  interpolated = 2,
  freeSlip     = 3,
  partialSlip  = 4,
};

namespace olb::parameters {

  struct TUNER_PARAM   : public descriptors::FIELD_BASE<1> { };
  struct COMPUTE_ERROR : public descriptors::TYPED_FIELD_BASE<bool,1> { };
  struct FLOW_TYPE     : public descriptors::TYPED_FIELD_BASE<FlowType, 1> { };
  struct BOUNDARY_TYPE : public descriptors::TYPED_FIELD_BASE<BoundaryType, 1> { };
  struct EOC           : public descriptors::TYPED_FIELD_BASE<bool,1> { };

};

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const T physLengthX   = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T physDeltaX    = parameters.get<parameters::DOMAIN_EXTENT>()[1]/parameters.get<parameters::RESOLUTION>();
  const T physLengthY   = parameters.get<parameters::DOMAIN_EXTENT>()[1];

  const Vector extent{physLengthX, physLengthY};
  const Vector origin{0, 0};
  IndicatorCuboid2D<T> cuboid(extent, origin);
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());

  if (parameters.get<parameters::FLOW_TYPE>() == FlowType::forced){
    mesh.getCuboidDecomposition().setPeriodicity({true, false});
  }
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

// Stores geometry information in form of material numbers
void prepareGeometry(MyCase& myCase)
{
  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();

  geometry.rename(0, 2);
  geometry.rename(2, 1, {1,1});

  const T physLengthX   = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T physLengthY   = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const T physDeltaX    = parameters.get<parameters::DOMAIN_EXTENT>()[1] / parameters.get<parameters::RESOLUTION>();

  if (parameters.get<parameters::FLOW_TYPE>() == FlowType::nonForced){
    Vector<T,2> extent(0.25*physDeltaX, physLengthY);
    Vector<T,2> origin;
    // Set material number for inflow
    IndicatorCuboid2D<T> inflow(extent, origin);
    geometry.rename(2, 3, 1, inflow);

    // Set material number for outflow
    origin[0] = physLengthX - physDeltaX / 4;
    IndicatorCuboid2D<T> outflow(extent, origin);
    geometry.rename(2, 4, 1, outflow);
  }
  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  // Removes all not needed boundary voxels inside the surface
  geometry.innerClean();
  geometry.checkForErrors();
  geometry.print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(MyCase& myCase)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(NavierStokes{});

  const T physLengthY   = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const int resolution  = parameters.get<parameters::RESOLUTION>();
  const T latticeRelaxationTime = parameters.get<parameters::LATTICE_RELAXATION_TIME>();
  const T maxVelocity   = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T Re            = parameters.get<parameters::REYNOLDS>();
  const T physDensity   = parameters.get<parameters::PHYS_CHAR_DENSITY>();
  const T tuner         = parameters.get<parameters::TUNER_PARAM>();

  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR>>(
    resolution,
    latticeRelaxationTime,
    physLengthY,
    maxVelocity,
    maxVelocity*physLengthY/Re,
    physDensity
  );

  const UnitConverter<T,DESCRIPTOR>& converter = lattice.getUnitConverter();
  converter.print();

  if (parameters.get<parameters::FLOW_TYPE>() == FlowType::forced){
    dynamics::set<ForcedBulkDynamics>(lattice, geometry, 1);
  } else {
    dynamics::set<BulkDynamics>(lattice, geometry, 1);
  }

  switch(parameters.get<parameters::BOUNDARY_TYPE>()){
    case BoundaryType::bounceBack:
      boundary::set<boundary::BounceBack>(lattice, geometry, 2);
      break;
    case BoundaryType::freeSlip:
      boundary::set<boundary::FullSlip>(lattice, geometry, 2);
      break;
    case BoundaryType::partialSlip:
      boundary::set<boundary::PartialSlip>(lattice, geometry, 2);
      lattice.template setParameter<descriptors::TUNER>(tuner);
      break;
    case BoundaryType::local:
      boundary::set<boundary::LocalVelocity>(lattice, geometry, 2);
      break;
    default:
      boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 2);
  }

  if (parameters.get<parameters::FLOW_TYPE>() == FlowType::nonForced){
    if (parameters.get<parameters::BOUNDARY_TYPE>() == BoundaryType::local){
      boundary::set<boundary::LocalVelocity>(lattice, geometry, 3);
      boundary::set<boundary::LocalPressure>(lattice, geometry, 4);
    } else {
      boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 3);
      boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 4);
    }
  }

  lattice.template setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase){
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(NavierStokes{});

  const UnitConverter<T,DESCRIPTOR>& converter = lattice.getUnitConverter();

  const T physLengthX   = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T physLengthY   = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const T physDeltaX    = parameters.get<parameters::DOMAIN_EXTENT>()[1] / parameters.get<parameters::RESOLUTION>();

  // Initial conditions
  T Lx = physLengthX;
  T Ly = physLengthY;

  if (parameters.get<parameters::BOUNDARY_TYPE>() == BoundaryType::bounceBack){
    Lx -= physDeltaX;
    Ly -= physDeltaX;
  }
  T p0 = 8. * converter.getPhysViscosity() * converter.getCharPhysVelocity() * Lx / (Ly*Ly);

  std::vector<T> poiseuilleForce(2, T());
  poiseuilleForce[0]= p0 / Lx / converter.getConversionFactorForce() * converter.getConversionFactorMass();
  AnalyticalConst2D<T,T> force(poiseuilleForce);

  if (parameters.get<parameters::FLOW_TYPE>() == FlowType::forced){
    // Initialize force
    fields::set<FORCE>(lattice, geometry.getMaterialIndicator({1,2}), force);
  } else {
    AnalyticalLinear2D<T,T> pressure(-p0/Lx, 0, p0);

    const T maxVelocity = converter.getCharPhysVelocity();

    bool bounceBack = false;
    if (parameters.get<parameters::BOUNDARY_TYPE>() == BoundaryType::bounceBack){
      bounceBack = true;
    }

    const T radius = (bounceBack) ? (T(0.5) * (physLengthY - physDeltaX)) : (T(0.5) * physLengthY);
    std::vector<T> axisPoint = {physLengthX/T{2}, physLengthY/T{2}};
    std::vector<T> axisDirection = {1, 0};
    Poiseuille2D<T> u(axisPoint, axisDirection, maxVelocity, radius);

    // Initialize all values of distribution functions to their local equilibrium
    momenta::setPressure(lattice, geometry.getMaterialIndicator({1,2,3,4}), pressure);
    momenta::setVelocity(lattice, geometry.getMaterialIndicator({1,2,3,4}), u);
  }

  // Make the lattice ready for simulation
  lattice.initialize();
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{ }

// Compute error norms
std::vector<MyCase::value_t> errorNorms(MyCase& myCase)
{
  OstreamManager clout( std::cout,"error" );
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry   = myCase.getGeometry();
  auto& lattice    = myCase.getLattice(NavierStokes{});

  const UnitConverter<T,DESCRIPTOR>& converter = lattice.getUnitConverter();
  lattice.setProcessingContext(ProcessingContext::Evaluation);

  const T maxVelocity = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T physLengthX = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T physLengthY = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const T physDeltaX  = parameters.get<parameters::DOMAIN_EXTENT>()[1] / parameters.get<parameters::RESOLUTION>();

  int tmp[]= { };
  T result[2]= { };
  std::vector<T> absErrors;

  bool bounceBack = false;
  if (parameters.get<parameters::BOUNDARY_TYPE>() == BoundaryType::bounceBack){
    bounceBack = true;
  }

  // velocity error
  const T radius = (bounceBack) ?
    T(0.5) * (physLengthY - physDeltaX) : T(0.5) * physLengthY;
  std::vector<T> axisPoint(2, T());
  axisPoint[0] = physLengthX/2.;
  axisPoint[1] = physLengthY/2.;
  std::vector<T> axisDirection(2, T());
  axisDirection[0] = 1;
  axisDirection[1] = 0;
  Poiseuille2D<T> uSol(axisPoint, axisDirection, maxVelocity, radius);
  SuperLatticePhysVelocity2D<T,DESCRIPTOR> u(lattice, converter);
  auto indicatorF = geometry.getMaterialIndicator(1);

  SuperAbsoluteErrorL1Norm2D<T> absVelocityErrorNormL1(u, uSol, indicatorF);
  absVelocityErrorNormL1(result, tmp);
  absErrors.push_back(result[0]);
  clout << "velocity-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm2D<T> relVelocityErrorNormL1(u, uSol, indicatorF);
  relVelocityErrorNormL1(result, tmp);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm2D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, tmp);
  absErrors.push_back(result[0]);
  clout << "velocity-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm2D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, tmp);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm2D<T> absVelocityErrorNormLinf(u, uSol, indicatorF);
  absVelocityErrorNormLinf(result, tmp);
  absErrors.push_back(result[0]);
  clout << "velocity-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm2D<T> relVelocityErrorNormLinf(u, uSol, indicatorF);
  relVelocityErrorNormLinf(result, tmp);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;

  // strainRate error
  PoiseuilleStrainRate2D<T,T,DESCRIPTOR> sSol(converter, T(2) * radius);
  SuperLatticePhysStrainRate2D<T,DESCRIPTOR> s(lattice, converter);

  SuperAbsoluteErrorL1Norm2D<T> absStrainRateErrorNormL1(s, sSol, indicatorF);
  absStrainRateErrorNormL1(result, tmp);
  absErrors.push_back(result[0]);
  clout << "strainRate-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm2D<T> relStrainRateErrorNormL1(s, sSol, indicatorF);
  relStrainRateErrorNormL1(result, tmp);
  clout << "; strainRate-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm2D<T> absStrainRateErrorNormL2(s, sSol, indicatorF);
  absStrainRateErrorNormL2(result, tmp);
  absErrors.push_back(result[0]);
  clout << "strainRate-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm2D<T> relStrainRateErrorNormL2(s, sSol, indicatorF);
  relStrainRateErrorNormL2(result, tmp);
  clout << "; strainRate-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm2D<T> absStrainRateErrorNormLinf(s, sSol, indicatorF);
  absStrainRateErrorNormLinf(result, tmp);
  absErrors.push_back(result[0]);
  clout << "strainRate-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm2D<T> relStrainRateErrorNormLinf(s, sSol, indicatorF);
  relStrainRateErrorNormLinf(result, tmp);
  clout << "; strainRate-Linf-error(rel)=" << result[0] << std::endl;


  if (parameters.get<parameters::FLOW_TYPE>() == FlowType::nonForced){
    // pressure error
    int Lx = converter.getLatticeLength(physLengthX);
    int Ly = converter.getLatticeLength(physLengthY);
    T p0 = 8.*converter.getLatticeViscosity()*converter.getCharLatticeVelocity()*Lx/T( Ly*Ly );
    AnalyticalLinear2D<T,T> pressureSol(-converter.getPhysPressure(p0)/physLengthX, 0, converter.getPhysPressure(p0));
    SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure(lattice, converter);

    SuperAbsoluteErrorL1Norm2D<T> absPressureErrorNormL1(pressure, pressureSol, indicatorF);
    absPressureErrorNormL1(result, tmp);
    absErrors.push_back(result[0]);
    clout << "pressure-L1-error(abs)=" << result[0];
    SuperRelativeErrorL1Norm2D<T> relPressureErrorNormL1(pressure, pressureSol, indicatorF);
    relPressureErrorNormL1(result, tmp);
    clout << "; pressure-L1-error(rel)=" << result[0] << std::endl;

    SuperAbsoluteErrorL2Norm2D<T> absPressureErrorNormL2(pressure, pressureSol, indicatorF);
    absPressureErrorNormL2(result, tmp);
    absErrors.push_back(result[0]);
    clout << "pressure-L2-error(abs)=" << result[0];
    SuperRelativeErrorL2Norm2D<T> relPressureErrorNormL2(pressure, pressureSol, indicatorF);
    relPressureErrorNormL2(result, tmp);
    clout << "; pressure-L2-error(rel)=" << result[0] << std::endl;

    SuperAbsoluteErrorLinfNorm2D<T> absPressureErrorNormLinf(pressure, pressureSol, indicatorF);
    absPressureErrorNormLinf(result, tmp);
    absErrors.push_back(result[0]);
    clout << "pressure-Linf-error(abs)=" << result[0];
    SuperRelativeErrorLinfNorm2D<T> relPressureErrorNormLinf(pressure, pressureSol, indicatorF);
    relPressureErrorNormLinf(result, tmp);
    clout << "; pressure-Linf-error(rel)=" << result[0] << std::endl;
  }
  return absErrors;
}

// Output to console and files
void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  OstreamManager clout(std::cout, "getResults");
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  auto& parameters = myCase.getParameters();
  auto& lattice    = myCase.getLattice(NavierStokes{});
  auto& geometry   = myCase.getGeometry();

  const UnitConverter<T,DESCRIPTOR>& converter = lattice.getUnitConverter();
  const T maxPhysT        = parameters.get<parameters::PHYS_MAX_T>();
  const T physLengthX     = parameters.get<parameters::DOMAIN_EXTENT>()[0];
  const T physLengthY     = parameters.get<parameters::DOMAIN_EXTENT>()[1];
  const T physDeltaX      = parameters.get<parameters::DOMAIN_EXTENT>()[1] / parameters.get<parameters::RESOLUTION>();
  const bool hasConverged = parameters.get<parameters::CONVERGED>();
  const bool computeError = parameters.get<parameters::COMPUTE_ERROR>();

  // absolute errors of the last evaluation, plotted over the resolution
  static std::vector<T> absErrors;

  const bool lastTimeStep = (hasConverged || (iT + 1 == converter.getLatticeTime(maxPhysT)));
  bool noslipBoundary = true;
  if (parameters.get<parameters::BOUNDARY_TYPE>() == BoundaryType::freeSlip ||
      parameters.get<parameters::BOUNDARY_TYPE>() == BoundaryType::partialSlip){
    noslipBoundary = false;
  }
  const int statIter = converter.getLatticeTime(parameters.get<parameters::PHYS_STAT_ITER_T>());

  // VTK and image output only if no EOC analysis
  if (parameters.get<parameters::VTK_ENABLED>()) {
    SuperVTMwriter2D<T> vtmWriter("poiseuille2d");
    const int vtmIter  = converter.getLatticeTime(parameters.get<parameters::PHYS_VTK_ITER_T>());

    SuperLatticePhysVelocity2D<T,DESCRIPTOR> velocity(lattice, converter);
    //geometryF<T,DESCRIPTOR> materials( geometry );
    vtmWriter.addFunctor(velocity);
    //vtmWriter.addFunctor( materials);

    SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure(lattice, converter);

    if (parameters.get<parameters::FLOW_TYPE>() == FlowType::nonForced){
      vtmWriter.addFunctor(pressure);
    }

    bool bounceBack = false;
    if (parameters.get<parameters::BOUNDARY_TYPE>() == BoundaryType::bounceBack){
      bounceBack = true;
    }

    const T maxVelocity = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
    const T radius = bounceBack ?
      T(0.5) * (physLengthY - physDeltaX) : T(0.5) * physLengthY;
    std::vector<T> axisPoint(2, T());
    axisPoint[0] = physLengthX/2.;
    axisPoint[1] = physLengthY/2.;
    std::vector<T> axisDirection(2, T());
    axisDirection[0] = 1;
    axisDirection[1] = 0;
    Poiseuille2D<T> analyticalVelocity(axisPoint, axisDirection, maxVelocity, radius);
    SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR> analyticalVelocityLattice(
      analyticalVelocity, lattice);
    analyticalVelocityLattice.getName() = "analytical solution";
    if (noslipBoundary) {
      vtmWriter.addFunctor(analyticalVelocityLattice);
    }

    if ( iT==0 ) {
      // Writes the geometry, cuboid no. and rank no. as vti file for visualization
      SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid(lattice);
      SuperLatticeRank2D<T, DESCRIPTOR> rank(lattice);
      SuperLatticeDiscreteNormal2D<T, DESCRIPTOR> discreteNormal(lattice, geometry, geometry.getMaterialIndicator({2, 3}));
      SuperLatticeDiscreteNormalType2D<T, DESCRIPTOR> discreteNormalType(
        lattice, geometry, geometry.getMaterialIndicator({2, 3, 4, 5}));
      vtmWriter.write(cuboid);
      vtmWriter.write(rank);
      vtmWriter.write(discreteNormal);
      vtmWriter.write(discreteNormalType);
      vtmWriter.createMasterFile();
    }

    // Writes the vtm files and profile text file
    if (iT%vtmIter==0 || lastTimeStep) {
      lattice.setProcessingContext(ProcessingContext::Evaluation);

      vtmWriter.write(iT);

      SuperEuklidNorm2D<T, DESCRIPTOR> normVel(velocity);
      BlockReduction2D2D<T> planeReduction(normVel, 600, BlockDataSyncMode::ReduceOnly);
      // write output as JPEG
      heatmap::write(planeReduction, iT);
    }
  }

  // Output on the console
  if (iT%statIter==0 || lastTimeStep) {
    // Timer console output
    timer.update(iT);
    timer.printStep();

    // Lattice statistics console output
    lattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Error norms
    if (noslipBoundary) {
      if ((computeError || lastTimeStep)){
        absErrors = errorNorms(myCase);
      }
    }
  }

  static Gnuplot<T> gplot("centerVelocity");
  static Gnuplot<T> gplotError("error", false, "set terminal png size 900,850",
                               Gnuplot<T>::LOGLOG, Gnuplot<T>::LINREG);
  // Gnuplot output
  if ((noslipBoundary) && (lastTimeStep) && parameters.get<parameters::GNUPLOT_ENABLED>()
      && parameters.get<parameters::EOC>()) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    std::vector<std::string> errorNames {
      "velocity L1 abs error", "velocity L2 abs error", "velocity Linf abs error",
      "strain rate L1 abs error", "strain rate L2 abs error", "strain rate Linf abs error",
      "pressure L1 abs error", "pressure L2 abs error", "pressure Linf abs error" };
    errorNames.resize(absErrors.size());
    gplotError.setLabel("resolution N", "absolute error");
    gplotError.setData(T(converter.getResolution()), absErrors, errorNames, "below offset 0,-1",
                       std::vector<char>(absErrors.size(), 'p'));
    gplotError.writePNG(-1, -1, "error");
  }

  if ((noslipBoundary) && (lastTimeStep) && parameters.get<parameters::GNUPLOT_ENABLED>()
      && !parameters.get<parameters::EOC>()) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    // plot velocity magnitude over line through the center of the simulation domain
    const T maxVelocity = converter.getPhysVelocity( converter.getCharLatticeVelocity() );
    T Ly = physLengthY / physDeltaX;
    bool bounceBack = false;
    if (parameters.get<parameters::BOUNDARY_TYPE>() == BoundaryType::bounceBack){
      bounceBack = true;
    }
    const T radius = (bounceBack) ?
      T(0.5) * (Ly - physDeltaX) : T(0.5) * physLengthY;
    std::vector<T> axisPoint{physLengthX/T(2), physLengthY/T(2)};
    std::vector<T> axisDirection{1, 0};
    Poiseuille2D<T> uSol(axisPoint, axisDirection, maxVelocity, radius);
    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(lattice, converter);
    AnalyticalFfromSuperF2D<T> intpolateVelocity(velocity, true);
    T point[2] { };
    point[0] = physLengthX/2.;
    T analytical[2] { };
    T numerical[2] { };
    for (int iY=0; iY<=Ly; ++iY) {
      point[1] = (T)iY/Ly;
      uSol(analytical, point);
      intpolateVelocity(numerical, point);
      gplot.setData(iY*physDeltaX, {analytical[0], numerical[0]}, {"analytical","numerical"});
    }
    gplot.writePNG(-1, -1, "centerVelocity");
  }
}

void simulate(MyCase& myCase){
  OstreamManager clout(std::cout, "simulate");
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();

  const T physMaxT = parameters.get<parameters::PHYS_MAX_T>();
  const std::size_t iTmax = myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(physMaxT);

  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  util::ValueTracer<T> converge(myCase.getLattice(NavierStokes{}).getUnitConverter().getLatticeTime(parameters.get<parameters::CONV_ITER>()),
                                parameters.get<parameters::CONVERGENCE_PRECISION>());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {

    if (converge.hasConverged()){
      parameters.set<parameters::CONVERGED>(true);
      getResults(myCase, timer, iT);
      clout << "Converged after " << iT << " iterations." << std::endl;
      break;
    }

    /// === Step 8.1: Update the Boundary Values and Fields at Times ===
    setTemporalValues(myCase, iT);

    /// === Step 8.2: Collide and Stream Execution ===
    myCase.getLattice(NavierStokes{}).collideAndStream();

    /// === Step 8.3: Computation and Output of the Results ===
    getResults(myCase, timer, iT);
    converge.takeValue(myCase.getLattice(NavierStokes{}).getStatistics().getMaxU(), false);
  }

  timer.stop();
  timer.printSummary();
}

void run(MyCase::ParametersD& myCaseParameters)
{
  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase);

  /// === Step 7: Set Initial Conditions ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);
}
