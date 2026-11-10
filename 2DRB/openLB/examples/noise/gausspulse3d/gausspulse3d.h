/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Philipp Spelten
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
#pragma once

#include <olb.h>
#include "../noiseauxiliary.h"

using namespace olb;

using MyCase = Case<
  names::NavierStokes, Lattice<double, descriptors::D3Q19<>>
>;

enum FarFieldType : int {
  ETERNAL       = 0,
  PERIODIC      = 1,
  LOCAL         = 2,
  INTERPOLATED  = 3,
  SPONGE        = 4,
  CBC_SPONGE    = 5,
  CBC           = 6,
  COARSE        = 7
};
// enum CollisionType : int {
//   BGK = 0,
//   HRR = 1,
//   HHRR = 2,
//   TRT = 3 ,
//   MRT = 4,
//   SMAGORINSKY_BGK = 5
// };
size_t defaultRes   = 51;
double defaultSigma = 0.01;
double defaultMa    = 0.1;
double defaultAmp   = 1e-3;
double defaultL     = 1.;
double defaultAlpha = log(2.) / (1. / 20. * 1. / 20.);


namespace olb::parameters {
  // struct COLLISION_TYPE       : public descriptors::TYPED_FIELD_BASE<CollisionType,1> { };
  // DOMAIN
  struct CORE_EXTENT          : public descriptors::FIELD_BASE<0,1> { };  // fluid domain size [m]
  // PULSE
  struct PULSE_AMPLITUDE      : public descriptors::FIELD_BASE<1> { };    // pressure amplitude of Gaussian pulse
  struct PULSE_PHYS_VELOCITY  : public descriptors::FIELD_BASE<0,1> { };    // background velocity of Gaussian pulse
  struct PULSE_ALPHA          : public descriptors::FIELD_BASE<1> { };    // distribution factor of Gaussian pulse
  // BOUNDARY CONDITIONS
  struct FAR_FIELD_TYPE       : public descriptors::TYPED_FIELD_BASE<FarFieldType,1> { };  // counter, see beginning of file
  struct ETERNALSCALE         : public descriptors::FIELD_BASE<1> { };               // amout to extend domain by for 'eternal' case
  struct SPONGE_DEPTH_LU      : public descriptors::TYPED_FIELD_BASE<size_t,1> { };  // number of points for sponge layer
  struct SPONGE_STRENGTH      : public descriptors::FIELD_BASE<1> { };               // maximum sponge strength of sponge layer
  struct CBC_SIGMA            : public descriptors::FIELD_BASE<1> { };               // factor for K1 in CBC
  // TIMING AND OUTPUTS
  struct INTERVAL_TABLE             : public descriptors::TYPED_FIELD_BASE<size_t,1> { };  // tabular outputs every INTERVAL_TABLE iterations (for line plot)
  struct DO_GRAPHICAL_OUTPUT  : public descriptors::TYPED_FIELD_BASE<bool,1> { };  // set to do images
  struct DO_L2_PLOT           : public descriptors::TYPED_FIELD_BASE<bool,1> { };  // set to do l2 plot
  struct WAIT_DEBUG           : public descriptors::TYPED_FIELD_BASE<bool,1> { };  // set to wait for gbd
  struct OUTDIR_SUFFIX        : public descriptors::TYPED_FIELD_BASE<std::string,1> { };  // set a custom suffix for the output directory
}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  OstreamManager clout(std::cout, "createMesh");
  clout << "createMesh ..." << std::endl;
  using T = MyCase::value_t;

  Vector        extentFluid     = parameters.get<parameters::CORE_EXTENT>();
  Vector        originFluid     = -0.5*extentFluid;
  FarFieldType  farFieldType    = parameters.get<parameters::FAR_FIELD_TYPE>();
  const T       physDeltaX      = parameters.get<parameters::PHYS_DELTA_X>();
  const T       spongeDepthPU   = physDeltaX * parameters.get<parameters::SPONGE_DEPTH_LU>();
  const T       eternalscale    = parameters.get<parameters::ETERNALSCALE>();

  Vector extentDomain = extentFluid;
  switch ( farFieldType ) {
    case ETERNAL:  // extend the domain
      extentDomain *= eternalscale;
      break;
    case PERIODIC:  // reference domain size
      break;
    case LOCAL:
    case CBC:
    case INTERPOLATED:  // add one layer in each direction
      extentDomain += 2 * physDeltaX;
      break;
    case SPONGE:
    case COARSE:  // add boundary layer
      extentDomain += 2 * spongeDepthPU;
      break;
    case CBC_SPONGE:  // add boundary layer and one layer for outlet
      extentDomain[0] += spongeDepthPU + physDeltaX;
      extentDomain[1] += 2 * spongeDepthPU;
      extentDomain[2] += 2 * spongeDepthPU;
      break;
  }
  Vector originDomain = -0.5*extentDomain;
  if ( farFieldType == CBC_SPONGE ) originDomain[0] = originFluid[0] - spongeDepthPU;

  clout << "Fluid Domain = " << extentFluid[0] << "^3; Simulation Domain = " << extentDomain[0] << "^3" << std::endl;
  parameters.set<parameters::DOMAIN_EXTENT>( {extentDomain[0],extentDomain[1],extentDomain[2]} );
  parameters.set<parameters::DOMAIN_ORIGIN>( {originDomain[0],originDomain[1],originDomain[2]} );

  IndicatorCuboid3D<T> cuboid( extentDomain, originDomain );
  size_t nC = ( farFieldType == COARSE ) ? 1 : singleton::mpi().getSize();
  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, nC);

  switch ( farFieldType ) {
    case LOCAL:
    case INTERPOLATED:
    case CBC:
      mesh.getCuboidDecomposition().setPeriodicity( {false, false, false} );
      break;
    case CBC_SPONGE:
      mesh.getCuboidDecomposition().setPeriodicity( {false, true, true} );
      break;
    case ETERNAL:
    case PERIODIC:
    case SPONGE:
    case COARSE:
      mesh.getCuboidDecomposition().setPeriodicity( {true, true, true} );
      break;
  }
  mesh.setOverlap( parameters.get<parameters::OVERLAP>() );

  clout << "createMesh ... OK" << std::endl;
  return mesh;
}

void prepareGeometry( MyCase& myCase ) {
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto&         geometry      = myCase.getGeometry();
  auto&         parameters    = myCase.getParameters();
  FarFieldType  farFieldType  = parameters.get<parameters::FAR_FIELD_TYPE>();
  const T       physDeltaX2   = parameters.get<parameters::PHYS_DELTA_X>() / 2.;
  const Vector  extentFluid   = parameters.get<parameters::CORE_EXTENT>();
  const Vector  originFluid   = -0.5*extentFluid;
  const T       spongeDepthPU = parameters.get<parameters::SPONGE_DEPTH_LU>() * physDeltaX2*2.;
  IndicatorCuboid3D<T> domainFluid( extentFluid, originFluid );

  // all nodes to temporary type
  geometry.rename(0, 2);
  Vector origin = geometry.getStatistics().getMinPhysR(2);
  Vector extent = geometry.getStatistics().getMaxPhysR(2) - geometry.getStatistics().getMinPhysR(2);

  switch ( farFieldType ) {
    case CBC:
    case LOCAL:
    case INTERPOLATED: {
      // set fluid material, keeping the outside layer on 2
      geometry.rename( 2, 1, {1, 1, 1} );

      // Set material number for outflow
      origin        -= physDeltaX2;
      extent        += 2 * physDeltaX2;
      origin[0] = geometry.getStatistics().getMaxPhysR(2)[0] - physDeltaX2;
      extent[0] = 2 * physDeltaX2;
      IndicatorCuboid3D<T> outflow(extent, origin);
      geometry.rename( 2, 6, 1, outflow );

      // Set material number for inflow (rest)
      geometry.rename( 2, 4 );
    } break;
    case SPONGE: {
      extent -= 6 * physDeltaX2;
      origin += 3 * physDeltaX2;
      IndicatorCuboid3D<T> spongeOutside(extent, origin);
      geometry.rename( 2, 3, spongeOutside );
      geometry.rename( 3, 1, domainFluid );
      geometry.rename( 2, 1 );
    } break;
    case CBC_SPONGE: {
      // rename all to fluid, except inflow/outflow
      geometry.rename( 2, 1, {1,0,0} );

      // inflow damping
      origin = geometry.getStatistics().getMinPhysR( 1 ) - physDeltaX2;
      extent = geometry.getStatistics().getMaxPhysR( 1 ) - geometry.getStatistics().getMinPhysR( 1 ) + 2*physDeltaX2;
      extent[0] = spongeDepthPU + 2*physDeltaX2;
      extent[0] += 2*physDeltaX2;  // one extra layer in x-direction to be renamed to outflow
      IndicatorCuboid3D<T> dampingIn( extent, origin );
      geometry.rename( 1, 3, dampingIn );

      // Set material number for inflow
      origin = geometry.getStatistics().getMinPhysR( 2 ) - physDeltaX2;
      extent = geometry.getStatistics().getMaxPhysR( 2 ) - geometry.getStatistics().getMinPhysR( 2 ) + 2*physDeltaX2;
      extent[0] = 2*physDeltaX2;
      IndicatorCuboid<T,3> inflow( extent, origin );
      geometry.rename( 2, 4, 3, inflow );

      // rename remaining to damping
      geometry.rename( 1, 3 );

      // rename core domain to fluid
      geometry.rename( 3, 1, domainFluid );

      // rename last two node layers before outflow to fluid (for normal calculation of outflow)
      origin = geometry.getStatistics().getMinPhysR( 3 ) - physDeltaX2;
      extent = geometry.getStatistics().getMaxPhysR( 3 ) - geometry.getStatistics().getMinPhysR( 3 ) + 2*physDeltaX2;
      origin[0] = geometry.getStatistics().getMaxPhysR( 3 )[0] - 3*physDeltaX2;
      extent[0] = 4*physDeltaX2;
      IndicatorCuboid3D<T> domainBeforeOutflow(extent, origin);
      geometry.rename( 3, 1, domainBeforeOutflow );

      // Set material number for outflow
      origin[0] = geometry.getStatistics().getMaxPhysR(2)[0] - 2*physDeltaX2;
      IndicatorCuboid3D<T> outflow(extent, origin);
      geometry.rename( 2, 6, 1, outflow );
    } break;
    case ETERNAL:
    case PERIODIC:
    case COARSE: {
      geometry.rename( 2, 1 );
    } break;
  }

  geometry.checkForErrors();
  geometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( MyCase& myCase ) {
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "prepareLattice ..." << std::endl;

  using T           = MyCase::value_t;
  using DESCRIPTOR  = MyCase::descriptor_t;
  auto& geometry    = myCase.getGeometry();
  auto& parameters  = myCase.getParameters();
  auto& lattice     = myCase.getLattice(names::NavierStokes{});

  const FarFieldType farFieldType = parameters.get<parameters::FAR_FIELD_TYPE>();
  const T physExtentX       = parameters.get<parameters::CORE_EXTENT>()[0];
  const T physCharLength    = physExtentX;
  const T physDeltaX        = parameters.get<parameters::PHYS_DELTA_X>();
  const T physDeltaT        = parameters.get<parameters::PHYS_DELTA_T>();
  const T physCharVelocity  = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T physCharViscosity = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T physCharDensity   = parameters.get<parameters::PHYS_CHAR_DENSITY>();
  const Vector uInfty       = parameters.get<parameters::PULSE_PHYS_VELOCITY>();
  Vector  domainLengths     = parameters.get<parameters::DOMAIN_EXTENT>();
  const T spongeStrength    = parameters.get<parameters::SPONGE_STRENGTH>();
  const T spongeDepthPU     = parameters.get<parameters::SPONGE_DEPTH_LU>() * physDeltaX;

  // Set up a unit converter with the characteristic physical units
  lattice.setUnitConverter<UnitConverter<T,DESCRIPTOR>>(
    physDeltaX,         // physDeltaX: spacing between two lattice cells in [m]
    physDeltaT,         // physDeltaT: time step in [s]
    physCharLength,     // charPhysLength: reference length of simulation geometry in [m]
    physCharVelocity,   // physCharVelocity: highest expected velocity during simulation in [m/s]
    physCharViscosity,  // physCharViscosity: physical kinematic viscosity in [m^2/s]
    physCharDensity     // physCharDensity: physical density [kg/m^3]
  );
  auto& converter = lattice.getUnitConverter();
  converter.print();

  using BULK = BGKdynamics<T,DESCRIPTOR>;
  // auto bulkIndicator = geometry.getMaterialIndicator( {1, 4, 6} );
  // lattice.defineDynamics<BULK>( bulkIndicator );
  lattice.defineDynamics<BULK>( geometry, 1 );
  switch ( farFieldType ) {
    case LOCAL:
      boundary::set<boundary::LocalVelocity>(lattice, geometry, 4);  // inflow
      boundary::set<boundary::LocalPressure>(lattice, geometry, 6);  // outflow
      break;
    case INTERPOLATED:
      boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 4);  // inflow
      boundary::set<boundary::InterpolatedPressure>(lattice, geometry, 6);  // outflow
      break;
    case CBC_SPONGE:
      domainLengths[0] += spongeDepthPU;  // "extend" the domain, such that DampingTerm sets sigma from the center of CORE_EXTENT
      // inflow
      boundary::set<T,DESCRIPTOR,boundary::LocalVelocity<T,DESCRIPTOR>>(
        lattice,
        geometry.getMaterialIndicator(4),  // boundaryI
        geometry.getMaterialIndicator(3),  // fluidI (here: sponge)
        geometry.getMaterialIndicator(0)   // outsideI
      );
      // outflow
      boundary::setCBC<T,DESCRIPTOR,boundary::Characteristic<T,DESCRIPTOR>>(
        lattice,
        geometry.getMaterialIndicator(6),  // boundaryI
        geometry.getMaterialIndicator(1),  // fluidI
        geometry.getMaterialIndicator(0)   // outsideI
      );
    case SPONGE:
      // sponge
      lattice.defineDynamics<
        SpongeLayerDynamics<T,DESCRIPTOR>
          >( geometry, 3 );
      break;
    case CBC:
      boundary::setCBC<T,DESCRIPTOR,boundary::Characteristic<T,DESCRIPTOR>>(
        lattice,
        geometry.getMaterialIndicator({4,6}),  // boundaryI
        geometry.getMaterialIndicator(1),  // fluidI
        geometry.getMaterialIndicator(0)   // outsideI
      );
      break;
    case ETERNAL:
    case PERIODIC:
    case COARSE:
      break;
  }

  if ( farFieldType == CBC || farFieldType == CBC_SPONGE ) {
    // Setting CBC fields
    lattice.setParameter<descriptors::cbc::CBC_SIGMA>(  parameters.get<parameters::CBC_SIGMA>() );
    lattice.setParameter<descriptors::cbc::CBC_L>(      converter.getCharPhysLength() );
    lattice.setParameter<descriptors::cbc::CBC_MA>(     converter.getCharLatticeVelocity() );
    lattice.setParameter<descriptors::cbc::RHO_INFTY>( 1. );
    lattice.setParameter<descriptors::cbc::U_INFTY>( {  converter.getLatticeVelocity(uInfty[0]),
                                                        converter.getLatticeVelocity(uInfty[1]),
                                                        converter.getLatticeVelocity(uInfty[2]) } );
    lattice.setParameter<descriptors::cbc::FLOW_DIRECTION>( 0 );
    lattice.setParameter<descriptors::cbc::FLOW_ORIENTATION>( 1 );
  }

  if ( farFieldType == SPONGE || farFieldType == CBC_SPONGE ) {
    // Setting sponge fields
    auto spongeIndicator = geometry.getMaterialIndicator({3});
    AnalyticalConst3D<T,T>  uFar( converter.getLatticeVelocity(uInfty[0]),
                                  converter.getLatticeVelocity(uInfty[1]),
                                  converter.getLatticeVelocity(uInfty[2]) );
    AnalyticalConst3D<T,T>  uxFar( converter.getLatticeVelocity(uInfty[0]) );
    AnalyticalConst3D<T,T>  uyFar( converter.getLatticeVelocity(uInfty[1]) );
    AnalyticalConst3D<T,T>  uzFar( converter.getLatticeVelocity(uInfty[2]) );
    AnalyticalConst3D<T,T>  rhoFar(1.);
    lattice.defineField<descriptors::UX>( spongeIndicator, uxFar );
    lattice.defineField<descriptors::UY>( spongeIndicator, uyFar );
    lattice.defineField<descriptors::UZ>( spongeIndicator, uzFar );
    lattice.defineField<descriptors::DENSITY>( spongeIndicator, rhoFar );
    // define sponge layer scaling
    Vector xMin = geometry.getStatistics().getMinPhysR( 3 );
    Vector xMax = geometry.getStatistics().getMaxPhysR( 3 );
    AnalyticalConst3D<T,T> spongeStrengthF( spongeStrength );
    if ( farFieldType == SPONGE ) {  // not outlet sponge for CBC_SPONGE
      T xBoundaries[6] = {xMin[0], xMax[0], xMin[1], xMax[1], xMin[2], xMax[2]};
      size_t directions[6] = {0,0,1,1,2,2};
      DampingTerm<3,T,6> sigma( directions, xBoundaries, spongeDepthPU );
      lattice.defineField<descriptors::DAMPING>(
        spongeIndicator, sigma*spongeStrengthF );
    } else if ( farFieldType == CBC_SPONGE ) {
      T xBoundaries[5] = {xMin[0], xMin[1], xMax[1], xMin[2], xMax[2]};
      size_t directions[5] = {0,1,1,2,2};
      DampingTerm<3,T,5> sigma( directions, xBoundaries, spongeDepthPU );
      lattice.defineField<descriptors::DAMPING>(
        spongeIndicator, sigma*spongeStrengthF );
    }
    // === alternatively, define sigma as a global constant
    // AnalyticalConst<3,T,T> sigma( spongeStrength );
    // lattice.defineField<descriptors::DAMPING>(
    //   spongeIndicator, sigma;
  }

  // CAREFUL: OMEGA must be set AFTER defining dynamics and boundary conditions
  lattice.setParameter<descriptors::OMEGA>( converter.getLatticeRelaxationFrequency() );
  clout << "prepareLattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  OstreamManager clout(std::cout, "setInitialValues");
  clout << std::endl << "setInitialValues ..." << std::endl;
  using T = MyCase::value_t;
  auto& lattice     = myCase.getLattice(names::NavierStokes{});
  auto& geometry    = myCase.getGeometry();
  auto& parameters  = myCase.getParameters();
  auto& converter   = lattice.getUnitConverter();

  const T amplitude   = parameters.get<parameters::PULSE_AMPLITUDE>();
  const T alpha       = parameters.get<parameters::PULSE_ALPHA>();
  const Vector uInfty = parameters.get<parameters::PULSE_PHYS_VELOCITY>();

  // Initial velocity and density
  AnalyticalConst3D<T,T> u( converter.getLatticeVelocity(uInfty[0]),
                            converter.getLatticeVelocity(uInfty[1]),
                            converter.getLatticeVelocity(uInfty[2]) );
  AcousticPulse<3,T>      densityProfile(1., amplitude, alpha);
  size_t res  = converter.getResolution();
  T physDx    = converter.getPhysDeltaX();
  linePlot<3,T>( densityProfile, res, physDx, "pulse_diag", "density [LU]", diagonal2d);
  linePlot<3,T>( densityProfile, res, physDx, "pulse_hline", "density [LU]", horizontal);

  /// Initialize populations to equilibrium state
  momenta::setVelocity(lattice, geometry.getMaterialIndicator({0,1,2,3,4,5,6}), u);
  momenta::setDensity(lattice, geometry.getMaterialIndicator({0,1,2,3,4,5,6}), densityProfile);
  auto bulkIndicator = geometry.getMaterialIndicator({0,1,2,3,4,5,6});
  lattice.defineRhoU(     bulkIndicator, densityProfile, u);
  lattice.iniEquilibrium( bulkIndicator, densityProfile, u);
  lattice.initialize();

  clout << "setInitialValues ... OK" << std::endl;
}

// write data to file system
void setPlotData( MyCase& myCase,
                  std::size_t iT,
                  Gnuplot<MyCase::value_t>& gplot_l2_abs,
                  MyCase::value_t Lp0 ) {
  OstreamManager clout(std::cout, "plot");
  using T = MyCase::value_t;
  auto& lattice     = myCase.getLattice(names::NavierStokes{});
  auto& geometry    = myCase.getGeometry();
  auto& parameters  = myCase.getParameters();
  auto& converter   = lattice.getUnitConverter();

  const size_t iTout  = parameters.get<parameters::INTERVAL_TABLE>();
  const size_t  iTlog       = std::min(
    parameters.get<parameters::LATTICE_VTK_ITER_T>(),
    parameters.get<parameters::LATTICE_STAT_ITER_T>());
  const Vector extent = parameters.get<parameters::CORE_EXTENT>();
  const Vector origin = -0.5*extent;

  if (iT % iTout == 0) {
    SuperIndicatorFfromIndicatorF3D<T> domainFluid( new IndicatorCuboid3D<T>( extent, origin ), geometry );
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    T Lpi = L2Norm<3, T, MyCase::descriptor_t>(lattice, converter, domainFluid);
    gplot_l2_abs.setData(T(iT), Lpi / Lp0);
    if ( iT % iTlog == 0 ) {
      clout << "step=" << iT << std::setw(4) << "; Lp2 = " << Lpi
        << "; Lp2/Lp2_0 = " << Lpi/Lp0 << std::endl;
    }
  }
}

void getGraphicalResults( MyCase& myCase,
                          std::size_t iT ) {

  /// Write vtk plots every 0.3 seconds (of phys. simulation time)
  using T = MyCase::value_t;
  auto&             parameters  = myCase.getParameters();
  auto&             lattice     = myCase.getLattice(names::NavierStokes{});
  auto&             converter   = lattice.getUnitConverter();
  const T           amplitude   = parameters.get<parameters::PULSE_AMPLITUDE>();
  const std::string suffix      = parameters.get<parameters::OUTDIR_SUFFIX>();
  const std::string name        = "gausspulse3d" + suffix;
  const size_t      iTgraph     = parameters.get<parameters::LATTICE_VTK_ITER_T>();

  SuperVTMwriter3D<T> vtmWriter(name);
  if (iT == 0) {
    vtmWriter.createMasterFile();
  }

  if ( iT % iTgraph == 0 ) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);
    lattice.scheduleBackgroundOutputVTK([&, name, iT](auto task) {
      SuperVTMwriter3D<T>        vtmWriter(name);
      SuperLatticePhysVelocity3D velocityF(lattice, converter);
      SuperLatticePhysPressure3D pressureF(lattice, converter);
      vtmWriter.addFunctor(velocityF);
      vtmWriter.addFunctor(pressureF);
      task(vtmWriter, iT);
    });

    // output pressure image
    SuperLatticePhysPressure3D pressure(lattice, converter);
    BlockReduction3D2D<T>      pressureReduction(pressure, Vector<T,3>({0, 0, 1}));
    heatmap::plotParam<T>      jpeg_ParamP;
    jpeg_ParamP.maxValue       = converter.getPhysPressure(+amplitude / 200);
    jpeg_ParamP.minValue       = converter.getPhysPressure(-amplitude / 200);
    jpeg_ParamP.colour         = "rainbow";
    jpeg_ParamP.fullScreenPlot = true;
    jpeg_ParamP.name           = name + "_pressure_" + std::to_string(iT) + ".jpeg";
    heatmap::write(pressureReduction, iT, jpeg_ParamP);

    std::stringstream ss;
    ss << std::setw(4) << std::setfill('0') << iT << "_" << name;
    T dist        = converter.getPhysDeltaX();
    T ndatapoints = converter.getResolution(); // number of data points along line plot
    AnalyticalFfromSuperF3D<T> pressure_interpolation(pressure, true, true);
    T                          pmin(converter.getPhysPressure(-amplitude / 50));
    T                          pmax(converter.getPhysPressure(+amplitude / 50));
    linePlot<3,T>(pressure_interpolation, ndatapoints, dist,
                  "pressure_hline_" + ss.str(), "pressure [PU]",
                  horizontal, false, true, pmin, pmax);
    linePlot<3,T>(pressure_interpolation, ndatapoints, dist,
                  "pressure_vline_" + ss.str(), "pressure [PU]",
                  vertical, false, true, pmin, pmax);
    linePlot<3,T>(pressure_interpolation, ndatapoints, dist,
                  "pressure_diagonal_" + ss.str(), "pressure [PU]",
                  diagonal2d, false, true, pmin, pmax);
  }
}

void simulate(MyCase& myCase) {
  OstreamManager clout(std::cout, "simulate");
  clout << "simulate ..." << std::endl;
  using T = MyCase::value_t;
  auto&         parameters  = myCase.getParameters();
  auto&         geometry    = myCase.getGeometry();
  auto&         lattice     = myCase.getLattice(names::NavierStokes{});
  auto&         converter   = lattice.getUnitConverter();
  size_t        iTmaxT      = converter.getLatticeTime(parameters.get<parameters::PHYS_MAX_T>());
  size_t        iTmax       = parameters.get<parameters::LATTICE_MAX_T>()+1;
  iTmax = ( iTmax > 0 ) ? std::min(iTmax, iTmaxT) : iTmaxT;
  const size_t  iTlog       = std::min(parameters.get<parameters::LATTICE_VTK_ITER_T>(), parameters.get<parameters::LATTICE_STAT_ITER_T>());

  // === Initialize pressure L2 norm plot
  const Vector extent = parameters.get<parameters::CORE_EXTENT>();
  const Vector origin = -0.5*extent;  // === calculate output intervals
  SuperIndicatorFfromIndicatorF3D<T> domainFluid( new IndicatorCuboid3D<T>( extent, origin ), geometry );
  Gnuplot<T> gplot_l2_abs("l2_absolute");
  gplot_l2_abs.setLabel("time []", "absolute L2 norm []");
  T Lp0 = L2Norm<3, T, MyCase::descriptor_t>( lattice, converter, domainFluid );
  const bool doGraphicalOutput = parameters.get<parameters::DO_GRAPHICAL_OUTPUT>();
  const bool doL2Plot = parameters.get<parameters::DO_L2_PLOT>();

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    if ( doL2Plot ) setPlotData(myCase, iT, gplot_l2_abs, Lp0);
    if ( doGraphicalOutput ) getGraphicalResults(myCase, iT);

    lattice.collideAndStream();

    if ( iT%iTlog == 0 ) {
      lattice.getStatistics().print(iT, converter.getPhysTime(iT));
      timer.print(iT);
    }
  }

  gplot_l2_abs.setYrange(1e-3, 1);
  gplot_l2_abs.setLogScale(2);
  gplot_l2_abs.writePNG(-1, -1, "gplot_l2_abs");

  timer.stop();
  timer.printSummary();
  clout << "simulate ... OK" << std::endl;
}

void setGetParameters( MyCase::ParametersD& myCaseParameters, int& argc, char** argv ) {

  using namespace olb::parameters;
  // DOMAIN
  myCaseParameters.set<LATTICE_CHAR_VELOCITY>(  defaultMa );
  myCaseParameters.set<CORE_EXTENT          >( {defaultL,defaultL,defaultL} );
  // PULSE
  myCaseParameters.set<PULSE_AMPLITUDE    >( defaultAmp );
  myCaseParameters.set<PULSE_PHYS_VELOCITY>( {1.,0.,0.} );
  myCaseParameters.set<PULSE_ALPHA        >( defaultAlpha );
  // BOUNDARY CONDITIONS
  myCaseParameters.set<FAR_FIELD_TYPE     >(        CBC );
  myCaseParameters.set<SPONGE_STRENGTH    >(         1. );
  myCaseParameters.set<SPONGE_DEPTH_LU    >(         20 );
  myCaseParameters.set<ETERNALSCALE       >(         3. );
  myCaseParameters.set<CBC_SIGMA          >( defaultSigma );
  // TIMING AND OUTPUTS
  myCaseParameters.set<DO_GRAPHICAL_OUTPUT>(      false );
  myCaseParameters.set<DO_L2_PLOT         >(       true );
  myCaseParameters.set<WAIT_DEBUG         >(      false );
  myCaseParameters.set<LATTICE_VTK_ITER_T >(         25 );
  myCaseParameters.set<LATTICE_STAT_ITER_T>(         25 );
  myCaseParameters.set<INTERVAL_TABLE     >(          5 );
  myCaseParameters.set<LATTICE_MAX_T      >(        150 );
  // STANDARD PARAMETERS
  myCaseParameters.set<PHYS_MAX_T         >(         .5 );
  myCaseParameters.set<PHYS_CHAR_VELOCITY >(        1.0 );
  myCaseParameters.set<PHYS_CHAR_VISCOSITY>(      0.001 );
  myCaseParameters.set<PHYS_CHAR_DENSITY  >(        1.0 );
  myCaseParameters.set<RESOLUTION         >( defaultRes );
  myCaseParameters.set<PHYS_DELTA_T       >([](MyCase::ParametersD& p) {
    return p.get<LATTICE_CHAR_VELOCITY>() / p.get<PHYS_CHAR_VELOCITY>() * p.get<PHYS_DELTA_X>();
  } );
  myCaseParameters.set<PHYS_CHAR_LENGTH>([](MyCase::ParametersD& p) {
    return p.get<CORE_EXTENT>()[0];
  });
  myCaseParameters.set<PHYS_DELTA_X>([](MyCase::ParametersD& p) {
    return p.get<PHYS_CHAR_LENGTH>() / p.get<RESOLUTION>();
  });
}

void setOutDir( MyCase::ParametersD& myCaseParameters ) {
  OstreamManager clout(std::cout, "setOutDir");
  using T = MyCase::value_t;
  FarFieldType  farFieldType    = myCaseParameters.get<parameters::FAR_FIELD_TYPE>();
  T             eternalscale    = myCaseParameters.get<parameters::ETERNALSCALE>();
  T             amplitude       = myCaseParameters.get<parameters::PULSE_AMPLITUDE>();
  Vector        coreExtent      = myCaseParameters.get<parameters::CORE_EXTENT>();
  size_t        res             = myCaseParameters.get<parameters::RESOLUTION>();
  size_t        spongeDepthLU   = myCaseParameters.get<parameters::SPONGE_DEPTH_LU>();
  T             spongeStrength  = myCaseParameters.get<parameters::SPONGE_STRENGTH>();
  T             cbcSigma        = myCaseParameters.get<parameters::CBC_SIGMA>();
  T             Ma              = myCaseParameters.get<parameters::LATTICE_CHAR_VELOCITY>();
  Vector        uInfty          = myCaseParameters.get<parameters::PULSE_PHYS_VELOCITY>();
  std::string   outDir          = myCaseParameters.get<parameters::OUTPUT_DIRECTORY>();
  std::string   suffix          = myCaseParameters.get<parameters::OUTDIR_SUFFIX>();

  std::stringstream outDirMod, caseSuffix;
  if ( outDir != "" ) outDirMod << outDir;
  else                outDirMod << "./tmp";
  switch ( farFieldType ) {
    case ETERNAL:
      clout << "Far field boundary condition is solved by just extending the domain to " << eternalscale << " times" << std::endl;
      caseSuffix << "_eternal"; break;
    case PERIODIC:
      clout << "Far field boundary condition type specified to periodic." << std::endl;
      caseSuffix << "_periodic"; break;
    case LOCAL:
      clout << "Far field boundary condition type specified to local." << std::endl;
      caseSuffix << "_local"; break;
    case INTERPOLATED:
      clout << "Far field boundary condition type specified to interpolated." << std::endl;
      caseSuffix << "_interpolated"; break;
    case SPONGE:
      clout << "Far field boundary condition type specified to sponge." << std::endl;
      caseSuffix << "_sponge"; break;
    case CBC_SPONGE:
      clout << "Far field boundary condition type specified to sponge with characteristic outlet." << std::endl;
      caseSuffix << "_cbcSponge"; break;
    case CBC:
      clout << "Far field boundary condition type specified to characteristic." << std::endl;
      caseSuffix << "_cbc"; break;
    case COARSE:
      clout << "Far field boundary condition type specified to coarser grid." << std::endl;
      caseSuffix << "_coarse"; break;
  }

  if ( res != defaultRes && farFieldType != COARSE ) caseSuffix << "_res" << res;
  if ( res != size_t(defaultRes/2.) && farFieldType == COARSE ) caseSuffix << "_res" << res;
  if ( cbcSigma != defaultSigma ) caseSuffix << "_sigma" << cbcSigma;
  if ( Ma != defaultMa )                caseSuffix << "_Ma" << Ma;
  if ( uInfty[0] != 1. || uInfty[1] != 0. || uInfty[2] != 0. ) caseSuffix << "_u" << uInfty[0] << "-" << uInfty[1] << "-" << uInfty[2];
  if ( amplitude != defaultAmp ) caseSuffix << "_a" << amplitude;
  if ( coreExtent[0] != defaultL )  caseSuffix << "_l" << coreExtent[0];
  if ( coreExtent[1] != coreExtent[0] || coreExtent[2] != coreExtent[0] ) caseSuffix << "x" << coreExtent[1] << "x" << coreExtent[2];
  if (farFieldType == ETERNAL) caseSuffix << "_scale" << eternalscale;
  if (farFieldType == SPONGE)   caseSuffix << "_bd" << spongeDepthLU << "x" << spongeStrength;
  if ( suffix != "" ) caseSuffix << "_" << suffix;
  myCaseParameters.set<parameters::OUTDIR_SUFFIX>(caseSuffix.str());
  outDirMod << caseSuffix.str();

  singleton::directories().setOutputDir(outDirMod.str() + "/");
  clout << "Output directory set to " << outDirMod.str() << '\n';
}
