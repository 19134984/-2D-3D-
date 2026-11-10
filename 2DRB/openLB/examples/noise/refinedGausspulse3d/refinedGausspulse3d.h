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

#include "olb.h"

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
void setFineInitialValues(MyCase::ParametersD& parameters,
                          SuperLattice<T, DESCRIPTOR>& lattice,
                          SuperGeometry<T, DESCRIPTOR::d>& geometry
) {
  OstreamManager clout(std::cout, "setInitialValues");
  clout << std::endl << "setInitialValues ..." << std::endl;
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
  Vector originDomain = parameters.get<parameters::DOMAIN_ORIGIN>() - 3.*physDx;
  Vector extentDomain = parameters.get<parameters::DOMAIN_EXTENT>() + 6.*physDx;
  SuperIndicatorFfromIndicatorF3D<T> domain( new IndicatorCuboid3D<T>( extentDomain, originDomain ), geometry);
  lattice.defineRhoU(     domain, densityProfile, u);
  lattice.iniEquilibrium( domain, densityProfile, u);
  lattice.initialize();

  clout << "setInitialValues ... OK" << std::endl;
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
void getInitialSetup( SuperLattice<T, DESCRIPTOR>& lattice,
                      SuperGeometry<T, DESCRIPTOR::d>& geometry,
                      std::string suffix ) {
  OstreamManager clout(std::cout, "getInitialSetup");
  clout << "getInitialSetup ..." << std::endl;

  auto&       converter = lattice.getUnitConverter();
  std::string name = "initialization" + suffix;
  size_t iT = 0;
  SuperVTMwriter3D<T> vtmWriter(name);
  vtmWriter.createMasterFile();

  lattice.scheduleBackgroundOutputVTK([&, name, iT](auto task) {
    SuperVTMwriter3D<T>        vtmWriter(name);
    SuperLatticePhysVelocity3D velocityF(lattice, converter);
    SuperLatticePhysPressure3D pressureF(lattice, converter);
    // leave all in LU (1. factor)
    SuperLatticePhysField<T, DESCRIPTOR, descriptors::DENSITY>  density( lattice, 1. );
    SuperGeometryF<T,3> materials( geometry );
    density.getName() = "densityField";
    materials.getName() = "geometry";
    vtmWriter.addFunctor( density );
    vtmWriter.addFunctor(velocityF);
    vtmWriter.addFunctor(pressureF);
    vtmWriter.addFunctor(materials);
    task(vtmWriter, iT);
  });

  clout << "getInitialSetup ... OK" << std::endl;
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
void setFinePlotData( MyCase& myCase,
                      SuperLattice<T, DESCRIPTOR>& lattice,
                      SuperGeometry<T, DESCRIPTOR::d>& geometry,
                      std::size_t iT,
                      Gnuplot<T>& gplot_l2_abs,
                      T Lp0 ) {
  auto& parameters    = myCase.getParameters();
  auto& converter     = lattice.getUnitConverter();
  const Vector extent = parameters.get<parameters::CORE_EXTENT>();
  const Vector origin = -0.5*extent;

  SuperIndicatorFfromIndicatorF3D<T> domainFluid( new IndicatorCuboid3D<T>( extent, origin ), geometry );
  lattice.setProcessingContext(ProcessingContext::Evaluation);
  T Lpi = L2Norm<3, T, DESCRIPTOR>(lattice, converter, domainFluid);
  gplot_l2_abs.setData(T(iT), Lpi / Lp0);
}

template<concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
void getFineGraphicalResults( MyCase& myCase,
                              std::size_t iT,
                              SuperLattice<T, DESCRIPTOR>& lattice ) {
  auto&             converter   = lattice.getUnitConverter();
  auto&             parameters  = myCase.getParameters();
  const T           amplitude   = parameters.get<parameters::PULSE_AMPLITUDE>();
  const size_t      iTgraph     = parameters.get<parameters::LATTICE_VTK_ITER_T>();
  const std::string suffix      = parameters.get<parameters::OUTDIR_SUFFIX>();
  const std::string name        = "gausspulse3d" + suffix + "_core";

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
    T                          dist        = converter.getPhysDeltaX();
    T                          ndatapoints = converter.getResolution(); // number of data points on line
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
