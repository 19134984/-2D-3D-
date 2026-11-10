/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Tim Bingert, Michael Rennick
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

/** flatInterfaceErrors2d.h
 * Provides error computations.
 */

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;

using V = FLOATING_POINT_TYPE;

// variables for eoc analysis
V orderParameterL1RelError = 0;
V orderParameterL2RelError = 0;
V orderParameterLinfRelError = 0;
V pressureL1RelError = 0;
V pressureL2RelError = 0;
V pressureLinfRelError = 0;
V velocityLinfAbsError = 0;

template <typename MyCase>
void errorFlatInterface( MyCase& myCase, bool Cahn_const )
{
  OstreamManager clout( std::cout,"error" );

  using T = MyCase::value_t;
  using NSDESCRIPTOR = typename MyCase::template descriptor_t_of<NavierStokes>;
  using PFDESCRIPTOR = typename MyCase::template descriptor_t_of<Component1>;
  auto& geometry = myCase.getGeometry();
  auto& params = myCase.getParameters();

  auto& sLatticeNS = myCase.getLattice(NavierStokes{});
  auto& sLatticePF = myCase.getLattice(Component1{});

  const auto& converter = sLatticeNS.getUnitConverter();

  Vector extent = params.template get<parameters::DOMAIN_EXTENT>();
  const int N = params.template get<parameters::RESOLUTION>();
  const int phaseLength = N/2;
  const T dx = extent[1] / N;
  const T phys_pressure = params.template get<parameters::PHYS_CHAR_PRESSURE>();
  const T w = params.template get<parameters::INTERFACE_WIDTH>();

  int tmp[]= { };
  T result[2]= { };

  AnalyticalConst2D<T,T> p_phys ( phys_pressure );
  AnalyticalConst2D<T,T> zero ( 0. );
  AnalyticalConst2D<T,T> one ( 1. );
  Vector<T,2> u(0., 0.);
  AnalyticalConst2D<T,T> zeroVelocity( u );
  SmoothIndicatorFactoredCuboid2D<T,T> interface_diffuse( {extent[0]/2., (N/2.+0.5)*dx}, 0, phaseLength*dx, w*dx/2, 0, {0,0}, 0, -1. );
  AnalyticalIdentity2D<T,T> phi0_diffuse( one + interface_diffuse );
  Vector<T,2> extend( 2./100.*extent[1], 0.25*extent[1] );
  Vector<T,2> origin1( 0., 0.);
  Vector<T,2> origin2( 0., 3./4.*extent[1]);
  IndicatorCuboid2D<T> ind1( extend, origin1 );
  IndicatorCuboid2D<T> ind2( extend, origin2 );
  SmoothIndicatorCuboid2D<T,T> interface_sharp1( ind1, T(0) );
  SmoothIndicatorCuboid2D<T,T> interface_sharp2( ind2, T(0) );
  AnalyticalIdentity2D<T,T> phi0_sharp( zero + interface_sharp1 + interface_sharp2 );

  SuperLatticeDensity2D<T, PFDESCRIPTOR> phi( sLatticePF );
  SuperLatticePhysIncPressure2D<T, NSDESCRIPTOR> p_hydro( sLatticeNS, converter );
  SuperLatticePhysVelocity2D<T, NSDESCRIPTOR> velocity( sLatticeNS, converter );

  auto indicatorF = geometry.getMaterialIndicator(1);

  //for order parameter
  if (Cahn_const) {
    SuperRelativeErrorL1Norm2D<T> relPhiErrorNormL1(phi, phi0_diffuse, indicatorF);
    relPhiErrorNormL1(result, tmp);
    clout << "phi-L1-error(rel)=" << result[0];
    orderParameterL1RelError = result[0];

    SuperRelativeErrorL2Norm2D<T> relPhiErrorNormL2(phi, phi0_diffuse, indicatorF);
    relPhiErrorNormL2(result, tmp);
    clout << "; phi-L2-error(rel)=" << result[0];
    orderParameterL2RelError = result[0];

    SuperRelativeErrorLinfNorm2D<T> relPhiErrorNormLinf(phi, phi0_diffuse, indicatorF);
    relPhiErrorNormLinf(result, tmp);
    clout << "; phi-Linf-error(rel)=" << result[0] << std::endl;
    orderParameterLinfRelError = result[0];
  }
  else {
    SuperRelativeErrorL1Norm2D<T> relPhiErrorNormL1(phi, phi0_sharp, indicatorF);
    relPhiErrorNormL1(result, tmp);
    clout << "phi-L1-error(rel)=" << result[0];
    orderParameterL1RelError = result[0];

    SuperAbsoluteErrorL2Norm2D<T> relPhiErrorNormL2(phi, phi0_sharp, indicatorF);
    relPhiErrorNormL2(result, tmp);
    clout << "; phi-L2-error(rel)=" << result[0];
    orderParameterL2RelError = result[0];

    SuperAbsoluteErrorLinfNorm2D<T> relPhiErrorNormLinf(phi, phi0_sharp, indicatorF);
    relPhiErrorNormLinf(result, tmp);
    clout << "; phi-Linf-error(rel)=" << result[0] << std::endl;
    orderParameterLinfRelError = result[0];
  }
  //for pressure
  SuperRelativeErrorL1Norm2D<T> relPErrorNormL1(p_hydro, p_phys, indicatorF);
  relPErrorNormL1(result, tmp);
  clout << "p-L1-error(rel)=" << result[0];
  pressureL1RelError = result[0];

  SuperRelativeErrorL2Norm2D<T> relPErrorNormL2(p_hydro, p_phys, indicatorF);
  relPErrorNormL2(result, tmp);
  clout << "; p-L2-error(rel)=" << result[0];
  pressureL2RelError = result[0];

  SuperRelativeErrorLinfNorm2D<T> relPErrorNormLinf(p_hydro, p_phys, indicatorF);
  relPErrorNormLinf(result, tmp);
  clout << "; p-Linf-error(rel)=" << result[0] << std::endl;
  pressureLinfRelError = result[0];

  //for velocity
  SuperAbsoluteErrorLinfNorm2D<T> absUErrorNormLinf(velocity, zeroVelocity, indicatorF);
  absUErrorNormLinf(result, tmp);
  clout << "u-Linf-error(abs)=" << result[0] << std::endl;
  velocityLinfAbsError = result[0];
}
