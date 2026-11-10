/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006, 2007, 2012 Jonas Latt, Mathias J. Krause,
 *  Louis Kronberg, Christian Vorwerk, Bastian Schäffauer
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

/* Spalding.cpp:
 * Test environment for the approximation of Spaldings wall function in turbulentWallModelPostProcessor3D
 * with two separate variations of Newtons method to calculate u+ and u_tau.
 * First, the approximation of u+ is validated by computing f(u+) = u+ ... - y+ = 0
 * with Newton's method for y+ = 1,...,1000.
 * Then, validation of approximating u_tau_correct is done by computing f(u_tau) = u/u_tau + ... - y/nu*u_tau = 0,
 * where u = u_tau_correct * u+ and y = nu*y+/u_tau_correct, for y+ = 1,...,1000.
 * The approximated values, number of iterations and errors e+ = |f(u+)| and e_tau = |u_tau - u_tau_correct| are
 * plotted for every y+ = 1,...,1000.
 */

#include <olb.h>
#include <cmath>
#include <iostream>

using namespace std;
using namespace olb;
using namespace olb::descriptors;
using namespace olb::util;
using namespace olb::graphics;

using T = double;

// Spalding wall function written in form of f(u+) = u+ + ... - y+ for spaldingOpt /= 2 and
// Spalding wall function written in form of f(u_tau) = u/u_tau + ... - y/nu*u_tau for spaldingOpt = 2.
// vector = {T u_tau, T u, T y, T nu, T u+, T y+, T spaldingOpt}
template <typename T, typename Par_f>
T spaldingWallFunction1(Par_f & params) any_platform{
  if(params[6] == 2.){
    params[4] = params[1]/params[0];
    params[5] = params[2]/params[3]*params[0];
  }
  T function = params[4] + 0.1108*(util::exp(0.4*params[4]) - 1. - 0.4*params[4] - 0.5*util::pow(0.4*params[4],2) - (1./6.)*util::pow(0.4*params[4],3)) - params[5];
  return function;
};

// Spalding wall function derivatives f'(u_tau), when k=2, and f'(u+), when k=1, for both cases of Newton's method.
// vector = {T x, T u, T y, T nu, T k}
template <typename T, typename Par_d>
T spaldingDerivative1(Par_d & param) any_platform{
  T derivative = 0.;
  if(param[4] == 2.){
    derivative = (-param[1]/(util::pow(param[0],2))) + 0.1108*( (util::exp(0.4*param[1]/param[0]))*(-0.4*param[1]/(util::pow(param[0],2))) + 0.4*param[1]/(util::pow(param[0],2)) + (util::pow(0.4*param[1]/param[0],2))/param[0] + 0.5*(util::pow(0.4*param[1]/param[0],3))/param[0] ) - (param[2]/param[3]);
  }
  if(param[4] == 1.){
    derivative = 1. + 0.1108*(util::exp(0.4*param[0])*0.4 - 0.4 - 0.4*param[0]*0.4 - 0.2*util::pow(0.4*param[0],2));
  }
  return derivative;
};

// Newton's method to approximate u+ and u_tau with respective function f and derivative d.
template <typename T, typename Funcf,  typename Funcd, typename Par_f, typename Par_d>
T newtonSpalding1(Funcf f, Funcd d, Par_f & argsf, Par_d & argsd, int optArg, T convergeCrit, Gnuplot<T>& gplot, T yplus) any_platform{
  OstreamManager clout( std::cout,"newton" );
  T j = 0;
  T x = argsf[optArg];
  T x_old = argsf[optArg];
  for(int i = 0; i<100; i++){
    j += 1;
    //if(i==99) clout << "i=99" << std::endl;
    x_old = x;
    argsf[optArg] = x;
    argsd[0] = x;
    x -= f(argsf)/d(argsd);
    if(util::abs((x_old - x)/x) <= convergeCrit) break;
    if(std::isnan(x)) break;
  }
  //clout << "j = " << j << std::endl;
  gplot.setData( yplus, { j }, { "j"}, "bottom right", { 'l' } ); // Plot number of iterations for every y+.
  return x;
};

// Main
int main( int argc, char* argv[] ) {

  // Plots:
  initialize( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  static Gnuplot<T> gplotjplus( "plotjplus", Gnuplot<T>::LINEAR );
  gplotjplus.setLabel("y_{plus}", "j");
  static Gnuplot<T> gplotjtau( "plotjtau", Gnuplot<T>::LINEAR );
  gplotjtau.setLabel("y_{plus}", "j");

  static Gnuplot<T> gplotu( "plotu", Gnuplot<T>::LOGLOG );
  gplotu.setLabel("y_{plus}", "u_{plus}");
  static Gnuplot<T> gploterr( "ploterr", Gnuplot<T>::LINEAR );
  gploterr.setLabel("y_{plus}", "error");

  static Gnuplot<T> gplotutau( "plottau", Gnuplot<T>::LINEAR );
  gplotutau.setLabel("y_{plus}", "u_{tau}");
  static Gnuplot<T> gploterrtau( "ploterrtau", Gnuplot<T>::LINEAR );
  gploterrtau.setLabel("y_{plus}", "error");

  // Initialisation of inital guesses for both Newton methods
  T u_plus1 = 0;
  T u_tau = 1.5;

  // Test case 1:
  //T nu = 0.00005;
  //T utau_corr = 0.0500256;

  // Test case 2:
  //T nu = 0.0008;
  //T utau_corr = 0.8005;

  // Test case 3:
  T nu = 0.002;
  T utau_corr = 2.005;

  T err_max = 0.;

  // loop over y_plus
  for(T y_plus = 1; y_plus < 1001; y_plus += 1){

    // Run the simulations

    // Approximation of u+
    u_plus1 = 20.;  // Initial guess of u+.
    Vector<T,7> params_u1 = {0., 0., 0., 0., u_plus1, y_plus, 1.};   // Vector for spalding wall function
    Vector<T,5> params_d1 = {u_plus1, 0., 0., 0., 1.};  // Vector for spalding derivative
    u_plus1 = newtonSpalding1( [](Vector<T,7> par1) { return spaldingWallFunction1<T>(par1); }, [](Vector<T,5> par2) { return spaldingDerivative1<T>(par2); }, params_u1, params_d1, 4, T{0.0001}, gplotjplus, y_plus);
    clout << "u_plus1= " << u_plus1 << std::endl;
    T error_plus = util::abs(u_plus1 + 0.1108*( util::exp(0.4*u_plus1) - 1. - 0.4*u_plus1 - 0.5*util::pow(0.4*u_plus1,2) - (1./6.)*util::pow(0.4*u_plus1,3) ) - y_plus);
    clout << "error uplus= " << error_plus << std::endl;
    if (err_max < error_plus) err_max = error_plus;

    // Approximation of u_tau
    T u_2 = u_plus1*utau_corr; // Define approx. u+ as "exact".
    T y = nu*y_plus/utau_corr;
    Vector<T,7> params_u2 = {u_tau, u_2, y, nu, u_plus1, y*u_tau/nu, 2.};  // Vector for spalding wall function
    Vector<T,5> params_d2 = {u_tau, u_2, y, nu, 2.};  // Vector for spalding derivative
    u_tau = newtonSpalding1( [](Vector<T,7> par1) { return spaldingWallFunction1<T>(par1); }, [](Vector<T,5> par2) { return spaldingDerivative1<T>(par2); }, params_u2, params_d2, 0, T{0.0001}, gplotjtau, y_plus);
    clout << "u_tau = " << u_tau << std::endl;
    T error_tau = util::abs(u_tau - utau_corr);
    clout << "error utau= " << error_tau << std::endl;

    // Gnuplot output
    gplotu.setData( y_plus, { u_plus1 }, { "u_{plus1}"}, "bottom right", { 'l' } );
    gploterr.setData( y_plus, { error_plus }, { "error"}, "bottom right", { 'l' } );

    gplotutau.setData( y_plus, { u_tau }, { "u_{tau}"}, "bottom right", { 'l' } );
    gploterrtau.setData( y_plus, { error_tau }, { "error"}, "bottom right", { 'l' } );
  }

  clout << "u+ err_max= " << err_max << std::endl;

  clout << "Done " << std::endl;
}
