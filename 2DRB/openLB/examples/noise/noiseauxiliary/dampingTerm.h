/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Philipp Spelten, Sylvain Franiatte
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

/* dampingTerm.h:
 * The newtonian cooling term suggested in M. Israeli, S.A. Orszag,
 * Approximation of radiation boundary conditions,
 * J. Comput. Phys. 41 (1981) 115–135.
*/

#ifndef DAMPING_TERM_H
#define DAMPING_TERM_H

namespace olb {

template <unsigned ndim, typename T, unsigned N>
class DampingTerm : public AnalyticalF<ndim, T, T> {
protected:
  T _A = 4;         // amplitude factor
  T _n = 4;         // order of polynomial
  size_t _direction[N];  // outside boundary position in direction
  T _xBoundary[N];  // outside boundary position in direction
  T _bdPU;          // boundary depth in PU

public:
  DampingTerm(  size_t direction[N],  // directions in which boundaries are pointed
                T xBoundary[N],       // outer positions along directions
                T boundaryDepthPU     // depth of damping layer
              ) : AnalyticalF<ndim, T, T>(1)
      , _bdPU(boundaryDepthPU)
  {
    for ( size_t nBoundary = 0; nBoundary < N; ++nBoundary ) {
      _xBoundary[nBoundary] = xBoundary[nBoundary];
      _direction[nBoundary] = direction[nBoundary];
    }
  };

  bool operator()(T output[], const T input[]) override
  {
    T x=1;
    for ( size_t nBoundary = 0; nBoundary < N; ++nBoundary ) {
      // normalize x to [0,1]
      T xi = std::abs(input[_direction[nBoundary]] - _xBoundary[nBoundary]) / _bdPU;
      if ( x < 1 ) {
        x = std::min( x, xi );
      }
    }
    if ( x < 1 ) {
      output[0] = -_A * (((std::pow(x-1, _n)) * (-x) * (_n + 1) * (_n + 2)) / std::pow(-_bdPU, _n + 2));
      // normalize sigma to [0,1]
      output[0] /= 9.8304 / std::pow(_bdPU, _n + 2); // calculated as max of function (at X=0.8) for A=n=4
      return true;
    }
    output[0] = 0;
    return true;
  };
};

} // namespace olb

#endif
