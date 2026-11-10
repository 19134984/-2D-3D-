/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Florian Kaiser
 *                2026 Johannes Grafen
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

/** \file Manufactured solution of the 2D linear elasticity problem on the
 *  combined ellipse domain, shared by the forward example
 *  `structuralmechanics/linearElasticity/quasistationary/combinedEllipse`
 *  and this parameter identification example.
 */
#pragma once

namespace olb {

/// Geometry of the combined ellipse domain
/**
 * `ellipse1` bounds the solid domain, `ellipse2` and `ellipse3` cut two holes
 * out of it. All extents are scaled by the characteristic physical length.
 */
template <typename T>
class ellipse2D
{
private:
  Vector<T,2> center1;
  Vector<T,2> center2;
  Vector<T,2> center3;

public:
  ellipse2D(T charL) :
    center1(Vector<T,2>(T(0.75) * charL, T(0.75) * charL)),
    center2(Vector<T,2>(T(0.95) * charL, T(0.6 ) * charL)),
    center3(Vector<T,2>(T(0.5 ) * charL, T(0.9 ) * charL)),
    ellipse1(IndicatorEllipse2D<T>(center1, T(0.693) * charL, T(0.548) * charL, T(-20))),
    ellipse2(IndicatorEllipse2D<T>(center2, T(0.1  ) * charL, T(0.134) * charL)),
    ellipse3(IndicatorEllipse2D<T>(center3, T(0.187) * charL, T(0.1  ) * charL))
  { }

  IndicatorEllipse2D<T> ellipse1;
  IndicatorEllipse2D<T> ellipse2;
  IndicatorEllipse2D<T> ellipse3;
};

/// Body force belonging to the manufactured displacement solution
template <typename T>
class ForceField2D : public AnalyticalF2D<T, T>
{
private:
  T _mu;
  T _lambda;
  T _conv;

public:
  ForceField2D(T mu, T lambda, T conv) : AnalyticalF2D<T, T>(2),
    _mu(mu),
    _lambda(lambda),
    _conv(conv) {}

  bool operator()(T output[], const T input[]) override
  {
    T x   = input[0];
    T y   = input[1];
    T pi  = T(M_PI);
    T pi2 = pi * pi;

    output[0] = (-(_lambda + _mu) * (-T(9) / T(2500) * pi2 * (y * y) * util::sin(T(2) * pi * x * y) - T(7) / T(2500) * pi * util::sin(T(2) * pi * y) * x) - _mu * (-T(9) / T(2500) * pi2 * (y * y) * util::sin(T(2) * pi * x * y) - T(9) / T(2500) * pi2 * (x * x) * util::sin(T(2) * pi * x * y))) / _conv;

    output[1] = (-(_lambda + _mu) * (T(9) / T(5000) * pi * util::cos(T(2) * pi * x * y) - T(9) / T(2500) * pi2 * y * x * util::sin(T(2) * pi * x * y) - T(7) / T(2500) * pi2 * util::cos(T(2) * pi * y) * (x * x + T(1))) - _mu * (T(7) / T(5000) * util::cos(T(2) * pi * y) - T(7) / T(2500) * pi2 * util::cos(T(2) * pi * y) * (x * x + T(1)))) / _conv;

    return true;
  }
};

/// Manufactured displacement solution
template <typename T>
class ManufacturedSolutionU2D : public AnalyticalF2D<T, T>
{
protected:
  bool _do;

public:
  ManufacturedSolutionU2D(bool run) : AnalyticalF2D<T, T>(2),
    _do(run) {}

  bool operator()(T output[], const T input[]) override
  {
    T x  = input[0];
    T y  = input[1];
    T pi = T(M_PI);

    output[0] = T(9) / T(10000) * util::sin(T(2) * pi * x * y);
    output[1] = T(7) / T(10000) * util::cos(T(2) * pi * y) * (x * x + T(1));

    return true;
  }
};

/// Stress belonging to the manufactured displacement solution
template <typename T>
class ManufacturedSolutionStress2D : public AnalyticalF2D<T, T>
{
private:
  T _mu;
  T _lambda;

public:
  ManufacturedSolutionStress2D(T mu, T lambda) : AnalyticalF2D<T, T>(4),
    _mu(mu),
    _lambda(lambda) {}

  bool operator()(T output[], const T input[]) override
  {
    T x  = input[0];
    T y  = input[1];
    T pi = T(M_PI);

    output[0] = (_lambda + _mu) * (T(9) / T(5000) * pi * y * util::cos(T(2) * pi * x * y) - T(7) / T(5000) * pi * util::sin(T(2) * pi * y) * (x * x + T(1))) + _mu * (T(9) / T(5000) * pi * y * util::cos(T(2) * pi * x * y) + T(7) / T(5000) * pi * util::sin(T(2) * pi * y) * (x * x + T(1)));

    output[1] = _mu * (T(7) / T(5000) * util::cos(T(2) * pi * y) * x + T(9) / T(5000) * pi * x * util::cos(T(2) * pi * x * y));

    output[2] = _mu * (T(7) / T(5000) * util::cos(T(2) * pi * y) * x + T(9) / T(5000) * pi * x * util::cos(T(2) * pi * x * y));

    output[3] = (_lambda + _mu) * (T(9) / T(5000) * pi * y * util::cos(T(2) * pi * x * y) - T(7) / T(5000) * pi * util::sin(T(2) * pi * y) * (x * x + T(1))) + _mu * (-T(7) / T(5000) * pi * util::sin(T(2) * pi * y) * (x * x + T(1)) - T(9) / T(5000) * pi * y * util::cos(T(2) * pi * x * y));

    return true;
  }
};

}  // namespace olb
