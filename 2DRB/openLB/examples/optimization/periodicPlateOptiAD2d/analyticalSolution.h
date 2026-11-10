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

//Manufactured solution for 2D linear elasticity problem peridic plate
namespace olb {

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
    T x = input[0];
    T y = input[1];
    T pi = T(M_PI);
    T pi2 = pi * pi;

    output[0] = (T(4) / T(625) * (_lambda + _mu) * pi2 * util::cos(T(2) * pi * x) * util::sin(T(2) * pi * y) + T(9) / T(1250) * _mu * pi2 * util::cos(T(2) * pi * x) * util::sin(T(2) * pi * y)) / _conv;

    output[1] = (T(4) / T(625) * (_lambda + _mu) * pi2 * util::sin(T(2) * pi * x) * util::cos(T(2) * pi * y) + T(7) / T(1250) * _mu * pi2 * util::sin(T(2) * pi * x) * util::cos(T(2) * pi * y)) / _conv;

    return true;
  }
};

template <typename T>
class ManufacturedSolutionU2D : public AnalyticalF2D<T, T>
{
  private:
    bool _do;

  public:
    ManufacturedSolutionU2D(bool run) : AnalyticalF2D<T, T>(2),
    _do(run) {}

    bool operator()(T output[], const T input[]) override
    {
      T x = input[0];
      T y = input[1];
      T pi = T(M_PI);

      output[0] = T(9) / T(10000) * util::cos(T(2) * pi * x) * util::sin(T(2) * pi * y);
      output[1] = T(7) / T(10000) * util::sin(T(2) * pi * x) * util::cos(T(2) * pi * y);

      return true;
    }
};

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
      T x = input[0];
      T y = input[1];
      T pi = T(M_PI);

      output[0] = -T(2) / T(625) * (_lambda + _mu) * pi * util::sin(T(2) * pi * x) * util::sin(T(2) * pi * y) - _mu * pi * util::sin(T(2) * pi * x) * util::sin(T(2) * pi * y) / T(2500);

      output[1] = T(2) / T(625) * _mu * util::cos(T(2) * pi * x) * pi * util::cos(T(2) * pi * y);

      output[2] = T(2) / T(625) * _mu * util::cos(T(2) * pi * x) * pi * util::cos(T(2) * pi * y);

      output[3] = -T(2) / T(625) * (_lambda + _mu) * pi * util::sin(T(2) * pi * x) * util::sin(T(2) * pi * y) + _mu * pi * util::sin(T(2) * pi * x) * util::sin(T(2) * pi * y) / T(2500);

      return true;
    }
};

}  // namespace olb
