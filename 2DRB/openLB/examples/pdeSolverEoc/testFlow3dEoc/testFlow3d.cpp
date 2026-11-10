/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006-2025 Mathias J. Krause, Fabian Klemens,
 *  Julius Jessberger, Shota Ito
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

#include "testFlow3d.h"

using T = MyCase::value_t;

/// Initialize gnuplot
static Gnuplot<T> gplot(
  "eoc",
  false,
  "set terminal png size 720, 720 font 'Arial,10'",
  Gnuplot<T>::LOGLOG,
  Gnuplot<T>::LINREG);

int main(int argc, char* argv[])
{
  /// @li Step 2: Set Parameters
  initialize(&argc, &argv);
  MyCase::ParametersD myCaseParametersD;
  {
    using namespace parameters;
    myCaseParametersD.set<RESOLUTION_START            >(                   11);
    myCaseParametersD.set<NUMBER_OF_RUNS              >(                    3);
    myCaseParametersD.set<DOMAIN_EXTENT               >(      {1.0, 1.0, 1.0});
    myCaseParametersD.set<PHYS_CHAR_VELOCITY          >(                  1.0);
    myCaseParametersD.set<PHYS_CHAR_VISCOSITY         >(                 0.01);
    myCaseParametersD.set<PHYS_CHAR_DENSITY           >(                    1);
    myCaseParametersD.set<PHYS_MAX_T                  >(                 50.0);
    myCaseParametersD.set<PHYS_BOUNDARY_VALUE_UPDATE_T>(                0.001);
    myCaseParametersD.set<RESOLUTION                  >(                   11);
    myCaseParametersD.set<LATTICE_CHAR_VELOCITY       >(                 0.07);
    myCaseParametersD.set<ERROR_VELOCITY_L1           >(                   0.);
    myCaseParametersD.set<ERROR_VELOCITY_L2           >(                   0.);
    myCaseParametersD.set<ERROR_VELOCITY_LINF         >(                   0.);
    myCaseParametersD.set<ERROR_PRESSURE_L1           >(                   0.);
    myCaseParametersD.set<ERROR_PRESSURE_L2           >(                   0.);
    myCaseParametersD.set<ERROR_PRESSURE_LINF         >(                   0.);
    myCaseParametersD.set<ERROR_STRAIN_RATE_L1        >(                   0.);
    myCaseParametersD.set<ERROR_STRAIN_RATE_L2        >(                   0.);
    myCaseParametersD.set<ERROR_STRAIN_RATE_LINF      >(                   0.);
    myCaseParametersD.set<ERROR_DISSIPATION_L1        >(                   0.);
    myCaseParametersD.set<ERROR_DISSIPATION_L2        >(                   0.);
    myCaseParametersD.set<ERROR_DISSIPATION_LINF      >(                   0.);
    myCaseParametersD.set<ERROR_STRESS_L1             >(                   0.);
    myCaseParametersD.set<ERROR_STRESS_L2             >(                   0.);
    myCaseParametersD.set<ERROR_STRESS_LINF           >(                   0.);
    myCaseParametersD.set<GEOMETRY_TYPE               >( GeometryType::sphere);
    myCaseParametersD.set<BOUNDARY_TYPE               >(BoundaryType::bouzidi);
  }

  T _errorVelL1[myCaseParametersD.get<parameters::NUMBER_OF_RUNS>()];
  T _errorVelL2[myCaseParametersD.get<parameters::NUMBER_OF_RUNS>()];
  T _errorVelLInf[myCaseParametersD.get<parameters::NUMBER_OF_RUNS>()];
  T _errorPreL1[myCaseParametersD.get<parameters::NUMBER_OF_RUNS>()];
  T _errorPreL2[myCaseParametersD.get<parameters::NUMBER_OF_RUNS>()];
  T _errorPreLInf[myCaseParametersD.get<parameters::NUMBER_OF_RUNS>()];
  T _errorStrL1[myCaseParametersD.get<parameters::NUMBER_OF_RUNS>()];
  T _errorStrL2[myCaseParametersD.get<parameters::NUMBER_OF_RUNS>()];
  T _errorStrLInf[myCaseParametersD.get<parameters::NUMBER_OF_RUNS>()];
  T _errorDisL1[myCaseParametersD.get<parameters::NUMBER_OF_RUNS>()];
  T _errorDisL2[myCaseParametersD.get<parameters::NUMBER_OF_RUNS>()];
  T _errorDisLInf[myCaseParametersD.get<parameters::NUMBER_OF_RUNS>()];
  T _errorStressL1[myCaseParametersD.get<parameters::NUMBER_OF_RUNS>()];
  T _errorStressL2[myCaseParametersD.get<parameters::NUMBER_OF_RUNS>()];
  T _errorStressLInf[myCaseParametersD.get<parameters::NUMBER_OF_RUNS>()];

  for(int i = 0; i < myCaseParametersD.get<parameters::NUMBER_OF_RUNS>(); ++i) {
    myCaseParametersD.set<parameters::RESOLUTION>(
      myCaseParametersD.get<parameters::RESOLUTION_START>() + i*10);
    myCaseParametersD.set<parameters::LATTICE_CHAR_VELOCITY>(
      myCaseParametersD.get<parameters::DOMAIN_EXTENT>()[0] /
      myCaseParametersD.get<parameters::RESOLUTION>());

    /// @li Step 3: Create Mesh
    Mesh mesh = createMesh(myCaseParametersD);

    /// @li Step 4: Create Case
    MyCase myCase(myCaseParametersD, mesh);

    /// @li Step 5: Prepare Geometry
    prepareGeometry(myCase);

    /// @li Step 6: Prepare Lattice
    prepareLattice(myCase);

    /// @li Step 7: Definition of Initial and Boundary Values and Fields
    setInitialValues(myCase);

    /// @li Step 8: Simulate
    simulate(myCase);

    // Save computed error norms
    _errorVelL1[i] = myCaseParametersD.get<parameters::ERROR_VELOCITY_L1>();
    _errorVelL2[i] = myCaseParametersD.get<parameters::ERROR_VELOCITY_L2>();
    _errorVelLInf[i] = myCaseParametersD.get<parameters::ERROR_VELOCITY_LINF>();
    _errorPreL1[i] = myCaseParametersD.get<parameters::ERROR_PRESSURE_L1>();
    _errorPreL2[i] = myCaseParametersD.get<parameters::ERROR_PRESSURE_L2>();
    _errorPreLInf[i] = myCaseParametersD.get<parameters::ERROR_PRESSURE_LINF>();
    _errorStrL1[i] = myCaseParametersD.get<parameters::ERROR_STRAIN_RATE_L1>();
    _errorStrL2[i] = myCaseParametersD.get<parameters::ERROR_STRAIN_RATE_L2>();
    _errorStrLInf[i] = myCaseParametersD.get<parameters::ERROR_STRAIN_RATE_LINF>();
    _errorDisL1[i] = myCaseParametersD.get<parameters::ERROR_DISSIPATION_L1>();
    _errorDisL2[i] = myCaseParametersD.get<parameters::ERROR_DISSIPATION_L2>();
    _errorDisLInf[i] = myCaseParametersD.get<parameters::ERROR_DISSIPATION_LINF>();
    _errorStressL1[i] = myCaseParametersD.get<parameters::ERROR_STRESS_L1>();
    _errorStressL2[i] = myCaseParametersD.get<parameters::ERROR_STRESS_L2>();
    _errorStressLInf[i] = myCaseParametersD.get<parameters::ERROR_STRESS_LINF>();
  }

  for(int i = 0; i < myCaseParametersD.get<parameters::NUMBER_OF_RUNS>(); ++i){
    int N = myCaseParametersD.get<parameters::RESOLUTION_START>() + 10 * i;
    gplot.setData(T(N),
                 {_errorVelL1[i], _errorVelL2[i], _errorVelLInf[i],
                  _errorPreL1[i], _errorPreL2[i], _errorPreLInf[i],
                  _errorStrL1[i], _errorStrL2[i], _errorStrLInf[i],
                  _errorDisL1[i], _errorDisL2[i], _errorDisLInf[i],
                  _errorStressL1[i], _errorStressL2[i], _errorStressLInf[i]},
                 {"Velocity L1", "Velocity L2", "Velocity LInf",
                  "Pressure L1", "Pressure L2", "Pressure LInf",
                  "Strain rate L1", "Strain rate L2", "Strain rate LInf",
                  "Dissipation L1", "Dissipation L2", "Dissipation LInf",
                  "Stress L1", "Stress L2", "Stress LInf"},
         "top right",
         {'p', 'p', 'p', 'p', 'p', 'p', 'p', 'p', 'p', 'p', 'p', 'p', 'p', 'p', 'p'});
  }
  gplot.writePNG();
}
