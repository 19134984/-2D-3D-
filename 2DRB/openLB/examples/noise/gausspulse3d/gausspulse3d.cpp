/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
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

/* gausspulse3d.cpp:
 * Benchmark from C. Bogey and C. Bailly, “Three-dimensional non-reﬂective
 * boundary conditions for acoustic simulations: far ﬁeld formulation and
 * validation test cases,” vol. 88, 2002.
 * Method from H. Xu and P. Sagaut, “Analysis of the absorbing layers for
 * the weakly-compressible lattice Boltzmann methods,” Journal of Computational Physics, vol. 245, pp. 14–42, Jul. 2013, doi: 10.1016/j.jcp.2013.02.051.
 * Use the option --BOUNDARY_CONDITION to vary how the far field is modeled:
 *  0=eternal (3x domain size to capture the vanishing pulse as ideal reference)
 *  1=periodic,
 *  2=local equilibrium BC,
 *  3=damping (with periodic boundaries)
 */

#include "gausspulse3d.h"

int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");
  clout << "Starting gausspulse3d ..." << std::endl;

  /// === Step 2a: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  // === Step 2b: set boundarytype depending on input values ===
  setGetParameters(myCaseParameters, argc, argv);
  // === Step 2c: set output directory depending on input values ===
  setOutDir(myCaseParameters);

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Lattice ===
  prepareLattice(myCase);

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValues(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);
}
