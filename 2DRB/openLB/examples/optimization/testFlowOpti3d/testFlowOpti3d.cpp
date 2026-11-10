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

#include "case.h"

/// Steps 2-8: Set up and run a simulation
int main(int argc, char* argv[])
{
  // Step 2: Set parameters
  initialize(&argc, &argv);

  MyCase::ParametersD myCaseParameters;
  {
    using namespace parameters;
    myCaseParameters.set<DOMAIN_EXTENT          >({1.0, 1.0, 1.0});
    myCaseParameters.set<PHYS_CHAR_VELOCITY     >(            1.0);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY    >(            0.1);
    myCaseParameters.set<PHYS_CHAR_DENSITY      >(            1.0);
    myCaseParameters.set<PHYS_MAX_T             >(            6.0);
    myCaseParameters.set<RESOLUTION             >(             11);
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(            0.8);
    myCaseParameters.set<ADJOINT_TYPE           >( AdjointType::continuous );
  }

  OptimizerParametersD<MyCase::value_t> optiParameters;
  {
    using namespace parameters;
    optiParameters.set<optimizer::EPS>(1e-1);
    optiParameters.set<optimizer::MAX_OPTIMIZATION_STEPS>(10);
    optiParameters.set<linesearch::MAX_STEP_ATTEMPTS>(20);
    optiParameters.set<linesearch::STEP_CONDITION>("StrongWolfe");
    optiParameters.set<optimizer::VERBOSE>(true);
  }

  /// Step 3: Create Mesh
  Mesh mesh = createMesh(myCaseParameters);

  /// Step B: Create Cases and prepare
  // Prepare controlled case
  MyCase myCase(myCaseParameters, mesh);
  prepareGeometry(myCase);
  prepareLattice(myCase);

  // Create adjoint case for gradient computation
  MyCase adjointCase(myCaseParameters, mesh);
  prepareGeometry(adjointCase);

  // Compute solution for the objective functional
  MyCase referenceCase(myCaseParameters, mesh);
  prepareGeometry(referenceCase);
  prepareLattice(referenceCase);
  setInitialValues(referenceCase);
  simulate(referenceCase);

  /// Step C: Create OptiCase and set cases
  MyOptiCase optiCase;
  optiCase.setCase<Controlled>(myCase);
  optiCase.setCase<Adjoint>(adjointCase);
  optiCase.setCase<Reference>(referenceCase);

  /// Step D: Set initial control
  setInitialControl(optiCase);

  /// Step E: Define objective routine
  optiCase.setObjective(computeObjective);

  /// Step F: Define derivative routine
  optiCase.setDerivative(computeDerivative);

  /// Step G: Create an Optimizer
  OptimizerLBFGS<MyOptiCase::value_t,std::vector<MyOptiCase::value_t>> optimizer(
    optiCase.getController().size(), optiParameters);
  optimizer.optimize(optiCase);
}
