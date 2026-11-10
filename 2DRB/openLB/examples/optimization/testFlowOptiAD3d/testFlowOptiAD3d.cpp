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
  /// @par Contents
  /// @li Step 2: Initialization
  initialize(&argc, &argv);
  MyCase::ParametersD myCaseParametersD;
  {
    using namespace parameters;
    myCaseParametersD.template set<PHYS_CHAR_VELOCITY     >(        1.0);
    myCaseParametersD.template set<PHYS_CHAR_LENGTH       >(        1.0);
    myCaseParametersD.template set<PHYS_CHAR_VISCOSITY    >(        0.1);
    myCaseParametersD.template set<PHYS_CHAR_DENSITY      >(        1.0);
    myCaseParametersD.template set<PHYS_MAX_T             >(        6.0);
    myCaseParametersD.template set<RESOLUTION             >(         11);
    myCaseParametersD.template set<LATTICE_RELAXATION_TIME>(        0.8);
    myCaseParametersD.template set<OBJECTIVE_TYPE>(ObjectiveType::velocity);
  }
  auto myADfCaseParametersD = myCaseParametersD.copyAs<MyADfCase::value_t>();

  OptimizerParametersD<MyCase::value_t> optiParameters;
  {
    using namespace parameters;
    optiParameters.set<optimizer::EPS>(1e-5);
    optiParameters.set<optimizer::MAX_OPTIMIZATION_STEPS>(10);
    optiParameters.set<linesearch::MAX_STEP_ATTEMPTS>(20);
    optiParameters.set<linesearch::STEP_CONDITION>("Wolfe");
    optiParameters.set<optimizer::VERBOSE>(true);
  }

  /// @li Step 3: Create Mesh
  Mesh mesh = createMesh(myCaseParametersD);
  Mesh aDfMesh = createMesh(myADfCaseParametersD);

  // ==== BELOW HERE OPTIMIZATION SPECIFIC ====
  /// @li Step B: Create Cases and prepare
  // Prepare controlled case
  MyCase myCase(myCaseParametersD, mesh);
  prepareGeometry(myCase);

  // Create ADf-typed case for gradient computation
  MyADfCase aDfCase(myADfCaseParametersD, aDfMesh);
  prepareGeometry(aDfCase);

  // Compute solution for the objective functional
  MyCase referenceCase(myCaseParametersD, mesh);
  prepareGeometry(referenceCase);
  prepareLattice(referenceCase);
  setInitialValues(referenceCase);
  simulate(referenceCase);

  /// @li Step C: Create OptiCase and set cases
  MyOptiCase optiCase;
  optiCase.setCase<Controlled>(myCase);
  optiCase.setCase<Derivatives>(aDfCase);
  optiCase.setCase<Reference>(referenceCase);

  /// @li Step D: Set initial control
  setInitialControl(optiCase);

  /// @li Step E: Define objective routine
  optiCase.setObjective([](auto& optiCase, auto name) {
    return computeObjective(optiCase, name);
  });

  /// @li Step G: Create an Optimizer
  OptimizerLBFGS<MyCase::value_t,std::vector<MyCase::value_t>> optimizer(
    optiCase.getController().size(), optiParameters);

  /// @li Step H: Optimize
  optimizer.optimize(optiCase);
}
