/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2026 Johannes Grafen
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

/// Setup and run the parameter identification
int main(int argc, char* argv[])
{
  /// === Step 1: Initialization ===
  OstreamManager clout(std::cout, "main");

  initialize(&argc, &argv);

  singleton::directories().setOutputDir("./tmp/");

  constexpr double targetE         = 0.1;
  constexpr double initialControlE = 0.15;

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<SOLID_CASE           >(  "continua"      );
    myCaseParameters.set<RESOLUTION           >(  40              );
    myCaseParameters.set<PHYS_CHAR_LENGTH     >(  1.0             );
    myCaseParameters.set<DOMAIN_EXTENT        >(  { 1.5, 1.5 }    );
    myCaseParameters.set<YOUNGS_MODULUS       >(  initialControlE );
    myCaseParameters.set<PHYS_DENSITY         >(  1.0             );
    myCaseParameters.set<POISSON_RATIO        >(  0.7             );
    myCaseParameters.set<DAMPING_FACTOR       >(  1.0             );
    myCaseParameters.set<PHYS_MAX_T           >(  60              );
    myCaseParameters.set<LATTICE_LOG_ITER_T   >(  100             );
    myCaseParameters.set<LATTICE_VTK_ITER_T   >(  1000            );
    myCaseParameters.set<CONVERGENCE_PRECISION>(  1e-7            );
    myCaseParameters.set<CONV_ITER            >(  10              );
    myCaseParameters.set<CONV_VALUE           >(  1               );
    myCaseParameters.set<OPTI_VTK_ENABLED     >(  true            );
    myCaseParameters.set<CONVERGED            >(  false           );
    myCaseParameters.set<INITIAL_CONTROL      >(  initialControlE );
    myCaseParameters.set<ANALYTICAL_REFERENCE >(  true            );
    myCaseParameters.set<TARGET_YOUNGS_MODULUS>(  targetE         );
    myCaseParameters.set<FIXED_TIME_STEP      >(  false           );
    myCaseParameters.set<FORWARD_ONLY         >(  false           );
  }

  // Compute additional derived parameters
  myCaseParameters.set<parameters::PHYS_DELTA_X>(
    myCaseParameters.get<parameters::DOMAIN_EXTENT>()[0] /
    myCaseParameters.get<parameters::RESOLUTION>());
  myCaseParameters.set<parameters::LATTICE_MAX_T>(static_cast<std::size_t>(
    myCaseParameters.get<parameters::PHYS_MAX_T>() /
    (myCaseParameters.get<parameters::PHYS_DELTA_X>() *
     myCaseParameters.get<parameters::PHYS_DELTA_X>()) + 0.5));
  myCaseParameters.set<parameters::ORIGIN>(
    Vector<MyCase::value_t,MyCase::d>{
      myCaseParameters.get<parameters::PHYS_DELTA_X>(),
      myCaseParameters.get<parameters::PHYS_DELTA_X>()});
  myCaseParameters.set<parameters::PHYS_CHAR_DISPLACEMENT>(
    myCaseParameters.get<parameters::PHYS_DELTA_X>());

  // ADf-typed copy for the gradient computation
  auto myADfCaseParameters = myCaseParameters.copyAs<MyADfCase::value_t>();

  // The reference case is the same setup, evaluated at the sought-after modulus
  auto myReferenceParameters = myCaseParameters.copyAs<MyCase::value_t>();
  myReferenceParameters.set<parameters::YOUNGS_MODULUS>(targetE);

  OptimizerParametersD<MyCase::value_t> optiParameters;
  {
    using namespace parameters;
    optiParameters.set<optimizer::EPS>(1e-10);
    optiParameters.set<optimizer::MAX_OPTIMIZATION_STEPS>(50);
    optiParameters.set<linesearch::MAX_STEP_ATTEMPTS>(100);
    optiParameters.set<optimizer::VERBOSE>(true);
  }

  /// === Step 3: Create Mesh ===
  Mesh mesh    = createMesh(myCaseParameters);
  Mesh adfMesh = createMesh(myADfCaseParameters);
  Mesh refMesh = createMesh(myReferenceParameters);

  // ==== BELOW HERE OPTIMIZATION SPECIFIC ====
  /// @li Step A: Create Cases and prepare
  // Prepare controlled case
  MyCase myCase(myCaseParameters, mesh);
  prepareGeometry(myCase);

  // Create ADf-typed case for gradient computation
  MyADfCase aDfCase(myADfCaseParameters, adfMesh);
  prepareGeometry(aDfCase);

  // Compute solution for the objective functional
  MyCase referenceCase(myReferenceParameters, refMesh);
  prepareGeometry(referenceCase);
  prepareLattice(referenceCase);
  setInitialValues(referenceCase);

  // Forward problem at the target modulus only
  if (myCaseParameters.get<parameters::FORWARD_ONLY>()) {
    simulate(referenceCase);
    writeErrors(referenceCase);
    return 0;
  }

  /// @li Step B: Run the numerical reference case
  // For comparison with a numerical reference solution, run a forward simulation
  // with the reference case at the target Young's modulus
  if (not myCaseParameters.get<parameters::ANALYTICAL_REFERENCE>()) {
    simulate(referenceCase);
  }

  /// @li Step C: Create OptiCase and set cases
  MyOptiCase optiCase;
  optiCase.setCase<Controlled >(myCase);
  optiCase.setCase<Derivatives>(aDfCase);
  optiCase.setCase<Reference  >(referenceCase);

  /// @li Step D: Set initial control
  setInitialControls(optiCase);

  /// @li Step E: Define objective routine
  optiCase.setObjective([](auto& optiCase, auto name) {
    return computeObjective(optiCase, name);
  });

  using T = MyCase::value_t;

  /// @li Step F: Create an Optimizer
  OptimizerLBFGS<T,std::vector<T>> optimizer(
    optiCase.getController().size(),
    optiParameters,
    std::vector<olb::opti::OptimizerLogType>(
      {OptimizerLogType::control, OptimizerLogType::objective,
       OptimizerLogType::norm_derivative, OptimizerLogType::rel_error}
    )
  );

  // Set a reference control to evaluate the relative error of the control update
  clout << "Initial Control: " << optiCase.getController().get()[0] << std::endl;
  optimizer.setReferenceControl(std::vector<T>{util::sqrt(targetE)});

  /// @li Step G: Optimize
  optimizer.optimize(optiCase);
}
