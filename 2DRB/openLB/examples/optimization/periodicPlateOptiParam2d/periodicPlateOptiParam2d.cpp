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

/// Setup and run a simulation
int main(int argc, char* argv[])
{
  OstreamManager clout(std::cout, "main");
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<SOLID_CASE>("continua");
    myCaseParameters.set<RESOLUTION>(40);
    myCaseParameters.set<YOUNGS_MODULUS>(0.11);
    myCaseParameters.set<PHYS_DENSITY>(1);
    myCaseParameters.set<POISSON_RATIO>(0.8);
    myCaseParameters.set<DAMPING_FACTOR>(1.);
    myCaseParameters.set<PHYS_CHAR_LENGTH>(1.0);
    myCaseParameters.set<DOMAIN_EXTENT>({1.0, 1.0});
    myCaseParameters.set<ORIGIN>({0.0, 0.0});
    myCaseParameters.set<PHYS_MAX_T>(100);
    myCaseParameters.set<LATTICE_LOG_ITER_T>(1000);
    myCaseParameters.set<LATTICE_VTK_ITER_T>(1000);
    myCaseParameters.set<CONV_ITER>(10);
    myCaseParameters.set<CONV_VALUE>(1);
    myCaseParameters.set<CONVERGENCE_PRECISION>(1e-7);
    myCaseParameters.set<CONVERGENCE_PRECISION_ADJOINT>(1e-7);
    myCaseParameters.set<VTK_ENABLED>(true);
    myCaseParameters.set<ANALYTICAL_REFERENCE>(true);
    myCaseParameters.set<CONVERGED_IT>(0);
    myCaseParameters.set<INITIAL_CONTROL>({-1.0, -1.0});
    myCaseParameters.set<SET_ADJOINT_MAX_IT>(false);
  }
  //Compute additional derived parameters
  myCaseParameters.set<parameters::PHYS_DELTA_X>(myCaseParameters.get<parameters::DOMAIN_EXTENT>()[0] /
                                                 myCaseParameters.get<parameters::RESOLUTION>());
  myCaseParameters.set<parameters::PHYS_DELTA_T>(myCaseParameters.get<parameters::PHYS_DELTA_X>() *
                                                 myCaseParameters.get<parameters::PHYS_DELTA_X>());
  myCaseParameters.set<parameters::PHYS_CHAR_DISPLACEMENT>(myCaseParameters.get<parameters::PHYS_DELTA_X>());

  /// Step 2.b: Set optimizer parameters
  OptimizerParametersD<MyCase::value_t> optiParameters;
  {
    using namespace parameters;
    optiParameters.set<optimizer::EPS>(1e-10);
    optiParameters.set<optimizer::MAX_OPTIMIZATION_STEPS>(50);
    optiParameters.set<linesearch::MAX_STEP_ATTEMPTS>(100);
    optiParameters.set<optimizer::VERBOSE>(true);
    optiParameters.set<LBFGS::HISTORY_SIZE>(20);
    optiParameters.set<linesearch::STEP_CONDITION>("StrongWolfe");
  }

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);
  prepareGeometry(myCase);
  prepareLattice(myCase);

  MyCase adjointCase(myCaseParameters, mesh);
  prepareGeometry(adjointCase);

  MyCase referenceCase(myCaseParameters, mesh);
  prepareGeometry(referenceCase);
  prepareLattice(referenceCase);
  setInitialValues(referenceCase);

  // For comparison with numerical reference solution, run a forward simulation with the reference case with the analytical force field
  if(not myCaseParameters.get<parameters::ANALYTICAL_REFERENCE>()) {
    clout << "Running reference simulation ..." << std::endl;
    simulate(referenceCase, myCaseParameters.get<parameters::CONVERGENCE_PRECISION>());
    clout << "Finished reference simulation." << std::endl;
  }
  MyOptiCase optiCase;
  optiCase.setCase<Controlled>(myCase);
  optiCase.setCase<Adjoint>(adjointCase);
  optiCase.setCase<Reference>(referenceCase);

  PeriodicPlate optimalitySystem;
  optiCase.set<Optimality>(optimalitySystem);
  setInitialControls(optiCase);

  optiCase.setObjective(computeObjective);

  optiCase.setDerivative(computeDerivative);

  /// @li Step G: Create an Optimizer
  using T = MyCase::value_t;
  OptimizerLBFGS<T, std::vector<T>> optimizer(
      optiCase.getController().size(),
      optiParameters,
      std::vector<olb::opti::OptimizerLogType>(
        {control, objective, norm_derivative, rel_error}
      )
  );

  // The reference control alpha* = (1 1) reproduces the analytical force field
  optimizer.setReferenceControl(std::vector<T>(optiCase.getController().size(), T(1)));

  /// @li Step H: Optimize
  optimizer.optimize(optiCase);
}
