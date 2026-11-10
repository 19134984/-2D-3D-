// Disable SIMD for ADf
/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Johannes Grafen
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

#include "../../examples/optimization/testFlowOpti3d/case.h"

#include "gtest/gtest.h"
#include "../tests/unit/helper/framework.h"

/// ensure olbInit is called for MPI
::testing::Environment* const foo_env = ::testing::AddGlobalTestEnvironment(new OlbInitTestEnvironment);

namespace olb::parameters {
struct OBJECTIVE_L2_NORM : public descriptors::FIELD_BASE<1> {};
struct OBJECTIVE_DERIVATIVE_L2_NORM : public descriptors::FIELD_BASE<1> {};
struct NUM_ITERATIONS : public descriptors::FIELD_BASE<1> {};
}; // namespace olb::parameters

namespace OptimizerOptions {
constexpr int         maxIt {10};
constexpr double      lambda {1.0};
constexpr int         l {20};
constexpr double      startCoefH {1e-4};
constexpr double      eps {6e-5};
constexpr int         maxStepAttempts {20};
constexpr std::string stepCondition {"Wolfe"};

constexpr double      controlEps {std::numeric_limits<double>::epsilon()};
constexpr bool        vectorBounds {false};
constexpr bool        verboseOn {true};
constexpr std::string fname {};
constexpr std::string logFileName {};
constexpr bool        withUpperBound {false};
constexpr double      upperBound {};
constexpr bool        withLowerBound {false};
constexpr double      lowerBound {};
constexpr bool        failOnMaxIter {false};
}; // namespace OptimizerOptions

constexpr double PRECISION = 1e-13;

void runAdjointOptiCase(MyCase::ParametersD& myCaseParameters)
{
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
  OptimizerLBFGS<MyOptiCase::value_t, std::vector<MyOptiCase::value_t>> optimizer(
      optiCase.getController().size(), OptimizerOptions::eps, OptimizerOptions::maxIt, OptimizerOptions::lambda,
      OptimizerOptions::maxStepAttempts, OptimizerOptions::stepCondition, OptimizerOptions::l,
      OptimizerOptions::startCoefH, OptimizerOptions::verboseOn, OptimizerOptions::fname, OptimizerOptions::logFileName,
      OptimizerOptions::withUpperBound, OptimizerOptions::upperBound, OptimizerOptions::withLowerBound,
      OptimizerOptions::lowerBound, OptimizerOptions::vectorBounds, OptimizerOptions::controlEps,
      OptimizerOptions::failOnMaxIter);
  /// @li Step H: Optimize
  optimizer.optimize(optiCase);

  myCaseParameters.template set<parameters::OBJECTIVE_L2_NORM>(optimizer.getObjective());
  myCaseParameters.template set<parameters::OBJECTIVE_DERIVATIVE_L2_NORM>(util::euklidN(optimizer.getDerivative()));
  myCaseParameters.template set<parameters::NUM_ITERATIONS>(optimizer.getIteration());
}

TEST(TestFlowOptiAdjointContinuous, TestVelocityObjective)
{
  MyCase::ParametersD myCaseParameters;
  {
    using namespace parameters;
    myCaseParameters.set<DOMAIN_EXTENT>({1.0, 1.0, 1.0});
    myCaseParameters.set<PHYS_CHAR_VELOCITY>(1.0);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(0.1);
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1.0);
    myCaseParameters.set<PHYS_MAX_T>(6.0);
    myCaseParameters.set<RESOLUTION>(11);
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.8);
    myCaseParameters.set<ADJOINT_TYPE>(AdjointType::continuous);
  }

  runAdjointOptiCase(myCaseParameters);

  EXPECT_NEAR(myCaseParameters.template get<parameters::OBJECTIVE_L2_NORM>(), 0.001158781596648804, PRECISION);
  EXPECT_NEAR(myCaseParameters.template get<parameters::OBJECTIVE_DERIVATIVE_L2_NORM>(), 0.05943654746696577, PRECISION);
  ASSERT_EQ(myCaseParameters.template get<parameters::NUM_ITERATIONS>(), 10);
}

TEST(TestFlowOptiAdjointDiscrete, TestVelocityObjective)
{
  MyCase::ParametersD myCaseParameters;
  {
    using namespace parameters;
    myCaseParameters.set<DOMAIN_EXTENT>({1.0, 1.0, 1.0});
    myCaseParameters.set<PHYS_CHAR_VELOCITY>(1.0);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(0.1);
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1.0);
    myCaseParameters.set<PHYS_MAX_T>(6.0);
    myCaseParameters.set<RESOLUTION>(11);
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(0.8);
    myCaseParameters.set<ADJOINT_TYPE>(AdjointType::discrete);
  }

  runAdjointOptiCase(myCaseParameters);

  EXPECT_NEAR(myCaseParameters.template get<parameters::OBJECTIVE_L2_NORM>(), 0.001161108335396035, PRECISION);
  EXPECT_NEAR(myCaseParameters.template get<parameters::OBJECTIVE_DERIVATIVE_L2_NORM>(), 0.06003185488536415, PRECISION);
  ASSERT_EQ(myCaseParameters.template get<parameters::NUM_ITERATIONS>(), 10);
}
