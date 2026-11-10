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

// Disable SIMD for ADf
#undef PLATFORM_CPU_SIMD

#include "../../examples/optimization/testFlowOptiAD3d/case.h"

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
constexpr double      eps {1e-5};
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

constexpr double PRECISION = 1e-14;

void runADfOptiCase(MyCase::ParametersD& myCaseParametersD, MyADfCase::ParametersD& myADfCaseParametersD)
{
  // ==== BELOW HERE OPTIMIZATION SPECIFIC ====
  /// @li Step B: Create Cases and prepare
  // Prepare controlled case
  Mesh   mesh = createMesh(myCaseParametersD);
  MyCase myCase(myCaseParametersD, mesh);
  prepareGeometry(myCase);

  // Create ADf-typed case for gradient computation
  Mesh      aDfMesh = createMesh(myADfCaseParametersD);
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

  OptimizerLBFGS<MyCase::value_t, std::vector<MyCase::value_t>> optimizer(
      optiCase.getController().size(), OptimizerOptions::eps, OptimizerOptions::maxIt, OptimizerOptions::lambda,
      OptimizerOptions::maxStepAttempts, OptimizerOptions::stepCondition, OptimizerOptions::l,
      OptimizerOptions::startCoefH, OptimizerOptions::verboseOn, OptimizerOptions::fname, OptimizerOptions::logFileName,
      OptimizerOptions::withUpperBound, OptimizerOptions::upperBound, OptimizerOptions::withLowerBound,
      OptimizerOptions::lowerBound, OptimizerOptions::vectorBounds, OptimizerOptions::controlEps,
      OptimizerOptions::failOnMaxIter);
  /// @li Step H: Optimize
  optimizer.optimize(optiCase);

  myCaseParametersD.template set<parameters::OBJECTIVE_L2_NORM>(optimizer.getObjective());
  myCaseParametersD.template set<parameters::OBJECTIVE_DERIVATIVE_L2_NORM>(util::euklidN(optimizer.getDerivative()));
  myCaseParametersD.template set<parameters::NUM_ITERATIONS>(optimizer.getIteration());
}

////////////////////////////////////////////////////////////////////////////////
TEST(TestFlowOptiAD, testVelocityObjective)
{
  MyCase::ParametersD myCaseParametersD;
  {
    using namespace parameters;
    myCaseParametersD.template set<PHYS_CHAR_VELOCITY>(1.0);
    myCaseParametersD.template set<PHYS_CHAR_LENGTH>(1.0);
    myCaseParametersD.template set<PHYS_CHAR_VISCOSITY>(0.1);
    myCaseParametersD.template set<PHYS_CHAR_DENSITY>(1.0);
    myCaseParametersD.template set<PHYS_MAX_T>(6.0);
    myCaseParametersD.template set<RESOLUTION>(11);
    myCaseParametersD.template set<LATTICE_RELAXATION_TIME>(0.8);
    myCaseParametersD.template set<OBJECTIVE_TYPE>(ObjectiveType::velocity);
  }
  auto myADfCaseParametersD = myCaseParametersD.copyAs<MyADfCase::value_t>();

  runADfOptiCase(myCaseParametersD, myADfCaseParametersD);

  EXPECT_NEAR(myCaseParametersD.template get<parameters::OBJECTIVE_L2_NORM>(), 2.67970616393e-14, PRECISION);
  EXPECT_NEAR(myCaseParametersD.template get<parameters::OBJECTIVE_DERIVATIVE_L2_NORM>(), 5.88862297352e-07, PRECISION);
  ASSERT_EQ(myCaseParametersD.template get<parameters::NUM_ITERATIONS>(), 10);
}

////////////////////////////////////////////////////////////////////////////////
TEST(TestFlowOptiAD, testDissipationObjective)
{
  MyCase::ParametersD myCaseParametersD;
  {
    using namespace parameters;
    myCaseParametersD.template set<PHYS_CHAR_VELOCITY>(1.0);
    myCaseParametersD.template set<PHYS_CHAR_LENGTH>(1.0);
    myCaseParametersD.template set<PHYS_CHAR_VISCOSITY>(0.1);
    myCaseParametersD.template set<PHYS_CHAR_DENSITY>(1.0);
    myCaseParametersD.template set<PHYS_MAX_T>(6.0);
    myCaseParametersD.template set<RESOLUTION>(11);
    myCaseParametersD.template set<LATTICE_RELAXATION_TIME>(0.8);
    myCaseParametersD.template set<OBJECTIVE_TYPE>(ObjectiveType::dissipation);
  }
  auto myADfCaseParametersD = myCaseParametersD.copyAs<MyADfCase::value_t>();

  runADfOptiCase(myCaseParametersD, myADfCaseParametersD);

  EXPECT_NEAR(myCaseParametersD.template get<parameters::OBJECTIVE_L2_NORM>(), 3.18124246448e-08, PRECISION);
  EXPECT_NEAR(myCaseParametersD.template get<parameters::OBJECTIVE_DERIVATIVE_L2_NORM>(), 0.00237311884109, PRECISION);
  ASSERT_EQ(myCaseParametersD.template get<parameters::NUM_ITERATIONS>(), 10);
}
