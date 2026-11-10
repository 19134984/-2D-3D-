// Disable SIMD for ADf
/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Shota Ito, Johannes Grafen
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
  struct SAMPLE_POS : public descriptors::FIELD_BASE<0,1> { };
  struct CONTINUOUS_ADJOINT_GRADIENTS : public descriptors::FIELD_BASE<0,1> { };
  struct DISCRETE_ADJOINT_GRADIENTS : public descriptors::FIELD_BASE<0,1> { };
  struct FDQ_GRADIENTS : public descriptors::FIELD_BASE<0,1> { };
} // namespace olb::parameters

constexpr double PRECISION_ADJOINTS = 1e-4;
constexpr double PRECISION_ADJOINT_FDQ = 1e-1;

// Introduce a shift in the control for a single cell position
auto getShiftedControl(MyOptiCase& optiCase, MyCase::value_t shift, unsigned dim) {
  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& lattice = controlledCase.getLattice(NavierStokes{});
  auto& geometry = controlledCase.getGeometry();
  auto& control = optiCase.getController();
  auto& parameters = controlledCase.getParameters();

  setInitialControl(optiCase);
  Vector<MyCase::value_t,3> samplePhysPos(parameters.get<parameters::SAMPLE_POS>());

  if ( auto samplePos = controlledCase.getMesh().getCuboidDecomposition().getLatticeR(samplePhysPos)) {
    if (controlledCase.getMesh().getLoadBalancer().isLocal((*samplePos)[0])) {
      auto cell = lattice.get(*samplePos);
      auto force = cell.template getFieldPointer<ControlledField>();
      force[dim] += shift;
    }
    control.template set<ControlledField>(geometry, 1, lattice);
  } else {
    throw std::runtime_error("Invalid sample position.");
  }
  return control.get();
}

void sampleDerivativesAndPrint(MyCase& adjointCase) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;

  auto& sLattice = adjointCase.getLattice(NavierStokes{});
  auto& parameters = adjointCase.getParameters();
  SuperLatticeField3D<T,DESCRIPTOR,opti::SENSITIVITY<ControlledField>> gradF(sLattice);
  AnalyticalFfromSuperF3D<T,T> grad(gradF);

  T output[3], input[3];
  input[0] = parameters.get<parameters::SAMPLE_POS>()[0];
  input[1] = parameters.get<parameters::SAMPLE_POS>()[1];
  input[2] = parameters.get<parameters::SAMPLE_POS>()[2];
  grad(output, input);
  const std::string type = (parameters.get<parameters::ADJOINT_TYPE>() == AdjointType::continuous ) ? "Conti. " : "Dis.   ";
  std::cout << type + "Gradient at " << parameters.get<parameters::SAMPLE_POS>()
                                    << ": (" << output[0]
                                    << ", "  << output[1]
                                    << ", "  << output[2]
                                    << ")" << std::endl;
  if(parameters.get<parameters::ADJOINT_TYPE>() == AdjointType::continuous) {
    parameters.set<parameters::CONTINUOUS_ADJOINT_GRADIENTS>({output[0], output[1], output[2]});
  } else {
    parameters.set<parameters::DISCRETE_ADJOINT_GRADIENTS>({output[0], output[1], output[2]});
  }
}

std::vector<MyCase::value_t> computeDerivativeWithSampling(MyOptiCase& optiCase) {
  auto& adjointCase = optiCase.getCase(Adjoint{});
  // Reset prior simulation lattice
  adjointCase.resetLattices();

  // Prepare new lattices
  prepareAdjointLattice(optiCase);
  setAdjointInitialValues(optiCase);

  // Execute simulation
  simulate(adjointCase);

  // Compute Derivative value from simulation results
  auto result = derivativeF(optiCase);
  // Evaluate adjoint-based gradients
  sampleDerivativesAndPrint(adjointCase);

  return result;
}

void computeDerivativesWithFDQ(MyOptiCase& optiCase, MyCase::value_t objective, MyCase::value_t shift = 1e-6) {
  Vector<MyCase::value_t,3> gradients({0.0, 0.0, 0.0});
  for (std::size_t dim = 0; dim < 3; ++dim) {
    auto shiftedControl = getShiftedControl(optiCase, shift, dim);
    auto shiftedObjective = optiCase.computeObjective(shiftedControl, 0);
    gradients[dim] = (shiftedObjective - objective) /
                     (shift * optiCase.getCase(Controlled{}).getLattice(NavierStokes{}).getUnitConverter().getConversionFactorForce());
  }
  std::cout << "FDQ.   Gradient at " << optiCase.getCase(Controlled{}).getParameters().get<parameters::SAMPLE_POS>()
                                      << ": " << gradients << std::endl;
  optiCase.getCase(Controlled{}).getParameters().set<parameters::FDQ_GRADIENTS>(gradients);
}

void runTestFlowGradientComparison(MyCase::ParametersD& myCaseParameters)
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
  optiCase.setDerivative(computeDerivativeWithSampling);

  auto& controls = optiCase.getController().get();
  std::vector<MyOptiCase::value_t> derivatives;
  derivatives.resize(controls.size(), 0);

  // Compute the Objective functional
  auto objective = optiCase.computeObjective(controls, 0);
  std::cout << "Objective: " << objective << std::endl;

  // Evaluate continuous-adjoint-based gradients
  derivatives = optiCase.computeDerivatives(controls, 0);

  // Evaluate discrete-adjoint-based gradients
  myCaseParameters.set<parameters::ADJOINT_TYPE>(AdjointType::discrete);
  derivatives = optiCase.computeDerivatives(controls, 0);

  // Evaluate FDQ-based gradients
  computeDerivativesWithFDQ(optiCase, objective);

  std::cout << "Summary:" << std::endl;
  std::cout << "Conti. Gradient at " << optiCase.getCase(Controlled{}).getParameters().get<parameters::SAMPLE_POS>()
                                     << ": " << optiCase.getCase(Controlled{}).getParameters().get<parameters::CONTINUOUS_ADJOINT_GRADIENTS>()
                                     << std::endl;
  std::cout << "Dis.   Gradient at " << optiCase.getCase(Controlled{}).getParameters().get<parameters::SAMPLE_POS>()
                                     << ": " << optiCase.getCase(Controlled{}).getParameters().get<parameters::DISCRETE_ADJOINT_GRADIENTS>()
                                     << std::endl;
  std::cout << "FDQ.   Gradient at " << optiCase.getCase(Controlled{}).getParameters().get<parameters::SAMPLE_POS>()
                                     << ": " << optiCase.getCase(Controlled{}).getParameters().get<parameters::FDQ_GRADIENTS>()
                                     << std::endl;
}

TEST(TestFlowOptiSensitivity, TestSensitivityComparison)
{
  MyCase::ParametersD myCaseParameters;
  {
    using namespace parameters;
    myCaseParameters.set<DOMAIN_EXTENT          >(        {1.0, 1.0, 1.0});
    myCaseParameters.set<PHYS_CHAR_VELOCITY     >(                    1.0);
    myCaseParameters.set<PHYS_CHAR_VISCOSITY    >(                    0.1);
    myCaseParameters.set<PHYS_CHAR_DENSITY      >(                    1.0);
    myCaseParameters.set<PHYS_MAX_T             >(                    6.0);
    myCaseParameters.set<RESOLUTION             >(                     21);
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(                    0.8);
    myCaseParameters.set<ADJOINT_TYPE           >(AdjointType::continuous);
    myCaseParameters.set<SAMPLE_POS             >(        {0.5, 0.5, 0.5});
  }

  runTestFlowGradientComparison(myCaseParameters);

  for (int d = 0; d < MyCase::descriptor_t::d; ++d) {
    EXPECT_NEAR(myCaseParameters.template get<parameters::CONTINUOUS_ADJOINT_GRADIENTS>()[d],
                myCaseParameters.template get<parameters::DISCRETE_ADJOINT_GRADIENTS>()[d],
                PRECISION_ADJOINTS);
    EXPECT_NEAR(myCaseParameters.template get<parameters::CONTINUOUS_ADJOINT_GRADIENTS>()[d],
                myCaseParameters.template get<parameters::FDQ_GRADIENTS>()[d],
                PRECISION_ADJOINT_FDQ);
    EXPECT_NEAR(myCaseParameters.template get<parameters::DISCRETE_ADJOINT_GRADIENTS>()[d],
                myCaseParameters.template get<parameters::FDQ_GRADIENTS>()[d],
                PRECISION_ADJOINT_FDQ);
  }
}
