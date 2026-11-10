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

namespace olb::parameters {
  struct ACTIVE : public fields::SCALAR<false> { };
  struct FDQ_SHIFT : public fields::SCALAR<1e-5> { };
  struct SAMPLE_POS : public descriptors::FIELD_BASE<0,1> { };
  struct PLOT_TYPE : public fields::STRING<"diagonal"> { };
} // namespace olb::parameters

std::vector<Vector<MyCase::value_t, MyCase::descriptor_t::d>> setSamplingPositions(MyCase& Case)
{
  using T = MyCase::value_t;
  auto& parameters = Case.getParameters();
  auto& converter = Case.getLattice(NavierCauchy{}).getUnitConverter();
  T const dx = converter.getPhysDeltaX();
  size_t const N = parameters.get<parameters::RESOLUTION>();
  size_t stepIdx{};

  // If N >= 20, sample every 20th node, otherwise sample every node
  if (N >= 20) {
    stepIdx = N/20;
  } else {
    stepIdx = 1;
  }
  // Sample every at every node from x=stepIdx*dx to x=1-stepIdx*dx (domain_extent - dx)
  size_t const numSamples = (N / stepIdx) - 1;

  std::vector<Vector<MyCase::value_t, MyCase::descriptor_t::d>> samplePhysPos(numSamples);
  for (size_t i = 0; i < numSamples; i++) {
    auto& pos = samplePhysPos[i];
    //Sample along a straight line at y=SAMPLE_POSE[1] or along the diagonal from (0,0) to (1,1)
    if ( parameters["gradient_comparison"].get<parameters::PLOT_TYPE>() == "straight") {
      pos[0] = (i+1) * stepIdx * dx;
      pos[1] = parameters["gradient_comparison"].get<parameters::SAMPLE_POS>()[1];
    } else if (parameters["gradient_comparison"].get<parameters::PLOT_TYPE>() == "diagonal") {
      pos[0] = (i+1) * stepIdx * dx;
      pos[1] = (i+1) * stepIdx * dx;
    }
  }

  return samplePhysPos;
}

//
std::vector<Vector<MyCase::value_t, MyCase::descriptor_t::d>> sampleDerivativesAtPoints(
  MyCase& adjointCase,
  const std::vector<Vector<MyCase::value_t, MyCase::descriptor_t::d>>& samplePos
) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;
  using ControlledField = PeriodicPlate::controls_t;
  constexpr unsigned dim = MyCase::d;

  auto& sLattice = adjointCase.getLattice(NavierCauchy{});
  SuperLatticeField2D<T,DESCRIPTOR,opti::SENSITIVITY<ControlledField>> gradF(sLattice);
  AnalyticalFfromSuperF2D<T,T> grad(gradF, true);

  std::vector<Vector<T, dim>> gradients(samplePos.size());
  for (size_t i = 0; i < samplePos.size(); i++) {
    T input[]{samplePos[i][0], samplePos[i][1]};
    T output[dim];
    grad(output, input);
    gradients[i] = output;
  }

  return gradients;
}

// Introduce a shift in the control for a single cell position
auto getShiftedControl(MyOptiCase& optiCase, MyCase::value_t shift, unsigned dim, Vector<MyCase::value_t, MyCase::descriptor_t::d> const& samplePhysPos) {
  using ControlledField = PeriodicPlate::controls_t;

  auto& controlledCase = optiCase.getCase(Controlled{});
  auto& lattice = controlledCase.getLattice(NavierCauchy{});
  auto& geometry = controlledCase.getGeometry();
  const auto& converter = lattice.getUnitConverter();

  auto& controller = optiCase.getController();

  setInitialControls(optiCase);
  auto& cuboidDecomposition = controlledCase.getMesh().getCuboidDecomposition();
  auto const minR = cuboidDecomposition.getMinPhysR();
  auto const maxR = cuboidDecomposition.getMaxPhysR() - cuboidDecomposition.getDeltaR();
  auto const tol = 0.5 * cuboidDecomposition.getDeltaR();
  if (samplePhysPos[0] < minR[0] - tol || samplePhysPos[0] > maxR[0] + tol ||
      samplePhysPos[1] < minR[1] - tol || samplePhysPos[1] > maxR[1] + tol) {
    std::ostringstream ss;
    ss << "Sample position is outside the domain: " << samplePhysPos << std::endl;
    throw std::runtime_error(ss.str());
  }
  auto samplePos = cuboidDecomposition.getLatticeR(samplePhysPos);
  if (!samplePos) {
    std::ostringstream ss;
    ss << "Invalid sample position: " << samplePhysPos << std::endl;
    throw std::runtime_error(ss.str());
  }

  if (controlledCase.getMesh().getLoadBalancer().isLocal((*samplePos)[0])) {
    auto cell = lattice.get(*samplePos);
    auto control = cell.template getFieldPointer<ControlledField>();
    control[dim] += shift * 1. / converter.getConversionFactorSolidForce();
  }
  controller.template set<ControlledField>(geometry, 1, lattice);
  return controller.get();
}

std::vector<Vector<MyCase::value_t, MyCase::descriptor_t::d>> computeDerivativesWithFDQAtPoints(
  MyOptiCase& optiCase,
  MyCase::value_t objective,
  std::vector<Vector<MyCase::value_t, MyCase::descriptor_t::d>> const& samplePos
) {
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;

  T const shift = optiCase.getCase(Controlled{}).getParameters()["gradient_comparison"].get<parameters::FDQ_SHIFT>();

  std::vector<Vector<T, DESCRIPTOR::d>> gradients(samplePos.size());
  for (size_t i = 0; i < samplePos.size(); ++i) {
    for (std::size_t dim = 0; dim < DESCRIPTOR::d; ++dim) {
      auto shiftedControl = getShiftedControl(optiCase, shift, dim, samplePos[i]);
      auto shiftedObjective = optiCase.computeObjective(shiftedControl, 0);
      gradients[i][dim] = (shiftedObjective - objective) /
                          (shift);
    }
  }

  return gradients;
}

void runGradientComparison(MyOptiCase& optiCase)
{
  OstreamManager clout("runGradientComparison");
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t;

  auto const samplingPositions = setSamplingPositions(optiCase.getCase(Controlled{}));
  clout << "Evaluating gradients at " << samplingPositions.size() << " sample positions." << std::endl;

  auto& controls = optiCase.getController().get();

  // Compute the Objective functional
  auto objective = optiCase.computeObjective(controls, 0);
  // Evaluate discrete-adjoint-based gradients
  optiCase.computeDerivatives(controls, 0);
  auto const adjointGrads = sampleDerivativesAtPoints(optiCase.getCase(Adjoint{}), samplingPositions);

  // Compute FDQ-based gradients at the same sample positions
  auto const fdqGrads = computeDerivativesWithFDQAtPoints(optiCase, objective, samplingPositions);

  // Compute absolute and relative errors between adjoint and FDQ gradients
  std::vector<Vector<T,DESCRIPTOR::d>> relErrors(adjointGrads.size());
  std::vector<Vector<T,DESCRIPTOR::d>> absErrors(adjointGrads.size());
  for(std::size_t i = 0; i < adjointGrads.size(); i++) {
    for(std::size_t d = 0; d < DESCRIPTOR::d; d++){
      absErrors[i][d] = util::abs(adjointGrads[i][d] - fdqGrads[i][d]);
      relErrors[i][d] = absErrors[i][d] / util::abs(fdqGrads[i][d]);
    }
  }

  Vector<T, DESCRIPTOR::d> totalL2error{};
  for(std::size_t d=0; d < DESCRIPTOR::d; d++) {
    for(std::size_t i=0; i < adjointGrads.size(); i++) {
      if(util::abs(adjointGrads[i][d]) > 1e-10)
        totalL2error[d] += util::pow(adjointGrads[i][d] - fdqGrads[i][d], 2);
    }
    totalL2error[d] = util::sqrt(totalL2error[d]);
  }

  clout << "total L2 errors (X,Y): " << totalL2error[0] << " " << totalL2error[1] << std::endl;


  CSV<MyCase::value_t> csvWriter("gradient_comparison",
    ';',
    {"Index", "SamplePosX", "SamplePosY", "AdjointGradX",
      "AdjointGradY", "FDQGradX", "FDQGradY", "ABSERRORX",
      "ABSERRORY", "RELERRORX", "RELERRORY"},
    ".csv");
  for (size_t i = 0; i < samplingPositions.size(); ++i) {
    csvWriter.writeDataFile(i, {samplingPositions[i][0], samplingPositions[i][1],
                                adjointGrads[i][0], adjointGrads[i][1],
                                fdqGrads[i][0], fdqGrads[i][1],
                                absErrors[i][0], absErrors[i][1],
                                relErrors[i][0], relErrors[i][1]
                              }
                            );
  }
}
