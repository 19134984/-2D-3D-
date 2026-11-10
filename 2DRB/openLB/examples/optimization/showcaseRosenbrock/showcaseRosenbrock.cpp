/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022-2025 Julius Jeßberger, Shota Ito
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

/** \file
 * This example illustrates the application of numerical optimization
 * algorithms on the (simple) example unconstrained minimization of the
 * Rosenbrock function, cf.
 * https://en.wikipedia.org/wiki/Rosenbrock_function. For a more complex
 * application with flow simulation, we refer to the example
 * domainIdenfitication3d.
 * The script is structured as a tutorial that can be read sequentially.
 */


#include <olb.h>

using namespace olb;
using namespace olb::opti;

using T = FLOATING_POINT_TYPE;
using C = std::vector<T>;

/** ---------------------------------------------------------------------------
 * Part 0: Preliminary remarks
 * ----------------------------------------------------------------------------
 * The n-dimensional Rosenbrock function has a global minimum at (1,1,1,...).
 * For some n, there is additionally another local minimum. The optimization
 * algorithms have to identify one of these minima.
 * Here, we use only gradient-based schemes which are fast but do not consider
 * the global shape of a function as they try to decrease the function locally.
 */

/** ---------------------------------------------------------------------------
 * Part 1a: Definition of an objective function
 * ----------------------------------------------------------------------------
 * We define the function which is to be minimized. We can use any c++-function
 * type to do so that satisfies the signature T (const std::vector<T>&) for
 * some type T.
 */

const double b (100);
const double eps (1.e-7);
constexpr int N (4);

// Implementation of the Rosenbrock function
template<typename T>
std::function<T (const std::vector<T>&)> rosenbrock
 = [](const std::vector<T>& control) -> T {
    T result(0);
    for (std::size_t i = 0; i < control.size(); i += 2) {
      T c1 = (control[i + 1] - control[i] * control[i]);
      T c2 = 1.0 - control[i];
      result += b * c1 * c1 + c2 * c2;
    }
    return result;
  };

/** ---------------------------------------------------------------------------
 * Part 1b: Derivative of the objective function
 * ----------------------------------------------------------------------------
 * Any gradient-based optimization scheme requires computation of the gradient
 * of the objective function. Several ways to compute this are implemented in
 * the OptiCase-classes.
 */

// For simple examples, we can just write down the analytical formula for the
// derivative
template<typename T>
std::vector<T> rosenbrockDerivAnalytical(const std::vector<T>& control)
{
  std::vector<T> derivatives(control.size());
  for (std::size_t i = 0; i < control.size(); i += 2) {
    T c1 = (control[i + 1] - control[i] * control[i]);
    T c2 = 1.0 - control[i];
    derivatives[i + 1] =  2.0 * b * c1;
    derivatives[i]     = -4.0 * b * c1 * control[i] - 2.0 * c2;
  }
  return derivatives;
}

/* We can then instantiate the OptiCaseAnalytical, which takes the objective
 * function and its derivative as arguments
 */
OptiCaseAnalytical<T,derivatives::Manual> ocAnalytical;

/* For complex examples, it is no longer feasible to write down the derivative
 * by hand. We can e.g. use finite difference quotients (forward or central),
 * where we have to tell the step width at construction.
 * For forward DQ, half of machine precision is a suitable (optimal) choice to
 * reduce cancellation errors, for central DQ, one third of machine precision
 * is good.
 */
OptiCaseAnalytical<T,derivatives::FDQ> ocForwardDQ;
OptiCaseAnalytical<T,derivatives::CDQ> ocCentralDQ;

/* For more accuracy, we can also use automatic differentiation.
 * To do so, we have to pass an ADf-typed version of the objective function.
 * Therefore, the objective function has to be templatized w.r.t. the
 * arithmetic data type.
 */
using U = util::ADf<T,N>;  // pass size of control vector as second argument
OptiCaseAnalytical<U,derivatives::ADf> ocAD;

// All of these opti-cases can later be used for optimization.


/** ---------------------------------------------------------------------------
 * Part 2: Optimization
 * ----------------------------------------------------------------------------
 * We select and call an optimization algorithm
 */

void optimizationRosenbrock(OptimizerParametersD<T>& optiParameters){

  /* Select optimization algorithm. Steepest descent is the slowest but most
   * stable, LBFGS and Barzilai-Borwein are considerably faster
   */
  OptimizerSteepestDescent<T,C> optimizerSD(N, optiParameters);
  OptimizerLBFGS<T,C> optimizerLBFGS(N, optiParameters);
  OptimizerBarzilaiBorwein<T,C> optimizerBB(N, optiParameters);

  // The last argument tells which data shall be written and plotted by gnuplot
  optimizerLBFGS.addLogging({OptimizerLogType::control, OptimizerLogType::derivative});

  // Set objective functional and derivative
  ocAnalytical.setObjective(rosenbrock<T>);
  ocForwardDQ.setObjective(rosenbrock<T>);
  ocCentralDQ.setObjective(rosenbrock<T>);
  ocAD.setObjective(rosenbrock<T>, rosenbrock<U>);

  ocAnalytical.setDerivative(rosenbrockDerivAnalytical<T>);

  // Define start value for minimum search
  ocAnalytical.getController().set({3.0, 0.0, 2.0, 0.0});
  ocForwardDQ.getController().set({3.0, 0.0, 2.0, 0.0});
  ocCentralDQ.getController().set({3.0, 0.0, 2.0, 0.0});
  ocAD.getController().set({3.0, 0.0, 2.0, 0.0});

  // Execute optimization
  optimizerLBFGS.optimize(ocAnalytical);
  //optimizerLBFGS.optimize(ocForwardDQ);
  //optimizerLBFGS.optimize(ocCentralDQ);
  //optimizerLBFGS.optimize(ocAD);

  // The result (minimizer) can now be accessed at optimizerLBFGS.getControl().
  // The selected data for output and their plot can be seen in files
  // tmp/gnuplotData/data/log.dat and tmp/gnuplotData/log.png, respectively.

  // We could also have used any of the other optimization algorithms and opti-
  // cases here. Feel free to try it out. Some variants are slower/ less stable,
  // so e.g. convergence criteria might be adjusted.

  // Evaluation: compare result with exact solution
  const std::vector<T> exactSolution(N, 1);  // exact minimum (1,1,1,...)
  const T errorSolution = util::euklidDistance(
    optimizerLBFGS.getControl().data(), exactSolution.data(), N);
  const T normDerivative = util::euklidN(
    optimizerLBFGS.getDerivative().data(), N);

  // objective value at numerical solution shall be near zero
  assert (optimizerLBFGS.getObjective() < 1.0e-16);
  // error between numerical and exact minimizer shall be low
  assert (errorSolution < 1.0e-12);
  // derivative at the numerical minimizer shall be low
  assert (normDerivative < 1.0e-7);
}


////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  initialize(&argc, &argv);

  // set parameters for optimization algorithms (mamimal number of steps etc.)
  OptimizerParametersD<T> optiParameters;
  {
    using namespace parameters;
    optiParameters.set<optimizer::EPS>(1.e-7);
    optiParameters.set<optimizer::MAX_OPTIMIZATION_STEPS>(1000);
    optiParameters.set<linesearch::MAX_STEP_ATTEMPTS>(50);
    optiParameters.set<linesearch::STEP_CONDITION>("Wolfe");
    optiParameters.set<optimizer::VERBOSE>(true);
  }

  optimizationRosenbrock(optiParameters);

  return 0;
}
