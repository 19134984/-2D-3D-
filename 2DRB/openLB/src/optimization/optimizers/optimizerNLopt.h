/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Shota Ito
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
 * Interface to external optimization libraries; here for NLopt
 * Link: https://nlopt.readthedocs.io/en/latest/
 */

#ifndef OPTIMIZER_EXTERNAL_H
#define OPTIMIZER_EXTERNAL_H

#include "optimizer.h"

#ifdef FEATURE_NLOPT
#include "nlopt.hpp"
#endif

namespace olb {

namespace opti {

template <typename T>
struct AbstractDifferentiableFunction {
  AbstractDifferentiableFunction() = default;

  virtual T evaluate(const std::vector<T>& controls);
  virtual std::vector<T> computeDerivatives(const std::vector<T>& controls);
};


/// Gradients are just passed as a function (and not computed by an own routine)
template <typename T>
class DifferentiableFunction : public AbstractDifferentiableFunction<T> {
private:
  std::function<T(const std::vector<T>&)> _function;
  std::function<std::vector<T>(const std::vector<T>&)> _derivative;

public:
  DifferentiableFunction(
    std::function<T(const std::vector<T>&)> function,
    std::function<std::vector<T>(const std::vector<T>&)> derivative)
   : _function(function), _derivative(derivative)
  { }

  T evaluate(const std::vector<T>& controls) override {
    return _function(controls);
  }

  std::vector<T> computeDerivatives(const std::vector<T>& controls) override {
    return _derivative(controls);
  }
};

/// Gradient computed through ADf
template<typename T, typename ADf>
class DifferentiableFunctionByADf : public AbstractDifferentiableFunction<T>{
private:
  std::function<T(const std::vector<T>&)> _function;
  std::function<ADf(const std::vector<ADf>&)> _adFunction;

public:
  DifferentiableFunctionByADf (
    std::function<T(const std::vector<T>&)> function,
    std::function<ADf(const std::vector<ADf>&)> adFunction)
   : _function(function), _adFunction(adFunction)
  { }

  T evaluate(const std::vector<T>& controls) override
  {
    return _function(controls);
  }

  std::vector<T> computeDerivatives(const std::vector<T>& controls) override
  {
    auto adControls = util::iniAD<ADf::dim,T,std::vector<T>>(controls);

    const ADf adResult = _adFunction(adControls);

    std::vector<T> derivatives(controls.size());
    for(std::size_t iDim = 0; iDim < controls.size(); ++iDim){
      derivatives[iDim] = adResult.d(iDim);
    }
    return derivatives;
  }
};

#ifdef FEATURE_NLOPT

/// Interface optimizer for NLopt; a open-source library for nonlinear optimization
/// Link: https://nlopt.readthedocs.io/en/latest/
template<typename T>
class OptimizerNLopt : public Optimizer<T,std::vector<T>> {
private:
  mutable OstreamManager clout;
  // NLopt optimizer
  nlopt::opt _optimizer;
  // Context passed to NLopt callbacks
  struct NLoptWrapper {
    OptimizerNLopt<T>* optimizer;
    AbstractDifferentiableFunction<T>* function;

    NLoptWrapper(OptimizerNLopt<T>* opt, AbstractDifferentiableFunction<T>* func)
    : optimizer(opt), function(func) {}
  };

  std::size_t _iT = 0;
  std::vector<std::unique_ptr<NLoptWrapper>> _wrappers;
  std::unique_ptr<AbstractDifferentiableFunction<T>> _objective;

  auto createNLoptFunction() {
    return [](const std::vector<T>& controls, std::vector<T>& derivatives, void* data) -> T {

      // Retrieve wrapper containing optimizer and function pointer
      auto* wrapper = reinterpret_cast<NLoptWrapper*>(data);
      auto* opt = wrapper->optimizer;
      auto* func = wrapper->function;

      // Evaluate objective or constraint value
      T value = func->evaluate(controls);
      opt->_value = value;

      // Compute gradient only if requested by NLopt
      if (derivatives.size() == controls.size()) {
        derivatives = func->computeDerivatives(controls);
        for (int i=0; i<derivatives.size(); ++i) {
          opt->_derivative[i] = derivatives[i];
        }
      }

      // Increment iteration counter only for objective evaluations
      if (func == opt->_objective.get()) {
        opt->_iT++;
        opt->print(opt->_iT);
      }

      return value;
    };
  }

public:
  /// Construction of an OptimizerConstrainedBFGS
  OptimizerNLopt(int dimCtrl, T eps, int maxIt,
                 int maxStepAttempts,
                 bool withUpperBound=false, T upperBound=T(),
                 bool withLowerBound=false, T lowerBound=T(),
                 std::vector<OptimizerLogType> gplotAnalysis = {})
    : Optimizer<T,std::vector<T>>(dimCtrl, eps, maxIt, true, "", "", withUpperBound, upperBound, withLowerBound,
                                  lowerBound, false, std::numeric_limits<T>::epsilon(), true, gplotAnalysis),
                                  clout(std::cout,"OptimizerNLopt")
  {
    static_assert(std::is_same<T, double>(), "NLopt optimizers require double precision.");

    // Default algorithm is set to SLSQP
    // https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#slsqp
    _optimizer = nlopt::opt("LD_SLSQP", (unsigned)(this->_dimCtrl));

    // Bound controls
    if( withLowerBound ){
      clout << "Set lower bound for optimizer: " << lowerBound << std::endl;
      _optimizer.set_lower_bounds(lowerBound);
    }
    if( withUpperBound ){
      clout << "Set upper bound for optimizer: " << upperBound << std::endl;
      _optimizer.set_upper_bounds(upperBound);
    }
    // Necessary stopping criteria
    _optimizer.set_xtol_rel(this->_controlEps);
    _optimizer.set_ftol_rel(this->_eps);
  };

  // Define which NLopt algorithm should be used
  void setNLoptAlgorithm(const char* algorithm) {
    _optimizer = nlopt::opt(algorithm, (unsigned)(this->_dimCtrl));
  }

  // Parameterize NLopt optimizer
  void setParameter(const char* name, T value) {
    _optimizer.set_param(name, value);
  }

  // Inequality constraint
  void addInequalityConstraint(AbstractDifferentiableFunction<T>& f) {
    _wrappers.emplace_back(std::make_unique<NLoptWrapper>(this, &f));
    _optimizer.add_inequality_constraint(createNLoptFunction(), _wrappers.back().get());
  }
  // Equality constraint
  void addEqualityConstraint(AbstractDifferentiableFunction<T>& f) {
    _wrappers.emplace_back(std::make_unique<NLoptWrapper>(this, &f));
    _optimizer.add_equality_constraint(createNLoptFunction(), _wrappers.back().get());
  }

  // Additional NLopt stopping criteria
  void setMaxEvaluations(int eval) {
    _optimizer.set_maxeval(eval);
  }

  void optimize(AbstractOptiCase<T>& optiCase) override {
    this->_optiCase = &optiCase;
    _iT = 0;
    _optimizer.set_maxeval(this->_maxIt);
    _objective = std::make_unique<DifferentiableFunction<T>>(
      [&](const std::vector<T>& controls) -> T {return optiCase.computeObjective(controls, this->_iT);},
      [&](const std::vector<T>& controls) -> std::vector<T> {return optiCase.computeDerivatives(controls, this->_iT);}
    );

    _wrappers.emplace_back(std::make_unique<NLoptWrapper>(this, _objective.get()));
    _optimizer.set_min_objective(createNLoptFunction(), _wrappers.back().get());

    try{
      clout << "Starting NLopt optimization..." << std::endl;
      T finalObjective = 0.;
      this->setControl(this->_optiCase->getControlVector());
      _optimizer.optimize(this->_control, finalObjective);

      clout << "Found minimum with objective value = " << finalObjective << std::endl;
    }
    catch(std::exception &e) {
      std::cerr << "nlopt failed: " << e.what() << std::endl;
    }
  }

  void optimizationStep() override { };

  int getNumEvaluations() {
    return _optimizer.get_numevals();
  }

  int getDimension() {
    return _optimizer.get_dimension();
  }

  unsigned getNumParams() {
    return _optimizer.num_params();
  }

};

#endif

}

}

#endif
