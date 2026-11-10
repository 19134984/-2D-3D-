/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2016 Mathias J. Krause, Benjamin Förster
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
 * The description of optimization algorithms -- header file.
 */


#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "optimization/core/optiCase.h"
#include "communication/mpiManager.h"
#include "io/gnuplotWriter.h"
#include "io/xmlReader.h"
#include "core/singleton.h"
#include "case/parameters.h"

namespace olb {

namespace parameters::optimizer {
  struct EPS                    : public olb::fields::TYPED_SCALAR<double,1e-10> { };
  struct FILE_NAME              : public olb::fields::STRING<"control.dat"> { };
  struct MAX_OPTIMIZATION_STEPS : public olb::fields::SCALAR<100> { };
  struct VERBOSE                : public olb::fields::TYPED_SCALAR<bool, true> { };
  struct LOG_FILE_NAME          : public olb::fields::STRING<"optiLog"> { };
  struct WITH_UPPER_BOUND       : public olb::fields::TYPED_SCALAR<bool, false> { };
  struct UPPER_BOUND            : public olb::fields::SCALAR<0.> { };
  struct WITH_LOWER_BOUND       : public olb::fields::TYPED_SCALAR<bool, false> { };
  struct LOWER_BOUND            : public olb::fields::SCALAR<0.> { };
  struct VECTOR_BOUNDS          : public olb::fields::TYPED_SCALAR<bool, false> { };
  struct CONTROL_EPS            : public olb::fields::SCALAR<std::numeric_limits<double>::epsilon()> { };
  struct FAIL_ON_MAX_ITER       : public olb::fields::TYPED_SCALAR<bool, true> { };
}

/// \namespace opti Optimization Code.
namespace opti {

template<typename> struct AbstractOptiCase;

enum OptimizerLogType {control, derivative, objective, error, rel_error, norm_derivative, weighted_norm_derivative};

template <typename T>
using OptimizerParametersD = ParametersD<T, descriptors::SPATIAL_DESCRIPTOR<0>>;

/// Interface for the use of various optimization algorithms.
/**
 * This class is intended to be derived from.
 * S is the underlying arithmetic data type
 * C is the container type that is used to store the vectors (e.g. std::vector<S>)
 * C is expected to allow construction via ContainerCreator class
 */
template<typename S, typename C>
class Optimizer {

private:
  mutable OstreamManager clout;

protected:
  /// Number of controlled variables
  int                                 _dimCtrl;
  /// Vector of controlled variables (size _dimCtrl)
  C                                   _control;
  /// Value of the objective functional evaluated
  /// for controlled variables saved in _control
  S                                   _value;
  /// Vector of derivatives of the object functional
  /// with respect to the controlled variables
  C                                   _derivative;
  /// Current iteration no.
  int                                 _it;
  /// Maximal number of iteration
  int                                 _maxIt;
  /// Fail when max number of iteration reached
  bool                                _failOnMaxIter;
  /// Optimizer stops if |_derivatives| < _eps
  S                                   _eps;
  /// Verbose
  bool                                _verboseOn;

  /// Bounded versions
  bool                                _withUpperBound;
  bool                                _withLowerBound;
  bool                                _vectorBounds;
  C                                   _boundedControl;
  C                                   _upperBound;
  C                                   _lowerBound;

  /// For setting tolerance of controls
  bool _controlsConverged;
  S _controlEps;

  /// Provides the Optimizer with methods to evaluate the
  /// value of an object functional and compute derivatives
  AbstractOptiCase<S>*                      _optiCase;

  /// control vector to compare with (for numerical evaluation)
  C                                   _referenceControl;
  S                                   _dx{1.0};

public:
  Gnuplot<S>                          gplot;
  /// For defining what kind of gnuplot analysis is wanted, if empty vector - no analysis,
  /// value, control and derivative are the possible options
  std::vector<OptimizerLogType>       _gplotAnalysis;

  Optimizer(int dimCtrl, S eps, int maxIt,
        bool verboseOn=true, const std::string fname="",
        const std::string logFileName="",
        bool withUpperBound=false, S upperBound=S(),
        bool withLowerBound=false, S lowerBound=S(),
        bool vectorBounds=false,
        S controlEps=S(std::numeric_limits<double>::epsilon() ),
        bool failOnMaxIter = true,
        std::vector<OptimizerLogType> gplotAnalysis = {OptimizerLogType::objective});

  Optimizer(int dimCtrl, OptimizerParametersD<S>& optimizerParameters,
    std::vector<OptimizerLogType> gplotAnalysis = {OptimizerLogType::objective}) :
    Optimizer(
      dimCtrl,
      optimizerParameters.template get<parameters::optimizer::EPS>(),
      optimizerParameters.template get<parameters::optimizer::MAX_OPTIMIZATION_STEPS>(),
      optimizerParameters.template get<parameters::optimizer::VERBOSE>(),
      optimizerParameters.template get<parameters::optimizer::FILE_NAME>(),
      optimizerParameters.template get<parameters::optimizer::LOG_FILE_NAME>(),
      optimizerParameters.template get<parameters::optimizer::WITH_UPPER_BOUND>(),
      optimizerParameters.template get<parameters::optimizer::UPPER_BOUND>(),
      optimizerParameters.template get<parameters::optimizer::WITH_LOWER_BOUND>(),
      optimizerParameters.template get<parameters::optimizer::LOWER_BOUND>(),
      optimizerParameters.template get<parameters::optimizer::VECTOR_BOUNDS>(),
      optimizerParameters.template get<parameters::optimizer::CONTROL_EPS>(),
      optimizerParameters.template get<parameters::optimizer::FAIL_ON_MAX_ITER>(),
      gplotAnalysis)
  {
    OstreamManager clout(std::cout, "Optimizer");
    clout << io::text::section("Optimizer parameters") << std::endl;
    optimizerParameters.print(clout);
    clout << std::endl;
  };

  virtual ~Optimizer() { };

  virtual void optimizationStep() = 0;

  void maxIterationReached();

  virtual void optimize();

  virtual void optimize(AbstractOptiCase<S>& optiCase) {
    _optiCase = &optiCase;
    setControl(_optiCase->getControlVector());
    optimize();
  }

  void simulate() {
    evaluateObjective(_control, _value);
  }

  void simulate(AbstractOptiCase<S>& optiCase) {
    _optiCase = &optiCase;
    simulate();
  }

  void evaluateObjective(const C& control, S& result) {
    result = _optiCase->computeObjective(control, _it);
  }

  void computeDerivatives(const C& control, C& derivatives) {
    derivatives = _optiCase->computeDerivatives(control, _it);
  }

  /// Prints information of the current optimization step it
  void print(int it);

  void setControl(C& control) {
    _control = std::move(control);
  }

  const C& getControl() const {
    return _control;
  }

  const C& getDerivative() const {
    return _derivative;
  }

  const S& getObjective() const {
    return _value;
  }

  int getIteration() const {
    return _it;
  }

  /// Enable logging regarding various parameters
  void addLogging(std::vector<OptimizerLogType> types) {
    _gplotAnalysis.insert(_gplotAnalysis.end(), types.begin(), types.end());
  };

  /// Writes the current control variables linewise into file fname
  void writeControlToFile(const std::string fname="control.dat");

  /// Reads the latest control variables from file fname
  void readControlFromFile(const std::string fname="control.dat");

  void setStartValue(S startValue);

  void setStartValue(C startValues, S factor=S(1));

  AbstractOptiCase<S>* getOptiCase() { return _optiCase; }

  void setOptiCase(AbstractOptiCase<S>* optiCase) { _optiCase = optiCase; }

  void setGnuplotData();

  /// set the reference value for the control vector (exact solution)
  void setReferenceControl(C result) {
    _referenceControl = result;
  }

  void setMeshSpacing(S meshSpacing) {
    _dx = meshSpacing;
  }
};

/// the gplotAnalysisString is gained from the xml file and the function than separates and prepares it to be used in the constructor
void getGnuplotTagsFromString(std::string gplotAnalysisString, std::vector<OptimizerLogType>& gplotAnalysis );

} // namespace opti

} // namespace olb

#endif
