/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Adrian Kummerlaender
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

#ifndef CASE_REPORTER_H
#define CASE_REPORTER_H

#include <cstddef>
#include <string>
#include <utility>

namespace olb {

/// Base of any case reporter, driven by Case::report
struct Reporter {
  virtual ~Reporter() = default;

  /// iT is the reference-level step; a refined level N advances 2^N substeps below it
  virtual void report(std::size_t iT) = 0;
  /// Ungated end-of-run report
  virtual void finalize(std::size_t iT) { }
};

/// Provides REPORTER with the concrete values of its declared parameters
template <typename CASE, typename REPORTER>
class ConcreteReporter final : public Reporter {
private:
  using value_t      = typename CASE::value_t;
  using descriptor_t = typename CASE::ParametersD::descriptor_t;
  using parameters_t = typename StaticParametersD<value_t,descriptor_t>
                         ::template include<typename REPORTER::parameters>;

  CASE&        _case;
  std::string  _group;
  REPORTER     _reporter;
  parameters_t _parameters;
  std::size_t  _version = 0;
  bool         _extracted = false;

  parameters_t& parameters() {
    auto& source = _case.getParameters()["reporter"][_group];
    if (!_extracted || source.version() != _version) {
      _parameters = source.template extract<typename REPORTER::parameters>();
      _version    = source.version();
      _extracted  = true;
    }
    return _parameters;
  }

public:
  template <typename... ARGS>
  ConcreteReporter(CASE& myCase, std::string group, ARGS&&... args)
    : _case(myCase)
    , _group(std::move(group))
    , _reporter(std::forward<ARGS>(args)...)
  {
    parameters();
  }

  void report(std::size_t iT) override {
    _reporter.report(_case, iT, parameters());
  }

  void finalize(std::size_t iT) override {
    if constexpr (requires { _reporter.finalize(_case, iT, parameters()); }) {
      _reporter.finalize(_case, iT, parameters());
    }
  }
};

}

#endif
