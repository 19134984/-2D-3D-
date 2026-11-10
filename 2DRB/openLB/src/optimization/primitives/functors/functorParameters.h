/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
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

#ifndef FUNCTOR_PARAMETERS_H
#define FUNCTOR_PARAMETERS_H

#include "case/parametersD.h"

namespace olb {

namespace functors {

template <concepts::PrimitiveFunctor FUNCTOR, typename T, typename DESCRIPTOR>
class FunctorParameters {
private:
  using type = FUNCTOR;

  ParametersD<T,DESCRIPTOR> _parameters;
  std::unordered_set<std::type_index> _missingParameters; // list of unset parameters
public:
  using parameters_t = typename FUNCTOR::parameters_t;

  FunctorParameters() {
    parameters_t::for_each([&](auto id) {
      _missingParameters.insert(id);
    });
  };

  template <typename KEY, typename VALUE>
    requires(parameters_t::template contains<KEY>())
  void set(VALUE&& value) {
    _missingParameters.erase(std::type_index(typeid(KEY)));
    return _parameters.template set<KEY>(std::forward<VALUE>(value));
  }

  // get value of parameter by KEY
  template <typename KEY>
  auto get() {
    return _parameters.template get<KEY>();
  };

  bool allParametersAreSet() {
    return _missingParameters.empty();
  }
};

// Helper function which sets all parameters in f to the passed instance
template <typename INSTANCE, typename FUNCTOR_PARAMETERS>
void setAllParametersIn(INSTANCE& instance, FUNCTOR_PARAMETERS& f) {
  if (!f.allParametersAreSet()) {
    throw std::runtime_error("setAllParametersIn: f is not fully parameterized");
  }
  FUNCTOR_PARAMETERS::parameters_t::for_each([&](auto id) {
    using KEY = typename decltype(id)::type;
    instance.template setParameter<KEY>(f.template get<KEY>());
  });
}

}

}

#endif
