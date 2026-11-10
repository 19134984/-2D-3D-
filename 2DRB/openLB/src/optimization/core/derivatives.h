/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Shota Ito
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

#ifndef DERIVATIVES_H
#define DERIVATIVES_H

namespace olb {

namespace opti {

namespace derivatives {

/// @brief   Logic to compute objective derivatives by manual definition
// Derivative needs to be passed via a callable (e.g. currently for adjoint optimization)
struct Manual {
  using parameters_t = meta::list<>;

  template <typename OBJECTIVE, typename DERIVATIVE, typename T=OBJECTIVE::result_type>
  static DERIVATIVE getDerivativeF(OBJECTIVE& objective) {
    return [&](const std::vector<T>& controls) -> std::vector<T> {
      throw std::logic_error("Logic for derivative");
    };
  }
};

/// @brief   Logic to compute objective derivatives using forward difference quotients
struct FDQ {
  using parameters_t = meta::list<parameters::FDQ_STENCIL_WIDTH>;

  template <typename OBJECTIVE, typename DERIVATIVE, typename PARAMETERS, typename T=OBJECTIVE::result_type>
  static DERIVATIVE getDerivativeF(OBJECTIVE& objective, PARAMETERS& parameters) {
    // return a lambda function defining the FDQ logic
    return [&](const std::vector<T>& controls) -> std::vector<T> {
      const T stencilWidth = parameters.template get<parameters::FDQ_STENCIL_WIDTH>();
      std::vector<T> derivatives(controls.size());

      std::vector<T> shiftedControl = controls;
      const T result = objective(controls);
      for (std::size_t d=0; d < controls.size(); ++d) {
        shiftedControl[d] += stencilWidth;
        const T shiftedResult = objective(shiftedControl);
        // FDQ stencil
        derivatives[d] = (shiftedResult - result) / stencilWidth;
        shiftedControl[d] = controls[d];
      }
      return derivatives;
    };
  }
};

/// @brief   Logic to compute objective derivatives using central difference quotients
struct CDQ {
  using parameters_t = meta::list<parameters::CDQ_STENCIL_WIDTH>;

  template <typename OBJECTIVE, typename DERIVATIVE, typename PARAMETERS, typename T=OBJECTIVE::result_type>
  // return a lambda function defining the CDQ logic
  static DERIVATIVE getDerivativeF(OBJECTIVE& objective, PARAMETERS& parameters) {
    return [&](const std::vector<T>& controls) -> std::vector<T> {
      const T stencilWidth = parameters.template get<parameters::CDQ_STENCIL_WIDTH>();
      std::vector<T> derivatives(controls.size());

      std::vector<T> shiftedControl = controls;
      for (std::size_t d=0; d<controls.size(); ++d) {
        shiftedControl[d] += stencilWidth;
        const T shiftedObjective_plus = objective(shiftedControl);
        shiftedControl[d] = controls[d] - stencilWidth;
        const T shiftedObjective_minus = objective(shiftedControl);
        // CDQ stencil
        derivatives[d] = T(0.5) * (shiftedObjective_plus - shiftedObjective_minus) / stencilWidth;
        shiftedControl[d] = controls[d];
      }
      return derivatives;
    };
  }
};

/// @brief   Logic to compute objective derivatives by automatic differentiation forward mode
struct ADf {
  using parameters_t = meta::list<>;

  // Define how to compute the derivative using the OptiCase instance
  template <typename OBJECTIVE, typename DERIVATIVE, typename T>
  static DERIVATIVE getDerivativeF(OBJECTIVE& objective) {
    return [&](const std::vector<T>& controls) -> std::vector<T> {
      using U = typename OBJECTIVE::result_type;
      std::vector<U> tmp = util::copyAs<U>(controls);
      util::iniDiagonal(tmp);
      U derivative = objective(tmp);
      return util::extractDerivatives(derivative);
    };
  }

};

} // namespace derivatives

} // namespace opti

} // namepsace olb

#endif
