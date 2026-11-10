/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Adrian Kummerlaender, Shota Ito, Julius Jeßberger, Mathias J. Krause
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

#ifndef CASE_CASE_H
#define CASE_CASE_H

#include "utilities/log/channels.h"
#include "mesh.h"
#include "parametersD.h"

#include "core/superLattice.h"

#include "utilities/typeMap.h"
#include "utilities/typeIndexedContainers.h"

#include "setters.h"
#include "reporter.h"
#include "statisticsReporter.h"
#include "timerReporter.h"
#include "refinement/vtkReporter.h"

#include "refinement/mesh.h"
#include "refinement/geometry.h"
#include "refinement/lattice.h"

namespace olb {

template <typename MAP>
class ConcreteCase {
public:
  /// Named components of this case
  using map_t = MAP;

  // Select first descriptor as reference for value type and dimension (for now)
  using reference_component_t = typename MAP::values_t::template get<0>;

  // Expose reference value type to use in apps with single lattice
  using value_t = typename reference_component_t::value_t;
  // Expose reference descriptor type to use in apps with single lattice
  using descriptor_t = typename reference_component_t::descriptor_t;

  template <typename VALUE_TYPE>
  using exchange_value_t = ConcreteCase<typename MAP::template map_values<
    meta::exchange_value_type<VALUE_TYPE>::template type
  >>;

  // Check first lattice to determine (un)refined state
  static constexpr bool is_refined = is_refined_lattice<reference_component_t>::value;

  using mesh_t = std::conditional_t<is_refined,
                                    RefinedMesh<value_t, reference_component_t::descriptor_t::d>,
                                    Mesh<value_t, reference_component_t::descriptor_t::d>>;

  using geometry_t = std::conditional_t<is_refined,
                                        RefinedGeometry<value_t, reference_component_t::descriptor_t::d>,
                                       SuperGeometry<value_t, reference_component_t::descriptor_t::d>>;

  using ParametersD = olb::ParametersD<typename reference_component_t::value_t,
                                       typename reference_component_t::descriptor_t>;

private:
  /// Reference to case parameters
  ParametersD& _parameters;

  /// Reference to mesh (shared between multiple cases)
  mesh_t& _mesh;

  /// Case-specific geometry
  std::unique_ptr<geometry_t> _geometry;

  /// Case-specific named components (e.g. lattices, ...)
  utilities::TypeIndexedTuple<typename MAP::template map_values<meta::unique_ptr_to>> _components;

  /// Arbitrary case-specific operators
  std::unordered_map<std::string, std::unique_ptr<ApplicableO>> _operators;

  /// Built-in reporters, keyed by name
  std::unordered_map<std::string, std::unique_ptr<Reporter>> _reporters;

  std::string _name{"Case"};

public:
  static constexpr unsigned d = reference_component_t::descriptor_t::d;

  template <typename NAME>
  using value_t_of = typename MAP::template value<NAME>::value_t;
  template <typename NAME>
  using descriptor_t_of = typename MAP::template value<NAME>::descriptor_t;

  /// (Re)construct all named components (stamping each with its case name) and reset operators
  void constructComponentsFromMesh();

  ConcreteCase(ParametersD& parameters, mesh_t& mesh);

  /// Print the component-to-lattice map
  void print(OstreamManager clout={std::cout, "Case"});

  /// Drive the registered reporters
  void report(std::size_t iT);
  /// Drive one end-of-run report per reporter
  void finalize(std::size_t iT);

  /// Register REPORTER under name; its parameters live in the like-named subgroup
  template <typename REPORTER, typename... ARGS>
  auto& setReporter(std::string name, ARGS&&... args) {
    using concrete_t = ConcreteReporter<ConcreteCase<MAP>, REPORTER>;
    auto ptr = std::make_unique<concrete_t>(*this, name, std::forward<ARGS>(args)...);
    auto& concrete = *ptr;
    _reporters[name] = std::move(ptr);
    return concrete;
  }

  Reporter& getReporter(std::string name) {
    return *_reporters.at(name);
  }

  /// Returns whether a reporter is registered under name
  bool hasReporter(std::string name) const {
    return _reporters.find(name) != _reporters.end();
  }

  /// Unregister the reporter under name; returns whether one was removed
  bool removeReporter(std::string name) {
    return _reporters.erase(name) > 0;
  }

public:

  ParametersD& getParameters() {
    return _parameters;
  }

  auto& getMesh() {
    return _mesh;
  }

  auto& getGeometry() {
    return *_geometry;
  }

  /// Returns reference to component of NAME
  template <typename NAME>
  auto& get(NAME) {
    return *_components.template get<NAME>();
  }

  /// Returns reference to lattice of NAME
  template <typename NAME>
  auto& getLattice(NAME name) {
    using T = value_t_of<NAME>;
    using DESCRIPTOR = descriptor_t_of<NAME>;
    using COMPONENT = typename MAP::template value<NAME>;
    if constexpr (is_refined_lattice<COMPONENT>::value) {
      return static_cast<RefinedLattice<T, DESCRIPTOR>&>(get(name));
    }
    else {
      return static_cast<Lattice<T, DESCRIPTOR>&>(get(name));
    }
  }

  ApplicableO& getOperator(std::string name) {
    return *_operators.at(name);
  }

  template <typename... ARGS>
  auto& setCouplingOperator(std::string name, ARGS&&... args) {
    auto ptr = constructUniqueCoupling(std::forward<ARGS&&>(args)...);
    auto& concrete = *ptr;
    _operators[name].reset(ptr.release());
    return concrete;
  }

  void resetGeometry() {
    _geometry.reset(new geometry_t(_mesh));
  }

  void resetLattices() {
    constructComponentsFromMesh();
  }

  void setName(std::string name) {
    _name = name;
  }

  template <typename NAME>
  bool hasName(NAME) {
    return  (NAME{}.name == _name);
  }
};

template <typename... MAP>
using Case = ConcreteCase<meta::map<MAP...>>;

template <typename T>
struct is_case {
  static constexpr bool value = false;
};

template <typename... MAP>
struct is_case<ConcreteCase<MAP...>> {
  static constexpr bool value = true;
};

template <typename T>
static constexpr bool is_case_v = is_case<T>::value;

}

#endif
