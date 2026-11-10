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

#ifndef CASE_CASE_HH
#define CASE_CASE_HH

#include "case.h"

#include "io/textTable.h"

namespace olb {

template <typename MAP>
void ConcreteCase<MAP>::constructComponentsFromMesh() {
  /// Necessary due to impending invalidation
  _operators.clear();
  MAP::keys_t::for_each([&](auto name) {
    using NAME = typename decltype(name)::type;
    using COMPONENT = typename MAP::template value<NAME>;
    OLB_LOG(logging::Setup, Debug, "constructing ", fields::name<COMPONENT>());
    _components.set(name, std::make_unique<COMPONENT>(_mesh));
    // Stamp the component with its case name. Each lattice forwards this internally
    // to its statistics identity (clash-free CSV / plot naming).
    if constexpr (requires (COMPONENT& c) { c.setName(std::string{}); }) {
      get(NAME{}).setName(std::string{fields::name<NAME>()});
    }
    OLB_LOG(logging::Setup, Debug, "constructed ", fields::name<COMPONENT>());
  });
}

template <typename MAP>
ConcreteCase<MAP>::ConcreteCase(ParametersD& parameters, mesh_t& mesh)
  : _parameters{parameters}
  , _mesh{mesh}
  , _geometry{}
  , _components{}
{
  setReporter<reporter::Statistics>("statistics");
  setReporter<reporter::Timer>("timer");
  _parameters["reporter"].template declare<parameters::enabled>();
  // Answer --help while the parameters are complete but nothing large exists yet
  if (singleton::input().helpRequested()) {
    singleton::input().printHelp();
    OstreamManager clout(std::cout, "help");
    clout << std::endl << io::text::section("Case Parameters") << std::endl;
    _parameters.print(clout);
    std::exit(0);
  }
  _geometry.reset(new geometry_t(mesh));

#ifndef FEATURE_CASE_SILENT
  OstreamManager clout(std::cout, "Case");
  clout << std::endl;
  clout << io::text::section("Mesh") << std::endl;
  _mesh.print(clout);
  clout << std::endl;
#endif
  constructComponentsFromMesh();
#ifndef FEATURE_CASE_SILENT
  clout << io::text::section("Components") << std::endl;
  print(clout);
  clout << std::endl;
  clout << io::text::section("Parameters") << std::endl;
  _parameters.print(clout);
  clout << std::endl;
#endif
  if (singleton::mpi().isMainProcessor()
      && !singleton::input().resolveBuiltin<parameters::no_introspection>(false)) {
    std::ofstream fout(singleton::directories().getIntrospectionOutDir()
                       + io::text::filenameSafe(_name) + "_parameters.log",
                       std::ios::trunc);
    _parameters.print(fout);
  }
}

template <typename MAP>
void ConcreteCase<MAP>::print(OstreamManager clout) {
  io::text::KeyValueList components;
  MAP::keys_t::for_each([&](auto name) {
    using NAME = typename decltype(name)::type;
    using COMPONENT = typename MAP::template value<NAME>;
    components.item(std::string{fields::name<NAME>()}, std::string{fields::name<COMPONENT>()});
  });
  components.render(clout, "");
}

template <typename MAP>
void ConcreteCase<MAP>::report(std::size_t iT) {
  logging::servicePendingFlush();
  if (iT == 0) {
    singleton::input().warnUnknownArgs();
  }
  if (!_parameters["reporter"].template get<parameters::enabled>()) {
    return;
  }
  for (auto& [_, reporter] : _reporters) {
    reporter->report(iT);
  }
}

template <typename MAP>
void ConcreteCase<MAP>::finalize(std::size_t iT) {
  if (!_parameters["reporter"].template get<parameters::enabled>()) {
    return;
  }
  for (auto& [_, reporter] : _reporters) {
    reporter->finalize(iT);
  }
}
}

#endif
