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

#ifndef IO_INPUT_H
#define IO_INPUT_H

#include "overrideSource.h"
#include "cliReader.h"
#include "ostreamManager.h"
#include "io/textTable.h"
#include "core/parameters.h"
#include "core/vector.h"
#include "utilities/log/log.h"

#include <cstdlib>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <cctype>

namespace olb {

/// The override frontends (CLI, XML, ...) that ParametersD resolves against
class InputRegistry {
private:
  /// Sources in descending precedence
  std::vector<std::unique_ptr<OverrideSource>> _sources;
  std::optional<std::string> _caseFile;
  std::map<std::string,OverrideSource::Query> _builtins;

  /// Query addressing a scalar built-in field (token / env var derived from the field name)
  template <typename FIELD>
  static OverrideSource::Query builtinQuery() {
    const std::string name = parameters::name<FIELD>();
    const bool isNative = parameters::isNativeName(name);
    const bool isFlag = isNative
      && std::is_same_v<typename FIELD::template value_type<double>, bool>;
    return OverrideSource::Query{{}, name, isNative, isFlag};
  }

  /// Document FIELD under --help, deriving its token and arity from the field name
  template <typename FIELD>
  void declareBuiltin(const std::string& description, const std::string& defaultValue = "") {
    const auto query = builtinQuery<FIELD>();
    const auto token = CliOverrideSource::tokenOf(query);
    _builtins.insert_or_assign(token, query);
    singleton::cli().declare(token, description, query.isFlag, defaultValue);
  }

  /// Document the OpenLB built-in options (listed by --help)
  void registerBuiltins() {
    declareBuiltin<parameters::no_introspection>(
      "Skip the per-component introspection dump written at lattice initialization", "off");
    singleton::cli().declare("--case-file",
      "Path to the case XML", false, "case.xml");
    singleton::cli().declare("--text-style",
      "Status-line grammar: semicolon, keyvalue or labeled", false, "semicolon");
    singleton::cli().declare("--table-style",
      "Table layout: ruled, columns, grid, framed or markdown", false, "ruled");
    declareBuiltin<parameters::text_ascii>(
      "Use ASCII instead of box-drawing glyphs", "off");
    singleton::cli().declare("--detail",
      "Statistics detail: glance, compact, summary or full", false, "compact");
    singleton::cli().declare("--text-color",
      "ANSI color: never, auto or always", false, "never");
    declareBuiltin<parameters::clout_timestamps>(
      "Prefix stdout lines with a timestamp", "off");
#ifdef FEATURE_PROFILING
    declareBuiltin<parameters::profile>(
      "Record a performance trace", "off");
#endif
#ifdef PARALLEL_MODE_OMP
    declareBuiltin<parameters::omp_affinity>(
      "Auto-manage OMP thread count and affinity", "on");
#endif
    singleton::cli().declare("--log",
      "Global log level: off or debug", false, "off");
    const auto channels = logging::channelNames();
    const auto levels   = logging::channelLevels();
    for (std::size_t i = 0; i < channels.size(); ++i) {
      std::string level(logging::name(levels[i]));
      std::transform(level.begin(), level.end(), level.begin(),
                     [](unsigned char c){ return std::tolower(c); });
      singleton::cli().declare("--log." + channels[i],
        "Log level for the " + channels[i] + " channel", false, level);
    }
    singleton::cli().declare("--help",
      "Print the built-in options and exit", true, "off");
  }

public:
  /// Whether --help was given; answered by Case once its parameters are complete
  bool helpRequested() const {
    return singleton::cli().contains("--help");
  }

  /// Print the documented built-in options
  void printHelp() const {
    OstreamManager clout(std::cout, "help");
    clout << std::endl << io::text::section("Built-in options") << std::endl;
    io::text::Table table({{"Option"}, {"Default"}, {"Description"}});
    for (const auto& [token, option] : singleton::cli().known()) {
      table.row(token, option.defaultValue.empty() ? "-" : option.defaultValue, option.description);
    }
    table.render(clout);
    clout << std::endl;
    clout << "Set any parameter via `--NAME VALUE`, `OLB_NAME=VALUE` or a case XML."
          << std::endl;
    clout << "Precedence: CLI > env (OLB_*) > XML > C++ defaults" << std::endl;
  }

public:
  InputRegistry() {
    _sources.push_back(std::make_unique<CliOverrideSource>(singleton::cli()));
    _sources.push_back(std::make_unique<EnvOverrideSource>());
  }

  /// Drop all sources but the CLI and environment ones
  void reset() {
    _sources.clear();
    _sources.push_back(std::make_unique<CliOverrideSource>(singleton::cli()));
    _sources.push_back(std::make_unique<EnvOverrideSource>());
    _caseFile.reset();
  }

  /// Resolve query against the sources in precedence order
  std::optional<std::string> resolve(const OverrideSource::Query& query) const {
    for (const auto& source : _sources) {
      if (auto value = source->lookup(query)) {
        return value;
      }
    }
    return std::nullopt;
  }

  /// Resolve a scalar built-in field, parsing the override through the field's own type
  template <typename FIELD>
  auto resolveBuiltin(typename FIELD::template value_type<double> fallback) const {
    using V = typename FIELD::template value_type<double>;
    if (auto raw = resolve(builtinQuery<FIELD>())) {
      if (auto parsed = Vector<V,1>::fromString(*raw)) {
        return (*parsed)[0];
      }
    }
    return fallback;
  }

  /// Render every declared option with its effective value
  void print(std::ostream& out, const std::string& indent = "  ") {
    io::text::KeyValueList list;
    for (const auto& [token, option] : singleton::cli().known()) {
      std::string value = option.defaultValue;
      if (auto builtin = _builtins.find(token); builtin != _builtins.end()) {
        if (auto resolved = resolve(builtin->second)) {
          value = *resolved;
        }
      } else if (singleton::cli().contains(token)) {
        value = singleton::cli().getValueOrFallback<std::string>(token, value);
      }
      list.item(token, value);
    }
    list.render(out, indent);
  }

  /// Append a case XML as an override source
  void loadCaseXml(const std::string& filename) {
    _caseFile = filename;
    _sources.push_back(std::make_unique<XmlOverrideSource>(filename));
  }

  /// The loaded case-file path, if any
  std::optional<std::string> getCaseFile() const {
    return _caseFile;
  }

  /// Document an option under --help
  void declare(const std::string& token, const std::string& description = "", bool isFlag = false,
               const std::string& defaultValue = "") {
    singleton::cli().declare(token, description, isFlag, defaultValue);
  }

  /// Capture the CLI, load an optional case XML and handle --help
  void load(int& argc, char** argv) {
    reset();
    singleton::cli().store(argc, argv);
    registerBuiltins();

    const std::string caseFile =
      singleton::cli().getValueOrFallback<std::string>("--case-file", "case.xml");
    if (std::filesystem::exists(caseFile)) {
      loadCaseXml(caseFile);
    }
  }

  /// Warn about --tokens no source consumed; only accurate once parameters resolved
  void warnUnknownArgs() const {
    for (const auto& arg : singleton::cli().getUnconsumed()) {
      OstreamManager clout(std::cout, "input");
      clout << "Warning: unknown parameter \"" << arg << "\" provided via CLI" << std::endl;
    }
  }
};

namespace singleton {

/// Process-global input authority
inline InputRegistry& input() {
  static InputRegistry instance;
  return instance;
}

}

}

#endif
