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

#ifndef IO_OVERRIDE_SOURCE_H
#define IO_OVERRIDE_SOURCE_H

#include "cliReader.h"
#include "xmlReader.h"

#include <optional>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cctype>

namespace olb {

/// A frontend that can override a parameter's C++ default (CLI, XML, ...)
class OverrideSource {
public:
  /// Names a parameter without implying any frontend's addressing
  struct Query {
    /// Enclosing subgroups, outermost first ("fluid","wall"); empty at the root
    std::vector<std::string> group;
    std::string name;
    /// Lowercase native option (import_mesh -> --import-mesh)
    bool isNative = false;
    /// Nullary flag: presence alone means true
    bool isFlag = false;
  };

  virtual ~OverrideSource() = default;

  /// Override value for query, or nullopt if this source does not specify it
  virtual std::optional<std::string> lookup(const Query& query) const = 0;
};

/// Overrides from the command line
class CliOverrideSource final : public OverrideSource {
private:
  CLIreader& _cli;

  /// --<group.path.>NAME, natives hyphenated
  static std::string token(const Query& query) {
    std::string name = query.name;
    if (query.isNative) {
      std::replace(name.begin(), name.end(), '_', '-');
    }
    std::string path;
    for (const auto& segment : query.group) {
      path += segment + ".";
    }
    return "--" + path + name;
  }

public:
  explicit CliOverrideSource(CLIreader& cli)
    : _cli(cli)
  { }

  std::optional<std::string> lookup(const Query& query) const override {
    const std::string key = token(query);
    if (query.isFlag) {
      return _cli.flagValue(key);
    }
    if (!_cli.contains(key)) {
      return std::nullopt;
    }
    return _cli.getValueOrFallback<std::string>(key, "[]");
  }

  /// CLI token matched for query
  static std::string tokenOf(const Query& query) {
    return token(query);
  }
};

/// Overrides from OLB_-prefixed environment variables
class EnvOverrideSource final : public OverrideSource {
private:
  /// OLB_<GROUP_>NAME, uppercased (import_mesh -> OLB_IMPORT_MESH)
  static std::string variable(const Query& query) {
    std::string name = "OLB_";
    for (const auto& segment : query.group) {
      name += segment + "_";
    }
    name += query.name;
    std::transform(name.begin(), name.end(), name.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    return name;
  }

public:
  std::optional<std::string> lookup(const Query& query) const override {
    const char* value = std::getenv(variable(query).c_str());
    if (!value && query.name == "text_color" && std::getenv("NO_COLOR")) {
      return std::string{"never"};
    }
    if (!value) {
      return std::nullopt;
    }
    if (query.isFlag) {
      const std::string text{value};
      return text.empty() ? std::optional<std::string>{"true"} : text;
    }
    return std::string{value};
  }
};

/// Overrides from a case XML file's nested <Parameters> tree
class XmlOverrideSource final : public OverrideSource {
private:
  std::unique_ptr<XMLreader> _xml;

public:
  explicit XmlOverrideSource(const std::string& filename)
    : _xml(std::make_unique<XMLreader>(filename))
  { }

  std::optional<std::string> lookup(const Query& query) const override {
    _xml->setWarningsOn(false);
    XMLreader const* node = &(*_xml)["Parameters"];
    for (const auto& segment : query.group) {
      node->setWarningsOn(false);
      node = &(*node)[segment];
    }
    node->setWarningsOn(false);
    std::string value;
    const bool found = (*node)[query.name].read(value, false, false);
    node->setWarningsOn(true);
    _xml->setWarningsOn(true);
    if (found) {
      return value;
    }
    return std::nullopt;
  }
};

}

#endif
