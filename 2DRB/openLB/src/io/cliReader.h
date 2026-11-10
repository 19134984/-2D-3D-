/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2023 Adrian Kummerlaender
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net
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

#ifndef IO_CLI_READER_H_
#define IO_CLI_READER_H_

#include <string>
#include <vector>
#include <map>
#include <charconv>
#include <optional>
#include <algorithm>

namespace olb {

/// Very simple CLI argument parser
class CLIreader {
public:
  /// A CLI option the runner knows about (built-in native or documented flag).
  /// Registered up front so it is listed by --help and not mis-reported as unknown.
  struct KnownOption {
    std::string description;
    bool isFlag = false;
    std::string defaultValue;
  };

private:
  std::vector<std::string> _tokens;
  std::vector<bool> _consumed;
  /// Documented options keyed by their full "--token"
  std::map<std::string, KnownOption> _known;

public:
  CLIreader() = default;

  CLIreader(int& argc, char** argv) {
    store(argc, argv);
  }

  /// (Re)populate the token list from argc/argv; used to fill the singleton at initialize()
  void store(int argc, char** argv) {
    _tokens.clear();
    for (int i=1; i < argc; ++i) {
      _tokens.emplace_back(argv[i]);
    }
    _consumed.assign(_tokens.size(), false);
    _known.clear();
  }

  /// Register a documented option (full "--token") so it feeds --help and is not
  /// reported as unknown even when resolved lazily (or never touched this run).
  void declare(const std::string& token, const std::string& description = "", bool isFlag = false,
               const std::string& defaultValue = "") {
    _known[token] = KnownOption{description, isFlag, defaultValue};
  }

  /// Whether token was registered via declare()
  bool isKnown(const std::string& token) const {
    return _known.find(token) != _known.end();
  }

  /// The documented options, for --help
  const std::map<std::string, KnownOption>& known() const {
    return _known;
  }

  /// Returns true iff name is specified
  bool contains(const std::string& name) {
    auto iter = std::find(_tokens.begin(), _tokens.end(), name);
    if (iter != _tokens.end()) {
      _consumed[std::distance(_tokens.begin(), iter)] = true;
      return true;
    }
    return false;
  }

  /// Returns value of token after name (i.e. by convention the value assigned to name)
  std::string operator[](const std::string& name) {
    auto iter = std::find(_tokens.begin(), _tokens.end(), name);
    if (iter != _tokens.end()) {
      _consumed[std::distance(_tokens.begin(), iter)] = true;
      if (++iter != _tokens.end()) {
        _consumed[std::distance(_tokens.begin(), iter)] = true;
        return *iter;
      }
    }
    return std::string{};
  }

  /// Flag lookup: nullopt if absent, its explicit value if the next token is not an
  /// option, else "true" for bare presence
  std::optional<std::string> flagValue(const std::string& name) {
    auto iter = std::find(_tokens.begin(), _tokens.end(), name);
    if (iter == _tokens.end()) {
      return std::nullopt;
    }
    _consumed[std::distance(_tokens.begin(), iter)] = true;
    auto next = iter + 1;
    if (next != _tokens.end() && next->rfind("--", 0) != 0) {
      _consumed[std::distance(_tokens.begin(), next)] = true;
      return *next;
    }
    return std::string{"true"};
  }

  /// Return value of name as TYPE or fallback if not provided
  template <typename TYPE>
  TYPE getValueOrFallback(const std::string& name, TYPE fallback) {
    if (contains(name)) {
      const std::string str = operator[](name);
      TYPE value{};
      // Handling required due to inconsistent support for std::from_chars
             if constexpr (std::is_integral_v<TYPE>) {
        auto result = std::from_chars(str.data(),
                                      str.data() + str.size(),
                                      value);
        if (result.ec != std::errc()) {
          return fallback;
        }
        return value;
      } else if constexpr (std::is_floating_point_v<TYPE>) {
        try {
          std::size_t nProcessed{};
          value = static_cast<TYPE>(std::stod(str, &nProcessed));
          if (nProcessed != str.size()) {
            return fallback;
          }
          return value;
        } catch (...) {
          return fallback;
        }
      } else {
        static_assert([](){ return false; }, "Unsupported type");
      }
    } else {
      return fallback;
    }
  }

  /// Return all tokens starting with -- that are neither consumed nor a documented option
  std::vector<std::string> getUnconsumed() const {
    std::vector<std::string> unconsumed;
    for (std::size_t i=0; i < _tokens.size(); ++i) {
      if (!_consumed[i] && _tokens[i].size() > 2 && _tokens[i].substr(0, 2) == "--"
          && !isKnown(_tokens[i])) {
        unconsumed.push_back(_tokens[i]);
      }
    }
    return unconsumed;
  }

};

template <>
inline std::string CLIreader::getValueOrFallback<std::string>(
  const std::string& name, std::string fallback) {
  if (contains(name)) {
    return operator[](name);
  } else {
    return fallback;
  }
}

namespace singleton {

/// Process-global CLI token store, populated once by olb::initialize(&argc,&argv)
/// (after MPI init, which may strip argv). The single CLI authority every
/// ParametersD resolves against.
inline CLIreader& cli() {
  static CLIreader instance;
  return instance;
}

}

}

#endif
