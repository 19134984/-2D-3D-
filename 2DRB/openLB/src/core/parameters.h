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

#ifndef CORE_PARAMETERS_H
#define CORE_PARAMETERS_H

// Framework built-in parameter tags and naming helpers; physics parameters live in case/parameters.h

#include "core/meta.h"

#include <string>
#include <string_view>
#include <cctype>

namespace olb {

namespace parameters {

/// Returns name of PARAMETER for human consumption (namespace-qualifier stripped)
template <typename PARAMETER>
std::string name() {
  std::string_view raw = meta::nice_name<PARAMETER>();
  // Look backwards for the last occurrence of "::"
  std::size_t pos = raw.rfind("::");

  if (pos != std::string_view::npos) {
    // Keep everything to the right of "::"
    return std::string(raw.substr(pos + 2));
  }

  // If "::" wasn't found (no namespace), just return the raw string
  return std::string(raw);
}

/// Naming convention: native (runner) options are lowercase, user physics
/// parameters uppercase. Used to derive the CLI token kind (hyphenated vs as-is).
inline bool isNativeName(std::string_view name) {
  for (char c : name) {
    if (std::isupper(static_cast<unsigned char>(c))) {
      return false;
    }
  }
  return true;
}

// ----- Framework built-in parameters (runner options + framework-read fields) -----

struct OVERLAP : public descriptors::TYPED_FIELD_BASE<unsigned,1> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,1>{3};
  }
};

/// Skip the per-component introspection dump written at lattice initialization (native flag --no-introspection)
struct no_introspection : public descriptors::TYPED_FIELD_BASE<bool,1> { };

/// Use ASCII instead of box-drawing glyphs (native flag --text-ascii)
struct text_ascii : public descriptors::TYPED_FIELD_BASE<bool,1> { };

/// Prefix stdout lines with a timestamp (native flag --clout-timestamps)
struct clout_timestamps : public descriptors::TYPED_FIELD_BASE<bool,1> { };

/// Record a performance trace (native flag --profile)
struct profile : public descriptors::TYPED_FIELD_BASE<bool,1> { };

/// Auto-manage OMP thread count and affinity for hybrid builds (native flag --omp_affinity)
struct omp_affinity : public descriptors::TYPED_FIELD_BASE<bool,1> { };

/// Enable a parameter group; under "reporter" (--reporter.enabled 0) it turns off all reporters
struct enabled : public descriptors::TYPED_FIELD_BASE<bool,1> {
  template <typename T, typename DESCRIPTOR>
  static constexpr auto getInitialValue() {
    return Vector<value_type<T>,1>{true};
  }
};

}

}

#endif
