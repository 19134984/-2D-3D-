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

#ifndef OLB_UTILITIES_LOG_LEVEL_H
#define OLB_UTILITIES_LOG_LEVEL_H

#include <cstdint>
#include <string_view>

namespace olb {

namespace logging {

/// Severity of an observability event, ordered by increasing verbosity.
enum class Level : std::uint8_t {
  Off = 0, Debug
};

constexpr std::string_view name(Level l) {
  return l == Level::Debug ? "DEBUG" : "OFF";
}

/// Parse a level name (case-insensitive prefix accepted: e.g. "deb" -> Debug).
Level parseLevel(std::string_view s);

} // namespace logging

} // namespace olb

#endif
