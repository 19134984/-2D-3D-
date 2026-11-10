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

#ifndef OLB_UTILITIES_LOG_CHANNEL_H
#define OLB_UTILITIES_LOG_CHANNEL_H

#include "level.h"

#include <atomic>
#include <array>
#include <cstddef>

namespace olb {

namespace logging {

/// Global COMPILE ceiling — nothing more verbose than this survives translation.
#ifndef OLB_LOG_CEILING
#define OLB_LOG_CEILING ::olb::logging::Level::Debug
#endif

inline constexpr std::size_t max_channels = 64;

namespace detail {

/// Registers a channel name once and returns its stable small id. Applies any
/// pending env/API runtime configuration for that name.
std::size_t registerChannel(const char* name, Level runtimeDefault);

/// Per-channel runtime verbosity threshold (mutable at runtime via env/API).
extern std::array<std::atomic<std::uint8_t>, max_channels> runtimeLevel;

} // namespace detail

/// Tier 2 gate: is (Channel, L) retained at compile time? Console logging is
/// always available; the profiler (OLB_SCOPE) is gated separately by
/// FEATURE_PROFILING. Two compile-time kill switches remain:
///   - Channel::compileCeiling (per-channel, set OLB_CHANNEL ceiling to Off to strip)
///   - OLB_LOG_CEILING         (global)
template <typename Channel, Level L>
inline constexpr bool compiled =
     (static_cast<int>(L) <= static_cast<int>(Channel::compileCeiling))
  && (static_cast<int>(L) <= static_cast<int>(OLB_LOG_CEILING));

/// Tier 3 gate: is (Channel, L) enabled right now? Cheap relaxed atomic load.
template <typename Channel>
inline bool enabled(Level l) {
  return static_cast<std::uint8_t>(l)
      <= detail::runtimeLevel[Channel::id()].load(std::memory_order_relaxed);
}

} // namespace logging

} // namespace olb

/// Declare an observability channel: a tag type carrying its compile ceiling,
/// plus a lazily-assigned stable id and registered name + runtime default.
///   OLB_CHANNEL(Comm, Level::Debug, Level::Off);
#define OLB_CHANNEL(NAME, COMPILE_CEILING, RUNTIME_DEFAULT)                    \
  struct NAME {                                                               \
    static constexpr ::olb::logging::Level compileCeiling = (COMPILE_CEILING);    \
    /* inline static data member: initialized once at startup (not lazily),   \
       so id() is a plain const load with no per-call thread-safe guard. */   \
    inline static const std::size_t _olb_log_id =                            \
      ::olb::logging::detail::registerChannel(#NAME, (RUNTIME_DEFAULT));          \
    static std::size_t id() { return _olb_log_id; }                          \
  }

#endif
