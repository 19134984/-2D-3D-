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

#ifndef OLB_UTILITIES_LOG_LOG_H
#define OLB_UTILITIES_LOG_LOG_H

#include "level.h"
#include "channel.h"

#include <functional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace olb {

namespace logging {

namespace detail {

/// The single runtime fan-out point: routes a fully-formed event to the active
/// sinks (stdout writer / Recorder buffer) for its channel. Non-inline so that
/// enabled call sites lower to a plain call and disabled ones to nothing.
void dispatch(std::size_t channel, Level level, std::string&& message);

/// Monotonic seconds since the Recorder epoch.
double getTime();

/// Append a timed scope to the shared Recorder buffer (olb::singleton::profiling()).
void recordScope(std::size_t channel, std::uint64_t nameHash, const char* name,
                 double start, double end);

/// Compile-time FNV-1a (mirrors the profiler's name hashing).
constexpr std::uint64_t fnv1a(const char* s) {
  std::uint64_t h = 0xcbf29ce484222325ULL;
  while (*s) { h ^= static_cast<std::uint8_t>(*s++); h *= 0x100000001b3ULL; }
  return h;
}

} // namespace detail

/// RAII timed scope on a channel. Empty (and fully elided) when the channel is
/// compile-stripped; otherwise records a begin/end pair into the Recorder iff
/// the channel is at Debug+ at runtime.
template <typename Channel>
struct ScopeGuard {
  const char*   name = nullptr;
  std::uint64_t hash = 0;
  double        start = 0.0;
  ScopeGuard(const char* n, std::uint64_t h) {
    if constexpr (compiled<Channel, Level::Debug>) {
      if (enabled<Channel>(Level::Debug)) {
        name = n; hash = h; start = detail::getTime();
      }
    }
  }
  ~ScopeGuard() {
    if constexpr (compiled<Channel, Level::Debug>) {
      if (name) {
        detail::recordScope(Channel::id(), hash, name, start, detail::getTime());
      }
    }
  }
  ScopeGuard(const ScopeGuard&) = delete;
  ScopeGuard& operator=(const ScopeGuard&) = delete;
};

/// Persist the Recorder buffer as a dependency-free raw .olbtrace binary (read
/// by the post-hoc tooling). Lives in libolbcore, so it works on any build.
/// No-op if empty. ADIOS BP5 output is a separate, explicit writer:
/// writePerformanceProfileADIOS2 in io/performanceProfilerWriter.h.
void flushTrace(const std::string& path = "olb_trace");

/// Arrange for flushTrace(pathPrefix) to run automatically on normal program
/// exit and on SIGINT/SIGTERM (so an aborted run still leaves a trace). Idempotent;
/// the flush happens at most once. olb::initialize() wires this up when built with
/// FEATURE_PROFILING.
void enableAutoFlush(const std::string& pathPrefix);

/// Flush the trace now, for hard exits that bypass the atexit handler.
void flushTraceOnExit();

/// Service a pending SIGUSR1 request: if one arrived, dump a numbered trace
/// window to the auto-flush path and clear the buffer. Call from a safe point
/// (e.g. once per step); cheap no-op otherwise.
void servicePendingFlush();

/// Format and route one event. Only ever instantiated for (Channel, Level)
/// pairs that survive the compile gate, and only ever *called* when enabled.
template <typename Channel, Level L, typename... Args>
inline void emit(Args&&... args) {
  std::ostringstream oss;
  (oss << ... << args);
  detail::dispatch(Channel::id(), L, oss.str());
}

/// Map a level word (off/error/warn/info/debug/trace) to its Level
Level parseLevel(std::string_view s);

/// Set the runtime verbosity of a channel ("*" for all)
void setLevel(std::string_view channel, Level l);

/// Where a channel's events go: "stdout"/"cout", "buffer"/"file" (trace only,
/// no console), "both", or "off". Mirror of the OLB_LOG "@route" suffix.
void setRoute(std::string_view channel, std::string_view route);

void setTimestamps(bool on);

/// Whether stdout timestamps are enabled
/// Single source of truth: both olb::logging and OstreamManager honour it.
bool timestamps();

/// Set the "should this rank/thread print to stdout?" predicate. Keeps log.cpp
/// MPI-agnostic; olb::initialize() wires this to singleton::mpi().isMainProcessor().
void setMainPredicate(std::function<bool()> pred);

/// Registered channel names indexed by channel id (for labelling traces).
std::vector<std::string> channelNames();

/// Current runtime level of each channel, indexed by channel id
std::vector<Level> channelLevels();

} // namespace logging

} // namespace olb

// -- Emission macros. Console logging (OLB_LOG/OLB_LOGF) is always compiled and
// -- gated by the per-channel compile ceiling + runtime level, so --log works in
// -- every build. Profiler timing spans (OLB_SCOPE) need the Recorder and stay
// -- behind FEATURE_PROFILING.
#define OLB_LOG(CH, LVL, ...)                                                  \
  do {                                                                        \
    if constexpr (::olb::logging::compiled<CH, ::olb::logging::Level::LVL>) {          \
      if (::olb::logging::enabled<CH>(::olb::logging::Level::LVL)) {                   \
        ::olb::logging::emit<CH, ::olb::logging::Level::LVL>(__VA_ARGS__);             \
      }                                                                       \
    }                                                                         \
  } while (0)

/// Deferred payload: the lambda body (and its captures' work) is confined to
/// the gated branch, so it never runs when the channel/level is stripped/off.
#define OLB_LOGF(CH, LVL, FN)                                                  \
  do {                                                                        \
    if constexpr (::olb::logging::compiled<CH, ::olb::logging::Level::LVL>) {          \
      if (::olb::logging::enabled<CH>(::olb::logging::Level::LVL)) {                   \
        ::olb::logging::emit<CH, ::olb::logging::Level::LVL>((FN)());                  \
      }                                                                       \
    }                                                                         \
  } while (0)

#define OLB_LOG_CAT_(a, b) a##b
#define OLB_LOG_CAT(a, b)  OLB_LOG_CAT_(a, b)

/// RAII timed scope (name hashed at compile time), recorded into the profiler.
#ifdef FEATURE_PROFILING
#define OLB_SCOPE(CH, NAME)                                                    \
  ::olb::logging::ScopeGuard<CH> OLB_LOG_CAT(_olb_scope_, __LINE__)               \
    { NAME, ::olb::logging::detail::fnv1a(NAME) }
#else
#define OLB_SCOPE(CH, NAME)
#endif

#endif
