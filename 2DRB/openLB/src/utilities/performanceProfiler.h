/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Adrian Kummerlaender
 *
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

#ifndef PERFORMANCE_PROFILER_H
#define PERFORMANCE_PROFILER_H

#include <vector>
#include <chrono>
#include <string>
#include <string_view>
#include <cstdint>
#include <map>
#include <mutex>
#include <algorithm>

#include "core/meta.h"

namespace olb {

namespace logging {

class Recorder {
public:
  enum class State : std::uint8_t {
    Disabled = 0,
    Enabled = 1
  };

  /// Channel id stored on records that originate from OLB_PROFILE (rather than
  /// from a named olb::logging channel). Lanes profiler scopes in the trace.
  static constexpr std::uint8_t PROFILE_CHANNEL = 0xFF;
  /// Dedicated channel for per-timestep OLB_PROFILE_STEP scopes, so inspection
  /// tooling can aggregate (mean step time) and exclude them from scope rankings.
  static constexpr std::uint8_t STEP_CHANNEL = 0xFE;

  struct Record {
    std::uint8_t  channel;
    std::uint64_t nameHash;
    double start;
    double end;
  };

  struct MessageRecord {
    std::uint8_t channel;
    double timestamp;
    std::string message;
  };

private:
  std::vector<Record> _buffer;
  std::vector<MessageRecord> _messages;
  std::chrono::steady_clock::time_point _epoch;
  std::map<std::uint64_t, std::string> _nameMap;
  std::mutex _mapMutex;
  std::size_t _stepCounter = 0;
  State _state = State::Disabled;

  // FNV-1a hash (64-bit)
  static constexpr std::uint64_t fnv1a(const char* str) {
    std::uint64_t hash = 0xcbf29ce484222325ULL;
    while (*str != '\0') {
      hash ^= static_cast<std::uint64_t>(*str++);
      hash *= 0x100000001b3ULL;
    }
    return hash;
  }

public:
  Recorder() {
    _epoch = std::chrono::steady_clock::now();
    _buffer.reserve(2000000);
  }

  void setRecording(bool state) {
    _state = state ? State::Enabled : State::Disabled;
  }

  bool isRecording() const {
    return _state == State::Enabled;
  }

  // OLB_PROFILE records are State-gated (explicit setRecording). olb::logging
  // events pass force=true: they have already passed their own channel/level
  // gate and record into the shared buffer regardless of the profiler State.
  void record(std::uint64_t nameHash, const char* name, double start, double end,
              std::uint8_t channel = PROFILE_CHANNEL, bool force = false) {
    if (!force && _state != State::Enabled) {
      return;
    }

    std::lock_guard<std::mutex> lock(_mapMutex);
    _buffer.push_back({channel, nameHash, start, end});
    if (_nameMap.find(nameHash) == _nameMap.end()) {
      _nameMap[nameHash] = std::string(name);
    }
  }

  void recordMessage(const char* message, std::uint8_t channel = PROFILE_CHANNEL,
                     bool force = false) {
    if (!force && _state != State::Enabled) {
      return;
    }

    std::lock_guard<std::mutex> lock(_mapMutex);
    _messages.push_back({channel, getTime(), std::string(message)});
  }

  void resetTimer() {
    _epoch = std::chrono::steady_clock::now();
  }

  std::size_t incrementStepCounter() {
    return _stepCounter++;
  }

  double getTime() const {
    std::chrono::duration<double> dist = std::chrono::steady_clock::now() - _epoch;
    return dist.count();
  }

  const std::vector<Record>& getBuffer() const { return _buffer; }
  const std::vector<MessageRecord>& getMessages() const { return _messages; }
  const std::map<std::uint64_t, std::string>& getNameMap() const { return _nameMap; }

  void clear() {
    _buffer.clear();
    _messages.clear();
  }

  static constexpr std::uint64_t hash(const char* name) {
    return fnv1a(name);
  }

  static std::uint64_t dynamic_hash(const std::string& name) {
    return fnv1a(name.c_str());
  }
};

} // namespace logging

namespace singleton {

logging::Recorder& profiling();

} // namespace singleton

namespace logging {

template <std::uint64_t NAME_HASH, std::uint8_t CHANNEL = Recorder::PROFILE_CHANNEL>
struct ScopedProfiler {
  const char* _name;
  double _start = 0.0;
  bool _active;

  ScopedProfiler(const char* name) : _name(name) {
    _active = singleton::profiling().isRecording();
    if (_active) {
      _start = singleton::profiling().getTime();
    }
  }
  ~ScopedProfiler() {
    if (_active) {
      singleton::profiling().record(NAME_HASH, _name, _start, singleton::profiling().getTime(), CHANNEL);
    }
  }
};

struct DynamicScopedProfiler {
  std::string _name;
  std::uint64_t _hash;
  double _start = 0.0;
  bool _active;

  DynamicScopedProfiler(std::string name) : _name(name) {
    _active = singleton::profiling().isRecording();
    if (_active) {
      _hash = Recorder::dynamic_hash(_name);
      _start = singleton::profiling().getTime();
    }
  }

  ~DynamicScopedProfiler() {
    if (_active) {
      singleton::profiling().record(_hash, _name.c_str(), _start, singleton::profiling().getTime());
    }
  }
};

#ifdef FEATURE_PROFILING
#define OLB_PROFILE(NAME) \
  static constexpr std::uint64_t _olb_prof_hash_##__LINE__ = olb::logging::Recorder::hash(NAME); \
  olb::logging::ScopedProfiler<_olb_prof_hash_##__LINE__> _olb_prof_timer_##__LINE__(NAME);

#define OLB_PROFILE_EVENT(NAME) \
  { \
    static constexpr std::uint64_t _olb_prof_event_hash_##__LINE__ = olb::logging::Recorder::hash(NAME); \
    if (olb::singleton::profiling().isRecording()) { \
      const double _olb_prof_event_time = olb::singleton::profiling().getTime(); \
      olb::singleton::profiling().record(_olb_prof_event_hash_##__LINE__, NAME, _olb_prof_event_time, _olb_prof_event_time); \
    } \
  }

template <std::size_t N>
struct StaticString {
  char data[N]{};
  static constexpr std::size_t size = N;
};

template <typename... Ts>
struct ProfileNameSize {
  static constexpr std::size_t value = (0 + ... + (1 + meta::nice_name<Ts>().size()));
};

template <std::size_t PrefixSize, typename... Ts>
constexpr auto build_profile_name(std::string_view prefix) {
  constexpr std::size_t total_size = PrefixSize + ProfileNameSize<Ts...>::value;
  StaticString<total_size + 1> result;
  std::size_t pos = 0;
  for (char c : prefix) {
    result.data[pos++] = c;
  }
  constexpr std::string_view names[] = { meta::nice_name<Ts>()... };
  for (auto s : names) {
    result.data[pos++] = '|';
    for (char c : s) {
      result.data[pos++] = c;
    }
  }
  result.data[pos] = '\0';
  return result;
}

#define OLB_PROFILE_T(NAME, ...) \
  static constexpr auto _olb_prof_name_##__LINE__ = olb::logging::build_profile_name<sizeof(NAME)-1, __VA_ARGS__>(NAME); \
  static constexpr std::uint64_t _olb_prof_hash_##__LINE__ = olb::logging::Recorder::hash(_olb_prof_name_##__LINE__.data); \
  olb::logging::ScopedProfiler<_olb_prof_hash_##__LINE__> _olb_prof_timer_##__LINE__(_olb_prof_name_##__LINE__.data);

template <typename... Ts>
std::string build_dynamic_profile_name(std::string action) {
  constexpr std::string_view names[] = { meta::nice_name<Ts>()... };
  for (auto s : names) {
    action += "|";
    action += s;
  }
  return action;
}

#define OLB_PROFILE_DYNAMIC_T(NAME, ...) \
  olb::logging::DynamicScopedProfiler _olb_prof_timer_##__LINE__(olb::logging::build_dynamic_profile_name<__VA_ARGS__>(NAME));

// Per-timestep scope on the dedicated STEP_CHANNEL with a constant name, so the
// steps aggregate (count / mean step time) and tooling can exclude them from the
// scope ranking (instead of one unique "Step N" name per step).
#define OLB_PROFILE_STEP() \
  static constexpr std::uint64_t _olb_prof_step_hash_##__LINE__ = olb::logging::Recorder::hash("Step"); \
  olb::logging::ScopedProfiler<_olb_prof_step_hash_##__LINE__, olb::logging::Recorder::STEP_CHANNEL> _olb_prof_step_timer_##__LINE__("Step");

#else
#define OLB_PROFILE(NAME)
#define OLB_PROFILE_EVENT(NAME)
#define OLB_PROFILE_T(NAME, ...)
#define OLB_PROFILE_DYNAMIC_T(NAME, ...)
#define OLB_PROFILE_STEP()
#endif

} // namespace logging

} // namespace olb

#endif
