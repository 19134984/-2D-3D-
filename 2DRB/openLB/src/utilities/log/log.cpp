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

#include "log.h"

#include "utilities/performanceProfiler.h"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <iostream>
#include <map>
#include <mutex>

#if defined(__unix__) || defined(__linux__) || defined(__APPLE__)
#include <signal.h>     // sigaction: block interrupts while the trace flushes
#include <unistd.h>     // write(): async-signal-safe notice
#endif
#include <string>
#include <vector>

namespace olb {

namespace logging {

namespace {
std::string toLower(std::string_view s) {
  std::string r(s);
  std::transform(r.begin(), r.end(), r.begin(),
                 [](unsigned char c){ return std::tolower(c); });
  return r;
}
} // namespace

Level parseLevel(std::string_view s) {
  const std::string l = toLower(s);
  if (l.starts_with("deb")) { return Level::Debug; }
  return Level::Off;
}

namespace {

enum class Route : std::uint8_t { Off = 0, Stdout = 1, Buffer = 2, Both = 3 };

struct Config {
  std::mutex mutex;
  std::map<std::string, Level> pendingLevel;   // by name, incl. "*" wildcard
  std::map<std::string, Route> pendingRoute;   // by name, incl. "*" wildcard
  std::vector<std::string>     names;          // id -> name
  std::vector<Route>           route;          // id -> route
  bool timestamps = false;
  std::function<bool()> isMain = []{ return true; }; // set to MPI rank-0 on integration
};

Route parseRoute(std::string_view s) {
  const std::string r = toLower(s);
  if (r == "buffer" || r == "file" || r == "trace") { return Route::Buffer; }
  if (r == "both")                                  { return Route::Both;   }
  if (r == "off"  || r == "none"  || r == "mute")   { return Route::Off;    }
  return Route::Stdout;                               // "stdout"/"cout"/"console"/default
}

Config& cfg() { static Config c; return c; }

} // namespace

namespace detail {

std::array<std::atomic<std::uint8_t>, max_channels> runtimeLevel{};

std::size_t registerChannel(const char* name, Level runtimeDefault) {
  auto& c = cfg();
  std::lock_guard<std::mutex> lock(c.mutex);
  const std::size_t id = c.names.size();
  c.names.emplace_back(name);

  // Effective config: explicit name > wildcard > default (case-insensitive).
  const std::string key = toLower(name);
  Level lvl = runtimeDefault;
  if (auto it = c.pendingLevel.find("*"); it != c.pendingLevel.end()) { lvl = it->second; }
  if (auto it = c.pendingLevel.find(key); it != c.pendingLevel.end()) { lvl = it->second; }
  runtimeLevel[id].store(static_cast<std::uint8_t>(lvl), std::memory_order_relaxed);

  Route rt = Route::Stdout;
  if (auto it = c.pendingRoute.find("*"); it != c.pendingRoute.end()) { rt = it->second; }
  if (auto it = c.pendingRoute.find(key); it != c.pendingRoute.end()) { rt = it->second; }
  c.route.push_back(rt);
  return id;
}

double getTime() { return singleton::profiling().getTime(); }

void recordScope(std::size_t channel, std::uint64_t nameHash, const char* name,
                 double start, double end) {
  singleton::profiling().record(nameHash, name, start, end,
                                static_cast<std::uint8_t>(channel), /*force=*/true);
}

void dispatch(std::size_t channel, [[maybe_unused]] Level level, std::string&& message) {
  auto& c = cfg();
  const Route r = (channel < c.route.size()) ? c.route[channel] : Route::Stdout;

  if (r == Route::Stdout || r == Route::Both) {
    if (c.isMain()) {
      const char* chName = (channel < c.names.size()) ? c.names[channel].c_str() : "?";
      std::lock_guard<std::mutex> lock(c.mutex);           // single emission (OMP-safe)
      if (c.timestamps) {                                  // [channel:date] — mirrors OstreamManager
        std::time_t t = std::time(nullptr);
        char buf[20];
        std::strftime(buf, sizeof buf, "%Y-%m-%d %H:%M:%S", std::localtime(&t));
        std::cout << '[' << chName << ':' << buf << "] " << message << '\n';
      } else {
        std::cout << '[' << chName << "] " << message << '\n';
      }
    }
  }
  if (r == Route::Buffer || r == Route::Both) {
    singleton::profiling().recordMessage(message.c_str(),
                                         static_cast<std::uint8_t>(channel), /*force=*/true);
  }
}

} // namespace detail

void setLevel(std::string_view channel, Level l) {
  auto& c = cfg();
  std::lock_guard<std::mutex> lock(c.mutex);
  std::string key = toLower(channel);
  c.pendingLevel[key] = l;                                  // for channels not yet registered
  for (std::size_t id = 0; id < c.names.size(); ++id) {     // and those already registered
    if (key == "*" || key == toLower(c.names[id])) {
      detail::runtimeLevel[id].store(static_cast<std::uint8_t>(l), std::memory_order_relaxed);
    }
  }
}

void setRoute(std::string_view channel, std::string_view route) {
  auto& c = cfg();
  std::lock_guard<std::mutex> lock(c.mutex);
  std::string key = toLower(channel);
  Route rt = parseRoute(route);
  c.pendingRoute[key] = rt;                                 // for channels not yet registered
  for (std::size_t id = 0; id < c.names.size(); ++id) {     // and those already registered
    if (key == "*" || key == toLower(c.names[id])) {
      c.route[id] = rt;
    }
  }
}

void setTimestamps(bool on) { cfg().timestamps = on; }
bool timestamps() { return cfg().timestamps; }

void setMainPredicate(std::function<bool()> pred) {
  auto& c = cfg();
  std::lock_guard<std::mutex> lock(c.mutex);
  c.isMain = std::move(pred);
}

std::vector<std::string> channelNames() {
  auto& c = cfg();
  std::lock_guard<std::mutex> lock(c.mutex);
  return c.names;
}

std::vector<Level> channelLevels() {
  auto& c = cfg();
  std::lock_guard<std::mutex> lock(c.mutex);
  std::vector<Level> levels;
  for (std::size_t id = 0; id < c.names.size(); ++id) {
    levels.push_back(static_cast<Level>(detail::runtimeLevel[id].load(std::memory_order_relaxed)));
  }
  return levels;
}

// Dependency-free raw trace dump as ADIOS alternative
void flushTrace(const std::string& path) {
  auto& rec = singleton::profiling();
  if (rec.getBuffer().empty() && rec.getMessages().empty()) {
    return;
  }

  const auto& buffer   = rec.getBuffer();
  const auto& messages = rec.getMessages();
  const auto& nameMap  = rec.getNameMap();
  // raw binary (SoA), read by script/profiling/olbprof.py:
  //   magic[8]="OLBLOG\3",
  //   channels{u32 n, (u32 len,char[len])* }   (name; id = index),
  //   names{u32 n, (u64 hash,u32 len,char[len])*},
  //   scopes{u32 n, (u8 ch,u64 hash,f64 start,f64 end)*},
  //   messages{u32 n, (u8 ch,f64 ts,u32 len,char[len])*}
  std::FILE* f = std::fopen((path + ".olbtrace").c_str(), "wb");
  if (!f) {
    return;
  }
  auto put = [&](const void* p, std::size_t n){ std::fwrite(p, 1, n, f); };

  put("OLBLOG\3", 8);
  const auto chNames = channelNames();
  const std::uint32_t nC = static_cast<std::uint32_t>(chNames.size());
  put(&nC, 4);
  for (const auto& name : chNames) {
    const std::uint32_t len = static_cast<std::uint32_t>(name.size());
    put(&len, 4); put(name.data(), len);
  }
  const std::uint32_t nN = static_cast<std::uint32_t>(nameMap.size());
  put(&nN, 4);
  for (const auto& [hash, name] : nameMap) {
    const std::uint32_t len = static_cast<std::uint32_t>(name.size());
    put(&hash, 8); put(&len, 4); put(name.data(), len);
  }
  const std::uint32_t nS = static_cast<std::uint32_t>(buffer.size());
  put(&nS, 4);
  for (const auto& s : buffer) {
    put(&s.channel, 1); put(&s.nameHash, 8); put(&s.start, 8); put(&s.end, 8);
  }
  const std::uint32_t nM = static_cast<std::uint32_t>(messages.size());
  put(&nM, 4);
  for (const auto& m : messages) {
    const std::uint32_t len = static_cast<std::uint32_t>(m.message.size());
    put(&m.channel, 1); put(&m.timestamp, 8); put(&len, 4); put(m.message.data(), len);
  }
  std::fclose(f);
}

namespace {

std::string autoFlushPath;
std::atomic<bool> autoFlushArmed{false};
std::atomic<bool> usr1Requested{false};   // SIGUSR1 -> on-demand window dump
std::atomic<int>  flushWindow{0};
void onSignalUsr1(int) { usr1Requested.store(true, std::memory_order_relaxed); }
void doAutoFlush() {
  bool expected = true;
  if (autoFlushArmed.compare_exchange_strong(expected, false)) {
    flushTrace(autoFlushPath);                   // at most once (exit OR signal)
  }
}
void onSignal(int sig) {
  // SIGINT/SIGTERM are masked for the whole handler (see sa_mask below), so a
  // second, impatient Ctrl+C is held pending and cannot truncate the trace
  // mid-write. A short async-signal-safe notice tells the user to wait.
#if defined(__unix__) || defined(__linux__) || defined(__APPLE__)
  static const char msg[] = "\n[logging] interrupted — flushing trace, please wait...\n";
  [[maybe_unused]] ssize_t n = ::write(STDERR_FILENO, msg, sizeof(msg) - 1);
#endif
  doAutoFlush();
  std::signal(sig, SIG_DFL);                      // restore + re-raise so abort still terminates
  std::raise(sig);                                // pending until the mask lifts on return → dies
}

} // namespace

void flushTraceOnExit() { doAutoFlush(); }

void enableAutoFlush(const std::string& pathPrefix) {
  autoFlushPath = pathPrefix;
  autoFlushArmed.store(true);
  std::atexit(doAutoFlush);
#if defined(__unix__) || defined(__linux__) || defined(__APPLE__)
  // sigaction (not std::signal) so we control masking: block both interrupts for
  // the duration of the handler so the flush always completes.
  struct sigaction sa{};
  sa.sa_handler = onSignal;
  sigemptyset(&sa.sa_mask);
  sigaddset(&sa.sa_mask, SIGINT);
  sigaddset(&sa.sa_mask, SIGTERM);
  sa.sa_flags = 0;                                // no SA_NODEFER: handled signal stays masked too
  sigaction(SIGINT,  &sa, nullptr);
  sigaction(SIGTERM, &sa, nullptr);
  struct sigaction su{};
  su.sa_handler = onSignalUsr1;
  sigemptyset(&su.sa_mask);
  su.sa_flags = SA_RESTART;
  sigaction(SIGUSR1, &su, nullptr);
#else
  std::signal(SIGINT,  onSignal);
  std::signal(SIGTERM, onSignal);
  std::signal(SIGUSR1, onSignalUsr1);
#endif
}

void servicePendingFlush() {
  bool expected = true;
  if (usr1Requested.compare_exchange_strong(expected, false)) {
    // Dump a numbered window and reset the buffer so the next window starts clean
    flushTrace(autoFlushPath + "_win" + std::to_string(flushWindow.fetch_add(1) + 1));
    singleton::profiling().clear();
  }
}

} // namespace logging

} // namespace olb
