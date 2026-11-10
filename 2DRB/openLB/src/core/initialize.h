/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Adrian Kummerlaender
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

#ifndef CORE_INITIALIZE_H
#define CORE_INITIALIZE_H

#include "olbInit.h"
#include "utilities/log/log.h"
#include "utilities/performanceProfiler.h"
#include "io/input.h"
#include "io/textTable.h"

#include <cstdio>
#include <unistd.h>
#ifdef PARALLEL_MODE_OMP
#include <omp.h>
#include <sched.h>
#include <fstream>
#include <set>
#include <utility>
#include <string>
#include <cstdlib>
#endif

namespace olb {

/// Resolve output presentation (text style, stdout timestamps) through the input authority
inline void configureOutput()
{
  using Query = OverrideSource::Query;
  auto& input = singleton::input();
  auto& style = io::text::style();

  if (auto mode = input.resolve(Query{{}, "text_color", true, false})) {
    style.color = (*mode == "always") || (*mode == "auto" && isatty(fileno(stdout)));
  }
  if (input.resolveBuiltin<parameters::text_ascii>(false)) {
    style.unicode = false;
  }
  if (auto value = input.resolve(Query{{}, "text_style", true, false})) {
    if (*value == "labeled")   style.fields = io::text::Style::Fields::Labeled;
    if (*value == "keyvalue")  style.fields = io::text::Style::Fields::KeyValue;
    if (*value == "semicolon") style.fields = io::text::Style::Fields::Semicolon;
  }
  if (auto value = input.resolve(Query{{}, "table_style", true, false})) {
    if (auto layout = io::text::parseTableLayout(*value)) {
      style.table = *layout;
    }
  }
  if (input.resolveBuiltin<parameters::clout_timestamps>(false)) {
    logging::setTimestamps(true);
  }
}

/// Resolve global and per-channel log levels through the input authority
inline void configureLogging()
{
  using Query = OverrideSource::Query;
  auto& input = singleton::input();
  if (auto level = input.resolve(Query{{}, "log", false, false})) {
    logging::setLevel("*", logging::parseLevel(*level));
  }
  for (const auto& channel : logging::channelNames()) {
    if (auto level = input.resolve(Query{{"log"}, channel, false, false})) {
      logging::setLevel(channel, logging::parseLevel(*level));
    }
  }
}

/// Initialize OpenLB
/**
 * Sets up MPI, thread pool and verifies platform requirements.
 **/
#ifdef PARALLEL_MODE_OMP
inline int physicalCoreCount()
{
  std::set<std::pair<int,int>> cores;
  const long n = sysconf(_SC_NPROCESSORS_ONLN);
  for (long c = 0; c < n; ++c) {
    auto read = [c](const char* field) {
      std::ifstream in("/sys/devices/system/cpu/cpu" + std::to_string(c) + "/topology/" + field);
      int v = -1; in >> v; return v;
    };
    const int pkg = read("physical_package_id");
    const int core = read("core_id");
    if (pkg >= 0 && core >= 0) cores.insert({pkg, core});
  }
  return cores.empty() ? static_cast<int>(n) : static_cast<int>(cores.size());
}

inline void autoAffinity(bool verbose)
{
  if (std::getenv("OMP_NUM_THREADS") || std::getenv("OMP_PLACES")) {
    return;
  }
  const long n = sysconf(_SC_NPROCESSORS_CONF);
  cpu_set_t set; CPU_ZERO(&set);
  for (long c = 0; c < n && c < CPU_SETSIZE; ++c) {
    CPU_SET(c, &set);
  }
  sched_setaffinity(0, sizeof(set), &set);
  const int cores = physicalCoreCount();
  setenv("OMP_WAIT_POLICY", "passive", 0);
  omp_set_num_threads(cores);
  if (verbose) {
    OstreamManager clout(std::cout, "initialize");
    clout << "Auto-affinity: unpinned rank over its cores, OMP capped to "
          << cores << " physical" << std::endl;
  }
}
#endif

void initialize(int *argc, char ***argv, bool multiOutput=false, bool verbose=true)
{
  std::set_terminate(terminateHandler);

  singleton::input().load(*argc, *argv);

#if defined(PARALLEL_MODE_OMP) && (defined(PLATFORM_GPU_CUDA) || defined(PLATFORM_GPU_HIP) || defined(PLATFORM_GPU_SYCL))
  if (!std::getenv("OLB_REEXECED") && !std::getenv("OMP_WAIT_POLICY")
      && !std::getenv("OMPI_COMM_WORLD_SIZE") && !std::getenv("PMIX_RANK")
      && singleton::input().resolveBuiltin<parameters::omp_affinity>(true)) {
    setenv("OMP_WAIT_POLICY", "passive", 1);
    setenv("OLB_REEXECED", "1", 1);
    execvp((*argv)[0], *argv);
  }
#endif

#ifdef PARALLEL_MODE_MPI
  bool mpiInitialized = false;
  if (singleton::mpi().init(argc, argv)) {
    mpiInitialized = true;
  }
#endif

  // Set up the process-global input authority after MPI init (which may
  // strip/reorder argv): capture the CLI, load an optional --case-file XML, and
  // handle --help. All parameter resolution reads from singleton::input().
  singleton::input().load(*argc, *argv);
  configureOutput();

  // Observability: route stdout to the main rank only, then apply log config.
  logging::setMainPredicate([]{ return singleton::mpi().isMainProcessor(); });
  configureLogging();

#ifdef FEATURE_PROFILING
  // Off by default; enable with --profile (or OLB_PROFILE=1). The trace dumps to
  // profilingData/ on normal exit and on abort (SIGINT/SIGTERM).
  if (singleton::input().resolveBuiltin<parameters::profile>(false)) {
    singleton::profiling().setRecording(true);
    logging::enableAutoFlush(singleton::directories().getProfilingOutDir()
                             + "trace_rank" + std::to_string(singleton::mpi().getRank()));
  }
#endif

#ifndef FEATURE_HIDE_LOGO
  // create an OstreamManager object in order to enable multi output
  if (singleton::mpi().isMainProcessor()) {
  std::cout << R"(
    ┃
    ┃  ┏━━━━┓      ▁▁▁▁                   ▁▁    ▁▁▁▁
    ┃  ┃    ┃     ╱ ▁▁ ╲▁▁▁▁  ▁▁▁  ▁▁▁▁  ╱ ╱   ╱ ▁▁ ╲
 ┏━━╋━┓┃    ┃    ╱ ╱ ╱ ╱ ▁▁ ╲╱ ▁ ╲╱ ▁▁ ╲╱ ╱   ╱ ╱▁╱ ╱
 ┃  ┗━╋╋━━━━┻┓  ╱ ╱▁╱ ╱ ╱▁╱ ╱  ▁▁╱ ╱ ╱ ╱ ╱▁▁▁╱ ╱▁╱ ╱
 ┃    ┃┃     ┃  ╲▁▁▁▁╱ ▁▁▁▁╱╲▁▁▁╱▁╱ ╱▁╱▁▁▁▁▁╱▁▁▁▁▁╱
 ┗━━━━┛┃     ┃      ╱▁╱ ==========================>>
       ┗━━━━━┛
)" << std::endl;
  }
#endif

  olb::OstreamManager clout(std::cout, "initialize");

  std::string platform;
#ifdef PLATFORM_CPU_SISD
  platform += " CPU_SISD";
#endif
#ifdef PLATFORM_CPU_SIMD
  platform += " CPU_SIMD";
#endif
#ifdef PLATFORM_GPU_CUDA
  platform += " GPU_CUDA";
#endif
#ifdef PLATFORM_GPU_HIP
  platform += " GPU_HIP";
#endif
#ifdef PLATFORM_GPU_SYCL
  platform += " GPU_SYCL";
#endif
  std::string parallel;
#ifdef PARALLEL_MODE_MPI
  parallel += " MPI";
#endif
#ifdef PARALLEL_MODE_OMP
  parallel += " OMP";
#endif

  io::text::KeyValueList info;
  info.item("Version", std::string{OLB_VERSION});
  info.item("Platform", platform.empty() ? std::string{"-"} : platform.substr(1));
  info.item("Parallel", parallel.empty() ? std::string{"None"} : parallel.substr(1));
#ifdef PARALLEL_MODE_MPI
  if (mpiInitialized) {
    info.item("MPI", std::to_string(singleton::mpi().getSize()) + " ranks");
    if (singleton::mpi().getNodes() > 1) {
      info.item("Nodes", std::to_string(singleton::mpi().getNodes()));
    }
  } else {
    info.item("MPI", std::string{"Failed!"});
  }
#endif
#ifdef PARALLEL_MODE_OMP
  if (singleton::input().resolveBuiltin<parameters::omp_affinity>(true)) {
    autoAffinity(verbose);
  }
  if (singleton::omp().init(argc, argv)) {
    const char* waitPolicy = std::getenv("OMP_WAIT_POLICY");
    info.item("OMP", std::to_string(singleton::omp().getSize()) + " threads, "
                     + (waitPolicy ? waitPolicy : "default") + " wait, affinity "
                     + (singleton::input().resolveBuiltin<parameters::omp_affinity>(true)
                        ? "managed" : "inherited"));
  } else {
    info.item("OMP", std::string{"Failed!"});
  }
#endif
  info.render(clout, "");

  if (singleton::mpi().isMainProcessor()
      && !singleton::input().resolveBuiltin<parameters::no_introspection>(false)) {
    std::ofstream fout(singleton::directories().getIntrospectionOutDir() + "builtin.log",
                       std::ios::trunc);
    singleton::input().print(fout);
  }

  clout.setMultiOutput(multiOutput);

  /// Verify requirements for using all enabled platforms
#ifdef PLATFORM_GPU_CUDA
  checkPlatform<Platform::GPU_CUDA>();
#endif
#ifdef PLATFORM_GPU_HIP
  checkPlatform<Platform::GPU_HIP>();
#endif
#ifdef PLATFORM_GPU_SYCL
  checkPlatform<Platform::GPU_SYCL>();
#endif

  int nThreads = 1;
  if (const char* envOlbNumThreads = std::getenv("OLB_NUM_THREADS")) {
    nThreads = std::stoi(envOlbNumThreads);
  }
  singleton::pool().init(nThreads, verbose);

  #ifdef FEATURE_VDB
  openvdb::initialize();
  #endif
}

/// Initialize OpenLB
void initialize(int argc, char **argv, bool multiOutput=false, bool verbose=true)
{
  initialize(&argc, &argv, multiOutput, verbose);
}

}

#endif
