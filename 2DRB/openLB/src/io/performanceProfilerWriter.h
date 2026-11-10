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

#ifndef PERFORMANCE_PROFILER_WRITER_H
#define PERFORMANCE_PROFILER_WRITER_H

#include "utilities/performanceProfiler.h"
#include "utilities/log/log.h"
#include "core/singleton.h"

#include <string>
#include <vector>

#ifdef FEATURE_ADIOS2
#include <adios2.h>
#endif

namespace olb {

namespace logging {

inline void writePerformanceProfileADIOS2(const std::string& filename,
                                          Recorder& profiler) {
#ifdef FEATURE_ADIOS2
#ifdef PARALLEL_MODE_MPI
    static MPI_Comm profilerComm = []() {
      MPI_Comm comm;
      MPI_Comm_dup(MPI_COMM_WORLD, &comm);
      return comm;
    }();
    static adios2::ADIOS adios(profilerComm);
#else
    static adios2::ADIOS adios;
#endif
    static adios2::IO io;
    static adios2::Engine engine;
    static bool initialized = false;

    if (!initialized) {
      io = adios.DeclareIO("ProfilerIO");
      io.SetEngine("BP5");
      std::string path = singleton::directories().getProfilingOutDir() + filename;
      engine = io.Open(path, adios2::Mode::Write);
      // The engine streams one ADIOS step per call and stays open across calls;
      // BP5 only finalizes its metadata index on Close(), so without this the
      // file is left incomplete and unreadable at program exit. Close the
      // (function-static) engine exactly once, at static destruction.
      static struct EngineCloser {
        adios2::Engine* e;
        ~EngineCloser() {
          if (e && *e) {
            e->Close();
          }
        }
      } closer{&engine};
      initialized = true;
    }

    const auto& buffer = profiler.getBuffer();
    const auto& nameMap = profiler.getNameMap();
    const std::size_t nEvents = buffer.size();
    const std::size_t rank = static_cast<std::size_t>(singleton::mpi().getRank());
    const std::size_t size = static_cast<std::size_t>(singleton::mpi().getSize());

    // 1. Metadata: Mapping hashes to original strings
    for (auto const& [hash, name] : nameMap) {
      std::string attrName = "name_map/" + std::to_string(hash);
      if (!io.InquireAttribute<std::string>(attrName)) {
        io.DefineAttribute<std::string>(attrName, name);
      }
    }
    // Channel id -> name, so inspection tooling can label channels.
    const auto chNames = channelNames();
    for (std::size_t id = 0; id < chNames.size(); ++id) {
      std::string attrName = "channel_map/" + std::to_string(id);
      if (!io.InquireAttribute<std::string>(attrName)) {
        io.DefineAttribute<std::string>(attrName, chNames[id]);
      }
    }

    // 2. Event data
    auto var_hash = io.InquireVariable<std::uint64_t>("hash");
    if (!var_hash) {
      var_hash = io.DefineVariable<std::uint64_t>("hash", {size, adios2::UnknownDim}, {rank, 0}, {1, adios2::UnknownDim});
      var_hash.SetMemorySpace(adios2::MemorySpace::Host);
    }

    auto var_start = io.InquireVariable<double>("start");
    if (!var_start) {
      var_start = io.DefineVariable<double>("start", {size, adios2::UnknownDim}, {rank, 0}, {1, adios2::UnknownDim});
      var_start.SetMemorySpace(adios2::MemorySpace::Host);
    }

    auto var_end = io.InquireVariable<double>("end");
    if (!var_end) {
      var_end = io.DefineVariable<double>("end", {size, adios2::UnknownDim}, {rank, 0}, {1, adios2::UnknownDim});
      var_end.SetMemorySpace(adios2::MemorySpace::Host);
    }

    // Per-event channel id (PROFILE_CHANNEL for OLB_PROFILE, else olb::logging id).
    auto var_channel = io.InquireVariable<std::uint8_t>("channel");
    if (!var_channel) {
      var_channel = io.DefineVariable<std::uint8_t>("channel", {size, adios2::UnknownDim}, {rank, 0}, {1, adios2::UnknownDim});
      var_channel.SetMemorySpace(adios2::MemorySpace::Host);
    }

    static struct {
      std::vector<std::uint64_t> hashes;
      std::vector<double> starts;
      std::vector<double> ends;
      std::vector<std::uint8_t> channels;
      std::string concatenated;
    } persistent_buffer;

    persistent_buffer.hashes.resize(nEvents);
    persistent_buffer.starts.resize(nEvents);
    persistent_buffer.ends.resize(nEvents);
    persistent_buffer.channels.resize(nEvents);

    for (std::size_t i = 0; i < nEvents; ++i) {
      persistent_buffer.hashes[i] = buffer[i].nameHash;
      persistent_buffer.starts[i] = buffer[i].start;
      persistent_buffer.ends[i] = buffer[i].end;
      persistent_buffer.channels[i] = buffer[i].channel;
    }

    engine.BeginStep();

    var_hash.SetSelection({{rank, 0}, {1, nEvents}});
    var_start.SetSelection({{rank, 0}, {1, nEvents}});
    var_end.SetSelection({{rank, 0}, {1, nEvents}});
    var_channel.SetSelection({{rank, 0}, {1, nEvents}});

    engine.Put(var_hash, persistent_buffer.hashes.data());
    engine.Put(var_start, persistent_buffer.starts.data());
    engine.Put(var_end, persistent_buffer.ends.data());
    engine.Put(var_channel, persistent_buffer.channels.data());

    const auto& messages = profiler.getMessages();
    if (!messages.empty()) {
      persistent_buffer.concatenated.clear();
      for (const auto& m : messages) {
        persistent_buffer.concatenated += std::to_string(m.timestamp) + ": " + m.message + "\n";
      }

      auto var_msg = io.InquireVariable<std::string>("messages");
      if (!var_msg) {
        var_msg = io.DefineVariable<std::string>("messages", {adios2::LocalValueDim});
      }
      engine.Put(var_msg, persistent_buffer.concatenated);
    }

    engine.EndStep();
#endif
}

} // namespace logging

} // namespace olb

#endif
