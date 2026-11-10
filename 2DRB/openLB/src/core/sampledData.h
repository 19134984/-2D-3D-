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

#ifndef CORE_SAMPLED_DATA_H
#define CORE_SAMPLED_DATA_H

#include <span>
#include <string>
#include <vector>
#include <fstream>

#include "singleton.h"

namespace olb {

/// Zero-copy view of sampled point data, decoupling extraction from writers
/**
 * A lightweight, non-owning view over an extractor's compact local buffers. Because
 * the underlying field storage is structure-of-arrays, both the point coordinates and
 * each result field are exposed as one span per component (deviation from the v3
 * interleaved-single-span sketch, dictated by the actual column storage). The
 * extractor owns the memory; a view is valid only until the next update() on that
 * extractor. Re-acquire after every update(); never store one across a step.
 **/
template <typename T>
struct SampledData {
  /// One named result field: nComp component spans, each of length nLocal
  struct FieldView {
    std::string name;
    unsigned    nComp;
    std::vector<std::span<const T>> components;
  };

  // flat geometry: one coordinate span per spatial axis (each of length nLocal)
  std::vector<std::span<const T>> pointCoords;

  // named result fields (zero-copy component spans)
  std::vector<FieldView> fields;

  // distribution metadata
  std::span<const std::size_t> globalIndices; // local -> global
  std::size_t nLocal = 0;
  std::size_t nGlobal = 0;
  int rank = 0;
  int nRanks = 1;
};

/// Writer abstraction consuming a SampledData view
template <typename T>
struct SampleWriter {
  virtual ~SampleWriter() = default;
  virtual void write(const SampledData<T>& data, const std::string& path, std::size_t iT) = 0;
};

/// Reference writer: each rank streams its owned subset as CSV (distributed, no gather)
/**
 * Serves the flat form directly from the distribution metadata in the view -- the
 * exascale composition (distributed view + per-rank parallel writer) in miniature,
 * and a backend-agnostic milestone for the SampledData boundary.
 **/
template <typename T>
class CsvSampleWriter final : public SampleWriter<T> {
public:
  void write(const SampledData<T>& data, const std::string& path, std::size_t iT) override
  {
    std::string fileName = path + "_" + std::to_string(iT)
                         + "_rank" + std::to_string(data.rank) + ".csv";
    std::ofstream out(fileName);
    out << "iGlobal";
    for (std::size_t d=0; d < data.pointCoords.size(); ++d) {
      out << ",x" << d;
    }
    for (const auto& f : data.fields) {
      for (unsigned c=0; c < f.nComp; ++c) {
        out << ',' << f.name << c;
      }
    }
    out << '\n';
    for (std::size_t i=0; i < data.nLocal; ++i) {
      out << data.globalIndices[i];
      for (std::size_t d=0; d < data.pointCoords.size(); ++d) {
        out << ',' << data.pointCoords[d][i];
      }
      for (const auto& f : data.fields) {
        for (unsigned c=0; c < f.nComp; ++c) {
          out << ',' << f.components[c][i];
        }
      }
      out << '\n';
    }
  }
};

}

#endif
