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

#ifndef IO_BINARY_TIME_SERIES_WRITER_H
#define IO_BINARY_TIME_SERIES_WRITER_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iterator>

namespace olb {

template <typename T>
class BinaryTimeSeriesWriter {
private:
  std::string _baseName;
  std::ofstream _binFile;
  std::vector<T> _buffer;

  std::size_t _nSamples;
  std::size_t _nComponents;
  std::size_t _bufferSize;
  std::size_t _timestepsWritten;
  const T _deltaT;

  std::vector<Vector<T,3>> _sampleRs;

  void writeDescription() {
    std::ofstream meta(_baseName + ".meta");
    meta << "filename=" << _baseName << ".bin" << "\n";
    meta << "dtype=" << (sizeof(T) == 8 ? "float64" : "float32") << "\n";
    meta << "n_points=" << _nSamples << "\n";
    meta << "n_components=" << _nComponents << "\n";
    meta << "n_timesteps_written=" << _timestepsWritten << "\n";
    meta << "timestep=" << _deltaT << "\n";

    meta << "---COORDINATES[x,y,z]---" << "\n";
    for(const auto& p : _sampleRs) {
      meta << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
    meta.close();
  }

public:
  template <typename ITER>
  BinaryTimeSeriesWriter(std::string filename,
                         ITER first, ITER last,
                         std::size_t nComponents,
                         T deltaT,
                         std::size_t bufferSteps = 100)
    : _baseName(filename),
      _nSamples(std::distance(first, last)),
      _nComponents(nComponents),
      _bufferSize(bufferSteps * std::distance(first, last) * nComponents),
      _timestepsWritten(0),
      _deltaT(deltaT),
      _sampleRs(first, last)
  {
    _binFile.open(_baseName + ".bin", std::ios::binary | std::ios::app);
    _buffer.reserve(_bufferSize);
    writeDescription();
  }

  BinaryTimeSeriesWriter(std::string filename,
                         const std::vector<Vector<T,3>>& sampleRs,
                         std::size_t nComponents,
                         T deltaT,
                         std::size_t bufferSteps = 100)
    : BinaryTimeSeriesWriter(filename, sampleRs.begin(), sampleRs.end(), nComponents, deltaT, bufferSteps)
  {}

  ~BinaryTimeSeriesWriter() {
    flush();
    _binFile.close();
  }

  void bufferSample(const std::vector<T>& data) {
    if (data.size() != _nSamples * _nComponents) {
      std::cerr << "Size mismatch, expected "
                << _nSamples * _nComponents << " got " << data.size() << std::endl;
      return;
    }

    _buffer.insert(_buffer.end(), data.begin(), data.end());

    if (_buffer.size() >= _bufferSize) {
      flush();
    }
  }

  void flush() {
    writeDescription();

    if (_buffer.empty()) {
      return;
    }

    _binFile.write(reinterpret_cast<const char*>(_buffer.data()), _buffer.size() * sizeof(T));

    std::size_t stepsInBatch = _buffer.size() / (_nSamples * _nComponents);
    _timestepsWritten += stepsInBatch;

    _buffer.clear();
  }

};

}

#endif
