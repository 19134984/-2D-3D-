/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007, 2008 Jonas Latt
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

#ifndef LATTICE_STATISTICS_HH
#define LATTICE_STATISTICS_HH

#include "latticeStatistics.h"
#include "io/textTable.h"
#include "utilities/omath.h"
#include <limits>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <filesystem>
#include "util.h"

namespace olb {

////////////////////// Class LatticeStatistics /////////////////

template<typename T, typename DESCRIPTOR>
std::size_t LatticeStatistics<T,DESCRIPTOR>::_instanceCounter = 0;

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T,DESCRIPTOR>::LatticeStatistics()
  : clout(std::cout,"LatticeStatistics")
  , _id(_instanceCounter++)
{
  initialize();
}

template<typename T, typename DESCRIPTOR>
void LatticeStatistics<T,DESCRIPTOR>::initialize()
{
  _partial   = Aggregatable{};
  _avRho     = T(1);
  _avEnergy  = T();
  _maxU      = T();
  _minRho    = T(1);
  _maxRho    = T(1);
  _minU      = T();
  _meanU     = T();
  _numCells  = 0;
  _firstCall = true;
  resetTime();
}

template<typename T, typename DESCRIPTOR>
void LatticeStatistics<T,DESCRIPTOR>::reset()
{
  const auto& a = _partial;
  if (a.nCells == 0) {
    _avRho = _avEnergy = _maxU = _minRho = _maxRho = _minU = _meanU = T();
    _numCells  = 0;
    _firstCall = false;
  }
  else {
    // The average density is used by the "ConstRhoBGK" model. Before the
    // simulation starts it may hold a nonsensical value, so averages are
    // initialized to 1 on the first call.
    if (_firstCall) {
      _avRho    = T(1);
      _avEnergy = T(1);
      _firstCall = false;
    }
    else {
      _avRho    = a.avRho    / T(a.nCells);
      _avEnergy = a.avEnergy / T(a.nCells);
    }
    _avEnergy *= T(0.5);                  // energy is 0.5*uSqr
    _maxU   = util::sqrt(a.maxU);         // u is sqrt(uSqr)
    _minRho = a.minRho;
    _maxRho = a.maxRho;
    _minU   = util::sqrt(a.minUSqr);
    _meanU  = a.sumU / T(a.nCells);
    _numCells = a.nCells;
  }
  _partial = Aggregatable{};
}

template<typename T, typename DESCRIPTOR>
typename LatticeStatistics<T,DESCRIPTOR>::Result
LatticeStatistics<T,DESCRIPTOR>::getResult() const
{
  return Result{_avRho, _avEnergy, _maxU, _minRho, _maxRho, _minU, _meanU, _numCells};
}

template<typename T, typename DESCRIPTOR>
void LatticeStatistics<T,DESCRIPTOR>::setResult(const Result& r)
{
  _avRho    = r.avRho;
  _avEnergy = r.avEnergy;
  _maxU     = r.maxU;
  _minRho   = r.minRho;
  _maxRho   = r.maxRho;
  _minU     = r.minU;
  _meanU    = r.meanU;
  _numCells = r.nCells;
  _partial  = Aggregatable{};
}

template<typename T, typename DESCRIPTOR>
void LatticeStatistics<T,DESCRIPTOR>::incrementStats(T rho, T uSqr)
{
  _partial.increment(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
void LatticeStatistics<T,DESCRIPTOR>::incrementStats(Aggregatable& aggregatable)
{
  _partial += aggregatable;
}

template<typename T, typename DESCRIPTOR>
void LatticeStatistics<T,DESCRIPTOR>::incrementStats()
{
  _partial.nCells += 1;
}

template<typename T, typename DESCRIPTOR>
const typename LatticeStatistics<T,DESCRIPTOR>::RegisteredField*
LatticeStatistics<T,DESCRIPTOR>::findRegistered(const std::string& name) const
{
  for (const auto& rf : _registered) {
    if (rf.name == name) {
      return &rf;
    }
  }
  return nullptr;
}

template<typename T, typename DESCRIPTOR>
StatisticEntry<T>
LatticeStatistics<T,DESCRIPTOR>::entryFor(const RegisteredField& rf) const
{
  StatisticEntry<T> e;
  e.name       = rf.name;
  e.components = rf.components;
  e.min.assign(rf.components, T());
  e.max.assign(rf.components, T());
  e.mean.assign(rf.components, T());
  const auto it = _fieldResults.find(rf.id);
  if (it != _fieldResults.end() && it->second.count > 0) {
    const auto& r    = it->second;
    const T     norm = T(1) / T(r.count);
    for (unsigned c = 0; c < rf.components; ++c) {
      if (c < r.min.size()) e.min[c]  = r.min[c];
      if (c < r.max.size()) e.max[c]  = r.max[c];
      if (c < r.sum.size()) e.mean[c] = r.sum[c] * norm;
    }
  }
  return e;
}

template<typename T, typename DESCRIPTOR>
T LatticeStatistics<T,DESCRIPTOR>::getMin(const std::string& name, unsigned iComponent) const
{
  if (name == "rho") return _minRho;
  if (name == "u")   return _minU;
  if (const auto* rf = findRegistered(name)) {
    const auto e = entryFor(*rf);
    if (iComponent < e.min.size()) return e.min[iComponent];
  }
  return T();
}

template<typename T, typename DESCRIPTOR>
T LatticeStatistics<T,DESCRIPTOR>::getMax(const std::string& name, unsigned iComponent) const
{
  if (name == "rho") return _maxRho;
  if (name == "u")   return _maxU;
  if (const auto* rf = findRegistered(name)) {
    const auto e = entryFor(*rf);
    if (iComponent < e.max.size()) return e.max[iComponent];
  }
  return T();
}

template<typename T, typename DESCRIPTOR>
T LatticeStatistics<T,DESCRIPTOR>::getMean(const std::string& name, unsigned iComponent) const
{
  if (name == "rho") return _avRho;
  if (name == "u")   return _meanU;
  if (const auto* rf = findRegistered(name)) {
    const auto e = entryFor(*rf);
    if (iComponent < e.mean.size()) return e.mean[iComponent];
  }
  return T();
}

template<typename T, typename DESCRIPTOR>
std::vector<StatisticEntry<T>> LatticeStatistics<T,DESCRIPTOR>::getEntries() const
{
  std::vector<StatisticEntry<T>> entries;
  entries.push_back(StatisticEntry<T>{"rho", 1, {_minRho}, {_maxRho}, {_avRho}});
  entries.push_back(StatisticEntry<T>{"u",   1, {_minU},   {_maxU},   {_meanU}});
  entries.push_back(StatisticEntry<T>{"energy", 1,
    {T{0.5}*_minU*_minU}, {T{0.5}*_maxU*_maxU}, {_avEnergy}});
  for (const auto& rf : _registered) {
    entries.push_back(entryFor(rf));
  }
  return entries;
}

template<typename T, typename DESCRIPTOR>
void LatticeStatistics<T,DESCRIPTOR>::incrementTime()
{
  ++_latticeTime;
}

template<typename T, typename DESCRIPTOR>
void LatticeStatistics<T,DESCRIPTOR>::resetTime(size_t value)
{
  _latticeTime = value;
}

template<typename T, typename DESCRIPTOR>
size_t LatticeStatistics<T,DESCRIPTOR>::getTime() const
{
  return _latticeTime;
}

template<typename T, typename DESCRIPTOR>
void LatticeStatistics<T,DESCRIPTOR>::printWithoutPlot(int iterationStep, T physicalTime) const
{
  std::stringstream extra;
  for (const auto& rf : _registered) {
    const auto e = entryFor(rf);
    extra << "; " << e.name << "=[";
    for (unsigned iC = 0; iC < e.components; ++iC) {
      if (iC > 0) extra << ", ";
      extra << e.min[iC] << "/" << e.max[iC] << "/" << e.mean[iC];
    }
    extra << "]";
  }
  clout << "step="      << iterationStep      << "; "
        << "physT="     << physicalTime       << "; "
        << "maxCFL="    << getMaxU()          << "; "
        << "avgEnergy=" << getAverageEnergy() << "; "
        << "avgRho="    << getAverageRho()
        << extra.str()
        << std::endl;
}

template<typename T, typename DESCRIPTOR>
void LatticeStatistics<T,DESCRIPTOR>::print(int iterationStep, T physicalTime, const std::string& postfix)
{
  const std::string fallback = _name.empty() ? std::to_string(_id) : _name;
  std::string plotPostfix = postfix.empty() ? fallback : postfix;
  if (!gplotAvgEnergy) {
    gplotAvgEnergy = std::make_unique<Gnuplot<T>>("avgEnergy_" + plotPostfix, false);
    gplotMaxCFL = std::make_unique<Gnuplot<T>>("maxCFL_" + plotPostfix, false);
  }

  printWithoutPlot(iterationStep, physicalTime);

  // Gnuplot output (avgEnergy / maxCFL PNGs)
  gplotAvgEnergy->setData(physicalTime, getAverageEnergy(), "avgEnergy", "bottom right", 'l');
  gplotMaxCFL->setData(physicalTime, getMaxU(), "maxCFL", "bottom right", 'l');
  gplotAvgEnergy->writePNG(-1, -1, "avgEnergy_" + plotPostfix);
  gplotMaxCFL->writePNG(-1, -1, "maxCFL_" + plotPostfix);
}

template<typename T, typename DESCRIPTOR>
void LatticeStatistics<T,DESCRIPTOR>::printSummary(int iterationStep, T physicalTime) const
{
  const auto entries = getEntries();
  std::size_t width = 5;
  for (const auto& e : entries) {
    for (unsigned iC = 0; iC < e.components; ++iC) {
      std::size_t w = e.name.size() + (e.components > 1 ? 3 : 0);
      width = std::max(width, w);
    }
  }
  auto rule = [&](){
    std::stringstream s;
    s << std::string(width+2, '-') << "+" << std::string(11,'-') << "+"
      << std::string(11,'-') << "+" << std::string(11,'-');
    clout << s.str() << std::endl;
  };
  clout << "step " << iterationStep << ", physT " << physicalTime << " s" << std::endl;
  clout << " " << std::setw(width) << std::left << "field" << std::right
        << " | " << std::setw(9) << "min"
        << " | " << std::setw(9) << "max"
        << " | " << std::setw(9) << "mean" << std::endl;
  rule();
  for (const auto& e : entries) {
    for (unsigned iC = 0; iC < e.components; ++iC) {
      std::string label = e.name;
      if (e.components > 1) label += "[" + std::to_string(iC) + "]";
      clout << " " << std::setw(width) << std::left << label << std::right
            << std::scientific << std::setprecision(3)
            << " | " << std::setw(9) << e.min[iC]
            << " | " << std::setw(9) << e.max[iC]
            << " | " << std::setw(9) << e.mean[iC] << std::endl;
    }
  }
  rule();
}

template<typename T, typename DESCRIPTOR>
void LatticeStatistics<T,DESCRIPTOR>::writeCSV(int iterationStep, T physicalTime, const std::string& name) const
{
  if (!singleton::mpi().isMainProcessor()) {
    return;
  }
  const auto entries = getEntries();
  const std::string dir = singleton::directories().getStatisticsOutDir();
  std::error_code ec;
  std::filesystem::create_directories(dir, ec);
  const std::string fallback = _name.empty() ? ("lattice" + std::to_string(_id)) : _name;
  const std::string base = io::text::filenameSafe(name.empty() ? fallback : name);
  const std::string file = dir + base + "_statistics.csv";

  const bool exists = std::filesystem::exists(file);
  std::ofstream out(file, std::ios::app);
  if (!out) {
    return;
  }
  if (!exists) {
    out << "iT,physT";
    for (const auto& e : entries) {
      for (unsigned iC = 0; iC < e.components; ++iC) {
        std::string label = e.name;
        if (e.components > 1) label += std::to_string(iC);
        out << "," << label << "_min," << label << "_max," << label << "_mean";
      }
    }
    out << "\n";
  }
  out << std::scientific << std::setprecision(8) << iterationStep << "," << physicalTime;
  for (const auto& e : entries) {
    for (unsigned iC = 0; iC < e.components; ++iC) {
      out << "," << e.min[iC] << "," << e.max[iC] << "," << e.mean[iC];
    }
  }
  out << "\n";
}

}  // namespace olb

#endif
