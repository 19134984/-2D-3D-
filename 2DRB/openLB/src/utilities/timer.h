/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2011 Lukas Baron, Mathias J. Krause
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

#ifndef TIMER_H
#define TIMER_H

#include "io/ostreamManager.h"
#include "io/xmlReader.h"
#include "core/unitConverter.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <sys/time.h>

namespace olb {

namespace util {

template<typename T>
class Timer {
private:
  OstreamManager clout;

  std::clock_t cpuTimeStart{}, cpuTimeCur{}, cpuTimeEnd{};

  using Clock = std::chrono::high_resolution_clock;
  using TimePoint = Clock::time_point;

  TimePoint wallTimeStart{}, wallTimeCur{}, wallTimeEnd{}, wallTimeLast{};

  std::size_t deltaTS{};
  std::size_t curTS{};
  std::size_t maxTS{};
  std::size_t numFC{};

  int    ltPas{},   ltRem{},   ltTot{};
  double ctPas{},   ctRem{},   ctTot{};
  T      rtPas{},   rtRem{},   rtTot{};
  T      rtPasMs{}, rtRemMs{}, rtTotMs{};

  unsigned int _printModeSummary;

public:
  Timer(int maxTimeSteps, std::size_t numFluidCells=1, unsigned int printModeSummary = 0);

  T getMLUPs();
  T getMLUPps();
  T getTotalMLUPs();
  T getTotalMLUPps();

  /// Progress in percent, 0 when no step target was set
  int getPercentage() const {
    return maxTS > 0 ? static_cast<int>(100 * curTS / maxTS) : 0;
  }
  /// Estimated remaining real time in seconds
  T getRealTimeRemaining() const {
    return rtRemMs / 1000;
  }

  void start();
  void update(int currentTimeStep);
  void stop();
  void reset();

  double getTotalCpuTime();
  T getTotalRealTime();
  T getTotalRealTimeMs();

  void printStep(int printMode=0);
  void printStep(OstreamManager& fout, int printMode=0);
  void print(std::size_t currentTimeStep, int printMode=0);

  void printSummary();
  void printSummary(OstreamManager& fout);
  void printShortSummary();
  void printShortSummary(OstreamManager& fout);
};

// Factory function /////////////////////////////////
template<typename T, typename DESCRIPTOR>
Timer<T>* createTimer(XMLreader& param, const UnitConverter<T,DESCRIPTOR>& converter, size_t numLatticePoints);
/////////////////////////////////////////////////////

} // namespace util

} // namespace olb

#endif
