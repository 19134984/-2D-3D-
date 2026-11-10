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

#ifndef TIMER_HH
#define TIMER_HH

#include "timer.h"
#include "communication/mpiManager.h"
#include "communication/ompManager.h"
#include "utilities/performanceProfiler.h"
#include "io/textTable.h"

namespace olb {

namespace util {

template<typename T>
Timer<T>::Timer(int maxTimeSteps, size_t numFluidCells, unsigned int printModeSummary)
  : clout(std::cout,"Timer"),
    deltaTS(0), curTS(0), maxTS(maxTimeSteps),
    numFC(numFluidCells), rtRemMs(1),
    _printModeSummary(printModeSummary)
{}

template<typename T>
T Timer<T>::getMLUPs()
{
  double seconds = std::chrono::duration<double>(wallTimeCur - wallTimeLast).count();
  double updates = static_cast<double>(numFC) * static_cast<double>(deltaTS);
  return static_cast<T>(updates / (seconds * 1e6));
}

template<typename T>
T Timer<T>::getMLUPps()
{
  T mlupps = getMLUPs() / singleton::mpi().getSize();
#ifdef PARALLEL_MODE_OMP
  mlupps /= singleton::omp().getSize();
#endif
  return mlupps;
}

template<typename T>
T Timer<T>::getTotalMLUPs()
{
  double seconds = std::chrono::duration<double>(wallTimeEnd - wallTimeStart).count();
  double updates = static_cast<double>(numFC) * static_cast<double>(curTS);
  return static_cast<T>(updates / (seconds * 1e6));
}

template<typename T>
T Timer<T>::getTotalMLUPps()
{
  T tmlupps = getTotalMLUPs() / singleton::mpi().getSize();
#ifdef PARALLEL_MODE_OMP
  tmlupps /= singleton::omp().getSize();
#endif
  return tmlupps;
}

template<typename T>
void Timer<T>::start()
{
  wallTimeStart = Clock::now();
  wallTimeCur = wallTimeStart;
  wallTimeLast = wallTimeStart;
  cpuTimeStart = std::clock();
}

template<typename T>
void Timer<T>::update(int currentTimeStep)
{
  OLB_PROFILE_EVENT("Host:Timer:Update");
  cpuTimeCur = std::clock();
  wallTimeLast = wallTimeCur;
  wallTimeCur = Clock::now();

  deltaTS = currentTimeStep - curTS;
  curTS   = currentTimeStep;

  double wallSecs = std::chrono::duration<double>(wallTimeCur - wallTimeStart).count();
  double wallMs = std::chrono::duration<double, std::milli>(wallTimeCur - wallTimeStart).count();
  double cpuSecs = static_cast<double>(cpuTimeCur - cpuTimeStart) / CLOCKS_PER_SEC;

  rtPas   = static_cast<T>(wallSecs);
  rtTot   = static_cast<T>(wallSecs * maxTS / util::max<int>(curTS, 1));
  rtRem   = rtTot - rtPas;

  rtPasMs = static_cast<T>(wallMs);
  rtTotMs = static_cast<T>(wallMs * maxTS / util::max<int>(curTS, 1));
  rtRemMs = rtTotMs - rtPasMs;

  ctPas   = cpuSecs;
  ctTot   = ctPas * maxTS / util::max<int>(curTS, 1);
  ctRem   = ctTot - ctPas;
}

template<typename T>
void Timer<T>::stop()
{
  cpuTimeEnd = std::clock();
  wallTimeEnd = Clock::now();
}

template<typename T>
void Timer<T>::reset()
{
  deltaTS = 0;
  curTS = 0;
  rtRemMs = 1;

  ltPas = 0;
  ltRem = 0;
  ltTot = 0;

  ctPas = 0.0;
  ctRem = 0.0;
  ctTot = 0.0;

  rtPas = T{};
  rtRem = T{};
  rtTot = T{};

  rtPasMs = T{};
  rtTotMs = T{};
}

template<typename T>
double Timer<T>::getTotalCpuTime()
{
  return static_cast<double>(cpuTimeEnd - cpuTimeStart) / CLOCKS_PER_SEC;
}

template<typename T>
T Timer<T>::getTotalRealTime()
{
  return static_cast<T>(std::chrono::duration<double>(wallTimeEnd - wallTimeStart).count());
}

template<typename T>
T Timer<T>::getTotalRealTimeMs()
{
  return static_cast<T>(std::chrono::duration<double, std::milli>(wallTimeEnd - wallTimeStart).count());
}

template<typename T>
void Timer<T>::print(std::size_t currentTimeStep, int printMode)
{
  if (currentTimeStep != curTS) {
    update(currentTimeStep);
  }
  printStep(printMode);
}

template<typename T>
void Timer<T>::printStep(int printMode)
{
  printStep(clout, printMode);
}

template<typename T>
void Timer<T>::printStep(OstreamManager& fout, int printMode)
{
  switch (printMode) {
  case 0:
    io::text::StatusLine()
      .field("step", curTS)
      .field("percent", 100.0 * curTS / maxTS)
      .field("passedTime", static_cast<double>(rtPasMs) / 1000.0)
      .field("remTime", static_cast<double>(rtRemMs) / 1000.0)
      .field("MLUPs", getMLUPs())
      .render(fout);
    break;

  case 1:
    fout
        << "latticeTS: "
        << curTS << "/" << maxTS << " (" << 100 * curTS / maxTS << "%); "
        << "pas/totTime: "
        << std::setprecision(2) << std::fixed << static_cast<double>(rtPasMs) / 1000.0 << "/"
        << std::setprecision(1) << std::fixed << static_cast<double>(rtTotMs) / 1000.0 << "s; "
        << "remTime: "
        << std::setw(2) << static_cast<int>(static_cast<double>(rtRemMs) / 1000.0) / 60 << "m "
        << std::setfill('0') << std::setw(4) << static_cast<double>(static_cast<int>(static_cast<double>(rtRemMs) / 100.0) % 600) / 10.0 << "s; "
        << std::setfill(' ')
        << "MLUPs: " << getMLUPs()
        << std::endl;
    break;

  case 2:
    fout
        << std::setw(21) << std::left << "Lattice-Timesteps"
        << std::setw(17) << std::left << "| CPU time/estim"
        << std::setw(18) << std::left << "| REAL time/estim"
        << std::setw(6)  << std::left << "| ETA"
        << std::setw(6)  << std::left << "| MLUPs"
        << std::endl << std::right
        << std::setw(6) << std::setprecision(2) << std::fixed << curTS << "/" << std::setw(6) << maxTS << " (" << std::setw(3) << 100 * curTS / maxTS << "%) |"
        << std::setw(7) << ctPas << "/" << std::setw(7) << ctTot << " |"
        << std::setw(8) << static_cast<double>(rtPasMs) / 1000.0 << "/" << std::setw(7) << static_cast<double>(rtTotMs) / 1000.0 << " |"
        << std::setw(4) << static_cast<int>(rtRemMs) / 1000 + 1 << " |"
        << std::setw(6) << getMLUPs()
        << std::endl;
    break;

  case 3:
    fout
        << "step " << curTS << "; "
        << "MLUPs=" << std::setw(8) << getMLUPs() << ", MLUPps=" << std::setw(8) << getMLUPps() << std::endl;
    break;

  default:
    fout << "Error in function printStep in class_timer.h: printMode=" << printMode << " not found" << std::endl << std::flush;
  }
}

template<typename T>
void Timer<T>::printSummary()
{
  printSummary(clout);
}

template<typename T>
void Timer<T>::printSummary(OstreamManager& fout)
{
  fout << std::endl;
  fout << "----------------Summary:Timer----------------" << std::endl;
  fout << "measured time (rt) : " << static_cast<int>(getTotalRealTimeMs()) / 1000 << "." << static_cast<int>(getTotalRealTimeMs()) - static_cast<int>(getTotalRealTimeMs()) / 1000 * 1000 << "s" << std::endl;
  fout << "measured time (cpu): " << std::setprecision(3) << std::fixed << getTotalCpuTime() << "s" << std::endl;
  if (numFC > 0 && curTS > 0) {
    fout << "average MLUPs :       " << getTotalMLUPs()  << std::endl;
    fout << "average MLUPps:       " << getTotalMLUPps() << std::endl;
  }
  fout << "---------------------------------------------" << std::endl;
}

template<typename T>
void Timer<T>::printShortSummary()
{
  printShortSummary(clout);
}

template<typename T>
void Timer<T>::printShortSummary(OstreamManager& fout)
{
  fout << "realTime=" << static_cast<int>(getTotalRealTimeMs()) / 1000 << "." << static_cast<int>(getTotalRealTimeMs()) - static_cast<int>(getTotalRealTimeMs()) / 1000 * 1000
       << "; cpuTime=" <<  std::setprecision(3) << std::fixed << getTotalCpuTime() << std::endl;
}

template<typename T, typename DESCRIPTOR>
Timer<T>* createTimer(XMLreader& param, const UnitConverter<T,DESCRIPTOR>& converter, size_t numLatticePoints)
{
  OstreamManager clout(std::cout,"createTimer");

  T physMaxT = T();
  T physStartT = T();
  int printModeSummary = 0;

  if ( ! param["Application"]["PhysParameters"]["PhysMaxTime"].read(physMaxT) ) {
    if ( ! param["Application"]["PhysParam"]["MaxTime"].read(physStartT) ) {
      clout << "PhysMaxTime not found" << std::endl;
    }
    else {
      clout << "Application::PhysParam::MaxTime needs to be renamed to Application::PhysParameters::PhysMaxTime" << std::endl;
    }
  }

  param.readOrWarn<int>("Output","Timer","PrintModeSummary", printModeSummary,true,false,false);

  int maxT = converter.getLatticeTime(physMaxT) + converter.getLatticeTime(physStartT);

  return new Timer<T>(maxT, numLatticePoints, printModeSummary);
}

} // namespace util

} // namespace olb

#endif
