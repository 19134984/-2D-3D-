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

#ifndef CASE_TIMER_REPORTER_H
#define CASE_TIMER_REPORTER_H

#include "reporter.h"
#include "refinement/lattice.h"
#include "io/textTable.h"
#include "io/ostreamManager.h"
#include "utilities/timer.h"

#include <memory>
#include <fstream>
#include <iomanip>

namespace olb {

namespace reporter {

/// Progress and throughput, gated by PHYS_ITER_T; shows ETA when MAX_PHYS_T is set
class Timer {
private:
  std::unique_ptr<util::Timer<double>> _timer;
  bool          _writeCsv;
  std::ofstream _csv;

  template <typename CASE>
  void writeSample(CASE& myCase) {
    if (!_writeCsv || !singleton::mpi().isMainProcessor()) {
      return;
    }
    if (!_csv.is_open()) {
      _csv.open(singleton::directories().getProfilingOutDir() + "performance.csv");
      _csv << "size, steps, ranks, MLUPs\n";
    }
    _csv << myCase.getParameters().template get<olb::parameters::RESOLUTION>() << ", "
         << maxSteps(myCase) << ", "
         << singleton::mpi().getSize() << ", "
         << std::fixed << std::setprecision(2) << _timer->getMLUPs() << "\n";
    _csv.flush();
  }

  template <typename CASE>
  static decltype(auto) referenceLattice(CASE& myCase) {
    auto& component = myCase.get(typename CASE::map_t::keys_t::template get<0>{});
    if constexpr (CASE::is_refined) {
      return component.getLevel(0);
    } else {
      return component;
    }
  }

  /// Cell updates per coarse step (MLUPs denominator); refined levels weighted by 2^level substeps
  template <typename CASE>
  static std::size_t collideCells(CASE& myCase) {
    auto& geometry = myCase.getGeometry();
    if constexpr (CASE::is_refined) {
      std::size_t total = 0;
      for (std::size_t iLevel = 0; iLevel < geometry.getLevelCount(); ++iLevel) {
        total += geometry.getLevel(iLevel).getStatistics().getNvoxel() * (std::size_t{1} << iLevel);
      }
      return total;
    } else {
      return geometry.getStatistics().getNvoxel();
    }
  }

  /// Total steps for progress: MAX_PHYS_T if set, else the step-count TIME_STEPS
  template <typename CASE>
  static std::size_t maxSteps(CASE& myCase) {
    const double maxPhysT = myCase.getParameters().template get<olb::parameters::PHYS_MAX_T>();
    if (maxPhysT > 0) {
      return referenceLattice(myCase).getUnitConverter().getLatticeTime(maxPhysT);
    }
    return myCase.getParameters().template get<olb::parameters::TIME_STEPS>();
  }

  /// Report cadence in steps: physical PHYS_ITER_T if set, else step-count LATTICE_ITER_T
  template <typename CASE, typename PARAMETERS>
  static std::size_t interval(CASE& myCase, PARAMETERS& params) {
    const double physInterval = params.template get<olb::parameters::PHYS_ITER_T>();
    if (physInterval > 0) {
      return referenceLattice(myCase).getUnitConverter().getLatticeTime(physInterval);
    }
    return static_cast<std::size_t>(params.template get<olb::parameters::LATTICE_ITER_T>());
  }

public:
  using parameters = meta::list<olb::parameters::PHYS_ITER_T,
                                olb::parameters::LATTICE_ITER_T>;

  /// writeCsv emits a machine-readable performance row to profilingData at finalize
  explicit Timer(bool writeCsv = false) : _writeCsv(writeCsv) { }

  template <typename CASE, typename PARAMETERS>
  void report(CASE& myCase, std::size_t iT, PARAMETERS& params) {
    if (!_timer) {
      _timer = std::make_unique<util::Timer<double>>(maxSteps(myCase), collideCells(myCase));
      _timer->start();
    }
    const std::size_t every = interval(myCase, params);
    if (every == 0 || iT % every != 0) {
      return;
    }
    _timer->update(static_cast<int>(iT));
    const double physT = referenceLattice(myCase).getUnitConverter().getPhysTime(iT);
    io::text::StatusLine line;
    line.field("step", iT).field("t", io::text::fixed(physT) + "s");
    if (maxSteps(myCase) > 0) {
      line.field("progress", io::text::percent(_timer->getPercentage() / 100.0))
          .field("eta", io::text::fixed(_timer->getRealTimeRemaining(), 1) + "s");
    }
    line.token(io::text::lups(_timer->getMLUPs() * 1e6));
    OstreamManager clout(std::cout, "Timer");
    line.render(clout);

    if (iT > 0) {
      writeSample(myCase);
    }
  }

  template <typename CASE, typename PARAMETERS>
  void finalize(CASE& myCase, std::size_t iT, PARAMETERS& params) {
    if (!_timer) {
      return;
    }
    _timer->update(static_cast<int>(iT));
    _timer->stop();
    OstreamManager clout(std::cout, "Timer");
    io::text::StatusLine()
      .token("done")
      .field("steps", iT)
      .field("wall", io::text::fixed(_timer->getTotalRealTime(), 2) + "s")
      .token(io::text::lups(_timer->getTotalMLUPs() * 1e6))
      .render(clout);

    writeSample(myCase);
  }
};

}

}

#endif
