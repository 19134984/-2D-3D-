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

#ifndef REFINEMENT_VTK_REPORTER_H
#define REFINEMENT_VTK_REPORTER_H

#include "case/reporter.h"
#include "lattice.h"
#include "amrWriter.h"

#include <string>
#include <utility>

namespace olb {

namespace refinement {

/// Per-level AMR VTK output for the NAME refined lattice, gated by PHYS_ITER_T
template <typename NAME, typename HOOK>
class VTKreporter {
private:
  std::string _fileName;
  HOOK        _addFunctors;
  bool        _masterWritten = false;

public:
  using parameters = meta::list<olb::parameters::PHYS_ITER_T>;

  VTKreporter(std::string fileName, HOOK hook)
    : _fileName(std::move(fileName)), _addFunctors(std::move(hook))
  { }

  template <typename CASE, typename PARAMETERS>
  void report(CASE& myCase, std::size_t iT, PARAMETERS& params) {
    auto& refinedLattice = myCase.getLattice(NAME{});
    auto& conv0 = refinedLattice.getLevel(0).getUnitConverter();

    const double physInterval = params.template get<olb::parameters::PHYS_ITER_T>();
    if (physInterval <= 0) {
      return;
    }
    const std::size_t every = conv0.getLatticeTime(physInterval);
    if (every == 0 || iT % every != 0) {
      return;
    }

    RefinedAMRwriter writer(myCase.getMesh(), _fileName);
    if (!_masterWritten) {
      writer.createMasterFile();
      _masterWritten = true;
    }
    refinedLattice.synchronizeForPostProcessing();
    writer.write(iT, refinedLattice, [this, iT](auto& levelWriter, std::size_t iLevel, auto& sLattice) {
      _addFunctors(levelWriter, iLevel, sLattice, iT);
    });
  }
};

/// Register a VTK reporter under group, writing to fileName; hook adds per-level functors
template <typename CASE, typename NAME, typename HOOK>
auto& setVTKreporter(CASE& myCase, std::string group, NAME, std::string fileName, HOOK&& hook) {
  return myCase.template setReporter<VTKreporter<NAME, std::decay_t<HOOK>>>(
    std::move(group), std::move(fileName), std::forward<HOOK>(hook));
}

}

}

#endif
