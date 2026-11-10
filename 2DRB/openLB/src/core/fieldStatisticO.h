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

#ifndef CORE_FIELD_STATISTIC_O_H
#define CORE_FIELD_STATISTIC_O_H

#include <limits>

#include "latticeStatistics.h"

namespace olb {

/// Per-block min/max/mean reduction of a registered descriptor FIELD
/**
 * Runs at statistics time (stage::Statistics) and writes its per-block result
 * into the block's LatticeStatistics. The generic type reduces the core cells
 * on the host (the field is synced to the host first, a no-op on CPU); GPU
 * platforms are specialized for a device reduction in the platform statistics
 * headers. Installed by LatticeStatistics::registerField.
 **/
template<typename FIELD>
struct FieldStatisticO {
  static constexpr OperatorScope scope = OperatorScope::PerBlock;

  int getPriority() const {
    return 0;
  }

  template<typename BLOCK>
  struct type {
    void setup(BLOCK&) { }

    void apply(BLOCK& block) {
      if (!block.statisticsEnabled()) {
        return;
      }
      using T          = typename BLOCK::value_t;
      using DESCRIPTOR = typename BLOCK::descriptor_t;
      const unsigned d = DESCRIPTOR::template size<FIELD>();

      // Host-side reduction of the core cells. Ensure the field is host-resident
      // (no-op on CPU). GPU platforms are specialized below for a device reduction.
      block.template getData<Array<FIELD>>().setProcessingContext(ProcessingContext::Evaluation);

      typename LatticeStatistics<T,DESCRIPTOR>::FieldResult r;
      r.components = d;
      r.min.assign(d, std::numeric_limits<T>::max());
      r.max.assign(d, std::numeric_limits<T>::lowest());
      r.sum.assign(d, T{});
      r.count = 0;

      block.forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
        auto cell = block.get(loc);
        for (unsigned iD = 0; iD < d; ++iD) {
          const T v = cell.template getFieldComponent<FIELD>(iD);
          r.min[iD] = util::min(r.min[iD], v);
          r.max[iD] = util::max(r.max[iD], v);
          r.sum[iD] += v;
        }
        r.count += 1;
      });

      block.getStatistics().template setFieldResult<FIELD>(r);
    }
  };

  template<typename BLOCK>
  void setup(BLOCK& block) {
    type<BLOCK>{}.setup(block);
  }

  template<typename BLOCK>
  void apply(BLOCK& block) {
    type<BLOCK>{}.apply(block);
  }
};

}  // namespace olb

#endif
