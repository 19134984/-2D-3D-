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

#ifndef REFINEMENT_SAMPLING_H
#define REFINEMENT_SAMPLING_H

#include <vector>
#include <memory>
#include <type_traits>

#include "lattice.h"
#include "mesh.h"
#include "utilities/samplingGenerators.h"
#include "functors/analytical/indicator/indicatorF3D.h"
#include "functors/analytical/indicator/indicComb3D.h"
#include "functors/lattice/indicator/superIndicatorF3D.h"

namespace olb {

namespace refinement {

namespace sampling {

/// Cell-centre points on an axis-aligned hyperplane, each level at native resolution.
/// Optional filter keeps only points inside its interior (e.g. the fluid-domain IndicatorCuboid3D).
template <typename T, typename DESCRIPTOR>
std::vector<Vector<T,DESCRIPTOR::d>> slice(RefinedLattice<T,DESCRIPTOR>& lattice,
                                           RefinedMesh<T,DESCRIPTOR::d>& mesh,
                                           const Hyperplane3D<T>& hyperplane,
                                           FunctorPtr<IndicatorF3D<std::type_identity_t<T>>>&& filter = {})
{
  constexpr unsigned D = DESCRIPTOR::d;
  const std::size_t L = lattice.getLevelCount();

  unsigned normal = 0;
  for (unsigned k = 1; k < D; ++k) {
    const T ak = hyperplane.normal[k] < T{0} ? -hyperplane.normal[k] : hyperplane.normal[k];
    const T an = hyperplane.normal[normal] < T{0} ? -hyperplane.normal[normal] : hyperplane.normal[normal];
    if (ak > an) {
      normal = k;
    }
  }
  const T slicePos = hyperplane.origin[normal];

  std::vector<Vector<T,D>> points;
  for (std::size_t iLevel = 0; iLevel < L; ++iLevel) {
    auto& decomp = lattice.getLevel(iLevel).getCuboidDecomposition();
    for (int iC = 0; iC < decomp.size(); ++iC) {
      const auto& cuboid = decomp.get(iC);
      if (cuboid.isCoarseBorder()) {
        continue;
      }
      const auto  o  = cuboid.getOrigin();
      const auto  e  = cuboid.getExtent();
      const T     dx = cuboid.getDeltaR();

      const int kNormal = static_cast<int>(util::floor((slicePos - o[normal]) / dx + T{0.5}));
      if (kNormal < 0 || kNormal >= e[normal]) {
        continue;
      }

      Vector<int,D> cells = e;
      cells[normal] = 1;
      std::size_t total = 1;
      for (unsigned k = 0; k < D; ++k) {
        total *= static_cast<std::size_t>(cells[k]);
      }

      Vector<int,D> idx(0);
      for (std::size_t lin = 0; lin < total; ++lin) {
        Vector<T,D> p;
        for (unsigned k = 0; k < D; ++k) {
          p[k] = (k == normal) ? slicePos : o[k] + idx[k]*dx;
        }
        bool keep = true;
        if (filter) {
          bool inside = false;
          (*filter)(&inside, p.data());
          keep = inside;
        }
        if (keep) {
          points.push_back(p);
        }
        for (unsigned k = 0; k < D; ++k) {
          if (++idx[k] < cells[k]) {
            break;
          }
          idx[k] = 0;
        }
      }
    }
  }
  return points;
}

/// Restrict a communicator stage to the sample regions (plus a stencil margin), per level
template <typename FIELD, typename REFINED_LATTICE, typename REFINED_GEOMETRY, typename STAGE, typename T>
void requestSampledOverlap(REFINED_LATTICE& lattice, REFINED_GEOMETRY& geometry, STAGE stage,
                           const std::vector<olb::sampling::AABB<T,3>>& regions, T margin)
{
  std::shared_ptr<IndicatorF3D<T>> indicator;
  for (const auto& r : regions) {
    std::shared_ptr<IndicatorF3D<T>> box = std::make_shared<IndicatorCuboid3D<T>>(
      r.extent + Vector<T,3>(2*margin), r.origin - Vector<T,3>(margin));
    indicator = indicator ? indicator + box : box;
  }
  for (std::size_t iLevel = 0; iLevel < lattice.getLevelCount(); ++iLevel) {
    auto& comm = lattice.getLevel(iLevel).getCommunicator(stage);
    comm.template requestField<FIELD>();
    std::unique_ptr<SuperIndicatorF<T,3>> restriction =
      std::make_unique<SuperIndicatorFfromIndicatorF3D<T>>(indicator, geometry.getLevel(iLevel));
    comm.requestOverlap(1, std::move(restriction));
    comm.exchangeRequests();
  }
}

}

}

}

#endif
