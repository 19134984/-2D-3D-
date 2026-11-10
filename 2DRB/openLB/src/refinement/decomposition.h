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

#ifndef REFINEMENT_DECOMPOSITION_H
#define REFINEMENT_DECOMPOSITION_H

#include "geometry/cuboidDecomposition.h"

namespace olb::refinement {

template <typename T, unsigned D>
std::vector<CuboidDecomposition<T,D>> resolveMixedCuboidDecomposition(
  CuboidDecomposition<T,D>& mixedDecomposition)
{
  std::map<T, std::vector<Cuboid<T,D>>> blueprint;
  std::set<T> deltas;

  for (int i = 0; i < mixedDecomposition.size(); ++i) {
    auto& cuboid = mixedDecomposition.get(i);
    cuboid.setReferenceCuboid(i);
    blueprint[cuboid.getDeltaR()].push_back(cuboid);
    deltas.insert(cuboid.getDeltaR());
  }

  std::vector<T> sortedDeltas(deltas.begin(), deltas.end());
  std::sort(sortedDeltas.rbegin(), sortedDeltas.rend());

  std::size_t nLevels = sortedDeltas.size();

  std::vector<CuboidDecomposition<T,D>> overlappingLevels;

  for (std::size_t i = 0; i < nLevels; ++i) {
    T currDelta = sortedDeltas[i];
    const auto& currBlueprint = blueprint[currDelta];

    auto levelMother = mixedDecomposition.getMotherCuboid();
    for (int j=0; j < i; ++j) {
      levelMother.refineCellCentered(2);
    }
    CuboidDecomposition<T,D> newDecomp(levelMother, 1);
    newDecomp.cuboids().clear();

    std::vector<Cuboid<T,D>> cuboidsForThisLevel;

    for (const auto& cub : currBlueprint) {
      cuboidsForThisLevel.push_back(cub);
    }

    if (i < nLevels - 1) {
      const auto& finerBlueprint = blueprint[sortedDeltas[i+1]];
      int coarsenFactor = static_cast<int>(std::round(currDelta / sortedDeltas[i+1]));

      for (const auto& fineCuboid : finerBlueprint) {
        bool isNeighbor = false;
        for (const auto& currCuboid : currBlueprint) {
          if (currCuboid.intersects(fineCuboid, 1)) {
            isNeighbor = true;
            break;
          }
        }
        if (isNeighbor) {
          auto coarsenedNeighbor = fineCuboid;
          coarsenedNeighbor.coarsenCellCentered(coarsenFactor);
          if (std::find(cuboidsForThisLevel.begin(),
                        cuboidsForThisLevel.end(),
                        coarsenedNeighbor) == cuboidsForThisLevel.end()) {
            cuboidsForThisLevel.push_back(coarsenedNeighbor);
          }
        }
      }
    }

    newDecomp.cuboids() = std::move(cuboidsForThisLevel);
    overlappingLevels.push_back(std::move(newDecomp));
  }

  return overlappingLevels;
}

}

#endif
