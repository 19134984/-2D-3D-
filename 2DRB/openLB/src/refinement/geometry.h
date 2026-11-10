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

#ifndef REFINEMENT_GEOMETRY_H
#define REFINEMENT_GEOMETRY_H

#include "mesh.h"
#include "geometry/superGeometry.h"

#include "utilities/log/log.h"

namespace olb {

template <typename T, unsigned D>
class RefinedGeometry {
private:
  std::vector<std::unique_ptr<SuperGeometry<T, D>>> _geometries;

public:
  RefinedGeometry(RefinedMesh<T, D>& refinedMesh)
  {
    OLB_SCOPE(logging::Setup, "construct refined geometry");
    for (std::size_t iLevel = 0; iLevel < refinedMesh.getLevelCount(); ++iLevel) {
      auto& mesh = refinedMesh.getLevel(iLevel);
      _geometries.push_back(std::make_unique<SuperGeometry<T, D>>(
        mesh.getCuboidDecomposition(),
        mesh.getLoadBalancer(),
        mesh.getOverlap()
      ));
    }
  }

  SuperGeometry<T, D>& getLevel(std::size_t iLevel)
  {
    return *_geometries[iLevel];
  }

  const SuperGeometry<T, D>& getLevel(std::size_t iLevel) const
  {
    return *_geometries[iLevel];
  }

  std::size_t getLevelCount() const
  {
    return _geometries.size();
  }

  /// Reset coarse-border blocks to an inactive material on every level
  void markCoarseBorders(int material)
  {
    for (std::size_t iLevel = 0; iLevel < getLevelCount(); ++iLevel) {
      auto& geometry = getLevel(iLevel);
      auto& decomp = geometry.getCuboidDecomposition();
      for (int iC = 0; iC < geometry.getLoadBalancer().size(); ++iC) {
        if (decomp.get(geometry.getLoadBalancer().glob(iC)).isCoarseBorder()) {
          geometry.getBlock(iC).reset(material);
        }
      }
    }
  }
};

}

#endif
