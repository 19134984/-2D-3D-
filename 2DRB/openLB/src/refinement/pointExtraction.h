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

#ifndef REFINEMENT_POINT_EXTRACTION_H
#define REFINEMENT_POINT_EXTRACTION_H

#include <vector>
#include <memory>
#include <functional>

#include "core/superLatticePointExtraction.h"
#include "utilities/samplingGenerators.h"
#include "mesh.h"
#include "lattice.h"

namespace olb {

namespace refinement {

/// Distributed point extraction across a refined lattice hierarchy
template <typename T, typename DESCRIPTOR, typename FUNCTOR>
class PointExtraction {
public:
  using value_t      = T;
  using descriptor_t = DESCRIPTOR;
  using functor_t    = FUNCTOR;

private:
  using Extractor = SuperLatticePointExtraction<T,DESCRIPTOR,FUNCTOR>;
  using Result = std::decay_t<decltype(std::declval<Extractor>().getValue(0))>;
  using Point  = Vector<T,DESCRIPTOR::d>;
  using Generator = typename Extractor::PointGeneratorF;
  using Region = sampling::AABB<T,DESCRIPTOR::d>;

  std::vector<std::unique_ptr<Extractor>> _extractor;

  std::vector<int>         _levelOfLocal;
  std::vector<std::size_t> _extractorLocalOf;
  std::size_t              _globalSize = 0;

  void build(RefinedLattice<T,DESCRIPTOR>& lattice, RefinedMesh<T,DESCRIPTOR::d>& mesh,
             const Region& region, const Generator& generator)
  {
    const std::size_t L = lattice.getLevelCount();
    auto& base = mesh.getDecompositionBase();

    std::vector<int> baseToLevel(base.size(), -1);
    for (std::size_t iLevel = 0; iLevel < L; ++iLevel) {
      auto& decomp = lattice.getLevel(iLevel).getCuboidDecomposition();
      for (int iC = 0; iC < decomp.size(); ++iC) {
        if (auto ref = decomp.get(iC).getReferenceCuboid()) {
          baseToLevel[*ref] = static_cast<int>(iLevel);
        }
      }
    }

    // Active levels: those whose base cuboids the region overlaps (rank-consistent, no point walk)
    std::vector<bool> active(L, false);
    for (int iC = 0; iC < base.size(); ++iC) {
      const int level = baseToLevel[iC];
      if (level < 0) {
        continue;
      }
      const auto& cuboid = base.get(iC);
      const auto  origin = cuboid.getOrigin();
      const auto  cells  = cuboid.getExtent();
      const T     dx     = cuboid.getDeltaR();
      Point extent;
      for (unsigned k = 0; k < DESCRIPTOR::d; ++k) {
        extent[k] = cells[k] * dx;
      }
      if (region.overlaps(origin, extent)) {
        active[level] = true;
      }
    }

    _extractor.resize(L);
    for (std::size_t iLevel = 0; iLevel < L; ++iLevel) {
      if (!active[iLevel]) {
        continue;
      }
      _extractor[iLevel] = std::make_unique<Extractor>(lattice.getLevel(iLevel),
        [&,iLevel](typename Extractor::IndexedConsumerF emit) {
          std::size_t iG = 0;
          generator([&](const Point& physR) {
            if (auto iCbase = base.getC(physR)) {
              if (baseToLevel[*iCbase] == static_cast<int>(iLevel)) {
                emit(iG, physR);
              }
            }
            ++iG;
          });
        });
    }

    for (std::size_t iLevel = 0; iLevel < L; ++iLevel) {
      if (!_extractor[iLevel]) {
        continue;
      }
      _globalSize += _extractor[iLevel]->getGlobalSize();
      for (std::size_t i = 0; i < _extractor[iLevel]->getLocalSize(); ++i) {
        _levelOfLocal.push_back(static_cast<int>(iLevel));
        _extractorLocalOf.push_back(i);
      }
    }
  }

  static Region boundingBox(const std::vector<Point>& points)
  {
    Region region;
    if (points.empty()) {
      return region;
    }
    Point mn = points[0], mx = points[0];
    for (const auto& p : points) {
      for (unsigned k = 0; k < DESCRIPTOR::d; ++k) {
        mn[k] = util::min(mn[k], p[k]);
        mx[k] = util::max(mx[k], p[k]);
      }
    }
    region.origin = mn;
    region.extent = mx - mn;
    return region;
  }

public:
  static constexpr unsigned d = DESCRIPTOR::d;

  PointExtraction(RefinedLattice<T,DESCRIPTOR>& lattice,
                  RefinedMesh<T,DESCRIPTOR::d>& mesh,
                  const Region& region,
                  const Generator& generator)
  {
    build(lattice, mesh, region, generator);
  }

  PointExtraction(RefinedLattice<T,DESCRIPTOR>& lattice,
                  RefinedMesh<T,DESCRIPTOR::d>& mesh,
                  const std::vector<Point>& points)
  {
    build(lattice, mesh, boundingBox(points), [&points](std::function<void(const Point&)> emit) {
      for (const auto& p : points) {
        emit(p);
      }
    });
  }

  PointExtraction(RefinedLattice<T,DESCRIPTOR>& lattice,
                  RefinedMesh<T,DESCRIPTOR::d>& mesh,
                  const HyperplaneLattice3D<T>& plane)
  {
    build(lattice, mesh, sampling::boundingBox(plane), sampling::pointsOf(plane));
  }

  template <typename FIELD>
  void setParameter(FieldD<T,DESCRIPTOR,FIELD> value) {
    for (auto& extractor : _extractor) {
      if (extractor) {
        auto copy = value;
        extractor->template setParameter<FIELD>(std::move(copy));
      }
    }
  }

  void update() {
    for (auto& extractor : _extractor) {
      if (extractor) {
        extractor->update();
      }
    }
  }

  std::size_t getLocalSize()  const { return _levelOfLocal.size(); }
  std::size_t getGlobalSize() const { return _globalSize; }

  std::size_t getGlobalIndex(std::size_t iLocal) const {
    const int level = _levelOfLocal[iLocal];
    return _extractor[level]->getGlobalIndex(_extractorLocalOf[iLocal]);
  }

  Vector<T,DESCRIPTOR::d> getPhysR(std::size_t iLocal) const {
    const int level = _levelOfLocal[iLocal];
    return _extractor[level]->getPhysR(_extractorLocalOf[iLocal]);
  }

  Result getValue(std::size_t iLocal) const {
    const int level = _levelOfLocal[iLocal];
    return _extractor[level]->getValue(_extractorLocalOf[iLocal]);
  }

  int getPointLevel(std::size_t iLocal) const {
    return _levelOfLocal[iLocal];
  }

};

}

}

#endif
