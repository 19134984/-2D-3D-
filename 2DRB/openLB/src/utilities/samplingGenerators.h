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

#ifndef UTILITIES_SAMPLING_GENERATORS_H
#define UTILITIES_SAMPLING_GENERATORS_H

#include <vector>
#include <functional>

#include "core/vector.h"
#include "hyperplaneLattice3D.h"

namespace olb {

namespace sampling {

/// Streaming point generator: visit each sampled point once via a callback
template <typename T, unsigned D>
using PointGenerator = std::function<void(std::function<void(const Vector<T,D>&)>)>;

/// Axis-aligned bounding box of a sample region (a zero-extent axis is a plane)
template <typename T, unsigned D>
struct AABB {
  Vector<T,D> origin;
  Vector<T,D> extent;

  bool overlaps(Vector<T,D> boxOrigin, Vector<T,D> boxExtent) const {
    for (unsigned k=0; k<D; ++k) {
      if (origin[k] + extent[k] < boxOrigin[k] || boxOrigin[k] + boxExtent[k] < origin[k]) {
        return false;
      }
    }
    return true;
  }
};

/// Regular structured region (box/plane): origin + per-axis extent
/**
 * A zero extent on one axis degenerates the box to a plane on that axis. The node
 * count per axis is extent/(spacing*stride)+1 (>=1); axis 0 varies fastest.
 **/
template <typename T, unsigned D>
struct StructuredRegion {
  Vector<T,D> origin;
  Vector<T,D> extent;
  T           spacing;
  int         stride;

  Vector<int,D> nodeCount() const {
    Vector<int,D> n;
    const T step = spacing * stride;
    for (unsigned k=0; k<D; ++k) {
      n[k] = (extent[k] > T{0}) ? static_cast<int>(extent[k]/step + T{0.5}) + 1 : 1;
    }
    return n;
  }
};

/// Explicit point list (scattered / random probes)
template <typename T, unsigned D>
PointGenerator<T,D> points(const std::vector<Vector<T,D>>& pts) {
  return [pts](std::function<void(const Vector<T,D>&)> emit) {
    for (const auto& p : pts) {
      emit(p);
    }
  };
}

/// Single monitor point
template <typename T, unsigned D>
PointGenerator<T,D> point(const Vector<T,D>& p) {
  return [p](std::function<void(const Vector<T,D>&)> emit) {
    emit(p);
  };
}

/// Structured box / volume, optionally strided (coarsened)
template <typename T, unsigned D>
PointGenerator<T,D> box(const StructuredRegion<T,D>& region) {
  return [region](std::function<void(const Vector<T,D>&)> emit) {
    const Vector<int,D> n = region.nodeCount();
    const T step = region.spacing * region.stride;
    std::size_t total = 1;
    for (unsigned k=0; k<D; ++k) {
      total *= static_cast<std::size_t>(n[k]);
    }
    Vector<int,D> idx(0);
    for (std::size_t lin=0; lin<total; ++lin) {
      Vector<T,D> p = region.origin;
      for (unsigned k=0; k<D; ++k) {
        p[k] += idx[k] * step;
      }
      emit(p);
      for (unsigned k=0; k<D; ++k) {           // odometer, axis 0 fastest
        if (++idx[k] < n[k]) break;
        idx[k] = 0;
      }
    }
  };
}

/// Convenience: box from origin + extent (+ optional spacing/stride)
template <typename T, unsigned D>
PointGenerator<T,D> box(const Vector<T,D>& origin, const Vector<T,D>& extent,
                        T spacing, int stride = 1) {
  return box<T,D>(StructuredRegion<T,D>{origin, extent, spacing, stride});
}

/// Plane / slice: a structured region with one zero extent (no stride coarsening)
template <typename T, unsigned D>
PointGenerator<T,D> plane(const Vector<T,D>& origin, const Vector<T,D>& extent, T spacing) {
  return box<T,D>(StructuredRegion<T,D>{origin, extent, spacing, 1});
}

/// Cell-centred bounded plane lattice from a corner origin and two in-plane extent vectors
template <typename T>
HyperplaneLattice3D<T> planeLattice(Vector<T,3> origin, Vector<T,3> uExtent, Vector<T,3> vExtent, T dx) {
  const int nx = static_cast<int>(norm(uExtent)/dx + T{0.5});
  const int ny = static_cast<int>(norm(vExtent)/dx + T{0.5});
  const Vector<T,3> shift = (normalize(uExtent) + normalize(vExtent)) * (T{0.5}*dx);
  Hyperplane3D<T> hyperplane;
  hyperplane.originAt(origin + shift).spannedBy(uExtent, vExtent);
  return HyperplaneLattice3D<T>(hyperplane, dx, nx, ny);
}

/// Axis-aligned convenience: extent has one zero axis (the normal); the other two span u and v
template <typename T>
HyperplaneLattice3D<T> planeLattice(Vector<T,3> origin, Vector<T,3> extent, T dx) {
  Vector<T,3> uExtent{}, vExtent{};
  bool haveU = false;
  for (unsigned k=0; k<3; ++k) {
    if (extent[k] != T{0}) {
      (haveU ? vExtent : uExtent)[k] = extent[k];
      haveU = true;
    }
  }
  return planeLattice(origin, uExtent, vExtent, dx);
}

/// Streaming generator over a plane lattice's points (u fastest)
template <typename T>
PointGenerator<T,3> pointsOf(HyperplaneLattice3D<T> lattice) {
  return [lattice](std::function<void(const Vector<T,3>&)> emit) {
    for (int iy=0; iy < lattice.getNy(); ++iy) {
      for (int ix=0; ix < lattice.getNx(); ++ix) {
        emit(lattice.getPhysR(ix, iy));
      }
    }
  };
}

/// Axis-aligned bounding box of a plane lattice
template <typename T>
AABB<T,3> boundingBox(const HyperplaneLattice3D<T>& lattice) {
  const Vector<T,3> o = lattice.getPhysOrigin();
  const Vector<T,3> u = lattice.getVectorU() * T(lattice.getNx());
  const Vector<T,3> v = lattice.getVectorV() * T(lattice.getNy());
  Vector<T,3> mn = o, mx = o;
  for (const auto& c : {o+u, o+v, o+u+v}) {
    for (unsigned k=0; k<3; ++k) {
      if (c[k] < mn[k]) { mn[k] = c[k]; }
      if (c[k] > mx[k]) { mx[k] = c[k]; }
    }
  }
  return {mn, mx - mn};
}

}

}

#endif
