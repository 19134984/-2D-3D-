/*  This file is part of the OpenLB library
*
 *  Copyright (C) 2026 Kilian Mahler
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

#ifndef OLB_REFINEMENT_OCTREE_OCTANT_H
#define OLB_REFINEMENT_OCTREE_OCTANT_H

#include <cstdint>
#include <array>
#include <vector>
#include <optional>
#include <list>
#include <cmath>
#include <algorithm>
#include <bit>

#include "core/vector.h"
#include "core/util.h"

namespace olb {

static const std::array<Vector<int, 3>, 8> neighborLookup = {
    {{{-1, -1, -1}},
     {{0, -1, -1}},
     {{-1, 0, -1}},
     {{0, 0, -1}},
     {{-1, -1, 0}},
     {{0, -1, 0}},
     {{-1, 0, 0}},
     {{0, 0, 0}}}
};

// based on : "Dynamic Octree Load Balancing Using Space-Filling Curves", Campbell, Devine, Gervasio, Teresco, 2003
// just translated, so we use our child indices and start with different orientation
static const std::array< std::array< int , 8 >, 24 > hilbertOrderLookup { { { 0,1,3,2,7,6,4,5},
{ 0,7,3,4,1,6,2,5},
{ 0,1,7,6,3,2,4,5},
{ 2,1,5,6,3,0,4,7},
{ 4,7,3,0,5,6,2,1},
{ 4,5,3,2,7,6,0,1},
{ 6,1,5,2,7,0,4,3},
{ 0,3,7,4,1,2,6,5},
{ 2,3,5,4,1,0,6,7},
{ 6,7,5,4,1,0,2,3},
{ 0,3,1,2,7,4,6,5},
{ 2,5,3,4,1,6,0,7},
{ 4,3,5,2,7,0,6,1},
{ 4,3,7,0,5,2,6,1},
{ 6,5,7,4,1,2,0,3},
{ 0,7,1,6,3,4,2,5},
{ 2,5,1,6,3,4,0,7},
{ 6,5,1,2,7,4,0,3},
{ 2,3,1,0,5,4,6,7},
{ 4,5,7,6,3,2,0,1},
{ 4,7,5,6,3,0,2,1},
{ 6,7,1,0,5,4,2,3},
{ 2,1,3,0,5,6,4,7},
{ 6,1,7,0,5,2,4,3} } };


static const std::array< std::array< int , 8 >, 24 > hilbertOrientation { { {  1,  2,  0,  3,  4,  0,  5,  6 },
                                                  {  0,  7,  1,  8,  5,  1,  4,  9 },
                                                  { 15,  0,  2, 22, 20,  2, 19, 23 },
                                                  { 20,  6,  3, 23, 15,  3, 16, 22 },
                                                  { 22, 13,  4, 12, 11,  4,  1, 20 },
                                                  { 11, 19,  5, 20, 22,  5,  0, 12 },
                                                  {  9,  3,  6,  2, 21,  6, 17,  0 },
                                                  { 10,  1,  7, 11, 12,  7, 13, 14 },
                                                  { 12,  9,  8, 14, 10,  8, 18, 11 },
                                                  {  6,  8,  9,  7, 17,  9, 21,  1 },
                                                  {  7, 15, 10, 16, 13, 10, 12, 17 },
                                                  {  5, 14, 11,  9,  0, 11, 22,  8 },
                                                  {  8, 20, 12, 19, 18, 12, 10,  5 },
                                                  { 18,  4, 13,  5,  8, 13,  7, 19 },
                                                  { 17, 11, 14,  1,  6, 14, 23,  7 },
                                                  {  2, 10, 15, 18, 19, 15, 20, 21 },
                                                  { 19, 17, 16, 21,  2, 16,  3, 18 },
                                                  { 14, 16, 17, 15, 23, 17,  6, 10 },
                                                  { 13, 21, 18, 17,  7, 18,  8, 16 },
                                                  { 16,  5, 19,  4,  3, 19,  2, 13 },
                                                  {  3, 12, 20, 13, 16, 20, 15,  4 },
                                                  { 23, 18, 21, 10, 14, 21,  9, 15 },
                                                  {  4, 23, 22,  6,  1, 22, 11,  3 },
                                                  { 21, 22, 23,  0,  9, 23, 14,  2 } } };



static const std::array<std::array<int, 4>, 4> hilbertOrderLookup2D {{
  {0, 3, 1, 2},
  { 0, 1, 3, 2 },
  { 2, 3, 1, 0},
  { 2, 1, 3, 0 }
}};


static const std::array<std::array<int, 4>, 4> hilbertOrientation2D {{
  { 1, 0, 0, 2},
  { 0, 1, 1, 3 },
  { 3, 2, 2, 0 },
  { 2, 3, 3, 1 }
}};


template <typename T, unsigned D>
class Octant {
private:

  Vector<int, D> _origin;

  int _refinementLevel;

  int _maxRefinement;

  int _extent;

  std::optional<int> process;

  std::optional<int> _globalID;

  bool _coarseBorderOctant;

  std::optional<int> _domainLevel;

  std::optional<std::uint64_t> _weight;

public:
  Octant(Vector<std::uint64_t, D> origin, int refinementLevel, int maxRefinement, std::optional<std::uint64_t> globalID = std::nullopt )
  {
    _origin             = origin;
    _refinementLevel    = refinementLevel;
    _maxRefinement      = maxRefinement;
    _extent             = util::pow(2, maxRefinement - refinementLevel);
    _coarseBorderOctant = false;
    process = std::nullopt;
    _domainLevel = std::nullopt;
    _weight = std::nullopt;
    _globalID = globalID;
  }

  void incDomainLevel()
  {
    if (!_domainLevel.has_value()) {
      _domainLevel = std::make_optional(_refinementLevel - 1);
    }
    else {
      _domainLevel = std::make_optional(_domainLevel.value() - 1);
    }
  }


  void addCuboid(std::vector<std::vector<Cuboid<T, D>>>& result, Vector<T, D> rootOrigin, T coarseVoxelSize,
                 int voxelPerCuboid)
  {
    if (_domainLevel.has_value()) {
      addCoarsendCuboid(result, rootOrigin, coarseVoxelSize, voxelPerCuboid);
      return;
    }
    T   delta      = coarseVoxelSize;
    T   localVoxel = coarseVoxelSize;
    int mask       = 0b1;
    mask <<= _maxRefinement;
    T            cubeLength  = coarseVoxelSize * voxelPerCuboid;
    Vector<T,D> innerOrigin;
    for (int i = 0; i < D; i ++) {
      innerOrigin[i] = (_origin[i] & mask ? 1: 0) * cubeLength;
    }


    for (int i = 0; i < _refinementLevel; i++) {
      delta /= 2;
      localVoxel /= 2;
      innerOrigin -= delta / 2;
      mask >>= 1;
      cubeLength  = localVoxel * voxelPerCuboid;
      for (int i = 0; i < D; i ++) {
        innerOrigin[i] = innerOrigin[i] + (_origin[i] & mask ? 1: 0) * cubeLength;
      }

    }

    result[_refinementLevel].emplace_back(Cuboid<T, D>(rootOrigin + innerOrigin, localVoxel, voxelPerCuboid, process));
    result[_refinementLevel].back().setWeight(calcWeight());
    if (_coarseBorderOctant) {

      innerOrigin += delta / 2;
      localVoxel *= 2;
      result[_refinementLevel - 1].emplace_back(
          Cuboid<T, D>(rootOrigin + innerOrigin, localVoxel, voxelPerCuboid / 2, process, true));
    }
  }


  void addCoarsendCuboid(std::vector<std::vector<Cuboid<T, D>>>& result, Vector<T, D> rootOrigin, T coarseVoxelSize,
                    int voxelPerCuboid)
  {
    T   delta      = coarseVoxelSize;
    T   localVoxel = coarseVoxelSize;
    int mask       = 0b1;
    mask <<= _maxRefinement;
    T            cubeLength  = coarseVoxelSize * voxelPerCuboid;
    Vector<T,D> innerOrigin;
    for (int i = 0; i < D; i ++) {
      innerOrigin[i] = (_origin[i] & mask ? 1: 0) * cubeLength;
    }

    for (int i = 0; i < _refinementLevel; i++) {
      delta /= 2;
      localVoxel /= 2;
      innerOrigin -= delta / 2;
      mask >>= 1;
      cubeLength  = localVoxel * voxelPerCuboid;
      for (int i = 0; i < D; i ++) {
        innerOrigin[i] = innerOrigin[i] + (_origin[i] & mask ? 1: 0) * cubeLength;
      }
    }

    result[_refinementLevel].emplace_back(Cuboid<T, D>(
        rootOrigin + innerOrigin, localVoxel, voxelPerCuboid * util::pow(2, _refinementLevel - _domainLevel.value()), process));
    result[_refinementLevel].back().setWeight(calcWeight());
    if (_coarseBorderOctant) {
      innerOrigin += delta / 2;
      localVoxel *= 2;

      result[_refinementLevel - 1].emplace_back(Cuboid<T, D>(
          rootOrigin + innerOrigin, localVoxel, voxelPerCuboid * util::pow(2, _refinementLevel - 1 - _domainLevel.value()), process, true));
    }
  }

  void addCuboidFine(std::vector<std::vector<Cuboid<T, D>>>& result, Vector<T, D> rootOrigin, T coarseVoxelSize,
                     int voxelPerCuboid)
  {
    if (!_coarseBorderOctant) {
      return;
    }
    if (_domainLevel.has_value()) {
      T   delta      = coarseVoxelSize;
      T   localVoxel = coarseVoxelSize;
      int mask       = 0b1;
      mask <<= _maxRefinement;
      T            cubeLength  = coarseVoxelSize * voxelPerCuboid;
      Vector<T, D> innerOrigin;
      for (int i = 0; i < D; i++) {
        innerOrigin[i] = ((_origin[i] & mask) ? 1 : 0) * cubeLength;
      }
      for (int i = 0; i < _refinementLevel; i++) {
        delta /= 2;
        localVoxel /= 2;
        innerOrigin -= delta / 2;
        mask >>= 1;
        cubeLength  = localVoxel * voxelPerCuboid;
        for (int j = 0; j < D; j++) {
          innerOrigin[j] = innerOrigin[j] + ((_origin[j] & mask) ? 1 : 0) * cubeLength;
        }
      }
      innerOrigin += delta / 2;
      localVoxel *= 2;

      result[_refinementLevel].emplace_back(Cuboid<T, D>(
          rootOrigin + innerOrigin, localVoxel, voxelPerCuboid * util::pow(2, _refinementLevel - 1 - _domainLevel.value()), process));
      result[_refinementLevel].back().setWeight(calcWeight());
      return;
    }
    T   delta      = coarseVoxelSize;
    T   localVoxel = coarseVoxelSize;
    int mask       = 0b1;
    mask <<= _maxRefinement;
    T            cubeLength  = coarseVoxelSize * voxelPerCuboid;
    Vector<T, D> innerOrigin;
    for (int i = 0; i < D; i++) {
      innerOrigin[i] = ((_origin[i] & mask) ? 1 : 0) * cubeLength;
    }

    for (int i = 0; i < _refinementLevel; i++) {
      delta /= 2;
      localVoxel /= 2;
      innerOrigin -= delta / 2;
      mask >>= 1;
      cubeLength  = localVoxel * voxelPerCuboid;
      for (int j = 0; j < D; j++) {
        innerOrigin[j] = innerOrigin[j] + ((_origin[j] & mask) ? 1 : 0) * cubeLength;
      }
    }
    innerOrigin += delta / 2;
    localVoxel *= 2;

    result[_refinementLevel].emplace_back(
        Cuboid<T, D>(rootOrigin + innerOrigin, localVoxel, voxelPerCuboid / 2, process));
    result[_refinementLevel].back().setWeight(calcWeight());
  }


  bool intersect(IndicatorF<T, D>& indicatorF, Vector<T, D> rootOrigin, T coarseVoxelSize, int voxelPerCuboid,
                      int padding)
  {
    T   delta      = coarseVoxelSize;
    T   localVoxel = coarseVoxelSize;
    int mask       = 0b1;
    mask <<= _maxRefinement;
    T            cubeLength  = coarseVoxelSize * voxelPerCuboid;
    Vector<T,D> innerOrigin;
    for (int i = 0; i < D; i++) {
      innerOrigin[i] = (_origin[i] & mask ? 1:0) * cubeLength;
    }

    for (int i = 0; i < _refinementLevel; i++) {
      delta /= 2;
      localVoxel /= 2;
      innerOrigin -= delta / 2;
      mask >>= 1;
      cubeLength  = localVoxel * voxelPerCuboid;
      for (int i = 0; i < D; i++) {
        innerOrigin[i] = innerOrigin[i] + ((_origin[i] & mask) ? 1 : 0) * cubeLength;
      }
    }

    const Vector<T, D> base = rootOrigin + innerOrigin;

    bool inside = false;
    T    point[D];
    int  finalPadding = _refinementLevel + 1 == _maxRefinement ? padding : 1;
    for (int iX = 0 - finalPadding; iX < voxelPerCuboid + finalPadding; ++iX) {
      for (int iY = 0 - finalPadding; iY < voxelPerCuboid + finalPadding; ++iY) {
        if constexpr (D == 2) {
          point[0] = base[0] + iX * localVoxel;
          point[1] = base[1] + iY * localVoxel;

          indicatorF(&inside, point);
          if (inside)
            return true;
        }
        else {
          for (int iZ = 0 - finalPadding; iZ < voxelPerCuboid + finalPadding; ++iZ) {

            point[0] = base[0] + iX * localVoxel;
            point[1] = base[1] + iY * localVoxel;
            point[2] = base[2] + iZ * localVoxel;

            indicatorF(&inside, point);
            if (inside)
              return true;
          }
        }
      }
    }

    return false;
  }

  void reduceWeights(IndicatorF<T, D>& indicatorF, Vector<T, D> rootOrigin, T coarseVoxelSize, int voxelPerCuboid)  {
    T   delta      = coarseVoxelSize;
    T   localVoxel = coarseVoxelSize;
    int mask       = 0b1;
    mask <<= _maxRefinement;
    T            cubeLength  = coarseVoxelSize * voxelPerCuboid;
    Vector<T,D> innerOrigin;
    for (int i = 0; i < D; i++) {
      innerOrigin[i] = (_origin[i] & mask ? 1:0) * cubeLength;
    }

    for (int i = 0; i < _refinementLevel; i++) {
      delta /= 2;
      localVoxel /= 2;
      innerOrigin -= delta / 2;
      mask >>= 1;
      cubeLength  = localVoxel * voxelPerCuboid;
      for (int i = 0; i < D; i++) {
        innerOrigin[i] = innerOrigin[i] + ((_origin[i] & mask) ? 1 : 0) * cubeLength;
      }
    }

    const Vector<T, D> base = rootOrigin + innerOrigin;
    std::uint64_t countInside = 0;
    std::uint64_t total = util::pow(_extent,D);

    bool inside = false;
    T    point[D];
    for (int iX = 0; iX < voxelPerCuboid; ++iX) {
      for (int iY = 0; iY < voxelPerCuboid; ++iY) {
        if constexpr (D == 2) {
          point[0] = base[0] + iX * localVoxel;
          point[1] = base[1] + iY * localVoxel;

          indicatorF(&inside, point);
          if (inside) {
            countInside++;
          }
        }
        else {
          for (int iZ = 0; iZ < voxelPerCuboid; ++iZ) {

            point[0] = base[0] + iX * localVoxel;
            point[1] = base[1] + iY * localVoxel;
            point[2] = base[2] + iZ * localVoxel;

            indicatorF(&inside, point);
            if (inside) {
              countInside++;
            }
          }
        }
      }
    }
    if (countInside == total) {
      _weight = std::make_optional(0);
      return;
    }
    if (total/(total-countInside) == 0) {
      _weight = std::make_optional(0);
      return;
    }
    std::uint64_t weight = calcWeight();
    _weight = std::make_optional(weight / (total/(total-countInside)));
  }

  std::list<Octant<T, D>> refine()
  {
    if (_refinementLevel == _maxRefinement) {
      throw std::logic_error("Maximum Refinement Level reached");
    }
    std::list<Octant<T, D>> childs;
    int                     childExtend = _extent / 2;
    if constexpr (D== 2) {
      for (int iY = 0; iY < 2; iY++) {
        for (int iX = 0; iX < 2; iX++) {
          Vector<T, D> offset(childExtend * iX, childExtend * iY);
          childs.push_back(Octant<T, D>(_origin + offset, _refinementLevel + 1, _maxRefinement));
        }
      }
    }
    else {
      for (int iZ = 0; iZ < 2; iZ++) {
        for (int iY = 0; iY < 2; iY++) {
          for (int iX = 0; iX < 2; iX++) {
            Vector<T, D> offset(childExtend * iX, childExtend * iY, childExtend * iZ);
            childs.push_back(Octant<T, D>(_origin + offset, _refinementLevel + 1, _maxRefinement));
          }
        }
      }
    }
    return childs;
  }


  int getRefinementLevel() const { return _refinementLevel; }


  int getDomainLevel() const
  {
    if (!_domainLevel.has_value()) {
      return _refinementLevel;
    }
    return _domainLevel.value();
  }

  std::uint64_t getId() const { return calcIDFast(_origin); }

  std::uint64_t getGlobalID() const
  {
    return _globalID.value();
  }

  Vector<int, D> getOrigin() const { return _origin; }

  static Vector<int, D> calcOrigin(std::uint64_t iC)
  {
    int x = 0;
    int y = 0;
    int z = 0;
    int k = 0;
    while (iC) {
      x  = x | ((iC & 0b1) << k);
      iC = iC >> 1;
      y  = y | ((iC & 0b1) << k);
      iC = iC >> 1;
      z  = z | ((iC & 0b1) << k);
      iC = iC >> 1;
      k++;
    }

    return Vector<int, D>(x, y, z);
  }

  static Vector<int, D> calcOriginFast(std::uint64_t id)
  {
    if constexpr (D == 2) {
      int x = compactBy2(id);
      int y = compactBy2(id >> 1);
      return {x,y};
    }
    else {
      int x = compactBy3(id);
      int y = compactBy3(id >> 1);
      int z = compactBy3(id >> 2);
      return {x, y, z};
    }
  }

  static std::uint64_t compactBy2(std::uint64_t a)
  {
    a &= 0x5555555555555555ULL;
    a = (a ^ (a >> 1))  & 0x3333333333333333ULL;
    a = (a ^ (a >> 2))  & 0x0f0f0f0f0f0f0f0fULL;
    a = (a ^ (a >> 4))  & 0x00ff00ff00ff00ffULL;
    a = (a ^ (a >> 8))  & 0x0000ffff0000ffffULL;
    a = (a ^ (a >> 16)) & 0x00000000ffffffffULL;
    return a;
  }

  static std::uint64_t compactBy3(std::uint64_t a)
  {
    a &= 0x1249249249249249ULL;
    a = (a ^ (a >> 2))  & 0x10c30c30c30c30c3ULL;
    a = (a ^ (a >> 4))  & 0x100f00f00f00f00fULL;
    a = (a ^ (a >> 8))  & 0x1f0000ff0000ffULL;
    a = (a ^ (a >> 16)) & 0x1f00000000ffffULL;
    a = (a ^ (a >> 32)) & 0x1fffffULL;
    return a;
  }


  /// calculates the id of the octant in refence to sibling octants
  int calcChildID() const
  {
    int i    = 0;
    int mask = 0b1 << (_maxRefinement - _refinementLevel);
    i |= _origin[0] & mask ? 1 : 0;
    i |= _origin[1] & mask ? 2 : 0;
    if constexpr (D == 3) {
      i |= _origin[2] & mask ? 4 : 0;
    }
    return i;
  }

  int calcHilbertChildID() const
  {
    int i    = 0;
    int mask = 0b1 << (_maxRefinement - _refinementLevel);
    i |= _origin[0] & mask ? 0 : 1;
    i |= _origin[1] & mask ? 0 : 2;
    if constexpr (D == 3) {
      i |= _origin[2] & mask ? 0 : 4;
    }
    return i;
  }

  /// calculates the origin of the parent octant
  Vector<int, D> getParentOrigin() const
  {
    int mask = 0b1 << (_maxRefinement - _refinementLevel);
    mask     = ~mask;
    if constexpr (D == 2) {
      return Vector<int, D>(_origin[0] & mask, _origin[1] & mask);
    }
    else {
      return Vector<int, D>(_origin[0] & mask, _origin[1] & mask, _origin[2] & mask);
    }
  }


  Vector<int, D> getParentOriginDomainLevel() const
  {
    if (!_domainLevel.has_value()) {
      return getParentOrigin();
    }
    int mask = 0b1 << (_maxRefinement - _domainLevel.value());
    mask     = ~mask;
    if constexpr (D == 2) {
      return Vector<int, D>(_origin[0] & mask, _origin[1] & mask);
    }
    else {
    return Vector<int, D>(_origin[0] & mask, _origin[1] & mask, _origin[2] & mask);
    }
  }

  /// calculates id of an Origin
  static std::uint64_t calcIDFast(Vector<int, D> origin)
  {
    if constexpr (D == 2) {
      std::uint64_t x = origin[0];
      std::uint64_t y = origin[1];
      return splitBy2(y) << 1 | splitBy2(x);
    }
    else {
      std::uint64_t x = origin[0];
      std::uint64_t y = origin[1];
      std::uint64_t z = origin[2];
      return splitBy3(z) << 2 | splitBy3(y) << 1 | splitBy3(x);
    }
  }

  /// splits uint64 into uint64 with 2x 0 between each bit
  static std::uint64_t splitBy3(std::uint64_t a)
  {
    a &= 0x1fffff; // only low 21 bits
    a = (a | a << 32) & 0x1f00000000ffff;
    a = (a | a << 16) & 0x1f0000ff0000ff;
    a = (a | a << 8) & 0x100f00f00f00f00f;
    a = (a | a << 4) & 0x10c30c30c30c30c3;
    a = (a | a << 2) & 0x1249249249249249;
    return a;
  }

  static std::uint64_t splitBy2(std::uint64_t a)
  {
    a &= 0x00000000ffffffffULL; // only low 32 bits
    a = (a | a << 16) & 0x0000ffff0000ffffULL;
    a = (a | a << 8)  & 0x00ff00ff00ff00ffULL;
    a = (a | a << 4)  & 0x0f0f0f0f0f0f0f0fULL;
    a = (a | a << 2)  & 0x3333333333333333ULL;
    a = (a | a << 1)  & 0x5555555555555555ULL;
    return a;
  }


  void setRefinementBorder() { _coarseBorderOctant = true; }


  bool isCoarseBorderOctant() { return _coarseBorderOctant; }


  std::uint64_t calcWeight() const
  {
    if (_weight.has_value()) {
      return _weight.value();
    }
    if (_domainLevel.has_value()) {
      std::uint64_t w = 0;
      if (_coarseBorderOctant && _refinementLevel >= (int)D) {
        w += (0b1ULL << _refinementLevel) + (0b1ULL << (_refinementLevel - (int)D));
      }
      else {
        w += (0b1ULL << _refinementLevel);
      }
      for (int i = _domainLevel.value(); i < _refinementLevel; i++) {
        w = w << D;
      }
      return w;
    }

    if (_coarseBorderOctant) {
      if (_refinementLevel >= (int)D) {
        return (0b1ULL << _refinementLevel) + (0b1ULL << (_refinementLevel - (int)D));
      }
      return (0b1ULL << _refinementLevel);
    }
    return 0b1ULL << _refinementLevel;
  }

  std::uint64_t calcEdgeWeight(Octant<T,D>& other, bool discount, bool remoteRoot = false) const
  {
    int ref = std::min(other.getRefinementLevel(), _refinementLevel);

    if (!discount) {
      return 0b1ULL << ref;
    }
    if (remoteRoot)
      return 1;

    std::uint64_t diff = getId() ^ other.getId();;

    diff <<= (64 - _maxRefinement * D);
    if (diff == 0)
      return 0b1ULL << ref * 2;

    int leading = std::countl_zero(diff);
    int sharedLevels = leading / D;

    return 0b1ULL << (sharedLevels * 2);
  }

  void setProcess(int p) { process = std::make_optional(p); }

  static bool precludeEqual(const Octant<T, D>& first, const Octant<T, D>& second)
  {
    std::uint64_t mask = ~0ULL << ((first._maxRefinement + 1 - first.getRefinementLevel()) * D);
    return (first.getRefinementLevel() <= second.getRefinementLevel()) &&
           (first.getId() & mask) == (second.getId() & mask);
  }

  static bool precludeTrue(Octant<T, D> first, Octant<T, D> second)
  {
    std::uint64_t mask = ~0ULL << ((first._maxRefinement + 1 - first.getRefinementLevel()) * D);
    return (first.getRefinementLevel() < second.getRefinementLevel()) &&
           (first.getId() & mask) == (second.getId() & mask);
  }

  static bool precludeTotalOrder(const Octant<T, D>& lhs, const Octant<T, D>& rhs)
  {
    std::uint64_t mask = ~0ULL << ((lhs._maxRefinement + 1 - lhs.getRefinementLevel()) * D);
    return ((lhs.getRefinementLevel() < rhs.getRefinementLevel()) && (lhs.getId() & mask) == (rhs.getId() & mask)) ||
           lhs.getId() < rhs.getId();
  }

  static bool idComp(const Octant<T, D>& lhs, const Octant<T, D>& rhs) { return lhs.getId() < rhs.getId(); }

  std::vector<Octant<T, D>> getCoarseLocalNeighbors() const
  {
    if (_refinementLevel == 0) {
      return {};
    }
    std::vector<Octant<T, D>> neighbors;
    Vector<int, D>            parentO      = getParentOrigin();
    const int                 parentLength = 1 << (_maxRefinement - _refinementLevel + 1);
    const int                 maxLength    = 1 << _maxRefinement;

    int lowx  = parentO[0] == 0 ? 0 : -1;
    int highx = parentO[0] + parentLength == maxLength ? 0 : 1;
    int lowy  = parentO[1] == 0 ? 0 : -1;
    int highy = parentO[1] + parentLength == maxLength ? 0 : 1;
    if constexpr (D == 3) {
      int lowz  = parentO[2] == 0 ? 0 : -1;
      int highz = parentO[2] + parentLength == maxLength ? 0 : 1;

      for (int z = lowz; z <= highz; z++) {
        for (int y = lowy; y <= highy; y++) {
          for (int x = lowx; x <= highx; x++) {
            neighbors.emplace_back(
                Octant<T, D>(parentO + Vector<T, D>(x, y, z) * parentLength, _refinementLevel - 1, _maxRefinement));
          }
        }
      }
    }
    else {
        for (int y = lowy; y <= highy; y++) {
          for (int x = lowx; x <= highx; x++) {
            neighbors.emplace_back(
                Octant<T, D>(parentO + Vector<T, D>(x, y) * parentLength, _refinementLevel - 1, _maxRefinement));
          }
        }
    }
    return neighbors;
  }

  std::vector<std::pair<Octant<T, D>, int>> getCoarseGlobalNeighbors() const
  {
    if (_refinementLevel == 0) {
      return {};
    }
    std::vector<std::pair<Octant<T, D>, int>> neighbors;
    Vector<int, D>                            parentO       = getParentOrigin();
    const int                                 parentLength  = 1 << (_maxRefinement - _refinementLevel + 1);
    const int                                 maxLength     = 1 << _maxRefinement;
    int                                       xNeighbors[3] = {1, 1, 1};
    int                                       yNeighbors[3] = {3, 3, 3};
    int                                       zNeighbors[3] = {9, 9, 9};

    xNeighbors[0]                                           = (parentO[0] == 0) ? 0 : 1;
    yNeighbors[0]                                           = (parentO[1] == 0) ? 0 : 3;
    if constexpr (D == 3) {
      zNeighbors[0]                                           = (parentO[2] == 0) ? 0 : 9;
    }

    xNeighbors[2] = (parentO[0] + parentLength == maxLength) ? 2 : 1;
    yNeighbors[2] = (parentO[1] + parentLength == maxLength) ? 6 : 3;
    if constexpr (D == 3) {
      zNeighbors[2] = (parentO[2] + parentLength == maxLength) ? 18 : 9;
    }
    if constexpr (D == 3) {
      for (int z = -1; z <= 1; z++) {
        for (int y = -1; y <= 1; y++) {
          for (int x = -1; x <= 1; x++) {
            Vector<int, D> newOrigin = (Vector<int, D>(x, y, z) * parentLength) + parentO + maxLength;
            neighbors.emplace_back(
                std::make_pair(Octant<T, D>(newOrigin % maxLength, _refinementLevel - 1, _maxRefinement),
                               xNeighbors[x + 1] + yNeighbors[y + 1] + zNeighbors[z + 1]));
          }
        }
      }
    }
    else {
      for (int y = -1; y <= 1; y++) {
        for (int x = -1; x <= 1; x++) {
          Vector<int, D> newOrigin = (Vector<int, D>(x, y) * parentLength) + parentO + maxLength;
          neighbors.emplace_back(
              std::make_pair(Octant<T, D>(newOrigin % maxLength, _refinementLevel - 1, _maxRefinement),
                             xNeighbors[x + 1] + yNeighbors[y + 1]));
        }
      }
    }

    return neighbors;
  }

  std::vector<std::pair<Octant<T, D>, int>> getAllTouching() const
  {

    std::vector<std::pair<Octant<T, D>, int>> neighbors;

    auto fine = getFineGlobalTouching();
    auto same = getSameGlobalTouching();
    auto coarse = getCoarseGlobalTouching();

    neighbors.insert(
    neighbors.end(),
    std::make_move_iterator(fine.begin()),
    std::make_move_iterator(fine.end()));
    neighbors.insert(
    neighbors.end(),
    std::make_move_iterator(same.begin()),
    std::make_move_iterator(same.end()));
    neighbors.insert(
    neighbors.end(),
    std::make_move_iterator(coarse.begin()),
    std::make_move_iterator(coarse.end()));
    return neighbors;
  }

  std::vector<std::pair<Octant<T, D>, int>> getSameGlobalTouching() const
  {
    std::vector<std::pair<Octant<T, D>, int>> neighbors;
    const int                                 localLength  = 1 << (_maxRefinement - _refinementLevel);
    const int                                 maxLength     = 1 << _maxRefinement;
    int                                       xNeighbors[3] = {1, 1, 1};
    int                                       yNeighbors[3] = {3, 3, 3};
    int                                       zNeighbors[3] = {9, 9, 9};
    xNeighbors[0]                                           = (_origin[0] == 0) ? 0 : 1;
    yNeighbors[0]                                           = (_origin[1] == 0) ? 0 : 3;
    if constexpr (D == 3) {
      zNeighbors[0]                                           = (_origin[2] == 0) ? 0 : 9;
    }

    xNeighbors[2] = (_origin[0] + localLength == maxLength) ? 2 : 1;
    yNeighbors[2] = (_origin[1] + localLength == maxLength) ? 6 : 3;

    if constexpr (D == 3) {
      zNeighbors[2] = (_origin[2] + localLength == maxLength) ? 18 : 9;
      for (int z = -1; z <= 1; z++) {
        for (int y = -1; y <= 1; y++) {
          for (int x = -1; x <= 1; x++) {
            if (z == 0 && y == 0 && x == 0) {
              continue;
            }
            Vector<int,D> newOrigin = Vector<int,D>(x,y,z) * localLength + _origin + maxLength;
            neighbors.emplace_back(
                std::make_pair(Octant<T, D>(newOrigin % maxLength, _refinementLevel, _maxRefinement),
                               xNeighbors[x + 1] + yNeighbors[y + 1] + zNeighbors[z + 1]));
          }
        }
      }
    }
    else {
        for (int y = -1; y <= 1; y++) {
          for (int x = -1; x <= 1; x++) {
            if (y == 0 && x == 0) {
              continue;
            }
            Vector<int,D> newOrigin = Vector<int,D>(x,y) * localLength + _origin + maxLength;
            neighbors.emplace_back(
                std::make_pair(Octant<T, D>(newOrigin % maxLength, _refinementLevel, _maxRefinement),
                               xNeighbors[x + 1] + yNeighbors[y + 1]));
          }
        }

    }

    return neighbors;
  }

  std::vector<std::pair<Octant<T, D>, int>> getFineGlobalTouching() const
  {
    if (_maxRefinement == _refinementLevel) {
      return {};
    }
    std::vector<std::pair<Octant<T, D>, int>> neighbors;
    const int                                 localLength  = 1 << (_maxRefinement - _refinementLevel);
    const int                                 childLength  = localLength >> 1;
    const int                                 maxLength     = 1 << _maxRefinement;
    int                                       xNeighbors[3] = {1, 1, 1};
    int                                       yNeighbors[3] = {3, 3, 3};
    int                                       zNeighbors[3] = {9, 9, 9};
    xNeighbors[0]                                           = (_origin[0] == 0) ? 0 : 1;
    yNeighbors[0]                                           = (_origin[1] == 0) ? 0 : 3;
    if constexpr (D == 3) {
      zNeighbors[0]                                           = (_origin[2] == 0) ? 0 : 9;
    }

    xNeighbors[2] = (_origin[0] + localLength == maxLength) ? 2 : 1;
    yNeighbors[2] = (_origin[1] + localLength == maxLength) ? 6 : 3;
    if constexpr (D == 3) {
      zNeighbors[2] = (_origin[2] + localLength == maxLength) ? 18 : 9;
    }

    for (int c = 0; c < (int)util::pow(2,D); c++) {
      Vector<int,D> currentOffset;
      for (int i = 0; i < (int)D; i++) {
        currentOffset[i] = (c & (int) util::pow(2,i) ? 1 : 0) * childLength;
      }
      int lowX    = neighborLookup[c][0];
      int lowY    = neighborLookup[c][1];
      if constexpr (D == 3) {
        int lowZ    = neighborLookup[c][2];
        for (int z = lowZ; z <= lowZ + 1; z++) {
          for (int y = lowY; y <= lowY + 1; y++) {
            for (int x = lowX; x <= lowX + 1; x++) {
              Vector<int,D> newOrigin = Vector<int,D>(x,y,z)*childLength + _origin + currentOffset + maxLength;
              neighbors.emplace_back(
                  std::make_pair(Octant<T, D>(newOrigin % maxLength, _refinementLevel + 1, _maxRefinement),
                                 xNeighbors[x + 1] + yNeighbors[y + 1] + zNeighbors[z + 1]));
            }
          }
        }
      }
      else {
        for (int y = lowY; y <= lowY + 1; y++) {
          for (int x = lowX; x <= lowX + 1; x++) {
            Vector<int,D> newOrigin = Vector<int,D>(x,y)*childLength + _origin + currentOffset + maxLength;
            neighbors.emplace_back(
                std::make_pair(Octant<T, D>(newOrigin % maxLength, _refinementLevel + 1, _maxRefinement),
                               xNeighbors[x + 1] + yNeighbors[y + 1]));
          }
        }
      }

    }
    return neighbors;
  }


  std::vector<std::pair<Octant<T, D>, int>> getCoarseGlobalTouching() const
  {
    if (_refinementLevel == 0) {
      return {};
    }
    std::vector<std::pair<Octant<T, D>, int>> neighbors;
    Vector<int, D>                            parentO       = getParentOrigin();
    const int                                 parentLength  = 1 << (_maxRefinement - _refinementLevel + 1);
    const int                                 maxLength     = 1 << _maxRefinement;
    int                                       xNeighbors[3] = {1, 1, 1};
    int                                       yNeighbors[3] = {3, 3, 3};
    int                                       zNeighbors[3] = {9, 9, 9};
    xNeighbors[0]                                           = (parentO[0] == 0) ? 0 : 1;
    yNeighbors[0]                                           = (parentO[1] == 0) ? 0 : 3;
    if constexpr (D == 3) {
      zNeighbors[0]                                           = (parentO[2] == 0) ? 0 : 9;
    }

    xNeighbors[2] = (parentO[0] + parentLength == maxLength) ? 2 : 1;
    yNeighbors[2] = (parentO[1] + parentLength == maxLength) ? 6 : 3;
    if constexpr (D == 3) {
      zNeighbors[2] = (parentO[2] + parentLength == maxLength) ? 18 : 9;
    }

    int childID = calcChildID();
    int lowX    = neighborLookup[childID][0];
    int lowY    = neighborLookup[childID][1];

    if constexpr (D == 3) {
      int lowZ    = neighborLookup[childID][2];
      for (int z = lowZ; z <= lowZ + 1; z++) {
        for (int y = lowY; y <= lowY + 1; y++) {
          for (int x = lowX; x <= lowX + 1; x++) {
            Vector<int, D> newOrigin = (Vector<int, D>(x, y, z) * parentLength) + parentO + maxLength;
            neighbors.emplace_back(
                std::make_pair(Octant<T, D>(newOrigin % maxLength, _refinementLevel - 1, _maxRefinement),
                               xNeighbors[x + 1] + yNeighbors[y + 1] + zNeighbors[z + 1]));
          }
        }
      }
    }
    else {
        for (int y = lowY; y <= lowY + 1; y++) {
          for (int x = lowX; x <= lowX + 1; x++) {
            Vector<int, D> newOrigin = (Vector<int, D>(x, y) * parentLength) + parentO + maxLength;
            neighbors.emplace_back(
                std::make_pair(Octant<T, D>(newOrigin % maxLength, _refinementLevel - 1, _maxRefinement),
                               xNeighbors[x + 1] + yNeighbors[y + 1]));
          }
        }

    }

    return neighbors;
  }

  void setID(std::uint64_t id)
  {
    _globalID = std::make_optional((int)id);
  }

  std::optional<int> getProcess() const  { return process; }

  bool operator<(const Octant<T, D>& rhs) const
  {
    return getId() < rhs.getId() || (getId() == rhs.getId() && _refinementLevel < rhs._refinementLevel);
  }

  bool operator==(const Octant<T, D>& rhs) const
  {
    return getId() == rhs.getId() && _refinementLevel == rhs._refinementLevel;
  }

  bool isBorder()
  {
    const int localLength = 1 << (_maxRefinement - _refinementLevel);
    const int maxLength   = 1 << _maxRefinement;
    if constexpr (D == 2) {
      return _origin[0] == 0 || _origin[1] == 0 || _origin[0] + localLength == maxLength ||
       _origin[1] + localLength == maxLength;
    }
    else {
      return _origin[0] == 0 || _origin[1] == 0 || _origin[2] == 0 || _origin[0] + localLength == maxLength ||
       _origin[1] + localLength == maxLength || _origin[2] + localLength == maxLength;
    }
  }

  static bool hilbertLEQ(Octant<T, D>& first, Octant<T, D>& second) {
    return first.calcHilbertID() < second.calcHilbertID() || (first.calcHilbertID() == second.calcHilbertID() && first.getRefinementLevel() < second.getRefinementLevel());
  }

  std::uint64_t calcHilbertID(int o = 0) const
  {
    std::uint64_t mortonID = getId();
    std::uint64_t mask = 0;
    for (int i = 0; i < (int)D; i++) {
      mask |= 0b1ULL << i;
    }

    std::uint64_t currentCoordinate;
    std::uint64_t hilbertID = 0b0;
    int currentOrder;
    int currentOrientation = o;

    for (int i = _maxRefinement-1; i >= 0; i--)
    {
      currentCoordinate = mask & (mortonID >> (i*D));
      if constexpr (D == 2) {
        currentOrder = hilbertOrderLookup2D[currentOrientation][currentCoordinate];
        currentOrientation = hilbertOrientation2D[currentOrientation][currentOrder];
      }
      else {
        currentOrder = hilbertOrderLookup[currentOrientation][currentCoordinate];
        currentOrientation = hilbertOrientation[currentOrientation][currentOrder];
      }
      hilbertID <<= D;
      hilbertID |= currentOrder;
    }

    return hilbertID;
  }

  static std::uint64_t calcHilbertID(std::uint64_t mortonID, int maxRefinementLevel)
  {
    std::uint64_t mask = 0;
    for (int i = 0; i < (int)D; i++) {
      mask |= 0b1ULL << i;
    }
    std::uint64_t currentCoordinate;
    std::uint64_t hilbertID = 0b0;
    int currentOrder;
    int currentOrientation = 0;

    for (int i = maxRefinementLevel-1; i >= 0; i--)
    {
      currentCoordinate = mask & (mortonID >> (i*D));

      if constexpr (D == 2) {
        currentOrder = hilbertOrderLookup2D[currentOrientation][currentCoordinate];
        currentOrientation = hilbertOrientation2D[currentOrientation][currentOrder];
      }
      else {
        currentOrder = hilbertOrderLookup[currentOrientation][currentCoordinate];
        currentOrientation = hilbertOrientation[currentOrientation][currentOrder];
      }
      hilbertID <<= D;
      hilbertID |= currentOrder;
    }

    return hilbertID;
  }

  static std::uint64_t calcHilbertID(Vector<int,D> origin, int maxRefinementLevel)
  {
    std::uint64_t hilbertID = 0;
    int currentOrientation = 0;
    int currentOrder;

    for (int i = maxRefinementLevel-1; i >= 0; --i)
    {
      std::uint64_t currentCoordinate = 0;
      for (int dim = 0; dim < (int)D; ++dim) {
        currentCoordinate |= (origin[dim] >> i & 1ULL) << dim;
      }

      if constexpr (D == 2) {
        currentOrder = hilbertOrderLookup2D[currentOrientation][currentCoordinate];
        currentOrientation = hilbertOrientation2D[currentOrientation][currentOrder];
      } else {
        currentOrder = hilbertOrderLookup[currentOrientation][currentCoordinate];
        currentOrientation = hilbertOrientation[currentOrientation][currentOrder];
      }

      hilbertID <<= D;
      hilbertID |= currentOrder;
    }

    return hilbertID;
  }

  static std::uint64_t calcHilbertOrientation(std::uint64_t mortonID, int maxRefinementLevel)
  {
    std::uint64_t mask = 0;
    for (int i = 0; i < (int)D; i++) {
      mask |= 0b1ULL << i;
    }
    std::uint64_t currentCoordinate;
    int currentOrder;
    int currentOrientation = 0;

      for (int i = maxRefinementLevel-1; i >= 0; i--)
      {
      currentCoordinate = mask & (mortonID >> (i*D));
        if constexpr (D == 2) {
          currentOrder = hilbertOrderLookup2D[currentOrientation][currentCoordinate];
          currentOrientation = hilbertOrientation2D[currentOrientation][currentOrder];
        }
        else {
          currentOrder = hilbertOrderLookup[currentOrientation][currentCoordinate];
          currentOrientation = hilbertOrientation[currentOrientation][currentOrder];
        }
      }
    return currentOrientation;
  }


  static std::uint64_t calcMortonID(std::uint64_t hilbertID, int maxRefinementLevel)
  {
    constexpr std::uint64_t mask = (1ULL << D) - 1ULL;

    std::uint64_t mortonID = 0;
    int currentOrientation = 0;

    for (int level = maxRefinementLevel - 1; level >= 0; --level)
    {
      std::uint64_t currentOrder = (hilbertID >> (level * D)) & mask;

      std::uint64_t currentCoordinate = 0;
      if constexpr (D == 2)
      {
        for (std::uint64_t c = 0; c < (1ULL << D); ++c)
        {
          if (hilbertOrderLookup2D[currentOrientation][c] == currentOrder)
          {
            currentCoordinate = c;
            break;
          }
        }
        currentOrientation = hilbertOrientation2D[currentOrientation][currentOrder];
      }
      else
      {
        for (std::uint64_t c = 0; c < (1ULL << D); ++c)
        {
          if (hilbertOrderLookup[currentOrientation][c] == currentOrder)
          {
            currentCoordinate = c;
            break;
          }
        }
        currentOrientation = hilbertOrientation[currentOrientation][currentOrder];
      }
      mortonID <<= D;
      mortonID |= currentCoordinate;
    }
    return mortonID;
  }


};

} // namespace olb

#endif // OLB_REFINEMENT_OCTREE_OCTANT_H
