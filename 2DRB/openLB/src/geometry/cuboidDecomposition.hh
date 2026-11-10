/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007, 2014 Mathias J. Krause
 *                2025 Adrian Kummerlaender
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

/** \file
 * The description of a vector of 3D cuboid  -- generic implementation.
 */


#ifndef CUBOID_DECOMPOSITION_HH
#define CUBOID_DECOMPOSITION_HH

#include <iostream>
#include <algorithm>
#include <set>
#include <limits>

#include "cuboidDecomposition.h"

namespace olb {

template <typename T, unsigned D>
CuboidDecomposition<T,D>::CuboidDecomposition(Vector<T,D> origin, T deltaR, Vector<int,D> extent, int nC)
  : _motherCuboid(origin, deltaR, extent)
  , _cuboids{_motherCuboid}
  , _periodicityOn{}
{
  split(0, nC);
}

template <typename T, unsigned D>
CuboidDecomposition<T,D>::CuboidDecomposition(const Cuboid<T,D>& motherCuboid, int nC)
  : CuboidDecomposition(motherCuboid.getOrigin(), motherCuboid.getDeltaR(), motherCuboid.getExtent(), nC)
{ }

template <typename T, unsigned D>
CuboidDecomposition<T,D>::CuboidDecomposition(IndicatorF<T,D>& indicatorF, T voxelSize, int nC)
  : _motherCuboid(indicatorF.getMin(),
                  voxelSize,
                  (indicatorF.getMax() - indicatorF.getMin()) / voxelSize + 1.5)
  , _cuboids{_motherCuboid}
  , _periodicityOn{}
{
  _cuboids.reserve(nC+2);
  if (nC > 1) {
    split(0, nC);
  }
  shrink(indicatorF);
}

template <typename T, unsigned D>
CuboidDecomposition<T,D>::CuboidDecomposition(IndicatorF<T,D>& indicatorF, T voxelSize, int nC, std::string minimizeBy)
  : _motherCuboid(indicatorF.getMin(),
                  voxelSize,
                  (indicatorF.getMax() - indicatorF.getMin()) / voxelSize + 1.5)
  , _cuboids{_motherCuboid}
  , _periodicityOn{}
{
         if (minimizeBy == "splitOnly") {
    split(0, nC);
  } else if (minimizeBy == "volume") {
    minimizeByVolume(*this, indicatorF, nC);
  }
  else if (minimizeBy == "weight") {
    minimizeByWeight(*this, indicatorF, nC);
  }
}

template <typename T, unsigned D>
CuboidDecomposition<T,D>::CuboidDecomposition(std::vector<Cuboid<T,D>> cuboids, T voxelSize)
  :_cuboids(std::move(cuboids))
  , _periodicityOn{}
{
  if (_cuboids.empty()) {
    if constexpr (D==2) {
      _motherCuboid = Cuboid<T,D>({0,0}, voxelSize,{0,0});
    }
    else {
      _motherCuboid = Cuboid<T,D>({0,0,0}, voxelSize,{0,0,0});

    }
    return;
  }
  Vector<T,D> max = _cuboids[0].getOrigin();
  Vector<T,D> min = _cuboids[0].getOrigin();
  for (Cuboid<T,D>& cuboid: _cuboids) {
    Vector<T,D> curMin = cuboid.getOrigin();
    Vector<T,D> curMax =  voxelSize * cuboid.getExtent() + cuboid.getOrigin();

    for (int i = 0; i < D; i ++) {
      if (curMax[i] > max[i]) {
        max[i] = curMax[i];
      }
      if (curMin[i] < min[i]) {
        min[i] = curMin[i];
      }
    }
  }
  _motherCuboid = Cuboid<T,D>(min, voxelSize, (max - min) / voxelSize + 1.5);
}

template <typename T, unsigned D>
int CuboidDecomposition<T,D>::size() const {
  OLB_PRECONDITION(_cuboids.size() < std::numeric_limits<int>::max());
  return static_cast<int>(_cuboids.size());
}

template <typename T, unsigned D>
std::optional<int> CuboidDecomposition<T,D>::getC(Vector<T,D> physR, int padding) const {
  for (int iC = 0; iC < size(); ++iC) {
    if (get(iC).isInside(physR, padding)) {
      return iC;
    }
  }
  return std::nullopt;
}

template <typename T, unsigned D>
std::optional<int> CuboidDecomposition<T,D>::getC(Vector<T,D> physR,
                                                  int padding,
                                                  const std::vector<int>& candidates) const {
  for (int iC : candidates) {
    if (get(iC).isInside(physR, padding)) {
      return iC;
    }
  }
  return std::nullopt;
}

template <typename T, unsigned D>
std::optional<int> CuboidDecomposition<T,D>::getC(Vector<T,D> physR,
                                                  int padding,
                                                  const std::vector<int>& candidates,
                                                  int* lastC) const {
  if (lastC && *lastC != -1) {
    if (get(*lastC).isInside(physR, padding)) {
      return *lastC;
    }
  }

  for (int iC : candidates) {
    if (get(iC).isInside(physR, padding)) {
      if (lastC) {
        *lastC = iC;
      }
      return iC;
    }
  }
  return std::nullopt;
}

template <typename T, unsigned D>
std::optional<LatticeR<D+1>> CuboidDecomposition<T,D>::getLatticeR(Vector<T,D> physR) const {
  if (auto iC = getC(physR, 0)) {
    const auto& c = get(*iC);
    return LatticeR<D>{util::floor((physR - c.getOrigin()) / c.getDeltaR() + 0.5)}.withPrefix(*iC);
  }
  return std::nullopt;
}

template <typename T, unsigned D>
std::optional<LatticeR<D+1>> CuboidDecomposition<T,D>::getSurroundingLatticeR(Vector<T,D> physR,
                                                                              int sourceC,
                                                                              int* lastC) const {
  // Use the neighbor cache when it has been built; otherwise fall back to a full
  // search (correct, just slower) so the decomposition is usable even before
  // buildNeighborCache has run.
  auto iC = (sourceC >= 0 && sourceC < static_cast<int>(_neighborCache.size()))
            ? getC(physR, 0, _neighborCache[sourceC], lastC)
            : getC(physR, 0);
  // The neighbor cache excludes the cuboid itself (see getNeighborhood), so a
  // periodic overlap that wraps onto the source cuboid (e.g. single-cuboid
  // periodicity) is not found via the cache. Fall back to the full search.
  if (!iC) {
    iC = getC(physR, 0);
  }
  if (iC) {
    const auto& c = get(*iC);
    return LatticeR<D>{util::floor((physR - c.getOrigin()) / c.getDeltaR() + 0.5)}.withPrefix(*iC);
  }
  return std::nullopt;
}

template <typename T, unsigned D>
std::optional<LatticeR<D+1>> CuboidDecomposition<T,D>::getFloorLatticeR(Vector<T,D> physR) const {
  if (auto iC = getC(physR, 0)) {
    const auto& c = get(*iC);
    return LatticeR<D>{util::floor((physR - c.getOrigin()) / c.getDeltaR())}.withPrefix(*iC);
  }
  return std::nullopt;
}

template <typename T, unsigned D>
T CuboidDecomposition<T,D>::getDeltaR() const {
  return _motherCuboid.getDeltaR();
}

template <typename T, unsigned D>
bool CuboidDecomposition<T,D>::isInside(Vector<T,D> physR) const {
  Vector<T,D> wrappedPhysR = physR;
  const Vector<T,D> mOrigin = _motherCuboid.getOrigin();
  const T deltaR = _motherCuboid.getDeltaR();

  for (unsigned iD = 0; iD < D; ++iD) {
    if (_periodicityOn[iD]) {
      const T length = deltaR * _motherCuboid.getExtent()[iD];
      wrappedPhysR[iD] = std::remainder(physR[iD] - mOrigin[iD] + length, length);
      // Floating point correction for remainder near zero
      if (wrappedPhysR[iD] * wrappedPhysR[iD] < 1e-6 * deltaR * deltaR) {
        wrappedPhysR[iD] = T{};
      }
      if (wrappedPhysR[iD] < 0) {
        wrappedPhysR[iD] += length;
      }
      wrappedPhysR[iD] += mOrigin[iD];
    }
  }

  for (const auto& cuboid : _cuboids) {
    if (cuboid.isInside(wrappedPhysR, 1)) {
      return true;
    }
  }
  return false;
}

template <typename T, unsigned D>
bool CuboidDecomposition<T,D>::isInsideCached(Vector<T,D> physR, int& lastC) const {
  if (lastC >= 0 && lastC < size()) {
    if (get(lastC).isInside(physR, 1)) {
      return true;
    }
    if (!_neighborCache.empty()) {
      for (int iN : _neighborCache[lastC]) {
        if (get(iN).isInside(physR, 1)) {
          lastC = iN;
          return true;
        }
      }
    }
  }

  for (int i = 0; i < size(); ++i) {
    if (get(i).isInside(physR, 1)) {
      lastC = i;
      return true;
    }
  }
  return false;
}

template <typename T, unsigned D>
const Cuboid<T,D>& CuboidDecomposition<T,D>::get(int iC) const
{
  return _cuboids.at(iC);
}

template <typename T, unsigned D>
Cuboid<T,D>& CuboidDecomposition<T,D>::get(int iC)
{
  return _cuboids[iC];
}

template <typename T, unsigned D>
const Cuboid<T,D>& CuboidDecomposition<T,D>::getMotherCuboid() const
{
  return _motherCuboid;
}

template <typename T, unsigned D>
Cuboid<T,D>& CuboidDecomposition<T,D>::getMotherCuboid()
{
  return _motherCuboid;
}

template <typename T, unsigned D>
Vector<T,D> CuboidDecomposition<T,D>::getPhysR(LatticeR<D+1> latticeR) const
{
  auto physR = get(latticeR[0]).getPhysR(latticeR.withoutPrefix());
  for (unsigned iD=0; iD < D; ++iD) {
    if (_periodicityOn[iD]) {
      physR[iD] = remainder(physR[iD] - _motherCuboid.getOrigin()[iD]
                                      + _motherCuboid.getDeltaR() * (_motherCuboid.getExtent()[iD]),
                            _motherCuboid.getDeltaR() * (_motherCuboid.getExtent()[iD]));
      // solving the rounding error problem for double
      if (physR[iD]*physR[iD] < 0.001 * _motherCuboid.getDeltaR()*_motherCuboid.getDeltaR()) {
        physR[iD] = T{};
      }
      // make it to mod instead remainer
      if (physR[iD] < 0) {
        physR[iD] += _motherCuboid.getDeltaR() * (_motherCuboid.getExtent()[iD]);
      }
      physR[iD] += _motherCuboid.getOrigin()[iD];
    }
  }
  return physR;
}

template <typename T, unsigned D>
void CuboidDecomposition<T,D>::setPeriodicity(Vector<bool,D> periodicity)
{
  _periodicityOn = periodicity;
}

template <typename T, unsigned D>
std::set<int> CuboidDecomposition<T,D>::getNeighborhood(int cuboid, int overlap) const
{
  std::set<int> dummy;
  if constexpr (D == 3) {
    for (int iC = 0; iC < size(); iC++) {
      if (cuboid == iC) {
        continue;
      }
      T globX = get(iC).getOrigin()[0];
      T globY = get(iC).getOrigin()[1];
      T globZ = get(iC).getOrigin()[2];
      T nX = get(iC).getNx();
      T nY = get(iC).getNy();
      T nZ = get(iC).getNz();
      T deltaR = get(iC).getDeltaR();
      if (get(cuboid).intersects({globX - overlap * deltaR,
                                  globY - overlap * deltaR,
                                  globZ - overlap * deltaR},
                                 {globX + (nX + overlap - 1)*deltaR,
                                  globY + (nY + overlap - 1)*deltaR,
                                  globZ + (nZ + overlap - 1)*deltaR}, overlap)) {
        dummy.insert(iC);
      }

      if (_periodicityOn[0]) {
        if (get(cuboid).getOrigin()[0] + (get(cuboid).getNx() + overlap - 1)*get(cuboid).getDeltaR() > getMaxPhysR()[0]) {
          Cuboid<T,D> cub({get(cuboid).getOrigin()[0]-getMaxPhysR()[0],
                           get(cuboid).getOrigin()[1],
                           get(cuboid).getOrigin()[2]},
                           get(cuboid).getDeltaR(),
                          {get(cuboid).getNx(),
                           get(cuboid).getNy(),
                           get(cuboid).getNz()});
          if (cub.intersects({globX - overlap * deltaR,
                               globY - overlap * deltaR,
                               globZ - overlap * deltaR},
                              {globY + (nY + overlap - 1)*deltaR,
                               globX + (nX + overlap - 1)*deltaR,
                               globZ + (nZ + overlap - 1)*deltaR}, overlap)) {
            dummy.insert(iC);
          }
        }
        if (get(cuboid).getOrigin()[0] - overlap*get(cuboid).getDeltaR() < getMinPhysR()[0]) {
          Cuboid<T,D> cub({get(cuboid).getOrigin()[0]+getMaxPhysR()[0],
                           get(cuboid).getOrigin()[1],
                           get(cuboid).getOrigin()[2]},
                           get(cuboid).getDeltaR(),
                          {get(cuboid).getNx(),
                           get(cuboid).getNy(),
                           get(cuboid).getNz()});
          if (cub.intersects({globX - overlap * deltaR,
                               globY - overlap * deltaR,
                               globZ - overlap * deltaR},
                              {globX + (nX + overlap - 1)*deltaR,
                               globY + (nY + overlap - 1)*deltaR,
                               globZ + (nZ + overlap - 1)*deltaR}, overlap)) {
            dummy.insert(iC);
          }
        }
      }

      if (_periodicityOn[1]) {
        if (get(cuboid).getOrigin()[1] + (get(cuboid).getNy() + overlap - 1)*get(cuboid).getDeltaR() > getMaxPhysR()[1]) {
          Cuboid<T,D> cub({get(cuboid).getOrigin()[0],
                           get(cuboid).getOrigin()[1]-getMaxPhysR()[1],
                           get(cuboid).getOrigin()[2]},
                           get(cuboid).getDeltaR(),
                          {get(cuboid).getNx(),
                           get(cuboid).getNy(),
                           get(cuboid).getNz()});
          if (cub.intersects({globX - overlap * deltaR,
                               globY - overlap * deltaR,
                               globZ - overlap * deltaR},
                              {globX + (nX + overlap - 1)*deltaR,
                               globY + (nY + overlap - 1)*deltaR,
                               globZ + (nZ + overlap - 1)*deltaR}, overlap)) {
            dummy.insert(iC);
          }
        }
        if (get(cuboid).getOrigin()[1] - overlap*get(cuboid).getDeltaR() < getMinPhysR()[1]) {
          Cuboid<T,D> cub({get(cuboid).getOrigin()[0],
                           get(cuboid).getOrigin()[1]+getMaxPhysR()[1],
                           get(cuboid).getOrigin()[2]},
                           get(cuboid).getDeltaR(),
                          {get(cuboid).getNx(),
                           get(cuboid).getNy(),
                           get(cuboid).getNz()});
          if (cub.intersects({globX - overlap * deltaR,
                               globY - overlap * deltaR,
                               globZ - overlap * deltaR},
                              {globX + (nX + overlap - 1)*deltaR,
                               globY + (nY + overlap - 1)*deltaR,
                               globZ + (nZ + overlap - 1)*deltaR}, overlap)) {
            dummy.insert(iC);
          }
        }
      }

      if (_periodicityOn[2]) {
        if (get(cuboid).getOrigin()[2] + (get(cuboid).getNz() + overlap - 1)*get(cuboid).getDeltaR() > getMaxPhysR()[2]) {
          Cuboid<T,D> cub({get(cuboid).getOrigin()[0],
                           get(cuboid).getOrigin()[1],
                           get(cuboid).getOrigin()[2]-getMaxPhysR()[2]},
                           get(cuboid).getDeltaR(),
                          {get(cuboid).getNx(),
                           get(cuboid).getNy(),
                           get(cuboid).getNz()});
          if (cub.intersects({globX - overlap * deltaR,
                               globY - overlap * deltaR,
                               globZ - overlap * deltaR},
                              {globX + (nX + overlap - 1)*deltaR,
                               globY + (nY + overlap - 1)*deltaR,
                               globZ + (nZ + overlap - 1)*deltaR}, overlap)) {
            dummy.insert(iC);
          }
        }
        if (get(cuboid).getOrigin()[2] - overlap*get(cuboid).getDeltaR() < getMinPhysR()[2]) {
          Cuboid<T,D> cub({get(cuboid).getOrigin()[0],
                           get(cuboid).getOrigin()[1],
                           get(cuboid).getOrigin()[2]+getMaxPhysR()[2]},
                           get(cuboid).getDeltaR(),
                          {get(cuboid).getNx(),
                           get(cuboid).getNy(),
                           get(cuboid).getNz()});
          if (cub.intersects({globX - overlap * deltaR,
                               globY - overlap * deltaR,
                               globZ - overlap * deltaR},
                              {globX + (nX + overlap - 1)*deltaR,
                               globY + (nY + overlap - 1)*deltaR,
                               globZ + (nZ + overlap - 1)*deltaR}, overlap)) {
            dummy.insert(iC);
          }
        }
      }
    }
  } else if constexpr (D == 2) {
    for (int iC = 0; iC < size(); iC++) {
      if (cuboid == iC) {
        continue;
      }
      T globX = get(iC).getOrigin()[0];
      T globY = get(iC).getOrigin()[1];
      T nX = get(iC).getNx();
      T nY = get(iC).getNy();
      T deltaR = get(iC).getDeltaR();
      if (get(cuboid).intersects({globX,
                                   globY - overlap * deltaR},
                                  {globX + (nX + overlap - 1)*deltaR,
                                   globY + (nY + overlap - 1)*deltaR}, overlap)) {
        dummy.insert(iC);
      }
      if (_periodicityOn[0]) {
        if (get(cuboid).getOrigin()[0] + (get(cuboid).getNx() + overlap - 1)*get(cuboid).getDeltaR() > getMaxPhysR()[0]) {
          Cuboid2D<T> cub({get(cuboid).getOrigin()[0]-getMaxPhysR()[0],
                           get(cuboid).getOrigin()[1]},
                           get(cuboid).getDeltaR(),
                          {get(cuboid).getNx(),
                           get(cuboid).getNy()});
          if (cub.intersects({globX - overlap * deltaR,
                               globY - overlap * deltaR},
                              {globX + (nX + overlap - 1)*deltaR,
                               globY + (nY + overlap - 1)*deltaR}, overlap)) {
            dummy.insert(iC);
          }
        }
        if (get(cuboid).getOrigin()[0] - overlap*get(cuboid).getDeltaR() < getMinPhysR()[0]) {
          Cuboid2D<T> cub({get(cuboid).getOrigin()[0]+getMaxPhysR()[0],
                           get(cuboid).getOrigin()[1]},
                           get(cuboid).getDeltaR(),
                          {get(cuboid).getNx(),
                           get(cuboid).getNy()});
          if (cub.intersects({globX - overlap * deltaR,
                               globY - overlap * deltaR},
                              {globX + (nX + overlap - 1)*deltaR,
                               globY + (nY + overlap - 1)*deltaR}, overlap)) {
            dummy.insert(iC);
          }
        }
      }

      if (_periodicityOn[1]) {
        if (get(cuboid).getOrigin()[1] + (get(cuboid).getNy() + overlap - 1)*get(cuboid).getDeltaR() > getMaxPhysR()[1]) {
          Cuboid2D<T> cub({get(cuboid).getOrigin()[0],
                           get(cuboid).getOrigin()[1]-getMaxPhysR()[1]},
                           get(cuboid).getDeltaR(),
                          {get(cuboid).getNx(),
                           get(cuboid).getNy()});
          if (cub.intersects({globX - overlap * deltaR,
                               globY - overlap * deltaR},
                              {globX + (nX + overlap - 1)*deltaR,
                               globY + (nY + overlap - 1)*deltaR}, overlap)) {
            dummy.insert(iC);
          }
        }
        if (get(cuboid).getOrigin()[1] - overlap*get(cuboid).getDeltaR() < getMinPhysR()[1]) {
          Cuboid2D<T> cub({get(cuboid).getOrigin()[0],
                           get(cuboid).getOrigin()[1]+getMaxPhysR()[1]},
                           get(cuboid).getDeltaR(),
                          {get(cuboid).getNx(),
                           get(cuboid).getNy()});
          if (cub.intersects({globX - overlap * deltaR,
                               globY - overlap * deltaR},
                              {globX + (nX + overlap - 1)*deltaR,
                               globY + (nY + overlap - 1)*deltaR}, overlap)) {
            dummy.insert(iC);
          }
        }
      }
    }
  }
  return dummy;
}

template <typename T, unsigned D>
T CuboidDecomposition<T,D>::getMinRatio() const
{
  T minRatio = 1.;
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if ((T)_cuboids[i].getNx() / (T)_cuboids[i].getNy() < minRatio) {
      minRatio = (T)_cuboids[i].getNx() / (T)_cuboids[i].getNy();
    }
    if constexpr (D == 3) {
      if ((T)_cuboids[i].getNy() / (T)_cuboids[i].getNz() < minRatio) {
        minRatio = (T)_cuboids[i].getNy() / (T)_cuboids[i].getNz();
      }
      if ((T)_cuboids[i].getNz() / (T)_cuboids[i].getNx() < minRatio) {
        minRatio = (T)_cuboids[i].getNz() / (T)_cuboids[i].getNx();
      }
    }
  }
  return minRatio;
}

template <typename T, unsigned D>
T CuboidDecomposition<T,D>::getMaxRatio() const
{
  T maxRatio = 1.;
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if ((T)_cuboids[i].getNx() / (T)_cuboids[i].getNy() > maxRatio) {
      maxRatio = (T)_cuboids[i].getNx() / (T)_cuboids[i].getNy();
    }
    if constexpr (D == 3) {
      if ((T)_cuboids[i].getNy() / (T)_cuboids[i].getNz() > maxRatio) {
        maxRatio = (T)_cuboids[i].getNy() / (T)_cuboids[i].getNz();
      }
      if ((T)_cuboids[i].getNz() / (T)_cuboids[i].getNx() > maxRatio) {
        maxRatio = (T)_cuboids[i].getNz() / (T)_cuboids[i].getNx();
      }
    }
  }
  return maxRatio;
}

template <typename T, unsigned D>
Vector<T,D> CuboidDecomposition<T,D>::getMinPhysR() const
{
  Vector<T,D> output(_cuboids[0].getOrigin());
  for (int i = 0; i < size(); i++) {
    for (unsigned iD=0; iD < D; ++iD) {
      if (_cuboids[i].getOrigin()[iD] < output[iD]) {
        output[iD] = _cuboids[i].getOrigin()[iD];
      }
    }
  }
  return output;
}

template <typename T, unsigned D>
Vector<T,D> CuboidDecomposition<T,D>::getMaxPhysR() const
{
  Vector<T,D> output(_cuboids[0].getOrigin());
  for (unsigned iD=0; iD < D; ++iD) {
    output[iD] += _cuboids[0].getExtent()[iD]*_cuboids[0].getDeltaR();
  }
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    for (unsigned iD=0; iD < D; ++iD) {
      if (_cuboids[i].getOrigin()[iD] + _cuboids[i].getExtent()[iD]*_cuboids[i].getDeltaR() > output[iD]) {
        output[iD] = _cuboids[i].getOrigin()[iD] + _cuboids[i].getExtent()[iD]*_cuboids[i].getDeltaR();
      }
    }
  }
  return output;
}

template <typename T, unsigned D>
T CuboidDecomposition<T,D>::getMinPhysVolume() const
{
  T minVolume = _cuboids[0].getPhysVolume();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getPhysVolume() < minVolume) {
      minVolume = _cuboids[i].getPhysVolume();
    }
  }
  return minVolume;
}

template <typename T, unsigned D>
T CuboidDecomposition<T,D>::getMaxPhysVolume() const
{
  T maxVolume = _cuboids[0].getPhysVolume();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getPhysVolume() > maxVolume) {
      maxVolume = _cuboids[i].getPhysVolume();
    }
  }
  return maxVolume;
}

template <typename T, unsigned D>
std::size_t CuboidDecomposition<T,D>::getMinLatticeVolume() const
{
  std::size_t minNodes = _cuboids[0].getLatticeVolume();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getLatticeVolume() < minNodes) {
      minNodes = _cuboids[i].getLatticeVolume();
    }
  }
  return minNodes;
}

template <typename T, unsigned D>
std::size_t CuboidDecomposition<T,D>::getMaxLatticeVolume() const
{
  std::size_t maxNodes = _cuboids[0].getLatticeVolume();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getLatticeVolume() > maxNodes) {
      maxNodes = _cuboids[i].getLatticeVolume();
    }
  }
  return maxNodes;
}

template <typename T, unsigned D>
std::size_t CuboidDecomposition<T,D>::getLatticeVolume() const
{
  std::size_t sumNodes = 0;
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    sumNodes += _cuboids[i].getLatticeVolume();
  }
  return sumNodes;
}

template <typename T, unsigned D>
std::size_t CuboidDecomposition<T,D>::getNumNodes() const
{
  std::size_t numNodes = 0;
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    numNodes = _cuboids[i].getLatticeVolume();
  }
  return numNodes;
}

template <typename T, unsigned D>
std::size_t CuboidDecomposition<T,D>::getMinLatticeWeight() const
{
  std::size_t minNodes = _cuboids[0].getWeight();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getWeight() < minNodes) {
      minNodes = _cuboids[i].getWeight();
    }
  }
  return minNodes;
}

template <typename T, unsigned D>
std::size_t CuboidDecomposition<T,D>::getMaxLatticeWeight() const
{
  std::size_t maxNodes = _cuboids[0].getWeight();
  for (unsigned i = 0; i < _cuboids.size(); i++) {
    if (_cuboids[i].getWeight() > maxNodes) {
      maxNodes = _cuboids[i].getWeight();
    }
  }
  return maxNodes;
}

template <typename T, unsigned D>
void CuboidDecomposition<T,D>::remove(int iC) {
  _cuboids.erase(_cuboids.begin() + iC);
}


template <typename T, unsigned D>
void CuboidDecomposition<T,D>::remove(IndicatorF<T,D>& indicatorF)
{
  std::vector<bool> allZero{};
  LatticeR<D+1> latticeR{};
  for (int iC = 0; iC < size(); iC++) {
    latticeR[0] = iC;
    allZero.push_back(true);
    for (int iX = 0; iX < _cuboids[iC].getNx(); iX++) {
      latticeR[1] = iX;
      for (int iY = 0; iY < _cuboids[iC].getNy(); iY++) {
        latticeR[2] = iY;
        if constexpr (D == 3) {
          for (int iZ = 0; iZ < _cuboids[iC].getNz(); iZ++) {
            latticeR[3] = iZ;
            auto physR = getPhysR(latticeR);
            bool inside[1];
            indicatorF(inside, physR.data());
            if (inside[0]) {
              allZero[iC] = 0;
            }
          }
        } else {
          auto physR = getPhysR(latticeR);
          bool inside[1];
          indicatorF(inside, physR.data());
          if (inside[0]) {
            allZero[iC] = 0;
          }
        }
      }
    }
  }
  for (int iC = _cuboids.size() - 1; iC >= 0; iC--) {
    if (allZero[iC]) {
      remove(iC);
    }
  }
}

template <typename T, unsigned D>
void CuboidDecomposition<T,D>::removeContaining(IndicatorF<T,D>& indicatorF)
{
  std::vector<bool> allZero{};
  LatticeR<D+1> latticeR{};
  for (int iC = 0; iC < size(); iC++) {
    latticeR[0] = iC;
    allZero.push_back(true);
    for (int iX = 0; iX < _cuboids[iC].getNx(); iX++) {
      latticeR[1] = iX;
      for (int iY = 0; iY < _cuboids[iC].getNy(); iY++) {
        latticeR[2] = iY;
        if constexpr (D == 3) {
          for (int iZ = 0; iZ < _cuboids[iC].getNz(); iZ++) {
            latticeR[3] = iZ;
            auto physR = getPhysR(latticeR);
            bool inside[1];
            indicatorF(inside, physR.data());
            if (inside[0]) {
              allZero[iC] = 0;
            }
          }
        } else {
          auto physR = getPhysR(latticeR);
          bool inside[1];
          indicatorF(inside, physR.data());
          if (inside[0]) {
            allZero[iC] = 0;
          }
        }
      }
    }
  }
  for (int iC = _cuboids.size() - 1; iC >= 0; iC--) {
    if (!allZero[iC]) {
      remove(iC);
    }
  }
}

template <typename T, unsigned D>
void CuboidDecomposition<T,D>::removeByWeight()
{
  std::vector<bool> allZero(_cuboids.size(), false);
  for (unsigned iC = 0; iC < _cuboids.size(); iC++) {
    allZero[iC] = (_cuboids[iC].getWeight() == 0);
  }
  for (int iC = _cuboids.size() - 1; iC >= 0; iC--) {
    if (allZero[iC]) {
      remove(iC);
    }
  }
}

template <typename T, unsigned D>
void CuboidDecomposition<T,D>::shrink(int iC, IndicatorF<T,D>& indicatorF)
{
  if constexpr (D == 3) {
    LatticeR<D+1> latticeR{};
    bool inside[1] {};

    latticeR[0] = iC;
    std::size_t fullCells = 0;
    int xN = get(iC).getNx();
    int yN = get(iC).getNy();
    int zN = get(iC).getNz();
    int maxX = 0;
    int maxY = 0;
    int maxZ = 0;
    int newX = xN - 1;
    int newY = yN - 1;
    int newZ = zN - 1;
    for (int iX = 0; iX < xN; iX++) {
      for (int iY = 0; iY < yN; iY++) {
        for (int iZ = 0; iZ < zN; iZ++) {
          latticeR[1] = iX;
          latticeR[2] = iY;
          latticeR[3] = iZ;
          auto physR = getPhysR(latticeR);
          indicatorF(inside,physR.data());
          if (inside[0]) {
            fullCells++;
            maxX = util::max(maxX, iX);
            maxY = util::max(maxY, iY);
            maxZ = util::max(maxZ, iZ);
            newX = util::min(newX, iX);
            newY = util::min(newY, iY);
            newZ = util::min(newZ, iZ);
          }
        }
      }
    }
    if (fullCells > 0) {
      get(iC).setWeight(fullCells);
      _cuboids[iC].resize({newX, newY, newZ}, {maxX - newX + 1, maxY - newY + 1, maxZ - newZ + 1});
    }
    else {
      remove(iC);
    }
  } else {
    LatticeR<D+1> latticeR{};
    bool inside[1] {};

    latticeR[0] = iC;
    std::size_t fullCells = 0;
    int xN = get(iC).getNx();
    int yN = get(iC).getNy();
    int maxX = 0;
    int maxY = 0;
    int newX = xN - 1;
    int newY = yN - 1;
    for (int iX = 0; iX < xN; iX++) {
      for (int iY = 0; iY < yN; iY++) {
        latticeR[1] = iX;
        latticeR[2] = iY;
        auto physR = getPhysR(latticeR);
        indicatorF(inside,physR.data());
        if (inside[0]) {
          fullCells++;
          maxX = util::max(maxX, iX);
          maxY = util::max(maxY, iY);
          newX = util::min(newX, iX);
          newY = util::min(newY, iY);
        }
      }
    }
    if (fullCells > 0) {
      get(iC).setWeight(fullCells);
      _cuboids[iC].resize({newX, newY}, {maxX - newX + 1, maxY - newY + 1});
    }
    else {
      remove(iC);
    }
  }
}

template <typename T, unsigned D>
void CuboidDecomposition<T,D>::shrink(IndicatorF<T,D>& indicatorF)
{
  for (int iC = size() - 1; iC >= 0; iC--) {
    shrink(iC, indicatorF);
  }
  updateMotherCuboid();
}

template <typename T, unsigned D>
void CuboidDecomposition<T,D>::updateMotherCuboid()
{
  if (_cuboids.empty()) {
    return;
  }
  // shrink mother cuboid
  Vector<T,D> minPhysR = getMinPhysR();
  Vector<T,D> maxPhysR = getMaxPhysR();
  T deltaR = getDeltaR();
  _motherCuboid = Cuboid<T,D>(minPhysR, deltaR, (maxPhysR - minPhysR) / deltaR + 0.5);
}

template <typename T, unsigned D>
void CuboidDecomposition<T,D>::refine(int factor)
{
  _motherCuboid.refine(factor);
  for (auto& cuboid : _cuboids) {
    cuboid.refine(factor);
  }
}

template <typename T, unsigned D>
bool CuboidDecomposition<T,D>::tryRefineTo(T goalDeltaR)
{
  const T tolerance = std::numeric_limits<T>::epsilon();
  const T currDeltaR = _motherCuboid.getDeltaR();
  const int factor = std::ceil(currDeltaR / goalDeltaR);
  if (util::fabs(currDeltaR / factor - goalDeltaR) < tolerance) {
    refine(factor);
    return true;
  } else {
    return false;
  }
}

template <typename T, unsigned D>
void CuboidDecomposition<T,D>::refineCellCentered(int factor)
{
  _motherCuboid.refineCellCentered(factor);
  for (auto& cuboid : _cuboids) {
    cuboid.refineCellCentered(factor);
  }
}

template <typename T, unsigned D>
void CuboidDecomposition<T,D>::refineCellCentered(int iC, int factor)
{
  _cuboids[iC].refineCellCentered(factor);
}

template <typename T, unsigned D>
void CuboidDecomposition<T,D>::split(int iC, int p)
{
  Cuboid<T,D> temp(_cuboids[iC]);
  temp.divideP(p, _cuboids);
  remove(iC);
}

template <typename T, unsigned D>
void CuboidDecomposition<T,D>::splitRegular(int iC, int width)
{
  Cuboid<T,D> temp(_cuboids[iC]);
  const int p = std::max(1, temp.getNx() / width);
  const int q = std::max(1, temp.getNy() / width);
  const int r = std::max(1, temp.getNz() / width);
  temp.divide({p, q, r}, _cuboids);
  remove(iC);
}

template <typename T, unsigned D>
void CuboidDecomposition<T,D>::splitByWeight(int iC, int p, IndicatorF<T,D>& indicatorF) requires (D == 3)
{
  T averageWeight = get(iC).getWeight() / (T) p;
  Cuboid<T,D> temp(_cuboids[iC]);

  int latticeR[4];
  latticeR[0] = iC;
  int xN = get(iC).getNx();
  int yN = get(iC).getNy();
  int zN = get(iC).getNz();
  T deltaR = get(iC).getDeltaR();
  int fullCells = 0;

  Vector<T, 3> globPos_child = get(iC).getOrigin();
  std::vector<int> extend_child = {xN, yN, zN};
  int localPos_child = 0;

  // looking for largest extend, because halfing the cuboid by its largest extend will result in the smallest surface and therfore in the least comminication cells
  if ( get(iC).getNx() >= get(iC).getNy() && get(iC).getNx() >= get(iC).getNz()) {
    // clout << "Cut in x direction!" << std::endl;

    // for each child cuboid search for the optimal cutting plane
    for ( int iChild = 0; iChild < p - 1; iChild++) {
      fullCells = 0;
      int fullCells_minusOne = 0;

      for (int iX = localPos_child; iX < xN; iX++) {
        fullCells_minusOne = fullCells;
        for (int iY = 0; iY < yN; iY++) {
          for (int iZ = 0; iZ < zN; iZ++) {
            latticeR[1] = iX;
            latticeR[2] = iY;
            latticeR[3] = iZ;
            auto physR = getPhysR(latticeR);
            bool inside[1];
            indicatorF(inside,physR.data());
            if (inside[0]) {
              fullCells++;
            }
          }
        }
        // the optimal cutting plane is determined, so that the child cuboid's cells inside the indicator are the closest to the total cells inside the indicator per number of children
        if ( fullCells >= averageWeight ) {
          if ( (fullCells - averageWeight) > (averageWeight - fullCells_minusOne) ) {
            iX--;
          }
          // clout << "found optimal iX = " << iX << std::endl;
          extend_child[0] = iX - localPos_child + 1;

          Cuboid<T,D> child(globPos_child, deltaR, extend_child);
          _cuboids.push_back(child);

          globPos_child[0] += extend_child[0]*deltaR;
          localPos_child += extend_child[0] + 1;
          // clout << "added child " << iChild << " of " << p << std::endl;

          break;
        }
      }

    }

    extend_child[0] = xN - localPos_child + p - 1;

    Cuboid<T,D> child(globPos_child, deltaR, extend_child);
    _cuboids.push_back(child);

    // clout << "added last child of " << p << std::endl;

  }
  else if ( get(iC).getNy() >= get(iC).getNx() && get(iC).getNy() >= get(iC).getNz()) {
    // clout << "Cut in y direction!" << std::endl;

    for ( int iChild = 0; iChild < p - 1; iChild++) {
      fullCells = 0;
      int fullCells_minusOne = 0;

      for (int iY = localPos_child; iY < yN; iY++) {
        fullCells_minusOne = fullCells;
        for (int iX = 0; iX < xN; iX++) {
          for (int iZ = 0; iZ < zN; iZ++) {
            latticeR[1] = iX;
            latticeR[2] = iY;
            latticeR[3] = iZ;
            auto physR = getPhysR(latticeR);
            bool inside[1];
            indicatorF(inside,physR.data());
            if (inside[0]) {
              fullCells++;
            }
          }
        }
        if ( fullCells >= averageWeight ) {
          if ( (fullCells - averageWeight) > (averageWeight - fullCells_minusOne) ) {
            iY--;
          }
          // clout << "found optimal iY = " << iY << std::endl;
          extend_child[1] = iY - localPos_child + 1;

          Cuboid<T,D> child(globPos_child, deltaR, extend_child);
          _cuboids.push_back(child);

          globPos_child[1] += extend_child[1]*deltaR;
          localPos_child += extend_child[1] + 1;
          // clout << "added child " << iChild << " of " << p << std::endl;

          break;
        }
      }

    }

    extend_child[1] = yN - localPos_child + p - 1;

    Cuboid<T,D> child(globPos_child, deltaR, extend_child);
    _cuboids.push_back(child);

    // clout << "added last child of " << p << std::endl;
  }
  else {
    // clout << "Cut in z direction!" << std::endl;

    for ( int iChild = 0; iChild < p - 1; iChild++) {
      fullCells = 0;
      int fullCells_minusOne = 0;

      for (int iZ = localPos_child; iZ < zN; iZ++) {
        fullCells_minusOne = fullCells;
        for (int iY = 0; iY < yN; iY++) {
          for (int iX = 0; iX < xN; iX++) {
            latticeR[1] = iX;
            latticeR[2] = iY;
            latticeR[3] = iZ;
            auto physR = getPhysR(latticeR);
            bool inside[1];
            indicatorF(inside,physR.data());
            if (inside[0]) {
              fullCells++;
            }
          }
        }
        if ( fullCells >= averageWeight ) {
          if ( (fullCells - averageWeight) > (averageWeight - fullCells_minusOne) ) {
            iZ--;
          }
          // clout << "found optimal iZ = " << iZ << std::endl;
          extend_child[2] = iZ - localPos_child + 1;

          Cuboid<T,D> child(globPos_child, deltaR, extend_child);
          _cuboids.push_back(child);

          globPos_child[2] += extend_child[2]*deltaR;
          localPos_child += extend_child[2] + 1;
          // clout << "added child " << iChild << " of " << p << std::endl;

          break;
        }
      }

    }

    extend_child[2] = zN - localPos_child + p - 1;

    Cuboid<T,D> child(globPos_child, deltaR, extend_child);
    _cuboids.push_back(child);

    // clout << "added last child of " << p << std::endl;
  }
}

template <typename T, unsigned D>
void CuboidDecomposition<T,D>::splitFractional(int iC, int iD, std::vector<T> fractions)
{
  Cuboid<T,D> tmp = _cuboids[iC];
  tmp.divideFractional(iD, fractions, _cuboids);
  remove(iC);
}

template <typename T, unsigned D>
std::vector<int> CuboidDecomposition<T,D>::createNestedBox(
  int baseC,
  const Vector<T,D>& relOrigin,
  const Vector<T,D>& relSize)
{
  std::vector<int> neighborIndices;
  int currentInner = baseC;

  for (int iD = 0; iD < D; ++iD) {
    std::vector<T> fractions;
    std::size_t innerFractionIndex = 0;

    if (relOrigin[iD] > std::numeric_limits<T>::epsilon()) {
      fractions.emplace_back(relOrigin[iD]);
      innerFractionIndex = 1;
    }
    fractions.emplace_back(relSize[iD]);
    T remainder = 1.0 - relOrigin[iD] - relSize[iD];
    if (remainder > std::numeric_limits<T>::epsilon()) {
      fractions.emplace_back(remainder);
    }

    if (fractions.size() <= 1) {
      continue;
    }

    std::size_t nBefore = _cuboids.size();
    splitFractional(currentInner, iD, fractions);

    for (int& idx : neighborIndices) {
      if (idx > currentInner) {
        idx--;
      }
    }

    int newStartIndex = static_cast<int>(nBefore) - 1;
    int nextInner = -1;
    for (std::size_t j = 0; j < fractions.size(); ++j) {
      int newIdx = newStartIndex + static_cast<int>(j);
      if (j == innerFractionIndex) {
        nextInner = newIdx;
      } else {
        neighborIndices.emplace_back(newIdx);
      }
    }
    currentInner = nextInner;
  }

  std::vector<int> allCreatedIndices;
  allCreatedIndices.emplace_back(currentInner);
  allCreatedIndices.insert(allCreatedIndices.end(),
                           neighborIndices.begin(),
                           neighborIndices.end());
  return allCreatedIndices;
}

template <typename T, unsigned D>
void CuboidDecomposition<T,D>::setWeights(IndicatorF<T,D>& indicatorF) requires (D == 3)
{
  #ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for schedule(dynamic,1)
  #endif
  for (int iC=0; iC < size(); ++iC) {
    int latticeR[4] { iC, 0, 0, 0 };
    int xN = get(iC).getNx();
    int yN = get(iC).getNy();
    int zN = get(iC).getNz();
    size_t fullCells = 0;
    for (int iX = 0; iX < xN; iX++) {
      for (int iY = 0; iY < yN; iY++) {
        for (int iZ = 0; iZ < zN; iZ++) {
          latticeR[1] = iX;
          latticeR[2] = iY;
          latticeR[3] = iZ;
          auto physR = getPhysR(latticeR);
          bool inside[1];
          indicatorF(inside,physR.data());
          if (inside[0]) {
            fullCells++;
          }
        }
      }
    }
    get(iC).setWeight(fullCells);
  }
}

template <typename T, unsigned D>
void CuboidDecomposition<T,D>::print() const
{
  OstreamManager clout(std::cout, "CuboidDecomposition");
  auto range = [](auto lo, auto hi) {
    return io::text::toCell(lo) + " / " + io::text::toCell(hi);
  };
  clout << io::text::section("Cuboid Structure") << std::endl;
  io::text::KeyValueList()
    .item("cuboids", size())
    .item("delta", getDeltaR())
    .item("ratio (min/max)", range(getMinRatio(), getMaxRatio()))
    .item("nodes (min/max)", range(getMinLatticeVolume(), getMaxLatticeVolume()))
    .item("weight (min/max)", range(getMinLatticeWeight(), getMaxLatticeWeight()))
    .render(clout);
}

template <typename T, unsigned D>
void CuboidDecomposition<T,D>::printExtended()
{
  OstreamManager clout(std::cout, "CuboidDecomposition");
  clout << "Mothercuboid :" << std::endl;
  getMotherCuboid().print();

  for (int iC = 0; iC < size(); iC++) {
    clout << "Cuboid #" << iC << ": " << std::endl;
    get(iC).print();
  }
}

template <typename T, unsigned D>
void CuboidDecomposition<T,D>::writeToExistingFile(std::string completeFileName, LoadBalancer<T>& loadBalancer)
{
  OstreamManager clout(std::cout, "CuboidDecomposition");
  std::ofstream fout;
  fout.flags(std::ios::scientific);
  fout.precision (std::numeric_limits<double>::digits10 + 1);
  if ( singleton::mpi().isMainProcessor() ) {

    // Open File
    fout.open(completeFileName.c_str(), std::ios::app);
    if (!fout) {
      clout << "Error: could not open " << completeFileName << std::endl;
    }

    // --- Preamble --- //
    fout << "<CuboidDecomposition dimension=\"3\" ";
    getMotherCuboid().writeAsXML(fout);
    fout << ">\n";

    // TODO: Move Cuboid XML Serialization to Cuboid3D class
    for (int iC = 0; iC < size(); ++iC) {
      fout << "<Cuboid ";
      get(iC).writeAsXML(fout);
      fout << " />\n";
    }

    fout << "</CuboidDecomposition>\n";

    fout.close();
  }
}


template <typename T, unsigned D>
void CuboidDecomposition<T,D>::writeToFile(std::string fileName, LoadBalancer<T>& loadBalancer)
{
  std::string fname = singleton::directories().getLogOutDir() + fileName + ".xml";
  std::ofstream fout;
  if (singleton::mpi().isMainProcessor()) {
    fout.open(fname.c_str(), std::ios::trunc);
    fout << "<?xml version=\"1.0\"?>\n";
    fout << "<XMLContent>\n";
    fout.close();
    fout.clear();
  }

  writeToExistingFile(fname, loadBalancer);

  if (singleton::mpi().isMainProcessor()) {
    fout.open(fname.c_str(), std::ios::app);
    fout << "</XMLContent>\n";
    fout.close();
  }
}

template<typename T, unsigned D>
std::unique_ptr<CuboidDecomposition<T,D>> createCuboidDecomposition(const XMLreader& reader)
{
  std::vector<T> origin = getDataFromTag<T>(reader, "origin", D);
  std::vector<int> extent = getDataFromTag<int>(reader, "extent", D);
  T deltaR = getDataFromTag<T>(reader, "deltaR", 1)[0];
  std::size_t weight = getDataFromTag<size_t>(reader, "weight", 1)[0];

  auto cGeo = std::make_unique<CuboidDecomposition<T,D>>(Vector<T,D>(origin.data()), deltaR, Vector<int,D>(extent.data()));
  cGeo->getMotherCuboid().setWeight(weight);
  cGeo->cuboids().clear();

  if (reader.getAttribute("periodicity") != "Attribute not found.") {
    std::stringstream ss(reader.getAttribute("periodicity"));
    Vector<bool,D> periodicity{};
    for (unsigned iD = 0; iD < D; ++iD) { int b = 0; ss >> b; periodicity[iD] = (b != 0); }
    cGeo->setPeriodicity(periodicity);
  }

  for ( const XMLreader* cub: reader ) {
    if (std::string(cub->getName()) != "Cuboid") continue;
    origin = getDataFromTag<T>(*cub, "origin", D);
    extent = getDataFromTag<int>(*cub, "extent", D);
    deltaR = getDataFromTag<T>(*cub, "deltaR", 1)[0];
    weight = getDataFromTag<int>(*cub, "weight", 1)[0];

    Cuboid<T,D> cuboid(Vector<T,D>(origin.data()), deltaR, Vector<int,D>(extent.data()));
    cuboid.setWeight(weight);

    if (cub->getAttribute("process") != "Attribute not found.") {
      cuboid.setPreCalculatedProcess(std::stoi(cub->getAttribute("process")));
    }
    if (cub->getAttribute("isCoarseBorder") != "Attribute not found.") {
      cuboid.setCoarseBorder(std::stoi(cub->getAttribute("isCoarseBorder")));
    }

    cGeo->cuboids().emplace_back(cuboid);
  }

  return cGeo;
}

template<typename T, unsigned D>
std::unique_ptr<CuboidDecomposition<T,D>> createCuboidDecomposition(std::string fileName)
{
  XMLreader reader(fileName);
  return createCuboidDecomposition<T,D>(reader["CuboidDecomposition"]);
}

}

#endif
