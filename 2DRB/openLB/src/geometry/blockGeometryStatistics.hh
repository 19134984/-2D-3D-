/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2011, 2014 Mathias J. Krause, Simon Zimny, 2026 Adrian Kummerlaender
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
 * Representation of a statistic for a block geometry -- generic implementation.
 */

#ifndef BLOCK_GEOMETRY_STATISTICS_HH
#define BLOCK_GEOMETRY_STATISTICS_HH

#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include "io/textTable.h"
#include "utilities/omath.h"

#include "geometry/blockGeometry.h"
#include "geometry/blockGeometryStatistics.h"

namespace olb {

template<typename T, unsigned D>
BlockGeometryStatistics<T,D>::BlockGeometryStatistics(BlockGeometry<T,D>* blockGeometry)
  : _blockGeometry(blockGeometry),
    clout(std::cout,"BlockGeometryStatistics")
{
  _statisticsUpdateNeeded = true;
}

template<typename T, unsigned D>
bool& BlockGeometryStatistics<T,D>::getStatisticsStatus()
{
  return _statisticsUpdateNeeded;
}

template<typename T, unsigned D>
bool const & BlockGeometryStatistics<T,D>::getStatisticsStatus() const
{
  return _statisticsUpdateNeeded;
}


template<typename T, unsigned D>
void BlockGeometryStatistics<T,D>::update(bool verbose)
{
  const_this = const_cast<const BlockGeometryStatistics<T,D>*>(this);

  if (getStatisticsStatus() ) {
    _material2n.clear();
    _blockGeometry->BlockStructureD<D>::forCoreSpatialLocations([&](LatticeR<D> latticeR) {
      takeStatistics(latticeR);
    });

    _nMaterials = _material2n.size();

    if (verbose) {
      clout << "updated" << std::endl;
    }
    getStatisticsStatus() = false;
  }
}


template<typename T, unsigned D>
int BlockGeometryStatistics<T,D>::getNmaterials()
{
  update();
  return const_this->getNmaterials();
}

template<typename T, unsigned D>
int BlockGeometryStatistics<T,D>::getNmaterials() const
{
  return _nMaterials;
}

template<typename T, unsigned D>
std::size_t BlockGeometryStatistics<T,D>::getNvoxel(int material)
{
  update();
  return const_this->getNvoxel(material);
}

template<typename T, unsigned D>
std::size_t BlockGeometryStatistics<T,D>::getNvoxel(int material) const
{
  try {
    return _material2n.at(material);
  }
  catch (std::out_of_range& ex) {
    return 0;
  }
}

template<typename T, unsigned D>
std::map<int, std::size_t> BlockGeometryStatistics<T,D>::getMaterial2n()
{
  update();
  return const_this->getMaterial2n();
}

template<typename T, unsigned D>
std::map<int, std::size_t> BlockGeometryStatistics<T,D>::getMaterial2n() const
{
  return _material2n;
}

template<typename T, unsigned D>
std::size_t BlockGeometryStatistics<T,D>::getNvoxel()
{
  update();
  return const_this->getNvoxel();
}

template<typename T, unsigned D>
std::size_t BlockGeometryStatistics<T,D>::getNvoxel() const
{
  std::size_t total = 0;
  for (const auto& material : _material2n) {
    total += material.second;
  }
  return total;
}

template<typename T, unsigned D>
std::vector<int> BlockGeometryStatistics<T,D>::getMinLatticeR(int material)
{
  update();
  return const_this->getMinLatticeR(material);
}

template<typename T, unsigned D>
std::vector<int> BlockGeometryStatistics<T,D>::getMinLatticeR(int material) const
{
  try {
    return _material2min.at(material);
  }
  catch (std::out_of_range& ex) {
    std::vector<int> null;
    return null;
  }
}

template<typename T, unsigned D>
std::vector<int> BlockGeometryStatistics<T,D>::getMaxLatticeR(int material)
{
  update();
  return const_this->getMaxLatticeR(material);
}

template<typename T, unsigned D>
std::vector<int> BlockGeometryStatistics<T,D>::getMaxLatticeR(int material) const
{
  try {
    return _material2max.at(material);
  }
  catch (std::out_of_range& ex) {
    std::vector<int> null;
    return null;
  }
}

template<typename T, unsigned D>
std::vector<T> BlockGeometryStatistics<T,D>::getMinPhysR(int material) const
{
  std::vector<T> tmp(D,T());
  Vector<T,D> physR;
  _blockGeometry->getPhysR(physR, &(getMinLatticeR(material)[0]));
  for (unsigned iDim=0; iDim<D; iDim++) {
    tmp[iDim] = physR[iDim];
  }
  return tmp;
}

template<typename T, unsigned D>
std::vector<T> BlockGeometryStatistics<T,D>::getMaxPhysR(int material) const
{
  std::vector<T> tmp(D,T());
  Vector<T,D> physR;
  _blockGeometry->getPhysR(physR, &(getMaxLatticeR(material)[0]));
  for (unsigned iDim=0; iDim<D; iDim++) {
    tmp[iDim] = physR[iDim];
  }
  return tmp;
}

template<typename T, unsigned D>
std::vector<T> BlockGeometryStatistics<T,D>::getLatticeExtend(int material)
{
  update();
  return const_this->getLatticeExtend(material);
}

template<typename T, unsigned D>
std::vector<T> BlockGeometryStatistics<T,D>::getLatticeExtend(int material) const
{
  try {
    std::vector<T> extend;
    for (unsigned iDim = 0; iDim < D; iDim++) {
      extend.push_back(_material2max.at(material)[iDim] - _material2min.at(material)[iDim]);
    }
    return extend;
  }
  catch (std::out_of_range& ex) {
    std::vector<T> null;
    return null;
  }
}

template<typename T, unsigned D>
std::vector<T> BlockGeometryStatistics<T,D>::getPhysExtend(int material)
{
  update();
  return const_this->getPhysExtend(material);
}

template<typename T, unsigned D>
std::vector<T> BlockGeometryStatistics<T,D>::getPhysExtend(int material) const
{
  std::vector<T> extend;
  for (unsigned iDim = 0; iDim < D; iDim++) {
    extend.push_back(getMaxPhysR(material)[iDim] - getMinPhysR(material)[iDim]);
  }
  return extend;
}

template<typename T, unsigned D>
std::vector<T> BlockGeometryStatistics<T,D>::getPhysRadius(int material)
{
  update();
  return const_this->getPhysRadius(material);
}

template<typename T, unsigned D>
std::vector<T> BlockGeometryStatistics<T,D>::getPhysRadius(int material) const
{
  std::vector<T> radius;
  for (unsigned iDim=0; iDim<D; iDim++) {
    radius.push_back((getMaxPhysR(material)[iDim] - getMinPhysR(material)[iDim])/2.);
  }
  return radius;
}

template<typename T, unsigned D>
std::vector<T> BlockGeometryStatistics<T,D>::getCenterPhysR(int material)
{
  update();
  return const_this->getCenterPhysR(material);
}

template<typename T, unsigned D>
std::vector<T> BlockGeometryStatistics<T,D>::getCenterPhysR(int material) const
{
  std::vector<T> center;
  for (unsigned iDim=0; iDim<D; iDim++) {
    center.push_back(getMinPhysR(material)[iDim] + getPhysRadius(material)[iDim]);
  }
  return center;
}

template<typename T, unsigned D>
std::vector<int> BlockGeometryStatistics<T,D>::getType(const int* input) const
{
  return getType(Vector<int,D>(input));
}

template<typename T, unsigned D>
std::vector<int> BlockGeometryStatistics<T,D>::getType(Vector<int,D> latticeR,
                                                       BlockIndicatorF<T,D>& fluidI,
                                                       BlockIndicatorF<T,D>& outsideI) const
{
  const auto [normalType, normal] = computeBoundaryTypeAndNormal(fluidI, outsideI, latticeR);
  std::vector<int> discreteNormal;
  discreteNormal.push_back(static_cast<int>(normalType));
  for (unsigned iDim=0; iDim<D; iDim++) {
    discreteNormal.push_back(normal[iDim]);
  }
  return discreteNormal;
}

template<typename T, unsigned D>
std::vector<int> BlockGeometryStatistics<T,D>::getType(Vector<int,D> latticeR) const
{
  BlockIndicatorMaterial<T,D> fluidI(*_blockGeometry, 1);
  BlockIndicatorMaterial<T,D> outsideI(*_blockGeometry, 0);
  return getType(latticeR, fluidI, outsideI);
}

template<typename T, unsigned D>
std::vector<int> BlockGeometryStatistics<T,D>::getType(int iX, int iY) const requires (D == 2)
{
  return getType(Vector<int,2>{iX,iY});
}

template<typename T, unsigned D>
std::vector<int> BlockGeometryStatistics<T,D>::getType(int iX, int iY, int iZ) requires (D == 3)
{
  update();
  return const_this->getType(iX, iY, iZ);
}

template<typename T, unsigned D>
std::vector<int> BlockGeometryStatistics<T,D>::getType(int iX, int iY, int iZ) const requires (D == 3)
{
  return getType(Vector<int,3>{iX,iY,iZ});
}

template<typename T, unsigned D>
std::vector<int> BlockGeometryStatistics<T,D>::computeNormal(Vector<int,D> latticeR) const
{
  std::vector<int> normal (D,int(0));
  for (unsigned iDim=0; iDim<D; iDim++) {
    if (latticeR[iDim] != 0) {
      Vector<int,D> pos = latticeR;
      pos[iDim] -= 1;
      if (_blockGeometry->getMaterial(pos) == 1) {
        normal[iDim] = -1;
      }
    }
    if (latticeR[iDim] != _n[iDim] - 1) {
      Vector<int,D> pos = latticeR;
      pos[iDim] += 1;
      if (_blockGeometry->getMaterial(pos) == 1) {
        normal[iDim] = 1;
      }
    }
  }
  return normal;
}

template<typename T, unsigned D>
std::vector<T> BlockGeometryStatistics<T,D>::computeNormal(int material)
{
  update();
  return const_this->computeNormal(material);
}

template<typename T, unsigned D>
std::vector<T> BlockGeometryStatistics<T,D>::computeNormal(int material) const
{
  std::vector<T> normal (D,T());
  std::vector<int> minC = getMinLatticeR(material);
  std::vector<int> maxC = getMaxLatticeR(material);
  if (minC.size() == D && maxC.size() == D) {
    Vector<int,D> min(minC.data());
    Vector<int,D> max(maxC.data());
    Vector<int,D> pos = min;
    while (true) {
      if (_blockGeometry->getMaterial(pos) == material) {
        auto n = computeNormal(pos);
        for (unsigned iDim=0; iDim<D; iDim++) {
          normal[iDim] += n[iDim];
        }
      }
      unsigned iDim = 0;
      for (; iDim < D; iDim++) {
        if (++pos[iDim] <= max[iDim]) {
          break;
        }
        pos[iDim] = min[iDim];
      }
      if (iDim == D) {
        break;
      }
    }
  }
  T norm = T();
  for (unsigned iDim=0; iDim<D; iDim++) {
    norm += normal[iDim]*normal[iDim];
  }
  norm = util::sqrt(norm);
  if (norm>0.) {
    for (unsigned iDim=0; iDim<D; iDim++) {
      normal[iDim] /= norm;
    }
  }
  return normal;
}

template<typename T, unsigned D>
std::vector<int> BlockGeometryStatistics<T,D>::computeDiscreteNormal(int material, T maxNorm)
{
  update();
  return const_this->computeDiscreteNormal(material, maxNorm);
}

template<typename T, unsigned D>
std::vector<int> BlockGeometryStatistics<T,D>::computeDiscreteNormal(int material, T maxNorm) const
{
  std::vector<T> normal = computeNormal(material);
  std::vector<int> discreteNormal(D,int(0));

  T smallestAngle = T(0);
  // iterate the {-1,0,1}^D stencil with the last space direction varying fastest
  Vector<int,D> offset(-1);
  while (true) {
    T sqrNorm = T();
    T dot = T();
    for (unsigned iDim=0; iDim<D; iDim++) {
      sqrNorm += offset[iDim]*offset[iDim];
      dot += offset[iDim]*normal[iDim];
    }
    T norm = util::sqrt(sqrNorm);
    if (norm>0.&& norm<maxNorm) {
      T angle = dot/norm;
      if (angle>=smallestAngle) {
        smallestAngle=angle;
        for (unsigned iDim=0; iDim<D; iDim++) {
          discreteNormal[iDim] = offset[iDim];
        }
      }
    }
    unsigned iDim = D;
    while (iDim-- > 0) {
      if (++offset[iDim] <= 1) {
        break;
      }
      offset[iDim] = -1;
    }
    if (offset == Vector<int,D>(-1)) {
      break;
    }
  }
  return discreteNormal;
}

template<typename T, unsigned D>
bool BlockGeometryStatistics<T,D>::check(int material, int iX, int iY,
    unsigned offsetX, unsigned offsetY) requires (D == 2)
{
  return const_this->check(material, iX, iY, offsetX, offsetY);
}

template<typename T, unsigned D>
bool BlockGeometryStatistics<T,D>::check(int material, int iX, int iY,
    unsigned offsetX, unsigned offsetY) const requires (D == 2)
{
  bool found = true;
  for (int iOffsetX = -offsetX; iOffsetX <= (int) offsetX; ++iOffsetX) {
    for (int iOffsetY = -offsetY; iOffsetY <= (int) offsetY; ++iOffsetY) {
      if (_blockGeometry->getMaterial({iX + iOffsetX, iY + iOffsetY}) != material) {
        found = false;
      }
    }
  }
  return found;
}

template<typename T, unsigned D>
bool BlockGeometryStatistics<T,D>::check(int material, int iX, int iY, int iZ,
    unsigned offsetX, unsigned offsetY, unsigned offsetZ) requires (D == 3)
{
  update();
  return const_this->check(material, iX, iY, iZ, offsetX, offsetY, offsetZ);
}

template<typename T, unsigned D>
bool BlockGeometryStatistics<T,D>::check(int material, int iX, int iY, int iZ,
    unsigned offsetX, unsigned offsetY, unsigned offsetZ) const requires (D == 3)
{
  bool found = true;
  for (int iOffsetX = -offsetX; iOffsetX <= (int) offsetX; ++iOffsetX) {
    for (int iOffsetY = -offsetY; iOffsetY <= (int) offsetY; ++iOffsetY) {
      for (int iOffsetZ = -offsetZ; iOffsetZ <= (int) offsetZ; ++iOffsetZ) {
        if (_blockGeometry->getMaterial(LatticeR<3>{iX + iOffsetX, iY + iOffsetY,
          iZ + iOffsetZ}) != material) {
          found = false;
        }
      }
    }
  }
  return found;
}

template<typename T, unsigned D>
bool BlockGeometryStatistics<T,D>::find(int material, unsigned offsetX,
                                        unsigned offsetY, int& foundX, int& foundY) requires (D == 2)
{
  return const_this->find(material, offsetX, offsetY, foundX, foundY);
}

template<typename T, unsigned D>
bool BlockGeometryStatistics<T,D>::find(int material, unsigned offsetX,
                                        unsigned offsetY, int& foundX, int& foundY) const requires (D == 2)
{
  bool found = false;
  for (foundX = 0; foundX < _n[0]; foundX++) {
    for (foundY = 0; foundY < _n[1]; foundY++) {
      found = check(material, foundX, foundY, offsetX, offsetY);
      if (found) {
        return found;
      }
    }
  }
  return found;
}

template<typename T, unsigned D>
bool BlockGeometryStatistics<T,D>::find(int material, unsigned offsetX,
                                        unsigned offsetY, unsigned offsetZ, int& foundX, int& foundY,
                                        int& foundZ) requires (D == 3)
{
  update();
  return const_this->find(material, offsetX, offsetY, offsetZ, foundX, foundY, foundZ);
}

template<typename T, unsigned D>
bool BlockGeometryStatistics<T,D>::find(int material, unsigned offsetX,
                                        unsigned offsetY, unsigned offsetZ, int& foundX, int& foundY,
                                        int& foundZ) const requires (D == 3)
{
  bool found = false;
  for (foundX = 0; foundX < _n[0]; foundX++) {
    for (foundY = 0; foundY < _n[1]; foundY++) {
      for (foundZ = 0; foundZ < _n[2]; foundZ++) {
        found = check(material, foundX, foundY, foundZ, offsetX,
                      offsetY, offsetZ);
        if (found) {
          return found;
        }
      }
    }
  }
  return found;
}

template<typename T, unsigned D>
void BlockGeometryStatistics<T,D>::print()
{
  update();
  return const_this->print();
}

template<typename T, unsigned D>
void BlockGeometryStatistics<T,D>::print() const
{
  try {
    for (const auto& material : _material2n) {
      const auto& lo = _material2min.at(material.first);
      const auto& hi = _material2max.at(material.first);
      std::string minR = "(", maxR = "(";
      for (unsigned iDim=0; iDim<D; iDim++) {
        minR += (iDim>0 ? "," : "") + io::text::toCell(lo[iDim]);
        maxR += (iDim>0 ? "," : "") + io::text::toCell(hi[iDim]);
      }
      minR += ")"; maxR += ")";
      io::text::StatusLine()
        .field("material", material.first)
        .field("cells", io::text::shortCount(material.second))
        .field("minLatticeR", minR)
        .field("maxLatticeR", maxR)
        .render(clout);
    }
  }
  catch (std::out_of_range& ex)
  { }
}

template<typename T, unsigned D>
void BlockGeometryStatistics<T,D>::takeStatistics(LatticeR<D> latticeR)
{
  int type = _blockGeometry->getMaterial(latticeR);
  if (_material2n.count(type) == 0) {
    std::vector<int> co;
    for (unsigned iDim=0; iDim<D; iDim++) {
      co.push_back(latticeR[iDim]);
    }
    _material2n[type] = 1;
    _material2min[type] = co;
    _material2max[type] = co;
  }
  else {
    _material2n[type]++;
    for (unsigned iDim=0; iDim<D; iDim++) {
      if (latticeR[iDim] < _material2min[type][iDim]) {
        _material2min[type][iDim] = latticeR[iDim];
      }
      if (latticeR[iDim] > _material2max[type][iDim]) {
        _material2max[type][iDim] = latticeR[iDim];
      }
    }
  }
}

} // namespace olb

#endif
