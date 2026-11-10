/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013, 2014 Mathias J. Krause, 2026 Adrian Kummerlaender
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
 * Representation of a statistic for a parallel geometry -- generic implementation.
 */

#ifndef SUPER_GEOMETRY_STATISTICS_HH
#define SUPER_GEOMETRY_STATISTICS_HH

#include <iostream>
#include <math.h>
#include <sstream>
#include <iomanip>
#include <limits>
#include "io/textTable.h"
#include <cstdint>
#include <algorithm>
#include "utilities/omath.h"

#include "geometry/superGeometry.h"
#include "geometry/superGeometryStatistics.h"
#include "utilities/vectorHelpers.h"
#include "core/vector.h"
#include "core/olbDebug.h"


namespace olb {

template<typename T, unsigned D>
SuperGeometryStatistics<T,D>::SuperGeometryStatistics(SuperGeometry<T,D>* superGeometry)
  : _superGeometry(superGeometry), _statisticsUpdateNeeded(true),
    _nMaterials(0), _minOverMaterial(0), _maxOverMaterial(0),
    clout(std::cout,"SuperGeometryStatistics")
{
}

template<typename T, unsigned D>
SuperGeometryStatistics<T,D>::SuperGeometryStatistics(SuperGeometryStatistics const& rhs)
  : _superGeometry(rhs._superGeometry), _statisticsUpdateNeeded(true),
    _nMaterials(rhs._nMaterials),
    _minOverMaterial(rhs._minOverMaterial), _maxOverMaterial(rhs._maxOverMaterial),
    clout(std::cout,"SuperGeometryStatistics")
{
}

template<typename T, unsigned D>
SuperGeometryStatistics<T,D>& SuperGeometryStatistics<T,D>::operator=(SuperGeometryStatistics const& rhs)
{
  _superGeometry = rhs._superGeometry;
  _statisticsUpdateNeeded = true;
  _nMaterials = rhs._nMaterials;
  _minOverMaterial = rhs._minOverMaterial;
  _maxOverMaterial = rhs._maxOverMaterial;
  return *this;
}


template<typename T, unsigned D>
bool& SuperGeometryStatistics<T,D>::getStatisticsStatus()
{
  return _statisticsUpdateNeeded;
}

template<typename T, unsigned D>
bool const & SuperGeometryStatistics<T,D>::getStatisticsStatus() const
{
  return _statisticsUpdateNeeded;
}

template<typename T, unsigned D>
void SuperGeometryStatistics<T,D>::update(bool verbose)
{
  const_this = const_cast<const SuperGeometryStatistics<T,D>*>(this);

#ifdef PARALLEL_MODE_MPI
  int updateReallyNeededGlobal = 0;
  if (_statisticsUpdateNeeded) {
    updateReallyNeededGlobal = 1;
  }
  singleton::mpi().reduceAndBcast(updateReallyNeededGlobal, MPI_SUM);

  if (updateReallyNeededGlobal>0) {
    _statisticsUpdateNeeded = true;
  }
#endif

  // check if update is really needed
  if (_statisticsUpdateNeeded) {
    int updateReallyNeeded = 0;
    for (int iCloc=0; iCloc<_superGeometry->getLoadBalancer().size(); iCloc++) {
      if (_superGeometry->getBlockGeometry(iCloc).getStatistics().getStatisticsStatus()) {
        auto& blockGeometry = const_cast<BlockGeometry<T,D>&>(_superGeometry->getBlockGeometry(iCloc));
        blockGeometry.getStatistics().update(false);
        updateReallyNeeded++;
      }
    }

#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(updateReallyNeeded, MPI_SUM);
#endif

    if (updateReallyNeeded==0) {
      _statisticsUpdateNeeded = false;
      return;
    }

    // size the dense reduction buffers to the largest material number so that
    // the reduction stays a single fixed-size MPI_Allreduce without a fixed cap
    int maxMaterial = 0;
    for (int iCloc=0; iCloc<_superGeometry->getLoadBalancer().size(); iCloc++) {
      auto const& localMaterial2n = _superGeometry->getBlockGeometry(iCloc).getStatistics().getMaterial2n();
      if (!localMaterial2n.empty()) {
        maxMaterial = std::max(maxMaterial, localMaterial2n.rbegin()->first);
      }
    }
#ifdef PARALLEL_MODE_MPI
    singleton::mpi().reduceAndBcast(maxMaterial, MPI_MAX);
#endif

    constexpr int maxMaterialSanityBound = 4096;
    if (maxMaterial > maxMaterialSanityBound) {
      clout << "Warning: largest material number " << maxMaterial
            << " exceeds " << maxMaterialSanityBound
            << ", statistics reduction buffers will be correspondingly large" << std::endl;
    }

    int numMaterials = maxMaterial + 1;
    std::vector<std::uint64_t> materialCount(numMaterials, 0);
    std::vector<T> materialMin(D * numMaterials, std::numeric_limits<T>::max());
    std::vector<T> materialMax(D * numMaterials, std::numeric_limits<T>::lowest());

    for (int iCloc=0; iCloc<_superGeometry->getLoadBalancer().size(); iCloc++) {
      auto const& localMaterial2n = _superGeometry->getBlockGeometry(iCloc).getStatistics().getMaterial2n();
      for (auto const& [mat, count] : localMaterial2n) {
        if (count > 0) {
          materialCount[mat] += static_cast<std::uint64_t>(count);
          std::vector<T> minPhysR = _superGeometry->getBlockGeometry(iCloc).getStatistics().getMinPhysR(mat);
          std::vector<T> maxPhysR = _superGeometry->getBlockGeometry(iCloc).getStatistics().getMaxPhysR(mat);
          for (unsigned iDim=0; iDim<D; iDim++) {
            materialMin[D * mat + iDim] = std::min(materialMin[D * mat + iDim], minPhysR[iDim]);
            materialMax[D * mat + iDim] = std::max(materialMax[D * mat + iDim], maxPhysR[iDim]);
          }
        }
      }
    }

#ifdef PARALLEL_MODE_MPI
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Datatype typeT = (sizeof(T) == sizeof(double)) ? MPI_DOUBLE : MPI_FLOAT;

    MPI_Allreduce(MPI_IN_PLACE, materialCount.data(), numMaterials, MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm);
    MPI_Allreduce(MPI_IN_PLACE, materialMin.data(), D * numMaterials, typeT, MPI_MIN, comm);
    MPI_Allreduce(MPI_IN_PLACE, materialMax.data(), D * numMaterials, typeT, MPI_MAX, comm);
#endif

    _material2n.clear();
    _material2min.clear();
    _material2max.clear();
    _nMaterials = 0;

    for (int mat = 0; mat <= maxMaterial; mat++) {
      if (materialCount[mat] > 0) {
        _nMaterials++;
        _material2n[mat] = static_cast<std::size_t>(materialCount[mat]);
        Vector<T,D> minR, maxR;
        for (unsigned iDim=0; iDim<D; iDim++) {
          minR[iDim] = materialMin[D * mat + iDim];
          maxR[iDim] = materialMax[D * mat + iDim];
        }
        _material2min[mat] = minR;
        _material2max[mat] = maxR;
      }
    }

    // update _minOverMaterial and _maxOverMaterial
    for (unsigned iDim = 0; iDim < D; iDim++) {
      for (auto const& [mat, minVec] : _material2min) {
        if (mat != 0) {
          if (_minOverMaterial[iDim] > minVec[iDim]) {
            _minOverMaterial[iDim] = minVec[iDim];
          }
        }
      }
      for (auto const& [mat, maxVec] : _material2max) {
        if (mat != 0) {
          if (_maxOverMaterial[iDim] < maxVec[iDim]) {
            _maxOverMaterial[iDim] = maxVec[iDim];
          }
        }
      }
    }

    if (verbose) {
      clout << "updated" << std::endl;
    }
    _statisticsUpdateNeeded = false;
  }
}

template<typename T, unsigned D>
int SuperGeometryStatistics<T,D>::getNmaterials()
{
  update();
  return const_this->getNmaterials();
}

template<typename T, unsigned D>
int SuperGeometryStatistics<T,D>::getNmaterials() const
{
  return _nMaterials;
}

template<typename T, unsigned D>
std::size_t SuperGeometryStatistics<T,D>::getNvoxel(int material)
{
  update(true);
  return const_this->getNvoxel(material);
}

template<typename T, unsigned D>
std::size_t SuperGeometryStatistics<T,D>::getNvoxel(int material) const
{
  try {
    return _material2n.at(material);
  }
  catch (std::out_of_range& ex) {
    return 0;
  }
}

template<typename T, unsigned D>
std::size_t SuperGeometryStatistics<T,D>::getNvoxel()
{
  update();
  return const_this->getNvoxel();
}

template<typename T, unsigned D>
std::size_t SuperGeometryStatistics<T,D>::getNvoxel() const
{
  std::size_t total = 0;
  for (const auto& material : _material2n) {
    if (material.first!=0) {
      total += material.second;
    }
  }
  return total;
}

template<typename T, unsigned D>
Vector<T, D> SuperGeometryStatistics<T,D>::getMinPhysR(int material)
{
  update();
  return const_this->getMinPhysR(material);
}

template<typename T, unsigned D>
Vector<T, D> SuperGeometryStatistics<T,D>::getMinPhysR(int material) const
{
  try {
    return _material2min.at(material);
  }
  catch (std::out_of_range& ex) {
    return Vector<T,D>(T());
  }
}

template<typename T, unsigned D>
Vector<T, D> SuperGeometryStatistics<T,D>::getMinPhysR()
{
  update();
  return const_this->getMinPhysR();
}

template<typename T, unsigned D>
Vector<T, D> SuperGeometryStatistics<T,D>::getMinPhysR() const
{
  return _minOverMaterial;
}

template<typename T, unsigned D>
Vector<T, D> SuperGeometryStatistics<T,D>::getMaxPhysR(int material)
{
  update();
  return const_this->getMaxPhysR(material);
}

template<typename T, unsigned D>
Vector<T, D> SuperGeometryStatistics<T,D>::getMaxPhysR(int material) const
{
  try {
    return _material2max.at(material);
  }
  catch (std::out_of_range& ex) {
    return Vector<T,D>(T());
  }
}

template<typename T, unsigned D>
Vector<T, D> SuperGeometryStatistics<T,D>::getMaxPhysR()
{
  update();
  return const_this->getMaxPhysR();
}

template<typename T, unsigned D>
Vector<T, D> SuperGeometryStatistics<T,D>::getMaxPhysR() const
{
  return _maxOverMaterial;
}

template<typename T, unsigned D>
Vector<T, D> SuperGeometryStatistics<T,D>::getPhysExtend(int material)
{
  update();
  return const_this->getPhysExtend(material);
}

template<typename T, unsigned D>
Vector<T, D> SuperGeometryStatistics<T,D>::getPhysExtend(int material) const
{
  try {
    Vector<T,D> extend;
    for (unsigned iDim = 0; iDim < D; iDim++) {
      extend[iDim] = _material2max.at(material)[iDim] - _material2min.at(material)[iDim];
    }
    return extend;
  }
  catch (std::out_of_range& ex) {
    return Vector<T,D>(T());
  }
}

template<typename T, unsigned D>
Vector<T, D> SuperGeometryStatistics<T,D>::getPhysRadius(int material)
{
  update();
  return const_this->getPhysRadius(material);
}

template<typename T, unsigned D>
Vector<T, D> SuperGeometryStatistics<T,D>::getPhysRadius(int material) const
{
  Vector<T,D> radius;
  for (unsigned iDim=0; iDim<D; iDim++) {
    radius[iDim] = (getMaxPhysR(material)[iDim] - getMinPhysR(material)[iDim])/2.;
  }
  return radius;
}

template<typename T, unsigned D>
Vector<T, D> SuperGeometryStatistics<T,D>::getCenterPhysR(int material)
{
  update();
  return const_this->getCenterPhysR(material);
}

template<typename T, unsigned D>
Vector<T, D> SuperGeometryStatistics<T,D>::getCenterPhysR(int material) const
{
  Vector<T,D> center;
  for (unsigned iDim=0; iDim<D; iDim++) {
    center[iDim] = getMinPhysR(material)[iDim] + getPhysRadius(material)[iDim];
  }
  return center;
}

template<typename T, unsigned D>
Vector<int, D> SuperGeometryStatistics<T,D>::getType(int iC, int iX, int iY) requires (D == 2)
{
  return const_this->getType(iC, iX, iY);
}

template<typename T, unsigned D>
Vector<int, D> SuperGeometryStatistics<T,D>::getType(int iC, int iX, int iY) const requires (D == 2)
{
  int iCloc=_superGeometry->getLoadBalancer().loc(iC);
  Vector<int,D> discreteNormal = _superGeometry->getBlockGeometry(iCloc).getStatistics().getType(iX, iY);
  return discreteNormal;
}

template<typename T, unsigned D>
Vector<int, D> SuperGeometryStatistics<T,D>::getType(int iC, int iX, int iY, int iZ) requires (D == 3)
{
  update();
  return const_this->getType(iC, iX, iY, iZ);
}

template<typename T, unsigned D>
Vector<int, D> SuperGeometryStatistics<T,D>::getType(int iC, int iX, int iY, int iZ) const requires (D == 3)
{
  int iCloc=_superGeometry->getLoadBalancer().loc(iC);
  Vector<int,D> discreteNormal = _superGeometry->getBlockGeometry(iCloc).getStatistics().getType(iX, iY, iZ);
  return discreteNormal;
}

template<typename T, unsigned D>
Vector<T, D> SuperGeometryStatistics<T,D>::computeNormal(int material)
{
  update();
  return const_this->computeNormal(material);
}

template<typename T, unsigned D>
Vector<T, D> SuperGeometryStatistics<T,D>::computeNormal(int material) const
{
  std::vector<T> normal (D,T());
  for (int iCloc=0; iCloc<_superGeometry->getLoadBalancer().size(); iCloc++) {
    for (unsigned iDim=0; iDim<D; iDim++) {
      if (_superGeometry->getBlockGeometry(iCloc).getStatistics().getNvoxel(material)!=0) {
        normal[iDim] += _superGeometry->getBlockGeometry(iCloc).getStatistics().computeNormal(material)[iDim]*_superGeometry->getBlockGeometry(iCloc).getStatistics().getNvoxel(material);
      }
    }
  }

#ifdef PARALLEL_MODE_MPI
  for (unsigned iDim=0; iDim<D; iDim++) {
    singleton::mpi().reduceAndBcast((normal[iDim]), MPI_SUM);
  }
#endif

  std::size_t nVoxel = getNvoxel(material);
  if constexpr (D == 2) {
    if (nVoxel != 0) {
      for (unsigned iDim=0; iDim<D; iDim++) {
        normal[iDim] /= nVoxel;
      }
    }
    OLB_ASSERT(nVoxel || (normal[0] == 0 && normal[1] == 0), "if no voxels found we expect the normal to be zero");
  }
  else {
    if (nVoxel == 0) {
      throw std::runtime_error("Unknown material number: " + std::to_string(material));
    }
    for (unsigned iDim=0; iDim<D; iDim++) {
      normal[iDim] /= nVoxel;
    }
  }

  T norm = T();
  for (unsigned iDim=0; iDim<D; iDim++) {
    norm += normal[iDim]*normal[iDim];
  }
  norm = util::sqrt(norm);
  if (norm>0.) {
    for (unsigned iDim=0; iDim<D; iDim++) {
      normal[iDim]/=norm;
    }
  }
  return Vector<T,D>(normal);
}

template<typename T, unsigned D>
Vector<int, D> SuperGeometryStatistics<T,D>::computeDiscreteNormal(int material, T maxNorm)
{
  update();
  return const_this->computeDiscreteNormal(material, maxNorm);
}

template<typename T, unsigned D>
Vector<int, D> SuperGeometryStatistics<T,D>::computeDiscreteNormal(int material, T maxNorm) const
{
  Vector<T,D> normal = computeNormal(material);
  Vector<int, D> discreteNormal(0);

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
T SuperGeometryStatistics<T,D>::computeMaxPhysDistance( int material ) const
{
  Vector<T,D> vec;
  for (unsigned iDim=0; iDim<D; iDim++) {
    vec[iDim] = getMaxPhysR(material)[iDim] - getMinPhysR(material)[iDim];
  }
  return norm(vec);
}

template<typename T, unsigned D>
T SuperGeometryStatistics<T,D>::computeMaxPhysDistance() const
{
  Vector<T,D> vec;
  for (unsigned iDim=0; iDim<D; iDim++) {
    vec[iDim] = getMaxPhysR()[iDim] - getMinPhysR()[iDim];
  }
  return norm(vec);
}

template<typename T, unsigned D>
void SuperGeometryStatistics<T,D>::print()
{
  update();
  return const_this->print();
}

template<typename T, unsigned D>
void SuperGeometryStatistics<T,D>::print() const
{
  try {
    std::size_t nCells = 0;
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
        .field("minPhysR", minR)
        .field("maxPhysR", maxR)
        .render(clout);
      nCells += material.second;
    }
    io::text::StatusLine().field("total cells", io::text::shortCount(nCells)).render(clout);
  }
  catch (std::out_of_range& ex) {
  }
}

} // namespace olb

#endif
