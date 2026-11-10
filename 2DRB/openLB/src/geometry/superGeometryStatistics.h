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
 * Representation of a statistic for a parallel geometry -- header file.
 */

#ifndef SUPER_GEOMETRY_STATISTICS_H
#define SUPER_GEOMETRY_STATISTICS_H


#include <map>
#include <string>
#include <vector>

#include "communication/mpiManager.h"
#include "core/vector.h"
#include "geometry/superGeometry.h"
#include "io/ostreamManager.h"


// All OpenLB code is contained in this namespace.
namespace olb {

/// Representation of a statistic for a parallel geometry
/** A super geometry statistic computes different integral
 * values, like total number of different materials,
 * materials of any kind, min./max. physical position, of an
 * underlying super geometry.
 *
 * This class is not intended to be derived from.
 */

template<typename T, unsigned D>
class SuperGeometry;

template<typename T, unsigned D>
class SuperGeometryStatistics {
private:
  /// Points to the underlying data from which the statistics is taken
  const SuperGeometry<T,D>* _superGeometry;
  /// Specifies if an update is needed
  bool _statisticsUpdateNeeded;

  /// Number of different material numbers
  int _nMaterials;
  /// Mapping a material number to the number of this kind found in the super geometry
  std::map<int, std::size_t> _material2n;
  /// Mapping a material number to the min. physical position in each space direction
  std::map<int, Vector<T, D> > _material2min;
  /// Mapping a material number to the max. physical position in each space direction
  std::map<int, Vector<T, D> > _material2max;
  /// Componentwise min extension over all material numbers. Default is 0.
  Vector<T, D> _minOverMaterial;
  /// Componentwise maximal extension over all material numbers. Default is 0.
  Vector<T, D> _maxOverMaterial;

  /// class specific cout
  mutable OstreamManager clout;

  const SuperGeometryStatistics<T,D>* const_this;

public:

  /// Constructor
  SuperGeometryStatistics(SuperGeometry<T,D>* superGeometry);
  /// Copy constructor
  SuperGeometryStatistics(SuperGeometryStatistics const& rhs);
  /// Copy assignment
  SuperGeometryStatistics<T,D>& operator=(SuperGeometryStatistics const& rhs);

  /// Read and write access to a flag, which indicates if an uptate is needed (=true)
  bool& getStatisticsStatus();
  /// Read only access to a flag, which indicates if an uptate is needed (=true)
  bool const & getStatisticsStatus() const;

  /// Updates the statistics if it is really needed
  void update(bool verbose=false);

  /// Returns the number of different materials
  int getNmaterials();
  int getNmaterials() const;
  /// Returns the number of voxels for a given material number
  std::size_t getNvoxel(int material);
  std::size_t getNvoxel(int material) const;
  /// Returns the number of voxels with material!=0
  std::size_t getNvoxel();
  std::size_t getNvoxel() const;
  /// Returns the min. phys position in each direction corresponding to material number
  Vector<T, D> getMinPhysR(int material);
  Vector<T, D> getMinPhysR(int material) const;
  /// Returns the min. phys position in each direction corresponding to all non-zero material numbers
  Vector<T, D> getMinPhysR();
  Vector<T, D> getMinPhysR() const;
  /// Returns the max. phys position in each direction corresponding to material number
  Vector<T, D> getMaxPhysR(int material);
  Vector<T, D> getMaxPhysR(int material) const;
  /// Returns the max. phys position in each direction corresponding to all non-zero material numbers
  Vector<T, D> getMaxPhysR();
  Vector<T, D> getMaxPhysR() const;
  /// Returns the phys extend as length in each direction
  Vector<T, D> getPhysExtend(int material);
  Vector<T, D> getPhysExtend(int material) const;
  /// Returns the phys radius as length in each direction
  Vector<T, D> getPhysRadius(int material);
  Vector<T, D> getPhysRadius(int material) const;
  /// Returns the center position
  Vector<T, D> getCenterPhysR(int material);
  Vector<T, D> getCenterPhysR(int material) const;
  /// Returns the boundary type which is characterized by a discrete normal (c.f. Zimny)
  Vector<int, D> getType(int iC, int iX, int iY) requires (D == 2);
  Vector<int, D> getType(int iC, int iX, int iY) const requires (D == 2);
  Vector<int, D> getType(int iC, int iX, int iY, int iZ) requires (D == 3);
  Vector<int, D> getType(int iC, int iX, int iY, int iZ) const requires (D == 3);

  /// Returns normal that points into the fluid for paraxial surfaces
  Vector<T, D> computeNormal (int material);
  Vector<T, D> computeNormal (int material) const;
  /// Returns discrete normal with norm maxNorm that points into the fluid for paraxial surfaces
  /// maxNorm=1.1 implies only normals parallel to the axises
  Vector<int, D> computeDiscreteNormal (int material, T maxNorm = 1.1);
  Vector<int, D> computeDiscreteNormal (int material, T maxNorm = 1.1) const;
  /// Returns util::sqrt( maxX^2 + maxY^2 + maxZ^2 ) max over a certain material number
  T computeMaxPhysDistance(int material) const;
  /// Returns util::sqrt( maxX^2 + maxY^2 + maxZ^2 ) max over all material numbers
  T computeMaxPhysDistance() const;
  /// Prints some statistic information, i.e. the number of voxels and min. max. physical position for each different material
  void print();
  void print() const;
};

} // namespace olb

#endif
