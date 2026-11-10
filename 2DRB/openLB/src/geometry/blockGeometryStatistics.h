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
 * Representation of a statistic for a block geometry -- header file.
 */

#ifndef BLOCK_GEOMETRY_STATISTICS_H
#define BLOCK_GEOMETRY_STATISTICS_H

#include <vector>
#include <map>
#include <string>

#include "core/vector.h"
#include "io/ostreamManager.h"
#include "discreteNormals.h"

// All OpenLB code is contained in this namespace.
namespace olb {

/// Representation of a statistic for a block geometry
/** A block geomety statistic computes different integral
 * values, like total number of different materials,
 * materials of any kind, min./max. physical position, of an
 * underlying block geoemtry structure.
 *
 * This class is not intended to be derived from.
 */

template<typename T, unsigned D>
class BlockGeometry;

template<typename T, unsigned D>
class BlockGeometryStatistics {
private:
  /// Points to the underlying data from which the statistics is taken
  BlockGeometry<T,D>* _blockGeometry{};
  /// Specifies if an update is needed
  bool _statisticsUpdateNeeded{};
  /// Number of voxels in each direction (legacy, kept zero-initialized)
  Vector<int,D> _n{};

  /// Number of different material numbers
  int _nMaterials{};
  /// Mapping a material number to the number of this kind found in the super geometry
  std::map<int, std::size_t> _material2n{};
  /// Mapping a material number to the min. lattice position in each space direction
  std::map<int, std::vector<int> > _material2min{};
  /// Mapping a material number to the max. lattice position in each space direction
  std::map<int, std::vector<int> > _material2max{};

  /// class specific cout
  mutable OstreamManager clout;

  const BlockGeometryStatistics<T,D>* const_this{};

public:
  /// Constructor
  BlockGeometryStatistics(BlockGeometry<T,D>* blockGeometry);

  /// Read and write access to a flag, which indicates if an uptate is needed (=true)
  bool& getStatisticsStatus();
  /// Read only access to a flag, which indicates if an uptate is needed (=true)
  bool const & getStatisticsStatus() const;
  /// Returns the map with the numbers of voxels for each  material
  std::map<int, std::size_t> getMaterial2n();
  std::map<int, std::size_t> getMaterial2n() const;

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
  /// Returns the min. lattice position in each direction
  std::vector<int> getMinLatticeR(int material);
  std::vector<int> getMinLatticeR(int material) const;
  /// Returns the max. lattice position in each direction
  std::vector<int> getMaxLatticeR(int material);
  std::vector<int> getMaxLatticeR(int material) const;
  /// Returns the min. phys position in each direction
  std::vector<T> getMinPhysR(int material) const;
  /// Returns the max. phys position in each direction
  std::vector<T> getMaxPhysR(int material) const;
  /// Returns the lattice extend as length in each direction
  std::vector<T> getLatticeExtend(int material);
  std::vector<T> getLatticeExtend(int material) const;
  /// Returns the phys extend as length in each direction
  std::vector<T> getPhysExtend(int material);
  std::vector<T> getPhysExtend(int material) const;
  /// Returns the phys radius as length in each direction
  std::vector<T> getPhysRadius(int material);
  std::vector<T> getPhysRadius(int material) const;
  /// Returns the center position
  std::vector<T> getCenterPhysR(int material);
  std::vector<T> getCenterPhysR(int material) const;

  /* Returns the boundary type which is characterized by a discrete normal (c.f. Zimny)
   * 0 element:   flat, corner or edge
   * following elements: orientation -1, 0 or 1 per space direction
   */
  std::vector<int> getType(const int* input) const;
  std::vector<int> getType(Vector<int,D> latticeR, BlockIndicatorF<T,D>& fluidI, BlockIndicatorF<T,D>& outsideI) const;
  std::vector<int> getType(int iX, int iY) const requires (D == 2);
  std::vector<int> getType(int iX, int iY, int iZ) requires (D == 3);
  std::vector<int> getType(int iX, int iY, int iZ) const requires (D == 3);

  /// Returns normal that points into the fluid for paraxial surfaces
  std::vector<T> computeNormal (int material);
  std::vector<T> computeNormal (int material) const;
  /// Returns discrete normal with norm maxNorm that points into the fluid for paraxial surfaces
  /// maxNorm=1.1 implies only normals parallel to the axises
  std::vector<int> computeDiscreteNormal (int material, T maxNorm = 1.1);
  std::vector<int> computeDiscreteNormal (int material, T maxNorm = 1.1) const;

  /// Returns true if in a neighbourhood of size offset around latticeR only voxels of the given material are found
  bool check(int material, int iX, int iY, unsigned offsetX, unsigned offsetY) requires (D == 2);
  bool check(int material, int iX, int iY, unsigned offsetX, unsigned offsetY) const requires (D == 2);
  bool check(int material, int iX, int iY, int iZ, unsigned offsetX, unsigned offsetY, unsigned offsetZ) requires (D == 3);
  bool check(int material, int iX, int iY, int iZ, unsigned offsetX, unsigned offsetY, unsigned offsetZ) const requires (D == 3);
  /// Returns true and a position if there is a neighbourhood of size offset with only voxels of the given material
  bool find(int material, unsigned offsetX, unsigned offsetY, int& iX, int& iY) requires (D == 2);
  bool find(int material, unsigned offsetX, unsigned offsetY, int& iX, int& iY) const requires (D == 2);
  bool find(int material, unsigned offsetX, unsigned offsetY, unsigned offsetZ, int& iX, int& iY, int& iZ) requires (D == 3);
  bool find(int material, unsigned offsetX, unsigned offsetY, unsigned offsetZ, int& iX, int& iY, int& iZ) const requires (D == 3);

  /// Prints some statistic information, i.e. the number of voxels and min. max. physical position for each different material
  void print();
  void print() const;

private:
  /// Returns normal that points into the fluid for paraxial surfaces
  std::vector<int> computeNormal(Vector<int,D> latticeR) const;
  /// Returns the boundary type which is characterized by a discrete normal (c.f. Zimny)
  std::vector<int> getType(Vector<int,D> latticeR) const;
  /// Helper function to simplify the implementation
  void takeStatistics(LatticeR<D> latticeR);
};

} // namespace olb

#endif
