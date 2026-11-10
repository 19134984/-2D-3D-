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

#ifndef REFINEMENT_LATTICE_H
#define REFINEMENT_LATTICE_H

#include "mesh.h"
#include "geometry.h"
#include "io/textTable.h"

#include "fields.h"
#include "indicator.h"
#include "superLatticeRefinement.h"

#include "algorithm/rohde.h"

#include "utilities/log/log.h"

#include <iomanip>

namespace olb {

template <typename T, typename DESCRIPTOR>
class RefinedLattice {
public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

private:
  std::vector<std::unique_ptr<SuperLattice<T,DESCRIPTOR>>>           _lattices;
  std::vector<std::unique_ptr<SuperLatticeRefinement<T,DESCRIPTOR>>> _couplers;

public:
  RefinedLattice(RefinedMesh<T,DESCRIPTOR::d>& refinedMesh)
  {
    OLB_SCOPE(logging::Setup, "construct refined lattice");
    for (std::size_t iLevel = 0; iLevel < refinedMesh.getLevelCount(); ++iLevel) {
      auto& mesh = refinedMesh.getLevel(iLevel);
      _lattices.push_back(std::make_unique<SuperLattice<T,DESCRIPTOR>>(
        mesh.getCuboidDecomposition(),
        mesh.getLoadBalancer(),
        mesh.getOverlap()
      ));
    }
  }

  template <concepts::LatticeDescriptor _DESCRIPTOR>
  void setUnitConverter(const UnitConverter<T,_DESCRIPTOR>& converter)
  {
    _lattices[0]->setUnitConverter(converter);
    for (std::size_t iLevel=1; iLevel < _lattices.size(); ++iLevel) {
      _lattices[iLevel]->template setUnitConverter<DESCRIPTOR>(
        convectivelyRefineUnitConverter(_lattices[iLevel-1]->getUnitConverter(), 2));
    }

    OstreamManager clout(std::cout, "RefinedLattice");

    auto& conv0 = _lattices[0]->getUnitConverter();
    clout << io::text::section("Physical Parameters") << std::endl;
    io::text::KeyValueList physical;
    physical.item("Char. length [m]",     conv0.getCharPhysLength());
    physical.item("Char. velocity [m/s]", conv0.getCharPhysVelocity());
    physical.item("Viscosity [m^2/s]",    conv0.getPhysViscosity());
    physical.item("Density [kg/m^3]",     conv0.getPhysDensity());
    physical.item("Reynolds number",      conv0.getReynoldsNumber());
    physical.item("Mach number",          conv0.getMachNumber());
    physical.render(clout, "");

    clout << std::endl;
    clout << io::text::section("Levels") << std::endl;
    using Align = io::text::Table::Align;
    io::text::Table levels({
      {"Level",       Align::Right}, {"Delta X [m]", Align::Right},
      {"Delta T [s]", Align::Right}, {"Lattice U",   Align::Right},
      {"Tau",         Align::Right}, {"Omega",       Align::Right},
    });
    for (std::size_t iLevel=0; iLevel < _lattices.size(); ++iLevel) {
      auto& conv = _lattices[iLevel]->getUnitConverter();
      levels.row(iLevel, io::text::sci((double)conv.getPhysDeltaX(), 4),
                 io::text::sci((double)conv.getPhysDeltaT(), 4),
                 io::text::fixed((double)conv.getCharLatticeVelocity(), 4),
                 io::text::fixed((double)conv.getLatticeRelaxationTime(), 6),
                 io::text::fixed((double)conv.getLatticeRelaxationFrequency(), 6));
    }
    levels.render(clout);
    clout << std::endl;
  }

  SuperLattice<T,DESCRIPTOR>& getLevel(std::size_t iLevel)
  {
    if (iLevel >= _lattices.size()) {
       throw std::out_of_range("Lattice index out of range");
    }
    return *_lattices[iLevel];
  }

  SuperLatticeRefinement<T,DESCRIPTOR>& getCoupler(std::size_t iLevel)
  {
    return *_couplers[iLevel];
  }

  std::size_t getLevelCount() const
  {
    return _lattices.size();
  }

  /// Name each level's SuperLattice, suffixing _L<level> (forwards to statistics internally)
  void setName(const std::string& name)
  {
    for (std::size_t iLevel = 0; iLevel < getLevelCount(); ++iLevel) {
      _lattices[iLevel]->setName(name + "_L" + std::to_string(iLevel));
    }
  }

  void initialize(RefinedGeometry<T,DESCRIPTOR::d>& refinedGeometry,
                  std::vector<int> coarseBulkMaterials = {1, 11})
  {
    OLB_SCOPE(logging::Setup, "initialize refined lattice couplers");
    for (std::size_t iLevel = 0; iLevel < getLevelCount() - 1; ++iLevel) {
      _couplers.push_back(refinement::rohde::makeCoupler(
        *_lattices[iLevel],   refinedGeometry.getLevel(iLevel),
        *_lattices[iLevel+1], refinedGeometry.getLevel(iLevel + 1),
        refinedGeometry.getLevel(iLevel).getMaterialIndicator(
          std::vector<int>(coarseBulkMaterials))
      ));
      _lattices[iLevel]->setProcessingContext(ProcessingContext::Simulation);
    }
    _lattices[getLevelCount() - 1]->setProcessingContext(ProcessingContext::Simulation);
  }

  void setNoDynamics(RefinedGeometry<T,DESCRIPTOR::d>& refinedGeometry, int from, int to)
  {
    for (std::size_t iLevel = 0; iLevel < _lattices.size(); ++iLevel) {
      auto& geometry = refinedGeometry.getLevel(iLevel);
      CuboidDecomposition<T,DESCRIPTOR::d>& decomp = geometry.getCuboidDecomposition();
      for (Cuboid<T,DESCRIPTOR::d>& c : decomp.cuboids()) {
        if (c.isCoarseBorder()) {
          Vector<T,DESCRIPTOR::d> origin = c.getOrigin() - 1.5 * c.getDeltaR();
          Vector<T,DESCRIPTOR::d> extent = c.getExtent() * c.getDeltaR() + 2 * c.getDeltaR();
          IndicatorCuboid<T,DESCRIPTOR::d> cuboid(extent, origin);
          geometry.rename(from, to, cuboid);
        }
      }
      _lattices[iLevel]->template defineDynamics<NoDynamics<T,DESCRIPTOR>>(geometry, to);
    }
  }

  /// Fill fine-level halos from coarse before reading fields for output
  void synchronizeForPostProcessing()
  {
    for (std::size_t iLevel = 0; iLevel < getLevelCount() - 1; ++iLevel) {
      _lattices[iLevel]->getCommunicator(stage::Full{}).communicate();
    }
    for (std::size_t iLevel = 0; iLevel < getLevelCount() - 1; ++iLevel) {
      _couplers[iLevel]->apply(meta::id<refinement::rohde::CoarseToFineForPostProcessingO>{});
    }
  }
};

template <typename T>
struct is_refined_lattice {
  static constexpr bool value = false;
};

template <typename T, typename DESCRIPTOR>
struct is_refined_lattice<RefinedLattice<T, DESCRIPTOR>> {
  static constexpr bool value = true;
};

}

#endif
