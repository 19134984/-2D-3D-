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

#ifndef REFINEMENT_MESH_H
#define REFINEMENT_MESH_H

#include "case/mesh.h"
#include "io/textTable.h"
#include "refinement/loadBalancer.h"

#include "communication/heuristicLoadBalancer.h"
#include "refinement/simpleLoadBalancer.h"
#include "io/xmlReader.h"

#include "utilities/log/log.h"

#include <iomanip>

#include <fstream>
#include <map>

namespace olb {

// Forward-declared: the octree mesh generator lives in src/contrib (the
// production mesh path is external via script/meshing). The RefinedMesh
// blueprint constructor below only performs dependent member calls, so it is
// instantiated solely by apps that include the contrib octree.
template <typename T, unsigned D>
class OctreeRefinedMeshBlueprint;

template <typename T, unsigned D>
class RefinedMesh {
private:
  std::vector<std::unique_ptr<Mesh<T,D>>> _meshes;
  std::unique_ptr<CuboidDecomposition<T,D>> _cDecompositionBase;

  void buildDecompositionBase()
  {
    auto& decomp0 = _meshes[0]->getCuboidDecomposition();
    _cDecompositionBase = std::make_unique<CuboidDecomposition<T,D>>(decomp0.getMotherCuboid(), 1);
    _cDecompositionBase->cuboids().clear();

    for (std::size_t iLevel = 0; iLevel < getLevelCount(); ++iLevel) {
      auto& decomp = getLevel(iLevel).getCuboidDecomposition();
      for (int iC = 0; iC < decomp.size(); ++iC) {
        auto& cuboid = decomp.get(iC);
        if (!cuboid.isCoarseBorder()) {
          int refIndex = _cDecompositionBase->cuboids().size();
          _cDecompositionBase->cuboids().push_back(cuboid);
          cuboid.setReferenceCuboid(refIndex);
        }
      }
    }
  }

public:
  RefinedMesh(std::vector<std::unique_ptr<CuboidDecomposition<T,D>>> decompositions,
              std::vector<std::unique_ptr<LoadBalancer<T>>> balancers)
  {
    OLB_SCOPE(logging::Setup, "construct refined mesh");
    if (decompositions.size() != balancers.size()) {
        throw std::runtime_error("Mismatched number of decompositions and load balancers.");
    }
    for (std::size_t iLevel = 0; iLevel < decompositions.size(); ++iLevel) {
      decompositions[iLevel]->buildNeighborCache(3);

      _meshes.push_back(std::unique_ptr<Mesh<T,D>>(new Mesh<T,D>(
        std::move(decompositions[iLevel]),
        std::move(balancers[iLevel])
      )));
      _meshes[iLevel]->setOverlap(3);
    }
    buildDecompositionBase();
  }

  RefinedMesh(std::string xmlFilename)
  {
    OLB_SCOPE(logging::Setup, "load refined mesh from XML");
    XMLreader xml(xmlFilename);

    const XMLreader* contentXml = nullptr;

    if (xml.getName() == "OpenLB") {
      // The root element is <OpenLB>, so its child is <Mesh>, and its grandchild is <XMLContent>
      contentXml = &(xml["Mesh"]["XMLContent"]);
    } else if (xml.getName() == "XMLContent") {
      // Legacy format: the root element is <XMLContent>
      contentXml = &xml;
    } else {
      // fallback
      contentXml = &(xml["XMLContent"]);
    }

    int maxRef = -1;
    for (const XMLreader* child : *contentXml) {
      std::string name = child->getName();
      if (name.starts_with("Level")) {
        int level = std::stoi(name.substr(5));
        maxRef = std::max(maxRef, level);
      }
    }

    if (maxRef < 0) {
      throw std::runtime_error("Invalid XML: No Level tags found.");
    }

    for (int i = 0; i <= maxRef; ++i) {
      std::string levelTag = "Level" + std::to_string(i);
      auto& decompXml = (*contentXml)[levelTag]["CuboidDecomposition"];
      auto cGeo = createCuboidDecomposition<T,D>(decompXml);

      cGeo->updateMotherCuboid();
      cGeo->buildNeighborCache(3);

      std::unique_ptr<LoadBalancer<T>> balancer;
      if (i == 0) {
        balancer = std::make_unique<PredefinedLoadBalancer<T>>(*cGeo);
      } else {
        // RefinedLoadBalancer required by the coupler via dynamic_cast, uses presetRank=true
        balancer = std::make_unique<RefinedLoadBalancer<T,D>>(
          _meshes[i - 1]->getCuboidDecomposition(),
          _meshes[i - 1]->getLoadBalancer(),
          *cGeo,
          true
        );
      }

      _meshes.push_back(std::unique_ptr<Mesh<T,D>>(new Mesh<T,D>(
        std::move(cGeo),
        std::move(balancer)
      )));
      _meshes[i]->setOverlap(3);
    }

    int maxProcess = -1;
    for (auto& mesh : _meshes) {
      for (auto& cuboid : mesh->getCuboidDecomposition().cuboids()) {
        if (cuboid.getPreCalculatedProcess().has_value()) {
          maxProcess = std::max(maxProcess, cuboid.getPreCalculatedProcess().value());
        }
      }
    }
    int expectedRanks = maxProcess + 1;
    if (expectedRanks > 0 && expectedRanks != singleton::mpi().getSize()) {
      throw std::runtime_error("Mismatch between imported XML ranks (" + std::to_string(expectedRanks) +
                               ") and MPI ranks (" + std::to_string(singleton::mpi().getSize()) + ").");
    }
    buildDecompositionBase();
  }

  RefinedMesh(OctreeRefinedMeshBlueprint<T,D>& refinedDecomposition)
  {
    OLB_SCOPE(logging::Setup, "construct refined mesh from octree blueprint");
    auto decompositions = refinedDecomposition.createCuboidDecomposition(false);
    auto decomposition0 = std::make_unique<CuboidDecomposition<T,D>>(decompositions[0]);
    decomposition0->buildNeighborCache(3);
    std::unique_ptr<LoadBalancer<T>> loadbalancer0 = std::make_unique<PredefinedLoadBalancer<T>>(*decomposition0);

    _meshes.push_back(std::unique_ptr<Mesh<T,D>>(new Mesh<T,D>(
      std::move(decomposition0),
      std::move(loadbalancer0)
    )));
    _meshes[0]->setOverlap(3);

    for (std::size_t iLevel = 1; iLevel <= refinedDecomposition.getMaxRefinementLevel(); ++iLevel) {
      auto decomposition = std::make_unique<CuboidDecomposition<T,D>>(decompositions[iLevel]);
      decomposition->buildNeighborCache(3);

      std::unique_ptr<LoadBalancer<T>> balancer = std::make_unique<RefinedLoadBalancer<T,D>>(
        _meshes[iLevel - 1]->getCuboidDecomposition(),
        _meshes[iLevel - 1]->getLoadBalancer(),
        *decomposition,
        true
      );

      _meshes.push_back(std::unique_ptr<Mesh<T,D>>(new Mesh<T,D>(
        std::move(decomposition),
        std::move(balancer)
      )));
      _meshes[iLevel]->setOverlap(3);
    }
    buildDecompositionBase();
  }

  Mesh<T,D>& getLevel(std::size_t iLevel)
  {
    return *_meshes[iLevel];
  }

  const Mesh<T,D>& getLevel(std::size_t iLevel) const
  {
    return *_meshes[iLevel];
  }

  std::size_t getLevelCount() const
  {
    return _meshes.size();
  }

  CuboidDecomposition<T,D>& getDecompositionBase()
  {
    return *_cDecompositionBase;
  }

  const CuboidDecomposition<T,D>& getDecompositionBase() const
  {
    return *_cDecompositionBase;
  }

  void print(OstreamManager clout = {std::cout, "RefinedMesh"})
  {
    using Align = io::text::Table::Align;
    io::text::Table table({
      {"Level",      Align::Right}, {"Resolution", Align::Right},
      {"Cells",      Align::Right}, {"Cuboids",    Align::Right},
      {"Imbalance",  Align::Right},
    });

    std::vector<double> weightedLoadPerRank(singleton::mpi().getSize(), 0.0);
    double totalCellsAllLevels = 0;
    std::size_t totalCuboidsAllLevels = 0;

    for (std::size_t iLevel = 0; iLevel < getLevelCount(); ++iLevel) {
      auto& mesh = getLevel(iLevel);
      auto& decomp = mesh.getCuboidDecomposition();
      auto& balancer = mesh.getLoadBalancer();

      std::vector<std::size_t> cellsPerRank(singleton::mpi().getSize(), 0);
      std::vector<std::size_t> weightsPerRank(singleton::mpi().getSize(), 0);
      for (std::size_t iLoc=0; iLoc < balancer.size(); ++iLoc) {
        cellsPerRank[singleton::mpi().getRank()] += decomp.get(balancer.glob(iLoc)).getLatticeVolume();
        weightsPerRank[singleton::mpi().getRank()] += decomp.get(balancer.glob(iLoc)).getWeight();
      }
#ifdef PARALLEL_MODE_MPI
      singleton::mpi().allReduceVect(cellsPerRank, MPI_SUM);
      singleton::mpi().allReduceVect(weightsPerRank, MPI_SUM);
#endif

      std::size_t maxWeights = 0;
      std::size_t totalCells = 0;
      double totalWeights = 0;
      for (std::size_t iRank = 0; iRank < weightsPerRank.size(); ++iRank) {
        maxWeights = std::max(maxWeights, weightsPerRank[iRank]);
        totalCells += cellsPerRank[iRank];
        totalWeights += (double)weightsPerRank[iRank];
        weightedLoadPerRank[iRank] += (double)weightsPerRank[iRank];
      }
      T imbalance = (totalWeights > 0) ? (T)(maxWeights / (totalWeights / singleton::mpi().getSize())) : 1;

      table.row(iLevel, io::text::sci((double)decomp.getDeltaR(), 4),
                io::text::shortCount(totalCells), decomp.size(),
                io::text::fixed((double)imbalance, 2));

      totalCellsAllLevels += (double)totalCells;
      totalCuboidsAllLevels += decomp.size();
    }

    double maxWeighted = 0.0;
    double totalWeighted = 0.0;
    for (auto n : weightedLoadPerRank) {
      maxWeighted = std::max(maxWeighted, n);
      totalWeighted += n;
    }
    T totalImbalance = (totalWeighted > 0.0) ? (T)(maxWeighted / (totalWeighted / singleton::mpi().getSize())) : 1;

    table.rule();
    table.row("total", "", io::text::shortCount((std::size_t)totalCellsAllLevels),
              totalCuboidsAllLevels, io::text::fixed((double)totalImbalance, 2));
    table.render(clout);
  }
};

}

#endif
