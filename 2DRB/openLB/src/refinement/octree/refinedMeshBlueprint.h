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

#ifndef OLB_REFINEMENT_OCTREE_REFINEDMESHBLUEPRINT_H
#define OLB_REFINEMENT_OCTREE_REFINEDMESHBLUEPRINT_H

#include "octant.h"
#include "root.h"

#include <filesystem>
#include <cstdint>
#include <cmath>
#include <cstdlib>

namespace olb {

template <typename T, unsigned D>
class OctreeRefinedMeshBlueprint {
private:
  Cuboid<T, D> _motherCuboid;

  std::vector<OctreeRoot<T, D>> _roots;

#ifdef PARALLEL_MODE_MPI
  std::vector<int>                neighborIDs;
  std::vector<int>                neighborRanks;
  singleton::MpiNonBlockingHelper _mpiNbHelper;
#endif

  int _maxRefinementLevel;
  int _rootSize;

public:
  OctreeRefinedMeshBlueprint(IndicatorF<T, D>& indicatorF,
                             T voxelSize,
                             int rootSize,
                             int maxRefinementLevel,
                             bool hilbert = false)
  {
    _maxRefinementLevel = maxRefinementLevel;
    _rootSize = rootSize;

    if (rootSize % 2 != 0 && maxRefinementLevel != 0) {
      throw std::logic_error("rootSize needs to be a multiple of 2");
    }

    const T blockExtent = voxelSize * rootSize;
    Vector<T,D> domainMin = indicatorF.getMin();
    Vector<T,D> domainExtent = indicatorF.getMax() - domainMin;

    Vector<int,D> forrestExtents;
    Vector<T,D> minV;

    for (unsigned iD=0; iD < D; ++iD) {
      T ratio = domainExtent[iD] / blockExtent;
      if (util::abs(ratio - util::round(ratio)) < 1e-5) {
        ratio = util::round(ratio) + 0.1;
      }
      forrestExtents[iD] = util::ceil(ratio) + 2;
      T totalForestLength = forrestExtents[iD] * blockExtent;
      minV[iD] = domainMin[iD] - (totalForestLength - domainExtent[iD]) * 0.5;
    }
    _motherCuboid = Cuboid<T,D>(minV, voxelSize, forrestExtents * rootSize);

    int maxExtent = max(forrestExtents);
    int logicalRef  = 0;

    maxExtent--;
    while (maxExtent) {
      maxExtent >>= 1;
      logicalRef++;
    }

    int max_order;
    if constexpr (D == 2) {
      if (hilbert) {
        max_order = util::pow(4,logicalRef)-1;
      }
      else {
        max_order = Octant<T,D>::calcIDFast({forrestExtents[0]-1, forrestExtents[1]-1});
      }
    }
    else {
      if (hilbert) {
        max_order = util::pow(8,logicalRef)-1;
      }
      else {
        max_order = Octant<T,D>::calcIDFast({forrestExtents[0]-1, forrestExtents[1]-1, forrestExtents[2]-1});
      }
    }
    int low  = 0;
    int high = max_order;

#ifdef PARALLEL_MODE_MPI
    std::vector<int> ranges;
    int              rank    = singleton::mpi().getRank();
    int              size    = singleton::mpi().getSize();
    int              perRank = (max_order+1) / size;

    if (max_order+1 < singleton::mpi().getSize()) {
      for (int i = 0; i <= max_order; i++) {
        ranges.push_back(i);
      }
      if (rank > ranges.size()) {
        _roots ={};
        neighborRanks = {};
        neighborIDs = {};
        return;
      }
      low = singleton::mpi().getRank();
      high = singleton::mpi().getRank();
    }
    else {
      for (int i = 1; i <= size; i++) {
        ranges.push_back(perRank * (i)-1);
      }
      ranges.back() = high;
      low  = rank == 0 ? low : ranges[rank - 1] + 1;
      high = rank == (size - 1) ? high : ranges[rank];
    }

#endif
    Vector<int, D> origin;
    for (std::uint64_t i = low; i <= high; i++) {

      if (hilbert) {
        origin = Octant<T, D>::calcOriginFast(Octant<T,D>::calcMortonID(i, logicalRef));
      }
      else {
        origin = Octant<T, D>::calcOriginFast(i);
      }

      if constexpr (D==3) {
        if (origin[0] >= forrestExtents[0] || origin[1] >= forrestExtents[1] || origin[2] >= forrestExtents[2]) {
          continue;
        }
      }
      else {
        if (origin[0] >= forrestExtents[0] || origin[1] >= forrestExtents[1]) {
          continue;
        }
      }


      Vector<T, D> offset = origin * blockExtent;
      _roots.emplace_back(OctreeRoot<T, D>(minV + offset, voxelSize, rootSize, maxRefinementLevel, i, logicalRef, forrestExtents, hilbert));
    }
#ifndef PARALLEL_MODE_MPI
    std::vector<std::uint64_t> linearized;
    for (int i = _roots.size() - 1; i >= 0; i--) {
      if (!_roots[i].intersect(indicatorF)) {
        linearized.push_back(_roots[i].getID());
        _roots.erase(_roots.begin() + i);
      }
    }
#else

    for (auto& root : _roots) {
      for (auto neighbor : root.getNeighbors()) {
        if (!neighbor.has_value()) {
          continue;
        }
        auto lb = std::lower_bound(ranges.begin(), ranges.end(), neighbor.value());
        if (lb == ranges.end()) {
          throw std::logic_error("should not happen");
        }
        int  rank   = std::distance(ranges.begin(), lb);
        auto lbRank = std::lower_bound(neighborRanks.begin(), neighborRanks.end(), rank);
        if (lbRank == neighborRanks.end() || *lbRank != rank) {
          neighborRanks.insert(lbRank, rank);
        }
      }
    }
    std::ranges::sort(neighborRanks);
    for (auto x : neighborRanks) {
      neighborIDs.push_back(ranges[x]);
    }

    std::vector<std::vector<std::uint64_t>> toSend(neighborRanks.size());
    for (int i = _roots.size() - 1; i >= 0; i--) {
      if (!_roots[i].intersect(indicatorF)) {
        for (auto n : _roots[i].getNeighbors()) {
          if (n.has_value()) {
            int index =
                std::distance(neighborIDs.begin(), std::lower_bound(neighborIDs.begin(), neighborIDs.end(), n.value()));
            auto z = _roots[i].getID();
            if (!std::binary_search(toSend[index].begin(), toSend[index].end(), _roots[i].getID())) {
              toSend[index].push_back(_roots[i].getID());
            }

          }
        }
        _roots.erase(_roots.begin() + i);
      }
    }
    for (auto& x : toSend) {
      std::reverse(x.begin(), x.end());
    }

    std::vector<std::vector<std::uint64_t>> erased = rank2Rank<std::uint64_t>(std::move(toSend));

    //linearize using k way sort
    std::vector<std::vector<std::uint64_t>::iterator> heads;

    int totalSize = 0;
    for (auto x : erased) {
      totalSize += x.size();
    }
    std::vector<std::uint64_t> linearized;
    linearized.reserve(totalSize);

    std::multimap<std::uint64_t, int> kHeap;
    for (int i = 0; i < erased.size(); i++) {
      if (erased[i].size() != 0) {
        kHeap.insert({*erased[i].begin(), i});
      }
      heads.push_back(erased[i].begin());
    }

    while (linearized.size() != totalSize) {
      int index = kHeap.begin()->second;
      linearized.push_back(kHeap.begin()->first);
      ++heads[index];
      if (heads[index] != erased[index].end()) {
        kHeap.insert({*heads[index], index});
      }
      kHeap.erase(kHeap.begin());
    }

#endif
#ifdef PARALLEL_MODE_OMP
#pragma omp parallel for
#endif
    for (std::size_t i = 0; i < _roots.size(); ++i) {
      _roots[i].removeNeighbors(linearized);
    }
#ifdef PARALLEL_MODE_MPI
    if (_roots.empty()) {
      neighborRanks.clear();
      neighborIDs.clear();
      return;
    }
    std::vector<bool> checked(neighborRanks.size());
    for (int i = 0; i < _roots.size(); i++) {
      for (auto& neighbor: _roots[i].getNeighbors()) {
        if (neighbor.has_value()) {
          auto lb = std::lower_bound(neighborIDs.begin(), neighborIDs.end(), neighbor.value());
          checked[std::distance(neighborIDs.begin(), lb)] = true;
        }
      }
    }
    checked[std::distance(neighborRanks.begin(), std::find(neighborRanks.begin(), neighborRanks.end(), singleton::mpi().getRank()))] = true;
    for (int i = neighborRanks.size()-1; i >= 0; i--) {
      if (!checked[i]) {
        neighborRanks.erase(neighborRanks.begin() + i);
        neighborIDs.erase(neighborIDs.begin() + i);
      }
    }
#endif
  }

  void refine(IndicatorF<T, D>& indicatorF, int refinementLevel, int padding = 1)
  {
#ifdef PARALLEL_MODE_OMP
#pragma omp parallel for
#endif
    for (std::size_t i = 0; i < _roots.size(); ++i) {
      _roots[i].refine(indicatorF, refinementLevel, padding);
    }
  }

  void reduceWeights(IndicatorF<T, D>& indicatorF)
  {
#ifdef PARALLEL_MODE_OMP
#pragma omp parallel for
#endif
    for (std::size_t i = 0; i < _roots.size(); ++i) {
      _roots[i].reduceWeights(indicatorF);
    }
  }

  void balanceGraph(int amountOfProcesses,std::string unique, bool discount)
  {
    std::uint64_t localOctants = 0;
    std::uint64_t acc = 0;
    std::uint64_t totalOctants = 0;

    for (auto& root : _roots) {
      localOctants += root.getOctants().size();
    }

#ifdef PARALLEL_MODE_MPI
    int rank = singleton::mpi().getRank();
    int size = singleton::mpi().getSize();

    std::vector<std::uint64_t> allLocalBlocks(size);
    singleton::mpi().allGather(
        &localOctants, 1,
        allLocalBlocks.data(), 1
    );

    for (int i = 0; i < rank; ++i) {
      acc += allLocalBlocks[i];
    }
    for (int i = 0; i < size; ++i) {
      totalOctants += allLocalBlocks[i];
    }

#endif

    for (auto& root : _roots) {
      acc = root.number(acc);
    }

    std::map<int, std::vector<Octant<T, D>>> comms;
    std::string dirname = "./tmp/graphData/"+unique+"/";
    std::string filename = "0";

#ifdef PARALLEL_MODE_MPI
    filename = std::to_string(singleton::mpi().getRank());
    std::vector<std::vector<std::uint64_t>>       toSendID;
    std::vector<std::vector<int>>            toSendRef;
    std::vector<std::vector<int>>            toSendOriginalOctant;
    std::vector<std::vector<std::uint64_t>>       toSendGlobalID;
    std::vector<std::uint64_t>                    receivers;
    std::vector<std::vector<bool>>           checked(neighborRanks.size(), std::vector<bool>(_roots.empty() ? 0: _roots.back().getID() - _roots.begin()->getID() +1, false));
#endif


    for (auto& root : _roots) {

      comms[root.getID()] = root.getBorderOctant();
#ifdef PARALLEL_MODE_MPI
      for (auto neighbor : root.getNeighbors()) {
        if (neighbor.has_value()) {
          auto lb = std::lower_bound(neighborIDs.begin(), neighborIDs.end(), neighbor.value());

          if (neighborRanks[std::distance(neighborIDs.begin(), lb)] != singleton::mpi().getRank() &&
              !checked[std::distance(neighborIDs.begin(), lb)][root.getID() - _roots[0].getID()]) {

            std::vector<std::uint64_t> vec;
            std::vector<int>      ref;
            std::vector<int>      oo;
            std::vector<std::uint64_t> gvec;
            for (const auto& x : comms[root.getID()]) {
              vec.push_back(x.getId());
              ref.push_back(x.getRefinementLevel());
              oo.push_back(root.getID());
              gvec.push_back(x.getGlobalID());
            }
            toSendGlobalID.push_back(gvec);
            toSendID.push_back(vec);
            toSendRef.push_back(ref);
            toSendOriginalOctant.push_back(oo);
            receivers.push_back(neighbor.value());
            checked[std::distance(neighborIDs.begin(), lb)][root.getID() - _roots[0].getID()] = true;
              }
        }
      }
#endif
    }

#ifdef PARALLEL_MODE_MPI

    auto [receivedIDs, _]          = root2Root(toSendID, receivers);
    auto [receivedRefinements, __] = root2Root(toSendRef, receivers);
    auto [originalOctant, ___]     = root2Root(toSendOriginalOctant, receivers);
    auto [receivedGlobalID, ____]     = root2Root(toSendGlobalID, receivers);

    for (int i = 0; i < originalOctant.size(); i++) {
      for (int j = 0; j < originalOctant[i].size(); j++) {
        comms[originalOctant[i][j]].push_back(
            Octant<T, D>(Octant<T, D>::calcOriginFast(receivedIDs[i][j]), receivedRefinements[i][j], _maxRefinementLevel, std::make_optional(receivedGlobalID[i][j])));
      }
    }
    if (singleton::mpi().isMainProcessor()) {
      std::filesystem::create_directories("./tmp/graphData/"+unique+"/data");
    }
    singleton::mpi().barrier();
#else
    std::filesystem::create_directories("./tmp/graphData/");
#endif


    std::ofstream file(dirname+"data/" + filename + "graph.tmp");
    file.close();
    for (auto& root : _roots) {
      std::map<int, std::vector<Octant<T, D>>> local;
      auto                                     neighbors = root.getNeighbors();
      for (int i = 0; i < neighbors.size(); i++) {
        if (neighbors[i].has_value()) {
          local[i] = comms[neighbors[i].value()];
        }
        else {
          local[i] = {};
        }
      }
      root.balanceGraph(local, dirname+"data/" + filename + "graph.tmp", discount);
    }

#ifdef PARALLEL_MODE_MPI
    singleton::mpi().barrier();
    if (singleton::mpi().isMainProcessor()) {
      std::ofstream merged(dirname + "data/graph.tmp");
      for (int i = 0; i < singleton::mpi().getSize(); i++) {
        std::ifstream curFile(dirname + "data/" +std::to_string(i) + "graph.tmp");
        curFile.seekg(0, std::ios::end);
        std::streamsize size = curFile.tellg();
        curFile.seekg(0, std::ios::beg);
        if (size > 0) {
          merged << curFile.rdbuf();
        }
      }
      merged.close();
    }

    if (singleton::mpi().isMainProcessor()) {
#endif

      std::ifstream written(dirname + "data/graph.tmp");
      std::ofstream graph(dirname + "graph.txt", std::ios::out | std::ios::trunc);

      int         number_of_nodes = 0;
      int         whitespaces     = 0;
      std::string line;
      while (std::getline(written, line)) {
        ++number_of_nodes;
        for (char x : line) {
          if (x == ' ') {
            ++whitespaces;
          }
        }
      }

      graph << number_of_nodes << " " << whitespaces / 4 << " 11" << std::endl;
      written.clear();
      written.seekg(0, std::ios::beg);
      graph << written.rdbuf();

      written.close();
      graph.close();

      //fast eco strong
      std::string command = "kaffpa "+dirname + "graph.txt --k=" + std::to_string(amountOfProcesses) +
                            " --preconfiguration=strong --output_filename="+dirname+"partition.txt";
      std::cout << command << std::endl;
      std::system(command.c_str());
#ifdef PARALLEL_MODE_MPI
    }
    singleton::mpi().barrier();
#endif
    std::ifstream partitions(dirname+"partition.txt");
    std::string _lines;
    for (std::uint64_t i = 0; i < acc - localOctants; ++i) {
      std::getline(partitions, _lines);
    }

    for (auto& root : _roots) {
      root.assignProcess(partitions);
    }
    partitions.close();
  }

  void loadBalanceGreedyLevel(int amountOfProcesses)
  {
    std::vector<int>              processed;
    std::vector<std::vector<int>> weightPerProcess(_maxRefinementLevel + 1, std::vector<int>(amountOfProcesses));
    std::vector<std::vector<Octant<T, D>*>> levelOctants(_maxRefinementLevel + 1);
    for (auto& root : _roots) {
      for (auto& octant : root.getOctants()) {
        levelOctants[octant.getRefinementLevel()].push_back(&octant);
      }
    }
    for (auto& level : levelOctants) {
      std::sort(level.begin(), level.end(), [](const Octant<T, D>* a, const Octant<T, D>* b) {
        return a->calcWeight() > b->calcWeight();
      });
    }

    for (int i = 0; i <= _maxRefinementLevel; i++) {
      for (auto octant : levelOctants[i]) {
        int inx = std::distance(weightPerProcess[i].begin(),
                                std::min_element(weightPerProcess[i].begin(), weightPerProcess[i].end()));
        octant->setProcess(inx);
        weightPerProcess[i][inx] += octant->calcWeight();
      }
    }
  }

  void transformHilbert()
  {
#ifdef PARALLEL_MODE_OMP
#pragma omp parallel for
#endif
    for (std::size_t i = 0; i < _roots.size(); ++i) {
      _roots[i].transformHilbert();
    }
  }

  void transformMorton()
  {
#ifdef PARALLEL_MODE_OMP
#pragma omp parallel for
#endif
    for (std::size_t i = 0; i < _roots.size(); ++i) {
      _roots[i].transformMorton();
    }
  }

  void outputResult(std::string unique)
  {
    std::uint64_t localOctants = 0;
    std::uint64_t acc = 0;
    std::uint64_t totalOctants = 0;

    for (auto& root : _roots) {
      localOctants += root.getOctants().size();
    }

#ifdef PARALLEL_MODE_MPI
    int rank = singleton::mpi().getRank();
    int size = singleton::mpi().getSize();

    std::vector<std::uint64_t> allLocalBlocks(size);
    singleton::mpi().allGather(
        &localOctants, 1,
        allLocalBlocks.data(), 1
    );

    for (int i = 0; i < rank; ++i) {
      acc += allLocalBlocks[i];
    }
    for (int i = 0; i < size; ++i) {
      totalOctants += allLocalBlocks[i];
    }

#endif

    for (auto& root : _roots) {
      acc = root.number(acc);
    }

    std::map<int, std::vector<Octant<T, D>>> comms;
    std::string dirname = "./tmp/graphData/"+unique+"/";
    std::string filename = "0";

#ifdef PARALLEL_MODE_MPI
    filename = std::to_string(singleton::mpi().getRank());
    std::vector<std::vector<std::uint64_t>>       toSendID;
    std::vector<std::vector<int>>            toSendRef;
    std::vector<std::vector<int>>            toSendOriginalOctant;
    std::vector<std::vector<std::uint64_t>>       toSendGlobalID;
    std::vector<std::uint64_t>                    receivers;
    std::vector<std::vector<bool>>           checked(neighborRanks.size(), std::vector<bool>(_roots.empty() ? 0: _roots.back().getID() - _roots.begin()->getID() +1, false));
#endif


    for (auto& root : _roots) {

      comms[root.getID()] = root.getBorderOctant();
#ifdef PARALLEL_MODE_MPI
      for (auto neighbor : root.getNeighbors()) {
        if (neighbor.has_value()) {
          auto lb = std::lower_bound(neighborIDs.begin(), neighborIDs.end(), neighbor.value());

          if (neighborRanks[std::distance(neighborIDs.begin(), lb)] != singleton::mpi().getRank() &&
              !checked[std::distance(neighborIDs.begin(), lb)][root.getID() - _roots[0].getID()]) {

            std::vector<std::uint64_t> vec;
            std::vector<int>      ref;
            std::vector<int>      oo;
            std::vector<std::uint64_t> gvec;
            for (const auto& x : comms[root.getID()]) {
              vec.push_back(x.getId());
              ref.push_back(x.getRefinementLevel());
              oo.push_back(root.getID());
              gvec.push_back(x.getGlobalID());
            }
            toSendGlobalID.push_back(gvec);
            toSendID.push_back(vec);
            toSendRef.push_back(ref);
            toSendOriginalOctant.push_back(oo);
            receivers.push_back(neighbor.value());
            checked[std::distance(neighborIDs.begin(), lb)][root.getID() - _roots[0].getID()] = true;
              }
        }
      }
#endif
    }

#ifdef PARALLEL_MODE_MPI

    auto [receivedIDs, _]          = root2Root(toSendID, receivers);
    auto [receivedRefinements, __] = root2Root(toSendRef, receivers);
    auto [originalOctant, ___]     = root2Root(toSendOriginalOctant, receivers);
    auto [receivedGlobalID, ____]     = root2Root(toSendGlobalID, receivers);

    for (int i = 0; i < originalOctant.size(); i++) {
      for (int j = 0; j < originalOctant[i].size(); j++) {
        comms[originalOctant[i][j]].push_back(
            Octant<T, D>(Octant<T, D>::calcOriginFast(receivedIDs[i][j]), receivedRefinements[i][j], _maxRefinementLevel, std::make_optional(receivedGlobalID[i][j])));
      }
    }
    if (singleton::mpi().isMainProcessor()) {
      std::filesystem::create_directories("./tmp/graphData/"+unique+"/data");
    }
    singleton::mpi().barrier();
#else
    std::filesystem::create_directories("./tmp/graphData/");
#endif


    std::ofstream file(dirname+"data/" + filename + "graph.tmp");
    std::ofstream file2(dirname+"data/" + filename + "partition.tmp");
    file.close();
    file2.close();
    for (auto& root : _roots) {
      std::map<int, std::vector<Octant<T, D>>> local;
      auto                                     neighbors = root.getNeighbors();
      for (int i = 0; i < neighbors.size(); i++) {
        if (neighbors[i].has_value()) {
          local[i] = comms[neighbors[i].value()];
        }
        else {
          local[i] = {};
        }
      }

      root.balanceGraph(local, dirname+"data/" + filename + "graph.tmp", false);
      root.writePartition(dirname+"data/" + filename + "partition.tmp");
    }

#ifdef PARALLEL_MODE_MPI
    singleton::mpi().barrier();
    if (singleton::mpi().isMainProcessor()) {
      std::ofstream merged(dirname + "data/graph.tmp");
      std::ofstream mergedP(dirname + "partition.txt");

      for (int i = 0; i < singleton::mpi().getSize(); i++) {
        std::ifstream curFile(dirname + "data/" +std::to_string(i) + "graph.tmp");
        std::ifstream curFileP(dirname + "data/" + std::to_string(i) + "partition.tmp");

        curFile.seekg(0, std::ios::end);

        std::streamsize size = curFile.tellg();

        curFile.seekg(0, std::ios::beg);
        curFileP.seekg(0, std::ios::beg);

        if (size > 0) {
          merged << curFile.rdbuf();
          mergedP << curFileP.rdbuf();
        }
      }
      merged.close();
      mergedP.close();
    }
#endif
    if (singleton::mpi().isMainProcessor()) {
      std::ifstream written(dirname + "data/graph.tmp");
      std::ofstream graph(dirname + "graph.txt", std::ios::out | std::ios::trunc);

      int         number_of_nodes = 0;
      int         whitespaces     = 0;
      std::string line;
      while (std::getline(written, line)) {
        ++number_of_nodes;
        for (char x : line) {
          if (x == ' ') {
            ++whitespaces;
          }
        }
      }

      graph << number_of_nodes << " " << whitespaces / 4 << " 11" << std::endl;
      written.clear();
      written.seekg(0, std::ios::beg);
      graph << written.rdbuf();

      written.close();
      graph.close();
  }
  }

  std::vector<CuboidDecomposition<T, D>> createCuboidDecomposition(bool heuristic = false)
  {
    std::vector<std::vector<Cuboid<T, D>>> decompositions;
    std::vector<T>                         voxelSizes;
    for (int i = 0; i <= _maxRefinementLevel; i++) {
      decompositions.push_back({});
      voxelSizes.push_back(voxelSizes.size() ? voxelSizes[i - 1] / 2 : _motherCuboid.getDeltaR());
    }

    for (auto& root : _roots) {
      root.addCuboids(decompositions);
    }

    if (heuristic) { // remove morton order to mimic original decomp
      for (int i = 0; i < _maxRefinementLevel; i++) {
        std::shuffle(decompositions[i].begin(), decompositions[i].end(),
          std::default_random_engine(9));
      }
    }

#ifdef PARALLEL_MODE_MPI
    int localCount = 0;

    for (auto& decomp : decompositions) {
      localCount += decomp.size();
    }
    int size = singleton::mpi().getSize();

    std::vector<int> receiveCounts(size);
    std::vector<int> displs(size);
    singleton::mpi().allGather(&localCount, 1, receiveCounts.data(), 1);
    displs[0] = 0;
    for (int r = 1; r < size; ++r) {
      displs[r] = displs[r - 1] + receiveCounts[r - 1];
    }
    std::uint64_t total = displs.back() + receiveCounts.back();

    std::vector<T> xsSend(localCount);
    std::vector<T> xsRecv(total);
    std::vector<T> ysSend(localCount);
    std::vector<T> ysRecv(total);
    std::vector<T> zsSend;
    std::vector<T> zsRecv;

    if constexpr(D==3) {
      zsSend.resize(localCount);
      zsRecv.resize(total);
    }
    std::vector<T> borderSend(localCount);
    std::vector<T> borderRecv(total);
    std::vector<T> processSend(localCount);
    std::vector<T> processRecv(total);
    std::vector<T> extentSend(localCount);
    std::vector<T> extentRecv(total);
    std::vector<T> levelSend(localCount);
    std::vector<T> levelRecv(total);

    int i = 0;
    int j = 0;
    for (auto& decomp : decompositions) {
      for (auto& cuboid : decomp) {

        Vector<T, D> o           = cuboid.getOrigin();
        int          currExtent  = cuboid.getExtent()[0];
        T            currVoxel   = cuboid.getDeltaR();

        int          currProcess = cuboid.getPreCalculatedProcess().has_value()? cuboid.getPreCalculatedProcess().value() : singleton::mpi().getRank();
        xsSend[i]    = o[0];
        ysSend[i]    = o[1];
        if constexpr (D == 3) {
          zsSend[i] = o[2];
        }

        borderSend[i] = cuboid.isCoarseBorder();
        processSend[i] = currProcess;
        extentSend[i]  = currExtent;
        levelSend[i]   = j;
        i++;
      }
      j++;
    }

    singleton::mpi().allGatherv(xsSend.data(), localCount, xsRecv.data(), receiveCounts.data(), displs.data());
    singleton::mpi().allGatherv(ysSend.data(), localCount, ysRecv.data(), receiveCounts.data(), displs.data());
    if constexpr (D == 3) {
      singleton::mpi().allGatherv(zsSend.data(), localCount, zsRecv.data(), receiveCounts.data(), displs.data());
    }
    singleton::mpi().allGatherv(processSend.data(), localCount, processRecv.data(), receiveCounts.data(), displs.data());
    singleton::mpi().allGatherv(borderSend.data(), localCount, borderRecv.data(), receiveCounts.data(), displs.data());
    singleton::mpi().allGatherv(extentSend.data(), localCount, extentRecv.data(), receiveCounts.data(), displs.data());
    singleton::mpi().allGatherv(levelSend.data(), localCount, levelRecv.data(), receiveCounts.data(), displs.data());

    for (int i = 0; i < decompositions.size(); i++) {
      decompositions[i].clear();
    }

    for (int i = 0; i < total; ++i) {
      if constexpr (D == 3) {
        decompositions[levelRecv[i]].emplace_back(Cuboid<T,D>(Vector<T, D>(xsRecv[i], ysRecv[i], zsRecv[i]),
                                                  voxelSizes[levelRecv[i]], extentRecv[i],
                                                  std::make_optional(processRecv[i]),borderRecv[i]));
      }
      else {
        decompositions[levelRecv[i]].emplace_back(Cuboid<T,D>(Vector<T, D>(xsRecv[i], ysRecv[i]), voxelSizes[levelRecv[i]],
                                                  extentRecv[i], std::make_optional(processRecv[i]),borderRecv[i]));
      }
    }
#endif

    std::vector<CuboidDecomposition<T, D>> level;
    for (int i=0; i < decompositions.size(); i++) {
      level.emplace_back(CuboidDecomposition<T, D>(std::move(decompositions[i]), voxelSizes[i]));
    }

    return level;
  }

  std::vector<CuboidDecomposition<T, D>> createFineCuboidDecomposition()
  {
    std::vector<std::vector<Cuboid<T, D>>> decompositions;
    std::vector<T>                         voxelSizes;
    for (int i = 0; i <= _maxRefinementLevel; i++) {
      decompositions.push_back({});
      voxelSizes.push_back(voxelSizes.size() ? voxelSizes[i - 1] / 2 : _motherCuboid.getDeltaR());
    }
    for (auto& root : _roots) {
      root.addCuboidFine(decompositions);
    }

    std::vector<CuboidDecomposition<T, D>> level;
    for (int i = 0; i < decompositions.size(); i++) {
      level.push_back(CuboidDecomposition<T, D>(decompositions[i], voxelSizes[i]));
    }
    return level;
  }

  void loadBalance(int processes)
  {
    std::uint64_t localTotal = 0;
    std::uint64_t globalTotal = 0;
    std::uint64_t currentTotal = 0;
    for (auto& root : _roots) {
      localTotal += root.sumWeight();
    }

#ifdef PARALLEL_MODE_MPI
    int rank = singleton::mpi().getRank();
    int size = singleton::mpi().getSize();
    std::vector<std::uint64_t> allLocalTotals(size);

    singleton::mpi().allGather(
        &localTotal, 1,
        allLocalTotals.data(), 1
    );

    for (int i = 0; i < rank; ++i) {
      currentTotal += allLocalTotals[i];
    }

    for (int i = 0; i < size; ++i) {
      globalTotal += allLocalTotals[i];
    }
#else
    globalTotal = localTotal;
#endif

    std::uint64_t idealShare = globalTotal / processes;
    std::vector<std::uint64_t> load;

    for (auto& root : _roots) {
      currentTotal = root.assignProcess(currentTotal, idealShare, processes - 1, load);
    }

  }

  int getMaxRefinementLevel() { return _maxRefinementLevel; }

  void balance()
  {

    std::map<int, std::vector<Octant<T, D>>> comms;

    for (auto& root : _roots) {
      auto rComms = root.balanceLocal();
      for (auto& [id, octants] : rComms) {
        for (auto& octant : octants) {
          auto possiblePreclude =
              std::lower_bound(comms[id].begin(), comms[id].end(), octant, Octant<T, D>::precludeTotalOrder);
          if (possiblePreclude == comms[id].end() || !Octant<T, D>::precludeEqual(octant, *possiblePreclude)) {
            comms[id].insert(possiblePreclude, octant);
          }
        }
      }
    }
#ifdef PARALLEL_MODE_MPI

    std::vector<std::vector<std::uint64_t>> toSendID;
    std::vector<std::vector<int>>      toSendRef;

    std::vector<std::uint64_t> receivers;
    for (auto& [id, octants] : comms) {

      auto lb = std::lower_bound(neighborIDs.begin(), neighborIDs.end(), id);
      if (neighborRanks[std::distance(neighborIDs.begin(), lb)] != singleton::mpi().getRank()) {
        std::vector<std::uint64_t> vec;
        std::vector<int>      ref;
        for (auto x : octants) {
          vec.push_back(x.getId());
          ref.push_back(x.getRefinementLevel());
        }
        toSendID.push_back(vec);
        toSendRef.push_back(ref);
        receivers.push_back(id);
      }
    }

    auto [receivedIDs, _]          = root2Root(toSendID, receivers);
    auto [receivedRefinements, __] = root2Root(toSendRef, receivers);

    for (int i = 0; i < receivedIDs.size(); i++) {
      for (int j = 0; j < receivedIDs[i].size(); j++) {
        int currIndex = _roots[i].getID();
        Octant<T, D> octant =
                    Octant<T, D>(Octant<T, D>::calcOriginFast(receivedIDs[i][j]), receivedRefinements[i][j], _maxRefinementLevel);

        auto possiblePreclude = std::lower_bound(comms[currIndex].begin(), comms[currIndex].end(), octant, Octant<T, D>::precludeTotalOrder);
        if (possiblePreclude == comms[currIndex].end() || !Octant<T, D>::precludeEqual(octant, *possiblePreclude)) {
          comms[currIndex].insert(possiblePreclude, octant);
        }
      }
    }

#endif

#ifdef PARALLEL_MODE_OMP
#pragma omp parallel for
#endif
    for (std::size_t i = 0; i < _roots.size(); ++i) {
      _roots[i].balanceBorders(comms[_roots[i].getID()]);
    }
  }

  void coarse(bool preBalance = false)
  {
#ifdef PARALLEL_MODE_OMP
#pragma omp parallel for
#endif
    for (std::size_t i = 0; i < _roots.size(); ++i) {
      _roots[i].coarse(preBalance);
    }
  }

  void coarseBorders()
  {
    for (auto& root : _roots) {
      root.coarseBorders();
    }
  }

  void markBorders()
  {
    if (_roots.empty()) {
      return;
    }
    std::map<int, std::vector<Octant<T, D>>> comms;

#ifdef PARALLEL_MODE_MPI
    std::vector<std::vector<std::uint64_t>>       toSendID;
    std::vector<std::vector<int>>            toSendRef;
    std::vector<std::vector<int>>            toSendOriginalOctant;
    std::vector<std::uint64_t>                    receivers;
    std::vector<std::vector<bool>>           checked(neighborRanks.size(), std::vector<bool>(_roots.back().getID() - _roots.begin()->getID() +1, false));
#endif
    for (auto& root : _roots) {

      comms[root.getID()] = root.getBorderOctant();
#ifdef PARALLEL_MODE_MPI
      for (auto neighbor : root.getNeighbors()) {
        if (neighbor.has_value()) {
          auto lb = std::lower_bound(neighborIDs.begin(), neighborIDs.end(), neighbor.value());

          if (neighborRanks[std::distance(neighborIDs.begin(), lb)] != singleton::mpi().getRank() &&
              !checked[std::distance(neighborIDs.begin(), lb)][root.getID() - _roots[0].getID()]) {

            std::vector<std::uint64_t> vec;
            std::vector<int>      ref;
            std::vector<int>      oo;
            for (const auto& x : comms[root.getID()]) {
              vec.push_back(x.getId());
              ref.push_back(x.getRefinementLevel());
              oo.push_back(root.getID());
            }
            toSendID.push_back(vec);
            toSendRef.push_back(ref);
            toSendOriginalOctant.push_back(oo);
            receivers.push_back(neighbor.value());
            checked[std::distance(neighborIDs.begin(), lb)][root.getID() - _roots[0].getID()] = true;
          }
        }
      }
#endif
    }

#ifdef PARALLEL_MODE_MPI

    auto [receivedIDs, _]          = root2Root(toSendID, receivers);
    auto [receivedRefinements, __] = root2Root(toSendRef, receivers);
    auto [originalOctant, ___]     = root2Root(toSendOriginalOctant, receivers);

    for (int i = 0; i < originalOctant.size(); i++) {
      for (int j = 0; j < originalOctant[i].size(); j++) {
        comms[originalOctant[i][j]].push_back(
            Octant<T, D>(Octant<T, D>::calcOriginFast(receivedIDs[i][j]), receivedRefinements[i][j], _maxRefinementLevel));
      }
    }

#endif

    for (auto& root : _roots) {
      std::map<int, std::vector<Octant<T, D>>> local;
      auto                                     neighbors = root.getNeighbors();
      for (int i = 0; i < neighbors.size(); i++) {
        if (neighbors[i].has_value()) {
          local[i] = comms[neighbors[i].value()];
        }
        else {
          local[i] = {};
        }
      }
      root.markBorder(local);
    }
  }

  void writeUnmergedXML(const std::string& filename)
  {
    auto decompositions = createCuboidDecomposition(false);
    for (auto& decomp : decompositions) {
      if (decomp.size() == 0) continue;

      std::vector<Cuboid<T,D>> correctedCuboids;
      correctedCuboids.reserve(decomp.size());

      for (const auto& cuboid : decomp.cuboids()) {
        auto newCuboid = Cuboid<T,D>(
              cuboid.getOrigin(),
              cuboid.getDeltaR(),
              cuboid.getExtent(),
              0, // Force precalculated process to 0
              cuboid.isCoarseBorder()
          );
        newCuboid.setWeight(cuboid.getWeight());
        correctedCuboids.push_back(newCuboid);
      }
      decomp.cuboids().swap(correctedCuboids);
      decomp.updateMotherCuboid();
    }

    if (singleton::mpi().isMainProcessor()) {
      std::ofstream fout(filename);
      fout << std::setprecision(17);
      fout << "<?xml version=\"1.0\"?>\n<XMLContent>\n";
      for (size_t i = 0; i <= (size_t)_maxRefinementLevel; ++i) {
        fout << "  <Level" << i << ">\n";
        fout << "    <CuboidDecomposition dimension=\"" << D << "\" ";
        decompositions[i].getMotherCuboid().writeAsXML(fout);
        fout << ">\n";
        for (int iC = 0; iC < decompositions[i].size(); ++iC) {
          fout << "      <Cuboid ";
          // Ensure process is set to 0 for initial export to satisfy PredefinedLoadBalancer
          if (!decompositions[i].get(iC).getPreCalculatedProcess().has_value()) {
            decompositions[i].get(iC).setPreCalculatedProcess(0);
          }
          decompositions[i].get(iC).writeAsXML(fout);
          fout << " />\n";
        }
        fout << "    </CuboidDecomposition>\n";
        fout << "  </Level" << i << ">\n";
      }
      fout << "</XMLContent>\n";
    }
  }

  void compare(OctreeRefinedMeshBlueprint& other)
  {
    auto first  = _roots.begin();
    auto second = other._roots.begin();
    while (first != _roots.end()) {
      first->compare(*second);
      ++first;
      ++second;
    }
  }

#ifdef PARALLEL_MODE_MPI
  /// Takes vector with custom length of vectors with messages to root id at same index in receivers
  /// Returns vectors with lengths of _root, first containing linearized received responses, second with cutoffs to delinearize
  template <typename U>
  std::pair<std::vector<std::vector<U>>, std::vector<std::vector<std::uint64_t>>>
  root2Root(std::vector<std::vector<U>>& toSend, std::vector<std::uint64_t>& receivers)
  {
    if (_roots.empty()) {
      return {};
    }
    std::vector<std::vector<U>>        linearizedToSend(neighborRanks.size());
    std::vector<std::vector<std::uint64_t>> cutsToSend(neighborRanks.size());
    std::vector<std::vector<std::uint64_t>> idsToSend(neighborRanks.size());

    for (int i = 0; i < receivers.size(); i++) {
      auto x     = std::lower_bound(neighborIDs.begin(), neighborIDs.end(), receivers[i]);
      int  index = std::distance(neighborIDs.begin(), x);

      linearizedToSend[index].insert(linearizedToSend[index].end(), toSend[i].begin(), toSend[i].end());
      cutsToSend[index].push_back(toSend[i].size());
      idsToSend[index].push_back(receivers[i]);
    }

    std::vector<std::vector<U>> linearizedReceived = rank2Rank(std::move(linearizedToSend));
    std::vector<std::vector<std::uint64_t>> cutsReceived = rank2Rank(std::move(cutsToSend));
    std::vector<std::vector<std::uint64_t>> idsReceived = rank2Rank(std::move(idsToSend));

    std::vector<int> lookup(_roots.back().getID() - _roots.front().getID() + 1);
    std::uint64_t lowest = _roots.front().getID();
    for (int i = 0; i < _roots.size(); i++) {
      lookup[_roots[i].getID() - lowest] = i;
    }

    for (int i = 0; i < idsReceived.size(); i++) {
      for (int j = 0; j < idsReceived[i].size(); j++) {
        idsReceived[i][j] -= lowest;
      }
    }

    std::vector<std::vector<U>> rootWiseLinearized(_roots.size());
    std::vector<std::vector<std::uint64_t>> rootWiseCuts(_roots.size());

    for (int i = 0; i < neighborRanks.size(); i++) {

      auto head = linearizedReceived[i].begin();
      auto tail = linearizedReceived[i].begin();
      for (int j = 0; j < idsReceived[i].size(); j++) {

        head += cutsReceived[i][j];
        rootWiseCuts[lookup[idsReceived[i][j]]].push_back(cutsReceived[i][j]);
        rootWiseLinearized[lookup[idsReceived[i][j]]].insert(rootWiseLinearized[lookup[idsReceived[i][j]]].end(), tail, head);
        tail = head;
      }
    }

    return std::make_pair(rootWiseLinearized, rootWiseCuts);
  }

  /// takes vector of length neighborRanks and sends inner vectors to rank index
  /// returns vector of lenght neighborRanks which got send my rank index
  template <typename U>
  std::vector<std::vector<U>> rank2Rank(std::vector<std::vector<U>> toSend)
  {
    if (_roots.empty()) {
      return toSend;
    }
    if (toSend.size() != neighborRanks.size()) [[unlikely]] {
      throw std::logic_error("invalid Send Operation");
    }

    std::vector<std::vector<U>> received(neighborRanks.size());
    std::vector<std::uint64_t>       sizesToSend(neighborRanks.size());
    std::vector<std::uint64_t>       sizesToReceive(neighborRanks.size());
    _mpiNbHelper.free();
    _mpiNbHelper.allocate(neighborRanks.size() * 2 - 2);
    int req = 0;

    for (int i = 0; i < toSend.size(); i++) {
      sizesToSend[i] = toSend[i].size();
    }

    for (int i = 0; i < neighborRanks.size(); i++) {
      if (neighborRanks[i] == singleton::mpi().getRank()) {
        continue;
      }
      singleton::mpi().iSend(&sizesToSend[i], 1, neighborRanks[i], _mpiNbHelper.get_mpiRequest(req++), 0);
    }

    for (int i = 0; i < neighborRanks.size(); i++) {
      if (neighborRanks[i] == singleton::mpi().getRank()) {
        continue;
      }

      singleton::mpi().iRecv(&sizesToReceive[i], 1, neighborRanks[i], _mpiNbHelper.get_mpiRequest(req++), 0);
    }

    singleton::mpi().waitAll(_mpiNbHelper);

    _mpiNbHelper.free();
    _mpiNbHelper.allocate(neighborRanks.size() * 2 - 2);
    req = 0;

    for (int i = 0; i < neighborRanks.size(); i++) {
      if (neighborRanks[i] == singleton::mpi().getRank()) {
        received[i] = std::move(toSend[i]);
        continue;
      }
      received[i].resize(sizesToReceive[i]);
    }

    for (int i = 0; i < neighborRanks.size(); i++) {
      if (neighborRanks[i] == singleton::mpi().getRank()) {
        continue;
      }
      U* buf = toSend[i].empty() ? nullptr : toSend[i].data();
      singleton::mpi().iSend(buf, toSend[i].size(), neighborRanks[i], _mpiNbHelper.get_mpiRequest(req++), 0);
    }

    for (int i = 0; i < neighborRanks.size(); i++) {
      if (neighborRanks[i] == singleton::mpi().getRank()) {
        continue;
      }
      singleton::mpi().iRecv(received[i].data(), sizesToReceive[i], neighborRanks[i],
                             _mpiNbHelper.get_mpiRequest(req++), 0);
    }
    singleton::mpi().waitAll(_mpiNbHelper);

    return received;
  }
#endif
};

} // namespace olb

#endif // OLB_REFINEMENT_OCTREE_REFINEDMESHBLUEPRINT_H
