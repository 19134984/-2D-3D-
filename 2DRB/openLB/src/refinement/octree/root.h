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

#ifndef OLB_REFINEMENT_OCTREE_ROOT_H
#define OLB_REFINEMENT_OCTREE_ROOT_H

#include <fstream>
#include <list>
#include <vector>
#include <optional>
#include <stack>
#include <string>
#include <map>
#include <algorithm>

#include "octant.h"
#include "core/vector.h"
#include "core/util.h"

namespace olb {

template <typename T, unsigned D>
class OctreeRoot {
private:
  /// ids of neighbors octree roots in order
  std::list<Octant<T, D>> _octants;

  /// Origin of octant relative to octree root coordinate system
  Vector<T, D> _origin;

  /// size 27 in D = 3, 9 in D = 2;
  std::vector<std::optional<int>> _neighborID;

  T _voxelSize;

  int _extent;

  int _maxRefinement;

  int _globalLogicalRefinement;

  int _id;

public:
  OctreeRoot(Vector<T, D> origin, T voxelSize, int extent, int maxRefinementLevel, int id, std::uint64_t globalLogicalRefinement, Vector<int,D> domainExtent, bool hilbert = false)
  : _neighborID(util::pow(3, D))
  {
    _origin = origin;
    Vector<T,D> octantOrigin;
    for (int i = 0; i < (int)D; i++) {
      octantOrigin[i] = 0;
    }
    _octants.push_back(Octant<T, D>(octantOrigin, 0, maxRefinementLevel));
    _voxelSize     = voxelSize;
    _extent        = extent;
    _maxRefinement = maxRefinementLevel;
    _id            = id;
    _globalLogicalRefinement = globalLogicalRefinement;

    Vector<int,D> MortonOrigin;
    if (hilbert) {
      MortonOrigin = Octant<T,D>::calcOriginFast(Octant<T,D>::calcMortonID(_id, _globalLogicalRefinement));
    }
    else {
      MortonOrigin = Octant<T,D>::calcOriginFast(id);
    }

    int lowx = MortonOrigin[0] == 0 ? 0: -1;
    int lowy = MortonOrigin[1] == 0 ? 0: -1;
    int lowz = 0;
    if constexpr (D == 3) {
      lowz = MortonOrigin[2] == 0 ? 0: -1;
    }
    Vector<int,D> n;
    for (int x = lowx; x <= 1; x++) {
      for (int y = lowy; y <= 1; y++) {
        if constexpr (D == 2) {
          if (hilbert) {
            n = Octant<T,D>::calcOriginFast(Octant<T,D>::calcMortonID(_id, _globalLogicalRefinement)) + Vector<int,D>(x, y);
          }
          else {
            n = Octant<T,D>::calcOriginFast(_id) + Vector<int,D>(x, y);
          }
          if (n[0] < 0 || n[0] >= domainExtent[0] || n[1] < 0 || n[1] >= domainExtent[1])  {
            _neighborID[x + 1 + (y + 1 ) * 3] = std::nullopt;
          }
          else if (hilbert) {
            _neighborID[x + 1 + (y + 1 ) * 3] = std::make_optional((int)Octant<T,D>::calcHilbertID(n, _globalLogicalRefinement));
          }
          else {
            _neighborID[x + 1 + (y + 1 ) * 3] = std::make_optional((int)Octant<T,D>::calcIDFast(n));
          }
          }
        else {
          for (int z = lowz; z <= 1; z++) {
            if (hilbert) {
              n = Octant<T,D>::calcOriginFast(Octant<T,D>::calcMortonID(_id, _globalLogicalRefinement)) + Vector<int,D>(x, y, z);
            }
            else {
              n = Octant<T,D>::calcOriginFast(_id) + Vector<int,D>(x, y, z);
            }
            if (n[0] < 0 || n[0] >= domainExtent[0] || n[1] < 0 || n[1] >= domainExtent[1] || n[2] < 0 || n[2] >= domainExtent[2])  {
              _neighborID[x + 1 + (y + 1 ) * 3 + (z+1) * 9] = std::nullopt;
            }
            else if (hilbert) {
                _neighborID[x + 1 + (y + 1 ) * 3 + (z+1) * 9] = std::make_optional((int)Octant<T,D>::calcHilbertID(n, _globalLogicalRefinement));
            }
            else {
              _neighborID[x + 1 + (y + 1 ) * 3 + (z+1) * 9] = std::make_optional((int)Octant<T,D>::calcIDFast(n));
              }
            }
          }
        }
      }
    }

  void addCuboids(std::vector<std::vector<Cuboid<T, D>>>& result)
  {
    for (Octant<T, D>& octant : _octants) {
      octant.addCuboid(result, _origin, _voxelSize, _extent);
    }
  }

  void addCuboidFine(std::vector<std::vector<Cuboid<T, D>>>& result)
  {
    for (Octant<T, D> octant : _octants) {
      octant.addCuboidFine(result, _origin, _voxelSize, _extent);
    }
  }

  bool intersect(IndicatorF<T, D>& indicatorF)
  {
    Vector<T, D> rootMax = _origin + (T)_extent * _voxelSize;
    Vector<T, D> iMin = indicatorF.getMin();
    Vector<T, D> iMax = indicatorF.getMax();
    for (unsigned iD=0; iD < D; ++iD) {
      if (_origin[iD] > iMax[iD] || rootMax[iD] < iMin[iD]) {
        return false;
      }
    }
    return true;
  }

  std::list<Octant<T,D>>& getOctants()
  {
    return _octants;
  }

  void coarse(bool preBalance)
  {
    std::stack<typename std::list<Octant<T, D>>::iterator> its;
    std::stack<typename std::list<Octant<T, D>>::iterator> lookaheads;
    std::stack<int> aheads;
    std::stack<bool> borders;

    auto it = _octants.begin();
    auto lookahead = _octants.begin();

    int ahead = 0;
    int amountOfChilds = util::pow(2,D);
    bool border = false;

    while (it != _octants.end()) {
      if (ahead == amountOfChilds) {
        it->incDomainLevel();
        if (border) {
          it->setRefinementBorder();
        }
        _octants.erase(std::next(it), lookahead);
        if (its.empty()) {
          lookahead = it;
          ahead = 0;
          border = false;
        }
        else {
          it = its.top();
          its.pop();
          lookahead = lookaheads.top();
          lookaheads.pop();
          ahead = aheads.top();
          aheads.pop();
          border = borders.top();
          borders.pop();
        }
      }
      else {
        if (lookahead != _octants.end() && (preBalance ||  it->getProcess().value() == lookahead->getProcess().value())) {
          if (it->getRefinementLevel() == lookahead->getRefinementLevel() && it->getDomainLevel() == lookahead->getDomainLevel() && it->getParentOriginDomainLevel() == lookahead->getParentOriginDomainLevel()) {
            border |= lookahead->isCoarseBorderOctant();
            ++lookahead;
            ahead++;
          }
          else {
            its.push(it);
            lookaheads.push(lookahead);
            borders.push(border);
            aheads.push(ahead);
            it = lookahead;
            ahead = 0;
            border = false;
          }
        }
        else {
          its = std::stack<typename std::list<Octant<T, D>>::iterator>();
          lookaheads = std::stack<typename std::list<Octant<T, D>>::iterator>();
          aheads = std::stack<int>();
          borders = std::stack<bool>();

          it = lookahead;
          ahead = 0;
          border = false;

        }
      }
    }
  }

  void coarseBorders()
  {
    std::stack<typename std::list<Octant<T, D>>::iterator> its;
    std::stack<typename std::list<Octant<T, D>>::iterator> lookaheads;
    std::stack<int> aheads;

    auto it = _octants.begin();
    auto lookahead = _octants.begin();

    int ahead = 0;
    while (it != _octants.end()) {
      if (ahead == 4) {
        it->incDomainLevel();
        _octants.erase(std::next(it), lookahead);
        if (its.empty()) {
          lookahead = it;
          ahead = 0;
        }
        else {
          it = its.top();
          its.pop();
          lookahead = lookaheads.top();
          lookaheads.pop();
          ahead = aheads.top();
          aheads.pop();
        }
      }
      else {
        if (it->getProcess().value() == lookahead->getProcess().value()) {
          if (it->getOrigin()[2] == lookahead->getOrigin()[2]) {
            if (it->isCoarseBorderOctant() && lookahead->isCoarseBorderOctant()) {
              if (it->getRefinementLevel() == lookahead->getRefinementLevel() && it->getParentOrigin() == lookahead->getParentOrigin()) {
                ahead++;
                ++lookahead;
              }
              else {
                it = lookahead;
                ahead = 0;
              }
            }
            else {
              if (its.empty()) {
                lookahead = it;
                ahead = 0;
              }
              else {
                it = its.top();
                its.pop();
                lookaheads.pop();
                ++lookahead;
                ahead = aheads.top();
                aheads.pop();
              }
            }
          }
          else {
            if (lookahead->isCoarseBorderOctant()) {
              its.push(it);
              lookaheads.push(lookahead);
              aheads.push(ahead);
              it = lookahead;
              ahead = 0;
            }
            else {
              ++lookahead;
            }
          }
        }
        else {
          its = std::stack<typename std::list<Octant<T, D>>::iterator>();
          lookaheads = std::stack<typename std::list<Octant<T, D>>::iterator>();
          aheads = std::stack<int>();

          it = lookahead;
          ahead = 0;
        }
      }
    }
  }

  void reduceWeights(IndicatorF<T, D>& indicatorF)
  {
    for (auto& it: _octants) {
      it.reduceWeights(indicatorF, _origin, _voxelSize, _extent);
    }
  }

  void refine(IndicatorF<T, D>& indicatorF, int refinementLevel, int padding)
  {
    for (auto it = _octants.begin(); it != _octants.end();) {
      if (it->getRefinementLevel() < refinementLevel && it->intersect(indicatorF, _origin, _voxelSize, _extent, padding)) {
        std::list<Octant<T, D>> children = it->refine();

        auto insertPos  = _octants.erase(it);
        auto childrenIT = children.begin();

        _octants.splice(insertPos, children);
        it = childrenIT;
      }
      else {
        ++it;
      }
    }
  }

  std::uint64_t number(std::uint64_t current)
  {
    for (auto& octant: _octants) {
      octant.setID(current);
      current++;
    }
    return current;
  }

  void writePartition(std::string fileName)
  {
    std::ofstream file(fileName, std::ios::app);
    for (auto& octant: _octants) {
      file << octant.getProcess().value() << std::endl;
    }
    file.close();
  }

  void balanceGraph(std::map<int, std::vector<Octant<T, D>>>& neighbors, std::string fileName, bool discount)
  {
    std::vector<std::pair<Octant<T,D>, std::uint64_t>> callback;
    std::uint64_t first = _octants.front().getGlobalID();
    std::uint64_t last = _octants.back().getGlobalID();
    // (Edge Weight, Node ID)
    std::vector<std::vector<std::pair<std::uint64_t, std::uint64_t>>> buff(last - first + 1);
    auto buffIt = buff.begin();
    std::vector<std::uint64_t> weights;
    weights.reserve(last - first + 1);
    int center = 0;
    for (int i = 0; i < (int)D; i++) {
      center += util::pow(3,i);
    }

    for (auto& octant: _octants) {
      weights.push_back(octant.calcWeight());
      for (auto [neighbor, neighborID]: octant.getAllTouching()) {
        if (neighborID != center) {
          auto pos = std::lower_bound(neighbors[neighborID].begin(),neighbors[neighborID].end(), neighbor);
          if (pos != neighbors[neighborID].end() && *pos == neighbor) {
            buffIt->push_back(std::make_pair(octant.calcEdgeWeight(*pos, discount, true), pos->getGlobalID()));
          }
        }
        else {
          callback.emplace_back(std::make_pair(neighbor,octant.getGlobalID()));
        }
      }
      ++buffIt;
    }

    std::sort(callback.begin(), callback.end());
    callback.erase( std::unique( callback.begin(), callback.end() ), callback.end());
    auto it = _octants.begin();
    auto it2 = callback.begin();
    while (it != _octants.end() && it2 != callback.end()) {
      if (*it == it2->first) {
        buff[it2->second - first].push_back(std::make_pair(it->calcEdgeWeight(it2->first, discount),it->getGlobalID()));
        ++it2;
      }
      else if (it2->first.getId() <= it->getId()) {
        ++it2;
      }
      else {
        ++it;
      }
    }

    std::ofstream file(fileName, std::ios::app);

    int i = 0;
      for (auto it_vec: buff) {
          file << weights[i];
          for (auto pair: it_vec) {
            file << " "<< (pair.second + 1) << " " << pair.first;
          }
          file << std::endl;
        i++;
      }
    file.close();
  }

  void refine(std::vector<Octant<T, D>>& required)
  {
    auto                  itMarked = required.begin();
    auto                  it       = _octants.begin();
    std::uint64_t              mask     = 0;
    std::vector<std::uint64_t> masks    = {0};

    for (int i = _maxRefinement; i > 0; i--) {

      mask |= ~(~0ULL << D) ;
      masks.push_back(mask);
      mask <<= D;
    }

    while (itMarked != required.end()) {
      if (it == _octants.end()) {
        throw std::logic_error("somethings wrong");
      }
      if (itMarked->getId() == it->getId() && it->getRefinementLevel() >= itMarked->getRefinementLevel()) {
        ++itMarked;
        continue;
      }
      if (it->getRefinementLevel() >= itMarked->getRefinementLevel()) {
        ++it;
        continue;
      }

      if (it->getId() == (itMarked->getId() & ~masks[_maxRefinement - it->getRefinementLevel()])) {
        std::list<Octant<T, D>> children = it->refine();

        auto insertPos  = _octants.erase(it);
        auto childrenIT = children.begin();

        _octants.splice(insertPos, children);
        it = childrenIT;
      }
      else {
        ++it;
      }
    }
  }
  std::vector<std::optional<int>> getNeighbors()
  {
    return _neighborID;
  }

  void removeNeighbors(std::vector<std::uint64_t>& neighbors)
  {
    for (auto& neighbor: _neighborID) {
      if (neighbor.has_value() && std::binary_search(neighbors.begin(), neighbors.end(), (std::uint64_t)neighbor.value())) {
        neighbor = std::nullopt;
      }
    }
  }

  void countCuboids(std::vector<int>& count)
  {
    for (auto octant : _octants) {
      count[octant.getRefinementLevel()] += 1;
    }
  }

  int getID() { return _id; }

  std::uint64_t sumWeight()
  {
    std::uint64_t total = 0;
    for (auto octant : _octants) {
      total += octant.calcWeight();
    }
    return total;
  }

  void assignProcess(std::ifstream& partition)
  {
    std::string part;
    for (auto& octant: _octants) {
      std::getline(partition, part);
      octant.setProcess(std::stoi(part));
    }
  }

  int assignProcess(std::uint64_t processed, std::uint64_t idealShare, std::uint64_t max, std::vector<std::uint64_t>& load)
  {
    int inc = processed / idealShare;
    for (auto& octant : _octants) {

      octant.setProcess(inc > (int)max ? (int)max : inc);
      processed += octant.calcWeight();
      if (processed / idealShare > (std::uint64_t)inc) {
        load.push_back(processed);
        inc++;
      }
    }
    return processed;
  }

  std::vector<Octant<T, D>> Reduce()
  {
    std::vector<Octant<T, D>> reduced;
    reduced.push_back(*_octants.begin());
    int i = 0;
    for (auto it = _octants.begin(); it != _octants.end(); ++it) {
      Octant<T, D> sibling = Octant<T, D>(it->getParentOrigin(), it->getRefinementLevel(), _maxRefinement);
      if (Octant<T, D>::precludeTrue(reduced[i], sibling)) {
        reduced[i] = sibling;
      }
      else if (!Octant<T, D>::precludeEqual(sibling, reduced[i])) {
        reduced.push_back(sibling);
        i++;
      }
    }
    return reduced;
  }

  std::map<int, std::vector<Octant<T,D>>> balanceLocal()
  {
    std::vector<Octant<T, D>>                reduced = Reduce();
    std::map<int, std::vector<std::pair<Octant<T, D>, bool>>> workingSet;

    for (auto o : reduced) {
      workingSet[o.getRefinementLevel()].push_back(std::make_pair(o, true));
    }
    int max_ref = 0;
    for (auto it = workingSet.begin(); it != workingSet.end(); ++it) {
      max_ref = std::max(it->first, max_ref);
    }
    std::vector<Octant<T, D>> all;
    for (auto [a,b]: workingSet[max_ref]) {
      all.push_back(a);
    }

    std::vector<std::vector<Octant<T,D>>> neighborRoots(util::pow(3,D));
    int center = 0;
    for (int i = 0; i < (int)D; i++) {
      center += util::pow(3,i);
    }
    int i = max_ref;
    while (i > 1) {

      for (auto [octant, legacy] : workingSet[i]) {
        if (legacy && i != max_ref) {
          auto possiblePreclude = std::lower_bound(all.begin(), all.end(), octant, Octant<T, D>::precludeTotalOrder);
          if (possiblePreclude != all.end() && Octant<T,D>::precludeEqual(octant, *possiblePreclude)) {
            continue;
          }
        }
        for (auto [neighbor, neighborID]: octant.getCoarseGlobalNeighbors()) {

          Octant<T, D> sibling =
              Octant<T, D>(neighbor.getParentOrigin(), neighbor.getRefinementLevel(), _maxRefinement);

          if (neighborID != (center)) {
            auto possiblePreclude = std::lower_bound(neighborRoots[neighborID].begin(), neighborRoots[neighborID].end(), sibling, Octant<T, D>::precludeTotalOrder);
            if(possiblePreclude == neighborRoots[neighborID].end() || !Octant<T,D>::precludeEqual(sibling, *possiblePreclude)) {
              neighborRoots[neighborID].insert(possiblePreclude, sibling);
            }
            continue;
          }

          if (sibling.getId() != octant.getId() && std::binary_search(all.begin(), all.end(), sibling, Octant<T,D>::idComp)) {
            continue;
          }
          workingSet[sibling.getRefinementLevel()].push_back(std::make_pair(sibling, false));
          auto possiblePreclude = std::lower_bound(all.begin(), all.end(), sibling, Octant<T, D>::precludeTotalOrder);
          if (possiblePreclude == all.end() || !Octant<T,D>::precludeEqual(sibling, *possiblePreclude)) {
            all.insert(possiblePreclude, sibling);
          }
        }
      }
      i--;
    }

    _octants = {};

    for (auto x: all) {
      _octants.push_back(x);
    }

    std::map<int, std::vector<Octant<T,D>>> comms;
    for (int idx = 0; idx < (int)neighborRoots.size(); idx++) {
      if (_neighborID[idx].has_value()) {
        comms[_neighborID[idx].value()] = neighborRoots[idx];
      }
    }
    return comms;
  }

  ///borders may have duplicates
  void balanceBorders(std::vector<Octant<T,D>> borders) {
    std::map<int, std::vector<Octant<T,D>>> workingSet;
    std::map<int, std::vector<Octant<T,D>>> workingSetBorder;

    for (auto o : _octants) {
      workingSet[o.getRefinementLevel()].push_back(o);
    }
    for (auto o : borders) {
      workingSetBorder[o.getRefinementLevel()].push_back(o);
    }

    int max_ref = 0;
    for (auto it = workingSet.begin(); it != workingSet.end(); ++it) {
      max_ref = std::max(it->first, max_ref);
    }
    for (auto it = workingSetBorder.begin(); it != workingSetBorder.end(); ++it) {
      max_ref = std::max(it->first, max_ref);
    }

    std::vector<Octant<T, D>> all;
    for (auto a: workingSet[max_ref]) {
      all.push_back(a);
    }
    for (auto a: workingSetBorder[max_ref]) {
      all.push_back(a);
    }

    int i = max_ref;
    while (i > 1) {

      for (auto octant : workingSet[i]) {
        auto possiblePreclude = std::lower_bound(all.begin(), all.end(), octant, Octant<T, D>::precludeTotalOrder);
        if (possiblePreclude == all.end() || !Octant<T,D>::precludeEqual(octant, *possiblePreclude)) {
          all.insert(possiblePreclude, octant);
        }
      }

      for (auto octant : workingSetBorder[i]) {
        auto possiblePreclude = std::lower_bound(all.begin(), all.end(), octant, Octant<T, D>::precludeTotalOrder);
        if (possiblePreclude != all.end() && Octant<T,D>::precludeEqual(octant, *possiblePreclude)) {
          continue;
        }
        all.insert(possiblePreclude, octant);
        for (auto neighbor: octant.getCoarseLocalNeighbors()) {

          Octant<T, D> sibling = Octant<T, D>(neighbor.getParentOrigin(), neighbor.getRefinementLevel(), _maxRefinement);

          if (std::binary_search(all.begin(), all.end(), sibling, Octant<T,D>::idComp)) {
            continue;
          }
          workingSetBorder[sibling.getRefinementLevel()].push_back(sibling);
          auto possiblePreclude = std::lower_bound(all.begin(), all.end(), sibling, Octant<T, D>::precludeTotalOrder);
          if (possiblePreclude == all.end() || !Octant<T,D>::precludeEqual(sibling, *possiblePreclude)) {
            all.insert(possiblePreclude, sibling);
          }
        }
      }
      i--;
    }
    complete(all);
  }

  std::vector<Octant<T,D>> getBorderOctant()
  {
    std::vector<Octant<T,D>> borders;
    for (auto o: _octants) {
      if (o.isBorder()) {
        borders.push_back(o);
      }
    }
    return borders;
  }

  void markBorder(std::map<int,std::vector<Octant<T,D>>>& neighbors)
  {
    if (_octants.size() == 1) {
      return;
    }
    std::vector<std::pair<Octant<T,D>, std::uint64_t>> callback;
    int center = 0;
    for (int i = 0; i < (int)D; i++) {
      center += util::pow(3,i);
    }
    for (auto& octant: _octants) {

      for (auto [neighbor, neighborID]: octant.getCoarseGlobalTouching()) {
        if (neighborID != center) {
          if (std::binary_search(neighbors[neighborID].begin(),neighbors[neighborID].end(), neighbor)) {
            octant.setRefinementBorder();
            break;
          }
        }
        else {
          callback.emplace_back(std::make_pair(neighbor,octant.getId()));
        }
      }
    }

    std::sort(callback.begin(), callback.end());
    callback.erase( std::unique( callback.begin(), callback.end() ), callback.end() );
    std::vector<std::uint64_t> called;
    auto it = _octants.begin();
    auto it2 = callback.begin();
    while (it != _octants.end() && it2 != callback.end()) {
      if (*it == it2->first) {
        called.push_back(it2->second);
        ++it2;
      }
      else if (it2->first.getId() <= it->getId()) {
        ++it2;
      }
      else {
        ++it;
      }
    }

    it = _octants.begin();
    std::sort(called.begin(), called.end());
    called.erase( std::unique( called.begin(), called.end() ), called.end() );
    auto it3 = called.begin();
    while (it != _octants.end() && it3 != called.end()) {
      if (it->getId() == *it3) {
        it->setRefinementBorder();
        ++it3;
        ++it;
      }
      else if (it->getId() < *it3) {
        ++it;
      }
      else {
        ++it3;
      }
    }
  }


  void compare(OctreeRoot<T,D>& other)
  {
    auto j = _octants.begin();
    auto i = other._octants.begin();

    while (j != _octants.end() && i != other._octants.end()) {
      if (j->getId() == i->getId()) {
        ++j;
        ++i;
      }
      else if (j->getId() < i->getId()) {
        j->setProcess(20);
        ++j;
      }
      else if (j->getId() > i->getId()) {
        i->setProcess(20);
        ++i;
      }
    }

    while (j != _octants.end()) {
      j->setProcess(20);
      ++j;
    }
    while (i != other._octants.end()) {
      i->setProcess(20);
      ++i;
    }
  }



  void complete(std::vector<Octant<T, D>>& reduced)
  {
    _octants = {};
    if constexpr (D == 2) {
      _octants.emplace_back(Octant<T, D>({0, 0}, 0, _maxRefinement));
    }
    else {
      _octants.emplace_back(Octant<T, D>({0, 0, 0}, 0, _maxRefinement));
    }
    std::sort(reduced.begin(), reduced.end());

    refine(reduced);
  }

  void transformHilbert()
  {
    int o = Octant<T,D>::calcHilbertOrientation(Octant<T,D>::calcMortonID(_id, _maxRefinement), _globalLogicalRefinement);
    _octants.sort([o](const Octant<T,D>& a, const Octant<T,D>& b) {
    return a.calcHilbertID(o) < b.calcHilbertID(o);
    });
  }




  void transformMorton()
  {
    _octants.sort();
  }

  static bool hilbertLEQ(OctreeRoot<T, D>& first, OctreeRoot<T, D>& second) {
    return Octant<T,D>::calcHilbertID(first._id, first._globalLogicalRefinement) < Octant<T,D>::calcHilbertID(second._id, second._globalLogicalRefinement);
  }

  bool operator<(const OctreeRoot<T, D>& rhs) const
  {
    return _id < rhs._id;
  }
};

} // namespace olb

#endif // OLB_REFINEMENT_OCTREE_ROOT_H
