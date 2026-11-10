/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Adrian Kummerlaender
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

#ifndef REFINEMENT_LOAD_BALANCER_H
#define REFINEMENT_LOAD_BALANCER_H

#include "simpleLoadBalancer.h"

#include <optional>

#include "communication/loadBalancer.h"
#include "geometry/cuboidDecomposition.h"

namespace olb {

/// Load balancer for refined lattice hierarchies
/**
 * Ensures that fine cuboids are assigned to the same
 * rank and platform as their corresponding coarse parents.
 *
 * Coarse parent load balancer must be aware of additional
 * weight due to refinement for a good balancing.
 **/
template<typename T, unsigned D>
class RefinedLoadBalancer final : public LoadBalancer<T> {
private:
  std::map<int,int> _coarseLoc;

public:
  RefinedLoadBalancer(CuboidDecomposition<T,D>& coarseDecomp,
                      LoadBalancer<T>& coarseBalancer,
                      CuboidDecomposition<T,D>& fineDecomp, bool presetRanks = false)
    : LoadBalancer<T>(0)
  {
    //OstreamManager clout(std::cout, "RefinedLoadBalancer");
    //clout.setMultiOutput(true);
    std::vector<int> rankBuffer(fineDecomp.size(), 0);
    std::vector<int> locBuffer(fineDecomp.size(), 0);
    std::set<int> freeCs;
    std::map<int,int> nLoc;
    for (int iC=0; iC < fineDecomp.size(); ++iC) {
      const auto& fineCuboid = fineDecomp.get(iC);
      if (auto latticeR = coarseDecomp.getLatticeR(fineCuboid.getOrigin()+0.5*fineCuboid.getDeltaR()*fineCuboid.getExtent())) {
        const int coarseC = (*latticeR)[0];
        if (coarseBalancer.isLocal(coarseC)) {
          rankBuffer[iC] = coarseBalancer.rank(coarseC);
          locBuffer[iC] = nLoc[rankBuffer[iC]]++;
          this->_platform[iC] = coarseBalancer.platform(coarseBalancer.loc(coarseC));
          //clout << "Matching fineC=" << iC << " to coarseC=" << coarseC << " on rank=" << rankBuffer[iC] << std::endl;
          _coarseLoc[locBuffer[iC]] = coarseBalancer.loc(coarseC);
        }
      } else {
        freeCs.emplace(iC);
      }
    }

    if (!freeCs.empty()) {
      //clout << "Balancing n=" << freeCs.size() << " free fine cuboids" << std::endl;
      // Build temporary decomposition containg only free cuboids
      CuboidDecomposition<T,D> freeDecomp = fineDecomp;
      for (int iC=freeDecomp.size()-1; iC >= 0; --iC) {
        if (freeCs.find(iC) == freeCs.end()) {
          freeDecomp.remove(iC);
        }
      }
      // Load balance free cuboids, assume that coarse is well balanced
      std::unique_ptr<LoadBalancer<T>> freeBalancer;

      if (presetRanks) {
        freeBalancer = std::make_unique<PredefinedLoadBalancer<T>>(freeDecomp);
      } else {
        freeBalancer = std::make_unique<HeuristicLoadBalancer<T>>(freeDecomp);
      }

      for (int iC=0; iC < freeDecomp.size(); ++iC) {
        if (freeBalancer->isLocal(iC)) {
          const auto& freeCuboid = freeDecomp.get(iC);
          if (auto latticeR = fineDecomp.getLatticeR(freeCuboid.getOrigin())) {
            const int fineC = (*latticeR)[0];
            if (freeCs.find(fineC) == freeCs.end()) {
              throw std::logic_error("Temporary decomposition inconsistent");
            } else {
              rankBuffer[fineC] = freeBalancer->rank(iC);
              locBuffer[fineC] = nLoc[rankBuffer[fineC]]++;
              this->_platform[fineC] = freeBalancer->platform(iC);
              //clout << "Assigning free fineC=" << fineC << " to rank=" << rankBuffer[fineC] << std::endl;
            }
          }
        }
      }
    }

    #ifdef PARALLEL_MODE_MPI
    std::vector<int> globalRankBuffer(fineDecomp.size(), 0);
    std::vector<int> globalLocBuffer(fineDecomp.size(), 0);
    // Broadcast assignments to all processes
    singleton::mpi().allreduce(rankBuffer.data(), globalRankBuffer.data(), rankBuffer.size(), MPI_SUM);
    singleton::mpi().allreduce(locBuffer.data(), globalLocBuffer.data(), locBuffer.size(), MPI_SUM);
    #else
    auto& globalRankBuffer = rankBuffer;
    auto& globalLocBuffer = locBuffer;
    #endif

    // Update internal LoadBalancer structure to match given assignment
    for (int iC=0; iC < fineDecomp.size(); ++iC) {
      this->_rank[iC] = globalRankBuffer[iC];
      this->_loc[iC] = globalLocBuffer[iC];
      if (globalRankBuffer[iC] == singleton::mpi().getRank()) {
        this->_glob.resize(std::max(static_cast<int>(this->_glob.size()), this->_loc[iC]+1));
        this->_glob[this->_loc[iC]] = iC;
        this->_size = this->_glob.size();
      }
    }

    //clout.setMultiOutput(false);
  }


  /// Returns full mapping from fine locs to coarse locs
  const std::map<int,int>& clocs() const override {
    return _coarseLoc;
  }

  /// Returns coarse loc for given fine loc if there is one
  std::optional<int> cloc(int fineLoc) const override {
    auto iter = _coarseLoc.find(fineLoc);
    if (iter != _coarseLoc.end()) {
      auto [_, cloc] = *iter;
      return cloc;
    }
    return std::nullopt;
  }

};

}

#endif
