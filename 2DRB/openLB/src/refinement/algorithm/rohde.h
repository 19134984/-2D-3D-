/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Adrian Kummerlaender
 *
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

#ifndef REFINEMENT_ALGORITHM_ROHDE_H
#define REFINEMENT_ALGORITHM_ROHDE_H

#include "utilities/performanceProfiler.h"

namespace olb::refinement::rohde {

struct CoarseToFineO {
  static constexpr refinement::OperatorScope scope = refinement::OperatorScope::PerCellCenteredCoarseCell;

  using parameters = meta::list<>;

  using data = meta::list<fields::refinement::NORMAL_TYPE,
                          fields::refinement::NORMAL>;

  template <typename COARSE_CELL, typename FINE_CELL, typename DATA, typename PARAMETERS>
  void apply(COARSE_CELL& cCell, FINE_CELL& fMother, DATA& data, PARAMETERS& params) any_platform {
    using DESCRIPTOR = typename COARSE_CELL::descriptor_t;

    const DiscreteNormalType type = data.template getField<fields::refinement::NORMAL_TYPE>();
    const Vector<int,DESCRIPTOR::d> normal = data.template getField<fields::refinement::NORMAL>();

    for (unsigned iOrthant=0; iOrthant < (unsigned{1} << DESCRIPTOR::d); ++iOrthant) {
      const Vector<int,DESCRIPTOR::d> o_i([iOrthant](unsigned j) -> int {
        return 2*((iOrthant >> (DESCRIPTOR::d-1-j)) & 1) - 1;
      });
      auto fCell = fMother.child(o_i);

      for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
        const auto& c_i = descriptors::c<DESCRIPTOR>(iPop);

        bool isInward = (c_i * normal <= 0);
        if (isInward) {
          fCell[iPop] = cCell[iPop];
        }
      }
    }
  }
};

struct FineToCoarseO {
  static constexpr refinement::OperatorScope scope = refinement::OperatorScope::PerCellCenteredCoarseCell;

  using parameters = meta::list<>;

  using data = meta::list<fields::refinement::NORMAL_TYPE,
                          fields::refinement::NORMAL>;

  template <typename COARSE_CELL, typename FINE_CELL, typename DATA, typename PARAMETERS>
  void apply(COARSE_CELL& cCell, FINE_CELL& fMother, DATA& data, PARAMETERS& params) any_platform {
    using DESCRIPTOR = typename COARSE_CELL::descriptor_t;
    const DiscreteNormalType type = data.template getField<fields::refinement::NORMAL_TYPE>();
    const Vector<int,DESCRIPTOR::d> normal = data.template getField<fields::refinement::NORMAL>();

    for (unsigned iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
      const auto& c_i = descriptors::c<DESCRIPTOR>(iPop);

      bool isOutward = (c_i * normal > 0);
      if (isOutward) {
        cCell[iPop] = 0;
        for (unsigned iOrthant=0; iOrthant < (unsigned{1} << DESCRIPTOR::d); ++iOrthant) {
          const Vector<int,DESCRIPTOR::d> o_i([iOrthant](unsigned j) -> int {
            return 2*((iOrthant >> (DESCRIPTOR::d-1-j)) & 1) - 1;
          });
          auto fCell = fMother.child(o_i);
          cCell[iPop] += fCell[iPop];
        }
        cCell[iPop] /= (unsigned{1} << DESCRIPTOR::d);
      }
    }
  }
};

struct CoarseToFineForPostProcessingO {
  static constexpr refinement::OperatorScope scope = refinement::OperatorScope::PerCellCenteredCoarseCell;

  using parameters = meta::list<>;
  using data = meta::list<fields::refinement::NORMAL_TYPE,
                          fields::refinement::NORMAL>;

  template <typename COARSE_CELL, typename FINE_CELL, typename DATA, typename PARAMETERS>
  void apply(COARSE_CELL& cCell, FINE_CELL& fMother, DATA& data, PARAMETERS& params) any_platform {
    using V = typename COARSE_CELL::value_t;
    using DESCRIPTOR = typename COARSE_CELL::descriptor_t;

    const DiscreteNormalType type = data.template getField<fields::refinement::NORMAL_TYPE>();
    const Vector<int,DESCRIPTOR::d> normal = data.template getField<fields::refinement::NORMAL>();

    auto coarseNeighbor = cCell.neighbor(normal);

    V rhoC{};
    Vector<V,DESCRIPTOR::d> uC{};
    coarseNeighbor.computeRhoU(rhoC, uC.data());

    // Flat Face Interpolation
    if (type == DiscreteNormalType::Flat) {
      V avgRhoF = V{0};
      Vector<V,DESCRIPTOR::d> avgUF{};

      for (unsigned iOrthant=0; iOrthant < (unsigned{1} << DESCRIPTOR::d); ++iOrthant) {
        const Vector<int,DESCRIPTOR::d> o_i([iOrthant](unsigned j) -> int {
          return 2*((iOrthant >> (DESCRIPTOR::d-1-j)) & 1) - 1;
        });

        // Anchor sampling to the bulk fine domain along the normal direction
        Vector o_inner = o_i;
        for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
          if (normal[iD] != 0) {
            o_inner[iD] = -normal[iD];
          }
        }

        V rhoF{}; Vector<V,DESCRIPTOR::d> uF{};
        fMother.child(o_inner).neighbor(-1*normal).computeRhoU(rhoF, uF.data());

        avgRhoF += rhoF;
        avgUF += uF;
      }

      const V num = static_cast<V>(unsigned{1} << DESCRIPTOR::d);
      V gradRho = (rhoC - (avgRhoF / num)) / V{1.75};
      Vector<V,DESCRIPTOR::d> gradU{};
      for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
        gradU[iD] = (uC[iD] - (avgUF[iD] / num)) / V{1.75};
      }

      for (unsigned iOrthant=0; iOrthant < (unsigned{1} << DESCRIPTOR::d); ++iOrthant) {
        const Vector<int,DESCRIPTOR::d> o_i([iOrthant](unsigned j) -> int {
          return 2*((iOrthant >> (DESCRIPTOR::d-1-j)) & 1) - 1;
        });

        Vector o_inner = o_i;
        for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
          if (normal[iD] != 0) {
            o_inner[iD] = -normal[iD];
          }
        }

        V localRhoF{};
        Vector<V,DESCRIPTOR::d> localUF{};
        fMother.child(o_inner).neighbor(-1*normal).computeRhoU(localRhoF, localUF.data());

        // Simple projection of orthant onto normal for linear distance fraction
        V proj{};
        for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
          proj += o_i[iD] * normal[iD];
        }

        V scalar = (V{0.25} * proj + V{0.75});
        V interpolatedRho = localRhoF + gradRho * scalar;
        Vector interpolatedU = localUF + gradU * scalar;
        fMother.child(o_i).defineRhoU(interpolatedRho, interpolatedU.data());
      }
    }
    else {
      for (unsigned iOrthant=0; iOrthant < (unsigned{1} << DESCRIPTOR::d); ++iOrthant) {
        const Vector<int,DESCRIPTOR::d> o_i([iOrthant](unsigned j) -> int {
          return 2*((iOrthant >> (DESCRIPTOR::d-1-j)) & 1) - 1;
        });

        V accRho{};
        Vector<V,DESCRIPTOR::d> accU{};
        unsigned nAxis = 0;
        for (unsigned a=0; a < DESCRIPTOR::d; ++a) {
          if (normal[a] == 0) {
            continue;
          }
          Vector<int,DESCRIPTOR::d> n_a(0);
          n_a[a] = normal[a];

          V rhoCa{}; Vector<V,DESCRIPTOR::d> uCa{};
          cCell.neighbor(n_a).computeRhoU(rhoCa, uCa.data());

          Vector<int,DESCRIPTOR::d> o_inner = o_i;
          for (unsigned d=0; d < DESCRIPTOR::d; ++d) {
            if (normal[d] != 0) {
              o_inner[d] = -normal[d];
            }
          }
          V rhoFa{}; Vector<V,DESCRIPTOR::d> uFa{};
          fMother.child(o_inner).neighbor(-1*normal).computeRhoU(rhoFa, uFa.data());

          const V proj   = static_cast<V>(o_i[a] * normal[a]); // +/-1
          const V scalar = V{0.25} * proj + V{0.75};
          accRho += rhoFa + ((rhoCa - rhoFa) / V{1.75}) * scalar;
          for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
            accU[iD] += uFa[iD] + ((uCa[iD] - uFa[iD]) / V{1.75}) * scalar;
          }
          ++nAxis;
        }

        const V inv = V{1} / static_cast<V>(nAxis == 0 ? 1 : nAxis);
        V outRho = accRho * inv;
        Vector<V,DESCRIPTOR::d> outU = accU * inv;
        fMother.child(o_i).defineRhoU(outRho, outU.data());
      }
    }
  }
};

template <typename T, typename DESCRIPTOR>
std::unique_ptr<SuperLatticeRefinement<T,DESCRIPTOR>> makeCoupler(
  SuperLattice<T,DESCRIPTOR>& sLatticeCoarse,
  SuperGeometry<T,DESCRIPTOR::d>& sGeometryCoarse,
  SuperLattice<T,DESCRIPTOR>& sLatticeFine,
  SuperGeometry<T,DESCRIPTOR::d>& sGeometryFine,
  FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& coarseBulkI)
{
  auto& loadBalancerFine = dynamic_cast<RefinedLoadBalancer<T,DESCRIPTOR::d>&>(
    sLatticeFine.getLoadBalancer());
  auto& cDecompositionFine = sLatticeFine.getCuboidDecomposition();

  SuperIndicatorDomainFrontierDistanceF<T,DESCRIPTOR> frontierI(cDecompositionFine,
                                                                sGeometryFine,
                                                                0);
  SuperIndicatorDomainFrontierDistanceF<T,DESCRIPTOR> insideI(cDecompositionFine,
                                                              sGeometryFine,
                                                              1);
  SuperIndicatorDomainFrontierDistanceF<T,DESCRIPTOR> outsideI(cDecompositionFine,
                                                               sGeometryFine,
                                                               -1);

  auto coupler = std::make_unique<SuperLatticeRefinement<T,DESCRIPTOR>>(
    sLatticeCoarse, sLatticeFine, loadBalancerFine);

  for (auto [iC, coarseC] : loadBalancerFine.clocs()) {
    auto& cBlock = sLatticeCoarse.getBlock(coarseC);
    auto& fBlock = sLatticeFine.getBlock(iC);

    cBlock.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> coarseLatticeR) {
      if (   frontierI.getBlockIndicatorF(iC)(2*coarseLatticeR)
          && coarseBulkI->getBlockIndicatorF(coarseC)(coarseLatticeR)
          && cBlock.getNeighborhoodRadius(coarseLatticeR) >= 1) {
        coupler->getBlock(iC).addCellCentered(coarseLatticeR);

        // Bitshift (unsigned{1} << 2) results in 4 orthants for 2D case,
        // (unsigned{1} << 3) will results in 8 orthants for the 3D case
        unsigned int numberOrthants = unsigned{1} << DESCRIPTOR::d;

        for (unsigned iOrthant = 0; iOrthant < numberOrthants; ++iOrthant) {
          Vector<int, DESCRIPTOR::d> offset;
          for (int j = 0; j < DESCRIPTOR::d; ++j) {
            offset[j] = (iOrthant >> (DESCRIPTOR::d - 1 - j)) & 1;
          }
          sLatticeFine.getBlock(iC).defineDynamics(
            2 * coarseLatticeR + offset, meta::id<NoDynamicsWithFixedDensityAndVelocity<T,DESCRIPTOR>>());
        }
      }
    });

    cBlock.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> coarseLatticeR) {
      if (   frontierI.getBlockIndicatorF(iC)(2*coarseLatticeR)
          && coarseBulkI->getBlockIndicatorF(coarseC)(coarseLatticeR)
          && cBlock.getNeighborhoodRadius(coarseLatticeR) >= 1) {
        if (auto index = coupler->getBlock(iC).getCoarseDataIndex(coarseLatticeR)) {
          auto [type, normal] = computeBoundaryTypeAndNormal(insideI.getBlockIndicatorF(iC),
                                                             outsideI.getBlockIndicatorF(iC),
                                                             2*coarseLatticeR);
          coupler->getBlock(iC).getData()
                 .template getField<fields::refinement::NORMAL_TYPE>().set(*index, type);
          coupler->getBlock(iC).getData()
                 .template getField<fields::refinement::NORMAL>().set(*index, normal);

          auto neighborR = coarseLatticeR + normal;

          cBlock.template defineDynamics<NoDynamicsWithBulkMomenta>(coarseLatticeR);
          cBlock.template defineDynamics<NoDynamicsWithBulkMomenta>(neighborR);
        } else {
          throw std::logic_error("Invalid cell-centered coupling setup");
        }
      } else {
        if (cDecompositionFine.isInside(
              cDecompositionFine.getPhysR(
                (2*coarseLatticeR).withPrefix(loadBalancerFine.glob(iC))))) {
          cBlock.template defineDynamics<NoDynamicsWithBulkMomenta>(coarseLatticeR);
        }
      }
    });

    fBlock.forSpatialLocations([&](LatticeR<DESCRIPTOR::d> latticeR) {
      if (fBlock.isPadding(latticeR) && !cDecompositionFine.isInside(
          cDecompositionFine.getPhysR(latticeR.withPrefix(loadBalancerFine.glob(iC))))) {
        fBlock.template defineDynamics<NoDynamics>(latticeR);
      }
    });
  }

  coupler->setProcessingContext(ProcessingContext::Simulation);

  return coupler;
}

template <typename REFINED_LATTICE, typename PRE_STEP_HOOK>
void executeStep(REFINED_LATTICE& refinedLattice,
                 unsigned iLevel, std::size_t iT,
                 PRE_STEP_HOOK&& preStepHook)
{
  OLB_PROFILE("Refinement:Step");
  preStepHook(iLevel, iT);

  if (iLevel == (unsigned)refinedLattice.getLevelCount() - 1) {
    refinedLattice.getLevel(iLevel).collideAndStream();
  }
  else {
    refinedLattice.getLevel(iLevel).collide();
    refinedLattice.getCoupler(iLevel).apply(meta::id<CoarseToFineO>{});

    for (int subStep = 0; subStep < 2; ++subStep) {
      executeStep(refinedLattice, iLevel + 1, 2 * iT + subStep, preStepHook);
    }

    refinedLattice.getLevel(iLevel).stream();
    refinedLattice.getCoupler(iLevel).apply(meta::id<FineToCoarseO>{});
  }
}

}

#endif
