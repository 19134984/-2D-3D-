/*  This file is part of the OpenLB library
 *
 *  Performance/GPU-portability variant of the convective-outflow boundary
 *  (sc26 SC submission): equilibrium-based StraightConvection without PREV_CELL.
 *  The canonical PREV_CELL version remains in boundaryPostProcessors3D.{h,hh}.
*/
#ifndef OLB_BOUNDARY_OPTIMIZED_CONVECTION_3D_H
#define OLB_BOUNDARY_OPTIMIZED_CONVECTION_3D_H

#include "boundaryPostProcessors3D.h"
#include "setBoundary.h"

namespace olb {

namespace optimized {

template <typename T, typename DESCRIPTOR, int direction, int orientation>
class StraightConvectionBoundaryProcessor3D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  struct PREV_CELL : public descriptors::FIELD_BASE<
    util::populationsContributingToVelocity<DESCRIPTOR,direction,-orientation>().size()
  > { };

  int getPriority() const {
    return 0;
  }

  //template <concepts::DynamicCell CELL>
  //void initialize(CELL& cell) any_platform;

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform;

};

template <typename T, typename DESCRIPTOR, int direction, int orientation>
template <concepts::DynamicCell CELL>
void StraightConvectionBoundaryProcessor3D<T,DESCRIPTOR,direction,orientation>::apply(CELL& cell) any_platform
{
  using V = typename CELL::value_t;
  constexpr auto missing = util::populationsContributingToVelocity<DESCRIPTOR,direction,-orientation>();

  V rho = 1;
  Vector<V,DESCRIPTOR::d> u{};
  cell.neighbor({-1, 0, 0}).computeU(u.data());

  Vector<V,DESCRIPTOR::q> fEq{};
  equilibria::ThirdOrder::template type<DESCRIPTOR,momenta::BulkTuple>().compute(cell, rho, u, fEq);

  for (unsigned i=0; i < missing.size(); ++i) {
    cell[missing[i]] = fEq[missing[i]];
  }
}

/// Convective-outflow boundary using the equilibrium-based (no PREV_CELL) processor
template <concepts::BaseType T, concepts::LatticeDescriptor DESCRIPTOR>
requires (DESCRIPTOR::d == 3)
struct InterpolatedConvection {
  using value_t      = T;
  using descriptor_t = DESCRIPTOR;

  CellDistance getNeighborhoodRadius() { return 1; }

  std::optional<DynamicsPromise<T,DESCRIPTOR>> getDynamics(DiscreteNormalType type,
                                                           DiscreteNormal<DESCRIPTOR> n)
  {
    return std::nullopt;
  }

  std::optional<PostProcessorPromise<T,DESCRIPTOR>> getPostProcessor(DiscreteNormalType type,
                                                                     DiscreteNormal<DESCRIPTOR> n)
  {
    switch (type) {
      case DiscreteNormalType::Flat:
        return boundaryhelper::promisePostProcessorForDirectionOrientation<
            T, DESCRIPTOR, StraightConvectionBoundaryProcessor3D>(n);
      default:
        return std::nullopt;
    }
  }
};

}  // namespace optimized

}  // namespace olb

#endif
