/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Florian Kaiser
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

#ifndef STRUCTURALMECHANICS_STRUCTURALMECHANICS_H
#define STRUCTURALMECHANICS_STRUCTURALMECHANICS_H

namespace olb {
namespace structuralmechanics {

  // Structural Mechanics (General)
  struct LAME_MU                  : public descriptors::FIELD_BASE<1,  0, 0> { };
  struct LAME_LAMBDA              : public descriptors::FIELD_BASE<1,  0, 0> { };
  struct BULK_MODULUS             : public descriptors::FIELD_BASE<1,  0, 0> { };
  struct DISPLACEMENT             : public descriptors::FIELD_BASE<0,  1, 0> { };
  struct VELOCITY                 : public descriptors::FIELD_BASE<0,  1, 0> { };
  struct CAUCHY_STRESS            : public descriptors::FIELD_BASE_CUSTOM_SIZE {
    template <typename DESCRIPTOR>
    static constexpr unsigned size() {
      return (DESCRIPTOR::d * (DESCRIPTOR::d));
    }

    template <typename T, typename DESCRIPTOR>
    static constexpr auto getInitialValue() {
      return Vector<value_type<T>, size<DESCRIPTOR>()>{};
    }
  };

  // Boolakee Quasistationary
  struct MOMENT_VECTOR            : public descriptors::FIELD_BASE<0,  0, 1> { };
  struct DAMPING_FACTOR           : public descriptors::FIELD_BASE<1,  0, 0> { };
  struct BARED_MOMENT_VECTOR      : public descriptors::FIELD_BASE<0,  0, 1> { };
  struct PREVIOUS_CELL            : public descriptors::FIELD_BASE<0,  0, 1> { };
  struct CURRENT_CELL             : public descriptors::FIELD_BASE<0,  0, 1> { };
  struct OMEGA_SOLID              : public descriptors::FIELD_BASE<0,  0, 1> { };
  struct CHAR_U                   : public descriptors::FIELD_BASE<1,  0, 0> { };
  struct CHAR_SIGMA               : public descriptors::FIELD_BASE<1,  0, 0> { };

  // Boolakee Elastodynamics
  struct CHAR_V                   : public descriptors::FIELD_BASE<1,  0, 0> { };
  struct CHAR_L                   : public descriptors::FIELD_BASE<1,  0, 0> { };
  struct CHAR_T                   : public descriptors::FIELD_BASE<1,  0, 0> { };
  struct C_K                      : public descriptors::FIELD_BASE<1,  0, 0> { };
  struct C_MU                     : public descriptors::FIELD_BASE<1,  0, 0> { };
  struct LATTICE_SPEED            : public descriptors::FIELD_BASE<1,  0, 0> { };

  // Initial macroscopic state U0
  template <unsigned N>
  struct U_ZERO_N    : public descriptors::FIELD_BASE<N> {};
  template <unsigned N>
  struct DX_U_ZERO_N : public descriptors::FIELD_BASE<N> {};
  template <unsigned N>
  struct DY_U_ZERO_N : public descriptors::FIELD_BASE<N> {};

  template <typename DESCRIPTOR> using U_ZERO = U_ZERO_N<DESCRIPTOR::n>;
  template <typename DESCRIPTOR> using DX_U_ZERO = DX_U_ZERO_N<DESCRIPTOR::n>;
  template <typename DESCRIPTOR> using DY_U_ZERO = DY_U_ZERO_N<DESCRIPTOR::n>;

  // Faust Elastodynamics
  struct GRAD_RHO                 : public descriptors::FIELD_BASE<0,  1, 0> { };
  struct SOLID_SOURCE             : public descriptors::FIELD_BASE<0,  1, 0> { };
  struct FULL_SOURCE              : public descriptors::FIELD_BASE<0,  0, 1> { };
  struct J                        : public descriptors::FIELD_BASE<0,  1, 0> { };
  struct RHO0                     : public descriptors::FIELD_BASE<1,  0, 0> { };
  struct POISSON_STRESS           : public descriptors::FIELD_BASE_CUSTOM_SIZE {
    template <typename DESCRIPTOR>
    static constexpr unsigned size() {
      return DESCRIPTOR::d * DESCRIPTOR::d;
    }

    template <typename V, typename DESCRIPTOR>
    static constexpr V getInitialValue() {
      return V();
    }
  };

  struct P_PI                     : public descriptors::FIELD_BASE_CUSTOM_SIZE {
    template <typename DESCRIPTOR>
    static constexpr unsigned size() {
      return DESCRIPTOR::d * DESCRIPTOR::d;
    }

    template <typename V, typename DESCRIPTOR>
    static constexpr V getInitialValue() {
      return V();
    }
  };

  struct FIRST_PIOLA_KIRCHHOFF    : public descriptors::FIELD_BASE_CUSTOM_SIZE {
    template <typename DESCRIPTOR>
    static constexpr unsigned size() {
      return DESCRIPTOR::d * DESCRIPTOR::d;
    }

    template <typename V, typename DESCRIPTOR>
    static constexpr V getInitialValue() {
      return V();
    }
  };

  // Boundary Conditions
  // Vector LBM
  struct VLBM_SOLID_DISTANCE_FIELD : public descriptors::FIELD_BASE_CUSTOM_SIZE {

    template <typename DESCRIPTOR>
    requires requires {DESCRIPTOR::n; }
    static constexpr unsigned size() {
      return DESCRIPTOR::q / DESCRIPTOR::n;
    }

    template <typename V, typename DESCRIPTOR>
    static constexpr V getInitialValue() {
      return V();
    }
  };

  // Traditional LBM
  struct SOLID_DISTANCE_FIELD     : public descriptors::FIELD_BASE<0,  0, 1> { };
  struct NEUMANN_NORMAL_X         : public descriptors::FIELD_BASE<0,  0, 1> { };
  struct NEUMANN_NORMAL_Y         : public descriptors::FIELD_BASE<0,  0, 1> { };

  // Quasistationary
  struct NEUMANN_SOLID_C          : public descriptors::FIELD_BASE<3,  0, 0> { };


  // Validation
  struct BOUNDARY_COORDS_X        : public descriptors::FIELD_BASE<0,  0, 1> { };
  struct BOUNDARY_COORDS_Y        : public descriptors::FIELD_BASE<0,  0, 1> { };
  // Reference Lattice moduli of the material for optimization problems
  struct REFERENCE_LAME_MU        : public descriptors::FIELD_BASE<1,  0, 0> { };
  struct REFERENCE_BULK_MODULUS   : public descriptors::FIELD_BASE<1,  0, 0> { };
  struct MS_AMPLITUDE             : public descriptors::FIELD_BASE<1,  0, 0> { };
  struct MS_WAVENUMBER            : public descriptors::FIELD_BASE<1,  0, 0> { };

  template <typename L1, typename L2> struct concat_parameter_lists;
  template <typename... A, typename... B>
  struct concat_parameter_lists<meta::list<A...>, meta::list<B...>> {
    using type = meta::list<A..., B...>;
  };

} // namespace structuralmechanics

// Traditional LBM Version
// DESCRIPTOR::n must not exist because of traditional populations
template <typename V, typename DESCRIPTOR, typename OPERATOR>
requires (!requires { DESCRIPTOR::n; })
void setStructuralMechanicsBoundary(BlockLattice<V,DESCRIPTOR>& block,
                                    BlockGeometry<V,DESCRIPTOR::d>& blockGeometry,
                                    BlockIndicatorF<V,DESCRIPTOR::d>& boundaryIndicator,
                                    BlockIndicatorF<V,DESCRIPTOR::d>& bulkIndicator,
                                    IndicatorF<V,DESCRIPTOR::d>& indicatorAnalyticalBoundary,
                                    bool outward,
                                    bool verbose) {
  OstreamManager clout(std::cout, "StructuralMechanicsBoundaryIterator");
  clout.setMultiOutput(true);

  int oppo;
  Vector<V, DESCRIPTOR::d> c, oppoC, normal;
  V dist, qIpop, norm;
  bool hasBoundary;

  const V deltaR = blockGeometry.getDeltaR();
  // for each solid cell: all of its fluid neighbors need population updates
  block.forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> node) {
    // Check if cell is boundary cell
    if (boundaryIndicator(node)) {
      // iterate over all distinct velocities
      for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
        c = descriptors::c<DESCRIPTOR>(iPop);
        Vector<V,DESCRIPTOR::d> bulk(node + c);

        // check if neighbor is fluid cell
        if (bulkIndicator(bulk)) {
          // Get coordinates of neighbor node because we need to calculate q_ij
          Vector<V, DESCRIPTOR::d> physCoords;
          blockGeometry.getPhysR(physCoords, bulk);
          Vector<V,DESCRIPTOR::d> bulkPhysR = {physCoords[0], physCoords[1]};

          oppo = descriptors::opposite<DESCRIPTOR>(iPop);
          oppoC = descriptors::c<DESCRIPTOR>(oppo);
          dist = V(-1);    // distance to boundary
          qIpop = V(-1);   // normed distance (Bouzidi distance) to boundary
          norm = deltaR * util::norm<DESCRIPTOR::d>(c);
          hasBoundary = indicatorAnalyticalBoundary.distance(dist, bulkPhysR, oppoC, blockGeometry.getIcGlob());
          // Check if distance calculation was performed correctly
          if (hasBoundary) {
            qIpop = dist / norm;
            // if distance function returned a dist. not suitable for Bouzidi -> fall-back
            if ((qIpop < 0) || (qIpop > 1)) {
              if (verbose) {
                clout << "Error, non suitable dist " << dist << ", for q_ij " << qIpop << " at lattice: (" << bulk << "), physical: (" << blockGeometry.getPhysR(bulk) << "), direction " << iPop << ". Fall-back to bounce-back." << std::endl;
              }
              // fall-back: half-way bounce back
              qIpop = V(0.5);
            }
          }
          // if distance function couldn't compute any distance -> fall-back
          else {
            if (verbose){
              clout << "Error, no boundary found at lattice:(" << bulk << "), physical: (" << blockGeometry.getPhysR(bulk) << "), direction: " << iPop << ".Fall-back to bounce-back." << std::endl;
            }
            // fall-back: half-way bounce back
            qIpop = V(0.5);
          }
          block.get(bulk).template setFieldComponent<structuralmechanics::SOLID_DISTANCE_FIELD>(oppo, qIpop);

          block.get(node).template setField<descriptors::CELL_COORDS>(node * deltaR);

          Vector<V,2> boundary = {
            bulkPhysR[0] + oppoC[0] * qIpop * deltaR,
            bulkPhysR[1] + oppoC[1] * qIpop * deltaR,
          };

          block.get(bulk).template setFieldComponent<structuralmechanics::BOUNDARY_COORDS_X>(oppo, boundary[0]);
          block.get(bulk).template setFieldComponent<structuralmechanics::BOUNDARY_COORDS_Y>(oppo, boundary[1]);

          bool hasNormal = indicatorAnalyticalBoundary.normal(normal, boundary, oppoC);
          if (! hasNormal) {
            if (verbose) {
              clout << "Normal not found at lattice:(" << bulk << "), physical: (" << blockGeometry.getPhysR(bulk) << "), direction: " << iPop << ". Fall-back to normal [0 0]." << std::endl;
            }
            normal[0] = 0;
            normal[1] = 0;
          }
          if (! outward) {
            normal[0] = -1 * normal[0];
            normal[1] = -1 * normal[1];
          }
          block.get(bulk).template setFieldComponent<structuralmechanics::NEUMANN_NORMAL_X>(oppo, normal[0]);
          block.get(bulk).template setFieldComponent<structuralmechanics::NEUMANN_NORMAL_Y>(oppo, normal[1]);

        }
      }
      block.addPostProcessor(typeid(stage::PostCollide),
      node,
      meta::id<OPERATOR>{});
    }
  });
};

// VectorLBM Version
// Requires DESCRIPTOR::n to exist because of Vector Cell
template <typename V, typename DESCRIPTOR, typename OPERATOR>
requires requires { DESCRIPTOR::n; }
void setStructuralMechanicsBoundary(BlockLattice<V,DESCRIPTOR>& block,
                                    BlockGeometry<V,DESCRIPTOR::d>& blockGeometry,
                                    BlockIndicatorF<V,DESCRIPTOR::d>& boundaryIndicator,
                                    BlockIndicatorF<V,DESCRIPTOR::d>& bulkIndicator,
                                    IndicatorF<V,DESCRIPTOR::d>& indicatorAnalyticalBoundary,
                                    bool outward,
                                    bool verbose) {
  OstreamManager clout(std::cout, "StructuralMechanicsBoundaryIterator");
  clout.setMultiOutput(true);

  int oppo;
  Vector<V, DESCRIPTOR::d> c, oppoC, normal;
  V dist, qIpop, norm;
  bool hasBoundary;

  const V deltaR = blockGeometry.getDeltaR();
  // for each solid cell: all of its fluid neighbors need population updates
  block.forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> node) {
    // Check if cell is boundary cell
    if (boundaryIndicator(node)) {
      // iterate over all distinct velocities
      for (int iPopV = 0; iPopV < DESCRIPTOR::base_descriptor_t::q; ++iPopV) {
        c = descriptors::cVector<DESCRIPTOR>(iPopV);
        Vector<V,DESCRIPTOR::d> bulk(node + c);

        // check if neighbor is fluid cell
        if (bulkIndicator(bulk)) {
          // Get coordinates of neighbor node because we need to calculate q_ij
          Vector<V, DESCRIPTOR::d> physCoords;
          blockGeometry.getPhysR(physCoords, bulk);
          Vector<V,DESCRIPTOR::d> bulkPhysR = {physCoords[0], physCoords[1]};

          oppo = descriptors::oppositeVector<DESCRIPTOR>(iPopV);
          oppoC = descriptors::cVector<DESCRIPTOR>(oppo);
          dist = V(-1);    // distance to boundary
          qIpop = V(-1);   // normed distance (Bouzidi distance) to boundary
          norm = deltaR * util::norm<DESCRIPTOR::d>(c);
          hasBoundary = indicatorAnalyticalBoundary.distance(dist, bulkPhysR, oppoC, blockGeometry.getIcGlob());
          // Check if distance calculation was performed correctly
          if (hasBoundary) {
            qIpop = dist / norm;
            // if distance function returned a dist. not suitable for Bouzidi -> fall-back
            if ((qIpop < 0) || (qIpop > 1)) {
              if (verbose) {
                clout << "Error, non suitable dist " << dist << ", for q_ij " << qIpop << " at lattice: (" << bulk << "), physical: (" << blockGeometry.getPhysR(bulk) << "), direction " << iPopV << ". Fall-back to bounce-back." << std::endl;
              }
              // fall-back: half-way bounce back
              qIpop = V(0.5);
            }
          }
          // if distance function couldn't compute any distance -> fall-back
          else {
            if (verbose) {
              clout << "Error, no boundary found at lattice:(" << bulk << "), physical: (" << blockGeometry.getPhysR(bulk) << "), direction: " << iPopV << ".Fall-back to bounce-back." << std::endl;
            }
            // fall-back: half-way bounce back
            qIpop = V(0.5);
          }
          block.get(bulk).template setFieldComponent<structuralmechanics::VLBM_SOLID_DISTANCE_FIELD>(oppo, qIpop);

          block.get(node).template setField<descriptors::CELL_COORDS>(node * deltaR);

          Vector<V,2> boundary = {
            bulkPhysR[0] + oppoC[0] * qIpop * deltaR,
            bulkPhysR[1] + oppoC[1] * qIpop * deltaR,
          };

          block.get(bulk).template setFieldComponent<structuralmechanics::BOUNDARY_COORDS_X>(oppo, boundary[0]);
          block.get(bulk).template setFieldComponent<structuralmechanics::BOUNDARY_COORDS_Y>(oppo, boundary[1]);

          bool hasNormal = indicatorAnalyticalBoundary.normal(normal, boundary, oppoC);
          if (! hasNormal) {
            if (verbose) {
              clout << "Normal not found at lattice:(" << bulk << "), physical: (" << blockGeometry.getPhysR(bulk) << "), direction: " << iPopV << ". Fall-back to normal [0 0]." << std::endl;
            }
            normal[0] = 0;
            normal[1] = 0;
          }
          if (! outward) {
            normal[0] = -1 * normal[0];
            normal[1] = -1 * normal[1];
          }
          block.get(bulk).template setFieldComponent<structuralmechanics::NEUMANN_NORMAL_X>(oppo, normal[0]);
          block.get(bulk).template setFieldComponent<structuralmechanics::NEUMANN_NORMAL_Y>(oppo, normal[1]);

        }
      }
      block.addPostProcessor(typeid(stage::PostCollide),
      node,
      meta::id<OPERATOR>{});
    }
  });
};

template<typename V, typename DESCRIPTOR, typename OPERATOR>
void setStructuralMechanicsBoundary(SuperLattice<V, DESCRIPTOR>& sLattice,
                                    FunctorPtr<SuperIndicatorF<V,DESCRIPTOR::d>>&& boundaryIndicator,
                                    FunctorPtr<SuperIndicatorF<V,DESCRIPTOR::d>>&& bulkIndicator,
                                    IndicatorF<V,DESCRIPTOR::d>& indicatorAnalyticalBoundary,
                                    bool outward = true,
                                    bool verbose = true) {
  OstreamManager clout(std::cout, "StructuralMechanicsBoundarySetter");
  auto& load = sLattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    setStructuralMechanicsBoundary<V,DESCRIPTOR, OPERATOR>(sLattice.getBlock(iC),
                                                           (bulkIndicator->getBlockIndicatorF(iC)).getBlockGeometry(),
                                                           boundaryIndicator->getBlockIndicatorF(iC),
                                                           bulkIndicator->getBlockIndicatorF(iC),
                                                           indicatorAnalyticalBoundary,
                                                           outward,
                                                           verbose);
  }
  clout.setMultiOutput(false);
}

}   // namespace olb

#include "elastodynamics/boolakee_elastodynamics.h"
#include "quasistationary/boolakee_quasistationary.h"

#endif
