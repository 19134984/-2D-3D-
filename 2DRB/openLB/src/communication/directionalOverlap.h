/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Adrian Kummerlaender
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

#ifndef COMMUNICATION_DIRECTIONAL_OVERLAP_H
#define COMMUNICATION_DIRECTIONAL_OVERLAP_H

#include <cstdint>
#include <optional>
#include <vector>
#include <typeindex>

#include "core/vector.h"
#include "descriptor/descriptor.h"

namespace olb {

/// Population-set math for directional (outgoing/incoming-only) overlap exchange
/**
 * The streaming overlap exchange only needs to transmit, per padding cell p, the
 * populations that post-shift propagation actually reads from p:
 *
 *   R(p) = { i : p + c_i ∈ core }
 *
 * For a cell one layer outside the core with padding signature d ∈ {-1,0,1}^D
 * (per axis: -1 below the core, +1 above, 0 inside) this is the canonical
 *
 *   M(d) = { i : ∀a with d_a ≠ 0 : c_i[a] = -d_a }
 *
 * e.g. 5/19 populations on a D3Q19 face, 1/19 on an edge and none on a
 * corner. Cells deeper than one layer are never read by propagation
 * (|c_i| ≤ 1) and drop out entirely.
 *
 * All functions are negotiation-time helpers (not hot-path) returning
 * population subsets as uint64 bitmasks.
 **/
namespace overlap {

/// Field components carried per streaming direction: 1 for a scalar lattice,
/// DESCRIPTOR::n for a vector lattice packing n components per velocity. The
/// directional mask carries one bit per direction, so a set bit selects the
/// direction's contiguous component block [k*n, k*n+n) = descriptors::q(k,·).
template <typename DESCRIPTOR>
constexpr unsigned componentsPerDirection()
{
  if constexpr (requires { DESCRIPTOR::n; }) {
    return DESCRIPTOR::n;
  } else {
    return 1;
  }
}

/// Number of streaming directions of DESCRIPTOR (q, or q/n for a vector lattice)
template <typename DESCRIPTOR>
constexpr unsigned directionCount()
{
  return DESCRIPTOR::q / componentsPerDirection<DESCRIPTOR>();
}

/// Masked-segment element size and components per mask bit of one field
struct MaskedFieldTraits {
  std::size_t elementSize;
  unsigned componentsPerDirection;
};

/// Resolve `field`'s masked traits: a mask bit selects a streaming direction
/// of a propagatable field resp. one component otherwise (component count
/// must fit the 64-bit mask; wider fields stay dense)
template <typename T, typename DESCRIPTOR>
std::optional<MaskedFieldTraits> maskedFieldTraits(std::type_index field)
{
  std::optional<MaskedFieldTraits> traits;
  DESCRIPTOR::fields_t::for_each([&](auto f) {
    using FIELD = typename decltype(f)::type;
    constexpr bool propagatable = descriptors::is_propagatable_field<FIELD>::value;
    if constexpr (propagatable || DESCRIPTOR::template size<FIELD>() <= 64) {
      if (field == std::type_index(typeid(FIELD))) {
        traits = MaskedFieldTraits{
          sizeof(precision_scalar_t<typename FIELD::template value_type<T>, precisionOf<DESCRIPTOR,FIELD>()>),
          propagatable ? componentsPerDirection<DESCRIPTOR>() : 1u
        };
      }
    }
  });
  return traits;
}

/// Full direction mask for DESCRIPTOR
template <typename DESCRIPTOR>
constexpr std::uint64_t densePopulations()
{
  constexpr unsigned Ndir = directionCount<DESCRIPTOR>();
  static_assert(Ndir <= 64, "directional mask is a uint64 bitmask");
  if constexpr (Ndir == 64) {
    return ~std::uint64_t{0};
  } else {
    return (std::uint64_t{1} << Ndir) - 1;
  }
}

/// Padding signature of latticeR w.r.t. the core box [0, core)
/**
 * Per axis: -1 if below the core, +1 if above, 0 if inside. The signature
 * does not encode depth — use requiredPopulations(latticeR, core) for the
 * exact per-cell set.
 **/
template <unsigned D>
constexpr Vector<int,D> paddingSignature(LatticeR<D> latticeR, LatticeR<D> core)
{
  Vector<int,D> d{};
  for (unsigned a=0; a < D; ++a) {
    d[a] = (latticeR[a] < 0) ? -1 : (latticeR[a] >= core[a] ? 1 : 0);
  }
  return d;
}

/// Canonical directional mask M(d) for a depth-1 padding signature d
/**
 * Upper bound on the exact per-cell set: R(p) ⊆ M(d(p)) for depth-1 cells,
 * with equality only when p is laterally interior — at cells near the core's
 * edges, populations with tangential components stream past the core and
 * drop out of R(p). Negotiation must therefore use the exact
 * requiredPopulations(latticeR, core) per cell; M(d) serves as documentation
 * and as a superset invariant. M(0) is the dense mask.
 **/
template <typename DESCRIPTOR>
constexpr std::uint64_t directionalPopulations(Vector<int,DESCRIPTOR::d> d)
{
  constexpr unsigned N = componentsPerDirection<DESCRIPTOR>();
  constexpr unsigned Ndir = DESCRIPTOR::q / N;
  static_assert(Ndir <= 64, "directional mask is a uint64 bitmask");
  std::uint64_t mask{};
  for (unsigned iDir=0; iDir < Ndir; ++iDir) {
    bool required = true;
    for (unsigned a=0; a < DESCRIPTOR::d; ++a) {
      if (d[a] != 0 && descriptors::c<DESCRIPTOR>(iDir*N,a) != -d[a]) {
        required = false;
        break;
      }
    }
    if (required) {
      mask |= std::uint64_t{1} << iDir;
    }
  }
  return mask;
}

/// Exact required set R(p) = { i : p + c_i ∈ core } for any padding cell
/**
 * This is the defining property: population i of cell p is only consumed by
 * post-shift propagation if it streams into a core cell. Equals
 * directionalPopulations(paddingSignature(p)) for depth-1 cells and the
 * empty set for deeper ones.
 **/
template <typename DESCRIPTOR>
constexpr std::uint64_t requiredPopulations(LatticeR<DESCRIPTOR::d> latticeR,
                                            LatticeR<DESCRIPTOR::d> core)
{
  constexpr unsigned N = componentsPerDirection<DESCRIPTOR>();
  constexpr unsigned Ndir = DESCRIPTOR::q / N;
  static_assert(Ndir <= 64, "directional mask is a uint64 bitmask");
  std::uint64_t mask{};
  for (unsigned iDir=0; iDir < Ndir; ++iDir) {
    bool inside = true;
    for (unsigned a=0; a < DESCRIPTOR::d; ++a) {
      const int target = latticeR[a] + descriptors::c<DESCRIPTOR>(iDir*N,a);
      if (target < 0 || target >= core[a]) {
        inside = false;
        break;
      }
    }
    if (inside) {
      mask |= std::uint64_t{1} << iDir;
    }
  }
  return mask;
}

/// Population count of a mask
constexpr unsigned populationCount(std::uint64_t mask)
{
  unsigned count = 0;
  while (mask) {
    count += mask & 1;
    mask >>= 1;
  }
  return count;
}

/// Type indices of DESCRIPTOR's propagatable fields (the shared streaming set)
/**
 * Every propagatable field shares the base descriptor stencil, hence one
 * R(p) and one negotiated mask. The comm layer nonetheless iterates this
 * set so the byte layout carries a segment per propagatable field.
 **/
template <typename DESCRIPTOR>
std::vector<std::type_index> propagatableFields()
{
  std::vector<std::type_index> fields;
  DESCRIPTOR::template filter<descriptors::is_propagatable_field>
            ::for_each([&](auto field) {
    fields.emplace_back(typeid(typename decltype(field)::type));
  });
  return fields;
}

/// Propagatable-field type indices of BLOCK, empty if BLOCK carries no descriptor
template <typename BLOCK>
std::vector<std::type_index> propagatableFieldsOf()
{
  if constexpr (requires { typename BLOCK::descriptor_t; }) {
    return propagatableFields<typename BLOCK::descriptor_t>();
  } else {
    return {};
  }
}

/// Component selector for directional streaming exchange of POPULATION
/**
 * For use with BlockCommunicationNeighborhood::requestField: restricts the
 * exchanged populations of each overlap cell to the exact set consumed by
 * post-shift propagation, R(p). The neighborhood supplies its block's core
 * extent as the second argument.
 **/
template <typename DESCRIPTOR>
auto streamingMask()
{
  return [](LatticeR<DESCRIPTOR::d> latticeR, LatticeR<DESCRIPTOR::d> core) -> std::uint64_t {
    return requiredPopulations<DESCRIPTOR>(latticeR, core);
  };
}

}

}

#endif
