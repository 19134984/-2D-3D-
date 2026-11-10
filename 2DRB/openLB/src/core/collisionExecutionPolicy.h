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

#ifndef CORE_COLLISION_EXECUTION_POLICY_H
#define CORE_COLLISION_EXECUTION_POLICY_H

#include <cstdint>

namespace olb {

/// Whether collision is split into outer / inner sub-domains for
/// communication-computation overlap.
/**
 * The split lets the halo exchange of the outer (boundary ring) cells
 * proceed while the inner (bulk) collision still runs. It only pays off when
 * the exchange frees the GPU, i.e. an off-node transfer with the host
 * blocked. For on-GPU (local copy) or intra-node (NVLink) halos the transfer
 * is itself GPU work that contends with the bulk collision, so splitting is
 * pure overhead.
 **/
enum struct CollisionOverlap : std::uint8_t {
  Auto,   ///< resolved from the decomposition topology (split iff off-node halo)
  Split,  ///< always split (comm/compute overlap)
  Fused,  ///< never split (single collision over the whole sub-domain)
};

/// How collision iterates the cells it applies to.
/**
 * Mask iterates the sub-domain bitmask with coalesced, indirection-free
 * addressing (fast when the operator covers most of the domain). List iterates
 * an explicit cell-id list (fewer launched threads but indirect addressing,
 * better when the operator is sparse).
 **/
enum struct CollisionDispatch : std::uint8_t {
  Auto,  ///< resolved per operator from its mask density
  Mask,  ///< iterate the sub-domain bitmask
  List,  ///< iterate an explicit cell-id list
};

/// Per-cuboid collision execution configuration.
/**
 * Owned per local cuboid by the LoadBalancer (the surface on which the user
 * expresses parallel intent) and consulted by the platform collision
 * operators. Extend by adding further independent axes. Each defaults to Auto
 * so that unset axes keep resolving heuristically.
 **/
struct CollisionExecutionPolicy {
  /// Outer/inner split for comm overlap (resolved from topology)
  CollisionOverlap  overlap  = CollisionOverlap::Auto;
  /// Dispatch of the bulk (inner if split, else whole sub-domain)
  CollisionDispatch bulk     = CollisionDispatch::Auto;
  /// Dispatch of the boundary ring (only meaningful when split)
  CollisionDispatch boundary = CollisionDispatch::Auto;
};

}

#endif
