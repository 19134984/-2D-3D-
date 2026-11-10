/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Shota Ito
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

/*  This case validates the generated discrete-adjoint based collision kernel
 *  using forward mode AD. For an arbitrary primal collision kernel the partial
 *  derivatives of a single post-collision population regarding all pre-
 *  collision populations are computed. The gradients are compared against
 *  the primal collision executed with the forward AD mode.
 *
 *  This case should be used for validation of generated adjoint-kernels.
 *  Per default, DISABLE_CSE is off and (if it exists) the CSE optimized
 *  kernel is used in this test.
 */

//#define DISABLE_CSE
#include <olb.h>
#include "gtest/gtest.h"
#include "../tests/unit/helper/framework.h"

/// ensure olbInit is called for MPI
::testing::Environment* const foo_env = ::testing::AddGlobalTestEnvironment(new OlbInitTestEnvironment);

using namespace olb;

using T = double;

using HBGKSmagorinskyDynamics = dynamics::Tuple<
  T, descriptors::D3Q19<>,
  momenta::Porous<momenta::BulkTuple>,
  equilibria::SecondOrder,
  collision::SmagorinskyEffectiveOmega<collision::BGK>
>;
using HRRSmagorinskyDynamics = dynamics::Tuple<
  T, descriptors::D3Q19<>,
  momenta::Porous<momenta::BulkTuple>,
  equilibria::ThirdOrder,
  collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB>
>;
template <int orientation, int direction>
using ZouHePressure = ZouHeDynamics<T, descriptors::D3Q19<>,
  dynamics::Tuple<T, descriptors::D3Q19<>,
    momenta::BulkTuple,
    equilibria::SecondOrder,
    collision::BGK,
    dynamics::DefaultCombination>,
  momenta::Tuple<momenta::FixedDensity,
    momenta::FixedPressureMomentum<orientation, direction>,
    momenta::BulkStress,
    momenta::DefineSeparately>,
  orientation, direction>;

using StationaryLinearElasticityMRT =
  structuralmechanics::StationaryLinearElasticityMRT<T, descriptors::D2Q8<>>;

// Extend the list of Dynamics here
using DYNAMICS_LIST = meta::list<
  ForcedBGKdynamics<T,descriptors::D3Q19<>>,
  PorousBGKdynamics<T,descriptors::D3Q19<>>,
  HBGKSmagorinskyDynamics,
  HRRSmagorinskyDynamics,
  ZouHePressure<1,1>,
  StationaryLinearElasticityMRT
>;

constexpr T PRECISION = 1e-14;

template <typename DUAL_DYNAMICS, typename LATTICE>
void initializeADR(LATTICE& lattice) {
  using V = typename DUAL_DYNAMICS::value_t;
  using DYNAMICS = typename DUAL_DYNAMICS::primalDynamics;
  using DESCRIPTOR = typename LATTICE::descriptor_t;
  // Set collision model
  lattice.setDynamics(0,meta::id<DUAL_DYNAMICS>{});
  // Initialize population field and perform seeding
  auto cell = lattice.get(LatticeR<DESCRIPTOR::d>(0));
  cell.template getFieldPointer<descriptors::POPULATION>()[1] = V(1.);          // compute sensitivity regarding post-collision f[1]
  for (int iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
    cell.template getFieldPointer<opti::F>()[iPop] = V(1.);                     // evaluate the jacobian matrix at this position
  }
  // Initialize accessed fields
  DYNAMICS::fields::for_each([&](auto id) {
    using FIELD = typename decltype(id)::type;
    for (std::size_t iDim=0; iDim<FIELD::template size<DESCRIPTOR>(); ++iDim) {
      cell.template getFieldPointer<FIELD>()[iDim] = V(0.1);
    }
  });
  // Initialize used parameters
  DYNAMICS::parameters::for_each([&](auto id) {
    using PARAMETER = typename decltype(id)::type;
    for (std::size_t iDim=0; iDim<PARAMETER::template size<DESCRIPTOR>(); ++iDim) {
      lattice.template setParameter<PARAMETER>(V(0.1));
    }
  });
}

template <concepts::DifferentiableDynamics DYNAMICS, typename LATTICE>
void initializeADF(LATTICE& lattice) {
  using V = typename DYNAMICS::value_t;
  using DESCRIPTOR = typename LATTICE::descriptor_t;
  // Set collision model
  lattice.setDynamics(0,meta::id<DYNAMICS>{});
  // Initialize population field and perform seeding
  auto cell = lattice.get(LatticeR<DESCRIPTOR::d>(0));
  for (int iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
    cell.template getFieldPointer<descriptors::POPULATION>()[iPop] = V(1.);                 // Compute primal collision at this position
    cell.template getFieldPointer<descriptors::POPULATION>()[iPop].setDiffVariable(iPop);   // Compute sensitivity regarding all pre-collision f
  }
  // Initialize accessed fields
  DYNAMICS::fields::for_each([&](auto id) {
    using FIELD = typename decltype(id)::type;
    for (std::size_t iDim=0; iDim<FIELD::template size<DESCRIPTOR>(); ++iDim) {
      cell.template getFieldPointer<FIELD>()[iDim] = V(0.1);
    }
  });
  // Initialize used parameters
  DYNAMICS::parameters::for_each([&](auto id) {
    using PARAMETER = typename decltype(id)::type;
    for (std::size_t iDim=0; iDim<PARAMETER::template size<DESCRIPTOR>(); ++iDim) {
      lattice.template setParameter<PARAMETER>(V(0.1));
    }
  });
}

template <typename LATTICE>
void dumpDataADR(LATTICE& lattice) {
  using DESCRIPTOR = typename LATTICE::descriptor_t;
  auto cell = lattice.get(LatticeR<DESCRIPTOR::d>(0));
  for (int iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
    std::cout << "ADR: f[" << iPop << "] = " << cell[iPop] << std::endl; // Print all pre-collision sensitivities
  }
  std::cout << "\n";
}

template <typename LATTICE>
void dumpDataADF(LATTICE& lattice) {
  using DESCRIPTOR = typename LATTICE::descriptor_t;
  auto cell = lattice.get(LatticeR<DESCRIPTOR::d>(0));
  for (int iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
    std::cout << "ADF: f[1] = " << cell[1].d()[iPop] << std::endl;       // Print sensitivity of post-collision f[1]
  }
  std::cout << "\n";
}

TEST(AdjointCollisionKernels, testForcedBGKdynamics) {
  DYNAMICS_LIST::for_each([&](auto id) {
    using DYNAMICS = typename decltype(id)::type;
    using DESCRIPTOR = typename DYNAMICS::descriptor_t;
    using V    = util::ADf<T,DESCRIPTOR::q>;
    using DYNAMICS_ADR = Dual<DYNAMICS>;
    using DYNAMICS_ADF = DYNAMICS::template exchange_value_type<V>;
    Vector<int,DESCRIPTOR::d> size(1);

    // Create a BlockLattice with a single cell to run a dummy simulation
    ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD> lattice(size, 0);
    ConcreteBlockLattice<V,DESCRIPTOR,Platform::CPU_SISD> latticeADF(size, 0);
    lattice.setStatisticsEnabled(false);
    latticeADF.setStatisticsEnabled(false);
    // Initialize simulation for forward and reverse AD collision
    initializeADR<DYNAMICS_ADR>(lattice);
    initializeADF<DYNAMICS_ADF>(latticeADF);
    // Print population field pre collide
    std::cout << "PreCollide:" << std::endl;
    dumpDataADR(lattice);
    dumpDataADF(latticeADF);
    // Perform single collision
    lattice.collide();
    latticeADF.collide();
    // Print population field post collide
    std::cout << "PostCollide:" << std::endl;
    dumpDataADR(lattice);
    dumpDataADF(latticeADF);
    // Check if results are identical
    const auto origin = LatticeR<DESCRIPTOR::d>(0);
    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      EXPECT_NEAR(lattice.get(origin)[iPop], latticeADF.get(origin)[1].d()[iPop], PRECISION);
    }
  });
}
