/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Mathias J. Krause
 *                2024 Julius Jessberger, Adrian Kummerlaender
 *                2026 Shota Ito, Ecem Öz
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

#ifndef DISCRETE_DUAL_DYNAMICS_H
#define DISCRETE_DUAL_DYNAMICS_H

#include "utilities/vectorHelpers.h"
#include "utilities/aDiffTape.h"

#include "dynamics/interface.h"

// All OpenLB code is contained in this namespace.
namespace olb {

namespace collision {

template <concepts::DifferentiableDynamics PRIMAL>
struct Dual {
template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {

#if defined(PLATFORM_CPU_SIMD) || defined(PLATFORM_GPU_CUDA) || defined(PLAFORM_GPU_SYCL) || defined(PLATFORM_GPU_HIP)
  static_assert(true, "collision::Dual is not supported on platforms beside CPU_SISD");
#endif

  if(checkIfStaticFieldListHasMissingElement<DESCRIPTOR>()) {
    throw std::runtime_error("PRIMAL::fields meta-list does not match the list obtained from introspection.");
  }

  using Tape = util::ADtape::Tape<V>;
  using ActiveType = util::ADtape::ActiveType<Tape>;
  using FieldsAccessedByPrimalDynamics = typename PRIMAL::fields;
  auto& tape = ActiveType::getTape();

  Vector<int,DESCRIPTOR::d> size(1);
  ConcreteBlockLattice<ActiveType,DESCRIPTOR,Platform::CPU_SISD> adrLattice(size,0);
  adrLattice.setStatisticsEnabled(false);
  cpu::Cell<ActiveType,DESCRIPTOR,Platform::CPU_SISD> adFcell(adrLattice);

  // Primal populations
  adFcell.template setField<descriptors::POPULATION>(cell.template getField<opti::F>());

  // Make sure reverseAD cell includes all accessed primal fields
  FieldsAccessedByPrimalDynamics::for_each([&](auto field) {
    using FIELD = typename decltype(field)::type;
    FieldD<ActiveType,DESCRIPTOR,FIELD> value{cell.template getField<FIELD>()};
    adFcell.template setField<FIELD>(value);
  });

  auto adParams = parameters.template copyAs<ActiveType>();

  int idFstar[DESCRIPTOR::q];

  for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    tape.registerInput(adFcell[iPop]);
    idFstar[iPop] = adFcell[iPop].getIdentifier();
  }

  PRIMAL{}.collide(adFcell, adParams);

  for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    tape.registerOutput(adFcell[iPop]);
    adFcell[iPop].gradient() = cell[iPop];
  }

  tape.evaluate();

  for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = tape.gradient(idFstar[iPop]) - cell.template getFieldComponent<opti::DJDF>(iPop);
  }

  tape.reset();

  V phi2 {};
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    phi2 += cell[iPop]*cell[iPop];
  }
  return {V(1) + cell[0], phi2};
};

template <typename DESCRIPTOR>
bool checkIfStaticFieldListHasMissingElement() {
  using V = double;
  using FieldsAccessedByPrimalDynamics = typename PRIMAL::fields;
  bool staticFieldListHasMissingElement = false;

  // Check first if fields are provided by the primal dynamics
  DynamicsPromise<V,DESCRIPTOR> promise(meta::id<typename PRIMAL::template exchange_value_type<V>>{});
  auto dynamicFields = promise.accessedFields(Platform::CPU_SISD);
  std::set<std::string> staticFieldNames;
  FieldsAccessedByPrimalDynamics::for_each([&](auto field) {
    using FIELD = typename decltype(field)::type;
    staticFieldNames.emplace(FieldTypePromise<V,DESCRIPTOR>{meta::id<Array<FIELD>>{}}.name());
  });
  for (auto fields : dynamicFields) {
    if((!staticFieldNames.contains(fields.name())) && (fields.checkIfFieldIsAnArray())) {
      std::cout << "Missing: " << fields.name() << std::endl;
      staticFieldListHasMissingElement = true;
    }
  }
  return staticFieldListHasMissingElement;
};

};

}

template <typename DUAL>
concept IsDualDynamics = requires { typename DUAL::primalDynamics; };

template <concepts::DifferentiableDynamics PRIMAL,
          typename T=typename PRIMAL::value_t, typename DESCRIPTOR=typename PRIMAL::descriptor_t>
struct Dual final : public dynamics::CustomCollision<T,DESCRIPTOR,momenta::BulkTuple> {
  // TODO: Decide which momenta to expose for dual dynamics
  using MomentaF = typename momenta::BulkTuple::template type<DESCRIPTOR>;
  using EquilibriumF = typename PRIMAL::EquilibriumF;
  using parameters = typename PRIMAL::parameters;

  // expose PRIMAL dynamics type
  using primalDynamics = PRIMAL;

  template <typename NEW_T>
  using exchange_value_type = Dual<typename PRIMAL::template exchange_value_type<NEW_T>,
                                   NEW_T,DESCRIPTOR>;

  template <typename NEW_DESCRIPTOR>
  using exchange_descriptor = Dual<typename PRIMAL::template exchange_descriptor<NEW_DESCRIPTOR>,
                                   T,NEW_DESCRIPTOR>;

  template <typename M>
  using exchange_momenta = Dual<typename PRIMAL::template exchange_momenta<M>,
                                T,DESCRIPTOR>;

  std::type_index id() override {
    return typeid(Dual);
  }

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<Dual>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) {
    return collision::Dual<PRIMAL>().apply(cell, parameters);
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    typename PRIMAL::EquilibriumF().compute(cell, rho, u, fEq);
  };

  std::string getName() const override {
    return "Dual<" + PRIMAL().getName() + ">";
  };

};

template <typename T, typename DESCRIPTOR>
struct ZeroPopulationDynamics final : public dynamics::CustomCollision<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
    momenta::ZeroMomentum,
    momenta::ZeroStress,
    momenta::DefineSeparately
  >
> {
  using EquilibriumF = typename equilibria::None::template type<DESCRIPTOR,momenta::Tuple<
    momenta::BulkDensity,
    momenta::ZeroMomentum,
    momenta::ZeroStress,
    momenta::DefineSeparately
  >>;

  using parameters = meta::list<>;

  template <typename NEW_T>
  using exchange_value_type = ZeroPopulationDynamics<NEW_T,DESCRIPTOR>;

  template <typename NEW_DESCRIPTOR>
  using exchange_descriptor = ZeroPopulationDynamics<T,NEW_DESCRIPTOR>;

  template <typename M>
  using exchange_momenta = ZeroPopulationDynamics<T,DESCRIPTOR>;

  std::type_index id() override {
    return typeid(ZeroPopulationDynamics);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<ZeroPopulationDynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = 0.0;
    }
    return {-1, -1};
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    for ( int iPop = 0; iPop < DESCRIPTOR::q; iPop++ ) {
      fEq[iPop] = 0;
    }
  };

  std::string getName() const override {
    return "ZeroPopulationDynamics";
  };
};

// Manually derived discrete adjoint collision kernels

namespace collision {

/// Adjoint collision operator for Boussinesq-coupled NSE-ADE simulation
struct AdjointADEBoussinesqBGK {
  using parameters = typename meta::list<descriptors::OMEGA,
                                         NavierStokesAdvectionDiffusionCoupling::T0,
                                         NavierStokesAdvectionDiffusionCoupling::FORCE_PREFACTOR>;

  static std::string getName() {
    return "AdjointADEBoussinesqBGK";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;

    template <concepts::Cell CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      // Preparation
      auto dJdG = cell.template getFieldPointer<opti::DJDG>();
      auto dCnsedGxphi = cell.template getFieldPointer<opti::DCnseDGxphi>();
      auto dCadedGxlambda= cell.template getFieldPointer<opti::DCadeDGxlambda>();

      // Adjoint collision
      for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop){
        cell[jPop] = -dJdG[jPop] + dCnsedGxphi[jPop] + dCadedGxlambda[jPop];
      }

      // Statistics
      V norm = V();
      for(int iPop = 0; iPop < DESCRIPTOR::q; ++iPop){
        norm += cell[iPop] * cell[iPop];
      }
      return {norm};
    }
  };
};

/// Adjoint collision operator for Boussinesq-coupled NSE-ADE simulation
struct AdjointNSEBoussinesqBGK {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static std::string getName() {
    return "AdjointNSEBoussinesqBGK";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    template <concepts::Cell CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      // Preparation
      auto dCnsedFxphi = cell.template getFieldPointer<opti::DCnseDFxphi>();
      auto dJdF = cell.template getFieldPointer<opti::DJDF>();
      auto dCadedFxlambda = cell.template getFieldPointer<opti::DCadeDFxlambda>();

      // Adjoint collision
      for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop){
        cell[jPop] = -dJdF[jPop] + dCnsedFxphi[jPop] + dCadedFxlambda [jPop];
      }

      // Statistics
      V norm = V();
      for(int iPop = 0; iPop < DESCRIPTOR::q; ++iPop){
        norm += cell[iPop] * cell[iPop];
      }
      return {norm};
    }
  };
};

} // namespace olb::collision

template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using DualAdvectionDiffusionBGKDynamics = dynamics::Tuple<
  T,
  DESCRIPTOR,
  MOMENTA,
  equilibria::FirstOrder,
  collision::AdjointADEBoussinesqBGK
>;

template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using DualNavierStokesBGKDynamics = dynamics::Tuple<
  T,
  DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::AdjointNSEBoussinesqBGK
>;

template <typename T, typename DESCRIPTOR, int orientation, int direction>
using DualLocalVelocityBoundary = Dual<
  CombinedRLBdynamics<T,DESCRIPTOR,
                      dynamics::Tuple<T, DESCRIPTOR,
                                      momenta::BulkTuple,
                                      equilibria::SecondOrder,
                                      collision::RLB,
                                      dynamics::DefaultCombination>,
                      momenta::Tuple<momenta::VelocityBoundaryDensity<orientation, direction>,
                                     momenta::FixedVelocityMomentumGeneric,
                                     momenta::RegularizedBoundaryStress<orientation, direction>,
                                     momenta::DefineSeparately>>>;

template <typename T, typename DESCRIPTOR, int orientation, int direction>
using DualLocalPressureBoundary = Dual<
  CombinedRLBdynamics<T,DESCRIPTOR,
                      dynamics::Tuple<T, DESCRIPTOR,
                                      momenta::BulkTuple,
                                      equilibria::SecondOrder,
                                      collision::RLB,
                                      dynamics::DefaultCombination>,
                      momenta::Tuple<momenta::FixedDensity,
                                     momenta::FixedPressureMomentum<orientation, direction>,
                                     momenta::RegularizedBoundaryStress<orientation, direction>,
                                     momenta::DefineSeparately>>>;

template <typename T, typename DESCRIPTOR, int orientation, int direction>
using DualZouHeVelocityBoundary = Dual<
  ZouHeDynamics<T, DESCRIPTOR,
                dynamics::Tuple<T, DESCRIPTOR,
                                momenta::BulkTuple,
                                equilibria::SecondOrder,
                                collision::BGK,
                                dynamics::DefaultCombination>,
                momenta::Tuple<momenta::VelocityBoundaryDensity<orientation, direction>,
                               momenta::FixedVelocityMomentumGeneric,
                               momenta::BulkStress,
                               momenta::DefineSeparately>,
                orientation, direction>>;

template <typename T, typename DESCRIPTOR, int orientation, int direction>
using DualZouHePressureBoundary = Dual<
  ZouHeDynamics<T, DESCRIPTOR,
                dynamics::Tuple<T, DESCRIPTOR,
                                momenta::BulkTuple,
                                equilibria::SecondOrder,
                                collision::BGK,
                                dynamics::DefaultCombination>,
                momenta::Tuple<momenta::FixedDensity,
                               momenta::FixedPressureMomentum<orientation, direction>,
                               momenta::BulkStress,
                               momenta::DefineSeparately>,
                orientation, direction>>;

} // namespace olb

#endif
