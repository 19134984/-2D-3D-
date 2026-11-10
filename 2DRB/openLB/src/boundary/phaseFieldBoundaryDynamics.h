/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Tim Bingert, Luiz Eduardo Czelusniak
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

#ifndef PHASE_FIELD_BOUNDARY_DYNAMICS_H
#define PHASE_FIELD_BOUNDARY_DYNAMICS_H

#include "dynamics/dynamics.h"

namespace olb {

/**
* Implementation of regularized Dirichlet boundary condition for the order parameter.
*/
template <typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA>
class CombinedPhaseFieldRLBdynamics final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
private:
  // Use same MOMENTA in combined and nested (boundary) dynamics
  using CORRECTED_DYNAMICS = typename DYNAMICS::template exchange_momenta<MOMENTA>;

public:
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename DYNAMICS::EquilibriumF;

  using parameters = typename DYNAMICS::parameters;

  template <typename NEW_T>
  using exchange_value_type = CombinedPhaseFieldRLBdynamics<
    NEW_T, DESCRIPTOR,
    typename DYNAMICS::template exchange_value_type<NEW_T>,
    MOMENTA
  >;

  template <typename NEW_DESCRIPTOR>
  using exchange_descriptor = CombinedPhaseFieldRLBdynamics<
    T, NEW_DESCRIPTOR,
    typename DYNAMICS::template exchange_descriptor<NEW_DESCRIPTOR>,
    MOMENTA
  >;

  template <typename M>
  using exchange_momenta = CombinedPhaseFieldRLBdynamics<T,DESCRIPTOR,DYNAMICS,M>;

  std::type_index id() override {
    return typeid(CombinedPhaseFieldRLBdynamics);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<CombinedPhaseFieldRLBdynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    V jNeq[DESCRIPTOR::d] { };

    const V phi = MomentaF().computeRho(cell);
    const auto u = cell.template getField<descriptors::VELOCITY>();
    MomentaF().computeJ(cell, jNeq);

    V fEq[DESCRIPTOR::q] { };
    EquilibriumF().compute(cell, phi, u, fEq);
    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = fEq[iPop] + equilibrium<DESCRIPTOR>::template fromJneqToFneq<V>(iPop, jNeq);
    }

    return typename CORRECTED_DYNAMICS::CollisionO().apply(cell, parameters);
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    for ( int iPop = 0; iPop < DESCRIPTOR::q; iPop++ ) {
      fEq[iPop] = equilibrium<DESCRIPTOR>::firstOrder(iPop, rho, u);
    }
  };

  std::string getName() const override {
    return "CombinedPhaseFieldRLBdynamics<" + MomentaF().getName() + ">";
  };

};

/**
* Implementation of convective boundary condition for the order parameter.
*/
template<typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA, int direction, int orientation>
struct PhaseFieldConvectiveOutletDynamics : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename DYNAMICS::EquilibriumF;

  //using parameters = typename DYNAMICS::parameters;
  using parameters = meta::list<descriptors::MAX_VELOCITY>;

  std::type_index id() override {
    return typeid(PhaseFieldConvectiveOutletDynamics);
  }

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<PhaseFieldConvectiveOutletDynamics>>();
  }

  constexpr static bool is_vectorizable = false;

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    V phi = MomentaF().computeRho(cell);
    V u[DESCRIPTOR::d];
    MomentaF().computeU(cell, u);
    V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
    return {phi, uSqr};
  }

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    EquilibriumF().compute(cell, rho, u, fEq);
  };

  std::string getName() const override {
    return "PhaseFieldConvectiveOutletDynamics<" + DYNAMICS().getName() + ">";
  };

};

/// Regularized incompressible (multiphase) pops followed by any other Dynamics
/**
 * Commonly used to model boundaries consisting of a specific boundary
 * momenta and standard dynamics::Tuple-declared dynamics.
 **/
template <typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA>
class CombinedIRLBdynamics final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
private:
  // Use same MOMENTA in combined and nested (boundary) dynamics
  using CORRECTED_DYNAMICS = typename DYNAMICS::template exchange_momenta<MOMENTA>;

public:
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename CORRECTED_DYNAMICS::EquilibriumF;

  using parameters = typename CORRECTED_DYNAMICS::parameters;

  template <typename NEW_T>
  using exchange_value_type = CombinedIRLBdynamics<
    NEW_T, DESCRIPTOR,
    typename DYNAMICS::template exchange_value_type<NEW_T>,
    MOMENTA
  >;

  template <typename NEW_DESCRIPTOR>
  using exchange_descriptor = CombinedIRLBdynamics<
    T, NEW_DESCRIPTOR,
    typename DYNAMICS::template exchange_descriptor<NEW_DESCRIPTOR>,
    MOMENTA
  >;

  template <typename M>
  using exchange_momenta = CombinedIRLBdynamics<T,DESCRIPTOR,DYNAMICS,M>;

  std::type_index id() override {
    return typeid(CombinedIRLBdynamics);
  }

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<CombinedIRLBdynamics>>();
  }

  template <concepts::Cell CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    V p, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR>::n];
    MomentaF().computeAllMomenta(cell,p,u,pi);
    V fEq[DESCRIPTOR::q] { };
    EquilibriumF().compute(cell, p, u, fEq);
    auto force = cell.template getField<descriptors::FORCE>();
    auto rho = cell.template getField<descriptors::RHO>();
    auto gradRho = cell.template getField<descriptors::NABLARHO>();
    auto laplaceRho = cell.template getField<descriptors::LAPLACERHO>();

    V u_gradRho = 0;
    for (int iD=0; iD < DESCRIPTOR::d; iD++) {
      u_gradRho += u[iD]*gradRho[iD];
    }
    for (int iAlpha=0; iAlpha < DESCRIPTOR::d; ++iAlpha) {
      for (int iBeta=iAlpha; iBeta < DESCRIPTOR::d; ++iBeta) {
        pi[ static_cast<int>( iAlpha * (DESCRIPTOR::d - ( iAlpha - 1 ) * 0.5) ) ] += V{0.5}/descriptors::invCs2<V,DESCRIPTOR>()*u_gradRho*(iAlpha == iBeta);
      }
    }
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = fEq[iPop];
      V c_F = 0;
      for (int iD=0; iD < DESCRIPTOR::d; iD++) {
        c_F += descriptors::c<DESCRIPTOR>(iPop,iD)*force[iD]*(rho + laplaceRho/4.);
      }
      V fNeq = equilibrium<DESCRIPTOR>::template fromPiToFneq<V>(iPop, pi) -
                descriptors::t<V,DESCRIPTOR>(iPop) * V{0.5} * (u_gradRho +
                descriptors::invCs2<V,DESCRIPTOR>()*c_F);
      cell[iPop] += fNeq;
    }

    return CORRECTED_DYNAMICS().collide(cell, parameters);
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    EquilibriumF().compute(cell, rho, u, fEq);
  };

  std::string getName() const override {
    return "CombinedIRLBdynamics<" + CORRECTED_DYNAMICS().getName() + ">";
  };

};


/// Equilibrium incompressible (multiphase) pops followed by any other Dynamics
/**
 * Used for geometrically complicated boundaries where no regularized formula is available.
 **/
template <typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA>
class CombinedEquilibriumDynamics final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
private:
  // Use same MOMENTA in combined and nested (boundary) dynamics
  using CORRECTED_DYNAMICS = typename DYNAMICS::template exchange_momenta<MOMENTA>;

public:
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename CORRECTED_DYNAMICS::EquilibriumF;

  using parameters = typename CORRECTED_DYNAMICS::parameters;

  template <typename NEW_T>
  using exchange_value_type = CombinedEquilibriumDynamics<
    NEW_T, DESCRIPTOR,
    typename DYNAMICS::template exchange_value_type<NEW_T>,
    MOMENTA
  >;

  template <typename NEW_DESCRIPTOR>
  using exchange_descriptor = CombinedEquilibriumDynamics<
    T, NEW_DESCRIPTOR,
    typename DYNAMICS::template exchange_descriptor<NEW_DESCRIPTOR>,
    MOMENTA
  >;

  template <typename M>
  using exchange_momenta = CombinedEquilibriumDynamics<T,DESCRIPTOR,DYNAMICS,M>;

  std::type_index id() override {
    return typeid(CombinedEquilibriumDynamics);
  }

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<OperatorParameters<CombinedEquilibriumDynamics>>();
  }

  template <concepts::Cell CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
    V p, u[DESCRIPTOR::d];
    MomentaF().computeRhoU(cell,p,u);
    V fEq[DESCRIPTOR::q] { };
    EquilibriumF().compute(cell, p, u, fEq);

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = fEq[iPop];
    }

    return CORRECTED_DYNAMICS().collide(cell, parameters);
  };

  void computeEquilibrium(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T fEq[DESCRIPTOR::q]) const override {
    EquilibriumF().compute(cell, rho, u, fEq);
  };

  std::string getName() const override {
    return "CombinedEquilibriumDynamics<" + CORRECTED_DYNAMICS().getName() + ">";
  };

};

template <typename T, typename DESCRIPTOR>
using NoCollideDynamicsExternalVelocity = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::ExternalVelocityTuple,
  equilibria::FirstOrder,
  collision::None
>;

}

#endif
