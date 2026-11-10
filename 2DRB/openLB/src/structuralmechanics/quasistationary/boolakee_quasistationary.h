/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Florian Kaiser, Oliver Boolakee
 *                2026 Florian Kaiser
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

/** \file
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- header file.
 */
#ifndef BOOLAKEE_QUASISTATIONARY_H
#define BOOLAKEE_QUASISTATIONARY_H

/**
 * @file
 * @brief Quasistationary Multiple-Relaxation-Time (MRT) Lattice Boltzmann Method (LBM)
 *        implementation for linear elasticity based on the formulations by Oliver Boolakee (2023):
 *        - Periodic collision dynamics:
 *          Boolakee, O. (2023). "A new lattice Boltzmann scheme for linear elastic solids:
 *          periodic problems." Computer Methods in Applied Mechanics and Engineering, 404.
 *        - Non-periodic collision dynamics and boundary conditions:
 *          Boolakee, O. (2023). "Dirichlet and Neumann boundary conditions for a lattice
 *          Boltzmann scheme for linear elastic solids on arbitrary domains."
 *          Computer Methods in Applied Mechanics and Engineering, 415.
 *
 * Usage Overview:
 *
 * 1. Dynamics:
 *    - StationaryLinearElasticityMRT
 *      Use exclusively with periodic boundary conditions.
 *    - StationaryLinearElasticityMRTBoundary
 *      Use when applying Dirichlet or Neumann boundary conditions on arbitrary domains.
 *
 * 2. Boundary Conditions:
 *    - StationaryLinearElasticityMRTDirichletValidation
 *      Validation model: embeds the manufactured displacement solution to set populations.
 *    - StationaryLinearElasticityMRTDirichletFromDisplacement
 *      Application model: takes a given lattice displacement to set populations.
 *    - StationaryLinearElasticityMRTNeumannValidation
 *      Validation model: embeds the manufactured Cauchy stress solution to set populations.
 *    - StationaryLinearElasticityMRTNeumannValidationFixedMaterial
 *      As above, but the manufactured stress is evaluated with REFERENCE_LAME_MU and
 *      REFERENCE_BULK_MODULUS instead of the moduli of the simulated material. Needed
 *      when the material is the control of an inverse problem. Both hold lattice moduli,
 *      converted with the conversion factor of the simulated lattice.
 *    - StationaryLinearElasticityMRTNeumannFromTraction
 *      Application model: takes a given lattice traction to set populations.
 *
 * 3. Unit Conversion:
 *    - StationaryMRTUnitConverter
 *      Required unit converter for all classes in this module.
 */

namespace olb {
namespace structuralmechanics {

  // ==== Implementation of collision for quasistationary linear elasticity with periodic boundary conditions ====
  // following Boolakee, O. (2023).
  // A new lattice Boltzmann scheme for linear elastic solids: periodic problems.
  // Computer Methods in Applied Mechanics and Engineering, 404.
  template <typename T, typename DESCRIPTOR>
  struct StationaryLinearElasticityMRT final : public dynamics::CustomCollision<T, DESCRIPTOR, momenta::BulkTuple> {
    using MomentaF     = typename momenta::BulkTuple::template type<DESCRIPTOR>;
    using EquilibriumF = typename equilibria::None::template type<DESCRIPTOR, momenta::BulkTuple>;
    using parameters   = meta::list<structuralmechanics::OMEGA_SOLID>;
    // Declared so the fields are communicated/allocated for the VTK writer and the control
    using fields       = meta::list<descriptors::FORCE, opti::CONTROL_VECTOR,
                                    structuralmechanics::DISPLACEMENT, structuralmechanics::CAUCHY_STRESS,
                                    structuralmechanics::MOMENT_VECTOR, structuralmechanics::BARED_MOMENT_VECTOR>;

    template <typename M>
    using exchange_momenta = StationaryLinearElasticityMRT<T, DESCRIPTOR>;

    template <typename NEW_T>
    using exchange_value_type = StationaryLinearElasticityMRT<NEW_T, DESCRIPTOR>;

    template <typename NEW_DESCRIPTOR>
    using exchange_descriptor = StationaryLinearElasticityMRT<T, NEW_DESCRIPTOR>;

    std::type_index id() override { return typeid(StationaryLinearElasticityMRT); };

    AbstractParameters<T, DESCRIPTOR>& getParameters(BlockLattice<T, DESCRIPTOR>& block) override {
      return block.template getData<OperatorParameters<StationaryLinearElasticityMRT>>();
    }

    template <typename CELL, typename PARAMETERS, typename V=CELL::value_t>
    CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
      const V                        theta = V(1) / descriptors::invCs2<V, DESCRIPTOR>();
      const Vector<V, DESCRIPTOR::q> Omega = parameters.template get<structuralmechanics::OMEGA_SOLID>();

      const auto f1      = cell.template getField<descriptors::FORCE>();
      const auto control = cell.template getField<opti::CONTROL_VECTOR>();
      Vector<V, DESCRIPTOR::q> force = {
        control[0] * f1[0],
        control[1] * f1[1],
        V(0.), V(0.), V(0.), V(0.), V(0.), V(0.)
      };

      Vector<Vector<V, DESCRIPTOR::q>, DESCRIPTOR::q> M;
      M[0] = Vector<V, DESCRIPTOR::q>(1,  0, -1,  0, 1, -1, -1,  1);
      M[1] = Vector<V, DESCRIPTOR::q>(0,  1,  0, -1, 1,  1, -1, -1);
      M[2] = Vector<V, DESCRIPTOR::q>(0,  0,  0,  0, 1, -1,  1, -1);
      M[3] = Vector<V, DESCRIPTOR::q>(1,  1,  1,  1, 2,  2,  2,  2);
      M[4] = Vector<V, DESCRIPTOR::q>(1, -1,  1, -1, 0,  0,  0,  0);
      M[5] = Vector<V, DESCRIPTOR::q>(0,  0,  0,  0, 1, -1, -1,  1);
      M[6] = Vector<V, DESCRIPTOR::q>(0,  0,  0,  0, 1,  1, -1, -1);
      M[7] = Vector<V, DESCRIPTOR::q>(0,  0,  0,  0, 1,  1,  1,  1);

      Vector<V, DESCRIPTOR::q> cell_vec([&cell](unsigned i) { return cell[i]; });
      Vector<V, DESCRIPTOR::q> m = M * cell_vec + V(0.5) * force;

      Vector<Vector<V, DESCRIPTOR::d>, DESCRIPTOR::q> momentaFactors;
      momentaFactors[0] = Vector<V, DESCRIPTOR::d>(    1,     0);
      momentaFactors[1] = Vector<V, DESCRIPTOR::d>(    0,     1);
      momentaFactors[2] = Vector<V, DESCRIPTOR::d>(    0,     0);
      momentaFactors[3] = Vector<V, DESCRIPTOR::d>(    0,     0);
      momentaFactors[4] = Vector<V, DESCRIPTOR::d>(    0,     0);
      momentaFactors[5] = Vector<V, DESCRIPTOR::d>(theta,     0);
      momentaFactors[6] = Vector<V, DESCRIPTOR::d>(    0, theta);
      momentaFactors[7] = Vector<V, DESCRIPTOR::d>(    0,     0);

      Vector<V, 2> m_first_two = {
        m[0],
        m[1]
      };

      Vector<V, DESCRIPTOR::q> meq = momentaFactors * m_first_two;

      Vector<V, DESCRIPTOR::q> bracket{};
      for (unsigned iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
        bracket[iPop] = V(0.5) * force[iPop] + m[iPop] + Omega[iPop] * (meq[iPop] - m[iPop]);
      }

      Vector<Vector<V, DESCRIPTOR::q>, DESCRIPTOR::q> invM;
      invM[0] = Vector<V, DESCRIPTOR::q>( V(0.5), V( 0.0), V( 0.0 ), V(0.25), V( 0.25), V(-0.5 ), V( 0.0 ), V(-0.5 ));
      invM[1] = Vector<V, DESCRIPTOR::q>( V(0.0), V( 0.5), V( 0.0 ), V(0.25), V(-0.25), V( 0.0 ), V(-0.5 ), V(-0.5 ));
      invM[2] = Vector<V, DESCRIPTOR::q>(V(-0.5), V( 0.0), V( 0.0 ), V(0.25), V( 0.25), V( 0.5 ), V( 0.0 ), V(-0.5 ));
      invM[3] = Vector<V, DESCRIPTOR::q>( V(0.0), V(-0.5), V( 0.0 ), V(0.25), V(-0.25), V( 0.0 ), V( 0.5 ), V(-0.5 ));
      invM[4] = Vector<V, DESCRIPTOR::q>( V(0.0), V( 0.0), V( 0.25), V(0.0 ), V( 0.0 ), V( 0.25), V( 0.25), V( 0.25));
      invM[5] = Vector<V, DESCRIPTOR::q>( V(0.0), V( 0.0), V(-0.25), V(0.0 ), V( 0.0 ), V(-0.25), V( 0.25), V( 0.25));
      invM[6] = Vector<V, DESCRIPTOR::q>( V(0.0), V( 0.0), V( 0.25), V(0.0 ), V( 0.0 ), V(-0.25), V(-0.25), V( 0.25));
      invM[7] = Vector<V, DESCRIPTOR::q>( V(0.0), V( 0.0), V(-0.25), V(0.0 ), V( 0.0 ), V( 0.25), V(-0.25), V( 0.25));

      for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] = invM[iPop] * bracket;
      }

      Vector<V, DESCRIPTOR::q> m_bared = V(0.5) * (m + bracket);

      cell.template setField<structuralmechanics::BARED_MOMENT_VECTOR>(m_bared);
      cell.template setField<structuralmechanics::MOMENT_VECTOR>(bracket);

      Vector<V, DESCRIPTOR::d> u = {
        m_bared[0],
        m_bared[1]
      };
      cell.template setField<structuralmechanics::DISPLACEMENT>(u);

      Vector<V, DESCRIPTOR::d * DESCRIPTOR::d> stress = {
        V(-0.5) * (m_bared[3] + m_bared[4]),
        -m_bared[2],
        -m_bared[2],
        V(-0.5) * (m_bared[3] - m_bared[4])
      };
      cell.template setField<structuralmechanics::CAUCHY_STRESS>(stress);

      return {V(1), util::normSqr<V, DESCRIPTOR::d>(u)};
    };

    void computeEquilibrium(ConstCell<T, DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d],
                            T fEq[DESCRIPTOR::q]) const override {
      EquilibriumF().compute(cell, rho, u, fEq);
    }

    std::string getName() const override { return "StationaryLinearElasticityMRT<" + MomentaF().getName() + ">"; };
  };

  // ==== Implementation of collision for quasistationary linear elasticity with boundary conditions ====
  // following Boolakee, O. (2023).
  // Dirichlet and Neumann boundary conditions for a lattice Boltzmann scheme for linear elastic solids on arbitrary domains.
  // Computer Methods in Applied Mechanics and Engineering, 415.
  template <typename T, typename DESCRIPTOR>
  struct StationaryLinearElasticityMRTBoundary final : public dynamics::CustomCollision<T, DESCRIPTOR, momenta::BulkTuple> {
    using MomentaF     = typename momenta::BulkTuple::template type<DESCRIPTOR>;
    using EquilibriumF = typename equilibria::None::template type<DESCRIPTOR, momenta::BulkTuple>;
    using parameters   = meta::list<structuralmechanics::OMEGA_SOLID>;

    template <typename M>
    using exchange_momenta = StationaryLinearElasticityMRTBoundary<T, DESCRIPTOR>;

    template <typename NEW_T>
    using exchange_value_type = StationaryLinearElasticityMRTBoundary<NEW_T, DESCRIPTOR>;

    template <typename NEW_DESCRIPTOR>
    using exchange_descriptor = StationaryLinearElasticityMRTBoundary<T, NEW_DESCRIPTOR>;

    std::type_index id() override { return typeid(StationaryLinearElasticityMRTBoundary); };

    AbstractParameters<T, DESCRIPTOR>& getParameters(BlockLattice<T, DESCRIPTOR>& block) override {
      return block.template getData<OperatorParameters<StationaryLinearElasticityMRTBoundary>>();
    }

    template <typename CELL, typename PARAMETERS, typename V = typename CELL::value_t>
    CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
      const V                        theta    = V(1) / descriptors::invCs2<V, DESCRIPTOR>();
      const Vector<V, DESCRIPTOR::q> Omega    = parameters.template get<structuralmechanics::OMEGA_SOLID>();
      const V                        omega_s  = Omega[3];

      V tau_s = V(1.) / omega_s - V(0.5);
      V tau_f = V(0.5);

      V gamma = theta * tau_f / ((V(1.) + theta) * (tau_s - tau_f));

      Vector<Vector<V, DESCRIPTOR::q>, DESCRIPTOR::q> M;
      M[0] = Vector<V, DESCRIPTOR::q>(    1,     0,    -1,     0,                       1,                      -1,                      -1,  1 );
      M[1] = Vector<V, DESCRIPTOR::q>(    0,     1,     0,    -1,                       1,                       1,                      -1, -1 );
      M[2] = Vector<V, DESCRIPTOR::q>(    0,     0,     0,     0,                       1,                      -1,                       1, -1 );
      M[3] = Vector<V, DESCRIPTOR::q>(    1,     1,     1,     1,                       2,                       2,                       2,  2 );
      M[4] = Vector<V, DESCRIPTOR::q>(    1,    -1,     1,    -1,                       0,                       0,                       0,  0 );
      M[5] = Vector<V, DESCRIPTOR::q>(    0,     0,     0,     0,                       1,                      -1,                      -1,  1 );
      M[6] = Vector<V, DESCRIPTOR::q>(    0,     0,     0,     0,                       1,                       1,                      -1, -1 );
      M[7] = Vector<V, DESCRIPTOR::q>(gamma, gamma, gamma, gamma, (V(1.) + V(2.) * gamma), (V(1.) + V(2.) * gamma), (V(1.) + V(2.) * gamma), (V(1.) + V(2.) * gamma) );

      const Vector<V, DESCRIPTOR::d> f1    = cell.template getField<descriptors::FORCE>();
      const Vector<V, DESCRIPTOR::q> force = {
        f1[0],
        f1[1],
        V(0.), V(0.), V(0.), V(0.), V(0.), V(0.)
      };

      Vector<V, DESCRIPTOR::q> cell_vec([&cell](unsigned i) { return cell[i]; });
      Vector<V, DESCRIPTOR::q> m = M * cell_vec + V(0.5) * force;

      Vector<Vector<V, DESCRIPTOR::d>, DESCRIPTOR::q> momentaFactors;
      momentaFactors[0] = Vector<V, DESCRIPTOR::d>(    1,     0);
      momentaFactors[1] = Vector<V, DESCRIPTOR::d>(    0,     1);
      momentaFactors[2] = Vector<V, DESCRIPTOR::d>(    0,     0);
      momentaFactors[3] = Vector<V, DESCRIPTOR::d>(    0,     0);
      momentaFactors[4] = Vector<V, DESCRIPTOR::d>(    0,     0);
      momentaFactors[5] = Vector<V, DESCRIPTOR::d>(theta,     0);
      momentaFactors[6] = Vector<V, DESCRIPTOR::d>(    0, theta);
      momentaFactors[7] = Vector<V, DESCRIPTOR::d>(    0,     0);

      Vector<V, 2> m_first_two = {
        m[0],
        m[1]
      };

      Vector<V, DESCRIPTOR::q> meq = momentaFactors * m_first_two;

      Vector<V, DESCRIPTOR::q> bracket;
      for (unsigned iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
        bracket[iPop] = V(0.5) * force[iPop] + m[iPop] + Omega[iPop] * (meq[iPop] - m[iPop]);
      }

      Vector<V, DESCRIPTOR::q> m_bared = V(0.5) * (m + bracket);

      Vector<Vector<V, DESCRIPTOR::q>, DESCRIPTOR::q> invM;
      invM[0] = Vector<V, DESCRIPTOR::q>( V( 0.5), V( 0.0), V( 0.0 ), V( 0.5 ) * gamma + V(0.25), V( 0.25), V(-0.5 ), V( 0.0 ), V(-0.5));
      invM[1] = Vector<V, DESCRIPTOR::q>( V( 0.0), V( 0.5), V( 0.0 ), V( 0.5 ) * gamma + V(0.25), V(-0.25), V( 0.0 ), V(-0.5 ), V(-0.5));
      invM[2] = Vector<V, DESCRIPTOR::q>( V(-0.5), V( 0.0), V( 0.0 ), V( 0.5 ) * gamma + V(0.25), V( 0.25), V( 0.5 ), V( 0.0 ), V(-0.5));
      invM[3] = Vector<V, DESCRIPTOR::q>( V( 0.0), V(-0.5), V( 0.0 ), V( 0.5 ) * gamma + V(0.25), V(-0.25), V( 0.0 ), V( 0.5 ), V(-0.5));
      invM[4] = Vector<V, DESCRIPTOR::q>( V( 0.0), V( 0.0), V( 0.25), V(-0.25) * gamma,           V( 0.0 ), V( 0.25), V( 0.25), V(0.25));
      invM[5] = Vector<V, DESCRIPTOR::q>( V( 0.0), V( 0.0), V(-0.25), V(-0.25) * gamma,           V( 0.0 ), V(-0.25), V( 0.25), V(0.25));
      invM[6] = Vector<V, DESCRIPTOR::q>( V( 0.0), V( 0.0), V( 0.25), V(-0.25) * gamma,           V( 0.0 ), V(-0.25), V(-0.25), V(0.25));
      invM[7] = Vector<V, DESCRIPTOR::q>( V( 0.0), V( 0.0), V(-0.25), V(-0.25) * gamma,           V( 0.0 ), V( 0.25), V(-0.25), V(0.25));

      for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] = invM[iPop] * bracket;
      }

      auto previous_cell = cell.template getField<structuralmechanics::CURRENT_CELL>();
      cell.template setField<structuralmechanics::CURRENT_CELL>(
        Vector<V, DESCRIPTOR::q>(cell[0], cell[1], cell[2], cell[3], cell[4], cell[5], cell[6], cell[7])
      );
      cell.template setField<structuralmechanics::PREVIOUS_CELL>(
        Vector<V, DESCRIPTOR::q>(
          previous_cell[0],
          previous_cell[1],
          previous_cell[2],
          previous_cell[3],
          previous_cell[4],
          previous_cell[5],
          previous_cell[6],
          previous_cell[7]
        )
      );

      cell.template setField<structuralmechanics::BARED_MOMENT_VECTOR>(m_bared);
      cell.template setField<structuralmechanics::MOMENT_VECTOR>(bracket);
      Vector<V, DESCRIPTOR::d> u = {
        m_bared[0],
        m_bared[1]
      };
      cell.template setField<structuralmechanics::DISPLACEMENT>(u);

      Vector<V, DESCRIPTOR::d * DESCRIPTOR::d> stress = {
        V(-0.5) * (m_bared[3] + m_bared[4]),
        -m_bared[2],
        -m_bared[2],
        V(-0.5) * (m_bared[3] - m_bared[4])
      };
      cell.template setField<structuralmechanics::CAUCHY_STRESS>(stress);

      return {V(1), util::normSqr<V, DESCRIPTOR::d>(u)};
    };

    void computeEquilibrium(ConstCell<T, DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d],
                            T fEq[DESCRIPTOR::q]) const override {
      EquilibriumF().compute(cell, rho, u, fEq);
    }

    std::string getName() const override { return "StationaryLinearElasticityMRTBoundary<" + MomentaF().getName() + ">"; };
  };

  struct StationaryLinearElasticityManufacturedDisplacement {

    using parameters = meta::list<
      structuralmechanics::CHAR_U
    >;

    template <typename DESCRIPTOR, typename V, typename PARAMETERS>
    static Vector<V, DESCRIPTOR::d> displacement(V x, V y, PARAMETERS& parameters) any_platform {
      const V convU = parameters.template get<structuralmechanics::CHAR_U>();

      Vector<V, DESCRIPTOR::d> uD = {
        (V(9) / V(10000)) * util::sin(V(2.) * V(M_PI) * x * y),
        (V(7) / V(10000)) * util::cos(V(2.) * V(M_PI) * y) * (x * x + V(1.))
      };
      uD /= convU;

      return uD;
    }
  };


  // ==== Implementation of the Dirichlet boundary condition for quasistationary linear elasticity ====
  // following Boolakee, O. (2023).
  // Dirichlet and Neumann boundary conditions for a lattice Boltzmann scheme for linear elastic solids on arbitrary domains.
  // Computer Methods in Applied Mechanics and Engineering, 415.
  template<typename DISPLACEMENT>
  struct StationaryLinearElasticityMRTDirichletFromDisplacement {
    static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

    int getPriority() const {
      return -1;
    }

    using parameters = concat_parameter_lists<
      meta::list<
        structuralmechanics::LAME_MU,
        structuralmechanics::BULK_MODULUS
      >,
      typename DISPLACEMENT::parameters
    >::type;

    template <typename CELL, typename PARAMETERS,
              typename DESCRIPTOR = CELL::descriptor_t,
              typename V = typename CELL::value_t>
    void apply(CELL& cell, PARAMETERS& parameters) any_platform {

      const V theta     = V(1) / descriptors::invCs2<V, DESCRIPTOR>();
      const V mu        = parameters.template get<structuralmechanics::LAME_MU>();
      const V bulk      = parameters.template get<structuralmechanics::BULK_MODULUS>();

      for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
        auto c        = descriptors::c<DESCRIPTOR>(iPop);
        auto neighbor = cell.neighbor(c);
        int oppo      = descriptors::opposite<DESCRIPTOR>(iPop);

        V q_ij        = neighbor.template getField<structuralmechanics::SOLID_DISTANCE_FIELD>()[oppo];
        if (!(neighbor.template getField<structuralmechanics::SOLID_DISTANCE_FIELD>()[oppo] > V(0))) { continue; }

        int i         = c[0];
        int j         = c[1];

        auto prev_cell = neighbor.template getField<structuralmechanics::PREVIOUS_CELL>();

        V x = neighbor.template getField<structuralmechanics::BOUNDARY_COORDS_X>()[oppo];
        V y = neighbor.template getField<structuralmechanics::BOUNDARY_COORDS_Y>()[oppo];

        Vector<V, DESCRIPTOR::d> u_D           = DISPLACEMENT::template displacement<DESCRIPTOR>(x, y, parameters);

        Vector<V, DESCRIPTOR::q> bared_moments = neighbor.template getField<structuralmechanics::BARED_MOMENT_VECTOR>();

        V dx_ux = V(-.25) * ((V(1.) / bulk) * bared_moments[3] + (V(1.) / mu) * bared_moments[4]);
        V dx_uy = V(-.25) *   V(2.) / mu    * bared_moments[2];

        V dy_ux = V(-.25) *   V(2.) / mu    * bared_moments[2];
        V dy_uy = V(-.25) * ((V(1.) / bulk) * bared_moments[3] - (V(1.) / mu) * bared_moments[4]);

        bool mask_ij = (i == 0 || j == 0);
        V w = mask_ij * (V(1.) - theta) + (1 - mask_ij) * (theta * V(0.5));

        // Calculate source term according to eq. 103, 104, 108, 109
        V term_v1 = std::abs(i) * dx_ux + std::abs(j) * dy_uy;
        V term_v2 = dx_ux + dy_uy + (i * j) * (dx_uy + dy_ux);
        V first_order_term = mask_ij * term_v1 + (1 - mask_ij) * term_v2;

        V s_ij = w * (i * u_D[0] + j * u_D[1])
               + w * (q_ij - V(0.5)) * first_order_term;

        cell[iPop] = prev_cell[oppo] + s_ij;
      }
    };
  };

  using StationaryLinearElasticityMRTDirichletValidation = StationaryLinearElasticityMRTDirichletFromDisplacement<StationaryLinearElasticityManufacturedDisplacement>;

  // The moduli of the manufactured stress are read from the parameters MU and BULK and
  // have to be lattice moduli.
  template <typename MU, typename BULK>
  struct TractionFromManufacturedStress {

    using parameters = meta::list<MU,BULK>;

    template <typename DESCRIPTOR, typename V, typename PARAMETERS>
    static Vector<V, DESCRIPTOR::d * DESCRIPTOR::d> cauchyStress(V x, V y, PARAMETERS& parameters) any_platform
    {
      const V bulk      = parameters.template get<BULK>();
      const V mu        = parameters.template get<MU>();


      const V pi = V(M_PI);

      const V A = (V(9) / V(5000)) * pi * y * util::cos(V(2) * pi * x * y);            // eps_xx
      const V B = (V(7) / V(5000)) * pi * (x * x + V(1)) * util::sin(V(2) * pi * y);   // -eps_yy
      const V S = mu * ((V(7) / V(5000)) * x * util::cos(V(2) * pi * y)
                      + (V(9) / V(5000)) * pi * x * util::cos(V(2) * pi * x * y));     // sigma_xy

      Vector<V, DESCRIPTOR::d * DESCRIPTOR::d> sigma = {
        (bulk + mu) * A - (bulk - mu) * B,
        S,
        S,
        (bulk - mu) * A - (bulk + mu) * B
      };

      return sigma;
    }

    template <typename DESCRIPTOR, typename V, typename PARAMETERS>
    static Vector<V, DESCRIPTOR::d> traction(V x, V y, const Vector<V, DESCRIPTOR::d>& n,
                                             PARAMETERS& parameters) any_platform
    {
      constexpr int d = DESCRIPTOR::d;

      const Vector<V, d * d> sigma = cauchyStress<DESCRIPTOR>(x, y, parameters);

      Vector<V, d> TN = {
        sigma[0] * n[0] + sigma[1] * n[1],
        sigma[2] * n[0] + sigma[3] * n[1]
      };

      return TN;
    }
  };

  /// Manufactured traction built from the moduli of the simulated material
  using StationaryLinearElasticityManufacturedTraction = TractionFromManufacturedStress<
    structuralmechanics::LAME_MU, structuralmechanics::BULK_MODULUS>;

  /// Manufactured traction built from moduli that are independent of the simulated material,
  /// for inverse problems in which the material itself is the control
  using StationaryLinearElasticityManufacturedTractionFixedMaterial = TractionFromManufacturedStress<
    structuralmechanics::REFERENCE_LAME_MU, structuralmechanics::REFERENCE_BULK_MODULUS>;


  // ==== Implementation of the Neumann boundary condition for quasistationary linear elasticity ====
  // following Boolakee, O. (2023).
  // Dirichlet and Neumann boundary conditions for a lattice Boltzmann scheme for linear elastic solids on arbitrary domains.
  // Computer Methods in Applied Mechanics and Engineering, 415.
  template<typename TRACTION>
  struct StationaryLinearElasticityMRTNeumannFromTraction {
    static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

    using parameters = typename concat_parameter_lists<
      meta::list<
        structuralmechanics::NEUMANN_SOLID_C
      >,
      typename TRACTION::parameters
    >::type;

    int getPriority() const {
      return 0;
    }

    template <typename CELL, typename PARAMETERS,
              typename DESCRIPTOR = typename CELL::descriptor_t,
              typename V = typename CELL::value_t>
    void apply(CELL& cell, PARAMETERS& parameters) any_platform {

      Vector<V, 3> constants = parameters.template get<structuralmechanics::NEUMANN_SOLID_C>();

      for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
        auto c        = descriptors::c<DESCRIPTOR>(iPop);
        int oppo      = descriptors::opposite<DESCRIPTOR>(iPop);
        auto neighbor = cell.neighbor(c);

        if (!(neighbor.template getField<structuralmechanics::SOLID_DISTANCE_FIELD>()[oppo] > V(0))) { continue; }

        int i = c[0];
        int j = c[1];

        V x = neighbor.template getField<structuralmechanics::BOUNDARY_COORDS_X >()[oppo];
        V y = neighbor.template getField<structuralmechanics::BOUNDARY_COORDS_Y >()[oppo];
        Vector<V, DESCRIPTOR::d> n = {
          neighbor.template getField<structuralmechanics::NEUMANN_NORMAL_X  >()[oppo],
          neighbor.template getField<structuralmechanics::NEUMANN_NORMAL_Y  >()[oppo]
        };

        bool mask_zeta = (util::abs(n[1]) > util::abs(n[0]));
        int zeta = mask_zeta * V(1.) + (1 - mask_zeta) * V(-1.);

        Vector<V, DESCRIPTOR::d> traction = TRACTION::template traction<DESCRIPTOR>(x, y, n, parameters);

        V sum{};
        for (int k = -1; k <= 1; k += 1) {
          for (int l = -1; l <= 1; l += 1) {
            bool mask_kl = (k != 0 || l != 0);

            int dim = findDim<DESCRIPTOR>(k, l);

            int safe_dim = (mask_kl > 0) ? dim : 0;
            V weight = calcWeight<V, DESCRIPTOR>(i, j, k, l, zeta, n, constants);

            sum += mask_kl * (weight * neighbor[safe_dim]);
          }
        }

        bool mask_ij = (i == 0 || j == 0);
        V s_ij_v1 = firstOrderSourceNu1<V, DESCRIPTOR>(i, j, traction);
        V s_ij_v2 = firstOrderSourceNu2<V, DESCRIPTOR>(i, j, traction, zeta);

        V s_ij = mask_ij * s_ij_v1 + (V(1.) - mask_ij) * s_ij_v2;

        cell[iPop] = sum + s_ij;
      }
    };

    template <typename V, typename DESCRIPTOR>
    V calcWeight(const int i, const int j, const int k, const int l, const int zeta, const Vector<V, DESCRIPTOR::d>& n, Vector<V, 3>& constants) any_platform {
      bool theta_ij = (i == -k && j == -l);

      bool mask_v1 = (i == 0 || j == 0);
      bool mask_v2 = 1 - mask_v1;

      int absK = util::abs(k);
      int absL = util::abs(l);
      int kl = k * l;

      V weight_v1 = -absK * absL * (V(1.0) + i * n[0] + j * n[1]) * constants[0]
                    - kl * (i * n[1] + j * n[0]) * constants[1]
                    - (absK * (1 - absL) * (util::abs(V(i)) + i * n[0])
                    + absL * (1 - absK) * (util::abs(V(j)) + j * n[1])) * constants[2]
                    - theta_ij;

      V z_plus  = (V(1.0) + V(zeta)) * V(0.5);
      V z_minus = (V(1.0) - V(zeta)) * V(0.5);

      V weight_v2 = -V(0.5) * absK * absL * (z_plus * i * n[0] + z_minus * j * n[1]) * constants[0]
                    - V(0.5) * kl * (V(i * j) + z_plus * i * n[1] + z_minus * j * n[0]) * constants[1]
                    - V(0.5) * (absK * (1 - absL) * z_plus * i * n[0]
                    + absL * (1 - absK) * z_minus * j * n[1]) * constants[2]
                    - theta_ij;

      return mask_v1 * weight_v1 + mask_v2 * weight_v2;
    }

    // Eq. 147
    template <typename V, typename DESCRIPTOR>
    V firstOrderSourceNu1(const int i, const int j, const Vector<V, DESCRIPTOR::d>& Traction) any_platform {
      return (i * Traction[0] + j * Traction[1]);
    };


    // Eq. 148
    template <typename V, typename DESCRIPTOR>
    V firstOrderSourceNu2(const int i, const int j, const Vector<V, DESCRIPTOR::d>& Traction, const int zeta) any_platform {
      return V(0.25) * (i * Traction[0] * (1 + zeta) + j * Traction[1] * (1 - zeta));
    };


    template <typename DESCRIPTOR>
    int findDim(const int i, const int j) any_platform {
      int found_dim = -1;

      for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
        Vector c = descriptors::c<DESCRIPTOR>(iPop);

        bool mask = (i == c[0] && j == c[1]);

        found_dim = (mask * iPop) + ((1 - mask) * found_dim);
      }
      return found_dim;
    }
  };

  using StationaryLinearElasticityMRTNeumannValidation = StationaryLinearElasticityMRTNeumannFromTraction<StationaryLinearElasticityManufacturedTraction>;
  using StationaryLinearElasticityMRTNeumannValidationFixedMaterial = StationaryLinearElasticityMRTNeumannFromTraction<StationaryLinearElasticityManufacturedTractionFixedMaterial>;
} // namespace structuralmechanics
} // namespace olb

#endif
