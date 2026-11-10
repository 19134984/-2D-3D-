/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2025 Florian Kaiser, Oliver Boolakee
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
#ifndef BOOLAKEE_ELASTODYNAMICS_H
#define BOOLAKEE_ELASTODYNAMICS_H

/**
 * @file
 * @brief Linear Elastodynamics Vector Lattice Boltzmann Method (LBM) implementation
 *        based on the formulations by Oliver Boolakee et al. (2025):
 *        - Collision dynamics and Dirichlet conditions:
 *          Boolakee, O. et al. (2025). "Lattice Boltzmann for linear elastodynamics:
 *          Periodic problems and Dirichlet boundary conditions."
 *          Computer Methods in Applied Mechanics and Engineering, 433.
 *        - Neumann boundary conditions:
 *          Boolakee, O. (2025). "Lattice Boltzmann for Solids."
 *          Doctoral Thesis, ETH Zürich. https://doi.org/10.3929/ethz-b-000742210
 *
 * Usage Overview:
 *
 * 1. Dynamics:
 *    - ElastodynamicsVectorLBM
 *      Core dynamics solver implementing collision for linear elastodynamics.
 *
 * 2. Boundary Conditions:
 *    - ElastodynamicsVectorLBMDirichletValidation
 *      Validation model: embeds the manufactured velocity solution to set populations.
 *    - ElastodynamicsVectorLBMDirichletFromVelocity
 *      Application model: takes a given lattice velocity to set populations.
 *    - ElastodynamicsVectorLBMNeumannValidation
 *      Validation model: embeds the manufactured Cauchy stress solution to set populations.
 *    - ElastodynamicsVectorLBMNeumannFromTraction
 *      Application model: takes a given lattice traction to set populations.
 *
 * 3. Unit Conversion:
 *    - ElastodynamicsVectorLBMUnitConverter
 *      Required unit converter for this module.
 */

namespace olb {
namespace structuralmechanics {

  // ==== Implementation of collision for linear elastodynamics ====
  // following Boolakee, O. et al. (2025).
  // Lattice Boltzmann for linear elastodynamics: Periodic problems and Dirichlet boundary conditions
  // Computer Methods in Applied Mechanics and Engineering, 433.
  template <typename T, typename DESCRIPTOR>
  struct ElastodynamicsVectorLBM final : public dynamics::CustomCollision<T, DESCRIPTOR, momenta::BulkTuple> {
    using MomentaF     = typename momenta::BulkTuple::template type<DESCRIPTOR>;
    using EquilibriumF = typename equilibria::None::template type<DESCRIPTOR, momenta::BulkTuple>;

    using parameters = meta::list<
      descriptors::DT,
      descriptors::OMEGA,
      structuralmechanics::C_MU,
      structuralmechanics::C_K,
      structuralmechanics::LATTICE_SPEED
    >;

    template <typename M>
    using exchange_momenta = ElastodynamicsVectorLBM<T, DESCRIPTOR>;

    template <typename NEW_T>
    using exchange_value_type = ElastodynamicsVectorLBM<NEW_T, DESCRIPTOR>;

    template <typename NEW_DESCRIPTOR>
    using exchange_descriptor = ElastodynamicsVectorLBM<T,NEW_DESCRIPTOR>;

    std::type_index id() override { return typeid(ElastodynamicsVectorLBM); };

    AbstractParameters<T, DESCRIPTOR>& getParameters(BlockLattice<T, DESCRIPTOR>& block) override {
      return block.template getData<OperatorParameters<ElastodynamicsVectorLBM>>();
    }

    template <typename CELL, typename PARAMETERS, typename V=CELL::value_t>
    CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {

      VectorCell cellV(cell);

      const V omega    = parameters.template get<descriptors::OMEGA                 >();
      const V c_K      = parameters.template get<structuralmechanics::C_K           >();
      const V c_mu     = parameters.template get<structuralmechanics::C_MU          >();
      const V cLattice = parameters.template get<structuralmechanics::LATTICE_SPEED>();

      Vector<V, DESCRIPTOR::d> f1    = cellV.template getField<descriptors::FORCE>();
      Vector<V, DESCRIPTOR::n> force = {f1[0], f1[1], V(0), V(0), V(0)};

      Vector<V, DESCRIPTOR::n> u_num{};
      for (int iPopV = 0; iPopV < DESCRIPTOR::base_descriptor_t::q; iPopV++) {
        u_num += cellV[iPopV];
      }
      u_num += V(0.5) * force;

      Vector<V, DESCRIPTOR::n> phiX = {
        c_K * u_num[2] + c_mu * u_num[3],
        c_mu * u_num[4],
        c_K * u_num[0],
        c_mu * u_num[0],
        c_mu * u_num[1]
      };
      Vector<V, DESCRIPTOR::n> phiY = {
        c_mu * u_num[4],
        c_K * u_num[2] - c_mu * u_num[3],
        c_K * u_num[1],
        -c_mu * u_num[1],
        c_mu * u_num[0]
      };

      Vector<V, DESCRIPTOR::d> vel(u_num[0], u_num[1]);

      Vector<V, DESCRIPTOR::d> velOld = cellV.template getField<structuralmechanics::VELOCITY>();
      auto                     disp   = cellV.template getFieldPointer<structuralmechanics::DISPLACEMENT>();
      disp += V(0.5) * (velOld + vel);
      cellV.template setField<structuralmechanics::VELOCITY>(vel);

      Vector<V, DESCRIPTOR::d * DESCRIPTOR::d> stress = {
        c_K * u_num[2] + c_mu * u_num[3],
        c_mu * u_num[4],
        c_mu * u_num[4],
        c_K * u_num[2] - c_mu * u_num[3],
      };
      stress *= cLattice;
      cellV.template setField<structuralmechanics::CAUCHY_STRESS>(stress);

      for (int iPop = 0; iPop < DESCRIPTOR::base_descriptor_t::q; ++iPop) {
        Vector<V, DESCRIPTOR::n> f_eq =
            descriptors::t<V, DESCRIPTOR>(iPop) * (u_num + (descriptors::cVector<DESCRIPTOR>(iPop)[0] * phiX +
                                                            descriptors::cVector<DESCRIPTOR>(iPop)[1] * phiY) *
                                                               descriptors::invCs2<V, DESCRIPTOR>());

        cellV[iPop] =
            (omega * f_eq) + ((V(1.) - omega) * cellV[iPop]) +
            (V(2.) - omega) * descriptors::t<V, DESCRIPTOR>(iPop) * force / descriptors::invCs2<V, DESCRIPTOR>();
      }

      return {1, util::normSqr<V, DESCRIPTOR::d>(vel)};
    };

    void computeEquilibrium(ConstCell<T, DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d],
                            T fEq[DESCRIPTOR::q]) const override {
      EquilibriumF().compute(cell, rho, u, fEq);
    }

    std::string getName() const override { return "ElastodynamicsVectorLBM<" + MomentaF().getName() + ">"; };
  };

  template <typename V>
  struct ElastodynamicsFlux {

    /// Phi_x(U) = [ c_K j_s + c_mu j_d , c_mu j_xy , c_K v_x , c_mu v_x , c_mu v_y ]
    static void phiX(Vector<V, 5>& out, const Vector<V, 5>& U, V cK, V cMu) any_platform {
      out[0] =  cK  * U[2] + cMu * U[3];
      out[1] =  cMu * U[4];
      out[2] =  cK  * U[0];
      out[3] =  cMu * U[0];
      out[4] =  cMu * U[1];
    }

    /// Phi_y(U) = [ c_mu j_xy , c_K j_s - c_mu j_d , c_K v_y , -c_mu v_y , c_mu v_x ]
    static void phiY(Vector<V, 5>& out, const Vector<V, 5>& U, V cK, V cMu) any_platform {
      out[0] =  cMu * U[4];
      out[1] =  cK  * U[2] - cMu * U[3];
      out[2] =  cK  * U[1];
      out[3] = -cMu * U[1];
      out[4] =  cMu * U[0];
    }
  };


// ---------------------------------------------------------------------------
// Stage 1: macroscopic initial state U_0, Eq. (18)
// ---------------------------------------------------------------------------
template <typename V, typename DESCRIPTOR>
struct IniElastodynamicsUZero {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct C_K  : public descriptors::FIELD_BASE<1> { };
  struct C_MU : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<C_K, C_MU>;

  template <typename CELL>
  static bool isBulk(CELL& c) any_platform {
    return c.template getField<descriptors::POROSITY>() > V(5);
  }

  /// first-order derivative of the DISPLACEMENT field along e, second-order accurate
  template <typename CELL>
  static void diffDisp(Vector<V, DESCRIPTOR::d>&  deriv,
                       CELL&                      cell,
                       Vector<V, DESCRIPTOR::d>   e,
                       Vector<V, DESCRIPTOR::d>   eo) any_platform {
    auto p = cell.neighbor(e);
    auto m = cell.neighbor(eo);
    const Vector<V, DESCRIPTOR::d> s0 = cell.template getField<structuralmechanics::DISPLACEMENT>();

    if (!isBulk(p)) {
      auto m2 = m.neighbor(eo);
      deriv = (V(3) * s0
               - V(4) * m.template getField<structuralmechanics::DISPLACEMENT>()
               + m2.template getField<structuralmechanics::DISPLACEMENT>()) / V(2);
    }
    else if (!isBulk(m)) {
      auto p2 = p.neighbor(e);
      deriv = (-V(3) * s0
               + V(4) * p.template getField<structuralmechanics::DISPLACEMENT>()
               - p2.template getField<structuralmechanics::DISPLACEMENT>()) / V(2);
    }
    else {
      deriv = (p.template getField<structuralmechanics::DISPLACEMENT>()
               - m.template getField<structuralmechanics::DISPLACEMENT>()) / V(2);
    }
  }

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform {
    auto& cell = cells.template get<names::NavierStokes>();

    const V cK  = parameters.template get<C_K>();
    const V cMu = parameters.template get<C_MU>();

    Vector<V, DESCRIPTOR::d> x{1, 0}, xopp{-1, 0}, y{0, 1}, yopp{0, -1};

    Vector<V, DESCRIPTOR::d> dxu{}, dyu{};
    diffDisp(dxu, cell, x, xopp);           // dxu = [ d_x u_x , d_x u_y ]
    diffDisp(dyu, cell, y, yopp);           // dyu = [ d_y u_x , d_y u_y ]

    const Vector<V, DESCRIPTOR::d> v0 = cell.template getField<structuralmechanics::VELOCITY>();

    // Eqs. (14)-(16) evaluated on the initial displacement -> Eq. (18)
    Vector<V, 5> U0 = { v0[0],
                        v0[1],
                        -cK  * (dxu[0] + dyu[1]),
                        -cMu * (dxu[0] - dyu[1]),
                        -cMu * (dxu[1] + dyu[0]) };

    cell.template setField<structuralmechanics::U_ZERO<DESCRIPTOR>>(U0);
  }
};


// ---------------------------------------------------------------------------
// Stage 2: second-order consistent population initialization, Eq. (32)
// ---------------------------------------------------------------------------
template <typename V, typename DESCRIPTOR>
struct IniPopulationsElastodynamics {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct C_K           : public descriptors::FIELD_BASE<1> { };
  struct C_MU          : public descriptors::FIELD_BASE<1> { };
  struct LATTICE_SPEED : public descriptors::FIELD_BASE<1> { };
  struct SECOND_ORDER  : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<
    C_K,
    C_MU,
    LATTICE_SPEED,
    SECOND_ORDER
  >;

  using flux = ElastodynamicsFlux<V>;

  template <typename CELL>
  static bool isBulk(CELL& c) any_platform {
    return c.template getField<descriptors::POROSITY>() > V(5);
  }

  template <typename CELL>
  static Vector<V, 5> localState(CELL cell) any_platform {
    return cell.template getField<structuralmechanics::U_ZERO<DESCRIPTOR>>();
  }

  template <typename CELL>
  static void directionalDerivative(Vector<V, 5>&            deriv,
                                    CELL&                    cell,
                                    Vector<V, DESCRIPTOR::d> e,
                                    Vector<V, DESCRIPTOR::d> eo) any_platform {
    auto p = cell.neighbor(e);
    auto m = cell.neighbor(eo);
    const Vector<V, 5> s0 = localState(cell);

    if (!isBulk(p)) {
      auto m2 = m.neighbor(eo);
      deriv = (V(3) * s0 - V(4) * localState(m) + localState(m2)) / V(2);
    }
    else if (!isBulk(m)) {
      auto p2 = p.neighbor(e);
      deriv = (-V(3) * s0 + V(4) * localState(p) - localState(p2)) / V(2);
    }
    else {
      deriv = (localState(p) - localState(m)) / V(2);
    }
  }

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform {
    auto&      cell = cells.template get<names::NavierStokes>();
    VectorCell cellV(cell);

    const V cK     = parameters.template get<C_K          >();
    const V cMu    = parameters.template get<C_MU         >();
    const V c      = parameters.template get<LATTICE_SPEED>();
    const V order2 = parameters.template get<SECOND_ORDER >();

    const Vector<V, 5> U0 = cell.template getField<structuralmechanics::U_ZERO<DESCRIPTOR>>();

    Vector<V, 5> phiXU0{}, phiYU0{};
    flux::phiX(phiXU0, U0, cK, cMu);
    flux::phiY(phiYU0, U0, cK, cMu);

    const Vector<V, DESCRIPTOR::d> b = cell.template getField<descriptors::FORCE>();
    const Vector<V, 5> B = { b[0], b[1], V(0), V(0), V(0) };

    Vector<V, 5> phiXB{}, phiYB{};
    Vector<V, 5> dxU{},  dyU{};
    Vector<V, 5> AxdxU{}, AxAxdxU{}, AyAxdxU{};
    Vector<V, 5> AydyU{}, AxAydyU{}, AyAydyU{};

    if (order2 > V(0.5)) {
      flux::phiX(phiXB, B, cK, cMu);
      flux::phiY(phiYB, B, cK, cMu);

      Vector<V, DESCRIPTOR::d> x{1, 0}, xopp{-1, 0}, y{0, 1}, yopp{0, -1};
      directionalDerivative(dxU, cell, x, xopp);
      directionalDerivative(dyU, cell, y, yopp);

      flux::phiX(AxdxU,   dxU,   cK, cMu);
      flux::phiX(AxAxdxU, AxdxU, cK, cMu);
      flux::phiY(AyAxdxU, AxdxU, cK, cMu);

      flux::phiY(AydyU,   dyU,   cK, cMu);
      flux::phiX(AxAydyU, AydyU, cK, cMu);
      flux::phiY(AyAydyU, AydyU, cK, cMu);
    }

    for (int iPop = 0; iPop < DESCRIPTOR::base_descriptor_t::q; ++iPop) {
      const auto localC = descriptors::cVector<DESCRIPTOR>(iPop);
      const V    i = localC[0];
      const V    j = localC[1];

      const V w      = descriptors::t<V, DESCRIPTOR>(iPop);
      const V invCs2 = descriptors::invCs2<V, DESCRIPTOR>();

      Vector<V, 5> f = w * (U0 + invCs2 * (i * phiXU0 + j * phiYU0));

      if (order2 > V(0.5)) {
        const Vector<V, 5> T0 = B + invCs2 * (i * phiXB + j * phiYB);

        const Vector<V, 5> Tx = i * dxU
                              + (V(2) * i * i - V(1)) * AxdxU
                              - V(2) * (i * AxAxdxU + j * AyAxdxU);

        const Vector<V, 5> Ty = j * dyU
                              + (V(2) * j * j - V(1)) * AydyU
                              - V(2) * (i * AxAydyU + j * AyAydyU);

        f -= V(0.5) * w * (T0 + Tx + Ty);
      }

      cellV[iPop] = f;
    }

    const Vector<V, DESCRIPTOR::d> v0 = { U0[0], U0[1] };

    auto u = cellV.template getFieldPointer<structuralmechanics::DISPLACEMENT>();
    u -= v0;

    Vector<V, DESCRIPTOR::d * DESCRIPTOR::d> stress0 = {
      phiXU0[0],
      phiXU0[1],
      phiYU0[0],
      phiYU0[1]
    };
    stress0 *= c;
    cell.template setField<structuralmechanics::CAUCHY_STRESS>(stress0);
  }
};

  struct ElastodynamicsManufacturedVelocity {

    using parameters = meta::list<
      descriptors::MAGIC,
      structuralmechanics::CHAR_V,
      structuralmechanics::MS_AMPLITUDE,
      structuralmechanics::MS_WAVENUMBER
    >;

    template <typename DESCRIPTOR, typename V, typename PARAMETERS>
    static Vector<V, DESCRIPTOR::d> velocity(V x, V y, PARAMETERS& parameters) any_platform {
      const V convV = parameters.template get<structuralmechanics::CHAR_V       >();
      const V A     = parameters.template get<structuralmechanics::MS_AMPLITUDE >();
      const V k     = parameters.template get<structuralmechanics::MS_WAVENUMBER>();
      const V t     = parameters.template get<descriptors::MAGIC                >();

      Vector<V, DESCRIPTOR::d> vD = {V(0.80000)*A*k*util::sin(k*(-V(0.80000)*t + y))*util::sin(V(2)*k*(-V(0.30000)*t + x))*util::sin(V(2)*k*(t + V(-0.10000))) + V(2)*A*k*util::sin(V(2)*k*(-V(0.30000)*t + x))*util::cos(k*(-V(0.80000)*t + y))*util::cos(V(2)*k*(t + V(-0.10000))) - V(0.60000)*A*k*util::sin(V(2)*k*(t + V(-0.10000)))*util::cos(k*(-V(0.80000)*t + y))*util::cos(V(2)*k*(-V(0.30000)*t + x)),
        V(1.40000)*A*k*util::sin(V(2)*k*(-V(0.70000)*t + x))*util::sin(k*(-V(0.10000)*t + y))*util::sin(V(2)*k*(t + V(0.40000))) + V(2)*A*k*util::sin(k*(-V(0.10000)*t + y))*util::cos(V(2)*k*(-V(0.70000)*t + x))*util::cos(V(2)*k*(t + V(0.40000))) - V(0.10000)*A*k*util::sin(V(2)*k*(t + V(0.40000)))*util::cos(V(2)*k*(-V(0.70000)*t + x))*util::cos(k*(-V(0.10000)*t + y))
      };
      vD /= convV;

      return vD;
    }
  };

  // ==== Implementation of the Dirichlet Boundary condition for linear elastodynamics ====
  // following Boolakee, O. et al. (2025).
  // Lattice Boltzmann for linear elastodynamics: Periodic problems and Dirichlet boundary conditions
  // Computer Methods in Applied Mechanics and Engineering, 433.
  template <typename VELOCITY>
  struct ElastodynamicsVectorLBMDirichletFromVelocity {
    static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

    int getPriority() const { return -1; }

    using parameters = typename concat_parameter_lists<
      meta::list<
        structuralmechanics::C_K,
        structuralmechanics::C_MU
      >,
      typename VELOCITY::parameters
    >::type;

    template <typename CELL, typename PARAMETERS,
              typename DESCRIPTOR = typename CELL::descriptor_t,
              typename V          = typename CELL::value_t>
    void apply(CELL& cell, PARAMETERS& parameters) any_platform {
      constexpr int d = DESCRIPTOR::d;
      constexpr int n = DESCRIPTOR::n;

      VectorCell cellV(cell);

      const V cK  = parameters.template get<structuralmechanics::C_K >();
      const V cMu = parameters.template get<structuralmechanics::C_MU>();

      for (int iPopV = 0; iPopV < DESCRIPTOR::base_descriptor_t::q; ++iPopV) {
        const int          oppo = descriptors::oppositeVector<DESCRIPTOR>(iPopV);
        const Vector<V, d> c1   = descriptors::cVector<DESCRIPTOR>(iPopV);

        auto neighbor = cellV.neighbor(c1);
        if (!(neighbor.template getField<structuralmechanics::VLBM_SOLID_DISTANCE_FIELD>()[oppo] > V(0))) { continue; }

        const V x = neighbor.template getField<structuralmechanics::BOUNDARY_COORDS_X>()[oppo];
        const V y = neighbor.template getField<structuralmechanics::BOUNDARY_COORDS_Y>()[oppo];

        const Vector<V, d> vD = VELOCITY::template velocity<DESCRIPTOR>(x, y, parameters);

        const V i = c1[0];
        const V j = c1[1];

        const Vector<V, n> pop = neighbor[oppo];

        // D = diag(-1, -1, 1, 1, 1), Eq. (38)
        const Vector<V, n> Dpop = { -pop[0], -pop[1], pop[2], pop[3], pop[4] };

        // S_ij v_D, Eq. (39)
        const Vector<V, n> S = {
          V(0.5) * vD[0],
          V(0.5) * vD[1],
          i * cK  * vD[0] + j * cK  * vD[1],
          i * cMu * vD[0] - j * cMu * vD[1],
          j * cMu * vD[0] + i * cMu * vD[1]
        };

        cellV[iPopV] = Dpop + S;
      }
    }
  };

  using ElastodynamicsVectorLBMDirichletValidation = ElastodynamicsVectorLBMDirichletFromVelocity<ElastodynamicsManufacturedVelocity>;

  struct ElastodynamicsManufacturedTraction {

    using parameters = meta::list<
      descriptors::MAGIC,
      structuralmechanics::MS_AMPLITUDE,
      structuralmechanics::MS_WAVENUMBER,
      structuralmechanics::LAME_LAMBDA,
      structuralmechanics::LAME_MU
    >;

    template <typename DESCRIPTOR, typename V, typename PARAMETERS>
    static Vector<V, DESCRIPTOR::d * DESCRIPTOR::d> cauchyStress(V x, V y, PARAMETERS& parameters) any_platform {
      const V lambda = parameters.template get<structuralmechanics::LAME_LAMBDA  >();
      const V mu     = parameters.template get<structuralmechanics::LAME_MU      >();
      const V A      = parameters.template get<structuralmechanics::MS_AMPLITUDE >();
      const V k      = parameters.template get<structuralmechanics::MS_WAVENUMBER>();
      const V t      = parameters.template get<descriptors::MAGIC                >();

      const V sa1 = util::sin(V(2) * k * (x - V(0.3) * t)), ca1 = util::cos(V(2) * k * (x - V(0.3) * t));
      const V sb1 = util::sin(       k * (y - V(0.8) * t)), cb1 = util::cos(       k * (y - V(0.8) * t));
      const V sa2 = util::sin(V(2) * k * (x - V(0.7) * t)), ca2 = util::cos(V(2) * k * (x - V(0.7) * t));
      const V sb2 = util::sin(       k * (y - V(0.1) * t)), cb2 = util::cos(       k * (y - V(0.1) * t));
      const V T1  = util::sin(V(2) * k * (t - V(0.1)));
      const V T2  = util::sin(V(2) * k * (t + V(0.4)));

      const V exx  = V(2) * A * k * ca1 * cb1 * T1;                             // d_x u_x
      const V eyy  =        A * k * ca2 * cb2 * T2;                             // d_y u_y
      const V exy  = -V(0.5) * A * k * sa1 * sb1 * T1 - A * k * sa2 * sb2 * T2; // sym. shear strain
      const V divU = exx + eyy;

      return {lambda * divU + V(2) * mu * exx, V(2) * mu * exy,
              V(2) * mu * exy,                 lambda * divU + V(2) * mu * eyy};
    }

    template <typename DESCRIPTOR, typename V, typename PARAMETERS>
    static Vector<V, DESCRIPTOR::d> traction(V x, V y, const Vector<V, DESCRIPTOR::d>& n,
                                             PARAMETERS& parameters) any_platform {
      constexpr int d = DESCRIPTOR::d;
      const Vector<V, d * d> sigma = cauchyStress<DESCRIPTOR>(x, y, parameters);
      Vector<V, d> T{};
      for (int a = 0; a < d; ++a) {
        for (int b = 0; b < d; ++b) {
          T[a] += sigma[a * d + b] * n[b];
        }
      }
      return T;
    }
  };

  // ==== Implementation of the Neumann boundary condition for linear elastodynamics ====
  // following Boolakee, O. et al. (2025).
  // Lattice Boltzmann for Solids
  // Doctoral Thesis, doi.org/10.3929/ethz-b-000742210
  template<typename TRACTION>
  struct ElastodynamicsVectorLBMNeumannFromTraction {
    static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

    int getPriority() const {
      return -1;
    }

    using parameters = typename concat_parameter_lists<
      meta::list<
        structuralmechanics::C_K,
        structuralmechanics::C_MU
      >,
      typename TRACTION::parameters
    >::type;

    template <typename CELL, typename PARAMETERS,
              typename DESCRIPTOR = CELL::descriptor_t,
              typename V = typename CELL::value_t>
    void apply(CELL& cell, PARAMETERS& parameters) any_platform {
      VectorCell cellV(cell);

      V c_K       = parameters.template get<structuralmechanics::C_K>();
      V c_mu      = parameters.template get<structuralmechanics::C_MU>();

      V a1    = (c_K * c_K - c_mu * c_mu) / (c_K * c_K + c_mu * c_mu);
      V a2    = V(2.) * c_K * c_mu / (c_K * c_K + c_mu * c_mu);

      V b1    = V(0.5) * c_K  / (c_K * c_K + c_mu * c_mu);
      V b2    = V(0.5) * c_mu / (c_K * c_K + c_mu * c_mu);
      V b3    = V(0.5) / c_mu;

      for (int iPopV = 0; iPopV < DESCRIPTOR::base_descriptor_t::q; iPopV++) {
        Vector<V, DESCRIPTOR::d> localC   = descriptors::cVector<DESCRIPTOR>(iPopV);
        auto neighbor                     = cellV.neighbor(localC);

        int oppo                          = descriptors::oppositeVector<DESCRIPTOR>(iPopV);
        if (!(neighbor.template getField<structuralmechanics::VLBM_SOLID_DISTANCE_FIELD>()[oppo] > V(0))) { continue; }

        int i = localC[0];
        int j = localC[1];

        Vector<V, DESCRIPTOR::n> pop = neighbor[oppo];
        Vector<V, DESCRIPTOR::n> N_ij_pop = {
          pop[0],
          pop[1],
          -a1 * pop[2] + (j * j - i * i) * a2 * pop[3],
          (j * j - i * i) * a2 * pop[2] + a1 * pop[3],
          -pop[4]
        };

        const V x = neighbor.template getField<structuralmechanics::BOUNDARY_COORDS_X>()[oppo];
        const V y = neighbor.template getField<structuralmechanics::BOUNDARY_COORDS_Y>()[oppo];

        const Vector<V, DESCRIPTOR::d> n = {
          neighbor.template getField<structuralmechanics::NEUMANN_NORMAL_X>()[oppo],
          neighbor.template getField<structuralmechanics::NEUMANN_NORMAL_Y>()[oppo]
        };

        Vector<V, DESCRIPTOR::d> Traction = TRACTION::template traction<DESCRIPTOR>(x, y, n, parameters);

        Vector<V, DESCRIPTOR::n> S_ij_trac = {
          (i * i + j * j) * Traction[0],
          (i * i + j * j) * Traction[1],
          (i * b1 * Traction[0] + j * b1 * Traction[1]),
          (i * b2 * Traction[0] - j * b2 * Traction[1]),
          (j * b3 * Traction[0] + i * b3 * Traction[1])
        };

        cellV[iPopV] = N_ij_pop + S_ij_trac;
      }
    }
  };

  using ElastodynamicsVectorLBMNeumannValidation = ElastodynamicsVectorLBMNeumannFromTraction<ElastodynamicsManufacturedTraction>;

} // namespace structuralmechanics
} // namespace olb
#endif
