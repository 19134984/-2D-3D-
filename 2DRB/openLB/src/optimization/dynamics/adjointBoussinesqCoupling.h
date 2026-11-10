/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Ecem Öz, Shota Ito
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

#ifndef ADJOINT_BOUSSINESQ_COUPLING_H
#define ADJOINT_BOUSSINESQ_COUPLING_H

#include "core/operator.h"

namespace olb {

/// Adjoint coupling operator of primal NavierStokesAdvectionDiffusionCoupling implementation
// Computes the four Jacobian matrices: collision operators for NSE and ADE lattice differentiated
// regarding the populations of the NSE and ADE, respectively. To prevent storing matrices in fields
// this coupling operator stores the product of the Jacobians with the adjoint populations.
// TODO: Once code generation is extended for coupling operators, compute the four Jacobians using fAD.
struct AdjointBoussinesqCoupling {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct FORCE_PREFACTOR : public descriptors::FIELD_BASE<0,1> { };
  struct T0 : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_NSE : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_ADE : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<FORCE_PREFACTOR,T0,OMEGA_NSE,OMEGA_ADE>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR_NSE = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;
    using DESCRIPTOR_ADE = typename CELLS::template value_t<names::Temperature>::descriptor_t;

    auto& cellNSE = cells.template get<names::NavierStokes>();
    auto& cellADE = cells.template get<names::Temperature>();
    const V omega_NSE = parameters.template get<OMEGA_NSE>();
    const V omega_ADE = parameters.template get<OMEGA_ADE>();
    const V t0 = parameters.template get<T0>();
    auto f  = cellNSE.template getField<opti::F>(); //primalNSE population
    auto g  = cellADE.template getField<opti::G>(); //primalADE population

    // Preparation
    V rho = V();
    V u[DESCRIPTOR_NSE::d] { };
    V u_kin[DESCRIPTOR_NSE::d]{ };
    V T = V();
    auto force = cellNSE.template getFieldPointer<descriptors::FORCE>();
    auto forcePrefactor = parameters.template get<FORCE_PREFACTOR>();

    for(int iPop = 0 ; iPop < DESCRIPTOR_ADE::q; ++iPop){
      T += g[iPop];
    }
    T += V(1);
    V temperatureDifference = T - t0;
    for (unsigned iD = 0; iD < DESCRIPTOR_NSE::d; ++iD) {
      force[iD] = forcePrefactor[iD] * temperatureDifference;
    }
    for(int iPop = 0 ; iPop < DESCRIPTOR_NSE::q; ++iPop){
      rho += f[iPop];
      for(int iD = 0; iD < DESCRIPTOR_NSE::d; ++iD){
        u[iD] += descriptors::c<DESCRIPTOR_NSE>(iPop,iD) * f[iPop];
      }
    }
    rho += V(1);
    for(int iD = 0; iD < DESCRIPTOR_NSE::d; iD++){
      u_kin[iD] = u[iD] / rho;
    }
    for(int iD = 0; iD < DESCRIPTOR_NSE::d; iD++){
      u[iD] = (u[iD] / rho + V(0.5) * force[iD]);
    }

    // Dervative of NSE-lattice collision regarding NSE-lattice population multiplied by
    // adjoint NSE-lattice populations ((\partial C_{nse}/ \partial f)^T \phi)
    auto dCnsedFxphi = cellNSE.template getFieldPointer<opti::DCnseDFxphi>();
    for (int jPop = 0; jPop < DESCRIPTOR_NSE::q; ++jPop){
      V duDf[DESCRIPTOR_NSE::d] { };
      for(int iD = 0; iD < DESCRIPTOR_NSE::d; ++iD){
        duDf[iD] = (descriptors::c<DESCRIPTOR_NSE>(jPop,iD) - u_kin[iD]) / rho ;
      }
      dCnsedFxphi[jPop] = V();
      for (int iPop = 0; iPop < DESCRIPTOR_NSE::q; ++iPop){
        V delta = (iPop == jPop) ? V(1) : V(0);

        V wi = descriptors::t<V,DESCRIPTOR_NSE>(iPop);

        V xi_i_dot_u = V();
        V u2 = V();
        V xi_i_dot_xi_j_minus_u = V();
        V xi_i_dot_duDf = V();
        V u_dot_duDf = V();
        V u_dot_xi_j_minus_u = V();

        for(int iD = 0; iD < DESCRIPTOR_NSE::d; ++iD){
          V xi_i = descriptors::c<DESCRIPTOR_NSE>(iPop,iD);
          V xi_j_minus_u = descriptors::c<DESCRIPTOR_NSE>(jPop,iD) - u_kin[iD];

          xi_i_dot_u += xi_i * u[iD];
          u2 += u[iD] * u[iD];
          xi_i_dot_duDf +=  xi_i * duDf[iD];
          u_dot_duDf += u[iD] * duDf[iD];
          xi_i_dot_xi_j_minus_u += xi_i * xi_j_minus_u;
          u_dot_xi_j_minus_u += u[iD] * xi_j_minus_u;
        }
        V phi_i_byhand = V(1)
                        + xi_i_dot_u * descriptors::invCs2<V,DESCRIPTOR_NSE>()
                        + V(0.5) * xi_i_dot_u * xi_i_dot_u
                        * descriptors::invCs2<V,DESCRIPTOR_NSE>()
                        * descriptors::invCs2<V,DESCRIPTOR_NSE>()
                        - V(0.5) * u2 * descriptors::invCs2<V,DESCRIPTOR_NSE>();

        V eqDeriv_byhand = wi * phi_i_byhand + wi * rho * (xi_i_dot_duDf * descriptors::invCs2<V,DESCRIPTOR_NSE>()
                                                          + xi_i_dot_u * xi_i_dot_duDf * descriptors::invCs2<V,DESCRIPTOR_NSE>()
                                                          * descriptors::invCs2<V,DESCRIPTOR_NSE>()
                                                          - u_dot_duDf * descriptors::invCs2<V,DESCRIPTOR_NSE>() );

        V Aji_byhand = delta - omega_NSE * (delta - eqDeriv_byhand);

        V dotWithForce = V();
        for( int iD = 0; iD <DESCRIPTOR_NSE::d; ++iD){
          V xi_i = descriptors::c<DESCRIPTOR_NSE>(iPop,iD);
          V component = - duDf[iD] * descriptors::invCs2<V,DESCRIPTOR_NSE>()
                        + xi_i_dot_duDf * xi_i * descriptors::invCs2<V,DESCRIPTOR_NSE>() * descriptors::invCs2<V,DESCRIPTOR_NSE>();
          dotWithForce += component * force[iD];
        }
        V Fb_deriv_f = (V(1) - V(0.5) * omega_NSE) * wi * dotWithForce;
        dCnsedFxphi[jPop] += ( Aji_byhand + Fb_deriv_f) * cellNSE[iPop];
      }
    }

    // Dervative of ADE-lattice collision regarding NSE-lattice population multiplied by
    // adjoint ADE-lattice populations ((\partial C_{ade}/ \partial f)^T \lambda)
    auto dCadedFxlambda = cellNSE.template getFieldPointer<opti::DCadeDFxlambda>();
    for(int jPop = 0; jPop < DESCRIPTOR_NSE::q; ++jPop){
      V duDf[DESCRIPTOR_NSE::d] { };
      for(int iD = 0; iD < DESCRIPTOR_NSE::d; ++iD){
        duDf[iD] = (descriptors::c<DESCRIPTOR_NSE>(jPop,iD) - u_kin[iD]) / rho ;
      }
      dCadedFxlambda[jPop] = V();
      for ( int iPop = 0; iPop < DESCRIPTOR_ADE::q; ++iPop){
        V wi = descriptors::t<V,DESCRIPTOR_ADE>(iPop);
        V xi_dot_duDf = V();

        for( int iD = 0; iD < DESCRIPTOR_ADE::d; ++iD){
          V xi_i = descriptors::c<DESCRIPTOR_ADE>(iPop,iD);
          xi_dot_duDf += xi_i * duDf[iD];
        }

        V geq_deriv = wi * T *  xi_dot_duDf * descriptors::invCs2<V,DESCRIPTOR_ADE>();
        dCadedFxlambda[jPop] += omega_ADE * geq_deriv * cellADE[iPop];
      }
    }

    // Dervative of NSE-lattice collision regarding ADE-lattice population multiplied by
    // adjoint NSE-lattice populations ((\partial C_{nse}/ \partial g)^T \phi)
    auto dCnsedGxphi = cellADE.template getFieldPointer<opti::DCnseDGxphi>();
    V duDg[DESCRIPTOR_NSE::d] { };
    for(int iD = 0; iD < DESCRIPTOR_NSE::d; ++iD){
      duDg[iD] = V(0.5) * forcePrefactor[iD] ;
    }
    for (int jPop = 0 ; jPop < DESCRIPTOR_ADE::q ; ++jPop){
      dCnsedGxphi[jPop] = V();
      for( int iPop = 0; iPop < DESCRIPTOR_NSE::q ; ++iPop){
       V wi = descriptors::t<V,DESCRIPTOR_NSE>(iPop);
       V xi_i_dot_u = V();
       V xi_i_dot_duDg = V();
       V u_dot_duDg = V();

       for (int iD = 0; iD < DESCRIPTOR_NSE::d; ++iD){
         V xi_i = descriptors::c<DESCRIPTOR_NSE>(iPop,iD);
         xi_i_dot_u += xi_i * u[iD];
         xi_i_dot_duDg += xi_i * duDg[iD];
         u_dot_duDg += u[iD] * duDg[iD];
       }

       V dfeqDg = wi * rho * (xi_i_dot_duDg * descriptors::invCs2<V,DESCRIPTOR_NSE>()
                               + xi_i_dot_u *  xi_i_dot_duDg * descriptors::invCs2<V,DESCRIPTOR_NSE>()
                               * descriptors::invCs2<V,DESCRIPTOR_NSE>()
                               - u_dot_duDg * descriptors::invCs2<V,DESCRIPTOR_NSE>());

       V dotWithForce = V();
       for (int iD = 0; iD < DESCRIPTOR_NSE::d; iD++){
         V xi_i = descriptors::c<DESCRIPTOR_NSE>(iPop,iD);
         V xi_i_minus_u = xi_i - u[iD];

         V component_1 =  xi_i_minus_u
                       * descriptors::invCs2<V,DESCRIPTOR_NSE>()
                       + xi_i_dot_u
                       * xi_i
                       * descriptors::invCs2<V,DESCRIPTOR_NSE>()
                       * descriptors::invCs2<V,DESCRIPTOR_NSE>() ;
         V component_2 =  -  duDg[iD]
                         * descriptors::invCs2<V,DESCRIPTOR_NSE>()
                         + xi_i_dot_duDg
                         * xi_i
                         * descriptors::invCs2<V,DESCRIPTOR_NSE>()
                         * descriptors::invCs2<V,DESCRIPTOR_NSE>() ;

         dotWithForce += component_1 * forcePrefactor[iD] + component_2 * force[iD];
       }
       V Fb_deriv_g =  (V(1) - V(0.5) * omega_NSE) * wi * dotWithForce ;
       dCnsedGxphi[jPop] += (  omega_NSE * dfeqDg + Fb_deriv_g) * cellNSE[iPop];
     }
   }

    // Dervative of ADE-lattice collision regarding ADE-lattice population multiplied by
    // adjoint ADE-lattice populations ((\partial C_{nse}/ \partial g)^T \phi)
    auto dCadedGxlambda= cellADE.template getFieldPointer<opti::DCadeDGxlambda>();
    for(int jPop = 0; jPop < DESCRIPTOR_ADE::q; ++jPop){
      dCadedGxlambda[jPop] = V();
      for(int iPop = 0; iPop < DESCRIPTOR_ADE::q; ++iPop){
        V delta = (iPop == jPop) ? V(1) : V(0);
        V wi = descriptors::t<V,DESCRIPTOR_ADE>(iPop);
        V xi_i_dot_u = V();
        V xi_i_dot_duDg = V();

        for (int iD = 0; iD < DESCRIPTOR_ADE::d; ++iD){
          V xi_i = descriptors::c<DESCRIPTOR_ADE>(iPop,iD);
          xi_i_dot_u += xi_i * u[iD];
          xi_i_dot_duDg += xi_i * duDg[iD];
        }

        V g_eqDeriv = wi * (V(1)
                          + xi_i_dot_u * descriptors::invCs2<V,DESCRIPTOR_ADE>()
                          + T * xi_i_dot_duDg * descriptors::invCs2<V,DESCRIPTOR_ADE>() );

        V Aij = delta - omega_ADE * (delta - g_eqDeriv);
        dCadedGxlambda[jPop] += Aij * cellADE[iPop];
      }
    }
  }
};

/// Evaluate the optimality condition using primal and adjoint results for the
/// primal NavierStokesAdvectionDiffusionCoupling.
// TODO: Once functors::DerivativeF is extended for coupling operators, compute the Jacobians using fAD.
struct AdjointBoussinesqOptimalityO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct T0 : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_NSE : public descriptors::FIELD_BASE<1> { };
  struct OMEGA_ADE : public descriptors::FIELD_BASE<1> { };
  struct D_FORCE_PREFACTOR_D_ALPHA : public descriptors::FIELD_BASE<0,1> { };
  struct FORCE_PREFACTOR : public descriptors::FIELD_BASE<0,1> { };
  struct CONVERSION_VELOCITY : public descriptors::FIELD_BASE<1> { };
  struct NORMALIZE : public descriptors::FIELD_BASE<1> { };
  struct FACTOR : public descriptors::FIELD_BASE<1> { };

  struct REFERENCE_VELOCITY : public descriptors::FIELD_BASE<0,1> { };
  struct SENSITIVITY : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<T0, OMEGA_NSE, OMEGA_ADE, D_FORCE_PREFACTOR_D_ALPHA, FORCE_PREFACTOR, CONVERSION_VELOCITY, NORMALIZE, FACTOR>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR_NSE = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;
    using DESCRIPTOR_ADE = typename CELLS::template value_t<names::Temperature>::descriptor_t;

    auto& cellNSE = cells.template get<names::NavierStokes>();
    auto& cellADE = cells.template get<names::Temperature>();

    const V omega_NSE = parameters.template get<OMEGA_NSE>();
    const V omega_ADE = parameters.template get<OMEGA_ADE>();
    const V t0 = parameters.template get<T0>();
    const V conversion = parameters.template get<CONVERSION_VELOCITY>();
    const V normalize = parameters.template get<NORMALIZE>();
    const V factor = parameters.template get<FACTOR>();

    auto d_force_prefactor_d_alpha = parameters.template get<D_FORCE_PREFACTOR_D_ALPHA>();
    auto forcePrefactor = parameters.template get<FORCE_PREFACTOR>();

    auto f  = cellNSE.template getField<opti::F>(); //primalNSE population
    auto g  = cellADE.template getField<opti::G>(); //primalADE population
    auto u_ref = cellNSE.template getField<REFERENCE_VELOCITY>();

    // Preparation
    V T = V();
    for(int iPop = 0 ; iPop < DESCRIPTOR_ADE::q; ++iPop){
      T += g[iPop];
    }
    T += V(1);
    V temperatureDifference = T - t0;
    V force[DESCRIPTOR_NSE::d]{ };
    for(int iD = 0; iD < DESCRIPTOR_NSE::d; ++iD){
      force[iD] = forcePrefactor[iD] * temperatureDifference;
    }
    V rho = V();
    V u[DESCRIPTOR_NSE::d] { };
    for(int iPop = 0 ; iPop < DESCRIPTOR_NSE::q; ++iPop){
      rho += f[iPop];
      for(int iD = 0; iD < DESCRIPTOR_NSE::d; ++iD){
        u[iD] += descriptors::c<DESCRIPTOR_NSE>(iPop,iD) * f[iPop];
      }
    }
    rho += V(1);
    for(int iD = 0; iD < DESCRIPTOR_NSE::d; iD++){
      u[iD] = u[iD] / rho + V(0.5) * force[iD];
    }
    V physU[DESCRIPTOR_NSE::d] {};
    for(int iD = 0 ; iD <DESCRIPTOR_NSE::d; iD++){
        physU[iD] = u[iD] * conversion;
    }
    V duDalpha[DESCRIPTOR_NSE::d] { };
    for(int iD = 0 ; iD < DESCRIPTOR_NSE::d; ++iD){
      duDalpha[iD] =  V(0.5) * d_force_prefactor_d_alpha[iD] * temperatureDifference;
    }

    // Total derivative of objective regarding controls
    auto sensitivity = cellNSE.template getField<SENSITIVITY>();
    // Derivative contains those three contributions
    V dC_nse_dalpha = V();
    V dC_ade_dalpha = V();
    V dJ_dalpha = V();

    // Partial derivative of objective regarding the controls
    for(int iD = 0; iD < DESCRIPTOR_NSE::d; ++iD){
      dJ_dalpha += (physU[iD] - u_ref[iD]) *  conversion *  V(0.5) * d_force_prefactor_d_alpha[iD] * temperatureDifference;
    }
    dJ_dalpha *= factor  / normalize;

    // Partial derivative of NSE-lattice collision operator regarding controls
    for(int iPop = 0; iPop < DESCRIPTOR_NSE::q ; ++iPop){
      V wi = descriptors::t<V,DESCRIPTOR_NSE>(iPop);
      V xi_i_dot_u = V();
      V xi_i_dot_dudalpha = V();
      V u_dot_dudalpha = V();

      for (int iD = 0; iD < DESCRIPTOR_NSE::d; iD++){
        V xi_i = descriptors::c<DESCRIPTOR_NSE>(iPop,iD);
        xi_i_dot_u += xi_i * u[iD];
        xi_i_dot_dudalpha += xi_i * duDalpha[iD];
        u_dot_dudalpha +=  u[iD] * duDalpha[iD];
      }

      V df_eq_dalpha = wi * rho
                       * ( xi_i_dot_dudalpha * descriptors::invCs2<V,DESCRIPTOR_NSE>()
                       + xi_i_dot_u * xi_i_dot_dudalpha * descriptors::invCs2<V,DESCRIPTOR_NSE>() * descriptors::invCs2<V,DESCRIPTOR_NSE>()
                       - u_dot_dudalpha * descriptors::invCs2<V,DESCRIPTOR_NSE>());

      V dotWithForce = V();
      for(int iD = 0; iD < DESCRIPTOR_NSE::d; iD++){
        V xi_i = descriptors::c<DESCRIPTOR_NSE>(iPop,iD);
        V component_1 = - duDalpha[iD] *  descriptors::invCs2<V,DESCRIPTOR_NSE>()
                        + xi_i_dot_dudalpha * xi_i * descriptors::invCs2<V,DESCRIPTOR_NSE>() * descriptors::invCs2<V,DESCRIPTOR_NSE>();
        V component_2 =  (xi_i - u[iD]) * descriptors::invCs2<V,DESCRIPTOR_NSE>()
                        + xi_i_dot_u * xi_i * descriptors::invCs2<V,DESCRIPTOR_NSE>() * descriptors::invCs2<V,DESCRIPTOR_NSE>() ;
        dotWithForce += component_1 * force[iD] + component_2 * d_force_prefactor_d_alpha[iD] * temperatureDifference;
      }

      V dFbdalpha = (V(1) - V(0.5) * omega_NSE) * wi * dotWithForce;
      dC_nse_dalpha -= (omega_NSE * df_eq_dalpha + dFbdalpha) * cellNSE[iPop];
    }

    // Partial derivative of ADE-lattice collision operator regarding controls
    for( int iPop = 0; iPop < DESCRIPTOR_ADE::q ; ++iPop){
      V wi = descriptors::t<V,DESCRIPTOR_ADE>(iPop);
      V xi_i_dot_dudalpha = V();

      for (int iD = 0; iD < DESCRIPTOR_ADE::d; iD++){
        V xi_i = descriptors::c<DESCRIPTOR_ADE>(iPop,iD);
        xi_i_dot_dudalpha += xi_i * duDalpha[iD];
      }

      V dg_eq_dalpha = wi * T * xi_i_dot_dudalpha * descriptors::invCs2<V,DESCRIPTOR_ADE>();
      dC_ade_dalpha -= omega_ADE * dg_eq_dalpha * cellADE[iPop];
    }
    // Compute final total derivative
    sensitivity = dJ_dalpha + dC_nse_dalpha + dC_ade_dalpha;
    cellNSE.template setField<SENSITIVITY>(sensitivity);
  }
};

/// Computes the partial derivative of the objective regarding the ADE-lattice populations
/// for the functors::L2Distance<functors::VelocityF> objective with the
/// AdjointBoussinesqCoupling.
// TODO: Once functors::DerivativeF is extended for coupling operators, compute the Jacobians using fAD.
struct AdjointBoussinesqObjectiveDerivativeO{
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct T0 : public descriptors::FIELD_BASE<1> { };
  struct CONVERSION_VELOCITY : public descriptors::FIELD_BASE<1> { };
  struct NORMALIZE : public descriptors::FIELD_BASE<1> { };
  struct FACTOR : public descriptors::FIELD_BASE<1> { };
  struct FORCE_PREFACTOR : public descriptors::FIELD_BASE<0,1> { };

  struct REFERENCE_VELOCITY : public descriptors::FIELD_BASE<0,1> { };

  using parameters = meta::list<T0, CONVERSION_VELOCITY, NORMALIZE, FACTOR, FORCE_PREFACTOR>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::NavierStokes>::value_t;
    using DESCRIPTOR_NSE = typename CELLS::template value_t<names::NavierStokes>::descriptor_t;
    using DESCRIPTOR_ADE = typename CELLS::template value_t<names::Temperature>::descriptor_t;

    const V t0 = parameters.template get<T0>();
    const V conversion = parameters.template get<CONVERSION_VELOCITY>();
    const V normalize = parameters.template get<NORMALIZE>();
    const V factor = parameters.template get<FACTOR>();
    auto forcePrefactor = parameters.template get<FORCE_PREFACTOR>();
    auto& cellNSE = cells.template get<names::NavierStokes>();
    auto& cellADE = cells.template get<names::Temperature>();
    auto f  = cellNSE.template getField<opti::F>();
    auto g  = cellADE.template getField<opti::G>();
    auto u_ref = cellNSE.template getField<REFERENCE_VELOCITY>();

    // Preparation
    V T = V();
    for(int iPop = 0 ; iPop < DESCRIPTOR_ADE::q; ++iPop){
      T += g[iPop];
    }
    T += V(1);
    V temperatureDifference = T - t0;
    V force[DESCRIPTOR_NSE::d]{ };
    for(int iD = 0; iD < DESCRIPTOR_NSE::d; ++iD){
      force[iD] = forcePrefactor[iD] * temperatureDifference;
    }
    V rho = V();
    V u[DESCRIPTOR_NSE::d] { };
    for(int iPop = 0 ; iPop < DESCRIPTOR_NSE::q; ++iPop){
      rho += f[iPop];
      for(int iD = 0; iD < DESCRIPTOR_NSE::d; ++iD){
        u[iD] += descriptors::c<DESCRIPTOR_NSE>(iPop,iD) * f[iPop];
      }
    }
    rho += V(1);
    for(int iD = 0; iD < DESCRIPTOR_NSE::d; iD++){
      u[iD] = u[iD] / rho + V(0.5) * force[iD];
    }
    V physU[DESCRIPTOR_NSE::d] {};
    for(int iD = 0 ; iD <DESCRIPTOR_NSE::d; iD++){
        physU[iD] = u[iD] * conversion;
    }

    // Compute the partial derivative of objective reg. ADE-lattice populations
    auto dJdg = cellADE.template getFieldPointer<opti::DJDG>();
    V result_dot = V();
    for(int iD = 0; iD < DESCRIPTOR_NSE::d; ++iD){
      result_dot += (physU[iD] - u_ref[iD]) * forcePrefactor[iD];
    }
    const V derivative = factor * V(0.5) / normalize * result_dot * conversion;
    for(int jPop = 0; jPop < DESCRIPTOR_ADE::q ; ++jPop){
      dJdg[jPop] = derivative; // for each jPop it is constant as objective is linear in g
    }
  }
};

}

#endif
