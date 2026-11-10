/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2024 Fedor Bukreev, Adrian Kummerländer
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

#ifndef TURBULENT_WALL_MODEL_POST_PROCESSOR_OPTIMIZED_H
#define TURBULENT_WALL_MODEL_POST_PROCESSOR_OPTIMIZED_H

namespace olb {
namespace optimized {

  // interpolation of the velocity along the surface normal on the chosen distance from the boundary cell
  template <typename CELL, typename V = typename CELL::value_t, typename DESCRIPTOR = typename CELL::descriptor_t>
  void interpolateVelocity(CELL& cell, Vector<V,DESCRIPTOR::d> distance, V u_y2[DESCRIPTOR::d]) any_platform{
    Vector<V,DESCRIPTOR::d> floorV;
    for( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
      floorV[iD] = util::floor(distance[iD]);
    }
    Vector<Vector<V,DESCRIPTOR::d>,(DESCRIPTOR::d==2)*4+(DESCRIPTOR::d==3)*8> surroundingPoints(floorV);
    surroundingPoints[1][0] += 1.;
    surroundingPoints[2][1] += 1.;
    surroundingPoints[3][0] += 1.; surroundingPoints[3][1] += 1.;
    if( DESCRIPTOR::d == 3) {
      surroundingPoints[4][2] += 1.;
      surroundingPoints[5][0] += 1.; surroundingPoints[5][2] += 1.;
      surroundingPoints[6][1] += 1.; surroundingPoints[6][2] += 1.;
      surroundingPoints[7][0] += 1.; surroundingPoints[7][1] += 1.; surroundingPoints[7][2] += 1.;
    }

    for (auto point : surroundingPoints) {
      const Vector<V,DESCRIPTOR::d> dist = distance - point;
      V weight = V(1);
      for( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
        weight *=(V(1.) - util::abs(dist[iD]));
      }
      V uNP[DESCRIPTOR::d] = {0.};
      cell.neighbor(point).computeU(uNP);
      for( int iD = 0; iD < DESCRIPTOR::d; iD++ ) {
        u_y2[iD] += uNP[iD]*weight;
      }
    }
  };

  template <int MODE, typename V>
  V spaldingWallFunction(V x, V u, V y, V nu, V y_plus_fixed = V(0)) any_platform {
    const V k = V(0.4);
    const V b = V(5.5);

    const V u_plus = (MODE == 2) ? (u / x) : x;
    const V y_plus = (MODE == 2) ? (y / nu * x) : y_plus_fixed;

    const V kp4 = k * u_plus;

    return u_plus + util::exp(-k * b) * (
             util::exp(kp4) - V(1) - kp4 - V(0.5) * util::pow(kp4, V(2)) - V(1./6.) * util::pow(kp4, V(3))
           ) - y_plus;
  }

  template <int MODE, typename V>
  V spaldingDerivative(V x, V u, V y, V nu) any_platform {
    const V k = V(0.4);
    const V b = V(5.5);
    const V exp_kb = util::exp(-k * b);

    if constexpr (MODE == 2) {
      const V x_inv = V(1) / x;
      const V x_inv2 = x_inv * x_inv;
      const V kp1_x = k * u * x_inv;

      return (-u * x_inv2) + exp_kb * (
               util::exp(kp1_x) * (-k * u * x_inv2) +
               k * u * x_inv2 +
               util::pow(kp1_x, V(2)) * x_inv +
               V(0.5) * util::pow(kp1_x, V(3)) * x_inv
             ) - (y / nu);
    } else {
      return V(1) + exp_kb * (
               util::exp(k * x) * k - k - k * k * x - V(0.5) * k * util::pow(k * x, V(2))
             );
    }
  }

  template <unsigned ITERS, typename V, typename Funcf, typename Funcd>
  V newtonSpalding(Funcf f, Funcd d, V x0, V convergeCrit) any_platform {
    const V eps = std::numeric_limits<V>::epsilon();
    V x = x0;
    #pragma unroll
    for (unsigned i = 0; i < ITERS; ++i) {
      V f_val = f(x);
      V d_val = d(x);

      // ensure minimum step, prevent NaN
      V filtered_d = std::copysign(util::max(util::abs(d_val), eps), d_val);
      V x_next = x - f_val / filtered_d;

      bool converged = util::abs(f_val) < convergeCrit;
      bool invalid = util::isnan(x_next);

      x = !(converged || invalid) ? x_next : x;
    }
    return x;
  }

  //Power Law profile
  //output = {u_tau, u_plus(x1)}
  template <typename V>
  Vector<V,2> powerLawWallFunction(V nu, V u2, V y2, V y1) any_platform{
    Vector<V,2> output(V(0), V(0));
    V u_tau = util::sqrt(V(0.0246384) * util::pow(nu, V(0.25)) * util::pow(u2, V(1.75)) / util::pow(y2, V(0.25)));
    V y_plus = y1 * u_tau / nu;

    if (y_plus > V(28.178)) {
      output[0] = u_tau;
      output[1] = V(5.5) + V(1./0.399) * util::log(y_plus);
    }
    else if (y_plus < V(28.178)  && y_plus > V(5.)) {
      output[0] = u_tau;
      output[1] = V(-3.05) + V(5.) * util::log(y_plus);
    }
    else {
      output[0] = util::sqrt(V(2.) * u2 * nu / y2);
      y_plus = y1 * output[0] / nu;
      output[1] = y_plus;
    }
    return output;
  };



//======================================================================
// ======== Wall function velocity and vsiscosity PostProcessor ======//
//======================================================================

template <bool bodyForce, bool interpolateSampleVelocity, bool useVanDriest, int wallFunctionProfile, bool movingWall>
struct TurbulentWallModelPostProcessor {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::list<descriptors::OMEGA, descriptors::SAMPLING_DISTANCE>;

  int getPriority() const {
    return 1;
  }

  template <typename CELL, typename PARAMETERS, typename V = typename CELL::value_t>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform {
    using DESCRIPTOR = typename CELL::descriptor_t;
    const V eps = std::numeric_limits<V>::epsilon();

    auto y1 = cell.template getField<descriptors::Y1>();
    V y1Norm = util::norm<DESCRIPTOR::d>(y1);

    if (y1Norm > V(1e-3)) {
      Vector<V, DESCRIPTOR::d> normal;
      for (int iD = 0; iD < DESCRIPTOR::d; iD++) {
        normal[iD] = y1[iD] / y1Norm;
      }

      V samplingCellDistance = parameters.template get<descriptors::SAMPLING_DISTANCE>();
      Vector<V, DESCRIPTOR::d> y2 = samplingCellDistance * normal;
      V u_y2[DESCRIPTOR::d] = {0.};
      auto y1y2 = y2 - y1;

      if constexpr (interpolateSampleVelocity) {
        optimized::interpolateVelocity(cell, y1y2, u_y2);
      } else {
        Vector<V, DESCRIPTOR::d> roundedNormal;
        for (int iD = 0; iD < DESCRIPTOR::d; iD++) {
          roundedNormal[iD] = int(util::round(y1y2[iD]));
        }
        cell.neighbor(roundedNormal).computeU(u_y2);
      }

      Vector<V, DESCRIPTOR::d> uWall{};
      if constexpr (movingWall) {
        uWall = cell.template getField<descriptors::VELOCITY>();
        for (int iD = 0; iD < DESCRIPTOR::d; iD++) {
          u_y2[iD] -= uWall[iD];
        }
      }

      V u_y2Norm = util::norm<DESCRIPTOR::d>(u_y2);
      if (u_y2Norm >= V(1e-7)) {
        V u_1Actual[DESCRIPTOR::d] = {0.};
        cell.computeU(u_1Actual);
        if constexpr (movingWall) {
          for (int iD = 0; iD < DESCRIPTOR::d; iD++) {
            u_1Actual[iD] -= uWall[iD];
          }
        }

        V u_y2normal = V(0);
        for (int iD = 0; iD < DESCRIPTOR::d; iD++) {
          u_y2normal += u_y2[iD] * normal[iD];
        }

        auto u_t = u_y2 - u_y2normal * normal;
        V u_tNorm = util::norm<DESCRIPTOR::d>(u_t);

        // Use a safe denominator that is never zero
        auto tangent = u_t / util::max(u_tNorm, eps);

        const V omega = parameters.template get<descriptors::OMEGA>();
        V nu = (V(1) / omega - 0.5) / descriptors::invCs2<V, DESCRIPTOR>();
        V u_2 = V(0);
        for (int iD = 0; iD < DESCRIPTOR::d; iD++) {
          u_2 += u_y2[iD] * tangent[iD];
        }

        V y2Norm = util::norm<DESCRIPTOR::d>(y2);
        V u_tau = cell.template getField<descriptors::U_TAU>();
        if (u_tau == V(0)) {
          u_tau = util::sqrt(nu * util::abs(u_2) / y2Norm);
        }

        V u_plusF = V(0);
        V abs_u2 = util::abs(u_2);

        if (abs_u2 > eps) {
          if constexpr (wallFunctionProfile == 0) {
            u_tau = powerLawWallFunction(nu, abs_u2, y2Norm, y1Norm)[0];
          } else {
            u_tau = newtonSpalding<10>(
              [=](V x) { return spaldingWallFunction<2>(x, abs_u2, y2Norm, nu); },
              [=](V x) { return spaldingDerivative<2>(x, abs_u2, y2Norm, nu); },
              u_tau,
              V(10) * eps
            );
          }
          cell.template setField<descriptors::U_TAU>(u_tau);

          V y_plus = y1Norm / nu * u_tau;
          if (y2Norm / nu * u_tau > V(2.5)) {
            V u_1t = V(0);
            for (int iD = 0; iD < DESCRIPTOR::d; iD++) {
              u_1t += u_1Actual[iD] * tangent[iD];
            }
            V abs_u1t = util::abs(u_1t);

            if constexpr (wallFunctionProfile == 0) {
              u_plusF = powerLawWallFunction(nu, abs_u2, y2Norm, y1Norm)[1];
            } else {
              V u_plus_init = abs_u1t / util::max(u_tau, eps);
              u_plusF = newtonSpalding<10>(
                [=](V x) { return spaldingWallFunction<1>(x, abs_u1t, y1Norm, nu, y_plus); },
                [=](V x) { return spaldingDerivative<1>(x, abs_u1t, y1Norm, nu); },
                u_plus_init,
                V(10) * eps
              );
            }

            if constexpr (useVanDriest) {
              V dUdY = (u_2 - u_plusF * u_tau) / util::max(util::norm<DESCRIPTOR::d>(y1y2), eps);
              V nuTurb = util::pow((V(0.4) * y1Norm * (V(1) - util::exp(-y_plus / V(26)))), 2) * util::abs(dUdY);
              cell.template setField<descriptors::VISCOSITY>(nu + nuTurb);
            }
          }
        }

        V u_1 = u_plusF * u_tau;
        auto ub_current = cell.template getField<descriptors::WMVELOCITY>();
        auto ub = u_1 * tangent;

        if constexpr (movingWall) {
          ub += uWall;

          // 0.5 under relaxation
          auto delta = V(0.5) * (ub - ub_current);

          //// no under relaxation
          //auto delta = ub - ub_current;

          V deltaNorm = util::norm<DESCRIPTOR::d>(delta);
          // continuous clamping (max allowed change per timestep 0.05)
          V scale = util::fmin(V(1), V(0.05) / util::max(deltaNorm, eps));
          ub = ub_current + delta * scale;
        }

        // max latticeU clamping
        V ubNorm = util::norm<DESCRIPTOR::d>(ub);
        ub *= util::fmin(V(1), V(0.3) / util::max(ubNorm, eps));

        V mask = V(u_tNorm > eps) * V(u_1 > eps);
        cell.template setField<descriptors::WMVELOCITY>(ub * mask + ub_current * (V(1) - mask));
      } else {
        cell.template setField<descriptors::WMVELOCITY>(0);
      }
    } else {
      cell.template setField<descriptors::WMVELOCITY>(0);
    }
  }
};

//======================================================================
// PostProcessor for fNeq from finite difference method (FDM) stress tensor for wall modelling at porous domains
//======================================================================

template<int rhoMethod>
struct TurbulentWallModelPorousFneqFDMPostProcessor {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 2;
  }

  template <typename CELL, typename V = typename CELL::value_t, typename DESCRIPTOR = typename CELL::descriptor_t>
  void interpolateDensity(CELL& cell, Vector<V,DESCRIPTOR::d> distance, V rho[1]) any_platform {
    Vector<V,DESCRIPTOR::d> floorV;
    for (int iD = 0; iD < DESCRIPTOR::d; ++iD) {
      floorV[iD] = util::floor(distance[iD]);
    }

    constexpr int nPts = (DESCRIPTOR::d == 2) ? 4 : 8;

    #pragma unroll
    for (int i = 0; i < nPts; ++i) {
      Vector<V, DESCRIPTOR::d> targetPt;
      V weight = V(1);
      for (int iD = 0; iD < DESCRIPTOR::d; ++iD) {
        targetPt[iD] = floorV[iD] + V((i >> iD) & 1);
        weight *= (V(1) - util::abs(distance[iD] - targetPt[iD]));
      }
      rho[0] += cell.neighbor(targetPt).computeRho() * weight;
    }
  }

  template <typename CELL, typename V = typename CELL::value_t>
  void apply(CELL& cell) any_platform {
    using DESCRIPTOR = typename CELL::descriptor_t;
    auto y1 = cell.template getField<descriptors::Y1>();
    V y1Norm = util::norm<DESCRIPTOR::d>(y1);

    if (y1Norm > V(1e-3)) {
      auto u = cell.template getField<descriptors::WMVELOCITY>();
      V uNorm = util::norm<DESCRIPTOR::d>(u);

      if (uNorm < V(1e-7)) {
        cell.template setField<descriptors::WMPOROSITY>(V(1));
        cell.template setField<collision::HYBRID>(V(1));
        if constexpr (rhoMethod == 1 || rhoMethod == 2) {
          cell.template setField<collision::HYBRID_RHO>(V(1));
        }
      } else {
        cell.template setField<descriptors::WMPOROSITY>(V(0));
        cell.template setField<collision::HYBRID>(V(0));

        if constexpr (DESCRIPTOR::d == 2) {
          using namespace olb::util::tensorIndices2D;
          V uxDx = V(0), uyDx = V(0), uxDy = V(0), uyDy = V(0);

          //----------- FDM in X-direction --------------------------------------------------//
          if (y1[0] > V(0) || cell.neighbor({-1,0}).template getField<descriptors::POROSITY>() == V(0)) {
            V uXp[DESCRIPTOR::d] {};
            cell.neighbor({1,0}).computeU(uXp);
            uxDx = uXp[0] - u[0];
            uyDx = uXp[1] - u[1];
          }
          else if (y1[0] < V(0) || cell.neighbor({1,0}).template getField<descriptors::POROSITY>() == V(0)) {
            V uXm[DESCRIPTOR::d] {};
            cell.neighbor({-1,0}).computeU(uXm);
            uxDx = u[0] - uXm[0];
            uyDx = u[1] - uXm[1];
          }
          if (util::abs(y1[0]) < V(1.e-8) &&
              cell.neighbor({-1,0}).template getField<descriptors::POROSITY>() != V(0) &&
              cell.neighbor({1,0}).template getField<descriptors::POROSITY>() != V(0)) {
            V uXp[DESCRIPTOR::d], uXm[DESCRIPTOR::d] {};
            cell.neighbor({1,0}).computeU(uXp);
            cell.neighbor({-1,0}).computeU(uXm);
            uxDx = V(0.5)*(uXp[0] - uXm[0]);
            uyDx = V(0.5)*(uXp[1] - uXm[1]);
          }

          //----------- FDM in Y-direction --------------------------------------------------//
          if (y1[1] > V(0) || cell.neighbor({0,-1}).template getField<descriptors::POROSITY>() == V(0)) {
            V uYp[DESCRIPTOR::d] {};
            cell.neighbor({0,1}).computeU(uYp);
            uxDy = uYp[0] - u[0];
            uyDy = uYp[1] - u[1];
          }
          else if (y1[1] < V(0) || cell.neighbor({0,1}).template getField<descriptors::POROSITY>() == V(0)) {
            V uYm[DESCRIPTOR::d] {};
            cell.neighbor({0,-1}).computeU(uYm);
            uxDy = u[0] - uYm[0];
            uyDy = u[1] - uYm[1];
          }
          if (util::abs(y1[1]) < V(1.e-8) &&
              cell.neighbor({0,-1}).template getField<descriptors::POROSITY>() != V(0) &&
              cell.neighbor({0,1}).template getField<descriptors::POROSITY>() != V(0)) {
            V uYp[DESCRIPTOR::d], uYm[DESCRIPTOR::d] {};
            cell.neighbor({0,1}).computeU(uYp);
            cell.neighbor({0,-1}).computeU(uYm);
            uxDy = V(0.5)*(uYp[0] - uYm[0]);
            uyDy = V(0.5)*(uYp[1] - uYm[1]);
          }

          //----------- TENSOR tensor with FDM --------------------------------------------------//
          V pi[util::TensorVal<DESCRIPTOR>::n] {V(0)};
          pi[xx] = uxDx;
          pi[yy] = uyDy;
          pi[xy] = V(0.5)*(uxDy + uyDx);
          cell.template setField<descriptors::TENSOR>(pi);
        }
        else {
          using namespace olb::util::tensorIndices3D;
          V uxDx = 0, uyDx = 0, uzDx = 0, uxDy = 0, uyDy = 0, uzDy = 0, uxDz = 0, uyDz = 0, uzDz = 0;
          //----------- FDM in X-direction --------------------------------------------------//
          if (y1[0] > V(0) || cell.neighbor({-1,0,0}).template getField<descriptors::POROSITY>() == V(0)) {
            V uXp[DESCRIPTOR::d] {};
            cell.neighbor({1,0,0}).computeU(uXp);
            uxDx = uXp[0] - u[0];
            uyDx = uXp[1] - u[1];
            uzDx = uXp[2] - u[2];
          }
          else if (y1[0] < V(0) || cell.neighbor({1,0,0}).template getField<descriptors::POROSITY>() == V(0)) {
            V uXm[DESCRIPTOR::d] {};
            cell.neighbor({-1,0,0}).computeU(uXm);
            uxDx = u[0] - uXm[0];
            uyDx = u[1] - uXm[1];
            uzDx = u[2] - uXm[2];
          }
          if (util::abs(y1[0]) < V(1.e-8) &&
              cell.neighbor({-1,0,0}).template getField<descriptors::POROSITY>() != V(0) &&
              cell.neighbor({1,0,0}).template getField<descriptors::POROSITY>() != V(0)) {
            V uXp[DESCRIPTOR::d], uXm[DESCRIPTOR::d] {};
            cell.neighbor({1,0,0}).computeU(uXp);
            cell.neighbor({-1,0,0}).computeU(uXm);
            uxDx = V(0.5)*(uXp[0] - uXm[0]);
            uyDx = V(0.5)*(uXp[1] - uXm[1]);
            uzDx = V(0.5)*(uXp[2] - uXm[2]);
          }

          //----------- FDM in Y-direction --------------------------------------------------//
          if (y1[1] > V(0) || cell.neighbor({0,-1,0}).template getField<descriptors::POROSITY>() == V(0)) {
            V uYp[DESCRIPTOR::d] {};
            cell.neighbor({0,1,0}).computeU(uYp);
            uxDy = uYp[0] - u[0];
            uyDy = uYp[1] - u[1];
            uzDy = uYp[2] - u[2];
          }
          else if (y1[1] < V(0) || cell.neighbor({0,1,0}).template getField<descriptors::POROSITY>() == V(0)) {
            V uYm[DESCRIPTOR::d] {};
            cell.neighbor({0,-1,0}).computeU(uYm);
            uxDy = u[0] - uYm[0];
            uyDy = u[1] - uYm[1];
            uzDy = u[2] - uYm[2];
          }
          if (util::abs(y1[1]) < V(1.e-8) &&
              cell.neighbor({0,-1,0}).template getField<descriptors::POROSITY>() != V(0) &&
              cell.neighbor({0,1,0}).template getField<descriptors::POROSITY>() != V(0)) {
            V uYp[DESCRIPTOR::d], uYm[DESCRIPTOR::d] {};
            cell.neighbor({0,1,0}).computeU(uYp);
            cell.neighbor({0,-1,0}).computeU(uYm);
            uxDy = V(0.5)*(uYp[0] - uYm[0]);
            uyDy = V(0.5)*(uYp[1] - uYm[1]);
            uzDy = V(0.5)*(uYp[2] - uYm[2]);
          }

          //----------- FDM in Z-direction --------------------------------------------------//
          if (y1[2] > V(0) || cell.neighbor({0,0,-1}).template getField<descriptors::POROSITY>() == V(0)) {
            V uZp[DESCRIPTOR::d] {};
            cell.neighbor({0,0,1}).computeU(uZp);
            uxDz = uZp[0] - u[0];
            uyDz = uZp[1] - u[1];
            uzDz = uZp[2] - u[2];
          }
          else if (y1[2] < V(0) || cell.neighbor({0,0,1}).template getField<descriptors::POROSITY>() == V(0)) {
            V uZm[DESCRIPTOR::d] {};
            cell.neighbor({0,0,-1}).computeU(uZm);
            uxDz = u[0] - uZm[0];
            uyDz = u[1] - uZm[1];
            uzDz = u[2] - uZm[2];
          }
          if (util::abs(y1[2]) < V(1.e-8) &&
              cell.neighbor({0,0,-1}).template getField<descriptors::POROSITY>() != V(0) &&
              cell.neighbor({0,0,1}).template getField<descriptors::POROSITY>() != V(0)) {
            V uZp[DESCRIPTOR::d], uZm[DESCRIPTOR::d] {};
            cell.neighbor({0,0,1}).computeU(uZp);
            cell.neighbor({0,0,-1}).computeU(uZm);
            uxDz = V(0.5)*(uZp[0] - uZm[0]);
            uyDz = V(0.5)*(uZp[1] - uZm[1]);
            uzDz = V(0.5)*(uZp[2] - uZm[2]);
          }

          //----------- TENSOR tensor with FDM --------------------------------------------------//
          V pi[util::TensorVal<DESCRIPTOR>::n] {V(0)};
          pi[xx] = uxDx;
          pi[yy] = uyDy;
          pi[zz] = uzDz;
          pi[xy] = V(0.5)*(uxDy + uyDx);
          pi[xz] = V(0.5)*(uxDz + uzDx);
          pi[yz] = V(0.5)*(uyDz + uzDy);

          cell.template setField<descriptors::TENSOR>(pi);
        }

        if constexpr (rhoMethod == 1) {
          Vector<V,DESCRIPTOR::d> roundedNormal;
          for (int iD = 0; iD < DESCRIPTOR::d; iD++) {
            roundedNormal[iD] = V(int(util::round(y1[iD] / y1Norm)));
          }
          V extRho = cell.neighbor(roundedNormal).computeRho();
          cell.template setField<descriptors::DENSITY>(extRho);
          cell.template setField<collision::HYBRID_RHO>(V(0));
        }
        else if constexpr (rhoMethod == 2) {
          cell.template setField<descriptors::DENSITY>(V(1.));
          cell.template setField<collision::HYBRID_RHO>(V(0));
        }
      }
    } else {
      cell.template setField<descriptors::WMPOROSITY>(V(1));
      cell.template setField<collision::HYBRID>(V(1));
      if constexpr (rhoMethod == 1 || rhoMethod == 2) {
        cell.template setField<collision::HYBRID_RHO>(V(1));
      }
    }
  }
};

}

}  // namespace optimized
#endif
