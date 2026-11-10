/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  Generic version of the collision, which modifies the particle
 *  distribution functions, by Orestis Malaspinas.
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

#ifndef FD_BOUNDARIES_2D_H
#define FD_BOUNDARIES_2D_H

#include "core/postProcessing.h"

#include "core/operator.h"
#include "boundaryMomentsDemand.h"

namespace olb {

/**
* This class computes the skordos BC
* on a flat wall in 2D but with a limited number of terms added to the
* equilibrium distributions (i.e. only the Q_i : Pi term)
*/
template <typename T, typename DESCRIPTOR, int direction, int orientation>
class StraightFdBoundaryProcessor2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const { return 0; }

  static std::vector<MomentsDemand<DESCRIPTOR::d>> getMomentsDemand() {
    return planeFdMomentsDemand<DESCRIPTOR,direction,orientation>();
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform;

private:
  template <int deriveDirection, typename CELL, typename V=CELL::value_t>
  void interpolateGradients(CELL& blockLattice, V velDeriv[DESCRIPTOR::d]) const any_platform;
};

/// Typed variant of StraightFdBoundaryProcessor2D: obtains the center cell's
/// momenta/equilibrium from the statically known DYNAMICS and omega from the
/// operator parameters, so no device-side dynamics dispatch is required (GPU_SYCL).
template <typename T, typename DESCRIPTOR, typename DYNAMICS, int direction, int orientation>
class StraightFdBoundaryProcessorOnDynamics2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::list<descriptors::OMEGA>;

  int getPriority() const { return 0; }

  static std::vector<MomentsDemand<DESCRIPTOR::d>> getMomentsDemand() {
    return planeFdMomentsDemand<DESCRIPTOR,direction,orientation>();
  }

  template <concepts::DynamicCell CELL, concepts::Parameters PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform;

private:
  template <int deriveDirection, typename CELL, typename V=CELL::value_t>
  void interpolateGradients(CELL& blockLattice, V velDeriv[DESCRIPTOR::d]) const any_platform;
};

/// PostProcessors for the chemical potential boundary condition in the free energy model.
/// The chemical potentials on the boundary are set equal to the chemical potential on the
/// fluid cell normal to the boundary. This is necessary because the coupling between the
/// lattices requires the calculation of the gradient of the chemical potential.
///
/// It would be preferable if these were implemented as a lattice coupling that ran
/// between the chemical potential and force lattice couplings. However there is no
/// access to the discrete normals in lattice couplings.
template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y>
struct FreeEnergyInletMomentum2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform {

    cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y}).template getField<descriptors::CHEM_POTENTIAL>());

    T rho0 = cell.computeRho();
    T rho1 = cell.neighbor({-NORMAL_X,-NORMAL_Y}).computeRho();

    cell.template setField<descriptors::CHEM_POTENTIAL>(
      cell.template getField<descriptors::CHEM_POTENTIAL>() + (rho1 / rho0 - 1) / descriptors::invCs2<T,DESCRIPTOR>());

  }

};

template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y>
struct FreeEnergyInletOrderParameter2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform {

    cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y}).template getField<descriptors::CHEM_POTENTIAL>());

  }

};

template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y>
struct FreeEnergyOutletMomentum2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform {

    T rhoBoundaryNew, rhoBoundaryOld, rhoBulk, u[2];

    rhoBoundaryOld = cell.computeRho();

    cell.neighbor({-NORMAL_X,-NORMAL_Y}).computeRhoU(rhoBulk, u);

    T uPerp = 0;

    Vector<T,2> normalVec({NORMAL_X,NORMAL_Y});

    if (normalVec[0] == 0) {
      uPerp = normalVec[1] * u[1];
    } else if (normalVec[1] == 0) {
            uPerp = normalVec[0] * u[0];
    }

    rhoBoundaryNew = (rhoBoundaryOld + uPerp * rhoBulk) / (1. + uPerp);
    cell.defineRho(rhoBoundaryNew);

    cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y}).template getField<descriptors::CHEM_POTENTIAL>());

    cell.template setField<descriptors::CHEM_POTENTIAL>(
      cell.template getField<descriptors::CHEM_POTENTIAL>() + (rhoBulk / rhoBoundaryNew - 1) / descriptors::invCs2<T,DESCRIPTOR>());

  }

};

template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y>
struct FreeEnergyOutletOrderParameter2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform {

    T rhoBoundaryNew, rhoBoundaryOld, rhoBulk, u[2];

    rhoBoundaryOld = cell.computeRho();

    cell.neighbor({-NORMAL_X,-NORMAL_Y}).computeRhoU(rhoBulk, u);

    T uPerp = 0;

    Vector<T,2> normalVec({NORMAL_X,NORMAL_Y});

    if (normalVec[0] == 0) {
      uPerp = normalVec[1] * u[1];
    } else if (normalVec[1] == 0) {
            uPerp = normalVec[0] * u[0];
    }

    rhoBoundaryNew = (rhoBoundaryOld + uPerp * rhoBulk) / (1. + uPerp);
    cell.defineRho(rhoBoundaryNew);

    cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y}).template getField<descriptors::CHEM_POTENTIAL>());

  }

};


/// PostProcessor for the wetting boundary condition in the free energy model. This is
/// required to set rho on the boundary (using the denisty of the neighbouring cell in
/// direction of inwards facing normal at the boundary), as the coupling between the
/// lattices requires the calculation of a density gradient.
template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y>
struct FreeEnergyWallMomentumProcessor2D {

public:
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<olb::descriptors::ADDEND>;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform{

    auto addend = parameters.template get<descriptors::ADDEND>();

    T rhoBulk = cell.neighbor({-NORMAL_X,-NORMAL_Y}).computeRho();
    T rhoTmp = 0.;

    for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
      rhoTmp += cell[iPop];
    }

    T rhoBoundary = rhoBulk + addend;
    rhoBoundary -= rhoTmp;

    cell[0] = rhoBoundary - 1.;

    cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y}).template getField<descriptors::CHEM_POTENTIAL>());

    cell.template setField<descriptors::CHEM_POTENTIAL>(
      cell.template getField<descriptors::CHEM_POTENTIAL>() + (rhoBulk / rhoBoundary - 1) / descriptors::invCs2<T,DESCRIPTOR>());

  }

};

template<typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y>
struct FreeEnergyWallOrderParameterProcessor2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<olb::descriptors::ADDEND>;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform{

    auto addend = parameters.template get<descriptors::ADDEND>();

    T rhoBulk = cell.neighbor({-NORMAL_X,-NORMAL_Y}).computeRho();
    T rhoTmp = 0.;

    for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
      rhoTmp += cell[iPop];
    }

    T rhoBoundary = rhoBulk + addend;
    rhoBoundary -= rhoTmp;

    cell[0] = rhoBoundary - 1.;

    cell.template setField<descriptors::CHEM_POTENTIAL>(
    cell.neighbor({-NORMAL_X,-NORMAL_Y}).template getField<descriptors::CHEM_POTENTIAL>());

  }

};

/// PostProcessor for pressure / velocity outflow boundaries in the free energy model.
/// The density / order parameters are prescribed to the outflow nodes such that they
/// obey the local-velocity convective boundary condition given in Lou, Gou, Shi (2013).
template <typename T, typename DESCRIPTOR, int NORMAL_X, int NORMAL_Y>
class FreeEnergyConvectiveProcessor2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform;
};

/**
* This class computes a convection BC on a flat wall in 2D
*/
template <typename T, typename DESCRIPTOR, int direction, int orientation>
class StraightConvectionBoundaryProcessor2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  struct PREV_CELL
      : public descriptors::FIELD_BASE<
            util::populationsContributingToVelocity<DESCRIPTOR, direction,
                                                    -orientation>()
                .size()> {};

  int getPriority() const {
    return 0;
  }

  template <concepts::DynamicCell CELL>
  void initialize(CELL& cell) any_platform;

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform;
};

/**
* This class computes a slip BC in 2D
*/

template <typename T, typename DESCRIPTOR>
class SlipBoundaryProcessor2D : public LocalPostProcessor2D<T, DESCRIPTOR> {
public:
  SlipBoundaryProcessor2D(int x0_, int x1_, int y0_, int y1_,
                          int discreteNormalX_, int discreteNormalY_);
  int  extent() const override { return 0; }
  int  extent(int whichDirection) const override { return 0; }
  void process(BlockLattice<T, DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice<T, DESCRIPTOR>& blockLattice, int x0_,
                        int x1_, int y0_, int y1_) override;

private:
  int reflectionPop[DESCRIPTOR::q];
  int x0, x1, y0, y1;
};

template <typename T, typename DESCRIPTOR>
class SlipBoundaryProcessorGenerator2D
    : public PostProcessorGenerator2D<T, DESCRIPTOR> {
public:
  SlipBoundaryProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_,
                                   int discreteNormalX_, int discreteNormalY_);
  PostProcessor2D<T, DESCRIPTOR>*          generate() const override;
  PostProcessorGenerator2D<T, DESCRIPTOR>* clone() const override;

private:
  int discreteNormalX;
  int discreteNormalY;
};


/**
* This class computes the skordos BC in 2D on a convex
* corner but with a limited number of terms added to the
* equilibrium distributions (i.e. only the Q_i : Pi term)
*/
template <typename T, typename DESCRIPTOR, int xNormal, int yNormal>
class OuterVelocityCornerProcessor2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 1;
  }

  static std::vector<MomentsDemand<DESCRIPTOR::d>> getMomentsDemand() {
    return outerVelocityCornerMomentsDemand<DESCRIPTOR::d>({xNormal,yNormal});
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform;
};

/// Typed variant of OuterVelocityCornerProcessor2D: obtains the center cell's
/// momenta/equilibrium from the statically known DYNAMICS and omega from the
/// operator parameters, so no device-side dynamics dispatch is required (GPU_SYCL).
template <typename T, typename DESCRIPTOR, typename DYNAMICS, int xNormal, int yNormal>
class OuterVelocityCornerProcessorOnDynamics2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::list<descriptors::OMEGA>;

  int getPriority() const {
    return 1;
  }

  static std::vector<MomentsDemand<DESCRIPTOR::d>> getMomentsDemand() {
    return outerVelocityCornerMomentsDemand<DESCRIPTOR::d>({xNormal,yNormal});
  }

  template <concepts::DynamicCell CELL, concepts::Parameters PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform;
};

template <typename T, typename DESCRIPTOR, int xNormal, int yNormal>
class IncOuterVelocityCornerProcessor2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 1;
  }

  static std::vector<MomentsDemand<DESCRIPTOR::d>> getMomentsDemand() {
    return outerVelocityCornerMomentsDemand<DESCRIPTOR::d>({xNormal,yNormal});
  }

  template <concepts::DynamicCell CELL>
  void apply(CELL& cell) any_platform;
};

/**
* This class computes the convective boundary condition for
* the populations of a phase field solving LBE.
*/
template<typename T, typename DESCRIPTOR, int xNormal,int yNormal>
struct FlatConvectivePhaseFieldPostProcessorA2D {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::MAX_VELOCITY>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform{

    auto u_conv = parameters.template get<descriptors::MAX_VELOCITY>();
    Vector<int,DESCRIPTOR::d> normal{-xNormal, -yNormal};
    auto cell1 = cell.neighbor(normal);
    auto outlet_cell = cell.template getField<descriptors::CONV_POPS>();

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      outlet_cell[iPop] = (outlet_cell[iPop] + u_conv * cell1[iPop]) / ((T)1 + u_conv);
      cell[iPop] = outlet_cell[iPop];
    }
    cell.template setField<descriptors::CONV_POPS>(outlet_cell);
  }

};

/**
* This class computes the convective boundary condition for
* phi at ghost nodes for a phase field solving LBE.
*/
template<typename T, typename DESCRIPTOR, int xNormal,int yNormal>
struct FlatConvectivePhaseFieldPostProcessorB2D {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;
  using parameters = meta::list<descriptors::MAX_VELOCITY>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform{

    auto u_conv = parameters.template get<descriptors::MAX_VELOCITY>();
    Vector<int,DESCRIPTOR::d> normal{xNormal, yNormal};
    auto outlet_cell = cell.template getField<descriptors::CONV_POPS>();

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = outlet_cell[iPop];
    }

    auto phi = cell.template getField<descriptors::STATISTIC>();
    phi[0] = cell.computeRho();
    cell.template setField<descriptors::STATISTIC>(phi);
    auto phiGhost = cell.neighbor(normal).template getField<descriptors::STATISTIC>();
    phiGhost[0] = (phiGhost[0] + u_conv * phi[0]) / ((T)1 + u_conv);
    if (phiGhost[0] > 1.0001) phiGhost[0] = 1;
    if (phiGhost[0] < -0.0001) phiGhost[0] = 0;
    cell.neighbor(normal).template setField<descriptors::STATISTIC>(phiGhost);

    Vector<int,DESCRIPTOR::d> tangent{yNormal, -xNormal};
    if (cell.neighbor(tangent).template getField<descriptors::SCALAR>() == 2 ||
        cell.neighbor(tangent).template getField<descriptors::BOUNDARY>() == 2) {
      cell.neighbor(normal+tangent).template setField<descriptors::STATISTIC>(phiGhost);
    }
    tangent *= -1;
    if (cell.neighbor(tangent).template getField<descriptors::SCALAR>() == 2 ||
        cell.neighbor(tangent).template getField<descriptors::BOUNDARY>() == 2) {
      cell.neighbor(normal+tangent).template setField<descriptors::STATISTIC>(phiGhost);
    }
  }

};

/**
* This class computes the Neumann boundary condition for
* the populations of a phase field solving LBE.
*/
template<typename T, typename DESCRIPTOR, int xNormal,int yNormal>
struct FlatNeumannPhaseFieldPostProcessorA2D {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform{

    Vector<int,2> normal{xNormal, yNormal};
    Vector<int,2> oppositeNormal{-xNormal, -yNormal};
    auto phi = cell.template getField<descriptors::STATISTIC>();
    Vector<T,2> uVec = cell.template getField<descriptors::VELOCITY>();
    T dotProduct{};
    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      dotProduct += uVec[iD]*normal[iD];
    }

    if ( dotProduct > 0 ) {
      auto cell1 = cell.neighbor(oppositeNormal);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] = cell1[iPop];
      }
      phi[0] = cell.computeRho();
    } else {
      phi[0] = cell.neighbor(normal).template getField<descriptors::STATISTIC>()[0];
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] = equilibrium<DESCRIPTOR>::firstOrder(iPop, phi[0], uVec);
      }
    }

    if (phi[0] > 1.0001) phi[0] = T{1};
    if (phi[0] < -0.0001) phi[0] = T{0};
    cell.template setField<descriptors::STATISTIC>(phi);
  }

};

/**
* This class computes the Neumann boundary condition for
* phi at ghost nodes for a phase field solving LBE.
*/
template<typename T, typename DESCRIPTOR, int xNormal,int yNormal>
struct FlatNeumannPhaseFieldPostProcessorB2D {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform{

    Vector<int,DESCRIPTOR::d> normal{xNormal, yNormal};

    auto phi = cell.template getField<descriptors::STATISTIC>();
    auto phiGhost = cell.neighbor(normal).template getField<descriptors::STATISTIC>();
    phiGhost[0] = phi[0];
    if (phiGhost[0] > 1.0001) phiGhost[0] = T{1};
    if (phiGhost[0] < -0.0001) phiGhost[0] = T{0};
    cell.neighbor(normal).template setField<descriptors::STATISTIC>(phiGhost);

    for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
      T dotProduct{};
      for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
        dotProduct += descriptors::c<DESCRIPTOR>(iPop)[iD]*normal[iD];
      }
      if ( dotProduct == 0 &&
      (cell.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getField<descriptors::SCALAR>() >= 1.9 ||
       cell.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getField<descriptors::SCALAR>() <= 1e-6))
      {
        cell.neighbor(normal+descriptors::c<DESCRIPTOR>(iPop)).template setField<descriptors::STATISTIC>(phiGhost);
      }
    }
  }

};

/**
* Setting the ghost node phi statistic equal to the phi of the fluid node for a phase field solving LBE.
*/
template<typename T, typename DESCRIPTOR, int xNormal,int yNormal>
struct PhaseFieldGhostProcessor2D {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 1;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform{
    Vector<int,DESCRIPTOR::d> normal{xNormal, yNormal};

    auto phi = cell.template getField<descriptors::STATISTIC>();
    auto phiGhost = cell.neighbor(normal).template getField<descriptors::STATISTIC>();
    phiGhost[0] = phi[0];
    cell.neighbor(normal).template setField<descriptors::STATISTIC>(phiGhost);

    for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
      T dotProduct{};
      for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
        dotProduct += descriptors::c<DESCRIPTOR>(iPop)[iD]*normal[iD];
      }
      if ( dotProduct == 0 &&
      (cell.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getField<descriptors::SCALAR>() >= 1.9 ||
       cell.neighbor(descriptors::c<DESCRIPTOR>(iPop)).template getField<descriptors::SCALAR>() <= 1e-6))
      {
        cell.neighbor(normal+descriptors::c<DESCRIPTOR>(iPop)).template setField<descriptors::STATISTIC>(phiGhost);
      }
    }
  }

};

/**
* This class overcomes overwriting populations from streaming
* for convective boundary conditions.
*/
template<typename T, typename DESCRIPTOR>
struct ConvectivePostProcessor2D {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform{
    auto outlet_cell = cell.template getField<descriptors::CONV_POPS>();
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = outlet_cell[iPop];
    }
  }

};

/**
* Extrapolating momentum normal to the boundary from one node inside the fluid.
*/
template<typename T, typename DESCRIPTOR, int direction, int orientation>
struct MomentumExtrapolatePostProcessor2D {
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 0;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform{

    auto rho = cell.template getField<descriptors::RHO>();
    auto j = cell.template getField<descriptors::MOMENTUM>();

    Vector<int,DESCRIPTOR::d> insideNormal{0,0};
    insideNormal[direction] = -orientation;
    Vector<int,DESCRIPTOR::d> normal{0,0};
    normal[direction] = orientation;
    T u[DESCRIPTOR::d] {};
    cell.computeU(u);
    Vector<T,2> uVec{};
    T dotProduct{};
    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      uVec[iD] = u[iD];
      dotProduct += uVec[iD]*normal[iD];
    }

    if ( dotProduct > 0 ) {
      auto cell1 = cell.neighbor(insideNormal);
      T u1[DESCRIPTOR::d];
      cell1.computeU(u1);
      for (int iDim = 0; iDim < DESCRIPTOR::d; ++iDim) {
        j[iDim] = rho*u1[iDim];
      }
    } else {
      for (int iDim = 0; iDim < DESCRIPTOR::d; ++iDim) {
        j[iDim] = T(0.);
      }
      j[direction] = rho*uVec[direction];
    }

    cell.template setField<descriptors::MOMENTUM>(j);
  }

};


} // namespace olb

#endif
