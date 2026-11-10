/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Adrian Kummerlaender
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

#ifndef GPU_SYCL_DYNAMICS_HH
#define GPU_SYCL_DYNAMICS_HH

#include "dynamics.h"
#include "dynamics/lbm.h"
#include "dynamics/interface.h"

namespace olb {

namespace gpu {

namespace sycl {

template <typename T, typename DESCRIPTOR> class DeviceContext;
template <typename T, typename DESCRIPTOR> class DataOnlyCell;

/// Virtual interface for device-side dynamically-dispatched dynamics access
template <typename T, typename DESCRIPTOR>
struct Dynamics {
  //SYCL_EXTERNAL virtual CellStatistic<T> collide(DeviceContext<T,DESCRIPTOR> lattice, CellID iCell)  = 0;

  //SYCL_EXTERNAL T computeRho (DataOnlyCell<T,DESCRIPTOR>& cell) {
  //  return cell.template getField<fields::moments::RHO>();
  //}
  //SYCL_EXTERNAL void computeU(DataOnlyCell<T,DESCRIPTOR>& cell, T* u) {
  //  auto uD = cell.template getField<fields::moments::U>();
  //  for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
  //    u[iD] = uD[iD];
  //  }
  //}
  //SYCL_EXTERNAL void computeRhoU(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u) {
  //  rho = computeRho(cell);
  //  computeU(cell, u);
  //}

  //SYCL_EXTERNAL void computeEquilibrium(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* fEq) {
  ////  typename DYNAMICS::EquilibriumF().compute(cell, rho, u, fEq);
  //  typename ThirdOrderRLBdynamics<T,DESCRIPTOR>::EquilibriumF().compute(cell, rho, u, fEq);
  //}

  //SYCL_EXTERNAL T getOmegaOrFallback(T fallback) {
  //  return fallback;
  //};

  //void computeJ   (DataOnlyCell<T,DESCRIPTOR>& cell,         T* j) { }
  //void computeStress    (DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* pi) { }
  //void computeAllMomenta(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* pi) { }

  //SYCL_EXTERNAL virtual void defineRho       (DataOnlyCell<T,DESCRIPTOR>& cell, T& rho             )  = 0;
  //SYCL_EXTERNAL virtual void defineU         (DataOnlyCell<T,DESCRIPTOR>& cell,         T* u       )  = 0;
  //SYCL_EXTERNAL virtual void defineRhoU      (DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u       )  = 0;
  //SYCL_EXTERNAL virtual void defineAllMomenta(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* pi)  = 0;

  //SYCL_EXTERNAL virtual void computeEquilibrium(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* fEq)  = 0;

  //void iniEquilibrium(DataOnlyCell<T,DESCRIPTOR>& cell, T rho, T* u)  {
  //  T fEq[DESCRIPTOR::q] { };
  //  computeEquilibrium(cell, rho, u, fEq);
  //  for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
  //    cell[iPop] = fEq[iPop];
  //  }
  //}

  //void iniRegularized(DataOnlyCell<T,DESCRIPTOR>& cell,
  //                    T rho, T* u, T* pi)  {
  //  iniEquilibrium(cell, rho, u);
  //  for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
  //    cell[iPop] += equilibrium<DESCRIPTOR>::template fromPiToFneq<T>(iPop, pi);
  //  }
  //}

  //SYCL_EXTERNAL virtual void inverseShiftRhoU(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u)   = 0;

};


/// Implementation of gpu::sycl::Dynamics for concrete DYNAMICS
template <typename T, typename DESCRIPTOR, typename DYNAMICS>
class ConcreteDynamics final : public Dynamics<T,DESCRIPTOR> {
private:
  ParametersOfOperatorD<T,DESCRIPTOR,DYNAMICS>* _parameters;

public:
  ConcreteDynamics(ParametersOfOperatorD<T,DESCRIPTOR,DYNAMICS>* parameters) :
    _parameters{parameters} {
  }

  //SYCL_EXTERNAL T getOmegaOrFallback(T fallback) override {
  //  if constexpr (DYNAMICS::parameters::template contains<descriptors::OMEGA>()) {
  //    return _parameters->template get<descriptors::OMEGA>();
  //  } else {
  //    return fallback;
  //  }
  //}

  //SYCL_EXTERNAL CellStatistic<T> collide(DeviceContext<T,DESCRIPTOR> lattice, CellID iCell) override  {
  //  DataOnlyCell<T,DESCRIPTOR> cell(lattice, iCell);
  //  return DYNAMICS().collide(cell, *_parameters);
  //}

  //SYCL_EXTERNAL T computeRho(DataOnlyCell<T,DESCRIPTOR>& cell) override  {
  //  return typename DYNAMICS::MomentaF().computeRho(cell);
  //}
  //SYCL_EXTERNAL void computeU(DataOnlyCell<T,DESCRIPTOR>& cell, T* u) override  {
  //  typename DYNAMICS::MomentaF().computeU(cell, u);
  //}
  //SYCL_EXTERNAL void computeJ(DataOnlyCell<T,DESCRIPTOR>& cell, T* j) override  {
  //  typename DYNAMICS::MomentaF().computeJ(cell, j);
  //}
  //SYCL_EXTERNAL void computeRhoU(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u) override  {
  //  typename DYNAMICS::MomentaF().computeRhoU(cell, rho, u);
  //}
  //SYCL_EXTERNAL void computeStress(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* pi) override  {
  //  typename DYNAMICS::MomentaF().computeStress(cell, rho, u, pi);
  //}
  //SYCL_EXTERNAL void computeAllMomenta(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* pi) override  {
  //  typename DYNAMICS::MomentaF().computeAllMomenta(cell, rho, u, pi);
  //}

  //SYCL_EXTERNAL T getOmegaOrFallback(T fallback) override  {
  //  if constexpr (DYNAMICS::parameters::template contains<descriptors::OMEGA>()) {
  //    return _parameters->template get<descriptors::OMEGA>();
  //  } else {
  //    return fallback;
  //  }
  //}

  //SYCL_EXTERNAL void computeEquilibrium(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* fEq) override  {
  //  typename DYNAMICS::EquilibriumF().compute(cell, rho, u, fEq);
  //}

  //SYCL_EXTERNAL void defineRho(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho) override  {
  //  typename DYNAMICS::MomentaF().defineRho(cell, rho);
  //}

  //SYCL_EXTERNAL void defineU(DataOnlyCell<T,DESCRIPTOR>& cell, T* u) override  {
  //  typename DYNAMICS::MomentaF().defineU(cell, u);
  //}

  //SYCL_EXTERNAL void defineRhoU(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u) override  {
  //  typename DYNAMICS::MomentaF().defineRhoU(cell, rho, u);
  //}

  //SYCL_EXTERNAL void defineAllMomenta(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* pi) override  {
  //  typename DYNAMICS::MomentaF().defineAllMomenta(cell, rho, u, pi);
  //}

  //SYCL_EXTERNAL void inverseShiftRhoU(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u) override  {
  //  typename DYNAMICS::MomentaF().inverseShiftRhoU(cell, rho, u);
  //}
};


/// Last node in a MaskedDynamics chain in kernel::call_operators
struct DynamicDispatchCollision {
  template <typename T, typename DESCRIPTOR>
  bool operator()(DeviceContext<T,DESCRIPTOR> lattice, CellID iCell) const {
    //if (auto* collisionO = lattice.template getField<DYNAMICS<T,DESCRIPTOR>>()[0][iCell]) {
    //  collisionO->collide(lattice, iCell);
    //  return true;
    //}
    return false;
  }

  template <typename T, typename DESCRIPTOR>
  bool operator()(DeviceContext<T,DESCRIPTOR> lattice, CellID iCell, CellStatistic<T>& statistic) const {
    //if (auto* collisionO = lattice.template getField<DYNAMICS<T,DESCRIPTOR>>()[0][iCell]) {
    //  statistic = collisionO->collide(lattice, iCell);
    //  return true;
    //}
    return false;
  }

};


}

}


}

#endif
