/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef CORE_PLATFORM_CPU_CELL_H
#define CORE_PLATFORM_CPU_CELL_H

#include <utility>

#include "dynamics/lbm.h"
#include "core/platform/convertingRef.h"

namespace olb {

template <typename T, typename DESCRIPTOR, Platform PLATFORM>
class ConcreteBlockD;

/// Implementations of CPU specifics
namespace cpu {

template <typename T, typename DESCRIPTOR, Platform PLATFORM>
class Cell;

template <typename T, typename DESCRIPTOR, Platform PLATFORM>
class DataOnlyCell;

/// Virtual interface for dynamically-dispatched dynamics access on CPU targets
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
struct Dynamics {
  virtual ~Dynamics() { }

  virtual CellStatistic<T> collide(DataOnlyCell<T,DESCRIPTOR,PLATFORM>& cell) = 0;

  virtual T    computeRho      (Cell<T,DESCRIPTOR,PLATFORM>& cell              ) = 0;
  virtual void computeU        (Cell<T,DESCRIPTOR,PLATFORM>& cell,         T* u) = 0;
  virtual void computeJ        (Cell<T,DESCRIPTOR,PLATFORM>& cell,         T* j) = 0;
  virtual void computeRhoU     (Cell<T,DESCRIPTOR,PLATFORM>& cell, T& rho, T* u) = 0;
  virtual void defineRho       (Cell<T,DESCRIPTOR,PLATFORM>& cell, T& rho      ) = 0;
  virtual void defineU         (Cell<T,DESCRIPTOR,PLATFORM>& cell,         T* u) = 0;
  virtual void defineRhoU      (Cell<T,DESCRIPTOR,PLATFORM>& cell, T& rho, T* u) = 0;
  virtual void defineAllMomenta(Cell<T,DESCRIPTOR,PLATFORM>& cell, T& rho, T* u, T* pi) = 0;

  virtual void computeStress    (Cell<T,DESCRIPTOR,PLATFORM>& cell, T& rho, T* u, T* pi) = 0;
  virtual void computeAllMomenta(Cell<T,DESCRIPTOR,PLATFORM>& cell, T& rho, T* u, T* pi) = 0;

  virtual void computeEquilibrium(Cell<T,DESCRIPTOR,PLATFORM>& cell, T rho, T* u, T* fEq) = 0;

  virtual T getOmegaOrFallback(T fallback) = 0;

  void iniEquilibrium(Cell<T,DESCRIPTOR,PLATFORM>& cell, T rho, T* u) {
    T fEq[DESCRIPTOR::q] { };
    computeEquilibrium(cell, rho, u, fEq);
    for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = fEq[iPop];
    }
  };

  void iniRegularized(Cell<T,DESCRIPTOR,PLATFORM>& cell,
                      T rho, T* u, T* pi) {
    iniEquilibrium(cell, rho, u);
    for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] += equilibrium<DESCRIPTOR>::template fromPiToFneq<T>(iPop, pi);
    }
  }

  virtual void inverseShiftRhoU(Cell<T,DESCRIPTOR,PLATFORM>& cell, T& rho, T* u) = 0;

};

/// Row concept for concrete non-LBM lattices on CPU platforms
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
class Row {
private:
  ConcreteBlockD<T,DESCRIPTOR,PLATFORM>* _lattice;
  CellID _iCell;

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  Row():
    _lattice(nullptr),
    _iCell(0) { }

  Row(ConcreteBlockD<T,DESCRIPTOR,PLATFORM>& lattice, CellID iCell=0):
    _lattice(&lattice),
    _iCell(iCell) { }

  CellID getCellId() const {
    return _iCell;
  }

  LatticeR<DESCRIPTOR::d> getLatticeR() const {
    return _lattice->getLatticeR(_iCell);
  }

  bool isCore() const {
    return _lattice->isInsideCore(_lattice->getLatticeR(_iCell));
  }

  void setCellId(CellID iCell) {
    _iCell = iCell;
  }

  template <typename FIELD>
  auto getField() const {
    return _lattice->template getField<FIELD>().getRow(_iCell);
  }

  template <typename FIELD>
  void setField(const FieldD<T,DESCRIPTOR,FIELD>& v) {
    return _lattice->template getField<FIELD>().setRow(_iCell, v);
  }

  template <typename FIELD>
  auto getFieldPointer() {
    return _lattice->template getField<FIELD>().getRowPointer(_iCell);
  }

  template <typename FIELD>
  decltype(auto) getFieldComponent(unsigned iD) {
    if constexpr (detail::isSparseColumnBacked<T,DESCRIPTOR,PLATFORM,FIELD>()) {
      // Ensure the container exists (sparse allocation creates no tiles), then read
      // through the const column so a plain read never materialises a tile
      auto& field = _lattice->template getField<FIELD>();
      return std::as_const(field)[iD][_iCell];
    } else {
      return (_lattice->template getField<FIELD>()[iD][_iCell]);
    }
  }

  template <typename FIELD>
  void setFieldComponent(unsigned iD, typename FIELD::template value_type<T> value) {
    _lattice->template getField<FIELD>()[iD][_iCell] = value;
  }

};

/// CPU specific field mirroring BlockDynamicsMap
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
struct DYNAMICS : public descriptors::TYPED_FIELD_BASE<Dynamics<T,DESCRIPTOR,PLATFORM>*,1> { };

/// Restricted cell exposing only populations and associated fields
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
class DataOnlyCell {
private:
  ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>* _lattice;
  CellID _iCell;

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  DataOnlyCell():
    _lattice(nullptr),
    _iCell(0) { }

  DataOnlyCell(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& lattice, std::size_t iCell=0):
    _lattice(&lattice),
    _iCell(iCell) { }

  CellID getCellId() const {
    return _iCell;
  }

  void setCellId(CellID iCell) {
    _iCell = iCell;
  }

  decltype(auto) operator[](unsigned iPop) {
    return makeConvertingRef<T>(_lattice->template getField<descriptors::POPULATION>()[iPop][_iCell]);
  }

  template <typename FIELD>
  auto getField() const {
    return _lattice->template getField<FIELD>().getRow(_iCell);
  }

  template <typename FIELD>
  void setField(const FieldD<T,DESCRIPTOR,FIELD>& v) {
    return _lattice->template getField<FIELD>().setRow(_iCell, v);
  }

  template <typename FIELD>
  auto getFieldPointer() {
    return _lattice->template getField<FIELD>().getRowPointer(_iCell);
  }

  template <typename FIELD>
  decltype(auto) getFieldComponent(unsigned iD) {
    if constexpr (detail::isSparseColumnBacked<T,DESCRIPTOR,PLATFORM,FIELD>()) {
      auto& field = _lattice->template getField<FIELD>();
      return std::as_const(field)[iD][_iCell];
    } else {
      return (_lattice->template getField<FIELD>()[iD][_iCell]);
    }
  }

  template <typename FIELD>
  void setFieldComponent(unsigned iD, typename FIELD::template value_type<T> value) {
    _lattice->template getField<FIELD>()[iD][_iCell] = value;
  }

};

/// Cell concept for concrete block lattices on CPU platforms
/**
 * Used for generic operators and non-legacy post processors
 **/
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
class Cell {
private:
  ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>* _lattice;
  CellID _iCell;

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  Cell():
    _lattice(nullptr),
    _iCell(0) { }

  Cell(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& lattice, std::size_t iCell=0):
    _lattice(&lattice),
    _iCell(iCell) { }

  CellID getCellId() const {
    return _iCell;
  }

  void setCellId(CellID iCell) {
    _iCell = iCell;
  }

  void setLatticeR(LatticeR<DESCRIPTOR::d> latticeR) {
    _iCell = _lattice->getCellId(latticeR);
  }

  LatticeR<DESCRIPTOR::d> getLatticeR() const {
    return _lattice->getLatticeR(_iCell);
  }

  Vector<T,DESCRIPTOR::d> getPhysR() const {
    return _lattice->getPhysR(_iCell);
  }

  bool isCore() const {
    return _lattice->isInsideCore(_lattice->getLatticeR(_iCell));
  }

  bool isPadding() const {
    return _lattice->isPadding(_iCell);
  }

  decltype(auto) operator[](unsigned iPop) {
    return makeConvertingRef<T>(_lattice->template getField<descriptors::POPULATION>()[iPop][_iCell]);
  }

  template <typename FIELD>
  auto getField() const {
    return _lattice->template getField<FIELD>().getRow(_iCell);
  }

  template <typename FIELD>
  void setField(const FieldD<T,DESCRIPTOR,FIELD>& v) {
    return _lattice->template getField<FIELD>().setRow(_iCell, v);
  }

  template <typename FIELD>
  auto getFieldPointer() {
    return _lattice->template getField<FIELD>().getRowPointer(_iCell);
  }

  template <typename FIELD>
  decltype(auto) getFieldComponent(unsigned iD) {
    if constexpr (detail::isSparseColumnBacked<T,DESCRIPTOR,PLATFORM,FIELD>()) {
      // Ensure the container exists (sparse allocation creates no tiles), then read
      // through the const column so a plain read never materialises a tile
      auto& field = _lattice->template getField<FIELD>();
      return std::as_const(field)[iD][_iCell];
    } else {
      return (_lattice->template getField<FIELD>()[iD][_iCell]);
    }
  }

  template <typename FIELD>
  void setFieldComponent(unsigned iD, typename FIELD::template value_type<T> value) {
    _lattice->template getField<FIELD>()[iD][_iCell] = value;
  }

  Cell<T,DESCRIPTOR,PLATFORM> neighbor(LatticeR<DESCRIPTOR::d> offset) {
    return {*_lattice, static_cast<std::size_t>(_iCell + _lattice->getNeighborDistance(offset))};
  }

  Dynamics<T,DESCRIPTOR,PLATFORM>& getDynamics() {
    return *getFieldComponent<DYNAMICS<T,DESCRIPTOR,PLATFORM>>(0);
  }

  T computeRho() {
    if constexpr (codegenTargetFor(PLATFORM) == CodegenTarget::StoredMoments) {
      return getFieldComponent<fields::moments::RHO>(0);
    } else {
      return getDynamics().computeRho(*this);
    }
  }
  void computeU(T* u) {
    if constexpr (codegenTargetFor(PLATFORM) == CodegenTarget::StoredMoments) {
      for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
        u[iD] = getFieldComponent<fields::moments::U>(iD);
      }
    } else {
      getDynamics().computeU(*this, u);
    }
  }
  void computeJ(T* j) {
    if constexpr (codegenTargetFor(PLATFORM) == CodegenTarget::StoredMoments) {
      for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
        j[iD] = getFieldComponent<fields::moments::J>(iD);
      }
    } else {
      getDynamics().computeJ(*this, j);
    }
  }
  void computeRhoU(T& rho, T* u) {
    if constexpr (codegenTargetFor(PLATFORM) == CodegenTarget::StoredMoments) {
      rho = computeRho();
      computeU(u);
    } else {
      getDynamics().computeRhoU(*this, rho, u);
    }
  }
  void computeStress(T* pi) {
    if constexpr (codegenTargetFor(PLATFORM) == CodegenTarget::StoredMoments) {
      for (unsigned iPi=0; iPi < fields::moments::PI::template size<DESCRIPTOR>(); ++iPi) {
        pi[iPi] = getFieldComponent<fields::moments::PI>(iPi);
      }
    } else {
      T rho, u[DESCRIPTOR::d] { };
      getDynamics().computeRhoU(*this, rho, u);
      getDynamics().computeStress(*this, rho, u, pi);
    }
  }
  void computeAllMomenta(T& rho, T* u, T* pi) {
    if constexpr (codegenTargetFor(PLATFORM) == CodegenTarget::StoredMoments) {
      rho = computeRho();
      computeU(u);
      computeStress(pi);
    } else {
      getDynamics().computeAllMomenta(*this, rho, u, pi);
    }
  }

  void defineRho(T& rho) {
    getDynamics().defineRho(*this, rho);
  }
  void defineU(T* u) {
    getDynamics().defineU(*this, u);
  }
  void defineRhoU(T rho, T* u) {
    getDynamics().defineRhoU(*this, rho, u);
  }
  void defineAllMomenta(T rho, T* u, T* pi) {
    getDynamics().defineAllMomenta(*this, rho, u, pi);
  }
  void definePopulations(const T* f) {
    for (int iPop=0; iPop < descriptors::q<DESCRIPTOR>(); ++iPop) {
      operator[](iPop) = f[iPop];
    }
  }

  void iniEquilibrium(T rho, T* u) {
    getDynamics().iniEquilibrium(*this, rho, u);
  }
  void iniRegularized(T rho, T* u, T* pi) {
    getDynamics().iniRegularized(*this, rho, u, pi);
  }
  void inverseShiftRhoU(T& rho, T* u) {
    getDynamics().inverseShiftRhoU(*this, rho, u);
  }

  T getOmegaOrFallback(T fallback) {
    return getDynamics().getOmegaOrFallback(fallback);
  }

  void computeEquilibrium(T rho, T* u, T* fEq) {
    getDynamics().computeEquilibrium(*this, rho, u, fEq);
  }

};

template <typename T, typename DESCRIPTOR, Platform PLATFORM>
class PlainCell {
private:
  ConcreteBlockD<T,DESCRIPTOR,PLATFORM>* _lattice;
  CellID _iCell;

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  PlainCell():
    _lattice(nullptr),
    _iCell(0) { }

  PlainCell(ConcreteBlockD<T,DESCRIPTOR,PLATFORM>& lattice, std::size_t iCell=0):
    _lattice(&lattice),
    _iCell(iCell) { }

  CellID getCellId() const {
    return _iCell;
  }

  LatticeR<DESCRIPTOR::d> getLatticeR() const {
    return _lattice->getLatticeR(_iCell);
  }

  bool isCore() const {
    return _lattice->isInsideCore(_lattice->getLatticeR(_iCell));
  }

  void setCellId(CellID iCell) {
    _iCell = iCell;
  }

  decltype(auto) operator[](unsigned iPop) {
    return makeConvertingRef<T>(_lattice->template getField<descriptors::POPULATION>()[iPop][_iCell]);
  }

  template <typename FIELD>
  auto getField() const {
    return _lattice->template getField<FIELD>().getRow(_iCell);
  }

  template <typename FIELD>
  void setField(const FieldD<T,DESCRIPTOR,FIELD>& v) {
    return _lattice->template getField<FIELD>().setRow(_iCell, v);
  }

  template <typename FIELD>
  auto getFieldPointer() {
    return _lattice->template getField<FIELD>().getRowPointer(_iCell);
  }

  template <typename FIELD>
  decltype(auto) getFieldComponent(unsigned iD) {
    if constexpr (detail::isSparseColumnBacked<T,DESCRIPTOR,PLATFORM,FIELD>()) {
      // Ensure the container exists (sparse allocation creates no tiles), then read
      // through the const column so a plain read never materialises a tile
      auto& field = _lattice->template getField<FIELD>();
      return std::as_const(field)[iD][_iCell];
    } else {
      return (_lattice->template getField<FIELD>()[iD][_iCell]);
    }
  }

  Cell<T,DESCRIPTOR,PLATFORM> neighbor(LatticeR<DESCRIPTOR::d> offset) {
    return {*_lattice, _iCell + _lattice->getNeighborDistance(offset)};
  }

};

}

}

#endif
