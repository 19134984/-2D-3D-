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

#ifndef GPU_SYCL_CONTEXT_HH
#define GPU_SYCL_CONTEXT_HH

#include "context.h"
#include "registry.hh"
#include "dynamics.hh"

namespace olb {

namespace gpu {

namespace sycl {

namespace concepts {

template <typename CONTAINER>
concept ConcreteFieldsContainer = requires (CONTAINER& container) {
  requires olb::concepts::BaseType<typename CONTAINER::value_t>;
  requires olb::concepts::Descriptor<typename CONTAINER::descriptor_t>;
  requires (CONTAINER::platform == Platform::GPU_SYCL);
};

}

/// Device accessor type of FIELD's column (dense pointer wrapper or sparse tile view)
template <typename T, typename DESCRIPTOR, typename FIELD>
using DeviceAccessorFor = typename detail::field_column_t<T,DESCRIPTOR,Platform::GPU_SYCL,FIELD>::Accessor;

/// Return pointers to on-device data for FIELD on lattice
template <concepts::ConcreteFieldsContainer CONTAINER, typename FIELD>
static auto getDeviceFieldPointer(CONTAINER& lattice) {
  using T = typename CONTAINER::value_t;
  using DESCRIPTOR = typename CONTAINER::descriptor_t;
  auto& fieldArray = lattice.template getField<FIELD>();
  // Every column (dense or sparse) exposes the same accessor interface
  return Vector<DeviceAccessorFor<T,DESCRIPTOR,FIELD>,DESCRIPTOR::template size<FIELD>()>([&](unsigned iD) {
    return fieldArray[iD].getDeviceAccessor();
  });
}

/// Return tuple of pointers to on-device data for multiple FIELDS on lattice
template <concepts::ConcreteFieldsContainer CONTAINER, typename... FIELDS>
static auto getDeviceFieldPointers(typename meta::list<FIELDS...>,
                                   CONTAINER& lattice) {
  return std::make_tuple(getDeviceFieldPointer<CONTAINER,FIELDS>(lattice)...);
}

/// Structure for passing pointers to on-device data into SYCL kernels
template <typename T, typename DESCRIPTOR>
class DeviceContext {
private:
  const std::size_t _nCells;

  /// Pointers to descriptor-declared field data
  decltype(getDeviceFieldPointers(typename DESCRIPTOR::fields_t(),
                                  std::declval<ConcreteBlockD<T,DESCRIPTOR,Platform::GPU_SYCL>&>()))
  _staticFieldsD;

  /// Pointer to type-erased field data
  void*** _customFieldsD;

  Vector<typename CyclicColumn<T>::Accessor,DESCRIPTOR::q> _populations;

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  template <concepts::ConcreteFieldsContainer CONTAINER>
  DeviceContext(CONTAINER& lattice) :
    _nCells(lattice.getNcells()),
    _staticFieldsD(getDeviceFieldPointers(typename DESCRIPTOR::fields_t(), lattice)),
    _customFieldsD(lattice.getDataRegistry().deviceData()),
    _populations([&lattice](unsigned iPop) -> typename CyclicColumn<T>::Accessor {
      return lattice.template getField<descriptors::POPULATION>()[iPop].getDeviceAccessor();
    })
  {
    static_assert(std::is_same_v<typename CONTAINER::value_t, T>,
                  "CONTAINER must share value_t");
    static_assert(std::is_same_v<typename CONTAINER::descriptor_t, DESCRIPTOR>,
                  "CONTAINER must share descriptor_t");
  }

  DeviceContext(ConcreteBlockD<T,DESCRIPTOR,Platform::GPU_SYCL>& lattice) :
    _nCells(lattice.getNcells()),
    _staticFieldsD(getDeviceFieldPointers(typename DESCRIPTOR::fields_t(), lattice)),
    _customFieldsD(lattice.getDataRegistry().deviceData())
  { }

  std::size_t getNcells() const {
    return _nCells;
  }

  CyclicColumn<T>::Accessor getPopulation(unsigned iPop) const {
    return _populations[iPop];
  }

  template <typename FIELD>
  auto getField() const {
    // Both paths yield a per-dimension array of the field's device accessor (dense or sparse)
    if constexpr (DESCRIPTOR::template provides<FIELD>()) {
      return std::get<(DESCRIPTOR::fields_t::template index<FIELD>())>(_staticFieldsD).data();
    } else {
      return reinterpret_cast<DeviceAccessorFor<T,DESCRIPTOR,FIELD>*>(
        _customFieldsD[gpu::sycl::field_type_index<FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_SYCL>, Array<FIELD>>.get()]);
    }
  }

};

/// Pointer to row of a D-dimensional field
template <typename T, typename DESCRIPTOR, typename FIELD>
class FieldPtr : public ScalarVector<typename FIELD::template value_type<T>,
                                     DESCRIPTOR::template size<FIELD>(),
                                     FieldPtr<T,DESCRIPTOR,FIELD>> {
private:
  const DeviceContext<T,DESCRIPTOR>& _data;
  std::size_t _index;

  friend typename ScalarVector<typename FIELD::template value_type<T>,
                               DESCRIPTOR::template size<FIELD>(),
                               FieldPtr>::type;

public:
  FieldPtr(const DeviceContext<T,DESCRIPTOR>& data, std::size_t index) :
    _data(data),
    _index(index) { }

  FieldPtr(FieldPtr<T,DESCRIPTOR,FIELD>&& rhs) :
    _data(rhs._data),
    _index(rhs._index) { }

  const typename FIELD::template value_type<T>* getComponentPointer(unsigned iDim) const
  {
    return &(_data.template getField<FIELD>()[iDim][_index]);
  }

  typename FIELD::template value_type<T>* getComponentPointer(unsigned iDim)
  {
    return &(_data.template getField<FIELD>()[iDim][_index]);
  }

  template <typename U, typename IMPL>
  FieldPtr<T,DESCRIPTOR,FIELD>& operator=(const GenericVector<U,DESCRIPTOR::template size<FIELD>(),IMPL>& rhs)
  {
    for (unsigned iDim=0; iDim < DESCRIPTOR::template size<FIELD>(); ++iDim) {
      this->operator[](iDim) = rhs[iDim];
    }
    return *this;
  }

};

/// Device-side implementation of the data-only Cell concept for collision steps
template <typename T, typename DESCRIPTOR>
class DataOnlyCell {
protected:
  const DeviceContext<T,DESCRIPTOR>& _data;
  CellID _iCell;

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  DataOnlyCell(const DeviceContext<T,DESCRIPTOR>& data, CellID iCell) :
    _data(data),
    _iCell(iCell)
  { }

  void setCellId(std::size_t iCell) {
    _iCell = iCell;
  }

  CellID getCellId() const {
    return _iCell;
  }

  value_t& operator[](unsigned int iPop) {
    return _data.getPopulation(iPop)[_iCell];
  }

  const value_t& operator[](unsigned int iPop) const {
    return _data.getPopulation(iPop)[_iCell];
  }

  template <typename FIELD>
  auto getField() const {
    if constexpr (std::is_same_v<FIELD, descriptors::POPULATION>) {
      // SYCL has no VMM ring, so the raw field pointer is unshifted. Read
      // populations through the cyclic accessor.
      return FieldD<value_t,descriptor_t,FIELD>([&](unsigned iPop) {
        return _data.getPopulation(iPop)[_iCell];
      });
    } else {
      auto fieldArray = _data.template getField<FIELD>();
      if constexpr (descriptor_t::template size<FIELD>() == 1) {
        return fieldArray[0].get(_iCell);
      } else {
        return FieldD<value_t,descriptor_t,FIELD>([&](unsigned iD) {
          return fieldArray[iD].get(_iCell);
        });
      }
    }
  }

  template <typename FIELD>
  typename FIELD::template value_type<T> getFieldComponent(unsigned iD) const {
    auto fieldArray = _data.template getField<FIELD>();
    return fieldArray[iD].get(_iCell);
  }

  template <typename FIELD>
  void setFieldComponent(unsigned iD, typename FIELD::template value_type<T> value) {
    auto fieldArray = _data.template getField<FIELD>();
    fieldArray[iD][_iCell] = value;
  }

  template <typename FIELD>
  FieldPtr<T,DESCRIPTOR,FIELD> getFieldPointer() {
    return FieldPtr<T,DESCRIPTOR,FIELD>(_data, _iCell);
  }

  template <typename FIELD>
  void setField(FieldD<value_t,descriptor_t,FIELD>&& value) {
    auto fieldArray = _data.template getField<FIELD>();
    for (unsigned iD=0; iD < descriptor_t::template size<FIELD>(); ++iD) {
      fieldArray[iD][_iCell] = value[iD];
    }
  }

  template <typename FIELD>
  void setField(const FieldD<value_t,descriptor_t,FIELD>& value) {
    auto fieldArray = _data.template getField<FIELD>();
    for (unsigned iD=0; iD < descriptor_t::template size<FIELD>(); ++iD) {
      fieldArray[iD][_iCell] = value[iD];
    }
  }

};

template <typename T, typename DESCRIPTOR>
class Cell;

/// Device-side view of a block lattice
/**
 * Used for non-local operators such as post processors
 **/
template <typename T, typename DESCRIPTOR>
class DeviceBlockLattice final : public DeviceContext<T,DESCRIPTOR> {
private:
  LatticeR<DESCRIPTOR::d> _core;
  LatticeR<DESCRIPTOR::d> _projection;
  int _padding;
  Vector<T,DESCRIPTOR::d> _origin;
  T _deltaR;

public:
  DeviceBlockLattice(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_SYCL>& lattice) :
    DeviceContext<T,DESCRIPTOR>(lattice),
    _core(lattice.getExtent()),
    _padding{lattice.getPadding()},
    _origin(lattice.getOrigin()),
    _deltaR(lattice.getDeltaR())
  {
    auto size = lattice.getExtent() + 2*_padding;
    if constexpr (DESCRIPTOR::d == 3) {
      _projection = {size[1]*size[2], size[2], 1};
    } else {
      _projection = {size[1], 1};
    }
  }

  CellID getCellId(LatticeR<DESCRIPTOR::d> loc) const {
    return (loc + _padding) * _projection;
  }

  CellDistance getNeighborDistance(LatticeR<DESCRIPTOR::d> dir) const {
    return dir * _projection;
  }

  Cell<T,DESCRIPTOR> get(CellID iCell) {
    return Cell<T,DESCRIPTOR>(*this, iCell);
  }

  Cell<T,DESCRIPTOR> get(LatticeR<DESCRIPTOR::d> loc) {
    return get(getCellId(loc));
  }

  LatticeR<DESCRIPTOR::d> getLatticeR(CellID iCell) const {
    if constexpr (DESCRIPTOR::d == 3) {
      const std::int32_t iX = iCell / _projection[0];
      const std::int32_t remainder = iCell % _projection[0];
      const std::int32_t iY = remainder / _projection[1];
      const std::int32_t iZ = remainder % _projection[1];
      return LatticeR<DESCRIPTOR::d>{iX - _padding, iY - _padding, iZ - _padding};
    } else {
      const std::int32_t iX = iCell / _projection[0];
      const std::int32_t iY = iCell % _projection[0];
      return LatticeR<DESCRIPTOR::d>{iX - _padding, iY - _padding};
    }
  };

  Vector<T,DESCRIPTOR::d> getPhysR(CellID iCell) const {
    return _origin + getLatticeR(iCell)*_deltaR;
  }

  bool isInside(LatticeR<DESCRIPTOR::d> loc) const {
    return loc >= -_padding && loc < _core + _padding;
  }

  bool isPadding(LatticeR<DESCRIPTOR::d> loc) const {
    return isInside(loc) && !(loc >= 0 && loc < _core);
  }

  bool isPadding(CellID iCell) const {
    return isPadding(getLatticeR(iCell));
  }

};

/// Device-side implementation of the Cell concept for post processors
/**
 * Adds neighborhood and dynamically-dispatched momenta access to DataOnlyCell
 **/
template <typename T, typename DESCRIPTOR>
class Cell final : public DataOnlyCell<T,DESCRIPTOR> {
protected:
  DeviceBlockLattice<T,DESCRIPTOR>& getBlock() {
    return static_cast<DeviceBlockLattice<T,DESCRIPTOR>&>(const_cast<DeviceContext<T,DESCRIPTOR>&>(this->_data));
  }

  const DeviceBlockLattice<T,DESCRIPTOR>& getBlock() const {
    return static_cast<const DeviceBlockLattice<T,DESCRIPTOR>&>(this->_data);
  }

public:
  Cell(DeviceBlockLattice<T,DESCRIPTOR>& data, CellID iCell) :
    DataOnlyCell<T,DESCRIPTOR>(data, iCell)
  { }

  Cell<T,DESCRIPTOR> neighbor(LatticeR<DESCRIPTOR::d> offset) {
    return {getBlock(), CellID(this->_iCell + getBlock().getNeighborDistance(offset))};
  }

  //Dynamics<T,DESCRIPTOR>& getDynamics() {
  //  return *this->template getField<gpu::sycl::DYNAMICS<T,DESCRIPTOR>>();
  //}

  std::size_t getCellId() const {
    return this->_iCell;
  }

  void setCellId(std::size_t iCell) {
    this->_iCell = iCell;
  }

  void setLatticeR(LatticeR<DESCRIPTOR::d> latticeR) {
    this->_iCell = getBlock().getCellId(latticeR);
  }

  LatticeR<DESCRIPTOR::d> getLatticeR() const {
    return getBlock().getLatticeR(this->_iCell);
  }

  Vector<T,DESCRIPTOR::d> getPhysR() const {
    return getBlock().getPhysR(this->_iCell);
  }

  bool isPadding() const {
    return getBlock().isPadding(this->_iCell);
  }

  T computeRho() {
    return this->template getField<fields::moments::RHO>();
  }
  void computeU(T* u) {
    auto uD = this->template getField<fields::moments::U>();
    for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
      u[iD] = uD[iD];
    }
  }
  void computeRhoU(T& rho, T* u) {
    rho = computeRho();
    computeU(u);
  }
  void computeEquilibrium(T& rho, T* u, T* fEq) {
    typename ThirdOrderRLBdynamics<T,DESCRIPTOR>::EquilibriumF().compute(*this, rho, u, fEq);
  }
  T getOmegaOrFallback(T fallback) {
    return fallback;
  };

  void computeJ(T* j) {
    auto jD = this->template getField<fields::moments::J>();
    for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
      j[iD] = jD[iD];
    }
  }
  void computeStress(T* pi) {
    auto piD = this->template getField<fields::moments::PI>();
    for (unsigned iPi=0; iPi < fields::moments::PI::template size<DESCRIPTOR>(); ++iPi) {
      pi[iPi] = piD[iPi];
    }
  }
  void computeAllMomenta(T& rho, T* u, T* pi) {
    rho = computeRho();
    computeU(u);
    computeStress(pi);
  }

  //void defineRho(T& rho) {
  //  return getDynamics().defineRho(*this, rho);
  //}
  //void defineU(T* u) {
  //  return getDynamics().defineU(*this, u);
  //}
  //void defineRhoU(T& rho, T* u) {
  //  return getDynamics().defineRhoU(*this, rho, u);
  //}
  //void defineAllMomenta(T& rho, T* u, T* pi) {
  //  return getDynamics().defineAllMomenta(*this, rho, u, pi);
  //}

  void definePopulations(const T* f) {
    for (int iPop=0; iPop < descriptors::q<DESCRIPTOR>(); ++iPop) {
      this->operator[](iPop) = f[iPop];
    }
  }

  //void iniEquilibrium(T rho, T* u) {
  //  getDynamics().iniEquilibrium(*this, rho, u);
  //}
  //void iniRegularized(T rho, T* u, T* pi) {
  //  getDynamics().iniRegularized(*this, rho, u, pi);
  //}
  //void inverseShiftRhoU(T& rho, T* u) {
  //  getDynamics().inverseShiftRhoU(*this, rho, u);
  //}

};

}

}

template <typename T, typename DESCRIPTOR, typename PARAMETERS>
ConcreteParametersD<T,DESCRIPTOR,Platform::GPU_SYCL,PARAMETERS>::ConcreteParametersD(std::size_t):
  _deviceParameters{gpu::sycl::device::malloc<ParametersD>(1)},
  parameters{}
{
  gpu::sycl::device::copyToDevice(&parameters,
                                  _deviceParameters.get(),
                                  sizeof(ParametersD));
}

template <typename T, typename DESCRIPTOR, typename PARAMETERS>
void ConcreteParametersD<T,DESCRIPTOR,Platform::GPU_SYCL,PARAMETERS>::setProcessingContext(
  ProcessingContext context)
{
  if (context == ProcessingContext::Simulation) {
    gpu::sycl::device::copyToDevice(&parameters,
                                    _deviceParameters.get(),
                                    sizeof(ParametersD));
  }
}

template<typename T, typename DESCRIPTOR, typename PARAMETERS>
std::size_t ConcreteParametersD<T,DESCRIPTOR,Platform::GPU_SYCL,PARAMETERS>::getNblock() const
{
  return decltype(parameters)::fields_t::size;
}

template<typename T, typename DESCRIPTOR, typename PARAMETERS>
std::size_t ConcreteParametersD<T,DESCRIPTOR,Platform::GPU_SYCL,PARAMETERS>::getSerializableSize() const
{
  std::size_t size = 0;
  decltype(parameters)::fields_t::for_each([&size](auto field) {
    using field_t = typename decltype(field)::type;
    size += FieldD<T,DESCRIPTOR,field_t>{}.getSerializableSize();
  });
  return size;
}

template<typename T, typename DESCRIPTOR, typename PARAMETERS>
bool* ConcreteParametersD<T,DESCRIPTOR,Platform::GPU_SYCL,PARAMETERS>::getBlock(
  std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;
  decltype(parameters)::fields_t::for_each([&](auto field) {
    using field_t = typename decltype(field)::type;
    if constexpr (DESCRIPTOR::template size<field_t>() == 1) {
      registerVar(iBlock, sizeBlock, currentBlock, dataPtr,
                  parameters.template get<field_t>(), loadingMode);
    } else {
      registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr,
                                      parameters.template get<field_t>(), loadingMode);
    }
  });
  return dataPtr;
}

}

namespace sycl {

template <typename T, typename DESCRIPTOR>
struct is_device_copyable<olb::gpu::sycl::DeviceContext<T, DESCRIPTOR>> : std::true_type {};

template <typename T, typename DESCRIPTOR>
struct is_device_copyable<olb::gpu::sycl::DeviceBlockLattice<T, DESCRIPTOR>> : std::true_type {};

template <typename T, unsigned D>
struct is_device_copyable<olb::Vector<T, D>> : std::true_type {};

}

#endif
