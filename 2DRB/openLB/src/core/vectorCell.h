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

#ifndef CORE_VECTOR_CELL_H
#define CORE_VECTOR_CELL_H

#include "scalarVector.h"
#include "descriptor/descriptor.h"

namespace olb {

template <typename CELL> requires concepts::VectorLatticeDescriptor<typename CELL::descriptor_t>
class VectorCell {
private:
  CELL _cell;

  /// Proxy for providing array view of dispersed components
  class PopulationPtr : public ScalarVector<typename CELL::value_t,
                                            CELL::descriptor_t::n,
                                            PopulationPtr> {
  private:
    CELL& _cell;
    unsigned _iPopV;

    friend typename ScalarVector<typename CELL::value_t,
                                 CELL::descriptor_t::n,
                                 PopulationPtr>::type;

  public:
    PopulationPtr(CELL& cell, unsigned iPopV) any_platform:
      _cell(cell),
      _iPopV(iPopV) { }

    const typename CELL::value_t* getComponentPointer(unsigned iDim) const any_platform
    {
      return &_cell[descriptors::q<typename CELL::descriptor_t>(_iPopV, iDim)];
    }

    typename CELL::value_t* getComponentPointer(unsigned iDim) any_platform
    {
      return &_cell[descriptors::q<typename CELL::descriptor_t>(_iPopV, iDim)];
    }

    template <typename U, typename IMPL>
    PopulationPtr& operator=(const GenericVector<U,CELL::descriptor_t::n,IMPL>& rhs) any_platform
    {
      for (unsigned iDim=0; iDim < CELL::descriptor_t::n; ++iDim) {
        this->operator[](iDim) = rhs[iDim];
      }
      return *this;
    }

  };

public:
  using value_t = typename CELL::value_t;
  using descriptor_t = typename CELL::descriptor_t;

  VectorCell(CELL cell) any_platform:
    _cell{cell} { }

  VectorCell neighbor(Vector<int,descriptor_t::d> c_i) any_platform {
    return _cell.neighbor(c_i);
  }

  PopulationPtr operator[](unsigned iPopV) any_platform {
    return PopulationPtr{_cell, iPopV};
  }

  template <typename FIELD>
  auto getField() const any_platform {
    return _cell.template getField<FIELD>();
  }

  template <typename FIELD>
  void setField(const FieldD<value_t,descriptor_t,FIELD>& v) any_platform {
    _cell.template setField<FIELD>(v);
  }

  template <typename FIELD>
  auto getFieldPointer() any_platform {
    return _cell.template getFieldPointer<FIELD>();
  }

  template <typename FIELD>
  auto getFieldPointer() const any_platform {
    return _cell.template getFieldPointer<FIELD>();
  }

  template <typename FIELD>
  auto getFieldComponent(unsigned iD) const any_platform {
    return _cell.template getFieldComponent<FIELD>(iD);
  }

  template <typename FIELD, typename VALUE>
  void setFieldComponent(unsigned iD, VALUE value) any_platform {
    _cell.template setFieldComponent<FIELD>(iD, value);
  }

};

}

#endif
