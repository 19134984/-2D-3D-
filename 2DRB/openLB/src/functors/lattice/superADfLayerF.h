/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Johannes Grafen
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

#ifndef SUPER_ADF_LAYER_F_H
#define SUPER_ADF_LAYER_F_H

#include <vector>

#include "superBaseF2D.h"
#include "superBaseF3D.h"
#include "utilities/aDiff.h"
#include "utilities/aliases.h"

namespace olb {

/// Super functor exposing one layer of an ADf-valued super functor as plain floating point values
/**
 * Evaluates the wrapped functor and returns either the primal value (iD < 0) or the
 * derivative with respect to the iD-th differentiation variable, e.g. to write the
 * results of an ADf-typed case with SuperVTMwriter<T,D,float,BaseType<T>>.
 **/
template <unsigned D, typename T> requires util::is_adf_v<T>
class SuperADfLayerF final : public SuperF<D,T,BaseType<T>> {
private:
  using base_t = SuperF<D,T,BaseType<T>>;

  SuperF<D,T,T>& _f;
  const int      _iD;
public:
  /**
   * \param f   ADf-valued super functor
   * \param iD  Index of the derivative to return, the primal value is returned for iD < 0
   **/
  SuperADfLayerF(SuperF<D,T,T>& f, int iD = -1)
    : base_t(f.getSuperStructure(), f.getTargetDim()),
      _f(f), _iD(iD)
  {
    OLB_ASSERT(_iD < static_cast<int>(T::dim), "derivative index exceeds ADf dimension");
    this->getName() = f.getName();
  }

  bool operator()(BaseType<T> output[], const int input[]) override
  {
    std::vector<T> evaluated(_f.getTargetDim());
    const bool result = _f(evaluated.data(), input);
    for (int iDim = 0; iDim < _f.getTargetDim(); ++iDim) {
      output[iDim] = (_iD < 0) ? evaluated[iDim].v() : evaluated[iDim].d(_iD);
    }
    return result;
  }
};

}

#endif
