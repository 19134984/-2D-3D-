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

#ifndef STATISTICS_FIELD_REDUCER_H
#define STATISTICS_FIELD_REDUCER_H

#include <typeinfo>

#include "latticeStatistics.h"
#include "fieldStatisticO.h"
#include "superLattice.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
template<typename FIELD>
void LatticeStatistics<T,DESCRIPTOR>::registerField(std::string name)
{
  if (name.empty()) {
    if constexpr (requires { FIELD::getName(); }) {
      name = FIELD::getName();
    } else {
      name = typeid(FIELD).name();
    }
  }
  RegisteredField rf{
    typeid(FIELD),
    name,
    DESCRIPTOR::template size<FIELD>(),
    [](SuperLattice<T,DESCRIPTOR>& sLattice) {
      sLattice.template addPostProcessor<stage::Statistics>(
        PostProcessorPromise<T,DESCRIPTOR>{meta::id<FieldStatisticO<FIELD>>{}});
    }
  };
  _registered.push_back(std::move(rf));
}

}  // namespace olb

#endif
