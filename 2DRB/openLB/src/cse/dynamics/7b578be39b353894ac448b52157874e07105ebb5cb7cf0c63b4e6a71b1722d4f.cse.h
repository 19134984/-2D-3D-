/*  This file is part of the OpenLB library
 *
 *  Copyright (C) Adrian Kummerlaender, Shota Ito
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

/*  ========================================================
 *  ==  WARNING: This is an automatically generated file, ==
 *  ==                  do not modify.                    ==
 *  ========================================================
 */

#pragma once


namespace olb {

namespace dynamics {

template <typename T, typename... FIELDS>
struct CSE<SourcedAdvectionDiffusionBGKdynamics<T, descriptors::D3Q7<FIELDS...> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x15 = cell.template getFieldComponent<olb::descriptors::SOURCE>(0);
auto x16 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x17 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x18 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{1} - x19;
auto x1 = V{0.25};
auto x2 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + V{0.5}*x15 + V{1};
auto x3 = -x15*(V{0.5}*x19 + V{-1});
auto x4 = V{4};
auto x20 = x16*x4;
auto x22 = x17*x4;
auto x23 = x18*x4;
auto x5 = V{0.125};
auto x6 = x2*x5;
auto x21 = V{1}*x6;
auto x7 = x3*x5;
auto x8 = cell[0]*x0 + x1*x19*(x2 + V{-1}) + x1*x3;
auto x9 = cell[1]*x0 - x19*(x21*(x20 + V{-1}) + x5) + x7;
auto x10 = cell[2]*x0 - x19*(x21*(x22 + V{-1}) + x5) + x7;
auto x11 = cell[3]*x0 - x19*(x21*(x23 + V{-1}) + x5) + x7;
auto x12 = cell[4]*x0 - x19*(x5 - x6*(x20 + V{1})) + x7;
auto x13 = cell[5]*x0 - x19*(x5 - x6*(x22 + V{1})) + x7;
auto x14 = cell[6]*x0 - x19*(x5 - x6*(x23 + V{1})) + x7;
cell[0] = x8;
cell[1] = x9;
cell[2] = x10;
cell[3] = x11;
cell[4] = x12;
cell[5] = x13;
cell[6] = x14;
return { x2, ((x16)*(x16)) + ((x17)*(x17)) + ((x18)*(x18)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
