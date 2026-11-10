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
struct CSE<SourcedAdvectionDiffusionBGKdynamics<T, descriptors::D3Q7<FIELDS...> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x24 = cell.template getFieldComponent<olb::descriptors::SOURCE>(0);
auto x25 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x26 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x27 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x28 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{1} - x28;
auto x1 = V{0.25};
auto x2 = V{0.5}*x24;
auto x29 = x2 + V{1};
auto x3 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + x29;
auto x4 = -x24*(V{0.5}*x28 + V{-1});
auto x5 = cell[0]*x0 + x1*x28*(x3 + V{-1}) + x1*x4;
auto x6 = V{4};
auto x30 = x25*x6;
auto x32 = x26*x6;
auto x33 = x27*x6;
auto x7 = V{0.125};
auto x8 = x3*x7;
auto x31 = V{1}*x8;
auto x9 = x4*x7;
auto x10 = cell[1]*x0 - x28*(x31*(x30 + V{-1}) + x7) + x9;
auto x11 = cell[2]*x0 - x28*(x31*(x32 + V{-1}) + x7) + x9;
auto x12 = cell[3]*x0 - x28*(x31*(x33 + V{-1}) + x7) + x9;
auto x13 = cell[4]*x0 - x28*(x7 - x8*(x30 + V{1})) + x9;
auto x14 = cell[5]*x0 - x28*(x7 - x8*(x32 + V{1})) + x9;
auto x15 = cell[6]*x0 - x28*(x7 - x8*(x33 + V{1})) + x9;
auto x16 = x10 + x11 + x12 + x13 + x14 + x15 + x29 + x5;
auto x34 = V{1} / (x16);
cell[0] = x5;
cell[1] = x10;
cell[2] = x11;
cell[3] = x12;
cell[4] = x13;
cell[5] = x14;
cell[6] = x15;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x16;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = -x34*(V{1}*x10 - x13);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = -x34*(V{1}*x11 - x14);
cell.template getFieldPointer<olb::fields::moments::U>()[2] = -x34*(V{1}*x12 - x15);
return { x3, ((x25)*(x25)) + ((x26)*(x26)) + ((x27)*(x27)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
