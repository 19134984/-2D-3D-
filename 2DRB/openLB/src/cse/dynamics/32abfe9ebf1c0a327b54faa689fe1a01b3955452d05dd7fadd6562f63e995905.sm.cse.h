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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::None, collision::None, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x23 = cell[10] + cell[14] + cell[16];
auto x24 = cell[11] + cell[18] + cell[5];
auto x25 = cell[12] + cell[7] + cell[9];
auto x0 = cell[0] + cell[13] + cell[15] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + x23 + x24 + x25 + V{1};
auto x26 = V{1} / (x0);
auto x1 = -V{1}*cell[4];
auto x27 = cell[13] + x1;
auto x2 = -V{1}*cell[6];
auto x28 = cell[15] + x2;
auto x3 = -V{1}*cell[8];
auto x29 = cell[17] + x3;
auto x4 = cell[0];
auto x5 = cell[10];
auto x6 = cell[11];
auto x7 = cell[12];
auto x8 = cell[13];
auto x9 = cell[14];
auto x10 = cell[15];
auto x11 = cell[16];
auto x12 = cell[17];
auto x13 = cell[18];
auto x14 = cell[1];
auto x15 = cell[2];
auto x16 = cell[3];
auto x17 = cell[4];
auto x18 = cell[5];
auto x19 = cell[6];
auto x20 = cell[7];
auto x21 = cell[8];
auto x22 = cell[9];
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x0;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x26*(-V{1}*cell[1] - V{1}*cell[5] - V{1}*cell[7] + x23 + x27 + x28);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x26*(-V{1}*cell[14] - V{1}*cell[2] - V{1}*cell[9] + x24 + x27 + x29);
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x26*(-V{1}*cell[16] - V{1}*cell[18] - V{1}*cell[3] + x25 + x28 + x29);
cell[0] = x4;
cell[10] = x5;
cell[11] = x6;
cell[12] = x7;
cell[13] = x8;
cell[14] = x9;
cell[15] = x10;
cell[16] = x11;
cell[17] = x12;
cell[18] = x13;
cell[1] = x14;
cell[2] = x15;
cell[3] = x16;
cell[4] = x17;
cell[5] = x18;
cell[6] = x19;
cell[7] = x20;
cell[8] = x21;
cell[9] = x22;
return { V{-1}, V{-1} };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
