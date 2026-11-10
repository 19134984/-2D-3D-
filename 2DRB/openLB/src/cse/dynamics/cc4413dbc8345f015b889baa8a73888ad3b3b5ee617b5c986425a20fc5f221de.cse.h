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
struct CSE<dynamics::Tuple<T, descriptors::D3Q27<FIELDS...>, momenta::Tuple<momenta::OneDensity, momenta::ZeroMomentum, momenta::ZeroStress, momenta::DefineSeparately>, equilibria::None, collision::None, dynamics::DefaultCombination>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x0 = cell[0];
auto x1 = cell[10];
auto x2 = cell[11];
auto x3 = cell[12];
auto x4 = cell[13];
auto x5 = cell[14];
auto x6 = cell[15];
auto x7 = cell[16];
auto x8 = cell[17];
auto x9 = cell[18];
auto x10 = cell[19];
auto x11 = cell[1];
auto x12 = cell[20];
auto x13 = cell[21];
auto x14 = cell[22];
auto x15 = cell[23];
auto x16 = cell[24];
auto x17 = cell[25];
auto x18 = cell[26];
auto x19 = cell[2];
auto x20 = cell[3];
auto x21 = cell[4];
auto x22 = cell[5];
auto x23 = cell[6];
auto x24 = cell[7];
auto x25 = cell[8];
auto x26 = cell[9];
cell[0] = x0;
cell[10] = x1;
cell[11] = x2;
cell[12] = x3;
cell[13] = x4;
cell[14] = x5;
cell[15] = x6;
cell[16] = x7;
cell[17] = x8;
cell[18] = x9;
cell[19] = x10;
cell[1] = x11;
cell[20] = x12;
cell[21] = x13;
cell[22] = x14;
cell[23] = x15;
cell[24] = x16;
cell[25] = x17;
cell[26] = x18;
cell[2] = x19;
cell[3] = x20;
cell[4] = x21;
cell[5] = x22;
cell[6] = x23;
cell[7] = x24;
cell[8] = x25;
cell[9] = x26;
return { V{-1}, V{-1} };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
