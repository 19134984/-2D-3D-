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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::CBCDensity, momenta::CBCMomentum<2, 1>, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x59 = parameters.template get<descriptors::OMEGA>();
auto x51 = cell[10] + cell[14] + cell[16];
auto x52 = cell[11] + cell[18] + cell[5];
auto x53 = cell[12] + cell[7] + cell[9];
auto x0 = V{1} - x59;
auto x1 = V{3};
auto x2 = -V{1}*cell[4];
auto x55 = cell[13] + x2;
auto x3 = -V{1}*cell[6];
auto x56 = cell[15] + x3;
auto x4 = cell[0] + cell[13] + cell[15] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + x51 + x52 + x53 + V{1};
auto x54 = V{1} / (x4);
auto x5 = x54*(-V{1}*cell[1] - V{1}*cell[5] - V{1}*cell[7] + x51 + x55 + x56);
auto x58 = ((x5)*(x5));
auto x6 = -V{1}*cell[8];
auto x57 = cell[17] + x6;
auto x7 = x54*(-V{1}*cell[14] - V{1}*cell[2] - V{1}*cell[9] + x52 + x55 + x57);
auto x60 = ((x7)*(x7));
auto x8 = x54*(-V{1}*cell[16] - V{1}*cell[18] - V{1}*cell[3] + x53 + x56 + x57);
auto x61 = ((x8)*(x8));
auto x9 = x58 + x60 + x61;
auto x10 = V{0.5}*x1*x9;
auto x62 = V{1} - x10;
auto x11 = V{0.333333333333333};
auto x12 = -V{1}*x5;
auto x13 = V{0.5}*((x1)*(x1));
auto x14 = V{0.0555555555555556};
auto x15 = x14*x4;
auto x16 = -V{1}*x7;
auto x17 = -V{1}*x8;
auto x18 = x12 + x16;
auto x19 = V{0.0277777777777778};
auto x20 = x19*x4;
auto x21 = x12 + x7;
auto x22 = x12 + x17;
auto x23 = x12 + x8;
auto x24 = x16 + x17;
auto x25 = x16 + x8;
auto x26 = x5 + x7;
auto x27 = x16 + x5;
auto x28 = x5 + x8;
auto x29 = x17 + x5;
auto x30 = x7 + x8;
auto x31 = x17 + x7;
auto x32 = cell[0]*x0 - x11*x59*(x4*(x10 + V{-1}) + V{1});
auto x33 = cell[10]*x0 - x59*(x14 - x15*(x1*x5 + x13*x58 + x62));
auto x34 = cell[11]*x0 - x59*(x14 - x15*(x1*x7 + x13*x60 + x62));
auto x35 = cell[12]*x0 - x59*(x14 - x15*(x1*x8 + x13*x61 + x62));
auto x36 = cell[13]*x0 - x59*(x19 - x20*(x1*x26 + x13*((x26)*(x26)) + x62));
auto x37 = cell[14]*x0 - x59*(x19 - x20*(x1*x27 + x13*((x27)*(x27)) + x62));
auto x38 = cell[15]*x0 - x59*(x19 - x20*(x1*x28 + x13*((x28)*(x28)) + x62));
auto x39 = cell[16]*x0 - x59*(x19 - x20*(x1*x29 + x13*((x29)*(x29)) + x62));
auto x40 = cell[17]*x0 - x59*(x19 - x20*(x1*x30 + x13*((x30)*(x30)) + x62));
auto x41 = cell[18]*x0 - x59*(x19 - x20*(x1*x31 + x13*((x31)*(x31)) + x62));
auto x42 = cell[1]*x0 - x59*(x14 - x15*(x1*x12 + x13*((x12)*(x12)) + x62));
auto x43 = cell[2]*x0 - x59*(x14 - x15*(x1*x16 + x13*((x16)*(x16)) + x62));
auto x44 = cell[3]*x0 - x59*(x14 - x15*(x1*x17 + x13*((x17)*(x17)) + x62));
auto x45 = cell[4]*x0 - x59*(x19 - x20*(x1*x18 + x13*((x18)*(x18)) + x62));
auto x46 = cell[5]*x0 - x59*(x19 - x20*(x1*x21 + x13*((x21)*(x21)) + x62));
auto x47 = cell[6]*x0 - x59*(x19 - x20*(x1*x22 + x13*((x22)*(x22)) + x62));
auto x48 = cell[7]*x0 - x59*(x19 - x20*(x1*x23 + x13*((x23)*(x23)) + x62));
auto x49 = cell[8]*x0 - x59*(x19 - x20*(x1*x24 + x13*((x24)*(x24)) + x62));
auto x50 = cell[9]*x0 - x59*(x19 - x20*(x1*x25 + x13*((x25)*(x25)) + x62));
cell[0] = x32;
cell[10] = x33;
cell[11] = x34;
cell[12] = x35;
cell[13] = x36;
cell[14] = x37;
cell[15] = x38;
cell[16] = x39;
cell[17] = x40;
cell[18] = x41;
cell[1] = x42;
cell[2] = x43;
cell[3] = x44;
cell[4] = x45;
cell[5] = x46;
cell[6] = x47;
cell[7] = x48;
cell[8] = x49;
cell[9] = x50;
return { x4, x9 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
