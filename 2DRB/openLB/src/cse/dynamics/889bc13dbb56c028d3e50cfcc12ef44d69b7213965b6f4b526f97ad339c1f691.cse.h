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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::ConstRhoBGK, dynamics::DefaultCombination>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x72 = parameters.template get<descriptors::OMEGA>();
auto x73 = parameters.template get<statistics::AVERAGE_RHO>();
auto x74 = cell[10] + cell[14] + cell[16];
auto x75 = cell[11] + cell[18] + cell[5];
auto x76 = cell[12] + cell[7] + cell[9];
auto x0 = V{3};
auto x1 = -V{1}*cell[4];
auto x78 = cell[13] + x1;
auto x2 = -V{1}*cell[6];
auto x79 = cell[15] + x2;
auto x3 = cell[0] + cell[13] + cell[15] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + x74 + x75 + x76 + V{1};
auto x77 = V{1} / (x3);
auto x4 = x77*(-V{1}*cell[1] - V{1}*cell[5] - V{1}*cell[7] + x74 + x78 + x79);
auto x81 = ((x4)*(x4));
auto x5 = -V{1}*cell[8];
auto x80 = cell[17] + x5;
auto x6 = x77*(-V{1}*cell[14] - V{1}*cell[2] - V{1}*cell[9] + x75 + x78 + x80);
auto x82 = ((x6)*(x6));
auto x7 = x77*(-V{1}*cell[16] - V{1}*cell[18] - V{1}*cell[3] + x76 + x79 + x80);
auto x83 = ((x7)*(x7));
auto x8 = x81 + x82 + x83;
auto x9 = V{0.5}*x0*x8;
auto x85 = V{1} - x9;
auto x10 = V{0.333333333333333};
auto x11 = -x10*(x3*(x9 + V{-1}) + V{1});
auto x12 = V{1} - x72;
auto x13 = V{1} - x73;
auto x14 = x13*x77 + V{1};
auto x15 = -V{1}*x4;
auto x16 = V{0.5}*((x0)*(x0));
auto x17 = V{0.0555555555555556};
auto x84 = -x17;
auto x18 = x17*x3;
auto x19 = x18*(x0*x15 + x16*((x15)*(x15)) + x85) + x84;
auto x20 = -V{1}*x6;
auto x21 = x18*(x0*x20 + x16*((x20)*(x20)) + x85) + x84;
auto x22 = -V{1}*x7;
auto x23 = x18*(x0*x22 + x16*((x22)*(x22)) + x85) + x84;
auto x24 = x15 + x20;
auto x25 = V{0.0277777777777778};
auto x86 = -x25;
auto x26 = x25*x3;
auto x27 = x26*(x0*x24 + x16*((x24)*(x24)) + x85) + x86;
auto x28 = x15 + x6;
auto x29 = x26*(x0*x28 + x16*((x28)*(x28)) + x85) + x86;
auto x30 = x15 + x22;
auto x31 = x26*(x0*x30 + x16*((x30)*(x30)) + x85) + x86;
auto x32 = x15 + x7;
auto x33 = x26*(x0*x32 + x16*((x32)*(x32)) + x85) + x86;
auto x34 = x20 + x22;
auto x35 = x26*(x0*x34 + x16*((x34)*(x34)) + x85) + x86;
auto x36 = x20 + x7;
auto x37 = x26*(x0*x36 + x16*((x36)*(x36)) + x85) + x86;
auto x38 = x18*(x0*x4 + x16*x81 + x85) + x84;
auto x39 = x18*(x0*x6 + x16*x82 + x85) + x84;
auto x40 = x18*(x0*x7 + x16*x83 + x85) + x84;
auto x41 = x4 + x6;
auto x42 = x26*(x0*x41 + x16*((x41)*(x41)) + x85) + x86;
auto x43 = x20 + x4;
auto x44 = x26*(x0*x43 + x16*((x43)*(x43)) + x85) + x86;
auto x45 = x4 + x7;
auto x46 = x26*(x0*x45 + x16*((x45)*(x45)) + x85) + x86;
auto x47 = x22 + x4;
auto x48 = x26*(x0*x47 + x16*((x47)*(x47)) + x85) + x86;
auto x49 = x6 + x7;
auto x50 = x26*(x0*x49 + x16*((x49)*(x49)) + x85) + x86;
auto x51 = x22 + x6;
auto x52 = x26*(x0*x51 + x16*((x51)*(x51)) + x85) + x86;
auto x53 = -x10 + x12*(cell[0] - x11) + x14*(x10 + x11);
auto x54 = x12*(cell[10] - x38) + x14*(x17 + x38) + x84;
auto x55 = x12*(cell[11] - x39) + x14*(x17 + x39) + x84;
auto x56 = x12*(cell[12] - x40) + x14*(x17 + x40) + x84;
auto x57 = x12*(cell[13] - x42) + x14*(x25 + x42) + x86;
auto x58 = x12*(cell[14] - x44) + x14*(x25 + x44) + x86;
auto x59 = x12*(cell[15] - x46) + x14*(x25 + x46) + x86;
auto x60 = x12*(cell[16] - x48) + x14*(x25 + x48) + x86;
auto x61 = x12*(cell[17] - x50) + x14*(x25 + x50) + x86;
auto x62 = x12*(cell[18] - x52) + x14*(x25 + x52) + x86;
auto x63 = x12*(cell[1] - x19) + x14*(x17 + x19) + x84;
auto x64 = x12*(cell[2] - x21) + x14*(x17 + x21) + x84;
auto x65 = x12*(cell[3] - x23) + x14*(x17 + x23) + x84;
auto x66 = x12*(cell[4] - x27) + x14*(x25 + x27) + x86;
auto x67 = x12*(cell[5] - x29) + x14*(x25 + x29) + x86;
auto x68 = x12*(cell[6] - x31) + x14*(x25 + x31) + x86;
auto x69 = x12*(cell[7] - x33) + x14*(x25 + x33) + x86;
auto x70 = x12*(cell[8] - x35) + x14*(x25 + x35) + x86;
auto x71 = x12*(cell[9] - x37) + x14*(x25 + x37) + x86;
cell[0] = x53;
cell[10] = x54;
cell[11] = x55;
cell[12] = x56;
cell[13] = x57;
cell[14] = x58;
cell[15] = x59;
cell[16] = x60;
cell[17] = x61;
cell[18] = x62;
cell[1] = x63;
cell[2] = x64;
cell[3] = x65;
cell[4] = x66;
cell[5] = x67;
cell[6] = x68;
cell[7] = x69;
cell[8] = x70;
cell[9] = x71;
return { x13 + x3, x8 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
