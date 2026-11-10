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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerCornerDensity3D<1, 1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerCornerStress3D<1, 1, -1>, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x125 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x126 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x127 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x128 = parameters.template get<descriptors::OMEGA>();
auto x129 = ((x125)*(x125));
auto x130 = ((x126)*(x126));
auto x131 = ((x127)*(x127));
auto x132 = cell[3] + cell[6];
auto x134 = cell[13] + cell[14] + cell[4] + cell[5];
auto x135 = cell[10] + x134;
auto x137 = cell[11] + cell[2];
auto x138 = V{2}*cell[16] + x137;
auto x140 = V{2}*cell[13];
auto x142 = cell[17] + cell[18] + cell[8] + cell[9];
auto x143 = cell[12] + x142;
auto x0 = V{3};
auto x146 = ((x0)*(x0));
auto x1 = x129 + x130 + x131;
auto x2 = V{0.5}*x0*x1;
auto x147 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = cell[0] + cell[1];
auto x136 = V{2}*cell[18] + x5 + V{1};
auto x6 = x132;
auto x7 = -V{1}*x127;
auto x8 = x7 + V{1};
auto x133 = V{1} / (x8);
auto x9 = x125 + V{1};
auto x141 = V{1} / (x9);
auto x10 = x126 + V{1};
auto x139 = V{1} / (x10);
auto x11 = V{0.333333333333333}*x133*(V{2}*cell[8] + x135 + x136 + x138 + V{2}*x6) + V{0.333333333333333}*x139*(cell[10] + V{2}*cell[11] + cell[12] + cell[15] + cell[16] + V{2}*cell[17] + V{2}*cell[5] + cell[7] + x132 + x136 + x140) + V{0.333333333333333}*x141*(cell[0] + V{2}*cell[10] + V{2}*cell[14] + V{2}*cell[15] + cell[3] + x138 + x140 + x143 + V{1});
auto x12 = V{1}/x0;
auto x13 = -x12;
auto x14 = x11*x125;
auto x15 = x12*(x11 + V{-1});
auto x144 = -x15;
auto x145 = cell[15] + cell[16] + cell[7] + x144;
auto x16 = cell[1] + cell[6] - x125*x14 + x135 + x145;
auto x17 = x13*x16;
auto x18 = x11*x126;
auto x19 = -x126*x18 + x134 + x137 + x142 + x144;
auto x20 = x13*x19;
auto x21 = x17 + x20;
auto x22 = -x11*x131 + x143 + x145 + x6;
auto x23 = x13*x22;
auto x24 = x4*(x11*x3 + x146*(V{0.5}*x21 + V{0.5}*x23) + V{-1});
auto x25 = V{1} - x128;
auto x26 = V{0.0277777777777778};
auto x150 = -x26;
auto x27 = x11*x26;
auto x28 = -V{1}*x126;
auto x29 = x125 + x28;
auto x30 = V{0.5}*x146;
auto x31 = x0*x29 + x147 + x30*((x29)*(x29));
auto x32 = V{0.5}*x146*x26;
auto x33 = V{1} - x12;
auto x34 = x19*x33;
auto x35 = x16*x33;
auto x148 = x23 + x35;
auto x149 = x148 + x34;
auto x36 = V{-2};
auto x37 = cell[13] - V{1}*cell[14] + cell[4] - V{1}*cell[5] - x126*x14;
auto x38 = x32*(x149 + x36*x37);
auto x151 = x150 + x38;
auto x39 = x151 + x27*x31;
auto x40 = x125 + x126;
auto x41 = x0*x40 + x147 + x30*((x40)*(x40));
auto x42 = V{2};
auto x43 = x32*(x149 + x37*x42);
auto x152 = x150 + x43;
auto x44 = x152 + x27*x41;
auto x45 = V{0.0555555555555556};
auto x153 = -x45;
auto x46 = x11*x45;
auto x47 = x0*x126 + x130*x30 + x147;
auto x48 = V{0.5}*x146*x45;
auto x49 = x17 + x34;
auto x50 = x48*(x23 + x49);
auto x154 = x153 + x50;
auto x51 = x154 + x46*x47;
auto x52 = x0*x125 + x129*x30 + x147;
auto x53 = x48*(x148 + x20);
auto x155 = x153 + x53;
auto x54 = x155 + x46*x52;
auto x55 = -V{1}*x125;
auto x56 = x126 + x55;
auto x57 = x0*x56 + x147 + x30*((x56)*(x56));
auto x58 = x151 + x27*x57;
auto x59 = x28 + x55;
auto x60 = x0*x59 + x147 + x30*((x59)*(x59));
auto x61 = x152 + x27*x60;
auto x62 = x0*x28 + x147 + x30*((x28)*(x28));
auto x63 = x154 + x46*x62;
auto x64 = x0*x55 + x147 + x30*((x55)*(x55));
auto x65 = x155 + x46*x64;
auto x66 = x24 + x65;
auto x67 = x126 + x7;
auto x68 = x0*x67 + x147 + x30*((x67)*(x67));
auto x69 = x22*x33;
auto x156 = x49 + x69;
auto x158 = x20 + x35 + x69;
auto x70 = cell[17] - V{1}*cell[18] + cell[8] - V{1}*cell[9] - x127*x18;
auto x71 = x32*(x156 + x36*x70);
auto x157 = x150 + x71;
auto x72 = x157 + x27*x68;
auto x163 = x54 + x66 + V{2}*x72 + V{1};
auto x73 = x125 + x7;
auto x74 = x0*x73 + x147 + x30*((x73)*(x73));
auto x75 = cell[15] - V{1}*cell[16] + cell[6] - V{1}*cell[7] - x127*x14;
auto x76 = x32*(x158 + x36*x75);
auto x159 = x150 + x76;
auto x77 = x159 + x27*x74;
auto x165 = x51 + x63 + V{2}*x77;
auto x78 = x28 + x7;
auto x79 = x0*x78 + x147 + x30*((x78)*(x78));
auto x80 = x32*(x156 + x42*x70);
auto x160 = x150 + x80;
auto x81 = x160 + x27*x79;
auto x82 = x0*x7 + x147 + x30*((x7)*(x7));
auto x83 = x48*(x21 + x69);
auto x161 = x153 + x83;
auto x84 = x161 + x46*x82;
auto x85 = x55 + x7;
auto x86 = x0*x85 + x147 + x30*((x85)*(x85));
auto x87 = x32*(x158 + x42*x75);
auto x162 = x150 + x87;
auto x88 = x162 + x27*x86;
auto x89 = x126 + x127;
auto x90 = x0*x89 + x147 + x30*((x89)*(x89));
auto x91 = x160 + x27*x90;
auto x92 = x0*x127 + x131*x30 + x147;
auto x93 = x161 + x46*x92;
auto x164 = V{2}*x44 + x84 + x93;
auto x94 = x127 + x28;
auto x95 = x0*x94 + x147 + x30*((x94)*(x94));
auto x96 = x157 + x27*x95;
auto x97 = x125 + x127;
auto x98 = x0*x97 + x147 + x30*((x97)*(x97));
auto x99 = x162 + x27*x98;
auto x100 = x127 + x55;
auto x101 = x0*x100 + x147 + x30*((x100)*(x100));
auto x102 = x101*x27 + x159;
auto x103 = V{0.333333333333333}*x133*(x163 + x165 + x39 + x44 + x58 + x61 + V{2}*x81 + V{2}*x84 + V{2}*x88) + V{0.333333333333333}*x139*(x102 + x163 + x164 + V{2}*x51 + V{2}*x58 + x77 + x88 + V{2}*x91 + x99) + V{0.333333333333333}*x141*(x164 + x165 + x24 + V{2}*x39 + V{2}*x54 + x72 + x81 + x91 + x96 + V{2}*x99 + V{1});
auto x104 = x103*x45;
auto x105 = x103*x26;
cell[0] = x128*x4*(x103*x3 + V{-1}) + x24*x25;
cell[10] = x128*(x104*x52 + x153) + x25*x54;
cell[11] = x128*(x104*x47 + x153) + x25*x51;
cell[12] = x128*(x104*x92 + x153) + x25*x93;
cell[13] = x128*(x105*x41 + x150) + x25*x44;
cell[14] = x128*(x105*x31 + x150) + x25*x39;
cell[15] = x128*(x105*x98 + x150) + x25*x99;
cell[16] = x128*(x105*x74 + x150) + x25*x77;
cell[17] = x128*(x105*x90 + x150) + x25*x91;
cell[18] = x128*(x105*x68 + x150) + x25*x72;
cell[1] = x128*(x104*x64 + x153) + x25*x65;
cell[2] = x128*(x104*x62 + x153) + x25*x63;
cell[3] = x128*(x104*x82 + x153) + x25*x84;
cell[4] = x128*(x105*x60 + x150) + x25*x61;
cell[5] = x128*(x105*x57 + x150) + x25*x58;
cell[6] = x128*(x105*x86 + x150) + x25*x88;
cell[7] = x102*x25 + x128*(x101*x105 + x150);
cell[8] = x128*(x105*x79 + x150) + x25*x81;
cell[9] = x128*(x105*x95 + x150) + x25*x96;
return { x103, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
