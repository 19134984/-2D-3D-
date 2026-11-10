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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<1, -1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<1, -1, -1>, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x125 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x126 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x127 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x128 = parameters.template get<descriptors::OMEGA>();
auto x129 = ((x125)*(x125));
auto x130 = ((x126)*(x126));
auto x131 = ((x127)*(x127));
auto x132 = cell[4] + cell[5];
auto x133 = cell[1] + x132;
auto x135 = cell[0] + cell[11] + cell[2] + V{1};
auto x136 = cell[12] + cell[17] + cell[18] + cell[8] + cell[9];
auto x138 = cell[13] + cell[14];
auto x139 = cell[10] + x138;
auto x0 = V{3};
auto x140 = ((x0)*(x0));
auto x1 = x129 + x130 + x131;
auto x2 = V{0.5}*x0*x1;
auto x141 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = cell[3] + cell[6];
auto x6 = -V{1}*x127;
auto x7 = x6 + V{1};
auto x137 = V{1} / (x7);
auto x8 = cell[6] + cell[7] + x133;
auto x9 = -V{1}*x125;
auto x10 = x9 + V{1};
auto x134 = V{1} / (x10);
auto x11 = V{0.5}*x134*(cell[3] + x135 + x136 + V{2}*x8) + V{0.5}*x137*(V{2}*cell[16] + V{2}*cell[18] + V{2}*cell[8] + x133 + x135 + x139 + V{2}*x5);
auto x12 = V{1}/x0;
auto x13 = -x12;
auto x14 = x125 + x127;
auto x15 = V{0.5}*x140;
auto x16 = x0*x14 + x141 + x15*((x14)*(x14));
auto x17 = V{0.0277777777777778};
auto x142 = -x17;
auto x18 = x11*x17;
auto x19 = x142 + x16*x18;
auto x20 = x6 + x9;
auto x21 = x0*x20 + x141 + x15*((x20)*(x20));
auto x22 = -x17 + x18*x21;
auto x23 = cell[6] + x19 - x22;
auto x24 = x11*x125;
auto x25 = x12*(x11 + V{-1});
auto x143 = -x25;
auto x144 = cell[16] + x143 + x23;
auto x26 = -x125*x24 + x139 + x144 + x8;
auto x27 = x13*x26;
auto x28 = x11*x126;
auto x29 = cell[11] + cell[17] + cell[18] + cell[2] + cell[8] + cell[9] - x126*x28 + x132 + x138 + x143;
auto x30 = x13*x29;
auto x31 = x27 + x30;
auto x32 = cell[7] - x11*x131 + x136 + x144 + x5;
auto x33 = x13*x32;
auto x34 = x4*(x11*x3 + x140*(V{0.5}*x31 + V{0.5}*x33) + V{-1});
auto x35 = V{1} - x128;
auto x36 = -V{1}*x126;
auto x37 = x125 + x36;
auto x38 = x0*x37 + x141 + x15*((x37)*(x37));
auto x39 = V{0.5}*x140*x17;
auto x40 = V{1} - x12;
auto x41 = x29*x40;
auto x42 = x26*x40;
auto x145 = x33 + x42;
auto x146 = x145 + x41;
auto x43 = V{-2};
auto x44 = cell[13] - V{1}*cell[14] + cell[4] - V{1}*cell[5] - x126*x24;
auto x45 = x39*(x146 + x43*x44);
auto x147 = x142 + x45;
auto x46 = x147 + x18*x38;
auto x47 = x125 + x126;
auto x48 = x0*x47 + x141 + x15*((x47)*(x47));
auto x49 = V{2};
auto x50 = x39*(x146 + x44*x49);
auto x148 = x142 + x50;
auto x51 = x148 + x18*x48;
auto x52 = V{0.0555555555555556};
auto x149 = -x52;
auto x53 = x11*x52;
auto x54 = x0*x126 + x130*x15 + x141;
auto x55 = V{0.5}*x140*x52;
auto x56 = x27 + x41;
auto x57 = x55*(x33 + x56);
auto x150 = x149 + x57;
auto x58 = x150 + x53*x54;
auto x59 = x0*x125 + x129*x15 + x141;
auto x60 = x55*(x145 + x30);
auto x151 = x149 + x60;
auto x61 = x151 + x53*x59;
auto x62 = x126 + x9;
auto x63 = x0*x62 + x141 + x15*((x62)*(x62));
auto x64 = x147 + x18*x63;
auto x65 = x36 + x9;
auto x66 = x0*x65 + x141 + x15*((x65)*(x65));
auto x67 = x148 + x18*x66;
auto x68 = x0*x36 + x141 + x15*((x36)*(x36));
auto x69 = x150 + x53*x68;
auto x70 = x0*x9 + x141 + x15*((x9)*(x9));
auto x71 = x151 + x53*x70;
auto x72 = x126 + x6;
auto x73 = x0*x72 + x141 + x15*((x72)*(x72));
auto x74 = x32*x40;
auto x152 = x56 + x74;
auto x154 = x30 + x42 + x74;
auto x75 = cell[17] - V{1}*cell[18] + cell[8] - V{1}*cell[9] - x127*x28;
auto x76 = x39*(x152 + x43*x75);
auto x153 = x142 + x76;
auto x77 = x153 + x18*x73;
auto x78 = x125 + x6;
auto x79 = x0*x78 + x141 + x15*((x78)*(x78));
auto x80 = -V{1}*cell[16] + cell[6] - V{1}*cell[7] - x127*x24 + x23;
auto x81 = x39*(x154 + x43*x80);
auto x155 = x142 + x81;
auto x82 = x155 + x18*x79;
auto x83 = x36 + x6;
auto x84 = x0*x83 + x141 + x15*((x83)*(x83));
auto x85 = x39*(x152 + x49*x75);
auto x156 = x142 + x85;
auto x86 = x156 + x18*x84;
auto x87 = x0*x6 + x141 + x15*((x6)*(x6));
auto x88 = x55*(x31 + x74);
auto x157 = x149 + x88;
auto x89 = x157 + x53*x87;
auto x90 = x39*(x154 + x49*x80);
auto x91 = x22 + x90;
auto x158 = x34 + x58 + x69 + V{2}*x91 + V{1};
auto x92 = x126 + x127;
auto x93 = x0*x92 + x141 + x15*((x92)*(x92));
auto x94 = x156 + x18*x93;
auto x95 = x0*x127 + x131*x15 + x141;
auto x96 = x157 + x53*x95;
auto x97 = x127 + x36;
auto x98 = x0*x97 + x141 + x15*((x97)*(x97));
auto x99 = x153 + x18*x98;
auto x100 = x127 + x9;
auto x101 = x0*x100 + x141 + x15*((x100)*(x100));
auto x102 = x101*x18 + x155;
auto x103 = V{0.5}*x134*(V{2}*x102 + x158 + V{2}*x64 + V{2}*x67 + V{2}*x71 + x77 + x86 + x89 + x94 + x96 + x99) + V{0.5}*x137*(x158 + x46 + x51 + x61 + x64 + x67 + x71 + V{2}*x77 + V{2}*x82 + V{2}*x86 + V{2}*x89);
auto x104 = x103*x52;
auto x105 = x103*x17;
cell[0] = x128*x4*(x103*x3 + V{-1}) + x34*x35;
cell[10] = x128*(x104*x59 + x149) + x35*x61;
cell[11] = x128*(x104*x54 + x149) + x35*x58;
cell[12] = x128*(x104*x95 + x149) + x35*x96;
cell[13] = x128*(x105*x48 + x142) + x35*x51;
cell[14] = x128*(x105*x38 + x142) + x35*x46;
cell[15] = x128*(x105*x16 + x142) + x35*(x19 + x90);
cell[16] = x128*(x105*x79 + x142) + x35*x82;
cell[17] = x128*(x105*x93 + x142) + x35*x94;
cell[18] = x128*(x105*x73 + x142) + x35*x77;
cell[1] = x128*(x104*x70 + x149) + x35*x71;
cell[2] = x128*(x104*x68 + x149) + x35*x69;
cell[3] = x128*(x104*x87 + x149) + x35*x89;
cell[4] = x128*(x105*x66 + x142) + x35*x67;
cell[5] = x128*(x105*x63 + x142) + x35*x64;
cell[6] = x128*(x105*x21 + x142) + x35*x91;
cell[7] = x102*x35 + x128*(x101*x105 + x142);
cell[8] = x128*(x105*x84 + x142) + x35*x86;
cell[9] = x128*(x105*x98 + x142) + x35*x99;
return { x103, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
