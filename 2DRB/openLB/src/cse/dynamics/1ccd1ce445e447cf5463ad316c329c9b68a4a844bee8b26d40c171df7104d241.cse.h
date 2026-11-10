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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<0, 1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<0, 1, 1>, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x125 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x126 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x127 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x128 = parameters.template get<descriptors::OMEGA>();
auto x129 = ((x125)*(x125));
auto x130 = ((x126)*(x126));
auto x131 = ((x127)*(x127));
auto x133 = cell[15] + cell[16] + cell[6] + cell[7];
auto x134 = cell[10] + x133;
auto x135 = cell[12] + cell[3];
auto x138 = cell[11] + cell[2];
auto x0 = V{3};
auto x141 = ((x0)*(x0));
auto x1 = x129 + x130 + x131;
auto x2 = V{0.5}*x0*x1;
auto x142 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = cell[0] + cell[1];
auto x136 = V{2}*cell[17] + x5 + V{1};
auto x6 = x126 + V{1};
auto x132 = V{1} / (x6);
auto x7 = x127 + V{1};
auto x137 = V{1} / (x7);
auto x8 = V{0.5}*x132*(V{2}*cell[11] + V{2}*cell[13] + V{2}*cell[18] + V{2}*cell[5] + x134 + x135 + x136) + V{0.5}*x137*(cell[10] + V{2}*cell[12] + cell[13] + cell[14] + V{2}*cell[15] + cell[4] + cell[5] + V{2}*cell[7] + V{2}*cell[9] + x136 + x138);
auto x9 = V{1}/x0;
auto x10 = -x9;
auto x11 = x125*x8;
auto x12 = x9*(x8 + V{-1});
auto x139 = -x12;
auto x140 = cell[13] + cell[14] + cell[4] + cell[5] + x139;
auto x13 = cell[1] - x11*x125 + x134 + x140;
auto x14 = x10*x13;
auto x15 = -V{1}*x126;
auto x16 = -V{1}*x127;
auto x17 = x15 + x16;
auto x18 = V{0.5}*x141;
auto x19 = x0*x17 + x142 + x18*((x17)*(x17));
auto x20 = V{0.0277777777777778};
auto x143 = -x20;
auto x21 = x20*x8;
auto x22 = x143 + x19*x21;
auto x23 = x126 + x127;
auto x24 = x0*x23 + x142 + x18*((x23)*(x23));
auto x25 = -x20 + x21*x24;
auto x26 = cell[17] + x22 - x25;
auto x144 = cell[17] + cell[18] + cell[9] + x26;
auto x27 = x126*x8;
auto x28 = -x126*x27 + x138 + x140 + x144;
auto x29 = x10*x28;
auto x30 = x14 + x29;
auto x31 = -x131*x8 + x133 + x135 + x139 + x144;
auto x32 = x10*x31;
auto x33 = x4*(x141*(V{0.5}*x30 + V{0.5}*x32) + x3*x8 + V{-1});
auto x34 = V{1} - x128;
auto x35 = x125 + x16;
auto x36 = x0*x35 + x142 + x18*((x35)*(x35));
auto x37 = V{0.5}*x141*x20;
auto x38 = V{1} - x9;
auto x39 = x31*x38;
auto x40 = x13*x38;
auto x145 = x29 + x40;
auto x146 = x145 + x39;
auto x41 = V{-2};
auto x42 = cell[15] - V{1}*cell[16] + cell[6] - V{1}*cell[7] - x11*x127;
auto x43 = x37*(x146 + x41*x42);
auto x147 = x143 + x43;
auto x44 = x147 + x21*x36;
auto x45 = x125 + x127;
auto x46 = x0*x45 + x142 + x18*((x45)*(x45));
auto x47 = V{2};
auto x48 = x37*(x146 + x42*x47);
auto x148 = x143 + x48;
auto x49 = x148 + x21*x46;
auto x50 = V{0.0555555555555556};
auto x149 = -x50;
auto x51 = x50*x8;
auto x52 = x0*x127 + x131*x18 + x142;
auto x53 = V{0.5}*x141*x50;
auto x54 = x53*(x30 + x39);
auto x150 = x149 + x54;
auto x55 = x150 + x51*x52;
auto x56 = x0*x125 + x129*x18 + x142;
auto x57 = x53*(x145 + x32);
auto x151 = x149 + x57;
auto x58 = x151 + x51*x56;
auto x59 = -V{1}*x125;
auto x60 = x127 + x59;
auto x61 = x0*x60 + x142 + x18*((x60)*(x60));
auto x62 = x147 + x21*x61;
auto x63 = x16 + x59;
auto x64 = x0*x63 + x142 + x18*((x63)*(x63));
auto x65 = x148 + x21*x64;
auto x66 = x0*x16 + x142 + x18*((x16)*(x16));
auto x67 = x150 + x51*x66;
auto x68 = x0*x59 + x142 + x18*((x59)*(x59));
auto x69 = x151 + x51*x68;
auto x70 = x33 + x69;
auto x71 = x126 + x16;
auto x72 = x0*x71 + x142 + x18*((x71)*(x71));
auto x73 = x28*x38;
auto x154 = x32 + x40 + x73;
auto x74 = x14 + x73;
auto x152 = x39 + x74;
auto x75 = cell[17] - V{1}*cell[18] - V{1}*cell[9] - x127*x27 + x26;
auto x76 = x37*(x152 + x41*x75);
auto x153 = x143 + x76;
auto x77 = x153 + x21*x72;
auto x78 = x37*(x152 + x47*x75);
auto x79 = x25 + x78;
auto x158 = x58 + x70 + V{2}*x79 + V{1};
auto x80 = x125 + x126;
auto x81 = x0*x80 + x142 + x18*((x80)*(x80));
auto x82 = cell[13] - V{1}*cell[14] + cell[4] - V{1}*cell[5] - x11*x126;
auto x83 = x37*(x154 + x47*x82);
auto x155 = x143 + x83;
auto x84 = x155 + x21*x81;
auto x85 = x126 + x59;
auto x86 = x0*x85 + x142 + x18*((x85)*(x85));
auto x87 = x37*(x154 + x41*x82);
auto x156 = x143 + x87;
auto x88 = x156 + x21*x86;
auto x89 = x0*x126 + x130*x18 + x142;
auto x90 = x53*(x32 + x74);
auto x157 = x149 + x90;
auto x91 = x157 + x51*x89;
auto x92 = x125 + x15;
auto x93 = x0*x92 + x142 + x18*((x92)*(x92));
auto x94 = x156 + x21*x93;
auto x95 = x15 + x59;
auto x96 = x0*x95 + x142 + x18*((x95)*(x95));
auto x97 = x155 + x21*x96;
auto x98 = x0*x15 + x142 + x18*((x15)*(x15));
auto x99 = x157 + x51*x98;
auto x100 = x127 + x15;
auto x101 = x0*x100 + x142 + x18*((x100)*(x100));
auto x102 = x101*x21 + x153;
auto x103 = V{0.5}*x132*(x158 + x44 + x49 + x55 + x62 + x65 + x67 + V{2}*x77 + V{2}*x84 + V{2}*x88 + V{2}*x91) + V{0.5}*x137*(V{2}*x102 + x158 + V{2}*x49 + V{2}*x55 + V{2}*x62 + x84 + x88 + x91 + x94 + x97 + x99);
auto x104 = x103*x50;
auto x105 = x103*x20;
cell[0] = x128*x4*(x103*x3 + V{-1}) + x33*x34;
cell[10] = x128*(x104*x56 + x149) + x34*x58;
cell[11] = x128*(x104*x89 + x149) + x34*x91;
cell[12] = x128*(x104*x52 + x149) + x34*x55;
cell[13] = x128*(x105*x81 + x143) + x34*x84;
cell[14] = x128*(x105*x93 + x143) + x34*x94;
cell[15] = x128*(x105*x46 + x143) + x34*x49;
cell[16] = x128*(x105*x36 + x143) + x34*x44;
cell[17] = x128*(x105*x24 + x143) + x34*x79;
cell[18] = x128*(x105*x72 + x143) + x34*x77;
cell[1] = x128*(x104*x68 + x149) + x34*x69;
cell[2] = x128*(x104*x98 + x149) + x34*x99;
cell[3] = x128*(x104*x66 + x149) + x34*x67;
cell[4] = x128*(x105*x96 + x143) + x34*x97;
cell[5] = x128*(x105*x86 + x143) + x34*x88;
cell[6] = x128*(x105*x64 + x143) + x34*x65;
cell[7] = x128*(x105*x61 + x143) + x34*x62;
cell[8] = x128*(x105*x19 + x143) + x34*(x22 + x78);
cell[9] = x102*x34 + x128*(x101*x105 + x143);
return { x103, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
