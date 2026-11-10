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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SaveVelocity<collision::BGK>, dynamics::DefaultCombination>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<0, 1>, momenta::RegularizedBoundaryStress<0, 1>, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x151 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x153 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1);
auto x154 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<0, 1>::VELOCITY>(2);
auto x155 = parameters.template get<descriptors::OMEGA>();
auto x148 = V{1} / (x151);
auto x149 = V{2}*cell[10];
auto x150 = cell[11] + cell[2];
auto x152 = cell[12] + cell[3];
auto x156 = ((x153)*(x153));
auto x157 = ((x154)*(x154));
auto x0 = V{3};
auto x160 = ((x0)*(x0));
auto x1 = V{0.5}*x0;
auto x2 = x148*(cell[0] + V{2}*cell[13] + V{2}*cell[14] + V{2}*cell[15] + V{2}*cell[16] + cell[17] + cell[18] + cell[8] + cell[9] + x149 + x150 + x152 + V{1}) + V{-1};
auto x158 = ((x2)*(x2));
auto x3 = x156;
auto x4 = x157;
auto x159 = x3 + x4;
auto x5 = x1*(x158 + x159);
auto x162 = -x5;
auto x163 = x162 + V{1};
auto x6 = V{0.333333333333333};
auto x7 = x151*x6;
auto x8 = V{1}/x0;
auto x9 = -x8;
auto x10 = V{0.0555555555555556};
auto x161 = -x10;
auto x11 = x10*x151;
auto x12 = V{0.5}*x160;
auto x13 = x11*(x0*x2 + x12*x158 + x163) + x161;
auto x14 = V{0.0277777777777778};
auto x164 = -x14;
auto x15 = x14*x151;
auto x16 = x153 + x2;
auto x17 = x15*(x0*x16 + x12*((x16)*(x16)) + x163) + x164;
auto x18 = V{2}*cell[13] - V{2}*x17;
auto x19 = -V{1}*x153;
auto x20 = x19 + x2;
auto x21 = x15*(x0*x20 + x12*((x20)*(x20)) + x163) + x164;
auto x22 = cell[14] - x21;
auto x23 = x154 + x2;
auto x24 = x15*(x0*x23 + x12*((x23)*(x23)) + x163) + x164;
auto x25 = V{2}*cell[15] - V{2}*x24;
auto x26 = -V{1}*x154;
auto x27 = x2 + x26;
auto x28 = x15*(x0*x27 + x12*((x27)*(x27)) + x163) + x164;
auto x29 = cell[16] - x28;
auto x30 = -V{2}*x13 + x149 + x18 + V{2}*x22 + x25 + V{2}*x29;
auto x31 = x30*x9;
auto x32 = x153 + x26;
auto x33 = x0*x32 + x12*((x32)*(x32)) + V{1};
auto x34 = -x14 + x15*(x162 + x33);
auto x35 = cell[18] - x34;
auto x36 = x153 + x154;
auto x37 = x0*x36 + x12*((x36)*(x36)) + V{1};
auto x38 = -x14 + x15*(x162 + x37);
auto x39 = cell[17] - x38;
auto x40 = x0*x153 + x12*x156 + V{1};
auto x41 = -x10 + x11*(x162 + x40);
auto x42 = x154 + x19;
auto x43 = x0*x42 + x12*((x42)*(x42)) + V{1};
auto x44 = -x14 + x15*(x162 + x43);
auto x45 = cell[9] - x44;
auto x46 = x0*x19 + x12*((x19)*(x19)) + V{1};
auto x47 = -x10 + x11*(x162 + x46);
auto x48 = x19 + x26;
auto x49 = x0*x48 + x12*((x48)*(x48)) + V{1};
auto x50 = -x14 + x15*(x162 + x49);
auto x51 = cell[8] - x50;
auto x165 = x35 + x39 + x45 + x51;
auto x52 = V{-2};
auto x53 = -V{1}*x52;
auto x54 = x150 + x165 + x18 + x22*x53 - x41 - x47;
auto x55 = x54*x9;
auto x56 = x31 + x55;
auto x57 = x0*x154 + x12*x157 + V{1};
auto x166 = -x57;
auto x58 = -x10 - x11*(x166 + x5);
auto x59 = x0*x26 + x12*((x26)*(x26)) + V{1};
auto x167 = -x59;
auto x60 = -x10 - x11*(x167 + x5);
auto x61 = x152 + x165 + x25 + x29*x53 - x58 - x60;
auto x62 = x61*x9;
auto x63 = V{0.5}*x160*x6*(x56 + x62) - x6 - x7*(x5 + V{-1});
auto x64 = V{1} - x155;
auto x65 = V{0.5}*x14*x160;
auto x66 = V{1} - x8;
auto x67 = x61*x66;
auto x68 = x54*x66;
auto x69 = x31 + x68;
auto x168 = x67 + x69;
auto x70 = -V{1}*x35 + x39 - V{1}*x45 + x51;
auto x71 = x65*(x168 + x52*x70);
auto x72 = x34 + x71;
auto x73 = V{2};
auto x74 = x65*(x168 + x70*x73);
auto x75 = x38 + x74;
auto x76 = V{0.5}*x10*x160;
auto x77 = x76*(x56 + x67);
auto x78 = x58 + x77;
auto x79 = x76*(x62 + x69);
auto x80 = x41 + x79;
auto x81 = x44 + x71;
auto x82 = x50 + x74;
auto x83 = x60 + x77;
auto x84 = x47 + x79;
auto x85 = x30*x66;
auto x169 = x55 + x85;
auto x170 = x169 + x67;
auto x171 = x62 + x68 + x85;
auto x86 = x25 + x29*x52;
auto x87 = x65*(x170 + x52*x86);
auto x88 = x28 + x87;
auto x89 = x65*(x170 + x73*x86);
auto x90 = x24 + x89;
auto x91 = x18 + x22*x52;
auto x92 = x65*(x171 + x52*x91);
auto x93 = x21 + x92;
auto x94 = x76*(x169 + x62);
auto x95 = x13 + x94;
auto x96 = x65*(x171 + x73*x91);
auto x97 = x17 + x96;
auto x98 = x148*(x63 + x72 + x75 + x78 + x80 + x81 + x82 + x83 + x84 + V{2}*x88 + V{2}*x90 + V{2}*x93 + V{2}*x95 + V{2}*x97 + V{1}) + V{-1};
auto x172 = ((x98)*(x98));
auto x99 = x159 + x172;
auto x100 = x1*x99;
auto x173 = V{1} - x100;
auto x101 = -x155*(x6 - x7*(V{1} - x100)) + x63*x64;
auto x102 = -V{1}*x2;
auto x103 = -V{1}*x98;
auto x104 = -x155*(x10 - x11*(-x100 + x46)) + x64*x84;
auto x105 = -x155*(x10 - x11*(-x100 - x167)) + x64*x83;
auto x106 = x102 + x19;
auto x107 = x103 + x19;
auto x108 = x102 + x153;
auto x109 = x103 + x153;
auto x110 = x102 + x26;
auto x111 = x103 + x26;
auto x112 = x102 + x154;
auto x113 = x103 + x154;
auto x114 = -x155*(x14 - x15*(-x100 + x49)) + x64*x82;
auto x115 = -x155*(x14 - x15*(-x100 + x43)) + x64*x81;
auto x116 = -x155*(x10 - x11*(x0*x98 + x12*x172 + x173)) + x64*x95;
auto x117 = -x155*(x10 - x11*(-x100 + x40)) + x64*x80;
auto x118 = -x155*(x10 - x11*(-x100 - x166)) + x64*x78;
auto x119 = x153 + x98;
auto x120 = -x155*(x14 - x15*(x0*x119 + x12*((x119)*(x119)) + x173)) + x64*x97;
auto x121 = x19 + x98;
auto x122 = -x155*(x14 - x15*(x0*x121 + x12*((x121)*(x121)) + x173)) + x64*x93;
auto x123 = x154 + x98;
auto x124 = -x155*(x14 - x15*(x0*x123 + x12*((x123)*(x123)) + x173)) + x64*x90;
auto x125 = x26 + x98;
auto x126 = -x155*(x14 - x15*(x0*x125 + x12*((x125)*(x125)) + x173)) + x64*x88;
auto x127 = -x155*(x14 - x15*(-x100 + x37)) + x64*x75;
auto x128 = -x155*(x14 - x15*(-x100 + x33)) + x64*x72;
cell[0] = x101;
cell[10] = x116;
cell[11] = x117;
cell[12] = x118;
cell[13] = x120;
cell[14] = x122;
cell[15] = x124;
cell[16] = x126;
cell[17] = x127;
cell[18] = x128;
cell[1] = -x155*(x10 - x11*(x0*x103 + x12*((x103)*(x103)) + x173)) + x64*(x11*(x0*x102 + x12*((x102)*(x102)) + x163) + x161 + x94);
cell[2] = x104;
cell[3] = x105;
cell[4] = -x155*(x14 - x15*(x0*x107 + x12*((x107)*(x107)) + x173)) + x64*(x15*(x0*x106 + x12*((x106)*(x106)) + x163) + x164 + x96);
cell[5] = -x155*(x14 - x15*(x0*x109 + x12*((x109)*(x109)) + x173)) + x64*(x15*(x0*x108 + x12*((x108)*(x108)) + x163) + x164 + x92);
cell[6] = -x155*(x14 - x15*(x0*x111 + x12*((x111)*(x111)) + x173)) + x64*(x15*(x0*x110 + x12*((x110)*(x110)) + x163) + x164 + x89);
cell[7] = -x155*(x14 - x15*(x0*x113 + x12*((x113)*(x113)) + x173)) + x64*(x15*(x0*x112 + x12*((x112)*(x112)) + x163) + x164 + x87);
cell[8] = x114;
cell[9] = x115;
cell.template getFieldPointer<olb::descriptors::VELOCITY>()[0] = x148*(x101 + x104 + x105 + x114 + x115 + V{2}*x116 + x117 + x118 + V{2}*x120 + V{2}*x122 + V{2}*x124 + V{2}*x126 + x127 + x128 + V{1}) + V{-1};
cell.template getFieldPointer<olb::descriptors::VELOCITY>()[1] = x153;
cell.template getFieldPointer<olb::descriptors::VELOCITY>()[2] = x154;
return { x151, x99 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
