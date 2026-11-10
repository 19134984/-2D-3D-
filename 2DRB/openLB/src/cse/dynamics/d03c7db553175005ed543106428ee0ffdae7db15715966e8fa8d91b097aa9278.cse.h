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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerCornerDensity3D<-1, 1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerCornerStress3D<-1, 1, -1>, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x126 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x127 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x128 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x129 = parameters.template get<descriptors::OMEGA>();
auto x130 = ((x126)*(x126));
auto x131 = ((x127)*(x127));
auto x132 = ((x128)*(x128));
auto x133 = cell[3] + cell[6];
auto x134 = cell[4] + cell[5];
auto x136 = cell[12] + V{1};
auto x137 = cell[11] + cell[2];
auto x138 = cell[17] + cell[18] + cell[8] + cell[9];
auto x139 = x137 + x138;
auto x141 = cell[13] + cell[14];
auto x142 = cell[10] + x141;
auto x0 = V{3};
auto x147 = ((x0)*(x0));
auto x1 = x130 + x131 + x132;
auto x2 = V{0.5}*x0*x1;
auto x148 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = cell[0] + cell[1];
auto x143 = V{2}*cell[18] + x5;
auto x6 = x133;
auto x7 = -V{1}*x128;
auto x8 = x7 + V{1};
auto x140 = V{1} / (x8);
auto x9 = cell[1] + cell[6] + cell[7] + x134;
auto x10 = -V{1}*x126;
auto x11 = x10 + V{1};
auto x135 = V{1} / (x11);
auto x12 = x127 + V{1};
auto x144 = V{1} / (x12);
auto x13 = V{0.333333333333333}*x135*(cell[0] + cell[3] + x136 + x139 + V{2}*x9) + V{0.333333333333333}*x140*(V{2}*cell[16] + V{2}*cell[8] + x134 + x137 + x142 + x143 + V{2}*x6 + V{1}) + V{0.333333333333333}*x144*(cell[10] + V{2}*cell[11] + V{2}*cell[13] + cell[15] + cell[16] + V{2}*cell[17] + V{2}*cell[5] + cell[7] + x133 + x136 + x143);
auto x14 = V{1}/x0;
auto x15 = -x14;
auto x16 = x126*x13;
auto x17 = x14*(x13 + V{-1});
auto x145 = -x17;
auto x146 = cell[15] + cell[16] + x145;
auto x18 = -x126*x16 + x142 + x146 + x9;
auto x19 = x15*x18;
auto x20 = x127*x13;
auto x21 = -x127*x20 + x134 + x139 + x141 + x145;
auto x22 = x15*x21;
auto x23 = x19 + x22;
auto x24 = cell[12] + cell[7] - x13*x132 + x138 + x146 + x6;
auto x25 = x15*x24;
auto x26 = x4*(x13*x3 + x147*(V{0.5}*x23 + V{0.5}*x25) + V{-1});
auto x27 = V{1} - x129;
auto x28 = V{0.0277777777777778};
auto x151 = -x28;
auto x29 = x13*x28;
auto x30 = -V{1}*x127;
auto x31 = x126 + x30;
auto x32 = V{0.5}*x147;
auto x33 = x0*x31 + x148 + x32*((x31)*(x31));
auto x34 = V{0.5}*x147*x28;
auto x35 = V{1} - x14;
auto x36 = x21*x35;
auto x37 = x18*x35;
auto x149 = x25 + x37;
auto x150 = x149 + x36;
auto x38 = V{-2};
auto x39 = cell[13] - V{1}*cell[14] + cell[4] - V{1}*cell[5] - x127*x16;
auto x40 = x34*(x150 + x38*x39);
auto x152 = x151 + x40;
auto x41 = x152 + x29*x33;
auto x42 = x126 + x127;
auto x43 = x0*x42 + x148 + x32*((x42)*(x42));
auto x44 = V{2};
auto x45 = x34*(x150 + x39*x44);
auto x153 = x151 + x45;
auto x46 = x153 + x29*x43;
auto x47 = V{0.0555555555555556};
auto x154 = -x47;
auto x48 = x13*x47;
auto x49 = x0*x127 + x131*x32 + x148;
auto x50 = V{0.5}*x147*x47;
auto x51 = x19 + x36;
auto x52 = x50*(x25 + x51);
auto x155 = x154 + x52;
auto x53 = x155 + x48*x49;
auto x54 = x0*x126 + x130*x32 + x148;
auto x55 = x50*(x149 + x22);
auto x156 = x154 + x55;
auto x56 = x156 + x48*x54;
auto x57 = x10 + x127;
auto x58 = x0*x57 + x148 + x32*((x57)*(x57));
auto x59 = x152 + x29*x58;
auto x60 = x10 + x30;
auto x61 = x0*x60 + x148 + x32*((x60)*(x60));
auto x62 = x153 + x29*x61;
auto x63 = x0*x30 + x148 + x32*((x30)*(x30));
auto x64 = x155 + x48*x63;
auto x65 = x0*x10 + x148 + x32*((x10)*(x10));
auto x66 = x156 + x48*x65;
auto x67 = x26 + x66;
auto x68 = x127 + x7;
auto x69 = x0*x68 + x148 + x32*((x68)*(x68));
auto x70 = x24*x35;
auto x157 = x51 + x70;
auto x159 = x22 + x37 + x70;
auto x71 = cell[17] - V{1}*cell[18] + cell[8] - V{1}*cell[9] - x128*x20;
auto x72 = x34*(x157 + x38*x71);
auto x158 = x151 + x72;
auto x73 = x158 + x29*x69;
auto x164 = x56 + x67 + V{2}*x73 + V{1};
auto x74 = x126 + x7;
auto x75 = x0*x74 + x148 + x32*((x74)*(x74));
auto x76 = cell[15] - V{1}*cell[16] + cell[6] - V{1}*cell[7] - x128*x16;
auto x77 = x34*(x159 + x38*x76);
auto x160 = x151 + x77;
auto x78 = x160 + x29*x75;
auto x79 = x30 + x7;
auto x80 = x0*x79 + x148 + x32*((x79)*(x79));
auto x81 = x34*(x157 + x44*x71);
auto x161 = x151 + x81;
auto x82 = x161 + x29*x80;
auto x83 = x0*x7 + x148 + x32*((x7)*(x7));
auto x84 = x50*(x23 + x70);
auto x162 = x154 + x84;
auto x85 = x162 + x48*x83;
auto x86 = x10 + x7;
auto x87 = x0*x86 + x148 + x32*((x86)*(x86));
auto x88 = x34*(x159 + x44*x76);
auto x163 = x151 + x88;
auto x89 = x163 + x29*x87;
auto x166 = x53 + x64 + V{2}*x89;
auto x90 = x127 + x128;
auto x91 = x0*x90 + x148 + x32*((x90)*(x90));
auto x92 = x161 + x29*x91;
auto x93 = x0*x128 + x132*x32 + x148;
auto x94 = x162 + x48*x93;
auto x165 = V{2}*x59 + x85 + x94;
auto x95 = x128 + x30;
auto x96 = x0*x95 + x148 + x32*((x95)*(x95));
auto x97 = x158 + x29*x96;
auto x98 = x10 + x128;
auto x99 = x0*x98 + x148 + x32*((x98)*(x98));
auto x100 = x160 + x29*x99;
auto x101 = x126 + x128;
auto x102 = x0*x101 + x148 + x32*((x101)*(x101));
auto x103 = x102*x29 + x163;
auto x104 = V{0.333333333333333}*x135*(V{2}*x100 + x165 + x166 + x26 + V{2}*x62 + V{2}*x66 + x73 + x82 + x92 + x97 + V{1}) + V{0.333333333333333}*x140*(x164 + x166 + x41 + x46 + x59 + x62 + V{2}*x78 + V{2}*x82 + V{2}*x85) + V{0.333333333333333}*x144*(x100 + x103 + x164 + x165 + V{2}*x46 + V{2}*x53 + x78 + x89 + V{2}*x92);
auto x105 = x104*x47;
auto x106 = x104*x28;
cell[0] = x129*x4*(x104*x3 + V{-1}) + x26*x27;
cell[10] = x129*(x105*x54 + x154) + x27*x56;
cell[11] = x129*(x105*x49 + x154) + x27*x53;
cell[12] = x129*(x105*x93 + x154) + x27*x94;
cell[13] = x129*(x106*x43 + x151) + x27*x46;
cell[14] = x129*(x106*x33 + x151) + x27*x41;
cell[15] = x103*x27 + x129*(x102*x106 + x151);
cell[16] = x129*(x106*x75 + x151) + x27*x78;
cell[17] = x129*(x106*x91 + x151) + x27*x92;
cell[18] = x129*(x106*x69 + x151) + x27*x73;
cell[1] = x129*(x105*x65 + x154) + x27*x66;
cell[2] = x129*(x105*x63 + x154) + x27*x64;
cell[3] = x129*(x105*x83 + x154) + x27*x85;
cell[4] = x129*(x106*x61 + x151) + x27*x62;
cell[5] = x129*(x106*x58 + x151) + x27*x59;
cell[6] = x129*(x106*x87 + x151) + x27*x89;
cell[7] = x100*x27 + x129*(x106*x99 + x151);
cell[8] = x129*(x106*x80 + x151) + x27*x82;
cell[9] = x129*(x106*x96 + x151) + x27*x97;
return { x104, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
