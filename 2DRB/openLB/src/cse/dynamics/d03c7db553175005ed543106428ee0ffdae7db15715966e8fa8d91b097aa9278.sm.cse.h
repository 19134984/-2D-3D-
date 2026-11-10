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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerCornerDensity3D<-1, 1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerCornerStress3D<-1, 1, -1>, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x146 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x147 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x148 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x149 = parameters.template get<descriptors::OMEGA>();
auto x150 = ((x146)*(x146));
auto x151 = ((x147)*(x147));
auto x152 = ((x148)*(x148));
auto x153 = cell[3] + cell[6];
auto x154 = cell[4] + cell[5];
auto x156 = cell[12] + V{1};
auto x157 = cell[11] + cell[2];
auto x158 = cell[17] + cell[18] + cell[8] + cell[9];
auto x159 = x157 + x158;
auto x161 = cell[13] + cell[14];
auto x162 = cell[10] + x161;
auto x0 = V{3};
auto x167 = ((x0)*(x0));
auto x1 = x150 + x151 + x152;
auto x2 = V{0.5}*x0*x1;
auto x168 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = cell[0] + cell[1];
auto x163 = V{2}*cell[18] + x5;
auto x6 = x153;
auto x7 = -V{1}*x148;
auto x8 = x7 + V{1};
auto x160 = V{1} / (x8);
auto x9 = cell[1] + cell[6] + cell[7] + x154;
auto x10 = -V{1}*x146;
auto x11 = x10 + V{1};
auto x155 = V{1} / (x11);
auto x12 = x147 + V{1};
auto x164 = V{1} / (x12);
auto x13 = V{0.333333333333333}*x155*(cell[0] + cell[3] + x156 + x159 + V{2}*x9) + V{0.333333333333333}*x160*(V{2}*cell[16] + V{2}*cell[8] + x154 + x157 + x162 + x163 + V{2}*x6 + V{1}) + V{0.333333333333333}*x164*(cell[10] + V{2}*cell[11] + V{2}*cell[13] + cell[15] + cell[16] + V{2}*cell[17] + V{2}*cell[5] + cell[7] + x153 + x156 + x163);
auto x14 = V{1}/x0;
auto x15 = -x14;
auto x16 = x13*x146;
auto x17 = x14*(x13 + V{-1});
auto x165 = -x17;
auto x166 = cell[15] + cell[16] + x165;
auto x18 = -x146*x16 + x162 + x166 + x9;
auto x19 = x15*x18;
auto x20 = x13*x147;
auto x21 = -x147*x20 + x154 + x159 + x161 + x165;
auto x22 = x15*x21;
auto x23 = x19 + x22;
auto x24 = cell[12] + cell[7] - x13*x152 + x158 + x166 + x6;
auto x25 = x15*x24;
auto x26 = x4*(x13*x3 + x167*(V{0.5}*x23 + V{0.5}*x25) + V{-1});
auto x27 = V{1} - x149;
auto x28 = V{0.0277777777777778};
auto x171 = -x28;
auto x29 = x13*x28;
auto x30 = -V{1}*x147;
auto x31 = x146 + x30;
auto x32 = V{0.5}*x167;
auto x33 = x0*x31 + x168 + x32*((x31)*(x31));
auto x34 = V{0.5}*x167*x28;
auto x35 = V{1} - x14;
auto x36 = x21*x35;
auto x37 = x18*x35;
auto x169 = x25 + x37;
auto x170 = x169 + x36;
auto x38 = V{-2};
auto x39 = cell[13] - V{1}*cell[14] + cell[4] - V{1}*cell[5] - x147*x16;
auto x40 = x34*(x170 + x38*x39);
auto x172 = x171 + x40;
auto x41 = x172 + x29*x33;
auto x42 = x146 + x147;
auto x43 = x0*x42 + x168 + x32*((x42)*(x42));
auto x44 = V{2};
auto x45 = x34*(x170 + x39*x44);
auto x173 = x171 + x45;
auto x46 = x173 + x29*x43;
auto x47 = V{0.0555555555555556};
auto x174 = -x47;
auto x48 = x13*x47;
auto x49 = x0*x147 + x151*x32 + x168;
auto x50 = V{0.5}*x167*x47;
auto x51 = x19 + x36;
auto x52 = x50*(x25 + x51);
auto x175 = x174 + x52;
auto x53 = x175 + x48*x49;
auto x54 = x0*x146 + x150*x32 + x168;
auto x55 = x50*(x169 + x22);
auto x176 = x174 + x55;
auto x56 = x176 + x48*x54;
auto x57 = x10 + x147;
auto x58 = x0*x57 + x168 + x32*((x57)*(x57));
auto x59 = x172 + x29*x58;
auto x60 = x10 + x30;
auto x61 = x0*x60 + x168 + x32*((x60)*(x60));
auto x62 = x173 + x29*x61;
auto x63 = x0*x30 + x168 + x32*((x30)*(x30));
auto x64 = x175 + x48*x63;
auto x65 = x0*x10 + x168 + x32*((x10)*(x10));
auto x66 = x176 + x48*x65;
auto x67 = x26 + x66;
auto x68 = x147 + x7;
auto x69 = x0*x68 + x168 + x32*((x68)*(x68));
auto x70 = x24*x35;
auto x177 = x51 + x70;
auto x179 = x22 + x37 + x70;
auto x71 = cell[17] - V{1}*cell[18] + cell[8] - V{1}*cell[9] - x148*x20;
auto x72 = x34*(x177 + x38*x71);
auto x178 = x171 + x72;
auto x73 = x178 + x29*x69;
auto x184 = x56 + x67 + V{2}*x73 + V{1};
auto x74 = x146 + x7;
auto x75 = x0*x74 + x168 + x32*((x74)*(x74));
auto x76 = cell[15] - V{1}*cell[16] + cell[6] - V{1}*cell[7] - x148*x16;
auto x77 = x34*(x179 + x38*x76);
auto x180 = x171 + x77;
auto x78 = x180 + x29*x75;
auto x79 = x30 + x7;
auto x80 = x0*x79 + x168 + x32*((x79)*(x79));
auto x81 = x34*(x177 + x44*x71);
auto x181 = x171 + x81;
auto x82 = x181 + x29*x80;
auto x83 = x0*x7 + x168 + x32*((x7)*(x7));
auto x84 = x50*(x23 + x70);
auto x182 = x174 + x84;
auto x85 = x182 + x48*x83;
auto x86 = x10 + x7;
auto x87 = x0*x86 + x168 + x32*((x86)*(x86));
auto x88 = x34*(x179 + x44*x76);
auto x183 = x171 + x88;
auto x89 = x183 + x29*x87;
auto x186 = x53 + x64 + V{2}*x89;
auto x90 = x147 + x148;
auto x91 = x0*x90 + x168 + x32*((x90)*(x90));
auto x92 = x181 + x29*x91;
auto x93 = x0*x148 + x152*x32 + x168;
auto x94 = x182 + x48*x93;
auto x185 = V{2}*x59 + x85 + x94;
auto x95 = x148 + x30;
auto x96 = x0*x95 + x168 + x32*((x95)*(x95));
auto x97 = x178 + x29*x96;
auto x98 = x10 + x148;
auto x99 = x0*x98 + x168 + x32*((x98)*(x98));
auto x100 = x180 + x29*x99;
auto x101 = x146 + x148;
auto x102 = x0*x101 + x168 + x32*((x101)*(x101));
auto x103 = x102*x29 + x183;
auto x104 = V{0.333333333333333}*x155*(V{2}*x100 + x185 + x186 + x26 + V{2}*x62 + V{2}*x66 + x73 + x82 + x92 + x97 + V{1}) + V{0.333333333333333}*x160*(x184 + x186 + x41 + x46 + x59 + x62 + V{2}*x78 + V{2}*x82 + V{2}*x85) + V{0.333333333333333}*x164*(x100 + x103 + x184 + x185 + V{2}*x46 + V{2}*x53 + x78 + x89 + V{2}*x92);
auto x105 = x149*x4*(x104*x3 + V{-1}) + x26*x27;
auto x106 = x104*x47;
auto x107 = x149*(x106*x65 + x174) + x27*x66;
auto x108 = x149*(x106*x63 + x174) + x27*x64;
auto x109 = x149*(x106*x83 + x174) + x27*x85;
auto x110 = x104*x28;
auto x111 = x149*(x110*x61 + x171) + x27*x62;
auto x112 = x149*(x110*x58 + x171) + x27*x59;
auto x113 = x149*(x110*x87 + x171) + x27*x89;
auto x114 = x100*x27 + x149*(x110*x99 + x171);
auto x115 = x149*(x110*x80 + x171) + x27*x82;
auto x116 = x149*(x110*x96 + x171) + x27*x97;
auto x117 = x149*(x106*x54 + x174) + x27*x56;
auto x118 = x149*(x106*x49 + x174) + x27*x53;
auto x189 = x108 + V{2}*x113 + x118;
auto x119 = x149*(x106*x93 + x174) + x27*x94;
auto x188 = x109 + V{2}*x112 + x119;
auto x120 = x149*(x110*x43 + x171) + x27*x46;
auto x121 = x149*(x110*x33 + x171) + x27*x41;
auto x122 = x103*x27 + x149*(x102*x110 + x171);
auto x123 = x149*(x110*x75 + x171) + x27*x78;
auto x124 = x149*(x110*x91 + x171) + x27*x92;
auto x125 = x149*(x110*x69 + x171) + x27*x73;
auto x126 = x105 + x107;
auto x187 = x117 + V{2}*x125 + x126 + V{1};
cell[0] = x105;
cell[10] = x117;
cell[11] = x118;
cell[12] = x119;
cell[13] = x120;
cell[14] = x121;
cell[15] = x122;
cell[16] = x123;
cell[17] = x124;
cell[18] = x125;
cell[1] = x107;
cell[2] = x108;
cell[3] = x109;
cell[4] = x111;
cell[5] = x112;
cell[6] = x113;
cell[7] = x114;
cell[8] = x115;
cell[9] = x116;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = V{0.333333333333333}*x155*(x105 + V{2}*x107 + V{2}*x111 + V{2}*x114 + x115 + x116 + x124 + x125 + x188 + x189 + V{1}) + V{0.333333333333333}*x160*(V{2}*x109 + x111 + x112 + V{2}*x115 + x120 + x121 + V{2}*x123 + x187 + x189) + V{0.333333333333333}*x164*(x113 + x114 + V{2}*x118 + V{2}*x120 + x122 + x123 + V{2}*x124 + x187 + x188);
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x146;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x147;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x148;
return { x104, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
