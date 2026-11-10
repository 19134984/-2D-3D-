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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<1, -1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<1, -1, -1>, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x143 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x144 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x145 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x146 = parameters.template get<descriptors::OMEGA>();
auto x147 = ((x143)*(x143));
auto x148 = ((x144)*(x144));
auto x149 = ((x145)*(x145));
auto x150 = cell[4] + cell[5];
auto x151 = cell[1] + x150;
auto x153 = cell[0] + cell[11] + cell[2] + V{1};
auto x154 = cell[12] + cell[17] + cell[18] + cell[8] + cell[9];
auto x156 = cell[13] + cell[14];
auto x157 = cell[10] + x156;
auto x0 = V{3};
auto x158 = ((x0)*(x0));
auto x1 = x147 + x148 + x149;
auto x2 = V{0.5}*x0*x1;
auto x159 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = cell[3] + cell[6];
auto x6 = -V{1}*x145;
auto x7 = x6 + V{1};
auto x155 = V{1} / (x7);
auto x8 = cell[6] + cell[7] + x151;
auto x9 = -V{1}*x143;
auto x10 = x9 + V{1};
auto x152 = V{1} / (x10);
auto x11 = V{0.5}*x152*(cell[3] + x153 + x154 + V{2}*x8) + V{0.5}*x155*(V{2}*cell[16] + V{2}*cell[18] + V{2}*cell[8] + x151 + x153 + x157 + V{2}*x5);
auto x12 = V{1}/x0;
auto x13 = -x12;
auto x14 = x143 + x145;
auto x15 = V{0.5}*x158;
auto x16 = x0*x14 + x15*((x14)*(x14)) + x159;
auto x17 = V{0.0277777777777778};
auto x160 = -x17;
auto x18 = x11*x17;
auto x19 = x16*x18 + x160;
auto x20 = x6 + x9;
auto x21 = x0*x20 + x15*((x20)*(x20)) + x159;
auto x22 = -x17 + x18*x21;
auto x23 = cell[6] + x19 - x22;
auto x24 = x11*x143;
auto x25 = x12*(x11 + V{-1});
auto x161 = -x25;
auto x162 = cell[16] + x161 + x23;
auto x26 = -x143*x24 + x157 + x162 + x8;
auto x27 = x13*x26;
auto x28 = x11*x144;
auto x29 = cell[11] + cell[17] + cell[18] + cell[2] + cell[8] + cell[9] - x144*x28 + x150 + x156 + x161;
auto x30 = x13*x29;
auto x31 = x27 + x30;
auto x32 = cell[7] - x11*x149 + x154 + x162 + x5;
auto x33 = x13*x32;
auto x34 = x4*(x11*x3 + x158*(V{0.5}*x31 + V{0.5}*x33) + V{-1});
auto x35 = V{1} - x146;
auto x36 = -V{1}*x144;
auto x37 = x143 + x36;
auto x38 = x0*x37 + x15*((x37)*(x37)) + x159;
auto x39 = V{0.5}*x158*x17;
auto x40 = V{1} - x12;
auto x41 = x29*x40;
auto x42 = x26*x40;
auto x163 = x33 + x42;
auto x164 = x163 + x41;
auto x43 = V{-2};
auto x44 = cell[13] - V{1}*cell[14] + cell[4] - V{1}*cell[5] - x144*x24;
auto x45 = x39*(x164 + x43*x44);
auto x165 = x160 + x45;
auto x46 = x165 + x18*x38;
auto x47 = x143 + x144;
auto x48 = x0*x47 + x15*((x47)*(x47)) + x159;
auto x49 = V{2};
auto x50 = x39*(x164 + x44*x49);
auto x166 = x160 + x50;
auto x51 = x166 + x18*x48;
auto x52 = V{0.0555555555555556};
auto x167 = -x52;
auto x53 = x11*x52;
auto x54 = x0*x144 + x148*x15 + x159;
auto x55 = V{0.5}*x158*x52;
auto x56 = x27 + x41;
auto x57 = x55*(x33 + x56);
auto x168 = x167 + x57;
auto x58 = x168 + x53*x54;
auto x59 = x0*x143 + x147*x15 + x159;
auto x60 = x55*(x163 + x30);
auto x169 = x167 + x60;
auto x61 = x169 + x53*x59;
auto x62 = x144 + x9;
auto x63 = x0*x62 + x15*((x62)*(x62)) + x159;
auto x64 = x165 + x18*x63;
auto x65 = x36 + x9;
auto x66 = x0*x65 + x15*((x65)*(x65)) + x159;
auto x67 = x166 + x18*x66;
auto x68 = x0*x36 + x15*((x36)*(x36)) + x159;
auto x69 = x168 + x53*x68;
auto x70 = x0*x9 + x15*((x9)*(x9)) + x159;
auto x71 = x169 + x53*x70;
auto x72 = x144 + x6;
auto x73 = x0*x72 + x15*((x72)*(x72)) + x159;
auto x74 = x32*x40;
auto x170 = x56 + x74;
auto x172 = x30 + x42 + x74;
auto x75 = cell[17] - V{1}*cell[18] + cell[8] - V{1}*cell[9] - x145*x28;
auto x76 = x39*(x170 + x43*x75);
auto x171 = x160 + x76;
auto x77 = x171 + x18*x73;
auto x78 = x143 + x6;
auto x79 = x0*x78 + x15*((x78)*(x78)) + x159;
auto x80 = -V{1}*cell[16] + cell[6] - V{1}*cell[7] - x145*x24 + x23;
auto x81 = x39*(x172 + x43*x80);
auto x173 = x160 + x81;
auto x82 = x173 + x18*x79;
auto x83 = x36 + x6;
auto x84 = x0*x83 + x15*((x83)*(x83)) + x159;
auto x85 = x39*(x170 + x49*x75);
auto x174 = x160 + x85;
auto x86 = x174 + x18*x84;
auto x87 = x0*x6 + x15*((x6)*(x6)) + x159;
auto x88 = x55*(x31 + x74);
auto x175 = x167 + x88;
auto x89 = x175 + x53*x87;
auto x90 = x39*(x172 + x49*x80);
auto x91 = x22 + x90;
auto x176 = x34 + x58 + x69 + V{2}*x91 + V{1};
auto x92 = x144 + x145;
auto x93 = x0*x92 + x15*((x92)*(x92)) + x159;
auto x94 = x174 + x18*x93;
auto x95 = x0*x145 + x149*x15 + x159;
auto x96 = x175 + x53*x95;
auto x97 = x145 + x36;
auto x98 = x0*x97 + x15*((x97)*(x97)) + x159;
auto x99 = x171 + x18*x98;
auto x100 = x145 + x9;
auto x101 = x0*x100 + x15*((x100)*(x100)) + x159;
auto x102 = x101*x18 + x173;
auto x103 = V{0.5}*x152*(V{2}*x102 + x176 + V{2}*x64 + V{2}*x67 + V{2}*x71 + x77 + x86 + x89 + x94 + x96 + x99) + V{0.5}*x155*(x176 + x46 + x51 + x61 + x64 + x67 + x71 + V{2}*x77 + V{2}*x82 + V{2}*x86 + V{2}*x89);
auto x104 = x146*x4*(x103*x3 + V{-1}) + x34*x35;
auto x105 = x103*x52;
auto x106 = x146*(x105*x70 + x167) + x35*x71;
auto x107 = x146*(x105*x68 + x167) + x35*x69;
auto x108 = x146*(x105*x87 + x167) + x35*x89;
auto x109 = x103*x17;
auto x110 = x146*(x109*x66 + x160) + x35*x67;
auto x111 = x146*(x109*x63 + x160) + x35*x64;
auto x112 = x146*(x109*x21 + x160) + x35*x91;
auto x113 = x102*x35 + x146*(x101*x109 + x160);
auto x114 = x146*(x109*x84 + x160) + x35*x86;
auto x115 = x146*(x109*x98 + x160) + x35*x99;
auto x116 = x146*(x105*x59 + x167) + x35*x61;
auto x117 = x146*(x105*x54 + x167) + x35*x58;
auto x177 = x104 + x107 + V{2}*x112 + x117 + V{1};
auto x118 = x146*(x105*x95 + x167) + x35*x96;
auto x119 = x146*(x109*x48 + x160) + x35*x51;
auto x120 = x146*(x109*x38 + x160) + x35*x46;
auto x121 = x146*(x109*x79 + x160) + x35*x82;
auto x122 = x146*(x109*x93 + x160) + x35*x94;
auto x123 = x146*(x109*x73 + x160) + x35*x77;
cell[0] = x104;
cell[10] = x116;
cell[11] = x117;
cell[12] = x118;
cell[13] = x119;
cell[14] = x120;
cell[15] = x146*(x109*x16 + x160) + x35*(x19 + x90);
cell[16] = x121;
cell[17] = x122;
cell[18] = x123;
cell[1] = x106;
cell[2] = x107;
cell[3] = x108;
cell[4] = x110;
cell[5] = x111;
cell[6] = x112;
cell[7] = x113;
cell[8] = x114;
cell[9] = x115;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = V{0.5}*x152*(V{2}*x106 + x108 + V{2}*x110 + V{2}*x111 + V{2}*x113 + x114 + x115 + x118 + x122 + x123 + x177) + V{0.5}*x155*(x106 + V{2}*x108 + x110 + x111 + V{2}*x114 + x116 + x119 + x120 + V{2}*x121 + V{2}*x123 + x177);
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x143;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x144;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x145;
return { x103, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
