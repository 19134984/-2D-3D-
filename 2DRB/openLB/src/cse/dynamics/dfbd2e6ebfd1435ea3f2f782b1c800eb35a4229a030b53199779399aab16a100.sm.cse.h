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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<0, -1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<0, -1, -1>, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x146 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x147 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x148 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x149 = parameters.template get<descriptors::OMEGA>();
auto x150 = ((x146)*(x146));
auto x151 = ((x147)*(x147));
auto x152 = ((x148)*(x148));
auto x153 = cell[2] + cell[4];
auto x154 = cell[3] + cell[6];
auto x158 = cell[15] + cell[16] + cell[7];
auto x159 = cell[12] + x158;
auto x0 = V{3};
auto x162 = ((x0)*(x0));
auto x1 = x150 + x151 + x152;
auto x2 = V{0.5}*x0*x1;
auto x163 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = cell[0] + cell[1];
auto x156 = cell[10] + V{2}*cell[8] + x5 + V{1};
auto x6 = x153;
auto x7 = -V{1}*x147;
auto x8 = x7 + V{1};
auto x157 = V{1} / (x8);
auto x9 = x154;
auto x10 = -V{1}*x148;
auto x11 = x10 + V{1};
auto x155 = V{1} / (x11);
auto x12 = V{0.5}*x155*(cell[11] + cell[13] + cell[14] + V{2}*cell[16] + V{2}*cell[18] + cell[5] + x153 + x156 + V{2}*x9) + V{0.5}*x157*(V{2}*cell[14] + V{2}*cell[9] + x154 + x156 + x159 + V{2}*x6);
auto x13 = V{1}/x0;
auto x14 = -x13;
auto x15 = x12*x146;
auto x16 = x13*(x12 + V{-1});
auto x160 = -x16;
auto x161 = cell[13] + cell[14] + cell[5] + x160;
auto x17 = cell[10] + cell[1] + cell[4] + cell[6] - x146*x15 + x158 + x161;
auto x18 = x14*x17;
auto x19 = x147 + x148;
auto x20 = V{0.5}*x162;
auto x21 = x0*x19 + x163 + x20*((x19)*(x19));
auto x22 = V{0.0277777777777778};
auto x164 = -x22;
auto x23 = x12*x22;
auto x24 = x164 + x21*x23;
auto x25 = x10 + x7;
auto x26 = x0*x25 + x163 + x20*((x25)*(x25));
auto x27 = -x22 + x23*x26;
auto x28 = cell[8] + x24 - x27;
auto x165 = cell[18] + cell[8] + cell[9] + x28;
auto x29 = x12*x147;
auto x30 = cell[11] - x147*x29 + x161 + x165 + x6;
auto x31 = x14*x30;
auto x32 = x18 + x31;
auto x33 = -x12*x152 + x159 + x160 + x165 + x9;
auto x34 = x14*x33;
auto x35 = x4*(x12*x3 + x162*(V{0.5}*x32 + V{0.5}*x34) + V{-1});
auto x36 = V{1} - x149;
auto x37 = x10 + x146;
auto x38 = x0*x37 + x163 + x20*((x37)*(x37));
auto x39 = V{0.5}*x162*x22;
auto x40 = V{1} - x13;
auto x41 = x33*x40;
auto x42 = x17*x40;
auto x166 = x31 + x42;
auto x167 = x166 + x41;
auto x43 = V{-2};
auto x44 = cell[15] - V{1}*cell[16] + cell[6] - V{1}*cell[7] - x148*x15;
auto x45 = x39*(x167 + x43*x44);
auto x168 = x164 + x45;
auto x46 = x168 + x23*x38;
auto x47 = x146 + x148;
auto x48 = x0*x47 + x163 + x20*((x47)*(x47));
auto x49 = V{2};
auto x50 = x39*(x167 + x44*x49);
auto x169 = x164 + x50;
auto x51 = x169 + x23*x48;
auto x52 = V{0.0555555555555556};
auto x170 = -x52;
auto x53 = x12*x52;
auto x54 = x0*x148 + x152*x20 + x163;
auto x55 = V{0.5}*x162*x52;
auto x56 = x55*(x32 + x41);
auto x171 = x170 + x56;
auto x57 = x171 + x53*x54;
auto x58 = x0*x146 + x150*x20 + x163;
auto x59 = x55*(x166 + x34);
auto x172 = x170 + x59;
auto x60 = x172 + x53*x58;
auto x61 = -V{1}*x146;
auto x62 = x148 + x61;
auto x63 = x0*x62 + x163 + x20*((x62)*(x62));
auto x64 = x168 + x23*x63;
auto x65 = x10 + x61;
auto x66 = x0*x65 + x163 + x20*((x65)*(x65));
auto x67 = x169 + x23*x66;
auto x68 = x0*x10 + x163 + x20*((x10)*(x10));
auto x69 = x171 + x53*x68;
auto x70 = x0*x61 + x163 + x20*((x61)*(x61));
auto x71 = x172 + x53*x70;
auto x72 = x35 + x71;
auto x73 = x146 + x7;
auto x74 = x0*x73 + x163 + x20*((x73)*(x73));
auto x75 = x30*x40;
auto x173 = x34 + x42 + x75;
auto x76 = cell[13] - V{1}*cell[14] + cell[4] - V{1}*cell[5] - x147*x15;
auto x77 = x39*(x173 + x43*x76);
auto x174 = x164 + x77;
auto x78 = x174 + x23*x74;
auto x79 = x148 + x7;
auto x80 = x0*x79 + x163 + x20*((x79)*(x79));
auto x81 = x18 + x75;
auto x175 = x41 + x81;
auto x82 = -V{1}*cell[18] + cell[8] - V{1}*cell[9] - x148*x29 + x28;
auto x83 = x39*(x175 + x43*x82);
auto x176 = x164 + x83;
auto x84 = x176 + x23*x80;
auto x85 = x39*(x175 + x49*x82);
auto x86 = x27 + x85;
auto x179 = x60 + x72 + V{2}*x86 + V{1};
auto x87 = x0*x7 + x163 + x20*((x7)*(x7));
auto x88 = x55*(x34 + x81);
auto x177 = x170 + x88;
auto x89 = x177 + x53*x87;
auto x90 = x61 + x7;
auto x91 = x0*x90 + x163 + x20*((x90)*(x90));
auto x92 = x39*(x173 + x49*x76);
auto x178 = x164 + x92;
auto x93 = x178 + x23*x91;
auto x94 = x146 + x147;
auto x95 = x0*x94 + x163 + x20*((x94)*(x94));
auto x96 = x178 + x23*x95;
auto x97 = x0*x147 + x151*x20 + x163;
auto x98 = x177 + x53*x97;
auto x99 = x147 + x61;
auto x100 = x0*x99 + x163 + x20*((x99)*(x99));
auto x101 = x100*x23 + x174;
auto x102 = x10 + x147;
auto x103 = x0*x102 + x163 + x20*((x102)*(x102));
auto x104 = x103*x23 + x176;
auto x105 = V{0.5}*x155*(x101 + V{2}*x104 + x179 + V{2}*x46 + V{2}*x67 + V{2}*x69 + x78 + x89 + x93 + x96 + x98) + V{0.5}*x157*(x179 + x46 + x51 + x57 + x64 + x67 + x69 + V{2}*x78 + V{2}*x84 + V{2}*x89 + V{2}*x93);
auto x106 = x149*x4*(x105*x3 + V{-1}) + x35*x36;
auto x107 = x105*x52;
auto x108 = x149*(x107*x70 + x170) + x36*x71;
auto x109 = x149*(x107*x87 + x170) + x36*x89;
auto x110 = x149*(x107*x68 + x170) + x36*x69;
auto x111 = x105*x22;
auto x112 = x149*(x111*x91 + x164) + x36*x93;
auto x113 = x101*x36 + x149*(x100*x111 + x164);
auto x114 = x149*(x111*x66 + x164) + x36*x67;
auto x115 = x149*(x111*x63 + x164) + x36*x64;
auto x116 = x149*(x111*x26 + x164) + x36*x86;
auto x117 = x149*(x111*x80 + x164) + x36*x84;
auto x118 = x149*(x107*x58 + x170) + x36*x60;
auto x119 = x149*(x107*x97 + x170) + x36*x98;
auto x120 = x149*(x107*x54 + x170) + x36*x57;
auto x121 = x149*(x111*x95 + x164) + x36*x96;
auto x122 = x149*(x111*x74 + x164) + x36*x78;
auto x123 = x149*(x111*x48 + x164) + x36*x51;
auto x124 = x149*(x111*x38 + x164) + x36*x46;
auto x125 = x104*x36 + x149*(x103*x111 + x164);
auto x126 = x106 + x108;
auto x180 = V{2}*x116 + x118 + x126 + V{1};
cell[0] = x106;
cell[10] = x118;
cell[11] = x119;
cell[12] = x120;
cell[13] = x121;
cell[14] = x122;
cell[15] = x123;
cell[16] = x124;
cell[17] = x149*(x111*x21 + x164) + x36*(x24 + x85);
cell[18] = x125;
cell[1] = x108;
cell[2] = x109;
cell[3] = x110;
cell[4] = x112;
cell[5] = x113;
cell[6] = x114;
cell[7] = x115;
cell[8] = x116;
cell[9] = x117;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = V{0.5}*x155*(x109 + V{2}*x110 + x112 + x113 + V{2}*x114 + x119 + x121 + x122 + V{2}*x124 + V{2}*x125 + x180) + V{0.5}*x157*(V{2}*x109 + x110 + V{2}*x112 + x114 + x115 + V{2}*x117 + x120 + V{2}*x122 + x123 + x124 + x180);
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x146;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x147;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x148;
return { x105, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
