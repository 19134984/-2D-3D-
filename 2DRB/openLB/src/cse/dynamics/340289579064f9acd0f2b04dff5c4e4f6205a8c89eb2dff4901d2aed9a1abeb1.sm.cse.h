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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerCornerDensity3D<1, 1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerCornerStress3D<1, 1, -1>, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x145 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x146 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x147 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x148 = parameters.template get<descriptors::OMEGA>();
auto x149 = ((x145)*(x145));
auto x150 = ((x146)*(x146));
auto x151 = ((x147)*(x147));
auto x152 = cell[3] + cell[6];
auto x154 = cell[13] + cell[14] + cell[4] + cell[5];
auto x155 = cell[10] + x154;
auto x157 = cell[11] + cell[2];
auto x158 = V{2}*cell[16] + x157;
auto x160 = V{2}*cell[13];
auto x162 = cell[17] + cell[18] + cell[8] + cell[9];
auto x163 = cell[12] + x162;
auto x0 = V{3};
auto x166 = ((x0)*(x0));
auto x1 = x149 + x150 + x151;
auto x2 = V{0.5}*x0*x1;
auto x167 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = cell[0] + cell[1];
auto x156 = V{2}*cell[18] + x5 + V{1};
auto x6 = x152;
auto x7 = -V{1}*x147;
auto x8 = x7 + V{1};
auto x153 = V{1} / (x8);
auto x9 = x145 + V{1};
auto x161 = V{1} / (x9);
auto x10 = x146 + V{1};
auto x159 = V{1} / (x10);
auto x11 = V{0.333333333333333}*x153*(V{2}*cell[8] + x155 + x156 + x158 + V{2}*x6) + V{0.333333333333333}*x159*(cell[10] + V{2}*cell[11] + cell[12] + cell[15] + cell[16] + V{2}*cell[17] + V{2}*cell[5] + cell[7] + x152 + x156 + x160) + V{0.333333333333333}*x161*(cell[0] + V{2}*cell[10] + V{2}*cell[14] + V{2}*cell[15] + cell[3] + x158 + x160 + x163 + V{1});
auto x12 = V{1}/x0;
auto x13 = -x12;
auto x14 = x11*x145;
auto x15 = x12*(x11 + V{-1});
auto x164 = -x15;
auto x165 = cell[15] + cell[16] + cell[7] + x164;
auto x16 = cell[1] + cell[6] - x14*x145 + x155 + x165;
auto x17 = x13*x16;
auto x18 = x11*x146;
auto x19 = -x146*x18 + x154 + x157 + x162 + x164;
auto x20 = x13*x19;
auto x21 = x17 + x20;
auto x22 = -x11*x151 + x163 + x165 + x6;
auto x23 = x13*x22;
auto x24 = x4*(x11*x3 + x166*(V{0.5}*x21 + V{0.5}*x23) + V{-1});
auto x25 = V{1} - x148;
auto x26 = V{0.0277777777777778};
auto x170 = -x26;
auto x27 = x11*x26;
auto x28 = -V{1}*x146;
auto x29 = x145 + x28;
auto x30 = V{0.5}*x166;
auto x31 = x0*x29 + x167 + x30*((x29)*(x29));
auto x32 = V{0.5}*x166*x26;
auto x33 = V{1} - x12;
auto x34 = x19*x33;
auto x35 = x16*x33;
auto x168 = x23 + x35;
auto x169 = x168 + x34;
auto x36 = V{-2};
auto x37 = cell[13] - V{1}*cell[14] + cell[4] - V{1}*cell[5] - x14*x146;
auto x38 = x32*(x169 + x36*x37);
auto x171 = x170 + x38;
auto x39 = x171 + x27*x31;
auto x40 = x145 + x146;
auto x41 = x0*x40 + x167 + x30*((x40)*(x40));
auto x42 = V{2};
auto x43 = x32*(x169 + x37*x42);
auto x172 = x170 + x43;
auto x44 = x172 + x27*x41;
auto x45 = V{0.0555555555555556};
auto x173 = -x45;
auto x46 = x11*x45;
auto x47 = x0*x146 + x150*x30 + x167;
auto x48 = V{0.5}*x166*x45;
auto x49 = x17 + x34;
auto x50 = x48*(x23 + x49);
auto x174 = x173 + x50;
auto x51 = x174 + x46*x47;
auto x52 = x0*x145 + x149*x30 + x167;
auto x53 = x48*(x168 + x20);
auto x175 = x173 + x53;
auto x54 = x175 + x46*x52;
auto x55 = -V{1}*x145;
auto x56 = x146 + x55;
auto x57 = x0*x56 + x167 + x30*((x56)*(x56));
auto x58 = x171 + x27*x57;
auto x59 = x28 + x55;
auto x60 = x0*x59 + x167 + x30*((x59)*(x59));
auto x61 = x172 + x27*x60;
auto x62 = x0*x28 + x167 + x30*((x28)*(x28));
auto x63 = x174 + x46*x62;
auto x64 = x0*x55 + x167 + x30*((x55)*(x55));
auto x65 = x175 + x46*x64;
auto x66 = x24 + x65;
auto x67 = x146 + x7;
auto x68 = x0*x67 + x167 + x30*((x67)*(x67));
auto x69 = x22*x33;
auto x176 = x49 + x69;
auto x178 = x20 + x35 + x69;
auto x70 = cell[17] - V{1}*cell[18] + cell[8] - V{1}*cell[9] - x147*x18;
auto x71 = x32*(x176 + x36*x70);
auto x177 = x170 + x71;
auto x72 = x177 + x27*x68;
auto x183 = x54 + x66 + V{2}*x72 + V{1};
auto x73 = x145 + x7;
auto x74 = x0*x73 + x167 + x30*((x73)*(x73));
auto x75 = cell[15] - V{1}*cell[16] + cell[6] - V{1}*cell[7] - x14*x147;
auto x76 = x32*(x178 + x36*x75);
auto x179 = x170 + x76;
auto x77 = x179 + x27*x74;
auto x185 = x51 + x63 + V{2}*x77;
auto x78 = x28 + x7;
auto x79 = x0*x78 + x167 + x30*((x78)*(x78));
auto x80 = x32*(x176 + x42*x70);
auto x180 = x170 + x80;
auto x81 = x180 + x27*x79;
auto x82 = x0*x7 + x167 + x30*((x7)*(x7));
auto x83 = x48*(x21 + x69);
auto x181 = x173 + x83;
auto x84 = x181 + x46*x82;
auto x85 = x55 + x7;
auto x86 = x0*x85 + x167 + x30*((x85)*(x85));
auto x87 = x32*(x178 + x42*x75);
auto x182 = x170 + x87;
auto x88 = x182 + x27*x86;
auto x89 = x146 + x147;
auto x90 = x0*x89 + x167 + x30*((x89)*(x89));
auto x91 = x180 + x27*x90;
auto x92 = x0*x147 + x151*x30 + x167;
auto x93 = x181 + x46*x92;
auto x184 = V{2}*x44 + x84 + x93;
auto x94 = x147 + x28;
auto x95 = x0*x94 + x167 + x30*((x94)*(x94));
auto x96 = x177 + x27*x95;
auto x97 = x145 + x147;
auto x98 = x0*x97 + x167 + x30*((x97)*(x97));
auto x99 = x182 + x27*x98;
auto x100 = x147 + x55;
auto x101 = x0*x100 + x167 + x30*((x100)*(x100));
auto x102 = x101*x27 + x179;
auto x103 = V{0.333333333333333}*x153*(x183 + x185 + x39 + x44 + x58 + x61 + V{2}*x81 + V{2}*x84 + V{2}*x88) + V{0.333333333333333}*x159*(x102 + x183 + x184 + V{2}*x51 + V{2}*x58 + x77 + x88 + V{2}*x91 + x99) + V{0.333333333333333}*x161*(x184 + x185 + x24 + V{2}*x39 + V{2}*x54 + x72 + x81 + x91 + x96 + V{2}*x99 + V{1});
auto x104 = x148*x4*(x103*x3 + V{-1}) + x24*x25;
auto x105 = x103*x45;
auto x106 = x148*(x105*x64 + x173) + x25*x65;
auto x107 = x148*(x105*x62 + x173) + x25*x63;
auto x108 = x148*(x105*x82 + x173) + x25*x84;
auto x109 = x103*x26;
auto x110 = x148*(x109*x60 + x170) + x25*x61;
auto x111 = x148*(x109*x57 + x170) + x25*x58;
auto x112 = x148*(x109*x86 + x170) + x25*x88;
auto x113 = x102*x25 + x148*(x101*x109 + x170);
auto x114 = x148*(x109*x79 + x170) + x25*x81;
auto x115 = x148*(x109*x95 + x170) + x25*x96;
auto x116 = x148*(x105*x52 + x173) + x25*x54;
auto x117 = x148*(x105*x47 + x173) + x25*x51;
auto x118 = x148*(x105*x92 + x173) + x25*x93;
auto x119 = x148*(x109*x41 + x170) + x25*x44;
auto x187 = x108 + x118 + V{2}*x119;
auto x120 = x148*(x109*x31 + x170) + x25*x39;
auto x121 = x148*(x109*x98 + x170) + x25*x99;
auto x122 = x148*(x109*x74 + x170) + x25*x77;
auto x188 = x107 + x117 + V{2}*x122;
auto x123 = x148*(x109*x90 + x170) + x25*x91;
auto x124 = x148*(x109*x68 + x170) + x25*x72;
auto x125 = x104 + x106;
auto x186 = x116 + V{2}*x124 + x125 + V{1};
cell[0] = x104;
cell[10] = x116;
cell[11] = x117;
cell[12] = x118;
cell[13] = x119;
cell[14] = x120;
cell[15] = x121;
cell[16] = x122;
cell[17] = x123;
cell[18] = x124;
cell[1] = x106;
cell[2] = x107;
cell[3] = x108;
cell[4] = x110;
cell[5] = x111;
cell[6] = x112;
cell[7] = x113;
cell[8] = x114;
cell[9] = x115;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = V{0.333333333333333}*x153*(V{2}*x108 + x110 + x111 + V{2}*x112 + V{2}*x114 + x119 + x120 + x186 + x188) + V{0.333333333333333}*x159*(V{2}*x111 + x112 + x113 + V{2}*x117 + x121 + x122 + V{2}*x123 + x186 + x187) + V{0.333333333333333}*x161*(x104 + x114 + x115 + V{2}*x116 + V{2}*x120 + V{2}*x121 + x123 + x124 + x187 + x188 + V{1});
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x145;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x146;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x147;
return { x103, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
