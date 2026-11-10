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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<0, 1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<0, 1, 1>, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x144 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x145 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x146 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x147 = parameters.template get<descriptors::OMEGA>();
auto x148 = ((x144)*(x144));
auto x149 = ((x145)*(x145));
auto x150 = ((x146)*(x146));
auto x152 = cell[15] + cell[16] + cell[6] + cell[7];
auto x153 = cell[10] + x152;
auto x154 = cell[12] + cell[3];
auto x157 = cell[11] + cell[2];
auto x0 = V{3};
auto x160 = ((x0)*(x0));
auto x1 = x148 + x149 + x150;
auto x2 = V{0.5}*x0*x1;
auto x161 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = cell[0] + cell[1];
auto x155 = V{2}*cell[17] + x5 + V{1};
auto x6 = x145 + V{1};
auto x151 = V{1} / (x6);
auto x7 = x146 + V{1};
auto x156 = V{1} / (x7);
auto x8 = V{0.5}*x151*(V{2}*cell[11] + V{2}*cell[13] + V{2}*cell[18] + V{2}*cell[5] + x153 + x154 + x155) + V{0.5}*x156*(cell[10] + V{2}*cell[12] + cell[13] + cell[14] + V{2}*cell[15] + cell[4] + cell[5] + V{2}*cell[7] + V{2}*cell[9] + x155 + x157);
auto x9 = V{1}/x0;
auto x10 = -x9;
auto x11 = x144*x8;
auto x12 = x9*(x8 + V{-1});
auto x158 = -x12;
auto x159 = cell[13] + cell[14] + cell[4] + cell[5] + x158;
auto x13 = cell[1] - x11*x144 + x153 + x159;
auto x14 = x10*x13;
auto x15 = -V{1}*x145;
auto x16 = -V{1}*x146;
auto x17 = x15 + x16;
auto x18 = V{0.5}*x160;
auto x19 = x0*x17 + x161 + x18*((x17)*(x17));
auto x20 = V{0.0277777777777778};
auto x162 = -x20;
auto x21 = x20*x8;
auto x22 = x162 + x19*x21;
auto x23 = x145 + x146;
auto x24 = x0*x23 + x161 + x18*((x23)*(x23));
auto x25 = -x20 + x21*x24;
auto x26 = cell[17] + x22 - x25;
auto x163 = cell[17] + cell[18] + cell[9] + x26;
auto x27 = x145*x8;
auto x28 = -x145*x27 + x157 + x159 + x163;
auto x29 = x10*x28;
auto x30 = x14 + x29;
auto x31 = -x150*x8 + x152 + x154 + x158 + x163;
auto x32 = x10*x31;
auto x33 = x4*(x160*(V{0.5}*x30 + V{0.5}*x32) + x3*x8 + V{-1});
auto x34 = V{1} - x147;
auto x35 = x144 + x16;
auto x36 = x0*x35 + x161 + x18*((x35)*(x35));
auto x37 = V{0.5}*x160*x20;
auto x38 = V{1} - x9;
auto x39 = x31*x38;
auto x40 = x13*x38;
auto x164 = x29 + x40;
auto x165 = x164 + x39;
auto x41 = V{-2};
auto x42 = cell[15] - V{1}*cell[16] + cell[6] - V{1}*cell[7] - x11*x146;
auto x43 = x37*(x165 + x41*x42);
auto x166 = x162 + x43;
auto x44 = x166 + x21*x36;
auto x45 = x144 + x146;
auto x46 = x0*x45 + x161 + x18*((x45)*(x45));
auto x47 = V{2};
auto x48 = x37*(x165 + x42*x47);
auto x167 = x162 + x48;
auto x49 = x167 + x21*x46;
auto x50 = V{0.0555555555555556};
auto x168 = -x50;
auto x51 = x50*x8;
auto x52 = x0*x146 + x150*x18 + x161;
auto x53 = V{0.5}*x160*x50;
auto x54 = x53*(x30 + x39);
auto x169 = x168 + x54;
auto x55 = x169 + x51*x52;
auto x56 = x0*x144 + x148*x18 + x161;
auto x57 = x53*(x164 + x32);
auto x170 = x168 + x57;
auto x58 = x170 + x51*x56;
auto x59 = -V{1}*x144;
auto x60 = x146 + x59;
auto x61 = x0*x60 + x161 + x18*((x60)*(x60));
auto x62 = x166 + x21*x61;
auto x63 = x16 + x59;
auto x64 = x0*x63 + x161 + x18*((x63)*(x63));
auto x65 = x167 + x21*x64;
auto x66 = x0*x16 + x161 + x18*((x16)*(x16));
auto x67 = x169 + x51*x66;
auto x68 = x0*x59 + x161 + x18*((x59)*(x59));
auto x69 = x170 + x51*x68;
auto x70 = x33 + x69;
auto x71 = x145 + x16;
auto x72 = x0*x71 + x161 + x18*((x71)*(x71));
auto x73 = x28*x38;
auto x173 = x32 + x40 + x73;
auto x74 = x14 + x73;
auto x171 = x39 + x74;
auto x75 = cell[17] - V{1}*cell[18] - V{1}*cell[9] - x146*x27 + x26;
auto x76 = x37*(x171 + x41*x75);
auto x172 = x162 + x76;
auto x77 = x172 + x21*x72;
auto x78 = x37*(x171 + x47*x75);
auto x79 = x25 + x78;
auto x177 = x58 + x70 + V{2}*x79 + V{1};
auto x80 = x144 + x145;
auto x81 = x0*x80 + x161 + x18*((x80)*(x80));
auto x82 = cell[13] - V{1}*cell[14] + cell[4] - V{1}*cell[5] - x11*x145;
auto x83 = x37*(x173 + x47*x82);
auto x174 = x162 + x83;
auto x84 = x174 + x21*x81;
auto x85 = x145 + x59;
auto x86 = x0*x85 + x161 + x18*((x85)*(x85));
auto x87 = x37*(x173 + x41*x82);
auto x175 = x162 + x87;
auto x88 = x175 + x21*x86;
auto x89 = x0*x145 + x149*x18 + x161;
auto x90 = x53*(x32 + x74);
auto x176 = x168 + x90;
auto x91 = x176 + x51*x89;
auto x92 = x144 + x15;
auto x93 = x0*x92 + x161 + x18*((x92)*(x92));
auto x94 = x175 + x21*x93;
auto x95 = x15 + x59;
auto x96 = x0*x95 + x161 + x18*((x95)*(x95));
auto x97 = x174 + x21*x96;
auto x98 = x0*x15 + x161 + x18*((x15)*(x15));
auto x99 = x176 + x51*x98;
auto x100 = x146 + x15;
auto x101 = x0*x100 + x161 + x18*((x100)*(x100));
auto x102 = x101*x21 + x172;
auto x103 = V{0.5}*x151*(x177 + x44 + x49 + x55 + x62 + x65 + x67 + V{2}*x77 + V{2}*x84 + V{2}*x88 + V{2}*x91) + V{0.5}*x156*(V{2}*x102 + x177 + V{2}*x49 + V{2}*x55 + V{2}*x62 + x84 + x88 + x91 + x94 + x97 + x99);
auto x104 = x147*x4*(x103*x3 + V{-1}) + x33*x34;
auto x105 = x103*x50;
auto x106 = x147*(x105*x68 + x168) + x34*x69;
auto x107 = x147*(x105*x98 + x168) + x34*x99;
auto x108 = x147*(x105*x66 + x168) + x34*x67;
auto x109 = x103*x20;
auto x110 = x147*(x109*x96 + x162) + x34*x97;
auto x111 = x147*(x109*x86 + x162) + x34*x88;
auto x112 = x147*(x109*x64 + x162) + x34*x65;
auto x113 = x147*(x109*x61 + x162) + x34*x62;
auto x114 = x102*x34 + x147*(x101*x109 + x162);
auto x115 = x147*(x105*x56 + x168) + x34*x58;
auto x116 = x147*(x105*x89 + x168) + x34*x91;
auto x117 = x147*(x105*x52 + x168) + x34*x55;
auto x118 = x147*(x109*x81 + x162) + x34*x84;
auto x119 = x147*(x109*x93 + x162) + x34*x94;
auto x120 = x147*(x109*x46 + x162) + x34*x49;
auto x121 = x147*(x109*x36 + x162) + x34*x44;
auto x122 = x147*(x109*x24 + x162) + x34*x79;
auto x123 = x147*(x109*x72 + x162) + x34*x77;
auto x124 = x104 + x106;
auto x178 = x115 + V{2}*x122 + x124 + V{1};
cell[0] = x104;
cell[10] = x115;
cell[11] = x116;
cell[12] = x117;
cell[13] = x118;
cell[14] = x119;
cell[15] = x120;
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
cell[8] = x147*(x109*x19 + x162) + x34*(x22 + x78);
cell[9] = x114;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = V{0.5}*x151*(x108 + V{2}*x111 + x112 + x113 + V{2}*x116 + x117 + V{2}*x118 + x120 + x121 + V{2}*x123 + x178) + V{0.5}*x156*(x107 + x110 + x111 + V{2}*x113 + V{2}*x114 + x116 + V{2}*x117 + x118 + x119 + V{2}*x120 + x178);
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x144;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x145;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x146;
return { x103, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
