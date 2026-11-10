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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<0, -1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<0, -1, 1>, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x145 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x146 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x147 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x148 = parameters.template get<descriptors::OMEGA>();
auto x149 = ((x145)*(x145));
auto x150 = ((x146)*(x146));
auto x151 = ((x147)*(x147));
auto x152 = cell[2] + cell[4];
auto x154 = cell[15] + cell[16] + cell[6] + cell[7];
auto x155 = cell[10] + x154;
auto x156 = cell[12] + cell[3];
auto x0 = V{3};
auto x161 = ((x0)*(x0));
auto x1 = x149 + x150 + x151;
auto x2 = V{0.5}*x0*x1;
auto x162 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = cell[0] + cell[1];
auto x157 = V{2}*cell[9] + x5 + V{1};
auto x6 = x152;
auto x7 = -V{1}*x146;
auto x8 = x7 + V{1};
auto x153 = V{1} / (x8);
auto x9 = x147 + V{1};
auto x158 = V{1} / (x9);
auto x10 = V{0.5}*x153*(V{2}*cell[14] + V{2}*cell[8] + x155 + x156 + x157 + V{2}*x6) + V{0.5}*x158*(cell[10] + cell[11] + V{2}*cell[12] + cell[13] + cell[14] + V{2}*cell[15] + V{2}*cell[17] + cell[5] + V{2}*cell[7] + x152 + x157);
auto x11 = V{1}/x0;
auto x12 = -x11;
auto x13 = x10*x145;
auto x14 = x11*(x10 + V{-1});
auto x159 = -x14;
auto x160 = cell[13] + cell[14] + cell[5] + x159;
auto x15 = cell[1] + cell[4] - x13*x145 + x155 + x160;
auto x16 = x12*x15;
auto x17 = -V{1}*x147;
auto x18 = x146 + x17;
auto x19 = V{0.5}*x161;
auto x20 = x0*x18 + x162 + x19*((x18)*(x18));
auto x21 = V{0.0277777777777778};
auto x163 = -x21;
auto x22 = x10*x21;
auto x23 = x163 + x20*x22;
auto x24 = x147 + x7;
auto x25 = x0*x24 + x162 + x19*((x24)*(x24));
auto x26 = -x21 + x22*x25;
auto x27 = cell[9] + x23 - x26;
auto x164 = cell[17] + cell[8] + cell[9] + x27;
auto x28 = x10*x146;
auto x29 = cell[11] - x146*x28 + x160 + x164 + x6;
auto x30 = x12*x29;
auto x31 = x16 + x30;
auto x32 = -x10*x151 + x154 + x156 + x159 + x164;
auto x33 = x12*x32;
auto x34 = x4*(x10*x3 + x161*(V{0.5}*x31 + V{0.5}*x33) + V{-1});
auto x35 = V{1} - x148;
auto x36 = x145 + x17;
auto x37 = x0*x36 + x162 + x19*((x36)*(x36));
auto x38 = V{0.5}*x161*x21;
auto x39 = V{1} - x11;
auto x40 = x32*x39;
auto x41 = x15*x39;
auto x165 = x30 + x41;
auto x166 = x165 + x40;
auto x42 = V{-2};
auto x43 = cell[15] - V{1}*cell[16] + cell[6] - V{1}*cell[7] - x13*x147;
auto x44 = x38*(x166 + x42*x43);
auto x167 = x163 + x44;
auto x45 = x167 + x22*x37;
auto x46 = x145 + x147;
auto x47 = x0*x46 + x162 + x19*((x46)*(x46));
auto x48 = V{2};
auto x49 = x38*(x166 + x43*x48);
auto x168 = x163 + x49;
auto x50 = x168 + x22*x47;
auto x51 = V{0.0555555555555556};
auto x169 = -x51;
auto x52 = x10*x51;
auto x53 = x0*x147 + x151*x19 + x162;
auto x54 = V{0.5}*x161*x51;
auto x55 = x54*(x31 + x40);
auto x170 = x169 + x55;
auto x56 = x170 + x52*x53;
auto x57 = x0*x145 + x149*x19 + x162;
auto x58 = x54*(x165 + x33);
auto x171 = x169 + x58;
auto x59 = x171 + x52*x57;
auto x60 = -V{1}*x145;
auto x61 = x147 + x60;
auto x62 = x0*x61 + x162 + x19*((x61)*(x61));
auto x63 = x167 + x22*x62;
auto x64 = x17 + x60;
auto x65 = x0*x64 + x162 + x19*((x64)*(x64));
auto x66 = x168 + x22*x65;
auto x67 = x0*x17 + x162 + x19*((x17)*(x17));
auto x68 = x170 + x52*x67;
auto x69 = x0*x60 + x162 + x19*((x60)*(x60));
auto x70 = x171 + x52*x69;
auto x71 = x34 + x70;
auto x72 = x145 + x7;
auto x73 = x0*x72 + x162 + x19*((x72)*(x72));
auto x74 = x29*x39;
auto x172 = x33 + x41 + x74;
auto x75 = cell[13] - V{1}*cell[14] + cell[4] - V{1}*cell[5] - x13*x146;
auto x76 = x38*(x172 + x42*x75);
auto x173 = x163 + x76;
auto x77 = x173 + x22*x73;
auto x78 = x16 + x74;
auto x174 = x40 + x78;
auto x79 = cell[17] + cell[8] - V{1}*cell[9] - x147*x28 - V{1}*x27;
auto x80 = x38*(x174 + x42*x79);
auto x81 = x26 + x80;
auto x178 = x59 + x71 + V{2}*x81 + V{1};
auto x82 = x17 + x7;
auto x83 = x0*x82 + x162 + x19*((x82)*(x82));
auto x84 = x38*(x174 + x48*x79);
auto x175 = x163 + x84;
auto x85 = x175 + x22*x83;
auto x86 = x0*x7 + x162 + x19*((x7)*(x7));
auto x87 = x54*(x33 + x78);
auto x176 = x169 + x87;
auto x88 = x176 + x52*x86;
auto x89 = x60 + x7;
auto x90 = x0*x89 + x162 + x19*((x89)*(x89));
auto x91 = x38*(x172 + x48*x75);
auto x177 = x163 + x91;
auto x92 = x177 + x22*x90;
auto x93 = x145 + x146;
auto x94 = x0*x93 + x162 + x19*((x93)*(x93));
auto x95 = x177 + x22*x94;
auto x96 = x0*x146 + x150*x19 + x162;
auto x97 = x176 + x52*x96;
auto x98 = x146 + x60;
auto x99 = x0*x98 + x162 + x19*((x98)*(x98));
auto x100 = x173 + x22*x99;
auto x101 = x146 + x147;
auto x102 = x0*x101 + x162 + x19*((x101)*(x101));
auto x103 = x102*x22 + x175;
auto x104 = V{0.5}*x153*(x178 + x45 + x50 + x56 + x63 + x66 + x68 + V{2}*x77 + V{2}*x85 + V{2}*x88 + V{2}*x92) + V{0.5}*x158*(x100 + V{2}*x103 + x178 + V{2}*x50 + V{2}*x56 + V{2}*x63 + x77 + x88 + x92 + x95 + x97);
auto x105 = x148*x4*(x104*x3 + V{-1}) + x34*x35;
auto x106 = x104*x51;
auto x107 = x148*(x106*x69 + x169) + x35*x70;
auto x108 = x148*(x106*x86 + x169) + x35*x88;
auto x109 = x148*(x106*x67 + x169) + x35*x68;
auto x110 = x104*x21;
auto x111 = x148*(x110*x90 + x163) + x35*x92;
auto x112 = x100*x35 + x148*(x110*x99 + x163);
auto x113 = x148*(x110*x65 + x163) + x35*x66;
auto x114 = x148*(x110*x62 + x163) + x35*x63;
auto x115 = x148*(x110*x83 + x163) + x35*x85;
auto x116 = x148*(x110*x25 + x163) + x35*x81;
auto x117 = x148*(x106*x57 + x169) + x35*x59;
auto x118 = x148*(x106*x96 + x169) + x35*x97;
auto x119 = x148*(x106*x53 + x169) + x35*x56;
auto x120 = x148*(x110*x94 + x163) + x35*x95;
auto x121 = x148*(x110*x73 + x163) + x35*x77;
auto x122 = x148*(x110*x47 + x163) + x35*x50;
auto x123 = x148*(x110*x37 + x163) + x35*x45;
auto x124 = x103*x35 + x148*(x102*x110 + x163);
auto x125 = x105 + x107;
auto x179 = V{2}*x116 + x117 + x125 + V{1};
cell[0] = x105;
cell[10] = x117;
cell[11] = x118;
cell[12] = x119;
cell[13] = x120;
cell[14] = x121;
cell[15] = x122;
cell[16] = x123;
cell[17] = x124;
cell[18] = x148*(x110*x20 + x163) + x35*(x23 + x80);
cell[1] = x107;
cell[2] = x108;
cell[3] = x109;
cell[4] = x111;
cell[5] = x112;
cell[6] = x113;
cell[7] = x114;
cell[8] = x115;
cell[9] = x116;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = V{0.5}*x153*(V{2}*x108 + x109 + V{2}*x111 + x113 + x114 + V{2}*x115 + x119 + V{2}*x121 + x122 + x123 + x179) + V{0.5}*x158*(x108 + x111 + x112 + V{2}*x114 + x118 + V{2}*x119 + x120 + x121 + V{2}*x122 + V{2}*x124 + x179);
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x145;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x146;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x147;
return { x104, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
