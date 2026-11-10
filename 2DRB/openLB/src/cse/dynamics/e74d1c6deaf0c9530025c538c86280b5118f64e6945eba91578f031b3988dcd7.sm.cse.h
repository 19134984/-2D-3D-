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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<1, -1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<1, -1, 1>, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x142 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x143 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x144 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x145 = parameters.template get<descriptors::OMEGA>();
auto x146 = ((x142)*(x142));
auto x147 = ((x143)*(x143));
auto x148 = ((x144)*(x144));
auto x150 = cell[13] + cell[14] + cell[4] + cell[5];
auto x151 = cell[10] + cell[1] + x150;
auto x152 = cell[11] + cell[2];
auto x153 = cell[0] + V{2}*cell[16] + x152 + V{1};
auto x155 = cell[17] + cell[18] + cell[8] + cell[9];
auto x156 = cell[12] + x155;
auto x0 = V{3};
auto x157 = ((x0)*(x0));
auto x1 = x146 + x147 + x148;
auto x2 = V{0.5}*x0*x1;
auto x158 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = cell[3] + cell[6];
auto x6 = -V{1}*x144;
auto x7 = x6 + V{1};
auto x149 = V{1} / (x7);
auto x8 = x142 + V{1};
auto x154 = V{1} / (x8);
auto x9 = V{0.5}*x149*(V{2}*cell[18] + V{2}*cell[8] + x151 + x153 + V{2}*x5) + V{0.5}*x154*(V{2}*cell[10] + V{2}*cell[13] + V{2}*cell[14] + V{2}*cell[15] + cell[3] + x153 + x156);
auto x10 = V{1}/x0;
auto x11 = -x10;
auto x12 = -V{1}*x142;
auto x13 = x12 + x144;
auto x14 = V{0.5}*x157;
auto x15 = x0*x13 + x14*((x13)*(x13)) + x158;
auto x16 = V{0.0277777777777778};
auto x159 = -x16;
auto x17 = x16*x9;
auto x18 = x15*x17 + x159;
auto x19 = x142 + x6;
auto x20 = x0*x19 + x14*((x19)*(x19)) + x158;
auto x21 = -x16 + x17*x20;
auto x22 = cell[16] + x18 - x21;
auto x23 = x142*x9;
auto x24 = x10*(x9 + V{-1});
auto x160 = -x24;
auto x161 = cell[15] + cell[16] + x160 + x22;
auto x25 = cell[6] - x142*x23 + x151 + x161;
auto x26 = x11*x25;
auto x27 = x143*x9;
auto x28 = -x143*x27 + x150 + x152 + x155 + x160;
auto x29 = x11*x28;
auto x30 = x26 + x29;
auto x31 = -x148*x9 + x156 + x161 + x5;
auto x32 = x11*x31;
auto x33 = x4*(x157*(V{0.5}*x30 + V{0.5}*x32) + x3*x9 + V{-1});
auto x34 = V{1} - x145;
auto x35 = -V{1}*x143;
auto x36 = x142 + x35;
auto x37 = x0*x36 + x14*((x36)*(x36)) + x158;
auto x38 = V{0.5}*x157*x16;
auto x39 = V{1} - x10;
auto x40 = x28*x39;
auto x41 = x25*x39;
auto x162 = x32 + x41;
auto x163 = x162 + x40;
auto x42 = V{-2};
auto x43 = cell[13] - V{1}*cell[14] + cell[4] - V{1}*cell[5] - x143*x23;
auto x44 = x38*(x163 + x42*x43);
auto x164 = x159 + x44;
auto x45 = x164 + x17*x37;
auto x46 = x142 + x143;
auto x47 = x0*x46 + x14*((x46)*(x46)) + x158;
auto x48 = V{2};
auto x49 = x38*(x163 + x43*x48);
auto x165 = x159 + x49;
auto x50 = x165 + x17*x47;
auto x51 = V{0.0555555555555556};
auto x166 = -x51;
auto x52 = x51*x9;
auto x53 = x0*x143 + x14*x147 + x158;
auto x54 = V{0.5}*x157*x51;
auto x55 = x26 + x40;
auto x56 = x54*(x32 + x55);
auto x167 = x166 + x56;
auto x57 = x167 + x52*x53;
auto x58 = x0*x142 + x14*x146 + x158;
auto x59 = x54*(x162 + x29);
auto x168 = x166 + x59;
auto x60 = x168 + x52*x58;
auto x61 = x12 + x143;
auto x62 = x0*x61 + x14*((x61)*(x61)) + x158;
auto x63 = x164 + x17*x62;
auto x64 = x12 + x35;
auto x65 = x0*x64 + x14*((x64)*(x64)) + x158;
auto x66 = x165 + x17*x65;
auto x67 = x0*x35 + x14*((x35)*(x35)) + x158;
auto x68 = x167 + x52*x67;
auto x69 = x0*x12 + x14*((x12)*(x12)) + x158;
auto x70 = x168 + x52*x69;
auto x71 = x143 + x6;
auto x72 = x0*x71 + x14*((x71)*(x71)) + x158;
auto x73 = x31*x39;
auto x169 = x55 + x73;
auto x171 = x29 + x41 + x73;
auto x74 = cell[17] - V{1}*cell[18] + cell[8] - V{1}*cell[9] - x144*x27;
auto x75 = x38*(x169 + x42*x74);
auto x170 = x159 + x75;
auto x76 = x17*x72 + x170;
auto x77 = cell[15] - V{1}*cell[16] + cell[6] - x144*x23 - V{1}*x22;
auto x78 = x38*(x171 + x42*x77);
auto x79 = x21 + x78;
auto x175 = x33 + x57 + x68 + V{2}*x79 + V{1};
auto x80 = x35 + x6;
auto x81 = x0*x80 + x14*((x80)*(x80)) + x158;
auto x82 = x38*(x169 + x48*x74);
auto x172 = x159 + x82;
auto x83 = x17*x81 + x172;
auto x84 = x0*x6 + x14*((x6)*(x6)) + x158;
auto x85 = x54*(x30 + x73);
auto x173 = x166 + x85;
auto x86 = x173 + x52*x84;
auto x87 = x12 + x6;
auto x88 = x0*x87 + x14*((x87)*(x87)) + x158;
auto x89 = x38*(x171 + x48*x77);
auto x174 = x159 + x89;
auto x90 = x17*x88 + x174;
auto x91 = x143 + x144;
auto x92 = x0*x91 + x14*((x91)*(x91)) + x158;
auto x93 = x17*x92 + x172;
auto x94 = x0*x144 + x14*x148 + x158;
auto x95 = x173 + x52*x94;
auto x96 = x144 + x35;
auto x97 = x0*x96 + x14*((x96)*(x96)) + x158;
auto x98 = x17*x97 + x170;
auto x99 = x142 + x144;
auto x100 = x0*x99 + x14*((x99)*(x99)) + x158;
auto x101 = x100*x17 + x174;
auto x102 = V{0.5}*x149*(x175 + x45 + x50 + x60 + x63 + x66 + x70 + V{2}*x76 + V{2}*x83 + V{2}*x86 + V{2}*x90) + V{0.5}*x154*(V{2}*x101 + x175 + V{2}*x45 + V{2}*x50 + V{2}*x60 + x76 + x83 + x86 + x93 + x95 + x98);
auto x103 = x145*x4*(x102*x3 + V{-1}) + x33*x34;
auto x104 = x102*x51;
auto x105 = x145*(x104*x69 + x166) + x34*x70;
auto x106 = x145*(x104*x67 + x166) + x34*x68;
auto x107 = x145*(x104*x84 + x166) + x34*x86;
auto x108 = x102*x16;
auto x109 = x145*(x108*x65 + x159) + x34*x66;
auto x110 = x145*(x108*x62 + x159) + x34*x63;
auto x111 = x145*(x108*x88 + x159) + x34*x90;
auto x112 = x145*(x108*x81 + x159) + x34*x83;
auto x113 = x145*(x108*x97 + x159) + x34*x98;
auto x114 = x145*(x104*x58 + x166) + x34*x60;
auto x115 = x145*(x104*x53 + x166) + x34*x57;
auto x116 = x145*(x104*x94 + x166) + x34*x95;
auto x117 = x145*(x108*x47 + x159) + x34*x50;
auto x118 = x145*(x108*x37 + x159) + x34*x45;
auto x119 = x101*x34 + x145*(x100*x108 + x159);
auto x120 = x145*(x108*x20 + x159) + x34*x79;
auto x176 = x103 + x106 + x115 + V{2}*x120 + V{1};
auto x121 = x145*(x108*x92 + x159) + x34*x93;
auto x122 = x145*(x108*x72 + x159) + x34*x76;
cell[0] = x103;
cell[10] = x114;
cell[11] = x115;
cell[12] = x116;
cell[13] = x117;
cell[14] = x118;
cell[15] = x119;
cell[16] = x120;
cell[17] = x121;
cell[18] = x122;
cell[1] = x105;
cell[2] = x106;
cell[3] = x107;
cell[4] = x109;
cell[5] = x110;
cell[6] = x111;
cell[7] = x145*(x108*x15 + x159) + x34*(x18 + x78);
cell[8] = x112;
cell[9] = x113;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = V{0.5}*x149*(x105 + V{2}*x107 + x109 + x110 + V{2}*x111 + V{2}*x112 + x114 + x117 + x118 + V{2}*x122 + x176) + V{0.5}*x154*(x107 + x112 + x113 + V{2}*x114 + x116 + V{2}*x117 + V{2}*x118 + V{2}*x119 + x121 + x122 + x176);
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x142;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x143;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x144;
return { x102, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
