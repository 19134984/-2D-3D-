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
struct CSE<dynamics::Tuple<T, descriptors::D3Q27<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, forcing::Guo<momenta::ForcedWithStress> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x142 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x143 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x144 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x145 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x146 = parameters.template get<descriptors::OMEGA>();
auto x147 = cell[12] + cell[25];
auto x148 = cell[6] + x147;
auto x149 = cell[11] + cell[24];
auto x150 = cell[4] + x149;
auto x151 = cell[19] + cell[20];
auto x152 = cell[14] + cell[18] + x151;
auto x153 = cell[1] + x150 + x152;
auto x154 = cell[23] + cell[26];
auto x155 = cell[10] + x154;
auto x156 = cell[8] + x155;
auto x157 = cell[16] + cell[7] + cell[9];
auto x158 = cell[3] + x156 + x157;
auto x159 = cell[13] + cell[17] + cell[5];
auto x160 = cell[21] + cell[22];
auto x161 = cell[15] + x159 + x160;
auto x162 = cell[2] + x161;
auto x0 = V{1}/x146;
auto x1 = -V{1}*cell[26];
auto x2 = -V{1}*cell[24];
auto x3 = -V{1}*cell[20];
auto x164 = cell[19] + x3;
auto x4 = -V{1}*cell[13];
auto x170 = cell[10] + cell[23] + x1 + x4;
auto x5 = -V{1}*cell[11];
auto x171 = x2 + x5;
auto x6 = -V{1}*cell[7];
auto x7 = -V{1}*cell[6];
auto x8 = -V{1}*cell[8];
auto x167 = cell[23] + x1 + x8;
auto x9 = -V{1}*cell[10];
auto x10 = -V{1}*cell[12];
auto x166 = cell[25] + x10 + x7 + x9;
auto x11 = -V{1}*cell[22];
auto x165 = cell[13] + cell[21] + x11;
auto x12 = cell[0] + x148 + x153 + x158 + x162 + V{1};
auto x163 = V{1} / (x12);
auto x13 = V{0.5}*x144;
auto x14 = x13 + x163*(cell[11] - V{1}*cell[3] + x157 + x164 + x165 + x166 + x167 + x2);
auto x177 = ((x14)*(x14));
auto x15 = -x13 + x14;
auto x16 = -V{1}*cell[4];
auto x169 = cell[24] + x16 + x5;
auto x17 = -V{1}*cell[5];
auto x168 = cell[17] + x17;
auto x18 = V{0.5}*x142;
auto x19 = x163*(-V{1}*cell[1] + x152 + x154 + x166 + x168 + x169 + x4 + x6) + x18;
auto x178 = ((x19)*(x19));
auto x20 = -x18 + x19;
auto x21 = x12*x20;
auto x22 = V{0.5}*x12;
auto x173 = V{2}*x22;
auto x23 = x148 - x15*x21 + x164 + x170 + x171 + x22*(x142*x15 + x144*x20) + x6;
auto x24 = ((x23)*(x23));
auto x25 = -V{1}*cell[25];
auto x172 = x10 + x25;
auto x26 = -V{1}*cell[18];
auto x27 = -V{1}*cell[9];
auto x28 = V{0.5}*x143;
auto x29 = x163*(cell[12] - V{1}*cell[2] + x161 + x167 + x169 + x25 + x26 + x27 + x9) + x28;
auto x179 = ((x29)*(x29));
auto x30 = -x28 + x29;
auto x31 = x150 + x168 + x170 + x172 - x21*x30 + x22*(x142*x30 + x143*x20) + x26;
auto x32 = V{3};
auto x176 = ((x32)*(x32));
auto x180 = x19*x32;
auto x181 = x29*x32;
auto x182 = x14*x32;
auto x33 = V{1}*(x12 + V{-1})/x32;
auto x174 = -x33;
auto x175 = x148 + x174;
auto x34 = x142*x20;
auto x35 = cell[7] + x153 + x155 + x159 + x173*x34 + x175 - x20*x21;
auto x36 = x144*x15;
auto x37 = cell[13] - x12*((x15)*(x15)) + x149 + x151 + x158 + x160 + x173*x36 + x175;
auto x38 = x12*x30;
auto x39 = -x15*x38 + x156 + x165 + x171 + x172 + x22*(x143*x15 + x144*x30) + x27;
auto x40 = x143*x30;
auto x41 = cell[18] + cell[9] + x147 + x150 + x156 + x162 + x173*x40 + x174 - x30*x38;
auto x42 = V{2}*V{1} / (x0 + util::sqrt(V{2.828428}*x163*x176*((x145)*(x145))*util::sqrt(V{2}*x24 + V{2}*((x31)*(x31)) + ((x35)*(x35)) + ((x37)*(x37)) + V{2}*((x39)*(x39)) + ((x41)*(x41))) + ((x0)*(x0))));
auto x43 = V{1} - x42;
auto x44 = x177 + x178 + x179;
auto x45 = V{0.5}*x32*x44;
auto x184 = V{1} - x45;
auto x46 = V{0.296296296296296};
auto x47 = -x142*x180;
auto x48 = -x143*x181;
auto x49 = x47 + x48;
auto x50 = -x144*x182;
auto x185 = x48 + x50;
auto x186 = x47 + x50;
auto x51 = -x12*(V{0.5}*x42 + V{-1});
auto x52 = -V{1}*x19;
auto x53 = x176;
auto x183 = V{1}*x53;
auto x54 = V{0.5}*x53;
auto x55 = V{0.0740740740740741};
auto x56 = x12*x55;
auto x57 = -x32*(x19 + V{1});
auto x58 = x51*x55;
auto x59 = -V{1}*x29;
auto x60 = -x32*(x29 + V{1});
auto x188 = -x60;
auto x61 = -V{1}*x14;
auto x62 = -x32*(x14 + V{1});
auto x187 = -x62;
auto x63 = x52 + x59;
auto x64 = V{0.0185185185185185};
auto x65 = x12*x64;
auto x66 = -x183*x63;
auto x67 = x51*x64;
auto x68 = x29 + x52;
auto x69 = x53*x68;
auto x70 = -x32*(x29 + V{-1});
auto x71 = x52 + x61;
auto x72 = -x183*x71;
auto x73 = x14 + x52;
auto x74 = x53*x73;
auto x75 = -x32*(x14 + V{-1});
auto x76 = x59 + x61;
auto x77 = -x183*x76;
auto x78 = x14 + x59;
auto x79 = x53*x78;
auto x80 = x61 + x63;
auto x81 = V{0.00462962962962963};
auto x82 = x12*x81;
auto x83 = -x183*x80;
auto x84 = x51*x81;
auto x85 = x14 + x63;
auto x86 = x53*x85;
auto x87 = -V{1}*x86;
auto x88 = x61 + x68;
auto x89 = x53*x88;
auto x90 = -V{1}*x89;
auto x91 = x14 + x68;
auto x92 = x53*x91;
auto x93 = -x32*(x19 + V{-1});
auto x94 = x19 + x29;
auto x95 = x53*x94;
auto x96 = x19 + x59;
auto x97 = x53*x96;
auto x98 = x14 + x19;
auto x99 = x53*x98;
auto x100 = x19 + x61;
auto x101 = x100*x53;
auto x102 = x14 + x29;
auto x103 = x102*x53;
auto x104 = x29 + x61;
auto x105 = x104*x53;
auto x106 = x14 + x94;
auto x107 = x106*x53;
auto x108 = x61 + x94;
auto x109 = x108*x53;
auto x110 = x14 + x96;
auto x111 = x110*x53;
auto x112 = x61 + x96;
auto x113 = x112*x53;
auto x114 = -V{1}*x113;
auto x115 = cell[0]*x43 - x42*x46*(x12*(x45 + V{-1}) + V{1}) + x46*x51*(x49 + x50);
auto x116 = cell[10]*x43 - x42*(x81 - x82*(x184 + x32*x80 + x54*((x80)*(x80)))) + x84*(x142*(x57 + x83) + x143*(x60 + x83) + x144*(x62 + x83));
auto x117 = cell[11]*x43 - x42*(x81 - x82*(x184 + x32*x85 + x54*((x85)*(x85)))) + x84*(x142*(x57 + x87) + x143*(x60 + x87) + x144*(x75 + x86));
auto x118 = cell[12]*x43 - x42*(x81 - x82*(x184 + x32*x88 + x54*((x88)*(x88)))) + x84*(x142*(x57 + x90) + x143*(x70 + x89) + x144*(x62 + x90));
auto x119 = cell[13]*x43 - x42*(x81 - x82*(x184 + x32*x91 + x54*((x91)*(x91)))) + x84*(x142*(x57 - V{1}*x92) + x143*(x70 + x92) + x144*(x75 + x92));
auto x120 = cell[14]*x43 - x42*(x55 - x56*(x178*x54 + x180 + x184)) + x58*(x142*(x19*x53 + x93) + x185);
auto x121 = cell[15]*x43 - x42*(x55 - x56*(x179*x54 + x181 + x184)) + x58*(x143*(x29*x53 + x70) + x186);
auto x122 = cell[16]*x43 - x42*(x55 - x56*(x177*x54 + x182 + x184)) + x58*(x144*(x14*x53 + x75) + x49);
auto x123 = cell[17]*x43 - x42*(x64 - x65*(x184 + x32*x94 + x54*((x94)*(x94)))) + x67*(x142*(x93 + x95) + x143*(x70 + x95) + x50);
auto x124 = cell[18]*x43 - x42*(x64 - x65*(x184 + x32*x96 + x54*((x96)*(x96)))) + x67*(x142*(x93 + x97) + x143*(x60 - V{1}*x97) + x50);
auto x125 = cell[19]*x43 - x42*(x64 - x65*(x184 + x32*x98 + x54*((x98)*(x98)))) + x67*(x142*(x93 + x99) + x144*(x75 + x99) + x48);
auto x126 = cell[1]*x43 - x42*(x55 - x56*(x184 + x32*x52 + x54*((x52)*(x52)))) + x58*(-x142*(x183*x52 - x57) + x185);
auto x127 = cell[20]*x43 - x42*(x64 - x65*(x100*x32 + x184 + x54*((x100)*(x100)))) + x67*(x142*(x101 + x93) - x144*(V{1}*x101 + x187) + x48);
auto x128 = cell[21]*x43 - x42*(x64 - x65*(x102*x32 + x184 + x54*((x102)*(x102)))) + x67*(x143*(x103 + x70) + x144*(x103 + x75) + x47);
auto x129 = cell[22]*x43 - x42*(x64 - x65*(x104*x32 + x184 + x54*((x104)*(x104)))) + x67*(x143*(x105 + x70) - x144*(V{1}*x105 + x187) + x47);
auto x130 = cell[23]*x43 - x42*(x81 - x82*(x106*x32 + x184 + x54*((x106)*(x106)))) + x84*(x142*(x107 + x93) + x143*(x107 + x70) + x144*(x107 + x75));
auto x131 = cell[24]*x43 - x42*(x81 - x82*(x108*x32 + x184 + x54*((x108)*(x108)))) + x84*(x142*(x109 + x93) + x143*(x109 + x70) - x144*(V{1}*x109 + x187));
auto x132 = cell[25]*x43 - x42*(x81 - x82*(x110*x32 + x184 + x54*((x110)*(x110)))) + x84*(x142*(x111 + x93) - x143*(V{1}*x111 + x188) + x144*(x111 + x75));
auto x133 = cell[26]*x43 - x42*(x81 - x82*(x112*x32 + x184 + x54*((x112)*(x112)))) + x84*(x142*(x113 + x93) + x143*(x114 + x60) + x144*(x114 + x62));
auto x134 = cell[2]*x43 - x42*(x55 - x56*(x184 + x32*x59 + x54*((x59)*(x59)))) + x58*(-x143*(x183*x59 + x188) + x186);
auto x135 = -(-cell[3]*x43 + x42*(x55 - x56*(x184 + x32*x61 + x54*((x61)*(x61)))) + x58*(x144*(x183*x61 + x187) - x49));
auto x136 = cell[4]*x43 - x42*(x64 - x65*(x184 + x32*x63 + x54*((x63)*(x63)))) + x67*(x142*(x57 + x66) + x143*(x60 + x66) + x50);
auto x137 = cell[5]*x43 - x42*(x64 - x65*(x184 + x32*x68 + x54*((x68)*(x68)))) + x67*(x142*(x57 - V{1}*x69) + x143*(x69 + x70) + x50);
auto x138 = cell[6]*x43 - x42*(x64 - x65*(x184 + x32*x71 + x54*((x71)*(x71)))) + x67*(x142*(x57 + x72) + x144*(x62 + x72) + x48);
auto x139 = cell[7]*x43 - x42*(x64 - x65*(x184 + x32*x73 + x54*((x73)*(x73)))) + x67*(x142*(x57 - V{1}*x74) + x144*(x74 + x75) + x48);
auto x140 = cell[8]*x43 - x42*(x64 - x65*(x184 + x32*x76 + x54*((x76)*(x76)))) + x67*(x143*(x60 + x77) + x144*(x62 + x77) + x47);
auto x141 = cell[9]*x43 - x42*(x64 - x65*(x184 + x32*x78 + x54*((x78)*(x78)))) + x67*(x143*(x60 - V{1}*x79) + x144*(x75 + x79) + x47);
cell[0] = x115;
cell[10] = x116;
cell[11] = x117;
cell[12] = x118;
cell[13] = x119;
cell[14] = x120;
cell[15] = x121;
cell[16] = x122;
cell[17] = x123;
cell[18] = x124;
cell[19] = x125;
cell[1] = x126;
cell[20] = x127;
cell[21] = x128;
cell[22] = x129;
cell[23] = x130;
cell[24] = x131;
cell[25] = x132;
cell[26] = x133;
cell[2] = x134;
cell[3] = x135;
cell[4] = x136;
cell[5] = x137;
cell[6] = x138;
cell[7] = x139;
cell[8] = x140;
cell[9] = x141;
return { x12, x44 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
