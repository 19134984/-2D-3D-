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
struct CSE<dynamics::Tuple<T, descriptors::D3Q27<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, forcing::Guo<momenta::ForcedWithStress> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x177 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x178 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x179 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x180 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x181 = parameters.template get<descriptors::OMEGA>();
auto x182 = cell[12] + cell[25];
auto x183 = cell[6] + x182;
auto x184 = cell[11] + cell[24];
auto x185 = cell[4] + x184;
auto x186 = cell[19] + cell[20];
auto x187 = cell[14] + cell[18] + x186;
auto x188 = cell[1] + x185 + x187;
auto x189 = cell[23] + cell[26];
auto x190 = cell[10] + x189;
auto x191 = cell[8] + x190;
auto x192 = cell[16] + cell[7] + cell[9];
auto x193 = cell[3] + x191 + x192;
auto x194 = cell[13] + cell[17] + cell[5];
auto x195 = cell[21] + cell[22];
auto x196 = cell[15] + x194 + x195;
auto x197 = cell[2] + x196;
auto x0 = V{1}/x181;
auto x1 = -V{1}*cell[26];
auto x2 = -V{1}*cell[24];
auto x3 = -V{1}*cell[20];
auto x199 = cell[19] + x3;
auto x4 = -V{1}*cell[13];
auto x205 = cell[10] + cell[23] + x1 + x4;
auto x5 = -V{1}*cell[11];
auto x206 = x2 + x5;
auto x6 = -V{1}*cell[7];
auto x7 = -V{1}*cell[6];
auto x8 = -V{1}*cell[8];
auto x202 = cell[23] + x1 + x8;
auto x9 = -V{1}*cell[10];
auto x10 = -V{1}*cell[12];
auto x201 = cell[25] + x10 + x7 + x9;
auto x11 = -V{1}*cell[22];
auto x200 = cell[13] + cell[21] + x11;
auto x12 = cell[0] + x183 + x188 + x193 + x197 + V{1};
auto x198 = V{1} / (x12);
auto x13 = V{0.5}*x179;
auto x14 = x13 + x198*(cell[11] - V{1}*cell[3] + x192 + x199 + x2 + x200 + x201 + x202);
auto x212 = ((x14)*(x14));
auto x15 = -x13 + x14;
auto x16 = -V{1}*cell[4];
auto x204 = cell[24] + x16 + x5;
auto x17 = -V{1}*cell[5];
auto x203 = cell[17] + x17;
auto x18 = V{0.5}*x177;
auto x19 = x18 + x198*(-V{1}*cell[1] + x187 + x189 + x201 + x203 + x204 + x4 + x6);
auto x213 = ((x19)*(x19));
auto x20 = -x18 + x19;
auto x21 = x12*x20;
auto x22 = V{0.5}*x12;
auto x208 = V{2}*x22;
auto x23 = -x15*x21 + x183 + x199 + x205 + x206 + x22*(x15*x177 + x179*x20) + x6;
auto x24 = ((x23)*(x23));
auto x25 = -V{1}*cell[25];
auto x207 = x10 + x25;
auto x26 = -V{1}*cell[18];
auto x27 = -V{1}*cell[9];
auto x28 = V{0.5}*x178;
auto x29 = x198*(cell[12] - V{1}*cell[2] + x196 + x202 + x204 + x25 + x26 + x27 + x9) + x28;
auto x214 = ((x29)*(x29));
auto x30 = -x28 + x29;
auto x31 = x185 + x203 + x205 + x207 - x21*x30 + x22*(x177*x30 + x178*x20) + x26;
auto x32 = V{3};
auto x211 = ((x32)*(x32));
auto x215 = x19*x32;
auto x216 = x29*x32;
auto x217 = x14*x32;
auto x33 = V{1}*(x12 + V{-1})/x32;
auto x209 = -x33;
auto x210 = x183 + x209;
auto x34 = x177*x20;
auto x35 = cell[7] + x188 + x190 + x194 - x20*x21 + x208*x34 + x210;
auto x36 = x15*x179;
auto x37 = cell[13] - x12*((x15)*(x15)) + x184 + x186 + x193 + x195 + x208*x36 + x210;
auto x38 = x12*x30;
auto x39 = -x15*x38 + x191 + x200 + x206 + x207 + x22*(x15*x178 + x179*x30) + x27;
auto x40 = x178*x30;
auto x41 = cell[18] + cell[9] + x182 + x185 + x191 + x197 + x208*x40 + x209 - x30*x38;
auto x42 = V{2}*V{1} / (x0 + util::sqrt(V{2.828428}*x198*x211*((x180)*(x180))*util::sqrt(V{2}*x24 + V{2}*((x31)*(x31)) + ((x35)*(x35)) + ((x37)*(x37)) + V{2}*((x39)*(x39)) + ((x41)*(x41))) + ((x0)*(x0))));
auto x43 = V{1} - x42;
auto x44 = x212 + x213 + x214;
auto x45 = V{0.5}*x32*x44;
auto x221 = V{1} - x45;
auto x46 = V{0.296296296296296};
auto x47 = -x177*x215;
auto x48 = -x178*x216;
auto x49 = x47 + x48;
auto x50 = -x179*x217;
auto x220 = x48 + x50;
auto x223 = x47 + x50;
auto x51 = -x12*(V{0.5}*x42 + V{-1});
auto x52 = cell[0]*x43 - x42*x46*(x12*(x45 + V{-1}) + V{1}) + x46*x51*(x49 + x50);
auto x53 = -V{1}*x19;
auto x54 = x211;
auto x219 = V{1}*x54;
auto x55 = V{0.5}*x54;
auto x56 = V{0.0740740740740741};
auto x57 = x12*x56;
auto x58 = -x32*(x19 + V{1});
auto x218 = -x58;
auto x59 = x51*x56;
auto x60 = cell[1]*x43 - x42*(x56 - x57*(x221 + x32*x53 + x55*((x53)*(x53)))) + x59*(-x177*(x218 + x219*x53) + x220);
auto x61 = -V{1}*x29;
auto x62 = -x32*(x29 + V{1});
auto x222 = -x62;
auto x63 = cell[2]*x43 - x42*(x56 - x57*(x221 + x32*x61 + x55*((x61)*(x61)))) + x59*(-x178*(x219*x61 + x222) + x223);
auto x64 = -V{1}*x14;
auto x65 = -x32*(x14 + V{1});
auto x224 = -x65;
auto x66 = -(-cell[3]*x43 + x42*(x56 - x57*(x221 + x32*x64 + x55*((x64)*(x64)))) + x59*(x179*(x219*x64 + x224) - x49));
auto x67 = x53 + x61;
auto x68 = V{0.0185185185185185};
auto x69 = x12*x68;
auto x70 = -x219*x67;
auto x71 = x51*x68;
auto x72 = cell[4]*x43 - x42*(x68 - x69*(x221 + x32*x67 + x55*((x67)*(x67)))) + x71*(x177*(x58 + x70) + x178*(x62 + x70) + x50);
auto x73 = x29 + x53;
auto x74 = x54*x73;
auto x75 = -x32*(x29 + V{-1});
auto x76 = cell[5]*x43 - x42*(x68 - x69*(x221 + x32*x73 + x55*((x73)*(x73)))) + x71*(x177*(x58 - V{1}*x74) + x178*(x74 + x75) + x50);
auto x77 = x53 + x64;
auto x78 = -x219*x77;
auto x79 = cell[6]*x43 - x42*(x68 - x69*(x221 + x32*x77 + x55*((x77)*(x77)))) + x71*(x177*(x58 + x78) + x179*(x65 + x78) + x48);
auto x80 = x14 + x53;
auto x81 = x54*x80;
auto x82 = -x32*(x14 + V{-1});
auto x83 = cell[7]*x43 - x42*(x68 - x69*(x221 + x32*x80 + x55*((x80)*(x80)))) + x71*(x177*(x58 - V{1}*x81) + x179*(x81 + x82) + x48);
auto x84 = x61 + x64;
auto x85 = -x219*x84;
auto x86 = cell[8]*x43 - x42*(x68 - x69*(x221 + x32*x84 + x55*((x84)*(x84)))) + x71*(x178*(x62 + x85) + x179*(x65 + x85) + x47);
auto x87 = x14 + x61;
auto x88 = x54*x87;
auto x89 = cell[9]*x43 - x42*(x68 - x69*(x221 + x32*x87 + x55*((x87)*(x87)))) + x71*(x178*(x62 - V{1}*x88) + x179*(x82 + x88) + x47);
auto x90 = x64 + x67;
auto x91 = V{0.00462962962962963};
auto x92 = x12*x91;
auto x93 = -x219*x90;
auto x94 = x51*x91;
auto x95 = cell[10]*x43 - x42*(x91 - x92*(x221 + x32*x90 + x55*((x90)*(x90)))) + x94*(x177*(x58 + x93) + x178*(x62 + x93) + x179*(x65 + x93));
auto x96 = x14 + x67;
auto x97 = x54*x96;
auto x98 = -V{1}*x97;
auto x99 = cell[11]*x43 - x42*(x91 - x92*(x221 + x32*x96 + x55*((x96)*(x96)))) + x94*(x177*(x58 + x98) + x178*(x62 + x98) + x179*(x82 + x97));
auto x100 = x64 + x73;
auto x101 = x100*x54;
auto x102 = -V{1}*x101;
auto x103 = cell[12]*x43 - x42*(x91 - x92*(x100*x32 + x221 + x55*((x100)*(x100)))) + x94*(x177*(x102 + x58) + x178*(x101 + x75) + x179*(x102 + x65));
auto x104 = x14 + x73;
auto x105 = x104*x54;
auto x106 = cell[13]*x43 - x42*(x91 - x92*(x104*x32 + x221 + x55*((x104)*(x104)))) + x94*(-x177*(V{1}*x105 + x218) + x178*(x105 + x75) + x179*(x105 + x82));
auto x107 = -x32*(x19 + V{-1});
auto x108 = cell[14]*x43 - x42*(x56 - x57*(x213*x55 + x215 + x221)) + x59*(x177*(x107 + x19*x54) + x220);
auto x109 = cell[15]*x43 - x42*(x56 - x57*(x214*x55 + x216 + x221)) + x59*(x178*(x29*x54 + x75) + x223);
auto x110 = cell[16]*x43 - x42*(x56 - x57*(x212*x55 + x217 + x221)) + x59*(x179*(x14*x54 + x82) + x49);
auto x111 = x19 + x29;
auto x112 = x111*x54;
auto x113 = cell[17]*x43 - x42*(x68 - x69*(x111*x32 + x221 + x55*((x111)*(x111)))) + x71*(x177*(x107 + x112) + x178*(x112 + x75) + x50);
auto x114 = x19 + x61;
auto x115 = x114*x54;
auto x116 = cell[18]*x43 - x42*(x68 - x69*(x114*x32 + x221 + x55*((x114)*(x114)))) + x71*(x177*(x107 + x115) - x178*(V{1}*x115 + x222) + x50);
auto x117 = x14 + x19;
auto x118 = x117*x54;
auto x119 = cell[19]*x43 - x42*(x68 - x69*(x117*x32 + x221 + x55*((x117)*(x117)))) + x71*(x177*(x107 + x118) + x179*(x118 + x82) + x48);
auto x120 = x19 + x64;
auto x121 = x120*x54;
auto x122 = cell[20]*x43 - x42*(x68 - x69*(x120*x32 + x221 + x55*((x120)*(x120)))) + x71*(x177*(x107 + x121) - x179*(V{1}*x121 + x224) + x48);
auto x123 = x14 + x29;
auto x124 = x123*x54;
auto x125 = cell[21]*x43 - x42*(x68 - x69*(x123*x32 + x221 + x55*((x123)*(x123)))) + x71*(x178*(x124 + x75) + x179*(x124 + x82) + x47);
auto x126 = x29 + x64;
auto x127 = x126*x54;
auto x128 = cell[22]*x43 - x42*(x68 - x69*(x126*x32 + x221 + x55*((x126)*(x126)))) + x71*(x178*(x127 + x75) - x179*(V{1}*x127 + x224) + x47);
auto x226 = x103 + x109 + x128 + x76;
auto x129 = x111 + x14;
auto x130 = x129*x54;
auto x131 = cell[23]*x43 - x42*(x91 - x92*(x129*x32 + x221 + x55*((x129)*(x129)))) + x94*(x177*(x107 + x130) + x178*(x130 + x75) + x179*(x130 + x82));
auto x227 = x110 + x131 + x83 + x89 + x99;
auto x132 = x111 + x64;
auto x133 = x132*x54;
auto x134 = cell[24]*x43 - x42*(x91 - x92*(x132*x32 + x221 + x55*((x132)*(x132)))) + x94*(x177*(x107 + x133) + x178*(x133 + x75) - x179*(V{1}*x133 + x224));
auto x135 = x114 + x14;
auto x136 = x135*x54;
auto x137 = cell[25]*x43 - x42*(x91 - x92*(x135*x32 + x221 + x55*((x135)*(x135)))) + x94*(x177*(x107 + x136) - x178*(V{1}*x136 + x222) + x179*(x136 + x82));
auto x138 = x114 + x64;
auto x139 = x138*x54;
auto x140 = -V{1}*x139;
auto x141 = cell[26]*x43 - x42*(x91 - x92*(x138*x32 + x221 + x55*((x138)*(x138)))) + x94*(x177*(x107 + x139) + x178*(x140 + x62) + x179*(x140 + x65));
auto x225 = x108 + x116 + x122 + x141;
auto x142 = x106 + x113 + x119 + x125 + x134 + x137 + x225 + x226 + x227 + x52 + x60 + x63 + x66 + x72 + x79 + x86 + x95 + V{1};
auto x228 = V{1} / (x142);
auto x143 = -V{1}*x72;
auto x144 = -V{1}*x79;
auto x145 = -V{1}*x95;
auto x146 = -V{1}*x99;
auto x229 = x113 + x131 + x134 + x143 + x145 + x146;
auto x147 = -V{1}*x103;
auto x230 = x119 + x137 + x144 + x147;
auto x148 = -V{1}*x86;
auto x149 = -V{1}*x141;
auto x231 = x106 + x125 + x148 + x149;
cell[0] = x52;
cell[10] = x95;
cell[11] = x99;
cell[12] = x103;
cell[13] = x106;
cell[14] = x108;
cell[15] = x109;
cell[16] = x110;
cell[17] = x113;
cell[18] = x116;
cell[19] = x119;
cell[1] = x60;
cell[20] = x122;
cell[21] = x125;
cell[22] = x128;
cell[23] = x131;
cell[24] = x134;
cell[25] = x137;
cell[26] = x141;
cell[2] = x63;
cell[3] = x66;
cell[4] = x72;
cell[5] = x76;
cell[6] = x79;
cell[7] = x83;
cell[8] = x86;
cell[9] = x89;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x142;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x18 + x228*(-V{1}*x106 + x225 + x229 + x230 - V{1}*x60 - V{1}*x76 - V{1}*x83);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x228*(-V{1}*x116 - V{1}*x137 + x226 + x229 + x231 - V{1}*x63 - V{1}*x89) + x28;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x13 + x228*(-V{1}*x122 - V{1}*x128 - V{1}*x134 + x145 + x227 + x230 + x231 - V{1}*x66);
return { x12, x44 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
