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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB>, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x208 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x209 = parameters.template get<descriptors::OMEGA>();
auto x210 = cell[17] + cell[9];
auto x211 = cell[15] + cell[7];
auto x212 = cell[12] + x210 + x211;
auto x213 = cell[3] + cell[6] + x212;
auto x214 = cell[13] + cell[5];
auto x215 = cell[11] + cell[18] + x214;
auto x216 = cell[2] + cell[4] + x215;
auto x217 = cell[10] + cell[14] + cell[16];
auto x218 = cell[1] + x217;
auto x0 = V{3};
auto x226 = ((x0)*(x0));
auto x227 = V{1}/x0;
auto x1 = -V{1}*cell[4];
auto x2 = -V{1}*cell[5];
auto x221 = cell[13] + x2;
auto x3 = -V{1}*cell[6];
auto x4 = -V{1}*cell[7];
auto x220 = cell[15] + x4;
auto x5 = cell[0] + cell[8] + x213 + x216 + x218 + V{1};
auto x219 = V{1} / (x5);
auto x6 = x219*(-V{1}*cell[1] + x1 + x217 + x220 + x221 + x3);
auto x222 = ((x6)*(x6));
auto x7 = x222;
auto x8 = -V{1}*cell[8];
auto x9 = -V{1}*cell[9];
auto x223 = cell[17] + x9;
auto x10 = -V{1}*cell[14];
auto x11 = x219*(-V{1}*cell[2] + x1 + x10 + x215 + x223 + x8);
auto x224 = ((x11)*(x11));
auto x12 = x224;
auto x13 = -V{1}*cell[16];
auto x14 = -V{1}*cell[18];
auto x15 = x219*(-V{1}*cell[3] + x13 + x14 + x212 + x3 + x8);
auto x225 = ((x15)*(x15));
auto x16 = x12 + x225 + x7;
auto x17 = V{0.5}*x0*x16;
auto x235 = V{1} - x17;
auto x18 = V{0.333333333333333};
auto x19 = x226;
auto x20 = V{0.5}*x19;
auto x21 = x227;
auto x232 = V{1} - x21;
auto x22 = -x21;
auto x23 = x21*(x5 + V{-1});
auto x228 = -x23;
auto x229 = cell[8] + x228;
auto x24 = cell[16] + cell[18] + x213 - x225*x5 + x229;
auto x25 = x22*x24;
auto x26 = x11*x5;
auto x27 = cell[14] - x11*x26 + x210 + x216 + x229;
auto x28 = x22*x27;
auto x29 = x5*x6;
auto x30 = cell[4] + cell[6] + x211 + x214 + x218 + x228 - x29*x6;
auto x31 = x22*x30;
auto x32 = x28 + x31;
auto x33 = V{1}/x209;
auto x34 = cell[6] + x13 - x15*x29 + x220;
auto x254 = x19*x34;
auto x35 = ((x34)*(x34));
auto x36 = cell[4] + x10 - x11*x29 + x221;
auto x241 = x19*x36;
auto x37 = cell[8] + x14 - x15*x26 + x223;
auto x261 = x19*x37;
auto x38 = V{1} - V{2}*V{1} / (x33 + util::sqrt(V{2.828428}*x219*x226*((x208)*(x208))*util::sqrt(V{2}*x35 + ((x24)*(x24)) + ((x27)*(x27)) + ((x30)*(x30)) + V{2}*((x36)*(x36)) + V{2}*((x37)*(x37))) + ((x33)*(x33))));
auto x39 = x18*(x20*x38*(x25 + x32) - x5*(x17 + V{-1}) + V{-1});
auto x40 = x0*x19;
auto x41 = V{0.166667}*x40;
auto x266 = x22*x41;
auto x42 = -x227;
auto x43 = -x42;
auto x230 = x41*x43;
auto x44 = -x230;
auto x45 = x11*x7;
auto x46 = x224*x6;
auto x47 = x45 - x46;
auto x237 = x41*x47;
auto x48 = V{0.5}*x40;
auto x231 = x43*x48;
auto x267 = x22*x48;
auto x49 = x231;
auto x50 = x45 + x46;
auto x238 = x48*x50;
auto x51 = x44*x47 + x49*x50;
auto x52 = x225*x6;
auto x53 = x11*x225;
auto x54 = x52 - x53;
auto x250 = x41*x54;
auto x251 = x250*(x22 - x43);
auto x55 = x52 + x53;
auto x247 = x48*x55;
auto x56 = x49*x55;
auto x236 = x235 + x56;
auto x57 = -V{1}*x6;
auto x58 = V{0.0555555555555556};
auto x233 = -x58;
auto x234 = x38*x58;
auto x59 = x5*x58;
auto x60 = -V{1}*x22;
auto x61 = 0;
auto x62 = V{2}*x15;
auto x63 = x24*x6 + x34*x62;
auto x64 = V{2}*x11;
auto x65 = x27*x6 + x36*x64;
auto x66 = x63 - x65;
auto x67 = V{2}*x48*x60;
auto x68 = x63 + x65;
auto x69 = x232;
auto x248 = x41*(x22 - x69);
auto x271 = x250*x69;
auto x70 = x30*x69;
auto x71 = x28 + x70;
auto x72 = x20*(x25 + x71);
auto x73 = x233 + x234*(x61*x66 + x67*x68 + x72) + x59*(x0*x57 + x20*((x57)*(x57)) + x230*x54 + x236 + x51);
auto x74 = x60*(x237 + x238);
auto x75 = -V{1}*x11;
auto x76 = V{2}*x6;
auto x77 = x11*x30 + x36*x76;
auto x78 = x11*x24 + x37*x62;
auto x79 = x77 - x78;
auto x80 = x77 + x78;
auto x81 = x27*x69;
auto x82 = x31 + x81;
auto x83 = x20*(x25 + x82);
auto x84 = x233 + x234*(x61*x79 + x67*x80 + x83) + x59*(x0*x75 + x20*((x75)*(x75)) + x236 + x44*x54 + x74);
auto x85 = x12*x15;
auto x86 = x15*x7;
auto x87 = x85 - x86;
auto x88 = x85 + x86;
auto x89 = -V{1}*x15;
auto x90 = x15*x27 + x37*x64;
auto x91 = x15*x30 + x34*x76;
auto x92 = x90 - x91;
auto x93 = x90 + x91;
auto x94 = x24*x69;
auto x95 = x20*(x32 + x94);
auto x96 = x233 + x234*(x61*x92 + x67*x93 + x95) + x59*(x0*x89 + x20*((x89)*(x89)) + x235 + x61*x87 + x67*x88);
auto x97 = -V{1}*x69;
auto x239 = -x97;
auto x240 = x41*(x239 + x60);
auto x98 = -x42 + V{-1};
auto x246 = -x98;
auto x253 = x250*x98;
auto x99 = x57 + x75;
auto x100 = V{0.0277777777777778};
auto x242 = -x100;
auto x245 = x100*x38;
auto x101 = x100*x5;
auto x102 = x240;
auto x103 = x102*x66;
auto x104 = -x240;
auto x105 = x104*x79;
auto x106 = x48*(x60 + x97);
auto x107 = x106*x68;
auto x243 = x103 + x107;
auto x108 = x106*x80;
auto x244 = x105 + x108;
auto x109 = x20*(x25 + x70 + x81);
auto x110 = x109 + x241;
auto x111 = x101*(x0*x99 + x20*((x99)*(x99)) + V{2}*x231*x55 + x235 + x237*(x246 + x97) + x238*(x97 + x98)) + x242 + x245*(x110 + x243 + x244);
auto x112 = x247*(x22 + x43);
auto x252 = x112 + x235;
auto x113 = x11 + x57;
auto x114 = -x248;
auto x115 = x114*x79;
auto x116 = x48*(x22 + x69);
auto x117 = x116*x80;
auto x249 = x115 + x117;
auto x118 = x109 - V{1}*x241;
auto x119 = x101*(x0*x113 + x20*((x113)*(x113)) + x237*(x246 + x69) + x238*(x69 + x98) - x251 + x252) + x242 + x245*(x118 + x243 + x249);
auto x120 = x102*x87;
auto x121 = x253;
auto x258 = x121 + x51;
auto x122 = x106*x88;
auto x123 = x247*x98;
auto x124 = x122 + x123;
auto x257 = x124 + x235;
auto x125 = x57 + x89;
auto x126 = x102*x92;
auto x127 = x104*x66;
auto x255 = x107 + x127;
auto x128 = x106*x93;
auto x256 = x126 + x128;
auto x129 = x20*(x71 + x94);
auto x130 = x129 + x254;
auto x131 = x101*(x0*x125 + x120 + x20*((x125)*(x125)) + x257 + x258) + x242 + x245*(x130 + x255 + x256);
auto x132 = x248;
auto x133 = x132*x87;
auto x134 = x116*x88;
auto x135 = x123 + x134;
auto x260 = x135 + x235;
auto x136 = x15 + x57;
auto x137 = x132*x92;
auto x138 = x116*x93;
auto x259 = x137 + x138;
auto x139 = x129 - V{1}*x254;
auto x140 = x101*(x0*x136 + x133 + x20*((x136)*(x136)) + x258 + x260) + x242 + x245*(x139 + x255 + x259);
auto x141 = x104*x87;
auto x142 = -x253;
auto x264 = x142 + x74;
auto x143 = x75 + x89;
auto x144 = x104*x92;
auto x263 = x128 + x144;
auto x145 = x102*x79;
auto x262 = x108 + x145;
auto x146 = x20*(x82 + x94);
auto x147 = x146 + x261;
auto x148 = x101*(x0*x143 + x141 + x20*((x143)*(x143)) + x257 + x264) + x242 + x245*(x147 + x262 + x263);
auto x149 = x114*x87;
auto x150 = x15 + x75;
auto x151 = x114*x92;
auto x265 = x138 + x151;
auto x152 = x146 - V{1}*x261;
auto x153 = x101*(x0*x150 + x149 + x20*((x150)*(x150)) + x260 + x264) + x242 + x245*(x152 + x262 + x265);
auto x154 = -x266;
auto x155 = x267;
auto x156 = x155*x50;
auto x157 = x154*x47 + x156;
auto x272 = x157 + x235;
auto x158 = x266;
auto x159 = x155*x55;
auto x268 = x159 + x235;
auto x160 = 0;
auto x161 = V{2}*x267;
auto x162 = x233 + x234*(x160*x66 + x161*x68 + x72) + x59*(x0*x6 + x157 + x158*x54 + x20*x222 + x268);
auto x163 = x156 + x158*x47;
auto x275 = x163 + x235;
auto x164 = x233 + x234*(x160*x79 + x161*x80 + x83) + x59*(x0*x11 + x154*x54 + x163 + x20*x224 + x268);
auto x165 = x233 + x234*(x160*x92 + x161*x93 + x95) + x59*(x0*x15 + x160*x87 + x161*x88 + x20*x225 + x235);
auto x278 = x140 + x153 + x165;
auto x166 = x11 + x6;
auto x167 = x132*x66;
auto x168 = x116*x68;
auto x269 = x167 + x168;
auto x169 = x101*(x0*x166 + x160*x54 + x161*x55 + x20*((x166)*(x166)) + x235 + V{2}*x238*x69) + x242 + x245*(x110 + x249 + x269);
auto x170 = x232;
auto x273 = x170*x250;
auto x171 = x6 + x75;
auto x172 = x101*(x0*x171 + x20*((x171)*(x171)) - x237*(x170 + x239) + x238*(x170 + x97) + x251 + x252) + x242 + x245*(x118 + x244 + x269);
auto x173 = x134 + x247*x69;
auto x174 = x15 + x6;
auto x175 = x114*x66;
auto x270 = x168 + x175;
auto x176 = x101*(x0*x174 + x133 + x173 + x20*((x174)*(x174)) + x271 + x272) + x242 + x245*(x130 + x259 + x270);
auto x177 = x122 + x170*x247;
auto x178 = x6 + x89;
auto x179 = x101*(x0*x178 + x120 + x177 + x20*((x178)*(x178)) + x272 + x273) + x242 + x245*(x139 + x256 + x270);
auto x276 = x162 + x172 + x179;
auto x180 = x11 + x15;
auto x181 = x132*x79;
auto x274 = x117 + x181;
auto x182 = x101*(x0*x180 + x149 + x173 + x20*((x180)*(x180)) - x271 + x275) + x242 + x245*(x147 + x265 + x274);
auto x183 = x11 + x89;
auto x184 = x101*(x0*x183 + x141 + x177 + x20*((x183)*(x183)) - x273 + x275) + x242 + x245*(x152 + x263 + x274);
auto x277 = x119 + x164 + x184;
auto x185 = x111 + x131 + x148 + x169 + x176 + x182 + x276 + x277 + x278 + x39 + x73 + x84 + x96 + V{1};
auto x279 = V{1} / (x185);
auto x186 = -V{1}*x111;
auto x280 = x169 + x186;
auto x187 = -V{1}*x131;
auto x281 = x176 + x187;
auto x188 = -V{1}*x148;
auto x282 = x182 + x188;
cell[0] = x39;
cell[10] = x162;
cell[11] = x164;
cell[12] = x165;
cell[13] = x169;
cell[14] = x172;
cell[15] = x176;
cell[16] = x179;
cell[17] = x182;
cell[18] = x184;
cell[1] = x73;
cell[2] = x84;
cell[3] = x96;
cell[4] = x111;
cell[5] = x119;
cell[6] = x131;
cell[7] = x140;
cell[8] = x148;
cell[9] = x153;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x185;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x279*(-V{1}*x119 - V{1}*x140 + x276 + x280 + x281 - V{1}*x73);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x279*(-V{1}*x153 - V{1}*x172 + x277 + x280 + x282 - V{1}*x84);
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x279*(-V{1}*x179 - V{1}*x184 + x278 + x281 + x282 - V{1}*x96);
return { x5, x16 };
}
};

}

}

// Generation Info: commit=61323934e549f5f16850d6c6280b91bc057cb6c3
