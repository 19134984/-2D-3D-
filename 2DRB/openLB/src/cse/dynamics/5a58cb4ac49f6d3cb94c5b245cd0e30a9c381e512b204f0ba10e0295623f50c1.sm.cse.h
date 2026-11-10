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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::ParameterFromCell<collision::LES::SMAGORINSKY, collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB> >, forcing::Guo<momenta::ForcedWithStress> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x244 = cell.template getFieldComponent<olb::collision::LES::SMAGORINSKY>(0);
auto x245 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x246 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x247 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x249 = parameters.template get<descriptors::OMEGA>();
auto x248 = cell[17] + cell[9];
auto x250 = cell[15] + cell[7];
auto x251 = cell[12] + x248 + x250;
auto x252 = cell[3] + cell[6] + x251;
auto x253 = cell[13] + cell[5];
auto x254 = cell[11] + cell[18] + x253;
auto x255 = cell[2] + cell[4] + x254;
auto x256 = cell[10] + cell[14] + cell[16];
auto x257 = cell[1] + x256;
auto x0 = V{3};
auto x265 = ((x0)*(x0));
auto x266 = V{1}/x0;
auto x1 = -V{1}*cell[4];
auto x2 = -V{1}*cell[5];
auto x260 = cell[13] + x2;
auto x3 = -V{1}*cell[6];
auto x4 = -V{1}*cell[7];
auto x259 = cell[15] + x4;
auto x5 = cell[0] + cell[8] + x252 + x255 + x257 + V{1};
auto x258 = V{1} / (x5);
auto x6 = V{0.5}*x245;
auto x7 = x258*(-V{1}*cell[1] + x1 + x256 + x259 + x260 + x3) + x6;
auto x261 = ((x7)*(x7));
auto x270 = x0*x7;
auto x8 = x261;
auto x9 = -V{1}*cell[8];
auto x10 = -V{1}*cell[9];
auto x262 = cell[17] + x10;
auto x11 = -V{1}*cell[14];
auto x12 = V{0.5}*x246;
auto x13 = x12 + x258*(-V{1}*cell[2] + x1 + x11 + x254 + x262 + x9);
auto x263 = ((x13)*(x13));
auto x271 = x0*x13;
auto x14 = x263;
auto x15 = -V{1}*cell[16];
auto x16 = -V{1}*cell[18];
auto x17 = V{0.5}*x247;
auto x18 = x17 + x258*(-V{1}*cell[3] + x15 + x16 + x251 + x3 + x9);
auto x264 = ((x18)*(x18));
auto x272 = x0*x18;
auto x19 = x14 + x264 + x8;
auto x20 = V{0.5}*x0*x19;
auto x281 = V{1} - x20;
auto x21 = V{0.333333333333333};
auto x22 = x265;
auto x279 = V{1}*x22;
auto x23 = V{0.5}*x22;
auto x24 = x266;
auto x275 = V{1} - x24;
auto x25 = -x24;
auto x26 = -x17 + x18;
auto x27 = x24*(x5 + V{-1});
auto x268 = -x27;
auto x269 = cell[8] + x268;
auto x28 = V{0.5}*x5;
auto x267 = V{2}*x28;
auto x29 = x247*x26;
auto x30 = cell[16] + cell[18] + x252 + x267*x29 + x269 - x5*((x26)*(x26));
auto x31 = x25*x30;
auto x32 = -x12 + x13;
auto x33 = x32*x5;
auto x34 = x246*x32;
auto x35 = cell[14] + x248 + x255 + x267*x34 + x269 - x32*x33;
auto x36 = x25*x35;
auto x37 = -x6 + x7;
auto x38 = x37*x5;
auto x39 = x245*x37;
auto x40 = cell[4] + cell[6] + x250 + x253 + x257 + x267*x39 + x268 - x37*x38;
auto x41 = x25*x40;
auto x42 = x36 + x41;
auto x43 = V{1}/x249;
auto x44 = cell[6] + x15 + x259 - x26*x38 + x28*(x245*x26 + x247*x37);
auto x300 = x22*x44;
auto x45 = ((x44)*(x44));
auto x46 = cell[4] + x11 + x260 + x28*(x245*x32 + x246*x37) - x32*x38;
auto x287 = x22*x46;
auto x47 = cell[8] + x16 - x26*x33 + x262 + x28*(x246*x26 + x247*x32);
auto x307 = x22*x47;
auto x48 = V{2}*V{1} / (x43 + util::sqrt(V{2.828428}*x258*x265*((x244)*(x244))*util::sqrt(V{2}*x45 + ((x30)*(x30)) + ((x35)*(x35)) + ((x40)*(x40)) + V{2}*((x46)*(x46)) + V{2}*((x47)*(x47))) + ((x43)*(x43))));
auto x49 = V{1} - x48;
auto x50 = -x245*x270;
auto x51 = -x246*x271;
auto x52 = x50 + x51;
auto x53 = -x247*x272;
auto x280 = x51 + x53;
auto x285 = x50 + x53;
auto x54 = -x5*(V{0.5}*x48 + V{-1});
auto x55 = x21*(x23*x49*(x31 + x42) - x5*(x20 + V{-1}) + x54*(x52 + x53) + V{-1});
auto x56 = x0*x22;
auto x57 = V{0.166667}*x56;
auto x312 = x25*x57;
auto x58 = -x266;
auto x59 = -x58;
auto x273 = x57*x59;
auto x60 = -x273;
auto x61 = x13*x8;
auto x62 = x263*x7;
auto x63 = x61 - x62;
auto x283 = x57*x63;
auto x64 = V{0.5}*x56;
auto x274 = x59*x64;
auto x313 = x25*x64;
auto x65 = x274;
auto x66 = x61 + x62;
auto x284 = x64*x66;
auto x67 = x60*x63 + x65*x66;
auto x68 = x264*x7;
auto x69 = x13*x264;
auto x70 = x68 - x69;
auto x296 = x57*x70;
auto x297 = x296*(x25 - x59);
auto x71 = x68 + x69;
auto x292 = x64*x71;
auto x72 = x65*x71;
auto x282 = x281 + x72;
auto x73 = -V{1}*x7;
auto x74 = V{0.0555555555555556};
auto x276 = -x74;
auto x277 = x49*x74;
auto x75 = x5*x74;
auto x76 = -V{1}*x25;
auto x77 = 0;
auto x78 = V{2}*x18;
auto x79 = x30*x7 + x44*x78;
auto x80 = V{2}*x13;
auto x81 = x35*x7 + x46*x80;
auto x82 = x79 - x81;
auto x83 = V{2}*x64*x76;
auto x84 = x79 + x81;
auto x85 = x275;
auto x293 = -x85;
auto x294 = x57*(x25 + x293);
auto x317 = x296*x85;
auto x86 = x40*x85;
auto x87 = x36 + x86;
auto x88 = x23*(x31 + x87);
auto x89 = -x0*(x7 + V{1});
auto x278 = -x89;
auto x90 = x54*x74;
auto x91 = x276 + x277*(x77*x82 + x83*x84 + x88) + x75*(x0*x73 + x23*((x73)*(x73)) + x273*x70 + x282 + x67) + x90*(-x245*(x278 + x279*x73) + x280);
auto x92 = x76*(x283 + x284);
auto x93 = -V{1}*x13;
auto x94 = V{2}*x7;
auto x95 = x13*x40 + x46*x94;
auto x96 = x13*x30 + x47*x78;
auto x97 = x95 - x96;
auto x98 = x95 + x96;
auto x99 = x35*x85;
auto x100 = x41 + x99;
auto x101 = x23*(x100 + x31);
auto x102 = -x0*(x13 + V{1});
auto x103 = x276 + x277*(x101 + x77*x97 + x83*x98) + x75*(x0*x93 + x23*((x93)*(x93)) + x282 + x60*x70 + x92) + x90*(x246*(x102 - x279*x93) + x285);
auto x104 = x14*x18;
auto x105 = x18*x8;
auto x106 = x104 - x105;
auto x107 = x104 + x105;
auto x108 = -V{1}*x18;
auto x109 = x18*x35 + x47*x80;
auto x110 = x18*x40 + x44*x94;
auto x111 = x109 - x110;
auto x112 = x109 + x110;
auto x113 = x30*x85;
auto x114 = x23*(x113 + x42);
auto x115 = -x0*(x18 + V{1});
auto x116 = x276 + x277*(x111*x77 + x112*x83 + x114) + x75*(x0*x108 + x106*x77 + x107*x83 + x23*((x108)*(x108)) + x281) - x90*(x247*(x108*x279 - x115) - x52);
auto x117 = -V{1}*x85;
auto x286 = x57*(x117 - x76);
auto x118 = -x58 + V{-1};
auto x299 = x118*x296;
auto x119 = x73 + x93;
auto x120 = V{0.0277777777777778};
auto x288 = -x120;
auto x291 = x120*x49;
auto x121 = x120*x5;
auto x122 = -x286;
auto x123 = x122*x82;
auto x124 = x286;
auto x125 = x124*x97;
auto x126 = x64*(x117 + x76);
auto x127 = x126*x84;
auto x289 = x123 + x127;
auto x128 = x126*x98;
auto x290 = x125 + x128;
auto x129 = x23*(x31 + x86 + x99);
auto x130 = x129 + x287;
auto x131 = -x119*x279;
auto x132 = x120*x54;
auto x133 = x121*(x0*x119 + x23*((x119)*(x119)) + V{2}*x274*x71 + x281 + x283*(x117 - x118) + x284*(x117 + x118)) + x132*(x245*(x131 + x89) + x246*(x102 + x131) + x53) + x288 + x291*(x130 + x289 + x290);
auto x134 = x292*(x25 + x59);
auto x298 = x134 + x281;
auto x135 = x13 + x73;
auto x136 = -x294;
auto x137 = x136*x97;
auto x138 = x64*(x25 + x85);
auto x139 = x138*x98;
auto x295 = x137 + x139;
auto x140 = x129 - V{1}*x287;
auto x141 = x135*x22;
auto x142 = -x0*(x13 + V{-1});
auto x143 = x121*(x0*x135 + x23*((x135)*(x135)) - x283*(x118 + x293) + x284*(x118 + x85) - x297 + x298) + x132*(-x245*(V{1}*x141 + x278) + x246*(x141 + x142) + x53) + x288 + x291*(x140 + x289 + x295);
auto x144 = x106*x122;
auto x145 = x299;
auto x304 = x145 + x67;
auto x146 = x107*x126;
auto x147 = x118*x292;
auto x148 = x146 + x147;
auto x303 = x148 + x281;
auto x149 = x108 + x73;
auto x150 = x111*x122;
auto x151 = x124*x82;
auto x301 = x127 + x151;
auto x152 = x112*x126;
auto x302 = x150 + x152;
auto x153 = x23*(x113 + x87);
auto x154 = x153 + x300;
auto x155 = -x149*x279;
auto x156 = x121*(x0*x149 + x144 + x23*((x149)*(x149)) + x303 + x304) + x132*(x245*(x155 + x89) + x247*(x115 + x155) + x51) + x288 + x291*(x154 + x301 + x302);
auto x157 = x294;
auto x158 = x106*x157;
auto x159 = x107*x138;
auto x160 = x147 + x159;
auto x306 = x160 + x281;
auto x161 = x18 + x73;
auto x162 = x111*x157;
auto x163 = x112*x138;
auto x305 = x162 + x163;
auto x164 = x153 - V{1}*x300;
auto x165 = x161*x22;
auto x166 = -x0*(x18 + V{-1});
auto x167 = x121*(x0*x161 + x158 + x23*((x161)*(x161)) + x304 + x306) + x132*(-x245*(V{1}*x165 + x278) + x247*(x165 + x166) + x51) + x288 + x291*(x164 + x301 + x305);
auto x168 = x106*x124;
auto x169 = -x299;
auto x310 = x169 + x92;
auto x170 = x108 + x93;
auto x171 = x111*x124;
auto x309 = x152 + x171;
auto x172 = x122*x97;
auto x308 = x128 + x172;
auto x173 = x23*(x100 + x113);
auto x174 = x173 + x307;
auto x175 = -x170*x279;
auto x176 = x121*(x0*x170 + x168 + x23*((x170)*(x170)) + x303 + x310) + x132*(x246*(x102 + x175) + x247*(x115 + x175) + x50) + x288 + x291*(x174 + x308 + x309);
auto x177 = x106*x136;
auto x178 = x18 + x93;
auto x179 = x111*x136;
auto x311 = x163 + x179;
auto x180 = x173 - V{1}*x307;
auto x181 = x178*x22;
auto x182 = x121*(x0*x178 + x177 + x23*((x178)*(x178)) + x306 + x310) + x132*(x246*(x102 - V{1}*x181) + x247*(x166 + x181) + x50) + x288 + x291*(x180 + x308 + x311);
auto x183 = -x312;
auto x184 = x313;
auto x185 = x184*x66;
auto x186 = x183*x63 + x185;
auto x318 = x186 + x281;
auto x187 = x312;
auto x188 = x184*x71;
auto x314 = x188 + x281;
auto x189 = 0;
auto x190 = V{2}*x313;
auto x191 = -x0*(x7 + V{-1});
auto x192 = x276 + x277*(x189*x82 + x190*x84 + x88) + x75*(x186 + x187*x70 + x23*x261 + x270 + x314) + x90*(x245*(x191 + x22*x7) + x280);
auto x193 = x185 + x187*x63;
auto x321 = x193 + x281;
auto x194 = x276 + x277*(x101 + x189*x97 + x190*x98) + x75*(x183*x70 + x193 + x23*x263 + x271 + x314) + x90*(x246*(x13*x22 + x142) + x285);
auto x195 = x276 + x277*(x111*x189 + x112*x190 + x114) + x75*(x106*x189 + x107*x190 + x23*x264 + x272 + x281) + x90*(x247*(x166 + x18*x22) + x52);
auto x324 = x167 + x182 + x195;
auto x196 = x13 + x7;
auto x197 = x157*x82;
auto x198 = x138*x84;
auto x315 = x197 + x198;
auto x199 = x196*x22;
auto x200 = x121*(x0*x196 + x189*x70 + x190*x71 + x23*((x196)*(x196)) + x281 + V{2}*x284*x85) + x132*(x245*(x191 + x199) + x246*(x142 + x199) + x53) + x288 + x291*(x130 + x295 + x315);
auto x201 = x275;
auto x319 = x201*x296;
auto x202 = x7 + x93;
auto x203 = x202*x22;
auto x204 = x121*(x0*x202 + x23*((x202)*(x202)) + x283*(x117 - x201) + x284*(x117 + x201) + x297 + x298) + x132*(x245*(x191 + x203) + x246*(x102 - V{1}*x203) + x53) + x288 + x291*(x140 + x290 + x315);
auto x205 = x159 + x292*x85;
auto x206 = x18 + x7;
auto x207 = x136*x82;
auto x316 = x198 + x207;
auto x208 = x206*x22;
auto x209 = x121*(x0*x206 + x158 + x205 + x23*((x206)*(x206)) + x317 + x318) + x132*(x245*(x191 + x208) + x247*(x166 + x208) + x51) + x288 + x291*(x154 + x305 + x316);
auto x210 = x146 + x201*x292;
auto x211 = x108 + x7;
auto x212 = x211*x22;
auto x213 = x121*(x0*x211 + x144 + x210 + x23*((x211)*(x211)) + x318 + x319) + x132*(x245*(x191 + x212) + x247*(x115 - V{1}*x212) + x51) + x288 + x291*(x164 + x302 + x316);
auto x322 = x192 + x204 + x213;
auto x214 = x13 + x18;
auto x215 = x157*x97;
auto x320 = x139 + x215;
auto x216 = x214*x22;
auto x217 = x121*(x0*x214 + x177 + x205 + x23*((x214)*(x214)) - x317 + x321) + x132*(x246*(x142 + x216) + x247*(x166 + x216) + x50) + x288 + x291*(x174 + x311 + x320);
auto x218 = x108 + x13;
auto x219 = x218*x22;
auto x220 = x121*(x0*x218 + x168 + x210 + x23*((x218)*(x218)) - x319 + x321) + x132*(x246*(x142 + x219) + x247*(x115 - V{1}*x219) + x50) + x288 + x291*(x180 + x309 + x320);
auto x323 = x143 + x194 + x220;
auto x221 = x103 + x116 + x133 + x156 + x176 + x200 + x209 + x217 + x322 + x323 + x324 + x55 + x91 + V{1};
auto x325 = V{1} / (x221);
auto x222 = -V{1}*x133;
auto x326 = x200 + x222;
auto x223 = -V{1}*x156;
auto x327 = x209 + x223;
auto x224 = -V{1}*x176;
auto x328 = x217 + x224;
cell[0] = x55;
cell[10] = x192;
cell[11] = x194;
cell[12] = x195;
cell[13] = x200;
cell[14] = x204;
cell[15] = x209;
cell[16] = x213;
cell[17] = x217;
cell[18] = x220;
cell[1] = x91;
cell[2] = x103;
cell[3] = x116;
cell[4] = x133;
cell[5] = x143;
cell[6] = x156;
cell[7] = x167;
cell[8] = x176;
cell[9] = x182;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x221;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x325*(-V{1}*x143 - V{1}*x167 + x322 + x326 + x327 - V{1}*x91) + x6;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x12 + x325*(-V{1}*x103 - V{1}*x182 - V{1}*x204 + x323 + x326 + x328);
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x17 + x325*(-V{1}*x116 - V{1}*x213 - V{1}*x220 + x324 + x327 + x328);
return { x5, x19 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
