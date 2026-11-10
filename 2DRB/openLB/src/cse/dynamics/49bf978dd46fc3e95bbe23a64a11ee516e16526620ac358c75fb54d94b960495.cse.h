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
struct CSE<Dual<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, dynamics::DefaultCombination> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x254 = cell.template getFieldComponent<olb::opti::DJDF>(11);
auto x251 = cell.template getFieldComponent<olb::opti::DJDF>(8);
auto x268 = cell.template getFieldComponent<olb::opti::F>(6);
auto x280 = cell.template getFieldComponent<olb::opti::F>(18);
auto x263 = cell.template getFieldComponent<olb::opti::F>(1);
auto x269 = cell.template getFieldComponent<olb::opti::F>(7);
auto x258 = cell.template getFieldComponent<olb::opti::DJDF>(15);
auto x264 = cell.template getFieldComponent<olb::opti::F>(2);
auto x260 = cell.template getFieldComponent<olb::opti::DJDF>(17);
auto x266 = cell.template getFieldComponent<olb::opti::F>(4);
auto x273 = cell.template getFieldComponent<olb::opti::F>(11);
auto x247 = cell.template getFieldComponent<olb::opti::DJDF>(4);
auto x255 = cell.template getFieldComponent<olb::opti::DJDF>(12);
auto x243 = cell.template getFieldComponent<olb::opti::DJDF>(0);
auto x250 = cell.template getFieldComponent<olb::opti::DJDF>(7);
auto x267 = cell.template getFieldComponent<olb::opti::F>(5);
auto x244 = cell.template getFieldComponent<olb::opti::DJDF>(1);
auto x282 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x257 = cell.template getFieldComponent<olb::opti::DJDF>(14);
auto x270 = cell.template getFieldComponent<olb::opti::F>(8);
auto x281 = parameters.template get<descriptors::OMEGA>();
auto x278 = cell.template getFieldComponent<olb::opti::F>(16);
auto x276 = cell.template getFieldComponent<olb::opti::F>(14);
auto x249 = cell.template getFieldComponent<olb::opti::DJDF>(6);
auto x242 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x259 = cell.template getFieldComponent<olb::opti::DJDF>(16);
auto x274 = cell.template getFieldComponent<olb::opti::F>(12);
auto x272 = cell.template getFieldComponent<olb::opti::F>(10);
auto x275 = cell.template getFieldComponent<olb::opti::F>(13);
auto x261 = cell.template getFieldComponent<olb::opti::DJDF>(18);
auto x253 = cell.template getFieldComponent<olb::opti::DJDF>(10);
auto x271 = cell.template getFieldComponent<olb::opti::F>(9);
auto x262 = cell.template getFieldComponent<olb::opti::F>(0);
auto x265 = cell.template getFieldComponent<olb::opti::F>(3);
auto x279 = cell.template getFieldComponent<olb::opti::F>(17);
auto x252 = cell.template getFieldComponent<olb::opti::DJDF>(9);
auto x277 = cell.template getFieldComponent<olb::opti::F>(15);
auto x256 = cell.template getFieldComponent<olb::opti::DJDF>(13);
auto x248 = cell.template getFieldComponent<olb::opti::DJDF>(5);
auto x246 = cell.template getFieldComponent<olb::opti::DJDF>(3);
auto x245 = cell.template getFieldComponent<olb::opti::DJDF>(2);
auto x283 = V{1}*x268;
auto x284 = V{1}*x270;
auto x285 = V{1}*x278;
auto x286 = V{1}*x280;
auto x287 = V{1}*x265;
auto x290 = x269 + x271 + x274;
auto x291 = x265 + x270 + x290;
auto x292 = x277 + x278;
auto x293 = x272 + x276 + x292;
auto x294 = x263 + x266 + x293;
auto x295 = x279 + x280;
auto x296 = x267 + x275;
auto x297 = x273 + x295 + x296;
auto x298 = x264 + x297;
auto x303 = V{1}*x269;
auto x304 = V{1}*x266;
auto x305 = V{1}*x267;
auto x306 = V{1}*x263;
auto x308 = V{1}*x276;
auto x309 = V{1}*x271;
auto x310 = V{1}*x264;
auto x0 = V{3};
auto x313 = ((x0)*(x0));
auto x1 = V{1}/x0;
auto x2 = -x283;
auto x3 = -x284;
auto x4 = -x285;
auto x288 = x277 + x4;
auto x5 = -x286;
auto x289 = x279 + x5;
auto x6 = x2 - x287 + x288 + x289 + x290 + x3;
auto x7 = x262 + x268 + x291 + x294 + x298 + V{1};
auto x299 = V{1} / (x7);
auto x300 = x242*x299;
auto x333 = V{1}/((x7)*(x7));
auto x8 = x300*x6;
auto x9 = x7*x8;
auto x10 = x1*(x7 + V{-1});
auto x301 = -x10;
auto x302 = x268 + x301;
auto x11 = x291 + x292 + x295 + x302 - x8*x9;
auto x12 = -x303;
auto x13 = -x304;
auto x14 = -x305;
auto x307 = x14 + x275;
auto x15 = x12 + x13 + x2 + x293 - x306 + x307;
auto x16 = x15*x300;
auto x17 = x16*x7;
auto x18 = x12 - x17*x8 + x268 + x288;
auto x19 = ((x18)*(x18));
auto x20 = -x308;
auto x21 = -x309;
auto x22 = x13 + x20 + x21 + x297 + x3 - x310;
auto x23 = x22*x300;
auto x24 = -x17*x23 + x20 + x266 + x307;
auto x311 = V{2}*x24;
auto x25 = x311;
auto x26 = -x16*x17 + x269 + x294 + x296 + x302;
auto x27 = x23*x7;
auto x28 = x21 - x27*x8 + x270 + x289;
auto x312 = V{2}*x28;
auto x29 = x312;
auto x30 = -x23*x27 + x266 + x270 + x271 + x276 + x298 + x301;
auto x31 = util::sqrt(V{2}*x19 + x24*x25 + x28*x29 + ((x11)*(x11)) + ((x26)*(x26)) + ((x30)*(x30)));
auto x32 = V{2.828428}*x313*((x282)*(x282));
auto x33 = x299*x32;
auto x34 = V{1}/x281;
auto x35 = util::sqrt(x31*x33 + ((x34)*(x34)));
auto x36 = V{0.5}*x34 + V{0.5}*x35;
auto x37 = V{0.0277777777777778};
auto x38 = x37*x7;
auto x39 = -V{1}*x8;
auto x40 = x23 + x39;
auto x41 = V{0.5}*x313;
auto x314 = x40*x41;
auto x321 = x41*x8;
auto x322 = x23*x41;
auto x323 = x16*x41;
auto x330 = x39*x41;
auto x42 = x314;
auto x43 = V{0.5}*x0;
auto x44 = x43*(((x16)*(x16)) + ((x23)*(x23)) + ((x8)*(x8)));
auto x315 = V{1} - x44;
auto x45 = x0*x40 + x315 + x40*x42;
auto x46 = x23 + x8;
auto x316 = x41*x46;
auto x47 = x316;
auto x48 = x0*x46 + x315 + x46*x47;
auto x49 = x16 + x39;
auto x317 = x41*x49;
auto x50 = x317;
auto x51 = x0*x49 + x315 + x49*x50;
auto x52 = x16 + x8;
auto x318 = x41*x52;
auto x53 = x318;
auto x54 = x0*x52 + x315 + x52*x53;
auto x55 = -V{1}*x23;
auto x331 = x41*x55;
auto x56 = x16 + x55;
auto x319 = x41*x56;
auto x57 = x319;
auto x58 = x0*x56 + x315 + x56*x57;
auto x59 = x16 + x23;
auto x320 = x41*x59;
auto x60 = x320;
auto x61 = x0*x59 + x315 + x59*x60;
auto x62 = V{0.0555555555555556};
auto x63 = x62*x7;
auto x64 = x321;
auto x65 = x0*x8 + x315 + x64*x8;
auto x66 = x322;
auto x67 = x0*x23 + x23*x66 + x315;
auto x68 = x323;
auto x69 = x0*x16 + x16*x68 + x315;
auto x70 = x55 + x8;
auto x324 = x41*x70;
auto x71 = x324;
auto x72 = x0*x70 + x315 + x70*x71;
auto x73 = x39 + x55;
auto x325 = x41*x73;
auto x74 = x325;
auto x75 = x0*x73 + x315 + x73*x74;
auto x76 = -V{1}*x16;
auto x332 = x41*x76;
auto x77 = x76 + x8;
auto x326 = x41*x77;
auto x78 = x326;
auto x79 = x0*x77 + x315 + x77*x78;
auto x80 = x39 + x76;
auto x327 = x41*x80;
auto x81 = x327;
auto x82 = x0*x80 + x315 + x80*x81;
auto x83 = x23 + x76;
auto x328 = x41*x83;
auto x84 = x328;
auto x85 = x0*x83 + x315 + x83*x84;
auto x86 = x55 + x76;
auto x329 = x41*x86;
auto x87 = x329;
auto x88 = x0*x86 + x315 + x86*x87;
auto x89 = x330;
auto x90 = x0*x39 + x315 + x39*x89;
auto x91 = x331;
auto x92 = x0*x55 + x315 + x55*x91;
auto x93 = x332;
auto x94 = x0*x76 + x315 + x76*x93;
auto x95 = V{1} - x44;
auto x96 = V{0.333333333333333};
auto x97 = x7*x96;
auto x98 = V{0.25}*(V{1}*cell[0]*x262 - cell[0]*(x95*x97 - x96) + V{1}*cell[10]*x272 + cell[10]*(x62 - x63*x69) + V{1}*cell[11]*x273 + cell[11]*(x62 - x63*x67) + V{1}*cell[12]*x274 + cell[12]*(x62 - x63*x65) + V{1}*cell[13]*x275 + cell[13]*(x37 - x38*x61) + cell[14]*x308 + cell[14]*(x37 - x38*x58) + V{1}*cell[15]*x277 + cell[15]*(x37 - x38*x54) + cell[16]*x285 + cell[16]*(x37 - x38*x51) + V{1}*cell[17]*x279 + cell[17]*(x37 - x38*x48) + cell[18]*x286 + cell[18]*(x37 - x38*x45) + cell[1]*x306 + cell[1]*(x62 - x63*x94) + cell[2]*x310 + cell[2]*(x62 - x63*x92) + cell[3]*x287 + cell[3]*(x62 - x63*x90) + cell[4]*x304 + cell[4]*(x37 - x38*x88) + cell[5]*x305 + cell[5]*(x37 - x38*x85) + cell[6]*x283 + cell[6]*(x37 - x38*x82) + cell[7]*x303 + cell[7]*(x37 - x38*x79) + cell[8]*x284 + cell[8]*(x37 - x38*x75) + cell[9]*x309 + cell[9]*(x37 - x38*x72))/(x35*((x36)*(x36)));
auto x99 = V{0.5}*x33*x98/x31;
auto x100 = x11*x99;
auto x101 = V{2}*x100;
auto x102 = -V{1}*x101;
auto x103 = x102*x8;
auto x104 = x99*(x29 + x312);
auto x105 = -V{1}*x104;
auto x106 = x105*x8;
auto x107 = x30*x99;
auto x108 = V{2}*x107;
auto x109 = -V{1}*x108;
auto x110 = x109*x23;
auto x111 = x18*x99;
auto x112 = V{4}*x111;
auto x113 = -V{1}*x112;
auto x114 = x113*x8;
auto x115 = x99*(x25 + x311);
auto x116 = -V{1}*x115;
auto x117 = x116*x23;
auto x118 = x26*x99;
auto x119 = V{2}*x118;
auto x120 = -V{1}*x119;
auto x121 = x120*x16;
auto x122 = x1*x102 + x1*x109 + x1*x120 + x103*x8 + x106*x23 + x110*x23 + x114*x16 + x117*x16 + x121*x16;
auto x123 = V{1}/x36;
auto x124 = cell[18]*x123;
auto x125 = cell[17]*x123;
auto x126 = cell[16]*x123;
auto x127 = cell[15]*x123;
auto x128 = cell[14]*x123;
auto x129 = cell[13]*x123;
auto x130 = cell[12]*x123;
auto x131 = cell[11]*x123;
auto x132 = cell[10]*x123;
auto x133 = cell[9]*x123;
auto x134 = cell[8]*x123;
auto x135 = cell[7]*x123;
auto x136 = cell[6]*x123;
auto x137 = cell[5]*x123;
auto x138 = cell[4]*x123;
auto x139 = cell[3]*x123;
auto x140 = cell[2]*x123;
auto x141 = cell[1]*x123;
auto x142 = cell[0]*x123;
auto x143 = x124*x37*x45 + x125*x37*x48 + x126*x37*x51 + x127*x37*x54 + x128*x37*x58 + x129*x37*x61 + x130*x62*x65 + x131*x62*x67 + x132*x62*x69 + x133*x37*x72 + x134*x37*x75 + x135*x37*x79 + x136*x37*x82 + x137*x37*x85 + x138*x37*x88 + x139*x62*x90 + x140*x62*x92 + x141*x62*x94 + x142*x95*x96;
auto x144 = V{1} - x123;
auto x145 = -x333*x6;
auto x146 = x124*x38;
auto x147 = x125*x38;
auto x148 = x126*x38;
auto x149 = x127*x38;
auto x150 = x128*x38;
auto x151 = x129*x38;
auto x152 = x130*x63;
auto x153 = x131*x63;
auto x154 = x132*x63;
auto x155 = x133*x38;
auto x156 = x134*x38;
auto x157 = x135*x38;
auto x158 = x136*x38;
auto x159 = x137*x38;
auto x160 = x138*x38;
auto x161 = x139*x63;
auto x334 = V{1}*x161;
auto x162 = x140*x63;
auto x337 = V{1}*x162;
auto x163 = x141*x63;
auto x339 = V{1}*x163;
auto x164 = -V{1}*x43*(x142*x97 + x146 + x147 + x148 + x149 + x150 + x151 + x152 + x153 + x154 + x155 + x156 + x157 + x158 + x159 + x160 + x161 + x162 + x163);
auto x165 = x164*x8;
auto x166 = x146*(x0 + x314 + x42);
auto x167 = x147*(x0 + x316 + x47);
auto x168 = x148*(x0 + x317 + x50);
auto x169 = x149*(x0 + x318 + x53);
auto x170 = x155*(x0 + x324 + x71);
auto x171 = -V{1}*x156*(x0 + x325 + x74);
auto x335 = x167 + x171;
auto x172 = x157*(x0 + x326 + x78);
auto x173 = -V{1}*x158*(x0 + x327 + x81);
auto x336 = x169 + x173;
auto x174 = x242*(x0*x152 - x0*x334 + x152*x321 + x152*x64 + V{2}*x165 - V{1}*x166 - V{1}*x168 + x170 + x172 - x330*x334 - x334*x89 + x335 + x336);
auto x175 = -x22*x333;
auto x176 = x164*x23;
auto x177 = x150*(x0 + x319 + x57);
auto x178 = x151*(x0 + x320 + x60);
auto x179 = x159*(x0 + x328 + x84);
auto x180 = -V{1}*x160*(x0 + x329 + x87);
auto x338 = x178 + x180;
auto x181 = x242*(x0*x153 - x0*x337 + x153*x322 + x153*x66 + x166 - V{1}*x170 + V{2}*x176 - V{1}*x177 + x179 - x331*x337 + x335 - x337*x91 + x338);
auto x182 = -x15*x333;
auto x183 = x16*x164;
auto x184 = x242*(x0*x154 - x0*x339 + x154*x323 + x154*x68 + x168 - V{1}*x172 + x177 - V{1}*x179 + V{2}*x183 - x332*x339 + x336 + x338 - x339*x93);
auto x185 = x145*x174 + x175*x181 + x182*x184;
auto x186 = -x31*x32*x333*x98;
auto x187 = x242*(x102*x9 + x103*x7 + x105*x27 + x113*x17);
auto x188 = x242*(x106*x7 + x109*x27 + x110*x7 + x116*x17);
auto x189 = x242*(x114*x7 + x117*x7 + x120*x17 + x121*x7);
auto x190 = x145*x187 + x175*x188 + x182*x189;
auto x340 = x122 + x143 + x185 + x186 + x190;
auto x341 = x119 + x340;
auto x343 = x108 + x340;
auto x346 = x101 + x340;
auto x191 = cell[0]*x144 - x243 + x340;
auto x192 = V{1}*x299;
auto x193 = x184*x192;
auto x194 = -V{1}*x193;
auto x195 = x189*x192;
auto x353 = x193 + x195 + x341;
auto x355 = x108 + x353;
auto x356 = x101 + x353;
auto x196 = -V{1}*x195;
auto x342 = x194 + x196 + x341;
auto x348 = x108 + x342;
auto x350 = x101 + x342;
auto x197 = cell[1]*x144 - x244 + x342;
auto x198 = x181*x192;
auto x199 = -V{1}*x198;
auto x200 = x188*x192;
auto x349 = x198 + x200;
auto x354 = x343 + x349;
auto x357 = x101 + x354;
auto x201 = -V{1}*x200;
auto x344 = x199 + x201;
auto x345 = x343 + x344;
auto x352 = x101 + x345;
auto x202 = cell[2]*x144 - x245 + x345;
auto x203 = x174*x192;
auto x204 = -V{1}*x203;
auto x205 = x187*x192;
auto x351 = x203 + x205;
auto x206 = -V{1}*x205;
auto x347 = x204 + x206;
auto x207 = cell[3]*x144 - x246 + x346 + x347;
auto x208 = cell[4]*x144 + x115 - x247 + x344 + x348;
auto x209 = cell[5]*x144 + x116 - x248 + x348 + x349;
auto x210 = cell[6]*x144 + x112 - x249 + x347 + x350;
auto x211 = cell[7]*x144 + x113 - x250 + x350 + x351;
auto x212 = cell[8]*x144 + x104 - x251 + x347 + x352;
auto x213 = cell[9]*x144 + x105 - x252 + x351 + x352;
auto x214 = cell[10]*x144 - x253 + x353;
auto x215 = cell[11]*x144 - x254 + x354;
auto x216 = cell[12]*x144 - x255 + x346 + x351;
auto x217 = cell[13]*x144 + x115 - x256 + x349 + x355;
auto x218 = cell[14]*x144 + x116 - x257 + x344 + x355;
auto x219 = cell[15]*x144 + x112 - x258 + x351 + x356;
auto x220 = cell[16]*x144 + x113 - x259 + x347 + x356;
auto x221 = cell[17]*x144 + x104 - x260 + x351 + x357;
auto x222 = cell[18]*x144 + x105 - x261 + x347 + x357;
cell[0] = x191;
cell[1] = x197;
cell[2] = x202;
cell[3] = x207;
cell[4] = x208;
cell[5] = x209;
cell[6] = x210;
cell[7] = x211;
cell[8] = x212;
cell[9] = x213;
cell[10] = x214;
cell[11] = x215;
cell[12] = x216;
cell[13] = x217;
cell[14] = x218;
cell[15] = x219;
cell[16] = x220;
cell[17] = x221;
cell[18] = x222;
return { x191 + V{1}, ((x191)*(x191)) + ((x197)*(x197)) + ((x202)*(x202)) + ((x207)*(x207)) + ((x208)*(x208)) + ((x209)*(x209)) + ((x210)*(x210)) + ((x211)*(x211)) + ((x212)*(x212)) + ((x213)*(x213)) + ((x214)*(x214)) + ((x215)*(x215)) + ((x216)*(x216)) + ((x217)*(x217)) + ((x218)*(x218)) + ((x219)*(x219)) + ((x220)*(x220)) + ((x221)*(x221)) + ((x222)*(x222)) };
}
};

}

}

// Generation Info: commit=4d59a7668c481a232fca09b0478b2faf6a58e0e4
