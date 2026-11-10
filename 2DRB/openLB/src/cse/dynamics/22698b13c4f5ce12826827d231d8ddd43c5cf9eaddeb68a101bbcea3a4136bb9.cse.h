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
struct CSE<Dual<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::PlainGuo> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x280 = cell.template getFieldComponent<olb::opti::F>(0);
auto x261 = cell.template getFieldComponent<olb::opti::DJDF>(0);
auto x272 = cell.template getFieldComponent<olb::opti::DJDF>(11);
auto x260 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x268 = cell.template getFieldComponent<olb::opti::DJDF>(7);
auto x263 = cell.template getFieldComponent<olb::opti::DJDF>(2);
auto x275 = cell.template getFieldComponent<olb::opti::DJDF>(14);
auto x269 = cell.template getFieldComponent<olb::opti::DJDF>(8);
auto x299 = parameters.template get<descriptors::OMEGA>();
auto x270 = cell.template getFieldComponent<olb::opti::DJDF>(9);
auto x271 = cell.template getFieldComponent<olb::opti::DJDF>(10);
auto x297 = cell.template getFieldComponent<olb::opti::F>(17);
auto x276 = cell.template getFieldComponent<olb::opti::DJDF>(15);
auto x295 = cell.template getFieldComponent<olb::opti::F>(15);
auto x287 = cell.template getFieldComponent<olb::opti::F>(7);
auto x265 = cell.template getFieldComponent<olb::opti::DJDF>(4);
auto x258 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x279 = cell.template getFieldComponent<olb::opti::DJDF>(18);
auto x289 = cell.template getFieldComponent<olb::opti::F>(9);
auto x288 = cell.template getFieldComponent<olb::opti::F>(8);
auto x282 = cell.template getFieldComponent<olb::opti::F>(2);
auto x259 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x277 = cell.template getFieldComponent<olb::opti::DJDF>(16);
auto x283 = cell.template getFieldComponent<olb::opti::F>(3);
auto x264 = cell.template getFieldComponent<olb::opti::DJDF>(3);
auto x274 = cell.template getFieldComponent<olb::opti::DJDF>(13);
auto x266 = cell.template getFieldComponent<olb::opti::DJDF>(5);
auto x296 = cell.template getFieldComponent<olb::opti::F>(16);
auto x290 = cell.template getFieldComponent<olb::opti::F>(10);
auto x285 = cell.template getFieldComponent<olb::opti::F>(5);
auto x284 = cell.template getFieldComponent<olb::opti::F>(4);
auto x278 = cell.template getFieldComponent<olb::opti::DJDF>(17);
auto x294 = cell.template getFieldComponent<olb::opti::F>(14);
auto x298 = cell.template getFieldComponent<olb::opti::F>(18);
auto x286 = cell.template getFieldComponent<olb::opti::F>(6);
auto x273 = cell.template getFieldComponent<olb::opti::DJDF>(12);
auto x267 = cell.template getFieldComponent<olb::opti::DJDF>(6);
auto x291 = cell.template getFieldComponent<olb::opti::F>(11);
auto x292 = cell.template getFieldComponent<olb::opti::F>(12);
auto x262 = cell.template getFieldComponent<olb::opti::DJDF>(1);
auto x281 = cell.template getFieldComponent<olb::opti::F>(1);
auto x293 = cell.template getFieldComponent<olb::opti::F>(13);
auto x302 = x290 + x294 + x296;
auto x303 = x285 + x291 + x298;
auto x304 = x287 + x289 + x292;
auto x0 = V{1} - V{0.5}*x299;
auto x1 = V{0.0277777777777778};
auto x2 = V{3};
auto x347 = V{1}*x2;
auto x348 = x260*x347;
auto x353 = x259*x347;
auto x355 = x258*x347;
auto x3 = -V{1}*x284;
auto x300 = x293 + x3;
auto x4 = -V{1}*x286;
auto x301 = x295 + x4;
auto x5 = -V{1}*x281 - V{1}*x285 - V{1}*x287 + x300 + x301 + x302;
auto x6 = x280 + x281 + x282 + x283 + x284 + x286 + x288 + x293 + x295 + x297 + x302 + x303 + x304 + V{1};
auto x305 = V{1} / (x6);
auto x336 = V{1}/((x6)*(x6));
auto x7 = V{0.5}*x258 + x305*x5;
auto x306 = x2*x7;
auto x8 = -x258*x306;
auto x9 = -V{1}*x288;
auto x307 = x297 + x9;
auto x10 = -V{1}*x282 - V{1}*x289 - V{1}*x294 + x300 + x303 + x307;
auto x11 = x10*x305 + V{0.5}*x259;
auto x308 = x11*x2;
auto x12 = -x2*(x11 + V{-1});
auto x13 = ((x2)*(x2));
auto x310 = V{1}*x13;
auto x14 = -V{1}*x283 - V{1}*x296 - V{1}*x298 + x301 + x304 + x307;
auto x15 = x14*x305 + V{0.5}*x260;
auto x309 = x15*x2;
auto x16 = -V{1}*x15;
auto x17 = x11 + x16;
auto x18 = x13*x17;
auto x19 = -x2*(x15 + V{1});
auto x315 = -x19;
auto x20 = x11 + x15;
auto x21 = x13*x20;
auto x22 = -x2*(x15 + V{-1});
auto x23 = -x259*x308;
auto x24 = -x2*(x7 + V{-1});
auto x25 = x16 + x7;
auto x26 = x13*x25;
auto x27 = x15 + x7;
auto x28 = x13*x27;
auto x29 = -x260*x309;
auto x312 = x23 + x29;
auto x313 = x29 + x8;
auto x30 = -V{1}*x11;
auto x31 = x30 + x7;
auto x32 = x13*x31;
auto x33 = -x2*(x11 + V{1});
auto x314 = -x33;
auto x34 = x11 + x7;
auto x35 = x13*x34;
auto x36 = V{0.0555555555555556};
auto x37 = x23 + x8;
auto x38 = x15 + x30;
auto x39 = x13*x38;
auto x40 = x16 + x30;
auto x41 = -x310*x40;
auto x42 = -x2*(x7 + V{1});
auto x43 = -V{1}*x7;
auto x44 = x15 + x43;
auto x45 = x13*x44;
auto x46 = x16 + x43;
auto x47 = -x310*x46;
auto x48 = x11 + x43;
auto x49 = x13*x48;
auto x50 = x30 + x43;
auto x51 = -x310*x50;
auto x52 = V{0.333333333333333};
auto x311 = cell[0]*x52;
auto x53 = x0*(cell[10]*x36*(x258*(x13*x7 + x24) + x312) + cell[11]*x36*(x259*(x11*x13 + x12) + x313) + cell[12]*x36*(x260*(x13*x15 + x22) + x37) + cell[13]*x1*(x258*(x24 + x35) + x259*(x12 + x35) + x29) + cell[14]*x1*(x258*(x24 + x32) - x259*(x314 + V{1}*x32) + x29) + cell[15]*x1*(x23 + x258*(x24 + x28) + x260*(x22 + x28)) + cell[16]*x1*(x23 + x258*(x24 + x26) + x260*(x19 - V{1}*x26)) + cell[17]*x1*(x259*(x12 + x21) + x260*(x21 + x22) + x8) + cell[18]*x1*(x259*(x12 + x18) - x260*(V{1}*x18 + x315) + x8) + cell[1]*x36*(-x258*(x310*x43 - x42) + x312) + cell[2]*x36*(-x259*(x30*x310 + x314) + x313) - cell[3]*x36*(x260*(x16*x310 + x315) - x37) + cell[4]*x1*(x258*(x42 + x51) + x259*(x33 + x51) + x29) + cell[5]*x1*(x258*(x42 - V{1}*x49) + x259*(x12 + x49) + x29) + cell[6]*x1*(x23 + x258*(x42 + x47) + x260*(x19 + x47)) + cell[7]*x1*(x23 + x258*(x42 - V{1}*x45) + x260*(x22 + x45)) + cell[8]*x1*(x259*(x33 + x41) + x260*(x19 + x41) + x8) + cell[9]*x1*(x259*(x33 - V{1}*x39) + x260*(x22 + x39) + x8) + x311*(x29 + x37));
auto x54 = V{0.5}*x13;
auto x316 = x17*x54;
auto x317 = x20*x54;
auto x318 = x25*x54;
auto x319 = x27*x54;
auto x320 = x31*x54;
auto x321 = x34*x54;
auto x322 = x15*x54;
auto x323 = x11*x54;
auto x324 = x54*x7;
auto x325 = x38*x54;
auto x326 = x40*x54;
auto x327 = x44*x54;
auto x328 = x46*x54;
auto x329 = x48*x54;
auto x330 = x50*x54;
auto x331 = x16*x54;
auto x332 = x30*x54;
auto x333 = x43*x54;
auto x55 = x316;
auto x56 = V{0.5}*x2;
auto x57 = x56*(((x11)*(x11)) + ((x15)*(x15)) + ((x7)*(x7)));
auto x334 = V{1} - x57;
auto x58 = cell[18]*x299;
auto x59 = x317;
auto x60 = cell[17]*x299;
auto x61 = x318;
auto x62 = cell[16]*x299;
auto x63 = x319;
auto x64 = cell[15]*x299;
auto x65 = x320;
auto x66 = cell[14]*x299;
auto x67 = x321;
auto x68 = cell[13]*x299;
auto x69 = x322;
auto x70 = cell[12]*x299;
auto x71 = x323;
auto x72 = cell[11]*x299;
auto x73 = x324;
auto x74 = cell[10]*x299;
auto x75 = x325;
auto x76 = cell[9]*x299;
auto x77 = x326;
auto x78 = cell[8]*x299;
auto x79 = x327;
auto x80 = cell[7]*x299;
auto x81 = x328;
auto x82 = cell[6]*x299;
auto x83 = x329;
auto x84 = cell[5]*x299;
auto x85 = x330;
auto x86 = cell[4]*x299;
auto x87 = x331;
auto x88 = cell[3]*x299;
auto x89 = x332;
auto x90 = cell[2]*x299;
auto x91 = x333;
auto x92 = cell[1]*x299;
auto x93 = cell[0]*x299;
auto x335 = x52*x93;
auto x94 = x1*x58*(x17*x2 + x17*x55 + x334) + x1*x60*(x2*x20 + x20*x59 + x334) + x1*x62*(x2*x25 + x25*x61 + x334) + x1*x64*(x2*x27 + x27*x63 + x334) + x1*x66*(x2*x31 + x31*x65 + x334) + x1*x68*(x2*x34 + x334 + x34*x67) + x1*x76*(x2*x38 + x334 + x38*x75) + x1*x78*(x2*x40 + x334 + x40*x77) + x1*x80*(x2*x44 + x334 + x44*x79) + x1*x82*(x2*x46 + x334 + x46*x81) + x1*x84*(x2*x48 + x334 + x48*x83) + x1*x86*(x2*x50 + x334 + x50*x85) - x335*(x57 + V{-1}) + x36*x70*(x15*x69 + x309 + x334) + x36*x72*(x11*x71 + x308 + x334) + x36*x74*(x306 + x334 + x7*x73) + x36*x88*(x16*x2 + x16*x87 + x334) + x36*x90*(x2*x30 + x30*x89 + x334) + x36*x92*(x2*x43 + x334 + x43*x91);
auto x95 = V{1} - x299;
auto x96 = -x14*x336;
auto x97 = x1*x6;
auto x98 = x58*x97;
auto x99 = x60*x97;
auto x100 = x62*x97;
auto x101 = x64*x97;
auto x102 = x66*x97;
auto x103 = x68*x97;
auto x104 = x36*x6;
auto x105 = x104*x70;
auto x106 = x104*x72;
auto x107 = x104*x74;
auto x108 = x76*x97;
auto x109 = x78*x97;
auto x110 = x80*x97;
auto x111 = x82*x97;
auto x112 = x84*x97;
auto x113 = x86*x97;
auto x114 = x104*x88;
auto x337 = V{1}*x114;
auto x115 = x104*x90;
auto x340 = V{1}*x115;
auto x116 = x104*x92;
auto x342 = V{1}*x116;
auto x117 = -V{1}*x56*(x100 + x101 + x102 + x103 + x105 + x106 + x107 + x108 + x109 + x110 + x111 + x112 + x113 + x114 + x115 + x116 + x335*x6 + x98 + x99);
auto x118 = x117*x15;
auto x119 = x98*(x2 + x316 + x55);
auto x120 = x99*(x2 + x317 + x59);
auto x121 = x100*(x2 + x318 + x61);
auto x122 = x101*(x2 + x319 + x63);
auto x123 = x108*(x2 + x325 + x75);
auto x124 = -V{1}*x109*(x2 + x326 + x77);
auto x338 = x120 + x124;
auto x125 = x110*(x2 + x327 + x79);
auto x126 = -V{1}*x111*(x2 + x328 + x81);
auto x339 = x122 + x126;
auto x127 = x105*x2 + x105*x322 + x105*x69 + V{2}*x118 - V{1}*x119 - V{1}*x121 + x123 + x125 - x2*x337 - x331*x337 - x337*x87 + x338 + x339;
auto x128 = -x10*x336;
auto x129 = x11*x117;
auto x130 = x102*(x2 + x320 + x65);
auto x131 = x103*(x2 + x321 + x67);
auto x132 = x112*(x2 + x329 + x83);
auto x133 = -V{1}*x113*(x2 + x330 + x85);
auto x341 = x131 + x133;
auto x134 = x106*x2 + x106*x323 + x106*x71 + x119 - V{1}*x123 + V{2}*x129 - V{1}*x130 + x132 - x2*x340 - x332*x340 + x338 - x340*x89 + x341;
auto x135 = -x336*x5;
auto x136 = x117*x7;
auto x137 = x107*x2 + x107*x324 + x107*x73 + x121 - V{1}*x125 + x130 - V{1}*x132 + V{2}*x136 - x2*x342 - x333*x342 + x339 + x341 - x342*x91;
auto x138 = x127*x96 + x128*x134 + x135*x137;
auto x139 = x0*x6;
auto x140 = x1*x139;
auto x141 = cell[18]*x140;
auto x142 = x141*x260;
auto x343 = V{1}*x142;
auto x143 = x141*x259;
auto x144 = -x13*(-x143 + x343);
auto x145 = cell[17]*x140;
auto x146 = x145*x260;
auto x147 = x145*x259;
auto x148 = x13*(x146 + x147);
auto x149 = cell[16]*x140;
auto x150 = x149*x260;
auto x344 = V{1}*x150;
auto x151 = x149*x258;
auto x152 = -x13*(-x151 + x344);
auto x153 = cell[15]*x140;
auto x154 = x153*x260;
auto x155 = x153*x258;
auto x156 = x13*(x154 + x155);
auto x157 = cell[14]*x140;
auto x158 = cell[13]*x140;
auto x159 = x139*x36;
auto x160 = cell[12]*x159;
auto x161 = x160*x260;
auto x162 = cell[11]*x159;
auto x163 = cell[10]*x159;
auto x164 = cell[9]*x140;
auto x165 = x164*x260;
auto x166 = x164*x259;
auto x345 = V{1}*x166;
auto x167 = x13*(x165 - x345);
auto x168 = cell[8]*x140;
auto x169 = x168*x260;
auto x170 = x168*x259;
auto x171 = x310*(x169 + x170);
auto x349 = -x148 - x171;
auto x172 = cell[7]*x140;
auto x173 = x172*x260;
auto x174 = x172*x258;
auto x346 = V{1}*x174;
auto x175 = x13*(x173 - x346);
auto x176 = cell[6]*x140;
auto x177 = x176*x260;
auto x178 = x176*x258;
auto x179 = x310*(x177 + x178);
auto x350 = -x156 - x179;
auto x180 = cell[5]*x140;
auto x181 = cell[4]*x140;
auto x182 = cell[3]*x159;
auto x183 = x182*x260;
auto x184 = cell[2]*x159;
auto x185 = cell[1]*x159;
auto x186 = x139*x311;
auto x187 = x13*x161 + V{1}*x13*x183 - V{1}*x144 - x146*x347 - V{1}*x152 - x154*x347 - x157*x348 - x158*x348 - x161*x347 - x162*x348 - x163*x348 - x165*x347 + x167 - x169*x347 - x173*x347 + x175 - x177*x347 - x180*x348 - x181*x348 - x183*x347 - x184*x348 - x185*x348 - x186*x348 - x2*x343 - x2*x344 - x349 - x350;
auto x188 = x157*x259;
auto x351 = V{1}*x188;
auto x189 = x157*x258;
auto x190 = -x13*(-x189 + x351);
auto x191 = x158*x259;
auto x192 = x158*x258;
auto x193 = x13*(x191 + x192);
auto x194 = x162*x259;
auto x195 = x180*x259;
auto x196 = x180*x258;
auto x352 = V{1}*x196;
auto x197 = x13*(x195 - x352);
auto x198 = x181*x259;
auto x199 = x181*x258;
auto x200 = x310*(x198 + x199);
auto x354 = -x193 - x200;
auto x201 = x184*x259;
auto x202 = x13*x194 + V{1}*x13*x201 - x143*x347 + x144 - x147*x347 - x149*x353 - x153*x353 - x160*x353 - x163*x353 - V{1}*x167 - x170*x347 - x172*x353 - x176*x353 - x182*x353 - x185*x353 - x186*x353 - V{1}*x190 - x191*x347 - x194*x347 - x195*x347 + x197 - x198*x347 - x2*x345 - x2*x351 - x201*x347 - x349 - x354;
auto x203 = x163*x258;
auto x204 = x185*x258;
auto x205 = x13*x203 + V{1}*x13*x204 - x141*x355 - x145*x355 - x151*x347 + x152 - x155*x347 - x160*x355 - x162*x355 - x164*x355 - x168*x355 - V{1}*x175 - x178*x347 - x182*x355 - x184*x355 - x186*x355 - x189*x347 + x190 - x192*x347 - V{1}*x197 - x199*x347 - x2*x346 - x2*x352 - x203*x347 - x204*x347 - x350 - x354;
auto x206 = x128*x202 + x135*x205 + x187*x96;
auto x356 = x138 + x206 + x53 + x94;
auto x207 = cell[0]*x95 - x261 + x356;
auto x208 = V{1}*x305;
auto x209 = x137*x208;
auto x210 = -V{1}*x209;
auto x211 = x205*x208;
auto x364 = x209 + x211;
auto x365 = x356 + x364;
auto x212 = -V{1}*x211;
auto x357 = x210 + x212 + x356;
auto x213 = cell[1]*x95 - x262 + x357;
auto x214 = x134*x208;
auto x215 = -V{1}*x214;
auto x216 = x202*x208;
auto x362 = x214 + x216;
auto x366 = x356 + x362;
auto x217 = -V{1}*x216;
auto x358 = x215 + x217;
auto x359 = x356 + x358;
auto x218 = cell[2]*x95 - x263 + x359;
auto x219 = x127*x208;
auto x220 = -V{1}*x219;
auto x221 = x187*x208;
auto x363 = x219 + x221;
auto x222 = -V{1}*x221;
auto x360 = x220 + x222;
auto x361 = x356 + x360;
auto x223 = cell[3]*x95 - x264 + x361;
auto x224 = cell[4]*x95 - x265 + x357 + x358;
auto x225 = cell[5]*x95 - x266 + x357 + x362;
auto x226 = cell[6]*x95 - x267 + x357 + x360;
auto x227 = cell[7]*x95 - x268 + x357 + x363;
auto x228 = cell[8]*x95 - x269 + x359 + x360;
auto x229 = cell[9]*x95 - x270 + x359 + x363;
auto x230 = cell[10]*x95 - x271 + x365;
auto x231 = cell[11]*x95 - x272 + x366;
auto x232 = cell[12]*x95 - x273 + x356 + x363;
auto x233 = cell[13]*x95 - x274 + x362 + x365;
auto x234 = cell[14]*x95 - x275 + x359 + x364;
auto x235 = cell[15]*x95 - x276 + x363 + x365;
auto x236 = cell[16]*x95 - x277 + x361 + x364;
auto x237 = cell[17]*x95 - x278 + x363 + x366;
auto x238 = cell[18]*x95 - x279 + x361 + x362;
cell[0] = x207;
cell[1] = x213;
cell[2] = x218;
cell[3] = x223;
cell[4] = x224;
cell[5] = x225;
cell[6] = x226;
cell[7] = x227;
cell[8] = x228;
cell[9] = x229;
cell[10] = x230;
cell[11] = x231;
cell[12] = x232;
cell[13] = x233;
cell[14] = x234;
cell[15] = x235;
cell[16] = x236;
cell[17] = x237;
cell[18] = x238;
return { x207 + V{1}, ((x207)*(x207)) + ((x213)*(x213)) + ((x218)*(x218)) + ((x223)*(x223)) + ((x224)*(x224)) + ((x225)*(x225)) + ((x226)*(x226)) + ((x227)*(x227)) + ((x228)*(x228)) + ((x229)*(x229)) + ((x230)*(x230)) + ((x231)*(x231)) + ((x232)*(x232)) + ((x233)*(x233)) + ((x234)*(x234)) + ((x235)*(x235)) + ((x236)*(x236)) + ((x237)*(x237)) + ((x238)*(x238)) };
}
};

}

}

// Generation Info: commit=4ea71c45fdf920fe461e9926f71dc63f7fa1062d
