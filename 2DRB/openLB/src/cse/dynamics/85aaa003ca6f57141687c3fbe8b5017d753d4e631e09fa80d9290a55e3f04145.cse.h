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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::ThirdOrder, collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB>, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x21 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x20 = parameters.template get<descriptors::OMEGA>();
auto x19 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x22 = cell[15] + cell[17];
auto x23 = cell[12] + x22;
auto x24 = cell[11] + cell[18];
auto x25 = cell[10] + cell[14] + cell[16];
auto x26 = cell[2] + cell[8] + cell[9];
auto x27 = cell[13] + cell[3];
auto x28 = cell[0] + cell[1] + cell[4] + cell[5] + cell[6] + cell[7] + x23 + x24 + x25 + x26 + x27;
auto x29 = x28 + V{1};
auto x30 = V{1} / (x29);
auto x31 = V{0.333333333333333}*cell[13];
auto x32 = V{0.333333333333333}*cell[14];
auto x33 = V{0.333333333333333}*cell[4];
auto x34 = V{0.333333333333333}*cell[5];
auto x35 = V{0.666666666666667}*cell[12];
auto x36 = V{0.666666666666667}*cell[3];
auto x37 = V{1}*x30;
auto x38 = ((x19)*(x19));
auto x39 = -cell[18];
auto x40 = -cell[3];
auto x41 = -cell[8];
auto x42 = cell[9] + x41;
auto x43 = x39 + x40 + x42;
auto x44 = -cell[6];
auto x45 = cell[7] + x44;
auto x46 = -cell[16] + x45;
auto x47 = x23 + x43 + x46;
auto x48 = ((x47)*(x47));
auto x49 = x38*x48;
auto x50 = x37*x49;
auto x51 = V{0.333333333333333}*cell[0];
auto x52 = V{0.333333333333333}*cell[10];
auto x53 = V{0.333333333333333}*cell[1];
auto x54 = V{0.666666666666667}*cell[17];
auto x55 = V{0.666666666666667}*cell[18];
auto x56 = V{0.666666666666667}*cell[8];
auto x57 = V{0.666666666666667}*cell[9];
auto x58 = x51 + x52 + x53 - x54 - x55 - x56 - x57;
auto x59 = V{0.333333333333333}*cell[11];
auto x60 = V{0.333333333333333}*cell[2];
auto x61 = V{0.666666666666667}*cell[15];
auto x62 = V{0.666666666666667}*cell[16];
auto x63 = V{0.666666666666667}*cell[6];
auto x64 = V{0.666666666666667}*cell[7];
auto x65 = x59 + x60 - x61 - x62 - x63 - x64;
auto x66 = x31 + x32 + x33 + x34 - x35 - x36 + x50 + x58 + x65;
auto x67 = V{0.333333333333333}*cell[15];
auto x68 = V{0.333333333333333}*cell[16];
auto x69 = V{0.333333333333333}*cell[6];
auto x70 = V{0.333333333333333}*cell[7];
auto x71 = V{0.666666666666667}*cell[11];
auto x72 = V{0.666666666666667}*cell[2];
auto x73 = cell[13] + cell[17];
auto x74 = -cell[2];
auto x75 = -cell[9];
auto x76 = x24 + x41 + x74 + x75;
auto x77 = -cell[4];
auto x78 = cell[5] + x77;
auto x79 = -cell[14] + x78;
auto x80 = x73 + x76 + x79;
auto x81 = ((x80)*(x80));
auto x82 = x38*x81;
auto x83 = x37*x82;
auto x84 = V{0.333333333333333}*cell[12];
auto x85 = V{0.333333333333333}*cell[3];
auto x86 = V{0.666666666666667}*cell[13];
auto x87 = V{0.666666666666667}*cell[14];
auto x88 = V{0.666666666666667}*cell[4];
auto x89 = V{0.666666666666667}*cell[5];
auto x90 = x84 + x85 - x86 - x87 - x88 - x89;
auto x91 = x58 + x67 + x68 + x69 + x70 - x71 - x72 + x83 + x90;
auto x92 = V{0.333333333333333}*cell[17];
auto x93 = V{0.333333333333333}*cell[18];
auto x94 = V{0.333333333333333}*cell[8];
auto x95 = V{0.333333333333333}*cell[9];
auto x96 = V{0.666666666666667}*cell[10];
auto x97 = V{0.666666666666667}*cell[1];
auto x98 = cell[13] + cell[15];
auto x99 = -cell[1];
auto x100 = -cell[7];
auto x101 = x100 + x44 + x99;
auto x102 = -cell[5] + x77;
auto x103 = x101 + x102 + x25 + x98;
auto x104 = ((x103)*(x103));
auto x105 = x104*x38;
auto x106 = x105*x37;
auto x107 = x106 + x51 + x65 + x90 + x92 + x93 + x94 + x95 - x96 - x97;
auto x108 = x103*x30;
auto x109 = x108*x80;
auto x110 = -cell[13] + cell[14] + x109*x38 + x78;
auto x111 = x108*x47;
auto x112 = -cell[15] + cell[16];
auto x113 = x111*x38 + x112 + x45;
auto x114 = x30*x80;
auto x115 = x114*x47;
auto x116 = -cell[17];
auto x117 = cell[18] + x116;
auto x118 = x115*x38 + x117 + x42;
auto x119 = V{1} - V{1} / (V{3.00000046417339}*util::sqrt(x30*((x21)*(x21))*util::sqrt(V{0.5}*((x107)*(x107)) + ((x110)*(x110)) + ((x113)*(x113)) + ((x118)*(x118)) + V{0.5}*((x66)*(x66)) + V{0.5}*((x91)*(x91))) + V{0.0277777691819762}/((x20)*(x20))) + V{0.5}/x20);
auto x120 = V{0.5}*x30;
auto x121 = V{1} / ((x29)*(x29));
auto x122 = V{1.5}*x121;
auto x123 = x122*x49;
auto x124 = x105*x122;
auto x125 = x122*x82;
auto x126 = x124 + x125 + V{-1};
auto x127 = x123 + x126;
auto x128 = x120*x19;
auto x129 = x103*x128;
auto x130 = x19*x37;
auto x131 = x130*x80;
auto x132 = x130*x47;
auto x133 = V{0.166666666666667}*x30;
auto x134 = V{6.93889390390723e-18}*cell[0];
auto x135 = V{0.0833333333333333}*cell[12];
auto x136 = V{0.0833333333333333}*cell[3];
auto x137 = V{0.0833333333333333}*x30;
auto x138 = -V{0.0833333333333333}*cell[13] - V{0.0833333333333333}*cell[14] - V{0.0833333333333333}*cell[4] - V{0.0833333333333333}*cell[5] + x134 + x135 + x136 - x137*x49;
auto x139 = V{0.0833333333333333}*cell[11];
auto x140 = V{0.0833333333333333}*cell[2];
auto x141 = -V{0.0833333333333333}*cell[15] - V{0.0833333333333333}*cell[16] - V{0.0833333333333333}*cell[6] - V{0.0833333333333333}*cell[7] - x137*x82 + x139 + x140;
auto x142 = -V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] - V{0.166666666666667}*cell[1] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x105*x133 + x138 + x141;
auto x143 = x28 + V{1};
auto x144 = V{3}*x121;
auto x145 = x105*x144;
auto x146 = -x123;
auto x147 = V{1} - x125;
auto x148 = util::pow(x19, 3)/util::pow(x29, 3);
auto x149 = V{6.000003}*x148;
auto x150 = x103*x149;
auto x151 = V{2.999997}*x148;
auto x152 = x151*x80;
auto x153 = V{3}*x19;
auto x154 = x108*x153;
auto x155 = -x154;
auto x156 = x150*x81;
auto x157 = x104*x152;
auto x158 = x155 + x156 + x157;
auto x159 = x150*x48 + x152*x48 + x158;
auto x160 = x128*x80;
auto x161 = x103*x130;
auto x162 = V{0.0833333333333333}*cell[10];
auto x163 = V{0.0833333333333333}*cell[1];
auto x164 = -V{0.0833333333333333}*cell[17] - V{0.0833333333333333}*cell[18] - V{0.0833333333333333}*cell[8] - V{0.0833333333333333}*cell[9] - x105*x137 + x162 + x163;
auto x165 = -V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] - V{0.166666666666667}*cell[2] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] + x133*x82 + x138 + x164;
auto x166 = -x124;
auto x167 = x144*x82;
auto x168 = x149*x80;
auto x169 = x103*x151;
auto x170 = x114*x153;
auto x171 = -x170;
auto x172 = x104*x168;
auto x173 = x169*x81;
auto x174 = x171 + x172 + x173;
auto x175 = x168*x48 + x169*x48 + x174;
auto x176 = x128*x47;
auto x177 = -V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] - V{0.166666666666667}*cell[3] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] + x133*x49 + x134 + x141 + x164;
auto x178 = x144*x49;
auto x179 = x30*x47;
auto x180 = x153*x179;
auto x181 = -x180;
auto x182 = V{9}*x148;
auto x183 = x182*x47;
auto x184 = x104*x183 + x181 + x183*x81;
auto x185 = x147 + x166;
auto x186 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x187 = V{4.5}*x121;
auto x188 = cell[10] + cell[16] + x101;
auto x189 = x187*x38*((2*cell[13] - 2*cell[4] + x188 + x22 + x76)*(2*cell[13] - 2*cell[4] + x188 + x22 + x76));
auto x190 = x127 + x154;
auto x191 = V{18}*x148;
auto x192 = -x103*x182*x48 + x103*x191*x81 + x104*x191*x80 + x170 - x182*x48*x80;
auto x193 = -x135;
auto x194 = -x136;
auto x195 = V{0.0416666666666667}*x30;
auto x196 = x195*x49;
auto x197 = V{0.25000025}*x19;
auto x198 = x108*x197;
auto x199 = x198*x91;
auto x200 = V{0.5000005}*x19;
auto x201 = x114*x200;
auto x202 = x110*x201;
auto x203 = V{2.49999999985601e-07}*x19;
auto x204 = x108*x203;
auto x205 = x204*x66;
auto x206 = V{4.99999999971202e-07}*x19;
auto x207 = x179*x206;
auto x208 = x113*x207;
auto x209 = -V{0.0416666666666667}*cell[0];
auto x210 = V{0.0833333333333333}*x30;
auto x211 = V{0.0416666666666667}*cell[10] + V{6.93889390390723e-18}*cell[17] + V{6.93889390390723e-18}*cell[18] + V{0.0416666666666667}*cell[1] + V{6.93889390390723e-18}*cell[8] + V{6.93889390390723e-18}*cell[9] - x105*x210 + x209;
auto x212 = V{0.0416666666666667}*cell[11] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[16] + V{0.0416666666666667}*cell[2] + V{6.93889390390723e-18}*cell[6] + V{6.93889390390723e-18}*cell[7] - x210*x82;
auto x213 = x193 + x194 + x196 + x199 + x202 - x205 - x208 + x211 + x212;
auto x214 = x114*x197;
auto x215 = x107*x214;
auto x216 = x108*x200;
auto x217 = x110*x216;
auto x218 = x114*x203;
auto x219 = x218*x66;
auto x220 = x118*x207;
auto x221 = x215 + x217 - x219 - x220;
auto x222 = V{0.25}*x109*x38;
auto x223 = V{0.375}*cell[13] - V{0.125}*cell[14] + V{0.375}*cell[4] - V{0.125}*cell[5] - x222;
auto x224 = -cell[11] + V{2}*cell[14] + cell[15] - V{2}*cell[5] + x116 + x188 + x26 + x39;
auto x225 = V{3.000006}*x148;
auto x226 = x225*x48*x80;
auto x227 = V{6.000012}*x148;
auto x228 = x103*x227*x81;
auto x229 = x103*x225*x48;
auto x230 = x104*x227*x80;
auto x231 = -x215 - x217 + x219 + x220;
auto x232 = -V{0.125}*cell[13] + V{0.375}*cell[14] - V{0.125}*cell[4] + V{0.375}*cell[5] + x222;
auto x233 = cell[10] + cell[14] + x102 + x99;
auto x234 = x187*x38*((cell[12] + 2*cell[15] - 2*cell[6] + x233 + x43 + x73)*(cell[12] + 2*cell[15] - 2*cell[6] + x233 + x43 + x73));
auto x235 = x148*x47;
auto x236 = V{9.000009}*x235;
auto x237 = x104*x236;
auto x238 = V{8.99999999948164e-06}*x235;
auto x239 = x238*x81;
auto x240 = x180 + x237 - x239;
auto x241 = V{5.999994}*x148;
auto x242 = x241*x48*x80;
auto x243 = V{12.000006}*x148;
auto x244 = x103*x243*x48;
auto x245 = -x156;
auto x246 = -x157;
auto x247 = x190 + x242 + x244 + x245 + x246;
auto x248 = -x139;
auto x249 = -x140;
auto x250 = x195*x82;
auto x251 = x198*x66;
auto x252 = x179*x200;
auto x253 = x113*x252;
auto x254 = x204*x91;
auto x255 = x114*x206;
auto x256 = x110*x255;
auto x257 = V{0.0416666666666667}*cell[12] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[14] + V{0.0416666666666667}*cell[3] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[5] - x210*x49;
auto x258 = x211 + x248 + x249 + x250 + x251 + x253 - x254 - x256 + x257;
auto x259 = x179*x197;
auto x260 = x107*x259;
auto x261 = x113*x216;
auto x262 = x179*x203;
auto x263 = x262*x91;
auto x264 = x118*x255;
auto x265 = x260 + x261 - x263 - x264;
auto x266 = V{0.25}*x111*x38;
auto x267 = V{0.375}*cell[15] - V{0.125}*cell[16] + V{0.375}*cell[6] - V{0.125}*cell[7] - x266;
auto x268 = -cell[12] + x27;
auto x269 = V{2}*cell[16] - V{2}*cell[7] + cell[8] + x117 + x233 + x268 + x75;
auto x270 = -x260 - x261 + x263 + x264;
auto x271 = -V{0.125}*cell[15] + V{0.375}*cell[16] - V{0.125}*cell[6] + V{0.375}*cell[7] + x266;
auto x272 = cell[11] + x74 + x79;
auto x273 = x187*x38*((cell[12] + 2*cell[17] - 2*cell[8] + x272 + x40 + x46 + x98)*(cell[12] + 2*cell[17] - 2*cell[8] + x272 + x40 + x46 + x98));
auto x274 = x236*x81;
auto x275 = x104*x238;
auto x276 = x180 + x274 - x275;
auto x277 = x103*x241*x48;
auto x278 = x243*x48*x80;
auto x279 = -x172;
auto x280 = -x173;
auto x281 = x127 + x170;
auto x282 = x277 + x278 + x279 + x280 + x281;
auto x283 = -x162;
auto x284 = -x163;
auto x285 = x105*x195;
auto x286 = x214*x66;
auto x287 = x118*x252;
auto x288 = x107*x218;
auto x289 = x108*x206;
auto x290 = x110*x289;
auto x291 = x209 + x212 + x257 + x283 + x284 + x285 + x286 + x287 - x288 - x290;
auto x292 = x259*x91;
auto x293 = x118*x201;
auto x294 = x107*x262;
auto x295 = x113*x289;
auto x296 = x292 + x293 - x294 - x295;
auto x297 = V{0.25}*x115*x38;
auto x298 = V{0.375}*cell[17] - V{0.125}*cell[18] + V{0.375}*cell[8] - V{0.125}*cell[9] - x297;
auto x299 = V{2}*cell[18] + cell[6] - V{2}*cell[9] + x100 + x112 + x268 + x272;
auto x300 = -x292 - x293 + x294 + x295;
auto x301 = -V{0.125}*cell[17] + V{0.375}*cell[18] - V{0.125}*cell[8] + V{0.375}*cell[9] + x297;
auto x302 = -x51;
auto x303 = x302 - x52 - x53 + x54 + x55 + x56 + x57;
auto x304 = -x84 - x85 + x86 + x87 + x88 + x89;
auto x305 = x303 + x304 - x67 - x68 - x69 - x70 + x71 + x72 - x83;
auto x306 = -x59 - x60 + x61 + x62 + x63 + x64;
auto x307 = x303 + x306 - x31 - x32 - x33 - x34 + x35 + x36 - x50;
auto x308 = -x110;
auto x309 = -x113;
auto x310 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x311 = x123 + V{-1};
auto x312 = -x106 + x302 + x304 + x306 - x92 - x93 - x94 - x95 + x96 + x97;
auto x313 = -x118;
auto x314 = x193 + x194 + x196 - x199 - x202 + x205 + x208 + x211 + x212;
auto x315 = x146 + x185;
auto x316 = x154 + x315;
auto x317 = x211 + x248 + x249 + x250 - x251 - x253 + x254 + x256 + x257;
auto x318 = x209 + x212 + x257 + x283 + x284 + x285 - x286 - x287 + x288 + x290;
auto x0 = -V{1}*x119*(-V{0.5}*cell[0] + V{4.16333634234434e-17}*cell[10] + V{4.16333634234434e-17}*cell[11] + V{4.16333634234434e-17}*cell[12] + V{0.5}*cell[13] + V{0.5}*cell[14] + V{0.5}*cell[15] + V{0.5}*cell[16] + V{0.5}*cell[17] + V{0.5}*cell[18] + V{4.16333634234434e-17}*cell[1] + V{4.16333634234434e-17}*cell[2] + V{4.16333634234434e-17}*cell[3] + V{0.5}*cell[4] + V{0.5}*cell[5] + V{0.5}*cell[6] + V{0.5}*cell[7] + V{0.5}*cell[8] + V{0.5}*cell[9] - x105*x120 - x120*x49 - x120*x82) - x127*(x31 + x32 + x33 + x34 + x51 + x52 + x53 + x59 + x60 + x67 + x68 + x69 + x70 + x84 + x85 + x92 + x93 + x94 + x95 + V{0.333333333333333}) + V{-0.333333333333333};
auto x1 = -V{1}*x119*(x110*x131 + x113*x132 + x129*x66 + x129*x91 + x142) + V{0.0555555555555556}*x143*(x145 + x146 + x147 + x159) + V{-0.0555555555555556};
auto x2 = -V{1}*x119*(x107*x160 + x110*x161 + x118*x132 + x160*x66 + x165) + V{0.0555555555555556}*x143*(x146 + x166 + x167 + x175 + V{1}) + V{-0.0555555555555556};
auto x3 = -V{1}*x119*(x107*x176 + x113*x161 + x118*x131 + x176*x91 + x177) + V{0.0555555555555556}*x143*(x178 + x184 + x185) + V{-0.0555555555555556};
auto x4 = V{1}*x119*(x213 + x221 + x223) - x186*(-x189 + x190 + x192) + V{-0.0277777777777778};
auto x5 = -(-V{1}*x119*(x213 + x231 + x232) + x186*(x171 - x187*x38*((x224)*(x224)) + x190 + x226 + x228 - x229 - x230) + V{0.0277777777777778});
auto x6 = V{1}*x119*(x258 + x265 + x267) - x186*(-x234 + x240 + x247) + V{-0.0277777777777778};
auto x7 = -(-V{1}*x119*(x258 + x270 + x271) + x186*(x181 - x187*x38*((x269)*(x269)) - x237 + x239 + x247) + V{0.0277777777777778});
auto x8 = V{1}*x119*(x291 + x296 + x298) - x186*(-x273 + x276 + x282) + V{-0.0277777777777778};
auto x9 = -(-V{1}*x119*(x291 + x300 + x301) + x186*(x181 - x187*x38*((x299)*(x299)) - x274 + x275 + x282) + V{0.0277777777777778});
auto x10 = -V{1}*x119*(x129*x305 + x129*x307 + x131*x308 + x132*x309 + x142) - x310*(x125 - x145 + x159 + x311) + V{-0.0555555555555556};
auto x11 = -V{1}*x119*(x132*x313 + x160*x307 + x160*x312 + x161*x308 + x165) - x310*(x124 - x167 + x175 + x311) + V{-0.0555555555555556};
auto x12 = -V{1}*x119*(x131*x313 + x161*x309 + x176*x305 + x176*x312 + x177) - x310*(x126 - x178 + x184) + V{-0.0555555555555556};
auto x13 = V{1}*x119*(x223 + x231 + x314) + x186*(x189 + x192 + x316) + V{-0.0277777777777778};
auto x14 = -(-V{1}*x119*(x221 + x232 + x314) + x186*(x155 - x187*x38*((x224)*(x224)) - x226 - x228 + x229 + x230 + x281) + V{0.0277777777777778});
auto x15 = V{1}*x119*(x267 + x270 + x317) + x186*(x234 + x240 + x242 + x244 + x245 + x246 + x316) + V{-0.0277777777777778};
auto x16 = -(-V{1}*x119*(x265 + x271 + x317) + x186*(x127 + x158 - x187*x38*((x269)*(x269)) + x240 - x242 - x244) + V{0.0277777777777778});
auto x17 = V{1}*x119*(x298 + x300 + x318) + x186*(x170 + x273 + x276 + x277 + x278 + x279 + x280 + x315) + V{-0.0277777777777778};
auto x18 = -(-V{1}*x119*(x296 + x301 + x318) + x186*(x127 + x174 - x187*x38*((x299)*(x299)) + x276 - x277 - x278) + V{0.0277777777777778});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
cell[9] = x9;
cell[10] = x10;
cell[11] = x11;
cell[12] = x12;
cell[13] = x13;
cell[14] = x14;
cell[15] = x15;
cell[16] = x16;
cell[17] = x17;
cell[18] = x18;
return { x29, V{1}*x121*x38*(x104 + x48 + x81) };
}
};

}

}

// Generation Info: commit=d03b8c680a234227b40a2e0c9e29e570326b7e51
