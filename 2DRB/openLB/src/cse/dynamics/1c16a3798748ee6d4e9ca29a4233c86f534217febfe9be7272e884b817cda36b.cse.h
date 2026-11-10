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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::ThirdOrder, collision::FrozenOmegaFromCell<collision::ThirdOrderRLB>, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x21 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x20 = cell.template getFieldComponent<olb::descriptors::OMEGA>(0);
auto x19 = cell.template getFieldComponent<olb::descriptors::AVERAGE_TURBULENT_TAU>(0);
auto x22 = V{1} / (x19 + V{1}/x20);
auto x23 = V{1} - V{1}*x22;
auto x24 = cell[13] + cell[15];
auto x25 = cell[10] + x24;
auto x26 = cell[11] + cell[5];
auto x27 = cell[12] + cell[7] + cell[9];
auto x28 = cell[2] + cell[8];
auto x29 = cell[18] + cell[3];
auto x30 = cell[16] + cell[6];
auto x31 = cell[0] + cell[14] + cell[17] + cell[1] + cell[4] + x25 + x26 + x27 + x28 + x29 + x30;
auto x32 = x31 + V{1};
auto x33 = V{1} / (x32);
auto x34 = V{0.5}*x33;
auto x35 = ((x21)*(x21));
auto x36 = -cell[1];
auto x37 = -cell[7];
auto x38 = -cell[6];
auto x39 = cell[16] + x38;
auto x40 = x36 + x37 + x39;
auto x41 = -cell[4];
auto x42 = cell[14] + x41;
auto x43 = -cell[5] + x42;
auto x44 = x25 + x40 + x43;
auto x45 = ((x44)*(x44));
auto x46 = x35*x45;
auto x47 = cell[13] + cell[17];
auto x48 = -cell[2];
auto x49 = -cell[9];
auto x50 = -cell[8];
auto x51 = cell[18] + x50;
auto x52 = x48 + x49 + x51;
auto x53 = -cell[14] + x26 + x41;
auto x54 = x47 + x52 + x53;
auto x55 = ((x54)*(x54));
auto x56 = x35*x55;
auto x57 = cell[15] + cell[17];
auto x58 = -cell[18];
auto x59 = -cell[3];
auto x60 = x50 + x58 + x59;
auto x61 = -cell[16] + x38;
auto x62 = x27 + x57 + x60 + x61;
auto x63 = ((x62)*(x62));
auto x64 = x35*x63;
auto x65 = V{0.333333333333333}*cell[0];
auto x66 = V{0.333333333333333}*cell[10];
auto x67 = V{0.333333333333333}*cell[11];
auto x68 = V{0.333333333333333}*cell[12];
auto x69 = V{0.333333333333333}*cell[13];
auto x70 = V{0.333333333333333}*cell[14];
auto x71 = V{0.333333333333333}*cell[15];
auto x72 = V{0.333333333333333}*cell[16];
auto x73 = V{0.333333333333333}*cell[17];
auto x74 = V{0.333333333333333}*cell[18];
auto x75 = V{0.333333333333333}*cell[1];
auto x76 = V{0.333333333333333}*cell[2];
auto x77 = V{0.333333333333333}*cell[3];
auto x78 = V{0.333333333333333}*cell[4];
auto x79 = V{0.333333333333333}*cell[5];
auto x80 = V{0.333333333333333}*cell[6];
auto x81 = V{0.333333333333333}*cell[7];
auto x82 = V{0.333333333333333}*cell[8];
auto x83 = V{0.333333333333333}*cell[9];
auto x84 = V{1} / ((x32)*(x32));
auto x85 = V{1.5}*x84;
auto x86 = x64*x85;
auto x87 = x46*x85;
auto x88 = x56*x85;
auto x89 = x87 + x88 + V{-1};
auto x90 = x86 + x89;
auto x91 = x21*x34;
auto x92 = x44*x91;
auto x93 = V{0.666666666666667}*cell[12];
auto x94 = V{0.666666666666667}*cell[3];
auto x95 = V{1}*x33;
auto x96 = x64*x95;
auto x97 = V{0.666666666666667}*cell[17];
auto x98 = V{0.666666666666667}*cell[18];
auto x99 = V{0.666666666666667}*cell[8];
auto x100 = V{0.666666666666667}*cell[9];
auto x101 = -x100 + x65 + x66 + x75 - x97 - x98 - x99;
auto x102 = V{0.666666666666667}*cell[15];
auto x103 = V{0.666666666666667}*cell[16];
auto x104 = V{0.666666666666667}*cell[6];
auto x105 = V{0.666666666666667}*cell[7];
auto x106 = -x102 - x103 - x104 - x105 + x67 + x76;
auto x107 = x101 + x106 + x69 + x70 + x78 + x79 - x93 - x94 + x96;
auto x108 = V{0.666666666666667}*cell[11];
auto x109 = V{0.666666666666667}*cell[2];
auto x110 = x56*x95;
auto x111 = V{0.666666666666667}*cell[13];
auto x112 = V{0.666666666666667}*cell[14];
auto x113 = V{0.666666666666667}*cell[4];
auto x114 = V{0.666666666666667}*cell[5];
auto x115 = -x111 - x112 - x113 - x114 + x68 + x77;
auto x116 = x101 - x108 - x109 + x110 + x115 + x71 + x72 + x80 + x81;
auto x117 = x21*x95;
auto x118 = x117*x54;
auto x119 = x33*x44;
auto x120 = x119*x54;
auto x121 = -cell[13] + cell[5] + x120*x35 + x42;
auto x122 = x117*x62;
auto x123 = -cell[15];
auto x124 = x119*x62;
auto x125 = cell[7] + x123 + x124*x35 + x39;
auto x126 = V{0.166666666666667}*x33;
auto x127 = V{6.93889390390723e-18}*cell[0];
auto x128 = V{0.0833333333333333}*cell[12];
auto x129 = V{0.0833333333333333}*cell[3];
auto x130 = V{0.0833333333333333}*x33;
auto x131 = -V{0.0833333333333333}*cell[13] - V{0.0833333333333333}*cell[14] - V{0.0833333333333333}*cell[4] - V{0.0833333333333333}*cell[5] + x127 + x128 + x129 - x130*x64;
auto x132 = V{0.0833333333333333}*cell[11];
auto x133 = V{0.0833333333333333}*cell[2];
auto x134 = -V{0.0833333333333333}*cell[15] - V{0.0833333333333333}*cell[16] - V{0.0833333333333333}*cell[6] - V{0.0833333333333333}*cell[7] - x130*x56 + x132 + x133;
auto x135 = -V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] - V{0.166666666666667}*cell[1] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x126*x46 + x131 + x134;
auto x136 = x31 + V{1};
auto x137 = V{3}*x84;
auto x138 = x137*x46;
auto x139 = -x86;
auto x140 = V{1} - x88;
auto x141 = util::pow(x21, 3)/util::pow(x32, 3);
auto x142 = V{6.000003}*x141;
auto x143 = x142*x44;
auto x144 = V{2.999997}*x141;
auto x145 = x144*x54;
auto x146 = V{3}*x21;
auto x147 = x119*x146;
auto x148 = -x147;
auto x149 = x143*x55;
auto x150 = x145*x45;
auto x151 = x148 + x149 + x150;
auto x152 = x143*x63 + x145*x63 + x151;
auto x153 = x54*x91;
auto x154 = V{0.666666666666667}*cell[10];
auto x155 = V{0.666666666666667}*cell[1];
auto x156 = x46*x95;
auto x157 = x106 + x115 - x154 - x155 + x156 + x65 + x73 + x74 + x82 + x83;
auto x158 = x117*x44;
auto x159 = x33*x54;
auto x160 = x159*x62;
auto x161 = -cell[17];
auto x162 = cell[9] + x161;
auto x163 = x160*x35 + x162 + x51;
auto x164 = V{0.0833333333333333}*cell[10];
auto x165 = V{0.0833333333333333}*cell[1];
auto x166 = -V{0.0833333333333333}*cell[17] - V{0.0833333333333333}*cell[18] - V{0.0833333333333333}*cell[8] - V{0.0833333333333333}*cell[9] - x130*x46 + x164 + x165;
auto x167 = -V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] - V{0.166666666666667}*cell[2] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] + x126*x56 + x131 + x166;
auto x168 = -x87;
auto x169 = x137*x56;
auto x170 = x142*x54;
auto x171 = x144*x44;
auto x172 = x146*x159;
auto x173 = -x172;
auto x174 = x170*x45;
auto x175 = x171*x55;
auto x176 = x173 + x174 + x175;
auto x177 = x170*x63 + x171*x63 + x176;
auto x178 = x62*x91;
auto x179 = -V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] - V{0.166666666666667}*cell[3] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] + x126*x64 + x127 + x134 + x166;
auto x180 = x137*x64;
auto x181 = x33*x62;
auto x182 = x146*x181;
auto x183 = -x182;
auto x184 = V{9}*x141;
auto x185 = x184*x62;
auto x186 = x183 + x185*x45 + x185*x55;
auto x187 = x140 + x168;
auto x188 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x189 = V{4.5}*x84;
auto x190 = cell[10] + x40;
auto x191 = x189*x35*((cell[11] + 2*cell[13] - 2*cell[4] + x190 + x52 + x57)*(cell[11] + 2*cell[13] - 2*cell[4] + x190 + x52 + x57));
auto x192 = x147 + x90;
auto x193 = V{18}*x141;
auto x194 = x172 - x184*x44*x63 - x184*x54*x63 + x193*x44*x55 + x193*x45*x54;
auto x195 = V{1} - x22;
auto x196 = -x128;
auto x197 = -x129;
auto x198 = V{0.0416666666666667}*x33;
auto x199 = x198*x64;
auto x200 = V{0.25000025}*x21;
auto x201 = x119*x200;
auto x202 = x116*x201;
auto x203 = V{0.5000005}*x21;
auto x204 = x159*x203;
auto x205 = x121*x204;
auto x206 = V{2.49999999985601e-07}*x21;
auto x207 = x119*x206;
auto x208 = x107*x207;
auto x209 = V{4.99999999971202e-07}*x21;
auto x210 = x181*x209;
auto x211 = x125*x210;
auto x212 = -V{0.0416666666666667}*cell[0];
auto x213 = V{0.0833333333333333}*x33;
auto x214 = V{0.0416666666666667}*cell[10] + V{6.93889390390723e-18}*cell[17] + V{6.93889390390723e-18}*cell[18] + V{0.0416666666666667}*cell[1] + V{6.93889390390723e-18}*cell[8] + V{6.93889390390723e-18}*cell[9] + x212 - x213*x46;
auto x215 = V{0.0416666666666667}*cell[11] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[16] + V{0.0416666666666667}*cell[2] + V{6.93889390390723e-18}*cell[6] + V{6.93889390390723e-18}*cell[7] - x213*x56;
auto x216 = x196 + x197 + x199 + x202 + x205 - x208 - x211 + x214 + x215;
auto x217 = x159*x200;
auto x218 = x157*x217;
auto x219 = x119*x203;
auto x220 = x121*x219;
auto x221 = x159*x206;
auto x222 = x107*x221;
auto x223 = x163*x210;
auto x224 = x218 + x220 - x222 - x223;
auto x225 = V{0.25}*x120*x35;
auto x226 = V{0.375}*cell[13] - V{0.125}*cell[14] + V{0.375}*cell[4] - V{0.125}*cell[5] - x225;
auto x227 = -cell[11] + V{2}*cell[14] + cell[15] - V{2}*cell[5] + x162 + x190 + x28 + x58;
auto x228 = V{3.000006}*x141;
auto x229 = x228*x54*x63;
auto x230 = V{6.000012}*x141;
auto x231 = x230*x44*x55;
auto x232 = x228*x44*x63;
auto x233 = x230*x45*x54;
auto x234 = -x218 - x220 + x222 + x223;
auto x235 = -V{0.125}*cell[13] + V{0.375}*cell[14] - V{0.125}*cell[4] + V{0.375}*cell[5] + x225;
auto x236 = cell[10] + x36 + x43;
auto x237 = x189*x35*((cell[12] + 2*cell[15] - 2*cell[6] + cell[9] + x236 + x47 + x60)*(cell[12] + 2*cell[15] - 2*cell[6] + cell[9] + x236 + x47 + x60));
auto x238 = x141*x62;
auto x239 = V{9.000009}*x238;
auto x240 = x239*x45;
auto x241 = V{8.99999999948164e-06}*x238;
auto x242 = x241*x55;
auto x243 = x182 + x240 - x242;
auto x244 = V{5.999994}*x141;
auto x245 = x244*x54*x63;
auto x246 = V{12.000006}*x141;
auto x247 = x246*x44*x63;
auto x248 = -x149;
auto x249 = -x150;
auto x250 = x192 + x245 + x247 + x248 + x249;
auto x251 = -x132;
auto x252 = -x133;
auto x253 = x198*x56;
auto x254 = x107*x201;
auto x255 = x181*x203;
auto x256 = x125*x255;
auto x257 = x116*x207;
auto x258 = x159*x209;
auto x259 = x121*x258;
auto x260 = V{0.0416666666666667}*cell[12] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[14] + V{0.0416666666666667}*cell[3] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[5] - x213*x64;
auto x261 = x214 + x251 + x252 + x253 + x254 + x256 - x257 - x259 + x260;
auto x262 = x181*x200;
auto x263 = x157*x262;
auto x264 = x125*x219;
auto x265 = x181*x206;
auto x266 = x116*x265;
auto x267 = x163*x258;
auto x268 = x263 + x264 - x266 - x267;
auto x269 = V{0.25}*x124*x35;
auto x270 = V{0.375}*cell[15] - V{0.125}*cell[16] + V{0.375}*cell[6] - V{0.125}*cell[7] - x269;
auto x271 = -cell[12] + cell[13];
auto x272 = V{2}*cell[16] - V{2}*cell[7] + cell[8] + x161 + x236 + x271 + x29 + x49;
auto x273 = -x263 - x264 + x266 + x267;
auto x274 = -V{0.125}*cell[15] + V{0.375}*cell[16] - V{0.125}*cell[6] + V{0.375}*cell[7] + x269;
auto x275 = x48 + x53;
auto x276 = x189*x35*((cell[12] + 2*cell[17] + cell[7] - 2*cell[8] + x24 + x275 + x59 + x61)*(cell[12] + 2*cell[17] + cell[7] - 2*cell[8] + x24 + x275 + x59 + x61));
auto x277 = x239*x55;
auto x278 = x241*x45;
auto x279 = x182 + x277 - x278;
auto x280 = x244*x44*x63;
auto x281 = x246*x54*x63;
auto x282 = -x174;
auto x283 = -x175;
auto x284 = x172 + x90;
auto x285 = x280 + x281 + x282 + x283 + x284;
auto x286 = -x164;
auto x287 = -x165;
auto x288 = x198*x46;
auto x289 = x107*x217;
auto x290 = x163*x255;
auto x291 = x157*x221;
auto x292 = x119*x209;
auto x293 = x121*x292;
auto x294 = x212 + x215 + x260 + x286 + x287 + x288 + x289 + x290 - x291 - x293;
auto x295 = x116*x262;
auto x296 = x163*x204;
auto x297 = x157*x265;
auto x298 = x125*x292;
auto x299 = x295 + x296 - x297 - x298;
auto x300 = V{0.25}*x160*x35;
auto x301 = V{0.375}*cell[17] - V{0.125}*cell[18] + V{0.375}*cell[8] - V{0.125}*cell[9] - x300;
auto x302 = V{2}*cell[18] + cell[3] - V{2}*cell[9] + x123 + x271 + x275 + x30 + x37;
auto x303 = -x295 - x296 + x297 + x298;
auto x304 = -V{0.125}*cell[17] + V{0.375}*cell[18] - V{0.125}*cell[8] + V{0.375}*cell[9] + x300;
auto x305 = -x65;
auto x306 = x100 + x305 - x66 - x75 + x97 + x98 + x99;
auto x307 = x111 + x112 + x113 + x114 - x68 - x77;
auto x308 = x108 + x109 - x110 + x306 + x307 - x71 - x72 - x80 - x81;
auto x309 = x102 + x103 + x104 + x105 - x67 - x76;
auto x310 = x306 + x309 - x69 - x70 - x78 - x79 + x93 + x94 - x96;
auto x311 = -x121;
auto x312 = -x125;
auto x313 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x314 = x86 + V{-1};
auto x315 = x154 + x155 - x156 + x305 + x307 + x309 - x73 - x74 - x82 - x83;
auto x316 = -x163;
auto x317 = x196 + x197 + x199 - x202 - x205 + x208 + x211 + x214 + x215;
auto x318 = x139 + x187;
auto x319 = x147 + x318;
auto x320 = x214 + x251 + x252 + x253 - x254 - x256 + x257 + x259 + x260;
auto x321 = x212 + x215 + x260 + x286 + x287 + x288 - x289 - x290 + x291 + x293;
auto x0 = -x23*(-V{0.5}*cell[0] + V{4.16333634234434e-17}*cell[10] + V{4.16333634234434e-17}*cell[11] + V{4.16333634234434e-17}*cell[12] + V{0.5}*cell[13] + V{0.5}*cell[14] + V{0.5}*cell[15] + V{0.5}*cell[16] + V{0.5}*cell[17] + V{0.5}*cell[18] + V{4.16333634234434e-17}*cell[1] + V{4.16333634234434e-17}*cell[2] + V{4.16333634234434e-17}*cell[3] + V{0.5}*cell[4] + V{0.5}*cell[5] + V{0.5}*cell[6] + V{0.5}*cell[7] + V{0.5}*cell[8] + V{0.5}*cell[9] - x34*x46 - x34*x56 - x34*x64) - x90*(x65 + x66 + x67 + x68 + x69 + x70 + x71 + x72 + x73 + x74 + x75 + x76 + x77 + x78 + x79 + x80 + x81 + x82 + x83 + V{0.333333333333333}) + V{-0.333333333333333};
auto x1 = V{0.0555555555555556}*x136*(x138 + x139 + x140 + x152) - x23*(x107*x92 + x116*x92 + x118*x121 + x122*x125 + x135) + V{-0.0555555555555556};
auto x2 = V{0.0555555555555556}*x136*(x139 + x168 + x169 + x177 + V{1}) - x23*(x107*x153 + x121*x158 + x122*x163 + x153*x157 + x167) + V{-0.0555555555555556};
auto x3 = V{0.0555555555555556}*x136*(x180 + x186 + x187) - x23*(x116*x178 + x118*x163 + x125*x158 + x157*x178 + x179) + V{-0.0555555555555556};
auto x4 = -x188*(-x191 + x192 + x194) + V{1}*x195*(x216 + x224 + x226) + V{-0.0277777777777778};
auto x5 = -(x188*(x173 - x189*x35*((x227)*(x227)) + x192 + x229 + x231 - x232 - x233) - V{1}*x195*(x216 + x234 + x235) + V{0.0277777777777778});
auto x6 = -x188*(-x237 + x243 + x250) + V{1}*x195*(x261 + x268 + x270) + V{-0.0277777777777778};
auto x7 = -(x188*(x183 - x189*x35*((x272)*(x272)) - x240 + x242 + x250) - V{1}*x195*(x261 + x273 + x274) + V{0.0277777777777778});
auto x8 = -x188*(-x276 + x279 + x285) + V{1}*x195*(x294 + x299 + x301) + V{-0.0277777777777778};
auto x9 = -(x188*(x183 - x189*x35*((x302)*(x302)) - x277 + x278 + x285) - V{1}*x195*(x294 + x303 + x304) + V{0.0277777777777778});
auto x10 = -x23*(x118*x311 + x122*x312 + x135 + x308*x92 + x310*x92) - x313*(-x138 + x152 + x314 + x88) + V{-0.0555555555555556};
auto x11 = -x23*(x122*x316 + x153*x310 + x153*x315 + x158*x311 + x167) - x313*(-x169 + x177 + x314 + x87) + V{-0.0555555555555556};
auto x12 = -x23*(x118*x316 + x158*x312 + x178*x308 + x178*x315 + x179) - x313*(-x180 + x186 + x89) + V{-0.0555555555555556};
auto x13 = x188*(x191 + x194 + x319) + x23*(x226 + x234 + x317) + V{-0.0277777777777778};
auto x14 = -(x188*(x148 - x189*x35*((x227)*(x227)) - x229 - x231 + x232 + x233 + x284) - V{1}*x195*(x224 + x235 + x317) + V{0.0277777777777778});
auto x15 = x188*(x237 + x243 + x245 + x247 + x248 + x249 + x319) + x23*(x270 + x273 + x320) + V{-0.0277777777777778};
auto x16 = -(x188*(x151 - x189*x35*((x272)*(x272)) + x243 - x245 - x247 + x90) - V{1}*x195*(x268 + x274 + x320) + V{0.0277777777777778});
auto x17 = x188*(x172 + x276 + x279 + x280 + x281 + x282 + x283 + x318) + x23*(x301 + x303 + x321) + V{-0.0277777777777778};
auto x18 = -(x188*(x176 - x189*x35*((x302)*(x302)) + x279 - x280 - x281 + x90) - V{1}*x195*(x299 + x304 + x321) + V{0.0277777777777778});
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
return { x32, V{1}*x35*x84*(x45 + x55 + x63) };
}
};

}

}

// Generation Info: commit=b53c7b1cd5bd508c995f735278effcb428fdd534
