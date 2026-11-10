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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::ThirdOrder, collision::ParameterFromCell<descriptors::OMEGA, collision::ThirdOrderRLB>, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = cell.template getFieldComponent<olb::descriptors::OMEGA>(0);
auto x20 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x21 = x19 + V{-1};
auto x22 = cell[13] + cell[15];
auto x23 = cell[10] + x22;
auto x24 = cell[11] + cell[5];
auto x25 = cell[12] + cell[7] + cell[9];
auto x26 = cell[2] + cell[8];
auto x27 = cell[18] + cell[3];
auto x28 = cell[16] + cell[6];
auto x29 = cell[0] + cell[14] + cell[17] + cell[1] + cell[4] + x23 + x24 + x25 + x26 + x27 + x28;
auto x30 = x29 + V{1};
auto x31 = V{1} / (x30);
auto x32 = V{0.5}*x31;
auto x33 = ((x20)*(x20));
auto x34 = -cell[1];
auto x35 = -cell[7];
auto x36 = -cell[6];
auto x37 = cell[16] + x36;
auto x38 = x34 + x35 + x37;
auto x39 = -cell[4];
auto x40 = cell[14] + x39;
auto x41 = -cell[5] + x40;
auto x42 = x23 + x38 + x41;
auto x43 = ((x42)*(x42));
auto x44 = x33*x43;
auto x45 = cell[13] + cell[17];
auto x46 = -cell[2];
auto x47 = -cell[9];
auto x48 = -cell[8];
auto x49 = cell[18] + x48;
auto x50 = x46 + x47 + x49;
auto x51 = -cell[14] + x24 + x39;
auto x52 = x45 + x50 + x51;
auto x53 = ((x52)*(x52));
auto x54 = x33*x53;
auto x55 = cell[15] + cell[17];
auto x56 = -cell[18];
auto x57 = -cell[3];
auto x58 = x48 + x56 + x57;
auto x59 = -cell[16] + x36;
auto x60 = x25 + x55 + x58 + x59;
auto x61 = ((x60)*(x60));
auto x62 = x33*x61;
auto x63 = V{0.333333333333333}*cell[0];
auto x64 = V{0.333333333333333}*cell[10];
auto x65 = V{0.333333333333333}*cell[11];
auto x66 = V{0.333333333333333}*cell[12];
auto x67 = V{0.333333333333333}*cell[13];
auto x68 = V{0.333333333333333}*cell[14];
auto x69 = V{0.333333333333333}*cell[15];
auto x70 = V{0.333333333333333}*cell[16];
auto x71 = V{0.333333333333333}*cell[17];
auto x72 = V{0.333333333333333}*cell[18];
auto x73 = V{0.333333333333333}*cell[1];
auto x74 = V{0.333333333333333}*cell[2];
auto x75 = V{0.333333333333333}*cell[3];
auto x76 = V{0.333333333333333}*cell[4];
auto x77 = V{0.333333333333333}*cell[5];
auto x78 = V{0.333333333333333}*cell[6];
auto x79 = V{0.333333333333333}*cell[7];
auto x80 = V{0.333333333333333}*cell[8];
auto x81 = V{0.333333333333333}*cell[9];
auto x82 = V{1} / ((x30)*(x30));
auto x83 = V{1.5}*x82;
auto x84 = x62*x83;
auto x85 = x44*x83;
auto x86 = x54*x83;
auto x87 = x85 + x86 + V{-1};
auto x88 = x84 + x87;
auto x89 = x20*x32;
auto x90 = x42*x89;
auto x91 = V{0.666666666666667}*cell[12];
auto x92 = V{0.666666666666667}*cell[3];
auto x93 = V{1}*x31;
auto x94 = x62*x93;
auto x95 = V{0.666666666666667}*cell[17];
auto x96 = V{0.666666666666667}*cell[18];
auto x97 = V{0.666666666666667}*cell[8];
auto x98 = V{0.666666666666667}*cell[9];
auto x99 = x63 + x64 + x73 - x95 - x96 - x97 - x98;
auto x100 = V{0.666666666666667}*cell[15];
auto x101 = V{0.666666666666667}*cell[16];
auto x102 = V{0.666666666666667}*cell[6];
auto x103 = V{0.666666666666667}*cell[7];
auto x104 = -x100 - x101 - x102 - x103 + x65 + x74;
auto x105 = x104 + x67 + x68 + x76 + x77 - x91 - x92 + x94 + x99;
auto x106 = V{0.666666666666667}*cell[11];
auto x107 = V{0.666666666666667}*cell[2];
auto x108 = x54*x93;
auto x109 = V{0.666666666666667}*cell[13];
auto x110 = V{0.666666666666667}*cell[14];
auto x111 = V{0.666666666666667}*cell[4];
auto x112 = V{0.666666666666667}*cell[5];
auto x113 = -x109 - x110 - x111 - x112 + x66 + x75;
auto x114 = -x106 - x107 + x108 + x113 + x69 + x70 + x78 + x79 + x99;
auto x115 = x20*x93;
auto x116 = x115*x52;
auto x117 = x31*x42;
auto x118 = x117*x52;
auto x119 = -cell[13] + cell[5] + x118*x33 + x40;
auto x120 = x115*x60;
auto x121 = -cell[15];
auto x122 = x117*x60;
auto x123 = cell[7] + x121 + x122*x33 + x37;
auto x124 = V{0.166666666666667}*x31;
auto x125 = V{6.93889390390723e-18}*cell[0];
auto x126 = V{0.0833333333333333}*cell[12];
auto x127 = V{0.0833333333333333}*cell[3];
auto x128 = V{0.0833333333333333}*x31;
auto x129 = -V{0.0833333333333333}*cell[13] - V{0.0833333333333333}*cell[14] - V{0.0833333333333333}*cell[4] - V{0.0833333333333333}*cell[5] + x125 + x126 + x127 - x128*x62;
auto x130 = V{0.0833333333333333}*cell[11];
auto x131 = V{0.0833333333333333}*cell[2];
auto x132 = -V{0.0833333333333333}*cell[15] - V{0.0833333333333333}*cell[16] - V{0.0833333333333333}*cell[6] - V{0.0833333333333333}*cell[7] - x128*x54 + x130 + x131;
auto x133 = -V{0.166666666666667}*cell[10] + V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] - V{0.166666666666667}*cell[1] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x124*x44 + x129 + x132;
auto x134 = x29 + V{1};
auto x135 = V{3}*x82;
auto x136 = x135*x44;
auto x137 = -x84;
auto x138 = V{1} - x86;
auto x139 = util::pow(x20, 3)/util::pow(x30, 3);
auto x140 = V{6.000003}*x139;
auto x141 = x140*x42;
auto x142 = V{2.999997}*x139;
auto x143 = x142*x52;
auto x144 = V{3}*x20;
auto x145 = x117*x144;
auto x146 = -x145;
auto x147 = x141*x53;
auto x148 = x143*x43;
auto x149 = x146 + x147 + x148;
auto x150 = x141*x61 + x143*x61 + x149;
auto x151 = x52*x89;
auto x152 = V{0.666666666666667}*cell[10];
auto x153 = V{0.666666666666667}*cell[1];
auto x154 = x44*x93;
auto x155 = x104 + x113 - x152 - x153 + x154 + x63 + x71 + x72 + x80 + x81;
auto x156 = x115*x42;
auto x157 = x31*x52;
auto x158 = x157*x60;
auto x159 = -cell[17];
auto x160 = cell[9] + x159;
auto x161 = x158*x33 + x160 + x49;
auto x162 = V{0.0833333333333333}*cell[10];
auto x163 = V{0.0833333333333333}*cell[1];
auto x164 = -V{0.0833333333333333}*cell[17] - V{0.0833333333333333}*cell[18] - V{0.0833333333333333}*cell[8] - V{0.0833333333333333}*cell[9] - x128*x44 + x162 + x163;
auto x165 = -V{0.166666666666667}*cell[11] + V{0.166666666666667}*cell[15] + V{0.166666666666667}*cell[16] - V{0.166666666666667}*cell[2] + V{0.166666666666667}*cell[6] + V{0.166666666666667}*cell[7] + x124*x54 + x129 + x164;
auto x166 = -x85;
auto x167 = x135*x54;
auto x168 = x140*x52;
auto x169 = x142*x42;
auto x170 = x144*x157;
auto x171 = -x170;
auto x172 = x168*x43;
auto x173 = x169*x53;
auto x174 = x171 + x172 + x173;
auto x175 = x168*x61 + x169*x61 + x174;
auto x176 = x60*x89;
auto x177 = -V{0.166666666666667}*cell[12] + V{0.166666666666667}*cell[13] + V{0.166666666666667}*cell[14] - V{0.166666666666667}*cell[3] + V{0.166666666666667}*cell[4] + V{0.166666666666667}*cell[5] + x124*x62 + x125 + x132 + x164;
auto x178 = x135*x62;
auto x179 = x31*x60;
auto x180 = x144*x179;
auto x181 = -x180;
auto x182 = V{9}*x139;
auto x183 = x182*x60;
auto x184 = x181 + x183*x43 + x183*x53;
auto x185 = x138 + x166;
auto x186 = -x126;
auto x187 = -x127;
auto x188 = V{0.0416666666666667}*x31;
auto x189 = x188*x62;
auto x190 = V{0.25000025}*x20;
auto x191 = x117*x190;
auto x192 = x114*x191;
auto x193 = V{0.5000005}*x20;
auto x194 = x157*x193;
auto x195 = x119*x194;
auto x196 = V{2.49999999985601e-07}*x20;
auto x197 = x117*x196;
auto x198 = x105*x197;
auto x199 = V{4.99999999971202e-07}*x20;
auto x200 = x179*x199;
auto x201 = x123*x200;
auto x202 = -V{0.0416666666666667}*cell[0];
auto x203 = V{0.0833333333333333}*x31;
auto x204 = V{0.0416666666666667}*cell[10] + V{6.93889390390723e-18}*cell[17] + V{6.93889390390723e-18}*cell[18] + V{0.0416666666666667}*cell[1] + V{6.93889390390723e-18}*cell[8] + V{6.93889390390723e-18}*cell[9] + x202 - x203*x44;
auto x205 = V{0.0416666666666667}*cell[11] + V{6.93889390390723e-18}*cell[15] + V{6.93889390390723e-18}*cell[16] + V{0.0416666666666667}*cell[2] + V{6.93889390390723e-18}*cell[6] + V{6.93889390390723e-18}*cell[7] - x203*x54;
auto x206 = x186 + x187 + x189 + x192 + x195 - x198 - x201 + x204 + x205;
auto x207 = x157*x190;
auto x208 = x155*x207;
auto x209 = x117*x193;
auto x210 = x119*x209;
auto x211 = x157*x196;
auto x212 = x105*x211;
auto x213 = x161*x200;
auto x214 = x208 + x210 - x212 - x213;
auto x215 = V{0.25}*x118*x33;
auto x216 = V{0.375}*cell[13] - V{0.125}*cell[14] + V{0.375}*cell[4] - V{0.125}*cell[5] - x215;
auto x217 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x218 = V{4.5}*x82;
auto x219 = cell[10] + x38;
auto x220 = x218*x33*((cell[11] + 2*cell[13] - 2*cell[4] + x219 + x50 + x55)*(cell[11] + 2*cell[13] - 2*cell[4] + x219 + x50 + x55));
auto x221 = x145 + x88;
auto x222 = V{18}*x139;
auto x223 = x170 - x182*x42*x61 - x182*x52*x61 + x222*x42*x53 + x222*x43*x52;
auto x224 = -x208 - x210 + x212 + x213;
auto x225 = -V{0.125}*cell[13] + V{0.375}*cell[14] - V{0.125}*cell[4] + V{0.375}*cell[5] + x215;
auto x226 = -cell[11] + V{2}*cell[14] + cell[15] - V{2}*cell[5] + x160 + x219 + x26 + x56;
auto x227 = V{3.000006}*x139;
auto x228 = x227*x52*x61;
auto x229 = V{6.000012}*x139;
auto x230 = x229*x42*x53;
auto x231 = x227*x42*x61;
auto x232 = x229*x43*x52;
auto x233 = -x130;
auto x234 = -x131;
auto x235 = x188*x54;
auto x236 = x105*x191;
auto x237 = x179*x193;
auto x238 = x123*x237;
auto x239 = x114*x197;
auto x240 = x157*x199;
auto x241 = x119*x240;
auto x242 = V{0.0416666666666667}*cell[12] + V{6.93889390390723e-18}*cell[13] + V{6.93889390390723e-18}*cell[14] + V{0.0416666666666667}*cell[3] + V{6.93889390390723e-18}*cell[4] + V{6.93889390390723e-18}*cell[5] - x203*x62;
auto x243 = x204 + x233 + x234 + x235 + x236 + x238 - x239 - x241 + x242;
auto x244 = x179*x190;
auto x245 = x155*x244;
auto x246 = x123*x209;
auto x247 = x179*x196;
auto x248 = x114*x247;
auto x249 = x161*x240;
auto x250 = x245 + x246 - x248 - x249;
auto x251 = V{0.25}*x122*x33;
auto x252 = V{0.375}*cell[15] - V{0.125}*cell[16] + V{0.375}*cell[6] - V{0.125}*cell[7] - x251;
auto x253 = cell[10] + x34 + x41;
auto x254 = x218*x33*((cell[12] + 2*cell[15] - 2*cell[6] + cell[9] + x253 + x45 + x58)*(cell[12] + 2*cell[15] - 2*cell[6] + cell[9] + x253 + x45 + x58));
auto x255 = x139*x60;
auto x256 = V{9.000009}*x255;
auto x257 = x256*x43;
auto x258 = V{8.99999999948164e-06}*x255;
auto x259 = x258*x53;
auto x260 = x180 + x257 - x259;
auto x261 = V{5.999994}*x139;
auto x262 = x261*x52*x61;
auto x263 = V{12.000006}*x139;
auto x264 = x263*x42*x61;
auto x265 = -x147;
auto x266 = -x148;
auto x267 = x221 + x262 + x264 + x265 + x266;
auto x268 = -x245 - x246 + x248 + x249;
auto x269 = -V{0.125}*cell[15] + V{0.375}*cell[16] - V{0.125}*cell[6] + V{0.375}*cell[7] + x251;
auto x270 = -cell[12] + cell[13];
auto x271 = V{2}*cell[16] - V{2}*cell[7] + cell[8] + x159 + x253 + x27 + x270 + x47;
auto x272 = -x162;
auto x273 = -x163;
auto x274 = x188*x44;
auto x275 = x105*x207;
auto x276 = x161*x237;
auto x277 = x155*x211;
auto x278 = x117*x199;
auto x279 = x119*x278;
auto x280 = x202 + x205 + x242 + x272 + x273 + x274 + x275 + x276 - x277 - x279;
auto x281 = x114*x244;
auto x282 = x161*x194;
auto x283 = x155*x247;
auto x284 = x123*x278;
auto x285 = x281 + x282 - x283 - x284;
auto x286 = V{0.25}*x158*x33;
auto x287 = V{0.375}*cell[17] - V{0.125}*cell[18] + V{0.375}*cell[8] - V{0.125}*cell[9] - x286;
auto x288 = x46 + x51;
auto x289 = x218*x33*((cell[12] + 2*cell[17] + cell[7] - 2*cell[8] + x22 + x288 + x57 + x59)*(cell[12] + 2*cell[17] + cell[7] - 2*cell[8] + x22 + x288 + x57 + x59));
auto x290 = x256*x53;
auto x291 = x258*x43;
auto x292 = x180 + x290 - x291;
auto x293 = x261*x42*x61;
auto x294 = x263*x52*x61;
auto x295 = -x172;
auto x296 = -x173;
auto x297 = x170 + x88;
auto x298 = x293 + x294 + x295 + x296 + x297;
auto x299 = -x281 - x282 + x283 + x284;
auto x300 = -V{0.125}*cell[17] + V{0.375}*cell[18] - V{0.125}*cell[8] + V{0.375}*cell[9] + x286;
auto x301 = V{2}*cell[18] + cell[3] - V{2}*cell[9] + x121 + x270 + x28 + x288 + x35;
auto x302 = -x63;
auto x303 = x302 - x64 - x73 + x95 + x96 + x97 + x98;
auto x304 = x109 + x110 + x111 + x112 - x66 - x75;
auto x305 = x106 + x107 - x108 + x303 + x304 - x69 - x70 - x78 - x79;
auto x306 = x100 + x101 + x102 + x103 - x65 - x74;
auto x307 = x303 + x306 - x67 - x68 - x76 - x77 + x91 + x92 - x94;
auto x308 = -x119;
auto x309 = -x123;
auto x310 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x311 = x84 + V{-1};
auto x312 = x152 + x153 - x154 + x302 + x304 + x306 - x71 - x72 - x80 - x81;
auto x313 = -x161;
auto x314 = x186 + x187 + x189 - x192 - x195 + x198 + x201 + x204 + x205;
auto x315 = x137 + x185;
auto x316 = x145 + x315;
auto x317 = x204 + x233 + x234 + x235 - x236 - x238 + x239 + x241 + x242;
auto x318 = x202 + x205 + x242 + x272 + x273 + x274 - x275 - x276 + x277 + x279;
auto x0 = x21*(-V{0.5}*cell[0] + V{4.16333634234434e-17}*cell[10] + V{4.16333634234434e-17}*cell[11] + V{4.16333634234434e-17}*cell[12] + V{0.5}*cell[13] + V{0.5}*cell[14] + V{0.5}*cell[15] + V{0.5}*cell[16] + V{0.5}*cell[17] + V{0.5}*cell[18] + V{4.16333634234434e-17}*cell[1] + V{4.16333634234434e-17}*cell[2] + V{4.16333634234434e-17}*cell[3] + V{0.5}*cell[4] + V{0.5}*cell[5] + V{0.5}*cell[6] + V{0.5}*cell[7] + V{0.5}*cell[8] + V{0.5}*cell[9] - x32*x44 - x32*x54 - x32*x62) - x88*(x63 + x64 + x65 + x66 + x67 + x68 + x69 + x70 + x71 + x72 + x73 + x74 + x75 + x76 + x77 + x78 + x79 + x80 + x81 + V{0.333333333333333}) + V{-0.333333333333333};
auto x1 = V{0.0555555555555556}*x134*(x136 + x137 + x138 + x150) - x21*(-x105*x90 - x114*x90 - x116*x119 - x120*x123 - x133) + V{-0.0555555555555556};
auto x2 = V{0.0555555555555556}*x134*(x137 + x166 + x167 + x175 + V{1}) - x21*(-x105*x151 - x119*x156 - x120*x161 - x151*x155 - x165) + V{-0.0555555555555556};
auto x3 = V{0.0555555555555556}*x134*(x178 + x184 + x185) - x21*(-x114*x176 - x116*x161 - x123*x156 - x155*x176 - x177) + V{-0.0555555555555556};
auto x4 = -x21*(x206 + x214 + x216) - x217*(-x220 + x221 + x223) + V{-0.0277777777777778};
auto x5 = -(x21*(x206 + x224 + x225) + x217*(x171 - x218*x33*((x226)*(x226)) + x221 + x228 + x230 - x231 - x232) + V{0.0277777777777778});
auto x6 = -x21*(x243 + x250 + x252) - x217*(-x254 + x260 + x267) + V{-0.0277777777777778};
auto x7 = -(x21*(x243 + x268 + x269) + x217*(x181 - x218*x33*((x271)*(x271)) - x257 + x259 + x267) + V{0.0277777777777778});
auto x8 = -x21*(x280 + x285 + x287) - x217*(-x289 + x292 + x298) + V{-0.0277777777777778};
auto x9 = -(x21*(x280 + x299 + x300) + x217*(x181 - x218*x33*((x301)*(x301)) - x290 + x291 + x298) + V{0.0277777777777778});
auto x10 = x21*(x116*x308 + x120*x309 + x133 + x305*x90 + x307*x90) - x310*(-x136 + x150 + x311 + x86) + V{-0.0555555555555556};
auto x11 = x21*(x120*x313 + x151*x307 + x151*x312 + x156*x308 + x165) - x310*(-x167 + x175 + x311 + x85) + V{-0.0555555555555556};
auto x12 = x21*(x116*x313 + x156*x309 + x176*x305 + x176*x312 + x177) - x310*(-x178 + x184 + x87) + V{-0.0555555555555556};
auto x13 = V{0.0277777777777778}*x134*(x220 + x223 + x316) - x21*(x216 + x224 + x314) + V{-0.0277777777777778};
auto x14 = -(x21*(x214 + x225 + x314) + x217*(x146 - x218*x33*((x226)*(x226)) - x228 - x230 + x231 + x232 + x297) + V{0.0277777777777778});
auto x15 = V{0.0277777777777778}*x134*(x254 + x260 + x262 + x264 + x265 + x266 + x316) - x21*(x252 + x268 + x317) + V{-0.0277777777777778};
auto x16 = -(x21*(x250 + x269 + x317) + x217*(x149 - x218*x33*((x271)*(x271)) + x260 - x262 - x264 + x88) + V{0.0277777777777778});
auto x17 = V{0.0277777777777778}*x134*(x170 + x289 + x292 + x293 + x294 + x295 + x296 + x315) - x21*(x287 + x299 + x318) + V{-0.0277777777777778};
auto x18 = -(x21*(x285 + x300 + x318) + x217*(x174 - x218*x33*((x301)*(x301)) + x292 - x293 - x294 + x88) + V{0.0277777777777778});
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
return { x30, V{1}*x33*x82*(x43 + x53 + x61) };
}
};

}

}

// Generation Info: commit=cc571609e6a7d0f0c4e0722342de06c84dd07797
