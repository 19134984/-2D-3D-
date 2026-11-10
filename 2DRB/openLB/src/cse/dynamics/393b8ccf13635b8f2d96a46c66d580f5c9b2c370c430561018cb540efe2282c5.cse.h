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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::ParameterFromCell<collision::HYBRID, collision::LocalSmagorinskyEffectiveOmega<collision::HRR> > >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x196 = cell.template getFieldComponent<olb::collision::HYBRID>(0);
auto x198 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x199 = cell.template getFieldComponent<olb::descriptors::TENSOR>(0);
auto x200 = cell.template getFieldComponent<olb::descriptors::TENSOR>(1);
auto x201 = cell.template getFieldComponent<olb::descriptors::TENSOR>(2);
auto x202 = cell.template getFieldComponent<olb::descriptors::TENSOR>(3);
auto x203 = cell.template getFieldComponent<olb::descriptors::TENSOR>(4);
auto x204 = cell.template getFieldComponent<olb::descriptors::TENSOR>(5);
auto x205 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x206 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x207 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x208 = cell.template getFieldComponent<olb::descriptors::WMPOROSITY>(0);
auto x209 = cell.template getFieldComponent<olb::descriptors::WMVELOCITY>(0);
auto x210 = cell.template getFieldComponent<olb::descriptors::WMVELOCITY>(1);
auto x211 = cell.template getFieldComponent<olb::descriptors::WMVELOCITY>(2);
auto x213 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x214 = parameters.template get<descriptors::OMEGA>();
auto x197 = cell[10] + cell[14] + cell[16];
auto x212 = cell[1] + cell[4] + x197;
auto x215 = cell[15] + cell[7];
auto x216 = cell[17] + cell[9];
auto x217 = cell[12] + x215 + x216;
auto x218 = cell[3] + cell[8] + x217;
auto x219 = cell[13] + cell[5];
auto x220 = cell[11] + cell[18] + x219;
auto x221 = cell[2] + x220;
auto x0 = V{3};
auto x230 = ((x0)*(x0));
auto x231 = V{1} / (x0);
auto x232 = V{1}*x231;
auto x1 = V{1} - x208;
auto x2 = V{1} - x198;
auto x3 = -V{1}*cell[4];
auto x4 = -V{1}*cell[5];
auto x223 = cell[13] + x4;
auto x5 = -V{1}*cell[6];
auto x6 = -V{1}*cell[7];
auto x222 = cell[15] + x6;
auto x7 = cell[0] + cell[6] + x212 + x218 + x221 + V{1};
auto x224 = V{1} / (x7);
auto x225 = x198*x224;
auto x8 = x1*x209 + x208*(x2*x205 + x225*(-V{1}*cell[1] + x197 + x222 + x223 + x3 + x5));
auto x226 = ((x8)*(x8));
auto x9 = x226;
auto x10 = -V{1}*cell[8];
auto x11 = -V{1}*cell[9];
auto x227 = cell[17] + x11;
auto x12 = -V{1}*cell[14];
auto x13 = x1*x210 + x208*(x2*x206 + x225*(-V{1}*cell[2] + x10 + x12 + x220 + x227 + x3));
auto x228 = ((x13)*(x13));
auto x14 = x228;
auto x15 = -V{1}*cell[16];
auto x16 = -V{1}*cell[18];
auto x17 = x1*x211 + x208*(x2*x207 + x225*(-V{1}*cell[3] + x10 + x15 + x16 + x217 + x5));
auto x229 = ((x17)*(x17));
auto x18 = x14 + x229 + x9;
auto x19 = V{0.5}*x0*x18;
auto x254 = V{1} - x19;
auto x20 = V{0.333333333333333};
auto x21 = x230;
auto x22 = V{0.5}*x21;
auto x23 = x232;
auto x238 = V{1} - x23;
auto x24 = -x23;
auto x25 = V{1}/x214;
auto x26 = x7*x8;
auto x27 = cell[6] + x15 - x17*x26 + x222;
auto x28 = ((x27)*(x27));
auto x29 = cell[4] + x12 - x13*x26 + x223;
auto x30 = x23*(x7 + V{-1});
auto x233 = -x30;
auto x234 = cell[6] + x233;
auto x31 = x212 + x215 + x219 + x234 - x26*x8;
auto x32 = cell[16] + cell[18] + x218 - x229*x7 + x234;
auto x33 = x13*x7;
auto x34 = cell[8] + x16 - x17*x33 + x227;
auto x35 = cell[14] + cell[4] + cell[8] - x13*x33 + x216 + x221 + x233;
auto x36 = V{2}*V{1} / (x25 + util::sqrt(V{2.828428}*x224*x230*((x213)*(x213))*util::sqrt(V{2}*x28 + V{2}*((x29)*(x29)) + ((x31)*(x31)) + ((x32)*(x32)) + V{2}*((x34)*(x34)) + ((x35)*(x35))) + ((x25)*(x25))));
auto x235 = V{1} / (x36);
auto x38 = V{2}*x231*x235*x7*(x196 + V{-1});
auto x39 = x196*x32 + x204*x38;
auto x40 = x24*x39;
auto x41 = x196*x35 + x202*x38;
auto x42 = x24*x41;
auto x43 = x196*x31 + x199*x38;
auto x44 = x24*x43;
auto x45 = x42 + x44;
auto x46 = V{1} - x36;
auto x47 = x0*x21;
auto x48 = V{0.166667}*x47;
auto x250 = x24*x48;
auto x49 = -x232;
auto x50 = -x49;
auto x236 = x48*x50;
auto x51 = -x236;
auto x52 = x13*x9;
auto x53 = x228*x8;
auto x54 = x52 - x53;
auto x239 = x48*x54;
auto x55 = V{0.5}*x47;
auto x237 = x50*x55;
auto x251 = x24*x55;
auto x56 = x237;
auto x57 = x52 + x53;
auto x240 = x55*x57;
auto x58 = x51*x54 + x56*x57;
auto x59 = x229*x8;
auto x60 = x13*x229;
auto x61 = x59 - x60;
auto x246 = x48*x61;
auto x261 = x246*(x24 - x50);
auto x62 = x59 + x60;
auto x243 = x55*x62;
auto x63 = x56*x62;
auto x273 = x254 + x63;
auto x64 = -V{1}*x8;
auto x65 = V{0.0555555555555556};
auto x252 = -x65;
auto x253 = x46*x65;
auto x66 = x65*x7;
auto x67 = -V{1}*x24;
auto x68 = 0;
auto x69 = x196*x27 + x201*x38;
auto x248 = x21*x69;
auto x70 = V{2}*x17;
auto x71 = x39*x8 + x69*x70;
auto x72 = x196*x29 + x200*x38;
auto x242 = x21*x72;
auto x73 = V{2}*x13;
auto x74 = x41*x8 + x72*x73;
auto x75 = x71 - x74;
auto x76 = V{2}*x55*x67;
auto x77 = x71 + x74;
auto x78 = x238;
auto x244 = -x78;
auto x245 = x48*(x24 + x244);
auto x265 = x246*x78;
auto x79 = x43*x78;
auto x80 = x42 + x79;
auto x81 = x22*(x40 + x80);
auto x82 = x67*(x239 + x240);
auto x83 = -V{1}*x13;
auto x84 = V{2}*x8;
auto x85 = x13*x43 + x72*x84;
auto x86 = x196*x34 + x203*x38;
auto x249 = x21*x86;
auto x87 = x13*x39 + x70*x86;
auto x88 = x85 - x87;
auto x89 = x85 + x87;
auto x90 = x41*x78;
auto x91 = x44 + x90;
auto x92 = x22*(x40 + x91);
auto x93 = x14*x17;
auto x94 = x17*x9;
auto x95 = x93 - x94;
auto x96 = x93 + x94;
auto x97 = -V{1}*x17;
auto x98 = x17*x41 + x73*x86;
auto x99 = x17*x43 + x69*x84;
auto x100 = x98 - x99;
auto x101 = x98 + x99;
auto x102 = x39*x78;
auto x103 = x22*(x102 + x45);
auto x104 = -V{1}*x78;
auto x241 = x48*(x104 - x67);
auto x105 = -x49 + V{-1};
auto x247 = x105*x246;
auto x106 = x64 + x83;
auto x107 = V{0.0277777777777778};
auto x256 = -x107;
auto x259 = x107*x46;
auto x108 = x107*x7;
auto x109 = -x241;
auto x110 = x109*x75;
auto x111 = x241;
auto x112 = x111*x88;
auto x113 = x55*(x104 + x67);
auto x114 = x113*x77;
auto x274 = x110 + x114;
auto x115 = x113*x89;
auto x260 = x112 + x115;
auto x116 = x22*(x40 + x79 + x90);
auto x117 = x116 + x242;
auto x118 = x243*(x24 + x50);
auto x262 = x118 + x254;
auto x119 = x13 + x64;
auto x120 = -x245;
auto x121 = x120*x88;
auto x122 = x55*(x24 + x78);
auto x123 = x122*x89;
auto x258 = x121 + x123;
auto x124 = x116 - V{1}*x242;
auto x125 = x109*x95;
auto x126 = x247;
auto x277 = x126 + x58;
auto x127 = x113*x96;
auto x128 = x105*x243;
auto x129 = x127 + x128;
auto x276 = x129 + x254;
auto x130 = x64 + x97;
auto x131 = x100*x109;
auto x132 = x111*x75;
auto x275 = x114 + x132;
auto x133 = x101*x113;
auto x267 = x131 + x133;
auto x134 = x22*(x102 + x80);
auto x135 = x134 + x248;
auto x136 = x245;
auto x137 = x136*x95;
auto x138 = x122*x96;
auto x139 = x128 + x138;
auto x278 = x139 + x254;
auto x140 = x17 + x64;
auto x141 = x100*x136;
auto x142 = x101*x122;
auto x264 = x141 + x142;
auto x143 = x134 - V{1}*x248;
auto x144 = x111*x95;
auto x145 = -x247;
auto x280 = x145 + x82;
auto x146 = x83 + x97;
auto x147 = x100*x111;
auto x272 = x133 + x147;
auto x148 = x109*x88;
auto x279 = x115 + x148;
auto x149 = x22*(x102 + x91);
auto x150 = x149 + x249;
auto x151 = x120*x95;
auto x152 = x17 + x83;
auto x153 = x100*x120;
auto x270 = x142 + x153;
auto x154 = x149 - V{1}*x249;
auto x155 = -x250;
auto x156 = x251;
auto x157 = x156*x57;
auto x158 = x155*x54 + x157;
auto x266 = x158 + x254;
auto x159 = x250;
auto x160 = x156*x62;
auto x255 = x160 + x254;
auto x161 = 0;
auto x162 = V{2}*x251;
auto x163 = x157 + x159*x54;
auto x271 = x163 + x254;
auto x164 = x13 + x8;
auto x165 = x136*x75;
auto x166 = x122*x77;
auto x257 = x165 + x166;
auto x167 = x238;
auto x268 = x167*x246;
auto x168 = x8 + x83;
auto x169 = x138 + x243*x78;
auto x170 = x17 + x8;
auto x171 = x120*x75;
auto x263 = x166 + x171;
auto x172 = x127 + x167*x243;
auto x173 = x8 + x97;
auto x174 = x13 + x17;
auto x175 = x136*x88;
auto x269 = x123 + x175;
auto x176 = x13 + x97;
cell[0] = x20*(x22*x46*(x40 + x45) - x7*(x19 + V{-1}) + V{-1});
cell[10] = x252 + x253*(x161*x75 + x162*x77 + x81) + x66*(x0*x8 + x158 + x159*x61 + x22*x226 + x255);
cell[11] = x252 + x253*(x161*x88 + x162*x89 + x92) + x66*(x0*x13 + x155*x61 + x163 + x22*x228 + x255);
cell[12] = x252 + x253*(x100*x161 + x101*x162 + x103) + x66*(x0*x17 + x161*x95 + x162*x96 + x22*x229 + x254);
cell[13] = x108*(x0*x164 + x161*x61 + x162*x62 + x22*((x164)*(x164)) + V{2}*x240*x78 + x254) + x256 + x259*(x117 + x257 + x258);
cell[14] = x108*(x0*x168 + x22*((x168)*(x168)) + x239*(x104 - x167) + x240*(x104 + x167) + x261 + x262) + x256 + x259*(x124 + x257 + x260);
cell[15] = x108*(x0*x170 + x137 + x169 + x22*((x170)*(x170)) + x265 + x266) + x256 + x259*(x135 + x263 + x264);
cell[16] = x108*(x0*x173 + x125 + x172 + x22*((x173)*(x173)) + x266 + x268) + x256 + x259*(x143 + x263 + x267);
cell[17] = x108*(x0*x174 + x151 + x169 + x22*((x174)*(x174)) - x265 + x271) + x256 + x259*(x150 + x269 + x270);
cell[18] = x108*(x0*x176 + x144 + x172 + x22*((x176)*(x176)) - x268 + x271) + x256 + x259*(x154 + x269 + x272);
cell[1] = x252 + x253*(x68*x75 + x76*x77 + x81) + x66*(x0*x64 + x22*((x64)*(x64)) + x236*x61 + x273 + x58);
cell[2] = x252 + x253*(x68*x88 + x76*x89 + x92) + x66*(x0*x83 + x22*((x83)*(x83)) + x273 + x51*x61 + x82);
cell[3] = x252 + x253*(x100*x68 + x101*x76 + x103) + x66*(x0*x97 + x22*((x97)*(x97)) + x254 + x68*x95 + x76*x96);
cell[4] = x108*(x0*x106 + x22*((x106)*(x106)) + V{2}*x237*x62 + x239*(x104 - x105) + x240*(x104 + x105) + x254) + x256 + x259*(x117 + x260 + x274);
cell[5] = x108*(x0*x119 + x22*((x119)*(x119)) - x239*(x105 + x244) + x240*(x105 + x78) - x261 + x262) + x256 + x259*(x124 + x258 + x274);
cell[6] = x108*(x0*x130 + x125 + x22*((x130)*(x130)) + x276 + x277) + x256 + x259*(x135 + x267 + x275);
cell[7] = x108*(x0*x140 + x137 + x22*((x140)*(x140)) + x277 + x278) + x256 + x259*(x143 + x264 + x275);
cell[8] = x108*(x0*x146 + x144 + x22*((x146)*(x146)) + x276 + x280) + x256 + x259*(x150 + x272 + x279);
cell[9] = x108*(x0*x152 + x151 + x22*((x152)*(x152)) + x278 + x280) + x256 + x259*(x154 + x270 + x279);
cell.template getFieldPointer<olb::descriptors::OMEGA>()[0] = x36;
return { x7, x18 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
