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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::MovingPorousMomentumCombinationNoWM<momenta::BulkMomentum>, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x212 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x213 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x214 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x215 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x216 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x217 = parameters.template get<descriptors::OMEGA>();
auto x218 = cell[17] + cell[9];
auto x219 = cell[15] + cell[7];
auto x220 = cell[12] + x218 + x219;
auto x221 = cell[3] + cell[6] + x220;
auto x222 = cell[13] + cell[5];
auto x223 = cell[11] + cell[18] + x222;
auto x224 = cell[2] + cell[4] + x223;
auto x225 = cell[10] + cell[14] + cell[16];
auto x226 = cell[1] + x225;
auto x0 = V{3};
auto x235 = ((x0)*(x0));
auto x236 = V{1}/x0;
auto x1 = V{1} - x212;
auto x2 = x1*x213;
auto x3 = -V{1}*cell[4];
auto x4 = -V{1}*cell[5];
auto x228 = cell[13] + x4;
auto x5 = -V{1}*cell[6];
auto x6 = -V{1}*cell[7];
auto x227 = cell[15] + x6;
auto x7 = cell[0] + cell[8] + x221 + x224 + x226 + V{1};
auto x229 = V{1} / (x7);
auto x230 = x212*x229;
auto x8 = x2 + x230*(-V{1}*cell[1] + x225 + x227 + x228 + x3 + x5);
auto x231 = ((x8)*(x8));
auto x9 = x231;
auto x10 = x1*x214;
auto x11 = -V{1}*cell[8];
auto x12 = -V{1}*cell[9];
auto x232 = cell[17] + x12;
auto x13 = -V{1}*cell[14];
auto x14 = x10 + x230*(-V{1}*cell[2] + x11 + x13 + x223 + x232 + x3);
auto x233 = ((x14)*(x14));
auto x15 = x233;
auto x16 = x1*x215;
auto x17 = -V{1}*cell[16];
auto x18 = -V{1}*cell[18];
auto x19 = x16 + x230*(-V{1}*cell[3] + x11 + x17 + x18 + x220 + x5);
auto x234 = ((x19)*(x19));
auto x20 = x15 + x234 + x9;
auto x21 = V{0.5}*x0*x20;
auto x244 = V{1} - x21;
auto x22 = V{0.333333333333333};
auto x23 = x235;
auto x24 = V{0.5}*x23;
auto x25 = x236;
auto x241 = V{1} - x25;
auto x26 = -x25;
auto x27 = x25*(x7 + V{-1});
auto x237 = -x27;
auto x238 = cell[8] + x237;
auto x28 = cell[16] + cell[18] + x221 - x234*x7 + x238;
auto x29 = x26*x28;
auto x30 = x14*x7;
auto x31 = cell[14] - x14*x30 + x218 + x224 + x238;
auto x32 = x26*x31;
auto x33 = x7*x8;
auto x34 = cell[4] + cell[6] + x219 + x222 + x226 + x237 - x33*x8;
auto x35 = x26*x34;
auto x36 = x32 + x35;
auto x37 = V{1}/x217;
auto x38 = cell[6] + x17 - x19*x33 + x227;
auto x262 = x23*x38;
auto x39 = ((x38)*(x38));
auto x40 = cell[4] + x13 - x14*x33 + x228;
auto x249 = x23*x40;
auto x41 = cell[8] + x18 - x19*x30 + x232;
auto x269 = x23*x41;
auto x42 = V{1} - V{2}*V{1} / (x37 + util::sqrt(V{2.828428}*x229*x235*((x216)*(x216))*util::sqrt(V{2}*x39 + ((x28)*(x28)) + ((x31)*(x31)) + ((x34)*(x34)) + V{2}*((x40)*(x40)) + V{2}*((x41)*(x41))) + ((x37)*(x37))));
auto x43 = x22*(x24*x42*(x29 + x36) - x7*(x21 + V{-1}) + V{-1});
auto x44 = x0*x23;
auto x45 = V{0.166667}*x44;
auto x274 = x26*x45;
auto x46 = -x236;
auto x47 = -x46;
auto x239 = x45*x47;
auto x48 = -x239;
auto x49 = x14*x9;
auto x50 = x233*x8;
auto x51 = x49 - x50;
auto x246 = x45*x51;
auto x52 = V{0.5}*x44;
auto x240 = x47*x52;
auto x275 = x26*x52;
auto x53 = x240;
auto x54 = x49 + x50;
auto x247 = x52*x54;
auto x55 = x48*x51 + x53*x54;
auto x56 = x234*x8;
auto x57 = x14*x234;
auto x58 = x56 - x57;
auto x258 = x45*x58;
auto x259 = x258*(x26 - x47);
auto x59 = x56 + x57;
auto x254 = x52*x59;
auto x60 = x53*x59;
auto x245 = x244 + x60;
auto x61 = -V{1}*x8;
auto x62 = V{0.0555555555555556};
auto x242 = -x62;
auto x243 = x42*x62;
auto x63 = x62*x7;
auto x64 = -V{1}*x26;
auto x65 = 0;
auto x66 = V{2}*x19;
auto x67 = x28*x8 + x38*x66;
auto x68 = V{2}*x14;
auto x69 = x31*x8 + x40*x68;
auto x70 = x67 - x69;
auto x71 = V{2}*x52*x64;
auto x72 = x67 + x69;
auto x73 = x241;
auto x255 = -x73;
auto x256 = x45*(x255 + x26);
auto x279 = x258*x73;
auto x74 = x34*x73;
auto x75 = x32 + x74;
auto x76 = x24*(x29 + x75);
auto x77 = x242 + x243*(x65*x70 + x71*x72 + x76) + x63*(x0*x61 + x239*x58 + x24*((x61)*(x61)) + x245 + x55);
auto x78 = x64*(x246 + x247);
auto x79 = -V{1}*x14;
auto x80 = V{2}*x8;
auto x81 = x14*x34 + x40*x80;
auto x82 = x14*x28 + x41*x66;
auto x83 = x81 - x82;
auto x84 = x81 + x82;
auto x85 = x31*x73;
auto x86 = x35 + x85;
auto x87 = x24*(x29 + x86);
auto x88 = x242 + x243*(x65*x83 + x71*x84 + x87) + x63*(x0*x79 + x24*((x79)*(x79)) + x245 + x48*x58 + x78);
auto x89 = x15*x19;
auto x90 = x19*x9;
auto x91 = x89 - x90;
auto x92 = x89 + x90;
auto x93 = -V{1}*x19;
auto x94 = x19*x31 + x41*x68;
auto x95 = x19*x34 + x38*x80;
auto x96 = x94 - x95;
auto x97 = x94 + x95;
auto x98 = x28*x73;
auto x99 = x24*(x36 + x98);
auto x100 = x242 + x243*(x65*x96 + x71*x97 + x99) + x63*(x0*x93 + x24*((x93)*(x93)) + x244 + x65*x91 + x71*x92);
auto x101 = -V{1}*x73;
auto x248 = x45*(x101 - x64);
auto x102 = -x46 + V{-1};
auto x261 = x102*x258;
auto x103 = x61 + x79;
auto x104 = V{0.0277777777777778};
auto x250 = -x104;
auto x253 = x104*x42;
auto x105 = x104*x7;
auto x106 = -x248;
auto x107 = x106*x70;
auto x108 = x248;
auto x109 = x108*x83;
auto x110 = x52*(x101 + x64);
auto x111 = x110*x72;
auto x251 = x107 + x111;
auto x112 = x110*x84;
auto x252 = x109 + x112;
auto x113 = x24*(x29 + x74 + x85);
auto x114 = x113 + x249;
auto x115 = x105*(x0*x103 + x24*((x103)*(x103)) + V{2}*x240*x59 + x244 + x246*(x101 - x102) + x247*(x101 + x102)) + x250 + x253*(x114 + x251 + x252);
auto x116 = x254*(x26 + x47);
auto x260 = x116 + x244;
auto x117 = x14 + x61;
auto x118 = -x256;
auto x119 = x118*x83;
auto x120 = x52*(x26 + x73);
auto x121 = x120*x84;
auto x257 = x119 + x121;
auto x122 = x113 - V{1}*x249;
auto x123 = x105*(x0*x117 + x24*((x117)*(x117)) - x246*(x102 + x255) + x247*(x102 + x73) - x259 + x260) + x250 + x253*(x122 + x251 + x257);
auto x124 = x106*x91;
auto x125 = x261;
auto x266 = x125 + x55;
auto x126 = x110*x92;
auto x127 = x102*x254;
auto x128 = x126 + x127;
auto x265 = x128 + x244;
auto x129 = x61 + x93;
auto x130 = x106*x96;
auto x131 = x108*x70;
auto x263 = x111 + x131;
auto x132 = x110*x97;
auto x264 = x130 + x132;
auto x133 = x24*(x75 + x98);
auto x134 = x133 + x262;
auto x135 = x105*(x0*x129 + x124 + x24*((x129)*(x129)) + x265 + x266) + x250 + x253*(x134 + x263 + x264);
auto x136 = x256;
auto x137 = x136*x91;
auto x138 = x120*x92;
auto x139 = x127 + x138;
auto x268 = x139 + x244;
auto x140 = x19 + x61;
auto x141 = x136*x96;
auto x142 = x120*x97;
auto x267 = x141 + x142;
auto x143 = x133 - V{1}*x262;
auto x144 = x105*(x0*x140 + x137 + x24*((x140)*(x140)) + x266 + x268) + x250 + x253*(x143 + x263 + x267);
auto x145 = x108*x91;
auto x146 = -x261;
auto x272 = x146 + x78;
auto x147 = x79 + x93;
auto x148 = x108*x96;
auto x271 = x132 + x148;
auto x149 = x106*x83;
auto x270 = x112 + x149;
auto x150 = x24*(x86 + x98);
auto x151 = x150 + x269;
auto x152 = x105*(x0*x147 + x145 + x24*((x147)*(x147)) + x265 + x272) + x250 + x253*(x151 + x270 + x271);
auto x153 = x118*x91;
auto x154 = x19 + x79;
auto x155 = x118*x96;
auto x273 = x142 + x155;
auto x156 = x150 - V{1}*x269;
auto x157 = x105*(x0*x154 + x153 + x24*((x154)*(x154)) + x268 + x272) + x250 + x253*(x156 + x270 + x273);
auto x158 = -x274;
auto x159 = x275;
auto x160 = x159*x54;
auto x161 = x158*x51 + x160;
auto x280 = x161 + x244;
auto x162 = x274;
auto x163 = x159*x59;
auto x276 = x163 + x244;
auto x164 = 0;
auto x165 = V{2}*x275;
auto x166 = x242 + x243*(x164*x70 + x165*x72 + x76) + x63*(x0*x8 + x161 + x162*x58 + x231*x24 + x276);
auto x167 = x160 + x162*x51;
auto x283 = x167 + x244;
auto x168 = x242 + x243*(x164*x83 + x165*x84 + x87) + x63*(x0*x14 + x158*x58 + x167 + x233*x24 + x276);
auto x169 = x242 + x243*(x164*x96 + x165*x97 + x99) + x63*(x0*x19 + x164*x91 + x165*x92 + x234*x24 + x244);
auto x286 = x144 + x157 + x169;
auto x170 = x14 + x8;
auto x171 = x136*x70;
auto x172 = x120*x72;
auto x277 = x171 + x172;
auto x173 = x105*(x0*x170 + x164*x58 + x165*x59 + x24*((x170)*(x170)) + x244 + V{2}*x247*x73) + x250 + x253*(x114 + x257 + x277);
auto x174 = x241;
auto x281 = x174*x258;
auto x175 = x79 + x8;
auto x176 = x105*(x0*x175 + x24*((x175)*(x175)) + x246*(x101 - x174) + x247*(x101 + x174) + x259 + x260) + x250 + x253*(x122 + x252 + x277);
auto x177 = x138 + x254*x73;
auto x178 = x19 + x8;
auto x179 = x118*x70;
auto x278 = x172 + x179;
auto x180 = x105*(x0*x178 + x137 + x177 + x24*((x178)*(x178)) + x279 + x280) + x250 + x253*(x134 + x267 + x278);
auto x181 = x126 + x174*x254;
auto x182 = x8 + x93;
auto x183 = x105*(x0*x182 + x124 + x181 + x24*((x182)*(x182)) + x280 + x281) + x250 + x253*(x143 + x264 + x278);
auto x284 = x166 + x176 + x183;
auto x184 = x14 + x19;
auto x185 = x136*x83;
auto x282 = x121 + x185;
auto x186 = x105*(x0*x184 + x153 + x177 + x24*((x184)*(x184)) - x279 + x283) + x250 + x253*(x151 + x273 + x282);
auto x187 = x14 + x93;
auto x188 = x105*(x0*x187 + x145 + x181 + x24*((x187)*(x187)) - x281 + x283) + x250 + x253*(x156 + x271 + x282);
auto x285 = x123 + x168 + x188;
auto x189 = x100 + x115 + x135 + x152 + x173 + x180 + x186 + x284 + x285 + x286 + x43 + x77 + x88 + V{1};
auto x289 = x212/x189;
auto x190 = -V{1}*x115;
auto x287 = x173 + x190;
auto x191 = -V{1}*x135;
auto x288 = x180 + x191;
auto x192 = -V{1}*x152;
auto x290 = x186 + x192;
cell[0] = x43;
cell[10] = x166;
cell[11] = x168;
cell[12] = x169;
cell[13] = x173;
cell[14] = x176;
cell[15] = x180;
cell[16] = x183;
cell[17] = x186;
cell[18] = x188;
cell[1] = x77;
cell[2] = x88;
cell[3] = x100;
cell[4] = x115;
cell[5] = x123;
cell[6] = x135;
cell[7] = x144;
cell[8] = x152;
cell[9] = x157;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x189;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x2 + x289*(-V{1}*x123 - V{1}*x144 + x284 + x287 + x288 - V{1}*x77);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x10 + x289*(-V{1}*x157 - V{1}*x176 + x285 + x287 + x290 - V{1}*x88);
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x16 + x289*(-V{1}*x100 - V{1}*x183 - V{1}*x188 + x286 + x288 + x290);
return { x7, x20 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
