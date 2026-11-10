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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::ParameterFromCell<collision::LES::SMAGORINSKY, collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB> >, forcing::Guo<momenta::ForcedWithStress> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x221 = cell.template getFieldComponent<olb::collision::LES::SMAGORINSKY>(0);
auto x222 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x223 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x224 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x226 = parameters.template get<descriptors::OMEGA>();
auto x225 = cell[17] + cell[9];
auto x227 = cell[15] + cell[7];
auto x228 = cell[12] + x225 + x227;
auto x229 = cell[3] + cell[6] + x228;
auto x230 = cell[13] + cell[5];
auto x231 = cell[11] + cell[18] + x230;
auto x232 = cell[2] + cell[4] + x231;
auto x233 = cell[10] + cell[14] + cell[16];
auto x234 = cell[1] + x233;
auto x0 = V{3};
auto x242 = ((x0)*(x0));
auto x243 = V{1}/x0;
auto x1 = -V{1}*cell[4];
auto x2 = -V{1}*cell[5];
auto x237 = cell[13] + x2;
auto x3 = -V{1}*cell[6];
auto x4 = -V{1}*cell[7];
auto x236 = cell[15] + x4;
auto x5 = cell[0] + cell[8] + x229 + x232 + x234 + V{1};
auto x235 = V{1} / (x5);
auto x6 = V{0.5}*x222;
auto x7 = x235*(-V{1}*cell[1] + x1 + x233 + x236 + x237 + x3) + x6;
auto x238 = ((x7)*(x7));
auto x247 = x0*x7;
auto x8 = x238;
auto x9 = -V{1}*cell[8];
auto x10 = -V{1}*cell[9];
auto x239 = cell[17] + x10;
auto x11 = -V{1}*cell[14];
auto x12 = V{0.5}*x223;
auto x13 = x12 + x235*(-V{1}*cell[2] + x1 + x11 + x231 + x239 + x9);
auto x240 = ((x13)*(x13));
auto x248 = x0*x13;
auto x14 = x240;
auto x15 = -V{1}*cell[16];
auto x16 = -V{1}*cell[18];
auto x17 = V{0.5}*x224;
auto x18 = x17 + x235*(-V{1}*cell[3] + x15 + x16 + x228 + x3 + x9);
auto x241 = ((x18)*(x18));
auto x249 = x0*x18;
auto x19 = x14 + x241 + x8;
auto x20 = V{0.5}*x0*x19;
auto x270 = V{1} - x20;
auto x21 = V{0.333333333333333};
auto x22 = x242;
auto x257 = V{1}*x22;
auto x23 = V{0.5}*x22;
auto x24 = x243;
auto x252 = V{1} - x24;
auto x25 = -x24;
auto x26 = -x17 + x18;
auto x27 = x24*(x5 + V{-1});
auto x245 = -x27;
auto x246 = cell[8] + x245;
auto x28 = V{0.5}*x5;
auto x244 = V{2}*x28;
auto x29 = x224*x26;
auto x30 = cell[16] + cell[18] + x229 + x244*x29 + x246 - x5*((x26)*(x26));
auto x31 = x25*x30;
auto x32 = -x12 + x13;
auto x33 = x32*x5;
auto x34 = x223*x32;
auto x35 = cell[14] + x225 + x232 + x244*x34 + x246 - x32*x33;
auto x36 = x25*x35;
auto x37 = -x6 + x7;
auto x38 = x37*x5;
auto x39 = x222*x37;
auto x40 = cell[4] + cell[6] + x227 + x230 + x234 + x244*x39 + x245 - x37*x38;
auto x41 = x25*x40;
auto x42 = x36 + x41;
auto x43 = V{1}/x226;
auto x44 = cell[6] + x15 + x236 - x26*x38 + x28*(x222*x26 + x224*x37);
auto x263 = x22*x44;
auto x45 = ((x44)*(x44));
auto x46 = cell[4] + x11 + x237 + x28*(x222*x32 + x223*x37) - x32*x38;
auto x256 = x22*x46;
auto x47 = cell[8] + x16 + x239 - x26*x33 + x28*(x223*x26 + x224*x32);
auto x264 = x22*x47;
auto x48 = V{2}*V{1} / (x43 + util::sqrt(V{2.828428}*x235*x242*((x221)*(x221))*util::sqrt(V{2}*x45 + ((x30)*(x30)) + ((x35)*(x35)) + ((x40)*(x40)) + V{2}*((x46)*(x46)) + V{2}*((x47)*(x47))) + ((x43)*(x43))));
auto x49 = V{1} - x48;
auto x50 = -x222*x247;
auto x51 = -x223*x248;
auto x52 = x50 + x51;
auto x53 = -x224*x249;
auto x269 = x51 + x53;
auto x272 = x50 + x53;
auto x54 = -x5*(V{0.5}*x48 + V{-1});
auto x55 = x0*x22;
auto x56 = V{0.166667}*x55;
auto x265 = x25*x56;
auto x57 = -x243;
auto x58 = -x57;
auto x250 = x56*x58;
auto x59 = -x250;
auto x60 = x13*x8;
auto x61 = x240*x7;
auto x62 = x60 - x61;
auto x253 = x56*x62;
auto x63 = V{0.5}*x55;
auto x251 = x58*x63;
auto x266 = x25*x63;
auto x64 = x251;
auto x65 = x60 + x61;
auto x254 = x63*x65;
auto x66 = x59*x62 + x64*x65;
auto x67 = x241*x7;
auto x68 = x13*x241;
auto x69 = x67 - x68;
auto x261 = x56*x69;
auto x278 = x261*(x25 - x58);
auto x70 = x67 + x68;
auto x258 = x63*x70;
auto x71 = x64*x70;
auto x291 = x270 + x71;
auto x72 = -V{1}*x7;
auto x73 = V{0.0555555555555556};
auto x267 = -x73;
auto x268 = x49*x73;
auto x74 = x5*x73;
auto x75 = -V{1}*x25;
auto x76 = 0;
auto x77 = V{2}*x18;
auto x78 = x30*x7 + x44*x77;
auto x79 = V{2}*x13;
auto x80 = x35*x7 + x46*x79;
auto x81 = x78 - x80;
auto x82 = V{2}*x63*x75;
auto x83 = x78 + x80;
auto x84 = x252;
auto x259 = -x84;
auto x260 = x56*(x25 + x259);
auto x282 = x261*x84;
auto x85 = x40*x84;
auto x86 = x36 + x85;
auto x87 = x23*(x31 + x86);
auto x88 = -x0*(x7 + V{1});
auto x290 = -x88;
auto x89 = x54*x73;
auto x90 = x75*(x253 + x254);
auto x91 = -V{1}*x13;
auto x92 = V{2}*x7;
auto x93 = x13*x40 + x46*x92;
auto x94 = x13*x30 + x47*x77;
auto x95 = x93 - x94;
auto x96 = x93 + x94;
auto x97 = x35*x84;
auto x98 = x41 + x97;
auto x99 = x23*(x31 + x98);
auto x100 = -x0*(x13 + V{1});
auto x101 = x14*x18;
auto x102 = x18*x8;
auto x103 = x101 - x102;
auto x104 = x101 + x102;
auto x105 = -V{1}*x18;
auto x106 = x18*x35 + x47*x79;
auto x107 = x18*x40 + x44*x92;
auto x108 = x106 - x107;
auto x109 = x106 + x107;
auto x110 = x30*x84;
auto x111 = x23*(x110 + x42);
auto x112 = -x0*(x18 + V{1});
auto x113 = -V{1}*x84;
auto x255 = x56*(x113 - x75);
auto x114 = -x57 + V{-1};
auto x262 = x114*x261;
auto x115 = x72 + x91;
auto x116 = V{0.0277777777777778};
auto x273 = -x116;
auto x276 = x116*x49;
auto x117 = x116*x5;
auto x118 = -x255;
auto x119 = x118*x81;
auto x120 = x255;
auto x121 = x120*x95;
auto x122 = x63*(x113 + x75);
auto x123 = x122*x83;
auto x292 = x119 + x123;
auto x124 = x122*x96;
auto x277 = x121 + x124;
auto x125 = x23*(x31 + x85 + x97);
auto x126 = x125 + x256;
auto x127 = -x115*x257;
auto x128 = x116*x54;
auto x129 = x258*(x25 + x58);
auto x279 = x129 + x270;
auto x130 = x13 + x72;
auto x131 = -x260;
auto x132 = x131*x95;
auto x133 = x63*(x25 + x84);
auto x134 = x133*x96;
auto x275 = x132 + x134;
auto x135 = x125 - V{1}*x256;
auto x136 = x130*x22;
auto x137 = -x0*(x13 + V{-1});
auto x138 = x103*x118;
auto x139 = x262;
auto x295 = x139 + x66;
auto x140 = x104*x122;
auto x141 = x114*x258;
auto x142 = x140 + x141;
auto x294 = x142 + x270;
auto x143 = x105 + x72;
auto x144 = x108*x118;
auto x145 = x120*x81;
auto x293 = x123 + x145;
auto x146 = x109*x122;
auto x284 = x144 + x146;
auto x147 = x23*(x110 + x86);
auto x148 = x147 + x263;
auto x149 = -x143*x257;
auto x150 = x260;
auto x151 = x103*x150;
auto x152 = x104*x133;
auto x153 = x141 + x152;
auto x296 = x153 + x270;
auto x154 = x18 + x72;
auto x155 = x108*x150;
auto x156 = x109*x133;
auto x281 = x155 + x156;
auto x157 = x147 - V{1}*x263;
auto x158 = x154*x22;
auto x159 = -x0*(x18 + V{-1});
auto x160 = x103*x120;
auto x161 = -x262;
auto x298 = x161 + x90;
auto x162 = x105 + x91;
auto x163 = x108*x120;
auto x289 = x146 + x163;
auto x164 = x118*x95;
auto x297 = x124 + x164;
auto x165 = x23*(x110 + x98);
auto x166 = x165 + x264;
auto x167 = -x162*x257;
auto x168 = x103*x131;
auto x169 = x18 + x91;
auto x170 = x108*x131;
auto x287 = x156 + x170;
auto x171 = x165 - V{1}*x264;
auto x172 = x169*x22;
auto x173 = -x265;
auto x174 = x266;
auto x175 = x174*x65;
auto x176 = x173*x62 + x175;
auto x283 = x176 + x270;
auto x177 = x265;
auto x178 = x174*x70;
auto x271 = x178 + x270;
auto x179 = 0;
auto x180 = V{2}*x266;
auto x181 = -x0*(x7 + V{-1});
auto x182 = x175 + x177*x62;
auto x288 = x182 + x270;
auto x183 = x13 + x7;
auto x184 = x150*x81;
auto x185 = x133*x83;
auto x274 = x184 + x185;
auto x186 = x183*x22;
auto x187 = x252;
auto x285 = x187*x261;
auto x188 = x7 + x91;
auto x189 = x188*x22;
auto x190 = x152 + x258*x84;
auto x191 = x18 + x7;
auto x192 = x131*x81;
auto x280 = x185 + x192;
auto x193 = x191*x22;
auto x194 = x140 + x187*x258;
auto x195 = x105 + x7;
auto x196 = x195*x22;
auto x197 = x13 + x18;
auto x198 = x150*x95;
auto x286 = x134 + x198;
auto x199 = x197*x22;
auto x200 = x105 + x13;
auto x201 = x200*x22;
cell[0] = x21*(x23*x49*(x31 + x42) - x5*(x20 + V{-1}) + x54*(x52 + x53) + V{-1});
cell[10] = x267 + x268*(x179*x81 + x180*x83 + x87) + x74*(x176 + x177*x69 + x23*x238 + x247 + x271) + x89*(x222*(x181 + x22*x7) + x269);
cell[11] = x267 + x268*(x179*x95 + x180*x96 + x99) + x74*(x173*x69 + x182 + x23*x240 + x248 + x271) + x89*(x223*(x13*x22 + x137) + x272);
cell[12] = x267 + x268*(x108*x179 + x109*x180 + x111) + x74*(x103*x179 + x104*x180 + x23*x241 + x249 + x270) + x89*(x224*(x159 + x18*x22) + x52);
cell[13] = x117*(x0*x183 + x179*x69 + x180*x70 + x23*((x183)*(x183)) + V{2}*x254*x84 + x270) + x128*(x222*(x181 + x186) + x223*(x137 + x186) + x53) + x273 + x276*(x126 + x274 + x275);
cell[14] = x117*(x0*x188 + x23*((x188)*(x188)) + x253*(x113 - x187) + x254*(x113 + x187) + x278 + x279) + x128*(x222*(x181 + x189) + x223*(x100 - V{1}*x189) + x53) + x273 + x276*(x135 + x274 + x277);
cell[15] = x117*(x0*x191 + x151 + x190 + x23*((x191)*(x191)) + x282 + x283) + x128*(x222*(x181 + x193) + x224*(x159 + x193) + x51) + x273 + x276*(x148 + x280 + x281);
cell[16] = x117*(x0*x195 + x138 + x194 + x23*((x195)*(x195)) + x283 + x285) + x128*(x222*(x181 + x196) + x224*(x112 - V{1}*x196) + x51) + x273 + x276*(x157 + x280 + x284);
cell[17] = x117*(x0*x197 + x168 + x190 + x23*((x197)*(x197)) - x282 + x288) + x128*(x223*(x137 + x199) + x224*(x159 + x199) + x50) + x273 + x276*(x166 + x286 + x287);
cell[18] = x117*(x0*x200 + x160 + x194 + x23*((x200)*(x200)) - x285 + x288) + x128*(x223*(x137 + x201) + x224*(x112 - V{1}*x201) + x50) + x273 + x276*(x171 + x286 + x289);
cell[1] = x267 + x268*(x76*x81 + x82*x83 + x87) + x74*(x0*x72 + x23*((x72)*(x72)) + x250*x69 + x291 + x66) + x89*(-x222*(x257*x72 + x290) + x269);
cell[2] = x267 + x268*(x76*x95 + x82*x96 + x99) + x74*(x0*x91 + x23*((x91)*(x91)) + x291 + x59*x69 + x90) + x89*(x223*(x100 - x257*x91) + x272);
cell[3] = x267 + x268*(x108*x76 + x109*x82 + x111) + x74*(x0*x105 + x103*x76 + x104*x82 + x23*((x105)*(x105)) + x270) - x89*(x224*(x105*x257 - x112) - x52);
cell[4] = x117*(x0*x115 + x23*((x115)*(x115)) + V{2}*x251*x70 + x253*(x113 - x114) + x254*(x113 + x114) + x270) + x128*(x222*(x127 + x88) + x223*(x100 + x127) + x53) + x273 + x276*(x126 + x277 + x292);
cell[5] = x117*(x0*x130 + x23*((x130)*(x130)) - x253*(x114 + x259) + x254*(x114 + x84) - x278 + x279) + x128*(-x222*(V{1}*x136 + x290) + x223*(x136 + x137) + x53) + x273 + x276*(x135 + x275 + x292);
cell[6] = x117*(x0*x143 + x138 + x23*((x143)*(x143)) + x294 + x295) + x128*(x222*(x149 + x88) + x224*(x112 + x149) + x51) + x273 + x276*(x148 + x284 + x293);
cell[7] = x117*(x0*x154 + x151 + x23*((x154)*(x154)) + x295 + x296) + x128*(-x222*(V{1}*x158 + x290) + x224*(x158 + x159) + x51) + x273 + x276*(x157 + x281 + x293);
cell[8] = x117*(x0*x162 + x160 + x23*((x162)*(x162)) + x294 + x298) + x128*(x223*(x100 + x167) + x224*(x112 + x167) + x50) + x273 + x276*(x166 + x289 + x297);
cell[9] = x117*(x0*x169 + x168 + x23*((x169)*(x169)) + x296 + x298) + x128*(x223*(x100 - V{1}*x172) + x224*(x159 + x172) + x50) + x273 + x276*(x171 + x287 + x297);
return { x5, x19 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
