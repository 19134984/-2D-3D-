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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::ParameterFromCell<collision::LES::SMAGORINSKY, collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB> > >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x216 = cell.template getFieldComponent<olb::collision::LES::SMAGORINSKY>(0);
auto x217 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x218 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x219 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x220 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x221 = cell.template getFieldComponent<olb::descriptors::WMPOROSITY>(0);
auto x222 = cell.template getFieldComponent<olb::descriptors::WMVELOCITY>(0);
auto x223 = cell.template getFieldComponent<olb::descriptors::WMVELOCITY>(1);
auto x224 = cell.template getFieldComponent<olb::descriptors::WMVELOCITY>(2);
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
auto x243 = ((x0)*(x0));
auto x244 = V{1}/x0;
auto x1 = V{1} - x221;
auto x2 = x1*x222;
auto x3 = V{1} - x217;
auto x4 = x218*x3;
auto x5 = -V{1}*cell[4];
auto x6 = -V{1}*cell[5];
auto x236 = cell[13] + x6;
auto x7 = -V{1}*cell[6];
auto x8 = -V{1}*cell[7];
auto x235 = cell[15] + x8;
auto x9 = cell[0] + cell[8] + x229 + x232 + x234 + V{1};
auto x237 = V{1} / (x9);
auto x238 = x217*x237;
auto x10 = x2 + x221*(x238*(-V{1}*cell[1] + x233 + x235 + x236 + x5 + x7) + x4);
auto x239 = ((x10)*(x10));
auto x11 = x239;
auto x12 = x1*x223;
auto x13 = x219*x3;
auto x14 = -V{1}*cell[8];
auto x15 = -V{1}*cell[9];
auto x240 = cell[17] + x15;
auto x16 = -V{1}*cell[14];
auto x17 = x12 + x221*(x13 + x238*(-V{1}*cell[2] + x14 + x16 + x231 + x240 + x5));
auto x241 = ((x17)*(x17));
auto x18 = x241;
auto x19 = x1*x224;
auto x20 = x220*x3;
auto x21 = -V{1}*cell[16];
auto x22 = -V{1}*cell[18];
auto x23 = x19 + x221*(x20 + x238*(-V{1}*cell[3] + x14 + x21 + x22 + x228 + x7));
auto x242 = ((x23)*(x23));
auto x24 = x11 + x18 + x242;
auto x25 = V{0.5}*x0*x24;
auto x252 = V{1} - x25;
auto x26 = V{0.333333333333333};
auto x27 = x243;
auto x28 = V{0.5}*x27;
auto x29 = x244;
auto x249 = V{1} - x29;
auto x30 = -x29;
auto x31 = x29*(x9 + V{-1});
auto x245 = -x31;
auto x246 = cell[8] + x245;
auto x32 = cell[16] + cell[18] + x229 - x242*x9 + x246;
auto x33 = x30*x32;
auto x34 = x17*x9;
auto x35 = cell[14] - x17*x34 + x225 + x232 + x246;
auto x36 = x30*x35;
auto x37 = x10*x9;
auto x38 = cell[4] + cell[6] - x10*x37 + x227 + x230 + x234 + x245;
auto x39 = x30*x38;
auto x40 = x36 + x39;
auto x41 = V{1}/x226;
auto x42 = cell[6] + x21 - x23*x37 + x235;
auto x270 = x27*x42;
auto x43 = ((x42)*(x42));
auto x44 = cell[4] + x16 - x17*x37 + x236;
auto x257 = x27*x44;
auto x45 = cell[8] + x22 - x23*x34 + x240;
auto x277 = x27*x45;
auto x46 = V{1} - V{2}*V{1} / (x41 + util::sqrt(V{2.828428}*x237*x243*((x216)*(x216))*util::sqrt(V{2}*x43 + ((x32)*(x32)) + ((x35)*(x35)) + ((x38)*(x38)) + V{2}*((x44)*(x44)) + V{2}*((x45)*(x45))) + ((x41)*(x41))));
auto x47 = x26*(x28*x46*(x33 + x40) - x9*(x25 + V{-1}) + V{-1});
auto x48 = x0*x27;
auto x49 = V{0.166667}*x48;
auto x282 = x30*x49;
auto x50 = -x244;
auto x51 = -x50;
auto x247 = x49*x51;
auto x52 = -x247;
auto x53 = x11*x17;
auto x54 = x10*x241;
auto x55 = x53 - x54;
auto x254 = x49*x55;
auto x56 = V{0.5}*x48;
auto x248 = x51*x56;
auto x283 = x30*x56;
auto x57 = x248;
auto x58 = x53 + x54;
auto x255 = x56*x58;
auto x59 = x52*x55 + x57*x58;
auto x60 = x10*x242;
auto x61 = x17*x242;
auto x62 = x60 - x61;
auto x266 = x49*x62;
auto x267 = x266*(x30 - x51);
auto x63 = x60 + x61;
auto x262 = x56*x63;
auto x64 = x57*x63;
auto x253 = x252 + x64;
auto x65 = -V{1}*x10;
auto x66 = V{0.0555555555555556};
auto x250 = -x66;
auto x251 = x46*x66;
auto x67 = x66*x9;
auto x68 = -V{1}*x30;
auto x69 = 0;
auto x70 = V{2}*x23;
auto x71 = x10*x32 + x42*x70;
auto x72 = V{2}*x17;
auto x73 = x10*x35 + x44*x72;
auto x74 = x71 - x73;
auto x75 = V{2}*x56*x68;
auto x76 = x71 + x73;
auto x77 = x249;
auto x263 = -x77;
auto x264 = x49*(x263 + x30);
auto x287 = x266*x77;
auto x78 = x38*x77;
auto x79 = x36 + x78;
auto x80 = x28*(x33 + x79);
auto x81 = x250 + x251*(x69*x74 + x75*x76 + x80) + x67*(x0*x65 + x247*x62 + x253 + x28*((x65)*(x65)) + x59);
auto x82 = x68*(x254 + x255);
auto x83 = -V{1}*x17;
auto x84 = V{2}*x10;
auto x85 = x17*x38 + x44*x84;
auto x86 = x17*x32 + x45*x70;
auto x87 = x85 - x86;
auto x88 = x85 + x86;
auto x89 = x35*x77;
auto x90 = x39 + x89;
auto x91 = x28*(x33 + x90);
auto x92 = x250 + x251*(x69*x87 + x75*x88 + x91) + x67*(x0*x83 + x253 + x28*((x83)*(x83)) + x52*x62 + x82);
auto x93 = x18*x23;
auto x94 = x11*x23;
auto x95 = x93 - x94;
auto x96 = x93 + x94;
auto x97 = -V{1}*x23;
auto x98 = x23*x35 + x45*x72;
auto x99 = x23*x38 + x42*x84;
auto x100 = x98 - x99;
auto x101 = x98 + x99;
auto x102 = x32*x77;
auto x103 = x28*(x102 + x40);
auto x104 = x250 + x251*(x100*x69 + x101*x75 + x103) + x67*(x0*x97 + x252 + x28*((x97)*(x97)) + x69*x95 + x75*x96);
auto x105 = -V{1}*x77;
auto x256 = x49*(x105 - x68);
auto x106 = -x50 + V{-1};
auto x269 = x106*x266;
auto x107 = x65 + x83;
auto x108 = V{0.0277777777777778};
auto x258 = -x108;
auto x261 = x108*x46;
auto x109 = x108*x9;
auto x110 = -x256;
auto x111 = x110*x74;
auto x112 = x256;
auto x113 = x112*x87;
auto x114 = x56*(x105 + x68);
auto x115 = x114*x76;
auto x259 = x111 + x115;
auto x116 = x114*x88;
auto x260 = x113 + x116;
auto x117 = x28*(x33 + x78 + x89);
auto x118 = x117 + x257;
auto x119 = x109*(x0*x107 + V{2}*x248*x63 + x252 + x254*(x105 - x106) + x255*(x105 + x106) + x28*((x107)*(x107))) + x258 + x261*(x118 + x259 + x260);
auto x120 = x262*(x30 + x51);
auto x268 = x120 + x252;
auto x121 = x17 + x65;
auto x122 = -x264;
auto x123 = x122*x87;
auto x124 = x56*(x30 + x77);
auto x125 = x124*x88;
auto x265 = x123 + x125;
auto x126 = x117 - V{1}*x257;
auto x127 = x109*(x0*x121 - x254*(x106 + x263) + x255*(x106 + x77) - x267 + x268 + x28*((x121)*(x121))) + x258 + x261*(x126 + x259 + x265);
auto x128 = x110*x95;
auto x129 = x269;
auto x274 = x129 + x59;
auto x130 = x114*x96;
auto x131 = x106*x262;
auto x132 = x130 + x131;
auto x273 = x132 + x252;
auto x133 = x65 + x97;
auto x134 = x100*x110;
auto x135 = x112*x74;
auto x271 = x115 + x135;
auto x136 = x101*x114;
auto x272 = x134 + x136;
auto x137 = x28*(x102 + x79);
auto x138 = x137 + x270;
auto x139 = x109*(x0*x133 + x128 + x273 + x274 + x28*((x133)*(x133))) + x258 + x261*(x138 + x271 + x272);
auto x140 = x264;
auto x141 = x140*x95;
auto x142 = x124*x96;
auto x143 = x131 + x142;
auto x276 = x143 + x252;
auto x144 = x23 + x65;
auto x145 = x100*x140;
auto x146 = x101*x124;
auto x275 = x145 + x146;
auto x147 = x137 - V{1}*x270;
auto x148 = x109*(x0*x144 + x141 + x274 + x276 + x28*((x144)*(x144))) + x258 + x261*(x147 + x271 + x275);
auto x149 = x112*x95;
auto x150 = -x269;
auto x280 = x150 + x82;
auto x151 = x83 + x97;
auto x152 = x100*x112;
auto x279 = x136 + x152;
auto x153 = x110*x87;
auto x278 = x116 + x153;
auto x154 = x28*(x102 + x90);
auto x155 = x154 + x277;
auto x156 = x109*(x0*x151 + x149 + x273 + x28*((x151)*(x151)) + x280) + x258 + x261*(x155 + x278 + x279);
auto x157 = x122*x95;
auto x158 = x23 + x83;
auto x159 = x100*x122;
auto x281 = x146 + x159;
auto x160 = x154 - V{1}*x277;
auto x161 = x109*(x0*x158 + x157 + x276 + x28*((x158)*(x158)) + x280) + x258 + x261*(x160 + x278 + x281);
auto x162 = -x282;
auto x163 = x283;
auto x164 = x163*x58;
auto x165 = x162*x55 + x164;
auto x288 = x165 + x252;
auto x166 = x282;
auto x167 = x163*x63;
auto x284 = x167 + x252;
auto x168 = 0;
auto x169 = V{2}*x283;
auto x170 = x250 + x251*(x168*x74 + x169*x76 + x80) + x67*(x0*x10 + x165 + x166*x62 + x239*x28 + x284);
auto x171 = x164 + x166*x55;
auto x291 = x171 + x252;
auto x172 = x250 + x251*(x168*x87 + x169*x88 + x91) + x67*(x0*x17 + x162*x62 + x171 + x241*x28 + x284);
auto x173 = x250 + x251*(x100*x168 + x101*x169 + x103) + x67*(x0*x23 + x168*x95 + x169*x96 + x242*x28 + x252);
auto x294 = x148 + x161 + x173;
auto x174 = x10 + x17;
auto x175 = x140*x74;
auto x176 = x124*x76;
auto x285 = x175 + x176;
auto x177 = x109*(x0*x174 + x168*x62 + x169*x63 + x252 + V{2}*x255*x77 + x28*((x174)*(x174))) + x258 + x261*(x118 + x265 + x285);
auto x178 = x249;
auto x289 = x178*x266;
auto x179 = x10 + x83;
auto x180 = x109*(x0*x179 + x254*(x105 - x178) + x255*(x105 + x178) + x267 + x268 + x28*((x179)*(x179))) + x258 + x261*(x126 + x260 + x285);
auto x181 = x142 + x262*x77;
auto x182 = x10 + x23;
auto x183 = x122*x74;
auto x286 = x176 + x183;
auto x184 = x109*(x0*x182 + x141 + x181 + x28*((x182)*(x182)) + x287 + x288) + x258 + x261*(x138 + x275 + x286);
auto x185 = x130 + x178*x262;
auto x186 = x10 + x97;
auto x187 = x109*(x0*x186 + x128 + x185 + x28*((x186)*(x186)) + x288 + x289) + x258 + x261*(x147 + x272 + x286);
auto x292 = x170 + x180 + x187;
auto x188 = x17 + x23;
auto x189 = x140*x87;
auto x290 = x125 + x189;
auto x190 = x109*(x0*x188 + x157 + x181 + x28*((x188)*(x188)) - x287 + x291) + x258 + x261*(x155 + x281 + x290);
auto x191 = x17 + x97;
auto x192 = x109*(x0*x191 + x149 + x185 + x28*((x191)*(x191)) - x289 + x291) + x258 + x261*(x160 + x279 + x290);
auto x293 = x127 + x172 + x192;
auto x193 = x104 + x119 + x139 + x156 + x177 + x184 + x190 + x292 + x293 + x294 + x47 + x81 + x92 + V{1};
auto x297 = x217/x193;
auto x194 = -V{1}*x119;
auto x295 = x177 + x194;
auto x195 = -V{1}*x139;
auto x296 = x184 + x195;
auto x196 = -V{1}*x156;
auto x298 = x190 + x196;
cell[0] = x47;
cell[10] = x170;
cell[11] = x172;
cell[12] = x173;
cell[13] = x177;
cell[14] = x180;
cell[15] = x184;
cell[16] = x187;
cell[17] = x190;
cell[18] = x192;
cell[1] = x81;
cell[2] = x92;
cell[3] = x104;
cell[4] = x119;
cell[5] = x127;
cell[6] = x139;
cell[7] = x148;
cell[8] = x156;
cell[9] = x161;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x193;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x2 + x221*(x297*(-V{1}*x127 - V{1}*x148 + x292 + x295 + x296 - V{1}*x81) + x4);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x12 + x221*(x13 + x297*(-V{1}*x161 - V{1}*x180 + x293 + x295 + x298 - V{1}*x92));
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x19 + x221*(x20 + x297*(-V{1}*x104 - V{1}*x187 - V{1}*x192 + x294 + x296 + x298));
return { x9, x24 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
