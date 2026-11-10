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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::ParameterFromCell<collision::LES::SMAGORINSKY, collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB> > >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x187 = cell.template getFieldComponent<olb::collision::LES::SMAGORINSKY>(0);
auto x188 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x189 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x190 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x191 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x192 = cell.template getFieldComponent<olb::descriptors::WMPOROSITY>(0);
auto x193 = cell.template getFieldComponent<olb::descriptors::WMVELOCITY>(0);
auto x194 = cell.template getFieldComponent<olb::descriptors::WMVELOCITY>(1);
auto x195 = cell.template getFieldComponent<olb::descriptors::WMVELOCITY>(2);
auto x197 = parameters.template get<descriptors::OMEGA>();
auto x196 = cell[17] + cell[9];
auto x198 = cell[15] + cell[7];
auto x199 = cell[12] + x196 + x198;
auto x200 = cell[3] + cell[6] + x199;
auto x201 = cell[13] + cell[5];
auto x202 = cell[11] + cell[18] + x201;
auto x203 = cell[2] + cell[4] + x202;
auto x204 = cell[10] + cell[14] + cell[16];
auto x205 = cell[1] + x204;
auto x0 = V{3};
auto x214 = ((x0)*(x0));
auto x215 = V{1}/x0;
auto x1 = V{1} - x192;
auto x2 = V{1} - x188;
auto x3 = -V{1}*cell[4];
auto x4 = -V{1}*cell[5];
auto x207 = cell[13] + x4;
auto x5 = -V{1}*cell[6];
auto x6 = -V{1}*cell[7];
auto x206 = cell[15] + x6;
auto x7 = cell[0] + cell[8] + x200 + x203 + x205 + V{1};
auto x208 = V{1} / (x7);
auto x209 = x188*x208;
auto x8 = x1*x193 + x192*(x189*x2 + x209*(-V{1}*cell[1] + x204 + x206 + x207 + x3 + x5));
auto x210 = ((x8)*(x8));
auto x9 = x210;
auto x10 = -V{1}*cell[8];
auto x11 = -V{1}*cell[9];
auto x211 = cell[17] + x11;
auto x12 = -V{1}*cell[14];
auto x13 = x1*x194 + x192*(x190*x2 + x209*(-V{1}*cell[2] + x10 + x12 + x202 + x211 + x3));
auto x212 = ((x13)*(x13));
auto x14 = x212;
auto x15 = -V{1}*cell[16];
auto x16 = -V{1}*cell[18];
auto x17 = x1*x195 + x192*(x191*x2 + x209*(-V{1}*cell[3] + x10 + x15 + x16 + x199 + x5));
auto x213 = ((x17)*(x17));
auto x18 = x14 + x213 + x9;
auto x19 = V{0.5}*x0*x18;
auto x236 = V{1} - x19;
auto x20 = V{0.333333333333333};
auto x21 = x214;
auto x22 = V{0.5}*x21;
auto x23 = x215;
auto x220 = V{1} - x23;
auto x24 = -x23;
auto x25 = x23*(x7 + V{-1});
auto x216 = -x25;
auto x217 = cell[8] + x216;
auto x26 = cell[16] + cell[18] + x200 - x213*x7 + x217;
auto x27 = x24*x26;
auto x28 = x13*x7;
auto x29 = cell[14] - x13*x28 + x196 + x203 + x217;
auto x30 = x24*x29;
auto x31 = x7*x8;
auto x32 = cell[4] + cell[6] + x198 + x201 + x205 + x216 - x31*x8;
auto x33 = x24*x32;
auto x34 = x30 + x33;
auto x35 = V{1}/x197;
auto x36 = cell[6] + x15 - x17*x31 + x206;
auto x230 = x21*x36;
auto x37 = ((x36)*(x36));
auto x38 = cell[4] + x12 - x13*x31 + x207;
auto x225 = x21*x38;
auto x39 = cell[8] + x16 - x17*x28 + x211;
auto x231 = x21*x39;
auto x40 = V{1} - V{2}*V{1} / (x35 + util::sqrt(V{2.828428}*x208*x214*((x187)*(x187))*util::sqrt(V{2}*x37 + ((x26)*(x26)) + ((x29)*(x29)) + ((x32)*(x32)) + V{2}*((x38)*(x38)) + V{2}*((x39)*(x39))) + ((x35)*(x35))));
auto x41 = x0*x21;
auto x42 = V{0.166667}*x41;
auto x232 = x24*x42;
auto x43 = -x215;
auto x44 = -x43;
auto x218 = x42*x44;
auto x45 = -x218;
auto x46 = x13*x9;
auto x47 = x212*x8;
auto x48 = x46 - x47;
auto x221 = x42*x48;
auto x49 = V{0.5}*x41;
auto x219 = x44*x49;
auto x233 = x24*x49;
auto x50 = x219;
auto x51 = x46 + x47;
auto x222 = x49*x51;
auto x52 = x45*x48 + x50*x51;
auto x53 = x213*x8;
auto x54 = x13*x213;
auto x55 = x53 - x54;
auto x228 = x42*x55;
auto x243 = x228*(x24 - x44);
auto x56 = x53 + x54;
auto x226 = x49*x56;
auto x57 = x50*x56;
auto x255 = x236 + x57;
auto x58 = -V{1}*x8;
auto x59 = V{0.0555555555555556};
auto x234 = -x59;
auto x235 = x40*x59;
auto x60 = x59*x7;
auto x61 = -V{1}*x24;
auto x62 = 0;
auto x63 = V{2}*x17;
auto x64 = x26*x8 + x36*x63;
auto x65 = V{2}*x13;
auto x66 = x29*x8 + x38*x65;
auto x67 = x64 - x66;
auto x68 = V{2}*x49*x61;
auto x69 = x64 + x66;
auto x70 = x220;
auto x227 = x42*(x24 - x70);
auto x247 = x228*x70;
auto x71 = x32*x70;
auto x72 = x30 + x71;
auto x73 = x22*(x27 + x72);
auto x74 = x61*(x221 + x222);
auto x75 = -V{1}*x13;
auto x76 = V{2}*x8;
auto x77 = x13*x32 + x38*x76;
auto x78 = x13*x26 + x39*x63;
auto x79 = x77 - x78;
auto x80 = x77 + x78;
auto x81 = x29*x70;
auto x82 = x33 + x81;
auto x83 = x22*(x27 + x82);
auto x84 = x14*x17;
auto x85 = x17*x9;
auto x86 = x84 - x85;
auto x87 = x84 + x85;
auto x88 = -V{1}*x17;
auto x89 = x17*x29 + x39*x65;
auto x90 = x17*x32 + x36*x76;
auto x91 = x89 - x90;
auto x92 = x89 + x90;
auto x93 = x26*x70;
auto x94 = x22*(x34 + x93);
auto x95 = -V{1}*x70;
auto x223 = -x95;
auto x224 = x42*(x223 + x61);
auto x96 = -x43 + V{-1};
auto x229 = x228*x96;
auto x257 = -x96;
auto x97 = x58 + x75;
auto x98 = V{0.0277777777777778};
auto x238 = -x98;
auto x241 = x40*x98;
auto x99 = x7*x98;
auto x100 = x224;
auto x101 = x100*x67;
auto x102 = -x224;
auto x103 = x102*x79;
auto x104 = x49*(x61 + x95);
auto x105 = x104*x69;
auto x256 = x101 + x105;
auto x106 = x104*x80;
auto x242 = x103 + x106;
auto x107 = x22*(x27 + x71 + x81);
auto x108 = x107 + x225;
auto x109 = x226*(x24 + x44);
auto x244 = x109 + x236;
auto x110 = x13 + x58;
auto x111 = -x227;
auto x112 = x111*x79;
auto x113 = x49*(x24 + x70);
auto x114 = x113*x80;
auto x240 = x112 + x114;
auto x115 = x107 - V{1}*x225;
auto x116 = x100*x86;
auto x117 = x229;
auto x260 = x117 + x52;
auto x118 = x104*x87;
auto x119 = x226*x96;
auto x120 = x118 + x119;
auto x259 = x120 + x236;
auto x121 = x58 + x88;
auto x122 = x100*x91;
auto x123 = x102*x67;
auto x258 = x105 + x123;
auto x124 = x104*x92;
auto x249 = x122 + x124;
auto x125 = x22*(x72 + x93);
auto x126 = x125 + x230;
auto x127 = x227;
auto x128 = x127*x86;
auto x129 = x113*x87;
auto x130 = x119 + x129;
auto x261 = x130 + x236;
auto x131 = x17 + x58;
auto x132 = x127*x91;
auto x133 = x113*x92;
auto x246 = x132 + x133;
auto x134 = x125 - V{1}*x230;
auto x135 = x102*x86;
auto x136 = -x229;
auto x263 = x136 + x74;
auto x137 = x75 + x88;
auto x138 = x102*x91;
auto x254 = x124 + x138;
auto x139 = x100*x79;
auto x262 = x106 + x139;
auto x140 = x22*(x82 + x93);
auto x141 = x140 + x231;
auto x142 = x111*x86;
auto x143 = x17 + x75;
auto x144 = x111*x91;
auto x252 = x133 + x144;
auto x145 = x140 - V{1}*x231;
auto x146 = -x232;
auto x147 = x233;
auto x148 = x147*x51;
auto x149 = x146*x48 + x148;
auto x248 = x149 + x236;
auto x150 = x232;
auto x151 = x147*x56;
auto x237 = x151 + x236;
auto x152 = 0;
auto x153 = V{2}*x233;
auto x154 = x148 + x150*x48;
auto x253 = x154 + x236;
auto x155 = x13 + x8;
auto x156 = x127*x67;
auto x157 = x113*x69;
auto x239 = x156 + x157;
auto x158 = x220;
auto x250 = x158*x228;
auto x159 = x75 + x8;
auto x160 = x129 + x226*x70;
auto x161 = x17 + x8;
auto x162 = x111*x67;
auto x245 = x157 + x162;
auto x163 = x118 + x158*x226;
auto x164 = x8 + x88;
auto x165 = x13 + x17;
auto x166 = x127*x79;
auto x251 = x114 + x166;
auto x167 = x13 + x88;
cell[0] = x20*(x22*x40*(x27 + x34) - x7*(x19 + V{-1}) + V{-1});
cell[10] = x234 + x235*(x152*x67 + x153*x69 + x73) + x60*(x0*x8 + x149 + x150*x55 + x210*x22 + x237);
cell[11] = x234 + x235*(x152*x79 + x153*x80 + x83) + x60*(x0*x13 + x146*x55 + x154 + x212*x22 + x237);
cell[12] = x234 + x235*(x152*x91 + x153*x92 + x94) + x60*(x0*x17 + x152*x86 + x153*x87 + x213*x22 + x236);
cell[13] = x238 + x241*(x108 + x239 + x240) + x99*(x0*x155 + x152*x55 + x153*x56 + x22*((x155)*(x155)) + V{2}*x222*x70 + x236);
cell[14] = x238 + x241*(x115 + x239 + x242) + x99*(x0*x159 + x22*((x159)*(x159)) - x221*(x158 + x223) + x222*(x158 + x95) + x243 + x244);
cell[15] = x238 + x241*(x126 + x245 + x246) + x99*(x0*x161 + x128 + x160 + x22*((x161)*(x161)) + x247 + x248);
cell[16] = x238 + x241*(x134 + x245 + x249) + x99*(x0*x164 + x116 + x163 + x22*((x164)*(x164)) + x248 + x250);
cell[17] = x238 + x241*(x141 + x251 + x252) + x99*(x0*x165 + x142 + x160 + x22*((x165)*(x165)) - x247 + x253);
cell[18] = x238 + x241*(x145 + x251 + x254) + x99*(x0*x167 + x135 + x163 + x22*((x167)*(x167)) - x250 + x253);
cell[1] = x234 + x235*(x62*x67 + x68*x69 + x73) + x60*(x0*x58 + x218*x55 + x22*((x58)*(x58)) + x255 + x52);
cell[2] = x234 + x235*(x62*x79 + x68*x80 + x83) + x60*(x0*x75 + x22*((x75)*(x75)) + x255 + x45*x55 + x74);
cell[3] = x234 + x235*(x62*x91 + x68*x92 + x94) + x60*(x0*x88 + x22*((x88)*(x88)) + x236 + x62*x86 + x68*x87);
cell[4] = x238 + x241*(x108 + x242 + x256) + x99*(x0*x97 + V{2}*x219*x56 + x22*((x97)*(x97)) + x221*(x257 + x95) + x222*(x95 + x96) + x236);
cell[5] = x238 + x241*(x115 + x240 + x256) + x99*(x0*x110 + x22*((x110)*(x110)) + x221*(x257 + x70) + x222*(x70 + x96) - x243 + x244);
cell[6] = x238 + x241*(x126 + x249 + x258) + x99*(x0*x121 + x116 + x22*((x121)*(x121)) + x259 + x260);
cell[7] = x238 + x241*(x134 + x246 + x258) + x99*(x0*x131 + x128 + x22*((x131)*(x131)) + x260 + x261);
cell[8] = x238 + x241*(x141 + x254 + x262) + x99*(x0*x137 + x135 + x22*((x137)*(x137)) + x259 + x263);
cell[9] = x238 + x241*(x145 + x252 + x262) + x99*(x0*x143 + x142 + x22*((x143)*(x143)) + x261 + x263);
return { x7, x18 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
