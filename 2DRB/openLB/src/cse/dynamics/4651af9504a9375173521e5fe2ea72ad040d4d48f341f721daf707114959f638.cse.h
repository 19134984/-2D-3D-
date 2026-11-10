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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::MovingPorousMomentumCombinationNoWM<momenta::BulkMomentum>, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x186 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x187 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x188 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x189 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x190 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x191 = parameters.template get<descriptors::OMEGA>();
auto x192 = cell[17] + cell[9];
auto x193 = cell[15] + cell[7];
auto x194 = cell[12] + x192 + x193;
auto x195 = cell[3] + cell[6] + x194;
auto x196 = cell[13] + cell[5];
auto x197 = cell[11] + cell[18] + x196;
auto x198 = cell[2] + cell[4] + x197;
auto x199 = cell[10] + cell[14] + cell[16];
auto x200 = cell[1] + x199;
auto x0 = V{3};
auto x209 = ((x0)*(x0));
auto x210 = V{1}/x0;
auto x1 = V{1} - x186;
auto x2 = -V{1}*cell[4];
auto x3 = -V{1}*cell[5];
auto x202 = cell[13] + x3;
auto x4 = -V{1}*cell[6];
auto x5 = -V{1}*cell[7];
auto x201 = cell[15] + x5;
auto x6 = cell[0] + cell[8] + x195 + x198 + x200 + V{1};
auto x203 = V{1} / (x6);
auto x204 = x186*x203;
auto x7 = x1*x187 + x204*(-V{1}*cell[1] + x199 + x2 + x201 + x202 + x4);
auto x205 = ((x7)*(x7));
auto x8 = x205;
auto x9 = -V{1}*cell[8];
auto x10 = -V{1}*cell[9];
auto x206 = cell[17] + x10;
auto x11 = -V{1}*cell[14];
auto x12 = x1*x188 + x204*(-V{1}*cell[2] + x11 + x197 + x2 + x206 + x9);
auto x207 = ((x12)*(x12));
auto x13 = x207;
auto x14 = -V{1}*cell[16];
auto x15 = -V{1}*cell[18];
auto x16 = x1*x189 + x204*(-V{1}*cell[3] + x14 + x15 + x194 + x4 + x9);
auto x208 = ((x16)*(x16));
auto x17 = x13 + x208 + x8;
auto x18 = V{0.5}*x0*x17;
auto x231 = V{1} - x18;
auto x19 = V{0.333333333333333};
auto x20 = x209;
auto x21 = V{0.5}*x20;
auto x22 = x210;
auto x215 = V{1} - x22;
auto x23 = -x22;
auto x24 = x22*(x6 + V{-1});
auto x211 = -x24;
auto x212 = cell[8] + x211;
auto x25 = cell[16] + cell[18] + x195 - x208*x6 + x212;
auto x26 = x23*x25;
auto x27 = x12*x6;
auto x28 = cell[14] - x12*x27 + x192 + x198 + x212;
auto x29 = x23*x28;
auto x30 = x6*x7;
auto x31 = cell[4] + cell[6] + x193 + x196 + x200 + x211 - x30*x7;
auto x32 = x23*x31;
auto x33 = x29 + x32;
auto x34 = V{1}/x191;
auto x35 = cell[6] + x14 - x16*x30 + x201;
auto x225 = x20*x35;
auto x36 = ((x35)*(x35));
auto x37 = cell[4] + x11 - x12*x30 + x202;
auto x220 = x20*x37;
auto x38 = cell[8] + x15 - x16*x27 + x206;
auto x226 = x20*x38;
auto x39 = V{1} - V{2}*V{1} / (x34 + util::sqrt(V{2.828428}*x203*x209*((x190)*(x190))*util::sqrt(V{2}*x36 + ((x25)*(x25)) + ((x28)*(x28)) + ((x31)*(x31)) + V{2}*((x37)*(x37)) + V{2}*((x38)*(x38))) + ((x34)*(x34))));
auto x40 = x0*x20;
auto x41 = V{0.166667}*x40;
auto x227 = x23*x41;
auto x42 = -x210;
auto x43 = -x42;
auto x213 = x41*x43;
auto x44 = -x213;
auto x45 = x12*x8;
auto x46 = x207*x7;
auto x47 = x45 - x46;
auto x216 = x41*x47;
auto x48 = V{0.5}*x40;
auto x214 = x43*x48;
auto x228 = x23*x48;
auto x49 = x214;
auto x50 = x45 + x46;
auto x217 = x48*x50;
auto x51 = x44*x47 + x49*x50;
auto x52 = x208*x7;
auto x53 = x12*x208;
auto x54 = x52 - x53;
auto x223 = x41*x54;
auto x238 = x223*(x23 - x43);
auto x55 = x52 + x53;
auto x221 = x48*x55;
auto x56 = x49*x55;
auto x250 = x231 + x56;
auto x57 = -V{1}*x7;
auto x58 = V{0.0555555555555556};
auto x229 = -x58;
auto x230 = x39*x58;
auto x59 = x58*x6;
auto x60 = -V{1}*x23;
auto x61 = 0;
auto x62 = V{2}*x16;
auto x63 = x25*x7 + x35*x62;
auto x64 = V{2}*x12;
auto x65 = x28*x7 + x37*x64;
auto x66 = x63 - x65;
auto x67 = V{2}*x48*x60;
auto x68 = x63 + x65;
auto x69 = x215;
auto x222 = x41*(x23 - x69);
auto x242 = x223*x69;
auto x70 = x31*x69;
auto x71 = x29 + x70;
auto x72 = x21*(x26 + x71);
auto x73 = x60*(x216 + x217);
auto x74 = -V{1}*x12;
auto x75 = V{2}*x7;
auto x76 = x12*x31 + x37*x75;
auto x77 = x12*x25 + x38*x62;
auto x78 = x76 - x77;
auto x79 = x76 + x77;
auto x80 = x28*x69;
auto x81 = x32 + x80;
auto x82 = x21*(x26 + x81);
auto x83 = x13*x16;
auto x84 = x16*x8;
auto x85 = x83 - x84;
auto x86 = x83 + x84;
auto x87 = -V{1}*x16;
auto x88 = x16*x28 + x38*x64;
auto x89 = x16*x31 + x35*x75;
auto x90 = x88 - x89;
auto x91 = x88 + x89;
auto x92 = x25*x69;
auto x93 = x21*(x33 + x92);
auto x94 = -V{1}*x69;
auto x218 = -x94;
auto x219 = x41*(x218 + x60);
auto x95 = -x42 + V{-1};
auto x224 = x223*x95;
auto x252 = -x95;
auto x96 = x57 + x74;
auto x97 = V{0.0277777777777778};
auto x233 = -x97;
auto x236 = x39*x97;
auto x98 = x6*x97;
auto x99 = x219;
auto x100 = x66*x99;
auto x101 = -x219;
auto x102 = x101*x78;
auto x103 = x48*(x60 + x94);
auto x104 = x103*x68;
auto x251 = x100 + x104;
auto x105 = x103*x79;
auto x237 = x102 + x105;
auto x106 = x21*(x26 + x70 + x80);
auto x107 = x106 + x220;
auto x108 = x221*(x23 + x43);
auto x239 = x108 + x231;
auto x109 = x12 + x57;
auto x110 = -x222;
auto x111 = x110*x78;
auto x112 = x48*(x23 + x69);
auto x113 = x112*x79;
auto x235 = x111 + x113;
auto x114 = x106 - V{1}*x220;
auto x115 = x85*x99;
auto x116 = x224;
auto x255 = x116 + x51;
auto x117 = x103*x86;
auto x118 = x221*x95;
auto x119 = x117 + x118;
auto x254 = x119 + x231;
auto x120 = x57 + x87;
auto x121 = x90*x99;
auto x122 = x101*x66;
auto x253 = x104 + x122;
auto x123 = x103*x91;
auto x244 = x121 + x123;
auto x124 = x21*(x71 + x92);
auto x125 = x124 + x225;
auto x126 = x222;
auto x127 = x126*x85;
auto x128 = x112*x86;
auto x129 = x118 + x128;
auto x256 = x129 + x231;
auto x130 = x16 + x57;
auto x131 = x126*x90;
auto x132 = x112*x91;
auto x241 = x131 + x132;
auto x133 = x124 - V{1}*x225;
auto x134 = x101*x85;
auto x135 = -x224;
auto x258 = x135 + x73;
auto x136 = x74 + x87;
auto x137 = x101*x90;
auto x249 = x123 + x137;
auto x138 = x78*x99;
auto x257 = x105 + x138;
auto x139 = x21*(x81 + x92);
auto x140 = x139 + x226;
auto x141 = x110*x85;
auto x142 = x16 + x74;
auto x143 = x110*x90;
auto x247 = x132 + x143;
auto x144 = x139 - V{1}*x226;
auto x145 = -x227;
auto x146 = x228;
auto x147 = x146*x50;
auto x148 = x145*x47 + x147;
auto x243 = x148 + x231;
auto x149 = x227;
auto x150 = x146*x55;
auto x232 = x150 + x231;
auto x151 = 0;
auto x152 = V{2}*x228;
auto x153 = x147 + x149*x47;
auto x248 = x153 + x231;
auto x154 = x12 + x7;
auto x155 = x126*x66;
auto x156 = x112*x68;
auto x234 = x155 + x156;
auto x157 = x215;
auto x245 = x157*x223;
auto x158 = x7 + x74;
auto x159 = x128 + x221*x69;
auto x160 = x16 + x7;
auto x161 = x110*x66;
auto x240 = x156 + x161;
auto x162 = x117 + x157*x221;
auto x163 = x7 + x87;
auto x164 = x12 + x16;
auto x165 = x126*x78;
auto x246 = x113 + x165;
auto x166 = x12 + x87;
cell[0] = x19*(x21*x39*(x26 + x33) - x6*(x18 + V{-1}) + V{-1});
cell[10] = x229 + x230*(x151*x66 + x152*x68 + x72) + x59*(x0*x7 + x148 + x149*x54 + x205*x21 + x232);
cell[11] = x229 + x230*(x151*x78 + x152*x79 + x82) + x59*(x0*x12 + x145*x54 + x153 + x207*x21 + x232);
cell[12] = x229 + x230*(x151*x90 + x152*x91 + x93) + x59*(x0*x16 + x151*x85 + x152*x86 + x208*x21 + x231);
cell[13] = x233 + x236*(x107 + x234 + x235) + x98*(x0*x154 + x151*x54 + x152*x55 + x21*((x154)*(x154)) + V{2}*x217*x69 + x231);
cell[14] = x233 + x236*(x114 + x234 + x237) + x98*(x0*x158 + x21*((x158)*(x158)) - x216*(x157 + x218) + x217*(x157 + x94) + x238 + x239);
cell[15] = x233 + x236*(x125 + x240 + x241) + x98*(x0*x160 + x127 + x159 + x21*((x160)*(x160)) + x242 + x243);
cell[16] = x233 + x236*(x133 + x240 + x244) + x98*(x0*x163 + x115 + x162 + x21*((x163)*(x163)) + x243 + x245);
cell[17] = x233 + x236*(x140 + x246 + x247) + x98*(x0*x164 + x141 + x159 + x21*((x164)*(x164)) - x242 + x248);
cell[18] = x233 + x236*(x144 + x246 + x249) + x98*(x0*x166 + x134 + x162 + x21*((x166)*(x166)) - x245 + x248);
cell[1] = x229 + x230*(x61*x66 + x67*x68 + x72) + x59*(x0*x57 + x21*((x57)*(x57)) + x213*x54 + x250 + x51);
cell[2] = x229 + x230*(x61*x78 + x67*x79 + x82) + x59*(x0*x74 + x21*((x74)*(x74)) + x250 + x44*x54 + x73);
cell[3] = x229 + x230*(x61*x90 + x67*x91 + x93) + x59*(x0*x87 + x21*((x87)*(x87)) + x231 + x61*x85 + x67*x86);
cell[4] = x233 + x236*(x107 + x237 + x251) + x98*(x0*x96 + x21*((x96)*(x96)) + V{2}*x214*x55 + x216*(x252 + x94) + x217*(x94 + x95) + x231);
cell[5] = x233 + x236*(x114 + x235 + x251) + x98*(x0*x109 + x21*((x109)*(x109)) + x216*(x252 + x69) + x217*(x69 + x95) - x238 + x239);
cell[6] = x233 + x236*(x125 + x244 + x253) + x98*(x0*x120 + x115 + x21*((x120)*(x120)) + x254 + x255);
cell[7] = x233 + x236*(x133 + x241 + x253) + x98*(x0*x130 + x127 + x21*((x130)*(x130)) + x255 + x256);
cell[8] = x233 + x236*(x140 + x249 + x257) + x98*(x0*x136 + x134 + x21*((x136)*(x136)) + x254 + x258);
cell[9] = x233 + x236*(x144 + x247 + x257) + x98*(x0*x142 + x141 + x21*((x142)*(x142)) + x256 + x258);
return { x6, x17 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
