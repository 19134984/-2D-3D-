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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::ParameterFromCell<collision::LES::SMAGORINSKY, collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB> > >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x185 = cell.template getFieldComponent<olb::collision::LES::SMAGORINSKY>(0);
auto x187 = parameters.template get<descriptors::OMEGA>();
auto x186 = cell[17] + cell[9];
auto x188 = cell[15] + cell[7];
auto x189 = cell[12] + x186 + x188;
auto x190 = cell[3] + cell[6] + x189;
auto x191 = cell[13] + cell[5];
auto x192 = cell[11] + cell[18] + x191;
auto x193 = cell[2] + cell[4] + x192;
auto x194 = cell[10] + cell[14] + cell[16];
auto x195 = cell[1] + x194;
auto x0 = V{3};
auto x203 = ((x0)*(x0));
auto x204 = V{1}/x0;
auto x1 = -V{1}*cell[4];
auto x2 = -V{1}*cell[5];
auto x198 = cell[13] + x2;
auto x3 = -V{1}*cell[6];
auto x4 = -V{1}*cell[7];
auto x197 = cell[15] + x4;
auto x5 = cell[0] + cell[8] + x190 + x193 + x195 + V{1};
auto x196 = V{1} / (x5);
auto x6 = x196*(-V{1}*cell[1] + x1 + x194 + x197 + x198 + x3);
auto x199 = ((x6)*(x6));
auto x7 = x199;
auto x8 = -V{1}*cell[8];
auto x9 = -V{1}*cell[9];
auto x200 = cell[17] + x9;
auto x10 = -V{1}*cell[14];
auto x11 = x196*(-V{1}*cell[2] + x1 + x10 + x192 + x200 + x8);
auto x201 = ((x11)*(x11));
auto x12 = x201;
auto x13 = -V{1}*cell[16];
auto x14 = -V{1}*cell[18];
auto x15 = x196*(-V{1}*cell[3] + x13 + x14 + x189 + x3 + x8);
auto x202 = ((x15)*(x15));
auto x16 = x12 + x202 + x7;
auto x17 = V{0.5}*x0*x16;
auto x225 = V{1} - x17;
auto x18 = V{0.333333333333333};
auto x19 = x203;
auto x20 = V{0.5}*x19;
auto x21 = x204;
auto x209 = V{1} - x21;
auto x22 = -x21;
auto x23 = x21*(x5 + V{-1});
auto x205 = -x23;
auto x206 = cell[8] + x205;
auto x24 = cell[16] + cell[18] + x190 - x202*x5 + x206;
auto x25 = x22*x24;
auto x26 = x11*x5;
auto x27 = cell[14] - x11*x26 + x186 + x193 + x206;
auto x28 = x22*x27;
auto x29 = x5*x6;
auto x30 = cell[4] + cell[6] + x188 + x191 + x195 + x205 - x29*x6;
auto x31 = x22*x30;
auto x32 = x28 + x31;
auto x33 = V{1}/x187;
auto x34 = cell[6] + x13 - x15*x29 + x197;
auto x219 = x19*x34;
auto x35 = ((x34)*(x34));
auto x36 = cell[4] + x10 - x11*x29 + x198;
auto x214 = x19*x36;
auto x37 = cell[8] + x14 - x15*x26 + x200;
auto x220 = x19*x37;
auto x38 = V{1} - V{2}*V{1} / (x33 + util::sqrt(V{2.828428}*x196*x203*((x185)*(x185))*util::sqrt(V{2}*x35 + ((x24)*(x24)) + ((x27)*(x27)) + ((x30)*(x30)) + V{2}*((x36)*(x36)) + V{2}*((x37)*(x37))) + ((x33)*(x33))));
auto x39 = x0*x19;
auto x40 = V{0.166667}*x39;
auto x221 = x22*x40;
auto x41 = -x204;
auto x42 = -x41;
auto x207 = x40*x42;
auto x43 = -x207;
auto x44 = x11*x7;
auto x45 = x201*x6;
auto x46 = x44 - x45;
auto x210 = x40*x46;
auto x47 = V{0.5}*x39;
auto x208 = x42*x47;
auto x222 = x22*x47;
auto x48 = x208;
auto x49 = x44 + x45;
auto x211 = x47*x49;
auto x50 = x43*x46 + x48*x49;
auto x51 = x202*x6;
auto x52 = x11*x202;
auto x53 = x51 - x52;
auto x217 = x40*x53;
auto x232 = x217*(x22 - x42);
auto x54 = x51 + x52;
auto x215 = x47*x54;
auto x55 = x48*x54;
auto x244 = x225 + x55;
auto x56 = -V{1}*x6;
auto x57 = V{0.0555555555555556};
auto x223 = -x57;
auto x224 = x38*x57;
auto x58 = x5*x57;
auto x59 = -V{1}*x22;
auto x60 = 0;
auto x61 = V{2}*x15;
auto x62 = x24*x6 + x34*x61;
auto x63 = V{2}*x11;
auto x64 = x27*x6 + x36*x63;
auto x65 = x62 - x64;
auto x66 = V{2}*x47*x59;
auto x67 = x62 + x64;
auto x68 = x209;
auto x216 = x40*(x22 - x68);
auto x236 = x217*x68;
auto x69 = x30*x68;
auto x70 = x28 + x69;
auto x71 = x20*(x25 + x70);
auto x72 = x59*(x210 + x211);
auto x73 = -V{1}*x11;
auto x74 = V{2}*x6;
auto x75 = x11*x30 + x36*x74;
auto x76 = x11*x24 + x37*x61;
auto x77 = x75 - x76;
auto x78 = x75 + x76;
auto x79 = x27*x68;
auto x80 = x31 + x79;
auto x81 = x20*(x25 + x80);
auto x82 = x12*x15;
auto x83 = x15*x7;
auto x84 = x82 - x83;
auto x85 = x82 + x83;
auto x86 = -V{1}*x15;
auto x87 = x15*x27 + x37*x63;
auto x88 = x15*x30 + x34*x74;
auto x89 = x87 - x88;
auto x90 = x87 + x88;
auto x91 = x24*x68;
auto x92 = x20*(x32 + x91);
auto x93 = -V{1}*x68;
auto x212 = -x93;
auto x213 = x40*(x212 + x59);
auto x94 = -x41 + V{-1};
auto x218 = x217*x94;
auto x246 = -x94;
auto x95 = x56 + x73;
auto x96 = V{0.0277777777777778};
auto x227 = -x96;
auto x230 = x38*x96;
auto x97 = x5*x96;
auto x98 = x213;
auto x99 = x65*x98;
auto x100 = -x213;
auto x101 = x100*x77;
auto x102 = x47*(x59 + x93);
auto x103 = x102*x67;
auto x245 = x103 + x99;
auto x104 = x102*x78;
auto x231 = x101 + x104;
auto x105 = x20*(x25 + x69 + x79);
auto x106 = x105 + x214;
auto x107 = x215*(x22 + x42);
auto x233 = x107 + x225;
auto x108 = x11 + x56;
auto x109 = -x216;
auto x110 = x109*x77;
auto x111 = x47*(x22 + x68);
auto x112 = x111*x78;
auto x229 = x110 + x112;
auto x113 = x105 - V{1}*x214;
auto x114 = x84*x98;
auto x115 = x218;
auto x249 = x115 + x50;
auto x116 = x102*x85;
auto x117 = x215*x94;
auto x118 = x116 + x117;
auto x248 = x118 + x225;
auto x119 = x56 + x86;
auto x120 = x89*x98;
auto x121 = x100*x65;
auto x247 = x103 + x121;
auto x122 = x102*x90;
auto x238 = x120 + x122;
auto x123 = x20*(x70 + x91);
auto x124 = x123 + x219;
auto x125 = x216;
auto x126 = x125*x84;
auto x127 = x111*x85;
auto x128 = x117 + x127;
auto x250 = x128 + x225;
auto x129 = x15 + x56;
auto x130 = x125*x89;
auto x131 = x111*x90;
auto x235 = x130 + x131;
auto x132 = x123 - V{1}*x219;
auto x133 = x100*x84;
auto x134 = -x218;
auto x252 = x134 + x72;
auto x135 = x73 + x86;
auto x136 = x100*x89;
auto x243 = x122 + x136;
auto x137 = x77*x98;
auto x251 = x104 + x137;
auto x138 = x20*(x80 + x91);
auto x139 = x138 + x220;
auto x140 = x109*x84;
auto x141 = x15 + x73;
auto x142 = x109*x89;
auto x241 = x131 + x142;
auto x143 = x138 - V{1}*x220;
auto x144 = -x221;
auto x145 = x222;
auto x146 = x145*x49;
auto x147 = x144*x46 + x146;
auto x237 = x147 + x225;
auto x148 = x221;
auto x149 = x145*x54;
auto x226 = x149 + x225;
auto x150 = 0;
auto x151 = V{2}*x222;
auto x152 = x146 + x148*x46;
auto x242 = x152 + x225;
auto x153 = x11 + x6;
auto x154 = x125*x65;
auto x155 = x111*x67;
auto x228 = x154 + x155;
auto x156 = x209;
auto x239 = x156*x217;
auto x157 = x6 + x73;
auto x158 = x127 + x215*x68;
auto x159 = x15 + x6;
auto x160 = x109*x65;
auto x234 = x155 + x160;
auto x161 = x116 + x156*x215;
auto x162 = x6 + x86;
auto x163 = x11 + x15;
auto x164 = x125*x77;
auto x240 = x112 + x164;
auto x165 = x11 + x86;
cell[0] = x18*(x20*x38*(x25 + x32) - x5*(x17 + V{-1}) + V{-1});
cell[10] = x223 + x224*(x150*x65 + x151*x67 + x71) + x58*(x0*x6 + x147 + x148*x53 + x199*x20 + x226);
cell[11] = x223 + x224*(x150*x77 + x151*x78 + x81) + x58*(x0*x11 + x144*x53 + x152 + x20*x201 + x226);
cell[12] = x223 + x224*(x150*x89 + x151*x90 + x92) + x58*(x0*x15 + x150*x84 + x151*x85 + x20*x202 + x225);
cell[13] = x227 + x230*(x106 + x228 + x229) + x97*(x0*x153 + x150*x53 + x151*x54 + x20*((x153)*(x153)) + V{2}*x211*x68 + x225);
cell[14] = x227 + x230*(x113 + x228 + x231) + x97*(x0*x157 + x20*((x157)*(x157)) - x210*(x156 + x212) + x211*(x156 + x93) + x232 + x233);
cell[15] = x227 + x230*(x124 + x234 + x235) + x97*(x0*x159 + x126 + x158 + x20*((x159)*(x159)) + x236 + x237);
cell[16] = x227 + x230*(x132 + x234 + x238) + x97*(x0*x162 + x114 + x161 + x20*((x162)*(x162)) + x237 + x239);
cell[17] = x227 + x230*(x139 + x240 + x241) + x97*(x0*x163 + x140 + x158 + x20*((x163)*(x163)) - x236 + x242);
cell[18] = x227 + x230*(x143 + x240 + x243) + x97*(x0*x165 + x133 + x161 + x20*((x165)*(x165)) - x239 + x242);
cell[1] = x223 + x224*(x60*x65 + x66*x67 + x71) + x58*(x0*x56 + x20*((x56)*(x56)) + x207*x53 + x244 + x50);
cell[2] = x223 + x224*(x60*x77 + x66*x78 + x81) + x58*(x0*x73 + x20*((x73)*(x73)) + x244 + x43*x53 + x72);
cell[3] = x223 + x224*(x60*x89 + x66*x90 + x92) + x58*(x0*x86 + x20*((x86)*(x86)) + x225 + x60*x84 + x66*x85);
cell[4] = x227 + x230*(x106 + x231 + x245) + x97*(x0*x95 + x20*((x95)*(x95)) + V{2}*x208*x54 + x210*(x246 + x93) + x211*(x93 + x94) + x225);
cell[5] = x227 + x230*(x113 + x229 + x245) + x97*(x0*x108 + x20*((x108)*(x108)) + x210*(x246 + x68) + x211*(x68 + x94) - x232 + x233);
cell[6] = x227 + x230*(x124 + x238 + x247) + x97*(x0*x119 + x114 + x20*((x119)*(x119)) + x248 + x249);
cell[7] = x227 + x230*(x132 + x235 + x247) + x97*(x0*x129 + x126 + x20*((x129)*(x129)) + x249 + x250);
cell[8] = x227 + x230*(x139 + x243 + x251) + x97*(x0*x135 + x133 + x20*((x135)*(x135)) + x248 + x252);
cell[9] = x227 + x230*(x143 + x241 + x251) + x97*(x0*x141 + x140 + x20*((x141)*(x141)) + x250 + x252);
return { x5, x16 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
