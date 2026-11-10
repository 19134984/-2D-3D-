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
struct CSE<Dual<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x156 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x157 = cell.template getFieldComponent<olb::opti::DJDF>(0);
auto x158 = cell.template getFieldComponent<olb::opti::DJDF>(1);
auto x159 = cell.template getFieldComponent<olb::opti::DJDF>(10);
auto x160 = cell.template getFieldComponent<olb::opti::DJDF>(11);
auto x161 = cell.template getFieldComponent<olb::opti::DJDF>(12);
auto x162 = cell.template getFieldComponent<olb::opti::DJDF>(13);
auto x163 = cell.template getFieldComponent<olb::opti::DJDF>(14);
auto x164 = cell.template getFieldComponent<olb::opti::DJDF>(15);
auto x165 = cell.template getFieldComponent<olb::opti::DJDF>(16);
auto x166 = cell.template getFieldComponent<olb::opti::DJDF>(17);
auto x167 = cell.template getFieldComponent<olb::opti::DJDF>(18);
auto x168 = cell.template getFieldComponent<olb::opti::DJDF>(2);
auto x169 = cell.template getFieldComponent<olb::opti::DJDF>(3);
auto x170 = cell.template getFieldComponent<olb::opti::DJDF>(4);
auto x171 = cell.template getFieldComponent<olb::opti::DJDF>(5);
auto x172 = cell.template getFieldComponent<olb::opti::DJDF>(6);
auto x173 = cell.template getFieldComponent<olb::opti::DJDF>(7);
auto x174 = cell.template getFieldComponent<olb::opti::DJDF>(8);
auto x175 = cell.template getFieldComponent<olb::opti::DJDF>(9);
auto x176 = cell.template getFieldComponent<olb::opti::F>(0);
auto x177 = cell.template getFieldComponent<olb::opti::F>(1);
auto x178 = cell.template getFieldComponent<olb::opti::F>(10);
auto x179 = cell.template getFieldComponent<olb::opti::F>(11);
auto x180 = cell.template getFieldComponent<olb::opti::F>(12);
auto x181 = cell.template getFieldComponent<olb::opti::F>(13);
auto x182 = cell.template getFieldComponent<olb::opti::F>(14);
auto x183 = cell.template getFieldComponent<olb::opti::F>(15);
auto x184 = cell.template getFieldComponent<olb::opti::F>(16);
auto x185 = cell.template getFieldComponent<olb::opti::F>(17);
auto x186 = cell.template getFieldComponent<olb::opti::F>(18);
auto x187 = cell.template getFieldComponent<olb::opti::F>(2);
auto x188 = cell.template getFieldComponent<olb::opti::F>(3);
auto x189 = cell.template getFieldComponent<olb::opti::F>(4);
auto x190 = cell.template getFieldComponent<olb::opti::F>(5);
auto x191 = cell.template getFieldComponent<olb::opti::F>(6);
auto x192 = cell.template getFieldComponent<olb::opti::F>(7);
auto x193 = cell.template getFieldComponent<olb::opti::F>(8);
auto x194 = cell.template getFieldComponent<olb::opti::F>(9);
auto x195 = parameters.template get<descriptors::OMEGA>();
auto x198 = x179 + x186 + x190;
auto x199 = x180 + x192 + x194;
auto x200 = x178 + x182 + x184;
auto x0 = V{0.0277777777777778};
auto x1 = -V{1}*x189;
auto x197 = x1 + x181;
auto x2 = -V{1}*x193;
auto x196 = x185 + x2;
auto x3 = -V{1}*x182 - V{1}*x187 - V{1}*x194 + x196 + x197 + x198;
auto x4 = x176 + x177 + x181 + x183 + x185 + x187 + x188 + x189 + x191 + x193 + x198 + x199 + x200 + V{1};
auto x201 = V{1} / (x4);
auto x202 = x156*x201;
auto x5 = x202*x3;
auto x6 = -V{1}*x191;
auto x203 = x183 + x6;
auto x7 = -V{1}*x184 - V{1}*x186 - V{1}*x188 + x196 + x199 + x203;
auto x8 = x202*x7;
auto x9 = -V{1}*x8;
auto x10 = x5 + x9;
auto x11 = V{3};
auto x12 = V{0.5}*((x11)*(x11));
auto x204 = x10*x12;
auto x210 = x12*x8;
auto x211 = x12*x5;
auto x219 = x12*x9;
auto x13 = x204;
auto x14 = V{0.5}*x11;
auto x15 = -V{1}*x177 - V{1}*x190 - V{1}*x192 + x197 + x200 + x203;
auto x16 = x15*x202;
auto x212 = x12*x16;
auto x17 = x14*(((x16)*(x16)) + ((x5)*(x5)) + ((x8)*(x8)));
auto x223 = V{1} - x17;
auto x18 = cell[18]*x195;
auto x19 = x5 + x8;
auto x205 = x12*x19;
auto x20 = x205;
auto x21 = cell[17]*x195;
auto x22 = x16 + x9;
auto x206 = x12*x22;
auto x23 = x206;
auto x24 = cell[16]*x195;
auto x25 = x16 + x8;
auto x207 = x12*x25;
auto x26 = x207;
auto x27 = cell[15]*x195;
auto x28 = -V{1}*x5;
auto x220 = x12*x28;
auto x29 = x16 + x28;
auto x208 = x12*x29;
auto x30 = x208;
auto x31 = cell[14]*x195;
auto x32 = x16 + x5;
auto x209 = x12*x32;
auto x33 = x209;
auto x34 = cell[13]*x195;
auto x35 = V{0.0555555555555556};
auto x36 = x210;
auto x37 = cell[12]*x195;
auto x38 = x211;
auto x39 = cell[11]*x195;
auto x40 = x212;
auto x41 = cell[10]*x195;
auto x42 = x28 + x8;
auto x213 = x12*x42;
auto x43 = x213;
auto x44 = cell[9]*x195;
auto x45 = x28 + x9;
auto x214 = x12*x45;
auto x46 = x214;
auto x47 = cell[8]*x195;
auto x48 = -V{1}*x16;
auto x221 = x12*x48;
auto x49 = x48 + x8;
auto x215 = x12*x49;
auto x50 = x215;
auto x51 = cell[7]*x195;
auto x52 = x48 + x9;
auto x216 = x12*x52;
auto x53 = x216;
auto x54 = cell[6]*x195;
auto x55 = x48 + x5;
auto x217 = x12*x55;
auto x56 = x217;
auto x57 = cell[5]*x195;
auto x58 = x28 + x48;
auto x218 = x12*x58;
auto x59 = x218;
auto x60 = cell[4]*x195;
auto x61 = x219;
auto x62 = cell[3]*x195;
auto x63 = x220;
auto x64 = cell[2]*x195;
auto x65 = x221;
auto x66 = cell[1]*x195;
auto x67 = V{0.333333333333333};
auto x68 = cell[0]*x195;
auto x222 = x67*x68;
auto x69 = x0*x18*(x10*x11 + x10*x13 + x223) + x0*x21*(x11*x19 + x19*x20 + x223) + x0*x24*(x11*x22 + x22*x23 + x223) + x0*x27*(x11*x25 + x223 + x25*x26) + x0*x31*(x11*x29 + x223 + x29*x30) + x0*x34*(x11*x32 + x223 + x32*x33) + x0*x44*(x11*x42 + x223 + x42*x43) + x0*x47*(x11*x45 + x223 + x45*x46) + x0*x51*(x11*x49 + x223 + x49*x50) + x0*x54*(x11*x52 + x223 + x52*x53) + x0*x57*(x11*x55 + x223 + x55*x56) + x0*x60*(x11*x58 + x223 + x58*x59) - x222*(x17 + V{-1}) + x35*x37*(x11*x8 + x223 + x36*x8) + x35*x39*(x11*x5 + x223 + x38*x5) + x35*x41*(x11*x16 + x16*x40 + x223) + x35*x62*(x11*x9 + x223 + x61*x9) + x35*x64*(x11*x28 + x223 + x28*x63) + x35*x66*(x11*x48 + x223 + x48*x65);
auto x70 = V{1} - x195;
auto x71 = x0*x4;
auto x72 = x18*x71;
auto x73 = x21*x71;
auto x74 = x24*x71;
auto x75 = x27*x71;
auto x76 = x31*x71;
auto x77 = x34*x71;
auto x78 = x35*x4;
auto x79 = x37*x78;
auto x80 = x39*x78;
auto x81 = x41*x78;
auto x82 = x44*x71;
auto x83 = x47*x71;
auto x84 = x51*x71;
auto x85 = x54*x71;
auto x86 = x57*x71;
auto x87 = x60*x71;
auto x88 = x62*x78;
auto x224 = V{1}*x88;
auto x89 = x64*x78;
auto x227 = V{1}*x89;
auto x90 = x66*x78;
auto x229 = V{1}*x90;
auto x91 = -V{1}*x14*(x222*x4 + x72 + x73 + x74 + x75 + x76 + x77 + x79 + x80 + x81 + x82 + x83 + x84 + x85 + x86 + x87 + x88 + x89 + x90);
auto x92 = x8*x91;
auto x93 = x72*(x11 + x13 + x204);
auto x94 = x73*(x11 + x20 + x205);
auto x95 = x74*(x11 + x206 + x23);
auto x96 = x75*(x11 + x207 + x26);
auto x97 = x82*(x11 + x213 + x43);
auto x98 = -V{1}*x83*(x11 + x214 + x46);
auto x225 = x94 + x98;
auto x99 = x84*(x11 + x215 + x50);
auto x100 = -V{1}*x85*(x11 + x216 + x53);
auto x226 = x100 + x96;
auto x101 = x156*(-x11*x224 + x11*x79 + x210*x79 - x219*x224 - x224*x61 + x225 + x226 + x36*x79 + V{2}*x92 - V{1}*x93 - V{1}*x95 + x97 + x99);
auto x102 = x5*x91;
auto x103 = x76*(x11 + x208 + x30);
auto x104 = x77*(x11 + x209 + x33);
auto x105 = x86*(x11 + x217 + x56);
auto x106 = -V{1}*x87*(x11 + x218 + x59);
auto x228 = x104 + x106;
auto x107 = x156*(V{2}*x102 - V{1}*x103 + x105 - x11*x227 + x11*x80 + x211*x80 - x220*x227 + x225 - x227*x63 + x228 + x38*x80 + x93 - V{1}*x97);
auto x108 = x16*x91;
auto x109 = x156*(x103 - V{1}*x105 + V{2}*x108 - x11*x229 + x11*x81 + x212*x81 - x221*x229 + x226 + x228 - x229*x65 + x40*x81 + x95 - V{1}*x99);
auto x110 = -V{1}*(x101*x7 + x107*x3 + x109*x15)/((x4)*(x4));
auto x230 = x110 + x69;
auto x111 = cell[0]*x70 - x157 + x230;
auto x112 = V{1}*x201;
auto x113 = x109*x112;
auto x234 = x113 + x230;
auto x114 = -V{1}*x113;
auto x231 = x114 + x230;
auto x115 = cell[1]*x70 - x158 + x231;
auto x116 = x107*x112;
auto x235 = x116 + x230;
auto x117 = -V{1}*x116;
auto x232 = x117 + x230;
auto x118 = cell[2]*x70 - x168 + x232;
auto x119 = x101*x112;
auto x120 = -V{1}*x119;
auto x233 = x120 + x230;
auto x121 = cell[3]*x70 - x169 + x233;
auto x122 = cell[4]*x70 + x117 - x170 + x231;
auto x123 = cell[5]*x70 + x116 - x171 + x231;
auto x124 = cell[6]*x70 + x120 - x172 + x231;
auto x125 = cell[7]*x70 + x119 - x173 + x231;
auto x126 = cell[8]*x70 + x120 - x174 + x232;
auto x127 = cell[9]*x70 + x119 - x175 + x232;
auto x128 = cell[10]*x70 - x159 + x234;
auto x129 = cell[11]*x70 - x160 + x235;
auto x130 = cell[12]*x70 + x119 - x161 + x230;
auto x131 = cell[13]*x70 + x116 - x162 + x234;
auto x132 = cell[14]*x70 + x113 - x163 + x232;
auto x133 = cell[15]*x70 + x119 - x164 + x234;
auto x134 = cell[16]*x70 + x113 - x165 + x233;
auto x135 = cell[17]*x70 + x119 - x166 + x235;
auto x136 = cell[18]*x70 + x116 - x167 + x233;
cell[0] = x111;
cell[10] = x128;
cell[11] = x129;
cell[12] = x130;
cell[13] = x131;
cell[14] = x132;
cell[15] = x133;
cell[16] = x134;
cell[17] = x135;
cell[18] = x136;
cell[1] = x115;
cell[2] = x118;
cell[3] = x121;
cell[4] = x122;
cell[5] = x123;
cell[6] = x124;
cell[7] = x125;
cell[8] = x126;
cell[9] = x127;
return { x111 + V{1}, ((x111)*(x111)) + ((x115)*(x115)) + ((x118)*(x118)) + ((x121)*(x121)) + ((x122)*(x122)) + ((x123)*(x123)) + ((x124)*(x124)) + ((x125)*(x125)) + ((x126)*(x126)) + ((x127)*(x127)) + ((x128)*(x128)) + ((x129)*(x129)) + ((x130)*(x130)) + ((x131)*(x131)) + ((x132)*(x132)) + ((x133)*(x133)) + ((x134)*(x134)) + ((x135)*(x135)) + ((x136)*(x136)) };
}
};

}

}

// Generation Info: commit=a053d6683a8855587b7af312fff11115a5b3dd07
