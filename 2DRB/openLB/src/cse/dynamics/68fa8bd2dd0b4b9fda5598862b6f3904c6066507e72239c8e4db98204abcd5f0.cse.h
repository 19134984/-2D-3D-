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
struct CSE<Dual<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<1, 1>, momenta::RegularizedBoundaryStress<1, 1>, momenta::DefineSeparately> > >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x256 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<1, 1>::VELOCITY>(0);
auto x295 = cell.template getFieldComponent<olb::opti::F>(17);
auto x272 = cell.template getFieldComponent<olb::opti::DJDF>(13);
auto x264 = cell.template getFieldComponent<olb::opti::DJDF>(5);
auto x275 = cell.template getFieldComponent<olb::opti::DJDF>(16);
auto x293 = cell.template getFieldComponent<olb::opti::F>(15);
auto x296 = cell.template getFieldComponent<olb::opti::F>(18);
auto x279 = cell.template getFieldComponent<olb::opti::F>(1);
auto x281 = cell.template getFieldComponent<olb::opti::F>(3);
auto x284 = cell.template getFieldComponent<olb::opti::F>(6);
auto x263 = cell.template getFieldComponent<olb::opti::DJDF>(4);
auto x285 = cell.template getFieldComponent<olb::opti::F>(7);
auto x297 = parameters.template get<descriptors::OMEGA>();
auto x261 = cell.template getFieldComponent<olb::opti::DJDF>(2);
auto x262 = cell.template getFieldComponent<olb::opti::DJDF>(3);
auto x267 = cell.template getFieldComponent<olb::opti::DJDF>(8);
auto x268 = cell.template getFieldComponent<olb::opti::DJDF>(9);
auto x260 = cell.template getFieldComponent<olb::opti::DJDF>(1);
auto x283 = cell.template getFieldComponent<olb::opti::F>(5);
auto x255 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x258 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<1, 1>::VELOCITY>(2);
auto x289 = cell.template getFieldComponent<olb::opti::F>(11);
auto x278 = cell.template getFieldComponent<olb::opti::F>(0);
auto x291 = cell.template getFieldComponent<olb::opti::F>(13);
auto x276 = cell.template getFieldComponent<olb::opti::DJDF>(17);
auto x288 = cell.template getFieldComponent<olb::opti::F>(10);
auto x259 = cell.template getFieldComponent<olb::opti::DJDF>(0);
auto x271 = cell.template getFieldComponent<olb::opti::DJDF>(12);
auto x294 = cell.template getFieldComponent<olb::opti::F>(16);
auto x270 = cell.template getFieldComponent<olb::opti::DJDF>(11);
auto x290 = cell.template getFieldComponent<olb::opti::F>(12);
auto x265 = cell.template getFieldComponent<olb::opti::DJDF>(6);
auto x274 = cell.template getFieldComponent<olb::opti::DJDF>(15);
auto x277 = cell.template getFieldComponent<olb::opti::DJDF>(18);
auto x269 = cell.template getFieldComponent<olb::opti::DJDF>(10);
auto x266 = cell.template getFieldComponent<olb::opti::DJDF>(7);
auto x273 = cell.template getFieldComponent<olb::opti::DJDF>(14);
auto x257 = V{1} / (x255);
auto x280 = V{2}*x289;
auto x282 = x281 + x290;
auto x286 = x279 + x288;
auto x302 = ((x256)*(x256));
auto x303 = ((x258)*(x258));
auto x0 = V{1}*x257;
auto x1 = x257*(x278 + x280 + x282 + V{2}*x283 + x284 + x285 + x286 + V{2}*x291 + x293 + x294 + V{2}*x295 + V{2}*x296 + V{1}) + V{-1};
auto x2 = V{3};
auto x287 = ((x2)*(x2));
auto x325 = V{1}*x2;
auto x3 = V{0.5}*x2;
auto x333 = V{1}*x3;
auto x4 = V{0.0555555555555556};
auto x307 = -x4;
auto x5 = x255*x4;
auto x327 = V{1}*x5;
auto x328 = x3*x327;
auto x6 = V{1}/x2;
auto x7 = V{1} - x6;
auto x8 = V{0.5}*x287*x4;
auto x9 = -x6;
auto x10 = V{0.0277777777777778};
auto x304 = -x10;
auto x11 = V{0.5}*x10*x287;
auto x317 = V{2}*x11;
auto x12 = V{1} - x297;
auto x292 = x11*x12;
auto x298 = x12*x8;
auto x13 = x292*(cell[18] + cell[9]);
auto x14 = x292*(cell[17] + cell[8]);
auto x15 = x292*(cell[16] + cell[7]);
auto x16 = x15*x7;
auto x17 = x292*(cell[15] + cell[6]);
auto x18 = x17*x7;
auto x19 = x292*(cell[14] + cell[5]);
auto x20 = x19*x7;
auto x21 = x292*(cell[13] + cell[4]);
auto x22 = x21*x7;
auto x23 = x298*(cell[12] + cell[3]);
auto x24 = x23*x9;
auto x301 = x20 + x22 + x24;
auto x25 = x298*(cell[11] + cell[2]);
auto x26 = x25*x9;
auto x300 = x16 + x18 + x26;
auto x27 = x298*(cell[10] + cell[1]);
auto x28 = V{0.333333333333333};
auto x29 = V{0.5}*x28*x287;
auto x299 = x29*x9;
auto x30 = cell[0]*x12*x299;
auto x31 = x13*x9 + x14*x9 + x27*x7 + x30 + x300 + x301;
auto x32 = x10*x255;
auto x331 = V{1}*x32;
auto x332 = x3*x331;
auto x33 = -V{1}*x258;
auto x34 = x256 + x33;
auto x35 = V{0.5}*x287;
auto x312 = x1*x35;
auto x326 = V{2}*x35;
auto x36 = x3*(x302 + x303 + ((x1)*(x1)));
auto x305 = V{1} - x36;
auto x37 = x304 + x32*(x2*x34 + x305 + x35*((x34)*(x34)));
auto x38 = x1 + x258;
auto x306 = x35*x38;
auto x39 = x306;
auto x40 = -x10 + x32*(x2*x38 + x305 + x38*x39);
auto x41 = V{2}*x295 - V{2}*x40;
auto x42 = x294 - x37;
auto x43 = x256 + x258;
auto x44 = x304 + x32*(x2*x43 + x305 + x35*((x43)*(x43)));
auto x45 = x293 - x44;
auto x46 = x307 + x5*(x2*x258 + x303*x35 + x305);
auto x47 = -V{1}*x256;
auto x48 = x258 + x47;
auto x49 = x304 + x32*(x2*x48 + x305 + x35*((x48)*(x48)));
auto x50 = x285 - x49;
auto x51 = x307 + x5*(x2*x33 + x305 + x35*((x33)*(x33)));
auto x52 = x33 + x47;
auto x53 = x304 + x32*(x2*x52 + x305 + x35*((x52)*(x52)));
auto x54 = x284 - x53;
auto x309 = x42 + x45 + x50 + x54;
auto x55 = V{-2};
auto x56 = -V{1}*x55;
auto x57 = x1 + x33;
auto x308 = x35*x57;
auto x58 = x308;
auto x59 = -x10 + x32*(x2*x57 + x305 + x57*x58);
auto x60 = x296 - x59;
auto x61 = x282 + x309 + x41 - x46 - x51 + x56*x60;
auto x62 = x61*x7;
auto x63 = x1 + x256;
auto x310 = x35*x63;
auto x64 = x310;
auto x65 = -x10 + x32*(x2*x63 + x305 + x63*x64);
auto x66 = V{2}*x291 - V{2}*x65;
auto x67 = x1 + x47;
auto x311 = x35*x67;
auto x68 = x311;
auto x69 = -x10 + x32*(x2*x67 + x305 + x67*x68);
auto x70 = x283 - x69;
auto x71 = x312;
auto x72 = -x4 + x5*(x1*x2 + x1*x71 + x305);
auto x313 = V{2}*x72;
auto x73 = x280 - x313 + x41 + V{2}*x60 + x66 + V{2}*x70;
auto x74 = x73*x9;
auto x75 = x307 + x5*(x2*x256 + x302*x35 + x305);
auto x76 = x307 + x5*(x2*x47 + x305 + x35*((x47)*(x47)));
auto x77 = x286 + x309 + x56*x70 + x66 - x75 - x76;
auto x78 = x7*x77;
auto x314 = x74 + x78;
auto x315 = x314 + x62;
auto x79 = -V{1}*x42 + x45 - V{1}*x50 + x54;
auto x80 = x11*(x315 + x55*x79);
auto x81 = x11*(x315 + V{2}*x79);
auto x82 = x77*x9;
auto x83 = x74 + x82;
auto x84 = x8*(x62 + x83);
auto x85 = x61*x9;
auto x86 = x8*(x314 + x85);
auto x87 = x255*x28;
auto x88 = x7*x73;
auto x318 = x78 + x85 + x88;
auto x89 = x82 + x88;
auto x316 = x62 + x89;
auto x90 = x41 + x55*x60;
auto x91 = x55*x70 + x66;
auto x92 = x257*(-x28 + x29*(x83 + x85) + x313 + x317*(x316 + V{2}*x90) + x317*(x316 + x55*x90) + x317*(x318 + V{2}*x91) + x317*(x318 + x55*x91) + x37 + V{2}*x40 + x44 + x46 + x49 + x51 + x53 + V{2}*x59 + V{2}*x65 + V{2}*x69 + x75 + x76 + V{2}*x8*(x85 + x89) + V{2}*x80 + V{2}*x81 + V{2}*x84 + V{2}*x86 - x87*(x36 + V{-1}) + V{1}) + V{-1};
auto x319 = V{1}*x92;
auto x323 = x35*x92;
auto x93 = cell[18]*x32;
auto x94 = cell[17]*x32;
auto x95 = cell[14]*x32;
auto x96 = cell[13]*x32;
auto x97 = cell[11]*x5;
auto x98 = cell[9]*x32;
auto x99 = cell[8]*x32;
auto x100 = cell[5]*x32;
auto x101 = cell[4]*x32;
auto x102 = cell[2]*x5;
auto x103 = -x3*x319*(cell[0]*x87 + cell[10]*x5 + cell[12]*x5 + cell[15]*x32 + cell[16]*x32 + cell[1]*x5 + cell[3]*x5 + cell[6]*x32 + cell[7]*x32 + x100 + x101 + x102 + x93 + x94 + x95 + x96 + x97 + x98 + x99);
auto x104 = x33 + x92;
auto x320 = x104*x35;
auto x105 = x320;
auto x106 = x258 + x92;
auto x321 = x106*x35;
auto x107 = x321;
auto x108 = -x319;
auto x109 = x108 + x256;
auto x110 = x256 + x92;
auto x322 = x110*x35;
auto x111 = x322;
auto x112 = x323;
auto x113 = x108 + x258;
auto x114 = x108 + x33;
auto x115 = x47 + x92;
auto x324 = x115*x35;
auto x116 = x324;
auto x117 = x108 + x47;
auto x118 = x0*(x100*x116 + x100*x2 + x100*x324 - x101*x117*x326 - x101*x325 - x102*x108*x326 - x102*x325 + V{2}*x103 + x105*x93 + x107*x94 - x109*x326*x95 + x111*x96 + x112*x97 - x113*x326*x98 - x114*x326*x99 + x2*x93 + x2*x94 + x2*x96 + x2*x97 + x320*x93 + x321*x94 + x322*x96 + x323*x97 - x325*x95 - x325*x98 - x325*x99);
auto x119 = x31*x328;
auto x120 = x7*(x13 + x14);
auto x121 = x27*x9;
auto x329 = x120 + x121 + x30;
auto x122 = x19*x9 + x21*x9 + x23*x7 + x300 + x329;
auto x330 = x122 + x31;
auto x123 = x122*x328;
auto x124 = x15*x55 + V{2}*x17;
auto x125 = x124 + x330;
auto x126 = x125*x332;
auto x127 = -V{1}*x124 + x330;
auto x128 = x127*x332;
auto x129 = x13*x55 + V{2}*x14;
auto x130 = V{2}*x15*x9 + V{2}*x17*x9 + V{2}*x25*x7 + V{2}*x301 + V{2}*x329;
auto x131 = x122*x56 + x129*x55 + x130;
auto x132 = -x131*x331;
auto x133 = V{2}*x122 + V{2}*x129 + x130;
auto x134 = -x133*x331;
auto x135 = x19*x55 + V{2}*x21;
auto x136 = x130 + V{2}*x135 + V{2}*x31;
auto x137 = -x136*x331;
auto x138 = -x130*x327;
auto x139 = x130 + x135*x55 + x31*x56;
auto x140 = -x139*x331;
auto x141 = -x92*(-V{2}*x119 - V{2}*x123 - V{2}*x126 - V{2}*x128 + x132*x333 + x134*x333 + x137*x333 + x138*x333 + x140*x333);
auto x142 = x0*(x105*x132 + x107*x134 + x111*x137 + x112*x138 + x116*x140 + x132*x2 + x132*x320 + x134*x2 + x134*x321 + x137*x2 + x137*x322 + x138*x2 + x138*x323 + x140*x2 + x140*x324 + V{2}*x141);
auto x334 = x118 + x142;
auto x143 = x31 + x334;
auto x144 = x143*x8;
auto x145 = x144*x7;
auto x146 = x122 + x334;
auto x147 = x146*x8;
auto x148 = x147*x9;
auto x149 = x125 + x334;
auto x150 = x11*x149;
auto x151 = x150*x7;
auto x152 = x127 + x334;
auto x153 = x11*x152;
auto x154 = x153*x7;
auto x155 = V{2}*x118;
auto x156 = V{2}*x142;
auto x335 = x155 + x156;
auto x157 = x131 + x335;
auto x158 = x11*x157;
auto x159 = x133 + x335;
auto x160 = x11*x159;
auto x161 = x136 + x335;
auto x162 = x11*x161;
auto x163 = x162*x7;
auto x164 = x130 + x335;
auto x165 = x164*x8;
auto x166 = x165*x9;
auto x167 = x139 + x335;
auto x168 = x11*x167;
auto x169 = x168*x7;
auto x170 = x334;
auto x171 = x170*x299;
auto x336 = V{2}*x151 + V{2}*x154 + x166 + x171;
auto x172 = V{2}*x145 + V{2}*x148 + x158*x9 + x160*x9 + x163 + x169 + x336;
auto x173 = x172*x328;
auto x174 = x144*x9;
auto x175 = x147*x7;
auto x176 = x7*(x158 + x160);
auto x177 = x162*x9 + x168*x9 + V{2}*x174 + V{2}*x175 + x176 + x336;
auto x337 = x172 + x177;
auto x178 = x177*x328;
auto x179 = V{2}*x150;
auto x180 = x153*x55;
auto x181 = V{2}*x179 + V{2}*x180;
auto x182 = x181 + x337;
auto x183 = x182*x332;
auto x184 = -V{1}*x181 + x337;
auto x185 = x184*x332;
auto x186 = x158*x55 + V{2}*x160;
auto x187 = x150*x9;
auto x188 = x153*x9;
auto x189 = V{4}*x148 + V{2}*x163 + V{2}*x165*x7 + V{2}*x169 + V{2}*x171 + V{4}*x174 + V{2}*x176 + V{4}*x187 + V{4}*x188;
auto x190 = x177*x56 + x186*x55 + x189;
auto x191 = -x190*x331;
auto x192 = V{2}*x177 + V{2}*x186 + x189;
auto x193 = -x192*x331;
auto x194 = V{2}*x162 + x168*x55;
auto x195 = V{2}*x172 + x189 + V{2}*x194;
auto x196 = -x195*x331;
auto x197 = -x189*x327;
auto x198 = x172*x56 + x189 + x194*x55;
auto x199 = -x198*x331;
auto x200 = -x1*(-V{2}*x173 - V{2}*x178 - V{2}*x183 - V{2}*x185 + x191*x333 + x193*x333 + x196*x333 + x197*x333 + x199*x333);
auto x201 = -x143*x328;
auto x202 = -x146*x328;
auto x203 = -x149*x332;
auto x204 = -x152*x332;
auto x205 = x157*x32;
auto x206 = x159*x32;
auto x207 = x161*x32;
auto x208 = x164*x5;
auto x209 = x167*x32;
auto x210 = -x1*(x170*x333*x87 - V{2}*x201 - V{2}*x202 - V{2}*x203 - V{2}*x204 + x205*x333 + x206*x333 + x207*x333 + x208*x333 + x209*x333);
auto x211 = x0*(x191*x2 + x191*x308 + x191*x58 + x193*x2 + x193*x306 + x193*x39 + x196*x2 + x196*x310 + x196*x64 + x197*x2 + x197*x312 + x197*x71 + x199*x2 + x199*x311 + x199*x68 + x2*x205 + x2*x206 + x2*x207 + x2*x208 + x2*x209 + V{2}*x200 + x205*x308 + x205*x58 + x206*x306 + x206*x39 + x207*x310 + x207*x64 + x208*x312 + x208*x71 + x209*x311 + x209*x68 + V{2}*x210);
auto x212 = x211 - x259;
auto x213 = x172 + x211;
auto x214 = x213 - x260;
auto x215 = -x261;
auto x216 = x177 + x211;
auto x217 = x216 - x262;
auto x218 = -x263;
auto x219 = V{2}*x211;
auto x220 = x198 + x219 - x264;
auto x221 = x182 + x211;
auto x222 = x221 - x265;
auto x223 = x184 + x211;
auto x224 = x223 - x266;
auto x225 = -x267;
auto x226 = -x268;
auto x227 = x213 - x269;
auto x228 = x189 + x219 - x270;
auto x229 = x216 - x271;
auto x230 = x195 + x219 - x272;
auto x231 = -x273;
auto x232 = x221 - x274;
auto x233 = x223 - x275;
auto x234 = x192 + x219 - x276;
auto x235 = x190 + x219 - x277;
cell[0] = x212;
cell[1] = x214;
cell[2] = x215;
cell[3] = x217;
cell[4] = x218;
cell[5] = x220;
cell[6] = x222;
cell[7] = x224;
cell[8] = x225;
cell[9] = x226;
cell[10] = x227;
cell[11] = x228;
cell[12] = x229;
cell[13] = x230;
cell[14] = x231;
cell[15] = x232;
cell[16] = x233;
cell[17] = x234;
cell[18] = x235;
return { x212 + V{1}, ((x212)*(x212)) + ((x214)*(x214)) + ((x215)*(x215)) + ((x217)*(x217)) + ((x218)*(x218)) + ((x220)*(x220)) + ((x222)*(x222)) + ((x224)*(x224)) + ((x225)*(x225)) + ((x226)*(x226)) + ((x227)*(x227)) + ((x228)*(x228)) + ((x229)*(x229)) + ((x230)*(x230)) + ((x231)*(x231)) + ((x232)*(x232)) + ((x233)*(x233)) + ((x234)*(x234)) + ((x235)*(x235)) };
}
};

}

}

// Generation Info: commit=4ea71c45fdf920fe461e9926f71dc63f7fa1062d
