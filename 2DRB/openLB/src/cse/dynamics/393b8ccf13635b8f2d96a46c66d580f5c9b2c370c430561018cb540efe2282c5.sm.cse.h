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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::ParameterFromCell<collision::HYBRID, collision::LocalSmagorinskyEffectiveOmega<collision::HRR> > >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x225 = cell.template getFieldComponent<olb::collision::HYBRID>(0);
auto x227 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x228 = cell.template getFieldComponent<olb::descriptors::TENSOR>(0);
auto x229 = cell.template getFieldComponent<olb::descriptors::TENSOR>(1);
auto x230 = cell.template getFieldComponent<olb::descriptors::TENSOR>(2);
auto x231 = cell.template getFieldComponent<olb::descriptors::TENSOR>(3);
auto x232 = cell.template getFieldComponent<olb::descriptors::TENSOR>(4);
auto x233 = cell.template getFieldComponent<olb::descriptors::TENSOR>(5);
auto x234 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x235 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x236 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x237 = cell.template getFieldComponent<olb::descriptors::WMPOROSITY>(0);
auto x238 = cell.template getFieldComponent<olb::descriptors::WMVELOCITY>(0);
auto x239 = cell.template getFieldComponent<olb::descriptors::WMVELOCITY>(1);
auto x240 = cell.template getFieldComponent<olb::descriptors::WMVELOCITY>(2);
auto x242 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x243 = parameters.template get<descriptors::OMEGA>();
auto x226 = cell[10] + cell[14] + cell[16];
auto x241 = cell[1] + cell[4] + x226;
auto x244 = cell[15] + cell[7];
auto x245 = cell[17] + cell[9];
auto x246 = cell[12] + x244 + x245;
auto x247 = cell[3] + cell[8] + x246;
auto x248 = cell[13] + cell[5];
auto x249 = cell[11] + cell[18] + x248;
auto x250 = cell[2] + x249;
auto x0 = V{3};
auto x259 = ((x0)*(x0));
auto x260 = V{1} / (x0);
auto x261 = V{1}*x260;
auto x1 = V{1} - x237;
auto x2 = x1*x238;
auto x3 = V{1} - x227;
auto x4 = x234*x3;
auto x5 = -V{1}*cell[4];
auto x6 = -V{1}*cell[5];
auto x252 = cell[13] + x6;
auto x7 = -V{1}*cell[6];
auto x8 = -V{1}*cell[7];
auto x251 = cell[15] + x8;
auto x9 = cell[0] + cell[6] + x241 + x247 + x250 + V{1};
auto x253 = V{1} / (x9);
auto x254 = x227*x253;
auto x10 = x2 + x237*(x254*(-V{1}*cell[1] + x226 + x251 + x252 + x5 + x7) + x4);
auto x255 = ((x10)*(x10));
auto x11 = x255;
auto x12 = x1*x239;
auto x13 = x235*x3;
auto x14 = -V{1}*cell[8];
auto x15 = -V{1}*cell[9];
auto x256 = cell[17] + x15;
auto x16 = -V{1}*cell[14];
auto x17 = x12 + x237*(x13 + x254*(-V{1}*cell[2] + x14 + x16 + x249 + x256 + x5));
auto x257 = ((x17)*(x17));
auto x18 = x257;
auto x19 = x1*x240;
auto x20 = x236*x3;
auto x21 = -V{1}*cell[16];
auto x22 = -V{1}*cell[18];
auto x23 = x19 + x237*(x20 + x254*(-V{1}*cell[3] + x14 + x21 + x22 + x246 + x7));
auto x258 = ((x23)*(x23));
auto x24 = x11 + x18 + x258;
auto x25 = V{0.5}*x0*x24;
auto x270 = V{1} - x25;
auto x26 = V{0.333333333333333};
auto x27 = x259;
auto x28 = V{0.5}*x27;
auto x29 = x261;
auto x267 = V{1} - x29;
auto x30 = -x29;
auto x31 = V{1}/x243;
auto x32 = x10*x9;
auto x33 = cell[6] + x21 - x23*x32 + x251;
auto x34 = ((x33)*(x33));
auto x35 = cell[4] + x16 - x17*x32 + x252;
auto x36 = x29*(x9 + V{-1});
auto x262 = -x36;
auto x263 = cell[6] + x262;
auto x37 = -x10*x32 + x241 + x244 + x248 + x263;
auto x38 = cell[16] + cell[18] + x247 - x258*x9 + x263;
auto x39 = x17*x9;
auto x40 = cell[8] + x22 - x23*x39 + x256;
auto x41 = cell[14] + cell[4] + cell[8] - x17*x39 + x245 + x250 + x262;
auto x42 = V{2}*V{1} / (x31 + util::sqrt(V{2.828428}*x253*x259*((x242)*(x242))*util::sqrt(V{2}*x34 + V{2}*((x35)*(x35)) + ((x37)*(x37)) + ((x38)*(x38)) + V{2}*((x40)*(x40)) + ((x41)*(x41))) + ((x31)*(x31))));
auto x264 = V{1} / (x42);
auto x44 = V{2}*x260*x264*x9*(x225 + V{-1});
auto x45 = x225*x38 + x233*x44;
auto x46 = x30*x45;
auto x47 = x225*x41 + x231*x44;
auto x48 = x30*x47;
auto x49 = x225*x37 + x228*x44;
auto x50 = x30*x49;
auto x51 = x48 + x50;
auto x52 = V{1} - x42;
auto x53 = x26*(x28*x52*(x46 + x51) - x9*(x25 + V{-1}) + V{-1});
auto x54 = x0*x27;
auto x55 = V{0.166667}*x54;
auto x300 = x30*x55;
auto x56 = -x261;
auto x57 = -x56;
auto x265 = x55*x57;
auto x58 = -x265;
auto x59 = x11*x17;
auto x60 = x10*x257;
auto x61 = x59 - x60;
auto x272 = x55*x61;
auto x62 = V{0.5}*x54;
auto x266 = x57*x62;
auto x301 = x30*x62;
auto x63 = x266;
auto x64 = x59 + x60;
auto x273 = x62*x64;
auto x65 = x58*x61 + x63*x64;
auto x66 = x10*x258;
auto x67 = x17*x258;
auto x68 = x66 - x67;
auto x284 = x55*x68;
auto x285 = x284*(x30 - x57);
auto x69 = x66 + x67;
auto x280 = x62*x69;
auto x70 = x63*x69;
auto x271 = x270 + x70;
auto x71 = -V{1}*x10;
auto x72 = V{0.0555555555555556};
auto x268 = -x72;
auto x269 = x52*x72;
auto x73 = x72*x9;
auto x74 = -V{1}*x30;
auto x75 = 0;
auto x76 = x225*x33 + x230*x44;
auto x288 = x27*x76;
auto x77 = V{2}*x23;
auto x78 = x10*x45 + x76*x77;
auto x79 = x225*x35 + x229*x44;
auto x275 = x27*x79;
auto x80 = V{2}*x17;
auto x81 = x10*x47 + x79*x80;
auto x82 = x78 - x81;
auto x83 = V{2}*x62*x74;
auto x84 = x78 + x81;
auto x85 = x267;
auto x281 = -x85;
auto x282 = x55*(x281 + x30);
auto x305 = x284*x85;
auto x86 = x49*x85;
auto x87 = x48 + x86;
auto x88 = x28*(x46 + x87);
auto x89 = x268 + x269*(x75*x82 + x83*x84 + x88) + x73*(x0*x71 + x265*x68 + x271 + x28*((x71)*(x71)) + x65);
auto x90 = x74*(x272 + x273);
auto x91 = -V{1}*x17;
auto x92 = V{2}*x10;
auto x93 = x17*x49 + x79*x92;
auto x94 = x225*x40 + x232*x44;
auto x295 = x27*x94;
auto x95 = x17*x45 + x77*x94;
auto x96 = x93 - x95;
auto x97 = x93 + x95;
auto x98 = x47*x85;
auto x99 = x50 + x98;
auto x100 = x28*(x46 + x99);
auto x101 = x268 + x269*(x100 + x75*x96 + x83*x97) + x73*(x0*x91 + x271 + x28*((x91)*(x91)) + x58*x68 + x90);
auto x102 = x18*x23;
auto x103 = x11*x23;
auto x104 = x102 - x103;
auto x105 = x102 + x103;
auto x106 = -V{1}*x23;
auto x107 = x23*x47 + x80*x94;
auto x108 = x23*x49 + x76*x92;
auto x109 = x107 - x108;
auto x110 = x107 + x108;
auto x111 = x45*x85;
auto x112 = x28*(x111 + x51);
auto x113 = x268 + x269*(x109*x75 + x110*x83 + x112) + x73*(x0*x106 + x104*x75 + x105*x83 + x270 + x28*((x106)*(x106)));
auto x114 = -V{1}*x85;
auto x274 = x55*(x114 - x74);
auto x115 = -x56 + V{-1};
auto x287 = x115*x284;
auto x116 = x71 + x91;
auto x117 = V{0.0277777777777778};
auto x276 = -x117;
auto x279 = x117*x52;
auto x118 = x117*x9;
auto x119 = -x274;
auto x120 = x119*x82;
auto x121 = x274;
auto x122 = x121*x96;
auto x123 = x62*(x114 + x74);
auto x124 = x123*x84;
auto x277 = x120 + x124;
auto x125 = x123*x97;
auto x278 = x122 + x125;
auto x126 = x28*(x46 + x86 + x98);
auto x127 = x126 + x275;
auto x128 = x118*(x0*x116 + V{2}*x266*x69 + x270 + x272*(x114 - x115) + x273*(x114 + x115) + x28*((x116)*(x116))) + x276 + x279*(x127 + x277 + x278);
auto x129 = x280*(x30 + x57);
auto x286 = x129 + x270;
auto x130 = x17 + x71;
auto x131 = -x282;
auto x132 = x131*x96;
auto x133 = x62*(x30 + x85);
auto x134 = x133*x97;
auto x283 = x132 + x134;
auto x135 = x126 - V{1}*x275;
auto x136 = x118*(x0*x130 - x272*(x115 + x281) + x273*(x115 + x85) + x28*((x130)*(x130)) - x285 + x286) + x276 + x279*(x135 + x277 + x283);
auto x137 = x104*x119;
auto x138 = x287;
auto x292 = x138 + x65;
auto x139 = x105*x123;
auto x140 = x115*x280;
auto x141 = x139 + x140;
auto x291 = x141 + x270;
auto x142 = x106 + x71;
auto x143 = x109*x119;
auto x144 = x121*x82;
auto x289 = x124 + x144;
auto x145 = x110*x123;
auto x290 = x143 + x145;
auto x146 = x28*(x111 + x87);
auto x147 = x146 + x288;
auto x148 = x118*(x0*x142 + x137 + x28*((x142)*(x142)) + x291 + x292) + x276 + x279*(x147 + x289 + x290);
auto x149 = x282;
auto x150 = x104*x149;
auto x151 = x105*x133;
auto x152 = x140 + x151;
auto x294 = x152 + x270;
auto x153 = x23 + x71;
auto x154 = x109*x149;
auto x155 = x110*x133;
auto x293 = x154 + x155;
auto x156 = x146 - V{1}*x288;
auto x157 = x118*(x0*x153 + x150 + x28*((x153)*(x153)) + x292 + x294) + x276 + x279*(x156 + x289 + x293);
auto x158 = x104*x121;
auto x159 = -x287;
auto x298 = x159 + x90;
auto x160 = x106 + x91;
auto x161 = x109*x121;
auto x297 = x145 + x161;
auto x162 = x119*x96;
auto x296 = x125 + x162;
auto x163 = x28*(x111 + x99);
auto x164 = x163 + x295;
auto x165 = x118*(x0*x160 + x158 + x28*((x160)*(x160)) + x291 + x298) + x276 + x279*(x164 + x296 + x297);
auto x166 = x104*x131;
auto x167 = x23 + x91;
auto x168 = x109*x131;
auto x299 = x155 + x168;
auto x169 = x163 - V{1}*x295;
auto x170 = x118*(x0*x167 + x166 + x28*((x167)*(x167)) + x294 + x298) + x276 + x279*(x169 + x296 + x299);
auto x171 = -x300;
auto x172 = x301;
auto x173 = x172*x64;
auto x174 = x171*x61 + x173;
auto x306 = x174 + x270;
auto x175 = x300;
auto x176 = x172*x69;
auto x302 = x176 + x270;
auto x177 = 0;
auto x178 = V{2}*x301;
auto x179 = x268 + x269*(x177*x82 + x178*x84 + x88) + x73*(x0*x10 + x174 + x175*x68 + x255*x28 + x302);
auto x180 = x173 + x175*x61;
auto x309 = x180 + x270;
auto x181 = x268 + x269*(x100 + x177*x96 + x178*x97) + x73*(x0*x17 + x171*x68 + x180 + x257*x28 + x302);
auto x182 = x268 + x269*(x109*x177 + x110*x178 + x112) + x73*(x0*x23 + x104*x177 + x105*x178 + x258*x28 + x270);
auto x312 = x157 + x170 + x182;
auto x183 = x10 + x17;
auto x184 = x149*x82;
auto x185 = x133*x84;
auto x303 = x184 + x185;
auto x186 = x118*(x0*x183 + x177*x68 + x178*x69 + x270 + V{2}*x273*x85 + x28*((x183)*(x183))) + x276 + x279*(x127 + x283 + x303);
auto x187 = x267;
auto x307 = x187*x284;
auto x188 = x10 + x91;
auto x189 = x118*(x0*x188 + x272*(x114 - x187) + x273*(x114 + x187) + x28*((x188)*(x188)) + x285 + x286) + x276 + x279*(x135 + x278 + x303);
auto x190 = x151 + x280*x85;
auto x191 = x10 + x23;
auto x192 = x131*x82;
auto x304 = x185 + x192;
auto x193 = x118*(x0*x191 + x150 + x190 + x28*((x191)*(x191)) + x305 + x306) + x276 + x279*(x147 + x293 + x304);
auto x194 = x139 + x187*x280;
auto x195 = x10 + x106;
auto x196 = x118*(x0*x195 + x137 + x194 + x28*((x195)*(x195)) + x306 + x307) + x276 + x279*(x156 + x290 + x304);
auto x310 = x179 + x189 + x196;
auto x197 = x17 + x23;
auto x198 = x149*x96;
auto x308 = x134 + x198;
auto x199 = x118*(x0*x197 + x166 + x190 + x28*((x197)*(x197)) - x305 + x309) + x276 + x279*(x164 + x299 + x308);
auto x200 = x106 + x17;
auto x201 = x118*(x0*x200 + x158 + x194 + x28*((x200)*(x200)) - x307 + x309) + x276 + x279*(x169 + x297 + x308);
auto x311 = x136 + x181 + x201;
auto x202 = x101 + x113 + x128 + x148 + x165 + x186 + x193 + x199 + x310 + x311 + x312 + x53 + x89 + V{1};
auto x315 = x227/x202;
auto x203 = -V{1}*x128;
auto x313 = x186 + x203;
auto x204 = -V{1}*x148;
auto x314 = x193 + x204;
auto x205 = -V{1}*x165;
auto x316 = x199 + x205;
cell[0] = x53;
cell[10] = x179;
cell[11] = x181;
cell[12] = x182;
cell[13] = x186;
cell[14] = x189;
cell[15] = x193;
cell[16] = x196;
cell[17] = x199;
cell[18] = x201;
cell[1] = x89;
cell[2] = x101;
cell[3] = x113;
cell[4] = x128;
cell[5] = x136;
cell[6] = x148;
cell[7] = x157;
cell[8] = x165;
cell[9] = x170;
cell.template getFieldPointer<olb::descriptors::OMEGA>()[0] = x42;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x202;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x2 + x237*(x315*(-V{1}*x136 - V{1}*x157 + x310 + x313 + x314 - V{1}*x89) + x4);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x12 + x237*(x13 + x315*(-V{1}*x101 - V{1}*x170 - V{1}*x189 + x311 + x313 + x316));
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x19 + x237*(x20 + x315*(-V{1}*x113 - V{1}*x196 - V{1}*x201 + x312 + x314 + x316));
return { x9, x24 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
