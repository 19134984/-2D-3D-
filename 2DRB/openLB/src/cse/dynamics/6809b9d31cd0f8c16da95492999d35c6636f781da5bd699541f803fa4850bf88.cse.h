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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::WMBulkDensity, momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::ParameterFromCell<collision::HYBRID, collision::ParameterFromCell<collision::HYBRID_RHO, collision::LocalVanDriestSmagorinskyEffectiveOmega<collision::ExternalRhoHRR> > > >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x198 = cell.template getFieldComponent<olb::collision::HYBRID>(0);
auto x199 = cell.template getFieldComponent<olb::collision::HYBRID_RHO>(0);
auto x200 = cell.template getFieldComponent<olb::descriptors::DENSITY>(0);
auto x202 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x203 = cell.template getFieldComponent<olb::descriptors::TENSOR>(0);
auto x204 = cell.template getFieldComponent<olb::descriptors::TENSOR>(1);
auto x205 = cell.template getFieldComponent<olb::descriptors::TENSOR>(2);
auto x206 = cell.template getFieldComponent<olb::descriptors::TENSOR>(3);
auto x207 = cell.template getFieldComponent<olb::descriptors::TENSOR>(4);
auto x208 = cell.template getFieldComponent<olb::descriptors::TENSOR>(5);
auto x209 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x210 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x211 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x212 = cell.template getFieldComponent<olb::descriptors::VISCOSITY>(0);
auto x213 = cell.template getFieldComponent<olb::descriptors::WMPOROSITY>(0);
auto x214 = cell.template getFieldComponent<olb::descriptors::WMVELOCITY>(0);
auto x215 = cell.template getFieldComponent<olb::descriptors::WMVELOCITY>(1);
auto x216 = cell.template getFieldComponent<olb::descriptors::WMVELOCITY>(2);
auto x219 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x220 = parameters.template get<descriptors::OMEGA>();
auto x201 = cell[10] + cell[14] + cell[16];
auto x217 = cell[1] + cell[4] + x201;
auto x218 = cell[15] + cell[7];
auto x221 = cell[17] + cell[9];
auto x222 = cell[12] + x218 + x221;
auto x223 = cell[3] + cell[8] + x222;
auto x224 = cell[13] + cell[5];
auto x225 = cell[11] + cell[18] + x224;
auto x226 = cell[2] + x225;
auto x0 = V{3};
auto x235 = ((x0)*(x0));
auto x236 = V{1} / (x0);
auto x237 = V{1}*x236;
auto x1 = V{1} - x213;
auto x2 = V{1} - x202;
auto x3 = -V{1}*cell[4];
auto x4 = -V{1}*cell[5];
auto x228 = cell[13] + x4;
auto x5 = -V{1}*cell[6];
auto x6 = -V{1}*cell[7];
auto x227 = cell[15] + x6;
auto x7 = x1*x200 + x213*(cell[0] + cell[6] + x217 + x223 + x226 + V{1});
auto x229 = V{1} / (x7);
auto x230 = x202*x229;
auto x8 = x1*x214 + x213*(x2*x209 + x230*(-V{1}*cell[1] + x201 + x227 + x228 + x3 + x5));
auto x231 = ((x8)*(x8));
auto x9 = x231;
auto x10 = -V{1}*cell[8];
auto x11 = -V{1}*cell[9];
auto x232 = cell[17] + x11;
auto x12 = -V{1}*cell[14];
auto x13 = x1*x215 + x213*(x2*x210 + x230*(-V{1}*cell[2] + x10 + x12 + x225 + x232 + x3));
auto x233 = ((x13)*(x13));
auto x14 = x233;
auto x15 = -V{1}*cell[16];
auto x16 = -V{1}*cell[18];
auto x17 = x1*x216 + x213*(x2*x211 + x230*(-V{1}*cell[3] + x10 + x15 + x16 + x222 + x5));
auto x234 = ((x17)*(x17));
auto x18 = x14 + x234 + x9;
auto x19 = V{0.5}*x0*x18;
auto x259 = V{1} - x19;
auto x20 = V{0.333333333333333};
auto x21 = x199*x7 - x200*(x199 + V{-1});
auto x22 = x235;
auto x23 = V{0.5}*x22;
auto x24 = x237;
auto x243 = V{1} - x24;
auto x25 = -x24;
auto x26 = V{1}/x220;
auto x27 = x7*x8;
auto x28 = cell[6] + x15 - x17*x27 + x227;
auto x29 = ((x28)*(x28));
auto x30 = cell[4] + x12 - x13*x27 + x228;
auto x31 = x24*(x7 + V{-1});
auto x238 = -x31;
auto x239 = cell[6] + x238;
auto x32 = x217 + x218 + x224 + x239 - x27*x8;
auto x33 = cell[16] + cell[18] + x223 - x234*x7 + x239;
auto x34 = x13*x7;
auto x35 = cell[8] + x16 - x17*x34 + x232;
auto x36 = cell[14] + cell[4] + cell[8] - x13*x34 + x221 + x226 + x238;
auto x37 = V{1};
auto x38 = V{1}*V{1} / (x37*(x0*x212 + V{0.5}) - (V{0.5}*x37 + V{-0.5})*(x26 + util::sqrt(V{2.828428}*x229*x235*((x219)*(x219))*util::sqrt(V{2}*x29 + V{2}*((x30)*(x30)) + ((x32)*(x32)) + ((x33)*(x33)) + V{2}*((x35)*(x35)) + ((x36)*(x36))) + ((x26)*(x26)))));
auto x240 = V{1} / (x38);
auto x40 = V{2}*x21*x236*x240*(x198 + V{-1});
auto x41 = x198*x33 + x208*x40;
auto x42 = x25*x41;
auto x43 = x198*x36 + x206*x40;
auto x44 = x25*x43;
auto x45 = x198*x32 + x203*x40;
auto x46 = x25*x45;
auto x47 = x44 + x46;
auto x48 = V{1} - x38;
auto x49 = x0*x22;
auto x50 = V{0.166667}*x49;
auto x255 = x25*x50;
auto x51 = -x237;
auto x52 = -x51;
auto x241 = x50*x52;
auto x53 = -x241;
auto x54 = x13*x9;
auto x55 = x233*x8;
auto x56 = x54 - x55;
auto x244 = x50*x56;
auto x57 = V{0.5}*x49;
auto x242 = x52*x57;
auto x256 = x25*x57;
auto x58 = x242;
auto x59 = x54 + x55;
auto x245 = x57*x59;
auto x60 = x53*x56 + x58*x59;
auto x61 = x234*x8;
auto x62 = x13*x234;
auto x63 = x61 - x62;
auto x251 = x50*x63;
auto x266 = x251*(x25 - x52);
auto x64 = x61 + x62;
auto x248 = x57*x64;
auto x65 = x58*x64;
auto x278 = x259 + x65;
auto x66 = -V{1}*x8;
auto x67 = V{0.0555555555555556};
auto x257 = -x67;
auto x258 = x48*x67;
auto x68 = x21*x67;
auto x69 = -V{1}*x25;
auto x70 = 0;
auto x71 = x198*x28 + x205*x40;
auto x253 = x22*x71;
auto x72 = V{2}*x17;
auto x73 = x41*x8 + x71*x72;
auto x74 = x198*x30 + x204*x40;
auto x247 = x22*x74;
auto x75 = V{2}*x13;
auto x76 = x43*x8 + x74*x75;
auto x77 = x73 - x76;
auto x78 = V{2}*x57*x69;
auto x79 = x73 + x76;
auto x80 = x243;
auto x249 = -x80;
auto x250 = x50*(x249 + x25);
auto x270 = x251*x80;
auto x81 = x45*x80;
auto x82 = x44 + x81;
auto x83 = x23*(x42 + x82);
auto x84 = x69*(x244 + x245);
auto x85 = -V{1}*x13;
auto x86 = V{2}*x8;
auto x87 = x13*x45 + x74*x86;
auto x88 = x198*x35 + x207*x40;
auto x254 = x22*x88;
auto x89 = x13*x41 + x72*x88;
auto x90 = x87 - x89;
auto x91 = x87 + x89;
auto x92 = x43*x80;
auto x93 = x46 + x92;
auto x94 = x23*(x42 + x93);
auto x95 = x14*x17;
auto x96 = x17*x9;
auto x97 = x95 - x96;
auto x98 = x95 + x96;
auto x99 = -V{1}*x17;
auto x100 = x17*x43 + x75*x88;
auto x101 = x17*x45 + x71*x86;
auto x102 = x100 - x101;
auto x103 = x100 + x101;
auto x104 = x41*x80;
auto x105 = x23*(x104 + x47);
auto x106 = -V{1}*x80;
auto x246 = x50*(x106 - x69);
auto x107 = -x51 + V{-1};
auto x252 = x107*x251;
auto x108 = x66 + x85;
auto x109 = V{0.0277777777777778};
auto x261 = -x109;
auto x264 = x109*x48;
auto x110 = x109*x21;
auto x111 = -x246;
auto x112 = x111*x77;
auto x113 = x246;
auto x114 = x113*x90;
auto x115 = x57*(x106 + x69);
auto x116 = x115*x79;
auto x279 = x112 + x116;
auto x117 = x115*x91;
auto x265 = x114 + x117;
auto x118 = x23*(x42 + x81 + x92);
auto x119 = x118 + x247;
auto x120 = x248*(x25 + x52);
auto x267 = x120 + x259;
auto x121 = x13 + x66;
auto x122 = -x250;
auto x123 = x122*x90;
auto x124 = x57*(x25 + x80);
auto x125 = x124*x91;
auto x263 = x123 + x125;
auto x126 = x118 - V{1}*x247;
auto x127 = x111*x97;
auto x128 = x252;
auto x282 = x128 + x60;
auto x129 = x115*x98;
auto x130 = x107*x248;
auto x131 = x129 + x130;
auto x281 = x131 + x259;
auto x132 = x66 + x99;
auto x133 = x102*x111;
auto x134 = x113*x77;
auto x280 = x116 + x134;
auto x135 = x103*x115;
auto x272 = x133 + x135;
auto x136 = x23*(x104 + x82);
auto x137 = x136 + x253;
auto x138 = x250;
auto x139 = x138*x97;
auto x140 = x124*x98;
auto x141 = x130 + x140;
auto x283 = x141 + x259;
auto x142 = x17 + x66;
auto x143 = x102*x138;
auto x144 = x103*x124;
auto x269 = x143 + x144;
auto x145 = x136 - V{1}*x253;
auto x146 = x113*x97;
auto x147 = -x252;
auto x285 = x147 + x84;
auto x148 = x85 + x99;
auto x149 = x102*x113;
auto x277 = x135 + x149;
auto x150 = x111*x90;
auto x284 = x117 + x150;
auto x151 = x23*(x104 + x93);
auto x152 = x151 + x254;
auto x153 = x122*x97;
auto x154 = x17 + x85;
auto x155 = x102*x122;
auto x275 = x144 + x155;
auto x156 = x151 - V{1}*x254;
auto x157 = -x255;
auto x158 = x256;
auto x159 = x158*x59;
auto x160 = x157*x56 + x159;
auto x271 = x160 + x259;
auto x161 = x255;
auto x162 = x158*x64;
auto x260 = x162 + x259;
auto x163 = 0;
auto x164 = V{2}*x256;
auto x165 = x159 + x161*x56;
auto x276 = x165 + x259;
auto x166 = x13 + x8;
auto x167 = x138*x77;
auto x168 = x124*x79;
auto x262 = x167 + x168;
auto x169 = x243;
auto x273 = x169*x251;
auto x170 = x8 + x85;
auto x171 = x140 + x248*x80;
auto x172 = x17 + x8;
auto x173 = x122*x77;
auto x268 = x168 + x173;
auto x174 = x129 + x169*x248;
auto x175 = x8 + x99;
auto x176 = x13 + x17;
auto x177 = x138*x90;
auto x274 = x125 + x177;
auto x178 = x13 + x99;
cell[0] = x20*(-x21*(x19 + V{-1}) + x23*x48*(x42 + x47) + V{-1});
cell[10] = x257 + x258*(x163*x77 + x164*x79 + x83) + x68*(x0*x8 + x160 + x161*x63 + x23*x231 + x260);
cell[11] = x257 + x258*(x163*x90 + x164*x91 + x94) + x68*(x0*x13 + x157*x63 + x165 + x23*x233 + x260);
cell[12] = x257 + x258*(x102*x163 + x103*x164 + x105) + x68*(x0*x17 + x163*x97 + x164*x98 + x23*x234 + x259);
cell[13] = x110*(x0*x166 + x163*x63 + x164*x64 + x23*((x166)*(x166)) + V{2}*x245*x80 + x259) + x261 + x264*(x119 + x262 + x263);
cell[14] = x110*(x0*x170 + x23*((x170)*(x170)) + x244*(x106 - x169) + x245*(x106 + x169) + x266 + x267) + x261 + x264*(x126 + x262 + x265);
cell[15] = x110*(x0*x172 + x139 + x171 + x23*((x172)*(x172)) + x270 + x271) + x261 + x264*(x137 + x268 + x269);
cell[16] = x110*(x0*x175 + x127 + x174 + x23*((x175)*(x175)) + x271 + x273) + x261 + x264*(x145 + x268 + x272);
cell[17] = x110*(x0*x176 + x153 + x171 + x23*((x176)*(x176)) - x270 + x276) + x261 + x264*(x152 + x274 + x275);
cell[18] = x110*(x0*x178 + x146 + x174 + x23*((x178)*(x178)) - x273 + x276) + x261 + x264*(x156 + x274 + x277);
cell[1] = x257 + x258*(x70*x77 + x78*x79 + x83) + x68*(x0*x66 + x23*((x66)*(x66)) + x241*x63 + x278 + x60);
cell[2] = x257 + x258*(x70*x90 + x78*x91 + x94) + x68*(x0*x85 + x23*((x85)*(x85)) + x278 + x53*x63 + x84);
cell[3] = x257 + x258*(x102*x70 + x103*x78 + x105) + x68*(x0*x99 + x23*((x99)*(x99)) + x259 + x70*x97 + x78*x98);
cell[4] = x110*(x0*x108 + x23*((x108)*(x108)) + V{2}*x242*x64 + x244*(x106 - x107) + x245*(x106 + x107) + x259) + x261 + x264*(x119 + x265 + x279);
cell[5] = x110*(x0*x121 + x23*((x121)*(x121)) - x244*(x107 + x249) + x245*(x107 + x80) - x266 + x267) + x261 + x264*(x126 + x263 + x279);
cell[6] = x110*(x0*x132 + x127 + x23*((x132)*(x132)) + x281 + x282) + x261 + x264*(x137 + x272 + x280);
cell[7] = x110*(x0*x142 + x139 + x23*((x142)*(x142)) + x282 + x283) + x261 + x264*(x145 + x269 + x280);
cell[8] = x110*(x0*x148 + x146 + x23*((x148)*(x148)) + x281 + x285) + x261 + x264*(x152 + x277 + x284);
cell[9] = x110*(x0*x154 + x153 + x23*((x154)*(x154)) + x283 + x285) + x261 + x264*(x156 + x275 + x284);
cell.template getFieldPointer<olb::descriptors::OMEGA>()[0] = x38;
return { x21, x18 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
