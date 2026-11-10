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

namespace operators {

template <typename T, typename...FIELDS>
struct CSE_O<operators::TrackAverageFunctorO<functors::DerivativeF<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, descriptors::POPULATION, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination> >, opti::AVERAGE_DJDF>,descriptors::D3Q19<FIELDS...>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
void apply(CELL& cell, PARAMETERS& parameters) any_platform {
using V = typename CELL::value_t;
auto x187 = parameters.template get<parameters::OMEGA>();
auto x193 = parameters.template get<parameters::FACTOR>();
auto x194 = parameters.template get<descriptors::LATTICE_TIME_T>();
auto x189 = parameters.template get<parameters::LATTICE_VISCOSITY>();
auto x190 = parameters.template get<parameters::CONVERSION_VELOCITY>();
auto x188 = parameters.template get<parameters::DX>();
auto x192 = parameters.template get<parameters::DT>();
auto x191 = parameters.template get<parameters::PHYS_CHAR_VISCOSITY>();
auto x195 = cell[15] + cell[16];
auto x196 = cell[10] + cell[14] + x195;
auto x197 = cell[1] + cell[4] + x196;
auto x198 = cell[12] + cell[7] + cell[9];
auto x199 = cell[3] + cell[8] + x198;
auto x200 = cell[13] + cell[5];
auto x201 = cell[17] + cell[18];
auto x202 = cell[11] + x200 + x201;
auto x203 = cell[2] + x202;
auto x216 = V{1} / (x194);
auto x220 = ((x190)*(x190));
auto x0 = V{0.5};
auto x219 = x0*x191;
auto x1 = -V{1}*cell[16];
auto x204 = cell[15] + x1;
auto x2 = -V{1}*cell[7];
auto x3 = -V{1}*cell[6];
auto x4 = -V{1}*cell[8];
auto x5 = -V{1}*cell[18];
auto x205 = cell[17] + x5;
auto x6 = cell[0] + cell[6] + x197 + x199 + x203 + V{1};
auto x7 = V{1}/x6;
auto x212 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*x7;
auto x8 = x7*(-V{1}*cell[3] + x198 + x204 + x205 + x3 + x4);
auto x9 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*x8;
auto x10 = -V{1}*cell[4];
auto x11 = -V{1}*cell[5];
auto x206 = cell[13] + x11;
auto x12 = x7*(-V{1}*cell[1] + x10 + x196 + x2 + x206 + x3);
auto x13 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*x12;
auto x14 = x13*x6;
auto x15 = cell[6] - x14*x9 + x2 + x204;
auto x16 = ((x15)*(x15));
auto x17 = -V{1}*cell[14];
auto x18 = -V{1}*cell[9];
auto x19 = x7*(-V{1}*cell[2] + x10 + x17 + x18 + x202 + x4);
auto x20 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*x19;
auto x21 = cell[4] - x14*x20 + x17 + x206;
auto x207 = V{2}*x21;
auto x22 = x207;
auto x23 = V{3};
auto x24 = V{1}/x23;
auto x25 = x24*(x6 + V{-1});
auto x208 = -x25;
auto x209 = cell[6] + x208;
auto x26 = cell[7] - x13*x14 + x197 + x200 + x209;
auto x27 = x6*x9;
auto x28 = x195 + x199 + x201 + x209 - x27*x9;
auto x29 = x20*x6;
auto x30 = cell[8] + x18 + x205 - x29*x9;
auto x210 = V{2}*x30;
auto x31 = x210;
auto x213 = x210 + x31;
auto x32 = cell[14] + cell[4] + cell[8] + cell[9] - x20*x29 + x203 + x208;
auto x33 = x187*x23*x7;
auto x34 = V{1}/x192;
auto x211 = x33*x34;
auto x35 = x211;
auto x36 = -V{2}*x211*util::pow(x35, V{1})*x7*(V{2}*x16 + x21*x22 + x30*x31 + ((x26)*(x26)) + ((x28)*(x28)) + ((x32)*(x32)));
auto x37 = V{1}*x8;
auto x38 = -x212*x37;
auto x39 = x14*x38;
auto x40 = V{1}*x12;
auto x41 = -x212*x40;
auto x42 = V{1}*x13;
auto x43 = x41*x6 + x42;
auto x44 = x43*x9;
auto x45 = -x15*(x39 + x44);
auto x46 = V{2}*x45;
auto x47 = V{1}*x19;
auto x48 = -x212*x47;
auto x49 = x14*x48;
auto x50 = x20*x43;
auto x51 = -x49 - x50;
auto x52 = V{1}*x24;
auto x215 = x52 + V{-1};
auto x53 = -x26*(x13*x43 + x14*x41 + x52);
auto x54 = V{2}*x53;
auto x55 = x207*x51 + x22*x51 + x54;
auto x56 = V{1}*x9;
auto x57 = -x28*(x27*x38 + x52 + x9*(x38*x6 + x56));
auto x58 = V{2}*x57;
auto x217 = V{2}*x46 + x58;
auto x59 = x29*x38;
auto x60 = V{1}*x20;
auto x61 = x48*x6 + x60;
auto x62 = x61*x9;
auto x63 = -x59 - x62;
auto x64 = x213*x63;
auto x65 = -x32*(x20*x61 + x29*x48 + x52);
auto x66 = V{2}*x65;
auto x218 = x55 + x66;
auto x67 = util::pow(x35, V{2});
auto x68 = x13*x41;
auto x69 = V{2}*x68;
auto x70 = x20*x48;
auto x71 = V{2}*x70;
auto x72 = x69 + x71;
auto x73 = x38*x9;
auto x74 = V{2}*x73;
auto x75 = -V{1}*x187*x191*(cell.template getFieldComponent<olb::descriptors::POROSITY>(0) + V{-1})/(((x188)*(x188))*x189);
auto x76 = x194 + V{-1};
auto x77 = V{-1};
auto x214 = -x77;
auto x78 = -x212*(x214 + x40);
auto x79 = x42 + x6*x78;
auto x80 = x79*x9;
auto x81 = -x15*(x39 + x80);
auto x82 = V{2}*x81;
auto x222 = x58 + V{2}*x82;
auto x83 = x20*x79;
auto x84 = -x49 - x83;
auto x85 = -x26*(x13*x79 + x14*x78 + x215);
auto x86 = V{2}*x85;
auto x226 = x222 + x86;
auto x87 = x207*x84 + x22*x84 + x86;
auto x221 = x66 + x87;
auto x88 = x13*x78;
auto x89 = V{2}*x88;
auto x228 = x74 + x89;
auto x90 = x71 + x89;
auto x91 = -x212*(x214 + x47);
auto x92 = x14*x91;
auto x93 = -x50 - x92;
auto x94 = x207*x93 + x22*x93 + x54;
auto x95 = x6*x91 + x60;
auto x96 = x9*x95;
auto x97 = -x59 - x96;
auto x98 = x213*x97;
auto x99 = -x32*(x20*x95 + x215 + x29*x91);
auto x100 = V{2}*x99;
auto x223 = x100 + x94;
auto x227 = x100 + x98;
auto x101 = x20*x91;
auto x102 = V{2}*x101;
auto x103 = x102 + x69;
auto x104 = -x212*(x214 + x37);
auto x105 = x104*x14;
auto x106 = -x15*(x105 + x44);
auto x107 = V{2}*x106;
auto x224 = V{2}*x107;
auto x108 = -x28*(x104*x27 + x215 + x9*(x104*x6 + x56));
auto x109 = V{2}*x108;
auto x231 = x109 + x224;
auto x110 = x104*x29;
auto x111 = -x110 - x62;
auto x112 = x111*x213;
auto x225 = x109 + x112;
auto x113 = x104*x9;
auto x114 = V{2}*x113;
auto x115 = -x83 - x92 + V{1};
auto x116 = -x212*(x47 + V{-1});
auto x117 = x116*x14;
auto x118 = -x117 - x214 - x83;
auto x119 = x116*x6 + x60;
auto x120 = x119*x9;
auto x121 = -x120 - x59;
auto x122 = x121*x213;
auto x123 = -x32*(x116*x29 + x119*x20 + x215);
auto x124 = V{2}*x123;
auto x229 = x122 + x124;
auto x125 = x116*x20;
auto x126 = V{2}*x125;
auto x127 = -x15*(x105 + x80 + V{-1});
auto x128 = V{2}*x127;
auto x129 = -x212*(x37 + V{-1});
auto x130 = x129*x14;
auto x131 = -x15*(x130 + x214 + x80);
auto x132 = V{2}*x131;
auto x133 = -x28*(x129*x27 + x215 + x9*(x129*x6 + x56));
auto x134 = V{2}*x133;
auto x135 = x129*x29;
auto x136 = -x135 - x62;
auto x137 = x136*x213;
auto x230 = x134 + x137;
auto x138 = x129*x9;
auto x139 = V{2}*x138;
auto x140 = -x110 - x96 + V{1};
auto x141 = -x15*(x130 + x44);
auto x142 = V{2}*x141;
auto x232 = V{2}*x142;
auto x233 = x134 + x232;
auto x143 = -x135 - x214 - x96;
auto x144 = -x212*(x40 + V{-1});
auto x145 = x144*x6 + x42;
auto x146 = x145*x9;
auto x147 = -x15*(x146 + x39);
auto x148 = V{2}*x147;
auto x235 = V{2}*x148 + x58;
auto x149 = x145*x20;
auto x150 = -x149 - x49;
auto x151 = -x26*(x13*x145 + x14*x144 + x215);
auto x152 = V{2}*x151;
auto x237 = x152 + x235;
auto x153 = x150*x207 + x150*x22 + x152;
auto x234 = x153 + x66;
auto x154 = x13*x144;
auto x155 = V{2}*x154;
auto x238 = x155 + x74;
auto x156 = x155 + x71;
auto x157 = -x117 - x50;
auto x158 = x157*x207 + x157*x22 + x54;
auto x236 = x124 + x158;
auto x159 = x126 + x69;
auto x160 = -x117 - x149 + V{1};
auto x161 = -x149 - x214 - x92;
auto x162 = -x15*(x130 + x146 + V{-1});
auto x163 = V{2}*x162;
auto x164 = -x15*(x105 + x146 + x214);
auto x165 = V{2}*x164;
auto x166 = -x120 - x135 + V{1};
auto x167 = -x110 - x120 - x214;
auto x168 = cell[0];
auto x169 = cell[1];
auto x170 = cell[2];
auto x171 = cell[3];
auto x172 = cell[4];
auto x173 = cell[5];
auto x174 = cell[6];
auto x175 = cell[7];
auto x176 = cell[8];
auto x177 = cell[9];
auto x178 = cell[10];
auto x179 = cell[11];
auto x180 = cell[12];
auto x181 = cell[13];
auto x182 = cell[14];
auto x183 = cell[15];
auto x184 = cell[16];
auto x185 = cell[17];
auto x186 = cell[18];
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[0] = x216*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(0)*x76 + x193*(x219*(x36 + x67*(x217 + x218 + x64)) + x220*x75*(x72 + x74)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[1] = x216*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(1)*x76 + x193*(x219*(x36 + x67*(x221 + x222 + x64)) + x220*x75*(x74 + x90)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[2] = x216*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(2)*x76 + x193*(x219*(x36 + x67*(x217 + x223 + x98)) + x220*x75*(x103 + x74)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[3] = x216*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(3)*x76 + x193*(x219*(x36 + x67*(x218 + x224 + x225)) + x220*x75*(x114 + x72)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[4] = x216*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(4)*x76 + x193*(x219*(x36 + x67*(x115*x207 + x115*x22 + x226 + x227)) + x220*x75*(x102 + x228)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[5] = x216*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(5)*x76 + x193*(x219*(x36 + x67*(x118*x207 + x118*x22 + x226 + x229)) + x220*x75*(x126 + x228)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[6] = x216*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(6)*x76 + x193*(x219*(x36 + x67*(V{2}*x128 + x221 + x225)) + x220*x75*(x114 + x90)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[7] = x216*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(7)*x76 + x193*(x219*(x36 + x67*(V{2}*x132 + x221 + x230)) + x220*x75*(x139 + x90)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[8] = x216*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(8)*x76 + x193*(x219*(x36 + x67*(x140*x210 + x140*x31 + x223 + x231)) + x220*x75*(x103 + x114)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[9] = x216*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(9)*x76 + x193*(x219*(x36 + x67*(x143*x210 + x143*x31 + x223 + x233)) + x220*x75*(x103 + x139)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[10] = x216*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(10)*x76 + x193*(x219*(x36 + x67*(x234 + x235 + x64)) + x220*x75*(x156 + x74)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[11] = x216*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(11)*x76 + x193*(x219*(x36 + x67*(x122 + x217 + x236)) + x220*x75*(x159 + x74)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[12] = x216*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(12)*x76 + x193*(x219*(x36 + x67*(x218 + x230 + x232)) + x220*x75*(x139 + x72)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[13] = x216*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(13)*x76 + x193*(x219*(x36 + x67*(x160*x207 + x160*x22 + x229 + x237)) + x220*x75*(x126 + x238)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[14] = x216*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(14)*x76 + x193*(x219*(x36 + x67*(x161*x207 + x161*x22 + x227 + x237)) + x220*x75*(x102 + x238)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[15] = x216*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(15)*x76 + x193*(x219*(x36 + x67*(V{2}*x163 + x230 + x234)) + x220*x75*(x139 + x156)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[16] = x216*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(16)*x76 + x193*(x219*(x36 + x67*(V{2}*x165 + x225 + x234)) + x220*x75*(x114 + x156)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[17] = x216*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(17)*x76 + x193*(x219*(x36 + x67*(x166*x210 + x166*x31 + x233 + x236)) + x220*x75*(x139 + x159)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[18] = x216*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(18)*x76 + x193*(x219*(x36 + x67*(x167*x210 + x167*x31 + x231 + x236)) + x220*x75*(x114 + x159)));
cell[0] = x168;
cell[1] = x169;
cell[2] = x170;
cell[3] = x171;
cell[4] = x172;
cell[5] = x173;
cell[6] = x174;
cell[7] = x175;
cell[8] = x176;
cell[9] = x177;
cell[10] = x178;
cell[11] = x179;
cell[12] = x180;
cell[13] = x181;
cell[14] = x182;
cell[15] = x183;
cell[16] = x184;
cell[17] = x185;
cell[18] = x186;
}
};

}

}

// Generation Info: commit=21521319fcfc23b85b6fd87340ebc70ad50346d4
