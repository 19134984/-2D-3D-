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
struct CSE<Dual<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<1, 1>, momenta::DefineSeparately> > >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x185 = cell.template getFieldComponent<olb::opti::DJDF>(5);
auto x196 = cell.template getFieldComponent<olb::opti::DJDF>(16);
auto x177 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x192 = cell.template getFieldComponent<olb::opti::DJDF>(12);
auto x180 = cell.template getFieldComponent<olb::opti::DJDF>(0);
auto x186 = cell.template getFieldComponent<olb::opti::DJDF>(6);
auto x197 = cell.template getFieldComponent<olb::opti::DJDF>(17);
auto x182 = cell.template getFieldComponent<olb::opti::DJDF>(2);
auto x188 = cell.template getFieldComponent<olb::opti::DJDF>(8);
auto x194 = cell.template getFieldComponent<olb::opti::DJDF>(14);
auto x187 = cell.template getFieldComponent<olb::opti::DJDF>(7);
auto x191 = cell.template getFieldComponent<olb::opti::DJDF>(11);
auto x189 = cell.template getFieldComponent<olb::opti::DJDF>(9);
auto x193 = cell.template getFieldComponent<olb::opti::DJDF>(13);
auto x218 = parameters.template get<descriptors::OMEGA>();
auto x181 = cell.template getFieldComponent<olb::opti::DJDF>(1);
auto x184 = cell.template getFieldComponent<olb::opti::DJDF>(4);
auto x178 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x183 = cell.template getFieldComponent<olb::opti::DJDF>(3);
auto x179 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x198 = cell.template getFieldComponent<olb::opti::DJDF>(18);
auto x190 = cell.template getFieldComponent<olb::opti::DJDF>(10);
auto x195 = cell.template getFieldComponent<olb::opti::DJDF>(15);
auto x200 = ((x177)*(x177));
auto x201 = ((x178)*(x178));
auto x202 = ((x179)*(x179));
auto x0 = V{1}/(x178 + V{1});
auto x1 = V{0.0277777777777778};
auto x2 = -V{1}*x179;
auto x3 = x178 + x2;
auto x4 = V{3};
auto x199 = ((x4)*(x4));
auto x5 = V{0.5}*x199;
auto x6 = V{0.5}*x4*(x200 + x201 + x202);
auto x203 = V{1} - x6;
auto x7 = x203 + x3*x4 + x5*((x3)*(x3));
auto x212 = x1*x7;
auto x228 = V{1}*x212;
auto x8 = x178 + x179;
auto x9 = x203 + x4*x8 + x5*((x8)*(x8));
auto x211 = x1*x9;
auto x229 = V{1}*x211;
auto x10 = x177 + x2;
auto x11 = x10*x4 + x203 + x5*((x10)*(x10));
auto x210 = x1*x11;
auto x12 = x177 + x179;
auto x13 = x12*x4 + x203 + x5*((x12)*(x12));
auto x209 = x1*x13;
auto x14 = -V{1}*x178;
auto x15 = x14 + x177;
auto x16 = x177 + x178;
auto x17 = x16*x4 + x203 + x5*((x16)*(x16));
auto x208 = x1*x17;
auto x226 = V{1}*x208;
auto x18 = V{0.0555555555555556};
auto x19 = x179*x4 + x202*x5 + x203;
auto x207 = x18*x19;
auto x20 = x178*x4 + x201*x5 + x203;
auto x206 = x18*x20;
auto x230 = V{1}*x206;
auto x21 = x177*x4 + x200*x5 + x203;
auto x205 = x18*x21;
auto x22 = x14 + x179;
auto x23 = x14 + x2;
auto x24 = -V{1}*x177;
auto x25 = x179 + x24;
auto x26 = x203 + x25*x4 + x5*((x25)*(x25));
auto x217 = x1*x26;
auto x27 = x2 + x24;
auto x28 = x203 + x27*x4 + x5*((x27)*(x27));
auto x216 = x1*x28;
auto x29 = x178 + x24;
auto x30 = x203 + x29*x4 + x5*((x29)*(x29));
auto x215 = x1*x30;
auto x227 = V{1}*x215;
auto x31 = x14 + x24;
auto x32 = x2*x4 + x203 + x5*((x2)*(x2));
auto x214 = x18*x32;
auto x33 = x203 + x24*x4 + x5*((x24)*(x24));
auto x213 = x18*x33;
auto x34 = V{0.333333333333333};
auto x35 = V{1} - x6;
auto x204 = x34*x35;
auto x36 = x0*(cell[0]*x204 + cell[10]*x205 + cell[11]*x206 + cell[12]*x207 + cell[13]*x208 + cell[14]*x1*(x15*x4 + x203 + x5*((x15)*(x15))) + cell[15]*x209 + cell[16]*x210 + cell[17]*x211 + cell[18]*x212 + cell[1]*x213 + cell[2]*x18*(x14*x4 + x203 + x5*((x14)*(x14))) + cell[3]*x214 + cell[4]*x1*(x203 + x31*x4 + x5*((x31)*(x31))) + cell[5]*x215 + cell[6]*x216 + cell[7]*x217 + cell[8]*x1*(x203 + x23*x4 + x5*((x23)*(x23))) + cell[9]*x1*(x203 + x22*x4 + x5*((x22)*(x22))));
auto x37 = V{2}*x36;
auto x38 = V{-2};
auto x39 = -V{1}*x38;
auto x40 = V{1}/x4;
auto x41 = V{1} - x40;
auto x42 = V{0.5}*x1*x199;
auto x43 = V{1} - x218;
auto x219 = x42*x43;
auto x44 = x219*(cell[18] + cell[9]);
auto x45 = x219*(cell[17] + cell[8]);
auto x46 = x41*(x44 + x45);
auto x47 = x219*(cell[16] + cell[7]);
auto x48 = x41*x47;
auto x49 = x219*(cell[15] + cell[6]);
auto x50 = x41*x49;
auto x51 = -x40;
auto x52 = x219*(cell[14] + cell[5]);
auto x53 = x219*(cell[13] + cell[4]);
auto x54 = V{0.5}*x18*x199;
auto x220 = x43*x54;
auto x55 = x220*(cell[12] + cell[3]);
auto x56 = x220*(cell[11] + cell[2]);
auto x57 = x51*x56;
auto x223 = x48 + x50 + x57;
auto x58 = x220*(cell[10] + cell[1]);
auto x59 = x51*x58;
auto x60 = V{0.5}*x199*x34;
auto x221 = x51*x60;
auto x61 = cell[0]*x221*x43;
auto x222 = x46 + x59 + x61;
auto x62 = x222 + x223 + x41*x55 + x51*x52 + x51*x53;
auto x63 = x38*x44 + V{2}*x45;
auto x64 = x41*x52;
auto x65 = x41*x53;
auto x66 = x51*x55;
auto x224 = x64 + x65 + x66;
auto x67 = V{2}*x222 + V{2}*x224 + V{2}*x41*x56 + V{2}*x47*x51 + V{2}*x49*x51;
auto x68 = x38*x63 + x39*x62 + x67;
auto x69 = V{2}*x62 + V{2}*x63 + x67;
auto x70 = x38*x52 + V{2}*x53;
auto x71 = x223 + x224 + x41*x58 + x44*x51 + x45*x51 + x61;
auto x225 = x62 + x71;
auto x72 = x67 + V{2}*x70 + V{2}*x71;
auto x73 = x38*x70 + x39*x71 + x67;
auto x74 = x38*x47 + V{2}*x49;
auto x75 = x225 - V{1}*x74;
auto x76 = -V{1}*x75;
auto x77 = x225 + x74;
auto x78 = -V{1}*x77;
auto x79 = -V{1}*x62;
auto x80 = -V{1}*x71;
auto x81 = x0*(x205*x80 + x207*x79 + x209*x78 + x210*x76 + x213*x80 + x214*x79 + x216*x78 + x217*x76 - x226*x72 - x227*x73 - x228*x68 - x229*x69 - x230*x67);
auto x232 = x36 + x81;
auto x82 = V{2}*x81;
auto x231 = x37 + x82;
auto x83 = x231 + x68;
auto x84 = x231 + x69;
auto x85 = x232 + x75;
auto x86 = x232 + x77;
auto x87 = x231 + x72;
auto x88 = x232 + x62;
auto x89 = x231 + x67;
auto x90 = x232 + x71;
auto x91 = x231 + x73;
auto x92 = x232;
auto x93 = x54*x90;
auto x94 = x51*x93;
auto x95 = x54*x88;
auto x96 = x41*x95;
auto x97 = x42*x86;
auto x98 = x41*x97;
auto x99 = x42*x85;
auto x100 = x41*x99;
auto x101 = x42*x83;
auto x102 = x42*x84;
auto x103 = x41*(x101 + x102);
auto x104 = x42*x87;
auto x105 = x54*x89;
auto x106 = x105*x51;
auto x107 = x42*x91;
auto x108 = x221*x92;
auto x233 = V{2}*x100 + x106 + x108 + V{2}*x98;
auto x109 = x103 + x104*x51 + x107*x51 + x233 + V{2}*x94 + V{2}*x96;
auto x110 = x101*x38 + V{2}*x102;
auto x111 = x51*x95;
auto x112 = x51*x97;
auto x113 = x51*x99;
auto x114 = x104*x41;
auto x115 = x107*x41;
auto x116 = V{2}*x103 + V{2}*x105*x41 + V{2}*x108 + V{4}*x111 + V{4}*x112 + V{4}*x113 + V{2}*x114 + V{2}*x115 + V{4}*x94;
auto x117 = x109*x39 + x110*x38 + x116;
auto x118 = V{2}*x109 + V{2}*x110 + x116;
auto x119 = V{2}*x104 + x107*x38;
auto x120 = x41*x93;
auto x121 = x101*x51 + x102*x51 + V{2}*x111 + x114 + x115 + V{2}*x120 + x233;
auto x234 = x109 + x121;
auto x122 = x116 + V{2}*x119 + V{2}*x121;
auto x123 = x116 + x119*x38 + x121*x39;
auto x124 = V{2}*x97;
auto x125 = x38*x99;
auto x126 = V{2}*x124 + V{2}*x125;
auto x127 = -V{1}*x126 + x234;
auto x128 = -V{1}*x127;
auto x129 = x126 + x234;
auto x130 = -V{1}*x129;
auto x131 = -V{1}*x109;
auto x132 = -V{1}*x121;
auto x133 = x0*(-x116*x230 - x117*x228 - x118*x229 - x122*x226 - x123*x227 + x128*x210 + x128*x217 + x130*x209 + x130*x216 + x131*x207 + x131*x214 + x132*x205 + x132*x213 + x204*x92 + x205*x90 + x206*x89 + x207*x88 + x208*x87 + x209*x86 + x210*x85 + x211*x84 + x212*x83 + x213*x90 + x214*x88 + x215*x91 + x216*x86 + x217*x85);
auto x134 = x133 - x180;
auto x135 = x121 + x133;
auto x136 = x135 - x181;
auto x137 = -x182;
auto x138 = x109 + x133;
auto x139 = x138 - x183;
auto x140 = -x184;
auto x141 = V{2}*x133;
auto x142 = x123 + x141 - x185;
auto x143 = x129 + x133;
auto x144 = x143 - x186;
auto x145 = x127 + x133;
auto x146 = x145 - x187;
auto x147 = -x188;
auto x148 = -x189;
auto x149 = x135 - x190;
auto x150 = x116 + x141 - x191;
auto x151 = x138 - x192;
auto x152 = x122 + x141 - x193;
auto x153 = -x194;
auto x154 = x143 - x195;
auto x155 = x145 - x196;
auto x156 = x118 + x141 - x197;
auto x157 = x117 + x141 - x198;
cell[0] = x134;
cell[1] = x136;
cell[2] = x137;
cell[3] = x139;
cell[4] = x140;
cell[5] = x142;
cell[6] = x144;
cell[7] = x146;
cell[8] = x147;
cell[9] = x148;
cell[10] = x149;
cell[11] = x150;
cell[12] = x151;
cell[13] = x152;
cell[14] = x153;
cell[15] = x154;
cell[16] = x155;
cell[17] = x156;
cell[18] = x157;
return { x134 + V{1}, ((x134)*(x134)) + ((x136)*(x136)) + ((x137)*(x137)) + ((x139)*(x139)) + ((x140)*(x140)) + ((x142)*(x142)) + ((x144)*(x144)) + ((x146)*(x146)) + ((x147)*(x147)) + ((x148)*(x148)) + ((x149)*(x149)) + ((x150)*(x150)) + ((x151)*(x151)) + ((x152)*(x152)) + ((x153)*(x153)) + ((x154)*(x154)) + ((x155)*(x155)) + ((x156)*(x156)) + ((x157)*(x157)) };
}
};

}

}

// Generation Info: commit=4ea71c45fdf920fe461e9926f71dc63f7fa1062d
