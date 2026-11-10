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
struct CSE<Dual<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<1, -1>, momenta::DefineSeparately> > >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x196 = cell.template getFieldComponent<olb::opti::DJDF>(14);
auto x181 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x197 = cell.template getFieldComponent<olb::opti::DJDF>(15);
auto x180 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x195 = cell.template getFieldComponent<olb::opti::DJDF>(13);
auto x189 = cell.template getFieldComponent<olb::opti::DJDF>(7);
auto x194 = cell.template getFieldComponent<olb::opti::DJDF>(12);
auto x184 = cell.template getFieldComponent<olb::opti::DJDF>(2);
auto x200 = cell.template getFieldComponent<olb::opti::DJDF>(18);
auto x182 = cell.template getFieldComponent<olb::opti::DJDF>(0);
auto x198 = cell.template getFieldComponent<olb::opti::DJDF>(16);
auto x190 = cell.template getFieldComponent<olb::opti::DJDF>(8);
auto x188 = cell.template getFieldComponent<olb::opti::DJDF>(6);
auto x179 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x193 = cell.template getFieldComponent<olb::opti::DJDF>(11);
auto x199 = cell.template getFieldComponent<olb::opti::DJDF>(17);
auto x183 = cell.template getFieldComponent<olb::opti::DJDF>(1);
auto x192 = cell.template getFieldComponent<olb::opti::DJDF>(10);
auto x191 = cell.template getFieldComponent<olb::opti::DJDF>(9);
auto x186 = cell.template getFieldComponent<olb::opti::DJDF>(4);
auto x187 = cell.template getFieldComponent<olb::opti::DJDF>(5);
auto x220 = parameters.template get<descriptors::OMEGA>();
auto x185 = cell.template getFieldComponent<olb::opti::DJDF>(3);
auto x202 = ((x179)*(x179));
auto x203 = ((x180)*(x180));
auto x204 = ((x181)*(x181));
auto x0 = -V{1}*x180;
auto x1 = V{1}/(x0 + V{1});
auto x2 = V{0.0277777777777778};
auto x3 = -V{1}*x181;
auto x4 = x179 + x3;
auto x5 = V{3};
auto x201 = ((x5)*(x5));
auto x6 = V{0.5}*x201;
auto x7 = V{0.5}*x5*(x202 + x203 + x204);
auto x205 = V{1} - x7;
auto x8 = x205 + x4*x5 + x6*((x4)*(x4));
auto x211 = x2*x8;
auto x9 = x180 + x3;
auto x10 = x180 + x181;
auto x11 = x179 + x181;
auto x12 = x11*x5 + x205 + x6*((x11)*(x11));
auto x210 = x12*x2;
auto x13 = x0 + x179;
auto x14 = x13*x5 + x205 + x6*((x13)*(x13));
auto x209 = x14*x2;
auto x228 = V{1}*x209;
auto x15 = x179 + x180;
auto x16 = V{0.0555555555555556};
auto x17 = x181*x5 + x204*x6 + x205;
auto x208 = x16*x17;
auto x18 = x179*x5 + x202*x6 + x205;
auto x207 = x16*x18;
auto x19 = x0 + x181;
auto x20 = x19*x5 + x205 + x6*((x19)*(x19));
auto x219 = x2*x20;
auto x230 = V{1}*x219;
auto x21 = x0 + x3;
auto x22 = x205 + x21*x5 + x6*((x21)*(x21));
auto x218 = x2*x22;
auto x231 = V{1}*x218;
auto x23 = -V{1}*x179;
auto x24 = x181 + x23;
auto x25 = x205 + x24*x5 + x6*((x24)*(x24));
auto x217 = x2*x25;
auto x26 = x23 + x3;
auto x27 = x205 + x26*x5 + x6*((x26)*(x26));
auto x216 = x2*x27;
auto x28 = x180 + x23;
auto x29 = x0 + x23;
auto x30 = x205 + x29*x5 + x6*((x29)*(x29));
auto x215 = x2*x30;
auto x232 = V{1}*x215;
auto x31 = x205 + x3*x5 + x6*((x3)*(x3));
auto x214 = x16*x31;
auto x32 = x0*x5 + x205 + x6*((x0)*(x0));
auto x213 = x16*x32;
auto x229 = V{1}*x213;
auto x33 = x205 + x23*x5 + x6*((x23)*(x23));
auto x212 = x16*x33;
auto x34 = V{0.333333333333333};
auto x35 = V{1} - x7;
auto x206 = x34*x35;
auto x36 = x1*(cell[0]*x206 + cell[10]*x207 + cell[11]*x16*(x180*x5 + x203*x6 + x205) + cell[12]*x208 + cell[13]*x2*(x15*x5 + x205 + x6*((x15)*(x15))) + cell[14]*x209 + cell[15]*x210 + cell[16]*x211 + cell[17]*x2*(x10*x5 + x205 + x6*((x10)*(x10))) + cell[18]*x2*(x205 + x5*x9 + x6*((x9)*(x9))) + cell[1]*x212 + cell[2]*x213 + cell[3]*x214 + cell[4]*x215 + cell[5]*x2*(x205 + x28*x5 + x6*((x28)*(x28))) + cell[6]*x216 + cell[7]*x217 + cell[8]*x218 + cell[9]*x219);
auto x37 = V{1}/x5;
auto x38 = -x37;
auto x39 = V{0.5}*x2*x201;
auto x40 = V{1} - x220;
auto x221 = x39*x40;
auto x41 = x221*(cell[18] + cell[9]);
auto x42 = x221*(cell[17] + cell[8]);
auto x43 = V{1} - x37;
auto x44 = x221*(cell[16] + cell[7]);
auto x45 = x43*x44;
auto x46 = x221*(cell[15] + cell[6]);
auto x47 = x43*x46;
auto x48 = x221*(cell[14] + cell[5]);
auto x49 = x43*x48;
auto x50 = x221*(cell[13] + cell[4]);
auto x51 = x43*x50;
auto x52 = V{0.5}*x16*x201;
auto x222 = x40*x52;
auto x53 = x222*(cell[12] + cell[3]);
auto x54 = x38*x53;
auto x225 = x49 + x51 + x54;
auto x55 = x222*(cell[11] + cell[2]);
auto x56 = x38*x55;
auto x224 = x45 + x47 + x56;
auto x57 = x222*(cell[10] + cell[1]);
auto x58 = V{0.5}*x201*x34;
auto x223 = x38*x58;
auto x59 = cell[0]*x223*x40;
auto x60 = x224 + x225 + x38*x41 + x38*x42 + x43*x57 + x59;
auto x61 = x43*(x41 + x42);
auto x62 = x38*x57;
auto x226 = x59 + x61 + x62;
auto x63 = x224 + x226 + x38*x48 + x38*x50 + x43*x53;
auto x227 = x60 + x63;
auto x64 = V{-2};
auto x65 = x44*x64 + V{2}*x46;
auto x66 = x227 - V{1}*x65;
auto x67 = -V{1}*x64;
auto x68 = x67*(x225 + x226 + x38*x44 + x38*x46 + x43*x55);
auto x69 = x48*x64 + V{2}*x50;
auto x70 = V{2}*x60 + x64*x69 + x68;
auto x71 = x41*x64 + V{2}*x42;
auto x72 = V{2}*x63 + x64*x71 + x68;
auto x73 = x63*x67 + x67*x71 + x68;
auto x74 = x60*x67 + x67*x69 + x68;
auto x75 = -V{1}*x66;
auto x76 = x227 + x65;
auto x77 = -V{1}*x76;
auto x78 = -V{1}*x63;
auto x79 = -V{1}*x60;
auto x80 = x1*(x207*x79 + x208*x78 + x210*x77 + x211*x75 + x212*x79 + x214*x78 + x216*x77 + x217*x75 - x228*x70 - x229*x68 - x230*x72 - x231*x73 - x232*x74);
auto x233 = x36 + x80;
auto x81 = x233 + x66;
auto x82 = x233 + x76;
auto x83 = V{2}*x36;
auto x84 = V{2}*x80;
auto x234 = x83 + x84;
auto x85 = x234 + x70;
auto x86 = x233 + x63;
auto x87 = x233 + x60;
auto x88 = x234 + x72;
auto x89 = x234 + x73;
auto x90 = x234 + x74;
auto x91 = x234 + x68;
auto x92 = x233;
auto x93 = x52*x87;
auto x94 = x38*x93;
auto x95 = x52*x86;
auto x96 = x38*x95;
auto x97 = x39*x82;
auto x98 = x38*x97;
auto x99 = x39*x81;
auto x100 = x38*x99;
auto x101 = x39*x85;
auto x102 = x101*x43;
auto x103 = x39*x88;
auto x104 = x103*x43;
auto x105 = x39*x89;
auto x106 = x105*x43;
auto x236 = x104 + x106 + V{2}*x94;
auto x107 = x39*x90;
auto x108 = x107*x43;
auto x109 = x52*x91;
auto x110 = x223*x92;
auto x235 = x102 + x108 + x110 + V{2}*x96;
auto x111 = x67*(V{2}*x100 + x109*x43 + x235 + x236 + V{2}*x98);
auto x112 = x101*x64 + V{2}*x107;
auto x113 = x43*x93;
auto x114 = x43*x97;
auto x115 = x43*x99;
auto x237 = x114 + x115;
auto x116 = x237;
auto x117 = x109*x38;
auto x238 = x116 + x117 + x237;
auto x118 = x103*x38 + x105*x38 + V{2}*x113 + x235 + x238;
auto x119 = x111 + x112*x64 + V{2}*x118;
auto x120 = x43*x95;
auto x121 = x101*x38 + x107*x38 + x110 + V{2}*x120 + x236 + x238;
auto x239 = x118 + x121;
auto x122 = x103*x64 + V{2}*x105;
auto x123 = x111 + V{2}*x121 + x122*x64;
auto x124 = x111 + x121*x67 + x122*x67;
auto x125 = x111 + x112*x67 + x118*x67;
auto x126 = V{2}*x97;
auto x127 = x64*x99;
auto x128 = V{2}*x126 + V{2}*x127;
auto x129 = -V{1}*x128 + x239;
auto x130 = -V{1}*x129;
auto x131 = x128 + x239;
auto x132 = -V{1}*x131;
auto x133 = -V{1}*x121;
auto x134 = -V{1}*x118;
auto x135 = x1*(-x111*x229 - x119*x228 - x123*x230 - x124*x231 - x125*x232 + x130*x211 + x130*x217 + x132*x210 + x132*x216 + x133*x208 + x133*x214 + x134*x207 + x134*x212 + x206*x92 + x207*x87 + x208*x86 + x209*x85 + x210*x82 + x211*x81 + x212*x87 + x213*x91 + x214*x86 + x215*x90 + x216*x82 + x217*x81 + x218*x89 + x219*x88);
auto x136 = x135 - x182;
auto x137 = x118 + x135;
auto x138 = x137 - x183;
auto x139 = V{2}*x135;
auto x140 = x111 + x139 - x184;
auto x141 = x121 + x135;
auto x142 = x141 - x185;
auto x143 = x125 + x139 - x186;
auto x144 = -x187;
auto x145 = x131 + x135;
auto x146 = x145 - x188;
auto x147 = x129 + x135;
auto x148 = x147 - x189;
auto x149 = x124 + x139 - x190;
auto x150 = x123 + x139 - x191;
auto x151 = x137 - x192;
auto x152 = -x193;
auto x153 = x141 - x194;
auto x154 = -x195;
auto x155 = x119 + x139 - x196;
auto x156 = x145 - x197;
auto x157 = x147 - x198;
auto x158 = -x199;
auto x159 = -x200;
cell[0] = x136;
cell[1] = x138;
cell[2] = x140;
cell[3] = x142;
cell[4] = x143;
cell[5] = x144;
cell[6] = x146;
cell[7] = x148;
cell[8] = x149;
cell[9] = x150;
cell[10] = x151;
cell[11] = x152;
cell[12] = x153;
cell[13] = x154;
cell[14] = x155;
cell[15] = x156;
cell[16] = x157;
cell[17] = x158;
cell[18] = x159;
return { x136 + V{1}, ((x136)*(x136)) + ((x138)*(x138)) + ((x140)*(x140)) + ((x142)*(x142)) + ((x143)*(x143)) + ((x144)*(x144)) + ((x146)*(x146)) + ((x148)*(x148)) + ((x149)*(x149)) + ((x150)*(x150)) + ((x151)*(x151)) + ((x152)*(x152)) + ((x153)*(x153)) + ((x154)*(x154)) + ((x155)*(x155)) + ((x156)*(x156)) + ((x157)*(x157)) + ((x158)*(x158)) + ((x159)*(x159)) };
}
};

}

}

// Generation Info: commit=4ea71c45fdf920fe461e9926f71dc63f7fa1062d
