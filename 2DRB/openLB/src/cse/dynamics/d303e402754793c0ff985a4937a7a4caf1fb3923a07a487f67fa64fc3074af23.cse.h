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
struct CSE<Dual<ZouHeDynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<1, 1>, momenta::BulkStress, momenta::DefineSeparately>, 1, 1> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x149 = cell.template getFieldComponent<olb::opti::F>(1);
auto x132 = cell.template getFieldComponent<olb::opti::DJDF>(3);
auto x144 = cell.template getFieldComponent<olb::opti::DJDF>(15);
auto x163 = cell.template getFieldComponent<olb::opti::F>(15);
auto x164 = cell.template getFieldComponent<olb::opti::F>(16);
auto x165 = cell.template getFieldComponent<olb::opti::F>(17);
auto x136 = cell.template getFieldComponent<olb::opti::DJDF>(7);
auto x158 = cell.template getFieldComponent<olb::opti::F>(10);
auto x133 = cell.template getFieldComponent<olb::opti::DJDF>(4);
auto x145 = cell.template getFieldComponent<olb::opti::DJDF>(16);
auto x139 = cell.template getFieldComponent<olb::opti::DJDF>(10);
auto x151 = cell.template getFieldComponent<olb::opti::F>(3);
auto x141 = cell.template getFieldComponent<olb::opti::DJDF>(12);
auto x137 = cell.template getFieldComponent<olb::opti::DJDF>(8);
auto x166 = cell.template getFieldComponent<olb::opti::F>(18);
auto x143 = cell.template getFieldComponent<olb::opti::DJDF>(14);
auto x155 = cell.template getFieldComponent<olb::opti::F>(7);
auto x161 = cell.template getFieldComponent<olb::opti::F>(13);
auto x140 = cell.template getFieldComponent<olb::opti::DJDF>(11);
auto x160 = cell.template getFieldComponent<olb::opti::F>(12);
auto x135 = cell.template getFieldComponent<olb::opti::DJDF>(6);
auto x138 = cell.template getFieldComponent<olb::opti::DJDF>(9);
auto x131 = cell.template getFieldComponent<olb::opti::DJDF>(2);
auto x128 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<1, 1>::VELOCITY>(2);
auto x154 = cell.template getFieldComponent<olb::opti::F>(6);
auto x125 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x130 = cell.template getFieldComponent<olb::opti::DJDF>(1);
auto x147 = cell.template getFieldComponent<olb::opti::DJDF>(18);
auto x134 = cell.template getFieldComponent<olb::opti::DJDF>(5);
auto x148 = cell.template getFieldComponent<olb::opti::F>(0);
auto x142 = cell.template getFieldComponent<olb::opti::DJDF>(13);
auto x159 = cell.template getFieldComponent<olb::opti::F>(11);
auto x153 = cell.template getFieldComponent<olb::opti::F>(5);
auto x129 = cell.template getFieldComponent<olb::opti::DJDF>(0);
auto x167 = parameters.template get<descriptors::OMEGA>();
auto x146 = cell.template getFieldComponent<olb::opti::DJDF>(17);
auto x126 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<1, 1>::VELOCITY>(0);
auto x127 = V{1} / (x125);
auto x150 = x155 + x160;
auto x152 = x158 + x164;
auto x0 = V{1} - x167;
auto x1 = V{1}*x127;
auto x2 = x148 + x149;
auto x156 = x151 + x154 + x163 + x2 + V{1};
auto x3 = x127*(x150 + x152 + V{2}*x153 + x156 + V{2}*x159 + V{2}*x161 + V{2}*x165 + V{2}*x166) + V{-1};
auto x177 = V{1}*x3;
auto x4 = V{3};
auto x205 = V{1}*x4;
auto x5 = V{0.5}*x4;
auto x6 = V{0.0277777777777778};
auto x7 = x125*x6;
auto x157 = x167*x7;
auto x174 = V{1}*x157;
auto x8 = cell[18]*x157;
auto x9 = cell[17]*x157;
auto x10 = cell[14]*x157;
auto x168 = V{1}*x10;
auto x11 = cell[13]*x157;
auto x12 = V{0.0555555555555556};
auto x13 = x12*x125;
auto x162 = x13*x167;
auto x173 = V{1}*x162;
auto x14 = cell[11]*x162;
auto x15 = cell[9]*x157;
auto x169 = V{1}*x15;
auto x16 = cell[8]*x157;
auto x170 = V{1}*x16;
auto x17 = cell[5]*x157;
auto x18 = cell[4]*x157;
auto x171 = V{1}*x18;
auto x19 = cell[2]*x162;
auto x172 = V{1}*x19;
auto x20 = -x3*x5*(V{0.333333333333333}*cell[0]*x125*x167 + cell[10]*x173 + cell[12]*x173 + cell[15]*x174 + cell[16]*x174 + cell[1]*x173 + cell[3]*x173 + cell[6]*x174 + cell[7]*x174 + V{1}*x11 + V{1}*x14 + x168 + x169 + V{1}*x17 + x170 + x171 + x172 + V{1}*x8 + V{1}*x9);
auto x21 = V{0.5}*((x4)*(x4));
auto x180 = x21*x3;
auto x22 = -V{1}*x128;
auto x23 = x22 + x3;
auto x175 = x21*x23;
auto x24 = x175;
auto x25 = x128 + x3;
auto x176 = x21*x25;
auto x26 = x176;
auto x27 = -x177;
auto x185 = x21*x27;
auto x28 = x126 + x27;
auto x178 = x21*x28;
auto x29 = x178;
auto x30 = x126 + x3;
auto x179 = x21*x30;
auto x31 = x179;
auto x32 = x180;
auto x33 = x128 + x27;
auto x181 = x21*x33;
auto x34 = x181;
auto x35 = x22 + x27;
auto x182 = x21*x35;
auto x36 = x182;
auto x37 = -V{1}*x126;
auto x38 = x3 + x37;
auto x183 = x21*x38;
auto x39 = x183;
auto x40 = x27 + x37;
auto x184 = x21*x40;
auto x41 = x184;
auto x42 = x185;
auto x43 = x1*(x11*x179 + x11*x31 + x11*x4 + x14*x180 + x14*x32 + x14*x4 - x168*x178 - x168*x29 - x168*x4 - x169*x181 - x169*x34 - x169*x4 + x17*x183 + x17*x39 + x17*x4 - x170*x182 - x170*x36 - x170*x4 - x171*x184 - x171*x4 - x171*x41 - x172*x185 - x172*x4 - x172*x42 + x175*x8 + x176*x9 + V{2}*x20 + x24*x8 + x26*x9 + x4*x8 + x4*x9);
auto x44 = x5*(((x126)*(x126)) + ((x128)*(x128)) + ((x3)*(x3)));
auto x186 = V{1} - x44;
auto x45 = x166 - x7*(x186 + x23*x24 + x23*x4) + x7*(x186 + x33*x34 + x33*x4);
auto x188 = x150 + x165 + x45;
auto x46 = -V{1}*x154;
auto x187 = x163 + x46;
auto x47 = x165 - x7*(x186 + x25*x26 + x25*x4) + x7*(x186 + x35*x36 + x35*x4);
auto x48 = -V{1}*x151 - V{1}*x164 - V{1}*x166 + x187 + x188 - V{1}*x47;
auto x49 = x153 + x7*(x186 + x28*x29 + x28*x4) - x7*(x186 + x38*x39 + x38*x4);
auto x189 = x152 + x161 + x49;
auto x50 = x161 - x7*(x186 + x30*x31 + x30*x4) + x7*(x186 + x4*x40 + x40*x41);
auto x51 = x13*(x186 + x27*x4 + x27*x42) - x13*(x186 + x3*x32 + x3*x4) + x153 + x156 + V{2}*x159 + x166 + x188 + x189 + x47 + x50;
auto x191 = V{1} / (x51);
auto x192 = V{1}*x191;
auto x52 = V{0.25};
auto x190 = V{1}*x52;
auto x53 = cell[9]*x0;
auto x54 = cell[8]*x0;
auto x55 = -x190*(x53 - V{1}*x54);
auto x56 = -V{1}*x149 - V{1}*x153 - V{1}*x155 + x187 + x189 - V{1}*x50;
auto x57 = cell[14]*x0;
auto x58 = cell[4]*x0;
auto x59 = -x190*(x57 - V{1}*x58);
auto x60 = x51*x55;
auto x61 = x51*x59;
auto x62 = x191*(-x192*x48*x60 - x192*x56*x61 + x48*x55 + x56*x59);
auto x63 = x192;
auto x64 = x60*x63;
auto x193 = x62 + x64;
auto x65 = x193 + x53;
auto x194 = x65*x7;
auto x66 = -V{1}*x194;
auto x67 = -V{1}*x64;
auto x195 = x62 + x67;
auto x68 = x195 + x54;
auto x196 = x68*x7;
auto x69 = -V{1}*x196;
auto x70 = x61*x63;
auto x197 = x62 + x70;
auto x71 = x197 + x57;
auto x198 = x7*x71;
auto x72 = x198;
auto x202 = V{1}*x72;
auto x73 = -V{1}*x70;
auto x199 = x62 + x73;
auto x74 = x199 + x58;
auto x200 = x7*x74;
auto x75 = -V{1}*x200;
auto x76 = cell[2]*x0 + x62;
auto x201 = x13*x76;
auto x77 = -V{1}*x201;
auto x78 = x194;
auto x203 = V{1}*x78;
auto x79 = x196;
auto x204 = V{1}*x79;
auto x80 = -V{1}*x198;
auto x81 = x200;
auto x206 = V{1}*x81;
auto x82 = x201;
auto x207 = V{1}*x82;
auto x83 = -x177*x5*(x66 + x69 + x72 + x75 + x77 + x78 + x79 + x80 + x81 + x82);
auto x84 = x1*(x175*x66 + x176*x69 - x178*x202 + x179*x75 + x180*x77 - x181*x203 - x182*x204 + x183*x80 - x184*x206 - x185*x207 - x202*x29 - x202*x4 - x203*x34 - x203*x4 - x204*x36 - x204*x4 - x205*x81 - x205*x82 - x206*x41 - x207*x42 + x24*x66 + x26*x69 + x31*x75 + x32*x77 + x39*x80 + x4*x66 + x4*x69 + x4*x75 + x4*x77 + x4*x80 + V{2}*x83);
auto x208 = x43 + x84;
auto x209 = x195 + x208;
auto x211 = x193 + x208;
auto x85 = cell[0]*x0 - x129 + x208 + x62;
auto x86 = cell[1]*x0 - x130 + x199 + x208;
auto x87 = -x131;
auto x88 = cell[3]*x0 - x132 + x209;
auto x89 = -x133;
auto x90 = V{2}*x43;
auto x91 = V{2}*x84;
auto x210 = x90 + x91;
auto x92 = cell[5]*x0 - x134 + x199 + x210 + x71;
auto x93 = cell[6]*x0 - x135 + x209 + x73;
auto x94 = cell[7]*x0 - x136 + x211 + x73;
auto x95 = -x137;
auto x96 = -x138;
auto x97 = cell[10]*x0 - x139 + x197 + x208;
auto x98 = cell[11]*x0 - x140 + x210 + x62 + x76;
auto x99 = cell[12]*x0 - x141 + x211;
auto x100 = cell[13]*x0 - x142 + x197 + x210 + x74;
auto x101 = -x143;
auto x102 = cell[15]*x0 - x144 + x211 + x70;
auto x103 = cell[16]*x0 - x145 + x209 + x70;
auto x104 = cell[17]*x0 - x146 + x193 + x210 + x68;
auto x105 = cell[18]*x0 - x147 + x195 + x210 + x65;
cell[0] = x85;
cell[1] = x86;
cell[2] = x87;
cell[3] = x88;
cell[4] = x89;
cell[5] = x92;
cell[6] = x93;
cell[7] = x94;
cell[8] = x95;
cell[9] = x96;
cell[10] = x97;
cell[11] = x98;
cell[12] = x99;
cell[13] = x100;
cell[14] = x101;
cell[15] = x102;
cell[16] = x103;
cell[17] = x104;
cell[18] = x105;
return { x85 + V{1}, ((x100)*(x100)) + ((x101)*(x101)) + ((x102)*(x102)) + ((x103)*(x103)) + ((x104)*(x104)) + ((x105)*(x105)) + ((x85)*(x85)) + ((x86)*(x86)) + ((x87)*(x87)) + ((x88)*(x88)) + ((x89)*(x89)) + ((x92)*(x92)) + ((x93)*(x93)) + ((x94)*(x94)) + ((x95)*(x95)) + ((x96)*(x96)) + ((x97)*(x97)) + ((x98)*(x98)) + ((x99)*(x99)) };
}
};

}

}

// Generation Info: commit=4ea71c45fdf920fe461e9926f71dc63f7fa1062d
