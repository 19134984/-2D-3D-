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
struct CSE<Dual<ZouHeDynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, 1, -1> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x126 = cell.template getFieldComponent<olb::opti::DJDF>(14);
auto x119 = cell.template getFieldComponent<olb::opti::DJDF>(7);
auto x120 = cell.template getFieldComponent<olb::opti::DJDF>(8);
auto x127 = cell.template getFieldComponent<olb::opti::DJDF>(15);
auto x118 = cell.template getFieldComponent<olb::opti::DJDF>(6);
auto x139 = cell.template getFieldComponent<olb::opti::F>(8);
auto x131 = cell.template getFieldComponent<olb::opti::F>(0);
auto x109 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x124 = cell.template getFieldComponent<olb::opti::DJDF>(12);
auto x125 = cell.template getFieldComponent<olb::opti::DJDF>(13);
auto x146 = cell.template getFieldComponent<olb::opti::F>(15);
auto x138 = cell.template getFieldComponent<olb::opti::F>(7);
auto x143 = cell.template getFieldComponent<olb::opti::F>(12);
auto x150 = parameters.template get<descriptors::OMEGA>();
auto x133 = cell.template getFieldComponent<olb::opti::F>(2);
auto x123 = cell.template getFieldComponent<olb::opti::DJDF>(11);
auto x116 = cell.template getFieldComponent<olb::opti::DJDF>(4);
auto x145 = cell.template getFieldComponent<olb::opti::F>(14);
auto x111 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x130 = cell.template getFieldComponent<olb::opti::DJDF>(18);
auto x117 = cell.template getFieldComponent<olb::opti::DJDF>(5);
auto x128 = cell.template getFieldComponent<olb::opti::DJDF>(16);
auto x140 = cell.template getFieldComponent<olb::opti::F>(9);
auto x110 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x141 = cell.template getFieldComponent<olb::opti::F>(10);
auto x115 = cell.template getFieldComponent<olb::opti::DJDF>(3);
auto x135 = cell.template getFieldComponent<olb::opti::F>(4);
auto x112 = cell.template getFieldComponent<olb::opti::DJDF>(0);
auto x113 = cell.template getFieldComponent<olb::opti::DJDF>(1);
auto x114 = cell.template getFieldComponent<olb::opti::DJDF>(2);
auto x132 = cell.template getFieldComponent<olb::opti::F>(1);
auto x122 = cell.template getFieldComponent<olb::opti::DJDF>(10);
auto x137 = cell.template getFieldComponent<olb::opti::F>(6);
auto x134 = cell.template getFieldComponent<olb::opti::F>(3);
auto x121 = cell.template getFieldComponent<olb::opti::DJDF>(9);
auto x129 = cell.template getFieldComponent<olb::opti::DJDF>(17);
auto x147 = cell.template getFieldComponent<olb::opti::F>(16);
auto x142 = ((x109)*(x109));
auto x144 = ((x110)*(x110));
auto x148 = ((x111)*(x111));
auto x161 = x138 + x143;
auto x162 = x141 + x147;
auto x0 = V{1} - x150;
auto x1 = -V{1}*x110;
auto x2 = x1 + V{1};
auto x136 = V{1} / (x2);
auto x3 = V{1}*x136;
auto x4 = V{0.0277777777777778};
auto x5 = -V{1}*x111;
auto x6 = x110 + x5;
auto x7 = V{3};
auto x8 = V{0.5}*((x7)*(x7));
auto x9 = V{0.5}*x7*(x142 + x144 + x148);
auto x149 = V{1} - x9;
auto x10 = x149 + x6*x7 + x8*((x6)*(x6));
auto x155 = x10*x4;
auto x11 = x110 + x111;
auto x12 = x11*x7 + x149 + x8*((x11)*(x11));
auto x154 = x12*x4;
auto x13 = x109 + x5;
auto x14 = x109 + x111;
auto x15 = x1 + x109;
auto x16 = x149 + x15*x7 + x8*((x15)*(x15));
auto x153 = x16*x4;
auto x17 = x109 + x110;
auto x18 = x149 + x17*x7 + x8*((x17)*(x17));
auto x152 = x18*x4;
auto x19 = V{0.0555555555555556};
auto x20 = x110*x7 + x144*x8 + x149;
auto x151 = x19*x20;
auto x21 = x1 + x111;
auto x22 = x149 + x21*x7 + x8*((x21)*(x21));
auto x160 = x22*x4;
auto x23 = x1 + x5;
auto x24 = x149 + x23*x7 + x8*((x23)*(x23));
auto x159 = x24*x4;
auto x25 = -V{1}*x109;
auto x26 = x111 + x25;
auto x27 = x25 + x5;
auto x28 = x110 + x25;
auto x29 = x149 + x28*x7 + x8*((x28)*(x28));
auto x158 = x29*x4;
auto x30 = x1 + x25;
auto x31 = x149 + x30*x7 + x8*((x30)*(x30));
auto x157 = x31*x4;
auto x32 = x1*x7 + x149 + x8*((x1)*(x1));
auto x156 = x19*x32;
auto x33 = x150*x3*(-V{0.333333333333333}*cell[0]*(x9 + V{-1}) + cell[10]*x19*(x109*x7 + x142*x8 + x149) + cell[11]*x151 + cell[12]*x19*(x111*x7 + x148*x8 + x149) + cell[13]*x152 + cell[14]*x153 + cell[15]*x4*(x14*x7 + x149 + x8*((x14)*(x14))) + cell[16]*x4*(x13*x7 + x149 + x8*((x13)*(x13))) + cell[17]*x154 + cell[18]*x155 + cell[1]*x19*(x149 + x25*x7 + x8*((x25)*(x25))) + cell[2]*x156 + cell[3]*x19*(x149 + x5*x7 + x8*((x5)*(x5))) + cell[4]*x157 + cell[5]*x158 + cell[6]*x4*(x149 + x27*x7 + x8*((x27)*(x27))) + cell[7]*x4*(x149 + x26*x7 + x8*((x26)*(x26))) + cell[8]*x159 + cell[9]*x160);
auto x34 = x131 + x132;
auto x163 = x134 + x137 + x146 + x34 + V{1};
auto x35 = x136*(V{2}*x133 + V{2}*x135 + V{2}*x139 + V{2}*x140 + V{2}*x145 + x161 + x162 + x163);
auto x36 = x35*x4;
auto x37 = x12*x36 + x139 - x24*x36;
auto x165 = x140 + x161 + x37;
auto x38 = -V{1}*x137;
auto x164 = x146 + x38;
auto x39 = x10*x36 + x140 - x22*x36;
auto x40 = -V{1}*x134 - V{1}*x139 - V{1}*x147 + x164 + x165 - V{1}*x39;
auto x41 = x135 + x18*x36 - x31*x36;
auto x166 = x145 + x162 + x41;
auto x42 = x19*x35;
auto x43 = x145 - x16*x36 + x29*x36;
auto x44 = V{2}*x133 + x135 + x139 + x163 + x165 + x166 + x20*x42 - x32*x42 + x39 + x43;
auto x167 = V{1} / (x44);
auto x168 = V{1}*x167;
auto x45 = V{0.25};
auto x46 = cell[17]*x0;
auto x47 = cell[18]*x0;
auto x48 = x45*(x46 - V{1}*x47);
auto x49 = -V{1}*x48;
auto x50 = -V{1}*x132 - V{1}*x135 - V{1}*x138 + x164 + x166 - V{1}*x43;
auto x51 = cell[13]*x0;
auto x52 = cell[5]*x0;
auto x53 = x45*(x51 - V{1}*x52);
auto x54 = -V{1}*x53;
auto x55 = x44*x49;
auto x56 = x44*x54;
auto x57 = x167*(-x168*x40*x55 - x168*x50*x56 + x40*x49 + x50*x54);
auto x58 = x168;
auto x59 = x55*x58;
auto x170 = x57 + x59;
auto x60 = -V{1}*x59;
auto x169 = x57 + x60;
auto x61 = x169 + x47;
auto x62 = x170 + x46;
auto x63 = x56*x58;
auto x172 = x57 + x63;
auto x64 = -V{1}*x63;
auto x171 = x57 + x64;
auto x65 = x171 + x52;
auto x66 = x172 + x51;
auto x67 = cell[11]*x0 + x57;
auto x68 = x3*(x109*x53 + x111*x48 + x151*x67 + x152*x66 - V{1}*x153*x65 + x154*x62 + x155*x61 - V{1}*x156*x67 - V{1}*x157*x66 + x158*x65 - V{1}*x159*x62 - V{1}*x160*x61);
auto x173 = x33 + x68;
auto x175 = x169 + x173;
auto x176 = x170 + x173;
auto x69 = cell[0]*x0 - x112 + x173 + x57;
auto x70 = cell[1]*x0 - x113 + x171 + x173;
auto x71 = V{2}*x33;
auto x72 = V{2}*x68;
auto x174 = x71 + x72;
auto x73 = cell[2]*x0 - x114 + x174 + x57 + x67;
auto x74 = cell[3]*x0 - x115 + x175;
auto x75 = cell[4]*x0 - x116 + x171 + x174 + x66;
auto x76 = -x117;
auto x77 = cell[6]*x0 - x118 + x175 + x64;
auto x78 = cell[7]*x0 - x119 + x176 + x64;
auto x79 = cell[8]*x0 - x120 + x169 + x174 + x62;
auto x80 = cell[9]*x0 - x121 + x170 + x174 + x61;
auto x81 = cell[10]*x0 - x122 + x172 + x173;
auto x82 = -x123;
auto x83 = cell[12]*x0 - x124 + x176;
auto x84 = -x125;
auto x85 = cell[14]*x0 - x126 + x172 + x174 + x65;
auto x86 = cell[15]*x0 - x127 + x176 + x63;
auto x87 = cell[16]*x0 - x128 + x175 + x63;
auto x88 = -x129;
auto x89 = -x130;
cell[0] = x69;
cell[1] = x70;
cell[2] = x73;
cell[3] = x74;
cell[4] = x75;
cell[5] = x76;
cell[6] = x77;
cell[7] = x78;
cell[8] = x79;
cell[9] = x80;
cell[10] = x81;
cell[11] = x82;
cell[12] = x83;
cell[13] = x84;
cell[14] = x85;
cell[15] = x86;
cell[16] = x87;
cell[17] = x88;
cell[18] = x89;
return { x69 + V{1}, ((x69)*(x69)) + ((x70)*(x70)) + ((x73)*(x73)) + ((x74)*(x74)) + ((x75)*(x75)) + ((x76)*(x76)) + ((x77)*(x77)) + ((x78)*(x78)) + ((x79)*(x79)) + ((x80)*(x80)) + ((x81)*(x81)) + ((x82)*(x82)) + ((x83)*(x83)) + ((x84)*(x84)) + ((x85)*(x85)) + ((x86)*(x86)) + ((x87)*(x87)) + ((x88)*(x88)) + ((x89)*(x89)) };
}
};

}

}

// Generation Info: commit=4ea71c45fdf920fe461e9926f71dc63f7fa1062d
