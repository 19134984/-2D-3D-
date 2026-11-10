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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, -1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<0, -1>, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x182 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x183 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x184 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x185 = parameters.template get<descriptors::OMEGA>();
auto x186 = ((x182)*(x182));
auto x187 = ((x183)*(x183));
auto x188 = ((x184)*(x184));
auto x190 = cell[12] + cell[3];
auto x191 = cell[11] + cell[2];
auto x0 = V{3};
auto x192 = ((x0)*(x0));
auto x1 = x186 + x187 + x188;
auto x2 = V{0.5}*x0*x1;
auto x193 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x200 = -x4;
auto x201 = x3*x4;
auto x5 = V{0.0277777777777778};
auto x194 = -x5;
auto x6 = -V{1}*x182;
auto x7 = x6 + V{1};
auto x189 = V{1} / (x7);
auto x8 = x189*(cell[0] + cell[17] + cell[18] + V{2}*cell[1] + V{2}*cell[4] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7] + cell[8] + cell[9] + x190 + x191 + V{1});
auto x9 = x5*x8;
auto x10 = -V{1}*x184;
auto x11 = x10 + x183;
auto x12 = V{0.5}*x192;
auto x13 = x0*x11 + x12*((x11)*(x11)) + x193;
auto x14 = x13*x9 + x194;
auto x15 = V{0.5}*x192*x5;
auto x16 = V{1}/x0;
auto x17 = V{1} - x16;
auto x18 = V{-2};
auto x19 = -V{1}*x18;
auto x20 = x10 + x6;
auto x21 = x0*x20 + x12*((x20)*(x20)) + x193;
auto x22 = x194 + x21*x9;
auto x23 = x19*(cell[6] - x22);
auto x24 = cell[18] - x14;
auto x25 = x183 + x184;
auto x26 = x0*x25 + x12*((x25)*(x25)) + x193;
auto x27 = x194 + x26*x9;
auto x28 = cell[17] - x27;
auto x29 = V{0.0555555555555556};
auto x30 = x29*x8;
auto x31 = x0*x184 + x12*x188 + x193;
auto x32 = -x29 + x30*x31;
auto x33 = -V{1}*x183;
auto x34 = x184 + x33;
auto x35 = x0*x34 + x12*((x34)*(x34)) + x193;
auto x36 = x194 + x35*x9;
auto x37 = cell[9] - x36;
auto x38 = x0*x10 + x12*((x10)*(x10)) + x193;
auto x39 = -x29 + x30*x38;
auto x40 = x10 + x33;
auto x41 = x0*x40 + x12*((x40)*(x40)) + x193;
auto x42 = x194 + x41*x9;
auto x43 = cell[8] - x42;
auto x195 = x24 + x28 + x37 + x43;
auto x44 = x184 + x6;
auto x45 = x0*x44 + x12*((x44)*(x44)) + x193;
auto x46 = x194 + x45*x9;
auto x47 = cell[7] - x46;
auto x48 = x190 + x195 + x23 - x32 - x39 + V{2}*x47;
auto x49 = x17*x48;
auto x50 = -x16;
auto x51 = x0*x6 + x12*((x6)*(x6)) + x193;
auto x52 = -x29 + x30*x51;
auto x53 = x33 + x6;
auto x54 = x0*x53 + x12*((x53)*(x53)) + x193;
auto x55 = -x5 + x54*x9;
auto x56 = x19*(cell[4] - x55);
auto x57 = x183 + x6;
auto x58 = x0*x57 + x12*((x57)*(x57)) + x193;
auto x59 = -x5 + x58*x9;
auto x60 = cell[5] - x59;
auto x61 = x19*x47 + x19*x60 + x19*(cell[1] - x52) + x23 + x56;
auto x62 = x50*x61;
auto x63 = x0*x183 + x12*x187 + x193;
auto x64 = -x29 + x30*x63;
auto x65 = x0*x33 + x12*((x33)*(x33)) + x193;
auto x66 = -x29 + x30*x65;
auto x67 = x191 + x195 + x56 + V{2}*x60 - x64 - x66;
auto x68 = x17*x67;
auto x69 = x62 + x68;
auto x196 = x49 + x69;
auto x70 = -V{1}*x24 + x28 - V{1}*x37 + x43;
auto x71 = x15*(x18*x70 + x196);
auto x72 = x14 + x71;
auto x73 = V{2};
auto x74 = x15*(x196 + x70*x73);
auto x75 = x27 + x74;
auto x76 = V{0.5}*x192*x29;
auto x77 = x50*x67;
auto x78 = x62 + x77;
auto x79 = x76*(x49 + x78);
auto x80 = x32 + x79;
auto x81 = x48*x50;
auto x82 = x76*(x69 + x81);
auto x83 = x64 + x82;
auto x84 = x36 + x71;
auto x85 = x42 + x74;
auto x86 = x39 + x79;
auto x203 = x80 + x86;
auto x87 = V{0.5}*x192*x4;
auto x88 = x66 + x82;
auto x202 = x83 + x88;
auto x89 = x17*x61;
auto x197 = x77 + x89;
auto x198 = x197 + x49;
auto x199 = x68 + x81 + x89;
auto x90 = x18*x47 + x23;
auto x91 = x15*(x18*x90 + x198) + x46;
auto x92 = x15*(x198 + x73*x90) + x22;
auto x93 = x18*x60 + x56;
auto x94 = x15*(x18*x93 + x199) + x59;
auto x95 = x52 + x76*(x197 + x81);
auto x96 = x15*(x199 + x73*x93) + x55;
auto x97 = x189*(x200 + x201*x8 + x202 + x203 + x72 + x75 + x84 + x85 + x87*(x78 + x81) + V{2}*x91 + V{2}*x92 + V{2}*x94 + V{2}*x95 + V{2}*x96 + V{1});
auto x98 = x29*x97;
auto x99 = -x29 + x51*x98;
auto x100 = x5*x97;
auto x101 = x100*x54 + x194;
auto x102 = -x19*(x101 - x96);
auto x103 = x100*x58 + x194;
auto x104 = -x103 + x94;
auto x105 = x100*x21 + x194;
auto x106 = -x19*(x105 - x92);
auto x107 = x100*x45 + x194;
auto x108 = -x107 + x91;
auto x109 = x102 + x104*x19 + x106 + x108*x19 + x19*(x95 - x99);
auto x110 = x109*x50;
auto x111 = x100*x13 + x194;
auto x112 = -x111 + x72;
auto x113 = x100*x26 + x194;
auto x114 = -x113 + x75;
auto x115 = -x29 + x63*x98;
auto x116 = x100*x35 + x194;
auto x117 = -x116 + x84;
auto x118 = -x29 + x65*x98;
auto x119 = x100*x41 + x194;
auto x120 = -x119 + x85;
auto x204 = x112 + x114 + x117 + x120;
auto x121 = x102 + V{2}*x104 - x115 - x118 + x202 + x204;
auto x122 = x121*x50;
auto x123 = x110 + x122;
auto x124 = -x29 + x31*x98;
auto x125 = -x29 + x38*x98;
auto x126 = x106 + V{2}*x108 - x124 - x125 + x203 + x204;
auto x127 = x126*x50;
auto x128 = V{1} - x185;
auto x206 = x128*x76;
auto x208 = x128*x15;
auto x129 = x128*x87*(x123 + x127) + x200 + x201*x97;
auto x130 = x109*x17;
auto x205 = x127 + x130;
auto x131 = x206*(x122 + x205);
auto x132 = x131 + x99;
auto x133 = x121*x17;
auto x207 = x133 + x205;
auto x134 = x110 + x133;
auto x135 = x206*(x127 + x134);
auto x136 = x118 + x135;
auto x137 = x126*x17;
auto x209 = x122 + x130 + x137;
auto x210 = x134 + x137;
auto x138 = x206*(x123 + x137);
auto x139 = x125 + x138;
auto x140 = x102 + x104*x18;
auto x141 = x208*(x140*x73 + x207);
auto x142 = x101 + x141;
auto x143 = x208*(x140*x18 + x207);
auto x144 = x103 + x143;
auto x145 = x106 + x108*x18;
auto x146 = x208*(x145*x73 + x209);
auto x147 = x105 + x146;
auto x148 = x208*(x145*x18 + x209);
auto x149 = x107 + x148;
auto x150 = -V{1}*x112 + x114 - V{1}*x117 + x120;
auto x151 = x208*(x150*x73 + x210);
auto x152 = x119 + x151;
auto x153 = x208*(x150*x18 + x210);
auto x154 = x116 + x153;
auto x155 = x115 + x135;
auto x156 = x124 + x138;
auto x157 = x182 + x183;
auto x158 = x182 + x33;
auto x159 = x182 + x184;
auto x160 = x10 + x182;
auto x161 = x113 + x151;
auto x162 = x111 + x153;
cell[0] = x129;
cell[10] = x131 - x29 + x98*(x0*x182 + x12*x186 + x193);
cell[11] = x155;
cell[12] = x156;
cell[13] = x100*(x0*x157 + x12*((x157)*(x157)) + x193) + x141 + x194;
cell[14] = x100*(x0*x158 + x12*((x158)*(x158)) + x193) + x143 + x194;
cell[15] = x100*(x0*x159 + x12*((x159)*(x159)) + x193) + x146 + x194;
cell[16] = x100*(x0*x160 + x12*((x160)*(x160)) + x193) + x148 + x194;
cell[17] = x161;
cell[18] = x162;
cell[1] = x132;
cell[2] = x136;
cell[3] = x139;
cell[4] = x142;
cell[5] = x144;
cell[6] = x147;
cell[7] = x149;
cell[8] = x152;
cell[9] = x154;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x189*(x129 + V{2}*x132 + x136 + x139 + V{2}*x142 + V{2}*x144 + V{2}*x147 + V{2}*x149 + x152 + x154 + x155 + x156 + x161 + x162 + V{1});
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x182;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x183;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x184;
return { x97, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
