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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, -1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<0, -1>, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x168 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x169 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x170 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x171 = parameters.template get<descriptors::OMEGA>();
auto x172 = ((x168)*(x168));
auto x173 = ((x169)*(x169));
auto x174 = ((x170)*(x170));
auto x176 = cell[12] + cell[3];
auto x177 = cell[11] + cell[2];
auto x0 = V{3};
auto x178 = ((x0)*(x0));
auto x1 = x172 + x173 + x174;
auto x2 = V{0.5}*x0*x1;
auto x179 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x186 = -x4;
auto x187 = x3*x4;
auto x5 = V{0.0277777777777778};
auto x180 = -x5;
auto x6 = -V{1}*x168;
auto x7 = x6 + V{1};
auto x175 = V{1} / (x7);
auto x8 = x175*(cell[0] + cell[17] + cell[18] + V{2}*cell[1] + V{2}*cell[4] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7] + cell[8] + cell[9] + x176 + x177 + V{1});
auto x9 = x5*x8;
auto x10 = -V{1}*x170;
auto x11 = x10 + x169;
auto x12 = V{0.5}*x178;
auto x13 = x0*x11 + x12*((x11)*(x11)) + x179;
auto x14 = x13*x9 + x180;
auto x15 = V{0.5}*x178*x5;
auto x16 = V{1}/x0;
auto x17 = V{1} - x16;
auto x18 = V{-2};
auto x19 = -V{1}*x18;
auto x20 = x10 + x6;
auto x21 = x0*x20 + x12*((x20)*(x20)) + x179;
auto x22 = x180 + x21*x9;
auto x23 = x19*(cell[6] - x22);
auto x24 = cell[18] - x14;
auto x25 = x169 + x170;
auto x26 = x0*x25 + x12*((x25)*(x25)) + x179;
auto x27 = x180 + x26*x9;
auto x28 = cell[17] - x27;
auto x29 = V{0.0555555555555556};
auto x30 = x29*x8;
auto x31 = x0*x170 + x12*x174 + x179;
auto x32 = -x29 + x30*x31;
auto x33 = -V{1}*x169;
auto x34 = x170 + x33;
auto x35 = x0*x34 + x12*((x34)*(x34)) + x179;
auto x36 = x180 + x35*x9;
auto x37 = cell[9] - x36;
auto x38 = x0*x10 + x12*((x10)*(x10)) + x179;
auto x39 = -x29 + x30*x38;
auto x40 = x10 + x33;
auto x41 = x0*x40 + x12*((x40)*(x40)) + x179;
auto x42 = x180 + x41*x9;
auto x43 = cell[8] - x42;
auto x181 = x24 + x28 + x37 + x43;
auto x44 = x170 + x6;
auto x45 = x0*x44 + x12*((x44)*(x44)) + x179;
auto x46 = x180 + x45*x9;
auto x47 = cell[7] - x46;
auto x48 = x176 + x181 + x23 - x32 - x39 + V{2}*x47;
auto x49 = x17*x48;
auto x50 = -x16;
auto x51 = x0*x6 + x12*((x6)*(x6)) + x179;
auto x52 = -x29 + x30*x51;
auto x53 = x33 + x6;
auto x54 = x0*x53 + x12*((x53)*(x53)) + x179;
auto x55 = -x5 + x54*x9;
auto x56 = x19*(cell[4] - x55);
auto x57 = x169 + x6;
auto x58 = x0*x57 + x12*((x57)*(x57)) + x179;
auto x59 = -x5 + x58*x9;
auto x60 = cell[5] - x59;
auto x61 = x19*x47 + x19*x60 + x19*(cell[1] - x52) + x23 + x56;
auto x62 = x50*x61;
auto x63 = x0*x169 + x12*x173 + x179;
auto x64 = -x29 + x30*x63;
auto x65 = x0*x33 + x12*((x33)*(x33)) + x179;
auto x66 = -x29 + x30*x65;
auto x67 = x177 + x181 + x56 + V{2}*x60 - x64 - x66;
auto x68 = x17*x67;
auto x69 = x62 + x68;
auto x182 = x49 + x69;
auto x70 = -V{1}*x24 + x28 - V{1}*x37 + x43;
auto x71 = x15*(x18*x70 + x182);
auto x72 = x14 + x71;
auto x73 = V{2};
auto x74 = x15*(x182 + x70*x73);
auto x75 = x27 + x74;
auto x76 = V{0.5}*x178*x29;
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
auto x189 = x80 + x86;
auto x87 = V{0.5}*x178*x4;
auto x88 = x66 + x82;
auto x188 = x83 + x88;
auto x89 = x17*x61;
auto x183 = x77 + x89;
auto x184 = x183 + x49;
auto x185 = x68 + x81 + x89;
auto x90 = x18*x47 + x23;
auto x91 = x15*(x18*x90 + x184) + x46;
auto x92 = x15*(x184 + x73*x90) + x22;
auto x93 = x18*x60 + x56;
auto x94 = x15*(x18*x93 + x185) + x59;
auto x95 = x52 + x76*(x183 + x81);
auto x96 = x15*(x185 + x73*x93) + x55;
auto x97 = x175*(x186 + x187*x8 + x188 + x189 + x72 + x75 + x84 + x85 + x87*(x78 + x81) + V{2}*x91 + V{2}*x92 + V{2}*x94 + V{2}*x95 + V{2}*x96 + V{1});
auto x98 = x29*x97;
auto x99 = -x29 + x51*x98;
auto x100 = x5*x97;
auto x101 = x100*x54 + x180;
auto x102 = -x19*(x101 - x96);
auto x103 = x100*x58 + x180;
auto x104 = -x103 + x94;
auto x105 = x100*x21 + x180;
auto x106 = -x19*(x105 - x92);
auto x107 = x100*x45 + x180;
auto x108 = -x107 + x91;
auto x109 = x102 + x104*x19 + x106 + x108*x19 + x19*(x95 - x99);
auto x110 = x109*x50;
auto x111 = x100*x13 + x180;
auto x112 = -x111 + x72;
auto x113 = x100*x26 + x180;
auto x114 = -x113 + x75;
auto x115 = -x29 + x63*x98;
auto x116 = x100*x35 + x180;
auto x117 = -x116 + x84;
auto x118 = -x29 + x65*x98;
auto x119 = x100*x41 + x180;
auto x120 = -x119 + x85;
auto x190 = x112 + x114 + x117 + x120;
auto x121 = x102 + V{2}*x104 - x115 - x118 + x188 + x190;
auto x122 = x121*x50;
auto x123 = x110 + x122;
auto x124 = -x29 + x31*x98;
auto x125 = -x29 + x38*x98;
auto x126 = x106 + V{2}*x108 - x124 - x125 + x189 + x190;
auto x127 = x126*x50;
auto x128 = V{1} - x171;
auto x192 = x128*x76;
auto x194 = x128*x15;
auto x129 = x109*x17;
auto x191 = x127 + x129;
auto x130 = x192*(x122 + x191);
auto x131 = x121*x17;
auto x193 = x131 + x191;
auto x132 = x110 + x131;
auto x133 = x192*(x127 + x132);
auto x134 = x126*x17;
auto x195 = x122 + x129 + x134;
auto x196 = x132 + x134;
auto x135 = x192*(x123 + x134);
auto x136 = x102 + x104*x18;
auto x137 = x194*(x136*x73 + x193);
auto x138 = x194*(x136*x18 + x193);
auto x139 = x106 + x108*x18;
auto x140 = x194*(x139*x73 + x195);
auto x141 = x194*(x139*x18 + x195);
auto x142 = -V{1}*x112 + x114 - V{1}*x117 + x120;
auto x143 = x194*(x142*x73 + x196);
auto x144 = x194*(x142*x18 + x196);
auto x145 = x168 + x169;
auto x146 = x168 + x33;
auto x147 = x168 + x170;
auto x148 = x10 + x168;
cell[0] = x128*x87*(x123 + x127) + x186 + x187*x97;
cell[10] = x130 - x29 + x98*(x0*x168 + x12*x172 + x179);
cell[11] = x115 + x133;
cell[12] = x124 + x135;
cell[13] = x100*(x0*x145 + x12*((x145)*(x145)) + x179) + x137 + x180;
cell[14] = x100*(x0*x146 + x12*((x146)*(x146)) + x179) + x138 + x180;
cell[15] = x100*(x0*x147 + x12*((x147)*(x147)) + x179) + x140 + x180;
cell[16] = x100*(x0*x148 + x12*((x148)*(x148)) + x179) + x141 + x180;
cell[17] = x113 + x143;
cell[18] = x111 + x144;
cell[1] = x130 + x99;
cell[2] = x118 + x133;
cell[3] = x125 + x135;
cell[4] = x101 + x137;
cell[5] = x103 + x138;
cell[6] = x105 + x140;
cell[7] = x107 + x141;
cell[8] = x119 + x143;
cell[9] = x116 + x144;
return { x97, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
