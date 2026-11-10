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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkPressure<momenta::ForcedMomentum<momenta::CorrectedIncompressibleBulkMomentum> >, momenta::CorrectedIncompressibleBulkMomentum, momenta::IncompressibleBulkStress, momenta::DefineToNEq>, equilibria::CorrectedMPIncompressible, collision::OmegaFromCellTauEff<collision::IncompressibleSmagorinskyEffectiveOmega<collision::CIRLB> >, forcing::CorrectedLiang<momenta::Forced> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x124 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x125 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x126 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x127 = cell.template getFieldComponent<olb::descriptors::LAPLACERHO>(0);
auto x128 = cell.template getFieldComponent<olb::descriptors::NABLARHO>(0);
auto x129 = cell.template getFieldComponent<olb::descriptors::NABLARHO>(1);
auto x130 = cell.template getFieldComponent<olb::descriptors::NABLARHO>(2);
auto x131 = cell.template getFieldComponent<olb::descriptors::RHO>(0);
auto x132 = cell.template getFieldComponent<olb::descriptors::SCALAR>(0);
auto x134 = cell.template getFieldComponent<olb::descriptors::TAU_EFF>(0);
auto x135 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x136 = parameters.template get<collision::LES::SMAGORINSKY_F>();
auto x140 = cell[10] + cell[14] + cell[16];
auto x142 = cell[13] + cell[5];
auto x143 = cell[11] + cell[18] + x142;
auto x144 = cell[15] + cell[7];
auto x145 = cell[17] + cell[9];
auto x146 = cell[12] + x144 + x145;
auto x155 = cell[1] + cell[6] + x140;
auto x156 = cell[2] + cell[8] + x143;
auto x157 = cell[3] + x146;
auto x175 = V{1}*x129;
auto x180 = V{1}*x130;
auto x187 = V{1}*x128;
auto x0 = V{3};
auto x133 = V{1} / (x0);
auto x161 = V{0.5}*x0;
auto x1 = V{1}*x133;
auto x2 = V{0.333333333333333};
auto x153 = x0*x2;
auto x3 = -V{1}*cell[4];
auto x4 = -V{1}*cell[5];
auto x139 = cell[13] + x4;
auto x5 = -V{1}*cell[6];
auto x6 = -V{1}*cell[7];
auto x138 = cell[15] + x6;
auto x7 = V{0.25}*x127;
auto x8 = x131 + x7;
auto x137 = V{1} / (x8);
auto x163 = x124*x8;
auto x164 = x125*x8;
auto x165 = x126*x8;
auto x9 = V{0.5}*x124 + x137*(-V{1}*cell[1] + x138 + x139 + x140 + x3 + x5);
auto x149 = ((x9)*(x9));
auto x150 = x128*x9;
auto x10 = -V{1}*cell[8];
auto x11 = -V{1}*cell[9];
auto x141 = cell[17] + x11;
auto x12 = -V{1}*cell[14];
auto x13 = V{0.5}*x125 + x137*(-V{1}*cell[2] + x10 + x12 + x141 + x143 + x3);
auto x147 = ((x13)*(x13));
auto x151 = x129*x13;
auto x14 = -V{1}*cell[16];
auto x15 = -V{1}*cell[18];
auto x16 = V{0.5}*x126 + x137*(-V{1}*cell[3] + x10 + x14 + x146 + x15 + x5);
auto x148 = ((x16)*(x16));
auto x152 = x130*x16;
auto x17 = x147 + x148 + x149;
auto x154 = V{0.5}*x17;
auto x18 = x150;
auto x19 = x151;
auto x20 = x152;
auto x21 = x18 + x19 + x20;
auto x22 = -x1*(cell[4] - x131*x153*x154 + x155 + x156 + x157 + V{0.5}*x21)/(x2 + V{-1});
auto x158 = -x22;
auto x159 = cell[4] + x158;
auto x23 = x153;
auto x24 = x154;
auto x170 = -x24;
auto x25 = -x1;
auto x26 = x131*x9;
auto x27 = x142 + x144 + x155 + x159 - x26*x9;
auto x28 = V{0.5}*x133;
auto x160 = V{2}*x28;
auto x29 = x21*x28;
auto x30 = x27 + x29;
auto x31 = x25*x30;
auto x32 = x13*x131;
auto x33 = cell[14] - x13*x32 + x145 + x156 + x159;
auto x34 = x29 + x33;
auto x35 = x25*x34;
auto x36 = x31 + x35;
auto x37 = cell[16] + cell[18] + cell[6] + cell[8] - x131*x148 + x157 + x158;
auto x38 = x29 + x37;
auto x39 = x25*x38;
auto x40 = V{1}*x134;
auto x41 = cell[6] + x138 + x14 - x16*x26;
auto x42 = x28*(x128*x16 + x130*x9) + x41;
auto x43 = ((x42)*(x42));
auto x44 = cell[4] + x12 - x13*x26 + x139;
auto x45 = x28*(x128*x13 + x129*x9) + x44;
auto x46 = x150*x28 + x160*x18 + x27;
auto x47 = x152*x28 + x160*x20 + x37;
auto x48 = cell[8] + x141 + x15 - x16*x32;
auto x49 = x28*(x129*x16 + x13*x130) + x48;
auto x50 = x151*x28 + x160*x19 + x33;
auto x51 = x132*x135 - x136*(x132 + V{-1});
auto x52 = V{0.5}*(x40 + util::sqrt(((x40)*(x40)) + V{2.828428}*((x0)*(x0))*((x51)*(x51))*util::sqrt(V{2}*x43 + V{2}*((x45)*(x45)) + ((x46)*(x46)) + ((x47)*(x47)) + V{2}*((x49)*(x49)) + ((x50)*(x50)))/x131));
auto x53 = V{1}/x52;
auto x54 = V{1} - x53;
auto x55 = V{0.0555555555555556};
auto x56 = x0*x55;
auto x57 = -V{1}*x9;
auto x58 = x161;
auto x59 = x161*x56;
auto x60 = V{1} - x1;
auto x61 = x30*x60;
auto x162 = x39 + x61;
auto x62 = x59*(x162 + x35);
auto x63 = V{0.5}*x55;
auto x64 = -V{1}*x163;
auto x65 = V{1} - V{0.5}*x53;
auto x169 = x56*x65;
auto x66 = -V{1}*x13;
auto x67 = x34*x60;
auto x166 = x162 + x67;
auto x68 = x31 + x67;
auto x69 = x59*(x39 + x68);
auto x171 = -x69;
auto x70 = -V{1}*x164;
auto x188 = x64 + x70;
auto x71 = -V{1}*x16;
auto x72 = x38*x60;
auto x167 = x35 + x61 + x72;
auto x168 = x68 + x72;
auto x73 = x59*(x36 + x72);
auto x172 = -x73;
auto x74 = -V{1}*x165;
auto x190 = x64 + x74;
auto x192 = x70 + x74;
auto x75 = V{0.0277777777777778};
auto x76 = x0*x75;
auto x174 = x65*x76;
auto x77 = x57 + x66;
auto x78 = x161*x76;
auto x79 = V{2};
auto x80 = x78*(x166 + x44*x79);
auto x81 = V{0.5}*x75;
auto x82 = -V{1}*x77;
auto x83 = x13 + x57;
auto x84 = V{-2};
auto x85 = x78*(x166 + x44*x84);
auto x177 = -x85;
auto x86 = x164;
auto x185 = x74 + x86;
auto x189 = x64 + x86;
auto x87 = x57 + x71;
auto x88 = x78*(x167 + x41*x79);
auto x179 = -x88;
auto x89 = -V{1}*x87;
auto x90 = x16 + x57;
auto x91 = x78*(x167 + x41*x84);
auto x182 = -x91;
auto x92 = x165;
auto x183 = x86 + x92;
auto x191 = x64 + x92;
auto x193 = x70 + x92;
auto x93 = x66 + x71;
auto x94 = x78*(x168 + x48*x79);
auto x184 = -x94;
auto x95 = -V{1}*x93;
auto x96 = x16 + x66;
auto x97 = x78*(x168 + x48*x84);
auto x186 = -x97;
auto x98 = x163;
auto x173 = x86 + x98;
auto x176 = x70 + x98;
auto x178 = x92 + x98;
auto x181 = x74 + x98;
auto x99 = x13 + x9;
auto x100 = x66 + x9;
auto x101 = x16 + x9;
auto x102 = x71 + x9;
auto x103 = x13 + x16;
auto x104 = x13 + x71;
cell[0] = -V{1}*x0*x22 - x23*(x131*x24 + x158) + V{0.5}*x54*(x0*x23*(x36 + x39) - x2*x21);
cell[10] = x169*(x18 + x98) + x54*(x62 - x63*(x0*x98 + x21)) + x56*(x131*(x149*x58 + x170 + x9) + x22 + x7*x9);
cell[11] = x169*(x19 + x86) - x54*(x171 + x63*(x0*x86 + x21)) + x56*(x13*x7 + x131*(x13 + x147*x58 + x170) + x22);
cell[12] = x169*(x20 + x92) - x54*(x172 + x63*(x0*x92 + x21)) + x56*(x131*(x148*x58 + x16 + x170) + x16*x7 + x22);
cell[13] = x174*(x128*x99 + x129*x99 + x173) + x54*(x80 - x81*(x0*x173 + x21)) + x76*(x131*(x170 + x58*((x99)*(x99)) + x99) + x22 + x7*x99);
cell[14] = x174*(x100*x128 - x100*x175 + x176) - x54*(x177 + x81*(x0*x176 + x21)) + x76*(x100*x7 + x131*(x100 + x170 + x58*((x100)*(x100))) + x22);
cell[15] = x174*(x101*x128 + x101*x130 + x178) - x54*(x179 + x81*(x0*x178 + x21)) + x76*(x101*x7 + x131*(x101 + x170 + x58*((x101)*(x101))) + x22);
cell[16] = x174*(x102*x128 - x102*x180 + x181) - x54*(x182 + x81*(x0*x181 + x21)) + x76*(x102*x7 + x131*(x102 + x170 + x58*((x102)*(x102))) + x22);
cell[17] = x174*(x103*x129 + x103*x130 + x183) - x54*(x184 + x81*(x0*x183 + x21)) + x76*(x103*x7 + x131*(x103 + x170 + x58*((x103)*(x103))) + x22);
cell[18] = x174*(x104*x129 - x104*x180 + x185) - x54*(x186 + x81*(x0*x185 + x21)) + x76*(x104*x7 + x131*(x104 + x170 + x58*((x104)*(x104))) + x22);
cell[1] = -x169*(x187*x57 - x64) + x54*(x62 - x63*(x0*x64 + x21)) + x56*(x131*(x170 + x57 + x58*((x57)*(x57))) + x22 + x57*x7);
cell[2] = -(x169*(x175*x66 - x70) + x54*(x171 + x63*(x0*x70 + x21)) - x56*(x131*(x170 + x58*((x66)*(x66)) + x66) + x22 + x66*x7));
cell[3] = -(x169*(x180*x71 - x74) + x54*(x172 + x63*(x0*x74 + x21)) - x56*(x131*(x170 + x58*((x71)*(x71)) + x71) + x22 + x7*x71));
cell[4] = x174*(x128*x82 + x129*x82 + x188) + x54*(x80 - x81*(x0*x188 + x21)) + x76*(x131*(x170 + x58*((x77)*(x77)) + x77) + x22 + x7*x77);
cell[5] = x174*(x129*x83 - x187*x83 + x189) - x54*(x177 + x81*(x0*x189 + x21)) + x76*(x131*(x170 + x58*((x83)*(x83)) + x83) + x22 + x7*x83);
cell[6] = x174*(x128*x89 + x130*x89 + x190) - x54*(x179 + x81*(x0*x190 + x21)) + x76*(x131*(x170 + x58*((x87)*(x87)) + x87) + x22 + x7*x87);
cell[7] = x174*(x130*x90 - x187*x90 + x191) - x54*(x182 + x81*(x0*x191 + x21)) + x76*(x131*(x170 + x58*((x90)*(x90)) + x90) + x22 + x7*x90);
cell[8] = x174*(x129*x95 + x130*x95 + x192) - x54*(x184 + x81*(x0*x192 + x21)) + x76*(x131*(x170 + x58*((x93)*(x93)) + x93) + x22 + x7*x93);
cell[9] = x174*(x130*x96 - x175*x96 + x193) - x54*(x186 + x81*(x0*x193 + x21)) + x76*(x131*(x170 + x58*((x96)*(x96)) + x96) + x22 + x7*x96);
cell.template getFieldPointer<olb::descriptors::TAU>()[0] = x52;
return { x131, x17 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
