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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkPressure<momenta::ForcedMomentum<momenta::CorrectedIncompressibleBulkMomentum> >, momenta::CorrectedIncompressibleBulkMomentum, momenta::IncompressibleBulkStress, momenta::DefineToNEq>, equilibria::CorrectedMPIncompressible, collision::OmegaFromCellTauEff<collision::IncompressibleSmagorinskyEffectiveOmega<collision::CIRLB> >, forcing::CorrectedLiang<momenta::Forced> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x153 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x154 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x155 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x156 = cell.template getFieldComponent<olb::descriptors::LAPLACERHO>(0);
auto x157 = cell.template getFieldComponent<olb::descriptors::NABLARHO>(0);
auto x158 = cell.template getFieldComponent<olb::descriptors::NABLARHO>(1);
auto x159 = cell.template getFieldComponent<olb::descriptors::NABLARHO>(2);
auto x160 = cell.template getFieldComponent<olb::descriptors::RHO>(0);
auto x161 = cell.template getFieldComponent<olb::descriptors::SCALAR>(0);
auto x163 = cell.template getFieldComponent<olb::descriptors::TAU_EFF>(0);
auto x164 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x165 = parameters.template get<collision::LES::SMAGORINSKY_F>();
auto x170 = cell[10] + cell[14] + cell[16];
auto x172 = cell[13] + cell[5];
auto x173 = cell[11] + cell[18] + x172;
auto x174 = cell[15] + cell[7];
auto x175 = cell[17] + cell[9];
auto x176 = cell[12] + x174 + x175;
auto x184 = cell[1] + cell[6] + x170;
auto x185 = cell[2] + cell[8] + x173;
auto x186 = cell[3] + x176;
auto x193 = V{1}*x157;
auto x197 = V{1}*x158;
auto x200 = V{1}*x159;
auto x0 = V{3};
auto x162 = V{1} / (x0);
auto x190 = V{0.5}*x0;
auto x1 = V{1}*x162;
auto x2 = V{0.333333333333333};
auto x166 = x0*x2;
auto x3 = -x1/(x2 + V{-1});
auto x4 = -V{0.5}*x166;
auto x183 = x160*x4;
auto x5 = -V{1}*cell[4];
auto x6 = -V{1}*cell[5];
auto x169 = cell[13] + x6;
auto x7 = -V{1}*cell[6];
auto x8 = -V{1}*cell[7];
auto x168 = cell[15] + x8;
auto x9 = V{0.25}*x156;
auto x10 = x160 + x9;
auto x167 = V{1} / (x10);
auto x192 = x10*x153;
auto x196 = x10*x154;
auto x199 = x10*x155;
auto x11 = V{0.5}*x153;
auto x12 = x11 + x167*(-V{1}*cell[1] + x168 + x169 + x170 + x5 + x7);
auto x177 = ((x12)*(x12));
auto x180 = x12*x157;
auto x13 = -V{1}*cell[8];
auto x14 = -V{1}*cell[9];
auto x171 = cell[17] + x14;
auto x15 = -V{1}*cell[14];
auto x16 = V{0.5}*x154;
auto x17 = x16 + x167*(-V{1}*cell[2] + x13 + x15 + x171 + x173 + x5);
auto x178 = ((x17)*(x17));
auto x181 = x158*x17;
auto x18 = -V{1}*cell[16];
auto x19 = -V{1}*cell[18];
auto x20 = V{0.5}*x155;
auto x21 = x167*(-V{1}*cell[3] + x13 + x176 + x18 + x19 + x7) + x20;
auto x179 = ((x21)*(x21));
auto x182 = x159*x21;
auto x22 = x177 + x178 + x179;
auto x23 = x180;
auto x24 = x181;
auto x25 = x182;
auto x26 = x23 + x24 + x25;
auto x27 = x3*(cell[4] + x183*x22 + x184 + x185 + x186 + V{0.5}*x26);
auto x187 = -x27;
auto x188 = cell[4] + x187;
auto x28 = x166;
auto x29 = V{0.5}*x22;
auto x195 = -x29;
auto x30 = -x1;
auto x31 = x12*x160;
auto x32 = -x12*x31 + x172 + x174 + x184 + x188;
auto x33 = V{0.5}*x162;
auto x189 = V{2}*x33;
auto x34 = x26*x33;
auto x35 = x32 + x34;
auto x36 = x30*x35;
auto x37 = x160*x17;
auto x38 = cell[14] - x17*x37 + x175 + x185 + x188;
auto x39 = x34 + x38;
auto x40 = x30*x39;
auto x41 = x36 + x40;
auto x42 = cell[16] + cell[18] + cell[6] + cell[8] - x160*x179 + x186 + x187;
auto x43 = x34 + x42;
auto x44 = x30*x43;
auto x45 = V{1}*x163;
auto x46 = cell[6] + x168 + x18 - x21*x31;
auto x47 = x33*(x12*x159 + x157*x21) + x46;
auto x48 = ((x47)*(x47));
auto x49 = cell[4] + x15 + x169 - x17*x31;
auto x50 = x33*(x12*x158 + x157*x17) + x49;
auto x51 = x180*x33 + x189*x23 + x32;
auto x52 = x182*x33 + x189*x25 + x42;
auto x53 = cell[8] + x171 + x19 - x21*x37;
auto x54 = x33*(x158*x21 + x159*x17) + x53;
auto x55 = x181*x33 + x189*x24 + x38;
auto x56 = x161*x164 - x165*(x161 + V{-1});
auto x57 = V{0.5}*(x45 + util::sqrt(((x45)*(x45)) + V{2.828428}*((x0)*(x0))*((x56)*(x56))*util::sqrt(V{2}*x48 + V{2}*((x50)*(x50)) + ((x51)*(x51)) + ((x52)*(x52)) + V{2}*((x54)*(x54)) + ((x55)*(x55)))/x160));
auto x58 = V{1}/x57;
auto x59 = V{1} - x58;
auto x60 = V{0.0555555555555556};
auto x61 = x0*x60;
auto x62 = -V{1}*x12;
auto x63 = x190;
auto x64 = x190*x61;
auto x65 = V{1} - x1;
auto x66 = x35*x65;
auto x191 = x44 + x66;
auto x67 = x64*(x191 + x40);
auto x68 = V{0.5}*x60;
auto x69 = -V{1}*x192;
auto x70 = V{1} - V{0.5}*x58;
auto x194 = x61*x70;
auto x71 = -x194*(x193*x62 - x69) + x59*(x67 - x68*(x0*x69 + x26)) + x61*(x160*(x195 + x62 + x63*((x62)*(x62))) + x27 + x62*x9);
auto x72 = -V{1}*x17;
auto x73 = x39*x65;
auto x202 = x191 + x73;
auto x74 = x36 + x73;
auto x75 = x64*(x44 + x74);
auto x198 = -x75;
auto x76 = -V{1}*x196;
auto x203 = x69 + x76;
auto x77 = -(x194*(x197*x72 - x76) + x59*(x198 + x68*(x0*x76 + x26)) - x61*(x160*(x195 + x63*((x72)*(x72)) + x72) + x27 + x72*x9));
auto x78 = -V{1}*x21;
auto x79 = x43*x65;
auto x207 = x40 + x66 + x79;
auto x211 = x74 + x79;
auto x80 = x64*(x41 + x79);
auto x201 = -x80;
auto x81 = -V{1}*x199;
auto x208 = x69 + x81;
auto x212 = x76 + x81;
auto x82 = -(x194*(x200*x78 - x81) + x59*(x201 + x68*(x0*x81 + x26)) - x61*(x160*(x195 + x63*((x78)*(x78)) + x78) + x27 + x78*x9));
auto x83 = V{0.0277777777777778};
auto x84 = x0*x83;
auto x204 = x70*x84;
auto x85 = x62 + x72;
auto x86 = x190*x84;
auto x87 = V{2};
auto x88 = x86*(x202 + x49*x87);
auto x89 = V{0.5}*x83;
auto x90 = -V{1}*x85;
auto x91 = x204*(x157*x90 + x158*x90 + x203) + x59*(x88 - x89*(x0*x203 + x26)) + x84*(x160*(x195 + x63*((x85)*(x85)) + x85) + x27 + x85*x9);
auto x92 = x17 + x62;
auto x93 = V{-2};
auto x94 = x86*(x202 + x49*x93);
auto x206 = -x94;
auto x95 = x196;
auto x205 = x69 + x95;
auto x219 = x81 + x95;
auto x96 = x204*(x158*x92 - x193*x92 + x205) - x59*(x206 + x89*(x0*x205 + x26)) + x84*(x160*(x195 + x63*((x92)*(x92)) + x92) + x27 + x9*x92);
auto x97 = x62 + x78;
auto x98 = x86*(x207 + x46*x87);
auto x209 = -x98;
auto x99 = -V{1}*x97;
auto x100 = x204*(x157*x99 + x159*x99 + x208) - x59*(x209 + x89*(x0*x208 + x26)) + x84*(x160*(x195 + x63*((x97)*(x97)) + x97) + x27 + x9*x97);
auto x101 = x21 + x62;
auto x102 = x86*(x207 + x46*x93);
auto x103 = x199;
auto x210 = x103 + x69;
auto x213 = x103 + x76;
auto x218 = x103 + x95;
auto x104 = x204*(x101*x159 - x101*x193 + x210) + x59*(x102 - x89*(x0*x210 + x26)) + x84*(x101*x9 + x160*(x101 + x195 + x63*((x101)*(x101))) + x27);
auto x105 = x72 + x78;
auto x106 = x86*(x211 + x53*x87);
auto x107 = -V{1}*x105;
auto x108 = x204*(x107*x158 + x107*x159 + x212) + x59*(x106 - x89*(x0*x212 + x26)) + x84*(x105*x9 + x160*(x105 + x195 + x63*((x105)*(x105))) + x27);
auto x109 = x21 + x72;
auto x110 = x86*(x211 + x53*x93);
auto x111 = x204*(x109*x159 - x109*x197 + x213) + x59*(x110 - x89*(x0*x213 + x26)) + x84*(x109*x9 + x160*(x109 + x195 + x63*((x109)*(x109))) + x27);
auto x112 = x192;
auto x214 = x112 + x95;
auto x215 = x112 + x76;
auto x216 = x103 + x112;
auto x217 = x112 + x81;
auto x113 = x194*(x112 + x23) + x59*(x67 - x68*(x0*x112 + x26)) + x61*(x12*x9 + x160*(x12 + x177*x63 + x195) + x27);
auto x114 = x194*(x24 + x95) - x59*(x198 + x68*(x0*x95 + x26)) + x61*(x160*(x17 + x178*x63 + x195) + x17*x9 + x27);
auto x115 = x194*(x103 + x25) - x59*(x201 + x68*(x0*x103 + x26)) + x61*(x160*(x179*x63 + x195 + x21) + x21*x9 + x27);
auto x225 = x104 + x111 + x115;
auto x116 = x12 + x17;
auto x117 = x204*(x116*x157 + x116*x158 + x214) + x59*(x88 - x89*(x0*x214 + x26)) + x84*(x116*x9 + x160*(x116 + x195 + x63*((x116)*(x116))) + x27);
auto x118 = x12 + x72;
auto x119 = x204*(x118*x157 - x118*x197 + x215) - x59*(x206 + x89*(x0*x215 + x26)) + x84*(x118*x9 + x160*(x118 + x195 + x63*((x118)*(x118))) + x27);
auto x120 = x12 + x21;
auto x121 = x204*(x120*x157 + x120*x159 + x216) - x59*(x209 + x89*(x0*x216 + x26)) + x84*(x120*x9 + x160*(x120 + x195 + x63*((x120)*(x120))) + x27);
auto x122 = x12 + x78;
auto x123 = x204*(x122*x157 - x122*x200 + x217) + x59*(x102 - x89*(x0*x217 + x26)) + x84*(x122*x9 + x160*(x122 + x195 + x63*((x122)*(x122))) + x27);
auto x222 = x113 + x119 + x123;
auto x124 = x17 + x21;
auto x125 = x204*(x124*x158 + x124*x159 + x218) + x59*(x106 - x89*(x0*x218 + x26)) + x84*(x124*x9 + x160*(x124 + x195 + x63*((x124)*(x124))) + x27);
auto x126 = x17 + x78;
auto x127 = x204*(x126*x158 - x126*x200 + x219) + x59*(x110 - x89*(x0*x219 + x26)) + x84*(x126*x9 + x160*(x126 + x195 + x63*((x126)*(x126))) + x27);
auto x224 = x114 + x127 + x96;
auto x128 = -V{1}*x91;
auto x220 = x117 + x128;
auto x129 = -V{1}*x100;
auto x221 = x121 + x129;
auto x130 = x11 + x167*(-V{1}*x104 + x220 + x221 + x222 - V{1}*x71 - V{1}*x96);
auto x131 = -V{1}*x108;
auto x223 = x125 + x131;
auto x132 = x16 + x167*(-V{1}*x111 - V{1}*x119 + x220 + x223 + x224 - V{1}*x77);
auto x133 = x167*(-V{1}*x123 - V{1}*x127 + x221 + x223 + x225 - V{1}*x82) + x20;
cell[0] = -V{1}*x0*x27 - x28*(x160*x29 + x187) + V{0.5}*x59*(x0*x28*(x41 + x44) - x2*x26);
cell[10] = x113;
cell[11] = x114;
cell[12] = x115;
cell[13] = x117;
cell[14] = x119;
cell[15] = x121;
cell[16] = x123;
cell[17] = x125;
cell[18] = x127;
cell[1] = x71;
cell[2] = x77;
cell[3] = x82;
cell[4] = x91;
cell[5] = x96;
cell[6] = x100;
cell[7] = x104;
cell[8] = x108;
cell[9] = x111;
cell.template getFieldPointer<olb::descriptors::TAU>()[0] = x57;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x3*(x100 + x108 + x117 + x121 + x125 + V{0.5}*x130*x157 + V{0.5}*x132*x158 + V{0.5}*x133*x159 + x183*(((x130)*(x130)) + ((x132)*(x132)) + ((x133)*(x133))) + x222 + x224 + x225 + x71 + x77 + x82 + x91);
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x130;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x132;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x133;
return { x160, x22 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
