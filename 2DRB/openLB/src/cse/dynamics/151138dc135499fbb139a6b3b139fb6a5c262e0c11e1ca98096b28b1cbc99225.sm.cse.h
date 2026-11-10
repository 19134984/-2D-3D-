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
struct CSE<dynamics::Tuple<T, descriptors::D3Q27<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::ParameterFromCell<collision::LES::SMAGORINSKY, collision::SmagorinskyEffectiveOmega<collision::BGK> >, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x129 = cell.template getFieldComponent<olb::collision::LES::SMAGORINSKY>(0);
auto x131 = parameters.template get<descriptors::OMEGA>();
auto x130 = cell[12] + cell[25];
auto x132 = cell[6] + x130;
auto x133 = cell[11] + cell[24];
auto x134 = cell[4] + x133;
auto x135 = cell[19] + cell[20];
auto x136 = cell[14] + cell[18] + x135;
auto x137 = cell[1] + x134 + x136;
auto x138 = cell[23] + cell[26];
auto x139 = cell[10] + x138;
auto x140 = cell[8] + x139;
auto x141 = cell[16] + cell[7] + cell[9];
auto x142 = cell[3] + x140 + x141;
auto x143 = cell[13] + cell[17] + cell[5];
auto x144 = cell[21] + cell[22];
auto x145 = cell[15] + x143 + x144;
auto x146 = cell[2] + x145;
auto x0 = V{1}/x131;
auto x1 = -V{1}*cell[26];
auto x2 = -V{1}*cell[24];
auto x3 = -V{1}*cell[20];
auto x148 = cell[19] + x3;
auto x4 = -V{1}*cell[13];
auto x154 = cell[10] + cell[23] + x1 + x4;
auto x5 = -V{1}*cell[11];
auto x155 = x2 + x5;
auto x6 = -V{1}*cell[7];
auto x7 = -V{1}*cell[6];
auto x8 = -V{1}*cell[8];
auto x151 = cell[23] + x1 + x8;
auto x9 = -V{1}*cell[10];
auto x10 = -V{1}*cell[12];
auto x150 = cell[25] + x10 + x7 + x9;
auto x11 = -V{1}*cell[22];
auto x149 = cell[13] + cell[21] + x11;
auto x12 = cell[0] + x132 + x137 + x142 + x146 + V{1};
auto x147 = V{1} / (x12);
auto x13 = x147*(cell[11] - V{1}*cell[3] + x141 + x148 + x149 + x150 + x151 + x2);
auto x159 = ((x13)*(x13));
auto x14 = -V{1}*cell[4];
auto x153 = cell[24] + x14 + x5;
auto x15 = -V{1}*cell[5];
auto x152 = cell[17] + x15;
auto x16 = x147*(-V{1}*cell[1] + x136 + x138 + x150 + x152 + x153 + x4 + x6);
auto x161 = ((x16)*(x16));
auto x17 = x12*x16;
auto x18 = -x13*x17 + x132 + x148 + x154 + x155 + x6;
auto x19 = ((x18)*(x18));
auto x20 = -V{1}*cell[25];
auto x156 = x10 + x20;
auto x21 = -V{1}*cell[18];
auto x22 = -V{1}*cell[9];
auto x23 = x147*(cell[12] - V{1}*cell[2] + x145 + x151 + x153 + x20 + x21 + x22 + x9);
auto x162 = ((x23)*(x23));
auto x24 = x134 + x152 + x154 + x156 - x17*x23 + x21;
auto x25 = V{3};
auto x160 = ((x25)*(x25));
auto x26 = V{1}*(x12 + V{-1})/x25;
auto x157 = -x26;
auto x158 = x132 + x157;
auto x27 = cell[7] + x137 + x139 + x143 + x158 - x16*x17;
auto x28 = cell[13] - x12*x159 + x133 + x135 + x142 + x144 + x158;
auto x29 = x12*x23;
auto x30 = -x13*x29 + x140 + x149 + x155 + x156 + x22;
auto x31 = cell[18] + cell[9] + x130 + x134 + x140 + x146 + x157 - x23*x29;
auto x32 = V{2}*V{1} / (x0 + util::sqrt(V{2.828428}*x147*x160*((x129)*(x129))*util::sqrt(V{2}*x19 + V{2}*((x24)*(x24)) + ((x27)*(x27)) + ((x28)*(x28)) + V{2}*((x30)*(x30)) + ((x31)*(x31))) + ((x0)*(x0))));
auto x33 = V{1} - x32;
auto x34 = x159 + x161 + x162;
auto x35 = V{0.5}*x25*x34;
auto x163 = V{1} - x35;
auto x36 = V{0.296296296296296};
auto x37 = cell[0]*x33 - x32*x36*(x12*(x35 + V{-1}) + V{1});
auto x38 = -V{1}*x16;
auto x39 = V{0.5}*x160;
auto x40 = V{0.0740740740740741};
auto x41 = x12*x40;
auto x42 = cell[1]*x33 - x32*(x40 - x41*(x163 + x25*x38 + x39*((x38)*(x38))));
auto x43 = -V{1}*x23;
auto x44 = cell[2]*x33 - x32*(x40 - x41*(x163 + x25*x43 + x39*((x43)*(x43))));
auto x45 = -V{1}*x13;
auto x46 = cell[3]*x33 - x32*(x40 - x41*(x163 + x25*x45 + x39*((x45)*(x45))));
auto x47 = x38 + x43;
auto x48 = V{0.0185185185185185};
auto x49 = x12*x48;
auto x50 = cell[4]*x33 - x32*(x48 - x49*(x163 + x25*x47 + x39*((x47)*(x47))));
auto x51 = x23 + x38;
auto x52 = cell[5]*x33 - x32*(x48 - x49*(x163 + x25*x51 + x39*((x51)*(x51))));
auto x53 = x38 + x45;
auto x54 = cell[6]*x33 - x32*(x48 - x49*(x163 + x25*x53 + x39*((x53)*(x53))));
auto x55 = x13 + x38;
auto x56 = cell[7]*x33 - x32*(x48 - x49*(x163 + x25*x55 + x39*((x55)*(x55))));
auto x57 = x43 + x45;
auto x58 = cell[8]*x33 - x32*(x48 - x49*(x163 + x25*x57 + x39*((x57)*(x57))));
auto x59 = x13 + x43;
auto x60 = cell[9]*x33 - x32*(x48 - x49*(x163 + x25*x59 + x39*((x59)*(x59))));
auto x61 = x45 + x47;
auto x62 = V{0.00462962962962963};
auto x63 = x12*x62;
auto x64 = cell[10]*x33 - x32*(x62 - x63*(x163 + x25*x61 + x39*((x61)*(x61))));
auto x65 = x13 + x47;
auto x66 = cell[11]*x33 - x32*(x62 - x63*(x163 + x25*x65 + x39*((x65)*(x65))));
auto x67 = x45 + x51;
auto x68 = cell[12]*x33 - x32*(x62 - x63*(x163 + x25*x67 + x39*((x67)*(x67))));
auto x69 = x13 + x51;
auto x70 = cell[13]*x33 - x32*(x62 - x63*(x163 + x25*x69 + x39*((x69)*(x69))));
auto x71 = cell[14]*x33 - x32*(x40 - x41*(x16*x25 + x161*x39 + x163));
auto x72 = cell[15]*x33 - x32*(x40 - x41*(x162*x39 + x163 + x23*x25));
auto x73 = cell[16]*x33 - x32*(x40 - x41*(x13*x25 + x159*x39 + x163));
auto x74 = x16 + x23;
auto x75 = cell[17]*x33 - x32*(x48 - x49*(x163 + x25*x74 + x39*((x74)*(x74))));
auto x76 = x16 + x43;
auto x77 = cell[18]*x33 - x32*(x48 - x49*(x163 + x25*x76 + x39*((x76)*(x76))));
auto x78 = x13 + x16;
auto x79 = cell[19]*x33 - x32*(x48 - x49*(x163 + x25*x78 + x39*((x78)*(x78))));
auto x80 = x16 + x45;
auto x81 = cell[20]*x33 - x32*(x48 - x49*(x163 + x25*x80 + x39*((x80)*(x80))));
auto x82 = x13 + x23;
auto x83 = cell[21]*x33 - x32*(x48 - x49*(x163 + x25*x82 + x39*((x82)*(x82))));
auto x84 = x23 + x45;
auto x85 = cell[22]*x33 - x32*(x48 - x49*(x163 + x25*x84 + x39*((x84)*(x84))));
auto x165 = x52 + x68 + x72 + x85;
auto x86 = x13 + x74;
auto x87 = cell[23]*x33 - x32*(x62 - x63*(x163 + x25*x86 + x39*((x86)*(x86))));
auto x166 = x56 + x60 + x66 + x73 + x87;
auto x88 = x45 + x74;
auto x89 = cell[24]*x33 - x32*(x62 - x63*(x163 + x25*x88 + x39*((x88)*(x88))));
auto x90 = x13 + x76;
auto x91 = cell[25]*x33 - x32*(x62 - x63*(x163 + x25*x90 + x39*((x90)*(x90))));
auto x92 = x45 + x76;
auto x93 = cell[26]*x33 - x32*(x62 - x63*(x163 + x25*x92 + x39*((x92)*(x92))));
auto x164 = x71 + x77 + x81 + x93;
auto x94 = x164 + x165 + x166 + x37 + x42 + x44 + x46 + x50 + x54 + x58 + x64 + x70 + x75 + x79 + x83 + x89 + x91 + V{1};
auto x167 = V{1} / (x94);
auto x95 = -V{1}*x50;
auto x96 = -V{1}*x54;
auto x97 = -V{1}*x64;
auto x98 = -V{1}*x66;
auto x168 = x75 + x87 + x89 + x95 + x97 + x98;
auto x99 = -V{1}*x68;
auto x169 = x79 + x91 + x96 + x99;
auto x100 = -V{1}*x58;
auto x101 = -V{1}*x93;
auto x170 = x100 + x101 + x70 + x83;
cell[0] = x37;
cell[10] = x64;
cell[11] = x66;
cell[12] = x68;
cell[13] = x70;
cell[14] = x71;
cell[15] = x72;
cell[16] = x73;
cell[17] = x75;
cell[18] = x77;
cell[19] = x79;
cell[1] = x42;
cell[20] = x81;
cell[21] = x83;
cell[22] = x85;
cell[23] = x87;
cell[24] = x89;
cell[25] = x91;
cell[26] = x93;
cell[2] = x44;
cell[3] = x46;
cell[4] = x50;
cell[5] = x52;
cell[6] = x54;
cell[7] = x56;
cell[8] = x58;
cell[9] = x60;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x94;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x167*(x164 + x168 + x169 - V{1}*x42 - V{1}*x52 - V{1}*x56 - V{1}*x70);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x167*(x165 + x168 + x170 - V{1}*x44 - V{1}*x60 - V{1}*x77 - V{1}*x91);
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x167*(x166 + x169 + x170 - V{1}*x46 - V{1}*x81 - V{1}*x85 - V{1}*x89 + x97);
return { x12, x34 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
