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
struct CSE<dynamics::Tuple<T, descriptors::D3Q27<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::ParameterFromCell<collision::LES::SMAGORINSKY, collision::SmagorinskyEffectiveOmega<collision::BGK> >, dynamics::DefaultCombination>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x94 = cell.template getFieldComponent<olb::collision::LES::SMAGORINSKY>(0);
auto x96 = parameters.template get<descriptors::OMEGA>();
auto x95 = cell[12] + cell[25];
auto x97 = cell[6] + x95;
auto x98 = cell[11] + cell[24];
auto x99 = cell[4] + x98;
auto x100 = cell[19] + cell[20];
auto x101 = cell[14] + cell[18] + x100;
auto x102 = cell[1] + x101 + x99;
auto x103 = cell[23] + cell[26];
auto x104 = cell[10] + x103;
auto x105 = cell[8] + x104;
auto x106 = cell[16] + cell[7] + cell[9];
auto x107 = cell[3] + x105 + x106;
auto x108 = cell[13] + cell[17] + cell[5];
auto x109 = cell[21] + cell[22];
auto x110 = cell[15] + x108 + x109;
auto x111 = cell[2] + x110;
auto x0 = V{1}/x96;
auto x1 = -V{1}*cell[26];
auto x2 = -V{1}*cell[24];
auto x3 = -V{1}*cell[20];
auto x113 = cell[19] + x3;
auto x4 = -V{1}*cell[13];
auto x119 = cell[10] + cell[23] + x1 + x4;
auto x5 = -V{1}*cell[11];
auto x120 = x2 + x5;
auto x6 = -V{1}*cell[7];
auto x7 = -V{1}*cell[6];
auto x8 = -V{1}*cell[8];
auto x116 = cell[23] + x1 + x8;
auto x9 = -V{1}*cell[10];
auto x10 = -V{1}*cell[12];
auto x115 = cell[25] + x10 + x7 + x9;
auto x11 = -V{1}*cell[22];
auto x114 = cell[13] + cell[21] + x11;
auto x12 = cell[0] + x102 + x107 + x111 + x97 + V{1};
auto x112 = V{1} / (x12);
auto x13 = x112*(cell[11] - V{1}*cell[3] + x106 + x113 + x114 + x115 + x116 + x2);
auto x124 = ((x13)*(x13));
auto x14 = -V{1}*cell[4];
auto x118 = cell[24] + x14 + x5;
auto x15 = -V{1}*cell[5];
auto x117 = cell[17] + x15;
auto x16 = x112*(-V{1}*cell[1] + x101 + x103 + x115 + x117 + x118 + x4 + x6);
auto x126 = ((x16)*(x16));
auto x17 = x12*x16;
auto x18 = x113 + x119 + x120 - x13*x17 + x6 + x97;
auto x19 = ((x18)*(x18));
auto x20 = -V{1}*cell[25];
auto x121 = x10 + x20;
auto x21 = -V{1}*cell[18];
auto x22 = -V{1}*cell[9];
auto x23 = x112*(cell[12] - V{1}*cell[2] + x110 + x116 + x118 + x20 + x21 + x22 + x9);
auto x127 = ((x23)*(x23));
auto x24 = x117 + x119 + x121 - x17*x23 + x21 + x99;
auto x25 = V{3};
auto x125 = ((x25)*(x25));
auto x26 = V{1}*(x12 + V{-1})/x25;
auto x122 = -x26;
auto x123 = x122 + x97;
auto x27 = cell[7] + x102 + x104 + x108 + x123 - x16*x17;
auto x28 = cell[13] + x100 + x107 + x109 - x12*x124 + x123 + x98;
auto x29 = x12*x23;
auto x30 = x105 + x114 + x120 + x121 - x13*x29 + x22;
auto x31 = cell[18] + cell[9] + x105 + x111 + x122 - x23*x29 + x95 + x99;
auto x32 = V{2}*V{1} / (x0 + util::sqrt(V{2.828428}*x112*x125*((x94)*(x94))*util::sqrt(V{2}*x19 + V{2}*((x24)*(x24)) + ((x27)*(x27)) + ((x28)*(x28)) + V{2}*((x30)*(x30)) + ((x31)*(x31))) + ((x0)*(x0))));
auto x33 = V{1} - x32;
auto x34 = x124 + x126 + x127;
auto x35 = V{0.5}*x25*x34;
auto x128 = V{1} - x35;
auto x36 = V{0.296296296296296};
auto x37 = -V{1}*x16;
auto x38 = V{0.5}*x125;
auto x39 = V{0.0740740740740741};
auto x40 = x12*x39;
auto x41 = -V{1}*x23;
auto x42 = -V{1}*x13;
auto x43 = x37 + x41;
auto x44 = V{0.0185185185185185};
auto x45 = x12*x44;
auto x46 = x23 + x37;
auto x47 = x37 + x42;
auto x48 = x13 + x37;
auto x49 = x41 + x42;
auto x50 = x13 + x41;
auto x51 = x42 + x43;
auto x52 = V{0.00462962962962963};
auto x53 = x12*x52;
auto x54 = x13 + x43;
auto x55 = x42 + x46;
auto x56 = x13 + x46;
auto x57 = x16 + x23;
auto x58 = x16 + x41;
auto x59 = x13 + x16;
auto x60 = x16 + x42;
auto x61 = x13 + x23;
auto x62 = x23 + x42;
auto x63 = x13 + x57;
auto x64 = x42 + x57;
auto x65 = x13 + x58;
auto x66 = x42 + x58;
auto x67 = cell[0]*x33 - x32*x36*(x12*(x35 + V{-1}) + V{1});
auto x68 = cell[10]*x33 - x32*(x52 - x53*(x128 + x25*x51 + x38*((x51)*(x51))));
auto x69 = cell[11]*x33 - x32*(x52 - x53*(x128 + x25*x54 + x38*((x54)*(x54))));
auto x70 = cell[12]*x33 - x32*(x52 - x53*(x128 + x25*x55 + x38*((x55)*(x55))));
auto x71 = cell[13]*x33 - x32*(x52 - x53*(x128 + x25*x56 + x38*((x56)*(x56))));
auto x72 = cell[14]*x33 - x32*(x39 - x40*(x126*x38 + x128 + x16*x25));
auto x73 = cell[15]*x33 - x32*(x39 - x40*(x127*x38 + x128 + x23*x25));
auto x74 = cell[16]*x33 - x32*(x39 - x40*(x124*x38 + x128 + x13*x25));
auto x75 = cell[17]*x33 - x32*(x44 - x45*(x128 + x25*x57 + x38*((x57)*(x57))));
auto x76 = cell[18]*x33 - x32*(x44 - x45*(x128 + x25*x58 + x38*((x58)*(x58))));
auto x77 = cell[19]*x33 - x32*(x44 - x45*(x128 + x25*x59 + x38*((x59)*(x59))));
auto x78 = cell[1]*x33 - x32*(x39 - x40*(x128 + x25*x37 + x38*((x37)*(x37))));
auto x79 = cell[20]*x33 - x32*(x44 - x45*(x128 + x25*x60 + x38*((x60)*(x60))));
auto x80 = cell[21]*x33 - x32*(x44 - x45*(x128 + x25*x61 + x38*((x61)*(x61))));
auto x81 = cell[22]*x33 - x32*(x44 - x45*(x128 + x25*x62 + x38*((x62)*(x62))));
auto x82 = cell[23]*x33 - x32*(x52 - x53*(x128 + x25*x63 + x38*((x63)*(x63))));
auto x83 = cell[24]*x33 - x32*(x52 - x53*(x128 + x25*x64 + x38*((x64)*(x64))));
auto x84 = cell[25]*x33 - x32*(x52 - x53*(x128 + x25*x65 + x38*((x65)*(x65))));
auto x85 = cell[26]*x33 - x32*(x52 - x53*(x128 + x25*x66 + x38*((x66)*(x66))));
auto x86 = cell[2]*x33 - x32*(x39 - x40*(x128 + x25*x41 + x38*((x41)*(x41))));
auto x87 = cell[3]*x33 - x32*(x39 - x40*(x128 + x25*x42 + x38*((x42)*(x42))));
auto x88 = cell[4]*x33 - x32*(x44 - x45*(x128 + x25*x43 + x38*((x43)*(x43))));
auto x89 = cell[5]*x33 - x32*(x44 - x45*(x128 + x25*x46 + x38*((x46)*(x46))));
auto x90 = cell[6]*x33 - x32*(x44 - x45*(x128 + x25*x47 + x38*((x47)*(x47))));
auto x91 = cell[7]*x33 - x32*(x44 - x45*(x128 + x25*x48 + x38*((x48)*(x48))));
auto x92 = cell[8]*x33 - x32*(x44 - x45*(x128 + x25*x49 + x38*((x49)*(x49))));
auto x93 = cell[9]*x33 - x32*(x44 - x45*(x128 + x25*x50 + x38*((x50)*(x50))));
cell[0] = x67;
cell[10] = x68;
cell[11] = x69;
cell[12] = x70;
cell[13] = x71;
cell[14] = x72;
cell[15] = x73;
cell[16] = x74;
cell[17] = x75;
cell[18] = x76;
cell[19] = x77;
cell[1] = x78;
cell[20] = x79;
cell[21] = x80;
cell[22] = x81;
cell[23] = x82;
cell[24] = x83;
cell[25] = x84;
cell[26] = x85;
cell[2] = x86;
cell[3] = x87;
cell[4] = x88;
cell[5] = x89;
cell[6] = x90;
cell[7] = x91;
cell[8] = x92;
cell[9] = x93;
return { x12, x34 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
