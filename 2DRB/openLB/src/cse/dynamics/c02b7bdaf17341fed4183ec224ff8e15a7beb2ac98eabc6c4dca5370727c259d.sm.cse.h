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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x113 = parameters.template get<descriptors::OMEGA>();
auto x114 = cell[10] + cell[14] + cell[16];
auto x115 = cell[1] + cell[6] + x114;
auto x116 = cell[13] + cell[5];
auto x117 = cell[11] + cell[18] + x116;
auto x118 = cell[2] + cell[8] + x117;
auto x119 = cell[15] + cell[7];
auto x120 = cell[17] + cell[9];
auto x121 = cell[12] + x119 + x120;
auto x122 = cell[3] + x121;
auto x0 = V{3};
auto x132 = ((x0)*(x0));
auto x1 = -V{1}*cell[4];
auto x2 = -V{1}*cell[5];
auto x124 = cell[13] + x2;
auto x3 = -V{1}*cell[6];
auto x4 = -V{1}*cell[7];
auto x125 = cell[15] + x4;
auto x5 = cell[0] + cell[4] + x115 + x118 + x122 + V{1};
auto x123 = V{1} / (x5);
auto x6 = x123*(-V{1}*cell[1] + x1 + x114 + x124 + x125 + x3);
auto x129 = ((x6)*(x6));
auto x7 = -V{1}*cell[8];
auto x8 = -V{1}*cell[9];
auto x126 = cell[17] + x8;
auto x9 = -V{1}*cell[14];
auto x10 = x123*(-V{1}*cell[2] + x1 + x117 + x126 + x7 + x9);
auto x127 = ((x10)*(x10));
auto x11 = -V{1}*cell[16];
auto x12 = -V{1}*cell[18];
auto x13 = x123*(-V{1}*cell[3] + x11 + x12 + x121 + x3 + x7);
auto x128 = ((x13)*(x13));
auto x14 = x127 + x128 + x129;
auto x15 = V{0.5}*x0*x14;
auto x135 = V{1} - x15;
auto x16 = V{0.333333333333333};
auto x17 = V{1}/x0;
auto x18 = -x17;
auto x19 = x5*x6;
auto x20 = x17*(x5 + V{-1});
auto x130 = -x20;
auto x131 = cell[4] + x130;
auto x21 = x115 + x116 + x119 + x131 - x19*x6;
auto x22 = x18*x21;
auto x23 = x10*x5;
auto x24 = cell[14] - x10*x23 + x118 + x120 + x131;
auto x25 = x18*x24;
auto x26 = x22 + x25;
auto x27 = cell[16] + cell[18] + cell[6] + cell[8] + x122 - x128*x5 + x130;
auto x28 = x18*x27;
auto x29 = V{1} - x113;
auto x30 = x16*(V{0.5}*x132*x29*(x26 + x28) - x5*(x15 + V{-1}) + V{-1});
auto x31 = -V{1}*x6;
auto x32 = V{0.5}*x132;
auto x33 = V{0.0555555555555556};
auto x136 = -x33;
auto x34 = x33*x5;
auto x35 = V{0.5}*x132*x33;
auto x134 = x29*x35;
auto x36 = V{1} - x17;
auto x37 = x21*x36;
auto x133 = x28 + x37;
auto x38 = x134*(x133 + x25);
auto x137 = x136 + x38;
auto x39 = x137 + x34*(x0*x31 + x135 + x32*((x31)*(x31)));
auto x40 = -V{1}*x10;
auto x41 = x24*x36;
auto x140 = x133 + x41;
auto x42 = x22 + x41;
auto x43 = x134*(x28 + x42);
auto x138 = x136 + x43;
auto x44 = x138 + x34*(x0*x40 + x135 + x32*((x40)*(x40)));
auto x45 = -V{1}*x13;
auto x46 = x27*x36;
auto x145 = x25 + x37 + x46;
auto x148 = x42 + x46;
auto x47 = x134*(x26 + x46);
auto x139 = x136 + x47;
auto x48 = x139 + x34*(x0*x45 + x135 + x32*((x45)*(x45)));
auto x49 = x31 + x40;
auto x50 = V{0.0277777777777778};
auto x142 = -x50;
auto x51 = x5*x50;
auto x52 = V{0.5}*x132*x50;
auto x141 = x29*x52;
auto x53 = V{2};
auto x54 = cell[4] - x10*x19 + x124 + x9;
auto x55 = x141*(x140 + x53*x54);
auto x143 = x142 + x55;
auto x56 = x143 + x51*(x0*x49 + x135 + x32*((x49)*(x49)));
auto x57 = x10 + x31;
auto x58 = V{-2};
auto x59 = x141*(x140 + x54*x58);
auto x144 = x142 + x59;
auto x60 = x144 + x51*(x0*x57 + x135 + x32*((x57)*(x57)));
auto x61 = x31 + x45;
auto x62 = cell[6] + x11 + x125 - x13*x19;
auto x63 = x141*(x145 + x53*x62);
auto x146 = x142 + x63;
auto x64 = x146 + x51*(x0*x61 + x135 + x32*((x61)*(x61)));
auto x65 = x13 + x31;
auto x66 = x141*(x145 + x58*x62);
auto x147 = x142 + x66;
auto x67 = x147 + x51*(x0*x65 + x135 + x32*((x65)*(x65)));
auto x68 = x40 + x45;
auto x69 = cell[8] + x12 + x126 - x13*x23;
auto x70 = x141*(x148 + x53*x69);
auto x149 = x142 + x70;
auto x71 = x149 + x51*(x0*x68 + x135 + x32*((x68)*(x68)));
auto x72 = x13 + x40;
auto x73 = x141*(x148 + x58*x69);
auto x150 = x142 + x73;
auto x74 = x150 + x51*(x0*x72 + x135 + x32*((x72)*(x72)));
auto x75 = x137 + x34*(x0*x6 + x129*x32 + x135);
auto x76 = x138 + x34*(x0*x10 + x127*x32 + x135);
auto x77 = x139 + x34*(x0*x13 + x128*x32 + x135);
auto x153 = x67 + x74 + x77;
auto x78 = x10 + x6;
auto x79 = x143 + x51*(x0*x78 + x135 + x32*((x78)*(x78)));
auto x80 = x40 + x6;
auto x81 = x144 + x51*(x0*x80 + x135 + x32*((x80)*(x80)));
auto x82 = x13 + x6;
auto x83 = x146 + x51*(x0*x82 + x135 + x32*((x82)*(x82)));
auto x84 = x45 + x6;
auto x85 = x147 + x51*(x0*x84 + x135 + x32*((x84)*(x84)));
auto x151 = x75 + x81 + x85;
auto x86 = x10 + x13;
auto x87 = x149 + x51*(x0*x86 + x135 + x32*((x86)*(x86)));
auto x88 = x10 + x45;
auto x89 = x150 + x51*(x0*x88 + x135 + x32*((x88)*(x88)));
auto x152 = x60 + x76 + x89;
auto x90 = x151 + x152 + x153 + x30 + x39 + x44 + x48 + x56 + x64 + x71 + x79 + x83 + x87 + V{1};
auto x154 = V{1} / (x90);
auto x91 = -V{1}*x56;
auto x155 = x79 + x91;
auto x92 = -V{1}*x64;
auto x156 = x83 + x92;
auto x93 = -V{1}*x71;
auto x157 = x87 + x93;
cell[0] = x30;
cell[10] = x75;
cell[11] = x76;
cell[12] = x77;
cell[13] = x79;
cell[14] = x81;
cell[15] = x83;
cell[16] = x85;
cell[17] = x87;
cell[18] = x89;
cell[1] = x39;
cell[2] = x44;
cell[3] = x48;
cell[4] = x56;
cell[5] = x60;
cell[6] = x64;
cell[7] = x67;
cell[8] = x71;
cell[9] = x74;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x90;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x154*(x151 + x155 + x156 - V{1}*x39 - V{1}*x60 - V{1}*x67);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x154*(x152 + x155 + x157 - V{1}*x44 - V{1}*x74 - V{1}*x81);
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x154*(x153 + x156 + x157 - V{1}*x48 - V{1}*x85 - V{1}*x89);
return { x5, x14 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
