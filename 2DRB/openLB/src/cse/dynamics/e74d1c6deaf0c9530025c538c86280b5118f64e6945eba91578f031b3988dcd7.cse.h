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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<1, -1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<1, -1, 1>, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x124 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x125 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x126 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x127 = parameters.template get<descriptors::OMEGA>();
auto x128 = ((x124)*(x124));
auto x129 = ((x125)*(x125));
auto x130 = ((x126)*(x126));
auto x132 = cell[13] + cell[14] + cell[4] + cell[5];
auto x133 = cell[10] + cell[1] + x132;
auto x134 = cell[11] + cell[2];
auto x135 = cell[0] + V{2}*cell[16] + x134 + V{1};
auto x137 = cell[17] + cell[18] + cell[8] + cell[9];
auto x138 = cell[12] + x137;
auto x0 = V{3};
auto x139 = ((x0)*(x0));
auto x1 = x128 + x129 + x130;
auto x2 = V{0.5}*x0*x1;
auto x140 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = cell[3] + cell[6];
auto x6 = -V{1}*x126;
auto x7 = x6 + V{1};
auto x131 = V{1} / (x7);
auto x8 = x124 + V{1};
auto x136 = V{1} / (x8);
auto x9 = V{0.5}*x131*(V{2}*cell[18] + V{2}*cell[8] + x133 + x135 + V{2}*x5) + V{0.5}*x136*(V{2}*cell[10] + V{2}*cell[13] + V{2}*cell[14] + V{2}*cell[15] + cell[3] + x135 + x138);
auto x10 = V{1}/x0;
auto x11 = -x10;
auto x12 = -V{1}*x124;
auto x13 = x12 + x126;
auto x14 = V{0.5}*x139;
auto x15 = x0*x13 + x14*((x13)*(x13)) + x140;
auto x16 = V{0.0277777777777778};
auto x141 = -x16;
auto x17 = x16*x9;
auto x18 = x141 + x15*x17;
auto x19 = x124 + x6;
auto x20 = x0*x19 + x14*((x19)*(x19)) + x140;
auto x21 = -x16 + x17*x20;
auto x22 = cell[16] + x18 - x21;
auto x23 = x124*x9;
auto x24 = x10*(x9 + V{-1});
auto x142 = -x24;
auto x143 = cell[15] + cell[16] + x142 + x22;
auto x25 = cell[6] - x124*x23 + x133 + x143;
auto x26 = x11*x25;
auto x27 = x125*x9;
auto x28 = -x125*x27 + x132 + x134 + x137 + x142;
auto x29 = x11*x28;
auto x30 = x26 + x29;
auto x31 = -x130*x9 + x138 + x143 + x5;
auto x32 = x11*x31;
auto x33 = x4*(x139*(V{0.5}*x30 + V{0.5}*x32) + x3*x9 + V{-1});
auto x34 = V{1} - x127;
auto x35 = -V{1}*x125;
auto x36 = x124 + x35;
auto x37 = x0*x36 + x14*((x36)*(x36)) + x140;
auto x38 = V{0.5}*x139*x16;
auto x39 = V{1} - x10;
auto x40 = x28*x39;
auto x41 = x25*x39;
auto x144 = x32 + x41;
auto x145 = x144 + x40;
auto x42 = V{-2};
auto x43 = cell[13] - V{1}*cell[14] + cell[4] - V{1}*cell[5] - x125*x23;
auto x44 = x38*(x145 + x42*x43);
auto x146 = x141 + x44;
auto x45 = x146 + x17*x37;
auto x46 = x124 + x125;
auto x47 = x0*x46 + x14*((x46)*(x46)) + x140;
auto x48 = V{2};
auto x49 = x38*(x145 + x43*x48);
auto x147 = x141 + x49;
auto x50 = x147 + x17*x47;
auto x51 = V{0.0555555555555556};
auto x148 = -x51;
auto x52 = x51*x9;
auto x53 = x0*x125 + x129*x14 + x140;
auto x54 = V{0.5}*x139*x51;
auto x55 = x26 + x40;
auto x56 = x54*(x32 + x55);
auto x149 = x148 + x56;
auto x57 = x149 + x52*x53;
auto x58 = x0*x124 + x128*x14 + x140;
auto x59 = x54*(x144 + x29);
auto x150 = x148 + x59;
auto x60 = x150 + x52*x58;
auto x61 = x12 + x125;
auto x62 = x0*x61 + x14*((x61)*(x61)) + x140;
auto x63 = x146 + x17*x62;
auto x64 = x12 + x35;
auto x65 = x0*x64 + x14*((x64)*(x64)) + x140;
auto x66 = x147 + x17*x65;
auto x67 = x0*x35 + x14*((x35)*(x35)) + x140;
auto x68 = x149 + x52*x67;
auto x69 = x0*x12 + x14*((x12)*(x12)) + x140;
auto x70 = x150 + x52*x69;
auto x71 = x125 + x6;
auto x72 = x0*x71 + x14*((x71)*(x71)) + x140;
auto x73 = x31*x39;
auto x151 = x55 + x73;
auto x153 = x29 + x41 + x73;
auto x74 = cell[17] - V{1}*cell[18] + cell[8] - V{1}*cell[9] - x126*x27;
auto x75 = x38*(x151 + x42*x74);
auto x152 = x141 + x75;
auto x76 = x152 + x17*x72;
auto x77 = cell[15] - V{1}*cell[16] + cell[6] - x126*x23 - V{1}*x22;
auto x78 = x38*(x153 + x42*x77);
auto x79 = x21 + x78;
auto x157 = x33 + x57 + x68 + V{2}*x79 + V{1};
auto x80 = x35 + x6;
auto x81 = x0*x80 + x14*((x80)*(x80)) + x140;
auto x82 = x38*(x151 + x48*x74);
auto x154 = x141 + x82;
auto x83 = x154 + x17*x81;
auto x84 = x0*x6 + x14*((x6)*(x6)) + x140;
auto x85 = x54*(x30 + x73);
auto x155 = x148 + x85;
auto x86 = x155 + x52*x84;
auto x87 = x12 + x6;
auto x88 = x0*x87 + x14*((x87)*(x87)) + x140;
auto x89 = x38*(x153 + x48*x77);
auto x156 = x141 + x89;
auto x90 = x156 + x17*x88;
auto x91 = x125 + x126;
auto x92 = x0*x91 + x14*((x91)*(x91)) + x140;
auto x93 = x154 + x17*x92;
auto x94 = x0*x126 + x130*x14 + x140;
auto x95 = x155 + x52*x94;
auto x96 = x126 + x35;
auto x97 = x0*x96 + x14*((x96)*(x96)) + x140;
auto x98 = x152 + x17*x97;
auto x99 = x124 + x126;
auto x100 = x0*x99 + x14*((x99)*(x99)) + x140;
auto x101 = x100*x17 + x156;
auto x102 = V{0.5}*x131*(x157 + x45 + x50 + x60 + x63 + x66 + x70 + V{2}*x76 + V{2}*x83 + V{2}*x86 + V{2}*x90) + V{0.5}*x136*(V{2}*x101 + x157 + V{2}*x45 + V{2}*x50 + V{2}*x60 + x76 + x83 + x86 + x93 + x95 + x98);
auto x103 = x102*x51;
auto x104 = x102*x16;
cell[0] = x127*x4*(x102*x3 + V{-1}) + x33*x34;
cell[10] = x127*(x103*x58 + x148) + x34*x60;
cell[11] = x127*(x103*x53 + x148) + x34*x57;
cell[12] = x127*(x103*x94 + x148) + x34*x95;
cell[13] = x127*(x104*x47 + x141) + x34*x50;
cell[14] = x127*(x104*x37 + x141) + x34*x45;
cell[15] = x101*x34 + x127*(x100*x104 + x141);
cell[16] = x127*(x104*x20 + x141) + x34*x79;
cell[17] = x127*(x104*x92 + x141) + x34*x93;
cell[18] = x127*(x104*x72 + x141) + x34*x76;
cell[1] = x127*(x103*x69 + x148) + x34*x70;
cell[2] = x127*(x103*x67 + x148) + x34*x68;
cell[3] = x127*(x103*x84 + x148) + x34*x86;
cell[4] = x127*(x104*x65 + x141) + x34*x66;
cell[5] = x127*(x104*x62 + x141) + x34*x63;
cell[6] = x127*(x104*x88 + x141) + x34*x90;
cell[7] = x127*(x104*x15 + x141) + x34*(x18 + x78);
cell[8] = x127*(x104*x81 + x141) + x34*x83;
cell[9] = x127*(x104*x97 + x141) + x34*x98;
return { x102, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
