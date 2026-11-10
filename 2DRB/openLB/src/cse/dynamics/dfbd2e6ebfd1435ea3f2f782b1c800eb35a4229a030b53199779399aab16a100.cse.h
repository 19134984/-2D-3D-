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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<0, -1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<0, -1, -1>, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x127 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x128 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x129 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x130 = parameters.template get<descriptors::OMEGA>();
auto x131 = ((x127)*(x127));
auto x132 = ((x128)*(x128));
auto x133 = ((x129)*(x129));
auto x134 = cell[2] + cell[4];
auto x135 = cell[3] + cell[6];
auto x139 = cell[15] + cell[16] + cell[7];
auto x140 = cell[12] + x139;
auto x0 = V{3};
auto x143 = ((x0)*(x0));
auto x1 = x131 + x132 + x133;
auto x2 = V{0.5}*x0*x1;
auto x144 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = cell[0] + cell[1];
auto x137 = cell[10] + V{2}*cell[8] + x5 + V{1};
auto x6 = x134;
auto x7 = -V{1}*x128;
auto x8 = x7 + V{1};
auto x138 = V{1} / (x8);
auto x9 = x135;
auto x10 = -V{1}*x129;
auto x11 = x10 + V{1};
auto x136 = V{1} / (x11);
auto x12 = V{0.5}*x136*(cell[11] + cell[13] + cell[14] + V{2}*cell[16] + V{2}*cell[18] + cell[5] + x134 + x137 + V{2}*x9) + V{0.5}*x138*(V{2}*cell[14] + V{2}*cell[9] + x135 + x137 + x140 + V{2}*x6);
auto x13 = V{1}/x0;
auto x14 = -x13;
auto x15 = x12*x127;
auto x16 = x13*(x12 + V{-1});
auto x141 = -x16;
auto x142 = cell[13] + cell[14] + cell[5] + x141;
auto x17 = cell[10] + cell[1] + cell[4] + cell[6] - x127*x15 + x139 + x142;
auto x18 = x14*x17;
auto x19 = x128 + x129;
auto x20 = V{0.5}*x143;
auto x21 = x0*x19 + x144 + x20*((x19)*(x19));
auto x22 = V{0.0277777777777778};
auto x145 = -x22;
auto x23 = x12*x22;
auto x24 = x145 + x21*x23;
auto x25 = x10 + x7;
auto x26 = x0*x25 + x144 + x20*((x25)*(x25));
auto x27 = -x22 + x23*x26;
auto x28 = cell[8] + x24 - x27;
auto x146 = cell[18] + cell[8] + cell[9] + x28;
auto x29 = x12*x128;
auto x30 = cell[11] - x128*x29 + x142 + x146 + x6;
auto x31 = x14*x30;
auto x32 = x18 + x31;
auto x33 = -x12*x133 + x140 + x141 + x146 + x9;
auto x34 = x14*x33;
auto x35 = x4*(x12*x3 + x143*(V{0.5}*x32 + V{0.5}*x34) + V{-1});
auto x36 = V{1} - x130;
auto x37 = x10 + x127;
auto x38 = x0*x37 + x144 + x20*((x37)*(x37));
auto x39 = V{0.5}*x143*x22;
auto x40 = V{1} - x13;
auto x41 = x33*x40;
auto x42 = x17*x40;
auto x147 = x31 + x42;
auto x148 = x147 + x41;
auto x43 = V{-2};
auto x44 = cell[15] - V{1}*cell[16] + cell[6] - V{1}*cell[7] - x129*x15;
auto x45 = x39*(x148 + x43*x44);
auto x149 = x145 + x45;
auto x46 = x149 + x23*x38;
auto x47 = x127 + x129;
auto x48 = x0*x47 + x144 + x20*((x47)*(x47));
auto x49 = V{2};
auto x50 = x39*(x148 + x44*x49);
auto x150 = x145 + x50;
auto x51 = x150 + x23*x48;
auto x52 = V{0.0555555555555556};
auto x151 = -x52;
auto x53 = x12*x52;
auto x54 = x0*x129 + x133*x20 + x144;
auto x55 = V{0.5}*x143*x52;
auto x56 = x55*(x32 + x41);
auto x152 = x151 + x56;
auto x57 = x152 + x53*x54;
auto x58 = x0*x127 + x131*x20 + x144;
auto x59 = x55*(x147 + x34);
auto x153 = x151 + x59;
auto x60 = x153 + x53*x58;
auto x61 = -V{1}*x127;
auto x62 = x129 + x61;
auto x63 = x0*x62 + x144 + x20*((x62)*(x62));
auto x64 = x149 + x23*x63;
auto x65 = x10 + x61;
auto x66 = x0*x65 + x144 + x20*((x65)*(x65));
auto x67 = x150 + x23*x66;
auto x68 = x0*x10 + x144 + x20*((x10)*(x10));
auto x69 = x152 + x53*x68;
auto x70 = x0*x61 + x144 + x20*((x61)*(x61));
auto x71 = x153 + x53*x70;
auto x72 = x35 + x71;
auto x73 = x127 + x7;
auto x74 = x0*x73 + x144 + x20*((x73)*(x73));
auto x75 = x30*x40;
auto x154 = x34 + x42 + x75;
auto x76 = cell[13] - V{1}*cell[14] + cell[4] - V{1}*cell[5] - x128*x15;
auto x77 = x39*(x154 + x43*x76);
auto x155 = x145 + x77;
auto x78 = x155 + x23*x74;
auto x79 = x129 + x7;
auto x80 = x0*x79 + x144 + x20*((x79)*(x79));
auto x81 = x18 + x75;
auto x156 = x41 + x81;
auto x82 = -V{1}*cell[18] + cell[8] - V{1}*cell[9] - x129*x29 + x28;
auto x83 = x39*(x156 + x43*x82);
auto x157 = x145 + x83;
auto x84 = x157 + x23*x80;
auto x85 = x39*(x156 + x49*x82);
auto x86 = x27 + x85;
auto x160 = x60 + x72 + V{2}*x86 + V{1};
auto x87 = x0*x7 + x144 + x20*((x7)*(x7));
auto x88 = x55*(x34 + x81);
auto x158 = x151 + x88;
auto x89 = x158 + x53*x87;
auto x90 = x61 + x7;
auto x91 = x0*x90 + x144 + x20*((x90)*(x90));
auto x92 = x39*(x154 + x49*x76);
auto x159 = x145 + x92;
auto x93 = x159 + x23*x91;
auto x94 = x127 + x128;
auto x95 = x0*x94 + x144 + x20*((x94)*(x94));
auto x96 = x159 + x23*x95;
auto x97 = x0*x128 + x132*x20 + x144;
auto x98 = x158 + x53*x97;
auto x99 = x128 + x61;
auto x100 = x0*x99 + x144 + x20*((x99)*(x99));
auto x101 = x100*x23 + x155;
auto x102 = x10 + x128;
auto x103 = x0*x102 + x144 + x20*((x102)*(x102));
auto x104 = x103*x23 + x157;
auto x105 = V{0.5}*x136*(x101 + V{2}*x104 + x160 + V{2}*x46 + V{2}*x67 + V{2}*x69 + x78 + x89 + x93 + x96 + x98) + V{0.5}*x138*(x160 + x46 + x51 + x57 + x64 + x67 + x69 + V{2}*x78 + V{2}*x84 + V{2}*x89 + V{2}*x93);
auto x106 = x105*x52;
auto x107 = x105*x22;
cell[0] = x130*x4*(x105*x3 + V{-1}) + x35*x36;
cell[10] = x130*(x106*x58 + x151) + x36*x60;
cell[11] = x130*(x106*x97 + x151) + x36*x98;
cell[12] = x130*(x106*x54 + x151) + x36*x57;
cell[13] = x130*(x107*x95 + x145) + x36*x96;
cell[14] = x130*(x107*x74 + x145) + x36*x78;
cell[15] = x130*(x107*x48 + x145) + x36*x51;
cell[16] = x130*(x107*x38 + x145) + x36*x46;
cell[17] = x130*(x107*x21 + x145) + x36*(x24 + x85);
cell[18] = x104*x36 + x130*(x103*x107 + x145);
cell[1] = x130*(x106*x70 + x151) + x36*x71;
cell[2] = x130*(x106*x87 + x151) + x36*x89;
cell[3] = x130*(x106*x68 + x151) + x36*x69;
cell[4] = x130*(x107*x91 + x145) + x36*x93;
cell[5] = x101*x36 + x130*(x100*x107 + x145);
cell[6] = x130*(x107*x66 + x145) + x36*x67;
cell[7] = x130*(x107*x63 + x145) + x36*x64;
cell[8] = x130*(x107*x26 + x145) + x36*x86;
cell[9] = x130*(x107*x80 + x145) + x36*x84;
return { x105, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
