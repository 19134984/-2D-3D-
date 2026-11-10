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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::InnerEdgeDensity3D<0, 1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::InnerEdgeStress3D<0, 1, -1>, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x126 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x127 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x128 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x129 = parameters.template get<descriptors::OMEGA>();
auto x130 = ((x126)*(x126));
auto x131 = ((x127)*(x127));
auto x132 = ((x128)*(x128));
auto x133 = cell[3] + cell[6];
auto x135 = cell[13] + cell[14] + cell[4] + cell[5];
auto x136 = cell[10] + x135;
auto x137 = cell[11] + cell[2];
auto x0 = V{3};
auto x142 = ((x0)*(x0));
auto x1 = x130 + x131 + x132;
auto x2 = V{0.5}*x0*x1;
auto x143 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = cell[0] + cell[1];
auto x138 = V{2}*cell[18] + x5 + V{1};
auto x6 = x127 + V{1};
auto x139 = V{1} / (x6);
auto x7 = x133;
auto x8 = -V{1}*x128;
auto x9 = x8 + V{1};
auto x134 = V{1} / (x9);
auto x10 = V{0.5}*x134*(V{2}*cell[16] + V{2}*cell[8] + x136 + x137 + x138 + V{2}*x7) + V{0.5}*x139*(cell[10] + V{2}*cell[11] + cell[12] + V{2}*cell[13] + cell[15] + cell[16] + V{2}*cell[17] + V{2}*cell[5] + cell[7] + x133 + x138);
auto x11 = V{1}/x0;
auto x12 = -x11;
auto x13 = x10*x126;
auto x14 = x11*(x10 + V{-1});
auto x140 = -x14;
auto x141 = cell[15] + cell[16] + cell[7] + x140;
auto x15 = cell[1] + cell[6] - x126*x13 + x136 + x141;
auto x16 = x12*x15;
auto x17 = -V{1}*x127;
auto x18 = x128 + x17;
auto x19 = V{0.5}*x142;
auto x20 = x0*x18 + x143 + x19*((x18)*(x18));
auto x21 = V{0.0277777777777778};
auto x144 = -x21;
auto x22 = x10*x21;
auto x23 = x144 + x20*x22;
auto x24 = x127 + x8;
auto x25 = x0*x24 + x143 + x19*((x24)*(x24));
auto x26 = -x21 + x22*x25;
auto x27 = cell[18] + x23 - x26;
auto x145 = cell[17] + cell[18] + cell[8] + x27;
auto x28 = x10*x127;
auto x29 = -x127*x28 + x135 + x137 + x140 + x145;
auto x30 = x12*x29;
auto x31 = x16 + x30;
auto x32 = cell[12] - x10*x132 + x141 + x145 + x7;
auto x33 = x12*x32;
auto x34 = x4*(x10*x3 + x142*(V{0.5}*x31 + V{0.5}*x33) + V{-1});
auto x35 = V{1} - x129;
auto x36 = x126 + x8;
auto x37 = x0*x36 + x143 + x19*((x36)*(x36));
auto x38 = V{0.5}*x142*x21;
auto x39 = V{1} - x11;
auto x40 = x32*x39;
auto x41 = x15*x39;
auto x146 = x30 + x41;
auto x147 = x146 + x40;
auto x42 = V{-2};
auto x43 = cell[15] - V{1}*cell[16] + cell[6] - V{1}*cell[7] - x128*x13;
auto x44 = x38*(x147 + x42*x43);
auto x148 = x144 + x44;
auto x45 = x148 + x22*x37;
auto x46 = x126 + x128;
auto x47 = x0*x46 + x143 + x19*((x46)*(x46));
auto x48 = V{2};
auto x49 = x38*(x147 + x43*x48);
auto x149 = x144 + x49;
auto x50 = x149 + x22*x47;
auto x51 = V{0.0555555555555556};
auto x150 = -x51;
auto x52 = x10*x51;
auto x53 = x0*x128 + x132*x19 + x143;
auto x54 = V{0.5}*x142*x51;
auto x55 = x54*(x31 + x40);
auto x151 = x150 + x55;
auto x56 = x151 + x52*x53;
auto x57 = x0*x126 + x130*x19 + x143;
auto x58 = x54*(x146 + x33);
auto x152 = x150 + x58;
auto x59 = x152 + x52*x57;
auto x60 = -V{1}*x126;
auto x61 = x128 + x60;
auto x62 = x0*x61 + x143 + x19*((x61)*(x61));
auto x63 = x148 + x22*x62;
auto x64 = x60 + x8;
auto x65 = x0*x64 + x143 + x19*((x64)*(x64));
auto x66 = x149 + x22*x65;
auto x67 = x0*x8 + x143 + x19*((x8)*(x8));
auto x68 = x151 + x52*x67;
auto x69 = x0*x60 + x143 + x19*((x60)*(x60));
auto x70 = x152 + x52*x69;
auto x71 = x34 + x70;
auto x72 = x29*x39;
auto x155 = x33 + x41 + x72;
auto x73 = x16 + x72;
auto x153 = x40 + x73;
auto x74 = cell[17] - V{1}*cell[18] + cell[8] - x128*x28 - V{1}*x27;
auto x75 = x38*(x153 + x42*x74);
auto x76 = x26 + x75;
auto x159 = x59 + x71 + V{2}*x76 + V{1};
auto x77 = x127 + x128;
auto x78 = x0*x77 + x143 + x19*((x77)*(x77));
auto x79 = x38*(x153 + x48*x74);
auto x154 = x144 + x79;
auto x80 = x154 + x22*x78;
auto x81 = x126 + x127;
auto x82 = x0*x81 + x143 + x19*((x81)*(x81));
auto x83 = cell[13] - V{1}*cell[14] + cell[4] - V{1}*cell[5] - x127*x13;
auto x84 = x38*(x155 + x48*x83);
auto x156 = x144 + x84;
auto x85 = x156 + x22*x82;
auto x86 = x127 + x60;
auto x87 = x0*x86 + x143 + x19*((x86)*(x86));
auto x88 = x38*(x155 + x42*x83);
auto x157 = x144 + x88;
auto x89 = x157 + x22*x87;
auto x90 = x0*x127 + x131*x19 + x143;
auto x91 = x54*(x33 + x73);
auto x158 = x150 + x91;
auto x92 = x158 + x52*x90;
auto x93 = x126 + x17;
auto x94 = x0*x93 + x143 + x19*((x93)*(x93));
auto x95 = x157 + x22*x94;
auto x96 = x17 + x60;
auto x97 = x0*x96 + x143 + x19*((x96)*(x96));
auto x98 = x156 + x22*x97;
auto x99 = x0*x17 + x143 + x19*((x17)*(x17));
auto x100 = x158 + x52*x99;
auto x101 = x17 + x8;
auto x102 = x0*x101 + x143 + x19*((x101)*(x101));
auto x103 = x102*x22 + x154;
auto x104 = V{0.5}*x134*(x100 + V{2}*x103 + x159 + V{2}*x45 + V{2}*x66 + V{2}*x68 + x85 + x89 + x92 + x95 + x98) + V{0.5}*x139*(x159 + x45 + x50 + x56 + x63 + x66 + x68 + V{2}*x80 + V{2}*x85 + V{2}*x89 + V{2}*x92);
auto x105 = x104*x51;
auto x106 = x104*x21;
cell[0] = x129*x4*(x104*x3 + V{-1}) + x34*x35;
cell[10] = x129*(x105*x57 + x150) + x35*x59;
cell[11] = x129*(x105*x90 + x150) + x35*x92;
cell[12] = x129*(x105*x53 + x150) + x35*x56;
cell[13] = x129*(x106*x82 + x144) + x35*x85;
cell[14] = x129*(x106*x94 + x144) + x35*x95;
cell[15] = x129*(x106*x47 + x144) + x35*x50;
cell[16] = x129*(x106*x37 + x144) + x35*x45;
cell[17] = x129*(x106*x78 + x144) + x35*x80;
cell[18] = x129*(x106*x25 + x144) + x35*x76;
cell[1] = x129*(x105*x69 + x150) + x35*x70;
cell[2] = x100*x35 + x129*(x105*x99 + x150);
cell[3] = x129*(x105*x67 + x150) + x35*x68;
cell[4] = x129*(x106*x97 + x144) + x35*x98;
cell[5] = x129*(x106*x87 + x144) + x35*x89;
cell[6] = x129*(x106*x65 + x144) + x35*x66;
cell[7] = x129*(x106*x62 + x144) + x35*x63;
cell[8] = x103*x35 + x129*(x102*x106 + x144);
cell[9] = x129*(x106*x20 + x144) + x35*(x23 + x75);
return { x104, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
