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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::TRT>, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x22 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x19 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x21 = parameters.template get<collision::TRT::MAGIC>();
auto x20 = parameters.template get<descriptors::OMEGA>();
auto x23 = cell[13] + cell[15];
auto x24 = cell[10] + x23;
auto x25 = cell[11] + cell[5];
auto x26 = cell[12] + cell[7] + cell[9];
auto x27 = cell[2] + cell[8];
auto x28 = cell[18] + cell[3];
auto x29 = cell[16] + cell[6];
auto x30 = cell[0] + cell[14] + cell[17] + cell[1] + cell[4] + x24 + x25 + x26 + x27 + x28 + x29 + V{1};
auto x31 = V{1} / (x30);
auto x32 = -cell[13];
auto x33 = cell[13] + cell[17];
auto x34 = -cell[2];
auto x35 = -cell[9];
auto x36 = -cell[8];
auto x37 = cell[18] + x36;
auto x38 = x34 + x35 + x37;
auto x39 = -cell[14];
auto x40 = -cell[4];
auto x41 = x25 + x39 + x40;
auto x42 = x33 + x38 + x41;
auto x43 = -cell[1];
auto x44 = -cell[7];
auto x45 = -cell[6];
auto x46 = cell[16] + x45;
auto x47 = x43 + x44 + x46;
auto x48 = cell[14] + x40;
auto x49 = -cell[5] + x48;
auto x50 = x24 + x47 + x49;
auto x51 = x31*x50;
auto x52 = ((x19)*(x19));
auto x53 = -cell[15];
auto x54 = cell[15] + cell[17];
auto x55 = -cell[18];
auto x56 = -cell[3];
auto x57 = x36 + x55 + x56;
auto x58 = -cell[16];
auto x59 = x45 + x58;
auto x60 = x26 + x54 + x57 + x59;
auto x61 = x31*x42;
auto x62 = -cell[17];
auto x63 = cell[9] + x62;
auto x64 = V{1}*x31;
auto x65 = ((x60)*(x60));
auto x66 = x52*x65;
auto x67 = V{0.333333333333333}*cell[0];
auto x68 = V{0.333333333333333}*cell[10] - V{0.666666666666667}*cell[17] - V{0.666666666666667}*cell[18] + V{0.333333333333333}*cell[1] - V{0.666666666666667}*cell[8] - V{0.666666666666667}*cell[9] + x67;
auto x69 = V{0.333333333333333}*cell[11] - V{0.666666666666667}*cell[15] - V{0.666666666666667}*cell[16] + V{0.333333333333333}*cell[2] - V{0.666666666666667}*cell[6] - V{0.666666666666667}*cell[7];
auto x70 = ((x42)*(x42));
auto x71 = x52*x70;
auto x72 = V{0.333333333333333}*cell[12] - V{0.666666666666667}*cell[13] - V{0.666666666666667}*cell[14] + V{0.333333333333333}*cell[3] - V{0.666666666666667}*cell[4] - V{0.666666666666667}*cell[5];
auto x73 = ((x50)*(x50));
auto x74 = x52*x73;
auto x75 = util::sqrt(V{25.455852}*x31*((x22)*(x22))*util::sqrt(V{2}*((x37 + x52*x60*x61 + x63)*(x37 + x52*x60*x61 + x63)) + V{2}*((cell[5] + x32 + x42*x51*x52 + x48)*(cell[5] + x32 + x42*x51*x52 + x48)) + V{2}*((cell[7] + x46 + x51*x52*x60 + x53)*(cell[7] + x46 + x51*x52*x60 + x53)) + ((-V{0.666666666666667}*cell[11] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] - V{0.666666666666667}*cell[2] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + x64*x71 + x68 + x72)*(-V{0.666666666666667}*cell[11] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] - V{0.666666666666667}*cell[2] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + x64*x71 + x68 + x72)) + ((-V{0.666666666666667}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] - V{0.666666666666667}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + x64*x66 + x68 + x69)*(-V{0.666666666666667}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] - V{0.666666666666667}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + x64*x66 + x68 + x69)) + ((-V{0.666666666666667}*cell[10] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] - V{0.666666666666667}*cell[1] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] + x64*x74 + x67 + x69 + x72)*(-V{0.666666666666667}*cell[10] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] - V{0.666666666666667}*cell[1] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] + x64*x74 + x67 + x69 + x72))) + V{1}/((x20)*(x20))) + V{1} / (x20);
auto x76 = V{1} / (x75);
auto x77 = V{1} / ((x30)*(x30));
auto x78 = V{1.5}*x77;
auto x79 = x74*x78;
auto x80 = x71*x78;
auto x81 = x66*x78;
auto x82 = x80 + x81 + V{-1};
auto x83 = x79 + x82;
auto x84 = V{1} / (V{2}*x21/(x75 + V{-1}) + V{0.5});
auto x85 = V{0.5}*cell[1];
auto x86 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x87 = V{3}*x19;
auto x88 = x51*x87;
auto x89 = -V{3}*x52*x73*x77 + x82 + x88;
auto x90 = V{3}*x77;
auto x91 = -x80;
auto x92 = V{1} - x81;
auto x93 = x91 + x92;
auto x94 = x88 + x93;
auto x95 = x74*x90 + x94;
auto x96 = x86*x95;
auto x97 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x98 = x76*(V{1}*cell[10] + V{1}*cell[1] + x89*x97 - x95*x97 + V{0.111111111111111});
auto x99 = V{0.5}*cell[2];
auto x100 = x61*x87;
auto x101 = x79 + V{-1};
auto x102 = x100 + x101 - V{3}*x52*x70*x77 + x81;
auto x103 = -x79;
auto x104 = x100 + x103;
auto x105 = x104 + x71*x90 + x92;
auto x106 = x105*x86;
auto x107 = x76*(V{1}*cell[11] + V{1}*cell[2] + x102*x97 - x105*x97 + V{0.111111111111111});
auto x108 = V{0.5}*cell[3];
auto x109 = x31*x60*x87;
auto x110 = x101 + x109 - V{3}*x52*x65*x77 + x80;
auto x111 = x103 + x109;
auto x112 = x111 + x66*x90 + x91 + V{1};
auto x113 = x112*x86;
auto x114 = x76*(V{1}*cell[12] + V{1}*cell[3] + x110*x97 - x112*x97 + V{0.111111111111111});
auto x115 = V{0.5}*cell[4];
auto x116 = V{0.0138888888888889}*cell[0] + V{0.0138888888888889}*cell[10] + V{0.0138888888888889}*cell[11] + V{0.0138888888888889}*cell[12] + V{0.0138888888888889}*cell[13] + V{0.0138888888888889}*cell[14] + V{0.0138888888888889}*cell[15] + V{0.0138888888888889}*cell[16] + V{0.0138888888888889}*cell[17] + V{0.0138888888888889}*cell[18] + V{0.0138888888888889}*cell[1] + V{0.0138888888888889}*cell[2] + V{0.0138888888888889}*cell[3] + V{0.0138888888888889}*cell[4] + V{0.0138888888888889}*cell[5] + V{0.0138888888888889}*cell[6] + V{0.0138888888888889}*cell[7] + V{0.0138888888888889}*cell[8] + V{0.0138888888888889}*cell[9] + V{0.0138888888888889};
auto x117 = cell[10] + x47;
auto x118 = ((cell[11] + 2*cell[13] - 2*cell[4] + x117 + x38 + x54)*(cell[11] + 2*cell[13] - 2*cell[4] + x117 + x38 + x54));
auto x119 = x83 + x88;
auto x120 = x100 - V{4.5}*x118*x52*x77 + x119;
auto x121 = V{4.5}*x77;
auto x122 = x104 + x118*x121*x52 + x94;
auto x123 = x116*x122;
auto x124 = x76*(V{1}*cell[13] + V{1}*cell[4] + x120*x86 - x122*x86 + V{0.0555555555555556});
auto x125 = -x100;
auto x126 = -cell[11];
auto x127 = V{2}*cell[14] + cell[15] - V{2}*cell[5] + x117 + x126 + x27 + x55 + x63;
auto x128 = x119 + x125 - V{4.5}*x52*x77*((x127)*(x127));
auto x129 = x83 - x88;
auto x130 = x100 + x129 - V{4.5}*x52*x77*((x127)*(x127));
auto x131 = V{0.5}*cell[14] - V{0.5}*cell[5];
auto x132 = x76*(V{1}*cell[14] + V{1}*cell[5] + x128*x86 + x130*x86 + V{0.0555555555555556});
auto x133 = V{0.5}*cell[6];
auto x134 = cell[10] + x43 + x49;
auto x135 = ((cell[12] + 2*cell[15] - 2*cell[6] + cell[9] + x134 + x33 + x57)*(cell[12] + 2*cell[15] - 2*cell[6] + cell[9] + x134 + x33 + x57));
auto x136 = x109 + x119 - V{4.5}*x135*x52*x77;
auto x137 = x111 + x121*x135*x52 + x94;
auto x138 = x116*x137;
auto x139 = x76*(V{1}*cell[15] + V{1}*cell[6] + x136*x86 - x137*x86 + V{0.0555555555555556});
auto x140 = -x109;
auto x141 = -cell[12];
auto x142 = cell[13] + x141;
auto x143 = V{2}*cell[16] - V{2}*cell[7] + cell[8] + x134 + x142 + x28 + x35 + x62;
auto x144 = x119 + x140 - V{4.5}*x52*x77*((x143)*(x143));
auto x145 = x109 + x129 - V{4.5}*x52*x77*((x143)*(x143));
auto x146 = V{0.5}*cell[16] - V{0.5}*cell[7];
auto x147 = x76*(V{1}*cell[16] + V{1}*cell[7] + x144*x86 + x145*x86 + V{0.0555555555555556});
auto x148 = V{0.5}*cell[8];
auto x149 = x34 + x41;
auto x150 = ((cell[12] + 2*cell[17] + cell[7] - 2*cell[8] + x149 + x23 + x56 + x59)*(cell[12] + 2*cell[17] + cell[7] - 2*cell[8] + x149 + x23 + x56 + x59));
auto x151 = x100 + x83;
auto x152 = x109 - V{4.5}*x150*x52*x77 + x151;
auto x153 = x104 + x109 + x121*x150*x52 + x93;
auto x154 = x116*x153;
auto x155 = x76*(V{1}*cell[17] + V{1}*cell[8] + x152*x86 - x153*x86 + V{0.0555555555555556});
auto x156 = V{2}*cell[18] + cell[3] - V{2}*cell[9] + x142 + x149 + x29 + x44 + x53;
auto x157 = x140 + x151 - V{4.5}*x52*x77*((x156)*(x156));
auto x158 = x109 + x125 - V{4.5}*x52*x77*((x156)*(x156)) + x83;
auto x159 = V{0.5}*cell[18] - V{0.5}*cell[9];
auto x160 = x76*(V{1}*cell[18] + V{1}*cell[9] + x157*x86 + x158*x86 + V{0.0555555555555556});
auto x0 = cell[0] - x76*(V{2}*cell[0] + x83*(V{0.666666666666667}*cell[0] + V{0.666666666666667}*cell[10] + V{0.666666666666667}*cell[11] + V{0.666666666666667}*cell[12] + V{0.666666666666667}*cell[13] + V{0.666666666666667}*cell[14] + V{0.666666666666667}*cell[15] + V{0.666666666666667}*cell[16] + V{0.666666666666667}*cell[17] + V{0.666666666666667}*cell[18] + V{0.666666666666667}*cell[1] + V{0.666666666666667}*cell[2] + V{0.666666666666667}*cell[3] + V{0.666666666666667}*cell[4] + V{0.666666666666667}*cell[5] + V{0.666666666666667}*cell[6] + V{0.666666666666667}*cell[7] + V{0.666666666666667}*cell[8] + V{0.666666666666667}*cell[9] + V{0.666666666666667}) + V{0.666666666666667});
auto x1 = cell[1] + x84*(V{0.5}*cell[10] - x85 - x86*x89 - x96) - x98;
auto x2 = cell[2] - x107 + x84*(V{0.5}*cell[11] - x102*x86 - x106 - x99);
auto x3 = cell[3] - x114 + x84*(V{0.5}*cell[12] - x108 - x110*x86 - x113);
auto x4 = cell[4] - x124 + x84*(V{0.5}*cell[13] - x115 - x116*x120 - x123);
auto x5 = cell[5] - x132 + x84*(-x116*x128 + x116*x130 + x131);
auto x6 = cell[6] - x139 + x84*(V{0.5}*cell[15] - x116*x136 - x133 - x138);
auto x7 = cell[7] - x147 + x84*(-x116*x144 + x116*x145 + x146);
auto x8 = cell[8] - x155 + x84*(V{0.5}*cell[17] - x116*x152 - x148 - x154);
auto x9 = cell[9] - x160 + x84*(-x116*x157 + x116*x158 + x159);
auto x10 = cell[10] - x84*(V{0.5}*cell[10] - x85 - x86*x89 - x96) - x98;
auto x11 = -x107 - x126 - x84*(V{0.5}*cell[11] - x102*x86 - x106 - x99);
auto x12 = -x114 - x141 - x84*(V{0.5}*cell[12] - x108 - x110*x86 - x113);
auto x13 = -x124 - x32 - x84*(V{0.5}*cell[13] - x115 - x116*x120 - x123);
auto x14 = -x132 - x39 - x84*(-x116*x128 + x116*x130 + x131);
auto x15 = -x139 - x53 - x84*(V{0.5}*cell[15] - x116*x136 - x133 - x138);
auto x16 = -x147 - x58 - x84*(-x116*x144 + x116*x145 + x146);
auto x17 = -x155 - x62 - x84*(V{0.5}*cell[17] - x116*x152 - x148 - x154);
auto x18 = -x160 - x55 - x84*(-x116*x157 + x116*x158 + x159);
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
cell[9] = x9;
cell[10] = x10;
cell[11] = x11;
cell[12] = x12;
cell[13] = x13;
cell[14] = x14;
cell[15] = x15;
cell[16] = x16;
cell[17] = x17;
cell[18] = x18;
return { x30, V{1}*x52*x77*(x65 + x70 + x73) };
}
};

}

}

// Generation Info: commit=cc8b060277f106540151cd2376323351f6e6960e
