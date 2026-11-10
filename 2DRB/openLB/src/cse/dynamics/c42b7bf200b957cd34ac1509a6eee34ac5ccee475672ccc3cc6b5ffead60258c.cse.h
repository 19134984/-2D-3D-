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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::TrackAveragePopulation<collision::SmagorinskyEffectiveOmega<collision::BGK> >, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(0);
auto x40 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x41 = parameters.template get<descriptors::LATTICE_TIME_T>();
auto x23 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(4);
auto x26 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(7);
auto x32 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(13);
auto x21 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(2);
auto x28 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(9);
auto x34 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(15);
auto x38 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x30 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(11);
auto x20 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(1);
auto x33 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(14);
auto x24 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(5);
auto x25 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(6);
auto x27 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(8);
auto x39 = parameters.template get<descriptors::OMEGA>();
auto x37 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(18);
auto x22 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(3);
auto x36 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(17);
auto x29 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(10);
auto x35 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(16);
auto x31 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(12);
auto x42 = cell[15] + cell[17];
auto x43 = cell[12] + x42;
auto x44 = cell[11] + cell[18];
auto x45 = cell[10] + cell[14] + cell[16];
auto x46 = cell[2] + cell[8] + cell[9];
auto x47 = cell[13] + cell[3];
auto x48 = cell[0] + cell[1] + cell[4] + cell[5] + cell[6] + cell[7] + x43 + x44 + x45 + x46 + x47 + V{1};
auto x49 = V{1} / (x48);
auto x50 = V{0.333333333333333}*cell[13];
auto x51 = V{0.333333333333333}*cell[14];
auto x52 = V{0.333333333333333}*cell[4];
auto x53 = V{0.333333333333333}*cell[5];
auto x54 = V{1}*x49;
auto x55 = ((x38)*(x38));
auto x56 = -cell[18];
auto x57 = -cell[3];
auto x58 = -cell[8];
auto x59 = cell[9] + x58;
auto x60 = x56 + x57 + x59;
auto x61 = -cell[6];
auto x62 = cell[7] + x61;
auto x63 = -cell[16] + x62;
auto x64 = x43 + x60 + x63;
auto x65 = ((x64)*(x64));
auto x66 = x55*x65;
auto x67 = V{0.333333333333333}*cell[0];
auto x68 = V{0.333333333333333}*cell[10];
auto x69 = V{0.333333333333333}*cell[1];
auto x70 = -V{0.666666666666667}*cell[17] - V{0.666666666666667}*cell[18] - V{0.666666666666667}*cell[8] - V{0.666666666666667}*cell[9] + x67 + x68 + x69;
auto x71 = V{0.333333333333333}*cell[11];
auto x72 = V{0.333333333333333}*cell[2];
auto x73 = -V{0.666666666666667}*cell[15] - V{0.666666666666667}*cell[16] - V{0.666666666666667}*cell[6] - V{0.666666666666667}*cell[7] + x71 + x72;
auto x74 = V{0.333333333333333}*cell[15];
auto x75 = V{0.333333333333333}*cell[16];
auto x76 = V{0.333333333333333}*cell[6];
auto x77 = V{0.333333333333333}*cell[7];
auto x78 = cell[13] + cell[17];
auto x79 = -cell[2];
auto x80 = -cell[9];
auto x81 = x44 + x58 + x79 + x80;
auto x82 = -cell[4];
auto x83 = cell[5] + x82;
auto x84 = -cell[14] + x83;
auto x85 = x78 + x81 + x84;
auto x86 = ((x85)*(x85));
auto x87 = x55*x86;
auto x88 = V{0.333333333333333}*cell[12];
auto x89 = V{0.333333333333333}*cell[3];
auto x90 = -V{0.666666666666667}*cell[13] - V{0.666666666666667}*cell[14] - V{0.666666666666667}*cell[4] - V{0.666666666666667}*cell[5] + x88 + x89;
auto x91 = V{0.333333333333333}*cell[17];
auto x92 = V{0.333333333333333}*cell[18];
auto x93 = V{0.333333333333333}*cell[8];
auto x94 = V{0.333333333333333}*cell[9];
auto x95 = cell[13] + cell[15];
auto x96 = -cell[1];
auto x97 = -cell[7];
auto x98 = x61 + x96 + x97;
auto x99 = -cell[5] + x82;
auto x100 = x45 + x95 + x98 + x99;
auto x101 = ((x100)*(x100));
auto x102 = x101*x55;
auto x103 = x100*x49;
auto x104 = -cell[15] + cell[16];
auto x105 = x49*x85;
auto x106 = -cell[17];
auto x107 = cell[18] + x106;
auto x108 = V{1} / (V{3.00000046417339}*util::sqrt(x49*((x40)*(x40))*util::sqrt(((x103*x55*x64 + x104 + x62)*(x103*x55*x64 + x104 + x62)) + ((x105*x55*x64 + x107 + x59)*(x105*x55*x64 + x107 + x59)) + ((-cell[13] + cell[14] + x103*x55*x85 + x83)*(-cell[13] + cell[14] + x103*x55*x85 + x83)) + V{0.5}*((-V{0.666666666666667}*cell[11] - V{0.666666666666667}*cell[2] + x54*x87 + x70 + x74 + x75 + x76 + x77 + x90)*(-V{0.666666666666667}*cell[11] - V{0.666666666666667}*cell[2] + x54*x87 + x70 + x74 + x75 + x76 + x77 + x90)) + V{0.5}*((-V{0.666666666666667}*cell[12] - V{0.666666666666667}*cell[3] + x50 + x51 + x52 + x53 + x54*x66 + x70 + x73)*(-V{0.666666666666667}*cell[12] - V{0.666666666666667}*cell[3] + x50 + x51 + x52 + x53 + x54*x66 + x70 + x73)) + V{0.5}*((-V{0.666666666666667}*cell[10] - V{0.666666666666667}*cell[1] + x102*x54 + x67 + x73 + x90 + x91 + x92 + x93 + x94)*(-V{0.666666666666667}*cell[10] - V{0.666666666666667}*cell[1] + x102*x54 + x67 + x73 + x90 + x91 + x92 + x93 + x94))) + V{0.0277777691819762}/((x39)*(x39))) + V{0.5}/x39);
auto x109 = V{1} / ((x48)*(x48));
auto x110 = V{1.5}*x109;
auto x111 = x102*x110;
auto x112 = x110*x87;
auto x113 = x110*x66;
auto x114 = x112 + x113 + V{-1};
auto x115 = x111 + x114;
auto x116 = V{1} - x108;
auto x117 = V{1}*cell[0]*x116 - x108*(x115*(x50 + x51 + x52 + x53 + x67 + x68 + x69 + x71 + x72 + x74 + x75 + x76 + x77 + x88 + x89 + x91 + x92 + x93 + x94 + V{0.333333333333333}) + V{0.333333333333333});
auto x118 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x119 = V{3}*x38;
auto x120 = x103*x119;
auto x121 = V{3}*x109;
auto x122 = x102*x121;
auto x123 = V{1}*cell[1]*x116 - x108*(x118*(x114 + x120 - x122) + V{0.0555555555555556});
auto x124 = x105*x119;
auto x125 = x121*x87;
auto x126 = x111 + V{-1};
auto x127 = V{1}*cell[2]*x116 - x108*(x118*(x113 + x124 - x125 + x126) + V{0.0555555555555556});
auto x128 = x119*x49*x64;
auto x129 = x121*x66;
auto x130 = V{1}*cell[3]*x116 - x108*(x118*(x112 + x126 + x128 - x129) + V{0.0555555555555556});
auto x131 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x132 = V{4.5}*x109;
auto x133 = cell[10] + cell[16] + x98;
auto x134 = x132*x55*((2*cell[13] - 2*cell[4] + x133 + x42 + x81)*(2*cell[13] - 2*cell[4] + x133 + x42 + x81));
auto x135 = x115 + x120;
auto x136 = V{1}*cell[4]*x116 - x108*(x131*(x124 - x134 + x135) + V{0.0277777777777778});
auto x137 = -x124;
auto x138 = -cell[11] + V{2}*cell[14] + cell[15] - V{2}*cell[5] + x106 + x133 + x46 + x56;
auto x139 = V{1}*cell[5]*x116 - x108*(x131*(-x132*x55*((x138)*(x138)) + x135 + x137) + V{0.0277777777777778});
auto x140 = cell[10] + cell[14] + x96 + x99;
auto x141 = x132*x55*((cell[12] + 2*cell[15] - 2*cell[6] + x140 + x60 + x78)*(cell[12] + 2*cell[15] - 2*cell[6] + x140 + x60 + x78));
auto x142 = V{1}*cell[6]*x116 - x108*(x131*(x128 + x135 - x141) + V{0.0277777777777778});
auto x143 = -x128;
auto x144 = -cell[12] + x47;
auto x145 = V{2}*cell[16] - V{2}*cell[7] + cell[8] + x107 + x140 + x144 + x80;
auto x146 = V{1}*cell[7]*x116 - x108*(x131*(-x132*x55*((x145)*(x145)) + x135 + x143) + V{0.0277777777777778});
auto x147 = cell[11] + x79 + x84;
auto x148 = x132*x55*((cell[12] + 2*cell[17] - 2*cell[8] + x147 + x57 + x63 + x95)*(cell[12] + 2*cell[17] - 2*cell[8] + x147 + x57 + x63 + x95));
auto x149 = x115 + x124;
auto x150 = V{1}*cell[8]*x116 - x108*(x131*(x128 - x148 + x149) + V{0.0277777777777778});
auto x151 = V{2}*cell[18] + cell[6] - V{2}*cell[9] + x104 + x144 + x147 + x97;
auto x152 = V{1}*cell[9]*x116 - x108*(x131*(-x132*x55*((x151)*(x151)) + x143 + x149) + V{0.0277777777777778});
auto x153 = -x112;
auto x154 = V{1} - x113;
auto x155 = x153 + x154;
auto x156 = x120 + x155;
auto x157 = V{1}*cell[10]*x116 + x108*(x118*(x122 + x156) + V{-0.0555555555555556});
auto x158 = -x111;
auto x159 = x124 + x158;
auto x160 = V{1}*cell[11]*x116 + x108*(x118*(x125 + x154 + x159) + V{-0.0555555555555556});
auto x161 = x128 + x158;
auto x162 = V{1}*cell[12]*x116 + x108*(x118*(x129 + x153 + x161 + V{1}) + V{-0.0555555555555556});
auto x163 = V{1}*cell[13]*x116 + x108*(x131*(x134 + x156 + x159) + V{-0.0277777777777778});
auto x164 = -x120;
auto x165 = V{1}*cell[14]*x116 - x108*(x131*(-x132*x55*((x138)*(x138)) + x149 + x164) + V{0.0277777777777778});
auto x166 = V{1}*cell[15]*x116 + x108*(x131*(x141 + x156 + x161) + V{-0.0277777777777778});
auto x167 = x115 + x128;
auto x168 = V{1}*cell[16]*x116 - x108*(x131*(-x132*x55*((x145)*(x145)) + x164 + x167) + V{0.0277777777777778});
auto x169 = V{1}*cell[17]*x116 + x108*(x131*(x128 + x148 + x155 + x159) + V{-0.0277777777777778});
auto x170 = V{1}*cell[18]*x116 - x108*(x131*(-x132*x55*((x151)*(x151)) + x137 + x167) + V{0.0277777777777778});
auto x171 = V{1} / (x41);
auto x172 = x41 + V{-1};
cell[0] = x117;
cell[1] = x123;
cell[2] = x127;
cell[3] = x130;
cell[4] = x136;
cell[5] = x139;
cell[6] = x142;
cell[7] = x146;
cell[8] = x150;
cell[9] = x152;
cell[10] = x157;
cell[11] = x160;
cell[12] = x162;
cell[13] = x163;
cell[14] = x165;
cell[15] = x166;
cell[16] = x168;
cell[17] = x169;
cell[18] = x170;
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[0] = x171*(x117 + x172*x19);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[1] = x171*(x123 + x172*x20);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[2] = x171*(x127 + x172*x21);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[3] = x171*(x130 + x172*x22);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[4] = x171*(x136 + x172*x23);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[5] = x171*(x139 + x172*x24);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[6] = x171*(x142 + x172*x25);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[7] = x171*(x146 + x172*x26);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[8] = x171*(x150 + x172*x27);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[9] = x171*(x152 + x172*x28);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[10] = x171*(x157 + x172*x29);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[11] = x171*(x160 + x172*x30);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[12] = x171*(x162 + x172*x31);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[13] = x171*(x163 + x172*x32);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[14] = x171*(x165 + x172*x33);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[15] = x171*(x166 + x172*x34);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[16] = x171*(x168 + x172*x35);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[17] = x171*(x169 + x172*x36);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[18] = x171*(x170 + x172*x37);
return { x48, V{1}*x109*x55*(x101 + x65 + x86) };
}
};

}

}

// Generation Info: commit=557a856e999bdc7da25fdc690eab4179ea52272c
