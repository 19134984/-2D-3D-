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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::TrackAveragePopulation<collision::BGK>, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x24 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(5);
auto x34 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(15);
auto x30 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(11);
auto x28 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(9);
auto x35 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(16);
auto x27 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(8);
auto x37 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(18);
auto x32 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(13);
auto x20 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(1);
auto x31 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(12);
auto x21 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(2);
auto x23 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(4);
auto x26 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(7);
auto x33 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(14);
auto x40 = parameters.template get<descriptors::LATTICE_TIME_T>();
auto x22 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(3);
auto x25 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(6);
auto x36 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(17);
auto x38 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x39 = parameters.template get<descriptors::OMEGA>();
auto x19 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(0);
auto x29 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(10);
auto x41 = x39 + V{-1};
auto x42 = cell[10] + cell[14];
auto x43 = cell[12] + cell[7];
auto x44 = cell[0] + cell[11] + cell[13] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[8] + cell[9] + x42 + x43;
auto x45 = x44 + V{1};
auto x46 = x44 + V{1};
auto x47 = V{1} / ((x46)*(x46));
auto x48 = V{1.5}*x47;
auto x49 = ((x38)*(x38));
auto x50 = cell[13] - cell[4];
auto x51 = cell[15] - cell[6];
auto x52 = x50 + x51;
auto x53 = -cell[1];
auto x54 = cell[16] - cell[7];
auto x55 = x53 + x54;
auto x56 = -cell[5] + x42;
auto x57 = x52 + x55 + x56;
auto x58 = ((x57)*(x57));
auto x59 = x49*x58;
auto x60 = x48*x59;
auto x61 = cell[17] - cell[8];
auto x62 = x50 + x61;
auto x63 = cell[18] - cell[9];
auto x64 = -cell[2];
auto x65 = cell[11] - cell[14] + cell[5] + x64;
auto x66 = x62 + x63 + x65;
auto x67 = ((x66)*(x66));
auto x68 = x49*x67;
auto x69 = x48*x68;
auto x70 = x51 + x61;
auto x71 = -cell[3];
auto x72 = -cell[18] + cell[9];
auto x73 = x71 + x72;
auto x74 = -cell[16] + x43;
auto x75 = x70 + x73 + x74;
auto x76 = ((x75)*(x75));
auto x77 = x49*x76;
auto x78 = x48*x77;
auto x79 = x69 + x78 + V{-1};
auto x80 = x60 + x79;
auto x81 = cell[0]*x41 + V{0.333333333333333}*x39*(x45*x80 + V{1});
auto x82 = V{0.0555555555555556}*x39;
auto x83 = V{3}*x38/x46;
auto x84 = x57*x83;
auto x85 = V{3}*x47;
auto x86 = x59*x85;
auto x87 = cell[1]*x41 + x82*(x45*(x79 + x84 - x86) + V{1});
auto x88 = x66*x83;
auto x89 = x68*x85;
auto x90 = x60 + V{-1};
auto x91 = cell[2]*x41 + x82*(x45*(x78 + x88 - x89 + x90) + V{1});
auto x92 = x75*x83;
auto x93 = x77*x85;
auto x94 = cell[3]*x41 + x82*(x45*(x69 + x90 + x92 - x93) + V{1});
auto x95 = V{0.0277777777777778}*x39;
auto x96 = V{4.5}*x47;
auto x97 = cell[10] + x55;
auto x98 = x49*x96*((cell[11] + 2*cell[13] - 2*cell[4] + x63 + x64 + x70 + x97)*(cell[11] + 2*cell[13] - 2*cell[4] + x63 + x64 + x70 + x97));
auto x99 = x80 + x84;
auto x100 = cell[4]*x41 + x95*(x45*(x88 - x98 + x99) + V{1});
auto x101 = -x88;
auto x102 = -cell[17] + cell[8];
auto x103 = -cell[11] + V{2}*cell[14] + cell[2] - V{2}*cell[5] + x102 + x51 + x72 + x97;
auto x104 = cell[5]*x41 + x95*(x45*(x101 - x49*x96*((x103)*(x103)) + x99) + V{1});
auto x105 = x53 + x56;
auto x106 = x49*x96*((cell[12] + 2*cell[15] - 2*cell[6] + x105 + x62 + x73)*(cell[12] + 2*cell[15] - 2*cell[6] + x105 + x62 + x73));
auto x107 = cell[6]*x41 + x95*(x45*(-x106 + x92 + x99) + V{1});
auto x108 = -x92;
auto x109 = -cell[12] + cell[3] + x50;
auto x110 = V{2}*cell[16] - V{2}*cell[7] + x102 + x105 + x109 + x63;
auto x111 = cell[7]*x41 + x95*(x45*(x108 - x49*x96*((x110)*(x110)) + x99) + V{1});
auto x112 = x49*x96*((2*cell[17] - 2*cell[8] + x52 + x65 + x71 + x74)*(2*cell[17] - 2*cell[8] + x52 + x65 + x71 + x74));
auto x113 = x80 + x88;
auto x114 = cell[8]*x41 + x95*(x45*(-x112 + x113 + x92) + V{1});
auto x115 = -cell[15] + V{2}*cell[18] + cell[6] - V{2}*cell[9] + x109 + x54 + x65;
auto x116 = cell[9]*x41 + x95*(x45*(x108 + x113 - x49*x96*((x115)*(x115))) + V{1});
auto x117 = -x69;
auto x118 = V{1} - x78;
auto x119 = x117 + x118;
auto x120 = x119 + x84;
auto x121 = -cell[10]*x41 + x82*(x45*(x120 + x86) + V{-1});
auto x122 = -x60;
auto x123 = x122 + x88;
auto x124 = -cell[11]*x41 + x82*(x45*(x118 + x123 + x89) + V{-1});
auto x125 = x122 + x92;
auto x126 = -cell[12]*x41 + x82*(x45*(x117 + x125 + x93 + V{1}) + V{-1});
auto x127 = -cell[13]*x41 + x95*(x45*(x120 + x123 + x98) + V{-1});
auto x128 = -x84;
auto x129 = cell[14]*x41 + x95*(x45*(x113 + x128 - x49*x96*((x103)*(x103))) + V{1});
auto x130 = -cell[15]*x41 + x95*(x45*(x106 + x120 + x125) + V{-1});
auto x131 = x80 + x92;
auto x132 = cell[16]*x41 + x95*(x45*(x128 + x131 - x49*x96*((x110)*(x110))) + V{1});
auto x133 = -cell[17]*x41 + x95*(x45*(x112 + x119 + x123 + x92) + V{-1});
auto x134 = cell[18]*x41 + x95*(x45*(x101 + x131 - x49*x96*((x115)*(x115))) + V{1});
auto x135 = V{1} / (x40);
auto x136 = x40 + V{-1};
cell[0] = -x81;
cell[1] = -x87;
cell[2] = -x91;
cell[3] = -x94;
cell[4] = -x100;
cell[5] = -x104;
cell[6] = -x107;
cell[7] = -x111;
cell[8] = -x114;
cell[9] = -x116;
cell[10] = x121;
cell[11] = x124;
cell[12] = x126;
cell[13] = x127;
cell[14] = -x129;
cell[15] = x130;
cell[16] = -x132;
cell[17] = x133;
cell[18] = -x134;
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[0] = x135*(x136*x19 - x81);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[1] = x135*(x136*x20 - x87);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[2] = x135*(x136*x21 - x91);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[3] = x135*(x136*x22 - x94);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[4] = x135*(-x100 + x136*x23);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[5] = x135*(-x104 + x136*x24);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[6] = x135*(-x107 + x136*x25);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[7] = x135*(-x111 + x136*x26);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[8] = x135*(-x114 + x136*x27);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[9] = x135*(-x116 + x136*x28);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[10] = x135*(x121 + x136*x29);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[11] = x135*(x124 + x136*x30);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[12] = x135*(x126 + x136*x31);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[13] = x135*(x127 + x136*x32);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[14] = x135*(-x129 + x136*x33);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[15] = x135*(x130 + x136*x34);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[16] = x135*(-x132 + x136*x35);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[17] = x135*(x133 + x136*x36);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[18] = x135*(-x134 + x136*x37);
return { x46, V{1}*x47*x49*(x58 + x67 + x76) };
}
};

}

}

// Generation Info: commit=0b639c483499c32c8066b84773b2d9417797d6d8
