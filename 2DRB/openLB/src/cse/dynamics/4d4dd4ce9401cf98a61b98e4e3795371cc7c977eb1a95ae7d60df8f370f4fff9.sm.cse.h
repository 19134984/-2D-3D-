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
struct CSE<SpongeLayerDynamics<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x111 = cell.template getFieldComponent<olb::descriptors::DAMPING>(0);
auto x112 = cell.template getFieldComponent<olb::descriptors::DENSITY>(0);
auto x113 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x114 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x115 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x116 = parameters.template get<descriptors::OMEGA>();
auto x117 = ((x113)*(x113));
auto x118 = ((x114)*(x114));
auto x119 = ((x115)*(x115));
auto x120 = cell[10] + cell[14] + cell[16];
auto x121 = cell[11] + cell[18] + cell[5];
auto x122 = cell[12] + cell[7] + cell[9];
auto x0 = V{3};
auto x1 = V{0.5}*x0;
auto x2 = x1*(x117 + x118 + x119);
auto x132 = V{1} - x2;
auto x3 = V{0.333333333333333};
auto x4 = -V{1}*cell[4];
auto x124 = cell[13] + x4;
auto x5 = -V{1}*cell[6];
auto x125 = cell[15] + x5;
auto x6 = cell[0] + cell[13] + cell[15] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + x120 + x121 + x122 + V{1};
auto x123 = V{1} / (x6);
auto x7 = x123*(-V{1}*cell[1] - V{1}*cell[5] - V{1}*cell[7] + x120 + x124 + x125);
auto x128 = ((x7)*(x7));
auto x8 = -V{1}*cell[8];
auto x126 = cell[17] + x8;
auto x9 = x123*(-V{1}*cell[14] - V{1}*cell[2] - V{1}*cell[9] + x121 + x124 + x126);
auto x129 = ((x9)*(x9));
auto x10 = x123*(-V{1}*cell[16] - V{1}*cell[18] - V{1}*cell[3] + x122 + x125 + x126);
auto x127 = ((x10)*(x10));
auto x11 = x127 + x128 + x129;
auto x12 = x1*x11;
auto x131 = V{1} - x12;
auto x13 = -x3*(x6*(x12 + V{-1}) + V{1});
auto x14 = cell[0] - x111*(-x112*x3*(V{1} - x2) + x13 + x3) - x116*(cell[0] - x13);
auto x15 = -V{1}*x113;
auto x16 = V{0.5}*((x0)*(x0));
auto x17 = V{0.0555555555555556};
auto x130 = -x17;
auto x18 = x112*x17;
auto x19 = -V{1}*x7;
auto x20 = x17*x6;
auto x21 = x130 + x20*(x0*x19 + x131 + x16*((x19)*(x19)));
auto x22 = -(-cell[1] + x111*(x17 - x18*(x0*x15 + x132 + x16*((x15)*(x15))) + x21) + x116*(cell[1] - x21));
auto x23 = -V{1}*x114;
auto x24 = -V{1}*x9;
auto x25 = x130 + x20*(x0*x24 + x131 + x16*((x24)*(x24)));
auto x26 = -(-cell[2] + x111*(x17 - x18*(x0*x23 + x132 + x16*((x23)*(x23))) + x25) + x116*(cell[2] - x25));
auto x27 = -V{1}*x115;
auto x28 = -V{1}*x10;
auto x29 = x130 + x20*(x0*x28 + x131 + x16*((x28)*(x28)));
auto x30 = -(-cell[3] + x111*(x17 - x18*(x0*x27 + x132 + x16*((x27)*(x27))) + x29) + x116*(cell[3] - x29));
auto x31 = x15 + x23;
auto x32 = V{0.0277777777777778};
auto x133 = -x32;
auto x33 = x112*x32;
auto x34 = x19 + x24;
auto x35 = x32*x6;
auto x36 = x133 + x35*(x0*x34 + x131 + x16*((x34)*(x34)));
auto x37 = -(-cell[4] + x111*(x32 - x33*(x0*x31 + x132 + x16*((x31)*(x31))) + x36) + x116*(cell[4] - x36));
auto x38 = x114 + x15;
auto x39 = x19 + x9;
auto x40 = x133 + x35*(x0*x39 + x131 + x16*((x39)*(x39)));
auto x41 = -(-cell[5] + x111*(x32 - x33*(x0*x38 + x132 + x16*((x38)*(x38))) + x40) + x116*(cell[5] - x40));
auto x42 = x15 + x27;
auto x43 = x19 + x28;
auto x44 = x133 + x35*(x0*x43 + x131 + x16*((x43)*(x43)));
auto x45 = -(-cell[6] + x111*(x32 - x33*(x0*x42 + x132 + x16*((x42)*(x42))) + x44) + x116*(cell[6] - x44));
auto x46 = x115 + x15;
auto x47 = x10 + x19;
auto x48 = x133 + x35*(x0*x47 + x131 + x16*((x47)*(x47)));
auto x49 = -(-cell[7] + x111*(x32 - x33*(x0*x46 + x132 + x16*((x46)*(x46))) + x48) + x116*(cell[7] - x48));
auto x50 = x23 + x27;
auto x51 = x24 + x28;
auto x52 = x133 + x35*(x0*x51 + x131 + x16*((x51)*(x51)));
auto x53 = -(-cell[8] + x111*(x32 - x33*(x0*x50 + x132 + x16*((x50)*(x50))) + x52) + x116*(cell[8] - x52));
auto x54 = x115 + x23;
auto x55 = x10 + x24;
auto x56 = x133 + x35*(x0*x55 + x131 + x16*((x55)*(x55)));
auto x57 = -(-cell[9] + x111*(x32 - x33*(x0*x54 + x132 + x16*((x54)*(x54))) + x56) + x116*(cell[9] - x56));
auto x58 = x130 + x20*(x0*x7 + x128*x16 + x131);
auto x59 = cell[10] - x111*(x17 - x18*(x0*x113 + x117*x16 + x132) + x58) - x116*(cell[10] - x58);
auto x60 = x130 + x20*(x0*x9 + x129*x16 + x131);
auto x61 = cell[11] - x111*(x17 - x18*(x0*x114 + x118*x16 + x132) + x60) - x116*(cell[11] - x60);
auto x62 = x130 + x20*(x0*x10 + x127*x16 + x131);
auto x63 = cell[12] - x111*(x17 - x18*(x0*x115 + x119*x16 + x132) + x62) - x116*(cell[12] - x62);
auto x136 = x49 + x57 + x63;
auto x64 = x113 + x114;
auto x65 = x7 + x9;
auto x66 = x133 + x35*(x0*x65 + x131 + x16*((x65)*(x65)));
auto x67 = -(-cell[13] + x111*(x32 - x33*(x0*x64 + x132 + x16*((x64)*(x64))) + x66) + x116*(cell[13] - x66));
auto x68 = x113 + x23;
auto x69 = x24 + x7;
auto x70 = x133 + x35*(x0*x69 + x131 + x16*((x69)*(x69)));
auto x71 = -(-cell[14] + x111*(x32 - x33*(x0*x68 + x132 + x16*((x68)*(x68))) + x70) + x116*(cell[14] - x70));
auto x72 = x113 + x115;
auto x73 = x10 + x7;
auto x74 = x133 + x35*(x0*x73 + x131 + x16*((x73)*(x73)));
auto x75 = -(-cell[15] + x111*(x32 - x33*(x0*x72 + x132 + x16*((x72)*(x72))) + x74) + x116*(cell[15] - x74));
auto x76 = x113 + x27;
auto x77 = x28 + x7;
auto x78 = x133 + x35*(x0*x77 + x131 + x16*((x77)*(x77)));
auto x79 = -(-cell[16] + x111*(x32 - x33*(x0*x76 + x132 + x16*((x76)*(x76))) + x78) + x116*(cell[16] - x78));
auto x134 = x59 + x71 + x79;
auto x80 = x114 + x115;
auto x81 = x10 + x9;
auto x82 = x133 + x35*(x0*x81 + x131 + x16*((x81)*(x81)));
auto x83 = -(-cell[17] + x111*(x32 - x33*(x0*x80 + x132 + x16*((x80)*(x80))) + x82) + x116*(cell[17] - x82));
auto x84 = x114 + x27;
auto x85 = x28 + x9;
auto x86 = x133 + x35*(x0*x85 + x131 + x16*((x85)*(x85)));
auto x87 = -(-cell[18] + x111*(x32 - x33*(x0*x84 + x132 + x16*((x84)*(x84))) + x86) + x116*(cell[18] - x86));
auto x135 = x41 + x61 + x87;
auto x88 = x134 + x135 + x136 + x14 + x22 + x26 + x30 + x37 + x45 + x53 + x67 + x75 + x83 + V{1};
auto x137 = V{1} / (x88);
auto x89 = -V{1}*x37;
auto x138 = x67 + x89;
auto x90 = -V{1}*x45;
auto x139 = x75 + x90;
auto x91 = -V{1}*x53;
auto x140 = x83 + x91;
cell[0] = x14;
cell[10] = x59;
cell[11] = x61;
cell[12] = x63;
cell[13] = x67;
cell[14] = x71;
cell[15] = x75;
cell[16] = x79;
cell[17] = x83;
cell[18] = x87;
cell[1] = x22;
cell[2] = x26;
cell[3] = x30;
cell[4] = x37;
cell[5] = x41;
cell[6] = x45;
cell[7] = x49;
cell[8] = x53;
cell[9] = x57;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x88;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x137*(x134 + x138 + x139 - V{1}*x22 - V{1}*x41 - V{1}*x49);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x137*(x135 + x138 + x140 - V{1}*x26 - V{1}*x57 - V{1}*x71);
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x137*(x136 + x139 + x140 - V{1}*x30 - V{1}*x79 - V{1}*x87);
return { x6, x11 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
