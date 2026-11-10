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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::ConstRhoBGK, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x95 = parameters.template get<descriptors::OMEGA>();
auto x96 = parameters.template get<statistics::AVERAGE_RHO>();
auto x97 = cell[10] + cell[14] + cell[16];
auto x98 = cell[11] + cell[18] + cell[5];
auto x99 = cell[12] + cell[7] + cell[9];
auto x0 = V{3};
auto x1 = -V{1}*cell[4];
auto x101 = cell[13] + x1;
auto x2 = -V{1}*cell[6];
auto x102 = cell[15] + x2;
auto x3 = cell[0] + cell[13] + cell[15] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + x97 + x98 + x99 + V{1};
auto x100 = V{1} / (x3);
auto x4 = x100*(-V{1}*cell[1] - V{1}*cell[5] - V{1}*cell[7] + x101 + x102 + x97);
auto x104 = ((x4)*(x4));
auto x5 = -V{1}*cell[8];
auto x103 = cell[17] + x5;
auto x6 = x100*(-V{1}*cell[14] - V{1}*cell[2] - V{1}*cell[9] + x101 + x103 + x98);
auto x105 = ((x6)*(x6));
auto x7 = x100*(-V{1}*cell[16] - V{1}*cell[18] - V{1}*cell[3] + x102 + x103 + x99);
auto x106 = ((x7)*(x7));
auto x8 = x104 + x105 + x106;
auto x9 = V{0.5}*x0*x8;
auto x108 = V{1} - x9;
auto x10 = V{0.333333333333333};
auto x11 = -x10*(x3*(x9 + V{-1}) + V{1});
auto x12 = V{1} - x95;
auto x13 = V{1} - x96;
auto x14 = x100*x13 + V{1};
auto x15 = -x10 + x12*(cell[0] - x11) + x14*(x10 + x11);
auto x16 = -V{1}*x4;
auto x17 = V{0.5}*((x0)*(x0));
auto x18 = V{0.0555555555555556};
auto x107 = -x18;
auto x19 = x18*x3;
auto x20 = x107 + x19*(x0*x16 + x108 + x17*((x16)*(x16)));
auto x21 = x107 + x12*(cell[1] - x20) + x14*(x18 + x20);
auto x22 = -V{1}*x6;
auto x23 = x107 + x19*(x0*x22 + x108 + x17*((x22)*(x22)));
auto x24 = x107 + x12*(cell[2] - x23) + x14*(x18 + x23);
auto x25 = -V{1}*x7;
auto x26 = x107 + x19*(x0*x25 + x108 + x17*((x25)*(x25)));
auto x27 = x107 + x12*(cell[3] - x26) + x14*(x18 + x26);
auto x28 = x16 + x22;
auto x29 = V{0.0277777777777778};
auto x109 = -x29;
auto x30 = x29*x3;
auto x31 = x109 + x30*(x0*x28 + x108 + x17*((x28)*(x28)));
auto x32 = x109 + x12*(cell[4] - x31) + x14*(x29 + x31);
auto x33 = x16 + x6;
auto x34 = x109 + x30*(x0*x33 + x108 + x17*((x33)*(x33)));
auto x35 = x109 + x12*(cell[5] - x34) + x14*(x29 + x34);
auto x36 = x16 + x25;
auto x37 = x109 + x30*(x0*x36 + x108 + x17*((x36)*(x36)));
auto x38 = x109 + x12*(cell[6] - x37) + x14*(x29 + x37);
auto x39 = x16 + x7;
auto x40 = x109 + x30*(x0*x39 + x108 + x17*((x39)*(x39)));
auto x41 = x109 + x12*(cell[7] - x40) + x14*(x29 + x40);
auto x42 = x22 + x25;
auto x43 = x109 + x30*(x0*x42 + x108 + x17*((x42)*(x42)));
auto x44 = x109 + x12*(cell[8] - x43) + x14*(x29 + x43);
auto x45 = x22 + x7;
auto x46 = x109 + x30*(x0*x45 + x108 + x17*((x45)*(x45)));
auto x47 = x109 + x12*(cell[9] - x46) + x14*(x29 + x46);
auto x48 = x107 + x19*(x0*x4 + x104*x17 + x108);
auto x49 = x107 + x12*(cell[10] - x48) + x14*(x18 + x48);
auto x50 = x107 + x19*(x0*x6 + x105*x17 + x108);
auto x51 = x107 + x12*(cell[11] - x50) + x14*(x18 + x50);
auto x52 = x107 + x19*(x0*x7 + x106*x17 + x108);
auto x53 = x107 + x12*(cell[12] - x52) + x14*(x18 + x52);
auto x112 = x41 + x47 + x53;
auto x54 = x4 + x6;
auto x55 = x109 + x30*(x0*x54 + x108 + x17*((x54)*(x54)));
auto x56 = x109 + x12*(cell[13] - x55) + x14*(x29 + x55);
auto x57 = x22 + x4;
auto x58 = x109 + x30*(x0*x57 + x108 + x17*((x57)*(x57)));
auto x59 = x109 + x12*(cell[14] - x58) + x14*(x29 + x58);
auto x60 = x4 + x7;
auto x61 = x109 + x30*(x0*x60 + x108 + x17*((x60)*(x60)));
auto x62 = x109 + x12*(cell[15] - x61) + x14*(x29 + x61);
auto x63 = x25 + x4;
auto x64 = x109 + x30*(x0*x63 + x108 + x17*((x63)*(x63)));
auto x65 = x109 + x12*(cell[16] - x64) + x14*(x29 + x64);
auto x110 = x49 + x59 + x65;
auto x66 = x6 + x7;
auto x67 = x109 + x30*(x0*x66 + x108 + x17*((x66)*(x66)));
auto x68 = x109 + x12*(cell[17] - x67) + x14*(x29 + x67);
auto x69 = x25 + x6;
auto x70 = x109 + x30*(x0*x69 + x108 + x17*((x69)*(x69)));
auto x71 = x109 + x12*(cell[18] - x70) + x14*(x29 + x70);
auto x111 = x35 + x51 + x71;
auto x72 = x110 + x111 + x112 + x15 + x21 + x24 + x27 + x32 + x38 + x44 + x56 + x62 + x68 + V{1};
auto x113 = V{1} / (x72);
auto x73 = -V{1}*x32;
auto x114 = x56 + x73;
auto x74 = -V{1}*x38;
auto x115 = x62 + x74;
auto x75 = -V{1}*x44;
auto x116 = x68 + x75;
cell[0] = x15;
cell[10] = x49;
cell[11] = x51;
cell[12] = x53;
cell[13] = x56;
cell[14] = x59;
cell[15] = x62;
cell[16] = x65;
cell[17] = x68;
cell[18] = x71;
cell[1] = x21;
cell[2] = x24;
cell[3] = x27;
cell[4] = x32;
cell[5] = x35;
cell[6] = x38;
cell[7] = x41;
cell[8] = x44;
cell[9] = x47;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x72;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x113*(x110 + x114 + x115 - V{1}*x21 - V{1}*x35 - V{1}*x41);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x113*(x111 + x114 + x116 - V{1}*x24 - V{1}*x47 - V{1}*x59);
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x113*(x112 + x115 + x116 - V{1}*x27 - V{1}*x65 - V{1}*x71);
return { x13 + x3, x8 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
