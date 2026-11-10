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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::PlainGuo>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x103 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x104 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x105 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x106 = parameters.template get<descriptors::OMEGA>();
auto x107 = cell[10] + cell[14] + cell[16];
auto x108 = cell[11] + cell[18] + cell[5];
auto x109 = cell[12] + cell[7] + cell[9];
auto x0 = V{1} - x106;
auto x1 = V{3};
auto x2 = -V{1}*cell[4];
auto x111 = cell[13] + x2;
auto x3 = -V{1}*cell[6];
auto x112 = cell[15] + x3;
auto x4 = cell[0] + cell[13] + cell[15] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + x107 + x108 + x109 + V{1};
auto x110 = V{1} / (x4);
auto x5 = V{0.5}*x103;
auto x6 = x110*(-V{1}*cell[1] - V{1}*cell[5] - V{1}*cell[7] + x107 + x111 + x112) + x5;
auto x115 = ((x6)*(x6));
auto x117 = x1*x6;
auto x7 = -V{1}*cell[8];
auto x113 = cell[17] + x7;
auto x8 = V{0.5}*x104;
auto x9 = x110*(-V{1}*cell[14] - V{1}*cell[2] - V{1}*cell[9] + x108 + x111 + x113) + x8;
auto x116 = ((x9)*(x9));
auto x118 = x1*x9;
auto x10 = V{0.5}*x105;
auto x11 = x10 + x110*(-V{1}*cell[16] - V{1}*cell[18] - V{1}*cell[3] + x109 + x112 + x113);
auto x114 = ((x11)*(x11));
auto x119 = x1*x11;
auto x12 = x114 + x115 + x116;
auto x13 = V{0.5}*x1*x12;
auto x122 = V{1} - x13;
auto x14 = V{0.333333333333333};
auto x15 = -x103*x117;
auto x16 = -x104*x118;
auto x17 = x15 + x16;
auto x18 = -x105*x119;
auto x121 = x16 + x18;
auto x123 = x15 + x18;
auto x19 = -x4*(V{0.5}*x106 + V{-1});
auto x20 = cell[0]*x0 - x106*x14*(x4*(x13 + V{-1}) + V{1}) + x14*x19*(x17 + x18);
auto x21 = -V{1}*x6;
auto x22 = ((x1)*(x1));
auto x120 = V{1}*x22;
auto x23 = V{0.5}*x22;
auto x24 = V{0.0555555555555556};
auto x25 = x24*x4;
auto x26 = -x1*(x6 + V{1});
auto x27 = x19*x24;
auto x28 = cell[1]*x0 - x106*(x24 - x25*(x1*x21 + x122 + x23*((x21)*(x21)))) + x27*(-x103*(x120*x21 - x26) + x121);
auto x29 = -V{1}*x9;
auto x30 = -x1*(x9 + V{1});
auto x31 = cell[2]*x0 - x106*(x24 - x25*(x1*x29 + x122 + x23*((x29)*(x29)))) + x27*(-x104*(x120*x29 - x30) + x123);
auto x32 = -V{1}*x11;
auto x33 = -x1*(x11 + V{1});
auto x34 = -(-cell[3]*x0 + x106*(x24 - x25*(x1*x32 + x122 + x23*((x32)*(x32)))) + x27*(x105*(x120*x32 - x33) - x17));
auto x35 = x21 + x29;
auto x36 = V{0.0277777777777778};
auto x37 = x36*x4;
auto x38 = -x120*x35;
auto x39 = x19*x36;
auto x40 = cell[4]*x0 - x106*(x36 - x37*(x1*x35 + x122 + x23*((x35)*(x35)))) + x39*(x103*(x26 + x38) + x104*(x30 + x38) + x18);
auto x41 = x21 + x9;
auto x42 = x22*x41;
auto x43 = -x1*(x9 + V{-1});
auto x44 = cell[5]*x0 - x106*(x36 - x37*(x1*x41 + x122 + x23*((x41)*(x41)))) + x39*(x103*(x26 - V{1}*x42) + x104*(x42 + x43) + x18);
auto x45 = x21 + x32;
auto x46 = -x120*x45;
auto x47 = cell[6]*x0 - x106*(x36 - x37*(x1*x45 + x122 + x23*((x45)*(x45)))) + x39*(x103*(x26 + x46) + x105*(x33 + x46) + x16);
auto x48 = x11 + x21;
auto x49 = x22*x48;
auto x50 = -x1*(x11 + V{-1});
auto x51 = cell[7]*x0 - x106*(x36 - x37*(x1*x48 + x122 + x23*((x48)*(x48)))) + x39*(x103*(x26 - V{1}*x49) + x105*(x49 + x50) + x16);
auto x52 = x29 + x32;
auto x53 = -x120*x52;
auto x54 = cell[8]*x0 - x106*(x36 - x37*(x1*x52 + x122 + x23*((x52)*(x52)))) + x39*(x104*(x30 + x53) + x105*(x33 + x53) + x15);
auto x55 = x11 + x29;
auto x56 = x22*x55;
auto x57 = cell[9]*x0 - x106*(x36 - x37*(x1*x55 + x122 + x23*((x55)*(x55)))) + x39*(x104*(x30 - V{1}*x56) + x105*(x50 + x56) + x15);
auto x58 = -x1*(x6 + V{-1});
auto x59 = cell[10]*x0 - x106*(x24 - x25*(x115*x23 + x117 + x122)) + x27*(x103*(x22*x6 + x58) + x121);
auto x60 = cell[11]*x0 - x106*(x24 - x25*(x116*x23 + x118 + x122)) + x27*(x104*(x22*x9 + x43) + x123);
auto x61 = cell[12]*x0 - x106*(x24 - x25*(x114*x23 + x119 + x122)) + x27*(x105*(x11*x22 + x50) + x17);
auto x126 = x51 + x57 + x61;
auto x62 = x6 + x9;
auto x63 = x22*x62;
auto x64 = cell[13]*x0 - x106*(x36 - x37*(x1*x62 + x122 + x23*((x62)*(x62)))) + x39*(x103*(x58 + x63) + x104*(x43 + x63) + x18);
auto x65 = x29 + x6;
auto x66 = x22*x65;
auto x67 = cell[14]*x0 - x106*(x36 - x37*(x1*x65 + x122 + x23*((x65)*(x65)))) + x39*(x103*(x58 + x66) + x104*(x30 - V{1}*x66) + x18);
auto x68 = x11 + x6;
auto x69 = x22*x68;
auto x70 = cell[15]*x0 - x106*(x36 - x37*(x1*x68 + x122 + x23*((x68)*(x68)))) + x39*(x103*(x58 + x69) + x105*(x50 + x69) + x16);
auto x71 = x32 + x6;
auto x72 = x22*x71;
auto x73 = cell[16]*x0 - x106*(x36 - x37*(x1*x71 + x122 + x23*((x71)*(x71)))) + x39*(x103*(x58 + x72) + x105*(x33 - V{1}*x72) + x16);
auto x124 = x59 + x67 + x73;
auto x74 = x11 + x9;
auto x75 = x22*x74;
auto x76 = cell[17]*x0 - x106*(x36 - x37*(x1*x74 + x122 + x23*((x74)*(x74)))) + x39*(x104*(x43 + x75) + x105*(x50 + x75) + x15);
auto x77 = x32 + x9;
auto x78 = x22*x77;
auto x79 = cell[18]*x0 - x106*(x36 - x37*(x1*x77 + x122 + x23*((x77)*(x77)))) + x39*(x104*(x43 + x78) + x105*(x33 - V{1}*x78) + x15);
auto x125 = x44 + x60 + x79;
auto x80 = x124 + x125 + x126 + x20 + x28 + x31 + x34 + x40 + x47 + x54 + x64 + x70 + x76 + V{1};
auto x127 = V{1} / (x80);
auto x81 = -V{1}*x40;
auto x128 = x64 + x81;
auto x82 = -V{1}*x47;
auto x129 = x70 + x82;
auto x83 = -V{1}*x54;
auto x130 = x76 + x83;
cell[0] = x20;
cell[10] = x59;
cell[11] = x60;
cell[12] = x61;
cell[13] = x64;
cell[14] = x67;
cell[15] = x70;
cell[16] = x73;
cell[17] = x76;
cell[18] = x79;
cell[1] = x28;
cell[2] = x31;
cell[3] = x34;
cell[4] = x40;
cell[5] = x44;
cell[6] = x47;
cell[7] = x51;
cell[8] = x54;
cell[9] = x57;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x80;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x127*(x124 + x128 + x129 - V{1}*x28 - V{1}*x44 - V{1}*x51) + x5;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x127*(x125 + x128 + x130 - V{1}*x31 - V{1}*x57 - V{1}*x67) + x8;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x10 + x127*(x126 + x129 + x130 - V{1}*x34 - V{1}*x73 - V{1}*x79);
return { x4, x12 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
