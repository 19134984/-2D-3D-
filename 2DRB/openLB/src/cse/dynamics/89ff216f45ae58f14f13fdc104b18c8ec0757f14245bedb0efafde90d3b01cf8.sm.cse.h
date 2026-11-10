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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::ScaledPlainGuo>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x106 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x107 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x108 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x109 = cell.template getFieldComponent<olb::descriptors::SCALAR>(0);
auto x110 = parameters.template get<descriptors::OMEGA>();
auto x111 = cell[10] + cell[14] + cell[16];
auto x112 = cell[11] + cell[18] + cell[5];
auto x113 = cell[12] + cell[7] + cell[9];
auto x0 = V{1} - x110;
auto x1 = V{3};
auto x2 = -V{1}*cell[4];
auto x115 = cell[13] + x2;
auto x3 = -V{1}*cell[6];
auto x116 = cell[15] + x3;
auto x4 = cell[0] + cell[13] + cell[15] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + x111 + x112 + x113 + V{1};
auto x114 = V{1} / (x4);
auto x5 = V{0.5}*x106;
auto x6 = x114*(-V{1}*cell[1] - V{1}*cell[5] - V{1}*cell[7] + x111 + x115 + x116) + x5;
auto x119 = ((x6)*(x6));
auto x121 = x1*x6;
auto x7 = -V{1}*cell[8];
auto x117 = cell[17] + x7;
auto x8 = V{0.5}*x107;
auto x9 = x114*(-V{1}*cell[14] - V{1}*cell[2] - V{1}*cell[9] + x112 + x115 + x117) + x8;
auto x120 = ((x9)*(x9));
auto x122 = x1*x9;
auto x10 = V{0.5}*x108;
auto x11 = x10 + x114*(-V{1}*cell[16] - V{1}*cell[18] - V{1}*cell[3] + x113 + x116 + x117);
auto x118 = ((x11)*(x11));
auto x123 = x1*x11;
auto x12 = x118 + x119 + x120;
auto x13 = V{0.5}*x1*x12;
auto x126 = V{1} - x13;
auto x14 = V{0.333333333333333};
auto x15 = x106*x109;
auto x16 = -x121*x15;
auto x17 = x107*x109;
auto x18 = -x122*x17;
auto x19 = x16 + x18;
auto x20 = x108*x109;
auto x21 = -x123*x20;
auto x125 = x18 + x21;
auto x127 = x16 + x21;
auto x22 = -x4*(V{0.5}*x110 + V{-1});
auto x23 = cell[0]*x0 - x110*x14*(x4*(x13 + V{-1}) + V{1}) + x14*x22*(x19 + x21);
auto x24 = -V{1}*x6;
auto x25 = ((x1)*(x1));
auto x124 = V{1}*x25;
auto x26 = V{0.5}*x25;
auto x27 = V{0.0555555555555556};
auto x28 = x27*x4;
auto x29 = -x1*(x6 + V{1});
auto x30 = x22*x27;
auto x31 = cell[1]*x0 - x110*(x27 - x28*(x1*x24 + x126 + x26*((x24)*(x24)))) + x30*(x125 - x15*(x124*x24 - x29));
auto x32 = -V{1}*x9;
auto x33 = -x1*(x9 + V{1});
auto x34 = cell[2]*x0 - x110*(x27 - x28*(x1*x32 + x126 + x26*((x32)*(x32)))) + x30*(x127 - x17*(x124*x32 - x33));
auto x35 = -V{1}*x11;
auto x36 = -x1*(x11 + V{1});
auto x37 = cell[3]*x0 - x110*(x27 - x28*(x1*x35 + x126 + x26*((x35)*(x35)))) + x30*(x19 - x20*(x124*x35 - x36));
auto x38 = x24 + x32;
auto x39 = V{0.0277777777777778};
auto x40 = x39*x4;
auto x41 = -x124*x38;
auto x42 = x22*x39;
auto x43 = cell[4]*x0 - x110*(x39 - x40*(x1*x38 + x126 + x26*((x38)*(x38)))) + x42*(x15*(x29 + x41) + x17*(x33 + x41) + x21);
auto x44 = x24 + x9;
auto x45 = x25*x44;
auto x46 = -x1*(x9 + V{-1});
auto x47 = cell[5]*x0 - x110*(x39 - x40*(x1*x44 + x126 + x26*((x44)*(x44)))) + x42*(x15*(x29 - V{1}*x45) + x17*(x45 + x46) + x21);
auto x48 = x24 + x35;
auto x49 = -x124*x48;
auto x50 = cell[6]*x0 - x110*(x39 - x40*(x1*x48 + x126 + x26*((x48)*(x48)))) + x42*(x15*(x29 + x49) + x18 + x20*(x36 + x49));
auto x51 = x11 + x24;
auto x52 = x25*x51;
auto x53 = -x1*(x11 + V{-1});
auto x54 = cell[7]*x0 - x110*(x39 - x40*(x1*x51 + x126 + x26*((x51)*(x51)))) + x42*(x15*(x29 - V{1}*x52) + x18 + x20*(x52 + x53));
auto x55 = x32 + x35;
auto x56 = -x124*x55;
auto x57 = cell[8]*x0 - x110*(x39 - x40*(x1*x55 + x126 + x26*((x55)*(x55)))) + x42*(x16 + x17*(x33 + x56) + x20*(x36 + x56));
auto x58 = x11 + x32;
auto x59 = x25*x58;
auto x60 = cell[9]*x0 - x110*(x39 - x40*(x1*x58 + x126 + x26*((x58)*(x58)))) + x42*(x16 + x17*(x33 - V{1}*x59) + x20*(x53 + x59));
auto x61 = -x1*(x6 + V{-1});
auto x62 = cell[10]*x0 - x110*(x27 - x28*(x119*x26 + x121 + x126)) + x30*(x125 + x15*(x25*x6 + x61));
auto x63 = cell[11]*x0 - x110*(x27 - x28*(x120*x26 + x122 + x126)) + x30*(x127 + x17*(x25*x9 + x46));
auto x64 = cell[12]*x0 - x110*(x27 - x28*(x118*x26 + x123 + x126)) + x30*(x19 + x20*(x11*x25 + x53));
auto x130 = x54 + x60 + x64;
auto x65 = x6 + x9;
auto x66 = x25*x65;
auto x67 = cell[13]*x0 - x110*(x39 - x40*(x1*x65 + x126 + x26*((x65)*(x65)))) + x42*(x15*(x61 + x66) + x17*(x46 + x66) + x21);
auto x68 = x32 + x6;
auto x69 = x25*x68;
auto x70 = cell[14]*x0 - x110*(x39 - x40*(x1*x68 + x126 + x26*((x68)*(x68)))) + x42*(x15*(x61 + x69) + x17*(x33 - V{1}*x69) + x21);
auto x71 = x11 + x6;
auto x72 = x25*x71;
auto x73 = cell[15]*x0 - x110*(x39 - x40*(x1*x71 + x126 + x26*((x71)*(x71)))) + x42*(x15*(x61 + x72) + x18 + x20*(x53 + x72));
auto x74 = x35 + x6;
auto x75 = x25*x74;
auto x76 = cell[16]*x0 - x110*(x39 - x40*(x1*x74 + x126 + x26*((x74)*(x74)))) + x42*(x15*(x61 + x75) + x18 + x20*(x36 - V{1}*x75));
auto x128 = x62 + x70 + x76;
auto x77 = x11 + x9;
auto x78 = x25*x77;
auto x79 = cell[17]*x0 - x110*(x39 - x40*(x1*x77 + x126 + x26*((x77)*(x77)))) + x42*(x16 + x17*(x46 + x78) + x20*(x53 + x78));
auto x80 = x35 + x9;
auto x81 = x25*x80;
auto x82 = cell[18]*x0 - x110*(x39 - x40*(x1*x80 + x126 + x26*((x80)*(x80)))) + x42*(x16 + x17*(x46 + x81) + x20*(x36 - V{1}*x81));
auto x129 = x47 + x63 + x82;
auto x83 = x128 + x129 + x130 + x23 + x31 + x34 + x37 + x43 + x50 + x57 + x67 + x73 + x79 + V{1};
auto x131 = V{1} / (x83);
auto x84 = -V{1}*x43;
auto x132 = x67 + x84;
auto x85 = -V{1}*x50;
auto x133 = x73 + x85;
auto x86 = -V{1}*x57;
auto x134 = x79 + x86;
cell[0] = x23;
cell[10] = x62;
cell[11] = x63;
cell[12] = x64;
cell[13] = x67;
cell[14] = x70;
cell[15] = x73;
cell[16] = x76;
cell[17] = x79;
cell[18] = x82;
cell[1] = x31;
cell[2] = x34;
cell[3] = x37;
cell[4] = x43;
cell[5] = x47;
cell[6] = x50;
cell[7] = x54;
cell[8] = x57;
cell[9] = x60;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x83;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x131*(x128 + x132 + x133 - V{1}*x31 - V{1}*x47 - V{1}*x54) + x5;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x131*(x129 + x132 + x134 - V{1}*x34 - V{1}*x60 - V{1}*x70) + x8;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x10 + x131*(x130 + x133 + x134 - V{1}*x37 - V{1}*x76 - V{1}*x82);
return { x4, x12 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
