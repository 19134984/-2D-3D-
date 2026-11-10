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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineUSeparatelyTrace>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x97 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x98 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x99 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x100 = parameters.template get<descriptors::OMEGA>();
auto x101 = ((x97)*(x97));
auto x102 = ((x98)*(x98));
auto x103 = ((x99)*(x99));
auto x104 = cell[15] + cell[16] + cell[6] + cell[7];
auto x105 = cell[10] + cell[1] + x104;
auto x106 = cell[17] + cell[18] + cell[8] + cell[9];
auto x107 = cell[11] + cell[2] + x106;
auto x108 = cell[12] + cell[3];
auto x0 = V{3};
auto x111 = ((x0)*(x0));
auto x1 = x101 + x102 + x103;
auto x2 = V{0.5}*x0*x1;
auto x114 = V{1} - x2;
auto x3 = V{0.333333333333333};
auto x4 = cell[0] + cell[13] + cell[14] + cell[4] + cell[5] + x105 + x107 + x108 + V{1};
auto x5 = V{1}/x0;
auto x6 = -x5;
auto x7 = x4*x97;
auto x8 = x5*(x4 + V{-1});
auto x109 = -x8;
auto x110 = cell[13] + cell[14] + cell[4] + cell[5] + x109;
auto x9 = x105 + x110 - x7*x97;
auto x10 = x6*x9;
auto x11 = x4*x98;
auto x12 = x107 - x11*x98 + x110;
auto x13 = x12*x6;
auto x14 = x10 + x13;
auto x15 = -x103*x4 + x104 + x106 + x108 + x109;
auto x16 = x15*x6;
auto x17 = V{1} - x100;
auto x18 = x3*(V{0.5}*x111*x17*(x14 + x16) - x4*(x2 + V{-1}) + V{-1});
auto x19 = -V{1}*x97;
auto x20 = V{0.5}*x111;
auto x21 = V{0.0555555555555556};
auto x115 = -x21;
auto x22 = x21*x4;
auto x23 = V{0.5}*x111*x21;
auto x113 = x17*x23;
auto x24 = V{1} - x5;
auto x25 = x24*x9;
auto x112 = x16 + x25;
auto x26 = x113*(x112 + x13);
auto x116 = x115 + x26;
auto x27 = x116 + x22*(x0*x19 + x114 + x20*((x19)*(x19)));
auto x28 = -V{1}*x98;
auto x29 = x12*x24;
auto x119 = x112 + x29;
auto x30 = x10 + x29;
auto x31 = x113*(x16 + x30);
auto x117 = x115 + x31;
auto x32 = x117 + x22*(x0*x28 + x114 + x20*((x28)*(x28)));
auto x33 = -V{1}*x99;
auto x34 = x15*x24;
auto x124 = x13 + x25 + x34;
auto x127 = x30 + x34;
auto x35 = x113*(x14 + x34);
auto x118 = x115 + x35;
auto x36 = x118 + x22*(x0*x33 + x114 + x20*((x33)*(x33)));
auto x37 = x19 + x28;
auto x38 = V{0.0277777777777778};
auto x121 = -x38;
auto x39 = x38*x4;
auto x40 = V{0.5}*x111*x38;
auto x120 = x17*x40;
auto x41 = V{2};
auto x42 = cell[13] - V{1}*cell[14] + cell[4] - V{1}*cell[5] - x7*x98;
auto x43 = x120*(x119 + x41*x42);
auto x122 = x121 + x43;
auto x44 = x122 + x39*(x0*x37 + x114 + x20*((x37)*(x37)));
auto x45 = x19 + x98;
auto x46 = V{-2};
auto x47 = x120*(x119 + x42*x46);
auto x123 = x121 + x47;
auto x48 = x123 + x39*(x0*x45 + x114 + x20*((x45)*(x45)));
auto x49 = x19 + x33;
auto x50 = cell[15] - V{1}*cell[16] + cell[6] - V{1}*cell[7] - x7*x99;
auto x51 = x120*(x124 + x41*x50);
auto x125 = x121 + x51;
auto x52 = x125 + x39*(x0*x49 + x114 + x20*((x49)*(x49)));
auto x53 = x19 + x99;
auto x54 = x120*(x124 + x46*x50);
auto x126 = x121 + x54;
auto x55 = x126 + x39*(x0*x53 + x114 + x20*((x53)*(x53)));
auto x56 = x28 + x33;
auto x57 = cell[17] - V{1}*cell[18] + cell[8] - V{1}*cell[9] - x11*x99;
auto x58 = x120*(x127 + x41*x57);
auto x128 = x121 + x58;
auto x59 = x128 + x39*(x0*x56 + x114 + x20*((x56)*(x56)));
auto x60 = x28 + x99;
auto x61 = x120*(x127 + x46*x57);
auto x129 = x121 + x61;
auto x62 = x129 + x39*(x0*x60 + x114 + x20*((x60)*(x60)));
auto x63 = x116 + x22*(x0*x97 + x101*x20 + x114);
auto x64 = x117 + x22*(x0*x98 + x102*x20 + x114);
auto x65 = x118 + x22*(x0*x99 + x103*x20 + x114);
auto x66 = x97 + x98;
auto x67 = x122 + x39*(x0*x66 + x114 + x20*((x66)*(x66)));
auto x68 = x28 + x97;
auto x69 = x123 + x39*(x0*x68 + x114 + x20*((x68)*(x68)));
auto x70 = x97 + x99;
auto x71 = x125 + x39*(x0*x70 + x114 + x20*((x70)*(x70)));
auto x72 = x33 + x97;
auto x73 = x126 + x39*(x0*x72 + x114 + x20*((x72)*(x72)));
auto x74 = x98 + x99;
auto x75 = x128 + x39*(x0*x74 + x114 + x20*((x74)*(x74)));
auto x76 = x33 + x98;
auto x77 = x129 + x39*(x0*x76 + x114 + x20*((x76)*(x76)));
cell[0] = x18;
cell[10] = x63;
cell[11] = x64;
cell[12] = x65;
cell[13] = x67;
cell[14] = x69;
cell[15] = x71;
cell[16] = x73;
cell[17] = x75;
cell[18] = x77;
cell[1] = x27;
cell[2] = x32;
cell[3] = x36;
cell[4] = x44;
cell[5] = x48;
cell[6] = x52;
cell[7] = x55;
cell[8] = x59;
cell[9] = x62;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x18 + x27 + x32 + x36 + x44 + x48 + x52 + x55 + x59 + x62 + x63 + x64 + x65 + x67 + x69 + x71 + x73 + x75 + x77 + V{1};
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x97;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x98;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x99;
return { x4, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
