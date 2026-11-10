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
struct CSE<dynamics::Tuple<T, descriptors::D3Q15<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, multiphase::OmegaFromCell<collision::BGK>, forcing::Guo<momenta::Forced> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x97 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x98 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x99 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x102 = parameters.template get<multiphase::OMEGA_LIQUID>();
auto x103 = parameters.template get<multiphase::OMEGA_VAPOR>();
auto x104 = parameters.template get<multiphase::RHO_LIQUID>();
auto x105 = parameters.template get<multiphase::RHO_VAPOR>();
auto x100 = cell[14] + cell[8];
auto x101 = cell[6] + cell[9];
auto x106 = cell[10] + cell[11] + cell[5];
auto x107 = -x105;
auto x0 = cell[0] + cell[12] + cell[13] + cell[1] + cell[2] + cell[3] + cell[4] + cell[7] + x100 + x101 + x106 + V{1};
auto x108 = V{1} / (x0);
auto x1 = (x102*(x0 + x107) - x103*(x0 - x104))/(x104 + x107);
auto x2 = V{1} - x1;
auto x3 = V{3};
auto x4 = -V{1}*cell[4];
auto x5 = -V{1}*cell[5];
auto x109 = cell[11] + cell[12] + x4 + x5;
auto x6 = -V{1}*cell[6];
auto x110 = cell[13] + x6;
auto x7 = V{0.5}*x97;
auto x8 = x108*(-V{1}*cell[1] - V{1}*cell[7] + x100 + x109 + x110) + x7;
auto x114 = ((x8)*(x8));
auto x115 = x3*x8;
auto x9 = -V{1}*cell[14];
auto x111 = cell[7] + x9;
auto x10 = V{0.5}*x98;
auto x11 = x10 + x108*(-V{1}*cell[13] - V{1}*cell[2] + x101 + x109 + x111);
auto x112 = ((x11)*(x11));
auto x116 = x11*x3;
auto x12 = V{0.5}*x99;
auto x13 = x108*(-V{1}*cell[12] - V{1}*cell[3] + x106 + x110 + x111 + x4) + x12;
auto x113 = ((x13)*(x13));
auto x117 = x13*x3;
auto x14 = x112 + x113 + x114;
auto x15 = V{0.5}*x14*x3;
auto x120 = V{1} - x15;
auto x16 = V{0.222222222222222};
auto x17 = -x115*x97;
auto x18 = -x116*x98;
auto x19 = x17 + x18;
auto x20 = -x117*x99;
auto x119 = x18 + x20;
auto x121 = x17 + x20;
auto x21 = -x0*(V{0.5}*x1 + V{-1});
auto x22 = cell[0]*x2 - x1*x16*(x0*(x15 + V{-1}) + V{1}) + x16*x21*(x19 + x20);
auto x23 = -V{1}*x8;
auto x24 = ((x3)*(x3));
auto x118 = V{1}*x24;
auto x25 = V{0.5}*x24;
auto x26 = V{0.111111111111111};
auto x27 = x0*x26;
auto x28 = -x3*(x8 + V{1});
auto x29 = x21*x26;
auto x30 = cell[1]*x2 - x1*(x26 - x27*(x120 + x23*x3 + x25*((x23)*(x23)))) + x29*(x119 - x97*(x118*x23 - x28));
auto x31 = -V{1}*x11;
auto x32 = -x3*(x11 + V{1});
auto x33 = cell[2]*x2 - x1*(x26 - x27*(x120 + x25*((x31)*(x31)) + x3*x31)) + x29*(x121 - x98*(x118*x31 - x32));
auto x34 = -V{1}*x13;
auto x35 = -x3*(x13 + V{1});
auto x36 = cell[3]*x2 - x1*(x26 - x27*(x120 + x25*((x34)*(x34)) + x3*x34)) + x29*(x19 - x99*(x118*x34 - x35));
auto x37 = x23 + x31;
auto x38 = x34 + x37;
auto x39 = V{0.0138888888888889};
auto x40 = x0*x39;
auto x41 = -x118*x38;
auto x42 = x21*x39;
auto x43 = cell[4]*x2 - x1*(x39 - x40*(x120 + x25*((x38)*(x38)) + x3*x38)) + x42*(x97*(x28 + x41) + x98*(x32 + x41) + x99*(x35 + x41));
auto x44 = x13 + x37;
auto x45 = x24*x44;
auto x46 = -V{1}*x45;
auto x47 = -x3*(x13 + V{-1});
auto x48 = cell[5]*x2 - x1*(x39 - x40*(x120 + x25*((x44)*(x44)) + x3*x44)) + x42*(x97*(x28 + x46) + x98*(x32 + x46) + x99*(x45 + x47));
auto x49 = x11 + x23;
auto x50 = x34 + x49;
auto x51 = x24*x50;
auto x52 = -V{1}*x51;
auto x53 = -x3*(x11 + V{-1});
auto x54 = cell[6]*x2 - x1*(x39 - x40*(x120 + x25*((x50)*(x50)) + x3*x50)) + x42*(x97*(x28 + x52) + x98*(x51 + x53) + x99*(x35 + x52));
auto x55 = x13 + x49;
auto x56 = x24*x55;
auto x57 = cell[7]*x2 - x1*(x39 - x40*(x120 + x25*((x55)*(x55)) + x3*x55)) + x42*(x97*(x28 - V{1}*x56) + x98*(x53 + x56) + x99*(x47 + x56));
auto x58 = -x3*(x8 + V{-1});
auto x59 = cell[8]*x2 - x1*(x26 - x27*(x114*x25 + x115 + x120)) + x29*(x119 + x97*(x24*x8 + x58));
auto x60 = cell[9]*x2 - x1*(x26 - x27*(x112*x25 + x116 + x120)) + x29*(x121 + x98*(x11*x24 + x53));
auto x123 = x54 + x60;
auto x61 = cell[10]*x2 - x1*(x26 - x27*(x113*x25 + x117 + x120)) + x29*(x19 + x99*(x13*x24 + x47));
auto x62 = x11 + x8;
auto x63 = x13 + x62;
auto x64 = x24*x63;
auto x65 = cell[11]*x2 - x1*(x39 - x40*(x120 + x25*((x63)*(x63)) + x3*x63)) + x42*(x97*(x58 + x64) + x98*(x53 + x64) + x99*(x47 + x64));
auto x124 = x48 + x61 + x65;
auto x66 = x34 + x62;
auto x67 = x24*x66;
auto x68 = cell[12]*x2 - x1*(x39 - x40*(x120 + x25*((x66)*(x66)) + x3*x66)) + x42*(x97*(x58 + x67) + x98*(x53 + x67) + x99*(x35 - V{1}*x67));
auto x69 = x31 + x8;
auto x70 = x13 + x69;
auto x71 = x24*x70;
auto x72 = cell[13]*x2 - x1*(x39 - x40*(x120 + x25*((x70)*(x70)) + x3*x70)) + x42*(x97*(x58 + x71) + x98*(x32 - V{1}*x71) + x99*(x47 + x71));
auto x73 = x34 + x69;
auto x74 = x24*x73;
auto x75 = -V{1}*x74;
auto x76 = cell[14]*x2 - x1*(x39 - x40*(x120 + x25*((x73)*(x73)) + x3*x73)) + x42*(x97*(x58 + x74) + x98*(x32 + x75) + x99*(x35 + x75));
auto x122 = x59 + x76;
auto x77 = x122 + x123 + x124 + x22 + x30 + x33 + x36 + x43 + x57 + x68 + x72 + V{1};
auto x125 = V{1} / (x77);
auto x78 = -V{1}*x43;
auto x79 = -V{1}*x48;
auto x126 = x65 + x68 + x78 + x79;
auto x80 = -V{1}*x54;
auto x127 = x72 + x80;
auto x81 = -V{1}*x76;
auto x128 = x57 + x81;
cell[0] = x22;
cell[10] = x61;
cell[11] = x65;
cell[12] = x68;
cell[13] = x72;
cell[14] = x76;
cell[1] = x30;
cell[2] = x33;
cell[3] = x36;
cell[4] = x43;
cell[5] = x48;
cell[6] = x54;
cell[7] = x57;
cell[8] = x59;
cell[9] = x60;
cell.template getFieldPointer<olb::descriptors::OMEGA>()[0] = x1;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x77;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x125*(x122 + x126 + x127 - V{1}*x30 - V{1}*x57) + x7;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x10 + x125*(x123 + x126 + x128 - V{1}*x33 - V{1}*x72);
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x12 + x125*(x124 + x127 + x128 - V{1}*x36 - V{1}*x68 + x78);
return { x0, x14 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
