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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x92 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x93 = parameters.template get<descriptors::OMEGA>();
auto x94 = cell[15] + cell[16];
auto x95 = cell[10] + cell[14] + x94;
auto x96 = cell[1] + cell[4] + x95;
auto x97 = cell[12] + cell[7] + cell[9];
auto x98 = cell[3] + cell[8] + x97;
auto x99 = cell[13] + cell[5];
auto x100 = cell[17] + cell[18];
auto x101 = cell[11] + x100 + x99;
auto x102 = cell[2] + x101;
auto x0 = V{1}/x93;
auto x1 = -V{1}*cell[16];
auto x104 = cell[15] + x1;
auto x2 = -V{1}*cell[7];
auto x3 = -V{1}*cell[6];
auto x4 = -V{1}*cell[8];
auto x5 = -V{1}*cell[18];
auto x105 = cell[17] + x5;
auto x6 = cell[0] + cell[6] + x102 + x96 + x98 + V{1};
auto x103 = V{1} / (x6);
auto x7 = x103*(-V{1}*cell[3] + x104 + x105 + x3 + x4 + x97);
auto x109 = ((x7)*(x7));
auto x8 = -V{1}*cell[4];
auto x9 = -V{1}*cell[5];
auto x106 = cell[13] + x9;
auto x10 = x103*(-V{1}*cell[1] + x106 + x2 + x3 + x8 + x95);
auto x111 = ((x10)*(x10));
auto x11 = x10*x6;
auto x12 = cell[6] + x104 - x11*x7 + x2;
auto x13 = ((x12)*(x12));
auto x14 = -V{1}*cell[14];
auto x15 = -V{1}*cell[9];
auto x16 = x103*(-V{1}*cell[2] + x101 + x14 + x15 + x4 + x8);
auto x112 = ((x16)*(x16));
auto x17 = cell[4] + x106 - x11*x16 + x14;
auto x18 = V{3};
auto x110 = ((x18)*(x18));
auto x19 = V{1}*(x6 + V{-1})/x18;
auto x107 = -x19;
auto x108 = cell[6] + x107;
auto x20 = cell[7] - x10*x11 + x108 + x96 + x99;
auto x21 = x100 + x108 - x109*x6 + x94 + x98;
auto x22 = x16*x6;
auto x23 = cell[8] + x105 + x15 - x22*x7;
auto x24 = cell[14] + cell[4] + cell[8] + cell[9] + x102 + x107 - x16*x22;
auto x25 = V{2}*V{1} / (x0 + util::sqrt(V{2.828428}*x103*x110*((x92)*(x92))*util::sqrt(V{2}*x13 + V{2}*((x17)*(x17)) + ((x20)*(x20)) + ((x21)*(x21)) + V{2}*((x23)*(x23)) + ((x24)*(x24))) + ((x0)*(x0))));
auto x26 = V{1} - x25;
auto x27 = x109 + x111 + x112;
auto x28 = V{0.5}*x18*x27;
auto x113 = V{1} - x28;
auto x29 = V{0.333333333333333};
auto x30 = cell[0]*x26 - x25*x29*(x6*(x28 + V{-1}) + V{1});
auto x31 = -V{1}*x10;
auto x32 = V{0.5}*x110;
auto x33 = V{0.0555555555555556};
auto x34 = x33*x6;
auto x35 = cell[1]*x26 - x25*(x33 - x34*(x113 + x18*x31 + x32*((x31)*(x31))));
auto x36 = -V{1}*x16;
auto x37 = cell[2]*x26 - x25*(x33 - x34*(x113 + x18*x36 + x32*((x36)*(x36))));
auto x38 = -V{1}*x7;
auto x39 = cell[3]*x26 - x25*(x33 - x34*(x113 + x18*x38 + x32*((x38)*(x38))));
auto x40 = x31 + x36;
auto x41 = V{0.0277777777777778};
auto x42 = x41*x6;
auto x43 = cell[4]*x26 - x25*(x41 - x42*(x113 + x18*x40 + x32*((x40)*(x40))));
auto x44 = x16 + x31;
auto x45 = cell[5]*x26 - x25*(x41 - x42*(x113 + x18*x44 + x32*((x44)*(x44))));
auto x46 = x31 + x38;
auto x47 = cell[6]*x26 - x25*(x41 - x42*(x113 + x18*x46 + x32*((x46)*(x46))));
auto x48 = x31 + x7;
auto x49 = cell[7]*x26 - x25*(x41 - x42*(x113 + x18*x48 + x32*((x48)*(x48))));
auto x50 = x36 + x38;
auto x51 = cell[8]*x26 - x25*(x41 - x42*(x113 + x18*x50 + x32*((x50)*(x50))));
auto x52 = x36 + x7;
auto x53 = cell[9]*x26 - x25*(x41 - x42*(x113 + x18*x52 + x32*((x52)*(x52))));
auto x54 = cell[10]*x26 - x25*(x33 - x34*(x10*x18 + x111*x32 + x113));
auto x55 = cell[11]*x26 - x25*(x33 - x34*(x112*x32 + x113 + x16*x18));
auto x56 = cell[12]*x26 - x25*(x33 - x34*(x109*x32 + x113 + x18*x7));
auto x116 = x49 + x53 + x56;
auto x57 = x10 + x16;
auto x58 = cell[13]*x26 - x25*(x41 - x42*(x113 + x18*x57 + x32*((x57)*(x57))));
auto x59 = x10 + x36;
auto x60 = cell[14]*x26 - x25*(x41 - x42*(x113 + x18*x59 + x32*((x59)*(x59))));
auto x61 = x10 + x7;
auto x62 = cell[15]*x26 - x25*(x41 - x42*(x113 + x18*x61 + x32*((x61)*(x61))));
auto x63 = x10 + x38;
auto x64 = cell[16]*x26 - x25*(x41 - x42*(x113 + x18*x63 + x32*((x63)*(x63))));
auto x114 = x54 + x60 + x64;
auto x65 = x16 + x7;
auto x66 = cell[17]*x26 - x25*(x41 - x42*(x113 + x18*x65 + x32*((x65)*(x65))));
auto x67 = x16 + x38;
auto x68 = cell[18]*x26 - x25*(x41 - x42*(x113 + x18*x67 + x32*((x67)*(x67))));
auto x115 = x45 + x55 + x68;
auto x69 = x114 + x115 + x116 + x30 + x35 + x37 + x39 + x43 + x47 + x51 + x58 + x62 + x66 + V{1};
auto x117 = V{1} / (x69);
auto x70 = -V{1}*x43;
auto x118 = x58 + x70;
auto x71 = -V{1}*x47;
auto x119 = x62 + x71;
auto x72 = -V{1}*x51;
auto x120 = x66 + x72;
cell[0] = x30;
cell[10] = x54;
cell[11] = x55;
cell[12] = x56;
cell[13] = x58;
cell[14] = x60;
cell[15] = x62;
cell[16] = x64;
cell[17] = x66;
cell[18] = x68;
cell[1] = x35;
cell[2] = x37;
cell[3] = x39;
cell[4] = x43;
cell[5] = x45;
cell[6] = x47;
cell[7] = x49;
cell[8] = x51;
cell[9] = x53;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x69;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x117*(x114 + x118 + x119 - V{1}*x35 - V{1}*x45 - V{1}*x49);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x117*(x115 + x118 + x120 - V{1}*x37 - V{1}*x53 - V{1}*x60);
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x117*(x116 + x119 + x120 - V{1}*x39 - V{1}*x64 - V{1}*x68);
return { x6, x27 };
}
};

}

}

// Generation Info: commit=61323934e549f5f16850d6c6280b91bc057cb6c3
