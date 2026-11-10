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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SaveVelocity<collision::BGK>, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, -1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<0, -1>, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x137 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x138 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x139 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x140 = parameters.template get<descriptors::OMEGA>();
auto x134 = ((x137)*(x137));
auto x135 = ((x138)*(x138));
auto x136 = ((x139)*(x139));
auto x142 = cell[11] + cell[2];
auto x143 = cell[12] + cell[3];
auto x0 = V{3};
auto x144 = ((x0)*(x0));
auto x1 = x134 + x135 + x136;
auto x2 = V{0.5}*x0*x1;
auto x145 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = -V{1}*x137;
auto x6 = x5 + V{1};
auto x141 = V{1} / (x6);
auto x7 = x141*(cell[0] + cell[17] + cell[18] + V{2}*cell[1] + V{2}*cell[4] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7] + cell[8] + cell[9] + x142 + x143 + V{1});
auto x8 = V{1}/x0;
auto x9 = -x8;
auto x10 = V{-2};
auto x11 = -V{1}*x10;
auto x12 = V{0.0555555555555556};
auto x151 = -x12;
auto x13 = x12*x7;
auto x14 = V{0.5}*x144;
auto x15 = x0*x5 + x14*((x5)*(x5)) + x145;
auto x16 = -x12 + x13*x15;
auto x17 = V{0.0277777777777778};
auto x152 = -x17;
auto x18 = x17*x7;
auto x19 = -V{1}*x138;
auto x20 = x19 + x5;
auto x21 = x0*x20 + x14*((x20)*(x20)) + x145;
auto x22 = -x17 + x18*x21;
auto x23 = x11*(cell[4] - x22);
auto x24 = x138 + x5;
auto x25 = x0*x24 + x14*((x24)*(x24)) + x145;
auto x26 = -x17 + x18*x25;
auto x27 = cell[5] - x26;
auto x28 = -V{1}*x139;
auto x29 = x28 + x5;
auto x30 = x0*x29 + x14*((x29)*(x29)) + x145;
auto x31 = -x17 + x18*x30;
auto x32 = x11*(cell[6] - x31);
auto x33 = x139 + x5;
auto x34 = x0*x33 + x14*((x33)*(x33)) + x145;
auto x35 = -x17 + x18*x34;
auto x36 = cell[7] - x35;
auto x37 = x11*x27 + x11*x36 + x11*(cell[1] - x16) + x23 + x32;
auto x38 = x37*x9;
auto x39 = x138 + x28;
auto x40 = x0*x39 + x14*((x39)*(x39)) + x145;
auto x41 = -x17 + x18*x40;
auto x42 = cell[18] - x41;
auto x43 = x138 + x139;
auto x44 = x0*x43 + x14*((x43)*(x43)) + x145;
auto x45 = -x17 + x18*x44;
auto x46 = cell[17] - x45;
auto x47 = x0*x138 + x135*x14 + x145;
auto x48 = -x12 + x13*x47;
auto x49 = x139 + x19;
auto x50 = x0*x49 + x14*((x49)*(x49)) + x145;
auto x51 = -x17 + x18*x50;
auto x52 = cell[9] - x51;
auto x53 = x0*x19 + x14*((x19)*(x19)) + x145;
auto x54 = -x12 + x13*x53;
auto x55 = x19 + x28;
auto x56 = x0*x55 + x14*((x55)*(x55)) + x145;
auto x57 = -x17 + x18*x56;
auto x58 = cell[8] - x57;
auto x146 = x42 + x46 + x52 + x58;
auto x59 = x142 + x146 + x23 + V{2}*x27 - x48 - x54;
auto x60 = x59*x9;
auto x61 = x38 + x60;
auto x62 = x0*x139 + x136*x14 + x145;
auto x63 = -x12 + x13*x62;
auto x64 = x0*x28 + x14*((x28)*(x28)) + x145;
auto x65 = -x12 + x13*x64;
auto x66 = x143 + x146 + x32 + V{2}*x36 - x63 - x65;
auto x67 = x66*x9;
auto x68 = x4*(x144*(V{0.5}*x61 + V{0.5}*x67) + x3*x7 + V{-1});
auto x69 = V{1} - x140;
auto x70 = V{0.5}*x144*x17;
auto x71 = V{1} - x8;
auto x72 = x66*x71;
auto x73 = x59*x71;
auto x74 = x38 + x73;
auto x147 = x72 + x74;
auto x75 = -V{1}*x42 + x46 - V{1}*x52 + x58;
auto x76 = x70*(x10*x75 + x147);
auto x77 = x41 + x76;
auto x78 = V{2};
auto x79 = x70*(x147 + x75*x78);
auto x80 = x45 + x79;
auto x81 = V{0.5}*x12*x144;
auto x82 = x81*(x61 + x72);
auto x83 = x63 + x82;
auto x84 = x81*(x67 + x74);
auto x85 = x48 + x84;
auto x86 = x51 + x76;
auto x87 = x57 + x79;
auto x88 = x65 + x82;
auto x89 = x54 + x84;
auto x90 = x37*x71;
auto x148 = x60 + x90;
auto x149 = x148 + x72;
auto x150 = x67 + x73 + x90;
auto x91 = x10*x36 + x32;
auto x92 = x70*(x10*x91 + x149);
auto x93 = x35 + x92;
auto x94 = x70*(x149 + x78*x91);
auto x95 = x31 + x94;
auto x96 = x10*x27 + x23;
auto x97 = x70*(x10*x96 + x150);
auto x98 = x26 + x97;
auto x99 = x81*(x148 + x67);
auto x100 = x16 + x99;
auto x101 = x70*(x150 + x78*x96);
auto x102 = x101 + x22;
auto x103 = x141*(V{2}*x100 + V{2}*x102 + x68 + x77 + x80 + x83 + x85 + x86 + x87 + x88 + x89 + V{2}*x93 + V{2}*x95 + V{2}*x98 + V{1});
auto x104 = x103*x12;
auto x105 = x103*x17;
auto x106 = x0*x137 + x134*x14 + x145;
auto x107 = x137 + x138;
auto x108 = x0*x107 + x14*((x107)*(x107)) + x145;
auto x109 = x137 + x19;
auto x110 = x0*x109 + x14*((x109)*(x109)) + x145;
auto x111 = x137 + x139;
auto x112 = x0*x111 + x14*((x111)*(x111)) + x145;
auto x113 = x137 + x28;
auto x114 = x0*x113 + x14*((x113)*(x113)) + x145;
cell[0] = x140*x4*(x103*x3 + V{-1}) + x68*x69;
cell[10] = x140*(x104*x106 + x151) + x69*(x106*x13 + x151 + x99);
cell[11] = x140*(x104*x47 + x151) + x69*x85;
cell[12] = x140*(x104*x62 + x151) + x69*x83;
cell[13] = x140*(x105*x108 + x152) + x69*(x101 + x108*x18 + x152);
cell[14] = x140*(x105*x110 + x152) + x69*(x110*x18 + x152 + x97);
cell[15] = x140*(x105*x112 + x152) + x69*(x112*x18 + x152 + x94);
cell[16] = x140*(x105*x114 + x152) + x69*(x114*x18 + x152 + x92);
cell[17] = x140*(x105*x44 + x152) + x69*x80;
cell[18] = x140*(x105*x40 + x152) + x69*x77;
cell[1] = x100*x69 + x140*(x104*x15 + x151);
cell[2] = x140*(x104*x53 + x151) + x69*x89;
cell[3] = x140*(x104*x64 + x151) + x69*x88;
cell[4] = x102*x69 + x140*(x105*x21 + x152);
cell[5] = x140*(x105*x25 + x152) + x69*x98;
cell[6] = x140*(x105*x30 + x152) + x69*x95;
cell[7] = x140*(x105*x34 + x152) + x69*x93;
cell[8] = x140*(x105*x56 + x152) + x69*x87;
cell[9] = x140*(x105*x50 + x152) + x69*x86;
cell.template getFieldPointer<olb::descriptors::VELOCITY>()[0] = x137;
cell.template getFieldPointer<olb::descriptors::VELOCITY>()[1] = x138;
cell.template getFieldPointer<olb::descriptors::VELOCITY>()[2] = x139;
return { x103, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
