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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<1, -1>, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x148 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x149 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x150 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x151 = parameters.template get<descriptors::OMEGA>();
auto x152 = ((x148)*(x148));
auto x153 = ((x149)*(x149));
auto x154 = ((x150)*(x150));
auto x156 = cell[10] + cell[1];
auto x157 = cell[12] + cell[3];
auto x0 = V{3};
auto x158 = ((x0)*(x0));
auto x1 = x152 + x153 + x154;
auto x2 = V{0.5}*x0*x1;
auto x159 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = -V{1}*x149;
auto x6 = x5 + V{1};
auto x155 = V{1} / (x6);
auto x7 = x155*(cell[0] + V{2}*cell[14] + cell[15] + cell[16] + V{2}*cell[2] + V{2}*cell[4] + cell[6] + cell[7] + V{2}*cell[8] + V{2}*cell[9] + x156 + x157 + V{1});
auto x8 = V{1}/x0;
auto x9 = -x8;
auto x10 = V{0.0277777777777778};
auto x166 = -x10;
auto x11 = x10*x7;
auto x12 = -V{1}*x150;
auto x13 = x12 + x148;
auto x14 = V{0.5}*x158;
auto x15 = x0*x13 + x14*((x13)*(x13)) + x159;
auto x16 = -x10 + x11*x15;
auto x17 = cell[16] - x16;
auto x18 = x148 + x150;
auto x19 = x0*x18 + x14*((x18)*(x18)) + x159;
auto x20 = -x10 + x11*x19;
auto x21 = cell[15] - x20;
auto x22 = V{0.0555555555555556};
auto x165 = -x22;
auto x23 = x22*x7;
auto x24 = x0*x148 + x14*x152 + x159;
auto x25 = -x22 + x23*x24;
auto x26 = -V{1}*x148;
auto x27 = x150 + x26;
auto x28 = x0*x27 + x14*((x27)*(x27)) + x159;
auto x29 = -x10 + x11*x28;
auto x30 = cell[7] - x29;
auto x31 = x0*x26 + x14*((x26)*(x26)) + x159;
auto x32 = -x22 + x23*x31;
auto x33 = x12 + x26;
auto x34 = x0*x33 + x14*((x33)*(x33)) + x159;
auto x35 = -x10 + x11*x34;
auto x36 = cell[6] - x35;
auto x160 = x17 + x21 + x30 + x36;
auto x37 = V{-2};
auto x38 = -V{1}*x37;
auto x39 = x26 + x5;
auto x40 = x0*x39 + x14*((x39)*(x39)) + x159;
auto x41 = -x10 + x11*x40;
auto x42 = x38*(cell[4] - x41);
auto x43 = x148 + x5;
auto x44 = x0*x43 + x14*((x43)*(x43)) + x159;
auto x45 = -x10 + x11*x44;
auto x46 = cell[14] - x45;
auto x47 = x156 + x160 - x25 - x32 + x42 + V{2}*x46;
auto x48 = x47*x9;
auto x49 = x0*x5 + x14*((x5)*(x5)) + x159;
auto x50 = -x22 + x23*x49;
auto x51 = x12 + x5;
auto x52 = x0*x51 + x14*((x51)*(x51)) + x159;
auto x53 = -x10 + x11*x52;
auto x54 = x38*(cell[8] - x53);
auto x55 = x150 + x5;
auto x56 = x0*x55 + x14*((x55)*(x55)) + x159;
auto x57 = -x10 + x11*x56;
auto x58 = cell[9] - x57;
auto x59 = x38*x46 + x38*x58 + x38*(cell[2] - x50) + x42 + x54;
auto x60 = x59*x9;
auto x61 = x48 + x60;
auto x62 = x0*x150 + x14*x154 + x159;
auto x63 = -x22 + x23*x62;
auto x64 = x0*x12 + x14*((x12)*(x12)) + x159;
auto x65 = -x22 + x23*x64;
auto x66 = x157 + x160 + x54 + V{2}*x58 - x63 - x65;
auto x67 = x66*x9;
auto x68 = x4*(x158*(V{0.5}*x61 + V{0.5}*x67) + x3*x7 + V{-1});
auto x69 = V{1} - x151;
auto x70 = V{0.5}*x10*x158;
auto x71 = V{1} - x8;
auto x72 = x66*x71;
auto x73 = x47*x71;
auto x161 = x60 + x73;
auto x162 = x161 + x72;
auto x74 = -V{1}*x17 + x21 - V{1}*x30 + x36;
auto x75 = x70*(x162 + x37*x74);
auto x76 = x16 + x75;
auto x77 = V{2};
auto x78 = x70*(x162 + x74*x77);
auto x79 = x20 + x78;
auto x80 = V{0.5}*x158*x22;
auto x81 = x80*(x61 + x72);
auto x82 = x63 + x81;
auto x83 = x80*(x161 + x67);
auto x84 = x25 + x83;
auto x85 = x29 + x75;
auto x86 = x35 + x78;
auto x87 = x65 + x81;
auto x88 = x32 + x83;
auto x89 = x59*x71;
auto x163 = x67 + x73 + x89;
auto x90 = x37*x46 + x42;
auto x91 = x70*(x163 + x37*x90);
auto x92 = x45 + x91;
auto x93 = x48 + x89;
auto x164 = x72 + x93;
auto x94 = x37*x58 + x54;
auto x95 = x70*(x164 + x37*x94);
auto x96 = x57 + x95;
auto x97 = x70*(x164 + x77*x94);
auto x98 = x53 + x97;
auto x99 = x80*(x67 + x93);
auto x100 = x50 + x99;
auto x101 = x70*(x163 + x77*x90);
auto x102 = x101 + x41;
auto x103 = x155*(V{2}*x100 + V{2}*x102 + x68 + x76 + x79 + x82 + x84 + x85 + x86 + x87 + x88 + V{2}*x92 + V{2}*x96 + V{2}*x98 + V{1});
auto x104 = x151*x4*(x103*x3 + V{-1}) + x68*x69;
auto x105 = x103*x22;
auto x106 = x151*(x105*x31 + x165) + x69*x88;
auto x107 = x100*x69 + x151*(x105*x49 + x165);
auto x108 = x151*(x105*x64 + x165) + x69*x87;
auto x109 = x10*x103;
auto x110 = x102*x69 - x151*(x10 - x109*x40);
auto x111 = x149 + x26;
auto x112 = x0*x111 + x14*((x111)*(x111)) + x159;
auto x113 = -x151*(x10 - x109*x34) + x69*x86;
auto x114 = -x151*(x10 - x109*x28) + x69*x85;
auto x115 = -x151*(x10 - x109*x52) + x69*x98;
auto x116 = -x151*(x10 - x109*x56) + x69*x96;
auto x117 = x151*(x105*x24 + x165) + x69*x84;
auto x118 = x0*x149 + x14*x153 + x159;
auto x119 = x151*(x105*x62 + x165) + x69*x82;
auto x120 = x148 + x149;
auto x121 = x0*x120 + x14*((x120)*(x120)) + x159;
auto x122 = -x151*(x10 - x109*x44) + x69*x92;
auto x123 = -x151*(x10 - x109*x19) + x69*x79;
auto x124 = -x151*(x10 - x109*x15) + x69*x76;
auto x125 = x149 + x150;
auto x126 = x0*x125 + x14*((x125)*(x125)) + x159;
auto x127 = x12 + x149;
auto x128 = x0*x127 + x14*((x127)*(x127)) + x159;
cell[0] = x104;
cell[10] = x117;
cell[11] = x151*(x105*x118 + x165) + x69*(x118*x23 + x165 + x99);
cell[12] = x119;
cell[13] = -x151*(x10 - x109*x121) + x69*(x101 + x11*x121 + x166);
cell[14] = x122;
cell[15] = x123;
cell[16] = x124;
cell[17] = -x151*(x10 - x109*x126) + x69*(x11*x126 + x166 + x97);
cell[18] = -x151*(x10 - x109*x128) + x69*(x11*x128 + x166 + x95);
cell[1] = x106;
cell[2] = x107;
cell[3] = x108;
cell[4] = x110;
cell[5] = -x151*(x10 - x109*x112) + x69*(x11*x112 + x166 + x91);
cell[6] = x113;
cell[7] = x114;
cell[8] = x115;
cell[9] = x116;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x155*(x104 + x106 + V{2}*x107 + x108 + V{2}*x110 + x113 + x114 + V{2}*x115 + V{2}*x116 + x117 + x119 + V{2}*x122 + x123 + x124 + V{1});
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x148;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x149;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x150;
return { x103, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
