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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<0, 1>, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x148 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x149 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x150 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x151 = parameters.template get<descriptors::OMEGA>();
auto x152 = ((x148)*(x148));
auto x153 = ((x149)*(x149));
auto x154 = ((x150)*(x150));
auto x156 = V{2}*cell[10];
auto x157 = cell[11] + cell[2];
auto x158 = cell[12] + cell[3];
auto x0 = V{3};
auto x159 = ((x0)*(x0));
auto x1 = x152 + x153 + x154;
auto x2 = V{0.5}*x0*x1;
auto x160 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = x148 + V{1};
auto x155 = V{1} / (x5);
auto x6 = x155*(cell[0] + V{2}*cell[13] + V{2}*cell[14] + V{2}*cell[15] + V{2}*cell[16] + cell[17] + cell[18] + cell[8] + cell[9] + x156 + x157 + x158 + V{1});
auto x7 = V{1}/x0;
auto x8 = -x7;
auto x9 = V{0.0555555555555556};
auto x161 = -x9;
auto x10 = x6*x9;
auto x11 = V{0.5}*x159;
auto x12 = x0*x148 + x11*x152 + x160;
auto x13 = x10*x12 + x161;
auto x14 = V{0.0277777777777778};
auto x167 = -x14;
auto x15 = x14*x6;
auto x16 = x148 + x149;
auto x17 = x0*x16 + x11*((x16)*(x16)) + x160;
auto x18 = -x14 + x15*x17;
auto x19 = V{2}*cell[13] - V{2}*x18;
auto x20 = -V{1}*x149;
auto x21 = x148 + x20;
auto x22 = x0*x21 + x11*((x21)*(x21)) + x160;
auto x23 = -x14 + x15*x22;
auto x24 = cell[14] - x23;
auto x25 = x148 + x150;
auto x26 = x0*x25 + x11*((x25)*(x25)) + x160;
auto x27 = -x14 + x15*x26;
auto x28 = V{2}*cell[15] - V{2}*x27;
auto x29 = -V{1}*x150;
auto x30 = x148 + x29;
auto x31 = x0*x30 + x11*((x30)*(x30)) + x160;
auto x32 = -x14 + x15*x31;
auto x33 = cell[16] - x32;
auto x34 = -V{2}*x13 + x156 + x19 + V{2}*x24 + x28 + V{2}*x33;
auto x35 = x34*x8;
auto x36 = x149 + x29;
auto x37 = x0*x36 + x11*((x36)*(x36)) + x160;
auto x38 = -x14 + x15*x37;
auto x39 = cell[18] - x38;
auto x40 = x149 + x150;
auto x41 = x0*x40 + x11*((x40)*(x40)) + x160;
auto x42 = -x14 + x15*x41;
auto x43 = cell[17] - x42;
auto x44 = x0*x149 + x11*x153 + x160;
auto x45 = x10*x44 + x161;
auto x46 = x150 + x20;
auto x47 = x0*x46 + x11*((x46)*(x46)) + x160;
auto x48 = -x14 + x15*x47;
auto x49 = cell[9] - x48;
auto x50 = x0*x20 + x11*((x20)*(x20)) + x160;
auto x51 = x10*x50 + x161;
auto x52 = x20 + x29;
auto x53 = x0*x52 + x11*((x52)*(x52)) + x160;
auto x54 = -x14 + x15*x53;
auto x55 = cell[8] - x54;
auto x162 = x39 + x43 + x49 + x55;
auto x56 = V{-2};
auto x57 = -V{1}*x56;
auto x58 = x157 + x162 + x19 + x24*x57 - x45 - x51;
auto x59 = x58*x8;
auto x60 = x35 + x59;
auto x61 = x0*x150 + x11*x154 + x160;
auto x62 = x10*x61 + x161;
auto x63 = x0*x29 + x11*((x29)*(x29)) + x160;
auto x64 = x10*x63 + x161;
auto x65 = x158 + x162 + x28 + x33*x57 - x62 - x64;
auto x66 = x65*x8;
auto x67 = x4*(x159*(V{0.5}*x60 + V{0.5}*x66) + x3*x6 + V{-1});
auto x68 = V{1} - x151;
auto x69 = V{0.5}*x14*x159;
auto x70 = V{1} - x7;
auto x71 = x65*x70;
auto x72 = x58*x70;
auto x73 = x35 + x72;
auto x163 = x71 + x73;
auto x74 = -V{1}*x39 + x43 - V{1}*x49 + x55;
auto x75 = x69*(x163 + x56*x74);
auto x76 = x38 + x75;
auto x77 = V{2};
auto x78 = x69*(x163 + x74*x77);
auto x79 = x42 + x78;
auto x80 = V{0.5}*x159*x9;
auto x81 = x80*(x60 + x71);
auto x82 = x62 + x81;
auto x83 = x80*(x66 + x73);
auto x84 = x45 + x83;
auto x85 = x48 + x75;
auto x86 = x54 + x78;
auto x87 = x64 + x81;
auto x88 = x51 + x83;
auto x89 = x34*x70;
auto x164 = x59 + x89;
auto x165 = x164 + x71;
auto x166 = x66 + x72 + x89;
auto x90 = x28 + x33*x56;
auto x91 = x69*(x165 + x56*x90);
auto x92 = x32 + x91;
auto x93 = x69*(x165 + x77*x90);
auto x94 = x27 + x93;
auto x95 = x19 + x24*x56;
auto x96 = x69*(x166 + x56*x95);
auto x97 = x23 + x96;
auto x98 = x80*(x164 + x66);
auto x99 = x13 + x98;
auto x100 = x69*(x166 + x77*x95);
auto x101 = x100 + x18;
auto x102 = x155*(V{2}*x101 + x67 + x76 + x79 + x82 + x84 + x85 + x86 + x87 + x88 + V{2}*x92 + V{2}*x94 + V{2}*x97 + V{2}*x99 + V{1});
auto x103 = x151*x4*(x102*x3 + V{-1}) + x67*x68;
auto x104 = -V{1}*x148;
auto x105 = x0*x104 + x11*((x104)*(x104)) + x160;
auto x106 = x102*x9;
auto x107 = x151*(x106*x50 + x161) + x68*x88;
auto x108 = x151*(x106*x63 + x161) + x68*x87;
auto x109 = x104 + x20;
auto x110 = x0*x109 + x11*((x109)*(x109)) + x160;
auto x111 = x102*x14;
auto x112 = x104 + x149;
auto x113 = x0*x112 + x11*((x112)*(x112)) + x160;
auto x114 = x104 + x29;
auto x115 = x0*x114 + x11*((x114)*(x114)) + x160;
auto x116 = x104 + x150;
auto x117 = x0*x116 + x11*((x116)*(x116)) + x160;
auto x118 = x151*(x111*x53 + x167) + x68*x86;
auto x119 = x151*(x111*x47 + x167) + x68*x85;
auto x120 = x151*(x106*x12 + x161) + x68*x99;
auto x121 = x151*(x106*x44 + x161) + x68*x84;
auto x122 = x151*(x106*x61 + x161) + x68*x82;
auto x123 = x101*x68 + x151*(x111*x17 + x167);
auto x124 = x151*(x111*x22 + x167) + x68*x97;
auto x125 = x151*(x111*x26 + x167) + x68*x94;
auto x126 = x151*(x111*x31 + x167) + x68*x92;
auto x127 = x151*(x111*x41 + x167) + x68*x79;
auto x128 = x151*(x111*x37 + x167) + x68*x76;
cell[0] = x103;
cell[10] = x120;
cell[11] = x121;
cell[12] = x122;
cell[13] = x123;
cell[14] = x124;
cell[15] = x125;
cell[16] = x126;
cell[17] = x127;
cell[18] = x128;
cell[1] = x151*(x105*x106 + x161) + x68*(x10*x105 + x161 + x98);
cell[2] = x107;
cell[3] = x108;
cell[4] = x151*(x110*x111 + x167) + x68*(x100 + x110*x15 + x167);
cell[5] = x151*(x111*x113 + x167) + x68*(x113*x15 + x167 + x96);
cell[6] = x151*(x111*x115 + x167) + x68*(x115*x15 + x167 + x93);
cell[7] = x151*(x111*x117 + x167) + x68*(x117*x15 + x167 + x91);
cell[8] = x118;
cell[9] = x119;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x155*(x103 + x107 + x108 + x118 + x119 + V{2}*x120 + x121 + x122 + V{2}*x123 + V{2}*x124 + V{2}*x125 + V{2}*x126 + x127 + x128 + V{1});
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x148;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x149;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x150;
return { x102, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
