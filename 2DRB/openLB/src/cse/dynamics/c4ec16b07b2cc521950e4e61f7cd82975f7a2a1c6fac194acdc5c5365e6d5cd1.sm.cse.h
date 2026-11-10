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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<2, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<2, 1>, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x148 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x149 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x150 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x151 = parameters.template get<descriptors::OMEGA>();
auto x152 = ((x148)*(x148));
auto x153 = ((x149)*(x149));
auto x154 = ((x150)*(x150));
auto x156 = V{2}*cell[12];
auto x157 = cell[10] + cell[1];
auto x158 = cell[11] + cell[2];
auto x0 = V{3};
auto x159 = ((x0)*(x0));
auto x1 = x152 + x153 + x154;
auto x2 = V{0.5}*x0*x1;
auto x160 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = x150 + V{1};
auto x155 = V{1} / (x5);
auto x6 = x155*(cell[0] + cell[13] + cell[14] + V{2}*cell[15] + V{2}*cell[17] + cell[4] + cell[5] + V{2}*cell[7] + V{2}*cell[9] + x156 + x157 + x158 + V{1});
auto x7 = V{1}/x0;
auto x8 = -x7;
auto x9 = V{0.0277777777777778};
auto x161 = -x9;
auto x10 = x6*x9;
auto x11 = -V{1}*x149;
auto x12 = x11 + x148;
auto x13 = V{0.5}*x159;
auto x14 = x0*x12 + x13*((x12)*(x12)) + x160;
auto x15 = x10*x14 + x161;
auto x16 = cell[14] - x15;
auto x17 = x148 + x149;
auto x18 = x0*x17 + x13*((x17)*(x17)) + x160;
auto x19 = x10*x18 + x161;
auto x20 = cell[13] - x19;
auto x21 = V{0.0555555555555556};
auto x167 = -x21;
auto x22 = x21*x6;
auto x23 = x0*x148 + x13*x152 + x160;
auto x24 = -x21 + x22*x23;
auto x25 = -V{1}*x148;
auto x26 = x149 + x25;
auto x27 = x0*x26 + x13*((x26)*(x26)) + x160;
auto x28 = x10*x27 + x161;
auto x29 = cell[5] - x28;
auto x30 = x0*x25 + x13*((x25)*(x25)) + x160;
auto x31 = -x21 + x22*x30;
auto x32 = x11 + x25;
auto x33 = x0*x32 + x13*((x32)*(x32)) + x160;
auto x34 = x10*x33 + x161;
auto x35 = cell[4] - x34;
auto x162 = x16 + x20 + x29 + x35;
auto x36 = V{-2};
auto x37 = -V{1}*x36;
auto x38 = x150 + x25;
auto x39 = x0*x38 + x13*((x38)*(x38)) + x160;
auto x40 = x10*x39 + x161;
auto x41 = cell[7] - x40;
auto x42 = x148 + x150;
auto x43 = x0*x42 + x13*((x42)*(x42)) + x160;
auto x44 = x10*x43 + x161;
auto x45 = V{2}*cell[15] - V{2}*x44;
auto x46 = x157 + x162 - x24 - x31 + x37*x41 + x45;
auto x47 = x46*x8;
auto x48 = x0*x149 + x13*x153 + x160;
auto x49 = -x21 + x22*x48;
auto x50 = x0*x11 + x13*((x11)*(x11)) + x160;
auto x51 = -x21 + x22*x50;
auto x52 = x11 + x150;
auto x53 = x0*x52 + x13*((x52)*(x52)) + x160;
auto x54 = x10*x53 + x161;
auto x55 = cell[9] - x54;
auto x56 = x149 + x150;
auto x57 = x0*x56 + x13*((x56)*(x56)) + x160;
auto x58 = x10*x57 + x161;
auto x59 = V{2}*cell[17] - V{2}*x58;
auto x60 = x158 + x162 + x37*x55 - x49 - x51 + x59;
auto x61 = x60*x8;
auto x62 = x47 + x61;
auto x63 = x0*x150 + x13*x154 + x160;
auto x64 = -x21 + x22*x63;
auto x65 = x156 + V{2}*x41 + x45 + V{2}*x55 + x59 - V{2}*x64;
auto x66 = x65*x8;
auto x67 = x4*(x159*(V{0.5}*x62 + V{0.5}*x66) + x3*x6 + V{-1});
auto x68 = V{1} - x151;
auto x69 = V{0.5}*x159*x9;
auto x70 = V{1} - x7;
auto x71 = x60*x70;
auto x72 = x46*x70;
auto x163 = x66 + x72;
auto x164 = x163 + x71;
auto x73 = -V{1}*x16 + x20 - V{1}*x29 + x35;
auto x74 = x69*(x164 + x36*x73);
auto x75 = x15 + x74;
auto x76 = V{2};
auto x77 = x69*(x164 + x73*x76);
auto x78 = x19 + x77;
auto x79 = V{0.5}*x159*x21;
auto x80 = x47 + x71;
auto x81 = x79*(x66 + x80);
auto x82 = x49 + x81;
auto x83 = x79*(x163 + x61);
auto x84 = x24 + x83;
auto x85 = x28 + x74;
auto x86 = x34 + x77;
auto x87 = x51 + x81;
auto x88 = x31 + x83;
auto x89 = x65*x70;
auto x165 = x80 + x89;
auto x166 = x61 + x72 + x89;
auto x90 = x36*x55 + x59;
auto x91 = x69*(x165 + x76*x90);
auto x92 = x58 + x91;
auto x93 = x36*x41 + x45;
auto x94 = x69*(x166 + x76*x93);
auto x95 = x44 + x94;
auto x96 = x79*(x62 + x89);
auto x97 = x64 + x96;
auto x98 = x69*(x166 + x36*x93);
auto x99 = x40 + x98;
auto x100 = x69*(x165 + x36*x90);
auto x101 = x100 + x54;
auto x102 = x155*(V{2}*x101 + x67 + x75 + x78 + x82 + x84 + x85 + x86 + x87 + x88 + V{2}*x92 + V{2}*x95 + V{2}*x97 + V{2}*x99 + V{1});
auto x103 = x151*x4*(x102*x3 + V{-1}) + x67*x68;
auto x104 = x102*x21;
auto x105 = x151*(x104*x30 + x167) + x68*x88;
auto x106 = x151*(x104*x50 + x167) + x68*x87;
auto x107 = -V{1}*x150;
auto x108 = x0*x107 + x13*((x107)*(x107)) + x160;
auto x109 = x102*x9;
auto x110 = x151*(x109*x33 + x161) + x68*x86;
auto x111 = x151*(x109*x27 + x161) + x68*x85;
auto x112 = x107 + x25;
auto x113 = x0*x112 + x13*((x112)*(x112)) + x160;
auto x114 = x151*(x109*x39 + x161) + x68*x99;
auto x115 = x107 + x11;
auto x116 = x0*x115 + x13*((x115)*(x115)) + x160;
auto x117 = x101*x68 + x151*(x109*x53 + x161);
auto x118 = x151*(x104*x23 + x167) + x68*x84;
auto x119 = x151*(x104*x48 + x167) + x68*x82;
auto x120 = x151*(x104*x63 + x167) + x68*x97;
auto x121 = x151*(x109*x18 + x161) + x68*x78;
auto x122 = x151*(x109*x14 + x161) + x68*x75;
auto x123 = x151*(x109*x43 + x161) + x68*x95;
auto x124 = x107 + x148;
auto x125 = x0*x124 + x13*((x124)*(x124)) + x160;
auto x126 = x151*(x109*x57 + x161) + x68*x92;
auto x127 = x107 + x149;
auto x128 = x0*x127 + x13*((x127)*(x127)) + x160;
cell[0] = x103;
cell[10] = x118;
cell[11] = x119;
cell[12] = x120;
cell[13] = x121;
cell[14] = x122;
cell[15] = x123;
cell[16] = x151*(x109*x125 + x161) + x68*(x10*x125 + x161 + x98);
cell[17] = x126;
cell[18] = x151*(x109*x128 + x161) + x68*(x10*x128 + x100 + x161);
cell[1] = x105;
cell[2] = x106;
cell[3] = x151*(x104*x108 + x167) + x68*(x108*x22 + x167 + x96);
cell[4] = x110;
cell[5] = x111;
cell[6] = x151*(x109*x113 + x161) + x68*(x10*x113 + x161 + x94);
cell[7] = x114;
cell[8] = x151*(x109*x116 + x161) + x68*(x10*x116 + x161 + x91);
cell[9] = x117;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x155*(x103 + x105 + x106 + x110 + x111 + V{2}*x114 + V{2}*x117 + x118 + x119 + V{2}*x120 + x121 + x122 + V{2}*x123 + V{2}*x126 + V{1});
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x148;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x149;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x150;
return { x102, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
