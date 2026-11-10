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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<2, -1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<2, -1>, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x134 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x135 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x136 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x137 = parameters.template get<descriptors::OMEGA>();
auto x138 = ((x134)*(x134));
auto x139 = ((x135)*(x135));
auto x140 = ((x136)*(x136));
auto x142 = cell[10] + cell[1];
auto x143 = cell[11] + cell[2];
auto x0 = V{3};
auto x144 = ((x0)*(x0));
auto x1 = x138 + x139 + x140;
auto x2 = V{0.5}*x0*x1;
auto x145 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = -V{1}*x136;
auto x6 = x5 + V{1};
auto x141 = V{1} / (x6);
auto x7 = x141*(cell[0] + cell[13] + cell[14] + V{2}*cell[16] + V{2}*cell[18] + V{2}*cell[3] + cell[4] + cell[5] + V{2}*cell[6] + V{2}*cell[8] + x142 + x143 + V{1});
auto x8 = V{1}/x0;
auto x9 = -x8;
auto x10 = V{0.0277777777777778};
auto x152 = -x10;
auto x11 = x10*x7;
auto x12 = -V{1}*x135;
auto x13 = x12 + x134;
auto x14 = V{0.5}*x144;
auto x15 = x0*x13 + x14*((x13)*(x13)) + x145;
auto x16 = -x10 + x11*x15;
auto x17 = cell[14] - x16;
auto x18 = x134 + x135;
auto x19 = x0*x18 + x14*((x18)*(x18)) + x145;
auto x20 = -x10 + x11*x19;
auto x21 = cell[13] - x20;
auto x22 = V{0.0555555555555556};
auto x151 = -x22;
auto x23 = x22*x7;
auto x24 = x0*x134 + x138*x14 + x145;
auto x25 = -x22 + x23*x24;
auto x26 = -V{1}*x134;
auto x27 = x135 + x26;
auto x28 = x0*x27 + x14*((x27)*(x27)) + x145;
auto x29 = -x10 + x11*x28;
auto x30 = cell[5] - x29;
auto x31 = x0*x26 + x14*((x26)*(x26)) + x145;
auto x32 = -x22 + x23*x31;
auto x33 = x12 + x26;
auto x34 = x0*x33 + x14*((x33)*(x33)) + x145;
auto x35 = -x10 + x11*x34;
auto x36 = cell[4] - x35;
auto x146 = x17 + x21 + x30 + x36;
auto x37 = V{-2};
auto x38 = -V{1}*x37;
auto x39 = x26 + x5;
auto x40 = x0*x39 + x14*((x39)*(x39)) + x145;
auto x41 = -x10 + x11*x40;
auto x42 = x38*(cell[6] - x41);
auto x43 = x134 + x5;
auto x44 = x0*x43 + x14*((x43)*(x43)) + x145;
auto x45 = -x10 + x11*x44;
auto x46 = cell[16] - x45;
auto x47 = x142 + x146 - x25 - x32 + x42 + V{2}*x46;
auto x48 = x47*x9;
auto x49 = x0*x135 + x139*x14 + x145;
auto x50 = -x22 + x23*x49;
auto x51 = x0*x12 + x14*((x12)*(x12)) + x145;
auto x52 = -x22 + x23*x51;
auto x53 = x12 + x5;
auto x54 = x0*x53 + x14*((x53)*(x53)) + x145;
auto x55 = -x10 + x11*x54;
auto x56 = x38*(cell[8] - x55);
auto x57 = x135 + x5;
auto x58 = x0*x57 + x14*((x57)*(x57)) + x145;
auto x59 = -x10 + x11*x58;
auto x60 = cell[18] - x59;
auto x61 = x143 + x146 - x50 - x52 + x56 + V{2}*x60;
auto x62 = x61*x9;
auto x63 = x48 + x62;
auto x64 = x0*x5 + x14*((x5)*(x5)) + x145;
auto x65 = -x22 + x23*x64;
auto x66 = x38*x46 + x38*x60 + x38*(cell[3] - x65) + x42 + x56;
auto x67 = x66*x9;
auto x68 = x4*(x144*(V{0.5}*x63 + V{0.5}*x67) + x3*x7 + V{-1});
auto x69 = V{1} - x137;
auto x70 = V{0.5}*x10*x144;
auto x71 = V{1} - x8;
auto x72 = x61*x71;
auto x73 = x47*x71;
auto x147 = x67 + x73;
auto x148 = x147 + x72;
auto x74 = -V{1}*x17 + x21 - V{1}*x30 + x36;
auto x75 = x70*(x148 + x37*x74);
auto x76 = x16 + x75;
auto x77 = V{2};
auto x78 = x70*(x148 + x74*x77);
auto x79 = x20 + x78;
auto x80 = V{0.5}*x144*x22;
auto x81 = x48 + x72;
auto x82 = x80*(x67 + x81);
auto x83 = x50 + x82;
auto x84 = x80*(x147 + x62);
auto x85 = x25 + x84;
auto x86 = x29 + x75;
auto x87 = x35 + x78;
auto x88 = x52 + x82;
auto x89 = x32 + x84;
auto x90 = x66*x71;
auto x149 = x81 + x90;
auto x150 = x62 + x73 + x90;
auto x91 = x37*x60 + x56;
auto x92 = x70*(x149 + x37*x91);
auto x93 = x59 + x92;
auto x94 = x37*x46 + x42;
auto x95 = x70*(x150 + x37*x94);
auto x96 = x45 + x95;
auto x97 = x70*(x149 + x77*x91);
auto x98 = x55 + x97;
auto x99 = x80*(x63 + x90);
auto x100 = x65 + x99;
auto x101 = x70*(x150 + x77*x94);
auto x102 = x101 + x41;
auto x103 = x141*(V{2}*x100 + V{2}*x102 + x68 + x76 + x79 + x83 + x85 + x86 + x87 + x88 + x89 + V{2}*x93 + V{2}*x96 + V{2}*x98 + V{1});
auto x104 = x103*x22;
auto x105 = x10*x103;
auto x106 = x136 + x26;
auto x107 = x0*x106 + x14*((x106)*(x106)) + x145;
auto x108 = x12 + x136;
auto x109 = x0*x108 + x14*((x108)*(x108)) + x145;
auto x110 = x0*x136 + x14*x140 + x145;
auto x111 = x134 + x136;
auto x112 = x0*x111 + x14*((x111)*(x111)) + x145;
auto x113 = x135 + x136;
auto x114 = x0*x113 + x14*((x113)*(x113)) + x145;
cell[0] = x137*x4*(x103*x3 + V{-1}) + x68*x69;
cell[10] = x137*(x104*x24 + x151) + x69*x85;
cell[11] = x137*(x104*x49 + x151) + x69*x83;
cell[12] = x137*(x104*x110 + x151) + x69*(x110*x23 + x151 + x99);
cell[13] = -x137*(x10 - x105*x19) + x69*x79;
cell[14] = -x137*(x10 - x105*x15) + x69*x76;
cell[15] = -x137*(x10 - x105*x112) + x69*(x101 + x11*x112 + x152);
cell[16] = -x137*(x10 - x105*x44) + x69*x96;
cell[17] = -x137*(x10 - x105*x114) + x69*(x11*x114 + x152 + x97);
cell[18] = -x137*(x10 - x105*x58) + x69*x93;
cell[1] = x137*(x104*x31 + x151) + x69*x89;
cell[2] = x137*(x104*x51 + x151) + x69*x88;
cell[3] = x100*x69 + x137*(x104*x64 + x151);
cell[4] = -x137*(x10 - x105*x34) + x69*x87;
cell[5] = -x137*(x10 - x105*x28) + x69*x86;
cell[6] = x102*x69 - x137*(x10 - x105*x40);
cell[7] = -x137*(x10 - x105*x107) + x69*(x107*x11 + x152 + x95);
cell[8] = -x137*(x10 - x105*x54) + x69*x98;
cell[9] = -x137*(x10 - x105*x109) + x69*(x109*x11 + x152 + x92);
return { x103, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
