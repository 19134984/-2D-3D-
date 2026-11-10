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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<1, 1>, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x134 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x135 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x136 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x137 = parameters.template get<descriptors::OMEGA>();
auto x138 = ((x134)*(x134));
auto x139 = ((x135)*(x135));
auto x140 = ((x136)*(x136));
auto x142 = V{2}*cell[11];
auto x143 = cell[10] + cell[1];
auto x144 = cell[12] + cell[3];
auto x0 = V{3};
auto x145 = ((x0)*(x0));
auto x1 = x138 + x139 + x140;
auto x2 = V{0.5}*x0*x1;
auto x146 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = x135 + V{1};
auto x141 = V{1} / (x5);
auto x6 = x141*(cell[0] + V{2}*cell[13] + cell[15] + cell[16] + V{2}*cell[17] + V{2}*cell[18] + V{2}*cell[5] + cell[6] + cell[7] + x142 + x143 + x144 + V{1});
auto x7 = V{1}/x0;
auto x8 = -x7;
auto x9 = V{0.0277777777777778};
auto x147 = -x9;
auto x10 = x6*x9;
auto x11 = -V{1}*x136;
auto x12 = x11 + x134;
auto x13 = V{0.5}*x145;
auto x14 = x0*x12 + x13*((x12)*(x12)) + x146;
auto x15 = x10*x14 + x147;
auto x16 = cell[16] - x15;
auto x17 = x134 + x136;
auto x18 = x0*x17 + x13*((x17)*(x17)) + x146;
auto x19 = x10*x18 + x147;
auto x20 = cell[15] - x19;
auto x21 = V{0.0555555555555556};
auto x153 = -x21;
auto x22 = x21*x6;
auto x23 = x0*x134 + x13*x138 + x146;
auto x24 = -x21 + x22*x23;
auto x25 = -V{1}*x134;
auto x26 = x136 + x25;
auto x27 = x0*x26 + x13*((x26)*(x26)) + x146;
auto x28 = x10*x27 + x147;
auto x29 = cell[7] - x28;
auto x30 = x0*x25 + x13*((x25)*(x25)) + x146;
auto x31 = -x21 + x22*x30;
auto x32 = x11 + x25;
auto x33 = x0*x32 + x13*((x32)*(x32)) + x146;
auto x34 = x10*x33 + x147;
auto x35 = cell[6] - x34;
auto x148 = x16 + x20 + x29 + x35;
auto x36 = V{-2};
auto x37 = -V{1}*x36;
auto x38 = x135 + x25;
auto x39 = x0*x38 + x13*((x38)*(x38)) + x146;
auto x40 = x10*x39 + x147;
auto x41 = cell[5] - x40;
auto x42 = x134 + x135;
auto x43 = x0*x42 + x13*((x42)*(x42)) + x146;
auto x44 = x10*x43 + x147;
auto x45 = V{2}*cell[13] - V{2}*x44;
auto x46 = x143 + x148 - x24 - x31 + x37*x41 + x45;
auto x47 = x46*x8;
auto x48 = x0*x135 + x13*x139 + x146;
auto x49 = -x21 + x22*x48;
auto x50 = x135 + x136;
auto x51 = x0*x50 + x13*((x50)*(x50)) + x146;
auto x52 = x10*x51 + x147;
auto x53 = V{2}*cell[17] - V{2}*x52;
auto x54 = x11 + x135;
auto x55 = x0*x54 + x13*((x54)*(x54)) + x146;
auto x56 = x10*x55 + x147;
auto x57 = cell[18] - x56;
auto x58 = x142 + V{2}*x41 + x45 - V{2}*x49 + x53 + V{2}*x57;
auto x59 = x58*x8;
auto x60 = x47 + x59;
auto x61 = x0*x136 + x13*x140 + x146;
auto x62 = -x21 + x22*x61;
auto x63 = x0*x11 + x13*((x11)*(x11)) + x146;
auto x64 = -x21 + x22*x63;
auto x65 = x144 + x148 + x37*x57 + x53 - x62 - x64;
auto x66 = x65*x8;
auto x67 = x4*(x145*(V{0.5}*x60 + V{0.5}*x66) + x3*x6 + V{-1});
auto x68 = V{1} - x137;
auto x69 = V{0.5}*x145*x9;
auto x70 = V{1} - x7;
auto x71 = x65*x70;
auto x72 = x46*x70;
auto x149 = x59 + x72;
auto x150 = x149 + x71;
auto x73 = -V{1}*x16 + x20 - V{1}*x29 + x35;
auto x74 = x69*(x150 + x36*x73);
auto x75 = x15 + x74;
auto x76 = V{2};
auto x77 = x69*(x150 + x73*x76);
auto x78 = x19 + x77;
auto x79 = V{0.5}*x145*x21;
auto x80 = x79*(x60 + x71);
auto x81 = x62 + x80;
auto x82 = x79*(x149 + x66);
auto x83 = x24 + x82;
auto x84 = x28 + x74;
auto x85 = x34 + x77;
auto x86 = x64 + x80;
auto x87 = x31 + x82;
auto x88 = x58*x70;
auto x152 = x66 + x72 + x88;
auto x89 = x47 + x88;
auto x151 = x71 + x89;
auto x90 = x36*x57 + x53;
auto x91 = x69*(x151 + x36*x90);
auto x92 = x56 + x91;
auto x93 = x69*(x151 + x76*x90);
auto x94 = x52 + x93;
auto x95 = x36*x41 + x45;
auto x96 = x69*(x152 + x76*x95);
auto x97 = x44 + x96;
auto x98 = x69*(x152 + x36*x95);
auto x99 = x40 + x98;
auto x100 = x79*(x66 + x89);
auto x101 = x100 + x49;
auto x102 = x141*(V{2}*x101 + x67 + x75 + x78 + x81 + x83 + x84 + x85 + x86 + x87 + V{2}*x92 + V{2}*x94 + V{2}*x97 + V{2}*x99 + V{1});
auto x103 = x102*x21;
auto x104 = -V{1}*x135;
auto x105 = x0*x104 + x13*((x104)*(x104)) + x146;
auto x106 = x104 + x25;
auto x107 = x0*x106 + x13*((x106)*(x106)) + x146;
auto x108 = x102*x9;
auto x109 = x104 + x11;
auto x110 = x0*x109 + x13*((x109)*(x109)) + x146;
auto x111 = x104 + x136;
auto x112 = x0*x111 + x13*((x111)*(x111)) + x146;
auto x113 = x104 + x134;
auto x114 = x0*x113 + x13*((x113)*(x113)) + x146;
cell[0] = x137*x4*(x102*x3 + V{-1}) + x67*x68;
cell[10] = x137*(x103*x23 + x153) + x68*x83;
cell[11] = x101*x68 + x137*(x103*x48 + x153);
cell[12] = x137*(x103*x61 + x153) + x68*x81;
cell[13] = x137*(x108*x43 + x147) + x68*x97;
cell[14] = x137*(x108*x114 + x147) + x68*(x10*x114 + x147 + x98);
cell[15] = x137*(x108*x18 + x147) + x68*x78;
cell[16] = x137*(x108*x14 + x147) + x68*x75;
cell[17] = x137*(x108*x51 + x147) + x68*x94;
cell[18] = x137*(x108*x55 + x147) + x68*x92;
cell[1] = x137*(x103*x30 + x153) + x68*x87;
cell[2] = x137*(x103*x105 + x153) + x68*(x100 + x105*x22 + x153);
cell[3] = x137*(x103*x63 + x153) + x68*x86;
cell[4] = x137*(x107*x108 + x147) + x68*(x10*x107 + x147 + x96);
cell[5] = x137*(x108*x39 + x147) + x68*x99;
cell[6] = x137*(x108*x33 + x147) + x68*x85;
cell[7] = x137*(x108*x27 + x147) + x68*x84;
cell[8] = x137*(x108*x110 + x147) + x68*(x10*x110 + x147 + x93);
cell[9] = x137*(x108*x112 + x147) + x68*(x10*x112 + x147 + x91);
return { x102, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
