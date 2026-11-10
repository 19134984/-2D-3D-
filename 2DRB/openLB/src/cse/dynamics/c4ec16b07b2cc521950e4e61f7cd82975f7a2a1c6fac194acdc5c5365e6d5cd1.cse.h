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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<2, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<2, 1>, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x134 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x135 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x136 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x137 = parameters.template get<descriptors::OMEGA>();
auto x138 = ((x134)*(x134));
auto x139 = ((x135)*(x135));
auto x140 = ((x136)*(x136));
auto x142 = V{2}*cell[12];
auto x143 = cell[10] + cell[1];
auto x144 = cell[11] + cell[2];
auto x0 = V{3};
auto x145 = ((x0)*(x0));
auto x1 = x138 + x139 + x140;
auto x2 = V{0.5}*x0*x1;
auto x146 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = x136 + V{1};
auto x141 = V{1} / (x5);
auto x6 = x141*(cell[0] + cell[13] + cell[14] + V{2}*cell[15] + V{2}*cell[17] + cell[4] + cell[5] + V{2}*cell[7] + V{2}*cell[9] + x142 + x143 + x144 + V{1});
auto x7 = V{1}/x0;
auto x8 = -x7;
auto x9 = V{0.0277777777777778};
auto x147 = -x9;
auto x10 = x6*x9;
auto x11 = -V{1}*x135;
auto x12 = x11 + x134;
auto x13 = V{0.5}*x145;
auto x14 = x0*x12 + x13*((x12)*(x12)) + x146;
auto x15 = x10*x14 + x147;
auto x16 = cell[14] - x15;
auto x17 = x134 + x135;
auto x18 = x0*x17 + x13*((x17)*(x17)) + x146;
auto x19 = x10*x18 + x147;
auto x20 = cell[13] - x19;
auto x21 = V{0.0555555555555556};
auto x153 = -x21;
auto x22 = x21*x6;
auto x23 = x0*x134 + x13*x138 + x146;
auto x24 = -x21 + x22*x23;
auto x25 = -V{1}*x134;
auto x26 = x135 + x25;
auto x27 = x0*x26 + x13*((x26)*(x26)) + x146;
auto x28 = x10*x27 + x147;
auto x29 = cell[5] - x28;
auto x30 = x0*x25 + x13*((x25)*(x25)) + x146;
auto x31 = -x21 + x22*x30;
auto x32 = x11 + x25;
auto x33 = x0*x32 + x13*((x32)*(x32)) + x146;
auto x34 = x10*x33 + x147;
auto x35 = cell[4] - x34;
auto x148 = x16 + x20 + x29 + x35;
auto x36 = V{-2};
auto x37 = -V{1}*x36;
auto x38 = x136 + x25;
auto x39 = x0*x38 + x13*((x38)*(x38)) + x146;
auto x40 = x10*x39 + x147;
auto x41 = cell[7] - x40;
auto x42 = x134 + x136;
auto x43 = x0*x42 + x13*((x42)*(x42)) + x146;
auto x44 = x10*x43 + x147;
auto x45 = V{2}*cell[15] - V{2}*x44;
auto x46 = x143 + x148 - x24 - x31 + x37*x41 + x45;
auto x47 = x46*x8;
auto x48 = x0*x135 + x13*x139 + x146;
auto x49 = -x21 + x22*x48;
auto x50 = x0*x11 + x13*((x11)*(x11)) + x146;
auto x51 = -x21 + x22*x50;
auto x52 = x11 + x136;
auto x53 = x0*x52 + x13*((x52)*(x52)) + x146;
auto x54 = x10*x53 + x147;
auto x55 = cell[9] - x54;
auto x56 = x135 + x136;
auto x57 = x0*x56 + x13*((x56)*(x56)) + x146;
auto x58 = x10*x57 + x147;
auto x59 = V{2}*cell[17] - V{2}*x58;
auto x60 = x144 + x148 + x37*x55 - x49 - x51 + x59;
auto x61 = x60*x8;
auto x62 = x47 + x61;
auto x63 = x0*x136 + x13*x140 + x146;
auto x64 = -x21 + x22*x63;
auto x65 = x142 + V{2}*x41 + x45 + V{2}*x55 + x59 - V{2}*x64;
auto x66 = x65*x8;
auto x67 = x4*(x145*(V{0.5}*x62 + V{0.5}*x66) + x3*x6 + V{-1});
auto x68 = V{1} - x137;
auto x69 = V{0.5}*x145*x9;
auto x70 = V{1} - x7;
auto x71 = x60*x70;
auto x72 = x46*x70;
auto x149 = x66 + x72;
auto x150 = x149 + x71;
auto x73 = -V{1}*x16 + x20 - V{1}*x29 + x35;
auto x74 = x69*(x150 + x36*x73);
auto x75 = x15 + x74;
auto x76 = V{2};
auto x77 = x69*(x150 + x73*x76);
auto x78 = x19 + x77;
auto x79 = V{0.5}*x145*x21;
auto x80 = x47 + x71;
auto x81 = x79*(x66 + x80);
auto x82 = x49 + x81;
auto x83 = x79*(x149 + x61);
auto x84 = x24 + x83;
auto x85 = x28 + x74;
auto x86 = x34 + x77;
auto x87 = x51 + x81;
auto x88 = x31 + x83;
auto x89 = x65*x70;
auto x151 = x80 + x89;
auto x152 = x61 + x72 + x89;
auto x90 = x36*x55 + x59;
auto x91 = x69*(x151 + x76*x90);
auto x92 = x58 + x91;
auto x93 = x36*x41 + x45;
auto x94 = x69*(x152 + x76*x93);
auto x95 = x44 + x94;
auto x96 = x79*(x62 + x89);
auto x97 = x64 + x96;
auto x98 = x69*(x152 + x36*x93);
auto x99 = x40 + x98;
auto x100 = x69*(x151 + x36*x90);
auto x101 = x100 + x54;
auto x102 = x141*(V{2}*x101 + x67 + x75 + x78 + x82 + x84 + x85 + x86 + x87 + x88 + V{2}*x92 + V{2}*x95 + V{2}*x97 + V{2}*x99 + V{1});
auto x103 = x102*x21;
auto x104 = -V{1}*x136;
auto x105 = x0*x104 + x13*((x104)*(x104)) + x146;
auto x106 = x102*x9;
auto x107 = x104 + x25;
auto x108 = x0*x107 + x13*((x107)*(x107)) + x146;
auto x109 = x104 + x11;
auto x110 = x0*x109 + x13*((x109)*(x109)) + x146;
auto x111 = x104 + x134;
auto x112 = x0*x111 + x13*((x111)*(x111)) + x146;
auto x113 = x104 + x135;
auto x114 = x0*x113 + x13*((x113)*(x113)) + x146;
cell[0] = x137*x4*(x102*x3 + V{-1}) + x67*x68;
cell[10] = x137*(x103*x23 + x153) + x68*x84;
cell[11] = x137*(x103*x48 + x153) + x68*x82;
cell[12] = x137*(x103*x63 + x153) + x68*x97;
cell[13] = x137*(x106*x18 + x147) + x68*x78;
cell[14] = x137*(x106*x14 + x147) + x68*x75;
cell[15] = x137*(x106*x43 + x147) + x68*x95;
cell[16] = x137*(x106*x112 + x147) + x68*(x10*x112 + x147 + x98);
cell[17] = x137*(x106*x57 + x147) + x68*x92;
cell[18] = x137*(x106*x114 + x147) + x68*(x10*x114 + x100 + x147);
cell[1] = x137*(x103*x30 + x153) + x68*x88;
cell[2] = x137*(x103*x50 + x153) + x68*x87;
cell[3] = x137*(x103*x105 + x153) + x68*(x105*x22 + x153 + x96);
cell[4] = x137*(x106*x33 + x147) + x68*x86;
cell[5] = x137*(x106*x27 + x147) + x68*x85;
cell[6] = x137*(x106*x108 + x147) + x68*(x10*x108 + x147 + x94);
cell[7] = x137*(x106*x39 + x147) + x68*x99;
cell[8] = x137*(x106*x110 + x147) + x68*(x10*x110 + x147 + x91);
cell[9] = x101*x68 + x137*(x106*x53 + x147);
return { x102, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
