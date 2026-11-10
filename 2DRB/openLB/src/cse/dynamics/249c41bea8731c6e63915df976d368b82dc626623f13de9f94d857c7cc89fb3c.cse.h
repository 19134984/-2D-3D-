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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<0, 1>, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x134 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x135 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x136 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x137 = parameters.template get<descriptors::OMEGA>();
auto x138 = ((x134)*(x134));
auto x139 = ((x135)*(x135));
auto x140 = ((x136)*(x136));
auto x142 = V{2}*cell[10];
auto x143 = cell[11] + cell[2];
auto x144 = cell[12] + cell[3];
auto x0 = V{3};
auto x145 = ((x0)*(x0));
auto x1 = x138 + x139 + x140;
auto x2 = V{0.5}*x0*x1;
auto x146 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.333333333333333};
auto x5 = x134 + V{1};
auto x141 = V{1} / (x5);
auto x6 = x141*(cell[0] + V{2}*cell[13] + V{2}*cell[14] + V{2}*cell[15] + V{2}*cell[16] + cell[17] + cell[18] + cell[8] + cell[9] + x142 + x143 + x144 + V{1});
auto x7 = V{1}/x0;
auto x8 = -x7;
auto x9 = V{0.0555555555555556};
auto x147 = -x9;
auto x10 = x6*x9;
auto x11 = V{0.5}*x145;
auto x12 = x0*x134 + x11*x138 + x146;
auto x13 = x10*x12 + x147;
auto x14 = V{0.0277777777777778};
auto x153 = -x14;
auto x15 = x14*x6;
auto x16 = x134 + x135;
auto x17 = x0*x16 + x11*((x16)*(x16)) + x146;
auto x18 = -x14 + x15*x17;
auto x19 = V{2}*cell[13] - V{2}*x18;
auto x20 = -V{1}*x135;
auto x21 = x134 + x20;
auto x22 = x0*x21 + x11*((x21)*(x21)) + x146;
auto x23 = -x14 + x15*x22;
auto x24 = cell[14] - x23;
auto x25 = x134 + x136;
auto x26 = x0*x25 + x11*((x25)*(x25)) + x146;
auto x27 = -x14 + x15*x26;
auto x28 = V{2}*cell[15] - V{2}*x27;
auto x29 = -V{1}*x136;
auto x30 = x134 + x29;
auto x31 = x0*x30 + x11*((x30)*(x30)) + x146;
auto x32 = -x14 + x15*x31;
auto x33 = cell[16] - x32;
auto x34 = -V{2}*x13 + x142 + x19 + V{2}*x24 + x28 + V{2}*x33;
auto x35 = x34*x8;
auto x36 = x135 + x29;
auto x37 = x0*x36 + x11*((x36)*(x36)) + x146;
auto x38 = -x14 + x15*x37;
auto x39 = cell[18] - x38;
auto x40 = x135 + x136;
auto x41 = x0*x40 + x11*((x40)*(x40)) + x146;
auto x42 = -x14 + x15*x41;
auto x43 = cell[17] - x42;
auto x44 = x0*x135 + x11*x139 + x146;
auto x45 = x10*x44 + x147;
auto x46 = x136 + x20;
auto x47 = x0*x46 + x11*((x46)*(x46)) + x146;
auto x48 = -x14 + x15*x47;
auto x49 = cell[9] - x48;
auto x50 = x0*x20 + x11*((x20)*(x20)) + x146;
auto x51 = x10*x50 + x147;
auto x52 = x20 + x29;
auto x53 = x0*x52 + x11*((x52)*(x52)) + x146;
auto x54 = -x14 + x15*x53;
auto x55 = cell[8] - x54;
auto x148 = x39 + x43 + x49 + x55;
auto x56 = V{-2};
auto x57 = -V{1}*x56;
auto x58 = x143 + x148 + x19 + x24*x57 - x45 - x51;
auto x59 = x58*x8;
auto x60 = x35 + x59;
auto x61 = x0*x136 + x11*x140 + x146;
auto x62 = x10*x61 + x147;
auto x63 = x0*x29 + x11*((x29)*(x29)) + x146;
auto x64 = x10*x63 + x147;
auto x65 = x144 + x148 + x28 + x33*x57 - x62 - x64;
auto x66 = x65*x8;
auto x67 = x4*(x145*(V{0.5}*x60 + V{0.5}*x66) + x3*x6 + V{-1});
auto x68 = V{1} - x137;
auto x69 = V{0.5}*x14*x145;
auto x70 = V{1} - x7;
auto x71 = x65*x70;
auto x72 = x58*x70;
auto x73 = x35 + x72;
auto x149 = x71 + x73;
auto x74 = -V{1}*x39 + x43 - V{1}*x49 + x55;
auto x75 = x69*(x149 + x56*x74);
auto x76 = x38 + x75;
auto x77 = V{2};
auto x78 = x69*(x149 + x74*x77);
auto x79 = x42 + x78;
auto x80 = V{0.5}*x145*x9;
auto x81 = x80*(x60 + x71);
auto x82 = x62 + x81;
auto x83 = x80*(x66 + x73);
auto x84 = x45 + x83;
auto x85 = x48 + x75;
auto x86 = x54 + x78;
auto x87 = x64 + x81;
auto x88 = x51 + x83;
auto x89 = x34*x70;
auto x150 = x59 + x89;
auto x151 = x150 + x71;
auto x152 = x66 + x72 + x89;
auto x90 = x28 + x33*x56;
auto x91 = x69*(x151 + x56*x90);
auto x92 = x32 + x91;
auto x93 = x69*(x151 + x77*x90);
auto x94 = x27 + x93;
auto x95 = x19 + x24*x56;
auto x96 = x69*(x152 + x56*x95);
auto x97 = x23 + x96;
auto x98 = x80*(x150 + x66);
auto x99 = x13 + x98;
auto x100 = x69*(x152 + x77*x95);
auto x101 = x100 + x18;
auto x102 = x141*(V{2}*x101 + x67 + x76 + x79 + x82 + x84 + x85 + x86 + x87 + x88 + V{2}*x92 + V{2}*x94 + V{2}*x97 + V{2}*x99 + V{1});
auto x103 = -V{1}*x134;
auto x104 = x0*x103 + x11*((x103)*(x103)) + x146;
auto x105 = x102*x9;
auto x106 = x103 + x20;
auto x107 = x0*x106 + x11*((x106)*(x106)) + x146;
auto x108 = x102*x14;
auto x109 = x103 + x135;
auto x110 = x0*x109 + x11*((x109)*(x109)) + x146;
auto x111 = x103 + x29;
auto x112 = x0*x111 + x11*((x111)*(x111)) + x146;
auto x113 = x103 + x136;
auto x114 = x0*x113 + x11*((x113)*(x113)) + x146;
cell[0] = x137*x4*(x102*x3 + V{-1}) + x67*x68;
cell[10] = x137*(x105*x12 + x147) + x68*x99;
cell[11] = x137*(x105*x44 + x147) + x68*x84;
cell[12] = x137*(x105*x61 + x147) + x68*x82;
cell[13] = x101*x68 + x137*(x108*x17 + x153);
cell[14] = x137*(x108*x22 + x153) + x68*x97;
cell[15] = x137*(x108*x26 + x153) + x68*x94;
cell[16] = x137*(x108*x31 + x153) + x68*x92;
cell[17] = x137*(x108*x41 + x153) + x68*x79;
cell[18] = x137*(x108*x37 + x153) + x68*x76;
cell[1] = x137*(x104*x105 + x147) + x68*(x10*x104 + x147 + x98);
cell[2] = x137*(x105*x50 + x147) + x68*x88;
cell[3] = x137*(x105*x63 + x147) + x68*x87;
cell[4] = x137*(x107*x108 + x153) + x68*(x100 + x107*x15 + x153);
cell[5] = x137*(x108*x110 + x153) + x68*(x110*x15 + x153 + x96);
cell[6] = x137*(x108*x112 + x153) + x68*(x112*x15 + x153 + x93);
cell[7] = x137*(x108*x114 + x153) + x68*(x114*x15 + x153 + x91);
cell[8] = x137*(x108*x53 + x153) + x68*x86;
cell[9] = x137*(x108*x47 + x153) + x68*x85;
return { x102, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
