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
struct CSE<Dual<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<1, -1>, momenta::DefineSeparately> > >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x86 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x90 = cell.template getFieldComponent<olb::opti::DJDF>(2);
auto x94 = cell.template getFieldComponent<olb::opti::DJDF>(6);
auto x93 = cell.template getFieldComponent<olb::opti::DJDF>(5);
auto x96 = cell.template getFieldComponent<olb::opti::DJDF>(8);
auto x87 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x91 = cell.template getFieldComponent<olb::opti::DJDF>(3);
auto x88 = cell.template getFieldComponent<olb::opti::DJDF>(0);
auto x92 = cell.template getFieldComponent<olb::opti::DJDF>(4);
auto x106 = parameters.template get<descriptors::OMEGA>();
auto x89 = cell.template getFieldComponent<olb::opti::DJDF>(1);
auto x95 = cell.template getFieldComponent<olb::opti::DJDF>(7);
auto x98 = ((x86)*(x86));
auto x99 = ((x87)*(x87));
auto x0 = -V{1}*x87;
auto x1 = V{1}/(x0 + V{1});
auto x2 = V{0.111111111111111};
auto x3 = V{3};
auto x97 = ((x3)*(x3));
auto x4 = V{0.5}*x97;
auto x5 = V{0.5}*x3*(x98 + x99);
auto x100 = V{1} - x5;
auto x6 = x100 + x3*x86 + x4*x98;
auto x111 = x2*x6;
auto x7 = V{1}/x3;
auto x8 = -x7;
auto x9 = V{0.5}*x2*x97;
auto x10 = V{1} - x106;
auto x101 = x10*x9;
auto x11 = x101*(cell[4] + cell[8]);
auto x12 = V{1} - x7;
auto x13 = V{0.0277777777777778};
auto x14 = V{0.5}*x13*x97;
auto x102 = x10*x14;
auto x15 = x102*(cell[3] + cell[7]);
auto x16 = x12*x15;
auto x17 = x101*(cell[2] + cell[6]);
auto x18 = x102*(cell[1] + cell[5]);
auto x19 = x12*x18;
auto x20 = V{0.444444444444444};
auto x21 = V{0.5}*x20*x97;
auto x103 = x21*x8;
auto x22 = cell[0]*x10*x103;
auto x104 = x16 + x19 + x22;
auto x23 = x104 + x11*x8 + x12*x17;
auto x24 = x86 + x87;
auto x25 = x0 + x86;
auto x26 = x100 + x25*x3 + x4*((x25)*(x25));
auto x110 = x13*x26;
auto x112 = V{1}*x110;
auto x27 = x0*x3 + x100 + x4*((x0)*(x0));
auto x109 = x2*x27;
auto x114 = V{1}*x109;
auto x28 = -V{1}*x86;
auto x29 = x0 + x28;
auto x30 = x100 + x29*x3 + x4*((x29)*(x29));
auto x108 = x13*x30;
auto x113 = V{1}*x108;
auto x31 = x100 + x28*x3 + x4*((x28)*(x28));
auto x107 = x2*x31;
auto x32 = x28 + x87;
auto x33 = V{1} - x5;
auto x105 = x20*x33;
auto x34 = x1*(cell[0]*x105 + cell[1]*x13*(x100 + x3*x32 + x4*((x32)*(x32))) + cell[2]*x107 + cell[3]*x108 + cell[4]*x109 + cell[5]*x110 + cell[6]*x111 + cell[7]*x13*(x100 + x24*x3 + x4*((x24)*(x24))) + cell[8]*x2*(x100 + x3*x87 + x4*x99));
auto x35 = V{-2};
auto x36 = -V{1}*x35;
auto x37 = x36*(x104 + x11*x12 + x17*x8);
auto x38 = V{2}*x15 + x18*x35;
auto x39 = V{2}*x23 + x35*x38 + x37;
auto x40 = x23*x36 + x36*x38 + x37;
auto x41 = -V{1}*x23;
auto x42 = -x1*(-x107*x41 - x111*x41 + x112*x39 + x113*x40 + x114*x37);
auto x115 = x34 + x42;
auto x43 = x115 + x23;
auto x44 = V{2}*x34;
auto x45 = V{2}*x42;
auto x116 = x44 + x45;
auto x46 = x116 + x39;
auto x47 = x116 + x37;
auto x48 = x116 + x40;
auto x49 = x115;
auto x50 = x43*x9;
auto x51 = x50*x8;
auto x52 = x14*x46;
auto x53 = x12*x52;
auto x54 = x47*x9;
auto x55 = x14*x48;
auto x56 = x12*x55;
auto x57 = x103*x49;
auto x117 = x53 + x56 + x57;
auto x58 = x36*(x117 + x12*x54 + V{2}*x51);
auto x59 = x35*x52 + V{2}*x55;
auto x60 = x12*x50;
auto x61 = x117 + x54*x8 + V{2}*x60;
auto x62 = x35*x59 + x58 + V{2}*x61;
auto x63 = x36*x59 + x36*x61 + x58;
auto x64 = -V{1}*x61;
auto x65 = x1*(x105*x49 + x107*x43 + x107*x64 + x108*x48 + x109*x47 + x110*x46 + x111*x43 + x111*x64 - x112*x62 - x113*x63 - x114*x58);
auto x66 = x65 - x88;
auto x67 = -x89;
auto x68 = x61 + x65;
auto x69 = x68 - x90;
auto x70 = V{2}*x65;
auto x71 = x63 + x70 - x91;
auto x72 = x58 + x70 - x92;
auto x73 = x62 + x70 - x93;
auto x74 = x68 - x94;
auto x75 = -x95;
auto x76 = -x96;
cell[0] = x66;
cell[1] = x67;
cell[2] = x69;
cell[3] = x71;
cell[4] = x72;
cell[5] = x73;
cell[6] = x74;
cell[7] = x75;
cell[8] = x76;
return { x66 + V{1}, ((x66)*(x66)) + ((x67)*(x67)) + ((x69)*(x69)) + ((x71)*(x71)) + ((x72)*(x72)) + ((x73)*(x73)) + ((x74)*(x74)) + ((x75)*(x75)) + ((x76)*(x76)) };
}
};

}

}

// Generation Info: commit=4ea71c45fdf920fe461e9926f71dc63f7fa1062d
