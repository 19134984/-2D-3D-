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
struct CSE<Dual<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<1, 1>, momenta::DefineSeparately> > >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x87 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x106 = parameters.template get<descriptors::OMEGA>();
auto x96 = cell.template getFieldComponent<olb::opti::DJDF>(8);
auto x90 = cell.template getFieldComponent<olb::opti::DJDF>(2);
auto x91 = cell.template getFieldComponent<olb::opti::DJDF>(3);
auto x86 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x94 = cell.template getFieldComponent<olb::opti::DJDF>(6);
auto x92 = cell.template getFieldComponent<olb::opti::DJDF>(4);
auto x88 = cell.template getFieldComponent<olb::opti::DJDF>(0);
auto x95 = cell.template getFieldComponent<olb::opti::DJDF>(7);
auto x89 = cell.template getFieldComponent<olb::opti::DJDF>(1);
auto x93 = cell.template getFieldComponent<olb::opti::DJDF>(5);
auto x98 = ((x86)*(x86));
auto x99 = ((x87)*(x87));
auto x0 = V{1}/(x87 + V{1});
auto x1 = V{0.111111111111111};
auto x2 = V{3};
auto x97 = ((x2)*(x2));
auto x3 = V{0.5}*x97;
auto x4 = V{0.5}*x2*(x98 + x99);
auto x100 = V{1} - x4;
auto x5 = x100 + x2*x87 + x3*x99;
auto x107 = x1*x5;
auto x112 = V{1}*x107;
auto x6 = V{0.0277777777777778};
auto x7 = x86 + x87;
auto x8 = x100 + x2*x7 + x3*((x7)*(x7));
auto x105 = x6*x8;
auto x114 = V{1}*x105;
auto x9 = x100 + x2*x86 + x3*x98;
auto x104 = x1*x9;
auto x10 = -V{1}*x87;
auto x11 = x10 + x86;
auto x12 = -V{1}*x86;
auto x13 = x10 + x12;
auto x14 = x100 + x12*x2 + x3*((x12)*(x12));
auto x103 = x1*x14;
auto x15 = x12 + x87;
auto x16 = x100 + x15*x2 + x3*((x15)*(x15));
auto x102 = x16*x6;
auto x113 = V{1}*x102;
auto x17 = V{0.444444444444444};
auto x18 = V{1} - x4;
auto x101 = x17*x18;
auto x19 = x0*(cell[0]*x101 + cell[1]*x102 + cell[2]*x103 + cell[3]*x6*(x100 + x13*x2 + x3*((x13)*(x13))) + cell[4]*x1*(x10*x2 + x100 + x3*((x10)*(x10))) + cell[5]*x6*(x100 + x11*x2 + x3*((x11)*(x11))) + cell[6]*x104 + cell[7]*x105 + cell[8]*x107);
auto x20 = V{2}*x19;
auto x21 = V{1}/x2;
auto x22 = V{1} - x21;
auto x23 = V{0.5}*x1*x97;
auto x24 = V{1} - x106;
auto x108 = x23*x24;
auto x25 = x108*(cell[4] + cell[8]);
auto x26 = V{0.5}*x6*x97;
auto x109 = x24*x26;
auto x27 = x109*(cell[3] + cell[7]);
auto x28 = x22*x27;
auto x29 = -x21;
auto x30 = x108*(cell[2] + cell[6]);
auto x31 = x109*(cell[1] + cell[5]);
auto x32 = x22*x31;
auto x33 = V{0.5}*x17*x97;
auto x110 = x29*x33;
auto x34 = cell[0]*x110*x24;
auto x111 = x28 + x32 + x34;
auto x35 = V{2}*x111 + V{2}*x22*x25 + V{2}*x29*x30;
auto x36 = V{-2};
auto x37 = V{2}*x27 + x31*x36;
auto x38 = x111 + x22*x30 + x25*x29;
auto x39 = x35 + V{2}*x37 + V{2}*x38;
auto x40 = -V{1}*x36;
auto x41 = x35 + x36*x37 + x38*x40;
auto x42 = -V{1}*x38;
auto x43 = -x0*(-x103*x42 - x104*x42 + x112*x35 + x113*x41 + x114*x39);
auto x116 = x19 + x43;
auto x44 = V{2}*x43;
auto x115 = x20 + x44;
auto x45 = x115 + x35;
auto x46 = x115 + x39;
auto x47 = x116 + x38;
auto x48 = x115 + x41;
auto x49 = x116;
auto x50 = x23*x47;
auto x51 = x29*x50;
auto x52 = x23*x45;
auto x53 = x26*x46;
auto x54 = x22*x53;
auto x55 = x26*x48;
auto x56 = x22*x55;
auto x57 = x110*x49;
auto x58 = V{2}*x22*x52 + V{4}*x51 + V{2}*x54 + V{2}*x56 + V{2}*x57;
auto x59 = x36*x55 + V{2}*x53;
auto x60 = x22*x50;
auto x61 = x29*x52 + x54 + x56 + x57 + V{2}*x60;
auto x62 = x58 + V{2}*x59 + V{2}*x61;
auto x63 = x36*x59 + x40*x61 + x58;
auto x64 = -V{1}*x61;
auto x65 = x0*(x101*x49 + x102*x48 + x103*x47 + x103*x64 + x104*x47 + x104*x64 + x105*x46 + x107*x45 - x112*x58 - x113*x63 - x114*x62);
auto x66 = x65 - x88;
auto x67 = V{2}*x65;
auto x68 = x63 + x67 - x89;
auto x69 = x61 + x65;
auto x70 = x69 - x90;
auto x71 = -x91;
auto x72 = -x92;
auto x73 = -x93;
auto x74 = x69 - x94;
auto x75 = x62 + x67 - x95;
auto x76 = x58 + x67 - x96;
cell[0] = x66;
cell[1] = x68;
cell[2] = x70;
cell[3] = x71;
cell[4] = x72;
cell[5] = x73;
cell[6] = x74;
cell[7] = x75;
cell[8] = x76;
return { x66 + V{1}, ((x66)*(x66)) + ((x68)*(x68)) + ((x70)*(x70)) + ((x71)*(x71)) + ((x72)*(x72)) + ((x73)*(x73)) + ((x74)*(x74)) + ((x75)*(x75)) + ((x76)*(x76)) };
}
};

}

}

// Generation Info: commit=4ea71c45fdf920fe461e9926f71dc63f7fa1062d
