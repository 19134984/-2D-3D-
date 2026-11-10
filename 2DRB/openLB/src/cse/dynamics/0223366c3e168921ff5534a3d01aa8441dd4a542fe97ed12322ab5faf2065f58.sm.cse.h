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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::ParameterFromCell<descriptors::OMEGA, collision::BGK>, forcing::Wagner>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<0, 1>, momenta::RegularizedBoundaryStress<0, 1>, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x105 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x106 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x107 = cell.template getFieldComponent<olb::descriptors::OMEGA>(0);
auto x108 = cell.template getFieldComponent<olb::descriptors::SCALAR>(0);
auto x109 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x111 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1);
auto x110 = V{1} / (x109);
auto x112 = V{2}*cell[6];
auto x113 = cell[4] + cell[8];
auto x114 = ((x111)*(x111));
auto x0 = V{3};
auto x117 = ((x0)*(x0));
auto x124 = x0*x111;
auto x1 = V{0.5}*x0;
auto x2 = x110*(cell[0] + V{2}*cell[5] + V{2}*cell[7] + x112 + x113 + V{1}) + V{-1};
auto x115 = ((x2)*(x2));
auto x3 = x114;
auto x4 = x1*(x115 + x3);
auto x120 = -x4;
auto x121 = x120 + V{1};
auto x5 = V{0.444444444444444};
auto x116 = x109*x5;
auto x6 = x116;
auto x7 = V{1}/x0;
auto x8 = -x7;
auto x9 = -V{1}*x111;
auto x10 = x2 + x9;
auto x11 = x117;
auto x131 = V{1}*x11;
auto x12 = V{0.5}*x11;
auto x13 = V{0.0277777777777778};
auto x118 = x109*x13;
auto x119 = -x13;
auto x14 = x118;
auto x15 = x119 + x14*(x0*x10 + x12*((x10)*(x10)) + x121);
auto x16 = cell[5] - x15;
auto x17 = V{0.111111111111111};
auto x122 = x109*x17;
auto x123 = -x17;
auto x18 = x122;
auto x19 = x123 + x18*(x0*x2 + x115*x12 + x121);
auto x20 = x111 + x2;
auto x21 = x119 + x14*(x0*x20 + x12*((x20)*(x20)) + x121);
auto x22 = V{2}*cell[7] - V{2}*x21;
auto x23 = x112 + V{2}*x16 - V{2}*x19 + x22;
auto x24 = x23*x8;
auto x25 = x0*x9 + x12*((x9)*(x9)) + V{1};
auto x26 = -x17 + x18*(x120 + x25);
auto x27 = x114*x12 + x124 + V{1};
auto x28 = -x17 + x18*(x120 + x27);
auto x29 = V{-2};
auto x125 = x16*x29;
auto x30 = x113 - V{1}*x125 + x22 - x26 - x28;
auto x31 = x30*x8;
auto x32 = V{0.5}*x117*x5*(x24 + x31) - x5 - x6*(x4 + V{-1});
auto x33 = V{1} - x107;
auto x34 = V{0.5}*x117*x17;
auto x35 = V{1} - x7;
auto x36 = x30*x35;
auto x37 = x34*(x24 + x36);
auto x38 = x28 + x37;
auto x39 = x26 + x37;
auto x40 = V{0.5}*x117*x13;
auto x41 = x23*x35;
auto x126 = x36 + x41;
auto x42 = x125 + x22;
auto x43 = x40*(x126 + V{2}*x42);
auto x44 = x21 + x43;
auto x45 = x40*(x126 + x29*x42);
auto x46 = x15 + x45;
auto x47 = x34*(x31 + x41);
auto x48 = x19 + x47;
auto x49 = x110*(x32 + x38 + x39 + V{2}*x44 + V{2}*x46 + V{2}*x48 + V{1}) + V{-1};
auto x127 = ((x49)*(x49));
auto x128 = x0*x49;
auto x50 = x127 + x3;
auto x51 = x1*x50;
auto x132 = -x51;
auto x136 = x132 + V{1};
auto x52 = -x105*x128;
auto x53 = -x106*x124;
auto x54 = V{-1.4999985};
auto x55 = V{0.0833333333333333}*x107*x108*x110;
auto x56 = V{1} - V{0.25}*x107;
auto x57 = x105*x56;
auto x58 = x105*x57 + x55;
auto x59 = x54*x58;
auto x129 = x52 + x59;
auto x60 = x106*x56;
auto x61 = x106*x60 + x55;
auto x62 = x54*x61;
auto x130 = x53 + x62;
auto x63 = -x107*(x5 + x6*(x51 + V{-1})) + x116*(x129 + x130) + x32*x33;
auto x64 = -V{1}*x2;
auto x65 = x111 + x64;
auto x66 = -V{1}*x49;
auto x67 = x111 + x66;
auto x68 = x11*x67;
auto x69 = -x0*(x49 + V{1});
auto x139 = -x69;
auto x70 = -x0*(x111 + V{-1});
auto x71 = V{3.0000015};
auto x72 = x58*x71;
auto x137 = x130 + x72;
auto x73 = V{-4.5};
auto x74 = x106*x57;
auto x75 = x73*x74;
auto x76 = x105*x60;
auto x77 = x73*x76;
auto x78 = x61*x71;
auto x133 = x129 + x78;
auto x134 = x72 + x78;
auto x135 = x134 + x75 + x77;
auto x79 = x64 + x9;
auto x80 = x66 + x9;
auto x81 = -x131*x80;
auto x82 = -x0*(x111 + V{1});
auto x83 = V{4.5};
auto x84 = x74*x83;
auto x85 = x76*x83;
auto x138 = x134 + x84 + x85;
auto x86 = -x107*(x17 - x18*(x132 + x25)) + x122*(-x106*(x131*x9 - x82) + x133) + x33*x39;
auto x87 = x49 + x9;
auto x88 = x11*x87;
auto x89 = -x0*(x49 + V{-1});
auto x90 = -x107*(x13 - x14*(x0*x87 + x12*((x87)*(x87)) + x136)) + x118*(x105*(x88 + x89) + x106*(x82 - V{1}*x88) + x135) + x33*x46;
auto x91 = -x107*(x17 - x18*(x12*x127 + x128 + x136)) + x122*(x105*(x11*x49 + x89) + x137) + x33*x48;
auto x92 = x111 + x49;
auto x93 = x11*x92;
auto x94 = -x107*(x13 - x14*(x0*x92 + x12*((x92)*(x92)) + x136)) + x118*(x105*(x89 + x93) + x106*(x70 + x93) + x138) + x33*x44;
auto x95 = -x107*(x17 - x18*(x132 + x27)) + x122*(x106*(x11*x111 + x70) + x133) + x33*x38;
cell[0] = x63;
cell[1] = -x107*(x13 - x14*(x0*x67 + x12*((x67)*(x67)) + x136)) + x118*(-x105*(x139 + V{1}*x68) + x106*(x68 + x70) + x135) + x33*(x119 + x14*(x0*x65 + x12*((x65)*(x65)) + x121) + x45);
cell[2] = -x107*(x17 - x18*(x0*x66 + x12*((x66)*(x66)) + x136)) + x122*(-x105*(x131*x66 + x139) + x137) + x33*(x123 + x18*(x0*x64 + x12*((x64)*(x64)) + x121) + x47);
cell[3] = -x107*(x13 - x14*(x0*x80 + x12*((x80)*(x80)) + x136)) + x118*(x105*(x69 + x81) + x106*(x81 + x82) + x138) + x33*(x119 + x14*(x0*x79 + x12*((x79)*(x79)) + x121) + x43);
cell[4] = x86;
cell[5] = x90;
cell[6] = x91;
cell[7] = x94;
cell[8] = x95;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x109;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x110*(x63 + x86 + V{2}*x90 + V{2}*x91 + V{2}*x94 + x95 + V{1}) + V{-1};
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x111;
return { x109, x50 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
