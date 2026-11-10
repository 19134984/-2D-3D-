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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::ParameterFromCell<descriptors::OMEGA, collision::BGK>, forcing::Wagner>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<0, -1>, momenta::RegularizedBoundaryStress<0, -1>, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x106 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x107 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x108 = cell.template getFieldComponent<olb::descriptors::OMEGA>(0);
auto x109 = cell.template getFieldComponent<olb::descriptors::SCALAR>(0);
auto x110 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x112 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<0, -1>::VELOCITY>(1);
auto x111 = V{1} / (x110);
auto x113 = cell[4] + cell[8];
auto x114 = ((x112)*(x112));
auto x0 = V{3};
auto x117 = ((x0)*(x0));
auto x124 = x0*x112;
auto x1 = V{0.5}*x0;
auto x2 = -V{1}*x111*(cell[0] + V{2}*cell[1] + V{2}*cell[2] + V{2}*cell[3] + x113 + V{1}) + V{1};
auto x115 = ((x2)*(x2));
auto x3 = x114;
auto x4 = x1*(x115 + x3);
auto x120 = -x4;
auto x121 = x120 + V{1};
auto x5 = V{0.444444444444444};
auto x116 = x110*x5;
auto x6 = x116;
auto x7 = V{1}/x0;
auto x8 = -x7;
auto x9 = V{-2};
auto x10 = -V{1}*x9;
auto x11 = -V{1}*x2;
auto x12 = x11 + x112;
auto x13 = x117;
auto x135 = V{1}*x13;
auto x14 = V{0.5}*x13;
auto x15 = V{0.0277777777777778};
auto x118 = x110*x15;
auto x119 = -x15;
auto x16 = x118;
auto x17 = x119 + x16*(x0*x12 + x121 + x14*((x12)*(x12)));
auto x18 = cell[1] - x17;
auto x19 = V{0.111111111111111};
auto x122 = x110*x19;
auto x123 = -x19;
auto x20 = x122;
auto x21 = x123 + x20*(x0*x11 + x121 + x14*((x11)*(x11)));
auto x22 = -V{1}*x112;
auto x23 = x11 + x22;
auto x24 = x119 + x16*(x0*x23 + x121 + x14*((x23)*(x23)));
auto x25 = x10*(cell[3] - x24);
auto x26 = x10*x18 + x10*(cell[2] - x21) + x25;
auto x27 = x26*x8;
auto x28 = x0*x22 + x14*((x22)*(x22)) + V{1};
auto x29 = -x19 + x20*(x120 + x28);
auto x30 = x114*x14 + x124 + V{1};
auto x31 = -x19 + x20*(x120 + x30);
auto x32 = x113 + V{2}*x18 + x25 - x29 - x31;
auto x33 = x32*x8;
auto x34 = V{0.5}*x117*x5*(x27 + x33) - x5 - x6*(x4 + V{-1});
auto x35 = V{1} - x108;
auto x36 = V{0.5}*x117*x19;
auto x37 = V{1} - x7;
auto x38 = x32*x37;
auto x39 = x36*(x27 + x38);
auto x40 = x31 + x39;
auto x41 = x29 + x39;
auto x42 = V{0.5}*x117*x15;
auto x43 = x26*x37;
auto x125 = x38 + x43;
auto x44 = x18*x9 + x25;
auto x45 = x42*(x125 + V{2}*x44);
auto x46 = x24 + x45;
auto x47 = x42*(x125 + x44*x9);
auto x48 = x17 + x47;
auto x49 = x36*(x33 + x43);
auto x50 = x21 + x49;
auto x51 = -V{1}*x111*(x34 + x40 + x41 + V{2}*x46 + V{2}*x48 + V{2}*x50 + V{1}) + V{1};
auto x126 = ((x51)*(x51));
auto x127 = x0*x51;
auto x52 = x126 + x3;
auto x53 = x1*x52;
auto x130 = -x53;
auto x131 = x130 + V{1};
auto x54 = -x106*x127;
auto x55 = -x107*x124;
auto x56 = V{-1.4999985};
auto x57 = V{0.0833333333333333}*x108*x109*x111;
auto x58 = V{1} - V{0.25}*x108;
auto x59 = x106*x58;
auto x60 = x106*x59 + x57;
auto x61 = x56*x60;
auto x129 = x54 + x61;
auto x62 = x107*x58;
auto x63 = x107*x62 + x57;
auto x64 = x56*x63;
auto x128 = x55 + x64;
auto x65 = -x108*(x5 + x6*(x53 + V{-1})) + x116*(x128 + x129) + x34*x35;
auto x66 = -V{1}*x51;
auto x67 = x112 + x66;
auto x68 = x13*x67;
auto x69 = -x0*(x51 + V{1});
auto x132 = -x69;
auto x70 = -x0*(x112 + V{-1});
auto x71 = V{3.0000015};
auto x72 = x60*x71;
auto x136 = x128 + x72;
auto x73 = V{-4.5};
auto x74 = x107*x59;
auto x75 = x73*x74;
auto x76 = x106*x62;
auto x77 = x73*x76;
auto x78 = x63*x71;
auto x133 = x72 + x78;
auto x134 = x133 + x75 + x77;
auto x138 = x129 + x78;
auto x79 = -x108*(x15 - x16*(x0*x67 + x131 + x14*((x67)*(x67)))) + x118*(-x106*(x132 + V{1}*x68) + x107*(x68 + x70) + x134) + x35*x48;
auto x80 = -x108*(x19 - x20*(x0*x66 + x131 + x14*((x66)*(x66)))) + x122*(-x106*(x132 + x135*x66) + x136) + x35*x50;
auto x81 = x22 + x66;
auto x82 = -x135*x81;
auto x83 = -x0*(x112 + V{1});
auto x84 = V{4.5};
auto x85 = x74*x84;
auto x86 = x76*x84;
auto x137 = x133 + x85 + x86;
auto x87 = -x108*(x15 - x16*(x0*x81 + x131 + x14*((x81)*(x81)))) + x118*(x106*(x69 + x82) + x107*(x82 + x83) + x137) + x35*x46;
auto x88 = -x108*(x19 - x20*(x130 + x28)) + x122*(-x107*(x135*x22 - x83) + x138) + x35*x41;
auto x89 = x2 + x22;
auto x90 = x22 + x51;
auto x91 = x13*x90;
auto x92 = -x0*(x51 + V{-1});
auto x93 = x112 + x2;
auto x94 = x112 + x51;
auto x95 = x13*x94;
auto x96 = -x108*(x19 - x20*(x130 + x30)) + x122*(x107*(x112*x13 + x70) + x138) + x35*x40;
cell[0] = x65;
cell[1] = x79;
cell[2] = x80;
cell[3] = x87;
cell[4] = x88;
cell[5] = -x108*(x15 - x16*(x0*x90 + x131 + x14*((x90)*(x90)))) + x118*(x106*(x91 + x92) + x107*(x83 - V{1}*x91) + x134) + x35*(x119 + x16*(x0*x89 + x121 + x14*((x89)*(x89))) + x47);
cell[6] = -x108*(x19 - x20*(x126*x14 + x127 + x131)) + x122*(x106*(x13*x51 + x92) + x136) + x35*(x123 + x20*(x0*x2 + x115*x14 + x121) + x49);
cell[7] = -x108*(x15 - x16*(x0*x94 + x131 + x14*((x94)*(x94)))) + x118*(x106*(x92 + x95) + x107*(x70 + x95) + x137) + x35*(x119 + x16*(x0*x93 + x121 + x14*((x93)*(x93))) + x45);
cell[8] = x96;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x110;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = -V{1}*x111*(x65 + V{2}*x79 + V{2}*x80 + V{2}*x87 + x88 + x96 + V{1}) + V{1};
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x112;
return { x110, x52 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
