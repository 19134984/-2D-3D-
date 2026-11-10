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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::ParameterFromCell<descriptors::OMEGA, collision::BGK>, forcing::Wagner>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<0, 1>, momenta::RegularizedBoundaryStress<0, 1>, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x99 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x100 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x101 = cell.template getFieldComponent<olb::descriptors::OMEGA>(0);
auto x102 = cell.template getFieldComponent<olb::descriptors::SCALAR>(0);
auto x103 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x105 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1);
auto x104 = V{1} / (x103);
auto x106 = V{2}*cell[6];
auto x107 = cell[4] + cell[8];
auto x108 = ((x105)*(x105));
auto x0 = V{3};
auto x111 = ((x0)*(x0));
auto x118 = x0*x105;
auto x1 = V{0.5}*x0;
auto x2 = x104*(cell[0] + V{2}*cell[5] + V{2}*cell[7] + x106 + x107 + V{1}) + V{-1};
auto x109 = ((x2)*(x2));
auto x3 = x108;
auto x4 = x1*(x109 + x3);
auto x114 = -x4;
auto x115 = x114 + V{1};
auto x5 = V{0.444444444444444};
auto x110 = x103*x5;
auto x6 = x110;
auto x7 = V{1}/x0;
auto x8 = -x7;
auto x9 = -V{1}*x105;
auto x10 = x2 + x9;
auto x11 = x111;
auto x123 = V{1}*x11;
auto x12 = V{0.5}*x11;
auto x13 = V{0.0277777777777778};
auto x112 = x103*x13;
auto x113 = -x13;
auto x14 = x112;
auto x15 = x113 + x14*(x0*x10 + x115 + x12*((x10)*(x10)));
auto x16 = cell[5] - x15;
auto x17 = V{0.111111111111111};
auto x116 = x103*x17;
auto x117 = -x17;
auto x18 = x116;
auto x19 = x117 + x18*(x0*x2 + x109*x12 + x115);
auto x20 = x105 + x2;
auto x21 = x113 + x14*(x0*x20 + x115 + x12*((x20)*(x20)));
auto x22 = V{2}*cell[7] - V{2}*x21;
auto x23 = x106 + V{2}*x16 - V{2}*x19 + x22;
auto x24 = x23*x8;
auto x25 = x0*x9 + x12*((x9)*(x9)) + V{1};
auto x26 = -x17 + x18*(x114 + x25);
auto x27 = x108*x12 + x118 + V{1};
auto x28 = -x17 + x18*(x114 + x27);
auto x29 = V{-2};
auto x119 = x16*x29;
auto x30 = x107 - V{1}*x119 + x22 - x26 - x28;
auto x31 = x30*x8;
auto x32 = V{0.5}*x111*x5*(x24 + x31) - x5 - x6*(x4 + V{-1});
auto x33 = V{1} - x101;
auto x34 = V{0.5}*x111*x17;
auto x35 = V{1} - x7;
auto x36 = x30*x35;
auto x37 = x34*(x24 + x36);
auto x38 = x28 + x37;
auto x39 = x26 + x37;
auto x40 = V{0.5}*x111*x13;
auto x41 = x23*x35;
auto x120 = x36 + x41;
auto x42 = x119 + x22;
auto x43 = x40*(x120 + V{2}*x42);
auto x44 = x21 + x43;
auto x45 = x40*(x120 + x29*x42);
auto x46 = x15 + x45;
auto x47 = x34*(x31 + x41);
auto x48 = x19 + x47;
auto x49 = x104*(x32 + x38 + x39 + V{2}*x44 + V{2}*x46 + V{2}*x48 + V{1}) + V{-1};
auto x121 = ((x49)*(x49));
auto x122 = x0*x49;
auto x50 = x121 + x3;
auto x51 = x1*x50;
auto x126 = -x51;
auto x127 = x126 + V{1};
auto x52 = -x122*x99;
auto x53 = -x100*x118;
auto x54 = V{-1.4999985};
auto x55 = V{0.0833333333333333}*x101*x102*x104;
auto x56 = V{1} - V{0.25}*x101;
auto x57 = x56*x99;
auto x58 = x55 + x57*x99;
auto x59 = x54*x58;
auto x125 = x52 + x59;
auto x60 = x100*x56;
auto x61 = x100*x60 + x55;
auto x62 = x54*x61;
auto x124 = x53 + x62;
auto x63 = -V{1}*x2;
auto x64 = x105 + x63;
auto x65 = -V{1}*x49;
auto x66 = x105 + x65;
auto x67 = x11*x66;
auto x68 = -x0*(x49 + V{1});
auto x128 = -x68;
auto x69 = -x0*(x105 + V{-1});
auto x70 = V{3.0000015};
auto x71 = x58*x70;
auto x131 = x124 + x71;
auto x72 = V{-4.5};
auto x73 = x100*x57;
auto x74 = x72*x73;
auto x75 = x60*x99;
auto x76 = x72*x75;
auto x77 = x61*x70;
auto x129 = x71 + x77;
auto x130 = x129 + x74 + x76;
auto x133 = x125 + x77;
auto x78 = x63 + x9;
auto x79 = x65 + x9;
auto x80 = -x123*x79;
auto x81 = -x0*(x105 + V{1});
auto x82 = V{4.5};
auto x83 = x73*x82;
auto x84 = x75*x82;
auto x132 = x129 + x83 + x84;
auto x85 = x49 + x9;
auto x86 = x11*x85;
auto x87 = -x0*(x49 + V{-1});
auto x88 = x105 + x49;
auto x89 = x11*x88;
cell[0] = -x101*(x5 + x6*(x51 + V{-1})) + x110*(x124 + x125) + x32*x33;
cell[1] = -x101*(x13 - x14*(x0*x66 + x12*((x66)*(x66)) + x127)) + x112*(x100*(x67 + x69) + x130 - x99*(x128 + V{1}*x67)) + x33*(x113 + x14*(x0*x64 + x115 + x12*((x64)*(x64))) + x45);
cell[2] = -x101*(x17 - x18*(x0*x65 + x12*((x65)*(x65)) + x127)) + x116*(x131 - x99*(x123*x65 + x128)) + x33*(x117 + x18*(x0*x63 + x115 + x12*((x63)*(x63))) + x47);
cell[3] = -x101*(x13 - x14*(x0*x79 + x12*((x79)*(x79)) + x127)) + x112*(x100*(x80 + x81) + x132 + x99*(x68 + x80)) + x33*(x113 + x14*(x0*x78 + x115 + x12*((x78)*(x78))) + x43);
cell[4] = -x101*(x17 - x18*(x126 + x25)) + x116*(-x100*(x123*x9 - x81) + x133) + x33*x39;
cell[5] = -x101*(x13 - x14*(x0*x85 + x12*((x85)*(x85)) + x127)) + x112*(x100*(x81 - V{1}*x86) + x130 + x99*(x86 + x87)) + x33*x46;
cell[6] = -x101*(x17 - x18*(x12*x121 + x122 + x127)) + x116*(x131 + x99*(x11*x49 + x87)) + x33*x48;
cell[7] = -x101*(x13 - x14*(x0*x88 + x12*((x88)*(x88)) + x127)) + x112*(x100*(x69 + x89) + x132 + x99*(x87 + x89)) + x33*x44;
cell[8] = -x101*(x17 - x18*(x126 + x27)) + x116*(x100*(x105*x11 + x69) + x133) + x33*x38;
return { x103, x50 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
