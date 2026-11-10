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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::ParameterFromCell<descriptors::OMEGA, collision::BGK>, forcing::Wagner>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<0, -1>, momenta::RegularizedBoundaryStress<0, -1>, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x100 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x101 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x102 = cell.template getFieldComponent<olb::descriptors::OMEGA>(0);
auto x103 = cell.template getFieldComponent<olb::descriptors::SCALAR>(0);
auto x104 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x106 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<0, -1>::VELOCITY>(1);
auto x105 = V{1} / (x104);
auto x107 = cell[4] + cell[8];
auto x108 = ((x106)*(x106));
auto x0 = V{3};
auto x111 = ((x0)*(x0));
auto x118 = x0*x106;
auto x1 = V{0.5}*x0;
auto x2 = -V{1}*x105*(cell[0] + V{2}*cell[1] + V{2}*cell[2] + V{2}*cell[3] + x107 + V{1}) + V{1};
auto x109 = ((x2)*(x2));
auto x3 = x108;
auto x4 = x1*(x109 + x3);
auto x114 = -x4;
auto x115 = x114 + V{1};
auto x5 = V{0.444444444444444};
auto x110 = x104*x5;
auto x6 = x110;
auto x7 = V{1}/x0;
auto x8 = -x7;
auto x9 = V{-2};
auto x10 = -V{1}*x9;
auto x11 = -V{1}*x2;
auto x12 = x106 + x11;
auto x13 = x111;
auto x122 = V{1}*x13;
auto x14 = V{0.5}*x13;
auto x15 = V{0.0277777777777778};
auto x112 = x104*x15;
auto x113 = -x15;
auto x16 = x112;
auto x17 = x113 + x16*(x0*x12 + x115 + x14*((x12)*(x12)));
auto x18 = cell[1] - x17;
auto x19 = V{0.111111111111111};
auto x116 = x104*x19;
auto x117 = -x19;
auto x20 = x116;
auto x21 = x117 + x20*(x0*x11 + x115 + x14*((x11)*(x11)));
auto x22 = -V{1}*x106;
auto x23 = x11 + x22;
auto x24 = x113 + x16*(x0*x23 + x115 + x14*((x23)*(x23)));
auto x25 = x10*(cell[3] - x24);
auto x26 = x10*x18 + x10*(cell[2] - x21) + x25;
auto x27 = x26*x8;
auto x28 = x0*x22 + x14*((x22)*(x22)) + V{1};
auto x29 = -x19 + x20*(x114 + x28);
auto x30 = x108*x14 + x118 + V{1};
auto x31 = -x19 + x20*(x114 + x30);
auto x32 = x107 + V{2}*x18 + x25 - x29 - x31;
auto x33 = x32*x8;
auto x34 = V{0.5}*x111*x5*(x27 + x33) - x5 - x6*(x4 + V{-1});
auto x35 = V{1} - x102;
auto x36 = V{0.5}*x111*x19;
auto x37 = V{1} - x7;
auto x38 = x32*x37;
auto x39 = x36*(x27 + x38);
auto x40 = x31 + x39;
auto x41 = x29 + x39;
auto x42 = V{0.5}*x111*x15;
auto x43 = x26*x37;
auto x119 = x38 + x43;
auto x44 = x18*x9 + x25;
auto x45 = x42*(x119 + V{2}*x44);
auto x46 = x24 + x45;
auto x47 = x42*(x119 + x44*x9);
auto x48 = x17 + x47;
auto x49 = x36*(x33 + x43);
auto x50 = x21 + x49;
auto x51 = -V{1}*x105*(x34 + x40 + x41 + V{2}*x46 + V{2}*x48 + V{2}*x50 + V{1}) + V{1};
auto x120 = ((x51)*(x51));
auto x121 = x0*x51;
auto x52 = x120 + x3;
auto x53 = x1*x52;
auto x125 = -x53;
auto x126 = x125 + V{1};
auto x54 = -x100*x121;
auto x55 = -x101*x118;
auto x56 = V{-1.4999985};
auto x57 = V{0.0833333333333333}*x102*x103*x105;
auto x58 = V{1} - V{0.25}*x102;
auto x59 = x100*x58;
auto x60 = x100*x59 + x57;
auto x61 = x56*x60;
auto x124 = x54 + x61;
auto x62 = x101*x58;
auto x63 = x101*x62 + x57;
auto x64 = x56*x63;
auto x123 = x55 + x64;
auto x65 = -V{1}*x51;
auto x66 = x106 + x65;
auto x67 = x13*x66;
auto x68 = -x0*(x51 + V{1});
auto x127 = -x68;
auto x69 = -x0*(x106 + V{-1});
auto x70 = V{3.0000015};
auto x71 = x60*x70;
auto x130 = x123 + x71;
auto x72 = V{-4.5};
auto x73 = x101*x59;
auto x74 = x72*x73;
auto x75 = x100*x62;
auto x76 = x72*x75;
auto x77 = x63*x70;
auto x128 = x71 + x77;
auto x129 = x128 + x74 + x76;
auto x132 = x124 + x77;
auto x78 = x22 + x65;
auto x79 = -x122*x78;
auto x80 = -x0*(x106 + V{1});
auto x81 = V{4.5};
auto x82 = x73*x81;
auto x83 = x75*x81;
auto x131 = x128 + x82 + x83;
auto x84 = x2 + x22;
auto x85 = x22 + x51;
auto x86 = x13*x85;
auto x87 = -x0*(x51 + V{-1});
auto x88 = x106 + x2;
auto x89 = x106 + x51;
auto x90 = x13*x89;
cell[0] = -x102*(x5 + x6*(x53 + V{-1})) + x110*(x123 + x124) + x34*x35;
cell[1] = -x102*(x15 - x16*(x0*x66 + x126 + x14*((x66)*(x66)))) + x112*(-x100*(x127 + V{1}*x67) + x101*(x67 + x69) + x129) + x35*x48;
cell[2] = -x102*(x19 - x20*(x0*x65 + x126 + x14*((x65)*(x65)))) + x116*(-x100*(x122*x65 + x127) + x130) + x35*x50;
cell[3] = -x102*(x15 - x16*(x0*x78 + x126 + x14*((x78)*(x78)))) + x112*(x100*(x68 + x79) + x101*(x79 + x80) + x131) + x35*x46;
cell[4] = -x102*(x19 - x20*(x125 + x28)) + x116*(-x101*(x122*x22 - x80) + x132) + x35*x41;
cell[5] = -x102*(x15 - x16*(x0*x85 + x126 + x14*((x85)*(x85)))) + x112*(x100*(x86 + x87) + x101*(x80 - V{1}*x86) + x129) + x35*(x113 + x16*(x0*x84 + x115 + x14*((x84)*(x84))) + x47);
cell[6] = -x102*(x19 - x20*(x120*x14 + x121 + x126)) + x116*(x100*(x13*x51 + x87) + x130) + x35*(x117 + x20*(x0*x2 + x109*x14 + x115) + x49);
cell[7] = -x102*(x15 - x16*(x0*x89 + x126 + x14*((x89)*(x89)))) + x112*(x100*(x87 + x90) + x101*(x69 + x90) + x131) + x35*(x113 + x16*(x0*x88 + x115 + x14*((x88)*(x88))) + x45);
cell[8] = -x102*(x19 - x20*(x125 + x30)) + x116*(x101*(x106*x13 + x69) + x132) + x35*x40;
return { x104, x52 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
