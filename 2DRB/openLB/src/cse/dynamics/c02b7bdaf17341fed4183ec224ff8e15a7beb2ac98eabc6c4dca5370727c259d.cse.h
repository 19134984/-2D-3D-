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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x90 = parameters.template get<descriptors::OMEGA>();
auto x91 = cell[10] + cell[14] + cell[16];
auto x92 = cell[1] + cell[6] + x91;
auto x93 = cell[13] + cell[5];
auto x94 = cell[11] + cell[18] + x93;
auto x95 = cell[2] + cell[8] + x94;
auto x96 = cell[15] + cell[7];
auto x97 = cell[17] + cell[9];
auto x98 = cell[12] + x96 + x97;
auto x99 = cell[3] + x98;
auto x0 = V{3};
auto x109 = ((x0)*(x0));
auto x1 = -V{1}*cell[4];
auto x2 = -V{1}*cell[5];
auto x101 = cell[13] + x2;
auto x3 = -V{1}*cell[6];
auto x4 = -V{1}*cell[7];
auto x102 = cell[15] + x4;
auto x5 = cell[0] + cell[4] + x92 + x95 + x99 + V{1};
auto x100 = V{1} / (x5);
auto x6 = x100*(-V{1}*cell[1] + x1 + x101 + x102 + x3 + x91);
auto x106 = ((x6)*(x6));
auto x7 = -V{1}*cell[8];
auto x8 = -V{1}*cell[9];
auto x103 = cell[17] + x8;
auto x9 = -V{1}*cell[14];
auto x10 = x100*(-V{1}*cell[2] + x1 + x103 + x7 + x9 + x94);
auto x104 = ((x10)*(x10));
auto x11 = -V{1}*cell[16];
auto x12 = -V{1}*cell[18];
auto x13 = x100*(-V{1}*cell[3] + x11 + x12 + x3 + x7 + x98);
auto x105 = ((x13)*(x13));
auto x14 = x104 + x105 + x106;
auto x15 = V{0.5}*x0*x14;
auto x116 = V{1} - x15;
auto x16 = V{0.333333333333333};
auto x17 = V{1}/x0;
auto x18 = -x17;
auto x19 = x5*x6;
auto x20 = x17*(x5 + V{-1});
auto x107 = -x20;
auto x108 = cell[4] + x107;
auto x21 = x108 - x19*x6 + x92 + x93 + x96;
auto x22 = x18*x21;
auto x23 = x10*x5;
auto x24 = cell[14] - x10*x23 + x108 + x95 + x97;
auto x25 = x18*x24;
auto x26 = x22 + x25;
auto x27 = cell[16] + cell[18] + cell[6] + cell[8] - x105*x5 + x107 + x99;
auto x28 = x18*x27;
auto x29 = V{1} - x90;
auto x30 = -V{1}*x6;
auto x31 = V{0.5}*x109;
auto x32 = V{0.0555555555555556};
auto x117 = -x32;
auto x33 = x32*x5;
auto x34 = V{0.5}*x109*x32;
auto x111 = x29*x34;
auto x35 = V{1} - x17;
auto x36 = x21*x35;
auto x110 = x28 + x36;
auto x37 = x111*(x110 + x25);
auto x118 = x117 + x37;
auto x38 = -V{1}*x10;
auto x39 = x24*x35;
auto x112 = x110 + x39;
auto x40 = x22 + x39;
auto x41 = x111*(x28 + x40);
auto x119 = x117 + x41;
auto x42 = -V{1}*x13;
auto x43 = x27*x35;
auto x114 = x25 + x36 + x43;
auto x115 = x40 + x43;
auto x44 = x111*(x26 + x43);
auto x120 = x117 + x44;
auto x45 = x30 + x38;
auto x46 = V{0.0277777777777778};
auto x121 = -x46;
auto x47 = x46*x5;
auto x48 = V{0.5}*x109*x46;
auto x113 = x29*x48;
auto x49 = V{2};
auto x50 = cell[4] - x10*x19 + x101 + x9;
auto x51 = x113*(x112 + x49*x50);
auto x122 = x121 + x51;
auto x52 = x10 + x30;
auto x53 = V{-2};
auto x54 = x113*(x112 + x50*x53);
auto x123 = x121 + x54;
auto x55 = x30 + x42;
auto x56 = cell[6] + x102 + x11 - x13*x19;
auto x57 = x113*(x114 + x49*x56);
auto x124 = x121 + x57;
auto x58 = x13 + x30;
auto x59 = x113*(x114 + x53*x56);
auto x125 = x121 + x59;
auto x60 = x38 + x42;
auto x61 = cell[8] + x103 + x12 - x13*x23;
auto x62 = x113*(x115 + x49*x61);
auto x126 = x121 + x62;
auto x63 = x13 + x38;
auto x64 = x113*(x115 + x53*x61);
auto x127 = x121 + x64;
auto x65 = x10 + x6;
auto x66 = x38 + x6;
auto x67 = x13 + x6;
auto x68 = x42 + x6;
auto x69 = x10 + x13;
auto x70 = x10 + x42;
cell[0] = x16*(V{0.5}*x109*x29*(x26 + x28) - x5*(x15 + V{-1}) + V{-1});
cell[10] = x118 + x33*(x0*x6 + x106*x31 + x116);
cell[11] = x119 + x33*(x0*x10 + x104*x31 + x116);
cell[12] = x120 + x33*(x0*x13 + x105*x31 + x116);
cell[13] = x122 + x47*(x0*x65 + x116 + x31*((x65)*(x65)));
cell[14] = x123 + x47*(x0*x66 + x116 + x31*((x66)*(x66)));
cell[15] = x124 + x47*(x0*x67 + x116 + x31*((x67)*(x67)));
cell[16] = x125 + x47*(x0*x68 + x116 + x31*((x68)*(x68)));
cell[17] = x126 + x47*(x0*x69 + x116 + x31*((x69)*(x69)));
cell[18] = x127 + x47*(x0*x70 + x116 + x31*((x70)*(x70)));
cell[1] = x118 + x33*(x0*x30 + x116 + x31*((x30)*(x30)));
cell[2] = x119 + x33*(x0*x38 + x116 + x31*((x38)*(x38)));
cell[3] = x120 + x33*(x0*x42 + x116 + x31*((x42)*(x42)));
cell[4] = x122 + x47*(x0*x45 + x116 + x31*((x45)*(x45)));
cell[5] = x123 + x47*(x0*x52 + x116 + x31*((x52)*(x52)));
cell[6] = x124 + x47*(x0*x55 + x116 + x31*((x55)*(x55)));
cell[7] = x125 + x47*(x0*x58 + x116 + x31*((x58)*(x58)));
cell[8] = x126 + x47*(x0*x60 + x116 + x31*((x60)*(x60)));
cell[9] = x127 + x47*(x0*x63 + x116 + x31*((x63)*(x63)));
return { x5, x14 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
