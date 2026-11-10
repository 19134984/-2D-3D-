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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::FrozenOmegaFromCell<collision::BGK>, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = cell.template getFieldComponent<olb::descriptors::AVERAGE_TURBULENT_TAU>(0);
auto x21 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x20 = cell.template getFieldComponent<olb::descriptors::OMEGA>(0);
auto x22 = V{1} / (x19 + V{1}/x20);
auto x23 = cell[10] + cell[14];
auto x24 = cell[12] + cell[7];
auto x25 = cell[0] + cell[11] + cell[13] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[8] + cell[9] + x23 + x24 + V{1};
auto x26 = V{1} / ((x25)*(x25));
auto x27 = V{1.5}*x26;
auto x28 = ((x21)*(x21));
auto x29 = cell[13] - cell[4];
auto x30 = cell[15] - cell[6];
auto x31 = x29 + x30;
auto x32 = -cell[1];
auto x33 = cell[16] - cell[7];
auto x34 = x32 + x33;
auto x35 = -cell[5] + x23;
auto x36 = x31 + x34 + x35;
auto x37 = ((x36)*(x36));
auto x38 = x28*x37;
auto x39 = x27*x38;
auto x40 = cell[17] - cell[8];
auto x41 = x29 + x40;
auto x42 = cell[18] - cell[9];
auto x43 = -cell[2];
auto x44 = cell[11] - cell[14] + cell[5] + x43;
auto x45 = x41 + x42 + x44;
auto x46 = ((x45)*(x45));
auto x47 = x28*x46;
auto x48 = x27*x47;
auto x49 = x30 + x40;
auto x50 = -cell[3];
auto x51 = -cell[18] + cell[9];
auto x52 = x50 + x51;
auto x53 = -cell[16] + x24;
auto x54 = x49 + x52 + x53;
auto x55 = ((x54)*(x54));
auto x56 = x28*x55;
auto x57 = x27*x56;
auto x58 = x48 + x57 + V{-1};
auto x59 = x39 + x58;
auto x60 = V{1} - V{1}*x22;
auto x61 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x62 = V{3}*x21/x25;
auto x63 = x36*x62;
auto x64 = V{3}*x26;
auto x65 = x38*x64;
auto x66 = x45*x62;
auto x67 = x47*x64;
auto x68 = x39 + V{-1};
auto x69 = x54*x62;
auto x70 = x56*x64;
auto x71 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x72 = V{4.5}*x26;
auto x73 = cell[10] + x34;
auto x74 = x28*x72*((cell[11] + 2*cell[13] - 2*cell[4] + x42 + x43 + x49 + x73)*(cell[11] + 2*cell[13] - 2*cell[4] + x42 + x43 + x49 + x73));
auto x75 = x59 + x63;
auto x76 = -x66;
auto x77 = -cell[17] + cell[8];
auto x78 = -cell[11] + V{2}*cell[14] + cell[2] - V{2}*cell[5] + x30 + x51 + x73 + x77;
auto x79 = x32 + x35;
auto x80 = x28*x72*((cell[12] + 2*cell[15] - 2*cell[6] + x41 + x52 + x79)*(cell[12] + 2*cell[15] - 2*cell[6] + x41 + x52 + x79));
auto x81 = -x69;
auto x82 = -cell[12] + cell[3] + x29;
auto x83 = V{2}*cell[16] - V{2}*cell[7] + x42 + x77 + x79 + x82;
auto x84 = x28*x72*((2*cell[17] - 2*cell[8] + x31 + x44 + x50 + x53)*(2*cell[17] - 2*cell[8] + x31 + x44 + x50 + x53));
auto x85 = x59 + x66;
auto x86 = -cell[15] + V{2}*cell[18] + cell[6] - V{2}*cell[9] + x33 + x44 + x82;
auto x87 = -x48;
auto x88 = V{1} - x57;
auto x89 = x87 + x88;
auto x90 = x63 + x89;
auto x91 = -x39;
auto x92 = x66 + x91;
auto x93 = x69 + x91;
auto x94 = -x63;
auto x95 = x59 + x69;
auto x0 = cell[0]*x60 - x22*(x59*(V{0.333333333333333}*cell[0] + V{0.333333333333333}*cell[10] + V{0.333333333333333}*cell[11] + V{0.333333333333333}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] + V{0.333333333333333}*cell[1] + V{0.333333333333333}*cell[2] + V{0.333333333333333}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] + V{0.333333333333333}) + V{0.333333333333333});
auto x1 = cell[1]*x60 - x22*(x61*(x58 + x63 - x65) + V{0.0555555555555556});
auto x2 = cell[2]*x60 - x22*(x61*(x57 + x66 - x67 + x68) + V{0.0555555555555556});
auto x3 = cell[3]*x60 - x22*(x61*(x48 + x68 + x69 - x70) + V{0.0555555555555556});
auto x4 = cell[4]*x60 - x22*(x71*(x66 - x74 + x75) + V{0.0277777777777778});
auto x5 = cell[5]*x60 - x22*(x71*(-x28*x72*((x78)*(x78)) + x75 + x76) + V{0.0277777777777778});
auto x6 = cell[6]*x60 - x22*(x71*(x69 + x75 - x80) + V{0.0277777777777778});
auto x7 = cell[7]*x60 - x22*(x71*(-x28*x72*((x83)*(x83)) + x75 + x81) + V{0.0277777777777778});
auto x8 = cell[8]*x60 - x22*(x71*(x69 - x84 + x85) + V{0.0277777777777778});
auto x9 = cell[9]*x60 - x22*(x71*(-x28*x72*((x86)*(x86)) + x81 + x85) + V{0.0277777777777778});
auto x10 = cell[10]*x60 + x22*(x61*(x65 + x90) + V{-0.0555555555555556});
auto x11 = cell[11]*x60 + x22*(x61*(x67 + x88 + x92) + V{-0.0555555555555556});
auto x12 = cell[12]*x60 + x22*(x61*(x70 + x87 + x93 + V{1}) + V{-0.0555555555555556});
auto x13 = cell[13]*x60 + x22*(x71*(x74 + x90 + x92) + V{-0.0277777777777778});
auto x14 = cell[14]*x60 - x22*(x71*(-x28*x72*((x78)*(x78)) + x85 + x94) + V{0.0277777777777778});
auto x15 = cell[15]*x60 + x22*(x71*(x80 + x90 + x93) + V{-0.0277777777777778});
auto x16 = cell[16]*x60 - x22*(x71*(-x28*x72*((x83)*(x83)) + x94 + x95) + V{0.0277777777777778});
auto x17 = cell[17]*x60 + x22*(x71*(x69 + x84 + x89 + x92) + V{-0.0277777777777778});
auto x18 = cell[18]*x60 - x22*(x71*(-x28*x72*((x86)*(x86)) + x76 + x95) + V{0.0277777777777778});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
cell[9] = x9;
cell[10] = x10;
cell[11] = x11;
cell[12] = x12;
cell[13] = x13;
cell[14] = x14;
cell[15] = x15;
cell[16] = x16;
cell[17] = x17;
cell[18] = x18;
return { x25, V{1}*x26*x28*(x37 + x46 + x55) };
}
};

}

}

// Generation Info: commit=9eea65bdd916a6bed3e590d6f1e0706fbbe388fa
