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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x19 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x21 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x20 = parameters.template get<descriptors::OMEGA>();
auto x22 = cell[13] + cell[15];
auto x23 = cell[10] + x22;
auto x24 = cell[11] + cell[5];
auto x25 = cell[12] + cell[7] + cell[9];
auto x26 = cell[2] + cell[8];
auto x27 = cell[18] + cell[3];
auto x28 = cell[16] + cell[6];
auto x29 = cell[0] + cell[14] + cell[17] + cell[1] + cell[4] + x23 + x24 + x25 + x26 + x27 + x28 + V{1};
auto x30 = V{1} / (x29);
auto x31 = cell[13] + cell[17];
auto x32 = -cell[2];
auto x33 = -cell[9];
auto x34 = -cell[8];
auto x35 = cell[18] + x34;
auto x36 = x32 + x33 + x35;
auto x37 = -cell[4];
auto x38 = -cell[14] + x24 + x37;
auto x39 = x31 + x36 + x38;
auto x40 = -cell[1];
auto x41 = -cell[7];
auto x42 = -cell[6];
auto x43 = cell[16] + x42;
auto x44 = x40 + x41 + x43;
auto x45 = cell[14] + x37;
auto x46 = -cell[5] + x45;
auto x47 = x23 + x44 + x46;
auto x48 = x30*x47;
auto x49 = ((x19)*(x19));
auto x50 = -cell[15];
auto x51 = cell[15] + cell[17];
auto x52 = -cell[18];
auto x53 = -cell[3];
auto x54 = x34 + x52 + x53;
auto x55 = -cell[16] + x42;
auto x56 = x25 + x51 + x54 + x55;
auto x57 = x30*x39;
auto x58 = -cell[17];
auto x59 = cell[9] + x58;
auto x60 = V{1}*x30;
auto x61 = ((x56)*(x56));
auto x62 = x49*x61;
auto x63 = V{0.333333333333333}*cell[0];
auto x64 = V{0.333333333333333}*cell[10] - V{0.666666666666667}*cell[17] - V{0.666666666666667}*cell[18] + V{0.333333333333333}*cell[1] - V{0.666666666666667}*cell[8] - V{0.666666666666667}*cell[9] + x63;
auto x65 = V{0.333333333333333}*cell[11] - V{0.666666666666667}*cell[15] - V{0.666666666666667}*cell[16] + V{0.333333333333333}*cell[2] - V{0.666666666666667}*cell[6] - V{0.666666666666667}*cell[7];
auto x66 = ((x39)*(x39));
auto x67 = x49*x66;
auto x68 = V{0.333333333333333}*cell[12] - V{0.666666666666667}*cell[13] - V{0.666666666666667}*cell[14] + V{0.333333333333333}*cell[3] - V{0.666666666666667}*cell[4] - V{0.666666666666667}*cell[5];
auto x69 = ((x47)*(x47));
auto x70 = x49*x69;
auto x71 = V{1} / (util::sqrt(V{25.455852}*x30*((x21)*(x21))*util::sqrt(V{2}*((x35 + x49*x56*x57 + x59)*(x35 + x49*x56*x57 + x59)) + V{2}*((-cell[13] + cell[5] + x39*x48*x49 + x45)*(-cell[13] + cell[5] + x39*x48*x49 + x45)) + V{2}*((cell[7] + x43 + x48*x49*x56 + x50)*(cell[7] + x43 + x48*x49*x56 + x50)) + ((-V{0.666666666666667}*cell[11] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] - V{0.666666666666667}*cell[2] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + x60*x67 + x64 + x68)*(-V{0.666666666666667}*cell[11] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] - V{0.666666666666667}*cell[2] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + x60*x67 + x64 + x68)) + ((-V{0.666666666666667}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] - V{0.666666666666667}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + x60*x62 + x64 + x65)*(-V{0.666666666666667}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] - V{0.666666666666667}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + x60*x62 + x64 + x65)) + ((-V{0.666666666666667}*cell[10] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] - V{0.666666666666667}*cell[1] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] + x60*x70 + x63 + x65 + x68)*(-V{0.666666666666667}*cell[10] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] - V{0.666666666666667}*cell[1] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] + x60*x70 + x63 + x65 + x68))) + V{1}/((x20)*(x20))) + V{1} / (x20));
auto x72 = V{1} - V{2}*x71;
auto x73 = V{1} / ((x29)*(x29));
auto x74 = V{1.5}*x73;
auto x75 = x70*x74;
auto x76 = x67*x74;
auto x77 = x62*x74;
auto x78 = x76 + x77 + V{-1};
auto x79 = x75 + x78;
auto x80 = V{0.111111111111111}*cell[0] + V{0.111111111111111}*cell[10] + V{0.111111111111111}*cell[11] + V{0.111111111111111}*cell[12] + V{0.111111111111111}*cell[13] + V{0.111111111111111}*cell[14] + V{0.111111111111111}*cell[15] + V{0.111111111111111}*cell[16] + V{0.111111111111111}*cell[17] + V{0.111111111111111}*cell[18] + V{0.111111111111111}*cell[1] + V{0.111111111111111}*cell[2] + V{0.111111111111111}*cell[3] + V{0.111111111111111}*cell[4] + V{0.111111111111111}*cell[5] + V{0.111111111111111}*cell[6] + V{0.111111111111111}*cell[7] + V{0.111111111111111}*cell[8] + V{0.111111111111111}*cell[9] + V{0.111111111111111};
auto x81 = V{3}*x19;
auto x82 = x48*x81;
auto x83 = V{3}*x73;
auto x84 = x70*x83;
auto x85 = x57*x81;
auto x86 = x67*x83;
auto x87 = x75 + V{-1};
auto x88 = x30*x56*x81;
auto x89 = x62*x83;
auto x90 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x91 = V{4.5}*x73;
auto x92 = cell[10] + x44;
auto x93 = x49*x91*((cell[11] + 2*cell[13] - 2*cell[4] + x36 + x51 + x92)*(cell[11] + 2*cell[13] - 2*cell[4] + x36 + x51 + x92));
auto x94 = x79 + x82;
auto x95 = -x85;
auto x96 = -cell[11] + V{2}*cell[14] + cell[15] - V{2}*cell[5] + x26 + x52 + x59 + x92;
auto x97 = cell[10] + x40 + x46;
auto x98 = x49*x91*((cell[12] + 2*cell[15] - 2*cell[6] + cell[9] + x31 + x54 + x97)*(cell[12] + 2*cell[15] - 2*cell[6] + cell[9] + x31 + x54 + x97));
auto x99 = -x88;
auto x100 = -cell[12] + cell[13];
auto x101 = V{2}*cell[16] - V{2}*cell[7] + cell[8] + x100 + x27 + x33 + x58 + x97;
auto x102 = x32 + x38;
auto x103 = x49*x91*((cell[12] + 2*cell[17] + cell[7] - 2*cell[8] + x102 + x22 + x53 + x55)*(cell[12] + 2*cell[17] + cell[7] - 2*cell[8] + x102 + x22 + x53 + x55));
auto x104 = x79 + x85;
auto x105 = V{2}*cell[18] + cell[3] - V{2}*cell[9] + x100 + x102 + x28 + x41 + x50;
auto x106 = -x76;
auto x107 = V{1} - x77;
auto x108 = x106 + x107;
auto x109 = x108 + x82;
auto x110 = -x75;
auto x111 = x110 + x85;
auto x112 = x110 + x88;
auto x113 = -x82;
auto x114 = x79 + x88;
auto x0 = cell[0]*x72 - x71*(x79*(V{0.666666666666667}*cell[0] + V{0.666666666666667}*cell[10] + V{0.666666666666667}*cell[11] + V{0.666666666666667}*cell[12] + V{0.666666666666667}*cell[13] + V{0.666666666666667}*cell[14] + V{0.666666666666667}*cell[15] + V{0.666666666666667}*cell[16] + V{0.666666666666667}*cell[17] + V{0.666666666666667}*cell[18] + V{0.666666666666667}*cell[1] + V{0.666666666666667}*cell[2] + V{0.666666666666667}*cell[3] + V{0.666666666666667}*cell[4] + V{0.666666666666667}*cell[5] + V{0.666666666666667}*cell[6] + V{0.666666666666667}*cell[7] + V{0.666666666666667}*cell[8] + V{0.666666666666667}*cell[9] + V{0.666666666666667}) + V{0.666666666666667});
auto x1 = cell[1]*x72 - x71*(x80*(x78 + x82 - x84) + V{0.111111111111111});
auto x2 = cell[2]*x72 - x71*(x80*(x77 + x85 - x86 + x87) + V{0.111111111111111});
auto x3 = cell[3]*x72 - x71*(x80*(x76 + x87 + x88 - x89) + V{0.111111111111111});
auto x4 = cell[4]*x72 - x71*(x90*(x85 - x93 + x94) + V{0.0555555555555556});
auto x5 = cell[5]*x72 - x71*(x90*(-x49*x91*((x96)*(x96)) + x94 + x95) + V{0.0555555555555556});
auto x6 = cell[6]*x72 - x71*(x90*(x88 + x94 - x98) + V{0.0555555555555556});
auto x7 = cell[7]*x72 - x71*(x90*(-x49*x91*((x101)*(x101)) + x94 + x99) + V{0.0555555555555556});
auto x8 = cell[8]*x72 - x71*(x90*(-x103 + x104 + x88) + V{0.0555555555555556});
auto x9 = cell[9]*x72 - x71*(x90*(x104 - x49*x91*((x105)*(x105)) + x99) + V{0.0555555555555556});
auto x10 = cell[10]*x72 + x71*(x80*(x109 + x84) + V{-0.111111111111111});
auto x11 = cell[11]*x72 + x71*(x80*(x107 + x111 + x86) + V{-0.111111111111111});
auto x12 = cell[12]*x72 + x71*(x80*(x106 + x112 + x89 + V{1}) + V{-0.111111111111111});
auto x13 = cell[13]*x72 + x71*(x90*(x109 + x111 + x93) + V{-0.0555555555555556});
auto x14 = cell[14]*x72 - x71*(x90*(x104 + x113 - x49*x91*((x96)*(x96))) + V{0.0555555555555556});
auto x15 = cell[15]*x72 + x71*(x90*(x109 + x112 + x98) + V{-0.0555555555555556});
auto x16 = cell[16]*x72 - x71*(x90*(x113 + x114 - x49*x91*((x101)*(x101))) + V{0.0555555555555556});
auto x17 = cell[17]*x72 + x71*(x90*(x103 + x108 + x111 + x88) + V{-0.0555555555555556});
auto x18 = cell[18]*x72 - x71*(x90*(x114 - x49*x91*((x105)*(x105)) + x95) + V{0.0555555555555556});
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
return { x29, V{1}*x49*x73*(x61 + x66 + x69) };
}
};

}

}

// Generation Info: commit=f977196d59888f9c6890546e1607f7031d0712fe
