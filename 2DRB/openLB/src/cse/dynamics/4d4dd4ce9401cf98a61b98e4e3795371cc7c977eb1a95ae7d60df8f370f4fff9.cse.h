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
struct CSE<SpongeLayerDynamics<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x88 = cell.template getFieldComponent<olb::descriptors::DAMPING>(0);
auto x89 = cell.template getFieldComponent<olb::descriptors::DENSITY>(0);
auto x90 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x91 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x92 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x93 = parameters.template get<descriptors::OMEGA>();
auto x94 = ((x90)*(x90));
auto x95 = ((x91)*(x91));
auto x96 = ((x92)*(x92));
auto x97 = cell[10] + cell[14] + cell[16];
auto x98 = cell[11] + cell[18] + cell[5];
auto x99 = cell[12] + cell[7] + cell[9];
auto x0 = V{3};
auto x1 = V{0.5}*x0;
auto x2 = x1*(x94 + x95 + x96);
auto x110 = V{1} - x2;
auto x3 = V{0.333333333333333};
auto x4 = -V{1}*cell[4];
auto x101 = cell[13] + x4;
auto x5 = -V{1}*cell[6];
auto x102 = cell[15] + x5;
auto x6 = cell[0] + cell[13] + cell[15] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + x97 + x98 + x99 + V{1};
auto x100 = V{1} / (x6);
auto x7 = x100*(-V{1}*cell[1] - V{1}*cell[5] - V{1}*cell[7] + x101 + x102 + x97);
auto x105 = ((x7)*(x7));
auto x8 = -V{1}*cell[8];
auto x103 = cell[17] + x8;
auto x9 = x100*(-V{1}*cell[14] - V{1}*cell[2] - V{1}*cell[9] + x101 + x103 + x98);
auto x106 = ((x9)*(x9));
auto x10 = x100*(-V{1}*cell[16] - V{1}*cell[18] - V{1}*cell[3] + x102 + x103 + x99);
auto x104 = ((x10)*(x10));
auto x11 = x104 + x105 + x106;
auto x12 = x1*x11;
auto x108 = V{1} - x12;
auto x13 = -x3*(x6*(x12 + V{-1}) + V{1});
auto x14 = -V{1}*x90;
auto x15 = V{0.5}*((x0)*(x0));
auto x16 = V{0.0555555555555556};
auto x107 = -x16;
auto x17 = x16*x89;
auto x18 = -V{1}*x7;
auto x19 = x16*x6;
auto x20 = x107 + x19*(x0*x18 + x108 + x15*((x18)*(x18)));
auto x21 = -V{1}*x91;
auto x22 = -V{1}*x9;
auto x23 = x107 + x19*(x0*x22 + x108 + x15*((x22)*(x22)));
auto x24 = -V{1}*x92;
auto x25 = -V{1}*x10;
auto x26 = x107 + x19*(x0*x25 + x108 + x15*((x25)*(x25)));
auto x27 = x14 + x21;
auto x28 = V{0.0277777777777778};
auto x109 = -x28;
auto x29 = x28*x89;
auto x30 = x18 + x22;
auto x31 = x28*x6;
auto x32 = x109 + x31*(x0*x30 + x108 + x15*((x30)*(x30)));
auto x33 = x14 + x91;
auto x34 = x18 + x9;
auto x35 = x109 + x31*(x0*x34 + x108 + x15*((x34)*(x34)));
auto x36 = x14 + x24;
auto x37 = x18 + x25;
auto x38 = x109 + x31*(x0*x37 + x108 + x15*((x37)*(x37)));
auto x39 = x14 + x92;
auto x40 = x10 + x18;
auto x41 = x109 + x31*(x0*x40 + x108 + x15*((x40)*(x40)));
auto x42 = x21 + x24;
auto x43 = x22 + x25;
auto x44 = x109 + x31*(x0*x43 + x108 + x15*((x43)*(x43)));
auto x45 = x21 + x92;
auto x46 = x10 + x22;
auto x47 = x109 + x31*(x0*x46 + x108 + x15*((x46)*(x46)));
auto x48 = x107 + x19*(x0*x7 + x105*x15 + x108);
auto x49 = x107 + x19*(x0*x9 + x106*x15 + x108);
auto x50 = x107 + x19*(x0*x10 + x104*x15 + x108);
auto x51 = x90 + x91;
auto x52 = x7 + x9;
auto x53 = x109 + x31*(x0*x52 + x108 + x15*((x52)*(x52)));
auto x54 = x21 + x90;
auto x55 = x22 + x7;
auto x56 = x109 + x31*(x0*x55 + x108 + x15*((x55)*(x55)));
auto x57 = x90 + x92;
auto x58 = x10 + x7;
auto x59 = x109 + x31*(x0*x58 + x108 + x15*((x58)*(x58)));
auto x60 = x24 + x90;
auto x61 = x25 + x7;
auto x62 = x109 + x31*(x0*x61 + x108 + x15*((x61)*(x61)));
auto x63 = x91 + x92;
auto x64 = x10 + x9;
auto x65 = x109 + x31*(x0*x64 + x108 + x15*((x64)*(x64)));
auto x66 = x24 + x91;
auto x67 = x25 + x9;
auto x68 = x109 + x31*(x0*x67 + x108 + x15*((x67)*(x67)));
auto x69 = cell[0] - x88*(x13 - x3*x89*(V{1} - x2) + x3) - x93*(cell[0] - x13);
auto x70 = cell[10] - x88*(x16 - x17*(x0*x90 + x110 + x15*x94) + x48) - x93*(cell[10] - x48);
auto x71 = cell[11] - x88*(x16 - x17*(x0*x91 + x110 + x15*x95) + x49) - x93*(cell[11] - x49);
auto x72 = cell[12] - x88*(x16 - x17*(x0*x92 + x110 + x15*x96) + x50) - x93*(cell[12] - x50);
auto x73 = -(-cell[13] + x88*(x28 - x29*(x0*x51 + x110 + x15*((x51)*(x51))) + x53) + x93*(cell[13] - x53));
auto x74 = -(-cell[14] + x88*(x28 - x29*(x0*x54 + x110 + x15*((x54)*(x54))) + x56) + x93*(cell[14] - x56));
auto x75 = -(-cell[15] + x88*(x28 - x29*(x0*x57 + x110 + x15*((x57)*(x57))) + x59) + x93*(cell[15] - x59));
auto x76 = -(-cell[16] + x88*(x28 - x29*(x0*x60 + x110 + x15*((x60)*(x60))) + x62) + x93*(cell[16] - x62));
auto x77 = -(-cell[17] + x88*(x28 - x29*(x0*x63 + x110 + x15*((x63)*(x63))) + x65) + x93*(cell[17] - x65));
auto x78 = -(-cell[18] + x88*(x28 - x29*(x0*x66 + x110 + x15*((x66)*(x66))) + x68) + x93*(cell[18] - x68));
auto x79 = -(-cell[1] + x88*(x16 - x17*(x0*x14 + x110 + x15*((x14)*(x14))) + x20) + x93*(cell[1] - x20));
auto x80 = -(-cell[2] + x88*(x16 - x17*(x0*x21 + x110 + x15*((x21)*(x21))) + x23) + x93*(cell[2] - x23));
auto x81 = -(-cell[3] + x88*(x16 - x17*(x0*x24 + x110 + x15*((x24)*(x24))) + x26) + x93*(cell[3] - x26));
auto x82 = -(-cell[4] + x88*(x28 - x29*(x0*x27 + x110 + x15*((x27)*(x27))) + x32) + x93*(cell[4] - x32));
auto x83 = -(-cell[5] + x88*(x28 - x29*(x0*x33 + x110 + x15*((x33)*(x33))) + x35) + x93*(cell[5] - x35));
auto x84 = -(-cell[6] + x88*(x28 - x29*(x0*x36 + x110 + x15*((x36)*(x36))) + x38) + x93*(cell[6] - x38));
auto x85 = -(-cell[7] + x88*(x28 - x29*(x0*x39 + x110 + x15*((x39)*(x39))) + x41) + x93*(cell[7] - x41));
auto x86 = -(-cell[8] + x88*(x28 - x29*(x0*x42 + x110 + x15*((x42)*(x42))) + x44) + x93*(cell[8] - x44));
auto x87 = -(-cell[9] + x88*(x28 - x29*(x0*x45 + x110 + x15*((x45)*(x45))) + x47) + x93*(cell[9] - x47));
cell[0] = x69;
cell[10] = x70;
cell[11] = x71;
cell[12] = x72;
cell[13] = x73;
cell[14] = x74;
cell[15] = x75;
cell[16] = x76;
cell[17] = x77;
cell[18] = x78;
cell[1] = x79;
cell[2] = x80;
cell[3] = x81;
cell[4] = x82;
cell[5] = x83;
cell[6] = x84;
cell[7] = x85;
cell[8] = x86;
cell[9] = x87;
return { x6, x11 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
