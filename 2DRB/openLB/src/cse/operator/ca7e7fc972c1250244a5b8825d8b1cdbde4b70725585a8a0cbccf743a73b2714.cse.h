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

namespace operators {

template <typename T, typename...FIELDS>
struct CSE_O<operators::TrackAverageFunctorO<functors::ReverseDerivativeF<functors::DissipationF, descriptors::POPULATION, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, dynamics::DefaultCombination> >, opti::AVERAGE_DJDF>,descriptors::D3Q19<FIELDS...>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
void apply(CELL& cell, PARAMETERS& parameters) any_platform {
using V = typename CELL::value_t;
auto x19 = parameters.template get<parameters::OMEGA>();
auto x21 = parameters.template get<parameters::PHYS_CHAR_VISCOSITY>();
auto x20 = parameters.template get<parameters::DT>();
auto x22 = parameters.template get<parameters::FACTOR>();
auto x23 = parameters.template get<descriptors::LATTICE_TIME_T>();
auto x24 = V{1} / (x23);
auto x25 = x23 + V{-1};
auto x26 = x21*x22;
auto x27 = cell[12] + cell[15] + cell[17];
auto x28 = cell[11] + cell[13] + cell[18];
auto x29 = cell[10] + cell[14] + cell[16];
auto x30 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + x27 + x28 + x29 + V{1};
auto x31 = V{1} / (x30);
auto x32 = V{1} / (x20);
auto x33 = x19*x32;
auto x34 = x31*x33;
auto x35 = util::pow(x34, V{2});
auto x36 = V{3}*x35;
auto x37 = V{0.333333333333333}*cell[0];
auto x38 = V{0.333333333333333}*cell[17];
auto x39 = V{0.333333333333333}*cell[18];
auto x40 = V{0.333333333333333}*cell[8];
auto x41 = V{0.333333333333333}*cell[9];
auto x42 = V{0.666666666666667}*cell[10];
auto x43 = V{0.666666666666667}*cell[1];
auto x44 = V{1}*x31;
auto x45 = -cell[4];
auto x46 = -cell[6];
auto x47 = cell[13] + cell[15] - cell[1] - cell[5] - cell[7] + x29 + x45 + x46;
auto x48 = ((x47)*(x47));
auto x49 = x44*x48;
auto x50 = V{0.333333333333333}*cell[11];
auto x51 = V{0.333333333333333}*cell[2];
auto x52 = V{0.666666666666667}*cell[15];
auto x53 = V{0.666666666666667}*cell[16];
auto x54 = V{0.666666666666667}*cell[6];
auto x55 = V{0.666666666666667}*cell[7];
auto x56 = x50 + x51 - x52 - x53 - x54 - x55;
auto x57 = V{0.333333333333333}*cell[12];
auto x58 = V{0.333333333333333}*cell[3];
auto x59 = V{0.666666666666667}*cell[13];
auto x60 = V{0.666666666666667}*cell[14];
auto x61 = V{0.666666666666667}*cell[4];
auto x62 = V{0.666666666666667}*cell[5];
auto x63 = x57 + x58 - x59 - x60 - x61 - x62;
auto x64 = x37 + x38 + x39 + x40 + x41 - x42 - x43 + x49 + x56 + x63;
auto x65 = x36*x64;
auto x66 = -cell[8];
auto x67 = cell[5] + x45;
auto x68 = -cell[14] + cell[17] - cell[2] - cell[9] + x28 + x66 + x67;
auto x69 = ((x68)*(x68));
auto x70 = V{0.333333333333333}*cell[10];
auto x71 = V{0.333333333333333}*cell[1];
auto x72 = V{0.666666666666667}*cell[17];
auto x73 = V{0.666666666666667}*cell[18];
auto x74 = V{0.666666666666667}*cell[8];
auto x75 = V{0.666666666666667}*cell[9];
auto x76 = x37 + x70 + x71 - x72 - x73 - x74 - x75;
auto x77 = -V{0.666666666666667}*cell[11] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] - V{0.666666666666667}*cell[2] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + x44*x69 + x63 + x76;
auto x78 = x36*x77;
auto x79 = V{0.333333333333333}*cell[13];
auto x80 = V{0.333333333333333}*cell[14];
auto x81 = V{0.333333333333333}*cell[4];
auto x82 = V{0.333333333333333}*cell[5];
auto x83 = V{0.666666666666667}*cell[12];
auto x84 = V{0.666666666666667}*cell[3];
auto x85 = cell[9] + x66;
auto x86 = cell[7] + x46;
auto x87 = -cell[16] - cell[18] - cell[3] + x27 + x85 + x86;
auto x88 = ((x87)*(x87));
auto x89 = x44*x88;
auto x90 = x56 + x76 + x79 + x80 + x81 + x82 - x83 - x84 + x89;
auto x91 = x36*x90;
auto x92 = V{1} / ((x30)*(x30));
auto x93 = x35*x92;
auto x94 = V{9}*x93;
auto x95 = x48*x64*x94;
auto x96 = x69*x77*x94;
auto x97 = x88*x90*x94;
auto x98 = V{18}*x93;
auto x99 = x31*x68;
auto x100 = -cell[13] + cell[14] + x47*x99 + x67;
auto x101 = -x100;
auto x102 = x47*x87;
auto x103 = -cell[15] + cell[16] + x102*x31 + x86;
auto x104 = -x103;
auto x105 = -x37 - x50 - x51 + x52 + x53 + x54 + x55;
auto x106 = x101*x68 + x104*x87 + x47*(x105 - x38 - x39 - x40 - x41 + x42 + x43 - x49 - x57 - x58 + x59 + x60 + x61 + x62);
auto x107 = x106*x47*x98;
auto x108 = -cell[17] + cell[18] + x85 + x87*x99;
auto x109 = -x108;
auto x110 = x104*x47 + x109*x68 + x87*(x105 - x70 - x71 + x72 + x73 + x74 + x75 - x79 - x80 - x81 - x82 + x83 + x84 - x89);
auto x111 = x110*x87*x98;
auto x112 = x101*x47;
auto x113 = x109*x87;
auto x114 = -x112 - x113 + x68*x77;
auto x115 = -V{18}*x114*x35*x68*x92;
auto x116 = x102*x103*x98;
auto x117 = util::pow(x34, V{1});
auto x118 = V{2}*((x100)*(x100)) + V{2}*((x103)*(x103)) + V{2}*((x108)*(x108)) + ((x64)*(x64)) + ((x77)*(x77)) + ((x90)*(x90));
auto x119 = -V{9}*x117*x118*x19*x32*x92;
auto x120 = -V{18}*x101*x35*x47*x68*x92;
auto x121 = -V{18}*x109*x35*x68*x87*x92;
auto x122 = x107 + x111 + x115 + x116 + x119 + x120 + x121 + x91 + x95 + x96 + x97;
auto x123 = x122 + x78;
auto x124 = V{18}*x35;
auto x125 = x124*x31;
auto x126 = x106*x125;
auto x127 = -V{6}*x35*x64;
auto x128 = x123 + x127;
auto x129 = x114*x125;
auto x130 = -x129;
auto x131 = -V{6}*x35*x77;
auto x132 = x122 + x131;
auto x133 = x132 + x65;
auto x134 = x110*x125;
auto x135 = -V{6}*x35*x90;
auto x136 = x107 + x111 + x115 + x116 + x119 + x120 + x121 + x135 + x78 + x95 + x96 + x97;
auto x137 = x136 + x65;
auto x138 = -x126;
auto x139 = x100*x124;
auto x140 = x129 + x139;
auto x141 = V{6}*x35;
auto x142 = x141*x77;
auto x143 = x68*x98;
auto x144 = x114*x143;
auto x145 = -x95;
auto x146 = -x96;
auto x147 = -x97;
auto x148 = -x107;
auto x149 = -x111;
auto x150 = V{9}*x117*x118*x33*x92;
auto x151 = x112*x143;
auto x152 = x113*x143;
auto x153 = -x116;
auto x154 = x141*x64 + x144 + x145 + x146 + x147 + x148 + x149 + x150 + x151 + x152 + x153;
auto x155 = x142 + x154 - x91;
auto x156 = x126 + x127;
auto x157 = x103*x124;
auto x158 = x136 + x156;
auto x159 = -x134;
auto x160 = x108*x124;
auto x161 = x129 + x160;
auto x162 = x141*x90;
auto x163 = x142 + x144 + x145 + x146 + x147 + x148 + x149 + x150 + x151 + x152 + x153 + x162 - x65;
auto x164 = x134 + x163;
auto x165 = x126 + x155;
auto x166 = x134 + x157;
auto x0 = cell[0];
auto x1 = cell[1];
auto x2 = cell[2];
auto x3 = cell[3];
auto x4 = cell[4];
auto x5 = cell[5];
auto x6 = cell[6];
auto x7 = cell[7];
auto x8 = cell[8];
auto x9 = cell[9];
auto x10 = cell[10];
auto x11 = cell[11];
auto x12 = cell[12];
auto x13 = cell[13];
auto x14 = cell[14];
auto x15 = cell[15];
auto x16 = cell[16];
auto x17 = cell[17];
auto x18 = cell[18];
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[0] = x24*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(0)*x25 - x26*(-x123 - x65));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[1] = x24*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(1)*x25 - x26*(-x126 - x128));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[2] = x24*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(2)*x25 - x26*(-x130 - x133));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[3] = x24*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(3)*x25 - x26*(-x134 - x137));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[4] = x24*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(4)*x25 - x26*(x138 + x140 + x155));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[5] = x24*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(5)*x25 - x26*(-x132 - x140 - x156));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[6] = x24*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(6)*x25 - x26*(-x134 + x157 - x158));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[7] = x24*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(7)*x25 - x26*(-x157 - x158 - x159));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[8] = x24*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(8)*x25 - x26*(x159 + x161 + x163));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[9] = x24*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(9)*x25 - x26*(x129 - x160 + x164));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[10] = x24*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(10)*x25 - x26*(-x128 - x138));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[11] = x24*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(11)*x25 - x26*(-x129 - x133));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[12] = x24*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(12)*x25 - x26*(-x137 - x159));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[13] = x24*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(13)*x25 - x26*(x130 + x139 + x165));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[14] = x24*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(14)*x25 - x26*(x129 - x139 + x165));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[15] = x24*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(15)*x25 - x26*(x126 + x154 + x162 + x166 - x78));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[16] = x24*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(16)*x25 - x26*(-x127 - x136 - x138 - x166));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[17] = x24*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(17)*x25 - x26*(x130 + x160 + x164));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[18] = x24*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(18)*x25 - x26*(-x107 - x111 - x115 - x116 - x119 - x120 - x121 - x131 - x134 - x135 - x161 - x65 - x95 - x96 - x97));
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
}
};

}

}

// Generation Info: commit=65ad251e99d893670b9a5c7db92e58559362372c
