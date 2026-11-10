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
struct CSE_O<operators::TrackAverageFunctorO<functors::ReverseDerivativeF<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, descriptors::POPULATION, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination> >, opti::AVERAGE_DJDF>,descriptors::D3Q19<FIELDS...>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
void apply(CELL& cell, PARAMETERS& parameters) any_platform {
using V = typename CELL::value_t;
auto x110 = parameters.template get<parameters::DX>();
auto x114 = parameters.template get<parameters::DT>();
auto x115 = parameters.template get<parameters::FACTOR>();
auto x112 = parameters.template get<parameters::CONVERSION_VELOCITY>();
auto x116 = parameters.template get<descriptors::LATTICE_TIME_T>();
auto x109 = parameters.template get<parameters::OMEGA>();
auto x113 = parameters.template get<parameters::PHYS_CHAR_VISCOSITY>();
auto x111 = parameters.template get<parameters::LATTICE_VISCOSITY>();
auto x117 = cell[15] + cell[16];
auto x118 = cell[10] + cell[14] + x117;
auto x119 = cell[1] + cell[4] + x118;
auto x120 = cell[12] + cell[7] + cell[9];
auto x121 = cell[3] + cell[8] + x120;
auto x122 = cell[13] + cell[5];
auto x123 = cell[17] + cell[18];
auto x124 = cell[11] + x122 + x123;
auto x125 = cell[2] + x124;
auto x138 = V{2}*cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x139 = V{1} / (x116);
auto x0 = V{3};
auto x1 = x0*x109;
auto x127 = x1/x114;
auto x2 = cell[0] + cell[6] + x119 + x121 + x125 + V{1};
auto x126 = V{1} / (x2);
auto x130 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*x126;
auto x136 = V{1} / ((x2)*(x2));
auto x137 = V{1}*x136;
auto x3 = x126*x127;
auto x4 = -V{1}*cell[16];
auto x128 = cell[15] + x4;
auto x5 = -V{1}*cell[7];
auto x6 = -V{1}*cell[6];
auto x7 = -V{1}*cell[8];
auto x8 = -V{1}*cell[18];
auto x129 = cell[17] + x8;
auto x9 = -V{1}*cell[3] + x120 + x128 + x129 + x6 + x7;
auto x10 = x130*x9;
auto x11 = -V{1}*cell[4];
auto x12 = -V{1}*cell[5];
auto x131 = cell[13] + x12;
auto x13 = -V{1}*cell[1] + x11 + x118 + x131 + x5 + x6;
auto x14 = x13*x130;
auto x15 = x14*x2;
auto x16 = cell[6] - x10*x15 + x128 + x5;
auto x17 = ((x16)*(x16));
auto x18 = -V{1}*cell[14];
auto x19 = -V{1}*cell[9];
auto x20 = -V{1}*cell[2] + x11 + x124 + x18 + x19 + x7;
auto x21 = x130*x20;
auto x22 = cell[4] + x131 - x15*x21 + x18;
auto x132 = V{2}*x22;
auto x23 = x132;
auto x24 = V{1}/x0;
auto x25 = x24*(x2 + V{-1});
auto x133 = -x25;
auto x134 = cell[6] + x133;
auto x26 = cell[7] + x119 + x122 + x134 - x14*x15;
auto x27 = x10*x2;
auto x28 = -x10*x27 + x117 + x121 + x123 + x134;
auto x29 = x2*x21;
auto x30 = cell[8] - x10*x29 + x129 + x19;
auto x135 = V{2}*x30;
auto x31 = x135;
auto x32 = cell[14] + cell[4] + cell[8] + cell[9] + x125 + x133 - x21*x29;
auto x33 = V{0.5}*x113;
auto x34 = util::pow(x3, V{2})*x33;
auto x35 = x28*x34;
auto x36 = V{2}*x35;
auto x37 = -V{1}*x36;
auto x38 = x10*x37;
auto x39 = x34*(x135 + x31);
auto x155 = x36 + x39;
auto x40 = -V{1}*x39;
auto x157 = x36 + x40;
auto x41 = x10*x40;
auto x42 = x32*x34;
auto x43 = V{2}*x42;
auto x44 = -V{1}*x43;
auto x45 = x21*x44;
auto x46 = x16*x34;
auto x47 = V{4}*x46;
auto x48 = -V{1}*x47;
auto x49 = x10*x48;
auto x50 = x34*(x132 + x23);
auto x51 = -V{1}*x50;
auto x52 = x21*x51;
auto x53 = x26*x34;
auto x54 = V{2}*x53;
auto x55 = -V{1}*x54;
auto x56 = x14*x55;
auto x57 = x10*x38 - V{2}*x127*x136*util::pow(x3, V{1})*x33*(V{2}*x17 + x22*x23 + x30*x31 + ((x26)*(x26)) + ((x28)*(x28)) + ((x32)*(x32))) + x14*x49 + x14*x52 + x14*x56 + x21*x41 + x21*x45 + x24*x37 + x24*x44 + x24*x55;
auto x58 = -x137*x9;
auto x59 = -x109*x113*(cell.template getFieldComponent<olb::descriptors::POROSITY>(0) + V{-1})*((x112)*(x112))/(((x110)*(x110))*x111);
auto x60 = x10*x59;
auto x61 = x138*x60;
auto x62 = -x137*x20;
auto x63 = x21*x59;
auto x64 = x138*x63;
auto x65 = -x13*x137;
auto x66 = x14*x59;
auto x67 = x138*x66;
auto x68 = x58*x61 + x62*x64 + x65*x67;
auto x69 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*(x15*x48 + x2*x38 + x27*x37 + x29*x40);
auto x70 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*(x15*x51 + x2*x41 + x2*x45 + x29*x44);
auto x71 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*(x15*x55 + x2*x49 + x2*x52 + x2*x56);
auto x72 = x58*x69 + x62*x70 + x65*x71;
auto x140 = x57 + x72;
auto x141 = x140 + x68;
auto x142 = x141 + x54;
auto x144 = x140 + x43;
auto x147 = x140 + x36;
auto x73 = x116 + V{-1};
auto x74 = V{1}*x126;
auto x75 = x67*x74;
auto x76 = -V{1}*x75;
auto x77 = x71*x74;
auto x158 = x75 + x77;
auto x159 = x158 + x54;
auto x78 = -V{1}*x77;
auto x143 = x76 + x78;
auto x150 = x143 + x54;
auto x79 = x64*x74;
auto x80 = -V{1}*x79;
auto x81 = x68 + x80;
auto x82 = x70*x74;
auto x151 = x144 + x82;
auto x160 = x151 + x79;
auto x83 = -V{1}*x82;
auto x145 = x144 + x83;
auto x146 = x145 + x81;
auto x156 = x145 + x80;
auto x84 = x61*x74;
auto x85 = x68 - V{1}*x84;
auto x86 = x69*x74;
auto x87 = -V{1}*x86;
auto x148 = x85 + x87;
auto x149 = x147 + x148;
auto x88 = x68 + x79;
auto x152 = x151 + x88;
auto x89 = x68 + x84;
auto x153 = x86 + x89;
auto x154 = x147 + x153;
auto x90 = cell[0];
auto x91 = cell[1];
auto x92 = cell[2];
auto x93 = cell[3];
auto x94 = cell[4];
auto x95 = cell[5];
auto x96 = cell[6];
auto x97 = cell[7];
auto x98 = cell[8];
auto x99 = cell[9];
auto x100 = cell[10];
auto x101 = cell[11];
auto x102 = cell[12];
auto x103 = cell[13];
auto x104 = cell[14];
auto x105 = cell[15];
auto x106 = cell[16];
auto x107 = cell[17];
auto x108 = cell[18];
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[0] = x139*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(0)*x73 + x115*x141);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[1] = x139*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(1)*x73 + x115*(x142 + x143));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[2] = x139*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(2)*x73 + x115*x146);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[3] = x139*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(3)*x73 + x115*x149);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[4] = x139*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(4)*x73 + x115*(x146 + x150 + x50));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[5] = x139*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(5)*x73 + x115*(x150 + x152 + x51));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[6] = x139*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(6)*x73 + x115*(x149 + x150 + x47));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[7] = x139*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(7)*x73 + x115*(x150 + x154 + x48));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[8] = x139*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(8)*x73 + x115*(x148 + x155 + x156));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[9] = x139*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(9)*x73 + x115*(x153 + x156 + x157));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[10] = x139*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(10)*x73 + x115*(x142 + x158));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[11] = x139*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(11)*x73 + x115*x152);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[12] = x139*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(12)*x73 + x115*x154);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[13] = x139*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(13)*x73 + x115*(x152 + x159 + x50));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[14] = x139*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(14)*x73 + x115*(x146 + x159 + x51));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[15] = x139*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(15)*x73 + x115*(x154 + x159 + x47));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[16] = x139*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(16)*x73 + x115*(x149 + x159 + x48));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[17] = x139*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(17)*x73 + x115*(x153 + x155 + x160));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[18] = x139*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(18)*x73 + x115*(x148 + x157 + x160));
cell[0] = x90;
cell[1] = x91;
cell[2] = x92;
cell[3] = x93;
cell[4] = x94;
cell[5] = x95;
cell[6] = x96;
cell[7] = x97;
cell[8] = x98;
cell[9] = x99;
cell[10] = x100;
cell[11] = x101;
cell[12] = x102;
cell[13] = x103;
cell[14] = x104;
cell[15] = x105;
cell[16] = x106;
cell[17] = x107;
cell[18] = x108;
}
};

}

}

// Generation Info: commit=21521319fcfc23b85b6fd87340ebc70ad50346d4
