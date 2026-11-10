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
struct CSE_O<operators::TrackAverageFunctorO<functors::ReverseDerivativeF<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, descriptors::POPULATION, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, dynamics::DefaultCombination> >, opti::AVERAGE_DJDF>,descriptors::D3Q19<FIELDS...>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
void apply(CELL& cell, PARAMETERS& parameters) any_platform {
using V = typename CELL::value_t;
auto x19 = parameters.template get<parameters::OMEGA>();
auto x21 = parameters.template get<parameters::LATTICE_VISCOSITY>();
auto x20 = parameters.template get<parameters::DX>();
auto x22 = parameters.template get<parameters::CONVERSION_VELOCITY>();
auto x23 = parameters.template get<parameters::PHYS_CHAR_VISCOSITY>();
auto x24 = parameters.template get<parameters::DT>();
auto x26 = parameters.template get<descriptors::LATTICE_TIME_T>();
auto x25 = parameters.template get<parameters::FACTOR>();
auto x27 = V{1} / (x26);
auto x28 = x26 + V{-1};
auto x29 = x23*x25;
auto x30 = cell[12] + cell[15] + cell[17];
auto x31 = cell[11] + cell[13] + cell[18];
auto x32 = cell[10] + cell[14] + cell[16];
auto x33 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + x30 + x31 + x32 + V{1};
auto x34 = V{1} / (x33);
auto x35 = x19/x24;
auto x36 = x34*x35;
auto x37 = util::pow(x36, V{2});
auto x38 = V{3}*x37;
auto x39 = V{0.333333333333333}*cell[0];
auto x40 = V{0.333333333333333}*cell[17];
auto x41 = V{0.333333333333333}*cell[18];
auto x42 = V{0.333333333333333}*cell[8];
auto x43 = V{0.333333333333333}*cell[9];
auto x44 = V{0.666666666666667}*cell[10];
auto x45 = V{0.666666666666667}*cell[1];
auto x46 = V{1}*x34;
auto x47 = ((cell.template getFieldComponent<olb::descriptors::POROSITY>(0))*(cell.template getFieldComponent<olb::descriptors::POROSITY>(0)));
auto x48 = -cell[4];
auto x49 = -cell[6];
auto x50 = cell[13] + cell[15] - cell[1] - cell[5] - cell[7] + x32 + x48 + x49;
auto x51 = ((x50)*(x50));
auto x52 = x46*x47*x51;
auto x53 = V{0.333333333333333}*cell[11];
auto x54 = V{0.333333333333333}*cell[2];
auto x55 = V{0.666666666666667}*cell[15];
auto x56 = V{0.666666666666667}*cell[16];
auto x57 = V{0.666666666666667}*cell[6];
auto x58 = V{0.666666666666667}*cell[7];
auto x59 = x53 + x54 - x55 - x56 - x57 - x58;
auto x60 = V{0.333333333333333}*cell[12];
auto x61 = V{0.333333333333333}*cell[3];
auto x62 = V{0.666666666666667}*cell[13];
auto x63 = V{0.666666666666667}*cell[14];
auto x64 = V{0.666666666666667}*cell[4];
auto x65 = V{0.666666666666667}*cell[5];
auto x66 = x60 + x61 - x62 - x63 - x64 - x65;
auto x67 = x39 + x40 + x41 + x42 + x43 - x44 - x45 + x52 + x59 + x66;
auto x68 = x38*x67;
auto x69 = -cell[8];
auto x70 = cell[5] + x48;
auto x71 = -cell[14] + cell[17] - cell[2] - cell[9] + x31 + x69 + x70;
auto x72 = ((x71)*(x71));
auto x73 = V{0.333333333333333}*cell[10];
auto x74 = V{0.333333333333333}*cell[1];
auto x75 = V{0.666666666666667}*cell[17];
auto x76 = V{0.666666666666667}*cell[18];
auto x77 = V{0.666666666666667}*cell[8];
auto x78 = V{0.666666666666667}*cell[9];
auto x79 = x39 + x73 + x74 - x75 - x76 - x77 - x78;
auto x80 = -V{0.666666666666667}*cell[11] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] - V{0.666666666666667}*cell[2] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + x46*x47*x72 + x66 + x79;
auto x81 = x38*x80;
auto x82 = V{0.333333333333333}*cell[13];
auto x83 = V{0.333333333333333}*cell[14];
auto x84 = V{0.333333333333333}*cell[4];
auto x85 = V{0.333333333333333}*cell[5];
auto x86 = V{0.666666666666667}*cell[12];
auto x87 = V{0.666666666666667}*cell[3];
auto x88 = cell[7] + x49;
auto x89 = cell[9] + x69;
auto x90 = -cell[16] - cell[18] - cell[3] + x30 + x88 + x89;
auto x91 = ((x90)*(x90));
auto x92 = x46*x47*x91;
auto x93 = x59 + x79 + x82 + x83 + x84 + x85 - x86 - x87 + x92;
auto x94 = V{1} / ((x33)*(x33));
auto x95 = V{9}*x94;
auto x96 = x37*x95;
auto x97 = x47*x51*x67*x96;
auto x98 = x47*x72*x80*x96;
auto x99 = x47*x91*x93*x96;
auto x100 = V{18}*x37;
auto x101 = x100*x94;
auto x102 = x101*x50;
auto x103 = x34*x50;
auto x104 = -cell[15] + cell[16] + x103*x47*x90 + x88;
auto x105 = -cell[13] + cell[14] + x103*x47*x71 + x70;
auto x106 = -x105;
auto x107 = x47*(-x104*x90 + x106*x71 - x50*x67);
auto x108 = x102*x107;
auto x109 = x101*x71;
auto x110 = x71*x90;
auto x111 = -cell[17] + cell[18] + x110*x34*x47 + x89;
auto x112 = -x111;
auto x113 = x106*x50 + x112*x90 - x71*x80;
auto x114 = x113*x47;
auto x115 = x109*x114;
auto x116 = x101*x90;
auto x117 = x47*(-x104*x50 + x112*x71 - x90*x93);
auto x118 = x116*x117;
auto x119 = x35*util::pow(x36, V{1})*x95*(V{2}*((x104)*(x104)) + V{2}*((x105)*(x105)) + V{2}*((x111)*(x111)) + ((x67)*(x67)) + ((x80)*(x80)) + ((x93)*(x93)));
auto x120 = -x119;
auto x121 = x102*x71;
auto x122 = x105*x121*x47;
auto x123 = x102*x104*x47*x90;
auto x124 = x101*x110;
auto x125 = x111*x124*x47;
auto x126 = V{2}*x19*(cell.template getFieldComponent<olb::descriptors::POROSITY>(0) + V{-1})/(((x20)*(x20))*x21);
auto x127 = x126/util::pow(x33, 3);
auto x128 = ((x22)*(x22));
auto x129 = x127*x128*x47*x51;
auto x130 = x127*x128*x47*x72;
auto x131 = x127*x128*x47*x91;
auto x132 = x108 + x115 + x118 + x120 + x122 + x123 + x125 + x129 + x130 + x131 + x38*x93 + x97 + x98 + x99;
auto x133 = x132 + x81;
auto x134 = V{6}*x37;
auto x135 = x134*x67;
auto x136 = -x135;
auto x137 = x133 + x136;
auto x138 = x100*x34;
auto x139 = x107*x138;
auto x140 = x126*x94;
auto x141 = x128*x47;
auto x142 = x140*x141*x50;
auto x143 = x139 + x142;
auto x144 = -x134*x80;
auto x145 = x132 + x144;
auto x146 = x145 + x68;
auto x147 = x114*x138;
auto x148 = x140*x141*x71;
auto x149 = x147 + x148;
auto x150 = x134*x93;
auto x151 = -x150;
auto x152 = x108 + x115 + x118 + x120 + x122 + x123 + x125 + x129 + x130 + x131 + x151 + x81 + x97 + x98 + x99;
auto x153 = x152 + x68;
auto x154 = x117*x138;
auto x155 = x140*x141*x90;
auto x156 = x154 + x155;
auto x157 = x100*x105;
auto x158 = -x157;
auto x159 = x136 + x143;
auto x160 = x145 + x159;
auto x161 = -x147 - x148;
auto x162 = x100*x104;
auto x163 = x152 + x159;
auto x164 = -x154 - x155;
auto x165 = x100*x111;
auto x166 = -x165;
auto x167 = x108 + x115 + x118 + x120 + x122 + x123 + x125 + x129 + x130 + x131 + x144 + x149 + x151 + x68 + x97 + x98 + x99;
auto x168 = -x139 - x142;
auto x169 = x136 + x168;
auto x170 = x145 + x169;
auto x171 = -x104;
auto x172 = -x39 - x53 - x54 + x55 + x56 + x57 + x58;
auto x173 = x47*(x112*x71 + x171*x50 + x90*(x172 - x73 - x74 + x75 + x76 + x77 + x78 - x82 - x83 - x84 - x85 + x86 + x87 - x92));
auto x174 = x47*(x106*x71 + x171*x90 + x50*(x172 - x40 - x41 - x42 - x43 + x44 + x45 - x52 - x60 - x61 + x62 + x63 + x64 + x65));
auto x175 = x108 + x115 + x118 + x120 + x122 + x123 + x125 + x129 + x130 + x131 + x144 + x151 + x161 + x68 + x97 + x98 + x99;
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
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[0] = x27*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(0)*x28 + x29*(x133 + x68));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[1] = x27*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(1)*x28 + x29*(x137 + x143));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[2] = x27*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(2)*x28 + x29*(x146 + x149));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[3] = x27*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(3)*x28 + x29*(x153 + x156));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[4] = x27*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(4)*x28 + x29*(x149 + x158 + x160));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[5] = x27*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(5)*x28 + x29*(x157 + x160 + x161));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[6] = x27*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(6)*x28 + x29*(x156 - x162 + x163));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[7] = x27*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(7)*x28 + x29*(x162 + x163 + x164));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[8] = x27*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(8)*x28 + x29*(x156 + x166 + x167));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[9] = x27*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(9)*x28 + x29*(x164 + x165 + x167));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[10] = x27*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(10)*x28 + x29*(x137 + x168));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[11] = x27*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(11)*x28 + x29*(x146 + x161));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[12] = x27*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(12)*x28 + x29*(x153 + x164));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[13] = x27*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(13)*x28 + x29*(x158 + x161 + x170));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[14] = x27*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(14)*x28 + x29*(x149 + x157 + x170));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[15] = x27*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(15)*x28 - x29*(-x102*x174 + x106*x121*x47 - x109*x113*x47 + x112*x124*x47 - x116*x173 + x119 - x123 - x129 - x130 - x131 + x135 + x138*x173 + x138*x174 + x142 + x150 + x155 + x162 - x81 - x97 - x98 - x99));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[16] = x27*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(16)*x28 + x29*(x152 + x156 + x162 + x169));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[17] = x27*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(17)*x28 + x29*(x164 + x166 + x175));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[18] = x27*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(18)*x28 + x29*(x156 + x165 + x175));
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

// Generation Info: commit=557a856e999bdc7da25fdc690eab4179ea52272c
