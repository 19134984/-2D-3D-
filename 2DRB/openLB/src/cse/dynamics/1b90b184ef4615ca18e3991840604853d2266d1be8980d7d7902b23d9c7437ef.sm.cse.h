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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::ParameterFromCell<collision::HYBRID, collision::LocalSmagorinskyEffectiveOmega<collision::HRR> >, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x101 = cell.template getFieldComponent<olb::collision::HYBRID>(0);
auto x103 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x104 = cell.template getFieldComponent<olb::descriptors::TENSOR>(0);
auto x105 = cell.template getFieldComponent<olb::descriptors::TENSOR>(1);
auto x106 = cell.template getFieldComponent<olb::descriptors::TENSOR>(2);
auto x107 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x108 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x109 = cell.template getFieldComponent<olb::descriptors::WMPOROSITY>(0);
auto x110 = cell.template getFieldComponent<olb::descriptors::WMVELOCITY>(0);
auto x111 = cell.template getFieldComponent<olb::descriptors::WMVELOCITY>(1);
auto x113 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x114 = parameters.template get<descriptors::OMEGA>();
auto x102 = cell[1] + cell[7];
auto x112 = cell[8] + x102;
auto x115 = cell[4] + x112;
auto x116 = cell[5] + cell[6];
auto x117 = cell[2] + x116;
auto x0 = V{3};
auto x123 = ((x0)*(x0));
auto x124 = V{1} / (x0);
auto x125 = V{1}*x124;
auto x1 = V{1} - x109;
auto x2 = x1*x110;
auto x3 = V{1} - x103;
auto x4 = x107*x3;
auto x5 = -V{1}*cell[1];
auto x118 = cell[7] + x5;
auto x6 = -V{1}*cell[3];
auto x7 = cell[0] + cell[3] + x115 + x117 + V{1};
auto x119 = V{1} / (x7);
auto x120 = x103*x119;
auto x8 = x109*(x120*(-V{1}*cell[2] + x116 + x118 + x6) + x4) + x2;
auto x121 = ((x8)*(x8));
auto x9 = x121;
auto x10 = x1*x111;
auto x11 = x108*x3;
auto x12 = -V{1}*cell[5];
auto x13 = x10 + x109*(x11 + x120*(-V{1}*cell[4] + x112 + x12 + x6));
auto x122 = ((x13)*(x13));
auto x14 = x122 + x9;
auto x15 = V{0.5}*x0*x14;
auto x138 = V{1} - x15;
auto x16 = V{0.444444444444444};
auto x17 = x123;
auto x18 = V{0.5}*x17;
auto x19 = x125;
auto x128 = V{1} - x19;
auto x20 = -x19;
auto x21 = V{1}/x114;
auto x22 = x19*(x7 + V{-1});
auto x126 = cell[3] - x22;
auto x23 = cell[5] + x115 - x122*x7 + x126;
auto x24 = x7*x8;
auto x25 = cell[3] + x118 + x12 - x13*x24;
auto x26 = x102 + x117 + x126 - x24*x8;
auto x27 = V{2}*V{1} / (x21 + util::sqrt(V{2.828428}*x119*x123*((x113)*(x113))*util::sqrt(((x23)*(x23)) + V{2}*((x25)*(x25)) + ((x26)*(x26))) + ((x21)*(x21))));
auto x127 = V{1} / (x27);
auto x29 = V{2}*x124*x127*x7*(x101 + V{-1});
auto x30 = x101*x23 + x106*x29;
auto x31 = x20*x30;
auto x32 = x101*x26 + x104*x29;
auto x33 = x20*x32;
auto x34 = V{1} - x27;
auto x35 = x16*(x18*x34*(x31 + x33) - x7*(x15 + V{-1}) + V{-1});
auto x36 = x0*x17;
auto x37 = V{0.166667}*x36;
auto x142 = x20*x37;
auto x38 = x128;
auto x39 = -x125;
auto x40 = -x39 + V{-1};
auto x41 = x13*x9;
auto x42 = x122*x8;
auto x43 = x41 - x42;
auto x136 = x37*x43;
auto x44 = V{0.5}*x36;
auto x45 = x41 + x42;
auto x137 = x44*x45;
auto x46 = -V{1}*x8;
auto x47 = x13 + x46;
auto x48 = V{0.0277777777777778};
auto x132 = -x48;
auto x135 = x34*x48;
auto x49 = x48*x7;
auto x50 = -V{1}*x38;
auto x133 = -x50;
auto x51 = x101*x25 + x105*x29;
auto x129 = V{2}*x51;
auto x130 = x17*x51;
auto x52 = x129*x8 + x13*x32;
auto x53 = x129*x13 + x30*x8;
auto x54 = x52 - x53;
auto x134 = x37*x54*(x133 + x38);
auto x55 = x52 + x53;
auto x131 = x44*x55;
auto x56 = x30*x38;
auto x57 = x32*x38;
auto x58 = x18*(x56 + x57);
auto x59 = -V{1}*x130 + x131*(x38 + x50) + x58;
auto x60 = x132 + x135*(x134 + x59) + x49*(x0*x47 + x136*(x38 - x40) + x137*(x38 + x40) + x138 + x18*((x47)*(x47)));
auto x61 = -x39;
auto x62 = V{0.111111111111111};
auto x139 = -x62;
auto x141 = x34*x62;
auto x63 = x62*x7;
auto x64 = -V{1}*x20;
auto x140 = x37*x64;
auto x65 = x44*x64;
auto x66 = x55*x65;
auto x67 = x18*(x31 + x57);
auto x68 = x139 + x141*(-x140*x54 + x66 + x67) + x63*(x0*x46 - x136*x61 + x137*x61 + x138 + x18*((x46)*(x46)));
auto x69 = -V{1}*x13;
auto x70 = x46 + x69;
auto x71 = x130 + x58;
auto x72 = x132 + x135*(V{2}*x131*x50 + x71) + x49*(x0*x70 - x136*(x133 + x40) + x137*(x40 + x50) + x138 + x18*((x70)*(x70)));
auto x73 = x140;
auto x74 = x18*(x33 + x56);
auto x75 = x139 + x141*(x54*x73 + x66 + x74) + x63*(x0*x69 + x138 + x18*((x69)*(x69)) + x43*x73 + x45*x65);
auto x76 = x128;
auto x77 = x69 + x8;
auto x78 = -(x135*(x134 - x59) + x48 - x49*(x0*x77 + x136*(x50 - x76) + x137*(x50 + x76) + x138 + x18*((x77)*(x77))));
auto x79 = -x142;
auto x80 = x20*x44;
auto x81 = x45*x80;
auto x143 = x138 + x81;
auto x82 = x55*x80;
auto x83 = x139 + x141*(x54*x79 + x67 + x82) + x63*(x0*x8 + x121*x18 + x143 + x43*x79);
auto x144 = x78 + x83;
auto x84 = 0;
auto x85 = V{2}*x38*x44;
auto x86 = x13 + x8;
auto x87 = x132 + x135*(x54*x84 + x55*x85 + x71) + x49*(x0*x86 + x138 + x18*((x86)*(x86)) + x43*x84 + x45*x85);
auto x88 = x142;
auto x89 = x139 + x141*(x54*x88 + x74 + x82) + x63*(x0*x13 + x122*x18 + x143 + x43*x88);
auto x145 = x60 + x89;
auto x90 = x144 + x145 + x35 + x68 + x72 + x75 + x87 + V{1};
auto x147 = x103/x90;
auto x91 = -V{1}*x72;
auto x146 = x87 + x91;
cell[0] = x35;
cell[1] = x60;
cell[2] = x68;
cell[3] = x72;
cell[4] = x75;
cell[5] = x78;
cell[6] = x83;
cell[7] = x87;
cell[8] = x89;
cell.template getFieldPointer<olb::descriptors::OMEGA>()[0] = x27;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x90;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x109*(x147*(x144 + x146 - V{1}*x60 - V{1}*x68) + x4) + x2;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x10 + x109*(x11 + x147*(x145 + x146 - V{1}*x75 - V{1}*x78));
return { x7, x14 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
