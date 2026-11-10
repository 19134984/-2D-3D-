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
struct CSE<Dual<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x84 = cell.template getFieldComponent<olb::opti::DJDF>(2);
auto x87 = cell.template getFieldComponent<olb::opti::DJDF>(5);
auto x91 = cell.template getFieldComponent<olb::opti::F>(0);
auto x96 = cell.template getFieldComponent<olb::opti::F>(5);
auto x88 = cell.template getFieldComponent<olb::opti::DJDF>(6);
auto x82 = cell.template getFieldComponent<olb::opti::DJDF>(0);
auto x94 = cell.template getFieldComponent<olb::opti::F>(3);
auto x83 = cell.template getFieldComponent<olb::opti::DJDF>(1);
auto x81 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x97 = cell.template getFieldComponent<olb::opti::F>(6);
auto x89 = cell.template getFieldComponent<olb::opti::DJDF>(7);
auto x98 = cell.template getFieldComponent<olb::opti::F>(7);
auto x95 = cell.template getFieldComponent<olb::opti::F>(4);
auto x85 = cell.template getFieldComponent<olb::opti::DJDF>(3);
auto x93 = cell.template getFieldComponent<olb::opti::F>(2);
auto x99 = cell.template getFieldComponent<olb::opti::F>(8);
auto x92 = cell.template getFieldComponent<olb::opti::F>(1);
auto x100 = parameters.template get<descriptors::OMEGA>();
auto x90 = cell.template getFieldComponent<olb::opti::DJDF>(8);
auto x86 = cell.template getFieldComponent<olb::opti::DJDF>(4);
auto x102 = x92 + x99;
auto x103 = x96 + x97;
auto x0 = V{0.111111111111111};
auto x1 = -V{1}*x94;
auto x101 = x1 + x98;
auto x2 = x101 + x102 - V{1}*x95 - V{1}*x96;
auto x3 = x102 + x103 + x91 + x93 + x94 + x95 + x98 + V{1};
auto x104 = V{1} / (x3);
auto x105 = x104*x81;
auto x4 = x105*x2;
auto x5 = V{3};
auto x6 = V{0.5}*((x5)*(x5));
auto x106 = x4*x6;
auto x7 = x106;
auto x8 = V{0.5}*x5;
auto x9 = x101 + x103 - V{1}*x92 - V{1}*x93;
auto x10 = x105*x9;
auto x108 = x10*x6;
auto x11 = x8*(((x10)*(x10)) + ((x4)*(x4)));
auto x114 = V{1} - x11;
auto x12 = cell[8]*x100;
auto x13 = V{0.0277777777777778};
auto x14 = x10 + x4;
auto x107 = x14*x6;
auto x15 = x107;
auto x16 = cell[7]*x100;
auto x17 = x108;
auto x18 = cell[6]*x100;
auto x19 = -V{1}*x4;
auto x110 = x19*x6;
auto x20 = x10 + x19;
auto x109 = x20*x6;
auto x21 = x109;
auto x22 = cell[5]*x100;
auto x23 = x110;
auto x24 = cell[4]*x100;
auto x25 = -V{1}*x10;
auto x112 = x25*x6;
auto x26 = x19 + x25;
auto x111 = x26*x6;
auto x27 = x111;
auto x28 = cell[3]*x100;
auto x29 = x112;
auto x30 = cell[2]*x100;
auto x31 = x25 + x4;
auto x113 = x31*x6;
auto x32 = x113;
auto x33 = cell[1]*x100;
auto x34 = V{0.444444444444444};
auto x35 = cell[0]*x100;
auto x115 = x34*x35;
auto x36 = x0*x12*(x114 + x4*x5 + x4*x7) + x0*x18*(x10*x17 + x10*x5 + x114) + x0*x24*(x114 + x19*x23 + x19*x5) + x0*x30*(x114 + x25*x29 + x25*x5) - x115*(x11 + V{-1}) + x13*x16*(x114 + x14*x15 + x14*x5) + x13*x22*(x114 + x20*x21 + x20*x5) + x13*x28*(x114 + x26*x27 + x26*x5) + x13*x33*(x114 + x31*x32 + x31*x5);
auto x37 = V{1} - x100;
auto x38 = x0*x3;
auto x39 = x12*x38;
auto x40 = x13*x3;
auto x41 = x16*x40;
auto x42 = x18*x38;
auto x43 = x22*x40;
auto x44 = x24*x38;
auto x116 = V{1}*x44;
auto x45 = x28*x40;
auto x46 = x30*x38;
auto x118 = V{1}*x46;
auto x47 = x33*x40;
auto x48 = -V{1}*x8*(x115*x3 + x39 + x41 + x42 + x43 + x44 + x45 + x46 + x47);
auto x49 = x4*x48;
auto x50 = x41*(x107 + x15 + x5);
auto x51 = x43*(x109 + x21 + x5);
auto x52 = -V{1}*x45*(x111 + x27 + x5);
auto x117 = x50 + x52;
auto x53 = x47*(x113 + x32 + x5);
auto x54 = x81*(x106*x39 - x110*x116 - x116*x23 - x116*x5 + x117 + x39*x5 + x39*x7 + V{2}*x49 - V{1}*x51 + x53);
auto x55 = x10*x48;
auto x56 = x81*(x108*x42 - x112*x118 + x117 - x118*x29 - x118*x5 + x17*x42 + x42*x5 + x51 - V{1}*x53 + V{2}*x55);
auto x57 = -V{1}*(x2*x54 + x56*x9)/((x3)*(x3));
auto x119 = x36 + x57;
auto x58 = cell[0]*x37 + x119 - x82;
auto x59 = V{1}*x104;
auto x60 = x56*x59;
auto x122 = x119 + x60;
auto x61 = -V{1}*x60;
auto x120 = x119 + x61;
auto x62 = x54*x59;
auto x63 = cell[1]*x37 + x120 + x62 - x83;
auto x64 = cell[2]*x37 + x120 - x84;
auto x65 = -V{1}*x62;
auto x121 = x119 + x65;
auto x66 = cell[3]*x37 + x120 + x65 - x85;
auto x67 = cell[4]*x37 + x121 - x86;
auto x68 = cell[5]*x37 + x121 + x60 - x87;
auto x69 = cell[6]*x37 + x122 - x88;
auto x70 = cell[7]*x37 + x122 + x62 - x89;
auto x71 = cell[8]*x37 + x119 + x62 - x90;
cell[0] = x58;
cell[1] = x63;
cell[2] = x64;
cell[3] = x66;
cell[4] = x67;
cell[5] = x68;
cell[6] = x69;
cell[7] = x70;
cell[8] = x71;
return { x58 + V{1}, ((x58)*(x58)) + ((x63)*(x63)) + ((x64)*(x64)) + ((x66)*(x66)) + ((x67)*(x67)) + ((x68)*(x68)) + ((x69)*(x69)) + ((x70)*(x70)) + ((x71)*(x71)) };
}
};

}

}

// Generation Info: commit=4ea71c45fdf920fe461e9926f71dc63f7fa1062d
