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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::MovingPorousMomentumCombinationNoWM<momenta::BulkMomentum>, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x92 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x93 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x94 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x95 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x96 = parameters.template get<descriptors::OMEGA>();
auto x97 = cell[1] + cell[7];
auto x98 = cell[8] + x97;
auto x99 = cell[4] + x98;
auto x100 = cell[5] + cell[6];
auto x101 = cell[2] + x100;
auto x0 = V{3};
auto x107 = ((x0)*(x0));
auto x108 = V{1}/x0;
auto x1 = V{1} - x92;
auto x2 = x1*x93;
auto x3 = -V{1}*cell[1];
auto x102 = cell[7] + x3;
auto x4 = -V{1}*cell[3];
auto x5 = cell[0] + cell[3] + x101 + x99 + V{1};
auto x103 = V{1} / (x5);
auto x104 = x103*x92;
auto x6 = x104*(-V{1}*cell[2] + x100 + x102 + x4) + x2;
auto x105 = ((x6)*(x6));
auto x7 = x105;
auto x8 = x1*x94;
auto x9 = -V{1}*cell[5];
auto x10 = x104*(-V{1}*cell[4] + x4 + x9 + x98) + x8;
auto x106 = ((x10)*(x10));
auto x11 = x106 + x7;
auto x12 = V{0.5}*x0*x11;
auto x120 = V{1} - x12;
auto x13 = V{0.444444444444444};
auto x14 = x107;
auto x15 = V{0.5}*x14;
auto x16 = x108;
auto x110 = V{1} - x16;
auto x17 = -x16;
auto x18 = x16*(x5 + V{-1});
auto x109 = cell[3] - x18;
auto x19 = cell[5] - x106*x5 + x109 + x99;
auto x20 = x17*x19;
auto x21 = x5*x6;
auto x22 = x101 + x109 - x21*x6 + x97;
auto x23 = x17*x22;
auto x24 = V{1}/x96;
auto x25 = cell[3] - x10*x21 + x102 + x9;
auto x111 = V{2}*x25;
auto x112 = x14*x25;
auto x26 = V{1} - V{2}*V{1} / (x24 + util::sqrt(V{2.828428}*x103*x107*((x95)*(x95))*util::sqrt(((x19)*(x19)) + ((x22)*(x22)) + V{2}*((x25)*(x25))) + ((x24)*(x24))));
auto x27 = x13*(x15*x26*(x20 + x23) - x5*(x12 + V{-1}) + V{-1});
auto x28 = x0*x14;
auto x29 = V{0.166667}*x28;
auto x124 = x17*x29;
auto x30 = x110;
auto x31 = -x108;
auto x32 = -x31 + V{-1};
auto x33 = x10*x7;
auto x34 = x106*x6;
auto x35 = x33 - x34;
auto x118 = x29*x35;
auto x36 = V{0.5}*x28;
auto x37 = x33 + x34;
auto x119 = x36*x37;
auto x38 = -V{1}*x6;
auto x39 = x10 + x38;
auto x40 = V{0.0277777777777778};
auto x114 = -x40;
auto x117 = x26*x40;
auto x41 = x40*x5;
auto x42 = -V{1}*x30;
auto x115 = -x42;
auto x43 = x10*x22 + x111*x6;
auto x44 = x10*x111 + x19*x6;
auto x45 = x43 - x44;
auto x116 = x29*x45*(x115 + x30);
auto x46 = x43 + x44;
auto x113 = x36*x46;
auto x47 = x19*x30;
auto x48 = x22*x30;
auto x49 = x15*(x47 + x48);
auto x50 = -V{1}*x112 + x113*(x30 + x42) + x49;
auto x51 = x114 + x117*(x116 + x50) + x41*(x0*x39 + x118*(x30 - x32) + x119*(x30 + x32) + x120 + x15*((x39)*(x39)));
auto x52 = -x31;
auto x53 = V{0.111111111111111};
auto x121 = -x53;
auto x123 = x26*x53;
auto x54 = x5*x53;
auto x55 = -V{1}*x17;
auto x122 = x29*x55;
auto x56 = x36*x55;
auto x57 = x46*x56;
auto x58 = x15*(x20 + x48);
auto x59 = x121 + x123*(-x122*x45 + x57 + x58) + x54*(x0*x38 - x118*x52 + x119*x52 + x120 + x15*((x38)*(x38)));
auto x60 = -V{1}*x10;
auto x61 = x38 + x60;
auto x62 = x112 + x49;
auto x63 = x114 + x117*(V{2}*x113*x42 + x62) + x41*(x0*x61 - x118*(x115 + x32) + x119*(x32 + x42) + x120 + x15*((x61)*(x61)));
auto x64 = x122;
auto x65 = x15*(x23 + x47);
auto x66 = x121 + x123*(x45*x64 + x57 + x65) + x54*(x0*x60 + x120 + x15*((x60)*(x60)) + x35*x64 + x37*x56);
auto x67 = x110;
auto x68 = x6 + x60;
auto x69 = -(x117*(x116 - x50) + x40 - x41*(x0*x68 + x118*(x42 - x67) + x119*(x42 + x67) + x120 + x15*((x68)*(x68))));
auto x70 = -x124;
auto x71 = x17*x36;
auto x72 = x37*x71;
auto x125 = x120 + x72;
auto x73 = x46*x71;
auto x74 = x121 + x123*(x45*x70 + x58 + x73) + x54*(x0*x6 + x105*x15 + x125 + x35*x70);
auto x126 = x69 + x74;
auto x75 = 0;
auto x76 = V{2}*x30*x36;
auto x77 = x10 + x6;
auto x78 = x114 + x117*(x45*x75 + x46*x76 + x62) + x41*(x0*x77 + x120 + x15*((x77)*(x77)) + x35*x75 + x37*x76);
auto x79 = x124;
auto x80 = x121 + x123*(x45*x79 + x65 + x73) + x54*(x0*x10 + x106*x15 + x125 + x35*x79);
auto x127 = x51 + x80;
auto x81 = x126 + x127 + x27 + x59 + x63 + x66 + x78 + V{1};
auto x129 = x92/x81;
auto x82 = -V{1}*x63;
auto x128 = x78 + x82;
cell[0] = x27;
cell[1] = x51;
cell[2] = x59;
cell[3] = x63;
cell[4] = x66;
cell[5] = x69;
cell[6] = x74;
cell[7] = x78;
cell[8] = x80;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x81;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x129*(x126 + x128 - V{1}*x51 - V{1}*x59) + x2;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x129*(x127 + x128 - V{1}*x66 - V{1}*x69) + x8;
return { x5, x11 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
