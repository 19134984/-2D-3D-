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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::MovingPorousMomentumCombinationNoWM<momenta::BulkMomentum>, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x79 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x80 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x81 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x82 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x83 = parameters.template get<descriptors::OMEGA>();
auto x84 = cell[1] + cell[7];
auto x85 = cell[8] + x84;
auto x86 = cell[4] + x85;
auto x87 = cell[5] + cell[6];
auto x88 = cell[2] + x87;
auto x0 = V{3};
auto x94 = ((x0)*(x0));
auto x95 = V{1}/x0;
auto x1 = V{1} - x79;
auto x2 = -V{1}*cell[1];
auto x89 = cell[7] + x2;
auto x3 = -V{1}*cell[3];
auto x4 = cell[0] + cell[3] + x86 + x88 + V{1};
auto x90 = V{1} / (x4);
auto x91 = x79*x90;
auto x5 = x1*x80 + x91*(-V{1}*cell[2] + x3 + x87 + x89);
auto x92 = ((x5)*(x5));
auto x6 = x92;
auto x7 = -V{1}*cell[5];
auto x8 = x1*x81 + x91*(-V{1}*cell[4] + x3 + x7 + x85);
auto x93 = ((x8)*(x8));
auto x9 = x6 + x93;
auto x10 = V{0.5}*x0*x9;
auto x109 = V{1} - x10;
auto x11 = V{0.444444444444444};
auto x12 = x94;
auto x13 = V{0.5}*x12;
auto x14 = x95;
auto x97 = V{1} - x14;
auto x15 = -x14;
auto x16 = x14*(x4 + V{-1});
auto x96 = cell[3] - x16;
auto x17 = cell[5] - x4*x93 + x86 + x96;
auto x18 = x15*x17;
auto x19 = x4*x5;
auto x20 = -x19*x5 + x84 + x88 + x96;
auto x21 = x15*x20;
auto x22 = V{1}/x83;
auto x23 = cell[3] - x19*x8 + x7 + x89;
auto x98 = V{2}*x23;
auto x99 = x12*x23;
auto x24 = V{1} - V{2}*V{1} / (x22 + util::sqrt(V{2.828428}*x90*x94*((x82)*(x82))*util::sqrt(((x17)*(x17)) + ((x20)*(x20)) + V{2}*((x23)*(x23))) + ((x22)*(x22))));
auto x25 = x0*x12;
auto x26 = V{0.166667}*x25;
auto x102 = x15*x26;
auto x27 = x97;
auto x28 = -x95;
auto x29 = -x28 + V{-1};
auto x30 = x6*x8;
auto x31 = x5*x93;
auto x32 = x30 - x31;
auto x107 = x26*x32;
auto x33 = V{0.5}*x25;
auto x34 = x30 + x31;
auto x108 = x33*x34;
auto x35 = -V{1}*x5;
auto x36 = x35 + x8;
auto x37 = V{0.0277777777777778};
auto x103 = -x37;
auto x106 = x24*x37;
auto x38 = x37*x4;
auto x39 = -V{1}*x27;
auto x104 = -x39;
auto x40 = x20*x8 + x5*x98;
auto x41 = x17*x5 + x8*x98;
auto x42 = x40 - x41;
auto x105 = x26*x42*(x104 + x27);
auto x43 = x40 + x41;
auto x100 = x33*x43;
auto x44 = x17*x27;
auto x45 = x20*x27;
auto x46 = x13*(x44 + x45);
auto x47 = x100*(x27 + x39) + x46 - V{1}*x99;
auto x48 = -x28;
auto x49 = V{0.111111111111111};
auto x110 = -x49;
auto x111 = x24*x49;
auto x50 = x4*x49;
auto x51 = -V{1}*x15;
auto x101 = x26*x51;
auto x52 = x33*x51;
auto x53 = x43*x52;
auto x54 = x13*(x18 + x45);
auto x55 = -V{1}*x8;
auto x56 = x35 + x55;
auto x57 = x46 + x99;
auto x58 = x101;
auto x59 = x13*(x21 + x44);
auto x60 = x97;
auto x61 = x5 + x55;
auto x62 = -x102;
auto x63 = x15*x33;
auto x64 = x34*x63;
auto x112 = x109 + x64;
auto x65 = x43*x63;
auto x66 = 0;
auto x67 = V{2}*x27*x33;
auto x68 = x5 + x8;
auto x69 = x102;
cell[0] = x11*(x13*x24*(x18 + x21) - x4*(x10 + V{-1}) + V{-1});
cell[1] = x103 + x106*(x105 + x47) + x38*(x0*x36 + x107*(x27 - x29) + x108*(x27 + x29) + x109 + x13*((x36)*(x36)));
cell[2] = x110 + x111*(-x101*x42 + x53 + x54) + x50*(x0*x35 - x107*x48 + x108*x48 + x109 + x13*((x35)*(x35)));
cell[3] = x103 + x106*(V{2}*x100*x39 + x57) + x38*(x0*x56 - x107*(x104 + x29) + x108*(x29 + x39) + x109 + x13*((x56)*(x56)));
cell[4] = x110 + x111*(x42*x58 + x53 + x59) + x50*(x0*x55 + x109 + x13*((x55)*(x55)) + x32*x58 + x34*x52);
cell[5] = -(x106*(x105 - x47) + x37 - x38*(x0*x61 + x107*(x39 - x60) + x108*(x39 + x60) + x109 + x13*((x61)*(x61))));
cell[6] = x110 + x111*(x42*x62 + x54 + x65) + x50*(x0*x5 + x112 + x13*x92 + x32*x62);
cell[7] = x103 + x106*(x42*x66 + x43*x67 + x57) + x38*(x0*x68 + x109 + x13*((x68)*(x68)) + x32*x66 + x34*x67);
cell[8] = x110 + x111*(x42*x69 + x59 + x65) + x50*(x0*x8 + x112 + x13*x93 + x32*x69);
return { x4, x9 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
