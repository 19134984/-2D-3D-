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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::ParameterFromCell<collision::LES::SMAGORINSKY, collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB> > >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x78 = cell.template getFieldComponent<olb::collision::LES::SMAGORINSKY>(0);
auto x80 = parameters.template get<descriptors::OMEGA>();
auto x79 = cell[1] + cell[7];
auto x81 = cell[8] + x79;
auto x82 = cell[4] + x81;
auto x83 = cell[5] + cell[6];
auto x84 = cell[2] + x83;
auto x0 = V{3};
auto x89 = ((x0)*(x0));
auto x90 = V{1}/x0;
auto x1 = -V{1}*cell[1];
auto x86 = cell[7] + x1;
auto x2 = -V{1}*cell[3];
auto x3 = cell[0] + cell[3] + x82 + x84 + V{1};
auto x85 = V{1} / (x3);
auto x4 = x85*(-V{1}*cell[2] + x2 + x83 + x86);
auto x87 = ((x4)*(x4));
auto x5 = x87;
auto x6 = -V{1}*cell[5];
auto x7 = x85*(-V{1}*cell[4] + x2 + x6 + x81);
auto x88 = ((x7)*(x7));
auto x8 = x5 + x88;
auto x9 = V{0.5}*x0*x8;
auto x104 = V{1} - x9;
auto x10 = V{0.444444444444444};
auto x11 = x89;
auto x12 = V{0.5}*x11;
auto x13 = x90;
auto x92 = V{1} - x13;
auto x14 = -x13;
auto x15 = x13*(x3 + V{-1});
auto x91 = cell[3] - x15;
auto x16 = cell[5] - x3*x88 + x82 + x91;
auto x17 = x14*x16;
auto x18 = x3*x4;
auto x19 = -x18*x4 + x79 + x84 + x91;
auto x20 = x14*x19;
auto x21 = V{1}/x80;
auto x22 = cell[3] - x18*x7 + x6 + x86;
auto x93 = V{2}*x22;
auto x94 = x11*x22;
auto x23 = V{1} - V{2}*V{1} / (x21 + util::sqrt(V{2.828428}*x85*x89*((x78)*(x78))*util::sqrt(((x16)*(x16)) + ((x19)*(x19)) + V{2}*((x22)*(x22))) + ((x21)*(x21))));
auto x24 = x0*x11;
auto x25 = V{0.166667}*x24;
auto x97 = x14*x25;
auto x26 = x92;
auto x27 = -x90;
auto x28 = -x27 + V{-1};
auto x29 = x5*x7;
auto x30 = x4*x88;
auto x31 = x29 - x30;
auto x102 = x25*x31;
auto x32 = V{0.5}*x24;
auto x33 = x29 + x30;
auto x103 = x32*x33;
auto x34 = -V{1}*x4;
auto x35 = x34 + x7;
auto x36 = V{0.0277777777777778};
auto x98 = -x36;
auto x101 = x23*x36;
auto x37 = x3*x36;
auto x38 = -V{1}*x26;
auto x99 = -x38;
auto x39 = x19*x7 + x4*x93;
auto x40 = x16*x4 + x7*x93;
auto x41 = x39 - x40;
auto x100 = x25*x41*(x26 + x99);
auto x42 = x39 + x40;
auto x95 = x32*x42;
auto x43 = x16*x26;
auto x44 = x19*x26;
auto x45 = x12*(x43 + x44);
auto x46 = x45 - V{1}*x94 + x95*(x26 + x38);
auto x47 = -x27;
auto x48 = V{0.111111111111111};
auto x105 = -x48;
auto x106 = x23*x48;
auto x49 = x3*x48;
auto x50 = -V{1}*x14;
auto x96 = x25*x50;
auto x51 = x32*x50;
auto x52 = x42*x51;
auto x53 = x12*(x17 + x44);
auto x54 = -V{1}*x7;
auto x55 = x34 + x54;
auto x56 = x45 + x94;
auto x57 = x96;
auto x58 = x12*(x20 + x43);
auto x59 = x92;
auto x60 = x4 + x54;
auto x61 = -x97;
auto x62 = x14*x32;
auto x63 = x33*x62;
auto x107 = x104 + x63;
auto x64 = x42*x62;
auto x65 = 0;
auto x66 = V{2}*x26*x32;
auto x67 = x4 + x7;
auto x68 = x97;
cell[0] = x10*(x12*x23*(x17 + x20) - x3*(x9 + V{-1}) + V{-1});
cell[1] = x101*(x100 + x46) + x37*(x0*x35 + x102*(x26 - x28) + x103*(x26 + x28) + x104 + x12*((x35)*(x35))) + x98;
cell[2] = x105 + x106*(-x41*x96 + x52 + x53) + x49*(x0*x34 - x102*x47 + x103*x47 + x104 + x12*((x34)*(x34)));
cell[3] = x101*(V{2}*x38*x95 + x56) + x37*(x0*x55 - x102*(x28 + x99) + x103*(x28 + x38) + x104 + x12*((x55)*(x55))) + x98;
cell[4] = x105 + x106*(x41*x57 + x52 + x58) + x49*(x0*x54 + x104 + x12*((x54)*(x54)) + x31*x57 + x33*x51);
cell[5] = -(x101*(x100 - x46) + x36 - x37*(x0*x60 + x102*(x38 - x59) + x103*(x38 + x59) + x104 + x12*((x60)*(x60))));
cell[6] = x105 + x106*(x41*x61 + x53 + x64) + x49*(x0*x4 + x107 + x12*x87 + x31*x61);
cell[7] = x101*(x41*x65 + x42*x66 + x56) + x37*(x0*x67 + x104 + x12*((x67)*(x67)) + x31*x65 + x33*x66) + x98;
cell[8] = x105 + x106*(x41*x68 + x58 + x64) + x49*(x0*x7 + x107 + x12*x88 + x31*x68);
return { x3, x8 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
