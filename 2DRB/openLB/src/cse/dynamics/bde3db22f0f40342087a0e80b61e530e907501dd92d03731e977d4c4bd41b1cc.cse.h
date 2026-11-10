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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::fsi::HLBM>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x50 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x51 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x52 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x53 = parameters.template get<descriptors::OMEGA>();
auto x54 = cell[5] + cell[6];
auto x55 = cell[1] + cell[8];
auto x0 = V{1} - x53;
auto x1 = V{3};
auto x2 = V{0.5}*x1;
auto x3 = -V{1}*cell[3];
auto x57 = cell[7] + x3;
auto x4 = cell[0] + cell[2] + cell[3] + cell[4] + cell[7] + x54 + x55 + V{1};
auto x56 = V{1} / (x4);
auto x5 = x56*(-V{1}*cell[1] - V{1}*cell[2] + x54 + x57);
auto x58 = ((x5)*(x5));
auto x6 = x56*(-V{1}*cell[4] - V{1}*cell[5] + x55 + x57);
auto x59 = ((x6)*(x6));
auto x7 = x58 + x59;
auto x8 = x2*x7;
auto x63 = V{1} - x8;
auto x9 = V{0.444444444444444};
auto x10 = x4*x9;
auto x11 = -x10*(x8 + V{-1}) - x9;
auto x12 = V{1} - x50;
auto x13 = -x12*(x5 - x51) + x5;
auto x60 = ((x13)*(x13));
auto x14 = x12*(x52 - x6) + x6;
auto x61 = ((x14)*(x14));
auto x15 = x2*(x60 + x61);
auto x65 = V{1} - x15;
auto x16 = -V{1}*x5;
auto x17 = x16 + x6;
auto x18 = V{0.5}*((x1)*(x1));
auto x19 = V{0.0277777777777778};
auto x62 = -x19;
auto x20 = x19*x4;
auto x21 = x20*(x1*x17 + x18*((x17)*(x17)) + x63) + x62;
auto x22 = -V{1}*x13;
auto x23 = x14 + x22;
auto x24 = V{0.111111111111111};
auto x64 = -x24;
auto x25 = x24*x4;
auto x26 = x25*(x1*x16 + x18*((x16)*(x16)) + x63) + x64;
auto x27 = -V{1}*x6;
auto x28 = x16 + x27;
auto x29 = x20*(x1*x28 + x18*((x28)*(x28)) + x63) + x62;
auto x30 = -V{1}*x14;
auto x31 = x22 + x30;
auto x32 = x25*(x1*x27 + x18*((x27)*(x27)) + x63) + x64;
auto x33 = x27 + x5;
auto x34 = x20*(x1*x33 + x18*((x33)*(x33)) + x63) + x62;
auto x35 = x13 + x30;
auto x36 = x25*(x1*x5 + x18*x58 + x63) + x64;
auto x37 = x5 + x6;
auto x38 = x20*(x1*x37 + x18*((x37)*(x37)) + x63) + x62;
auto x39 = x13 + x14;
auto x40 = x25*(x1*x6 + x18*x59 + x63) + x64;
auto x41 = cell[0]*x0 - x10*(x15 + V{-1}) + x11*x53 - x11 - x9;
auto x42 = cell[1]*x0 + x20*(x1*x23 + x18*((x23)*(x23)) + x65) + x21*x53 - x21 + x62;
auto x43 = cell[2]*x0 + x25*(x1*x22 + x18*((x22)*(x22)) + x65) + x26*x53 - x26 + x64;
auto x44 = cell[3]*x0 + x20*(x1*x31 + x18*((x31)*(x31)) + x65) + x29*x53 - x29 + x62;
auto x45 = cell[4]*x0 + x25*(x1*x30 + x18*((x30)*(x30)) + x65) + x32*x53 - x32 + x64;
auto x46 = cell[5]*x0 + x20*(x1*x35 + x18*((x35)*(x35)) + x65) + x34*x53 - x34 + x62;
auto x47 = cell[6]*x0 + x25*(x1*x13 + x18*x60 + x65) + x36*x53 - x36 + x64;
auto x48 = cell[7]*x0 + x20*(x1*x39 + x18*((x39)*(x39)) + x65) + x38*x53 - x38 + x62;
auto x49 = cell[8]*x0 + x25*(x1*x14 + x18*x61 + x65) + x40*x53 - x40 + x64;
cell[0] = x41;
cell[1] = x42;
cell[2] = x43;
cell[3] = x44;
cell[4] = x45;
cell[5] = x46;
cell[6] = x47;
cell[7] = x48;
cell[8] = x49;
return { x4, x7 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
