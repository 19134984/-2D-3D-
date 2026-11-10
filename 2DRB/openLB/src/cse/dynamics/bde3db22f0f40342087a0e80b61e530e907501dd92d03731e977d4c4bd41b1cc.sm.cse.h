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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::fsi::HLBM>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x64 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x65 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x66 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x67 = parameters.template get<descriptors::OMEGA>();
auto x68 = cell[5] + cell[6];
auto x69 = cell[1] + cell[8];
auto x74 = -x65;
auto x75 = -x66;
auto x0 = V{1} - x67;
auto x1 = V{3};
auto x2 = V{0.5}*x1;
auto x3 = -V{1}*cell[3];
auto x71 = cell[7] + x3;
auto x4 = cell[0] + cell[2] + cell[3] + cell[4] + cell[7] + x68 + x69 + V{1};
auto x70 = V{1} / (x4);
auto x5 = x70*(-V{1}*cell[1] - V{1}*cell[2] + x68 + x71);
auto x72 = ((x5)*(x5));
auto x6 = x70*(-V{1}*cell[4] - V{1}*cell[5] + x69 + x71);
auto x73 = ((x6)*(x6));
auto x7 = x72 + x73;
auto x8 = x2*x7;
auto x79 = V{1} - x8;
auto x9 = V{0.444444444444444};
auto x10 = x4*x9;
auto x11 = -x10*(x8 + V{-1}) - x9;
auto x12 = V{1} - x64;
auto x13 = -x12*(x5 + x74) + x5;
auto x76 = ((x13)*(x13));
auto x14 = -x12*(x6 + x75) + x6;
auto x77 = ((x14)*(x14));
auto x15 = x2*(x76 + x77);
auto x80 = V{1} - x15;
auto x16 = cell[0]*x0 - x10*(x15 + V{-1}) + x11*x67 - x11 - x9;
auto x17 = -V{1}*x5;
auto x18 = x17 + x6;
auto x19 = V{0.5}*((x1)*(x1));
auto x20 = V{0.0277777777777778};
auto x78 = -x20;
auto x21 = x20*x4;
auto x22 = x21*(x1*x18 + x19*((x18)*(x18)) + x79) + x78;
auto x23 = -V{1}*x13;
auto x24 = x14 + x23;
auto x25 = cell[1]*x0 + x21*(x1*x24 + x19*((x24)*(x24)) + x80) + x22*x67 - x22 + x78;
auto x26 = V{0.111111111111111};
auto x81 = -x26;
auto x27 = x26*x4;
auto x28 = x27*(x1*x17 + x19*((x17)*(x17)) + x79) + x81;
auto x29 = cell[2]*x0 + x27*(x1*x23 + x19*((x23)*(x23)) + x80) + x28*x67 - x28 + x81;
auto x30 = -V{1}*x6;
auto x31 = x17 + x30;
auto x32 = x21*(x1*x31 + x19*((x31)*(x31)) + x79) + x78;
auto x33 = -V{1}*x14;
auto x34 = x23 + x33;
auto x35 = cell[3]*x0 + x21*(x1*x34 + x19*((x34)*(x34)) + x80) + x32*x67 - x32 + x78;
auto x36 = x27*(x1*x30 + x19*((x30)*(x30)) + x79) + x81;
auto x37 = cell[4]*x0 + x27*(x1*x33 + x19*((x33)*(x33)) + x80) + x36*x67 - x36 + x81;
auto x38 = x30 + x5;
auto x39 = x21*(x1*x38 + x19*((x38)*(x38)) + x79) + x78;
auto x40 = x13 + x33;
auto x41 = cell[5]*x0 + x21*(x1*x40 + x19*((x40)*(x40)) + x80) + x39*x67 - x39 + x78;
auto x42 = x27*(x1*x5 + x19*x72 + x79) + x81;
auto x43 = cell[6]*x0 + x27*(x1*x13 + x19*x76 + x80) + x42*x67 - x42 + x81;
auto x82 = x41 + x43;
auto x44 = x5 + x6;
auto x45 = x21*(x1*x44 + x19*((x44)*(x44)) + x79) + x78;
auto x46 = x13 + x14;
auto x47 = cell[7]*x0 + x21*(x1*x46 + x19*((x46)*(x46)) + x80) + x45*x67 - x45 + x78;
auto x48 = x27*(x1*x6 + x19*x73 + x79) + x81;
auto x49 = cell[8]*x0 + x27*(x1*x14 + x19*x77 + x80) + x48*x67 - x48 + x81;
auto x83 = x25 + x49;
auto x50 = x16 + x29 + x35 + x37 + x47 + x82 + x83 + V{1};
auto x84 = V{1} / (x50);
auto x51 = -V{1}*x35;
auto x85 = x47 + x51;
auto x52 = x84*(-V{1}*x25 - V{1}*x29 + x82 + x85);
auto x53 = V{0.5}*x12*x50;
auto x54 = x84*(-V{1}*x37 - V{1}*x41 + x83 + x85);
cell[0] = x16;
cell[1] = x25;
cell[2] = x29;
cell[3] = x35;
cell[4] = x37;
cell[5] = x41;
cell[6] = x43;
cell[7] = x47;
cell[8] = x49;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x50;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x52 - x53*(x52 + x74);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = -x53*(x54 + x75) + x54;
return { x4, x7 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
