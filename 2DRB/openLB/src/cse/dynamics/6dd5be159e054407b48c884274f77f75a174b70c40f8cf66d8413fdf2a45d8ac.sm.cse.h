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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::CBCoutsideDensity, momenta::CBCoutsideMomentum, momenta::BulkStress, momenta::DefineUSeparatelyTrace>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x64 = cell.template getFieldComponent<olb::fields::cbc::CBC_EDGE_VELOCITY>(0);
auto x65 = cell.template getFieldComponent<olb::fields::cbc::CBC_EDGE_VELOCITY>(1);
auto x66 = cell.template getFieldComponent<olb::fields::cbc::CBC_EDGE_VELOCITY>(2);
auto x75 = parameters.template get<descriptors::OMEGA>();
auto x67 = ((x64)*(x64));
auto x68 = ((x65)*(x65));
auto x69 = ((x66)*(x66));
auto x0 = V{1} - x75;
auto x1 = V{3};
auto x2 = x67 + x68 + x69;
auto x3 = V{0.5}*x1*x2;
auto x71 = V{1} - x3;
auto x4 = V{0.333333333333333};
auto x5 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{1};
auto x6 = cell[0]*x0 - x4*x75*(x5*(x3 + V{-1}) + V{1});
auto x7 = -V{1}*x64;
auto x8 = V{0.5}*((x1)*(x1));
auto x9 = V{0.0555555555555556};
auto x70 = -x9;
auto x10 = x5*x9;
auto x11 = cell[1]*x0 + x75*(x10*(x1*x7 + x71 + x8*((x7)*(x7))) + x70);
auto x12 = -V{1}*x65;
auto x13 = cell[2]*x0 + x75*(x10*(x1*x12 + x71 + x8*((x12)*(x12))) + x70);
auto x14 = -V{1}*x66;
auto x15 = cell[3]*x0 + x75*(x10*(x1*x14 + x71 + x8*((x14)*(x14))) + x70);
auto x16 = x12 + x7;
auto x17 = V{0.0277777777777778};
auto x18 = x17*x5;
auto x19 = cell[4]*x0 - x75*(x17 - x18*(x1*x16 + x71 + x8*((x16)*(x16))));
auto x20 = x65 + x7;
auto x21 = cell[5]*x0 - x75*(x17 - x18*(x1*x20 + x71 + x8*((x20)*(x20))));
auto x22 = x14 + x7;
auto x23 = cell[6]*x0 - x75*(x17 - x18*(x1*x22 + x71 + x8*((x22)*(x22))));
auto x24 = x66 + x7;
auto x25 = cell[7]*x0 - x75*(x17 - x18*(x1*x24 + x71 + x8*((x24)*(x24))));
auto x26 = x12 + x14;
auto x27 = cell[8]*x0 - x75*(x17 - x18*(x1*x26 + x71 + x8*((x26)*(x26))));
auto x28 = x12 + x66;
auto x29 = cell[9]*x0 - x75*(x17 - x18*(x1*x28 + x71 + x8*((x28)*(x28))));
auto x30 = cell[10]*x0 + x75*(x10*(x1*x64 + x67*x8 + x71) + x70);
auto x31 = cell[11]*x0 + x75*(x10*(x1*x65 + x68*x8 + x71) + x70);
auto x32 = cell[12]*x0 + x75*(x10*(x1*x66 + x69*x8 + x71) + x70);
auto x33 = x64 + x65;
auto x34 = cell[13]*x0 - x75*(x17 - x18*(x1*x33 + x71 + x8*((x33)*(x33))));
auto x35 = x12 + x64;
auto x36 = cell[14]*x0 - x75*(x17 - x18*(x1*x35 + x71 + x8*((x35)*(x35))));
auto x37 = x64 + x66;
auto x38 = cell[15]*x0 - x75*(x17 - x18*(x1*x37 + x71 + x8*((x37)*(x37))));
auto x39 = x14 + x64;
auto x40 = cell[16]*x0 - x75*(x17 - x18*(x1*x39 + x71 + x8*((x39)*(x39))));
auto x41 = x65 + x66;
auto x42 = cell[17]*x0 - x75*(x17 - x18*(x1*x41 + x71 + x8*((x41)*(x41))));
auto x43 = x14 + x65;
auto x44 = cell[18]*x0 - x75*(x17 - x18*(x1*x43 + x71 + x8*((x43)*(x43))));
cell[0] = x6;
cell[10] = x30;
cell[11] = x31;
cell[12] = x32;
cell[13] = x34;
cell[14] = x36;
cell[15] = x38;
cell[16] = x40;
cell[17] = x42;
cell[18] = x44;
cell[1] = x11;
cell[2] = x13;
cell[3] = x15;
cell[4] = x19;
cell[5] = x21;
cell[6] = x23;
cell[7] = x25;
cell[8] = x27;
cell[9] = x29;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x11 + x13 + x15 + x19 + x21 + x23 + x25 + x27 + x29 + x30 + x31 + x32 + x34 + x36 + x38 + x40 + x42 + x44 + x6 + V{1};
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x64;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x65;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x66;
return { x5, x2 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
