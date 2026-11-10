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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::FirstOrder, collision::BGK, forcing::ConservativePhaseField>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x33 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x34 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x35 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x36 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x37 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{1} - x37;
auto x1 = V{0.444444444444444};
auto x2 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + V{1};
auto x3 = cell[0]*x0 + x1*x37*(x2 + V{-1});
auto x4 = V{3};
auto x5 = -V{1}*x35;
auto x6 = V{0.0277777777777778};
auto x7 = x2*x6;
auto x8 = V{1} - V{0.5}*x37;
auto x9 = x4*x6;
auto x38 = x8*x9;
auto x10 = -V{1}*x33;
auto x11 = cell[1]*x0 - x37*(x6 - x7*(x4*(x36 + x5) + V{1})) + x38*(x10 + x34);
auto x12 = V{0.111111111111111};
auto x13 = x12*x2;
auto x14 = x12*x4;
auto x39 = x14*x8;
auto x15 = cell[2]*x0 + x10*x39 - x37*(x12 - x13*(x4*x5 + V{1}));
auto x16 = -V{1}*x36;
auto x17 = -V{1}*x34;
auto x18 = cell[3]*x0 - x37*(x6 - x7*(x4*(x16 + x5) + V{1})) + x38*(x10 + x17);
auto x19 = cell[4]*x0 + x17*x39 - x37*(x12 - x13*(x16*x4 + V{1}));
auto x20 = cell[5]*x0 - x37*(x6 - x7*(x4*(x16 + x35) + V{1})) + x38*(x17 + x33);
auto x21 = cell[6]*x0 + x33*x39 - x37*(x12 - x13*(x35*x4 + V{1}));
auto x22 = cell[7]*x0 - x37*(x6 - x7*(x4*(x35 + x36) + V{1})) + x38*(x33 + x34);
auto x23 = cell[8]*x0 + x34*x39 - x37*(x12 - x13*(x36*x4 + V{1}));
cell[0] = x3;
cell[1] = x11;
cell[2] = x15;
cell[3] = x18;
cell[4] = x19;
cell[5] = x20;
cell[6] = x21;
cell[7] = x22;
cell[8] = x23;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x11 + x15 + x18 + x19 + x20 + x21 + x22 + x23 + x3 + V{1};
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x35;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x36;
return { x2, ((x35)*(x35)) + ((x36)*(x36)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
