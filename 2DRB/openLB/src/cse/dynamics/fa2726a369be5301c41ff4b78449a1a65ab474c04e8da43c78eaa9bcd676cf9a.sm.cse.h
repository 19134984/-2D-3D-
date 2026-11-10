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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::CahnHilliardZerothOrder, collision::BGK, forcing::WellBalancedCahnHilliard>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x29 = cell.template getFieldComponent<olb::descriptors::CHEM_POTENTIAL>(0);
auto x30 = cell.template getFieldComponent<olb::descriptors::SOURCE>(0);
auto x31 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x32 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x33 = parameters.template get<descriptors::OMEGA>();
auto x34 = x33*(x29 + V{-1});
auto x0 = V{1} - x33;
auto x1 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + V{1};
auto x2 = V{0.444444444444444};
auto x3 = V{-2};
auto x4 = cell[0]*x0 + x2*x3*x30 - x33*(-x1 + x2 - x29*(x2 + V{-1}));
auto x5 = V{0.0277777777777778};
auto x6 = x34*x5;
auto x7 = V{1.5};
auto x8 = x30*x5*(x3 + V{2}*x7);
auto x35 = x6 + x8;
auto x9 = cell[1]*x0 + x35;
auto x10 = V{0.111111111111111};
auto x11 = x10*x34;
auto x12 = x10*x30*(x3 + x7);
auto x36 = x11 + x12;
auto x13 = cell[2]*x0 + x36;
auto x14 = cell[3]*x0 + x35;
auto x15 = cell[4]*x0 + x36;
auto x16 = cell[5]*x0 + x35;
auto x17 = cell[6]*x0 + x36;
auto x18 = cell[7]*x0 + x35;
auto x19 = cell[8]*x0 + x36;
cell[0] = x4;
cell[1] = x9;
cell[2] = x13;
cell[3] = x14;
cell[4] = x15;
cell[5] = x16;
cell[6] = x17;
cell[7] = x18;
cell[8] = x19;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x13 + x14 + x15 + x16 + x17 + x18 + x19 + x4 + x9 + V{1};
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x31;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x32;
return { x1, ((x31)*(x31)) + ((x32)*(x32)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
