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
struct CSE<AdvectionDiffusionCornerDynamics2D<T, descriptors::D2Q5<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q5<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, 1, -1>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x25 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x26 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x27 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x28 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x29 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x30 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{0.333333333333333};
auto x1 = V{0.166666666666667};
auto x32 = V{2}*x1;
auto x2 = x1*x27;
auto x3 = V{3};
auto x31 = x26*x3;
auto x33 = x25*x3;
auto x4 = x31 + V{1};
auto x5 = V{1} - V{1}*x31;
auto x6 = -cell[2] + x2*x4 + x2*x5 - x32;
auto x7 = V{1} - V{1}*x33;
auto x8 = x33 + V{1};
auto x9 = -cell[3] + x2*x7 + x2*x8 - x32;
auto x10 = cell[0] + cell[2] + cell[3] + x6 + x9 + V{1};
auto x11 = x1*x10;
auto x12 = -x1 + x11*x7;
auto x13 = V{1} - x30;
auto x14 = x1*x3;
auto x34 = x13*x14;
auto x15 = -x1 + x11*x8;
auto x16 = cell[3] + V{1}*x12 - x15 - V{1}*x9;
auto x35 = x16*x34;
auto x17 = -x1 + x11*x5;
auto x18 = -x1 + x11*x4;
auto x19 = -V{1}*cell[2] + V{1}*x17 - x18 + x6;
auto x36 = x19*x34;
cell[0] = x0*(x10 + V{-1});
cell[1] = x12 - V{1}*x35;
cell[2] = x17 - V{1}*x36;
cell[3] = x15 + x35;
cell[4] = x18 + x36;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x27;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x28;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x29;
return { x10, ((x25)*(x25)) + ((x26)*(x26)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
