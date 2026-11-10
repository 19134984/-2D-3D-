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
struct CSE<AdvectionDiffusionBoundariesDynamics<T, descriptors::D2Q5<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q5<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, 1, -1>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x17 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x18 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x19 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x20 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x21 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x22 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{0.333333333333333};
auto x1 = V{3};
auto x23 = x1*x17;
auto x25 = x1*x18;
auto x2 = V{0.166666666666667};
auto x3 = x19*x2;
auto x24 = V{1}*x3;
auto x4 = -x2 - x24*(x23 + V{-1});
auto x5 = V{1} - x22;
auto x6 = x1*x2;
auto x26 = x5*x6;
auto x7 = -x2 + x3*(x23 + V{1});
auto x8 = -V{1}*cell[1] + cell[3] + V{1}*x4 - x7;
auto x27 = x26*x8;
auto x9 = -x2 - x24*(x25 + V{-1});
auto x10 = -x2 + x3*(x25 + V{1});
auto x11 = -cell[0] - cell[1] - V{2}*cell[2] - cell[3] - x10 + x19 + V{1}*x9 + V{-1};
auto x28 = x11*x26;
cell[0] = x0*(x19 + V{-1});
cell[1] = -V{1}*x27 + x4;
cell[2] = -V{1}*x28 + x9;
cell[3] = x27 + x7;
cell[4] = x10 + x28;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x19;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x20;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x21;
return { x19, ((x17)*(x17)) + ((x18)*(x18)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
