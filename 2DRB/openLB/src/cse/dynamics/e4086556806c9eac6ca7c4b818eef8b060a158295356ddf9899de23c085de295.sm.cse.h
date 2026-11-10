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
struct CSE<dynamics::Tuple<T, descriptors::D2Q5<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::FirstOrder, collision::OmegaFromCellTauEff<collision::BGK>, AdvectionDiffusionExternalVelocityCollision>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x18 = cell.template getFieldComponent<olb::descriptors::TAU_EFF>(0);
auto x19 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x20 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x0 = V{1}/x18;
auto x1 = V{1} - x0;
auto x2 = V{0.333333333333333};
auto x3 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + V{1};
auto x4 = cell[0]*x1 + x0*x2*(x3 + V{-1});
auto x5 = V{3};
auto x21 = x19*x5;
auto x23 = x20*x5;
auto x6 = V{0.166666666666667};
auto x7 = x3*x6;
auto x22 = V{1}*x7;
auto x8 = cell[1]*x1 - x0*(x22*(x21 + V{-1}) + x6);
auto x9 = cell[2]*x1 - x0*(x22*(x23 + V{-1}) + x6);
auto x10 = cell[3]*x1 - x0*(x6 - x7*(x21 + V{1}));
auto x11 = cell[4]*x1 - x0*(x6 - x7*(x23 + V{1}));
auto x12 = x10 + x11 + x4 + x8 + x9 + V{1};
auto x24 = V{1} / (x12);
cell[0] = x4;
cell[1] = x8;
cell[2] = x9;
cell[3] = x10;
cell[4] = x11;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x12;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x24*(x10 - V{1}*x8);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x24*(x11 - V{1}*x9);
return { x3, ((x19)*(x19)) + ((x20)*(x20)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
