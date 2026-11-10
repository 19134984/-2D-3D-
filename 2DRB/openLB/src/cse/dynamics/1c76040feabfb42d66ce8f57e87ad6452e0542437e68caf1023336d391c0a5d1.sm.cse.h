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
struct CSE<dynamics::Tuple<T, descriptors::D3Q7<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::BGK, AdvectionDiffusionExternalVelocityCollision>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x21 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x22 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x23 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x24 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{1} - x24;
auto x1 = V{0.25};
auto x2 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + V{1};
auto x3 = cell[0]*x0 + x1*x24*(x2 + V{-1});
auto x4 = V{4};
auto x25 = x21*x4;
auto x27 = x22*x4;
auto x28 = x23*x4;
auto x5 = V{0.125};
auto x6 = x2*x5;
auto x26 = V{1}*x6;
auto x7 = cell[1]*x0 - x24*(x26*(x25 + V{-1}) + x5);
auto x8 = cell[2]*x0 - x24*(x26*(x27 + V{-1}) + x5);
auto x9 = cell[3]*x0 - x24*(x26*(x28 + V{-1}) + x5);
auto x10 = cell[4]*x0 - x24*(x5 - x6*(x25 + V{1}));
auto x11 = cell[5]*x0 - x24*(x5 - x6*(x27 + V{1}));
auto x12 = cell[6]*x0 - x24*(x5 - x6*(x28 + V{1}));
auto x13 = x10 + x11 + x12 + x3 + x7 + x8 + x9 + V{1};
auto x29 = V{1} / (x13);
cell[0] = x3;
cell[1] = x7;
cell[2] = x8;
cell[3] = x9;
cell[4] = x10;
cell[5] = x11;
cell[6] = x12;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x13;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x29*(x10 - V{1}*x7);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x29*(x11 - V{1}*x8);
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x29*(x12 - V{1}*x9);
return { x2, ((x21)*(x21)) + ((x22)*(x22)) + ((x23)*(x23)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
