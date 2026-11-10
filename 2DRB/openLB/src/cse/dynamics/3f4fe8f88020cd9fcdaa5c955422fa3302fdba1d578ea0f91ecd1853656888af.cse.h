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
struct CSE<dynamics::Tuple<T, descriptors::D2Q5<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::BGK, AdvectionDiffusionExternalVelocityCollision>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x11 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x12 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x13 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{1} - x13;
auto x1 = V{0.333333333333333};
auto x2 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + V{1};
auto x3 = V{3};
auto x14 = x11*x3;
auto x16 = x12*x3;
auto x4 = V{0.166666666666667};
auto x5 = x2*x4;
auto x15 = V{1}*x5;
auto x6 = cell[0]*x0 + x1*x13*(x2 + V{-1});
auto x7 = cell[1]*x0 - x13*(x15*(x14 + V{-1}) + x4);
auto x8 = cell[2]*x0 - x13*(x15*(x16 + V{-1}) + x4);
auto x9 = cell[3]*x0 - x13*(x4 - x5*(x14 + V{1}));
auto x10 = cell[4]*x0 - x13*(x4 - x5*(x16 + V{1}));
cell[0] = x6;
cell[1] = x7;
cell[2] = x8;
cell[3] = x9;
cell[4] = x10;
return { x2, ((x11)*(x11)) + ((x12)*(x12)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
