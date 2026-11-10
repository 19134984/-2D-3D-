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
struct CSE<AdvectionDiffusionBoundariesDynamics<T, descriptors::D3Q7<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q7<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, 0, 1>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x22 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x23 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x24 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x25 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x29 = parameters.template get<descriptors::OMEGA>();
auto x28 = V{1}*cell[2];
auto x30 = V{1}*cell[3];
auto x0 = V{0.25};
auto x1 = V{4};
auto x26 = x1*x22;
auto x31 = x1*x23;
auto x32 = x1*x24;
auto x2 = V{0.125};
auto x3 = x2*x25;
auto x27 = V{1}*x3;
auto x4 = -x2 - x27*(x26 + V{-1});
auto x5 = V{1} - x29;
auto x6 = x1*x2;
auto x33 = x5*x6;
auto x7 = -x2 + x3*(x26 + V{1});
auto x8 = V{1}*cell[0] + V{2}*cell[4] + V{1}*cell[5] + V{1}*cell[6] - V{1}*x25 + x28 + x30 + V{1}*x4 - x7 + V{1};
auto x34 = x33*x8;
auto x9 = -x2 - x27*(x31 + V{-1});
auto x10 = -x2 + x3*(x31 + V{1});
auto x11 = cell[5] - x10 - x28 + V{1}*x9;
auto x35 = x11*x33;
auto x12 = -x2 - x27*(x32 + V{-1});
auto x13 = -x2 + x3*(x32 + V{1});
auto x14 = cell[6] + V{1}*x12 - x13 - x30;
auto x36 = x14*x33;
cell[0] = x0*(x25 + V{-1});
cell[1] = -V{1}*x34 + x4;
cell[2] = -V{1}*x35 + x9;
cell[3] = x12 - V{1}*x36;
cell[4] = x34 + x7;
cell[5] = x10 + x35;
cell[6] = x13 + x36;
return { x25, ((x22)*(x22)) + ((x23)*(x23)) + ((x24)*(x24)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
