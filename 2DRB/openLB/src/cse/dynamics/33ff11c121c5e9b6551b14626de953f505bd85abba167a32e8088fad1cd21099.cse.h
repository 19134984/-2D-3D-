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
struct CSE<AdvectionDiffusionEdgesDynamics<T, descriptors::D3Q7<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q7<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, 2, -1, -1>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x24 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x25 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x26 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x27 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x31 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{0.25};
auto x1 = V{4};
auto x28 = x1*x24;
auto x30 = x1*x25;
auto x32 = x1*x26;
auto x2 = V{0.125};
auto x3 = x2*x27;
auto x29 = V{1}*x3;
auto x4 = -x2 - x29*(x28 + V{-1});
auto x5 = V{1} - x31;
auto x6 = x1*x2;
auto x33 = x5*x6;
auto x7 = cell[1] - x4;
auto x8 = -x2 + x3*(x28 + V{1});
auto x9 = -V{2}*x7;
auto x34 = x33*x9;
auto x10 = -x2 - x29*(x30 + V{-1});
auto x11 = cell[2] - x10;
auto x12 = -x2 + x3*(x30 + V{1});
auto x13 = -V{2}*x11;
auto x35 = x13*x33;
auto x14 = -x2 - x29*(x32 + V{-1});
auto x15 = -x2 + x3*(x32 + V{1});
auto x16 = -V{1}*cell[3] + cell[6] + V{1}*x14 - x15;
auto x36 = x16*x33;
cell[0] = x0*(x27 + V{-1});
cell[1] = -V{1}*x34 + x4;
cell[2] = x10 - V{1}*x35;
cell[3] = x14 - V{1}*x36;
cell[4] = x34 + x8;
cell[5] = x12 + x35;
cell[6] = x15 + x36;
return { x27, ((x24)*(x24)) + ((x25)*(x25)) + ((x26)*(x26)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
