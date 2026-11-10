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
struct CSE<dynamics::Tuple<T, descriptors::D3Q7<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::FirstOrder, collision::BGK, forcing::ConservativePhaseField>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x15 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x16 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x17 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x18 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x19 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x20 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x21 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{1} - x21;
auto x1 = V{0.25};
auto x2 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + V{1};
auto x3 = V{4};
auto x24 = x18*x3;
auto x27 = x19*x3;
auto x29 = x20*x3;
auto x4 = V{0.125};
auto x5 = x2*x4;
auto x25 = V{1}*x5;
auto x6 = V{1} - V{0.5}*x21;
auto x7 = x3*x4;
auto x22 = x6*x7;
auto x23 = x15*x22;
auto x26 = x16*x22;
auto x28 = x17*x22;
auto x8 = cell[0]*x0 + x1*x21*(x2 + V{-1});
auto x9 = cell[1]*x0 - x21*(-x25*(V{1} - x24) + x4) - V{1}*x23;
auto x10 = cell[2]*x0 - x21*(-x25*(V{1} - x27) + x4) - V{1}*x26;
auto x11 = cell[3]*x0 - x21*(-x25*(V{1} - x29) + x4) - V{1}*x28;
auto x12 = cell[4]*x0 - x21*(x4 - x5*(x24 + V{1})) + x23;
auto x13 = cell[5]*x0 - x21*(x4 - x5*(x27 + V{1})) + x26;
auto x14 = cell[6]*x0 - x21*(x4 - x5*(x29 + V{1})) + x28;
cell[0] = x8;
cell[1] = x9;
cell[2] = x10;
cell[3] = x11;
cell[4] = x12;
cell[5] = x13;
cell[6] = x14;
return { x2, ((x18)*(x18)) + ((x19)*(x19)) + ((x20)*(x20)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
