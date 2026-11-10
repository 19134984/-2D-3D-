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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::FirstOrder, collision::BGK, forcing::ConservativePhaseField>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x24 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x25 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x26 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x27 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x28 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{1} - x28;
auto x1 = V{0.444444444444444};
auto x2 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + V{1};
auto x3 = V{3};
auto x4 = -V{1}*x26;
auto x5 = V{0.0277777777777778};
auto x6 = x2*x5;
auto x7 = V{1} - V{0.5}*x28;
auto x8 = x3*x5;
auto x29 = x7*x8;
auto x9 = -V{1}*x24;
auto x10 = V{0.111111111111111};
auto x11 = x10*x2;
auto x12 = x10*x3;
auto x30 = x12*x7;
auto x13 = -V{1}*x27;
auto x14 = -V{1}*x25;
auto x15 = cell[0]*x0 + x1*x28*(x2 + V{-1});
auto x16 = cell[1]*x0 - x28*(x5 - x6*(x3*(x27 + x4) + V{1})) + x29*(x25 + x9);
auto x17 = cell[2]*x0 - x28*(x10 - x11*(x3*x4 + V{1})) + x30*x9;
auto x18 = cell[3]*x0 - x28*(x5 - x6*(x3*(x13 + x4) + V{1})) + x29*(x14 + x9);
auto x19 = cell[4]*x0 + x14*x30 - x28*(x10 - x11*(x13*x3 + V{1}));
auto x20 = cell[5]*x0 - x28*(x5 - x6*(x3*(x13 + x26) + V{1})) + x29*(x14 + x24);
auto x21 = cell[6]*x0 + x24*x30 - x28*(x10 - x11*(x26*x3 + V{1}));
auto x22 = cell[7]*x0 - x28*(x5 - x6*(x3*(x26 + x27) + V{1})) + x29*(x24 + x25);
auto x23 = cell[8]*x0 + x25*x30 - x28*(x10 - x11*(x27*x3 + V{1}));
cell[0] = x15;
cell[1] = x16;
cell[2] = x17;
cell[3] = x18;
cell[4] = x19;
cell[5] = x20;
cell[6] = x21;
cell[7] = x22;
cell[8] = x23;
return { x2, ((x26)*(x26)) + ((x27)*(x27)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
