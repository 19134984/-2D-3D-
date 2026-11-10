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
struct CSE<ZouHeDynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, 1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, 0, 1>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x10 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x11 = parameters.template get<descriptors::OMEGA>();
auto x9 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x12 = x11 + V{-1};
auto x13 = (cell[0] + cell[4] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7] + cell[8] + V{1})/(x9 + V{1});
auto x14 = ((x9)*(x9));
auto x15 = V{1.5}*x14;
auto x16 = ((x10)*(x10));
auto x17 = V{1.5}*x16;
auto x18 = x17 + V{-1};
auto x19 = x15 + x18;
auto x20 = V{0.0138888888888889}*x13;
auto x21 = V{3}*x9;
auto x22 = -x21;
auto x23 = x10 - x9;
auto x24 = V{4.5}*((x23)*(x23));
auto x25 = V{3}*x10;
auto x26 = x19 + x25;
auto x27 = V{4.5}*((x10 + x9)*(x10 + x9));
auto x28 = V{1} - x17;
auto x29 = -x15;
auto x30 = x25 + x29;
auto x31 = x28 + x30;
auto x32 = x21 + x27 + x31;
auto x33 = x20*x32;
auto x34 = x21 + x26 - x27;
auto x35 = -x25;
auto x36 = V{4.5}*((x23)*(x23));
auto x37 = V{0.5}*cell[4] - V{0.5}*cell[8] + V{0.5}*x10*x13;
auto x38 = V{0.0277777777777778}*x11;
auto x39 = x22 + x31 + x36;
auto x40 = V{0.111111111111111}*x13;
auto x41 = V{3}*x14;
auto x42 = x28 + x41;
auto x43 = x21 + x42;
auto x44 = V{0.111111111111111}*x11;
auto x45 = x29 + x35;
auto x46 = x21 + x24 + x28 + x45;
auto x47 = -x34;
auto x48 = V{3}*x16 + V{1};
auto x0 = -cell[0]*x12 - V{0.444444444444444}*x11*(x13*x19 + V{1});
auto x1 = -x12*(cell[5] - x20*x34 + x20*(x22 - x24 + x26) - x20*(x19 + x21 + x35 - x36) - x33 + x37) - x38*(-x13*x39 + V{1});
auto x2 = -x12*(cell[6] - x40*x43 - x40*(x18 + x21 - x41)) - x44*(-x13*(x22 + x42) + V{1});
auto x3 = x12*(-cell[7] - x20*x39 + x20*x46 - x20*x47 + x33 + x37) - x38*(-x13*x47 + V{1});
auto x4 = -cell[4]*x12 - x44*(-x13*(x45 + x48) + V{1});
auto x5 = -cell[5]*x12 - x38*(-x13*x46 + V{1});
auto x6 = -cell[6]*x12 - x44*(-x13*x43 + V{1});
auto x7 = -cell[7]*x12 - x38*(-x13*x32 + V{1});
auto x8 = -cell[8]*x12 - x44*(-x13*(x30 + x48) + V{1});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { V{1}*x13, x14 + x16 };
}
};

}

}

// Generation Info: commit=dedfcf8ee2d0fa6903db6f80111ed262a9148a4a
