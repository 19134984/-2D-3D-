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
struct CSE<ZouHeDynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, -1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, 0, -1>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x11 = parameters.template get<descriptors::OMEGA>();
auto x10 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x9 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x12 = x11 + V{-1};
auto x13 = x9 + V{-1};
auto x14 = cell[0] + V{2}*cell[1] + V{2}*cell[2] + V{2}*cell[3] + cell[4] + cell[8] + V{1};
auto x15 = -x14/x13;
auto x16 = ((x9)*(x9));
auto x17 = V{1.5}*x16;
auto x18 = ((x10)*(x10));
auto x19 = V{1.5}*x18;
auto x20 = x19 + V{-1};
auto x21 = x17 + x20;
auto x22 = V{0.0277777777777778}*x11;
auto x23 = V{3}*x9;
auto x24 = -x23;
auto x25 = x10 - x9;
auto x26 = V{4.5}*((x25)*(x25));
auto x27 = V{1} - x19;
auto x28 = -x17;
auto x29 = V{3}*x10;
auto x30 = x28 + x29;
auto x31 = x27 + x30;
auto x32 = V{0.111111111111111}*x11;
auto x33 = V{3}*x16;
auto x34 = x27 + x33;
auto x35 = ((x10 + x9)*(x10 + x9));
auto x36 = x21 + x29;
auto x37 = x23 - V{4.5}*x35 + x36;
auto x38 = -x29;
auto x39 = x28 + x38;
auto x40 = V{3}*x18 + V{1};
auto x41 = V{1} / (x13);
auto x42 = x14*x41;
auto x43 = V{0.0138888888888889}*x42;
auto x44 = V{4.5}*((x25)*(x25));
auto x45 = x24 + x36 - x44;
auto x46 = x23 + x31 + V{4.5}*x35;
auto x47 = x43*x46;
auto x48 = x37*x43;
auto x49 = x43*(x21 + x23 - x26 + x38);
auto x50 = -V{0.5}*cell[4] + V{0.5}*cell[8] + V{0.5}*x10*x42;
auto x51 = V{0.111111111111111}*x42;
auto x52 = x23 + x34;
auto x0 = -cell[0]*x12 - V{0.444444444444444}*x11*(x15*x21 + V{1});
auto x1 = -cell[1]*x12 - x22*(-x15*(x24 + x26 + x31) + V{1});
auto x2 = -cell[2]*x12 - x32*(-x15*(x24 + x34) + V{1});
auto x3 = -cell[3]*x12 - x22*(x15*x37 + V{1});
auto x4 = -cell[4]*x12 - x32*(-x15*(x39 + x40) + V{1});
auto x5 = -x12*(cell[1] + x43*x45 - x47 - x48 - x49 + x50) - x22*(-x15*(x23 + x27 + x39 + x44) + V{1});
auto x6 = -x12*(cell[2] - x51*x52 - x51*(x20 + x23 - x33)) - x32*(-x15*x52 + V{1});
auto x7 = -x12*(cell[3] + V{0.0138888888888889}*x14*x41*x45 - x47 - x48 - x49 - x50) - x22*(-x15*x46 + V{1});
auto x8 = -cell[8]*x12 - x32*(-x15*(x30 + x40) + V{1});
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
return { -V{1}*x42, x16 + x18 };
}
};

}

}

// Generation Info: commit=dedfcf8ee2d0fa6903db6f80111ed262a9148a4a
