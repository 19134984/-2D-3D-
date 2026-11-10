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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, -1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::ConstRhoBGK, dynamics::DefaultCombination>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x37 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x38 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x39 = parameters.template get<descriptors::OMEGA>();
auto x40 = parameters.template get<statistics::AVERAGE_RHO>();
auto x41 = ((x37)*(x37));
auto x42 = ((x38)*(x38));
auto x0 = V{3};
auto x1 = x41 + x42;
auto x2 = V{0.5}*x0*x1;
auto x44 = V{1} - x2;
auto x3 = V{0.444444444444444};
auto x4 = -V{1}*x37;
auto x5 = (cell[0] + V{2}*cell[1] + V{2}*cell[2] + V{2}*cell[3] + cell[4] + cell[8] + V{1})/(x4 + V{1});
auto x6 = -x3*(x5*(x2 + V{-1}) + V{1});
auto x7 = V{1} - x39;
auto x8 = V{1} - x40;
auto x9 = V{1} + x8/x5;
auto x10 = x38 + x4;
auto x11 = V{0.5}*((x0)*(x0));
auto x12 = V{0.0277777777777778};
auto x43 = -x12;
auto x13 = x12*x5;
auto x14 = x13*(x0*x10 + x11*((x10)*(x10)) + x44) + x43;
auto x15 = V{0.111111111111111};
auto x45 = -x15;
auto x16 = x15*x5;
auto x17 = x16*(x0*x4 + x11*((x4)*(x4)) + x44) + x45;
auto x18 = -V{1}*x38;
auto x19 = x18 + x4;
auto x20 = x13*(x0*x19 + x11*((x19)*(x19)) + x44) + x43;
auto x21 = x16*(x0*x18 + x11*((x18)*(x18)) + x44) + x45;
auto x22 = x18 + x37;
auto x23 = x13*(x0*x22 + x11*((x22)*(x22)) + x44) + x43;
auto x24 = x16*(x0*x37 + x11*x41 + x44) + x45;
auto x25 = x37 + x38;
auto x26 = x13*(x0*x25 + x11*((x25)*(x25)) + x44) + x43;
auto x27 = x16*(x0*x38 + x11*x42 + x44) + x45;
auto x28 = -x3 + x7*(cell[0] - x6) + x9*(x3 + x6);
auto x29 = x43 + x7*(cell[1] - x14) + x9*(x12 + x14);
auto x30 = x45 + x7*(cell[2] - x17) + x9*(x15 + x17);
auto x31 = x43 + x7*(cell[3] - x20) + x9*(x12 + x20);
auto x32 = x45 + x7*(cell[4] - x21) + x9*(x15 + x21);
auto x33 = x43 + x7*(cell[5] - x23) + x9*(x12 + x23);
auto x34 = x45 + x7*(cell[6] - x24) + x9*(x15 + x24);
auto x35 = x43 + x7*(cell[7] - x26) + x9*(x12 + x26);
auto x36 = x45 + x7*(cell[8] - x27) + x9*(x15 + x27);
cell[0] = x28;
cell[1] = x29;
cell[2] = x30;
cell[3] = x31;
cell[4] = x32;
cell[5] = x33;
cell[6] = x34;
cell[7] = x35;
cell[8] = x36;
return { x5 + x8, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
