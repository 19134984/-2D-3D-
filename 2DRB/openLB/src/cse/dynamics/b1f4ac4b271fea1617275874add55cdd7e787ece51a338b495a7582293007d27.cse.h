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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<1, 1>, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x26 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x27 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<1, 1>::VELOCITY>(0);
auto x29 = parameters.template get<descriptors::OMEGA>();
auto x30 = ((x27)*(x27));
auto x0 = V{1} - x29;
auto x1 = V{3};
auto x2 = V{-1} + (cell[0] + V{2}*cell[1] + cell[2] + cell[6] + V{2}*cell[7] + V{2}*cell[8] + V{1})/x26;
auto x28 = ((x2)*(x2));
auto x3 = x28 + x30;
auto x4 = V{0.5}*x1*x3;
auto x32 = V{1} - x4;
auto x5 = V{0.444444444444444};
auto x6 = -V{1}*x27;
auto x7 = x2 + x6;
auto x8 = V{0.5}*((x1)*(x1));
auto x9 = V{0.0277777777777778};
auto x31 = -x9;
auto x10 = x26*x9;
auto x11 = V{0.111111111111111};
auto x12 = x11*x26;
auto x13 = -V{1}*x2;
auto x14 = x13 + x6;
auto x15 = x13 + x27;
auto x16 = x2 + x27;
auto x17 = cell[0]*x0 - x29*x5*(x26*(x4 + V{-1}) + V{1});
auto x18 = cell[1]*x0 + x29*(x10*(x1*x7 + x32 + x8*((x7)*(x7))) + x31);
auto x19 = cell[2]*x0 - x29*(x11 - x12*(x1*x6 + x32 + x8*((x6)*(x6))));
auto x20 = cell[3]*x0 + x29*(x10*(x1*x14 + x32 + x8*((x14)*(x14))) + x31);
auto x21 = cell[4]*x0 - x29*(x11 - x12*(x1*x13 + x32 + x8*((x13)*(x13))));
auto x22 = cell[5]*x0 + x29*(x10*(x1*x15 + x32 + x8*((x15)*(x15))) + x31);
auto x23 = cell[6]*x0 - x29*(x11 - x12*(x1*x27 + x30*x8 + x32));
auto x24 = cell[7]*x0 + x29*(x10*(x1*x16 + x32 + x8*((x16)*(x16))) + x31);
auto x25 = cell[8]*x0 - x29*(x11 - x12*(x1*x2 + x28*x8 + x32));
cell[0] = x17;
cell[1] = x18;
cell[2] = x19;
cell[3] = x20;
cell[4] = x21;
cell[5] = x22;
cell[6] = x23;
cell[7] = x24;
cell[8] = x25;
return { x26, x3 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
