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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::ConstRhoBGK, dynamics::DefaultCombination>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x40 = parameters.template get<descriptors::OMEGA>();
auto x41 = parameters.template get<statistics::AVERAGE_RHO>();
auto x42 = cell[5] + cell[6];
auto x43 = cell[1] + cell[8];
auto x0 = V{3};
auto x1 = -V{1}*cell[3];
auto x45 = cell[7] + x1;
auto x2 = cell[0] + cell[2] + cell[3] + cell[4] + cell[7] + x42 + x43 + V{1};
auto x44 = V{1} / (x2);
auto x3 = x44*(-V{1}*cell[1] - V{1}*cell[2] + x42 + x45);
auto x46 = ((x3)*(x3));
auto x4 = x44*(-V{1}*cell[4] - V{1}*cell[5] + x43 + x45);
auto x47 = ((x4)*(x4));
auto x5 = x46 + x47;
auto x6 = V{0.5}*x0*x5;
auto x49 = V{1} - x6;
auto x7 = V{0.444444444444444};
auto x8 = -x7*(x2*(x6 + V{-1}) + V{1});
auto x9 = V{1} - x40;
auto x10 = V{1} - x41;
auto x11 = x10*x44 + V{1};
auto x12 = -V{1}*x3;
auto x13 = x12 + x4;
auto x14 = V{0.5}*((x0)*(x0));
auto x15 = V{0.0277777777777778};
auto x48 = -x15;
auto x16 = x15*x2;
auto x17 = x16*(x0*x13 + x14*((x13)*(x13)) + x49) + x48;
auto x18 = V{0.111111111111111};
auto x50 = -x18;
auto x19 = x18*x2;
auto x20 = x19*(x0*x12 + x14*((x12)*(x12)) + x49) + x50;
auto x21 = -V{1}*x4;
auto x22 = x12 + x21;
auto x23 = x16*(x0*x22 + x14*((x22)*(x22)) + x49) + x48;
auto x24 = x19*(x0*x21 + x14*((x21)*(x21)) + x49) + x50;
auto x25 = x21 + x3;
auto x26 = x16*(x0*x25 + x14*((x25)*(x25)) + x49) + x48;
auto x27 = x19*(x0*x3 + x14*x46 + x49) + x50;
auto x28 = x3 + x4;
auto x29 = x16*(x0*x28 + x14*((x28)*(x28)) + x49) + x48;
auto x30 = x19*(x0*x4 + x14*x47 + x49) + x50;
auto x31 = x11*(x7 + x8) - x7 + x9*(cell[0] - x8);
auto x32 = x11*(x15 + x17) + x48 + x9*(cell[1] - x17);
auto x33 = x11*(x18 + x20) + x50 + x9*(cell[2] - x20);
auto x34 = x11*(x15 + x23) + x48 + x9*(cell[3] - x23);
auto x35 = x11*(x18 + x24) + x50 + x9*(cell[4] - x24);
auto x36 = x11*(x15 + x26) + x48 + x9*(cell[5] - x26);
auto x37 = x11*(x18 + x27) + x50 + x9*(cell[6] - x27);
auto x38 = x11*(x15 + x29) + x48 + x9*(cell[7] - x29);
auto x39 = x11*(x18 + x30) + x50 + x9*(cell[8] - x30);
cell[0] = x31;
cell[1] = x32;
cell[2] = x33;
cell[3] = x34;
cell[4] = x35;
cell[5] = x36;
cell[6] = x37;
cell[7] = x38;
cell[8] = x39;
return { x10 + x2, x5 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
