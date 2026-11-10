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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x29 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x30 = parameters.template get<descriptors::OMEGA>();
auto x31 = cell[5] + cell[6];
auto x32 = cell[1] + cell[8];
auto x0 = V{1} - x30;
auto x1 = V{3};
auto x2 = -V{1}*cell[3];
auto x33 = cell[7] + x2;
auto x3 = cell[0] + cell[2] + cell[3] + cell[4] + cell[7] + x31 + x32 + V{1};
auto x34 = x29/x3;
auto x4 = x34*(-V{1}*cell[1] - V{1}*cell[2] + x31 + x33);
auto x35 = ((x4)*(x4));
auto x5 = x34*(-V{1}*cell[4] - V{1}*cell[5] + x32 + x33);
auto x36 = ((x5)*(x5));
auto x6 = x35 + x36;
auto x7 = V{0.5}*x1*x6;
auto x37 = V{1} - x7;
auto x8 = V{0.444444444444444};
auto x9 = -V{1}*x4;
auto x10 = x5 + x9;
auto x11 = V{0.5}*((x1)*(x1));
auto x12 = V{0.0277777777777778};
auto x13 = x12*x3;
auto x14 = V{0.111111111111111};
auto x15 = x14*x3;
auto x16 = -V{1}*x5;
auto x17 = x16 + x9;
auto x18 = x16 + x4;
auto x19 = x4 + x5;
auto x20 = cell[0]*x0 - x30*x8*(x3*(x7 + V{-1}) + V{1});
auto x21 = cell[1]*x0 - x30*(x12 - x13*(x1*x10 + x11*((x10)*(x10)) + x37));
auto x22 = cell[2]*x0 - x30*(x14 - x15*(x1*x9 + x11*((x9)*(x9)) + x37));
auto x23 = cell[3]*x0 - x30*(x12 - x13*(x1*x17 + x11*((x17)*(x17)) + x37));
auto x24 = cell[4]*x0 - x30*(x14 - x15*(x1*x16 + x11*((x16)*(x16)) + x37));
auto x25 = cell[5]*x0 - x30*(x12 - x13*(x1*x18 + x11*((x18)*(x18)) + x37));
auto x26 = cell[6]*x0 - x30*(x14 - x15*(x1*x4 + x11*x35 + x37));
auto x27 = cell[7]*x0 - x30*(x12 - x13*(x1*x19 + x11*((x19)*(x19)) + x37));
auto x28 = cell[8]*x0 - x30*(x14 - x15*(x1*x5 + x11*x36 + x37));
cell[0] = x20;
cell[1] = x21;
cell[2] = x22;
cell[3] = x23;
cell[4] = x24;
cell[5] = x25;
cell[6] = x26;
cell[7] = x27;
cell[8] = x28;
return { x3, x6 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
