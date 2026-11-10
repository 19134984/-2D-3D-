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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::OmegaFromCellTauEff<collision::BGK>, forcing::Guo<momenta::Forced> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x44 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x45 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x46 = cell.template getFieldComponent<olb::descriptors::TAU_EFF>(0);
auto x47 = cell[5] + cell[6];
auto x48 = cell[1] + cell[8];
auto x0 = V{1}/x46;
auto x1 = V{1} - x0;
auto x2 = V{3};
auto x3 = -V{1}*cell[3];
auto x50 = cell[7] + x3;
auto x4 = cell[0] + cell[2] + cell[3] + cell[4] + cell[7] + x47 + x48 + V{1};
auto x49 = V{1} / (x4);
auto x5 = V{0.5}*x44 + x49*(-V{1}*cell[1] - V{1}*cell[2] + x47 + x50);
auto x51 = ((x5)*(x5));
auto x53 = x2*x5;
auto x6 = V{0.5}*x45 + x49*(-V{1}*cell[4] - V{1}*cell[5] + x48 + x50);
auto x52 = ((x6)*(x6));
auto x54 = x2*x6;
auto x7 = x51 + x52;
auto x8 = V{0.5}*x2*x7;
auto x57 = V{1} - x8;
auto x9 = V{0.444444444444444};
auto x10 = -x44*x53;
auto x11 = -x45*x54;
auto x12 = -x4*(V{0.5}*x0 + V{-1});
auto x13 = -V{1}*x5;
auto x14 = x13 + x6;
auto x15 = ((x2)*(x2));
auto x55 = V{1}*x15;
auto x16 = V{0.5}*x15;
auto x17 = V{0.0277777777777778};
auto x18 = x17*x4;
auto x19 = x14*x15;
auto x20 = -x2*(x5 + V{1});
auto x56 = -x20;
auto x21 = -x2*(x6 + V{-1});
auto x22 = x12*x17;
auto x23 = V{0.111111111111111};
auto x24 = x23*x4;
auto x25 = x12*x23;
auto x26 = -V{1}*x6;
auto x27 = x13 + x26;
auto x28 = -x27*x55;
auto x29 = -x2*(x6 + V{1});
auto x30 = x26 + x5;
auto x31 = x15*x30;
auto x32 = -x2*(x5 + V{-1});
auto x33 = x5 + x6;
auto x34 = x15*x33;
auto x35 = cell[0]*x1 - x0*x9*(x4*(x8 + V{-1}) + V{1}) + x12*x9*(x10 + x11);
auto x36 = -(-cell[1]*x1 + x0*(x17 - x18*(x14*x2 + x16*((x14)*(x14)) + x57)) + x22*(x44*(V{1}*x19 + x56) - x45*(x19 + x21)));
auto x37 = cell[2]*x1 - x0*(x23 - x24*(x13*x2 + x16*((x13)*(x13)) + x57)) + x25*(x11 - x44*(x13*x55 + x56));
auto x38 = cell[3]*x1 - x0*(x17 - x18*(x16*((x27)*(x27)) + x2*x27 + x57)) + x22*(x44*(x20 + x28) + x45*(x28 + x29));
auto x39 = cell[4]*x1 - x0*(x23 - x24*(x16*((x26)*(x26)) + x2*x26 + x57)) + x25*(x10 - x45*(x26*x55 - x29));
auto x40 = cell[5]*x1 - x0*(x17 - x18*(x16*((x30)*(x30)) + x2*x30 + x57)) + x22*(x44*(x31 + x32) + x45*(x29 - V{1}*x31));
auto x41 = cell[6]*x1 - x0*(x23 - x24*(x16*x51 + x53 + x57)) + x25*(x11 + x44*(x15*x5 + x32));
auto x42 = cell[7]*x1 - x0*(x17 - x18*(x16*((x33)*(x33)) + x2*x33 + x57)) + x22*(x44*(x32 + x34) + x45*(x21 + x34));
auto x43 = cell[8]*x1 - x0*(x23 - x24*(x16*x52 + x54 + x57)) + x25*(x10 + x45*(x15*x6 + x21));
cell[0] = x35;
cell[1] = x36;
cell[2] = x37;
cell[3] = x38;
cell[4] = x39;
cell[5] = x40;
cell[6] = x41;
cell[7] = x42;
cell[8] = x43;
return { x4, x7 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
