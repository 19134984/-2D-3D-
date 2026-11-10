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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::ConstRhoBGK, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x51 = parameters.template get<descriptors::OMEGA>();
auto x52 = parameters.template get<statistics::AVERAGE_RHO>();
auto x53 = cell[5] + cell[6];
auto x54 = cell[1] + cell[8];
auto x0 = V{3};
auto x1 = -V{1}*cell[3];
auto x56 = cell[7] + x1;
auto x2 = cell[0] + cell[2] + cell[3] + cell[4] + cell[7] + x53 + x54 + V{1};
auto x55 = V{1} / (x2);
auto x3 = x55*(-V{1}*cell[1] - V{1}*cell[2] + x53 + x56);
auto x57 = ((x3)*(x3));
auto x4 = x55*(-V{1}*cell[4] - V{1}*cell[5] + x54 + x56);
auto x58 = ((x4)*(x4));
auto x5 = x57 + x58;
auto x6 = V{0.5}*x0*x5;
auto x60 = V{1} - x6;
auto x7 = V{0.444444444444444};
auto x8 = -x7*(x2*(x6 + V{-1}) + V{1});
auto x9 = V{1} - x51;
auto x10 = V{1} - x52;
auto x11 = x10*x55 + V{1};
auto x12 = x11*(x7 + x8) - x7 + x9*(cell[0] - x8);
auto x13 = -V{1}*x3;
auto x14 = x13 + x4;
auto x15 = V{0.5}*((x0)*(x0));
auto x16 = V{0.0277777777777778};
auto x59 = -x16;
auto x17 = x16*x2;
auto x18 = x17*(x0*x14 + x15*((x14)*(x14)) + x60) + x59;
auto x19 = x11*(x16 + x18) + x59 + x9*(cell[1] - x18);
auto x20 = V{0.111111111111111};
auto x61 = -x20;
auto x21 = x2*x20;
auto x22 = x21*(x0*x13 + x15*((x13)*(x13)) + x60) + x61;
auto x23 = x11*(x20 + x22) + x61 + x9*(cell[2] - x22);
auto x24 = -V{1}*x4;
auto x25 = x13 + x24;
auto x26 = x17*(x0*x25 + x15*((x25)*(x25)) + x60) + x59;
auto x27 = x11*(x16 + x26) + x59 + x9*(cell[3] - x26);
auto x28 = x21*(x0*x24 + x15*((x24)*(x24)) + x60) + x61;
auto x29 = x11*(x20 + x28) + x61 + x9*(cell[4] - x28);
auto x30 = x24 + x3;
auto x31 = x17*(x0*x30 + x15*((x30)*(x30)) + x60) + x59;
auto x32 = x11*(x16 + x31) + x59 + x9*(cell[5] - x31);
auto x33 = x21*(x0*x3 + x15*x57 + x60) + x61;
auto x34 = x11*(x20 + x33) + x61 + x9*(cell[6] - x33);
auto x62 = x32 + x34;
auto x35 = x3 + x4;
auto x36 = x17*(x0*x35 + x15*((x35)*(x35)) + x60) + x59;
auto x37 = x11*(x16 + x36) + x59 + x9*(cell[7] - x36);
auto x38 = x21*(x0*x4 + x15*x58 + x60) + x61;
auto x39 = x11*(x20 + x38) + x61 + x9*(cell[8] - x38);
auto x63 = x19 + x39;
auto x40 = x12 + x23 + x27 + x29 + x37 + x62 + x63 + V{1};
auto x64 = V{1} / (x40);
auto x41 = -V{1}*x27;
auto x65 = x37 + x41;
cell[0] = x12;
cell[1] = x19;
cell[2] = x23;
cell[3] = x27;
cell[4] = x29;
cell[5] = x32;
cell[6] = x34;
cell[7] = x37;
cell[8] = x39;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x40;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x64*(-V{1}*x19 - V{1}*x23 + x62 + x65);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x64*(-V{1}*x29 - V{1}*x32 + x63 + x65);
return { x10 + x2, x5 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
