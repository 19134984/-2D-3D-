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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::PlainGuo>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x43 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x44 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x45 = parameters.template get<descriptors::OMEGA>();
auto x46 = cell[5] + cell[6];
auto x47 = cell[1] + cell[8];
auto x0 = V{1} - x45;
auto x1 = V{3};
auto x2 = -V{1}*cell[3];
auto x49 = cell[7] + x2;
auto x3 = cell[0] + cell[2] + cell[3] + cell[4] + cell[7] + x46 + x47 + V{1};
auto x48 = V{1} / (x3);
auto x4 = V{0.5}*x43 + x48*(-V{1}*cell[1] - V{1}*cell[2] + x46 + x49);
auto x50 = ((x4)*(x4));
auto x52 = x1*x4;
auto x5 = V{0.5}*x44 + x48*(-V{1}*cell[4] - V{1}*cell[5] + x47 + x49);
auto x51 = ((x5)*(x5));
auto x53 = x1*x5;
auto x6 = x50 + x51;
auto x7 = V{0.5}*x1*x6;
auto x56 = V{1} - x7;
auto x8 = V{0.444444444444444};
auto x9 = -x43*x52;
auto x10 = -x44*x53;
auto x11 = -x3*(V{0.5}*x45 + V{-1});
auto x12 = -V{1}*x4;
auto x13 = x12 + x5;
auto x14 = ((x1)*(x1));
auto x54 = V{1}*x14;
auto x15 = V{0.5}*x14;
auto x16 = V{0.0277777777777778};
auto x17 = x16*x3;
auto x18 = x13*x14;
auto x19 = -x1*(x4 + V{1});
auto x55 = -x19;
auto x20 = -x1*(x5 + V{-1});
auto x21 = x11*x16;
auto x22 = V{0.111111111111111};
auto x23 = x22*x3;
auto x24 = x11*x22;
auto x25 = -V{1}*x5;
auto x26 = x12 + x25;
auto x27 = -x26*x54;
auto x28 = -x1*(x5 + V{1});
auto x29 = x25 + x4;
auto x30 = x14*x29;
auto x31 = -x1*(x4 + V{-1});
auto x32 = x4 + x5;
auto x33 = x14*x32;
auto x34 = cell[0]*x0 + x11*x8*(x10 + x9) - x45*x8*(x3*(x7 + V{-1}) + V{1});
auto x35 = -(-cell[1]*x0 + x21*(x43*(V{1}*x18 + x55) - x44*(x18 + x20)) + x45*(x16 - x17*(x1*x13 + x15*((x13)*(x13)) + x56)));
auto x36 = cell[2]*x0 + x24*(x10 - x43*(x12*x54 + x55)) - x45*(x22 - x23*(x1*x12 + x15*((x12)*(x12)) + x56));
auto x37 = cell[3]*x0 + x21*(x43*(x19 + x27) + x44*(x27 + x28)) - x45*(x16 - x17*(x1*x26 + x15*((x26)*(x26)) + x56));
auto x38 = -(-cell[4]*x0 + x24*(x44*(x25*x54 - x28) - x9) + x45*(x22 - x23*(x1*x25 + x15*((x25)*(x25)) + x56)));
auto x39 = cell[5]*x0 + x21*(x43*(x30 + x31) + x44*(x28 - V{1}*x30)) - x45*(x16 - x17*(x1*x29 + x15*((x29)*(x29)) + x56));
auto x40 = cell[6]*x0 + x24*(x10 + x43*(x14*x4 + x31)) - x45*(x22 - x23*(x15*x50 + x52 + x56));
auto x41 = cell[7]*x0 + x21*(x43*(x31 + x33) + x44*(x20 + x33)) - x45*(x16 - x17*(x1*x32 + x15*((x32)*(x32)) + x56));
auto x42 = cell[8]*x0 + x24*(x44*(x14*x5 + x20) + x9) - x45*(x22 - x23*(x15*x51 + x53 + x56));
cell[0] = x34;
cell[1] = x35;
cell[2] = x36;
cell[3] = x37;
cell[4] = x38;
cell[5] = x39;
cell[6] = x40;
cell[7] = x41;
cell[8] = x42;
return { x3, x6 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
