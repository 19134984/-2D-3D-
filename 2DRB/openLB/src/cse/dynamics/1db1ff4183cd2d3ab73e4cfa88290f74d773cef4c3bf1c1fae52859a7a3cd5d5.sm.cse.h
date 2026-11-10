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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::PlainGuo>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x56 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x57 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x58 = parameters.template get<descriptors::OMEGA>();
auto x59 = cell[5] + cell[6];
auto x60 = cell[1] + cell[8];
auto x0 = V{1} - x58;
auto x1 = V{3};
auto x2 = -V{1}*cell[3];
auto x62 = cell[7] + x2;
auto x3 = cell[0] + cell[2] + cell[3] + cell[4] + cell[7] + x59 + x60 + V{1};
auto x61 = V{1} / (x3);
auto x4 = V{0.5}*x56;
auto x5 = x4 + x61*(-V{1}*cell[1] - V{1}*cell[2] + x59 + x62);
auto x63 = ((x5)*(x5));
auto x65 = x1*x5;
auto x6 = V{0.5}*x57;
auto x7 = x6 + x61*(-V{1}*cell[4] - V{1}*cell[5] + x60 + x62);
auto x64 = ((x7)*(x7));
auto x66 = x1*x7;
auto x8 = x63 + x64;
auto x9 = V{0.5}*x1*x8;
auto x68 = V{1} - x9;
auto x10 = V{0.444444444444444};
auto x11 = -x56*x65;
auto x12 = -x57*x66;
auto x13 = -x3*(V{0.5}*x58 + V{-1});
auto x14 = cell[0]*x0 + x10*x13*(x11 + x12) - x10*x58*(x3*(x9 + V{-1}) + V{1});
auto x15 = -V{1}*x5;
auto x16 = x15 + x7;
auto x17 = ((x1)*(x1));
auto x69 = V{1}*x17;
auto x18 = V{0.5}*x17;
auto x19 = V{0.0277777777777778};
auto x20 = x19*x3;
auto x21 = x16*x17;
auto x22 = -x1*(x5 + V{1});
auto x67 = -x22;
auto x23 = -x1*(x7 + V{-1});
auto x24 = x13*x19;
auto x25 = -(-cell[1]*x0 + x24*(x56*(V{1}*x21 + x67) - x57*(x21 + x23)) + x58*(x19 - x20*(x1*x16 + x18*((x16)*(x16)) + x68)));
auto x26 = V{0.111111111111111};
auto x27 = x26*x3;
auto x28 = x13*x26;
auto x29 = cell[2]*x0 + x28*(x12 - x56*(x15*x69 + x67)) - x58*(x26 - x27*(x1*x15 + x18*((x15)*(x15)) + x68));
auto x30 = -V{1}*x7;
auto x31 = x15 + x30;
auto x32 = -x31*x69;
auto x33 = -x1*(x7 + V{1});
auto x34 = cell[3]*x0 + x24*(x56*(x22 + x32) + x57*(x32 + x33)) - x58*(x19 - x20*(x1*x31 + x18*((x31)*(x31)) + x68));
auto x35 = cell[4]*x0 + x28*(x11 - x57*(x30*x69 - x33)) - x58*(x26 - x27*(x1*x30 + x18*((x30)*(x30)) + x68));
auto x36 = x30 + x5;
auto x37 = x17*x36;
auto x38 = -x1*(x5 + V{-1});
auto x39 = cell[5]*x0 + x24*(x56*(x37 + x38) + x57*(x33 - V{1}*x37)) - x58*(x19 - x20*(x1*x36 + x18*((x36)*(x36)) + x68));
auto x40 = cell[6]*x0 + x28*(x12 + x56*(x17*x5 + x38)) - x58*(x26 - x27*(x18*x63 + x65 + x68));
auto x70 = x39 + x40;
auto x41 = x5 + x7;
auto x42 = x17*x41;
auto x43 = cell[7]*x0 + x24*(x56*(x38 + x42) + x57*(x23 + x42)) - x58*(x19 - x20*(x1*x41 + x18*((x41)*(x41)) + x68));
auto x44 = cell[8]*x0 + x28*(x11 + x57*(x17*x7 + x23)) - x58*(x26 - x27*(x18*x64 + x66 + x68));
auto x71 = x25 + x44;
auto x45 = x14 + x29 + x34 + x35 + x43 + x70 + x71 + V{1};
auto x72 = V{1} / (x45);
auto x46 = -V{1}*x34;
auto x73 = x43 + x46;
cell[0] = x14;
cell[1] = x25;
cell[2] = x29;
cell[3] = x34;
cell[4] = x35;
cell[5] = x39;
cell[6] = x40;
cell[7] = x43;
cell[8] = x44;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x45;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x4 + x72*(-V{1}*x25 - V{1}*x29 + x70 + x73);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x6 + x72*(-V{1}*x35 - V{1}*x39 + x71 + x73);
return { x3, x8 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
