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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x49 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x50 = parameters.template get<descriptors::OMEGA>();
auto x51 = cell[1] + cell[8];
auto x52 = cell[4] + x51;
auto x53 = cell[5] + cell[7];
auto x54 = cell[6] + x53;
auto x55 = cell[2] + x54;
auto x0 = V{1}/x50;
auto x1 = -V{1}*cell[3];
auto x2 = -V{1}*cell[5];
auto x57 = cell[7] + x2;
auto x3 = cell[0] + cell[3] + x52 + x55 + V{1};
auto x56 = V{1} / (x3);
auto x4 = x56*(-V{1}*cell[4] + x1 + x51 + x57);
auto x58 = ((x4)*(x4));
auto x5 = V{3};
auto x60 = ((x5)*(x5));
auto x6 = V{1}*(x3 + V{-1})/x5;
auto x59 = cell[3] - x6;
auto x7 = -x3*x58 + x52 + x53 + x59;
auto x8 = -V{1}*cell[1];
auto x9 = x56*(-V{1}*cell[2] + x1 + x54 + x8);
auto x61 = ((x9)*(x9));
auto x10 = x3*x9;
auto x11 = cell[3] - x10*x4 + x57 + x8;
auto x12 = cell[1] - x10*x9 + x55 + x59;
auto x13 = V{2}*V{1} / (x0 + util::sqrt(V{2.828428}*x56*x60*((x49)*(x49))*util::sqrt(V{2}*((x11)*(x11)) + ((x12)*(x12)) + ((x7)*(x7))) + ((x0)*(x0))));
auto x14 = V{1} - x13;
auto x15 = x58 + x61;
auto x16 = V{0.5}*x15*x5;
auto x62 = V{1} - x16;
auto x17 = V{0.444444444444444};
auto x18 = cell[0]*x14 - x13*x17*(x3*(x16 + V{-1}) + V{1});
auto x19 = -V{1}*x9;
auto x20 = x19 + x4;
auto x21 = V{0.5}*x60;
auto x22 = V{0.0277777777777778};
auto x23 = x22*x3;
auto x24 = cell[1]*x14 - x13*(x22 - x23*(x20*x5 + x21*((x20)*(x20)) + x62));
auto x25 = V{0.111111111111111};
auto x26 = x25*x3;
auto x27 = cell[2]*x14 - x13*(x25 - x26*(x19*x5 + x21*((x19)*(x19)) + x62));
auto x28 = -V{1}*x4;
auto x29 = x19 + x28;
auto x30 = cell[3]*x14 - x13*(x22 - x23*(x21*((x29)*(x29)) + x29*x5 + x62));
auto x31 = cell[4]*x14 - x13*(x25 - x26*(x21*((x28)*(x28)) + x28*x5 + x62));
auto x32 = x28 + x9;
auto x33 = cell[5]*x14 - x13*(x22 - x23*(x21*((x32)*(x32)) + x32*x5 + x62));
auto x34 = cell[6]*x14 - x13*(x25 - x26*(x21*x61 + x5*x9 + x62));
auto x63 = x33 + x34;
auto x35 = x4 + x9;
auto x36 = cell[7]*x14 - x13*(x22 - x23*(x21*((x35)*(x35)) + x35*x5 + x62));
auto x37 = cell[8]*x14 - x13*(x25 - x26*(x21*x58 + x4*x5 + x62));
auto x64 = x24 + x37;
auto x38 = x18 + x27 + x30 + x31 + x36 + x63 + x64 + V{1};
auto x65 = V{1} / (x38);
auto x39 = -V{1}*x30;
auto x66 = x36 + x39;
cell[0] = x18;
cell[1] = x24;
cell[2] = x27;
cell[3] = x30;
cell[4] = x31;
cell[5] = x33;
cell[6] = x34;
cell[7] = x36;
cell[8] = x37;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x38;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x65*(-V{1}*x24 - V{1}*x27 + x63 + x66);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x65*(-V{1}*x31 - V{1}*x33 + x64 + x66);
return { x3, x15 };
}
};

}

}

// Generation Info: commit=61323934e549f5f16850d6c6280b91bc057cb6c3
