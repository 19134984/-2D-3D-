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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::BGK, AdvectionDiffusionExternalVelocityCollision>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x53 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x54 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x55 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x56 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{1} - x56;
auto x1 = V{0.333333333333333};
auto x2 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{1};
auto x3 = cell[0]*x0 + x1*x56*(x2 + V{-1});
auto x4 = V{3};
auto x5 = -V{1}*x53;
auto x6 = V{0.0555555555555556};
auto x7 = x2*x6;
auto x8 = cell[1]*x0 - x56*(x6 - x7*(x4*x5 + V{1}));
auto x9 = -V{1}*x54;
auto x10 = cell[2]*x0 - x56*(x6 - x7*(x4*x9 + V{1}));
auto x11 = -V{1}*x55;
auto x12 = cell[3]*x0 - x56*(x6 - x7*(x11*x4 + V{1}));
auto x13 = V{0.0277777777777778};
auto x14 = x13*x2;
auto x15 = cell[4]*x0 - x56*(x13 - x14*(x4*(x5 + x9) + V{1}));
auto x16 = cell[5]*x0 - x56*(x13 - x14*(x4*(x5 + x54) + V{1}));
auto x17 = cell[6]*x0 - x56*(x13 - x14*(x4*(x11 + x5) + V{1}));
auto x18 = cell[7]*x0 - x56*(x13 - x14*(x4*(x5 + x55) + V{1}));
auto x19 = cell[8]*x0 - x56*(x13 - x14*(x4*(x11 + x9) + V{1}));
auto x20 = cell[9]*x0 - x56*(x13 - x14*(x4*(x55 + x9) + V{1}));
auto x21 = cell[10]*x0 - x56*(x6 - x7*(x4*x53 + V{1}));
auto x22 = cell[11]*x0 - x56*(x6 - x7*(x4*x54 + V{1}));
auto x23 = cell[12]*x0 - x56*(x6 - x7*(x4*x55 + V{1}));
auto x59 = x18 + x20 + x23;
auto x24 = cell[13]*x0 - x56*(x13 - x14*(x4*(x53 + x54) + V{1}));
auto x25 = cell[14]*x0 - x56*(x13 - x14*(x4*(x53 + x9) + V{1}));
auto x26 = cell[15]*x0 - x56*(x13 - x14*(x4*(x53 + x55) + V{1}));
auto x27 = cell[16]*x0 - x56*(x13 - x14*(x4*(x11 + x53) + V{1}));
auto x57 = x21 + x25 + x27;
auto x28 = cell[17]*x0 - x56*(x13 - x14*(x4*(x54 + x55) + V{1}));
auto x29 = cell[18]*x0 - x56*(x13 - x14*(x4*(x11 + x54) + V{1}));
auto x58 = x16 + x22 + x29;
auto x30 = x10 + x12 + x15 + x17 + x19 + x24 + x26 + x28 + x3 + x57 + x58 + x59 + x8 + V{1};
auto x60 = V{1} / (x30);
auto x31 = -V{1}*x15;
auto x61 = x24 + x31;
auto x32 = -V{1}*x17;
auto x62 = x26 + x32;
auto x33 = -V{1}*x19;
auto x63 = x28 + x33;
cell[0] = x3;
cell[10] = x21;
cell[11] = x22;
cell[12] = x23;
cell[13] = x24;
cell[14] = x25;
cell[15] = x26;
cell[16] = x27;
cell[17] = x28;
cell[18] = x29;
cell[1] = x8;
cell[2] = x10;
cell[3] = x12;
cell[4] = x15;
cell[5] = x16;
cell[6] = x17;
cell[7] = x18;
cell[8] = x19;
cell[9] = x20;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x30;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x60*(-V{1}*x16 - V{1}*x18 + x57 + x61 + x62 - V{1}*x8);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x60*(-V{1}*x10 - V{1}*x20 - V{1}*x25 + x58 + x61 + x63);
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x60*(-V{1}*x12 - V{1}*x27 - V{1}*x29 + x59 + x62 + x63);
return { x2, ((x53)*(x53)) + ((x54)*(x54)) + ((x55)*(x55)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
