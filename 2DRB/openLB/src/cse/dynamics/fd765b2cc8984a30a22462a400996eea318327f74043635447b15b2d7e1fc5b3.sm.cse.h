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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::CahnHilliardZerothOrder, collision::BGK, forcing::WellBalancedCahnHilliard>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x49 = cell.template getFieldComponent<olb::descriptors::CHEM_POTENTIAL>(0);
auto x50 = cell.template getFieldComponent<olb::descriptors::SOURCE>(0);
auto x51 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x52 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x53 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x54 = parameters.template get<descriptors::OMEGA>();
auto x55 = x54*(x49 + V{-1});
auto x0 = V{1} - x54;
auto x1 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{1};
auto x2 = V{0.333333333333333};
auto x3 = V{-2};
auto x4 = cell[0]*x0 + x2*x3*x50 - x54*(-x1 + x2 - x49*(x2 + V{-1}));
auto x5 = V{0.0555555555555556};
auto x6 = x5*x55;
auto x7 = V{1.5};
auto x8 = x5*x50*(x3 + x7);
auto x56 = x6 + x8;
auto x9 = cell[1]*x0 + x56;
auto x10 = cell[2]*x0 + x56;
auto x11 = cell[3]*x0 + x56;
auto x12 = V{0.0277777777777778};
auto x13 = x12*x55;
auto x14 = x12*x50*(x3 + V{2}*x7);
auto x57 = x13 + x14;
auto x15 = cell[4]*x0 + x57;
auto x16 = cell[5]*x0 + x57;
auto x17 = cell[6]*x0 + x57;
auto x18 = cell[7]*x0 + x57;
auto x19 = cell[8]*x0 + x57;
auto x20 = cell[9]*x0 + x57;
auto x21 = cell[10]*x0 + x56;
auto x22 = cell[11]*x0 + x56;
auto x23 = cell[12]*x0 + x56;
auto x24 = cell[13]*x0 + x57;
auto x25 = cell[14]*x0 + x57;
auto x26 = cell[15]*x0 + x57;
auto x27 = cell[16]*x0 + x57;
auto x28 = cell[17]*x0 + x57;
auto x29 = cell[18]*x0 + x57;
cell[0] = x4;
cell[10] = x21;
cell[11] = x22;
cell[12] = x23;
cell[13] = x24;
cell[14] = x25;
cell[15] = x26;
cell[16] = x27;
cell[17] = x28;
cell[18] = x29;
cell[1] = x9;
cell[2] = x10;
cell[3] = x11;
cell[4] = x15;
cell[5] = x16;
cell[6] = x17;
cell[7] = x18;
cell[8] = x19;
cell[9] = x20;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x10 + x11 + x15 + x16 + x17 + x18 + x19 + x20 + x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 + x29 + x4 + x9 + V{1};
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x51;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x52;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x53;
return { x1, ((x51)*(x51)) + ((x52)*(x52)) + ((x53)*(x53)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
