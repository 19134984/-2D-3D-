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
struct CSE<AdvectionDiffusionCornerDynamics3D<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, -1, 1, 1>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x70 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x71 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x72 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x73 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x74 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x75 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x76 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x77 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{0.333333333333333};
auto x1 = V{0.0555555555555556};
auto x2 = x1*x73;
auto x3 = V{3};
auto x4 = -V{1}*x70;
auto x5 = -x1 + x2*(x3*x4 + V{1});
auto x6 = V{1} - x77;
auto x7 = x1*x3;
auto x82 = x6*x7;
auto x8 = cell[1] - x5;
auto x9 = V{0.0277777777777778};
auto x78 = -x9;
auto x10 = x73*x9;
auto x11 = -V{1}*x71;
auto x12 = x10*(x3*(x11 + x4) + V{1}) + x78;
auto x13 = -V{1}*cell[4] + V{1}*x12;
auto x14 = x10*(x3*(x4 + x71) + V{1}) + x78;
auto x15 = cell[5] - x14;
auto x16 = -V{1}*x72;
auto x17 = x10*(x3*(x16 + x4) + V{1}) + x78;
auto x18 = -V{1}*cell[6] + V{1}*x17;
auto x19 = x10*(x3*(x4 + x72) + V{1}) + x78;
auto x20 = cell[7] - x19;
auto x21 = -x1 + x2*(x3*x70 + V{1});
auto x22 = x10*(x3*(x70 + x71) + V{1}) + x78;
auto x23 = cell[13] - x22;
auto x79 = x13 + x23;
auto x24 = x10*(x3*(x11 + x70) + V{1}) + x78;
auto x25 = -x15;
auto x26 = x10*(x3*(x70 + x72) + V{1}) + x78;
auto x27 = cell[15] - x26;
auto x80 = x18 + x27;
auto x28 = x10*(x3*(x16 + x70) + V{1}) + x78;
auto x29 = -x20;
auto x30 = -V{1}*x15 - V{1}*x20 + x25 + x29 + x79 - V{2}*x8 + x80;
auto x31 = -V{1}*x30;
auto x32 = -x1 + x2*(x11*x3 + V{1});
auto x33 = x10*(x3*(x71 + x72) + V{1}) + x78;
auto x34 = cell[17] - x33;
auto x35 = -x1 + x2*(x3*x71 + V{1});
auto x36 = cell[11] - x35;
auto x37 = x10*(x3*(x11 + x16) + V{1}) + x78;
auto x38 = V{1}*x34;
auto x81 = x34 + x38;
auto x39 = x10*(x3*(x11 + x72) + V{1}) + x78;
auto x40 = cell[9] - x39;
auto x41 = x10*(x3*(x16 + x71) + V{1}) + x78;
auto x42 = cell[18] - x41;
auto x43 = x15 - V{1}*x25 + V{2}*x36 - V{1}*x40 + x42 + x79 + x81;
auto x44 = -V{1}*x43;
auto x45 = -x1 + x2*(x16*x3 + V{1});
auto x46 = -x1 + x2*(x3*x72 + V{1});
auto x47 = cell[12] - x46;
auto x48 = x20 - V{1}*x29 + x40 - V{1}*x42 + V{2}*x47 + x80 + x81;
auto x49 = -V{1}*x48;
auto x50 = x3*x9;
auto x83 = x50*x6;
cell[0] = x0*(x73 + V{-1});
cell[10] = x21 + x30*x82;
cell[11] = x35 + x43*x82;
cell[12] = x46 + x48*x82;
cell[13] = x22 + x83*(x30 + x43);
cell[14] = x24 + x83*(x30 + x44);
cell[15] = x26 + x83*(x30 + x48);
cell[16] = x28 + x83*(x30 + x49);
cell[17] = x33 + x83*(x43 + x48);
cell[18] = x41 + x83*(x43 + x49);
cell[1] = x31*x82 + x5;
cell[2] = x32 + x44*x82;
cell[3] = x45 + x49*x82;
cell[4] = x12 + x83*(x31 + x44);
cell[5] = x14 + x83*(x31 + x43);
cell[6] = x17 + x83*(x31 + x49);
cell[7] = x19 + x83*(x31 + x48);
cell[8] = x37 + x83*(x44 + x49);
cell[9] = x39 + x83*(x44 + x48);
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x73;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x74;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x75;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x76;
return { x73, ((x70)*(x70)) + ((x71)*(x71)) + ((x72)*(x72)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
