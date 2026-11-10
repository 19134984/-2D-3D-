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
struct CSE<AdvectionDiffusionBoundariesDynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, 1, -1>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x27 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x28 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x29 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x30 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x0 = V{0.0277777777777778};
auto x31 = V{2}*x0;
auto x1 = x0*x30;
auto x2 = V{3};
auto x3 = -V{1}*x27;
auto x4 = -V{1}*x28;
auto x5 = V{0.0555555555555556};
auto x6 = x30*x5;
auto x7 = -V{1}*x29;
auto x8 = cell[0];
auto x9 = cell[10];
auto x10 = -cell[2] - V{2}*x5 + x6*(x2*x28 + V{1}) + x6*(x2*x4 + V{1});
auto x11 = cell[12];
auto x12 = -cell[4] + x1*(x2*(x27 + x28) + V{1}) + x1*(x2*(x3 + x4) + V{1}) - x31;
auto x13 = cell[14];
auto x14 = cell[15];
auto x15 = cell[16];
auto x16 = -cell[8] + x1*(x2*(x28 + x29) + V{1}) + x1*(x2*(x4 + x7) + V{1}) - x31;
auto x17 = -cell[9] + x1*(x2*(x28 + x7) + V{1}) + x1*(x2*(x29 + x4) + V{1}) - x31;
auto x18 = cell[1];
auto x19 = cell[2];
auto x20 = cell[3];
auto x21 = cell[4];
auto x22 = -cell[14] + x1*(x2*(x27 + x4) + V{1}) + x1*(x2*(x28 + x3) + V{1}) - x31;
auto x23 = cell[6];
auto x24 = cell[7];
auto x25 = cell[8];
auto x26 = cell[9];
cell[0] = x8;
cell[10] = x9;
cell[11] = x10;
cell[12] = x11;
cell[13] = x12;
cell[14] = x13;
cell[15] = x14;
cell[16] = x15;
cell[17] = x16;
cell[18] = x17;
cell[1] = x18;
cell[2] = x19;
cell[3] = x20;
cell[4] = x21;
cell[5] = x22;
cell[6] = x23;
cell[7] = x24;
cell[8] = x25;
cell[9] = x26;
return { V{-1}, V{-1} };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
