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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::BGK, AdvectionDiffusionExternalVelocityCollision>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x30 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x31 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x32 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x33 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{1} - x33;
auto x1 = V{0.333333333333333};
auto x2 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{1};
auto x3 = V{3};
auto x4 = -V{1}*x30;
auto x5 = V{0.0555555555555556};
auto x6 = x2*x5;
auto x7 = -V{1}*x31;
auto x8 = -V{1}*x32;
auto x9 = V{0.0277777777777778};
auto x34 = -x9;
auto x10 = x2*x9;
auto x11 = cell[0]*x0 + x1*x33*(x2 + V{-1});
auto x12 = cell[10]*x0 - x33*(x5 - x6*(x3*x30 + V{1}));
auto x13 = cell[11]*x0 - x33*(x5 - x6*(x3*x31 + V{1}));
auto x14 = cell[12]*x0 - x33*(x5 - x6*(x3*x32 + V{1}));
auto x15 = cell[13]*x0 + x33*(x10*(x3*(x30 + x31) + V{1}) + x34);
auto x16 = cell[14]*x0 + x33*(x10*(x3*(x30 + x7) + V{1}) + x34);
auto x17 = cell[15]*x0 + x33*(x10*(x3*(x30 + x32) + V{1}) + x34);
auto x18 = cell[16]*x0 + x33*(x10*(x3*(x30 + x8) + V{1}) + x34);
auto x19 = cell[17]*x0 + x33*(x10*(x3*(x31 + x32) + V{1}) + x34);
auto x20 = cell[18]*x0 + x33*(x10*(x3*(x31 + x8) + V{1}) + x34);
auto x21 = cell[1]*x0 - x33*(x5 - x6*(x3*x4 + V{1}));
auto x22 = cell[2]*x0 - x33*(x5 - x6*(x3*x7 + V{1}));
auto x23 = cell[3]*x0 - x33*(x5 - x6*(x3*x8 + V{1}));
auto x24 = cell[4]*x0 + x33*(x10*(x3*(x4 + x7) + V{1}) + x34);
auto x25 = cell[5]*x0 + x33*(x10*(x3*(x31 + x4) + V{1}) + x34);
auto x26 = cell[6]*x0 + x33*(x10*(x3*(x4 + x8) + V{1}) + x34);
auto x27 = cell[7]*x0 + x33*(x10*(x3*(x32 + x4) + V{1}) + x34);
auto x28 = cell[8]*x0 + x33*(x10*(x3*(x7 + x8) + V{1}) + x34);
auto x29 = cell[9]*x0 + x33*(x10*(x3*(x32 + x7) + V{1}) + x34);
cell[0] = x11;
cell[10] = x12;
cell[11] = x13;
cell[12] = x14;
cell[13] = x15;
cell[14] = x16;
cell[15] = x17;
cell[16] = x18;
cell[17] = x19;
cell[18] = x20;
cell[1] = x21;
cell[2] = x22;
cell[3] = x23;
cell[4] = x24;
cell[5] = x25;
cell[6] = x26;
cell[7] = x27;
cell[8] = x28;
cell[9] = x29;
return { x2, ((x30)*(x30)) + ((x31)*(x31)) + ((x32)*(x32)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
