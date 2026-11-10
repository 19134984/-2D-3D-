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
struct CSE<AdvectionDiffusionCornerDynamics3D<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, -1, -1, 1>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x71 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x72 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x73 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x74 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x78 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{0.333333333333333};
auto x1 = V{0.0555555555555556};
auto x2 = x1*x74;
auto x3 = V{3};
auto x4 = -V{1}*x71;
auto x5 = -x1 + x2*(x3*x4 + V{1});
auto x6 = V{1} - x78;
auto x7 = x1*x3;
auto x79 = x6*x7;
auto x8 = cell[1] - x5;
auto x9 = -V{1}*x72;
auto x10 = V{0.0277777777777778};
auto x11 = x10*x74;
auto x12 = -x10 + x11*(x3*(x4 + x9) + V{1});
auto x13 = cell[4] - x12;
auto x14 = -V{1}*x13;
auto x15 = -x10 + x11*(x3*(x4 + x72) + V{1});
auto x16 = cell[5] - x15;
auto x17 = -V{1}*x73;
auto x18 = -x10 + x11*(x3*(x17 + x4) + V{1});
auto x19 = -V{1}*cell[6] + V{1}*x18;
auto x20 = -x10 + x11*(x3*(x4 + x73) + V{1});
auto x21 = cell[7] - x20;
auto x22 = -x1 + x2*(x3*x71 + V{1});
auto x23 = -x10 + x11*(x3*(x71 + x72) + V{1});
auto x24 = -x13;
auto x75 = x14 + x24;
auto x25 = -x10 + x11*(x3*(x71 + x9) + V{1});
auto x26 = cell[14] - x25;
auto x27 = -x10 + x11*(x3*(x71 + x73) + V{1});
auto x28 = cell[15] - x27;
auto x76 = x19 + x28;
auto x29 = -x10 + x11*(x3*(x17 + x71) + V{1});
auto x30 = -x21;
auto x31 = -V{1}*x16 - V{1}*x21 + x26 + x30 + x75 + x76 - V{2}*x8;
auto x32 = -V{1}*x31;
auto x33 = -x1 + x2*(x3*x9 + V{1});
auto x34 = cell[2] - x33;
auto x35 = -x10 + x11*(x3*(x17 + x9) + V{1});
auto x36 = -V{1}*cell[8] + V{1}*x35;
auto x37 = -x10 + x11*(x3*(x73 + x9) + V{1});
auto x38 = cell[9] - x37;
auto x39 = -x1 + x2*(x3*x72 + V{1});
auto x40 = -x10 + x11*(x3*(x72 + x73) + V{1});
auto x41 = cell[17] - x40;
auto x77 = x36 + x41;
auto x42 = -x10 + x11*(x3*(x17 + x72) + V{1});
auto x43 = -x38;
auto x44 = x16 - V{1}*x26 - V{2}*x34 - V{1}*x38 + x43 + x75 + x77;
auto x45 = -V{1}*x44;
auto x46 = -x1 + x2*(x17*x3 + V{1});
auto x47 = -x1 + x2*(x3*x73 + V{1});
auto x48 = cell[12] - x47;
auto x49 = x21 - V{1}*x30 + x38 - V{1}*x43 + V{2}*x48 + x76 + x77;
auto x50 = -V{1}*x49;
auto x51 = x10*x3;
auto x80 = x51*x6;
cell[0] = x0*(x74 + V{-1});
cell[10] = x22 + x31*x79;
cell[11] = x39 + x44*x79;
cell[12] = x47 + x49*x79;
cell[13] = x23 + x80*(x31 + x44);
cell[14] = x25 + x80*(x31 + x45);
cell[15] = x27 + x80*(x31 + x49);
cell[16] = x29 + x80*(x31 + x50);
cell[17] = x40 + x80*(x44 + x49);
cell[18] = x42 + x80*(x44 + x50);
cell[1] = x32*x79 + x5;
cell[2] = x33 + x45*x79;
cell[3] = x46 + x50*x79;
cell[4] = x12 + x80*(x32 + x45);
cell[5] = x15 + x80*(x32 + x44);
cell[6] = x18 + x80*(x32 + x50);
cell[7] = x20 + x80*(x32 + x49);
cell[8] = x35 + x80*(x45 + x50);
cell[9] = x37 + x80*(x45 + x49);
return { x74, ((x71)*(x71)) + ((x72)*(x72)) + ((x73)*(x73)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
