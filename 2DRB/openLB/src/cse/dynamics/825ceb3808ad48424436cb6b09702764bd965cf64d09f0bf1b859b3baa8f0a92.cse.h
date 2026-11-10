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
struct CSE<AdvectionDiffusionEdgesDynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, 1, -1, 1>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x70 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x71 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x72 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x73 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x77 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{0.333333333333333};
auto x1 = V{0.0555555555555556};
auto x2 = x1*x73;
auto x3 = V{3};
auto x4 = -V{1}*x70;
auto x5 = -x1 + x2*(x3*x4 + V{1});
auto x6 = V{1} - x77;
auto x7 = x1*x3;
auto x79 = x6*x7;
auto x8 = cell[1] - x5;
auto x9 = V{0.0277777777777778};
auto x74 = -x9;
auto x10 = x73*x9;
auto x11 = -V{1}*x71;
auto x12 = x10*(x3*(x11 + x4) + V{1}) + x74;
auto x13 = cell[4] - x12;
auto x14 = -V{1}*x13;
auto x15 = x10*(x3*(x4 + x71) + V{1}) + x74;
auto x16 = cell[5] - x15;
auto x17 = -V{1}*x72;
auto x18 = x10*(x3*(x17 + x4) + V{1}) + x74;
auto x19 = -V{1}*cell[6] + V{1}*x18;
auto x20 = x10*(x3*(x4 + x72) + V{1}) + x74;
auto x21 = cell[7] - x20;
auto x22 = -x1 + x2*(x3*x70 + V{1});
auto x23 = x10*(x3*(x70 + x71) + V{1}) + x74;
auto x24 = -x13;
auto x76 = x14 + x24;
auto x25 = x10*(x3*(x11 + x70) + V{1}) + x74;
auto x26 = -x16;
auto x27 = x10*(x3*(x70 + x72) + V{1}) + x74;
auto x28 = cell[15] - x27;
auto x75 = x19 + x28;
auto x29 = x10*(x3*(x17 + x70) + V{1}) + x74;
auto x30 = -x21;
auto x31 = -V{1}*x16 - V{1}*x21 + x26 + x30 + x75 + x76 - V{2}*x8;
auto x32 = -V{1}*x31;
auto x33 = -x1 + x2*(x11*x3 + V{1});
auto x34 = x10*(x3*(x71 + x72) + V{1}) + x74;
auto x35 = cell[17] - x34;
auto x36 = x10*(x3*(x11 + x17) + V{1}) + x74;
auto x37 = V{1}*x35;
auto x78 = x35 + x37;
auto x38 = x10*(x3*(x11 + x72) + V{1}) + x74;
auto x39 = cell[9] - x38;
auto x40 = -x1 + x2*(x3*x71 + V{1});
auto x41 = x10*(x3*(x17 + x71) + V{1}) + x74;
auto x42 = -x39;
auto x43 = cell[11] - V{1}*cell[2] + x16 - V{1}*x26 + V{1}*x33 - V{1}*x39 - x40 + x42 + x76 + x78;
auto x44 = -V{1}*x43;
auto x45 = -x1 + x2*(x17*x3 + V{1});
auto x46 = -x1 + x2*(x3*x72 + V{1});
auto x47 = cell[12] - x46;
auto x48 = x21 - V{1}*x30 + x39 - V{1}*x42 + V{2}*x47 + x75 + x78;
auto x49 = -V{1}*x48;
auto x50 = x3*x9;
auto x80 = x50*x6;
cell[0] = x0*(x73 + V{-1});
cell[10] = x22 + x31*x79;
cell[11] = x40 + x43*x79;
cell[12] = x46 + x48*x79;
cell[13] = x23 + x80*(x31 + x43);
cell[14] = x25 + x80*(x31 + x44);
cell[15] = x27 + x80*(x31 + x48);
cell[16] = x29 + x80*(x31 + x49);
cell[17] = x34 + x80*(x43 + x48);
cell[18] = x41 + x80*(x43 + x49);
cell[1] = x32*x79 + x5;
cell[2] = x33 + x44*x79;
cell[3] = x45 + x49*x79;
cell[4] = x12 + x80*(x32 + x44);
cell[5] = x15 + x80*(x32 + x43);
cell[6] = x18 + x80*(x32 + x49);
cell[7] = x20 + x80*(x32 + x48);
cell[8] = x36 + x80*(x44 + x49);
cell[9] = x38 + x80*(x44 + x48);
return { x73, ((x70)*(x70)) + ((x71)*(x71)) + ((x72)*(x72)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
