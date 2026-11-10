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
struct CSE<AdvectionDiffusionEdgesDynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, 1, 1, 1>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x69 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x70 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x71 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x72 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x76 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{0.333333333333333};
auto x1 = V{0.0555555555555556};
auto x2 = x1*x72;
auto x3 = V{3};
auto x4 = -V{1}*x69;
auto x5 = -x1 + x2*(x3*x4 + V{1});
auto x6 = V{1} - x76;
auto x7 = x1*x3;
auto x77 = x6*x7;
auto x8 = V{0.0277777777777778};
auto x9 = x72*x8;
auto x10 = -x8 + x9*(x3*(x69 + x71) + V{1});
auto x11 = cell[15] - x10;
auto x12 = -V{1}*x70;
auto x13 = -x8 + x9*(x3*(x12 + x69) + V{1});
auto x14 = cell[14] - x13;
auto x15 = -x8 + x9*(x3*(x69 + x70) + V{1});
auto x16 = cell[13] - x15;
auto x17 = -x1 + x2*(x3*x69 + V{1});
auto x18 = cell[10] - x17;
auto x19 = -x8 + x9*(x3*(x12 + x4) + V{1});
auto x20 = V{1}*x16;
auto x74 = x16 + x20;
auto x21 = -x8 + x9*(x3*(x4 + x70) + V{1});
auto x22 = -x14;
auto x23 = -V{1}*x71;
auto x24 = -x8 + x9*(x3*(x23 + x4) + V{1});
auto x25 = V{1}*x11;
auto x73 = x11 + x25;
auto x26 = -x8 + x9*(x3*(x4 + x71) + V{1});
auto x27 = cell[7] - x26;
auto x28 = -x8 + x9*(x3*(x23 + x69) + V{1});
auto x29 = cell[16] - x28;
auto x30 = x14 + V{2}*x18 - V{1}*x22 - V{1}*x27 + x29 + x73 + x74;
auto x31 = -V{1}*x30;
auto x32 = -x1 + x2*(x12*x3 + V{1});
auto x33 = -x8 + x9*(x3*(x70 + x71) + V{1});
auto x34 = cell[17] - x33;
auto x35 = -x8 + x9*(x3*(x12 + x23) + V{1});
auto x36 = V{1}*x34;
auto x75 = x34 + x36;
auto x37 = -x8 + x9*(x3*(x12 + x71) + V{1});
auto x38 = cell[9] - x37;
auto x39 = -x1 + x2*(x3*x70 + V{1});
auto x40 = -x8 + x9*(x3*(x23 + x70) + V{1});
auto x41 = -x38;
auto x42 = cell[11] - V{1}*cell[2] - V{1}*x14 + x22 + V{1}*x32 - V{1}*x38 - x39 + x41 + x74 + x75;
auto x43 = -V{1}*x42;
auto x44 = -x1 + x2*(x23*x3 + V{1});
auto x45 = -x1 + x2*(x3*x71 + V{1});
auto x46 = cell[12] - x45;
auto x47 = x27 - V{1}*x29 + x38 - V{1}*x41 + V{2}*x46 + x73 + x75;
auto x48 = -V{1}*x47;
auto x49 = x3*x8;
auto x78 = x49*x6;
cell[0] = x0*(x72 + V{-1});
cell[10] = x17 + x30*x77;
cell[11] = x39 + x42*x77;
cell[12] = x45 + x47*x77;
cell[13] = x15 + x78*(x30 + x42);
cell[14] = x13 + x78*(x30 + x43);
cell[15] = x10 + x78*(x30 + x47);
cell[16] = x28 + x78*(x30 + x48);
cell[17] = x33 + x78*(x42 + x47);
cell[18] = x40 + x78*(x42 + x48);
cell[1] = x31*x77 + x5;
cell[2] = x32 + x43*x77;
cell[3] = x44 + x48*x77;
cell[4] = x19 + x78*(x31 + x43);
cell[5] = x21 + x78*(x31 + x42);
cell[6] = x24 + x78*(x31 + x48);
cell[7] = x26 + x78*(x31 + x47);
cell[8] = x35 + x78*(x43 + x48);
cell[9] = x37 + x78*(x43 + x47);
return { x72, ((x69)*(x69)) + ((x70)*(x70)) + ((x71)*(x71)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
