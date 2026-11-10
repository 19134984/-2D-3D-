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
struct CSE<AdvectionDiffusionEdgesDynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, 2, 1, 1>, olb::CodegenTarget::Generic> {
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
auto x10 = -V{1}*x71;
auto x11 = -x8 + x9*(x3*(x10 + x69) + V{1});
auto x12 = cell[16] - x11;
auto x13 = -x8 + x9*(x3*(x69 + x71) + V{1});
auto x14 = cell[15] - x13;
auto x15 = -x8 + x9*(x3*(x69 + x70) + V{1});
auto x16 = cell[13] - x15;
auto x17 = -x1 + x2*(x3*x69 + V{1});
auto x18 = cell[10] - x17;
auto x19 = -V{1}*x70;
auto x20 = -x8 + x9*(x3*(x19 + x4) + V{1});
auto x21 = V{1}*x16;
auto x73 = x16 + x21;
auto x22 = -x8 + x9*(x3*(x4 + x70) + V{1});
auto x23 = cell[5] - x22;
auto x24 = -x8 + x9*(x3*(x10 + x4) + V{1});
auto x25 = V{1}*x14;
auto x74 = x14 + x25;
auto x26 = -x8 + x9*(x3*(x4 + x71) + V{1});
auto x27 = -x12;
auto x28 = -x8 + x9*(x3*(x19 + x69) + V{1});
auto x29 = cell[14] - x28;
auto x30 = x12 + V{2}*x18 - V{1}*x23 - V{1}*x27 + x29 + x73 + x74;
auto x31 = -V{1}*x30;
auto x32 = -x1 + x2*(x19*x3 + V{1});
auto x33 = -x8 + x9*(x3*(x10 + x70) + V{1});
auto x34 = cell[18] - x33;
auto x35 = -x8 + x9*(x3*(x70 + x71) + V{1});
auto x36 = cell[17] - x35;
auto x37 = -x1 + x2*(x3*x70 + V{1});
auto x38 = cell[11] - x37;
auto x39 = -x8 + x9*(x3*(x10 + x19) + V{1});
auto x40 = V{1}*x36;
auto x75 = x36 + x40;
auto x41 = -x8 + x9*(x3*(x19 + x71) + V{1});
auto x42 = -x34;
auto x43 = x23 - V{1}*x29 + x34 + V{2}*x38 - V{1}*x42 + x73 + x75;
auto x44 = -V{1}*x43;
auto x45 = -x1 + x2*(x10*x3 + V{1});
auto x46 = -x1 + x2*(x3*x71 + V{1});
auto x47 = cell[12] - V{1}*cell[3] - V{1}*x12 + x27 - V{1}*x34 + x42 + V{1}*x45 - x46 + x74 + x75;
auto x48 = -V{1}*x47;
auto x49 = x3*x8;
auto x78 = x49*x6;
cell[0] = x0*(x72 + V{-1});
cell[10] = x17 + x30*x77;
cell[11] = x37 + x43*x77;
cell[12] = x46 + x47*x77;
cell[13] = x15 + x78*(x30 + x43);
cell[14] = x28 + x78*(x30 + x44);
cell[15] = x13 + x78*(x30 + x47);
cell[16] = x11 + x78*(x30 + x48);
cell[17] = x35 + x78*(x43 + x47);
cell[18] = x33 + x78*(x43 + x48);
cell[1] = x31*x77 + x5;
cell[2] = x32 + x44*x77;
cell[3] = x45 + x48*x77;
cell[4] = x20 + x78*(x31 + x44);
cell[5] = x22 + x78*(x31 + x43);
cell[6] = x24 + x78*(x31 + x48);
cell[7] = x26 + x78*(x31 + x47);
cell[8] = x39 + x78*(x44 + x48);
cell[9] = x41 + x78*(x44 + x47);
return { x72, ((x69)*(x69)) + ((x70)*(x70)) + ((x71)*(x71)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
