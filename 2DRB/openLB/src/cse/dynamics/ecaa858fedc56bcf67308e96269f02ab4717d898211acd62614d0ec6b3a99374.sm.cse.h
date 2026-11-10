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
struct CSE<AdvectionDiffusionEdgesDynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, 2, -1, -1>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x72 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x73 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x74 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x75 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x76 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x77 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x78 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x79 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{0.333333333333333};
auto x1 = V{3};
auto x2 = -V{1}*x72;
auto x3 = V{0.0555555555555556};
auto x4 = x3*x75;
auto x5 = -x3 + x4*(x1*x2 + V{1});
auto x6 = V{1} - x79;
auto x7 = x1*x3;
auto x83 = x6*x7;
auto x8 = cell[1] - x5;
auto x9 = -V{1}*x73;
auto x10 = V{0.0277777777777778};
auto x11 = x10*x75;
auto x12 = -x10 + x11*(x1*(x2 + x9) + V{1});
auto x13 = cell[4] - x12;
auto x14 = -V{1}*x13;
auto x15 = -x10 + x11*(x1*(x2 + x73) + V{1});
auto x16 = cell[5] - x15;
auto x17 = -V{1}*x74;
auto x18 = -x10 + x11*(x1*(x17 + x2) + V{1});
auto x19 = cell[6] - x18;
auto x20 = -V{1}*x19;
auto x21 = -x10 + x11*(x1*(x2 + x74) + V{1});
auto x22 = cell[7] - x21;
auto x23 = -x3 + x4*(x1*x72 + V{1});
auto x24 = -x10 + x11*(x1*(x72 + x73) + V{1});
auto x25 = -x13;
auto x80 = x14 + x25;
auto x26 = -x10 + x11*(x1*(x72 + x9) + V{1});
auto x27 = cell[14] - x26;
auto x28 = -x10 + x11*(x1*(x72 + x74) + V{1});
auto x29 = -x19;
auto x81 = x20 + x29;
auto x30 = -x10 + x11*(x1*(x17 + x72) + V{1});
auto x31 = -x22;
auto x32 = -V{1}*x16 - V{1}*x22 + x27 + x31 - V{2}*x8 + x80 + x81;
auto x33 = -V{1}*x32;
auto x34 = -x3 + x4*(x1*x9 + V{1});
auto x35 = cell[2] - x34;
auto x36 = -x10 + x11*(x1*(x17 + x9) + V{1});
auto x37 = cell[8] - x36;
auto x38 = -V{1}*x37;
auto x39 = -x10 + x11*(x1*(x74 + x9) + V{1});
auto x40 = cell[9] - x39;
auto x41 = -x3 + x4*(x1*x73 + V{1});
auto x42 = -x10 + x11*(x1*(x73 + x74) + V{1});
auto x43 = -x37;
auto x82 = x38 + x43;
auto x44 = -x10 + x11*(x1*(x17 + x73) + V{1});
auto x45 = -x40;
auto x46 = x16 - V{1}*x27 - V{2}*x35 - V{1}*x40 + x45 + x80 + x82;
auto x47 = -V{1}*x46;
auto x48 = -x3 + x4*(x1*x17 + V{1});
auto x49 = -x3 + x4*(x1*x74 + V{1});
auto x50 = cell[12] - V{1}*cell[3] + x22 - V{1}*x31 + x40 - V{1}*x45 + V{1}*x48 - x49 + x81 + x82;
auto x51 = -V{1}*x50;
auto x52 = x1*x10;
auto x84 = x52*x6;
cell[0] = x0*(x75 + V{-1});
cell[10] = x23 + x32*x83;
cell[11] = x41 + x46*x83;
cell[12] = x49 + x50*x83;
cell[13] = x24 + x84*(x32 + x46);
cell[14] = x26 + x84*(x32 + x47);
cell[15] = x28 + x84*(x32 + x50);
cell[16] = x30 + x84*(x32 + x51);
cell[17] = x42 + x84*(x46 + x50);
cell[18] = x44 + x84*(x46 + x51);
cell[1] = x33*x83 + x5;
cell[2] = x34 + x47*x83;
cell[3] = x48 + x51*x83;
cell[4] = x12 + x84*(x33 + x47);
cell[5] = x15 + x84*(x33 + x46);
cell[6] = x18 + x84*(x33 + x51);
cell[7] = x21 + x84*(x33 + x50);
cell[8] = x36 + x84*(x47 + x51);
cell[9] = x39 + x84*(x47 + x50);
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x75;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x76;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x77;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x78;
return { x75, ((x72)*(x72)) + ((x73)*(x73)) + ((x74)*(x74)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
