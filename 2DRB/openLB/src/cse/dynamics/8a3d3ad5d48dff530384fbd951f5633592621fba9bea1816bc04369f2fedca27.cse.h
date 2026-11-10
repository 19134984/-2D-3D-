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
struct CSE<AdvectionDiffusionEdgesDynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, 0, -1, -1>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x72 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x73 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x74 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x75 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x79 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{0.333333333333333};
auto x1 = V{0.0555555555555556};
auto x2 = x1*x75;
auto x3 = V{3};
auto x4 = -V{1}*x72;
auto x5 = -x1 + x2*(x3*x4 + V{1});
auto x6 = V{1} - x79;
auto x7 = x1*x3;
auto x80 = x6*x7;
auto x8 = V{0.0277777777777778};
auto x9 = x75*x8;
auto x10 = -V{1}*x74;
auto x11 = -x8 + x9*(x3*(x10 + x72) + V{1});
auto x12 = cell[16] - x11;
auto x13 = -V{1}*x73;
auto x14 = -x8 + x9*(x3*(x13 + x72) + V{1});
auto x15 = cell[14] - x14;
auto x16 = -x8 + x9*(x3*(x13 + x4) + V{1});
auto x17 = cell[4] - x16;
auto x18 = -V{1}*x17;
auto x19 = -x8 + x9*(x3*(x4 + x73) + V{1});
auto x20 = -x15;
auto x21 = -x8 + x9*(x3*(x10 + x4) + V{1});
auto x22 = cell[6] - x21;
auto x23 = -V{1}*x22;
auto x24 = -x8 + x9*(x3*(x4 + x74) + V{1});
auto x25 = -x12;
auto x26 = -x1 + x2*(x3*x72 + V{1});
auto x27 = -x8 + x9*(x3*(x72 + x73) + V{1});
auto x28 = -x17;
auto x76 = x18 + x28;
auto x29 = -x8 + x9*(x3*(x72 + x74) + V{1});
auto x30 = -x22;
auto x77 = x23 + x30;
auto x31 = cell[10] - V{1}*cell[1] + x12 + x15 - V{1}*x20 - V{1}*x25 - x26 + V{1}*x5 + x76 + x77;
auto x32 = -V{1}*x31;
auto x33 = -x1 + x2*(x13*x3 + V{1});
auto x34 = cell[2] - x33;
auto x35 = -x8 + x9*(x3*(x10 + x13) + V{1});
auto x36 = cell[8] - x35;
auto x37 = -V{1}*x36;
auto x38 = -x8 + x9*(x3*(x13 + x74) + V{1});
auto x39 = cell[9] - x38;
auto x40 = -x1 + x2*(x3*x73 + V{1});
auto x41 = -x8 + x9*(x3*(x73 + x74) + V{1});
auto x42 = -x36;
auto x78 = x37 + x42;
auto x43 = -x8 + x9*(x3*(x10 + x73) + V{1});
auto x44 = cell[18] - x43;
auto x45 = -V{1}*x15 + x20 - V{2}*x34 - V{1}*x39 + x44 + x76 + x78;
auto x46 = -V{1}*x45;
auto x47 = -x1 + x2*(x10*x3 + V{1});
auto x48 = cell[3] - x47;
auto x49 = -x1 + x2*(x3*x74 + V{1});
auto x50 = -V{1}*x12 + x25 + x39 - V{1}*x44 - V{2}*x48 + x77 + x78;
auto x51 = -V{1}*x50;
auto x52 = x3*x8;
auto x81 = x52*x6;
cell[0] = x0*(x75 + V{-1});
cell[10] = x26 + x31*x80;
cell[11] = x40 + x45*x80;
cell[12] = x49 + x50*x80;
cell[13] = x27 + x81*(x31 + x45);
cell[14] = x14 + x81*(x31 + x46);
cell[15] = x29 + x81*(x31 + x50);
cell[16] = x11 + x81*(x31 + x51);
cell[17] = x41 + x81*(x45 + x50);
cell[18] = x43 + x81*(x45 + x51);
cell[1] = x32*x80 + x5;
cell[2] = x33 + x46*x80;
cell[3] = x47 + x51*x80;
cell[4] = x16 + x81*(x32 + x46);
cell[5] = x19 + x81*(x32 + x45);
cell[6] = x21 + x81*(x32 + x51);
cell[7] = x24 + x81*(x32 + x50);
cell[8] = x35 + x81*(x46 + x51);
cell[9] = x38 + x81*(x46 + x50);
return { x75, ((x72)*(x72)) + ((x73)*(x73)) + ((x74)*(x74)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
