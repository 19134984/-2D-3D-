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
struct CSE<AdvectionDiffusionEdgesDynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, 0, 1, -1>, olb::CodegenTarget::Generic> {
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
auto x78 = x6*x7;
auto x8 = V{0.0277777777777778};
auto x9 = x73*x8;
auto x10 = -V{1}*x72;
auto x11 = -x8 + x9*(x3*(x10 + x70) + V{1});
auto x12 = cell[16] - x11;
auto x13 = -x8 + x9*(x3*(x70 + x71) + V{1});
auto x14 = cell[13] - x13;
auto x15 = -V{1}*x71;
auto x16 = -x8 + x9*(x3*(x15 + x4) + V{1});
auto x17 = V{1}*x14;
auto x74 = x14 + x17;
auto x18 = -x8 + x9*(x3*(x4 + x71) + V{1});
auto x19 = cell[5] - x18;
auto x20 = -x8 + x9*(x3*(x10 + x4) + V{1});
auto x21 = cell[6] - x20;
auto x22 = -V{1}*x21;
auto x23 = -x8 + x9*(x3*(x4 + x72) + V{1});
auto x24 = -x12;
auto x25 = -x1 + x2*(x3*x70 + V{1});
auto x26 = -x8 + x9*(x3*(x15 + x70) + V{1});
auto x27 = -x19;
auto x28 = -x8 + x9*(x3*(x70 + x72) + V{1});
auto x29 = -x21;
auto x75 = x22 + x29;
auto x30 = cell[10] - V{1}*cell[1] + x12 - V{1}*x19 - V{1}*x24 - x25 + x27 + V{1}*x5 + x74 + x75;
auto x31 = -V{1}*x30;
auto x32 = -x1 + x2*(x15*x3 + V{1});
auto x33 = -x8 + x9*(x3*(x10 + x71) + V{1});
auto x34 = cell[18] - x33;
auto x35 = -x1 + x2*(x3*x71 + V{1});
auto x36 = cell[11] - x35;
auto x37 = -x8 + x9*(x3*(x10 + x15) + V{1});
auto x38 = -V{1}*cell[8] + V{1}*x37;
auto x39 = -x8 + x9*(x3*(x15 + x72) + V{1});
auto x40 = -x34;
auto x41 = -x8 + x9*(x3*(x71 + x72) + V{1});
auto x42 = cell[17] - x41;
auto x76 = x38 + x42;
auto x43 = x19 - V{1}*x27 + x34 + V{2}*x36 - V{1}*x40 + x74 + x76;
auto x44 = -V{1}*x43;
auto x45 = -x1 + x2*(x10*x3 + V{1});
auto x46 = cell[3] - x45;
auto x47 = -x1 + x2*(x3*x72 + V{1});
auto x48 = -V{1}*x12 + x24 - V{1}*x34 + x40 - V{2}*x46 + x75 + x76;
auto x49 = -V{1}*x48;
auto x50 = x3*x8;
auto x79 = x50*x6;
cell[0] = x0*(x73 + V{-1});
cell[10] = x25 + x30*x78;
cell[11] = x35 + x43*x78;
cell[12] = x47 + x48*x78;
cell[13] = x13 + x79*(x30 + x43);
cell[14] = x26 + x79*(x30 + x44);
cell[15] = x28 + x79*(x30 + x48);
cell[16] = x11 + x79*(x30 + x49);
cell[17] = x41 + x79*(x43 + x48);
cell[18] = x33 + x79*(x43 + x49);
cell[1] = x31*x78 + x5;
cell[2] = x32 + x44*x78;
cell[3] = x45 + x49*x78;
cell[4] = x16 + x79*(x31 + x44);
cell[5] = x18 + x79*(x31 + x43);
cell[6] = x20 + x79*(x31 + x49);
cell[7] = x23 + x79*(x31 + x48);
cell[8] = x37 + x79*(x44 + x49);
cell[9] = x39 + x79*(x44 + x48);
return { x73, ((x70)*(x70)) + ((x71)*(x71)) + ((x72)*(x72)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
