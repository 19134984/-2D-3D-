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
struct CSE<CombinedAdvectionDiffusionRLBdynamics<T, descriptors::D3Q7<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q7<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::HeatFluxBoundaryDensity<1, 1>, momenta::FixedVelocityMomentumAD, momenta::NoStress, momenta::DefineToEq> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x51 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x52 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x53 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x54 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumAD::VELOCITY>(0);
auto x55 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumAD::VELOCITY>(1);
auto x56 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumAD::VELOCITY>(2);
auto x57 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{0.25};
auto x1 = V{0.125};
auto x60 = -x1;
auto x2 = x52 + V{1};
auto x58 = V{1} / (x2);
auto x3 = x55*x58;
auto x4 = -x3 + x58*(cell[0] + cell[1] + cell[3] + cell[4] + V{2}*cell[5] + cell[6] + V{1});
auto x5 = x1*x4;
auto x6 = V{4};
auto x59 = x53*x6;
auto x62 = x52*x6;
auto x64 = x51*x6;
auto x7 = x59 + V{1};
auto x8 = x1*x6;
auto x10 = x56;
auto x61 = x10*x8;
auto x11 = x5*x7 + x60 + x61;
auto x12 = x62 + V{1};
auto x14 = x55;
auto x63 = x14*x8;
auto x15 = x12*x5 + x60 + x63;
auto x16 = x64 + V{1};
auto x18 = x54;
auto x65 = x18*x8;
auto x19 = x16*x5 + x60 + x65;
auto x20 = V{1} - V{1}*x59;
auto x21 = -x1 + x20*x5 - V{1}*x61;
auto x22 = V{1} - V{1}*x62;
auto x23 = -x1 + x22*x5 - V{1}*x63;
auto x24 = V{1} - V{1}*x64;
auto x25 = -x1 + x24*x5 - V{1}*x65;
auto x26 = x0*x4 - x0 + x11 + x15 + x19 + x21 + x23 + x25 + V{1};
auto x27 = x0*(x26 + V{-1});
auto x28 = x1*x26;
auto x29 = -x1 + x24*x28;
auto x30 = V{1} - x57;
auto x66 = x30*x8;
auto x31 = -x1 + x16*x28;
auto x32 = x19 - V{1}*x25 + V{1}*x29 - x31;
auto x67 = x32*x66;
auto x33 = x29 - V{1}*x67;
auto x34 = -x1 + x22*x28;
auto x35 = -x1 + x12*x28;
auto x36 = x15 - V{1}*x23 + V{1}*x34 - x35;
auto x69 = x36*x66;
auto x37 = -x1 + x20*x28;
auto x38 = -x1 + x28*x7;
auto x39 = x11 - V{1}*x21 + V{1}*x37 - x38;
auto x68 = x39*x66;
auto x40 = x37 - V{1}*x68;
auto x41 = x31 + x67;
auto x42 = x35 + x69;
auto x43 = x38 + x68;
cell[0] = x27;
cell[1] = x33;
cell[2] = x34 - V{1}*x69;
cell[3] = x40;
cell[4] = x41;
cell[5] = x42;
cell[6] = x43;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = -x3 + x58*(x27 + x33 + x40 + x41 + V{2}*x42 + x43 + V{1});
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x54;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x55;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x56;
return { x26, ((x51)*(x51)) + ((x52)*(x52)) + ((x53)*(x53)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
