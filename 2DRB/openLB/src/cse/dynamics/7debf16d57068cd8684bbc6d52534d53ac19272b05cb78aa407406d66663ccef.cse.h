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
struct CSE<CombinedAdvectionDiffusionRLBdynamics<T, descriptors::D3Q7<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q7<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::HeatFluxBoundaryDensity<1, 1>, momenta::FixedVelocityMomentumAD, momenta::NoStress, momenta::DefineToEq> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x44 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x45 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x46 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x47 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumAD::VELOCITY>(0);
auto x48 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumAD::VELOCITY>(1);
auto x49 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumAD::VELOCITY>(2);
auto x50 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{0.25};
auto x1 = V{0.125};
auto x52 = -x1;
auto x2 = x45 + V{1};
auto x3 = (cell[0] + cell[1] + cell[3] + cell[4] + V{2}*cell[5] + cell[6] - x48 + V{1})/x2;
auto x4 = x1*x3;
auto x5 = V{4};
auto x51 = x46*x5;
auto x54 = x45*x5;
auto x56 = x44*x5;
auto x6 = x51 + V{1};
auto x7 = x1*x5;
auto x9 = x49;
auto x53 = x7*x9;
auto x10 = x4*x6 + x52 + x53;
auto x11 = x54 + V{1};
auto x13 = x48;
auto x55 = x13*x7;
auto x14 = x11*x4 + x52 + x55;
auto x15 = x56 + V{1};
auto x17 = x47;
auto x57 = x17*x7;
auto x18 = x15*x4 + x52 + x57;
auto x19 = V{1} - V{1}*x51;
auto x20 = -x1 + x19*x4 - V{1}*x53;
auto x21 = V{1} - V{1}*x54;
auto x22 = -x1 + x21*x4 - V{1}*x55;
auto x23 = V{1} - V{1}*x56;
auto x24 = -x1 + x23*x4 - V{1}*x57;
auto x25 = x0*x3 - x0 + x10 + x14 + x18 + x20 + x22 + x24 + V{1};
auto x26 = x1*x25;
auto x27 = -x1 + x23*x26;
auto x28 = V{1} - x50;
auto x58 = x28*x7;
auto x29 = -x1 + x15*x26;
auto x30 = x18 - V{1}*x24 + V{1}*x27 - x29;
auto x59 = x30*x58;
auto x31 = -x1 + x21*x26;
auto x32 = -x1 + x11*x26;
auto x33 = x14 - V{1}*x22 + V{1}*x31 - x32;
auto x60 = x33*x58;
auto x34 = -x1 + x19*x26;
auto x35 = -x1 + x26*x6;
auto x36 = x10 - V{1}*x20 + V{1}*x34 - x35;
auto x61 = x36*x58;
cell[0] = x0*(x25 + V{-1});
cell[1] = x27 - V{1}*x59;
cell[2] = x31 - V{1}*x60;
cell[3] = x34 - V{1}*x61;
cell[4] = x29 + x59;
cell[5] = x32 + x60;
cell[6] = x35 + x61;
return { x25, ((x44)*(x44)) + ((x45)*(x45)) + ((x46)*(x46)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
