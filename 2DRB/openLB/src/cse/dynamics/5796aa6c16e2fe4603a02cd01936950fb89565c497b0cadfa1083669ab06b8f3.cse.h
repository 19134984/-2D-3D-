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
struct CSE<CombinedAdvectionDiffusionRLBdynamics<T, descriptors::D2Q5<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q5<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::HeatFluxBoundaryDensity<1, 1>, momenta::FixedVelocityMomentumAD, momenta::NoStress, momenta::DefineToEq> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x33 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x34 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x35 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumAD::VELOCITY>(0);
auto x36 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumAD::VELOCITY>(1);
auto x37 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{0.333333333333333};
auto x1 = V{0.166666666666667};
auto x39 = -x1;
auto x2 = x34 + V{1};
auto x3 = (cell[0] + cell[1] + cell[3] + V{2}*cell[4] - x36 + V{1})/x2;
auto x4 = x1*x3;
auto x5 = V{3};
auto x38 = x34*x5;
auto x41 = x33*x5;
auto x6 = x38 + V{1};
auto x7 = x1*x5;
auto x9 = x36;
auto x40 = x7*x9;
auto x10 = x39 + x4*x6 + x40;
auto x11 = x41 + V{1};
auto x13 = x35;
auto x42 = x13*x7;
auto x14 = x11*x4 + x39 + x42;
auto x15 = V{1} - V{1}*x38;
auto x16 = -x1 + x15*x4 - V{1}*x40;
auto x17 = V{1} - V{1}*x41;
auto x18 = -x1 + x17*x4 - V{1}*x42;
auto x19 = x0*x3 - x0 + x10 + x14 + x16 + x18 + V{1};
auto x20 = x1*x19;
auto x21 = -x1 + x17*x20;
auto x22 = V{1} - x37;
auto x43 = x22*x7;
auto x23 = -x1 + x11*x20;
auto x24 = x14 - V{1}*x18 + V{1}*x21 - x23;
auto x44 = x24*x43;
auto x25 = -x1 + x15*x20;
auto x26 = -x1 + x20*x6;
auto x27 = x10 - V{1}*x16 + V{1}*x25 - x26;
auto x45 = x27*x43;
cell[0] = x0*(x19 + V{-1});
cell[1] = x21 - V{1}*x44;
cell[2] = x25 - V{1}*x45;
cell[3] = x23 + x44;
cell[4] = x26 + x45;
return { x19, ((x33)*(x33)) + ((x34)*(x34)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
