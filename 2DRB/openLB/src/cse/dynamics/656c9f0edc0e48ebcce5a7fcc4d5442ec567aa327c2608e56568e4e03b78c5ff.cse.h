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
struct CSE<CombinedAdvectionDiffusionRLBdynamics<T, descriptors::D2Q5<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q5<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedTemperatureMomentum<1, -1>, momenta::NoStress, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x34 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x35 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x36 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x37 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{0.333333333333333};
auto x1 = V{0.166666666666667};
auto x2 = x1*x36;
auto x3 = V{3};
auto x38 = x3*x35;
auto x41 = x3*x34;
auto x4 = x38 + V{1};
auto x5 = x1*x3;
auto x6 = V{1} - V{1}*x38;
auto x7 = -x1 + x2*x6;
auto x39 = -x7;
auto x9 = -V{2}*cell[2] - V{2}*x39;
auto x40 = x5*x9;
auto x10 = -x1 + x2*x4 + x40;
auto x11 = x41 + V{1};
auto x12 = -x1 + x11*x2;
auto x13 = V{1} - V{1}*x41;
auto x14 = -x1 + x13*x2;
auto x16 = -V{1}*cell[1] + cell[3] - x12 + V{1}*x14;
auto x42 = x16*x5;
auto x17 = x12 + x42;
auto x18 = -x39 - V{1}*x40;
auto x19 = x14 - V{1}*x42;
auto x20 = x0*x36 - x0 + x10 + x17 + x18 + x19 + V{1};
auto x21 = x1*x20;
auto x22 = -x1 + x13*x21;
auto x23 = V{1} - x37;
auto x43 = x23*x5;
auto x24 = -x1 + x11*x21;
auto x25 = x17 - V{1}*x19 + V{1}*x22 - x24;
auto x44 = x25*x43;
auto x26 = -x1 + x21*x6;
auto x27 = -x1 + x21*x4;
auto x28 = x10 - V{1}*x18 + V{1}*x26 - x27;
auto x45 = x28*x43;
cell[0] = x0*(x20 + V{-1});
cell[1] = x22 - V{1}*x44;
cell[2] = x26 - V{1}*x45;
cell[3] = x24 + x44;
cell[4] = x27 + x45;
return { x20, ((x34)*(x34)) + ((x35)*(x35)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
