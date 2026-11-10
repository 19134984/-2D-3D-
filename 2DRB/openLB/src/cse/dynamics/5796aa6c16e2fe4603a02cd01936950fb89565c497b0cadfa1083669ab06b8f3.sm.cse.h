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
struct CSE<CombinedAdvectionDiffusionRLBdynamics<T, descriptors::D2Q5<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q5<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::HeatFluxBoundaryDensity<1, 1>, momenta::FixedVelocityMomentumAD, momenta::NoStress, momenta::DefineToEq> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x38 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x39 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x40 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumAD::VELOCITY>(0);
auto x41 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumAD::VELOCITY>(1);
auto x42 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{0.333333333333333};
auto x1 = V{0.166666666666667};
auto x45 = -x1;
auto x2 = x39 + V{1};
auto x43 = V{1} / (x2);
auto x3 = x41*x43;
auto x4 = -x3 + x43*(cell[0] + cell[1] + cell[3] + V{2}*cell[4] + V{1});
auto x5 = x1*x4;
auto x6 = V{3};
auto x44 = x39*x6;
auto x47 = x38*x6;
auto x7 = x44 + V{1};
auto x8 = x1*x6;
auto x10 = x41;
auto x46 = x10*x8;
auto x11 = x45 + x46 + x5*x7;
auto x12 = x47 + V{1};
auto x14 = x40;
auto x48 = x14*x8;
auto x15 = x12*x5 + x45 + x48;
auto x16 = V{1} - V{1}*x44;
auto x17 = -x1 + x16*x5 - V{1}*x46;
auto x18 = V{1} - V{1}*x47;
auto x19 = -x1 + x18*x5 - V{1}*x48;
auto x20 = x0*x4 - x0 + x11 + x15 + x17 + x19 + V{1};
auto x21 = x0*(x20 + V{-1});
auto x22 = x1*x20;
auto x23 = -x1 + x18*x22;
auto x24 = V{1} - x42;
auto x49 = x24*x8;
auto x25 = -x1 + x12*x22;
auto x26 = x15 - V{1}*x19 + V{1}*x23 - x25;
auto x50 = x26*x49;
auto x27 = x23 - V{1}*x50;
auto x28 = -x1 + x16*x22;
auto x29 = -x1 + x22*x7;
auto x30 = x11 - V{1}*x17 + V{1}*x28 - x29;
auto x51 = x30*x49;
auto x31 = x25 + x50;
auto x32 = x29 + x51;
cell[0] = x21;
cell[1] = x27;
cell[2] = x28 - V{1}*x51;
cell[3] = x31;
cell[4] = x32;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = -x3 + x43*(x21 + x27 + x31 + V{2}*x32 + V{1});
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x40;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x41;
return { x20, ((x38)*(x38)) + ((x39)*(x39)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
