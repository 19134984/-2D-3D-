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
struct CSE<CombinedAdvectionDiffusionRLBdynamics<T, descriptors::D2Q5<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q5<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedTemperatureMomentum<0, -1>, momenta::NoStress, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x37 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x38 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x39 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x40 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{0.333333333333333};
auto x1 = V{0.166666666666667};
auto x2 = x1*x39;
auto x3 = V{3};
auto x41 = x3*x38;
auto x44 = x3*x37;
auto x4 = x41 + V{1};
auto x5 = -x1 + x2*x4;
auto x6 = x1*x3;
auto x7 = V{1} - V{1}*x41;
auto x8 = -x1 + x2*x7;
auto x42 = x5 - V{1}*x8;
auto x10 = -V{1}*cell[2] + cell[4] - x42;
auto x43 = x10*x6;
auto x11 = x43 + x5;
auto x12 = x44 + V{1};
auto x13 = V{1} - V{1}*x44;
auto x14 = -x1 + x13*x2;
auto x16 = -V{2}*cell[1] + V{2}*x14;
auto x45 = x16*x6;
auto x17 = -x1 + x12*x2 + x45;
auto x18 = -V{1}*x43 + x8;
auto x19 = x14 - V{1}*x45;
auto x20 = x0*x39 - x0 + x11 + x17 + x18 + x19 + V{1};
auto x21 = x1*x20;
auto x22 = -x1 + x13*x21;
auto x23 = V{1} - x40;
auto x46 = x23*x6;
auto x24 = -x1 + x12*x21;
auto x25 = x17 - V{1}*x19 + V{1}*x22 - x24;
auto x47 = x25*x46;
auto x26 = x22 - V{1}*x47;
auto x27 = -x1 + x21*x7;
auto x28 = -x1 + x21*x4;
auto x29 = x11 - V{1}*x18 + V{1}*x27 - x28;
auto x48 = x29*x46;
auto x30 = x27 - V{1}*x48;
auto x31 = x28 + x48;
cell[0] = x0*(x20 + V{-1});
cell[1] = x26;
cell[2] = x30;
cell[3] = x24 + x47;
cell[4] = x31;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x39;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = V{2}*x14 - V{2}*x26;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = -V{1}*x30 + x31 - x42;
return { x20, ((x37)*(x37)) + ((x38)*(x38)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
