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
struct CSE<CombinedAdvectionDiffusionRLBdynamics<T, descriptors::D3Q7<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q7<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedTemperatureMomentum<1, -1>, momenta::NoStress, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x52 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x53 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x54 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x55 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x56 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{0.25};
auto x1 = V{0.125};
auto x2 = x1*x55;
auto x3 = V{4};
auto x57 = x3*x54;
auto x60 = x3*x53;
auto x62 = x3*x52;
auto x4 = x57 + V{1};
auto x5 = -x1 + x2*x4;
auto x6 = x1*x3;
auto x7 = V{1} - V{1}*x57;
auto x8 = -x1 + x2*x7;
auto x58 = x5 - V{1}*x8;
auto x10 = -V{1}*cell[3] + cell[6] - x58;
auto x59 = x10*x6;
auto x11 = x5 + x59;
auto x12 = x60 + V{1};
auto x13 = V{1} - V{1}*x60;
auto x14 = -x1 + x13*x2;
auto x16 = -V{2}*cell[2] + V{2}*x14;
auto x61 = x16*x6;
auto x17 = -x1 + x12*x2 + x61;
auto x18 = x62 + V{1};
auto x19 = -x1 + x18*x2;
auto x20 = V{1} - V{1}*x62;
auto x21 = -x1 + x2*x20;
auto x63 = x19 - V{1}*x21;
auto x23 = -V{1}*cell[1] + cell[4] - x63;
auto x64 = x23*x6;
auto x24 = x19 + x64;
auto x25 = -V{1}*x59 + x8;
auto x26 = x14 - V{1}*x61;
auto x27 = x21 - V{1}*x64;
auto x28 = x0*x55 - x0 + x11 + x17 + x24 + x25 + x26 + x27 + V{1};
auto x29 = x1*x28;
auto x30 = -x1 + x20*x29;
auto x31 = V{1} - x56;
auto x65 = x31*x6;
auto x32 = -x1 + x18*x29;
auto x33 = x24 - V{1}*x27 + V{1}*x30 - x32;
auto x66 = x33*x65;
auto x34 = x30 - V{1}*x66;
auto x35 = -x1 + x13*x29;
auto x36 = -x1 + x12*x29;
auto x37 = x17 - V{1}*x26 + V{1}*x35 - x36;
auto x67 = x37*x65;
auto x38 = x35 - V{1}*x67;
auto x39 = -x1 + x29*x7;
auto x40 = -x1 + x29*x4;
auto x41 = x11 - V{1}*x25 + V{1}*x39 - x40;
auto x68 = x41*x65;
auto x42 = x39 - V{1}*x68;
auto x43 = x32 + x66;
auto x44 = x40 + x68;
cell[0] = x0*(x28 + V{-1});
cell[1] = x34;
cell[2] = x38;
cell[3] = x42;
cell[4] = x43;
cell[5] = x36 + x67;
cell[6] = x44;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x55;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = -V{1}*x34 + x43 - x63;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = V{2}*x14 - V{2}*x38;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = -V{1}*x42 + x44 - x58;
return { x28, ((x52)*(x52)) + ((x53)*(x53)) + ((x54)*(x54)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
