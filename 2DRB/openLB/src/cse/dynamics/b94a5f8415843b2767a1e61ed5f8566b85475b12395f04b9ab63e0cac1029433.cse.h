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
struct CSE<CombinedAdvectionDiffusionRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::NoStress, momenta::DefineToEq>, equilibria::FirstOrder, collision::AdvectionDiffusionRLB, AdvectionDiffusionExternalVelocityCollision>, momenta::Tuple<momenta::FixedDensity, momenta::FixedTemperatureMomentum<0, -1>, momenta::NoStress, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x68 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x69 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x70 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x71 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{0.444444444444444};
auto x1 = V{0.111111111111111};
auto x2 = x1*x70;
auto x3 = V{3};
auto x4 = x3*x69 + V{1};
auto x5 = -x1 + x2*x4;
auto x6 = x1*x3;
auto x7 = -V{1}*x69;
auto x8 = x3*x7 + V{1};
auto x9 = -x1 + x2*x8;
auto x10 = -V{1}*x68;
auto x11 = x3*(x10 + x69) + V{1};
auto x12 = V{0.0277777777777778};
auto x72 = -x12;
auto x13 = x12*x70;
auto x14 = x11*x13 + x72;
auto x15 = cell[1] - x14;
auto x16 = x3*(x10 + x7) + V{1};
auto x17 = -x12 + x13*x16;
auto x18 = -V{1}*cell[3] + V{1}*x17;
auto x20 = -V{1}*cell[4] + cell[8] + V{2}*x15 + V{2}*x18 - x5 + V{1}*x9;
auto x21 = x20*x6 + x5;
auto x22 = x3*(x68 + x69) + V{1};
auto x23 = x12*x3;
auto x24 = x10*x3 + V{1};
auto x25 = -x1 + x2*x24;
auto x27 = -V{2}*cell[2] - V{2}*x15 + V{2}*x18 + V{2}*x25;
auto x28 = x13*x22 + x23*(x20 + x27) + x72;
auto x29 = x3*x68 + V{1};
auto x30 = -x1 + x2*x29 + x27*x6;
auto x31 = x3*(x68 + x7) + V{1};
auto x32 = -V{1}*x20;
auto x33 = x13*x31 + x23*(x27 + x32) + x72;
auto x34 = x32*x6 + x9;
auto x35 = -V{1}*x27;
auto x36 = x17 + x23*(x32 + x35);
auto x37 = x25 + x35*x6;
auto x38 = x14 + x23*(x20 + x35);
auto x39 = x0*x70 - x0 + x21 + x28 + x30 + x33 + x34 + x36 + x37 + x38 + V{1};
auto x40 = x12*x39;
auto x41 = x11*x40 + x72;
auto x42 = V{1} - x71;
auto x74 = x23*x42;
auto x75 = x42*x6;
auto x43 = x38 - x41;
auto x44 = -x12 + x16*x40;
auto x45 = -V{1}*x36 + V{1}*x44;
auto x46 = x1*x39;
auto x47 = -x1 + x46*x8;
auto x48 = -x12 + x31*x40;
auto x49 = x33 - x48;
auto x50 = -x12 + x22*x40;
auto x51 = x28 - x50;
auto x73 = x45 + x51;
auto x52 = -x1 + x4*x46;
auto x53 = x21 - V{1}*x34 + x43 + V{1}*x47 - V{1}*x49 - x52 + x73;
auto x54 = -x1 + x24*x46;
auto x55 = -x1 + x29*x46;
auto x56 = x30 - V{1}*x37 - V{1}*x43 + x49 + V{1}*x54 - x55 + x73;
auto x57 = -V{1}*x56;
auto x58 = -V{1}*x53;
cell[0] = x0*(x39 + V{-1});
cell[1] = x41 + x74*(x53 + x57);
cell[2] = x54 + x57*x75;
cell[3] = x44 + x74*(x57 + x58);
cell[4] = x47 + x58*x75;
cell[5] = x48 + x74*(x56 + x58);
cell[6] = x55 + x56*x75;
cell[7] = x50 + x74*(x53 + x56);
cell[8] = x52 + x53*x75;
return { x39, ((x68)*(x68)) + ((x69)*(x69)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
