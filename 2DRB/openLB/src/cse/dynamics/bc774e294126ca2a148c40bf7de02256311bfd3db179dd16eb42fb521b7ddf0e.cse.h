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
struct CSE<structuralmechanics::StationaryLinearElasticityMRTBoundary<T, descriptors::D2Q8<FIELDS...> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x60 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x61 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x74 = cell.template getFieldComponent<olb::structuralmechanics::CURRENT_CELL>(0);
auto x75 = cell.template getFieldComponent<olb::structuralmechanics::CURRENT_CELL>(1);
auto x76 = cell.template getFieldComponent<olb::structuralmechanics::CURRENT_CELL>(2);
auto x77 = cell.template getFieldComponent<olb::structuralmechanics::CURRENT_CELL>(3);
auto x78 = cell.template getFieldComponent<olb::structuralmechanics::CURRENT_CELL>(4);
auto x79 = cell.template getFieldComponent<olb::structuralmechanics::CURRENT_CELL>(5);
auto x80 = cell.template getFieldComponent<olb::structuralmechanics::CURRENT_CELL>(6);
auto x81 = cell.template getFieldComponent<olb::structuralmechanics::CURRENT_CELL>(7);
auto x102 = parameters.template get<structuralmechanics::OMEGA_SOLID>()[2];
auto x103 = parameters.template get<structuralmechanics::OMEGA_SOLID>()[3];
auto x104 = parameters.template get<structuralmechanics::OMEGA_SOLID>()[4];
auto x105 = parameters.template get<structuralmechanics::OMEGA_SOLID>()[5];
auto x106 = parameters.template get<structuralmechanics::OMEGA_SOLID>()[6];
auto x107 = parameters.template get<structuralmechanics::OMEGA_SOLID>()[7];
auto x64 = cell[0] + cell[2];
auto x0 = V{0.5}*x60;
auto x1 = -V{1}*cell[5];
auto x62 = cell[4] + x1;
auto x2 = -V{1}*cell[6];
auto x63 = cell[7] + x2;
auto x3 = cell[0] - V{1}*cell[2] + x0 + x62 + x63;
auto x4 = x0 + x3;
auto x66 = V{0.5}*x4;
auto x5 = cell[1] + cell[3] + V{2}*cell[4] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7] + x64;
auto x6 = x5*(V{1} - x103);
auto x7 = V{0.333333333333333};
auto x8 = -V{0.5}*x7/((V{1} - 1/x103)*(x7 + V{1}));
auto x9 = x6*(V{0.5}*x8 + V{0.25});
auto x10 = -V{1}*cell[3];
auto x11 = -V{1}*cell[1] + x10 + x64;
auto x12 = x11*(V{1} - x104);
auto x65 = V{0.25}*x12;
auto x13 = x65;
auto x14 = x62;
auto x15 = x14 + x63;
auto x16 = -x105*(x15 - x3*x7) + x15;
auto x67 = V{0.5}*x16;
auto x83 = V{0.25}*x16;
auto x17 = V{2}*x8 + V{1};
auto x18 = cell[0]*x8 + cell[1]*x8 + cell[2]*x8 + cell[3]*x8 + cell[4]*x17 + cell[5]*x17 + cell[6]*x17 + cell[7]*x17;
auto x19 = x18*(V{1} - x107);
auto x20 = -V{0.5}*x19;
auto x68 = x20 + x9;
auto x69 = x13 + x68;
auto x21 = x66 - x67 + x69;
auto x22 = V{0.5}*x61;
auto x23 = -V{1}*cell[7];
auto x70 = cell[4] + cell[5] + x2 + x23;
auto x24 = cell[1] + x10 + x22 + x70;
auto x25 = x22 + x24;
auto x71 = V{0.5}*x25;
auto x26 = -x65;
auto x73 = x26 + x68;
auto x27 = x70;
auto x28 = x106*(x24*x7 - x27) + x27;
auto x72 = V{0.5}*x28;
auto x84 = V{0.25}*x28;
auto x29 = x71 - x72 + x73;
auto x30 = -x66 + x67 + x69;
auto x31 = -x71 + x72 + x73;
auto x32 = cell[6] + x14 + x23;
auto x33 = x32*(V{1} - x102);
auto x82 = V{0.25}*x33;
auto x34 = -V{0.25}*x6*x8;
auto x35 = x34 + x82;
auto x36 = x83;
auto x37 = x84;
auto x38 = V{0.25}*x19;
auto x85 = x37 + x38;
auto x39 = x35 + x36 + x85;
auto x40 = x34 - x82;
auto x41 = -x83;
auto x42 = x40 + x41 + x85;
auto x43 = -x84;
auto x86 = x38 + x43;
auto x44 = x35 + x41 + x86;
auto x45 = x36 + x40 + x86;
auto x46 = V{0.5}*x3 + V{0.5}*x4;
auto x47 = V{0.5}*x24 + V{0.5}*x25;
auto x48 = V{0.5}*x32 + V{0.5}*x33;
auto x49 = V{0.5}*x5 + V{0.5}*x6;
auto x50 = V{0.5}*x11 + V{0.5}*x12;
auto x51 = -V{1}*x48;
cell[0] = x21;
cell[1] = x29;
cell[2] = x30;
cell[3] = x31;
cell[4] = x39;
cell[5] = x42;
cell[6] = x44;
cell[7] = x45;
cell.template getFieldPointer<olb::structuralmechanics::BARED_MOMENT_VECTOR>()[0] = x46;
cell.template getFieldPointer<olb::structuralmechanics::BARED_MOMENT_VECTOR>()[1] = x47;
cell.template getFieldPointer<olb::structuralmechanics::BARED_MOMENT_VECTOR>()[2] = x48;
cell.template getFieldPointer<olb::structuralmechanics::BARED_MOMENT_VECTOR>()[3] = x49;
cell.template getFieldPointer<olb::structuralmechanics::BARED_MOMENT_VECTOR>()[4] = x50;
cell.template getFieldPointer<olb::structuralmechanics::BARED_MOMENT_VECTOR>()[5] = V{0.5}*x15 + V{0.5}*x16;
cell.template getFieldPointer<olb::structuralmechanics::BARED_MOMENT_VECTOR>()[6] = V{0.5}*x27 + V{0.5}*x28;
cell.template getFieldPointer<olb::structuralmechanics::BARED_MOMENT_VECTOR>()[7] = V{0.5}*x18 + V{0.5}*x19;
cell.template getFieldPointer<olb::structuralmechanics::CAUCHY_STRESS>()[0] = -V{0.5}*x49 - V{0.5}*x50;
cell.template getFieldPointer<olb::structuralmechanics::CAUCHY_STRESS>()[1] = x51;
cell.template getFieldPointer<olb::structuralmechanics::CAUCHY_STRESS>()[2] = x51;
cell.template getFieldPointer<olb::structuralmechanics::CAUCHY_STRESS>()[3] = -V{0.5}*x49 + V{0.5}*x50;
cell.template getFieldPointer<olb::structuralmechanics::CURRENT_CELL>()[0] = x21;
cell.template getFieldPointer<olb::structuralmechanics::CURRENT_CELL>()[1] = x29;
cell.template getFieldPointer<olb::structuralmechanics::CURRENT_CELL>()[2] = x30;
cell.template getFieldPointer<olb::structuralmechanics::CURRENT_CELL>()[3] = x31;
cell.template getFieldPointer<olb::structuralmechanics::CURRENT_CELL>()[4] = x39;
cell.template getFieldPointer<olb::structuralmechanics::CURRENT_CELL>()[5] = x42;
cell.template getFieldPointer<olb::structuralmechanics::CURRENT_CELL>()[6] = x44;
cell.template getFieldPointer<olb::structuralmechanics::CURRENT_CELL>()[7] = x45;
cell.template getFieldPointer<olb::structuralmechanics::DISPLACEMENT>()[0] = x46;
cell.template getFieldPointer<olb::structuralmechanics::DISPLACEMENT>()[1] = x47;
cell.template getFieldPointer<olb::structuralmechanics::MOMENT_VECTOR>()[0] = x4;
cell.template getFieldPointer<olb::structuralmechanics::MOMENT_VECTOR>()[1] = x25;
cell.template getFieldPointer<olb::structuralmechanics::MOMENT_VECTOR>()[2] = x33;
cell.template getFieldPointer<olb::structuralmechanics::MOMENT_VECTOR>()[3] = x6;
cell.template getFieldPointer<olb::structuralmechanics::MOMENT_VECTOR>()[4] = x12;
cell.template getFieldPointer<olb::structuralmechanics::MOMENT_VECTOR>()[5] = x16;
cell.template getFieldPointer<olb::structuralmechanics::MOMENT_VECTOR>()[6] = x28;
cell.template getFieldPointer<olb::structuralmechanics::MOMENT_VECTOR>()[7] = x19;
cell.template getFieldPointer<olb::structuralmechanics::PREVIOUS_CELL>()[0] = x74;
cell.template getFieldPointer<olb::structuralmechanics::PREVIOUS_CELL>()[1] = x75;
cell.template getFieldPointer<olb::structuralmechanics::PREVIOUS_CELL>()[2] = x76;
cell.template getFieldPointer<olb::structuralmechanics::PREVIOUS_CELL>()[3] = x77;
cell.template getFieldPointer<olb::structuralmechanics::PREVIOUS_CELL>()[4] = x78;
cell.template getFieldPointer<olb::structuralmechanics::PREVIOUS_CELL>()[5] = x79;
cell.template getFieldPointer<olb::structuralmechanics::PREVIOUS_CELL>()[6] = x80;
cell.template getFieldPointer<olb::structuralmechanics::PREVIOUS_CELL>()[7] = x81;
return { V{1}, ((x46)*(x46)) + ((x47)*(x47)) };
}
};

}

}

// Generation Info: commit=a03b2b7627cc7e10c490423e6846d9c24d39422e
