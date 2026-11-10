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
struct CSE<structuralmechanics::StationaryLinearElasticityMRT<T, descriptors::D2Q8<FIELDS...> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x50 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x51 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x52 = cell.template getFieldComponent<olb::opti::CONTROL_VECTOR>(0);
auto x53 = cell.template getFieldComponent<olb::opti::CONTROL_VECTOR>(1);
auto x78 = parameters.template get<structuralmechanics::OMEGA_SOLID>()[2];
auto x79 = parameters.template get<structuralmechanics::OMEGA_SOLID>()[3];
auto x80 = parameters.template get<structuralmechanics::OMEGA_SOLID>()[4];
auto x81 = parameters.template get<structuralmechanics::OMEGA_SOLID>()[5];
auto x82 = parameters.template get<structuralmechanics::OMEGA_SOLID>()[6];
auto x83 = parameters.template get<structuralmechanics::OMEGA_SOLID>()[7];
auto x56 = cell[0] + cell[2];
auto x58 = cell[4] + cell[5];
auto x0 = V{0.5}*x50*x52;
auto x1 = -V{1}*cell[5];
auto x54 = cell[4] + x1;
auto x2 = -V{1}*cell[6];
auto x55 = cell[7] + x2;
auto x3 = cell[0] - V{1}*cell[2] + x0 + x54 + x55;
auto x4 = x0 + x3;
auto x63 = V{0.5}*x4;
auto x5 = cell[1] + cell[3] + V{2}*cell[4] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7] + x56;
auto x6 = x5*(V{1} - x79);
auto x7 = V{0.25}*x6;
auto x8 = -V{1}*cell[3];
auto x9 = -V{1}*cell[1] + x56 + x8;
auto x10 = x9*(V{1} - x80);
auto x57 = V{0.25}*x10;
auto x11 = x57;
auto x12 = x54;
auto x13 = x12 + x55;
auto x14 = V{0.333333333333333};
auto x15 = x13 - x81*(x13 - x14*x3);
auto x61 = V{0.25}*x15;
auto x64 = V{0.5}*x15;
auto x16 = x58;
auto x17 = cell[6] + cell[7] + x16;
auto x18 = x17*(V{1} - x83);
auto x19 = -V{0.5}*x18;
auto x65 = x19 + x7;
auto x66 = x11 + x65;
auto x20 = V{0.5}*x51*x53;
auto x21 = -V{1}*cell[7];
auto x59 = x2 + x21;
auto x22 = cell[1] + x20 + x58 + x59 + x8;
auto x23 = x20 + x22;
auto x67 = V{0.5}*x23;
auto x24 = -x57;
auto x69 = x24 + x65;
auto x25 = x16 + x59;
auto x26 = x25 + x82*(x14*x22 - x25);
auto x62 = V{0.25}*x26;
auto x68 = V{0.5}*x26;
auto x27 = cell[6] + x12 + x21;
auto x28 = x27*(V{1} - x78);
auto x60 = V{0.25}*x28;
auto x29 = x60;
auto x30 = x61;
auto x31 = x62;
auto x32 = V{0.25}*x18;
auto x70 = x31 + x32;
auto x33 = -x60;
auto x34 = -x61;
auto x35 = -x62;
auto x71 = x32 + x35;
auto x36 = V{0.5}*x3 + V{0.5}*x4;
auto x37 = V{0.5}*x22 + V{0.5}*x23;
auto x38 = V{0.5}*x27 + V{0.5}*x28;
auto x39 = V{0.5}*x5 + V{0.5}*x6;
auto x40 = V{0.5}*x10 + V{0.5}*x9;
auto x41 = -V{1}*x38;
cell[0] = x63 - x64 + x66;
cell[1] = x67 - x68 + x69;
cell[2] = -x63 + x64 + x66;
cell[3] = -x67 + x68 + x69;
cell[4] = x29 + x30 + x70;
cell[5] = x33 + x34 + x70;
cell[6] = x29 + x34 + x71;
cell[7] = x30 + x33 + x71;
cell.template getFieldPointer<olb::structuralmechanics::BARED_MOMENT_VECTOR>()[0] = x36;
cell.template getFieldPointer<olb::structuralmechanics::BARED_MOMENT_VECTOR>()[1] = x37;
cell.template getFieldPointer<olb::structuralmechanics::BARED_MOMENT_VECTOR>()[2] = x38;
cell.template getFieldPointer<olb::structuralmechanics::BARED_MOMENT_VECTOR>()[3] = x39;
cell.template getFieldPointer<olb::structuralmechanics::BARED_MOMENT_VECTOR>()[4] = x40;
cell.template getFieldPointer<olb::structuralmechanics::BARED_MOMENT_VECTOR>()[5] = V{0.5}*x13 + V{0.5}*x15;
cell.template getFieldPointer<olb::structuralmechanics::BARED_MOMENT_VECTOR>()[6] = V{0.5}*x25 + V{0.5}*x26;
cell.template getFieldPointer<olb::structuralmechanics::BARED_MOMENT_VECTOR>()[7] = V{0.5}*x17 + V{0.5}*x18;
cell.template getFieldPointer<olb::structuralmechanics::CAUCHY_STRESS>()[0] = -V{0.5}*x39 - V{0.5}*x40;
cell.template getFieldPointer<olb::structuralmechanics::CAUCHY_STRESS>()[1] = x41;
cell.template getFieldPointer<olb::structuralmechanics::CAUCHY_STRESS>()[2] = x41;
cell.template getFieldPointer<olb::structuralmechanics::CAUCHY_STRESS>()[3] = -V{0.5}*x39 + V{0.5}*x40;
cell.template getFieldPointer<olb::structuralmechanics::DISPLACEMENT>()[0] = x36;
cell.template getFieldPointer<olb::structuralmechanics::DISPLACEMENT>()[1] = x37;
cell.template getFieldPointer<olb::structuralmechanics::MOMENT_VECTOR>()[0] = x4;
cell.template getFieldPointer<olb::structuralmechanics::MOMENT_VECTOR>()[1] = x23;
cell.template getFieldPointer<olb::structuralmechanics::MOMENT_VECTOR>()[2] = x28;
cell.template getFieldPointer<olb::structuralmechanics::MOMENT_VECTOR>()[3] = x6;
cell.template getFieldPointer<olb::structuralmechanics::MOMENT_VECTOR>()[4] = x10;
cell.template getFieldPointer<olb::structuralmechanics::MOMENT_VECTOR>()[5] = x15;
cell.template getFieldPointer<olb::structuralmechanics::MOMENT_VECTOR>()[6] = x26;
cell.template getFieldPointer<olb::structuralmechanics::MOMENT_VECTOR>()[7] = x18;
return { V{1}, ((x36)*(x36)) + ((x37)*(x37)) };
}
};

}

}

// Generation Info: commit=a03b2b7627cc7e10c490423e6846d9c24d39422e
