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
struct CSE<Dual<structuralmechanics::StationaryLinearElasticityMRT<T, descriptors::D2Q8<FIELDS...> > >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x63 = cell.template getFieldComponent<olb::opti::DJDF>(0);
auto x64 = cell.template getFieldComponent<olb::opti::DJDF>(1);
auto x65 = cell.template getFieldComponent<olb::opti::DJDF>(2);
auto x66 = cell.template getFieldComponent<olb::opti::DJDF>(3);
auto x67 = cell.template getFieldComponent<olb::opti::DJDF>(4);
auto x68 = cell.template getFieldComponent<olb::opti::DJDF>(5);
auto x69 = cell.template getFieldComponent<olb::opti::DJDF>(6);
auto x70 = cell.template getFieldComponent<olb::opti::DJDF>(7);
auto x103 = parameters.template get<structuralmechanics::OMEGA_SOLID>()[2];
auto x104 = parameters.template get<structuralmechanics::OMEGA_SOLID>()[3];
auto x105 = parameters.template get<structuralmechanics::OMEGA_SOLID>()[4];
auto x106 = parameters.template get<structuralmechanics::OMEGA_SOLID>()[5];
auto x107 = parameters.template get<structuralmechanics::OMEGA_SOLID>()[6];
auto x108 = parameters.template get<structuralmechanics::OMEGA_SOLID>()[7];
auto x59 = V{0.5}*cell[2];
auto x60 = V{0.5}*cell[0];
auto x61 = V{0.25}*cell[7];
auto x62 = V{0.25}*cell[6];
auto x71 = V{0.25}*cell[5];
auto x74 = V{0.25}*cell[1] + V{0.25}*cell[3];
auto x75 = V{0.5}*cell[3];
auto x76 = V{0.5}*cell[1];
auto x0 = -x59;
auto x1 = x0 + x60;
auto x3 = V{0.333333333333333};
auto x4 = x61;
auto x5 = -x62;
auto x6 = -x71;
auto x7 = V{0.25}*cell[4];
auto x72 = x6 + x7;
auto x8 = -x60;
auto x73 = x4 + x8;
auto x9 = x5 + x59 + x72 + x73;
auto x10 = x106*x9;
auto x11 = x1 + x10*x3;
auto x12 = V{0.25}*cell[2];
auto x13 = V{0.25}*cell[0];
auto x14 = x12 + x13 - x74;
auto x15 = x14*(V{1} - V{1}*x105);
auto x16 = x12 + x13 + x74;
auto x17 = x16*(V{1} - V{1}*x104);
auto x18 = x15 + x17;
auto x19 = x11 + x18 - x63;
auto x20 = -x75;
auto x21 = x20 + x76;
auto x23 = -x61;
auto x24 = x71;
auto x25 = -x76;
auto x77 = x24 + x25 + x7;
auto x26 = x23 + x5 + x75 + x77;
auto x27 = x107*x26;
auto x28 = x21 + x27*x3;
auto x29 = -V{1}*x15 + x17;
auto x30 = x28 + x29 - x64;
auto x31 = -V{1}*x11;
auto x32 = x18 + x31 - x65;
auto x33 = -V{1}*x28;
auto x34 = x29 + x33 - x66;
auto x35 = x62;
auto x36 = x23 + x35 + x72;
auto x37 = x36*(V{1} - V{1}*x103);
auto x38 = -V{1}*x10 + x9;
auto x80 = x11 + x38;
auto x39 = x0 + x20 + x35 + x73 + x77;
auto x40 = x39*(V{1} - V{1}*x108);
auto x41 = x26 - V{1}*x27;
auto x42 = x40 + x41;
auto x79 = x28 + x42;
auto x43 = V{2}*x17;
auto x78 = x37 + x43;
auto x44 = -x67 + x78 + x79 + x80;
auto x45 = -V{1}*x38;
auto x81 = x31 + x45;
auto x46 = -V{1}*x37;
auto x82 = x43 + x46;
auto x47 = -x68 + x79 + x81 + x82;
auto x48 = x40 - V{1}*x41;
auto x83 = x33 + x48;
auto x49 = -x69 + x78 + x81 + x83;
auto x50 = -x70 + x80 + x82 + x83;
cell[0] = x19;
cell[1] = x30;
cell[2] = x32;
cell[3] = x34;
cell[4] = x44;
cell[5] = x47;
cell[6] = x49;
cell[7] = x50;
return { x19 + V{1}, ((x19)*(x19)) + ((x30)*(x30)) + ((x32)*(x32)) + ((x34)*(x34)) + ((x44)*(x44)) + ((x47)*(x47)) + ((x49)*(x49)) + ((x50)*(x50)) };
}
};

}

}

// Generation Info: commit=a03b2b7627cc7e10c490423e6846d9c24d39422e
