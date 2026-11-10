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
struct CSE<structuralmechanics::ElastodynamicsVectorLBM<T, descriptors::VECTOR<descriptors::D2Q4<FIELDS...>, 5> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x47 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x48 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x53 = cell.template getFieldComponent<olb::structuralmechanics::DISPLACEMENT>(0);
auto x54 = cell.template getFieldComponent<olb::structuralmechanics::DISPLACEMENT>(1);
auto x55 = cell.template getFieldComponent<olb::structuralmechanics::VELOCITY>(0);
auto x56 = cell.template getFieldComponent<olb::structuralmechanics::VELOCITY>(1);
auto x58 = parameters.template get<descriptors::OMEGA>();
auto x59 = parameters.template get<structuralmechanics::C_K>();
auto x60 = parameters.template get<structuralmechanics::C_MU>();
auto x61 = parameters.template get<structuralmechanics::LATTICE_SPEED>();
auto x0 = V{0.25};
auto x62 = x0*x58;
auto x1 = -x0*(x58 + V{-2});
auto x2 = V{2};
auto x49 = x1/x2;
auto x3 = x47*x49;
auto x4 = V{1} - x58;
auto x5 = cell[0] + cell[10] + cell[15] + cell[5] + V{0.5}*x47;
auto x6 = cell[13] + cell[18] + cell[3] + cell[8];
auto x7 = x6*x60;
auto x8 = cell[12] + cell[17] + cell[2] + cell[7];
auto x63 = -x8;
auto x9 = x59*x8;
auto x10 = x7 + x9;
auto x57 = x10*x2;
auto x11 = x48*x49;
auto x12 = cell[11] + cell[16] + cell[1] + cell[6] + V{0.5}*x48;
auto x52 = x12*x60;
auto x13 = cell[14] + cell[19] + cell[4] + cell[9];
auto x14 = x13*x60;
auto x50 = x14*x2;
auto x15 = x50;
auto x16 = 0;
auto x17 = x5*x59;
auto x64 = x17*x2;
auto x18 = x5*x60;
auto x51 = x18*x2;
auto x19 = x51;
auto x20 = x52;
auto x65 = x2*x20;
auto x21 = -x7 + x9;
auto x66 = x2*x21;
auto x22 = x12*x59;
auto x67 = x2*x22;
auto x23 = -V{1}*x52;
auto x68 = x2*x23;
auto x24 = -V{1}*x50;
auto x25 = -V{1}*x51;
auto x26 = x14*x61;
auto x27 = cell[0]*x4 + x3 + x62*(x5 + x57);
auto x28 = cell[10]*x4 + x3 - x62*(-x5 + V{1}*x57);
auto x29 = cell[11]*x4 + x11 + x62*(x12 + x24);
auto x30 = cell[12]*x4 + x16 - x62*(x63 + V{1}*x64);
auto x31 = cell[13]*x4 + x16 + x62*(x25 + x6);
auto x32 = cell[14]*x4 + x16 + x62*(x13 - V{1}*x65);
auto x33 = cell[15]*x4 + x3 + x62*(x24 + x5);
auto x34 = cell[16]*x4 + x11 + x62*(x12 - V{1}*x66);
auto x35 = cell[17]*x4 + x16 - x62*(x63 + V{1}*x67);
auto x36 = cell[18]*x4 + x16 - x62*(-x6 + V{1}*x68);
auto x37 = cell[19]*x4 + x16 + x62*(x13 + x25);
auto x38 = cell[1]*x4 + x11 + x62*(x12 + x15);
auto x39 = cell[2]*x4 + x16 + x62*(x64 + x8);
auto x40 = cell[3]*x4 + x16 + x62*(x19 + x6);
auto x41 = cell[4]*x4 + x16 + x62*(x13 + x65);
auto x42 = cell[5]*x4 + x3 + x62*(x15 + x5);
auto x43 = cell[6]*x4 + x11 + x62*(x12 + x66);
auto x44 = cell[7]*x4 + x16 + x62*(x67 + x8);
auto x45 = cell[8]*x4 + x16 + x62*(x6 + x68);
auto x46 = cell[9]*x4 + x16 + x62*(x13 + x19);
cell.template getFieldPointer<olb::structuralmechanics::CAUCHY_STRESS>()[0] = x10*x61;
cell.template getFieldPointer<olb::structuralmechanics::CAUCHY_STRESS>()[1] = x26;
cell.template getFieldPointer<olb::structuralmechanics::CAUCHY_STRESS>()[2] = x26;
cell.template getFieldPointer<olb::structuralmechanics::CAUCHY_STRESS>()[3] = x21*x61;
cell.template getFieldPointer<olb::structuralmechanics::DISPLACEMENT>()[0] = V{0.5}*x5 + x53 + V{0.5}*x55;
cell.template getFieldPointer<olb::structuralmechanics::DISPLACEMENT>()[1] = V{0.5}*x12 + x54 + V{0.5}*x56;
cell.template getFieldPointer<olb::structuralmechanics::VELOCITY>()[0] = x5;
cell.template getFieldPointer<olb::structuralmechanics::VELOCITY>()[1] = x12;
cell[0] = x27;
cell[10] = x28;
cell[11] = x29;
cell[12] = x30;
cell[13] = x31;
cell[14] = x32;
cell[15] = x33;
cell[16] = x34;
cell[17] = x35;
cell[18] = x36;
cell[19] = x37;
cell[1] = x38;
cell[2] = x39;
cell[3] = x40;
cell[4] = x41;
cell[5] = x42;
cell[6] = x43;
cell[7] = x44;
cell[8] = x45;
cell[9] = x46;
return { V{1}, ((x12)*(x12)) + ((x5)*(x5)) };
}
};

}

}

// Generation Info: commit=a03b2b7627cc7e10c490423e6846d9c24d39422e
