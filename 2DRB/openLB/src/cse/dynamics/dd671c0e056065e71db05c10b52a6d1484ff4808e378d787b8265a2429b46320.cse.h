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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::ParameterFromCell<descriptors::OMEGA, collision::BGK>, forcing::Wagner>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x60 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x61 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x62 = cell.template getFieldComponent<olb::descriptors::OMEGA>(0);
auto x63 = cell.template getFieldComponent<olb::descriptors::SCALAR>(0);
auto x64 = cell[5] + cell[6];
auto x65 = cell[1] + cell[8];
auto x0 = V{1} - x62;
auto x1 = V{3};
auto x2 = -V{1}*cell[3];
auto x67 = cell[7] + x2;
auto x3 = cell[0] + cell[2] + cell[3] + cell[4] + cell[7] + x64 + x65 + V{1};
auto x66 = V{1} / (x3);
auto x4 = x66*(-V{1}*cell[1] - V{1}*cell[2] + x64 + x67);
auto x68 = ((x4)*(x4));
auto x70 = x1*x4;
auto x5 = x66*(-V{1}*cell[4] - V{1}*cell[5] + x65 + x67);
auto x69 = ((x5)*(x5));
auto x71 = x1*x5;
auto x6 = x68 + x69;
auto x7 = V{0.5}*x1*x6;
auto x77 = V{1} - x7;
auto x8 = V{0.444444444444444};
auto x9 = -x60*x70;
auto x10 = -x61*x71;
auto x11 = V{-1.4999985};
auto x12 = V{0.0833333333333333}*x62*x63*x66;
auto x13 = V{1} - V{0.25}*x62;
auto x14 = x13*x60;
auto x15 = x12 + x14*x60;
auto x16 = x11*x15;
auto x76 = x16 + x9;
auto x17 = x13*x61;
auto x18 = x12 + x17*x61;
auto x19 = x11*x18;
auto x75 = x10 + x19;
auto x20 = -V{1}*x4;
auto x21 = x20 + x5;
auto x22 = ((x1)*(x1));
auto x74 = V{1}*x22;
auto x23 = V{0.5}*x22;
auto x24 = V{0.0277777777777778};
auto x72 = x24*x3;
auto x25 = x72;
auto x26 = x21*x22;
auto x27 = -x1*(x4 + V{1});
auto x78 = -x27;
auto x28 = -x1*(x5 + V{-1});
auto x29 = V{3.0000015};
auto x30 = x15*x29;
auto x81 = x30 + x75;
auto x31 = V{-4.5};
auto x32 = x14*x61;
auto x33 = x31*x32;
auto x34 = x17*x60;
auto x35 = x31*x34;
auto x36 = x18*x29;
auto x79 = x30 + x36;
auto x80 = x33 + x35 + x79;
auto x83 = x36 + x76;
auto x37 = V{0.111111111111111};
auto x73 = x3*x37;
auto x38 = x73;
auto x39 = -V{1}*x5;
auto x40 = x20 + x39;
auto x41 = -x40*x74;
auto x42 = -x1*(x5 + V{1});
auto x43 = V{4.5};
auto x44 = x32*x43;
auto x45 = x34*x43;
auto x82 = x44 + x45 + x79;
auto x46 = x39 + x4;
auto x47 = x22*x46;
auto x48 = -x1*(x4 + V{-1});
auto x49 = x4 + x5;
auto x50 = x22*x49;
auto x51 = cell[0]*x0 + x3*x8*(x75 + x76) - x62*x8*(x3*(x7 + V{-1}) + V{1});
auto x52 = cell[1]*x0 - x62*(x24 - x25*(x1*x21 + x23*((x21)*(x21)) + x77)) + x72*(-x60*(V{1}*x26 + x78) + x61*(x26 + x28) + x80);
auto x53 = cell[2]*x0 - x62*(x37 - x38*(x1*x20 + x23*((x20)*(x20)) + x77)) + x73*(-x60*(x20*x74 + x78) + x81);
auto x54 = cell[3]*x0 - x62*(x24 - x25*(x1*x40 + x23*((x40)*(x40)) + x77)) + x72*(x60*(x27 + x41) + x61*(x41 + x42) + x82);
auto x55 = cell[4]*x0 - x62*(x37 - x38*(x1*x39 + x23*((x39)*(x39)) + x77)) + x73*(-x61*(x39*x74 - x42) + x83);
auto x56 = cell[5]*x0 - x62*(x24 - x25*(x1*x46 + x23*((x46)*(x46)) + x77)) + x72*(x60*(x47 + x48) + x61*(x42 - V{1}*x47) + x80);
auto x57 = cell[6]*x0 - x62*(x37 - x38*(x23*x68 + x70 + x77)) + x73*(x60*(x22*x4 + x48) + x81);
auto x58 = cell[7]*x0 - x62*(x24 - x25*(x1*x49 + x23*((x49)*(x49)) + x77)) + x72*(x60*(x48 + x50) + x61*(x28 + x50) + x82);
auto x59 = cell[8]*x0 - x62*(x37 - x38*(x23*x69 + x71 + x77)) + x73*(x61*(x22*x5 + x28) + x83);
cell[0] = x51;
cell[1] = x52;
cell[2] = x53;
cell[3] = x54;
cell[4] = x55;
cell[5] = x56;
cell[6] = x57;
cell[7] = x58;
cell[8] = x59;
return { x3, x6 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
