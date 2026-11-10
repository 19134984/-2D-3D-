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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::ParameterFromCell<descriptors::OMEGA, collision::BGK>, forcing::Wagner>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x71 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x72 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x73 = cell.template getFieldComponent<olb::descriptors::OMEGA>(0);
auto x74 = cell.template getFieldComponent<olb::descriptors::SCALAR>(0);
auto x75 = cell[5] + cell[6];
auto x76 = cell[1] + cell[8];
auto x0 = V{1} - x73;
auto x1 = V{3};
auto x2 = -V{1}*cell[3];
auto x78 = cell[7] + x2;
auto x3 = cell[0] + cell[2] + cell[3] + cell[4] + cell[7] + x75 + x76 + V{1};
auto x77 = V{1} / (x3);
auto x4 = x77*(-V{1}*cell[1] - V{1}*cell[2] + x75 + x78);
auto x79 = ((x4)*(x4));
auto x81 = x1*x4;
auto x5 = x77*(-V{1}*cell[4] - V{1}*cell[5] + x76 + x78);
auto x80 = ((x5)*(x5));
auto x82 = x1*x5;
auto x6 = x79 + x80;
auto x7 = V{0.5}*x1*x6;
auto x86 = V{1} - x7;
auto x8 = V{0.444444444444444};
auto x9 = -x71*x81;
auto x10 = -x72*x82;
auto x11 = V{-1.4999985};
auto x12 = V{0.0833333333333333}*x73*x74*x77;
auto x13 = V{1} - V{0.25}*x73;
auto x14 = x13*x71;
auto x15 = x12 + x14*x71;
auto x16 = x11*x15;
auto x84 = x16 + x9;
auto x17 = x13*x72;
auto x18 = x12 + x17*x72;
auto x19 = x11*x18;
auto x83 = x10 + x19;
auto x20 = cell[0]*x0 + x3*x8*(x83 + x84) - x73*x8*(x3*(x7 + V{-1}) + V{1});
auto x21 = -V{1}*x4;
auto x22 = x21 + x5;
auto x23 = ((x1)*(x1));
auto x91 = V{1}*x23;
auto x24 = V{0.5}*x23;
auto x25 = V{0.0277777777777778};
auto x85 = x25*x3;
auto x26 = x85;
auto x27 = x22*x23;
auto x28 = -x1*(x4 + V{1});
auto x87 = -x28;
auto x29 = -x1*(x5 + V{-1});
auto x30 = V{3.0000015};
auto x31 = x15*x30;
auto x92 = x31 + x83;
auto x32 = V{-4.5};
auto x33 = x14*x72;
auto x34 = x32*x33;
auto x35 = x17*x71;
auto x36 = x32*x35;
auto x37 = x18*x30;
auto x88 = x31 + x37;
auto x89 = x34 + x36 + x88;
auto x94 = x37 + x84;
auto x38 = cell[1]*x0 - x73*(x25 - x26*(x1*x22 + x24*((x22)*(x22)) + x86)) + x85*(-x71*(V{1}*x27 + x87) + x72*(x27 + x29) + x89);
auto x39 = V{0.111111111111111};
auto x90 = x3*x39;
auto x40 = x90;
auto x41 = cell[2]*x0 - x73*(x39 - x40*(x1*x21 + x24*((x21)*(x21)) + x86)) + x90*(-x71*(x21*x91 + x87) + x92);
auto x42 = -V{1}*x5;
auto x43 = x21 + x42;
auto x44 = -x43*x91;
auto x45 = -x1*(x5 + V{1});
auto x46 = V{4.5};
auto x47 = x33*x46;
auto x48 = x35*x46;
auto x93 = x47 + x48 + x88;
auto x49 = cell[3]*x0 - x73*(x25 - x26*(x1*x43 + x24*((x43)*(x43)) + x86)) + x85*(x71*(x28 + x44) + x72*(x44 + x45) + x93);
auto x50 = cell[4]*x0 - x73*(x39 - x40*(x1*x42 + x24*((x42)*(x42)) + x86)) + x90*(-x72*(x42*x91 - x45) + x94);
auto x51 = x4 + x42;
auto x52 = x23*x51;
auto x53 = -x1*(x4 + V{-1});
auto x54 = cell[5]*x0 - x73*(x25 - x26*(x1*x51 + x24*((x51)*(x51)) + x86)) + x85*(x71*(x52 + x53) + x72*(x45 - V{1}*x52) + x89);
auto x55 = cell[6]*x0 - x73*(x39 - x40*(x24*x79 + x81 + x86)) + x90*(x71*(x23*x4 + x53) + x92);
auto x95 = x54 + x55;
auto x56 = x4 + x5;
auto x57 = x23*x56;
auto x58 = cell[7]*x0 - x73*(x25 - x26*(x1*x56 + x24*((x56)*(x56)) + x86)) + x85*(x71*(x53 + x57) + x72*(x29 + x57) + x93);
auto x59 = cell[8]*x0 - x73*(x39 - x40*(x24*x80 + x82 + x86)) + x90*(x72*(x23*x5 + x29) + x94);
auto x96 = x38 + x59;
auto x60 = x20 + x41 + x49 + x50 + x58 + x95 + x96 + V{1};
auto x97 = V{1} / (x60);
auto x61 = -V{1}*x49;
auto x98 = x58 + x61;
cell[0] = x20;
cell[1] = x38;
cell[2] = x41;
cell[3] = x49;
cell[4] = x50;
cell[5] = x54;
cell[6] = x55;
cell[7] = x58;
cell[8] = x59;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x60;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x97*(-V{1}*x38 - V{1}*x41 + x95 + x98);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x97*(-V{1}*x50 - V{1}*x54 + x96 + x98);
return { x3, x6 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
