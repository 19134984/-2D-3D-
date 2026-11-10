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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<0, 1>, momenta::RegularizedBoundaryStress<0, 1>, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x88 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x90 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1);
auto x91 = parameters.template get<descriptors::OMEGA>();
auto x89 = ((x90)*(x90));
auto x93 = V{1} / (x88);
auto x94 = V{2}*cell[6];
auto x95 = cell[4] + cell[8];
auto x0 = V{0.444444444444444};
auto x1 = x0*x88;
auto x2 = V{3};
auto x92 = ((x2)*(x2));
auto x3 = V{0.5}*x2;
auto x4 = x89;
auto x5 = V{0.111111111111111};
auto x100 = -x5;
auto x6 = x5*x88;
auto x7 = V{0.5}*x92;
auto x8 = x2*x90 + x7*x89 + V{1};
auto x97 = -x8;
auto x9 = x93*(cell[0] + V{2}*cell[5] + V{2}*cell[7] + x94 + x95 + V{1}) + V{-1};
auto x96 = ((x9)*(x9));
auto x10 = x3*(x4 + x96);
auto x99 = V{1} - x10;
auto x11 = -x5 - x6*(x10 + x97);
auto x12 = V{0.5}*x5*x92;
auto x13 = V{1}/x2;
auto x14 = -x13;
auto x15 = -V{1}*x90;
auto x16 = x15 + x9;
auto x17 = V{0.0277777777777778};
auto x98 = -x17;
auto x18 = x17*x88;
auto x19 = x18*(x16*x2 + x7*((x16)*(x16)) + x99) + x98;
auto x20 = cell[5] - x19;
auto x21 = x100 + x6*(x2*x9 + x7*x96 + x99);
auto x22 = x9 + x90;
auto x23 = x18*(x2*x22 + x7*((x22)*(x22)) + x99) + x98;
auto x24 = V{2}*cell[7] - V{2}*x23;
auto x25 = V{2}*x20 - V{2}*x21 + x24 + x94;
auto x26 = x14*x25;
auto x27 = V{1} - x13;
auto x28 = x15*x2 + x7*((x15)*(x15)) + V{1};
auto x29 = -x5 - x6*(x10 - x28);
auto x30 = V{-2};
auto x31 = -V{1}*x30;
auto x32 = -x11 + x20*x31 + x24 - x29 + x95;
auto x33 = x27*x32;
auto x34 = x12*(x26 + x33);
auto x35 = x11 + x34;
auto x36 = V{0.5}*x0*x92;
auto x37 = x14*x32;
auto x38 = x29 + x34;
auto x103 = x35 + x38;
auto x39 = V{0.5}*x17*x92;
auto x40 = x25*x27;
auto x101 = x33 + x40;
auto x41 = V{2};
auto x42 = x20*x30 + x24;
auto x43 = x23 + x39*(x101 + x41*x42);
auto x44 = x19 + x39*(x101 + x30*x42);
auto x45 = x12*(x37 + x40) + x21;
auto x102 = V{2}*x45;
auto x46 = x93*(-x0 - x1*(x10 + V{-1}) + x102 + x103 + x36*(x26 + x37) + V{2}*x43 + V{2}*x44 + V{1}) + V{-1};
auto x104 = ((x46)*(x46));
auto x47 = x104 + x4;
auto x48 = x3*x47;
auto x105 = -x48;
auto x106 = x105 + V{1};
auto x49 = x15 + x46;
auto x50 = x18*(x106 + x2*x49 + x7*((x49)*(x49))) + x98;
auto x51 = x44 - x50;
auto x52 = x100 + x6*(x104*x7 + x106 + x2*x46);
auto x53 = x46 + x90;
auto x54 = x18*(x106 + x2*x53 + x7*((x53)*(x53))) + x98;
auto x55 = V{2}*x43 - V{2}*x54;
auto x56 = x102 + V{2}*x51 - V{2}*x52 + x55;
auto x57 = x14*x56;
auto x58 = -x5 + x6*(x105 + x28);
auto x59 = -x5 - x6*(x48 + x97);
auto x60 = x103 + x31*x51 + x55 - x58 - x59;
auto x61 = x14*x60;
auto x62 = V{1} - x91;
auto x108 = x39*x62;
auto x109 = x12*x62;
auto x63 = -x0 - x1*(x48 + V{-1}) + x36*x62*(x57 + x61);
auto x64 = -V{1}*x46;
auto x65 = x64 + x90;
auto x66 = x27*x56;
auto x67 = x30*x51 + x55;
auto x68 = x27*x60;
auto x107 = x66 + x68;
auto x69 = x108*(x107 + x30*x67);
auto x70 = x109*(x61 + x66);
auto x71 = x15 + x64;
auto x72 = x108*(x107 + x41*x67);
auto x73 = x109*(x57 + x68);
auto x74 = x58 + x73;
auto x75 = x50 + x69;
auto x76 = x52 + x70;
auto x77 = x54 + x72;
auto x78 = x59 + x73;
cell[0] = x63;
cell[1] = x18*(x106 + x2*x65 + x7*((x65)*(x65))) + x69 + x98;
cell[2] = x100 + x6*(x106 + x2*x64 + x7*((x64)*(x64))) + x70;
cell[3] = x18*(x106 + x2*x71 + x7*((x71)*(x71))) + x72 + x98;
cell[4] = x74;
cell[5] = x75;
cell[6] = x76;
cell[7] = x77;
cell[8] = x78;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x88;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x93*(x63 + x74 + V{2}*x75 + V{2}*x76 + V{2}*x77 + x78 + V{1}) + V{-1};
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x90;
return { x88, x47 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
