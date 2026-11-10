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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<1, 1>, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x89 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x90 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x91 = parameters.template get<descriptors::OMEGA>();
auto x92 = ((x89)*(x89));
auto x93 = ((x90)*(x90));
auto x95 = V{2}*cell[8];
auto x96 = cell[2] + cell[6];
auto x0 = V{3};
auto x97 = ((x0)*(x0));
auto x1 = x92 + x93;
auto x2 = V{0.5}*x0*x1;
auto x98 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.444444444444444};
auto x101 = -x4;
auto x103 = x3*x4;
auto x5 = V{0.111111111111111};
auto x99 = -x5;
auto x6 = x90 + V{1};
auto x94 = V{1} / (x6);
auto x7 = x94*(cell[0] + V{2}*cell[1] + V{2}*cell[7] + x95 + x96 + V{1});
auto x8 = x5*x7;
auto x9 = V{0.5}*x97;
auto x10 = x0*x89 + x9*x92 + x98;
auto x11 = x10*x8 + x99;
auto x12 = V{0.5}*x5*x97;
auto x13 = V{1}/x0;
auto x14 = -x13;
auto x15 = V{0.0277777777777778};
auto x108 = -x15;
auto x16 = x15*x7;
auto x17 = x89 + x90;
auto x18 = x0*x17 + x9*((x17)*(x17)) + x98;
auto x19 = -x15 + x16*x18;
auto x20 = V{2}*cell[7] - V{2}*x19;
auto x21 = -V{1}*x89;
auto x22 = x21 + x90;
auto x23 = x0*x22 + x9*((x22)*(x22)) + x98;
auto x24 = -x15 + x16*x23;
auto x25 = cell[1] - x24;
auto x26 = x0*x90 + x9*x93 + x98;
auto x27 = x26*x8 + x99;
auto x28 = x20 + V{2}*x25 - V{2}*x27 + x95;
auto x29 = x14*x28;
auto x30 = V{1} - x13;
auto x31 = x0*x21 + x9*((x21)*(x21)) + x98;
auto x32 = x31*x8 + x99;
auto x33 = V{-2};
auto x34 = -V{1}*x33;
auto x35 = -x11 + x20 + x25*x34 - x32 + x96;
auto x36 = x30*x35;
auto x37 = x12*(x29 + x36);
auto x38 = x11 + x37;
auto x39 = V{0.5}*x4*x97;
auto x40 = x14*x35;
auto x41 = x32 + x37;
auto x104 = x38 + x41;
auto x42 = x28*x30;
auto x100 = x36 + x42;
auto x43 = x12*(x40 + x42) + x27;
auto x102 = V{2}*x43;
auto x44 = V{0.5}*x15*x97;
auto x45 = x20 + x25*x33;
auto x46 = x24 + x44*(x100 + x33*x45);
auto x47 = V{2};
auto x48 = x19 + x44*(x100 + x45*x47);
auto x49 = x94*(x101 + x102 + x103*x7 + x104 + x39*(x29 + x40) + V{2}*x46 + V{2}*x48 + V{1});
auto x50 = x49*x5;
auto x51 = x31*x50 + x99;
auto x52 = x10*x50 + x99;
auto x53 = x15*x49;
auto x54 = -x15 + x23*x53;
auto x55 = x46 - x54;
auto x56 = -x15 + x18*x53;
auto x57 = V{2}*x48 - V{2}*x56;
auto x58 = x104 + x34*x55 - x51 - x52 + x57;
auto x59 = x14*x58;
auto x60 = x26*x50 + x99;
auto x61 = x102 + V{2}*x55 + x57 - V{2}*x60;
auto x62 = x14*x61;
auto x63 = V{1} - x91;
auto x106 = x44*x63;
auto x107 = x12*x63;
auto x64 = x101 + x103*x49 + x39*x63*(x59 + x62);
auto x65 = x30*x58;
auto x66 = x33*x55 + x57;
auto x67 = x30*x61;
auto x105 = x65 + x67;
auto x68 = x106*(x105 + x33*x66);
auto x69 = x54 + x68;
auto x70 = x107*(x62 + x65);
auto x71 = x51 + x70;
auto x72 = -V{1}*x90;
auto x73 = x21 + x72;
auto x74 = x106*(x105 + x47*x66);
auto x75 = x107*(x59 + x67);
auto x76 = x72 + x89;
auto x77 = x52 + x70;
auto x78 = x56 + x74;
auto x79 = x60 + x75;
cell[0] = x64;
cell[1] = x69;
cell[2] = x71;
cell[3] = x108 + x53*(x0*x73 + x9*((x73)*(x73)) + x98) + x74;
cell[4] = x50*(x0*x72 + x9*((x72)*(x72)) + x98) + x75 + x99;
cell[5] = x108 + x53*(x0*x76 + x9*((x76)*(x76)) + x98) + x68;
cell[6] = x77;
cell[7] = x78;
cell[8] = x79;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x94*(x64 + V{2}*x69 + x71 + x77 + V{2}*x78 + V{2}*x79 + V{1});
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x89;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x90;
return { x49, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
