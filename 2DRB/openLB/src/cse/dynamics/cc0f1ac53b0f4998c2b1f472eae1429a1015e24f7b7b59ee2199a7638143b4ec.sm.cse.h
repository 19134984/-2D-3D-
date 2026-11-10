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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, -1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<0, -1>, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x89 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x90 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x91 = parameters.template get<descriptors::OMEGA>();
auto x92 = ((x89)*(x89));
auto x93 = ((x90)*(x90));
auto x95 = cell[4] + cell[8];
auto x0 = V{3};
auto x96 = ((x0)*(x0));
auto x1 = x92 + x93;
auto x2 = V{0.5}*x0*x1;
auto x97 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.444444444444444};
auto x101 = -x4;
auto x102 = x3*x4;
auto x5 = V{0.111111111111111};
auto x98 = -x5;
auto x6 = -V{1}*x89;
auto x7 = x6 + V{1};
auto x94 = V{1} / (x7);
auto x8 = x94*(cell[0] + V{2}*cell[1] + V{2}*cell[2] + V{2}*cell[3] + x95 + V{1});
auto x9 = x5*x8;
auto x10 = V{0.5}*x96;
auto x11 = x0*x90 + x10*x93 + x97;
auto x12 = x11*x9 + x98;
auto x13 = V{0.5}*x5*x96;
auto x14 = V{1}/x0;
auto x15 = -x14;
auto x16 = V{-2};
auto x17 = -V{1}*x16;
auto x18 = x6 + x90;
auto x19 = x0*x18 + x10*((x18)*(x18)) + x97;
auto x20 = V{0.0277777777777778};
auto x99 = -x20;
auto x21 = x20*x8;
auto x22 = x19*x21 + x99;
auto x23 = cell[1] - x22;
auto x24 = x0*x6 + x10*((x6)*(x6)) + x97;
auto x25 = x24*x9 + x98;
auto x26 = -V{1}*x90;
auto x27 = x26 + x6;
auto x28 = x0*x27 + x10*((x27)*(x27)) + x97;
auto x29 = -x20 + x21*x28;
auto x30 = x17*(cell[3] - x29);
auto x31 = x17*x23 + x17*(cell[2] - x25) + x30;
auto x32 = x15*x31;
auto x33 = V{1} - x14;
auto x34 = x0*x26 + x10*((x26)*(x26)) + x97;
auto x35 = x34*x9 + x98;
auto x36 = -x12 + V{2}*x23 + x30 - x35 + x95;
auto x37 = x33*x36;
auto x38 = x13*(x32 + x37);
auto x39 = x12 + x38;
auto x40 = V{0.5}*x4*x96;
auto x41 = x15*x36;
auto x42 = x35 + x38;
auto x103 = x39 + x42;
auto x43 = V{0.5}*x20*x96;
auto x44 = x31*x33;
auto x100 = x37 + x44;
auto x45 = V{2};
auto x46 = x16*x23 + x30;
auto x47 = x29 + x43*(x100 + x45*x46);
auto x48 = x22 + x43*(x100 + x16*x46);
auto x49 = x13*(x41 + x44) + x25;
auto x50 = x94*(x101 + x102*x8 + x103 + x40*(x32 + x41) + V{2}*x47 + V{2}*x48 + V{2}*x49 + V{1});
auto x51 = x20*x50;
auto x52 = x19*x51 + x99;
auto x53 = x48 - x52;
auto x54 = x5*x50;
auto x55 = x24*x54 + x98;
auto x56 = -x20 + x28*x51;
auto x57 = x17*(x47 - x56);
auto x58 = x17*x53 + x17*(x49 - x55) + x57;
auto x59 = x15*x58;
auto x60 = x34*x54 + x98;
auto x61 = x11*x54 + x98;
auto x62 = x103 + V{2}*x53 + x57 - x60 - x61;
auto x63 = x15*x62;
auto x64 = V{1} - x91;
auto x105 = x43*x64;
auto x106 = x13*x64;
auto x65 = x101 + x102*x50 + x40*x64*(x59 + x63);
auto x66 = x33*x58;
auto x67 = x16*x53 + x57;
auto x68 = x33*x62;
auto x104 = x66 + x68;
auto x69 = x105*(x104 + x16*x67);
auto x70 = x52 + x69;
auto x71 = x106*(x63 + x66);
auto x72 = x55 + x71;
auto x73 = x105*(x104 + x45*x67);
auto x74 = x56 + x73;
auto x75 = x106*(x59 + x68);
auto x76 = x60 + x75;
auto x77 = x26 + x89;
auto x78 = x89 + x90;
auto x79 = x61 + x75;
cell[0] = x65;
cell[1] = x70;
cell[2] = x72;
cell[3] = x74;
cell[4] = x76;
cell[5] = x51*(x0*x77 + x10*((x77)*(x77)) + x97) + x69 + x99;
cell[6] = x54*(x0*x89 + x10*x92 + x97) + x71 + x98;
cell[7] = x51*(x0*x78 + x10*((x78)*(x78)) + x97) + x73 + x99;
cell[8] = x79;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x94*(x65 + V{2}*x70 + V{2}*x72 + V{2}*x74 + x76 + x79 + V{1});
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x89;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x90;
return { x50, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
