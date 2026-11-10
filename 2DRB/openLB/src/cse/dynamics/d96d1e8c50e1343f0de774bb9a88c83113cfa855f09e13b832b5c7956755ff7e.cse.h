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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<1, -1>, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x83 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x84 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x85 = parameters.template get<descriptors::OMEGA>();
auto x86 = ((x83)*(x83));
auto x87 = ((x84)*(x84));
auto x89 = cell[2] + cell[6];
auto x0 = V{3};
auto x90 = ((x0)*(x0));
auto x1 = x86 + x87;
auto x2 = V{0.5}*x0*x1;
auto x91 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.444444444444444};
auto x95 = -x4;
auto x96 = x3*x4;
auto x5 = V{0.111111111111111};
auto x92 = -x5;
auto x6 = -V{1}*x84;
auto x7 = x6 + V{1};
auto x88 = V{1} / (x7);
auto x8 = x88*(cell[0] + V{2}*cell[3] + V{2}*cell[4] + V{2}*cell[5] + x89 + V{1});
auto x9 = x5*x8;
auto x10 = V{0.5}*x90;
auto x11 = x0*x83 + x10*x86 + x91;
auto x12 = x11*x9 + x92;
auto x13 = V{0.5}*x5*x90;
auto x14 = V{1}/x0;
auto x15 = -x14;
auto x16 = V{-2};
auto x17 = -V{1}*x16;
auto x18 = -V{1}*x83;
auto x19 = x18 + x6;
auto x20 = x0*x19 + x10*((x19)*(x19)) + x91;
auto x21 = V{0.0277777777777778};
auto x93 = -x21;
auto x22 = x21*x8;
auto x23 = x20*x22 + x93;
auto x24 = x17*(cell[3] - x23);
auto x25 = x0*x6 + x10*((x6)*(x6)) + x91;
auto x26 = x25*x9 + x92;
auto x27 = x6 + x83;
auto x28 = x0*x27 + x10*((x27)*(x27)) + x91;
auto x29 = -x21 + x22*x28;
auto x30 = cell[5] - x29;
auto x31 = x17*x30 + x17*(cell[4] - x26) + x24;
auto x32 = x15*x31;
auto x33 = V{1} - x14;
auto x34 = x0*x18 + x10*((x18)*(x18)) + x91;
auto x35 = x34*x9 + x92;
auto x36 = -x12 + x24 + V{2}*x30 - x35 + x89;
auto x37 = x33*x36;
auto x38 = x13*(x32 + x37);
auto x39 = x12 + x38;
auto x40 = V{0.5}*x4*x90;
auto x41 = x15*x36;
auto x42 = x35 + x38;
auto x97 = x39 + x42;
auto x43 = V{0.5}*x21*x90;
auto x44 = x16*x30 + x24;
auto x45 = x31*x33;
auto x94 = x37 + x45;
auto x46 = x29 + x43*(x16*x44 + x94);
auto x47 = V{2};
auto x48 = x23 + x43*(x44*x47 + x94);
auto x49 = x13*(x41 + x45) + x26;
auto x50 = x88*(x40*(x32 + x41) + V{2}*x46 + V{2}*x48 + V{2}*x49 + x8*x96 + x95 + x97 + V{1});
auto x51 = x5*x50;
auto x52 = x34*x51 + x92;
auto x53 = x11*x51 + x92;
auto x54 = x21*x50;
auto x55 = x20*x54 + x93;
auto x56 = x17*(x48 - x55);
auto x57 = -x21 + x28*x54;
auto x58 = x46 - x57;
auto x59 = -x52 - x53 + x56 + V{2}*x58 + x97;
auto x60 = x15*x59;
auto x61 = x25*x51 + x92;
auto x62 = x17*x58 + x17*(x49 - x61) + x56;
auto x63 = x15*x62;
auto x64 = V{1} - x85;
auto x99 = x43*x64;
auto x100 = x13*x64;
auto x65 = x18 + x84;
auto x66 = x33*x59;
auto x67 = x16*x58 + x56;
auto x68 = x33*x62;
auto x98 = x66 + x68;
auto x69 = x99*(x16*x67 + x98);
auto x70 = x100*(x63 + x66);
auto x71 = x99*(x47*x67 + x98);
auto x72 = x100*(x60 + x68);
auto x73 = x83 + x84;
cell[0] = x40*x64*(x60 + x63) + x50*x96 + x95;
cell[1] = x54*(x0*x65 + x10*((x65)*(x65)) + x91) + x69 + x93;
cell[2] = x52 + x70;
cell[3] = x55 + x71;
cell[4] = x61 + x72;
cell[5] = x57 + x69;
cell[6] = x53 + x70;
cell[7] = x54*(x0*x73 + x10*((x73)*(x73)) + x91) + x71 + x93;
cell[8] = x51*(x0*x84 + x10*x87 + x91) + x72 + x92;
return { x50, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
