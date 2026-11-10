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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<0, 1>, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x89 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x90 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x91 = parameters.template get<descriptors::OMEGA>();
auto x92 = ((x89)*(x89));
auto x93 = ((x90)*(x90));
auto x95 = V{2}*cell[6];
auto x96 = cell[4] + cell[8];
auto x0 = V{3};
auto x97 = ((x0)*(x0));
auto x1 = x92 + x93;
auto x2 = V{0.5}*x0*x1;
auto x98 = V{1} - x2;
auto x3 = V{1} - x2;
auto x4 = V{0.444444444444444};
auto x102 = -x4;
auto x104 = x3*x4;
auto x5 = V{0.111111111111111};
auto x99 = -x5;
auto x6 = x89 + V{1};
auto x94 = V{1} / (x6);
auto x7 = x94*(cell[0] + V{2}*cell[5] + V{2}*cell[7] + x95 + x96 + V{1});
auto x8 = x5*x7;
auto x9 = V{0.5}*x97;
auto x10 = x0*x90 + x9*x93 + x98;
auto x11 = x10*x8 + x99;
auto x12 = V{0.5}*x5*x97;
auto x13 = V{1}/x0;
auto x14 = -x13;
auto x15 = -V{1}*x90;
auto x16 = x15 + x89;
auto x17 = x0*x16 + x9*((x16)*(x16)) + x98;
auto x18 = V{0.0277777777777778};
auto x100 = -x18;
auto x19 = x18*x7;
auto x20 = x100 + x17*x19;
auto x21 = cell[5] - x20;
auto x22 = x0*x89 + x9*x92 + x98;
auto x23 = x22*x8 + x99;
auto x24 = x89 + x90;
auto x25 = x0*x24 + x9*((x24)*(x24)) + x98;
auto x26 = -x18 + x19*x25;
auto x27 = V{2}*cell[7] - V{2}*x26;
auto x28 = V{2}*x21 - V{2}*x23 + x27 + x95;
auto x29 = x14*x28;
auto x30 = V{1} - x13;
auto x31 = x0*x15 + x9*((x15)*(x15)) + x98;
auto x32 = x31*x8 + x99;
auto x33 = V{-2};
auto x34 = -V{1}*x33;
auto x35 = -x11 + x21*x34 + x27 - x32 + x96;
auto x36 = x30*x35;
auto x37 = x12*(x29 + x36);
auto x38 = x11 + x37;
auto x39 = V{0.5}*x4*x97;
auto x40 = x14*x35;
auto x41 = x32 + x37;
auto x105 = x38 + x41;
auto x42 = V{0.5}*x18*x97;
auto x43 = x28*x30;
auto x101 = x36 + x43;
auto x44 = V{2};
auto x45 = x21*x33 + x27;
auto x46 = x26 + x42*(x101 + x44*x45);
auto x47 = x20 + x42*(x101 + x33*x45);
auto x48 = x12*(x40 + x43) + x23;
auto x103 = V{2}*x48;
auto x49 = x94*(x102 + x103 + x104*x7 + x105 + x39*(x29 + x40) + V{2}*x46 + V{2}*x47 + V{1});
auto x50 = x18*x49;
auto x51 = x100 + x17*x50;
auto x52 = x47 - x51;
auto x53 = x49*x5;
auto x54 = x22*x53 + x99;
auto x55 = -x18 + x25*x50;
auto x56 = V{2}*x46 - V{2}*x55;
auto x57 = x103 + V{2}*x52 - V{2}*x54 + x56;
auto x58 = x14*x57;
auto x59 = x31*x53 + x99;
auto x60 = x10*x53 + x99;
auto x61 = x105 + x34*x52 + x56 - x59 - x60;
auto x62 = x14*x61;
auto x63 = V{1} - x91;
auto x107 = x42*x63;
auto x108 = x12*x63;
auto x64 = x102 + x104*x49 + x39*x63*(x58 + x62);
auto x65 = -V{1}*x89;
auto x66 = x65 + x90;
auto x67 = x30*x57;
auto x68 = x33*x52 + x56;
auto x69 = x30*x61;
auto x106 = x67 + x69;
auto x70 = x107*(x106 + x33*x68);
auto x71 = x108*(x62 + x67);
auto x72 = x15 + x65;
auto x73 = x107*(x106 + x44*x68);
auto x74 = x108*(x58 + x69);
auto x75 = x59 + x74;
auto x76 = x51 + x70;
auto x77 = x54 + x71;
auto x78 = x55 + x73;
auto x79 = x60 + x74;
cell[0] = x64;
cell[1] = x100 + x50*(x0*x66 + x9*((x66)*(x66)) + x98) + x70;
cell[2] = x53*(x0*x65 + x9*((x65)*(x65)) + x98) + x71 + x99;
cell[3] = x100 + x50*(x0*x72 + x9*((x72)*(x72)) + x98) + x73;
cell[4] = x75;
cell[5] = x76;
cell[6] = x77;
cell[7] = x78;
cell[8] = x79;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x94*(x64 + x75 + V{2}*x76 + V{2}*x77 + V{2}*x78 + x79 + V{1});
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x89;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x90;
return { x49, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
