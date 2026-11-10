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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineUSeparatelyTrace>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x78 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x79 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x80 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x81 = parameters.template get<descriptors::OMEGA>();
auto x82 = ((x78)*(x78));
auto x83 = ((x79)*(x79));
auto x84 = ((x80)*(x80));
auto x85 = cell[15] + cell[16] + cell[6] + cell[7];
auto x86 = cell[10] + cell[1] + x85;
auto x87 = cell[17] + cell[18] + cell[8] + cell[9];
auto x88 = cell[11] + cell[2] + x87;
auto x89 = cell[12] + cell[3];
auto x0 = V{3};
auto x92 = ((x0)*(x0));
auto x1 = x82 + x83 + x84;
auto x2 = V{0.5}*x0*x1;
auto x99 = V{1} - x2;
auto x3 = V{0.333333333333333};
auto x4 = cell[0] + cell[13] + cell[14] + cell[4] + cell[5] + x86 + x88 + x89 + V{1};
auto x5 = V{1}/x0;
auto x6 = -x5;
auto x7 = x4*x78;
auto x8 = x5*(x4 + V{-1});
auto x90 = -x8;
auto x91 = cell[13] + cell[14] + cell[4] + cell[5] + x90;
auto x9 = -x7*x78 + x86 + x91;
auto x10 = x6*x9;
auto x11 = x4*x79;
auto x12 = -x11*x79 + x88 + x91;
auto x13 = x12*x6;
auto x14 = x10 + x13;
auto x15 = -x4*x84 + x85 + x87 + x89 + x90;
auto x16 = x15*x6;
auto x17 = V{1} - x81;
auto x18 = -V{1}*x78;
auto x19 = V{0.5}*x92;
auto x20 = V{0.0555555555555556};
auto x100 = -x20;
auto x21 = x20*x4;
auto x22 = V{0.5}*x20*x92;
auto x94 = x17*x22;
auto x23 = V{1} - x5;
auto x24 = x23*x9;
auto x93 = x16 + x24;
auto x25 = x94*(x13 + x93);
auto x101 = x100 + x25;
auto x26 = -V{1}*x79;
auto x27 = x12*x23;
auto x95 = x27 + x93;
auto x28 = x10 + x27;
auto x29 = x94*(x16 + x28);
auto x102 = x100 + x29;
auto x30 = -V{1}*x80;
auto x31 = x15*x23;
auto x97 = x13 + x24 + x31;
auto x98 = x28 + x31;
auto x32 = x94*(x14 + x31);
auto x103 = x100 + x32;
auto x33 = x18 + x26;
auto x34 = V{0.0277777777777778};
auto x104 = -x34;
auto x35 = x34*x4;
auto x36 = V{0.5}*x34*x92;
auto x96 = x17*x36;
auto x37 = V{2};
auto x38 = cell[13] - V{1}*cell[14] + cell[4] - V{1}*cell[5] - x7*x79;
auto x39 = x96*(x37*x38 + x95);
auto x105 = x104 + x39;
auto x40 = x18 + x79;
auto x41 = V{-2};
auto x42 = x96*(x38*x41 + x95);
auto x106 = x104 + x42;
auto x43 = x18 + x30;
auto x44 = cell[15] - V{1}*cell[16] + cell[6] - V{1}*cell[7] - x7*x80;
auto x45 = x96*(x37*x44 + x97);
auto x107 = x104 + x45;
auto x46 = x18 + x80;
auto x47 = x96*(x41*x44 + x97);
auto x108 = x104 + x47;
auto x48 = x26 + x30;
auto x49 = cell[17] - V{1}*cell[18] + cell[8] - V{1}*cell[9] - x11*x80;
auto x50 = x96*(x37*x49 + x98);
auto x109 = x104 + x50;
auto x51 = x26 + x80;
auto x52 = x96*(x41*x49 + x98);
auto x110 = x104 + x52;
auto x53 = x78 + x79;
auto x54 = x26 + x78;
auto x55 = x78 + x80;
auto x56 = x30 + x78;
auto x57 = x79 + x80;
auto x58 = x30 + x79;
cell[0] = x3*(V{0.5}*x17*x92*(x14 + x16) - x4*(x2 + V{-1}) + V{-1});
cell[10] = x101 + x21*(x0*x78 + x19*x82 + x99);
cell[11] = x102 + x21*(x0*x79 + x19*x83 + x99);
cell[12] = x103 + x21*(x0*x80 + x19*x84 + x99);
cell[13] = x105 + x35*(x0*x53 + x19*((x53)*(x53)) + x99);
cell[14] = x106 + x35*(x0*x54 + x19*((x54)*(x54)) + x99);
cell[15] = x107 + x35*(x0*x55 + x19*((x55)*(x55)) + x99);
cell[16] = x108 + x35*(x0*x56 + x19*((x56)*(x56)) + x99);
cell[17] = x109 + x35*(x0*x57 + x19*((x57)*(x57)) + x99);
cell[18] = x110 + x35*(x0*x58 + x19*((x58)*(x58)) + x99);
cell[1] = x101 + x21*(x0*x18 + x19*((x18)*(x18)) + x99);
cell[2] = x102 + x21*(x0*x26 + x19*((x26)*(x26)) + x99);
cell[3] = x103 + x21*(x0*x30 + x19*((x30)*(x30)) + x99);
cell[4] = x105 + x35*(x0*x33 + x19*((x33)*(x33)) + x99);
cell[5] = x106 + x35*(x0*x40 + x19*((x40)*(x40)) + x99);
cell[6] = x107 + x35*(x0*x43 + x19*((x43)*(x43)) + x99);
cell[7] = x108 + x35*(x0*x46 + x19*((x46)*(x46)) + x99);
cell[8] = x109 + x35*(x0*x48 + x19*((x48)*(x48)) + x99);
cell[9] = x110 + x35*(x0*x51 + x19*((x51)*(x51)) + x99);
return { x4, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
