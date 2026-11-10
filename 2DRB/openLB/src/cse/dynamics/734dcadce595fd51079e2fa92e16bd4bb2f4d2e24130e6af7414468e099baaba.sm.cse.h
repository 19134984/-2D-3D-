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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::CBCDensity, momenta::CBCMomentum<1, 1>, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x82 = parameters.template get<descriptors::OMEGA>();
auto x74 = cell[10] + cell[14] + cell[16];
auto x75 = cell[11] + cell[18] + cell[5];
auto x76 = cell[12] + cell[7] + cell[9];
auto x0 = V{1} - x82;
auto x1 = V{3};
auto x2 = -V{1}*cell[4];
auto x78 = cell[13] + x2;
auto x3 = -V{1}*cell[6];
auto x79 = cell[15] + x3;
auto x4 = cell[0] + cell[13] + cell[15] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + x74 + x75 + x76 + V{1};
auto x77 = V{1} / (x4);
auto x5 = x77*(-V{1}*cell[1] - V{1}*cell[5] - V{1}*cell[7] + x74 + x78 + x79);
auto x81 = ((x5)*(x5));
auto x6 = -V{1}*cell[8];
auto x80 = cell[17] + x6;
auto x7 = x77*(-V{1}*cell[14] - V{1}*cell[2] - V{1}*cell[9] + x75 + x78 + x80);
auto x83 = ((x7)*(x7));
auto x8 = x77*(-V{1}*cell[16] - V{1}*cell[18] - V{1}*cell[3] + x76 + x79 + x80);
auto x84 = ((x8)*(x8));
auto x9 = x81 + x83 + x84;
auto x10 = V{0.5}*x1*x9;
auto x85 = V{1} - x10;
auto x11 = V{0.333333333333333};
auto x12 = cell[0]*x0 - x11*x82*(x4*(x10 + V{-1}) + V{1});
auto x13 = -V{1}*x5;
auto x14 = V{0.5}*((x1)*(x1));
auto x15 = V{0.0555555555555556};
auto x16 = x15*x4;
auto x17 = cell[1]*x0 - x82*(x15 - x16*(x1*x13 + x14*((x13)*(x13)) + x85));
auto x18 = -V{1}*x7;
auto x19 = cell[2]*x0 - x82*(x15 - x16*(x1*x18 + x14*((x18)*(x18)) + x85));
auto x20 = -V{1}*x8;
auto x21 = cell[3]*x0 - x82*(x15 - x16*(x1*x20 + x14*((x20)*(x20)) + x85));
auto x22 = x13 + x18;
auto x23 = V{0.0277777777777778};
auto x24 = x23*x4;
auto x25 = cell[4]*x0 - x82*(x23 - x24*(x1*x22 + x14*((x22)*(x22)) + x85));
auto x26 = x13 + x7;
auto x27 = cell[5]*x0 - x82*(x23 - x24*(x1*x26 + x14*((x26)*(x26)) + x85));
auto x28 = x13 + x20;
auto x29 = cell[6]*x0 - x82*(x23 - x24*(x1*x28 + x14*((x28)*(x28)) + x85));
auto x30 = x13 + x8;
auto x31 = cell[7]*x0 - x82*(x23 - x24*(x1*x30 + x14*((x30)*(x30)) + x85));
auto x32 = x18 + x20;
auto x33 = cell[8]*x0 - x82*(x23 - x24*(x1*x32 + x14*((x32)*(x32)) + x85));
auto x34 = x18 + x8;
auto x35 = cell[9]*x0 - x82*(x23 - x24*(x1*x34 + x14*((x34)*(x34)) + x85));
auto x36 = cell[10]*x0 - x82*(x15 - x16*(x1*x5 + x14*x81 + x85));
auto x37 = cell[11]*x0 - x82*(x15 - x16*(x1*x7 + x14*x83 + x85));
auto x38 = cell[12]*x0 - x82*(x15 - x16*(x1*x8 + x14*x84 + x85));
auto x88 = x31 + x35 + x38;
auto x39 = x5 + x7;
auto x40 = cell[13]*x0 - x82*(x23 - x24*(x1*x39 + x14*((x39)*(x39)) + x85));
auto x41 = x18 + x5;
auto x42 = cell[14]*x0 - x82*(x23 - x24*(x1*x41 + x14*((x41)*(x41)) + x85));
auto x43 = x5 + x8;
auto x44 = cell[15]*x0 - x82*(x23 - x24*(x1*x43 + x14*((x43)*(x43)) + x85));
auto x45 = x20 + x5;
auto x46 = cell[16]*x0 - x82*(x23 - x24*(x1*x45 + x14*((x45)*(x45)) + x85));
auto x86 = x36 + x42 + x46;
auto x47 = x7 + x8;
auto x48 = cell[17]*x0 - x82*(x23 - x24*(x1*x47 + x14*((x47)*(x47)) + x85));
auto x49 = x20 + x7;
auto x50 = cell[18]*x0 - x82*(x23 - x24*(x1*x49 + x14*((x49)*(x49)) + x85));
auto x87 = x27 + x37 + x50;
auto x51 = x12 + x17 + x19 + x21 + x25 + x29 + x33 + x40 + x44 + x48 + x86 + x87 + x88 + V{1};
auto x89 = V{1} / (x51);
auto x52 = -V{1}*x25;
auto x90 = x40 + x52;
auto x53 = -V{1}*x29;
auto x91 = x44 + x53;
auto x54 = -V{1}*x33;
auto x92 = x48 + x54;
cell[0] = x12;
cell[10] = x36;
cell[11] = x37;
cell[12] = x38;
cell[13] = x40;
cell[14] = x42;
cell[15] = x44;
cell[16] = x46;
cell[17] = x48;
cell[18] = x50;
cell[1] = x17;
cell[2] = x19;
cell[3] = x21;
cell[4] = x25;
cell[5] = x27;
cell[6] = x29;
cell[7] = x31;
cell[8] = x33;
cell[9] = x35;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x51;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x89*(-V{1}*x17 - V{1}*x27 - V{1}*x31 + x86 + x90 + x91);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x89*(-V{1}*x19 - V{1}*x35 - V{1}*x42 + x87 + x90 + x92);
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x89*(-V{1}*x21 - V{1}*x46 - V{1}*x50 + x88 + x91 + x92);
return { x4, x9 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
