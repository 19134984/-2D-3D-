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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::SecondOrder, collision::BGK, forcing::ShanChen>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x69 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x70 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x71 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x72 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x73 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x74 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x75 = parameters.template get<descriptors::OMEGA>();
auto x76 = V{1} / (x75);
auto x0 = V{1} - x75;
auto x1 = V{3};
auto x2 = x69*x76 + x72;
auto x77 = ((x2)*(x2));
auto x3 = x70*x76 + x73;
auto x78 = ((x3)*(x3));
auto x4 = x71*x76 + x74;
auto x79 = ((x4)*(x4));
auto x5 = V{0.5}*x1*(x77 + x78 + x79);
auto x80 = V{1} - x5;
auto x6 = V{0.333333333333333};
auto x7 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{1};
auto x8 = cell[0]*x0 - x6*x75*(x7*(x5 + V{-1}) + V{1});
auto x9 = -V{1}*x2;
auto x10 = V{0.5}*((x1)*(x1));
auto x11 = V{0.0555555555555556};
auto x12 = x11*x7;
auto x13 = cell[1]*x0 - x75*(x11 - x12*(x1*x9 + x10*((x9)*(x9)) + x80));
auto x14 = -V{1}*x3;
auto x15 = cell[2]*x0 - x75*(x11 - x12*(x1*x14 + x10*((x14)*(x14)) + x80));
auto x16 = -V{1}*x4;
auto x17 = cell[3]*x0 - x75*(x11 - x12*(x1*x16 + x10*((x16)*(x16)) + x80));
auto x18 = x14 + x9;
auto x19 = V{0.0277777777777778};
auto x20 = x19*x7;
auto x21 = cell[4]*x0 - x75*(x19 - x20*(x1*x18 + x10*((x18)*(x18)) + x80));
auto x22 = x3 + x9;
auto x23 = cell[5]*x0 - x75*(x19 - x20*(x1*x22 + x10*((x22)*(x22)) + x80));
auto x24 = x16 + x9;
auto x25 = cell[6]*x0 - x75*(x19 - x20*(x1*x24 + x10*((x24)*(x24)) + x80));
auto x26 = x4 + x9;
auto x27 = cell[7]*x0 - x75*(x19 - x20*(x1*x26 + x10*((x26)*(x26)) + x80));
auto x28 = x14 + x16;
auto x29 = cell[8]*x0 - x75*(x19 - x20*(x1*x28 + x10*((x28)*(x28)) + x80));
auto x30 = x14 + x4;
auto x31 = cell[9]*x0 - x75*(x19 - x20*(x1*x30 + x10*((x30)*(x30)) + x80));
auto x32 = cell[10]*x0 - x75*(x11 - x12*(x1*x2 + x10*x77 + x80));
auto x33 = cell[11]*x0 - x75*(x11 - x12*(x1*x3 + x10*x78 + x80));
auto x34 = cell[12]*x0 - x75*(x11 - x12*(x1*x4 + x10*x79 + x80));
auto x35 = x2 + x3;
auto x36 = cell[13]*x0 - x75*(x19 - x20*(x1*x35 + x10*((x35)*(x35)) + x80));
auto x37 = x14 + x2;
auto x38 = cell[14]*x0 - x75*(x19 - x20*(x1*x37 + x10*((x37)*(x37)) + x80));
auto x39 = x2 + x4;
auto x40 = cell[15]*x0 - x75*(x19 - x20*(x1*x39 + x10*((x39)*(x39)) + x80));
auto x41 = x16 + x2;
auto x42 = cell[16]*x0 - x75*(x19 - x20*(x1*x41 + x10*((x41)*(x41)) + x80));
auto x43 = x3 + x4;
auto x44 = cell[17]*x0 - x75*(x19 - x20*(x1*x43 + x10*((x43)*(x43)) + x80));
auto x45 = x16 + x3;
auto x46 = cell[18]*x0 - x75*(x19 - x20*(x1*x45 + x10*((x45)*(x45)) + x80));
auto x47 = V{0.5}*x69 + x72;
auto x48 = V{0.5}*x70 + x73;
auto x49 = V{0.5}*x71 + x74;
cell[0] = x8;
cell[10] = x32;
cell[11] = x33;
cell[12] = x34;
cell[13] = x36;
cell[14] = x38;
cell[15] = x40;
cell[16] = x42;
cell[17] = x44;
cell[18] = x46;
cell[1] = x13;
cell[2] = x15;
cell[3] = x17;
cell[4] = x21;
cell[5] = x23;
cell[6] = x25;
cell[7] = x27;
cell[8] = x29;
cell[9] = x31;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x13 + x15 + x17 + x21 + x23 + x25 + x27 + x29 + x31 + x32 + x33 + x34 + x36 + x38 + x40 + x42 + x44 + x46 + x8 + V{1};
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x47;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x48;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x49;
return { x7, ((x47)*(x47)) + ((x48)*(x48)) + ((x49)*(x49)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
