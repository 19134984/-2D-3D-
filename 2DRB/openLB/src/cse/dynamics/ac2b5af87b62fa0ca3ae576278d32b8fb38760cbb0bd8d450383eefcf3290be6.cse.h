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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::SecondOrder, collision::BGK, forcing::ShanChen>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x50 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x51 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x52 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x53 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x54 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x55 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x56 = parameters.template get<descriptors::OMEGA>();
auto x57 = V{1} / (x56);
auto x0 = V{1} - x56;
auto x1 = V{3};
auto x2 = x50*x57 + x53;
auto x58 = ((x2)*(x2));
auto x3 = x51*x57 + x54;
auto x59 = ((x3)*(x3));
auto x4 = x52*x57 + x55;
auto x60 = ((x4)*(x4));
auto x5 = V{0.5}*x1*(x58 + x59 + x60);
auto x61 = V{1} - x5;
auto x6 = V{0.333333333333333};
auto x7 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{1};
auto x8 = -V{1}*x2;
auto x9 = V{0.5}*((x1)*(x1));
auto x10 = V{0.0555555555555556};
auto x11 = x10*x7;
auto x12 = -V{1}*x3;
auto x13 = -V{1}*x4;
auto x14 = x12 + x8;
auto x15 = V{0.0277777777777778};
auto x16 = x15*x7;
auto x17 = x3 + x8;
auto x18 = x13 + x8;
auto x19 = x4 + x8;
auto x20 = x12 + x13;
auto x21 = x12 + x4;
auto x22 = x2 + x3;
auto x23 = x12 + x2;
auto x24 = x2 + x4;
auto x25 = x13 + x2;
auto x26 = x3 + x4;
auto x27 = x13 + x3;
auto x28 = V{0.5}*x50 + x53;
auto x29 = V{0.5}*x51 + x54;
auto x30 = V{0.5}*x52 + x55;
auto x31 = cell[0]*x0 - x56*x6*(x7*(x5 + V{-1}) + V{1});
auto x32 = cell[10]*x0 - x56*(x10 - x11*(x1*x2 + x58*x9 + x61));
auto x33 = cell[11]*x0 - x56*(x10 - x11*(x1*x3 + x59*x9 + x61));
auto x34 = cell[12]*x0 - x56*(x10 - x11*(x1*x4 + x60*x9 + x61));
auto x35 = cell[13]*x0 - x56*(x15 - x16*(x1*x22 + x61 + x9*((x22)*(x22))));
auto x36 = cell[14]*x0 - x56*(x15 - x16*(x1*x23 + x61 + x9*((x23)*(x23))));
auto x37 = cell[15]*x0 - x56*(x15 - x16*(x1*x24 + x61 + x9*((x24)*(x24))));
auto x38 = cell[16]*x0 - x56*(x15 - x16*(x1*x25 + x61 + x9*((x25)*(x25))));
auto x39 = cell[17]*x0 - x56*(x15 - x16*(x1*x26 + x61 + x9*((x26)*(x26))));
auto x40 = cell[18]*x0 - x56*(x15 - x16*(x1*x27 + x61 + x9*((x27)*(x27))));
auto x41 = cell[1]*x0 - x56*(x10 - x11*(x1*x8 + x61 + x9*((x8)*(x8))));
auto x42 = cell[2]*x0 - x56*(x10 - x11*(x1*x12 + x61 + x9*((x12)*(x12))));
auto x43 = cell[3]*x0 - x56*(x10 - x11*(x1*x13 + x61 + x9*((x13)*(x13))));
auto x44 = cell[4]*x0 - x56*(x15 - x16*(x1*x14 + x61 + x9*((x14)*(x14))));
auto x45 = cell[5]*x0 - x56*(x15 - x16*(x1*x17 + x61 + x9*((x17)*(x17))));
auto x46 = cell[6]*x0 - x56*(x15 - x16*(x1*x18 + x61 + x9*((x18)*(x18))));
auto x47 = cell[7]*x0 - x56*(x15 - x16*(x1*x19 + x61 + x9*((x19)*(x19))));
auto x48 = cell[8]*x0 - x56*(x15 - x16*(x1*x20 + x61 + x9*((x20)*(x20))));
auto x49 = cell[9]*x0 - x56*(x15 - x16*(x1*x21 + x61 + x9*((x21)*(x21))));
cell[0] = x31;
cell[10] = x32;
cell[11] = x33;
cell[12] = x34;
cell[13] = x35;
cell[14] = x36;
cell[15] = x37;
cell[16] = x38;
cell[17] = x39;
cell[18] = x40;
cell[1] = x41;
cell[2] = x42;
cell[3] = x43;
cell[4] = x44;
cell[5] = x45;
cell[6] = x46;
cell[7] = x47;
cell[8] = x48;
cell[9] = x49;
return { x7, ((x28)*(x28)) + ((x29)*(x29)) + ((x30)*(x30)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
