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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SaveVelocity<collision::BGK>, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x44 = parameters.template get<descriptors::OMEGA>();
auto x42 = cell[5] + cell[6];
auto x43 = cell[1] + cell[8];
auto x0 = V{1} - x44;
auto x1 = V{3};
auto x2 = -V{1}*cell[3];
auto x46 = cell[7] + x2;
auto x3 = cell[0] + cell[2] + cell[3] + cell[4] + cell[7] + x42 + x43 + V{1};
auto x45 = V{1} / (x3);
auto x4 = x45*(-V{1}*cell[1] - V{1}*cell[2] + x42 + x46);
auto x47 = ((x4)*(x4));
auto x5 = x45*(-V{1}*cell[4] - V{1}*cell[5] + x43 + x46);
auto x48 = ((x5)*(x5));
auto x6 = x47 + x48;
auto x7 = V{0.5}*x1*x6;
auto x49 = V{1} - x7;
auto x8 = V{0.444444444444444};
auto x9 = cell[0]*x0 - x44*x8*(x3*(x7 + V{-1}) + V{1});
auto x10 = -V{1}*x4;
auto x11 = x10 + x5;
auto x12 = V{0.5}*((x1)*(x1));
auto x13 = V{0.0277777777777778};
auto x14 = x13*x3;
auto x15 = cell[1]*x0 - x44*(x13 - x14*(x1*x11 + x12*((x11)*(x11)) + x49));
auto x16 = V{0.111111111111111};
auto x17 = x16*x3;
auto x18 = cell[2]*x0 - x44*(x16 - x17*(x1*x10 + x12*((x10)*(x10)) + x49));
auto x19 = -V{1}*x5;
auto x20 = x10 + x19;
auto x21 = cell[3]*x0 - x44*(x13 - x14*(x1*x20 + x12*((x20)*(x20)) + x49));
auto x22 = cell[4]*x0 - x44*(x16 - x17*(x1*x19 + x12*((x19)*(x19)) + x49));
auto x23 = x19 + x4;
auto x24 = cell[5]*x0 - x44*(x13 - x14*(x1*x23 + x12*((x23)*(x23)) + x49));
auto x25 = cell[6]*x0 - x44*(x16 - x17*(x1*x4 + x12*x47 + x49));
auto x50 = x24 + x25;
auto x26 = x4 + x5;
auto x27 = cell[7]*x0 - x44*(x13 - x14*(x1*x26 + x12*((x26)*(x26)) + x49));
auto x28 = cell[8]*x0 - x44*(x16 - x17*(x1*x5 + x12*x48 + x49));
auto x51 = x15 + x28;
auto x29 = -V{1}*x21;
auto x53 = x27 + x29;
auto x30 = x18 + x21 + x22 + x27 + x50 + x51 + x9 + V{1};
auto x52 = V{1} / (x30);
auto x31 = x52*(-V{1}*x15 - V{1}*x18 + x50 + x53);
auto x32 = x52*(-V{1}*x22 - V{1}*x24 + x51 + x53);
cell[0] = x9;
cell[1] = x15;
cell[2] = x18;
cell[3] = x21;
cell[4] = x22;
cell[5] = x24;
cell[6] = x25;
cell[7] = x27;
cell[8] = x28;
cell.template getFieldPointer<olb::descriptors::VELOCITY>()[0] = x31;
cell.template getFieldPointer<olb::descriptors::VELOCITY>()[1] = x32;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x30;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x31;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x32;
return { x3, x6 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
