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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::SecondOrder, collision::BGK, forcing::ShanChen>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x38 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x39 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x40 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x41 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x42 = parameters.template get<descriptors::OMEGA>();
auto x43 = V{1} / (x42);
auto x0 = V{1} - x42;
auto x1 = V{3};
auto x2 = x38*x43 + x40;
auto x44 = ((x2)*(x2));
auto x3 = x39*x43 + x41;
auto x45 = ((x3)*(x3));
auto x4 = V{0.5}*x1*(x44 + x45);
auto x46 = V{1} - x4;
auto x5 = V{0.444444444444444};
auto x6 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + V{1};
auto x7 = cell[0]*x0 - x42*x5*(x6*(x4 + V{-1}) + V{1});
auto x8 = -V{1}*x2;
auto x9 = x3 + x8;
auto x10 = V{0.5}*((x1)*(x1));
auto x11 = V{0.0277777777777778};
auto x12 = x11*x6;
auto x13 = cell[1]*x0 - x42*(x11 - x12*(x1*x9 + x10*((x9)*(x9)) + x46));
auto x14 = V{0.111111111111111};
auto x15 = x14*x6;
auto x16 = cell[2]*x0 - x42*(x14 - x15*(x1*x8 + x10*((x8)*(x8)) + x46));
auto x17 = -V{1}*x3;
auto x18 = x17 + x8;
auto x19 = cell[3]*x0 - x42*(x11 - x12*(x1*x18 + x10*((x18)*(x18)) + x46));
auto x20 = cell[4]*x0 - x42*(x14 - x15*(x1*x17 + x10*((x17)*(x17)) + x46));
auto x21 = x17 + x2;
auto x22 = cell[5]*x0 - x42*(x11 - x12*(x1*x21 + x10*((x21)*(x21)) + x46));
auto x23 = cell[6]*x0 - x42*(x14 - x15*(x1*x2 + x10*x44 + x46));
auto x24 = x2 + x3;
auto x25 = cell[7]*x0 - x42*(x11 - x12*(x1*x24 + x10*((x24)*(x24)) + x46));
auto x26 = cell[8]*x0 - x42*(x14 - x15*(x1*x3 + x10*x45 + x46));
auto x27 = V{0.5}*x38 + x40;
auto x28 = V{0.5}*x39 + x41;
cell[0] = x7;
cell[1] = x13;
cell[2] = x16;
cell[3] = x19;
cell[4] = x20;
cell[5] = x22;
cell[6] = x23;
cell[7] = x25;
cell[8] = x26;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x13 + x16 + x19 + x20 + x22 + x23 + x25 + x26 + x7 + V{1};
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x27;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x28;
return { x6, ((x27)*(x27)) + ((x28)*(x28)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
