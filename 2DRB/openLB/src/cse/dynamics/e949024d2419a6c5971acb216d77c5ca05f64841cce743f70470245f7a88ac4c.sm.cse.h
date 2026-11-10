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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineUSeparatelyTrace>, equilibria::SecondOrder, collision::ConstRhoBGK, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x46 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x47 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x48 = parameters.template get<descriptors::OMEGA>();
auto x49 = parameters.template get<statistics::AVERAGE_RHO>();
auto x50 = ((x46)*(x46));
auto x51 = ((x47)*(x47));
auto x0 = V{3};
auto x1 = x50 + x51;
auto x2 = V{0.5}*x0*x1;
auto x53 = V{1} - x2;
auto x3 = V{0.444444444444444};
auto x4 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + V{1};
auto x5 = -x3*(x4*(x2 + V{-1}) + V{1});
auto x6 = V{1} - x48;
auto x7 = V{1} - x49;
auto x8 = V{1} + x7/x4;
auto x9 = -x3 + x6*(cell[0] - x5) + x8*(x3 + x5);
auto x10 = -V{1}*x46;
auto x11 = x10 + x47;
auto x12 = V{0.5}*((x0)*(x0));
auto x13 = V{0.0277777777777778};
auto x52 = -x13;
auto x14 = x13*x4;
auto x15 = x14*(x0*x11 + x12*((x11)*(x11)) + x53) + x52;
auto x16 = x52 + x6*(cell[1] - x15) + x8*(x13 + x15);
auto x17 = V{0.111111111111111};
auto x54 = -x17;
auto x18 = x17*x4;
auto x19 = x18*(x0*x10 + x12*((x10)*(x10)) + x53) + x54;
auto x20 = x54 + x6*(cell[2] - x19) + x8*(x17 + x19);
auto x21 = -V{1}*x47;
auto x22 = x10 + x21;
auto x23 = x14*(x0*x22 + x12*((x22)*(x22)) + x53) + x52;
auto x24 = x52 + x6*(cell[3] - x23) + x8*(x13 + x23);
auto x25 = x18*(x0*x21 + x12*((x21)*(x21)) + x53) + x54;
auto x26 = x54 + x6*(cell[4] - x25) + x8*(x17 + x25);
auto x27 = x21 + x46;
auto x28 = x14*(x0*x27 + x12*((x27)*(x27)) + x53) + x52;
auto x29 = x52 + x6*(cell[5] - x28) + x8*(x13 + x28);
auto x30 = x18*(x0*x46 + x12*x50 + x53) + x54;
auto x31 = x54 + x6*(cell[6] - x30) + x8*(x17 + x30);
auto x32 = x46 + x47;
auto x33 = x14*(x0*x32 + x12*((x32)*(x32)) + x53) + x52;
auto x34 = x52 + x6*(cell[7] - x33) + x8*(x13 + x33);
auto x35 = x18*(x0*x47 + x12*x51 + x53) + x54;
auto x36 = x54 + x6*(cell[8] - x35) + x8*(x17 + x35);
cell[0] = x9;
cell[1] = x16;
cell[2] = x20;
cell[3] = x24;
cell[4] = x26;
cell[5] = x29;
cell[6] = x31;
cell[7] = x34;
cell[8] = x36;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x16 + x20 + x24 + x26 + x29 + x31 + x34 + x36 + x9 + V{1};
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x46;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x47;
return { x4 + x7, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
