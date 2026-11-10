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
struct CSE<dynamics::Tuple<T, descriptors::D3Q27<FIELDS...>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, equilibria::FirstOrder, collision::FixedEquilibrium, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x42 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x43 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x44 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x45 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x0 = V{0.296296296296296};
auto x1 = V{3};
auto x2 = -V{1}*x43;
auto x3 = V{0.0740740740740741};
auto x4 = x3*x42;
auto x5 = -V{1}*x44;
auto x6 = -V{1}*x45;
auto x7 = x2 + x5;
auto x8 = V{0.0185185185185185};
auto x9 = x42*x8;
auto x10 = x2 + x44;
auto x11 = V{0.00462962962962963};
auto x12 = x11*x42;
auto x13 = x43 + x44;
auto x14 = x43 + x5;
cell[0] = x0*(x42 + V{-1});
cell[10] = -x11 + x12*(x1*(x6 + x7) + V{1});
cell[11] = -x11 + x12*(x1*(x45 + x7) + V{1});
cell[12] = -x11 + x12*(x1*(x10 + x6) + V{1});
cell[13] = -x11 + x12*(x1*(x10 + x45) + V{1});
cell[14] = -x3 + x4*(x1*x43 + V{1});
cell[15] = -x3 + x4*(x1*x44 + V{1});
cell[16] = -x3 + x4*(x1*x45 + V{1});
cell[17] = -x8 + x9*(x1*x13 + V{1});
cell[18] = -x8 + x9*(x1*x14 + V{1});
cell[19] = -x8 + x9*(x1*(x43 + x45) + V{1});
cell[1] = -x3 + x4*(x1*x2 + V{1});
cell[20] = -x8 + x9*(x1*(x43 + x6) + V{1});
cell[21] = -x8 + x9*(x1*(x44 + x45) + V{1});
cell[22] = -x8 + x9*(x1*(x44 + x6) + V{1});
cell[23] = -x11 + x12*(x1*(x13 + x45) + V{1});
cell[24] = -x11 + x12*(x1*(x13 + x6) + V{1});
cell[25] = -x11 + x12*(x1*(x14 + x45) + V{1});
cell[26] = -x11 + x12*(x1*(x14 + x6) + V{1});
cell[2] = -x3 + x4*(x1*x5 + V{1});
cell[3] = -x3 + x4*(x1*x6 + V{1});
cell[4] = -x8 + x9*(x1*x7 + V{1});
cell[5] = -x8 + x9*(x1*x10 + V{1});
cell[6] = -x8 + x9*(x1*(x2 + x6) + V{1});
cell[7] = -x8 + x9*(x1*(x2 + x45) + V{1});
cell[8] = -x8 + x9*(x1*(x5 + x6) + V{1});
cell[9] = -x8 + x9*(x1*(x45 + x5) + V{1});
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x42;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x43;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x44;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x45;
return { x42, ((x43)*(x43)) + ((x44)*(x44)) + ((x45)*(x45)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
