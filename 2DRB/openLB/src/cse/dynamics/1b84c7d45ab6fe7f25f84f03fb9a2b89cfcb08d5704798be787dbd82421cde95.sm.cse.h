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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::FixedDensity, momenta::FixedVelocityMomentumGeneric, momenta::ZeroStress, momenta::DefineSeparately>, equilibria::FirstOrder, collision::FixedEquilibrium, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x28 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x29 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x30 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x31 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x0 = V{0.333333333333333};
auto x1 = V{3};
auto x2 = -V{1}*x29;
auto x3 = V{0.0555555555555556};
auto x4 = x28*x3;
auto x5 = -V{1}*x30;
auto x6 = -V{1}*x31;
auto x7 = V{0.0277777777777778};
auto x8 = x28*x7;
cell[0] = x0*(x28 + V{-1});
cell[10] = -x3 + x4*(x1*x29 + V{1});
cell[11] = -x3 + x4*(x1*x30 + V{1});
cell[12] = -x3 + x4*(x1*x31 + V{1});
cell[13] = -x7 + x8*(x1*(x29 + x30) + V{1});
cell[14] = -x7 + x8*(x1*(x29 + x5) + V{1});
cell[15] = -x7 + x8*(x1*(x29 + x31) + V{1});
cell[16] = -x7 + x8*(x1*(x29 + x6) + V{1});
cell[17] = -x7 + x8*(x1*(x30 + x31) + V{1});
cell[18] = -x7 + x8*(x1*(x30 + x6) + V{1});
cell[1] = -x3 + x4*(x1*x2 + V{1});
cell[2] = -x3 + x4*(x1*x5 + V{1});
cell[3] = -x3 + x4*(x1*x6 + V{1});
cell[4] = -x7 + x8*(x1*(x2 + x5) + V{1});
cell[5] = -x7 + x8*(x1*(x2 + x30) + V{1});
cell[6] = -x7 + x8*(x1*(x2 + x6) + V{1});
cell[7] = -x7 + x8*(x1*(x2 + x31) + V{1});
cell[8] = -x7 + x8*(x1*(x5 + x6) + V{1});
cell[9] = -x7 + x8*(x1*(x31 + x5) + V{1});
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x28;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x29;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x30;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x31;
return { x28, ((x29)*(x29)) + ((x30)*(x30)) + ((x31)*(x31)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
