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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FreeEnergyMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::FreeEnergy, collision::FreeEnergy, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x38 = cell.template getFieldComponent<olb::descriptors::CHEM_POTENTIAL>(0);
auto x39 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x40 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x41 = parameters.template get<collision::FreeEnergy::GAMMA>();
auto x42 = parameters.template get<descriptors::OMEGA>();
auto x43 = ((x39)*(x39));
auto x44 = ((x40)*(x40));
auto x0 = V{1} - x42;
auto x1 = V{3};
auto x2 = x43 + x44;
auto x3 = V{0.5}*x1*x2;
auto x46 = V{1} - x3;
auto x4 = V{0.444444444444444};
auto x5 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + V{1};
auto x6 = -x1*x38*x41 + x5;
auto x45 = x42*x6;
auto x7 = cell[0]*x0 - x4*x42*(x5*(x3 + V{-1}) + V{1}) - x45*(x4 + V{-1});
auto x8 = -V{1}*x39;
auto x9 = x40 + x8;
auto x10 = V{0.5}*((x1)*(x1));
auto x11 = V{0.0277777777777778};
auto x12 = x11*x5;
auto x13 = x11*x45;
auto x14 = -(-cell[1]*x0 + x13 + x42*(x11 - x12*(x1*x9 + x10*((x9)*(x9)) + x46)));
auto x15 = V{0.111111111111111};
auto x16 = x15*x5;
auto x17 = x15*x45;
auto x18 = -(-cell[2]*x0 + x17 + x42*(x15 - x16*(x1*x8 + x10*((x8)*(x8)) + x46)));
auto x19 = -V{1}*x40;
auto x20 = x19 + x8;
auto x21 = -(-cell[3]*x0 + x13 + x42*(x11 - x12*(x1*x20 + x10*((x20)*(x20)) + x46)));
auto x22 = -(-cell[4]*x0 + x17 + x42*(x15 - x16*(x1*x19 + x10*((x19)*(x19)) + x46)));
auto x23 = x19 + x39;
auto x24 = -(-cell[5]*x0 + x13 + x42*(x11 - x12*(x1*x23 + x10*((x23)*(x23)) + x46)));
auto x25 = cell[6]*x0 - x17 - x42*(x15 - x16*(x1*x39 + x10*x43 + x46));
auto x26 = x39 + x40;
auto x27 = -(-cell[7]*x0 + x13 + x42*(x11 - x12*(x1*x26 + x10*((x26)*(x26)) + x46)));
auto x28 = cell[8]*x0 - x17 - x42*(x15 - x16*(x1*x40 + x10*x44 + x46));
cell[0] = x7;
cell[1] = x14;
cell[2] = x18;
cell[3] = x21;
cell[4] = x22;
cell[5] = x24;
cell[6] = x25;
cell[7] = x27;
cell[8] = x28;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x14 + x18 + x21 + x22 + x24 + x25 + x27 + x28 + x7 + V{1};
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x39;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x40;
return { x5, x2 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
