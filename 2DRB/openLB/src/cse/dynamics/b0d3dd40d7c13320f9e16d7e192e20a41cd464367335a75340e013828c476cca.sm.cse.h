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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineUSeparatelyTrace>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x35 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x36 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x37 = parameters.template get<descriptors::OMEGA>();
auto x38 = ((x35)*(x35));
auto x39 = ((x36)*(x36));
auto x0 = V{1} - x37;
auto x1 = V{3};
auto x2 = x38 + x39;
auto x3 = V{0.5}*x1*x2;
auto x40 = V{1} - x3;
auto x4 = V{0.444444444444444};
auto x5 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + V{1};
auto x6 = cell[0]*x0 - x37*x4*(x5*(x3 + V{-1}) + V{1});
auto x7 = -V{1}*x35;
auto x8 = x36 + x7;
auto x9 = V{0.5}*((x1)*(x1));
auto x10 = V{0.0277777777777778};
auto x11 = x10*x5;
auto x12 = cell[1]*x0 - x37*(x10 - x11*(x1*x8 + x40 + x9*((x8)*(x8))));
auto x13 = V{0.111111111111111};
auto x14 = x13*x5;
auto x15 = cell[2]*x0 - x37*(x13 - x14*(x1*x7 + x40 + x9*((x7)*(x7))));
auto x16 = -V{1}*x36;
auto x17 = x16 + x7;
auto x18 = cell[3]*x0 - x37*(x10 - x11*(x1*x17 + x40 + x9*((x17)*(x17))));
auto x19 = cell[4]*x0 - x37*(x13 - x14*(x1*x16 + x40 + x9*((x16)*(x16))));
auto x20 = x16 + x35;
auto x21 = cell[5]*x0 - x37*(x10 - x11*(x1*x20 + x40 + x9*((x20)*(x20))));
auto x22 = cell[6]*x0 - x37*(x13 - x14*(x1*x35 + x38*x9 + x40));
auto x23 = x35 + x36;
auto x24 = cell[7]*x0 - x37*(x10 - x11*(x1*x23 + x40 + x9*((x23)*(x23))));
auto x25 = cell[8]*x0 - x37*(x13 - x14*(x1*x36 + x39*x9 + x40));
cell[0] = x6;
cell[1] = x12;
cell[2] = x15;
cell[3] = x18;
cell[4] = x19;
cell[5] = x21;
cell[6] = x22;
cell[7] = x24;
cell[8] = x25;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x12 + x15 + x18 + x19 + x21 + x22 + x24 + x25 + x6 + V{1};
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x35;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x36;
return { x5, x2 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
