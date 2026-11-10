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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<1, 1>, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x32 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x33 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<1, 1>::VELOCITY>(0);
auto x35 = parameters.template get<descriptors::OMEGA>();
auto x34 = V{1} / (x32);
auto x37 = ((x33)*(x33));
auto x0 = V{1} - x35;
auto x1 = V{3};
auto x2 = x34*(cell[0] + V{2}*cell[1] + cell[2] + cell[6] + V{2}*cell[7] + V{2}*cell[8] + V{1}) + V{-1};
auto x36 = ((x2)*(x2));
auto x3 = x36 + x37;
auto x4 = V{0.5}*x1*x3;
auto x38 = V{1} - x4;
auto x5 = V{0.444444444444444};
auto x6 = cell[0]*x0 - x35*x5*(x32*(x4 + V{-1}) + V{1});
auto x7 = -V{1}*x33;
auto x8 = x2 + x7;
auto x9 = V{0.5}*((x1)*(x1));
auto x10 = V{0.0277777777777778};
auto x11 = x10*x32;
auto x12 = cell[1]*x0 - x35*(x10 - x11*(x1*x8 + x38 + x9*((x8)*(x8))));
auto x13 = V{0.111111111111111};
auto x14 = x13*x32;
auto x15 = cell[2]*x0 - x35*(x13 - x14*(x1*x7 + x38 + x9*((x7)*(x7))));
auto x16 = -V{1}*x2;
auto x17 = x16 + x7;
auto x18 = x16 + x33;
auto x19 = cell[6]*x0 - x35*(x13 - x14*(x1*x33 + x37*x9 + x38));
auto x20 = x2 + x33;
auto x21 = cell[7]*x0 - x35*(x10 - x11*(x1*x20 + x38 + x9*((x20)*(x20))));
auto x22 = cell[8]*x0 - x35*(x13 - x14*(x1*x2 + x36*x9 + x38));
cell[0] = x6;
cell[1] = x12;
cell[2] = x15;
cell[3] = cell[3]*x0 - x35*(x10 - x11*(x1*x17 + x38 + x9*((x17)*(x17))));
cell[4] = cell[4]*x0 - x35*(x13 - x14*(x1*x16 + x38 + x9*((x16)*(x16))));
cell[5] = cell[5]*x0 - x35*(x10 - x11*(x1*x18 + x38 + x9*((x18)*(x18))));
cell[6] = x19;
cell[7] = x21;
cell[8] = x22;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x32;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x33;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x34*(V{2}*x12 + x15 + x19 + V{2}*x21 + V{2}*x22 + x6 + V{1}) + V{-1};
return { x32, x3 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
