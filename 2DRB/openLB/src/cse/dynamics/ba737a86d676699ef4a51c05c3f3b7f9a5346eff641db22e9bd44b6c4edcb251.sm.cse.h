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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, -1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x33 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x34 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x35 = parameters.template get<descriptors::OMEGA>();
auto x36 = ((x33)*(x33));
auto x37 = ((x34)*(x34));
auto x0 = V{1} - x35;
auto x1 = V{3};
auto x2 = x36 + x37;
auto x3 = V{0.5}*x1*x2;
auto x39 = V{1} - x3;
auto x4 = V{0.444444444444444};
auto x5 = -V{1}*x33;
auto x6 = x5 + V{1};
auto x38 = V{1} / (x6);
auto x7 = x38*(cell[0] + V{2}*cell[1] + V{2}*cell[2] + V{2}*cell[3] + cell[4] + cell[8] + V{1});
auto x8 = cell[0]*x0 - x35*x4*(x7*(x3 + V{-1}) + V{1});
auto x9 = x34 + x5;
auto x10 = V{0.5}*((x1)*(x1));
auto x11 = V{0.0277777777777778};
auto x12 = x11*x7;
auto x13 = cell[1]*x0 - x35*(x11 - x12*(x1*x9 + x10*((x9)*(x9)) + x39));
auto x14 = V{0.111111111111111};
auto x15 = x14*x7;
auto x16 = cell[2]*x0 - x35*(x14 - x15*(x1*x5 + x10*((x5)*(x5)) + x39));
auto x17 = -V{1}*x34;
auto x18 = x17 + x5;
auto x19 = cell[3]*x0 - x35*(x11 - x12*(x1*x18 + x10*((x18)*(x18)) + x39));
auto x20 = cell[4]*x0 - x35*(x14 - x15*(x1*x17 + x10*((x17)*(x17)) + x39));
auto x21 = x17 + x33;
auto x22 = x33 + x34;
auto x23 = cell[8]*x0 - x35*(x14 - x15*(x1*x34 + x10*x37 + x39));
cell[0] = x8;
cell[1] = x13;
cell[2] = x16;
cell[3] = x19;
cell[4] = x20;
cell[5] = cell[5]*x0 - x35*(x11 - x12*(x1*x21 + x10*((x21)*(x21)) + x39));
cell[6] = cell[6]*x0 - x35*(x14 - x15*(x1*x33 + x10*x36 + x39));
cell[7] = cell[7]*x0 - x35*(x11 - x12*(x1*x22 + x10*((x22)*(x22)) + x39));
cell[8] = x23;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x38*(V{2}*x13 + V{2}*x16 + V{2}*x19 + x20 + x23 + x8 + V{1});
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x33;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x34;
return { x7, x2 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
