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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, 1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::ConstRhoBGK, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x44 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x45 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x46 = parameters.template get<descriptors::OMEGA>();
auto x47 = parameters.template get<statistics::AVERAGE_RHO>();
auto x48 = ((x44)*(x44));
auto x49 = ((x45)*(x45));
auto x0 = V{3};
auto x1 = x48 + x49;
auto x2 = V{0.5}*x0*x1;
auto x52 = V{1} - x2;
auto x3 = V{0.444444444444444};
auto x4 = x44 + V{1};
auto x50 = V{1} / (x4);
auto x5 = x50*(cell[0] + cell[4] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7] + cell[8] + V{1});
auto x6 = -x3*(x5*(x2 + V{-1}) + V{1});
auto x7 = V{1} - x46;
auto x8 = V{1} - x47;
auto x9 = V{1} + x8/x5;
auto x10 = -x3 + x7*(cell[0] - x6) + x9*(x3 + x6);
auto x11 = -V{1}*x44;
auto x12 = x11 + x45;
auto x13 = V{0.5}*((x0)*(x0));
auto x14 = V{0.0277777777777778};
auto x51 = -x14;
auto x15 = x14*x5;
auto x16 = x15*(x0*x12 + x13*((x12)*(x12)) + x52) + x51;
auto x17 = V{0.111111111111111};
auto x53 = -x17;
auto x18 = x17*x5;
auto x19 = x18*(x0*x11 + x13*((x11)*(x11)) + x52) + x53;
auto x20 = -V{1}*x45;
auto x21 = x11 + x20;
auto x22 = x15*(x0*x21 + x13*((x21)*(x21)) + x52) + x51;
auto x23 = x18*(x0*x20 + x13*((x20)*(x20)) + x52) + x53;
auto x24 = x53 + x7*(cell[4] - x23) + x9*(x17 + x23);
auto x25 = x20 + x44;
auto x26 = x15*(x0*x25 + x13*((x25)*(x25)) + x52) + x51;
auto x27 = x51 + x7*(cell[5] - x26) + x9*(x14 + x26);
auto x28 = x18*(x0*x44 + x13*x48 + x52) + x53;
auto x29 = x53 + x7*(cell[6] - x28) + x9*(x17 + x28);
auto x30 = x44 + x45;
auto x31 = x15*(x0*x30 + x13*((x30)*(x30)) + x52) + x51;
auto x32 = x51 + x7*(cell[7] - x31) + x9*(x14 + x31);
auto x33 = x18*(x0*x45 + x13*x49 + x52) + x53;
auto x34 = x53 + x7*(cell[8] - x33) + x9*(x17 + x33);
cell[0] = x10;
cell[1] = x51 + x7*(cell[1] - x16) + x9*(x14 + x16);
cell[2] = x53 + x7*(cell[2] - x19) + x9*(x17 + x19);
cell[3] = x51 + x7*(cell[3] - x22) + x9*(x14 + x22);
cell[4] = x24;
cell[5] = x27;
cell[6] = x29;
cell[7] = x32;
cell[8] = x34;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x50*(x10 + x24 + V{2}*x27 + V{2}*x29 + V{2}*x32 + x34 + V{1});
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x44;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x45;
return { x5 + x8, x1 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
