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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::VelocityBoundaryDensity<2, -1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x45 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x46 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x47 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x48 = parameters.template get<descriptors::OMEGA>();
auto x49 = ((x45)*(x45));
auto x50 = ((x46)*(x46));
auto x51 = ((x47)*(x47));
auto x0 = V{1} - x48;
auto x1 = V{3};
auto x2 = x49 + x50 + x51;
auto x3 = V{0.5}*x1*x2;
auto x53 = V{1} - x3;
auto x4 = V{0.333333333333333};
auto x5 = -V{1}*x47;
auto x6 = (cell[0] + cell[10] + cell[11] + cell[13] + cell[14] + V{2}*cell[16] + V{2}*cell[18] + cell[1] + cell[2] + V{2}*cell[3] + cell[4] + cell[5] + V{2}*cell[6] + V{2}*cell[8] + V{1})/(x5 + V{1});
auto x7 = -V{1}*x45;
auto x8 = V{0.5}*((x1)*(x1));
auto x9 = V{0.0555555555555556};
auto x52 = -x9;
auto x10 = x6*x9;
auto x11 = -V{1}*x46;
auto x12 = x11 + x7;
auto x13 = V{0.0277777777777778};
auto x14 = x13*x6;
auto x15 = x46 + x7;
auto x16 = x5 + x7;
auto x17 = x47 + x7;
auto x18 = x11 + x5;
auto x19 = x11 + x47;
auto x20 = x45 + x46;
auto x21 = x11 + x45;
auto x22 = x45 + x47;
auto x23 = x45 + x5;
auto x24 = x46 + x47;
auto x25 = x46 + x5;
auto x26 = cell[0]*x0 - x4*x48*(x6*(x3 + V{-1}) + V{1});
auto x27 = cell[10]*x0 + x48*(x10*(x1*x45 + x49*x8 + x53) + x52);
auto x28 = cell[11]*x0 + x48*(x10*(x1*x46 + x50*x8 + x53) + x52);
auto x29 = cell[12]*x0 + x48*(x10*(x1*x47 + x51*x8 + x53) + x52);
auto x30 = cell[13]*x0 - x48*(x13 - x14*(x1*x20 + x53 + x8*((x20)*(x20))));
auto x31 = cell[14]*x0 - x48*(x13 - x14*(x1*x21 + x53 + x8*((x21)*(x21))));
auto x32 = cell[15]*x0 - x48*(x13 - x14*(x1*x22 + x53 + x8*((x22)*(x22))));
auto x33 = cell[16]*x0 - x48*(x13 - x14*(x1*x23 + x53 + x8*((x23)*(x23))));
auto x34 = cell[17]*x0 - x48*(x13 - x14*(x1*x24 + x53 + x8*((x24)*(x24))));
auto x35 = cell[18]*x0 - x48*(x13 - x14*(x1*x25 + x53 + x8*((x25)*(x25))));
auto x36 = cell[1]*x0 + x48*(x10*(x1*x7 + x53 + x8*((x7)*(x7))) + x52);
auto x37 = cell[2]*x0 + x48*(x10*(x1*x11 + x53 + x8*((x11)*(x11))) + x52);
auto x38 = cell[3]*x0 + x48*(x10*(x1*x5 + x53 + x8*((x5)*(x5))) + x52);
auto x39 = cell[4]*x0 - x48*(x13 - x14*(x1*x12 + x53 + x8*((x12)*(x12))));
auto x40 = cell[5]*x0 - x48*(x13 - x14*(x1*x15 + x53 + x8*((x15)*(x15))));
auto x41 = cell[6]*x0 - x48*(x13 - x14*(x1*x16 + x53 + x8*((x16)*(x16))));
auto x42 = cell[7]*x0 - x48*(x13 - x14*(x1*x17 + x53 + x8*((x17)*(x17))));
auto x43 = cell[8]*x0 - x48*(x13 - x14*(x1*x18 + x53 + x8*((x18)*(x18))));
auto x44 = cell[9]*x0 - x48*(x13 - x14*(x1*x19 + x53 + x8*((x19)*(x19))));
cell[0] = x26;
cell[10] = x27;
cell[11] = x28;
cell[12] = x29;
cell[13] = x30;
cell[14] = x31;
cell[15] = x32;
cell[16] = x33;
cell[17] = x34;
cell[18] = x35;
cell[1] = x36;
cell[2] = x37;
cell[3] = x38;
cell[4] = x39;
cell[5] = x40;
cell[6] = x41;
cell[7] = x42;
cell[8] = x43;
cell[9] = x44;
return { x6, x2 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
