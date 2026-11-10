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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<0, 1>, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x59 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x61 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1);
auto x62 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<0, 1>::VELOCITY>(2);
auto x63 = parameters.template get<descriptors::OMEGA>();
auto x60 = V{1} / (x59);
auto x65 = ((x61)*(x61));
auto x66 = ((x62)*(x62));
auto x0 = V{1} - x63;
auto x1 = V{3};
auto x2 = x60*(cell[0] + V{2}*cell[10] + cell[11] + cell[12] + V{2}*cell[13] + V{2}*cell[14] + V{2}*cell[15] + V{2}*cell[16] + cell[17] + cell[18] + cell[2] + cell[3] + cell[8] + cell[9] + V{1}) + V{-1};
auto x64 = ((x2)*(x2));
auto x3 = x64 + x65 + x66;
auto x4 = V{0.5}*x1*x3;
auto x68 = V{1} - x4;
auto x5 = V{0.333333333333333};
auto x6 = cell[0]*x0 - x5*x63*(x59*(x4 + V{-1}) + V{1});
auto x7 = -V{1}*x2;
auto x8 = V{0.5}*((x1)*(x1));
auto x9 = V{0.0555555555555556};
auto x67 = -x9;
auto x10 = x59*x9;
auto x11 = -V{1}*x61;
auto x12 = cell[2]*x0 + x63*(x10*(x1*x11 + x68 + x8*((x11)*(x11))) + x67);
auto x13 = -V{1}*x62;
auto x14 = cell[3]*x0 + x63*(x10*(x1*x13 + x68 + x8*((x13)*(x13))) + x67);
auto x15 = x11 + x7;
auto x16 = V{0.0277777777777778};
auto x17 = x16*x59;
auto x18 = x61 + x7;
auto x19 = x13 + x7;
auto x20 = x62 + x7;
auto x21 = x11 + x13;
auto x22 = cell[8]*x0 - x63*(x16 - x17*(x1*x21 + x68 + x8*((x21)*(x21))));
auto x23 = x11 + x62;
auto x24 = cell[9]*x0 - x63*(x16 - x17*(x1*x23 + x68 + x8*((x23)*(x23))));
auto x25 = cell[10]*x0 + x63*(x10*(x1*x2 + x64*x8 + x68) + x67);
auto x26 = cell[11]*x0 + x63*(x10*(x1*x61 + x65*x8 + x68) + x67);
auto x27 = cell[12]*x0 + x63*(x10*(x1*x62 + x66*x8 + x68) + x67);
auto x28 = x2 + x61;
auto x29 = cell[13]*x0 - x63*(x16 - x17*(x1*x28 + x68 + x8*((x28)*(x28))));
auto x30 = x11 + x2;
auto x31 = cell[14]*x0 - x63*(x16 - x17*(x1*x30 + x68 + x8*((x30)*(x30))));
auto x32 = x2 + x62;
auto x33 = cell[15]*x0 - x63*(x16 - x17*(x1*x32 + x68 + x8*((x32)*(x32))));
auto x34 = x13 + x2;
auto x35 = cell[16]*x0 - x63*(x16 - x17*(x1*x34 + x68 + x8*((x34)*(x34))));
auto x36 = x61 + x62;
auto x37 = cell[17]*x0 - x63*(x16 - x17*(x1*x36 + x68 + x8*((x36)*(x36))));
auto x38 = x13 + x61;
auto x39 = cell[18]*x0 - x63*(x16 - x17*(x1*x38 + x68 + x8*((x38)*(x38))));
cell[0] = x6;
cell[10] = x25;
cell[11] = x26;
cell[12] = x27;
cell[13] = x29;
cell[14] = x31;
cell[15] = x33;
cell[16] = x35;
cell[17] = x37;
cell[18] = x39;
cell[1] = cell[1]*x0 + x63*(x10*(x1*x7 + x68 + x8*((x7)*(x7))) + x67);
cell[2] = x12;
cell[3] = x14;
cell[4] = cell[4]*x0 - x63*(x16 - x17*(x1*x15 + x68 + x8*((x15)*(x15))));
cell[5] = cell[5]*x0 - x63*(x16 - x17*(x1*x18 + x68 + x8*((x18)*(x18))));
cell[6] = cell[6]*x0 - x63*(x16 - x17*(x1*x19 + x68 + x8*((x19)*(x19))));
cell[7] = cell[7]*x0 - x63*(x16 - x17*(x1*x20 + x68 + x8*((x20)*(x20))));
cell[8] = x22;
cell[9] = x24;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x59;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x60*(x12 + x14 + x22 + x24 + V{2}*x25 + x26 + x27 + V{2}*x29 + V{2}*x31 + V{2}*x33 + V{2}*x35 + x37 + x39 + x6 + V{1}) + V{-1};
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x61;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x62;
return { x59, x3 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
