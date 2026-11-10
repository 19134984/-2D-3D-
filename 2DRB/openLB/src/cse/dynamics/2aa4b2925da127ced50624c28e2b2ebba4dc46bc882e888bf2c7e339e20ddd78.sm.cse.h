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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<2, 1>, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x59 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x60 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<2, 1>::VELOCITY>(0);
auto x61 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<2, 1>::VELOCITY>(1);
auto x63 = parameters.template get<descriptors::OMEGA>();
auto x62 = V{1} / (x59);
auto x65 = ((x60)*(x60));
auto x66 = ((x61)*(x61));
auto x0 = V{1} - x63;
auto x1 = V{3};
auto x2 = x62*(cell[0] + cell[10] + cell[11] + V{2}*cell[12] + cell[13] + cell[14] + V{2}*cell[15] + V{2}*cell[17] + cell[1] + cell[2] + cell[4] + cell[5] + V{2}*cell[7] + V{2}*cell[9] + V{1}) + V{-1};
auto x64 = ((x2)*(x2));
auto x3 = x64 + x65 + x66;
auto x4 = V{0.5}*x1*x3;
auto x68 = V{1} - x4;
auto x5 = V{0.333333333333333};
auto x6 = cell[0]*x0 - x5*x63*(x59*(x4 + V{-1}) + V{1});
auto x7 = -V{1}*x60;
auto x8 = V{0.5}*((x1)*(x1));
auto x9 = V{0.0555555555555556};
auto x67 = -x9;
auto x10 = x59*x9;
auto x11 = cell[1]*x0 + x63*(x10*(x1*x7 + x68 + x8*((x7)*(x7))) + x67);
auto x12 = -V{1}*x61;
auto x13 = cell[2]*x0 + x63*(x10*(x1*x12 + x68 + x8*((x12)*(x12))) + x67);
auto x14 = -V{1}*x2;
auto x15 = x12 + x7;
auto x16 = V{0.0277777777777778};
auto x17 = x16*x59;
auto x18 = cell[4]*x0 - x63*(x16 - x17*(x1*x15 + x68 + x8*((x15)*(x15))));
auto x19 = x61 + x7;
auto x20 = cell[5]*x0 - x63*(x16 - x17*(x1*x19 + x68 + x8*((x19)*(x19))));
auto x21 = x14 + x7;
auto x22 = x2 + x7;
auto x23 = cell[7]*x0 - x63*(x16 - x17*(x1*x22 + x68 + x8*((x22)*(x22))));
auto x24 = x12 + x14;
auto x25 = x12 + x2;
auto x26 = cell[9]*x0 - x63*(x16 - x17*(x1*x25 + x68 + x8*((x25)*(x25))));
auto x27 = cell[10]*x0 + x63*(x10*(x1*x60 + x65*x8 + x68) + x67);
auto x28 = cell[11]*x0 + x63*(x10*(x1*x61 + x66*x8 + x68) + x67);
auto x29 = cell[12]*x0 + x63*(x10*(x1*x2 + x64*x8 + x68) + x67);
auto x30 = x60 + x61;
auto x31 = cell[13]*x0 - x63*(x16 - x17*(x1*x30 + x68 + x8*((x30)*(x30))));
auto x32 = x12 + x60;
auto x33 = cell[14]*x0 - x63*(x16 - x17*(x1*x32 + x68 + x8*((x32)*(x32))));
auto x34 = x2 + x60;
auto x35 = cell[15]*x0 - x63*(x16 - x17*(x1*x34 + x68 + x8*((x34)*(x34))));
auto x36 = x14 + x60;
auto x37 = x2 + x61;
auto x38 = cell[17]*x0 - x63*(x16 - x17*(x1*x37 + x68 + x8*((x37)*(x37))));
auto x39 = x14 + x61;
cell[0] = x6;
cell[10] = x27;
cell[11] = x28;
cell[12] = x29;
cell[13] = x31;
cell[14] = x33;
cell[15] = x35;
cell[16] = cell[16]*x0 - x63*(x16 - x17*(x1*x36 + x68 + x8*((x36)*(x36))));
cell[17] = x38;
cell[18] = cell[18]*x0 - x63*(x16 - x17*(x1*x39 + x68 + x8*((x39)*(x39))));
cell[1] = x11;
cell[2] = x13;
cell[3] = cell[3]*x0 + x63*(x10*(x1*x14 + x68 + x8*((x14)*(x14))) + x67);
cell[4] = x18;
cell[5] = x20;
cell[6] = cell[6]*x0 - x63*(x16 - x17*(x1*x21 + x68 + x8*((x21)*(x21))));
cell[7] = x23;
cell[8] = cell[8]*x0 - x63*(x16 - x17*(x1*x24 + x68 + x8*((x24)*(x24))));
cell[9] = x26;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x59;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x60;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x61;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x62*(x11 + x13 + x18 + x20 + V{2}*x23 + V{2}*x26 + x27 + x28 + V{2}*x29 + x31 + x33 + V{2}*x35 + V{2}*x38 + x6 + V{1}) + V{-1};
return { x59, x3 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
