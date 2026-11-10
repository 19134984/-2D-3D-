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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, 1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x60 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x61 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x62 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x63 = parameters.template get<descriptors::OMEGA>();
auto x64 = ((x60)*(x60));
auto x65 = ((x61)*(x61));
auto x66 = ((x62)*(x62));
auto x0 = V{1} - x63;
auto x1 = V{3};
auto x2 = x64 + x65 + x66;
auto x3 = V{0.5}*x1*x2;
auto x68 = V{1} - x3;
auto x4 = V{0.333333333333333};
auto x5 = x60 + V{1};
auto x67 = V{1} / (x5);
auto x6 = x67*(cell[0] + V{2}*cell[10] + cell[11] + cell[12] + V{2}*cell[13] + V{2}*cell[14] + V{2}*cell[15] + V{2}*cell[16] + cell[17] + cell[18] + cell[2] + cell[3] + cell[8] + cell[9] + V{1});
auto x7 = cell[0]*x0 - x4*x63*(x6*(x3 + V{-1}) + V{1});
auto x8 = -V{1}*x60;
auto x9 = V{0.5}*((x1)*(x1));
auto x10 = V{0.0555555555555556};
auto x11 = x10*x6;
auto x12 = -V{1}*x61;
auto x13 = cell[2]*x0 - x63*(x10 - x11*(x1*x12 + x68 + x9*((x12)*(x12))));
auto x14 = -V{1}*x62;
auto x15 = cell[3]*x0 - x63*(x10 - x11*(x1*x14 + x68 + x9*((x14)*(x14))));
auto x16 = x12 + x8;
auto x17 = V{0.0277777777777778};
auto x18 = x17*x6;
auto x19 = x61 + x8;
auto x20 = x14 + x8;
auto x21 = x62 + x8;
auto x22 = x12 + x14;
auto x23 = cell[8]*x0 - x63*(x17 - x18*(x1*x22 + x68 + x9*((x22)*(x22))));
auto x24 = x12 + x62;
auto x25 = cell[9]*x0 - x63*(x17 - x18*(x1*x24 + x68 + x9*((x24)*(x24))));
auto x26 = cell[10]*x0 - x63*(x10 - x11*(x1*x60 + x64*x9 + x68));
auto x27 = cell[11]*x0 - x63*(x10 - x11*(x1*x61 + x65*x9 + x68));
auto x28 = cell[12]*x0 - x63*(x10 - x11*(x1*x62 + x66*x9 + x68));
auto x29 = x60 + x61;
auto x30 = cell[13]*x0 - x63*(x17 - x18*(x1*x29 + x68 + x9*((x29)*(x29))));
auto x31 = x12 + x60;
auto x32 = cell[14]*x0 - x63*(x17 - x18*(x1*x31 + x68 + x9*((x31)*(x31))));
auto x33 = x60 + x62;
auto x34 = cell[15]*x0 - x63*(x17 - x18*(x1*x33 + x68 + x9*((x33)*(x33))));
auto x35 = x14 + x60;
auto x36 = cell[16]*x0 - x63*(x17 - x18*(x1*x35 + x68 + x9*((x35)*(x35))));
auto x37 = x61 + x62;
auto x38 = cell[17]*x0 - x63*(x17 - x18*(x1*x37 + x68 + x9*((x37)*(x37))));
auto x39 = x14 + x61;
auto x40 = cell[18]*x0 - x63*(x17 - x18*(x1*x39 + x68 + x9*((x39)*(x39))));
cell[0] = x7;
cell[10] = x26;
cell[11] = x27;
cell[12] = x28;
cell[13] = x30;
cell[14] = x32;
cell[15] = x34;
cell[16] = x36;
cell[17] = x38;
cell[18] = x40;
cell[1] = cell[1]*x0 - x63*(x10 - x11*(x1*x8 + x68 + x9*((x8)*(x8))));
cell[2] = x13;
cell[3] = x15;
cell[4] = cell[4]*x0 - x63*(x17 - x18*(x1*x16 + x68 + x9*((x16)*(x16))));
cell[5] = cell[5]*x0 - x63*(x17 - x18*(x1*x19 + x68 + x9*((x19)*(x19))));
cell[6] = cell[6]*x0 - x63*(x17 - x18*(x1*x20 + x68 + x9*((x20)*(x20))));
cell[7] = cell[7]*x0 - x63*(x17 - x18*(x1*x21 + x68 + x9*((x21)*(x21))));
cell[8] = x23;
cell[9] = x25;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x67*(x13 + x15 + x23 + x25 + V{2}*x26 + x27 + x28 + V{2}*x30 + V{2}*x32 + V{2}*x34 + V{2}*x36 + x38 + x40 + x7 + V{1});
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x60;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x61;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x62;
return { x6, x2 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
