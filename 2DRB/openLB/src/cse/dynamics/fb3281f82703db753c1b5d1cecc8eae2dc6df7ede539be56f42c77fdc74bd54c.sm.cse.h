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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::VelocityBoundaryDensity<1, -1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
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
auto x5 = -V{1}*x61;
auto x6 = x5 + V{1};
auto x67 = V{1} / (x6);
auto x7 = x67*(cell[0] + cell[10] + cell[12] + V{2}*cell[14] + cell[15] + cell[16] + cell[1] + V{2}*cell[2] + cell[3] + V{2}*cell[4] + cell[6] + cell[7] + V{2}*cell[8] + V{2}*cell[9] + V{1});
auto x8 = cell[0]*x0 - x4*x63*(x7*(x3 + V{-1}) + V{1});
auto x9 = -V{1}*x60;
auto x10 = V{0.5}*((x1)*(x1));
auto x11 = V{0.0555555555555556};
auto x12 = x11*x7;
auto x13 = cell[1]*x0 - x63*(x11 - x12*(x1*x9 + x10*((x9)*(x9)) + x68));
auto x14 = cell[2]*x0 - x63*(x11 - x12*(x1*x5 + x10*((x5)*(x5)) + x68));
auto x15 = -V{1}*x62;
auto x16 = cell[3]*x0 - x63*(x11 - x12*(x1*x15 + x10*((x15)*(x15)) + x68));
auto x17 = x5 + x9;
auto x18 = V{0.0277777777777778};
auto x19 = x18*x7;
auto x20 = cell[4]*x0 - x63*(x18 - x19*(x1*x17 + x10*((x17)*(x17)) + x68));
auto x21 = x61 + x9;
auto x22 = x15 + x9;
auto x23 = cell[6]*x0 - x63*(x18 - x19*(x1*x22 + x10*((x22)*(x22)) + x68));
auto x24 = x62 + x9;
auto x25 = cell[7]*x0 - x63*(x18 - x19*(x1*x24 + x10*((x24)*(x24)) + x68));
auto x26 = x15 + x5;
auto x27 = cell[8]*x0 - x63*(x18 - x19*(x1*x26 + x10*((x26)*(x26)) + x68));
auto x28 = x5 + x62;
auto x29 = cell[9]*x0 - x63*(x18 - x19*(x1*x28 + x10*((x28)*(x28)) + x68));
auto x30 = cell[10]*x0 - x63*(x11 - x12*(x1*x60 + x10*x64 + x68));
auto x31 = cell[12]*x0 - x63*(x11 - x12*(x1*x62 + x10*x66 + x68));
auto x32 = x60 + x61;
auto x33 = x5 + x60;
auto x34 = cell[14]*x0 - x63*(x18 - x19*(x1*x33 + x10*((x33)*(x33)) + x68));
auto x35 = x60 + x62;
auto x36 = cell[15]*x0 - x63*(x18 - x19*(x1*x35 + x10*((x35)*(x35)) + x68));
auto x37 = x15 + x60;
auto x38 = cell[16]*x0 - x63*(x18 - x19*(x1*x37 + x10*((x37)*(x37)) + x68));
auto x39 = x61 + x62;
auto x40 = x15 + x61;
cell[0] = x8;
cell[10] = x30;
cell[11] = cell[11]*x0 - x63*(x11 - x12*(x1*x61 + x10*x65 + x68));
cell[12] = x31;
cell[13] = cell[13]*x0 - x63*(x18 - x19*(x1*x32 + x10*((x32)*(x32)) + x68));
cell[14] = x34;
cell[15] = x36;
cell[16] = x38;
cell[17] = cell[17]*x0 - x63*(x18 - x19*(x1*x39 + x10*((x39)*(x39)) + x68));
cell[18] = cell[18]*x0 - x63*(x18 - x19*(x1*x40 + x10*((x40)*(x40)) + x68));
cell[1] = x13;
cell[2] = x14;
cell[3] = x16;
cell[4] = x20;
cell[5] = cell[5]*x0 - x63*(x18 - x19*(x1*x21 + x10*((x21)*(x21)) + x68));
cell[6] = x23;
cell[7] = x25;
cell[8] = x27;
cell[9] = x29;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x67*(x13 + V{2}*x14 + x16 + V{2}*x20 + x23 + x25 + V{2}*x27 + V{2}*x29 + x30 + x31 + V{2}*x34 + x36 + x38 + x8 + V{1});
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x60;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x61;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x62;
return { x7, x2 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
