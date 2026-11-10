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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, -1>, momenta::FixedVelocityMomentumGeneric, momenta::BulkStress, momenta::DefineSeparately>, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
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
auto x5 = -V{1}*x60;
auto x6 = x5 + V{1};
auto x67 = V{1} / (x6);
auto x7 = x67*(cell[0] + cell[11] + cell[12] + cell[17] + cell[18] + V{2}*cell[1] + cell[2] + cell[3] + V{2}*cell[4] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7] + cell[8] + cell[9] + V{1});
auto x8 = cell[0]*x0 - x4*x63*(x7*(x3 + V{-1}) + V{1});
auto x9 = V{0.5}*((x1)*(x1));
auto x10 = V{0.0555555555555556};
auto x11 = x10*x7;
auto x12 = cell[1]*x0 - x63*(x10 - x11*(x1*x5 + x68 + x9*((x5)*(x5))));
auto x13 = -V{1}*x61;
auto x14 = cell[2]*x0 - x63*(x10 - x11*(x1*x13 + x68 + x9*((x13)*(x13))));
auto x15 = -V{1}*x62;
auto x16 = cell[3]*x0 - x63*(x10 - x11*(x1*x15 + x68 + x9*((x15)*(x15))));
auto x17 = x13 + x5;
auto x18 = V{0.0277777777777778};
auto x19 = x18*x7;
auto x20 = cell[4]*x0 - x63*(x18 - x19*(x1*x17 + x68 + x9*((x17)*(x17))));
auto x21 = x5 + x61;
auto x22 = cell[5]*x0 - x63*(x18 - x19*(x1*x21 + x68 + x9*((x21)*(x21))));
auto x23 = x15 + x5;
auto x24 = cell[6]*x0 - x63*(x18 - x19*(x1*x23 + x68 + x9*((x23)*(x23))));
auto x25 = x5 + x62;
auto x26 = cell[7]*x0 - x63*(x18 - x19*(x1*x25 + x68 + x9*((x25)*(x25))));
auto x27 = x13 + x15;
auto x28 = cell[8]*x0 - x63*(x18 - x19*(x1*x27 + x68 + x9*((x27)*(x27))));
auto x29 = x13 + x62;
auto x30 = cell[9]*x0 - x63*(x18 - x19*(x1*x29 + x68 + x9*((x29)*(x29))));
auto x31 = cell[11]*x0 - x63*(x10 - x11*(x1*x61 + x65*x9 + x68));
auto x32 = cell[12]*x0 - x63*(x10 - x11*(x1*x62 + x66*x9 + x68));
auto x33 = x60 + x61;
auto x34 = x13 + x60;
auto x35 = x60 + x62;
auto x36 = x15 + x60;
auto x37 = x61 + x62;
auto x38 = cell[17]*x0 - x63*(x18 - x19*(x1*x37 + x68 + x9*((x37)*(x37))));
auto x39 = x15 + x61;
auto x40 = cell[18]*x0 - x63*(x18 - x19*(x1*x39 + x68 + x9*((x39)*(x39))));
cell[0] = x8;
cell[10] = cell[10]*x0 - x63*(x10 - x11*(x1*x60 + x64*x9 + x68));
cell[11] = x31;
cell[12] = x32;
cell[13] = cell[13]*x0 - x63*(x18 - x19*(x1*x33 + x68 + x9*((x33)*(x33))));
cell[14] = cell[14]*x0 - x63*(x18 - x19*(x1*x34 + x68 + x9*((x34)*(x34))));
cell[15] = cell[15]*x0 - x63*(x18 - x19*(x1*x35 + x68 + x9*((x35)*(x35))));
cell[16] = cell[16]*x0 - x63*(x18 - x19*(x1*x36 + x68 + x9*((x36)*(x36))));
cell[17] = x38;
cell[18] = x40;
cell[1] = x12;
cell[2] = x14;
cell[3] = x16;
cell[4] = x20;
cell[5] = x22;
cell[6] = x24;
cell[7] = x26;
cell[8] = x28;
cell[9] = x30;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x67*(V{2}*x12 + x14 + x16 + V{2}*x20 + V{2}*x22 + V{2}*x24 + V{2}*x26 + x28 + x30 + x31 + x32 + x38 + x40 + x8 + V{1});
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x60;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x61;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x62;
return { x7, x2 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
