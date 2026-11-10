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
struct CSE<SourcedAdvectionDiffusionBGKdynamics<T, descriptors::D3Q19<FIELDS...> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x57 = cell.template getFieldComponent<olb::descriptors::SOURCE>(0);
auto x58 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x59 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x60 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x61 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{1} - x61;
auto x1 = V{0.333333333333333};
auto x2 = V{0.5}*x57;
auto x62 = x2 + V{1};
auto x3 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + x62;
auto x4 = -x57*(V{0.5}*x61 + V{-1});
auto x5 = cell[0]*x0 + x1*x4 + x1*x61*(x3 + V{-1});
auto x6 = V{3};
auto x7 = -V{1}*x58;
auto x8 = V{0.0555555555555556};
auto x9 = x3*x8;
auto x10 = x4*x8;
auto x11 = cell[1]*x0 + x10 - x61*(x8 - x9*(x6*x7 + V{1}));
auto x12 = -V{1}*x59;
auto x13 = cell[2]*x0 + x10 - x61*(x8 - x9*(x12*x6 + V{1}));
auto x14 = -V{1}*x60;
auto x15 = cell[3]*x0 + x10 - x61*(x8 - x9*(x14*x6 + V{1}));
auto x16 = V{0.0277777777777778};
auto x17 = x16*x3;
auto x18 = x16*x4;
auto x19 = cell[4]*x0 + x18 - x61*(x16 - x17*(x6*(x12 + x7) + V{1}));
auto x20 = cell[5]*x0 + x18 - x61*(x16 - x17*(x6*(x59 + x7) + V{1}));
auto x21 = cell[6]*x0 + x18 - x61*(x16 - x17*(x6*(x14 + x7) + V{1}));
auto x22 = cell[7]*x0 + x18 - x61*(x16 - x17*(x6*(x60 + x7) + V{1}));
auto x23 = cell[8]*x0 + x18 - x61*(x16 - x17*(x6*(x12 + x14) + V{1}));
auto x24 = cell[9]*x0 + x18 - x61*(x16 - x17*(x6*(x12 + x60) + V{1}));
auto x25 = cell[10]*x0 + x10 - x61*(x8 - x9*(x58*x6 + V{1}));
auto x26 = cell[11]*x0 + x10 - x61*(x8 - x9*(x59*x6 + V{1}));
auto x27 = cell[12]*x0 + x10 - x61*(x8 - x9*(x6*x60 + V{1}));
auto x65 = x22 + x24 + x27;
auto x28 = cell[13]*x0 + x18 - x61*(x16 - x17*(x6*(x58 + x59) + V{1}));
auto x29 = cell[14]*x0 + x18 - x61*(x16 - x17*(x6*(x12 + x58) + V{1}));
auto x30 = cell[15]*x0 + x18 - x61*(x16 - x17*(x6*(x58 + x60) + V{1}));
auto x31 = cell[16]*x0 + x18 - x61*(x16 - x17*(x6*(x14 + x58) + V{1}));
auto x63 = x25 + x29 + x31;
auto x32 = cell[17]*x0 + x18 - x61*(x16 - x17*(x6*(x59 + x60) + V{1}));
auto x33 = cell[18]*x0 + x18 - x61*(x16 - x17*(x6*(x14 + x59) + V{1}));
auto x64 = x20 + x26 + x33;
auto x34 = x11 + x13 + x15 + x19 + x21 + x23 + x28 + x30 + x32 + x5 + x62 + x63 + x64 + x65;
auto x66 = V{1} / (x34);
auto x35 = -V{1}*x19;
auto x67 = x28 + x35;
auto x36 = -V{1}*x21;
auto x68 = x30 + x36;
auto x37 = -V{1}*x23;
auto x69 = x32 + x37;
cell[0] = x5;
cell[10] = x25;
cell[11] = x26;
cell[12] = x27;
cell[13] = x28;
cell[14] = x29;
cell[15] = x30;
cell[16] = x31;
cell[17] = x32;
cell[18] = x33;
cell[1] = x11;
cell[2] = x13;
cell[3] = x15;
cell[4] = x19;
cell[5] = x20;
cell[6] = x21;
cell[7] = x22;
cell[8] = x23;
cell[9] = x24;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x34;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x66*(-V{1}*x11 - V{1}*x20 - V{1}*x22 + x63 + x67 + x68);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x66*(-V{1}*x13 - V{1}*x24 - V{1}*x29 + x64 + x67 + x69);
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x66*(-V{1}*x15 - V{1}*x31 - V{1}*x33 + x65 + x68 + x69);
return { x3, ((x58)*(x58)) + ((x59)*(x59)) + ((x60)*(x60)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
