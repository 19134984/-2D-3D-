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
struct CSE<SourcedAdvectionDiffusionBGKdynamics<T, descriptors::D3Q19<FIELDS...> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x33 = cell.template getFieldComponent<olb::descriptors::SOURCE>(0);
auto x34 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x35 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x36 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x37 = parameters.template get<descriptors::OMEGA>();
auto x0 = V{1} - x37;
auto x1 = V{0.333333333333333};
auto x2 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{0.5}*x33 + V{1};
auto x3 = -x33*(V{0.5}*x37 + V{-1});
auto x4 = V{3};
auto x5 = -V{1}*x34;
auto x6 = V{0.0555555555555556};
auto x7 = x2*x6;
auto x8 = x3*x6;
auto x9 = -V{1}*x35;
auto x10 = -V{1}*x36;
auto x11 = V{0.0277777777777778};
auto x12 = x11*x2;
auto x13 = x11*x3;
auto x14 = cell[0]*x0 + x1*x3 + x1*x37*(x2 + V{-1});
auto x15 = cell[10]*x0 - x37*(x6 - x7*(x34*x4 + V{1})) + x8;
auto x16 = cell[11]*x0 - x37*(x6 - x7*(x35*x4 + V{1})) + x8;
auto x17 = cell[12]*x0 - x37*(x6 - x7*(x36*x4 + V{1})) + x8;
auto x18 = cell[13]*x0 + x13 - x37*(x11 - x12*(x4*(x34 + x35) + V{1}));
auto x19 = cell[14]*x0 + x13 - x37*(x11 - x12*(x4*(x34 + x9) + V{1}));
auto x20 = cell[15]*x0 + x13 - x37*(x11 - x12*(x4*(x34 + x36) + V{1}));
auto x21 = cell[16]*x0 + x13 - x37*(x11 - x12*(x4*(x10 + x34) + V{1}));
auto x22 = cell[17]*x0 + x13 - x37*(x11 - x12*(x4*(x35 + x36) + V{1}));
auto x23 = cell[18]*x0 + x13 - x37*(x11 - x12*(x4*(x10 + x35) + V{1}));
auto x24 = cell[1]*x0 - x37*(x6 - x7*(x4*x5 + V{1})) + x8;
auto x25 = cell[2]*x0 - x37*(x6 - x7*(x4*x9 + V{1})) + x8;
auto x26 = cell[3]*x0 - x37*(x6 - x7*(x10*x4 + V{1})) + x8;
auto x27 = cell[4]*x0 + x13 - x37*(x11 - x12*(x4*(x5 + x9) + V{1}));
auto x28 = cell[5]*x0 + x13 - x37*(x11 - x12*(x4*(x35 + x5) + V{1}));
auto x29 = cell[6]*x0 + x13 - x37*(x11 - x12*(x4*(x10 + x5) + V{1}));
auto x30 = cell[7]*x0 + x13 - x37*(x11 - x12*(x4*(x36 + x5) + V{1}));
auto x31 = cell[8]*x0 + x13 - x37*(x11 - x12*(x4*(x10 + x9) + V{1}));
auto x32 = cell[9]*x0 + x13 - x37*(x11 - x12*(x4*(x36 + x9) + V{1}));
cell[0] = x14;
cell[10] = x15;
cell[11] = x16;
cell[12] = x17;
cell[13] = x18;
cell[14] = x19;
cell[15] = x20;
cell[16] = x21;
cell[17] = x22;
cell[18] = x23;
cell[1] = x24;
cell[2] = x25;
cell[3] = x26;
cell[4] = x27;
cell[5] = x28;
cell[6] = x29;
cell[7] = x30;
cell[8] = x31;
cell[9] = x32;
return { x2, ((x34)*(x34)) + ((x35)*(x35)) + ((x36)*(x36)) };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
