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

namespace operators {

template <typename T, typename... FIELDS>
struct CSE_O<StraightConvectionBoundaryProcessor3D<T, descriptors::D3Q19<FIELDS...>, 0, -1>,descriptors::D3Q19<FIELDS...>> {
template <concepts::Cell CELL>
void apply(CELL& cell) any_platform {
using V = typename CELL::value_t;
using DESCRIPTOR = typename CELL::descriptor_t;
V v_S0 = cell.computeRho();
V v_V0 [DESCRIPTOR::d]; cell.computeU(v_V0);
V v_S1 = cell.neighbor({1,0,0}).computeRho();
V v_V1 [DESCRIPTOR::d]; cell.neighbor({1,0,0}).computeU(v_V1);
V v_S2 = cell.neighbor({2,0,0}).computeRho();
V v_V2 [DESCRIPTOR::d]; cell.neighbor({2,0,0}).computeU(v_V2);
auto x0 = V{3};
auto x1 = V{3}*v_S0;
auto x2 = V{4}*v_S1;
auto x3 = -V{0.5}*v_S0*v_V0[0];
auto x4 = x3*(v_S2*v_V2[0] + v_V0[0]*x1 - v_V1[0]*x2);
auto x5 = V{0.0277777777777778}*x0;
auto x6 = x3*(v_S2*v_V2[1] + v_V0[1]*x1 - v_V1[1]*x2);
auto x7 = x3*(v_S2*v_V2[2] + v_V0[2]*x1 - v_V1[2]*x2);
auto x8 = cell[0];
auto x9 = cell[1];
auto x10 = cell[2];
auto x11 = cell[3];
auto x12 = cell[4];
auto x13 = cell[5];
auto x14 = cell[6];
auto x15 = cell[7];
auto x16 = cell[8];
auto x17 = cell[9];
auto x18 = cell.template getFieldComponent<olb::StraightConvectionBoundaryProcessor3D<olb::Expr, olb::descriptors::D3Q19<>, 0, -1>::PREV_CELL>(0);
auto x19 = cell[11];
auto x20 = cell[12];
auto x21 = cell.template getFieldComponent<olb::StraightConvectionBoundaryProcessor3D<olb::Expr, olb::descriptors::D3Q19<>, 0, -1>::PREV_CELL>(1);
auto x22 = cell.template getFieldComponent<olb::StraightConvectionBoundaryProcessor3D<olb::Expr, olb::descriptors::D3Q19<>, 0, -1>::PREV_CELL>(2);
auto x23 = cell.template getFieldComponent<olb::StraightConvectionBoundaryProcessor3D<olb::Expr, olb::descriptors::D3Q19<>, 0, -1>::PREV_CELL>(3);
auto x24 = cell.template getFieldComponent<olb::StraightConvectionBoundaryProcessor3D<olb::Expr, olb::descriptors::D3Q19<>, 0, -1>::PREV_CELL>(4);
auto x25 = cell[17];
auto x26 = cell[18];
cell.template getFieldPointer<olb::StraightConvectionBoundaryProcessor3D<olb::Expr, olb::descriptors::D3Q19<>, 0, -1>::PREV_CELL>()[0] = cell.template getFieldComponent<olb::StraightConvectionBoundaryProcessor3D<olb::Expr, olb::descriptors::D3Q19<>, 0, -1>::PREV_CELL>(0) + V{0.0555555555555556}*x0*x4;
cell.template getFieldPointer<olb::StraightConvectionBoundaryProcessor3D<olb::Expr, olb::descriptors::D3Q19<>, 0, -1>::PREV_CELL>()[1] = cell.template getFieldComponent<olb::StraightConvectionBoundaryProcessor3D<olb::Expr, olb::descriptors::D3Q19<>, 0, -1>::PREV_CELL>(1) + x5*(x4 + x6);
cell.template getFieldPointer<olb::StraightConvectionBoundaryProcessor3D<olb::Expr, olb::descriptors::D3Q19<>, 0, -1>::PREV_CELL>()[2] = cell.template getFieldComponent<olb::StraightConvectionBoundaryProcessor3D<olb::Expr, olb::descriptors::D3Q19<>, 0, -1>::PREV_CELL>(2) + x5*(x4 - V{1}*x6);
cell.template getFieldPointer<olb::StraightConvectionBoundaryProcessor3D<olb::Expr, olb::descriptors::D3Q19<>, 0, -1>::PREV_CELL>()[3] = cell.template getFieldComponent<olb::StraightConvectionBoundaryProcessor3D<olb::Expr, olb::descriptors::D3Q19<>, 0, -1>::PREV_CELL>(3) + x5*(x4 + x7);
cell.template getFieldPointer<olb::StraightConvectionBoundaryProcessor3D<olb::Expr, olb::descriptors::D3Q19<>, 0, -1>::PREV_CELL>()[4] = cell.template getFieldComponent<olb::StraightConvectionBoundaryProcessor3D<olb::Expr, olb::descriptors::D3Q19<>, 0, -1>::PREV_CELL>(4) + x5*(x4 - V{1}*x7);
cell[0] = x8;
cell[1] = x9;
cell[2] = x10;
cell[3] = x11;
cell[4] = x12;
cell[5] = x13;
cell[6] = x14;
cell[7] = x15;
cell[8] = x16;
cell[9] = x17;
cell[10] = x18;
cell[11] = x19;
cell[12] = x20;
cell[13] = x21;
cell[14] = x22;
cell[15] = x23;
cell[16] = x24;
cell[17] = x25;
cell[18] = x26;
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
