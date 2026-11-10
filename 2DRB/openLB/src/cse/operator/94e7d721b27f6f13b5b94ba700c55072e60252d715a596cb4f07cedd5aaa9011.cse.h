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
struct CSE_O<PlaneFdBoundaryProcessor3D<T, descriptors::D3Q19<FIELDS...>, 0, -1>,descriptors::D3Q19<FIELDS...>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
void apply(CELL& cell, PARAMETERS& parameters) any_platform {
using V = typename CELL::value_t;
using DESCRIPTOR = typename CELL::descriptor_t;
auto x53 = cell.getOmegaOrFallback(std::numeric_limits<V>::signaling_NaN());
V v_S0 = cell.computeRho();
V v_V0 [DESCRIPTOR::d]; cell.computeU(v_V0);
V v_V1 [DESCRIPTOR::d]; cell.computeU(v_V1);
V v_V2 [DESCRIPTOR::d]; cell.neighbor({1,0,0}).computeU(v_V2);
V v_V3 [DESCRIPTOR::d]; cell.neighbor({2,0,0}).computeU(v_V3);
V v_V4 [DESCRIPTOR::d]; cell.neighbor({0,1,0}).computeU(v_V4);
V v_V5 [DESCRIPTOR::d]; cell.neighbor({0,-1,0}).computeU(v_V5);
V v_V6 [DESCRIPTOR::d]; cell.neighbor({0,0,1}).computeU(v_V6);
V v_V7 [DESCRIPTOR::d]; cell.neighbor({0,0,-1}).computeU(v_V7);
V v_P0 [DESCRIPTOR::q]; cell.computeEquilibrium(v_S0, v_V0, v_P0);
auto x0 = V{3};
auto x52 = V{1}/x0;
auto x55 = ((x0)*(x0));
auto x1 = x52;
auto x2 = -x1;
auto x3 = -v_S0*x52/x53;
auto x54 = V{1}*x3;
auto x4 = V{-3};
auto x57 = V{0.5}*x4;
auto x5 = x3*(V{1}*v_V1[0]*x4 + V{4}*v_V2[0] - V{1}*v_V3[0]);
auto x6 = x2*x5;
auto x7 = x54*(v_V4[1] - v_V5[1]);
auto x8 = x2*x7;
auto x9 = x6 + x8;
auto x10 = x54*(v_V6[2] - v_V7[2]);
auto x11 = x10*x2;
auto x12 = V{0.0277777777777778}*x55;
auto x13 = V{1} - x1;
auto x14 = x13*x5;
auto x56 = x11 + x14;
auto x15 = x12*(x56 + x8);
auto x16 = x13*x7;
auto x58 = x16 + x56;
auto x17 = x16 + x6;
auto x18 = x12*(x11 + x17);
auto x19 = x10*x13;
auto x59 = x14 + x19 + x8;
auto x60 = x17 + x19;
auto x20 = x12*(x19 + x9);
auto x21 = V{0.0138888888888889}*x55;
auto x22 = V{2};
auto x23 = x3*(v_V1[1]*x57 + V{2}*v_V2[1] - V{0.5}*v_V3[1] + V{0.5}*v_V4[0] - V{0.5}*v_V5[0]);
auto x24 = x21*(x22*x23 + x58);
auto x25 = V{-2};
auto x26 = x21*(x23*x25 + x58);
auto x27 = x3*(v_V1[2]*x57 + V{2}*v_V2[2] - V{0.5}*v_V3[2] + V{0.5}*v_V6[0] - V{0.5}*v_V7[0]);
auto x28 = x21*(x22*x27 + x59);
auto x29 = x21*(x25*x27 + x59);
auto x30 = V{0.5}*x3*(v_V4[2] - v_V5[2] + v_V6[1] - v_V7[1]);
auto x31 = x21*(x22*x30 + x60);
auto x32 = x21*(x25*x30 + x60);
cell[0] = v_P0[0] + x55*(V{0.166666666666667}*x11 + V{0.166666666666667}*x9);
cell[1] = v_P0[1] + x15;
cell[2] = v_P0[2] + x18;
cell[3] = v_P0[3] + x20;
cell[4] = v_P0[4] + x24;
cell[5] = v_P0[5] + x26;
cell[6] = v_P0[6] + x28;
cell[7] = v_P0[7] + x29;
cell[8] = v_P0[8] + x31;
cell[9] = v_P0[9] + x32;
cell[10] = v_P0[10] + x15;
cell[11] = v_P0[11] + x18;
cell[12] = v_P0[12] + x20;
cell[13] = v_P0[13] + x24;
cell[14] = v_P0[14] + x26;
cell[15] = v_P0[15] + x28;
cell[16] = v_P0[16] + x29;
cell[17] = v_P0[17] + x31;
cell[18] = v_P0[18] + x32;
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
