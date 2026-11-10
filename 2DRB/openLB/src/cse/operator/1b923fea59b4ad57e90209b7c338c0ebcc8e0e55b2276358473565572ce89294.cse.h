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
struct CSE_O<OuterVelocityCornerProcessor2D<T, descriptors::D2Q9<FIELDS...>, 1, 1>,descriptors::D2Q9<FIELDS...>> {
template <concepts::Cell CELL>
void apply(CELL& cell) any_platform {
using V = typename CELL::value_t;
using DESCRIPTOR = typename CELL::descriptor_t;
auto x28 = cell.getOmegaOrFallback(std::numeric_limits<V>::signaling_NaN());
V v_S0 = cell.neighbor({-1,0}).computeRho();
V v_S1 = cell.neighbor({0,-1}).computeRho();
V v_S2 = cell.neighbor({-2,0}).computeRho();
V v_S3 = cell.neighbor({0,-2}).computeRho();
auto x29 = V{0.666666666666667}*v_S0 + V{0.666666666666667}*v_S1 - V{0.166666666666667}*v_S2 - V{0.166666666666667}*v_S3;
V v_V0 [DESCRIPTOR::d]; cell.computeU(v_V0);
V v_V1 [DESCRIPTOR::d]; cell.neighbor({-1,0}).computeU(v_V1);
V v_V2 [DESCRIPTOR::d]; cell.neighbor({-2,0}).computeU(v_V2);
V v_V3 [DESCRIPTOR::d]; cell.computeU(v_V3);
V v_V4 [DESCRIPTOR::d]; cell.neighbor({0,-1}).computeU(v_V4);
V v_V5 [DESCRIPTOR::d]; cell.neighbor({0,-2}).computeU(v_V5);
V v_V6 [DESCRIPTOR::d]; cell.computeU(v_V6);
V v_S4 = x29;
V v_P0 [DESCRIPTOR::q]; cell.computeEquilibrium(v_S4, v_V6, v_P0);
auto x0 = V{3};
auto x30 = V{1} / (x0);
auto x32 = ((x0)*(x0));
auto x1 = V{1}*x30;
auto x2 = -x1;
auto x3 = -x29*x30/x28;
auto x4 = V{-3};
auto x31 = V{1}*x4;
auto x33 = V{0.5}*x4;
auto x5 = -x3*(v_V0[0]*x31 + V{4}*v_V1[0] - V{1}*v_V2[0]);
auto x6 = x2*x5;
auto x7 = -x3*(v_V3[1]*x31 + V{4}*v_V4[1] - V{1}*v_V5[1]);
auto x8 = x2*x7;
auto x9 = V{0.0138888888888889}*x32;
auto x10 = V{1} - x1;
auto x11 = x10*x5;
auto x12 = -x3*(v_V0[1]*x33 + V{2}*v_V1[1] - V{0.5}*v_V2[1] + v_V3[0]*x33 + V{2}*v_V4[0] - V{0.5}*v_V5[0]);
auto x34 = V{2}*x12;
auto x13 = x10*x7;
auto x35 = x11 + x13;
auto x14 = x9*(-x34 + x35);
auto x15 = V{0.0555555555555556}*x32;
auto x16 = x15*(x11 + x8);
auto x17 = x9*(x34 + x35);
auto x18 = x15*(x13 + x6);
cell[0] = v_P0[0] + x32*(V{0.222222222222222}*x6 + V{0.222222222222222}*x8);
cell[1] = v_P0[1] + x14;
cell[2] = v_P0[2] + x16;
cell[3] = v_P0[3] + x17;
cell[4] = v_P0[4] + x18;
cell[5] = v_P0[5] + x14;
cell[6] = v_P0[6] + x16;
cell[7] = v_P0[7] + x17;
cell[8] = v_P0[8] + x18;
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
