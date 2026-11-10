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

template <typename... FIELDS>
struct CSE_O<operators::TrackAverageFunctorO<functors::FarFieldDragF<0>, opti::AVERAGE_J>,descriptors::D3Q19<FIELDS...>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
void apply(CELL& cell, PARAMETERS& parameters) any_platform {
using V = typename CELL::value_t;
using DESCRIPTOR = typename CELL::descriptor_t;
auto x22 = parameters.template get<parameters::LATTICE_REFERENCE_PRESSURE>();
auto x23 = parameters.template get<descriptors::LATTICE_TIME_T>();
auto x20 = parameters.template get<parameters::CONVERSION_FORCE>();
auto x21 = parameters.template get<parameters::LATTICE_REFERENCE_VELOCITY>();
V v_S0 = cell.computeRho();
V v_V0 [DESCRIPTOR::d]; cell.computeU(v_V0);
auto x0 = v_S0 + V{-1};
auto x1 = cell[0];
auto x2 = cell[1];
auto x3 = cell[2];
auto x4 = cell[3];
auto x5 = cell[4];
auto x6 = cell[5];
auto x7 = cell[6];
auto x8 = cell[7];
auto x9 = cell[8];
auto x10 = cell[9];
auto x11 = cell[10];
auto x12 = cell[11];
auto x13 = cell[12];
auto x14 = cell[13];
auto x15 = cell[14];
auto x16 = cell[15];
auto x17 = cell[16];
auto x18 = cell[17];
auto x19 = cell[18];
cell.template getFieldPointer<olb::opti::AVERAGE_J>()[0] = (cell.template getFieldComponent<olb::opti::AVERAGE_J>(0)*(x23 + V{-1}) - x20*(v_V0[0]*x0*(v_V0[0] - x21) + V{0.333333333333333}*x0 - x22))/x23;
cell[0] = x1;
cell[1] = x2;
cell[2] = x3;
cell[3] = x4;
cell[4] = x5;
cell[5] = x6;
cell[6] = x7;
cell[7] = x8;
cell[8] = x9;
cell[9] = x10;
cell[10] = x11;
cell[11] = x12;
cell[12] = x13;
cell[13] = x14;
cell[14] = x15;
cell[15] = x16;
cell[16] = x17;
cell[17] = x18;
cell[18] = x19;
}
};

}

}

// Generation Info: commit=6decf22927dac6a87a8720c954e63974d643a84b
