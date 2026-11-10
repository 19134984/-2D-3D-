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
struct CSE_O<operators::TrackAverageFunctorO<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, opti::AVERAGE_J>,descriptors::D3Q19<FIELDS...>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
void apply(CELL& cell, PARAMETERS& parameters) any_platform {
using V = typename CELL::value_t;
using DESCRIPTOR = typename CELL::descriptor_t;
auto x22 = parameters.template get<parameters::LATTICE_VISCOSITY>();
auto x21 = parameters.template get<parameters::DX>();
auto x26 = parameters.template get<descriptors::LATTICE_TIME_T>();
auto x25 = parameters.template get<parameters::DT>();
auto x24 = parameters.template get<parameters::PHYS_CHAR_VISCOSITY>();
auto x23 = parameters.template get<parameters::CONVERSION_VELOCITY>();
auto x20 = parameters.template get<parameters::OMEGA>();
V v_S0 = cell.computeRho();
V v_V0 [DESCRIPTOR::d]; cell.computeU(v_V0);
V v_T0 [util::TensorVal<DESCRIPTOR>::n]; cell.computeAllMomenta(v_S0, v_V0, v_T0);
V v_V1 [DESCRIPTOR::d]; cell.computeU(v_V1);
auto x0 = ((v_T0[2])*(v_T0[2]));
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
cell.template getFieldPointer<olb::opti::AVERAGE_J>()[0] = (cell.template getFieldComponent<olb::opti::AVERAGE_J>(0)*(x26 + V{-1}) - x20*x24*(cell.template getFieldComponent<olb::descriptors::POROSITY>(0) + V{-1})*((x23)*(x23))*(((v_V1[0])*(v_V1[0])) + ((v_V1[1])*(v_V1[1])) + ((v_V1[2])*(v_V1[2])))/(((x21)*(x21))*x22) + V{4.5}*x24*util::pow(x20/(v_S0*x25), V{2})*(V{2}*x0 + ((v_T0[0])*(v_T0[0])) + V{2}*((v_T0[1])*(v_T0[1])) + ((v_T0[3])*(v_T0[3])) + V{2}*((v_T0[4])*(v_T0[4])) + ((v_T0[5])*(v_T0[5]))))/x26;
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

// Generation Info: commit=21521319fcfc23b85b6fd87340ebc70ad50346d4
