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
struct CSE_O<operators::TrackAverageFunctorO<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, opti::AVERAGE_J>,descriptors::D2Q9<FIELDS...>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
void apply(CELL& cell, PARAMETERS& parameters) any_platform {
using V = typename CELL::value_t;
using DESCRIPTOR = typename CELL::descriptor_t;
auto x12 = parameters.template get<parameters::CONVERSION_VELOCITY>();
auto x9 = parameters.template get<parameters::OMEGA>();
auto x10 = parameters.template get<parameters::DX>();
auto x15 = parameters.template get<descriptors::LATTICE_TIME_T>();
auto x11 = parameters.template get<parameters::LATTICE_VISCOSITY>();
auto x13 = parameters.template get<parameters::PHYS_CHAR_VISCOSITY>();
auto x14 = parameters.template get<parameters::DT>();
V v_S0 = cell.computeRho();
V v_V0 [DESCRIPTOR::d]; cell.computeU(v_V0);
V v_T0 [util::TensorVal<DESCRIPTOR>::n]; cell.computeAllMomenta(v_S0, v_V0, v_T0);
V v_V1 [DESCRIPTOR::d]; cell.computeU(v_V1);
auto x0 = cell[0];
auto x1 = cell[1];
auto x2 = cell[2];
auto x3 = cell[3];
auto x4 = cell[4];
auto x5 = cell[5];
auto x6 = cell[6];
auto x7 = cell[7];
auto x8 = cell[8];
cell.template getFieldPointer<olb::opti::AVERAGE_J>()[0] = (cell.template getFieldComponent<olb::opti::AVERAGE_J>(0)*(x15 + V{-1}) + V{4.5}*x13*util::pow(x9/(v_S0*x14), V{2})*(((v_T0[0])*(v_T0[0])) + V{2}*((v_T0[1])*(v_T0[1])) + ((v_T0[2])*(v_T0[2]))) - x13*x9*(cell.template getFieldComponent<olb::descriptors::POROSITY>(0) + V{-1})*((x12)*(x12))*(((v_V1[0])*(v_V1[0])) + ((v_V1[1])*(v_V1[1])))/(((x10)*(x10))*x11))/x15;
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
}
};

}

}

// Generation Info: commit=b0dc2eff70ace778d8c09131cb5af79150b5dff9
