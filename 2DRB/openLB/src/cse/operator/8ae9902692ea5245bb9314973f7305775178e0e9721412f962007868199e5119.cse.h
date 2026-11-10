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
struct CSE_O<operators::TrackAverageFunctorO<functors::DerivativeF<functors::FarFieldDragF<0>, descriptors::POROSITY, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::ThirdOrder, collision::FrozenOmegaFromCell<collision::ThirdOrderRLB>, dynamics::DefaultCombination> >, opti::AVERAGE_DJDALPHA<descriptors::POROSITY> >,descriptors::D3Q19<FIELDS...>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
void apply(CELL& cell, PARAMETERS& parameters) any_platform {
using V = typename CELL::value_t;
auto x27 = parameters.template get<descriptors::LATTICE_TIME_T>();
auto x23 = parameters.template get<parameters::CONVERSION_FORCE>();
auto x26 = parameters.template get<parameters::FACTOR>();
auto x24 = parameters.template get<parameters::LATTICE_REFERENCE_VELOCITY>();
auto x25 = cell[10] + cell[13] + cell[14] + cell[15] + cell[16];
auto x0 = cell[0] + cell[11] + cell[12] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + x25 + V{1};
auto x1 = V{1}*(-V{1}*cell[1] - V{1}*cell[4] - V{1}*cell[5] - V{1}*cell[6] - V{1}*cell[7] + x25)/x0;
auto x2 = V{1}*x1;
auto x3 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*x1*x2;
auto x4 = cell[0];
auto x5 = cell[1];
auto x6 = cell[2];
auto x7 = cell[3];
auto x8 = cell[4];
auto x9 = cell[5];
auto x10 = cell[6];
auto x11 = cell[7];
auto x12 = cell[8];
auto x13 = cell[9];
auto x14 = cell[10];
auto x15 = cell[11];
auto x16 = cell[12];
auto x17 = cell[13];
auto x18 = cell[14];
auto x19 = cell[15];
auto x20 = cell[16];
auto x21 = cell[17];
auto x22 = cell[18];
cell.template getFieldPointer<olb::opti::AVERAGE_DJDALPHA<olb::descriptors::POROSITY> >()[0] = (cell.template getFieldComponent<olb::opti::AVERAGE_DJDALPHA<olb::descriptors::POROSITY> >(0)*(x27 + V{-1}) + x23*x26*(x0 + V{-1})*(x2*x24 - V{2}*x3))/x27;
cell[0] = x4;
cell[1] = x5;
cell[2] = x6;
cell[3] = x7;
cell[4] = x8;
cell[5] = x9;
cell[6] = x10;
cell[7] = x11;
cell[8] = x12;
cell[9] = x13;
cell[10] = x14;
cell[11] = x15;
cell[12] = x16;
cell[13] = x17;
cell[14] = x18;
cell[15] = x19;
cell[16] = x20;
cell[17] = x21;
cell[18] = x22;
}
};

}

}

// Generation Info: commit=7ed6d8a1694b31c30b112c2d6e5b03f3c572d55c
