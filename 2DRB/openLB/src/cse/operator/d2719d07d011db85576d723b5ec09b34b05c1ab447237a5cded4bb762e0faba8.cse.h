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

template <typename T, typename...FIELDS>
struct CSE_O<operators::TrackAverageFunctorO<functors::DerivativeF<functors::FarFieldDragF<0>, descriptors::POROSITY, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, dynamics::DefaultCombination> >, opti::AVERAGE_DJDALPHA<descriptors::POROSITY> >,descriptors::D3Q19<FIELDS...>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
void apply(CELL& cell, PARAMETERS& parameters) any_platform {
using V = typename CELL::value_t;
auto x20 = parameters.template get<parameters::LATTICE_REFERENCE_VELOCITY>();
auto x22 = parameters.template get<parameters::FACTOR>();
auto x23 = parameters.template get<descriptors::LATTICE_TIME_T>();
auto x19 = cell[10] + cell[13] + cell[14] + cell[15] + cell[16];
auto x21 = cell[0] + cell[11] + cell[12] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + x19;
auto x24 = V{1} / (x21 + V{1});
auto x0 = cell[0];
auto x1 = cell[1];
auto x2 = cell[2];
auto x3 = cell[3];
auto x4 = cell[4];
auto x5 = cell[5];
auto x6 = cell[6];
auto x7 = cell[7];
auto x8 = cell[8];
auto x9 = cell[9];
auto x10 = cell[10];
auto x11 = cell[11];
auto x12 = cell[12];
auto x13 = cell[13];
auto x14 = cell[14];
auto x15 = cell[15];
auto x16 = cell[16];
auto x17 = cell[17];
auto x18 = cell[18];
cell.template getFieldPointer<olb::opti::AVERAGE_DJDALPHA<olb::descriptors::POROSITY> >()[0] = (cell.template getFieldComponent<olb::opti::AVERAGE_DJDALPHA<olb::descriptors::POROSITY> >(0)*(x23 + V{-1}) + x21*x22*x24*(V{1}*x20 - x24*(V{2}*cell[10] + V{2}*cell[13] + V{2}*cell[14] + V{2}*cell[15] + V{2}*cell[16] - V{2}*cell[1] - V{2}*cell[4] - V{2}*cell[5] - V{2}*cell[6] - V{2}*cell[7]))*(-cell[1] - cell[4] - cell[5] - cell[6] - cell[7] + x19))/x23;
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
cell[9] = x9;
cell[10] = x10;
cell[11] = x11;
cell[12] = x12;
cell[13] = x13;
cell[14] = x14;
cell[15] = x15;
cell[16] = x16;
cell[17] = x17;
cell[18] = x18;
}
};

}

}

// Generation Info: commit=9fa5ee2ead7964b5746886a513054fe377560274
