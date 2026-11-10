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
struct CSE_O<operators::TrackAverageFunctorO<functors::DerivativeF<functors::DissipationF, descriptors::POROSITY, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, dynamics::DefaultCombination> >, opti::AVERAGE_DJDALPHA<descriptors::POROSITY> >,descriptors::D3Q19<FIELDS...>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
void apply(CELL& cell, PARAMETERS& parameters) any_platform {
using V = typename CELL::value_t;
auto x20 = parameters.template get<parameters::DT>();
auto x22 = parameters.template get<parameters::FACTOR>();
auto x19 = parameters.template get<parameters::OMEGA>();
auto x23 = parameters.template get<descriptors::LATTICE_TIME_T>();
auto x21 = parameters.template get<parameters::PHYS_CHAR_VISCOSITY>();
auto x24 = cell[10] + cell[13] + cell[15];
auto x25 = cell[11] + cell[17] + cell[5];
auto x26 = cell[12] + cell[7] + cell[9];
auto x27 = cell[14] + cell[16];
auto x28 = V{1} / (cell[0] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + x24 + x25 + x26 + x27 + V{1});
auto x29 = -V{1}*cell[4];
auto x30 = cell[14] + x29;
auto x31 = -V{1}*cell[6];
auto x32 = cell[16] + x31;
auto x33 = -cell[5] - cell[7] + x24;
auto x34 = -V{1}*cell[1] + x30 + x32 + x33;
auto x35 = ((x34)*(x34));
auto x36 = V{1}*x28;
auto x37 = V{0.333333}*cell[0];
auto x38 = V{0.333333}*cell[12] - V{0.666667}*cell[13] - V{0.666667}*cell[14] + V{0.333333}*cell[3] - V{0.666667}*cell[4] - V{0.666667}*cell[5] + x37;
auto x39 = V{0.333333}*cell[11] - V{0.666667}*cell[15] - V{0.666667}*cell[16] + V{0.333333}*cell[2] - V{0.666667}*cell[6] - V{0.666667}*cell[7];
auto x40 = -V{1}*cell[8];
auto x41 = cell[18] + x40;
auto x42 = cell[13] - cell[14] - cell[2] - cell[9] + x25 + x29 + x41;
auto x43 = ((x42)*(x42));
auto x44 = V{0.333333}*cell[10] - V{0.666667}*cell[17] - V{0.666667}*cell[18] + V{0.333333}*cell[1] - V{0.666667}*cell[8] - V{0.666667}*cell[9];
auto x45 = cell[15] - cell[16] + cell[17] - cell[18] - cell[3] + x26 + x31 + x40;
auto x46 = ((x45)*(x45));
auto x47 = x34*x36;
auto x48 = -V{8}*cell[1] - V{8}*cell[4] - V{8}*cell[6] + V{8}*x27 + V{8}*x33;
auto x49 = x42*x45;
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
cell.template getFieldPointer<olb::opti::AVERAGE_DJDALPHA<olb::descriptors::POROSITY> >()[0] = (cell.template getFieldComponent<olb::opti::AVERAGE_DJDALPHA<olb::descriptors::POROSITY> >(0)*(x23 + V{-1}) + V{4.5}*x21*x22*x28*util::pow(x19*x28/x20, V{2})*(V{4}*x35*(-V{0.666667}*cell[10] + V{0.333333}*cell[17] + V{0.333333}*cell[18] - V{0.666667}*cell[1] + V{0.333333}*cell[8] + V{0.333333}*cell[9] + x35*x36 + x38 + x39) + x42*x48*(-V{1}*cell[13] + cell[5] + x30 + x42*x47) + V{4}*x43*(-V{0.666667}*cell[11] + V{0.333333}*cell[15] + V{0.333333}*cell[16] - V{0.666667}*cell[2] + V{0.333333}*cell[6] + V{0.333333}*cell[7] + x36*x43 + x38 + x44) + x45*x48*(-V{1}*cell[15] + cell[7] + x32 + x45*x47) + V{4}*x46*(-V{0.666667}*cell[12] + V{0.333333}*cell[13] + V{0.333333}*cell[14] - V{0.666667}*cell[3] + V{0.333333}*cell[4] + V{0.333333}*cell[5] + x36*x46 + x37 + x39 + x44) + V{8}*x49*(-V{1}*cell[17] + cell[9] + x36*x49 + x41)))/x23;
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

// Generation Info: commit=5690cdacbcabc20c4032537b2189a882dc41a9bd
