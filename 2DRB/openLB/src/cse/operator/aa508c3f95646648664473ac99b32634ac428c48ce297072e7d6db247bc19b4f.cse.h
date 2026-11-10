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
struct CSE_O<operators::TrackAverageFunctorO<functors::DerivativeF<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, descriptors::POROSITY, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination> >, opti::AVERAGE_DJDALPHA<descriptors::POROSITY> >,descriptors::D3Q19<FIELDS...>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
void apply(CELL& cell, PARAMETERS& parameters) any_platform {
using V = typename CELL::value_t;
auto x63 = parameters.template get<parameters::PHYS_CHAR_VISCOSITY>();
auto x66 = parameters.template get<descriptors::LATTICE_TIME_T>();
auto x64 = parameters.template get<parameters::DT>();
auto x59 = parameters.template get<parameters::OMEGA>();
auto x60 = parameters.template get<parameters::DX>();
auto x61 = parameters.template get<parameters::LATTICE_VISCOSITY>();
auto x65 = parameters.template get<parameters::FACTOR>();
auto x62 = parameters.template get<parameters::CONVERSION_VELOCITY>();
auto x67 = cell[15] + cell[16];
auto x68 = cell[10] + cell[14] + x67;
auto x69 = cell[1] + cell[4] + x68;
auto x70 = cell[12] + cell[7] + cell[9];
auto x71 = cell[3] + cell[8] + x70;
auto x72 = cell[13] + cell[5];
auto x73 = cell[17] + cell[18];
auto x74 = cell[11] + x72 + x73;
auto x75 = cell[2] + x74;
auto x82 = ((x62)*(x62));
auto x0 = -V{1}*cell[16];
auto x76 = cell[15] + x0;
auto x1 = -V{1}*cell[7];
auto x2 = -V{1}*cell[6];
auto x3 = -V{1}*cell[8];
auto x4 = -V{1}*cell[18];
auto x77 = cell[17] + x4;
auto x5 = cell[0] + cell[6] + x69 + x71 + x75 + V{1};
auto x6 = V{1}/x5;
auto x7 = x6*(-V{1}*cell[3] + x2 + x3 + x70 + x76 + x77);
auto x8 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*x7;
auto x81 = x5*x8;
auto x9 = -V{1}*cell[4];
auto x10 = -V{1}*cell[5];
auto x78 = cell[13] + x10;
auto x11 = x6*(-V{1}*cell[1] + x1 + x2 + x68 + x78 + x9);
auto x12 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*x11;
auto x13 = x12*x5;
auto x14 = V{1}*x7;
auto x15 = V{1}*x11;
auto x16 = x15*x5;
auto x17 = -(x13*x14 + x16*x8)*(cell[6] + x1 - x13*x8 + x76);
auto x18 = V{2}*x17;
auto x19 = -V{1}*cell[9];
auto x20 = -V{1}*cell[14];
auto x21 = x6*(-V{1}*cell[2] + x19 + x20 + x3 + x74 + x9);
auto x22 = V{1}*x21;
auto x23 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*x21;
auto x24 = -x13*x22 - x16*x23;
auto x25 = cell[4] - x13*x23 + x20 + x78;
auto x26 = V{3};
auto x27 = V{1}*(x5 + V{-1})/x26;
auto x79 = -x27;
auto x80 = cell[6] + x79;
auto x28 = -(x12*x16 + x13*x15)*(cell[7] - x12*x13 + x69 + x72 + x80);
auto x29 = x81;
auto x30 = -x14*(x29 + x81)*(-x29*x8 + x67 + x71 + x73 + x80);
auto x31 = x23*x5;
auto x32 = x22*x5;
auto x33 = -x14*x31 - x32*x8;
auto x34 = cell[8] + x19 - x31*x8 + x77;
auto x35 = -(x22*x31 + x23*x32)*(cell[14] + cell[4] + cell[8] + cell[9] - x23*x31 + x75 + x79);
auto x36 = x12*x15;
auto x37 = x22*x23;
auto x38 = x14*x8;
auto x39 = V{1}*x59/(((x60)*(x60))*x61);
auto x40 = cell[0];
auto x41 = cell[1];
auto x42 = cell[2];
auto x43 = cell[3];
auto x44 = cell[4];
auto x45 = cell[5];
auto x46 = cell[6];
auto x47 = cell[7];
auto x48 = cell[8];
auto x49 = cell[9];
auto x50 = cell[10];
auto x51 = cell[11];
auto x52 = cell[12];
auto x53 = cell[13];
auto x54 = cell[14];
auto x55 = cell[15];
auto x56 = cell[16];
auto x57 = cell[17];
auto x58 = cell[18];
cell.template getFieldPointer<olb::opti::AVERAGE_DJDALPHA<olb::descriptors::POROSITY> >()[0] = (cell.template getFieldComponent<olb::opti::AVERAGE_DJDALPHA<olb::descriptors::POROSITY> >(0)*(x66 + V{-1}) - x63*x65*(V{2}*x39*x82*(cell.template getFieldComponent<olb::descriptors::POROSITY>(0) + V{-1})*(x36 + x37 + x38) + V{1}*x39*x82*(((x12)*(x12)) + ((x23)*(x23)) + ((x8)*(x8))) - V{0.5}*util::pow(x26*x59*x6/x64, V{2})*(V{2}*x18 + V{4}*x24*x25 + V{2}*x28 + V{2}*x30 + V{4}*x33*x34 + V{2}*x35)))/x66;
cell[0] = x40;
cell[1] = x41;
cell[2] = x42;
cell[3] = x43;
cell[4] = x44;
cell[5] = x45;
cell[6] = x46;
cell[7] = x47;
cell[8] = x48;
cell[9] = x49;
cell[10] = x50;
cell[11] = x51;
cell[12] = x52;
cell[13] = x53;
cell[14] = x54;
cell[15] = x55;
cell[16] = x56;
cell[17] = x57;
cell[18] = x58;
}
};

}

}

// Generation Info: commit=21521319fcfc23b85b6fd87340ebc70ad50346d4
