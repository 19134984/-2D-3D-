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
struct CSE_O<operators::TrackAverageFunctorO<functors::DerivativeF<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, descriptors::POROSITY, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination> >, opti::AVERAGE_DJDALPHA<descriptors::POROSITY> >,descriptors::D2Q9<FIELDS...>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
void apply(CELL& cell, PARAMETERS& parameters) any_platform {
using V = typename CELL::value_t;
auto x37 = parameters.template get<parameters::DT>();
auto x38 = parameters.template get<parameters::FACTOR>();
auto x34 = parameters.template get<parameters::LATTICE_VISCOSITY>();
auto x36 = parameters.template get<parameters::PHYS_CHAR_VISCOSITY>();
auto x39 = parameters.template get<descriptors::LATTICE_TIME_T>();
auto x33 = parameters.template get<parameters::DX>();
auto x32 = parameters.template get<parameters::OMEGA>();
auto x35 = parameters.template get<parameters::CONVERSION_VELOCITY>();
auto x40 = cell[1] + cell[8];
auto x41 = cell[4] + x40;
auto x42 = cell[5] + cell[7];
auto x43 = cell[6] + x42;
auto x44 = cell[2] + x43;
auto x48 = ((x35)*(x35));
auto x0 = -V{1}*cell[3];
auto x1 = -V{1}*cell[5];
auto x45 = cell[7] + x1;
auto x2 = cell[0] + cell[3] + x41 + x44 + V{1};
auto x3 = V{1}/x2;
auto x4 = x3*(-V{1}*cell[4] + x0 + x40 + x45);
auto x5 = V{1}*x4;
auto x6 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*x4;
auto x46 = x2*x6;
auto x7 = x46;
auto x8 = V{3};
auto x9 = V{1}*(x2 + V{-1})/x8;
auto x47 = cell[3] - x9;
auto x10 = -x5*(x46 + x7)*(x41 + x42 + x47 - x6*x7);
auto x11 = -V{1}*cell[1];
auto x12 = x3*(-V{1}*cell[2] + x0 + x11 + x43);
auto x13 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*x12;
auto x14 = x13*x2;
auto x15 = V{1}*x12;
auto x16 = x15*x2;
auto x17 = -x14*x5 - x16*x6;
auto x18 = cell[3] + x11 - x14*x6 + x45;
auto x19 = -(x13*x16 + x14*x15)*(cell[1] - x13*x14 + x44 + x47);
auto x20 = x13*x15;
auto x21 = x5*x6;
auto x22 = V{1}*x32/(((x33)*(x33))*x34);
auto x23 = cell[0];
auto x24 = cell[1];
auto x25 = cell[2];
auto x26 = cell[3];
auto x27 = cell[4];
auto x28 = cell[5];
auto x29 = cell[6];
auto x30 = cell[7];
auto x31 = cell[8];
cell.template getFieldPointer<olb::opti::AVERAGE_DJDALPHA<olb::descriptors::POROSITY> >()[0] = (cell.template getFieldComponent<olb::opti::AVERAGE_DJDALPHA<olb::descriptors::POROSITY> >(0)*(x39 + V{-1}) - x36*x38*(V{2}*x22*x48*(cell.template getFieldComponent<olb::descriptors::POROSITY>(0) + V{-1})*(x20 + x21) + V{1}*x22*x48*(((x13)*(x13)) + ((x6)*(x6))) - V{0.5}*util::pow(x3*x32*x8/x37, V{2})*(V{2}*x10 + V{4}*x17*x18 + V{2}*x19)))/x39;
cell[0] = x23;
cell[1] = x24;
cell[2] = x25;
cell[3] = x26;
cell[4] = x27;
cell[5] = x28;
cell[6] = x29;
cell[7] = x30;
cell[8] = x31;
}
};

}

}

// Generation Info: commit=21521319fcfc23b85b6fd87340ebc70ad50346d4
