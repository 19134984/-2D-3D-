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
struct CSE_O<operators::TrackAverageFunctorO<functors::ReverseDerivativeF<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, descriptors::POPULATION, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination> >, opti::AVERAGE_DJDF>,descriptors::D2Q9<FIELDS...>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
void apply(CELL& cell, PARAMETERS& parameters) any_platform {
using V = typename CELL::value_t;
auto x69 = parameters.template get<parameters::DT>();
auto x66 = parameters.template get<parameters::LATTICE_VISCOSITY>();
auto x65 = parameters.template get<parameters::DX>();
auto x67 = parameters.template get<parameters::CONVERSION_VELOCITY>();
auto x68 = parameters.template get<parameters::PHYS_CHAR_VISCOSITY>();
auto x64 = parameters.template get<parameters::OMEGA>();
auto x71 = parameters.template get<descriptors::LATTICE_TIME_T>();
auto x70 = parameters.template get<parameters::FACTOR>();
auto x72 = cell[1] + cell[8];
auto x73 = cell[4] + x72;
auto x74 = cell[5] + cell[7];
auto x75 = cell[6] + x74;
auto x76 = cell[2] + x75;
auto x85 = V{2}*cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x86 = V{1} / (x71);
auto x0 = V{3};
auto x1 = x0*x64;
auto x78 = x1/x69;
auto x2 = cell[0] + cell[3] + x73 + x76 + V{1};
auto x77 = V{1} / (x2);
auto x80 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*x77;
auto x83 = V{1} / ((x2)*(x2));
auto x84 = V{1}*x83;
auto x3 = x77*x78;
auto x4 = -V{1}*cell[3];
auto x5 = -V{1}*cell[5];
auto x79 = cell[7] + x5;
auto x6 = -V{1}*cell[4] + x4 + x72 + x79;
auto x7 = x6*x80;
auto x8 = x2*x7;
auto x9 = V{1}/x0;
auto x10 = x9*(x2 + V{-1});
auto x81 = cell[3] - x10;
auto x11 = -x7*x8 + x73 + x74 + x81;
auto x12 = -V{1}*cell[1];
auto x13 = -V{1}*cell[2] + x12 + x4 + x75;
auto x14 = x13*x80;
auto x15 = x14*x2;
auto x16 = cell[3] + x12 - x15*x7 + x79;
auto x82 = V{2}*x16;
auto x17 = x82;
auto x18 = cell[1] - x14*x15 + x76 + x81;
auto x19 = V{0.5}*x68;
auto x20 = x19*util::pow(x3, V{2});
auto x21 = x11*x20;
auto x22 = V{2}*x21;
auto x23 = -V{1}*x22;
auto x24 = x23*x7;
auto x25 = x20*(x17 + x82);
auto x26 = -V{1}*x25;
auto x27 = x26*x7;
auto x28 = x18*x20;
auto x29 = V{2}*x28;
auto x30 = -V{1}*x29;
auto x31 = x14*x30;
auto x32 = x14*x27 + x14*x31 - V{2}*x19*util::pow(x3, V{1})*x78*x83*(x16*x17 + ((x11)*(x11)) + ((x18)*(x18))) + x23*x9 + x24*x7 + x30*x9;
auto x33 = -x6*x84;
auto x34 = -x64*x68*(cell.template getFieldComponent<olb::descriptors::POROSITY>(0) + V{-1})*((x67)*(x67))/(((x65)*(x65))*x66);
auto x35 = x34*x7;
auto x36 = x35*x85;
auto x37 = -x13*x84;
auto x38 = x14*x34;
auto x39 = x38*x85;
auto x40 = x33*x36 + x37*x39;
auto x41 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*(x15*x26 + x2*x24 + x23*x8);
auto x42 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*(x15*x30 + x2*x27 + x2*x31);
auto x43 = x33*x41 + x37*x42;
auto x87 = x32 + x43;
auto x88 = x40 + x87;
auto x89 = x22 + x87;
auto x93 = x29 + x88;
auto x44 = x71 + V{-1};
auto x45 = V{1}*x77;
auto x46 = x42*x45;
auto x47 = -V{1}*x46;
auto x48 = x41*x45;
auto x49 = x36*x45;
auto x50 = x40 + x49;
auto x90 = x48 + x50 + x89;
auto x51 = x39*x45;
auto x95 = x46 + x51;
auto x96 = x29 + x95;
auto x52 = -V{1}*x51;
auto x91 = x47 + x52;
auto x92 = x29 + x91;
auto x53 = x40 - V{1}*x49;
auto x54 = -V{1}*x48;
auto x94 = x53 + x54 + x89;
auto x55 = cell[0];
auto x56 = cell[1];
auto x57 = cell[2];
auto x58 = cell[3];
auto x59 = cell[4];
auto x60 = cell[5];
auto x61 = cell[6];
auto x62 = cell[7];
auto x63 = cell[8];
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[0] = x86*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(0)*x44 + x70*x88);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[1] = x86*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(1)*x44 + x70*(x26 + x90 + x92));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[2] = x86*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(2)*x44 + x70*(x91 + x93));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[3] = x86*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(3)*x44 + x70*(x25 + x92 + x94));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[4] = x86*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(4)*x44 + x70*x94);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[5] = x86*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(5)*x44 + x70*(x26 + x94 + x96));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[6] = x86*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(6)*x44 + x70*(x93 + x95));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[7] = x86*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(7)*x44 + x70*(x25 + x90 + x96));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[8] = x86*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(8)*x44 + x70*x90);
cell[0] = x55;
cell[1] = x56;
cell[2] = x57;
cell[3] = x58;
cell[4] = x59;
cell[5] = x60;
cell[6] = x61;
cell[7] = x62;
cell[8] = x63;
}
};

}

}

// Generation Info: commit=21521319fcfc23b85b6fd87340ebc70ad50346d4
