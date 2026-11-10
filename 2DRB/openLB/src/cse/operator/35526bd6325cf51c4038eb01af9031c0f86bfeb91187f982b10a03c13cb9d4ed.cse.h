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
struct CSE_O<operators::TrackAverageFunctorO<functors::DerivativeF<functors::FarFieldDragF<0>, descriptors::POPULATION, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::ThirdOrder, collision::FrozenOmegaFromCell<collision::ThirdOrderRLB>, dynamics::DefaultCombination> >, opti::AVERAGE_DJDF>,descriptors::D3Q19<FIELDS...>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
void apply(CELL& cell, PARAMETERS& parameters) any_platform {
using V = typename CELL::value_t;
auto x40 = parameters.template get<parameters::FACTOR>();
auto x37 = parameters.template get<parameters::CONVERSION_FORCE>();
auto x41 = parameters.template get<descriptors::LATTICE_TIME_T>();
auto x38 = parameters.template get<parameters::LATTICE_REFERENCE_VELOCITY>();
auto x39 = cell[10] + cell[13] + cell[14] + cell[15] + cell[16];
auto x44 = x37*x40;
auto x45 = V{1} / (x41);
auto x0 = V{-0.333333333333333};
auto x1 = cell[0] + cell[11] + cell[12] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + x39 + V{1};
auto x2 = x1 + V{-1};
auto x3 = V{1}/x1;
auto x42 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*x3;
auto x4 = x3*(-V{1}*cell[1] - V{1}*cell[4] - V{1}*cell[5] - V{1}*cell[6] - V{1}*cell[7] + x39);
auto x5 = V{1}*x4;
auto x6 = -x42*x5;
auto x7 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*x4;
auto x8 = x6*x7;
auto x9 = V{1}*x7*(x38 - x7);
auto x43 = x0 + x9;
auto x10 = x44*(x2*(x38*x6 - V{2}*x8) + x43);
auto x11 = x41 + V{-1};
auto x12 = -x42*(x5 + V{1});
auto x13 = x12*x7;
auto x14 = x44*(x2*(x12*x38 - V{2}*x13) + x43);
auto x15 = -x42*(x5 + V{-1});
auto x16 = x15*x7;
auto x17 = x44*(x2*(x15*x38 - V{2}*x16) + x43);
auto x18 = cell[0];
auto x19 = cell[1];
auto x20 = cell[2];
auto x21 = cell[3];
auto x22 = cell[4];
auto x23 = cell[5];
auto x24 = cell[6];
auto x25 = cell[7];
auto x26 = cell[8];
auto x27 = cell[9];
auto x28 = cell[10];
auto x29 = cell[11];
auto x30 = cell[12];
auto x31 = cell[13];
auto x32 = cell[14];
auto x33 = cell[15];
auto x34 = cell[16];
auto x35 = cell[17];
auto x36 = cell[18];
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[0] = x45*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(0)*x11 + x10);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[1] = x45*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(1)*x11 + x14);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[2] = x45*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(2)*x11 + x10);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[3] = x45*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(3)*x11 + x10);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[4] = x45*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(4)*x11 + x14);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[5] = x45*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(5)*x11 + x14);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[6] = x45*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(6)*x11 + x14);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[7] = x45*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(7)*x11 + x14);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[8] = x45*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(8)*x11 + x10);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[9] = x45*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(9)*x11 + x10);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[10] = x45*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(10)*x11 + x17);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[11] = x45*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(11)*x11 + x10);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[12] = x45*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(12)*x11 + x10);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[13] = x45*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(13)*x11 + x17);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[14] = x45*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(14)*x11 + x17);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[15] = x45*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(15)*x11 + x17);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[16] = x45*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(16)*x11 + x17);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[17] = x45*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(17)*x11 + x10);
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[18] = x45*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(18)*x11 + x10);
cell[0] = x18;
cell[1] = x19;
cell[2] = x20;
cell[3] = x21;
cell[4] = x22;
cell[5] = x23;
cell[6] = x24;
cell[7] = x25;
cell[8] = x26;
cell[9] = x27;
cell[10] = x28;
cell[11] = x29;
cell[12] = x30;
cell[13] = x31;
cell[14] = x32;
cell[15] = x33;
cell[16] = x34;
cell[17] = x35;
cell[18] = x36;
}
};

}

}

// Generation Info: commit=7ed6d8a1694b31c30b112c2d6e5b03f3c572d55c
