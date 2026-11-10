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

namespace dynamics {

template <typename T, typename... FIELDS>
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::FreeEnergyInletOutletDensity, momenta::FreeEnergyInletOutletMomentum<0, -1>, momenta::RegularizedBoundaryStress<0, -1>, momenta::DefineSeparately>, equilibria::FreeEnergy, collision::FreeEnergyInletOutlet<0, -1>, dynamics::DefaultCombination>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x33 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x34 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x37 = parameters.template get<descriptors::OMEGA>();
auto x35 = x33*x37;
auto x0 = V{1} - x37;
auto x1 = V{3};
auto x2 = ((x34)*(x34));
auto x3 = V{0.5}*x1*x2;
auto x4 = V{1} - x3;
auto x5 = V{0.444444444444444};
auto x6 = cell[0]*x0 - x35*(x5 + V{-1}) + x37*x5*(x33*x4 + V{-1});
auto x7 = -V{1}*x34;
auto x8 = x1*x7 - x3 + V{1} + V{0.5}*((x1)*(x1))*((x7)*(x7));
auto x9 = V{0.0277777777777778};
auto x10 = x37*x9*(x33*x8 + V{-1});
auto x11 = x35*x9;
auto x36 = x10 - x11;
auto x12 = cell[1]*x0 + x36;
auto x13 = V{0.111111111111111};
auto x38 = x13*x33;
auto x14 = x38;
auto x15 = x37*x38;
auto x16 = cell[2]*x0 - x15 - x37*(x13 - x14*x8);
auto x17 = cell[3]*x0 + x36;
auto x18 = -x37*(x13 - x14*x4);
auto x39 = -x15 + x18;
auto x19 = cell[4]*x0 + x39;
auto x20 = cell[8]*x0 + x39;
auto x21 = -x12 - x16 - x17 - x19 - x20 + x33 - x6 + V{-1};
auto x22 = x13 + V{2}*x9;
auto x40 = x21/x22;
auto x23 = x40*x9;
cell[0] = x6;
cell[1] = x12;
cell[2] = x16;
cell[3] = x17;
cell[4] = x19;
cell[5] = x23;
cell[6] = x13*x40;
cell[7] = x23;
cell[8] = x20;
return { x33, x2 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
