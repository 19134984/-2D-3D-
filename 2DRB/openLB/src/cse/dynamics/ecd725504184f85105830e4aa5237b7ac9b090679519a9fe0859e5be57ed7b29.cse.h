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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::FreeEnergyInletOutletDensity, momenta::FreeEnergyInletOutletMomentum<0, 1>, momenta::RegularizedBoundaryStress<0, 1>, momenta::DefineSeparately>, equilibria::FreeEnergy, collision::FreeEnergyInletOutlet<0, 1>, dynamics::DefaultCombination>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x32 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x33 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x36 = parameters.template get<descriptors::OMEGA>();
auto x34 = ((x33)*(x33));
auto x35 = x32*x36;
auto x0 = V{1} - x36;
auto x1 = V{3};
auto x2 = x34;
auto x3 = V{0.5}*x1*x2;
auto x4 = V{1} - x3;
auto x5 = V{0.444444444444444};
auto x6 = cell[0]*x0 - x35*(x5 + V{-1}) + x36*x5*(x32*x4 + V{-1});
auto x7 = V{0.0277777777777778};
auto x8 = V{0.111111111111111};
auto x37 = x32*x8;
auto x38 = -x8;
auto x9 = x37;
auto x10 = x36*(x38 + x4*x9);
auto x11 = x36*x37;
auto x39 = -x11;
auto x40 = x10 + x39;
auto x12 = cell[4]*x0 + x40;
auto x13 = x1*x33 - x3 + V{0.5}*x34*((x1)*(x1)) + V{1};
auto x14 = x36*x7*(x13*x32 + V{-1});
auto x15 = x35*x7;
auto x41 = x14 - x15;
auto x16 = cell[5]*x0 + x41;
auto x17 = cell[6]*x0 + x36*(x13*x9 + x38) + x39;
auto x18 = cell[7]*x0 + x41;
auto x19 = cell[8]*x0 + x40;
auto x20 = -x12 - x16 - x17 - x18 - x19 + x32 - x6 + V{-1};
auto x21 = V{2}*x7 + x8;
auto x42 = x20/x21;
auto x22 = x42*x7;
cell[0] = x6;
cell[1] = x22;
cell[2] = x42*x8;
cell[3] = x22;
cell[4] = x12;
cell[5] = x16;
cell[6] = x17;
cell[7] = x18;
cell[8] = x19;
return { x32, x2 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
