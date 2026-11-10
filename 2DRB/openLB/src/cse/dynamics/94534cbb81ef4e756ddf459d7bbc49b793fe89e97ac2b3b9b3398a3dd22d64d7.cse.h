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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FreeEnergyMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::FreeEnergy, collision::FreeEnergy, dynamics::DefaultCombination>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x29 = cell.template getFieldComponent<olb::descriptors::CHEM_POTENTIAL>(0);
auto x30 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x31 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x32 = parameters.template get<collision::FreeEnergy::GAMMA>();
auto x33 = parameters.template get<descriptors::OMEGA>();
auto x34 = ((x30)*(x30));
auto x35 = ((x31)*(x31));
auto x0 = V{1} - x33;
auto x1 = V{3};
auto x2 = x34 + x35;
auto x3 = V{0.5}*x1*x2;
auto x37 = V{1} - x3;
auto x4 = V{0.444444444444444};
auto x5 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + V{1};
auto x6 = -x1*x29*x32 + x5;
auto x36 = x33*x6;
auto x7 = -V{1}*x30;
auto x8 = x31 + x7;
auto x9 = V{0.5}*((x1)*(x1));
auto x10 = V{0.0277777777777778};
auto x11 = x10*x5;
auto x12 = x10*x36;
auto x13 = V{0.111111111111111};
auto x14 = x13*x5;
auto x15 = x13*x36;
auto x16 = -V{1}*x31;
auto x17 = x16 + x7;
auto x18 = x16 + x30;
auto x19 = x30 + x31;
auto x20 = cell[0]*x0 - x33*x4*(x5*(x3 + V{-1}) + V{1}) - x36*(x4 + V{-1});
auto x21 = -(-cell[1]*x0 + x12 + x33*(x10 - x11*(x1*x8 + x37 + x9*((x8)*(x8)))));
auto x22 = -(-cell[2]*x0 + x15 + x33*(x13 - x14*(x1*x7 + x37 + x9*((x7)*(x7)))));
auto x23 = -(-cell[3]*x0 + x12 + x33*(x10 - x11*(x1*x17 + x37 + x9*((x17)*(x17)))));
auto x24 = -(-cell[4]*x0 + x15 + x33*(x13 - x14*(x1*x16 + x37 + x9*((x16)*(x16)))));
auto x25 = -(-cell[5]*x0 + x12 + x33*(x10 - x11*(x1*x18 + x37 + x9*((x18)*(x18)))));
auto x26 = cell[6]*x0 - x15 - x33*(x13 - x14*(x1*x30 + x34*x9 + x37));
auto x27 = -(-cell[7]*x0 + x12 + x33*(x10 - x11*(x1*x19 + x37 + x9*((x19)*(x19)))));
auto x28 = cell[8]*x0 - x15 - x33*(x13 - x14*(x1*x31 + x35*x9 + x37));
cell[0] = x20;
cell[1] = x21;
cell[2] = x22;
cell[3] = x23;
cell[4] = x24;
cell[5] = x25;
cell[6] = x26;
cell[7] = x27;
cell[8] = x28;
return { x5, x2 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
