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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, dynamics::DefaultCombination>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x38 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x39 = parameters.template get<descriptors::OMEGA>();
auto x40 = cell[1] + cell[8];
auto x41 = cell[4] + x40;
auto x42 = cell[5] + cell[7];
auto x43 = cell[6] + x42;
auto x44 = cell[2] + x43;
auto x0 = V{1}/x39;
auto x1 = -V{1}*cell[3];
auto x2 = -V{1}*cell[5];
auto x46 = cell[7] + x2;
auto x3 = cell[0] + cell[3] + x41 + x44 + V{1};
auto x45 = V{1} / (x3);
auto x4 = x45*(-V{1}*cell[4] + x1 + x40 + x46);
auto x47 = ((x4)*(x4));
auto x5 = V{3};
auto x49 = ((x5)*(x5));
auto x6 = V{1}*(x3 + V{-1})/x5;
auto x48 = cell[3] - x6;
auto x7 = -x3*x47 + x41 + x42 + x48;
auto x8 = -V{1}*cell[1];
auto x9 = x45*(-V{1}*cell[2] + x1 + x43 + x8);
auto x50 = ((x9)*(x9));
auto x10 = x3*x9;
auto x11 = cell[3] - x10*x4 + x46 + x8;
auto x12 = cell[1] - x10*x9 + x44 + x48;
auto x13 = V{2}*V{1} / (x0 + util::sqrt(V{2.828428}*x45*x49*((x38)*(x38))*util::sqrt(V{2}*((x11)*(x11)) + ((x12)*(x12)) + ((x7)*(x7))) + ((x0)*(x0))));
auto x14 = V{1} - x13;
auto x15 = x47 + x50;
auto x16 = V{0.5}*x15*x5;
auto x51 = V{1} - x16;
auto x17 = V{0.444444444444444};
auto x18 = -V{1}*x9;
auto x19 = x18 + x4;
auto x20 = V{0.5}*x49;
auto x21 = V{0.0277777777777778};
auto x22 = x21*x3;
auto x23 = V{0.111111111111111};
auto x24 = x23*x3;
auto x25 = -V{1}*x4;
auto x26 = x18 + x25;
auto x27 = x25 + x9;
auto x28 = x4 + x9;
auto x29 = cell[0]*x14 - x13*x17*(x3*(x16 + V{-1}) + V{1});
auto x30 = cell[1]*x14 - x13*(x21 - x22*(x19*x5 + x20*((x19)*(x19)) + x51));
auto x31 = cell[2]*x14 - x13*(x23 - x24*(x18*x5 + x20*((x18)*(x18)) + x51));
auto x32 = cell[3]*x14 - x13*(x21 - x22*(x20*((x26)*(x26)) + x26*x5 + x51));
auto x33 = cell[4]*x14 - x13*(x23 - x24*(x20*((x25)*(x25)) + x25*x5 + x51));
auto x34 = cell[5]*x14 - x13*(x21 - x22*(x20*((x27)*(x27)) + x27*x5 + x51));
auto x35 = cell[6]*x14 - x13*(x23 - x24*(x20*x50 + x5*x9 + x51));
auto x36 = cell[7]*x14 - x13*(x21 - x22*(x20*((x28)*(x28)) + x28*x5 + x51));
auto x37 = cell[8]*x14 - x13*(x23 - x24*(x20*x47 + x4*x5 + x51));
cell[0] = x29;
cell[1] = x30;
cell[2] = x31;
cell[3] = x32;
cell[4] = x33;
cell[5] = x34;
cell[6] = x35;
cell[7] = x36;
cell[8] = x37;
return { x3, x15 };
}
};

}

}

// Generation Info: commit=61323934e549f5f16850d6c6280b91bc057cb6c3
