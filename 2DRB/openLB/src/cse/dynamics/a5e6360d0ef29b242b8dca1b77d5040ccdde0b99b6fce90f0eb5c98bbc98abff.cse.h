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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::SecondOrder, collision::BGK, forcing::MCGuo<momenta::Identity> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x40 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x41 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x42 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x43 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x44 = parameters.template get<descriptors::OMEGA>();
auto x45 = ((x42)*(x42));
auto x46 = ((x43)*(x43));
auto x0 = V{1} - x44;
auto x1 = V{3};
auto x47 = x1*x42;
auto x48 = x1*x43;
auto x2 = x45 + x46;
auto x3 = V{0.5}*x1*x2;
auto x51 = V{1} - x3;
auto x4 = V{0.444444444444444};
auto x5 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + V{1};
auto x6 = -x40*x47;
auto x7 = -x41*x48;
auto x8 = -x5*(V{0.5}*x44 + V{-1});
auto x9 = -V{1}*x42;
auto x10 = x43 + x9;
auto x11 = ((x1)*(x1));
auto x49 = V{1}*x11;
auto x12 = V{0.5}*x11;
auto x13 = V{0.0277777777777778};
auto x14 = x13*x5;
auto x15 = x10*x11;
auto x16 = -x1*(x42 + V{1});
auto x50 = -x16;
auto x17 = -x1*(x43 + V{-1});
auto x18 = x13*x8;
auto x19 = V{0.111111111111111};
auto x20 = x19*x5;
auto x21 = x19*x8;
auto x22 = -V{1}*x43;
auto x23 = x22 + x9;
auto x24 = -x23*x49;
auto x25 = -x1*(x43 + V{1});
auto x26 = x22 + x42;
auto x27 = x11*x26;
auto x28 = -x1*(x42 + V{-1});
auto x29 = x42 + x43;
auto x30 = x11*x29;
auto x31 = cell[0]*x0 - x4*x44*(x5*(x3 + V{-1}) + V{1}) + x4*x8*(x6 + x7);
auto x32 = -(-cell[1]*x0 + x18*(x40*(V{1}*x15 + x50) - x41*(x15 + x17)) + x44*(x13 - x14*(x1*x10 + x12*((x10)*(x10)) + x51)));
auto x33 = -(-cell[2]*x0 + x21*(x40*(x49*x9 + x50) - x7) + x44*(x19 - x20*(x1*x9 + x12*((x9)*(x9)) + x51)));
auto x34 = cell[3]*x0 + x18*(x40*(x16 + x24) + x41*(x24 + x25)) - x44*(x13 - x14*(x1*x23 + x12*((x23)*(x23)) + x51));
auto x35 = -(-cell[4]*x0 + x21*(x41*(x22*x49 - x25) - x6) + x44*(x19 - x20*(x1*x22 + x12*((x22)*(x22)) + x51)));
auto x36 = cell[5]*x0 + x18*(x40*(x27 + x28) + x41*(x25 - V{1}*x27)) - x44*(x13 - x14*(x1*x26 + x12*((x26)*(x26)) + x51));
auto x37 = cell[6]*x0 + x21*(x40*(x11*x42 + x28) + x7) - x44*(x19 - x20*(x12*x45 + x47 + x51));
auto x38 = cell[7]*x0 + x18*(x40*(x28 + x30) + x41*(x17 + x30)) - x44*(x13 - x14*(x1*x29 + x12*((x29)*(x29)) + x51));
auto x39 = cell[8]*x0 + x21*(x41*(x11*x43 + x17) + x6) - x44*(x19 - x20*(x12*x46 + x48 + x51));
cell[0] = x31;
cell[1] = x32;
cell[2] = x33;
cell[3] = x34;
cell[4] = x35;
cell[5] = x36;
cell[6] = x37;
cell[7] = x38;
cell[8] = x39;
return { x5, x2 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
