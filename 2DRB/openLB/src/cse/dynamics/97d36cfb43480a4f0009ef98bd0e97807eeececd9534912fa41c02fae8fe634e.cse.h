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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, multiphase::OmegaFromCell<collision::BGK>, forcing::Guo<momenta::Forced> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x44 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x45 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x48 = parameters.template get<multiphase::OMEGA_LIQUID>();
auto x49 = parameters.template get<multiphase::OMEGA_VAPOR>();
auto x50 = parameters.template get<multiphase::RHO_LIQUID>();
auto x51 = parameters.template get<multiphase::RHO_VAPOR>();
auto x46 = cell[5] + cell[6];
auto x47 = cell[1] + cell[8];
auto x52 = -x51;
auto x0 = cell[0] + cell[2] + cell[3] + cell[4] + cell[7] + x46 + x47 + V{1};
auto x53 = V{1} / (x0);
auto x1 = (x48*(x0 + x52) - x49*(x0 - x50))/(x50 + x52);
auto x2 = V{1} - x1;
auto x3 = V{3};
auto x4 = -V{1}*cell[3];
auto x54 = cell[7] + x4;
auto x5 = V{0.5}*x44 + x53*(-V{1}*cell[1] - V{1}*cell[2] + x46 + x54);
auto x55 = ((x5)*(x5));
auto x57 = x3*x5;
auto x6 = V{0.5}*x45 + x53*(-V{1}*cell[4] - V{1}*cell[5] + x47 + x54);
auto x56 = ((x6)*(x6));
auto x58 = x3*x6;
auto x7 = x55 + x56;
auto x8 = V{0.5}*x3*x7;
auto x61 = V{1} - x8;
auto x9 = V{0.444444444444444};
auto x10 = -x44*x57;
auto x11 = -x45*x58;
auto x12 = -x0*(V{0.5}*x1 + V{-1});
auto x13 = -V{1}*x5;
auto x14 = x13 + x6;
auto x15 = ((x3)*(x3));
auto x59 = V{1}*x15;
auto x16 = V{0.5}*x15;
auto x17 = V{0.0277777777777778};
auto x18 = x0*x17;
auto x19 = x14*x15;
auto x20 = -x3*(x5 + V{1});
auto x60 = -x20;
auto x21 = -x3*(x6 + V{-1});
auto x22 = x12*x17;
auto x23 = V{0.111111111111111};
auto x24 = x0*x23;
auto x25 = x12*x23;
auto x26 = -V{1}*x6;
auto x27 = x13 + x26;
auto x28 = -x27*x59;
auto x29 = -x3*(x6 + V{1});
auto x30 = x26 + x5;
auto x31 = x15*x30;
auto x32 = -x3*(x5 + V{-1});
auto x33 = x5 + x6;
auto x34 = x15*x33;
auto x35 = cell[0]*x2 - x1*x9*(x0*(x8 + V{-1}) + V{1}) + x12*x9*(x10 + x11);
auto x36 = -(-cell[1]*x2 + x1*(x17 - x18*(x14*x3 + x16*((x14)*(x14)) + x61)) + x22*(x44*(V{1}*x19 + x60) - x45*(x19 + x21)));
auto x37 = cell[2]*x2 - x1*(x23 - x24*(x13*x3 + x16*((x13)*(x13)) + x61)) + x25*(x11 - x44*(x13*x59 + x60));
auto x38 = cell[3]*x2 - x1*(x17 - x18*(x16*((x27)*(x27)) + x27*x3 + x61)) + x22*(x44*(x20 + x28) + x45*(x28 + x29));
auto x39 = cell[4]*x2 - x1*(x23 - x24*(x16*((x26)*(x26)) + x26*x3 + x61)) + x25*(x10 - x45*(x26*x59 - x29));
auto x40 = cell[5]*x2 - x1*(x17 - x18*(x16*((x30)*(x30)) + x3*x30 + x61)) + x22*(x44*(x31 + x32) + x45*(x29 - V{1}*x31));
auto x41 = cell[6]*x2 - x1*(x23 - x24*(x16*x55 + x57 + x61)) + x25*(x11 + x44*(x15*x5 + x32));
auto x42 = cell[7]*x2 - x1*(x17 - x18*(x16*((x33)*(x33)) + x3*x33 + x61)) + x22*(x44*(x32 + x34) + x45*(x21 + x34));
auto x43 = cell[8]*x2 - x1*(x23 - x24*(x16*x56 + x58 + x61)) + x25*(x10 + x45*(x15*x6 + x21));
cell.template getFieldPointer<olb::descriptors::OMEGA>()[0] = x1;
cell[0] = x35;
cell[1] = x36;
cell[2] = x37;
cell[3] = x38;
cell[4] = x39;
cell[5] = x40;
cell[6] = x41;
cell[7] = x42;
cell[8] = x43;
return { x0, x7 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
