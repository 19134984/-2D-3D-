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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FreeEnergyMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::FreeEnergy, collision::FreeEnergy, dynamics::DefaultCombination>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x48 = cell.template getFieldComponent<olb::descriptors::CHEM_POTENTIAL>(0);
auto x49 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x50 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x51 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x52 = parameters.template get<collision::FreeEnergy::GAMMA>();
auto x53 = parameters.template get<descriptors::OMEGA>();
auto x54 = ((x49)*(x49));
auto x55 = ((x50)*(x50));
auto x56 = ((x51)*(x51));
auto x0 = V{1} - x53;
auto x1 = V{3};
auto x2 = x54 + x55 + x56;
auto x3 = V{0.5}*x1*x2;
auto x60 = V{1} - x3;
auto x4 = V{0.333333333333333};
auto x5 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{1};
auto x6 = -x1*x48*x52 + x5;
auto x57 = x53*x6;
auto x7 = -V{1}*x49;
auto x8 = V{0.5}*((x1)*(x1));
auto x9 = V{0.0555555555555556};
auto x59 = -x9;
auto x10 = x5*x9;
auto x11 = x57*x9;
auto x58 = -x11;
auto x12 = -V{1}*x50;
auto x13 = -V{1}*x51;
auto x14 = x12 + x7;
auto x15 = V{0.0277777777777778};
auto x16 = x15*x5;
auto x17 = x15*x57;
auto x18 = x50 + x7;
auto x19 = x13 + x7;
auto x20 = x51 + x7;
auto x21 = x12 + x13;
auto x22 = x12 + x51;
auto x23 = x49 + x50;
auto x24 = x12 + x49;
auto x25 = x49 + x51;
auto x26 = x13 + x49;
auto x27 = x50 + x51;
auto x28 = x13 + x50;
auto x29 = cell[0]*x0 - x4*x53*(x5*(x3 + V{-1}) + V{1}) - x57*(x4 + V{-1});
auto x30 = cell[10]*x0 + x53*(x10*(x1*x49 + x54*x8 + x60) + x59) + x58;
auto x31 = cell[11]*x0 + x53*(x10*(x1*x50 + x55*x8 + x60) + x59) + x58;
auto x32 = cell[12]*x0 + x53*(x10*(x1*x51 + x56*x8 + x60) + x59) + x58;
auto x33 = -(-cell[13]*x0 + x17 + x53*(x15 - x16*(x1*x23 + x60 + x8*((x23)*(x23)))));
auto x34 = -(-cell[14]*x0 + x17 + x53*(x15 - x16*(x1*x24 + x60 + x8*((x24)*(x24)))));
auto x35 = -(-cell[15]*x0 + x17 + x53*(x15 - x16*(x1*x25 + x60 + x8*((x25)*(x25)))));
auto x36 = -(-cell[16]*x0 + x17 + x53*(x15 - x16*(x1*x26 + x60 + x8*((x26)*(x26)))));
auto x37 = -(-cell[17]*x0 + x17 + x53*(x15 - x16*(x1*x27 + x60 + x8*((x27)*(x27)))));
auto x38 = -(-cell[18]*x0 + x17 + x53*(x15 - x16*(x1*x28 + x60 + x8*((x28)*(x28)))));
auto x39 = cell[1]*x0 + x53*(x10*(x1*x7 + x60 + x8*((x7)*(x7))) + x59) + x58;
auto x40 = cell[2]*x0 + x53*(x10*(x1*x12 + x60 + x8*((x12)*(x12))) + x59) + x58;
auto x41 = cell[3]*x0 + x53*(x10*(x1*x13 + x60 + x8*((x13)*(x13))) + x59) + x58;
auto x42 = -(-cell[4]*x0 + x17 + x53*(x15 - x16*(x1*x14 + x60 + x8*((x14)*(x14)))));
auto x43 = -(-cell[5]*x0 + x17 + x53*(x15 - x16*(x1*x18 + x60 + x8*((x18)*(x18)))));
auto x44 = -(-cell[6]*x0 + x17 + x53*(x15 - x16*(x1*x19 + x60 + x8*((x19)*(x19)))));
auto x45 = -(-cell[7]*x0 + x17 + x53*(x15 - x16*(x1*x20 + x60 + x8*((x20)*(x20)))));
auto x46 = -(-cell[8]*x0 + x17 + x53*(x15 - x16*(x1*x21 + x60 + x8*((x21)*(x21)))));
auto x47 = -(-cell[9]*x0 + x17 + x53*(x15 - x16*(x1*x22 + x60 + x8*((x22)*(x22)))));
cell[0] = x29;
cell[10] = x30;
cell[11] = x31;
cell[12] = x32;
cell[13] = x33;
cell[14] = x34;
cell[15] = x35;
cell[16] = x36;
cell[17] = x37;
cell[18] = x38;
cell[1] = x39;
cell[2] = x40;
cell[3] = x41;
cell[4] = x42;
cell[5] = x43;
cell[6] = x44;
cell[7] = x45;
cell[8] = x46;
cell[9] = x47;
return { x5, x2 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
