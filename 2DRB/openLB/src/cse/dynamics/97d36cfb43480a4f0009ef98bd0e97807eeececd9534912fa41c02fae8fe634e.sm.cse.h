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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, multiphase::OmegaFromCell<collision::BGK>, forcing::Guo<momenta::Forced> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x57 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x58 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x61 = parameters.template get<multiphase::OMEGA_LIQUID>();
auto x62 = parameters.template get<multiphase::OMEGA_VAPOR>();
auto x63 = parameters.template get<multiphase::RHO_LIQUID>();
auto x64 = parameters.template get<multiphase::RHO_VAPOR>();
auto x59 = cell[5] + cell[6];
auto x60 = cell[1] + cell[8];
auto x65 = -x64;
auto x0 = cell[0] + cell[2] + cell[3] + cell[4] + cell[7] + x59 + x60 + V{1};
auto x66 = V{1} / (x0);
auto x1 = (x61*(x0 + x65) - x62*(x0 - x63))/(x63 + x65);
auto x2 = V{1} - x1;
auto x3 = V{3};
auto x4 = -V{1}*cell[3];
auto x67 = cell[7] + x4;
auto x5 = V{0.5}*x57;
auto x6 = x5 + x66*(-V{1}*cell[1] - V{1}*cell[2] + x59 + x67);
auto x68 = ((x6)*(x6));
auto x70 = x3*x6;
auto x7 = V{0.5}*x58;
auto x8 = x66*(-V{1}*cell[4] - V{1}*cell[5] + x60 + x67) + x7;
auto x69 = ((x8)*(x8));
auto x71 = x3*x8;
auto x9 = x68 + x69;
auto x10 = V{0.5}*x3*x9;
auto x73 = V{1} - x10;
auto x11 = V{0.444444444444444};
auto x12 = -x57*x70;
auto x13 = -x58*x71;
auto x14 = -x0*(V{0.5}*x1 + V{-1});
auto x15 = cell[0]*x2 - x1*x11*(x0*(x10 + V{-1}) + V{1}) + x11*x14*(x12 + x13);
auto x16 = -V{1}*x6;
auto x17 = x16 + x8;
auto x18 = ((x3)*(x3));
auto x74 = V{1}*x18;
auto x19 = V{0.5}*x18;
auto x20 = V{0.0277777777777778};
auto x21 = x0*x20;
auto x22 = x17*x18;
auto x23 = -x3*(x6 + V{1});
auto x72 = -x23;
auto x24 = -x3*(x8 + V{-1});
auto x25 = x14*x20;
auto x26 = -(-cell[1]*x2 + x1*(x20 - x21*(x17*x3 + x19*((x17)*(x17)) + x73)) + x25*(x57*(V{1}*x22 + x72) - x58*(x22 + x24)));
auto x27 = V{0.111111111111111};
auto x28 = x0*x27;
auto x29 = x14*x27;
auto x30 = cell[2]*x2 - x1*(x27 - x28*(x16*x3 + x19*((x16)*(x16)) + x73)) + x29*(x13 - x57*(x16*x74 + x72));
auto x31 = -V{1}*x8;
auto x32 = x16 + x31;
auto x33 = -x32*x74;
auto x34 = -x3*(x8 + V{1});
auto x35 = cell[3]*x2 - x1*(x20 - x21*(x19*((x32)*(x32)) + x3*x32 + x73)) + x25*(x57*(x23 + x33) + x58*(x33 + x34));
auto x36 = cell[4]*x2 - x1*(x27 - x28*(x19*((x31)*(x31)) + x3*x31 + x73)) + x29*(x12 - x58*(x31*x74 - x34));
auto x37 = x31 + x6;
auto x38 = x18*x37;
auto x39 = -x3*(x6 + V{-1});
auto x40 = cell[5]*x2 - x1*(x20 - x21*(x19*((x37)*(x37)) + x3*x37 + x73)) + x25*(x57*(x38 + x39) + x58*(x34 - V{1}*x38));
auto x41 = cell[6]*x2 - x1*(x27 - x28*(x19*x68 + x70 + x73)) + x29*(x13 + x57*(x18*x6 + x39));
auto x75 = x40 + x41;
auto x42 = x6 + x8;
auto x43 = x18*x42;
auto x44 = cell[7]*x2 - x1*(x20 - x21*(x19*((x42)*(x42)) + x3*x42 + x73)) + x25*(x57*(x39 + x43) + x58*(x24 + x43));
auto x45 = cell[8]*x2 - x1*(x27 - x28*(x19*x69 + x71 + x73)) + x29*(x12 + x58*(x18*x8 + x24));
auto x76 = x26 + x45;
auto x46 = x15 + x30 + x35 + x36 + x44 + x75 + x76 + V{1};
auto x77 = V{1} / (x46);
auto x47 = -V{1}*x35;
auto x78 = x44 + x47;
cell[0] = x15;
cell[1] = x26;
cell[2] = x30;
cell[3] = x35;
cell[4] = x36;
cell[5] = x40;
cell[6] = x41;
cell[7] = x44;
cell[8] = x45;
cell.template getFieldPointer<olb::descriptors::OMEGA>()[0] = x1;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x46;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x5 + x77*(-V{1}*x26 - V{1}*x30 + x75 + x78);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x7 + x77*(-V{1}*x36 - V{1}*x40 + x76 + x78);
return { x0, x9 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
