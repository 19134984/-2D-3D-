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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, forcing::Guo<momenta::ForcedWithStress> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x59 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x60 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x61 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x62 = parameters.template get<descriptors::OMEGA>();
auto x63 = cell[3] + cell[7];
auto x64 = cell[1] + cell[8];
auto x65 = cell[4] + x64;
auto x66 = cell[5] + cell[6];
auto x67 = cell[2] + x66;
auto x0 = V{1}/x62;
auto x1 = -V{1}*cell[3];
auto x69 = cell[7] + x1;
auto x2 = -V{1}*cell[5];
auto x3 = cell[0] + x63 + x65 + x67 + V{1};
auto x68 = V{1} / (x3);
auto x4 = V{0.5}*x60;
auto x5 = x4 + x68*(-V{1}*cell[4] + x2 + x64 + x69);
auto x74 = ((x5)*(x5));
auto x6 = -x4 + x5;
auto x7 = V{3};
auto x72 = ((x7)*(x7));
auto x76 = x5*x7;
auto x8 = V{1}*(x3 + V{-1})/x7;
auto x71 = x63 - x8;
auto x9 = V{0.5}*x3;
auto x70 = V{2}*x9;
auto x10 = x6*x60;
auto x11 = cell[5] + x10*x70 - x3*((x6)*(x6)) + x65 + x71;
auto x12 = -V{1}*cell[1];
auto x13 = V{0.5}*x59;
auto x14 = x13 + x68*(-V{1}*cell[2] + x12 + x66 + x69);
auto x73 = ((x14)*(x14));
auto x75 = x14*x7;
auto x15 = -x13 + x14;
auto x16 = x15*x3;
auto x17 = x12 - x16*x6 + x2 + x63 + x9*(x15*x60 + x59*x6);
auto x18 = x15*x59;
auto x19 = cell[1] - x15*x16 + x18*x70 + x67 + x71;
auto x20 = V{2}*V{1} / (x0 + util::sqrt(V{2.828428}*x68*x72*((x61)*(x61))*util::sqrt(((x11)*(x11)) + V{2}*((x17)*(x17)) + ((x19)*(x19))) + ((x0)*(x0))));
auto x21 = V{1} - x20;
auto x22 = x73 + x74;
auto x23 = V{0.5}*x22*x7;
auto x79 = V{1} - x23;
auto x24 = V{0.444444444444444};
auto x25 = -x59*x75;
auto x26 = -x60*x76;
auto x27 = -x3*(V{0.5}*x20 + V{-1});
auto x28 = -V{1}*x14;
auto x29 = x28 + x5;
auto x30 = x72;
auto x77 = V{1}*x30;
auto x31 = V{0.5}*x30;
auto x32 = V{0.0277777777777778};
auto x33 = x3*x32;
auto x34 = x29*x30;
auto x35 = -x7*(x14 + V{1});
auto x78 = -x35;
auto x36 = -x7*(x5 + V{-1});
auto x37 = x27*x32;
auto x38 = V{0.111111111111111};
auto x39 = x3*x38;
auto x40 = x27*x38;
auto x41 = -V{1}*x5;
auto x42 = x28 + x41;
auto x43 = -x42*x77;
auto x44 = -x7*(x5 + V{1});
auto x45 = x14 + x41;
auto x46 = x30*x45;
auto x47 = -x7*(x14 + V{-1});
auto x48 = x14 + x5;
auto x49 = x30*x48;
auto x50 = cell[0]*x21 - x20*x24*(x3*(x23 + V{-1}) + V{1}) + x24*x27*(x25 + x26);
auto x51 = -(-cell[1]*x21 + x20*(x32 - x33*(x29*x7 + x31*((x29)*(x29)) + x79)) + x37*(x59*(V{1}*x34 + x78) - x60*(x34 + x36)));
auto x52 = cell[2]*x21 - x20*(x38 - x39*(x28*x7 + x31*((x28)*(x28)) + x79)) + x40*(x26 - x59*(x28*x77 + x78));
auto x53 = cell[3]*x21 - x20*(x32 - x33*(x31*((x42)*(x42)) + x42*x7 + x79)) + x37*(x59*(x35 + x43) + x60*(x43 + x44));
auto x54 = cell[4]*x21 - x20*(x38 - x39*(x31*((x41)*(x41)) + x41*x7 + x79)) + x40*(x25 - x60*(x41*x77 - x44));
auto x55 = cell[5]*x21 - x20*(x32 - x33*(x31*((x45)*(x45)) + x45*x7 + x79)) + x37*(x59*(x46 + x47) + x60*(x44 - V{1}*x46));
auto x56 = cell[6]*x21 - x20*(x38 - x39*(x31*x73 + x75 + x79)) + x40*(x26 + x59*(x14*x30 + x47));
auto x57 = cell[7]*x21 - x20*(x32 - x33*(x31*((x48)*(x48)) + x48*x7 + x79)) + x37*(x59*(x47 + x49) + x60*(x36 + x49));
auto x58 = cell[8]*x21 - x20*(x38 - x39*(x31*x74 + x76 + x79)) + x40*(x25 + x60*(x30*x5 + x36));
cell[0] = x50;
cell[1] = x51;
cell[2] = x52;
cell[3] = x53;
cell[4] = x54;
cell[5] = x55;
cell[6] = x56;
cell[7] = x57;
cell[8] = x58;
return { x3, x22 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
