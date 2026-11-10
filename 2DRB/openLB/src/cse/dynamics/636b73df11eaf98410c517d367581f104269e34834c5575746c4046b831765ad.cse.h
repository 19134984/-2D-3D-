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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::SecondOrder, collision::BGK, forcing::MCGuo<momenta::Identity> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x71 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x72 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x73 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x74 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x75 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x76 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x77 = parameters.template get<descriptors::OMEGA>();
auto x78 = ((x74)*(x74));
auto x79 = ((x75)*(x75));
auto x80 = ((x76)*(x76));
auto x0 = V{1} - x77;
auto x1 = V{3};
auto x81 = x1*x74;
auto x82 = x1*x75;
auto x83 = x1*x76;
auto x2 = x78 + x79 + x80;
auto x3 = V{0.5}*x1*x2;
auto x86 = V{1} - x3;
auto x4 = V{0.333333333333333};
auto x5 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{1};
auto x6 = -x71*x81;
auto x7 = -x72*x82;
auto x8 = x6 + x7;
auto x9 = -x73*x83;
auto x85 = x7 + x9;
auto x87 = x6 + x9;
auto x10 = -x5*(V{0.5}*x77 + V{-1});
auto x11 = -V{1}*x74;
auto x12 = ((x1)*(x1));
auto x84 = V{1}*x12;
auto x13 = V{0.5}*x12;
auto x14 = V{0.0555555555555556};
auto x15 = x14*x5;
auto x16 = -x1*(x74 + V{1});
auto x17 = x10*x14;
auto x18 = -V{1}*x75;
auto x19 = -x1*(x75 + V{1});
auto x20 = -V{1}*x76;
auto x21 = -x1*(x76 + V{1});
auto x22 = x11 + x18;
auto x23 = V{0.0277777777777778};
auto x24 = x23*x5;
auto x25 = -x22*x84;
auto x26 = x10*x23;
auto x27 = x11 + x75;
auto x28 = x12*x27;
auto x29 = -x1*(x75 + V{-1});
auto x30 = x11 + x20;
auto x31 = -x30*x84;
auto x32 = x11 + x76;
auto x33 = x12*x32;
auto x34 = -x1*(x76 + V{-1});
auto x35 = x18 + x20;
auto x36 = -x35*x84;
auto x37 = x18 + x76;
auto x38 = x12*x37;
auto x39 = -x1*(x74 + V{-1});
auto x40 = x74 + x75;
auto x41 = x12*x40;
auto x42 = x18 + x74;
auto x43 = x12*x42;
auto x44 = x74 + x76;
auto x45 = x12*x44;
auto x46 = x20 + x74;
auto x47 = x12*x46;
auto x48 = x75 + x76;
auto x49 = x12*x48;
auto x50 = x20 + x75;
auto x51 = x12*x50;
auto x52 = cell[0]*x0 + x10*x4*(x8 + x9) - x4*x77*(x5*(x3 + V{-1}) + V{1});
auto x53 = cell[10]*x0 + x17*(x71*(x12*x74 + x39) + x85) - x77*(x14 - x15*(x13*x78 + x81 + x86));
auto x54 = cell[11]*x0 + x17*(x72*(x12*x75 + x29) + x87) - x77*(x14 - x15*(x13*x79 + x82 + x86));
auto x55 = cell[12]*x0 + x17*(x73*(x12*x76 + x34) + x8) - x77*(x14 - x15*(x13*x80 + x83 + x86));
auto x56 = cell[13]*x0 + x26*(x71*(x39 + x41) + x72*(x29 + x41) + x9) - x77*(x23 - x24*(x1*x40 + x13*((x40)*(x40)) + x86));
auto x57 = cell[14]*x0 + x26*(x71*(x39 + x43) + x72*(x19 - V{1}*x43) + x9) - x77*(x23 - x24*(x1*x42 + x13*((x42)*(x42)) + x86));
auto x58 = cell[15]*x0 + x26*(x7 + x71*(x39 + x45) + x73*(x34 + x45)) - x77*(x23 - x24*(x1*x44 + x13*((x44)*(x44)) + x86));
auto x59 = cell[16]*x0 + x26*(x7 + x71*(x39 + x47) + x73*(x21 - V{1}*x47)) - x77*(x23 - x24*(x1*x46 + x13*((x46)*(x46)) + x86));
auto x60 = cell[17]*x0 + x26*(x6 + x72*(x29 + x49) + x73*(x34 + x49)) - x77*(x23 - x24*(x1*x48 + x13*((x48)*(x48)) + x86));
auto x61 = cell[18]*x0 + x26*(x6 + x72*(x29 + x51) + x73*(x21 - V{1}*x51)) - x77*(x23 - x24*(x1*x50 + x13*((x50)*(x50)) + x86));
auto x62 = cell[1]*x0 + x17*(-x71*(x11*x84 - x16) + x85) - x77*(x14 - x15*(x1*x11 + x13*((x11)*(x11)) + x86));
auto x63 = cell[2]*x0 + x17*(-x72*(x18*x84 - x19) + x87) - x77*(x14 - x15*(x1*x18 + x13*((x18)*(x18)) + x86));
auto x64 = -(-cell[3]*x0 + x17*(x73*(x20*x84 - x21) - x8) + x77*(x14 - x15*(x1*x20 + x13*((x20)*(x20)) + x86)));
auto x65 = cell[4]*x0 + x26*(x71*(x16 + x25) + x72*(x19 + x25) + x9) - x77*(x23 - x24*(x1*x22 + x13*((x22)*(x22)) + x86));
auto x66 = cell[5]*x0 + x26*(x71*(x16 - V{1}*x28) + x72*(x28 + x29) + x9) - x77*(x23 - x24*(x1*x27 + x13*((x27)*(x27)) + x86));
auto x67 = cell[6]*x0 + x26*(x7 + x71*(x16 + x31) + x73*(x21 + x31)) - x77*(x23 - x24*(x1*x30 + x13*((x30)*(x30)) + x86));
auto x68 = cell[7]*x0 + x26*(x7 + x71*(x16 - V{1}*x33) + x73*(x33 + x34)) - x77*(x23 - x24*(x1*x32 + x13*((x32)*(x32)) + x86));
auto x69 = cell[8]*x0 + x26*(x6 + x72*(x19 + x36) + x73*(x21 + x36)) - x77*(x23 - x24*(x1*x35 + x13*((x35)*(x35)) + x86));
auto x70 = cell[9]*x0 + x26*(x6 + x72*(x19 - V{1}*x38) + x73*(x34 + x38)) - x77*(x23 - x24*(x1*x37 + x13*((x37)*(x37)) + x86));
cell[0] = x52;
cell[10] = x53;
cell[11] = x54;
cell[12] = x55;
cell[13] = x56;
cell[14] = x57;
cell[15] = x58;
cell[16] = x59;
cell[17] = x60;
cell[18] = x61;
cell[1] = x62;
cell[2] = x63;
cell[3] = x64;
cell[4] = x65;
cell[5] = x66;
cell[6] = x67;
cell[7] = x68;
cell[8] = x69;
cell[9] = x70;
return { x5, x2 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
