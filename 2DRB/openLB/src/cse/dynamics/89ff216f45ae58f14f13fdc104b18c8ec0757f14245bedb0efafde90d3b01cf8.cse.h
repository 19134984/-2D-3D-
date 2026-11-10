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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::ScaledPlainGuo>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x80 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x81 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x82 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x83 = cell.template getFieldComponent<olb::descriptors::SCALAR>(0);
auto x84 = parameters.template get<descriptors::OMEGA>();
auto x85 = cell[10] + cell[14] + cell[16];
auto x86 = cell[11] + cell[18] + cell[5];
auto x87 = cell[12] + cell[7] + cell[9];
auto x0 = V{1} - x84;
auto x1 = V{3};
auto x2 = -V{1}*cell[4];
auto x89 = cell[13] + x2;
auto x3 = -V{1}*cell[6];
auto x90 = cell[15] + x3;
auto x4 = cell[0] + cell[13] + cell[15] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + x85 + x86 + x87 + V{1};
auto x88 = V{1} / (x4);
auto x5 = V{0.5}*x80 + x88*(-V{1}*cell[1] - V{1}*cell[5] - V{1}*cell[7] + x85 + x89 + x90);
auto x92 = ((x5)*(x5));
auto x95 = x1*x5;
auto x6 = -V{1}*cell[8];
auto x91 = cell[17] + x6;
auto x7 = V{0.5}*x81 + x88*(-V{1}*cell[14] - V{1}*cell[2] - V{1}*cell[9] + x86 + x89 + x91);
auto x93 = ((x7)*(x7));
auto x96 = x1*x7;
auto x8 = V{0.5}*x82 + x88*(-V{1}*cell[16] - V{1}*cell[18] - V{1}*cell[3] + x87 + x90 + x91);
auto x94 = ((x8)*(x8));
auto x97 = x1*x8;
auto x9 = x92 + x93 + x94;
auto x10 = V{0.5}*x1*x9;
auto x100 = V{1} - x10;
auto x11 = V{0.333333333333333};
auto x12 = x80*x83;
auto x13 = -x12*x95;
auto x14 = x81*x83;
auto x15 = -x14*x96;
auto x16 = x13 + x15;
auto x17 = x82*x83;
auto x18 = -x17*x97;
auto x99 = x15 + x18;
auto x101 = x13 + x18;
auto x19 = -x4*(V{0.5}*x84 + V{-1});
auto x20 = -V{1}*x5;
auto x21 = ((x1)*(x1));
auto x98 = V{1}*x21;
auto x22 = V{0.5}*x21;
auto x23 = V{0.0555555555555556};
auto x24 = x23*x4;
auto x25 = -x1*(x5 + V{1});
auto x26 = x19*x23;
auto x27 = -V{1}*x7;
auto x28 = -x1*(x7 + V{1});
auto x29 = -V{1}*x8;
auto x30 = -x1*(x8 + V{1});
auto x31 = x20 + x27;
auto x32 = V{0.0277777777777778};
auto x33 = x32*x4;
auto x34 = -x31*x98;
auto x35 = x19*x32;
auto x36 = x20 + x7;
auto x37 = x21*x36;
auto x38 = -x1*(x7 + V{-1});
auto x39 = x20 + x29;
auto x40 = -x39*x98;
auto x41 = x20 + x8;
auto x42 = x21*x41;
auto x43 = -x1*(x8 + V{-1});
auto x44 = x27 + x29;
auto x45 = -x44*x98;
auto x46 = x27 + x8;
auto x47 = x21*x46;
auto x48 = -x1*(x5 + V{-1});
auto x49 = x5 + x7;
auto x50 = x21*x49;
auto x51 = x27 + x5;
auto x52 = x21*x51;
auto x53 = x5 + x8;
auto x54 = x21*x53;
auto x55 = x29 + x5;
auto x56 = x21*x55;
auto x57 = x7 + x8;
auto x58 = x21*x57;
auto x59 = x29 + x7;
auto x60 = x21*x59;
auto x61 = cell[0]*x0 + x11*x19*(x16 + x18) - x11*x84*(x4*(x10 + V{-1}) + V{1});
auto x62 = cell[10]*x0 + x26*(x12*(x21*x5 + x48) + x99) - x84*(x23 - x24*(x100 + x22*x92 + x95));
auto x63 = cell[11]*x0 + x26*(x101 + x14*(x21*x7 + x38)) - x84*(x23 - x24*(x100 + x22*x93 + x96));
auto x64 = cell[12]*x0 + x26*(x16 + x17*(x21*x8 + x43)) - x84*(x23 - x24*(x100 + x22*x94 + x97));
auto x65 = cell[13]*x0 + x35*(x12*(x48 + x50) + x14*(x38 + x50) + x18) - x84*(x32 - x33*(x1*x49 + x100 + x22*((x49)*(x49))));
auto x66 = cell[14]*x0 + x35*(x12*(x48 + x52) + x14*(x28 - V{1}*x52) + x18) - x84*(x32 - x33*(x1*x51 + x100 + x22*((x51)*(x51))));
auto x67 = cell[15]*x0 + x35*(x12*(x48 + x54) + x15 + x17*(x43 + x54)) - x84*(x32 - x33*(x1*x53 + x100 + x22*((x53)*(x53))));
auto x68 = cell[16]*x0 + x35*(x12*(x48 + x56) + x15 + x17*(x30 - V{1}*x56)) - x84*(x32 - x33*(x1*x55 + x100 + x22*((x55)*(x55))));
auto x69 = cell[17]*x0 + x35*(x13 + x14*(x38 + x58) + x17*(x43 + x58)) - x84*(x32 - x33*(x1*x57 + x100 + x22*((x57)*(x57))));
auto x70 = cell[18]*x0 + x35*(x13 + x14*(x38 + x60) + x17*(x30 - V{1}*x60)) - x84*(x32 - x33*(x1*x59 + x100 + x22*((x59)*(x59))));
auto x71 = cell[1]*x0 + x26*(-x12*(x20*x98 - x25) + x99) - x84*(x23 - x24*(x1*x20 + x100 + x22*((x20)*(x20))));
auto x72 = cell[2]*x0 + x26*(x101 - x14*(x27*x98 - x28)) - x84*(x23 - x24*(x1*x27 + x100 + x22*((x27)*(x27))));
auto x73 = cell[3]*x0 + x26*(x16 - x17*(x29*x98 - x30)) - x84*(x23 - x24*(x1*x29 + x100 + x22*((x29)*(x29))));
auto x74 = cell[4]*x0 + x35*(x12*(x25 + x34) + x14*(x28 + x34) + x18) - x84*(x32 - x33*(x1*x31 + x100 + x22*((x31)*(x31))));
auto x75 = cell[5]*x0 + x35*(x12*(x25 - V{1}*x37) + x14*(x37 + x38) + x18) - x84*(x32 - x33*(x1*x36 + x100 + x22*((x36)*(x36))));
auto x76 = cell[6]*x0 + x35*(x12*(x25 + x40) + x15 + x17*(x30 + x40)) - x84*(x32 - x33*(x1*x39 + x100 + x22*((x39)*(x39))));
auto x77 = cell[7]*x0 + x35*(x12*(x25 - V{1}*x42) + x15 + x17*(x42 + x43)) - x84*(x32 - x33*(x1*x41 + x100 + x22*((x41)*(x41))));
auto x78 = cell[8]*x0 + x35*(x13 + x14*(x28 + x45) + x17*(x30 + x45)) - x84*(x32 - x33*(x1*x44 + x100 + x22*((x44)*(x44))));
auto x79 = cell[9]*x0 + x35*(x13 + x14*(x28 - V{1}*x47) + x17*(x43 + x47)) - x84*(x32 - x33*(x1*x46 + x100 + x22*((x46)*(x46))));
cell[0] = x61;
cell[10] = x62;
cell[11] = x63;
cell[12] = x64;
cell[13] = x65;
cell[14] = x66;
cell[15] = x67;
cell[16] = x68;
cell[17] = x69;
cell[18] = x70;
cell[1] = x71;
cell[2] = x72;
cell[3] = x73;
cell[4] = x74;
cell[5] = x75;
cell[6] = x76;
cell[7] = x77;
cell[8] = x78;
cell[9] = x79;
return { x4, x9 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
