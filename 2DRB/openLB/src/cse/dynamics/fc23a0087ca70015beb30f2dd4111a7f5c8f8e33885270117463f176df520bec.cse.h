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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::BGK, forcing::PlainGuo>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x77 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x78 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x79 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x80 = parameters.template get<descriptors::OMEGA>();
auto x81 = cell[10] + cell[14] + cell[16];
auto x82 = cell[11] + cell[18] + cell[5];
auto x83 = cell[12] + cell[7] + cell[9];
auto x0 = V{1} - x80;
auto x1 = V{3};
auto x2 = -V{1}*cell[4];
auto x85 = cell[13] + x2;
auto x3 = -V{1}*cell[6];
auto x86 = cell[15] + x3;
auto x4 = cell[0] + cell[13] + cell[15] + cell[17] + cell[1] + cell[2] + cell[3] + cell[4] + cell[6] + cell[8] + x81 + x82 + x83 + V{1};
auto x84 = V{1} / (x4);
auto x5 = V{0.5}*x77 + x84*(-V{1}*cell[1] - V{1}*cell[5] - V{1}*cell[7] + x81 + x85 + x86);
auto x88 = ((x5)*(x5));
auto x91 = x1*x5;
auto x6 = -V{1}*cell[8];
auto x87 = cell[17] + x6;
auto x7 = V{0.5}*x78 + x84*(-V{1}*cell[14] - V{1}*cell[2] - V{1}*cell[9] + x82 + x85 + x87);
auto x89 = ((x7)*(x7));
auto x92 = x1*x7;
auto x8 = V{0.5}*x79 + x84*(-V{1}*cell[16] - V{1}*cell[18] - V{1}*cell[3] + x83 + x86 + x87);
auto x90 = ((x8)*(x8));
auto x93 = x1*x8;
auto x9 = x88 + x89 + x90;
auto x10 = V{0.5}*x1*x9;
auto x96 = V{1} - x10;
auto x11 = V{0.333333333333333};
auto x12 = -x77*x91;
auto x13 = -x78*x92;
auto x14 = x12 + x13;
auto x15 = -x79*x93;
auto x95 = x13 + x15;
auto x97 = x12 + x15;
auto x16 = -x4*(V{0.5}*x80 + V{-1});
auto x17 = -V{1}*x5;
auto x18 = ((x1)*(x1));
auto x94 = V{1}*x18;
auto x19 = V{0.5}*x18;
auto x20 = V{0.0555555555555556};
auto x21 = x20*x4;
auto x22 = -x1*(x5 + V{1});
auto x23 = x16*x20;
auto x24 = -V{1}*x7;
auto x25 = -x1*(x7 + V{1});
auto x26 = -V{1}*x8;
auto x27 = -x1*(x8 + V{1});
auto x28 = x17 + x24;
auto x29 = V{0.0277777777777778};
auto x30 = x29*x4;
auto x31 = -x28*x94;
auto x32 = x16*x29;
auto x33 = x17 + x7;
auto x34 = x18*x33;
auto x35 = -x1*(x7 + V{-1});
auto x36 = x17 + x26;
auto x37 = -x36*x94;
auto x38 = x17 + x8;
auto x39 = x18*x38;
auto x40 = -x1*(x8 + V{-1});
auto x41 = x24 + x26;
auto x42 = -x41*x94;
auto x43 = x24 + x8;
auto x44 = x18*x43;
auto x45 = -x1*(x5 + V{-1});
auto x46 = x5 + x7;
auto x47 = x18*x46;
auto x48 = x24 + x5;
auto x49 = x18*x48;
auto x50 = x5 + x8;
auto x51 = x18*x50;
auto x52 = x26 + x5;
auto x53 = x18*x52;
auto x54 = x7 + x8;
auto x55 = x18*x54;
auto x56 = x26 + x7;
auto x57 = x18*x56;
auto x58 = cell[0]*x0 + x11*x16*(x14 + x15) - x11*x80*(x4*(x10 + V{-1}) + V{1});
auto x59 = cell[10]*x0 + x23*(x77*(x18*x5 + x45) + x95) - x80*(x20 - x21*(x19*x88 + x91 + x96));
auto x60 = cell[11]*x0 + x23*(x78*(x18*x7 + x35) + x97) - x80*(x20 - x21*(x19*x89 + x92 + x96));
auto x61 = cell[12]*x0 + x23*(x14 + x79*(x18*x8 + x40)) - x80*(x20 - x21*(x19*x90 + x93 + x96));
auto x62 = cell[13]*x0 + x32*(x15 + x77*(x45 + x47) + x78*(x35 + x47)) - x80*(x29 - x30*(x1*x46 + x19*((x46)*(x46)) + x96));
auto x63 = cell[14]*x0 + x32*(x15 + x77*(x45 + x49) + x78*(x25 - V{1}*x49)) - x80*(x29 - x30*(x1*x48 + x19*((x48)*(x48)) + x96));
auto x64 = cell[15]*x0 + x32*(x13 + x77*(x45 + x51) + x79*(x40 + x51)) - x80*(x29 - x30*(x1*x50 + x19*((x50)*(x50)) + x96));
auto x65 = cell[16]*x0 + x32*(x13 + x77*(x45 + x53) + x79*(x27 - V{1}*x53)) - x80*(x29 - x30*(x1*x52 + x19*((x52)*(x52)) + x96));
auto x66 = cell[17]*x0 + x32*(x12 + x78*(x35 + x55) + x79*(x40 + x55)) - x80*(x29 - x30*(x1*x54 + x19*((x54)*(x54)) + x96));
auto x67 = cell[18]*x0 + x32*(x12 + x78*(x35 + x57) + x79*(x27 - V{1}*x57)) - x80*(x29 - x30*(x1*x56 + x19*((x56)*(x56)) + x96));
auto x68 = cell[1]*x0 + x23*(-x77*(x17*x94 - x22) + x95) - x80*(x20 - x21*(x1*x17 + x19*((x17)*(x17)) + x96));
auto x69 = cell[2]*x0 + x23*(-x78*(x24*x94 - x25) + x97) - x80*(x20 - x21*(x1*x24 + x19*((x24)*(x24)) + x96));
auto x70 = cell[3]*x0 + x23*(x14 - x79*(x26*x94 - x27)) - x80*(x20 - x21*(x1*x26 + x19*((x26)*(x26)) + x96));
auto x71 = cell[4]*x0 + x32*(x15 + x77*(x22 + x31) + x78*(x25 + x31)) - x80*(x29 - x30*(x1*x28 + x19*((x28)*(x28)) + x96));
auto x72 = cell[5]*x0 + x32*(x15 + x77*(x22 - V{1}*x34) + x78*(x34 + x35)) - x80*(x29 - x30*(x1*x33 + x19*((x33)*(x33)) + x96));
auto x73 = cell[6]*x0 + x32*(x13 + x77*(x22 + x37) + x79*(x27 + x37)) - x80*(x29 - x30*(x1*x36 + x19*((x36)*(x36)) + x96));
auto x74 = cell[7]*x0 + x32*(x13 + x77*(x22 - V{1}*x39) + x79*(x39 + x40)) - x80*(x29 - x30*(x1*x38 + x19*((x38)*(x38)) + x96));
auto x75 = cell[8]*x0 + x32*(x12 + x78*(x25 + x42) + x79*(x27 + x42)) - x80*(x29 - x30*(x1*x41 + x19*((x41)*(x41)) + x96));
auto x76 = cell[9]*x0 + x32*(x12 + x78*(x25 - V{1}*x44) + x79*(x40 + x44)) - x80*(x29 - x30*(x1*x43 + x19*((x43)*(x43)) + x96));
cell[0] = x58;
cell[10] = x59;
cell[11] = x60;
cell[12] = x61;
cell[13] = x62;
cell[14] = x63;
cell[15] = x64;
cell[16] = x65;
cell[17] = x66;
cell[18] = x67;
cell[1] = x68;
cell[2] = x69;
cell[3] = x70;
cell[4] = x71;
cell[5] = x72;
cell[6] = x73;
cell[7] = x74;
cell[8] = x75;
cell[9] = x76;
return { x4, x9 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
