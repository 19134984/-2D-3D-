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
struct CSE<dynamics::Tuple<T, descriptors::D3Q15<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, multiphase::OmegaFromCell<collision::BGK>, forcing::Guo<momenta::Forced> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x74 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x75 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x76 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x79 = parameters.template get<multiphase::OMEGA_LIQUID>();
auto x80 = parameters.template get<multiphase::OMEGA_VAPOR>();
auto x81 = parameters.template get<multiphase::RHO_LIQUID>();
auto x82 = parameters.template get<multiphase::RHO_VAPOR>();
auto x77 = cell[14] + cell[8];
auto x78 = cell[6] + cell[9];
auto x83 = cell[10] + cell[11] + cell[5];
auto x84 = -x82;
auto x0 = cell[0] + cell[12] + cell[13] + cell[1] + cell[2] + cell[3] + cell[4] + cell[7] + x77 + x78 + x83 + V{1};
auto x85 = V{1} / (x0);
auto x1 = (x79*(x0 + x84) - x80*(x0 - x81))/(x81 + x84);
auto x2 = V{1} - x1;
auto x3 = V{3};
auto x4 = -V{1}*cell[4];
auto x5 = -V{1}*cell[5];
auto x86 = cell[11] + cell[12] + x4 + x5;
auto x6 = -V{1}*cell[6];
auto x87 = cell[13] + x6;
auto x7 = V{0.5}*x74 + x85*(-V{1}*cell[1] - V{1}*cell[7] + x77 + x86 + x87);
auto x90 = ((x7)*(x7));
auto x92 = x3*x7;
auto x8 = -V{1}*cell[14];
auto x88 = cell[7] + x8;
auto x9 = V{0.5}*x75 + x85*(-V{1}*cell[13] - V{1}*cell[2] + x78 + x86 + x88);
auto x91 = ((x9)*(x9));
auto x93 = x3*x9;
auto x10 = V{0.5}*x76 + x85*(-V{1}*cell[12] - V{1}*cell[3] + x4 + x83 + x87 + x88);
auto x89 = ((x10)*(x10));
auto x94 = x10*x3;
auto x11 = x89 + x90 + x91;
auto x12 = V{0.5}*x11*x3;
auto x96 = V{1} - x12;
auto x13 = V{0.222222222222222};
auto x14 = -x74*x92;
auto x15 = -x75*x93;
auto x16 = x14 + x15;
auto x17 = -x76*x94;
auto x97 = x15 + x17;
auto x98 = x14 + x17;
auto x18 = -x0*(V{0.5}*x1 + V{-1});
auto x19 = -V{1}*x7;
auto x20 = ((x3)*(x3));
auto x95 = V{1}*x20;
auto x21 = V{0.5}*x20;
auto x22 = V{0.111111111111111};
auto x23 = x0*x22;
auto x24 = -x3*(x7 + V{1});
auto x25 = x18*x22;
auto x26 = -V{1}*x9;
auto x27 = -x3*(x9 + V{1});
auto x28 = -V{1}*x10;
auto x29 = -x3*(x10 + V{1});
auto x30 = x19 + x26;
auto x31 = x28 + x30;
auto x32 = V{0.0138888888888889};
auto x33 = x0*x32;
auto x34 = -x31*x95;
auto x35 = x18*x32;
auto x36 = x10 + x30;
auto x37 = x20*x36;
auto x38 = -V{1}*x37;
auto x39 = -x3*(x10 + V{-1});
auto x40 = x19 + x9;
auto x41 = x28 + x40;
auto x42 = x20*x41;
auto x43 = -V{1}*x42;
auto x44 = -x3*(x9 + V{-1});
auto x45 = x10 + x40;
auto x46 = x20*x45;
auto x47 = -x3*(x7 + V{-1});
auto x48 = x7 + x9;
auto x49 = x10 + x48;
auto x50 = x20*x49;
auto x51 = x28 + x48;
auto x52 = x20*x51;
auto x53 = x26 + x7;
auto x54 = x10 + x53;
auto x55 = x20*x54;
auto x56 = x28 + x53;
auto x57 = x20*x56;
auto x58 = -V{1}*x57;
auto x59 = cell[0]*x2 - x1*x13*(x0*(x12 + V{-1}) + V{1}) + x13*x18*(x16 + x17);
auto x60 = cell[10]*x2 - x1*(x22 - x23*(x21*x89 + x94 + x96)) + x25*(x16 + x76*(x10*x20 + x39));
auto x61 = cell[11]*x2 - x1*(x32 - x33*(x21*((x49)*(x49)) + x3*x49 + x96)) + x35*(x74*(x47 + x50) + x75*(x44 + x50) + x76*(x39 + x50));
auto x62 = cell[12]*x2 - x1*(x32 - x33*(x21*((x51)*(x51)) + x3*x51 + x96)) + x35*(x74*(x47 + x52) + x75*(x44 + x52) + x76*(x29 - V{1}*x52));
auto x63 = cell[13]*x2 - x1*(x32 - x33*(x21*((x54)*(x54)) + x3*x54 + x96)) + x35*(x74*(x47 + x55) + x75*(x27 - V{1}*x55) + x76*(x39 + x55));
auto x64 = cell[14]*x2 - x1*(x32 - x33*(x21*((x56)*(x56)) + x3*x56 + x96)) + x35*(x74*(x47 + x57) + x75*(x27 + x58) + x76*(x29 + x58));
auto x65 = cell[1]*x2 - x1*(x22 - x23*(x19*x3 + x21*((x19)*(x19)) + x96)) + x25*(-x74*(x19*x95 - x24) + x97);
auto x66 = cell[2]*x2 - x1*(x22 - x23*(x21*((x26)*(x26)) + x26*x3 + x96)) + x25*(-x75*(x26*x95 - x27) + x98);
auto x67 = cell[3]*x2 - x1*(x22 - x23*(x21*((x28)*(x28)) + x28*x3 + x96)) + x25*(x16 - x76*(x28*x95 - x29));
auto x68 = cell[4]*x2 - x1*(x32 - x33*(x21*((x31)*(x31)) + x3*x31 + x96)) + x35*(x74*(x24 + x34) + x75*(x27 + x34) + x76*(x29 + x34));
auto x69 = cell[5]*x2 - x1*(x32 - x33*(x21*((x36)*(x36)) + x3*x36 + x96)) + x35*(x74*(x24 + x38) + x75*(x27 + x38) + x76*(x37 + x39));
auto x70 = cell[6]*x2 - x1*(x32 - x33*(x21*((x41)*(x41)) + x3*x41 + x96)) + x35*(x74*(x24 + x43) + x75*(x42 + x44) + x76*(x29 + x43));
auto x71 = cell[7]*x2 - x1*(x32 - x33*(x21*((x45)*(x45)) + x3*x45 + x96)) + x35*(x74*(x24 - V{1}*x46) + x75*(x44 + x46) + x76*(x39 + x46));
auto x72 = cell[8]*x2 - x1*(x22 - x23*(x21*x90 + x92 + x96)) + x25*(x74*(x20*x7 + x47) + x97);
auto x73 = cell[9]*x2 - x1*(x22 - x23*(x21*x91 + x93 + x96)) + x25*(x75*(x20*x9 + x44) + x98);
cell.template getFieldPointer<olb::descriptors::OMEGA>()[0] = x1;
cell[0] = x59;
cell[10] = x60;
cell[11] = x61;
cell[12] = x62;
cell[13] = x63;
cell[14] = x64;
cell[1] = x65;
cell[2] = x66;
cell[3] = x67;
cell[4] = x68;
cell[5] = x69;
cell[6] = x70;
cell[7] = x71;
cell[8] = x72;
cell[9] = x73;
return { x0, x11 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
