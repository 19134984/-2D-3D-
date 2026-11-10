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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, dynamics::DefaultCombination>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x69 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x70 = parameters.template get<descriptors::OMEGA>();
auto x71 = cell[15] + cell[16];
auto x72 = cell[10] + cell[14] + x71;
auto x73 = cell[1] + cell[4] + x72;
auto x74 = cell[12] + cell[7] + cell[9];
auto x75 = cell[3] + cell[8] + x74;
auto x76 = cell[13] + cell[5];
auto x77 = cell[17] + cell[18];
auto x78 = cell[11] + x76 + x77;
auto x79 = cell[2] + x78;
auto x0 = V{1}/x70;
auto x1 = -V{1}*cell[16];
auto x81 = cell[15] + x1;
auto x2 = -V{1}*cell[7];
auto x3 = -V{1}*cell[6];
auto x4 = -V{1}*cell[8];
auto x5 = -V{1}*cell[18];
auto x82 = cell[17] + x5;
auto x6 = cell[0] + cell[6] + x73 + x75 + x79 + V{1};
auto x80 = V{1} / (x6);
auto x7 = x80*(-V{1}*cell[3] + x3 + x4 + x74 + x81 + x82);
auto x86 = ((x7)*(x7));
auto x8 = -V{1}*cell[4];
auto x9 = -V{1}*cell[5];
auto x83 = cell[13] + x9;
auto x10 = x80*(-V{1}*cell[1] + x2 + x3 + x72 + x8 + x83);
auto x88 = ((x10)*(x10));
auto x11 = x10*x6;
auto x12 = cell[6] - x11*x7 + x2 + x81;
auto x13 = ((x12)*(x12));
auto x14 = -V{1}*cell[14];
auto x15 = -V{1}*cell[9];
auto x16 = x80*(-V{1}*cell[2] + x14 + x15 + x4 + x78 + x8);
auto x89 = ((x16)*(x16));
auto x17 = cell[4] - x11*x16 + x14 + x83;
auto x18 = V{3};
auto x87 = ((x18)*(x18));
auto x19 = V{1}*(x6 + V{-1})/x18;
auto x84 = -x19;
auto x85 = cell[6] + x84;
auto x20 = cell[7] - x10*x11 + x73 + x76 + x85;
auto x21 = -x6*x86 + x71 + x75 + x77 + x85;
auto x22 = x16*x6;
auto x23 = cell[8] + x15 - x22*x7 + x82;
auto x24 = cell[14] + cell[4] + cell[8] + cell[9] - x16*x22 + x79 + x84;
auto x25 = V{2}*V{1} / (x0 + util::sqrt(V{2.828428}*x80*x87*((x69)*(x69))*util::sqrt(V{2}*x13 + V{2}*((x17)*(x17)) + ((x20)*(x20)) + ((x21)*(x21)) + V{2}*((x23)*(x23)) + ((x24)*(x24))) + ((x0)*(x0))));
auto x26 = V{1} - x25;
auto x27 = x86 + x88 + x89;
auto x28 = V{0.5}*x18*x27;
auto x90 = V{1} - x28;
auto x29 = V{0.333333333333333};
auto x30 = -V{1}*x10;
auto x31 = V{0.5}*x87;
auto x32 = V{0.0555555555555556};
auto x33 = x32*x6;
auto x34 = -V{1}*x16;
auto x35 = -V{1}*x7;
auto x36 = x30 + x34;
auto x37 = V{0.0277777777777778};
auto x38 = x37*x6;
auto x39 = x16 + x30;
auto x40 = x30 + x35;
auto x41 = x30 + x7;
auto x42 = x34 + x35;
auto x43 = x34 + x7;
auto x44 = x10 + x16;
auto x45 = x10 + x34;
auto x46 = x10 + x7;
auto x47 = x10 + x35;
auto x48 = x16 + x7;
auto x49 = x16 + x35;
auto x50 = cell[0]*x26 - x25*x29*(x6*(x28 + V{-1}) + V{1});
auto x51 = cell[10]*x26 - x25*(x32 - x33*(x10*x18 + x31*x88 + x90));
auto x52 = cell[11]*x26 - x25*(x32 - x33*(x16*x18 + x31*x89 + x90));
auto x53 = cell[12]*x26 - x25*(x32 - x33*(x18*x7 + x31*x86 + x90));
auto x54 = cell[13]*x26 - x25*(x37 - x38*(x18*x44 + x31*((x44)*(x44)) + x90));
auto x55 = cell[14]*x26 - x25*(x37 - x38*(x18*x45 + x31*((x45)*(x45)) + x90));
auto x56 = cell[15]*x26 - x25*(x37 - x38*(x18*x46 + x31*((x46)*(x46)) + x90));
auto x57 = cell[16]*x26 - x25*(x37 - x38*(x18*x47 + x31*((x47)*(x47)) + x90));
auto x58 = cell[17]*x26 - x25*(x37 - x38*(x18*x48 + x31*((x48)*(x48)) + x90));
auto x59 = cell[18]*x26 - x25*(x37 - x38*(x18*x49 + x31*((x49)*(x49)) + x90));
auto x60 = cell[1]*x26 - x25*(x32 - x33*(x18*x30 + x31*((x30)*(x30)) + x90));
auto x61 = cell[2]*x26 - x25*(x32 - x33*(x18*x34 + x31*((x34)*(x34)) + x90));
auto x62 = cell[3]*x26 - x25*(x32 - x33*(x18*x35 + x31*((x35)*(x35)) + x90));
auto x63 = cell[4]*x26 - x25*(x37 - x38*(x18*x36 + x31*((x36)*(x36)) + x90));
auto x64 = cell[5]*x26 - x25*(x37 - x38*(x18*x39 + x31*((x39)*(x39)) + x90));
auto x65 = cell[6]*x26 - x25*(x37 - x38*(x18*x40 + x31*((x40)*(x40)) + x90));
auto x66 = cell[7]*x26 - x25*(x37 - x38*(x18*x41 + x31*((x41)*(x41)) + x90));
auto x67 = cell[8]*x26 - x25*(x37 - x38*(x18*x42 + x31*((x42)*(x42)) + x90));
auto x68 = cell[9]*x26 - x25*(x37 - x38*(x18*x43 + x31*((x43)*(x43)) + x90));
cell[0] = x50;
cell[10] = x51;
cell[11] = x52;
cell[12] = x53;
cell[13] = x54;
cell[14] = x55;
cell[15] = x56;
cell[16] = x57;
cell[17] = x58;
cell[18] = x59;
cell[1] = x60;
cell[2] = x61;
cell[3] = x62;
cell[4] = x63;
cell[5] = x64;
cell[6] = x65;
cell[7] = x66;
cell[8] = x67;
cell[9] = x68;
return { x6, x27 };
}
};

}

}

// Generation Info: commit=61323934e549f5f16850d6c6280b91bc057cb6c3
