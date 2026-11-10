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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::SecondOrder, collision::BGK, forcing::MCGuo<momenta::Identity> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x90 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x91 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x92 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x93 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x94 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x95 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(2);
auto x96 = parameters.template get<descriptors::OMEGA>();
auto x97 = ((x93)*(x93));
auto x98 = ((x94)*(x94));
auto x99 = ((x95)*(x95));
auto x0 = V{1} - x96;
auto x1 = V{3};
auto x100 = x1*x93;
auto x101 = x1*x94;
auto x102 = x1*x95;
auto x2 = x97 + x98 + x99;
auto x3 = V{0.5}*x1*x2;
auto x105 = V{1} - x3;
auto x4 = V{0.333333333333333};
auto x5 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{1};
auto x6 = -x100*x90;
auto x7 = -x101*x91;
auto x8 = x6 + x7;
auto x9 = -x102*x92;
auto x104 = x7 + x9;
auto x106 = x6 + x9;
auto x10 = -x5*(V{0.5}*x96 + V{-1});
auto x11 = cell[0]*x0 + x10*x4*(x8 + x9) - x4*x96*(x5*(x3 + V{-1}) + V{1});
auto x12 = -V{1}*x93;
auto x13 = ((x1)*(x1));
auto x103 = V{1}*x13;
auto x14 = V{0.5}*x13;
auto x15 = V{0.0555555555555556};
auto x16 = x15*x5;
auto x17 = -x1*(x93 + V{1});
auto x18 = x10*x15;
auto x19 = cell[1]*x0 + x18*(x104 - x90*(x103*x12 - x17)) - x96*(x15 - x16*(x1*x12 + x105 + x14*((x12)*(x12))));
auto x20 = -V{1}*x94;
auto x21 = -x1*(x94 + V{1});
auto x22 = cell[2]*x0 + x18*(x106 - x91*(x103*x20 - x21)) - x96*(x15 - x16*(x1*x20 + x105 + x14*((x20)*(x20))));
auto x23 = -V{1}*x95;
auto x24 = -x1*(x95 + V{1});
auto x25 = cell[3]*x0 + x18*(x8 - x92*(x103*x23 - x24)) - x96*(x15 - x16*(x1*x23 + x105 + x14*((x23)*(x23))));
auto x26 = x12 + x20;
auto x27 = V{0.0277777777777778};
auto x28 = x27*x5;
auto x29 = -x103*x26;
auto x30 = x10*x27;
auto x31 = cell[4]*x0 + x30*(x9 + x90*(x17 + x29) + x91*(x21 + x29)) - x96*(x27 - x28*(x1*x26 + x105 + x14*((x26)*(x26))));
auto x32 = x12 + x94;
auto x33 = x13*x32;
auto x34 = -x1*(x94 + V{-1});
auto x35 = cell[5]*x0 + x30*(x9 + x90*(x17 - V{1}*x33) + x91*(x33 + x34)) - x96*(x27 - x28*(x1*x32 + x105 + x14*((x32)*(x32))));
auto x36 = x12 + x23;
auto x37 = -x103*x36;
auto x38 = cell[6]*x0 + x30*(x7 + x90*(x17 + x37) + x92*(x24 + x37)) - x96*(x27 - x28*(x1*x36 + x105 + x14*((x36)*(x36))));
auto x39 = x12 + x95;
auto x40 = x13*x39;
auto x41 = -x1*(x95 + V{-1});
auto x42 = cell[7]*x0 + x30*(x7 + x90*(x17 - V{1}*x40) + x92*(x40 + x41)) - x96*(x27 - x28*(x1*x39 + x105 + x14*((x39)*(x39))));
auto x43 = x20 + x23;
auto x44 = -x103*x43;
auto x45 = cell[8]*x0 + x30*(x6 + x91*(x21 + x44) + x92*(x24 + x44)) - x96*(x27 - x28*(x1*x43 + x105 + x14*((x43)*(x43))));
auto x46 = x20 + x95;
auto x47 = x13*x46;
auto x48 = cell[9]*x0 + x30*(x6 + x91*(x21 - V{1}*x47) + x92*(x41 + x47)) - x96*(x27 - x28*(x1*x46 + x105 + x14*((x46)*(x46))));
auto x49 = -x1*(x93 + V{-1});
auto x50 = cell[10]*x0 + x18*(x104 + x90*(x13*x93 + x49)) - x96*(x15 - x16*(x100 + x105 + x14*x97));
auto x51 = cell[11]*x0 + x18*(x106 + x91*(x13*x94 + x34)) - x96*(x15 - x16*(x101 + x105 + x14*x98));
auto x52 = cell[12]*x0 + x18*(x8 + x92*(x13*x95 + x41)) - x96*(x15 - x16*(x102 + x105 + x14*x99));
auto x53 = x93 + x94;
auto x54 = x13*x53;
auto x55 = cell[13]*x0 + x30*(x9 + x90*(x49 + x54) + x91*(x34 + x54)) - x96*(x27 - x28*(x1*x53 + x105 + x14*((x53)*(x53))));
auto x56 = x20 + x93;
auto x57 = x13*x56;
auto x58 = cell[14]*x0 + x30*(x9 + x90*(x49 + x57) + x91*(x21 - V{1}*x57)) - x96*(x27 - x28*(x1*x56 + x105 + x14*((x56)*(x56))));
auto x59 = x93 + x95;
auto x60 = x13*x59;
auto x61 = cell[15]*x0 + x30*(x7 + x90*(x49 + x60) + x92*(x41 + x60)) - x96*(x27 - x28*(x1*x59 + x105 + x14*((x59)*(x59))));
auto x62 = x23 + x93;
auto x63 = x13*x62;
auto x64 = cell[16]*x0 + x30*(x7 + x90*(x49 + x63) + x92*(x24 - V{1}*x63)) - x96*(x27 - x28*(x1*x62 + x105 + x14*((x62)*(x62))));
auto x65 = x94 + x95;
auto x66 = x13*x65;
auto x67 = cell[17]*x0 + x30*(x6 + x91*(x34 + x66) + x92*(x41 + x66)) - x96*(x27 - x28*(x1*x65 + x105 + x14*((x65)*(x65))));
auto x68 = x23 + x94;
auto x69 = x13*x68;
auto x70 = cell[18]*x0 + x30*(x6 + x91*(x34 + x69) + x92*(x24 - V{1}*x69)) - x96*(x27 - x28*(x1*x68 + x105 + x14*((x68)*(x68))));
cell[0] = x11;
cell[10] = x50;
cell[11] = x51;
cell[12] = x52;
cell[13] = x55;
cell[14] = x58;
cell[15] = x61;
cell[16] = x64;
cell[17] = x67;
cell[18] = x70;
cell[1] = x19;
cell[2] = x22;
cell[3] = x25;
cell[4] = x31;
cell[5] = x35;
cell[6] = x38;
cell[7] = x42;
cell[8] = x45;
cell[9] = x48;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x11 + x19 + x22 + x25 + x31 + x35 + x38 + x42 + x45 + x48 + x50 + x51 + x52 + x55 + x58 + x61 + x64 + x67 + x70 + V{1};
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x93;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x94;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x95;
return { x5, x2 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
