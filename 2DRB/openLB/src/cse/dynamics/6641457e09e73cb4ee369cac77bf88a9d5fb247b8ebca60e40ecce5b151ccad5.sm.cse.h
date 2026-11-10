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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::ParameterFromCell<collision::LES::SMAGORINSKY, collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB> > >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x89 = cell.template getFieldComponent<olb::collision::LES::SMAGORINSKY>(0);
auto x91 = parameters.template get<descriptors::OMEGA>();
auto x90 = cell[1] + cell[7];
auto x92 = cell[8] + x90;
auto x93 = cell[4] + x92;
auto x94 = cell[5] + cell[6];
auto x95 = cell[2] + x94;
auto x0 = V{3};
auto x100 = ((x0)*(x0));
auto x101 = V{1}/x0;
auto x1 = -V{1}*cell[1];
auto x97 = cell[7] + x1;
auto x2 = -V{1}*cell[3];
auto x3 = cell[0] + cell[3] + x93 + x95 + V{1};
auto x96 = V{1} / (x3);
auto x4 = x96*(-V{1}*cell[2] + x2 + x94 + x97);
auto x98 = ((x4)*(x4));
auto x5 = x98;
auto x6 = -V{1}*cell[5];
auto x7 = x96*(-V{1}*cell[4] + x2 + x6 + x92);
auto x99 = ((x7)*(x7));
auto x8 = x5 + x99;
auto x9 = V{0.5}*x0*x8;
auto x113 = V{1} - x9;
auto x10 = V{0.444444444444444};
auto x11 = x100;
auto x12 = V{0.5}*x11;
auto x13 = x101;
auto x103 = V{1} - x13;
auto x14 = -x13;
auto x15 = x13*(x3 + V{-1});
auto x102 = cell[3] - x15;
auto x16 = cell[5] + x102 - x3*x99 + x93;
auto x17 = x14*x16;
auto x18 = x3*x4;
auto x19 = x102 - x18*x4 + x90 + x95;
auto x20 = x14*x19;
auto x21 = V{1}/x91;
auto x22 = cell[3] - x18*x7 + x6 + x97;
auto x104 = V{2}*x22;
auto x105 = x11*x22;
auto x23 = V{1} - V{2}*V{1} / (x21 + util::sqrt(V{2.828428}*x100*x96*((x89)*(x89))*util::sqrt(((x16)*(x16)) + ((x19)*(x19)) + V{2}*((x22)*(x22))) + ((x21)*(x21))));
auto x24 = x10*(x12*x23*(x17 + x20) - x3*(x9 + V{-1}) + V{-1});
auto x25 = x0*x11;
auto x26 = V{0.166667}*x25;
auto x117 = x14*x26;
auto x27 = x103;
auto x28 = -x101;
auto x29 = -x28 + V{-1};
auto x30 = x5*x7;
auto x31 = x4*x99;
auto x32 = x30 - x31;
auto x111 = x26*x32;
auto x33 = V{0.5}*x25;
auto x34 = x30 + x31;
auto x112 = x33*x34;
auto x35 = -V{1}*x4;
auto x36 = x35 + x7;
auto x37 = V{0.0277777777777778};
auto x107 = -x37;
auto x110 = x23*x37;
auto x38 = x3*x37;
auto x39 = -V{1}*x27;
auto x108 = -x39;
auto x40 = x104*x4 + x19*x7;
auto x41 = x104*x7 + x16*x4;
auto x42 = x40 - x41;
auto x109 = x26*x42*(x108 + x27);
auto x43 = x40 + x41;
auto x106 = x33*x43;
auto x44 = x16*x27;
auto x45 = x19*x27;
auto x46 = x12*(x44 + x45);
auto x47 = -V{1}*x105 + x106*(x27 + x39) + x46;
auto x48 = x107 + x110*(x109 + x47) + x38*(x0*x36 + x111*(x27 - x29) + x112*(x27 + x29) + x113 + x12*((x36)*(x36)));
auto x49 = -x28;
auto x50 = V{0.111111111111111};
auto x114 = -x50;
auto x116 = x23*x50;
auto x51 = x3*x50;
auto x52 = -V{1}*x14;
auto x115 = x26*x52;
auto x53 = x33*x52;
auto x54 = x43*x53;
auto x55 = x12*(x17 + x45);
auto x56 = x114 + x116*(-x115*x42 + x54 + x55) + x51*(x0*x35 - x111*x49 + x112*x49 + x113 + x12*((x35)*(x35)));
auto x57 = -V{1}*x7;
auto x58 = x35 + x57;
auto x59 = x105 + x46;
auto x60 = x107 + x110*(V{2}*x106*x39 + x59) + x38*(x0*x58 - x111*(x108 + x29) + x112*(x29 + x39) + x113 + x12*((x58)*(x58)));
auto x61 = x115;
auto x62 = x12*(x20 + x44);
auto x63 = x114 + x116*(x42*x61 + x54 + x62) + x51*(x0*x57 + x113 + x12*((x57)*(x57)) + x32*x61 + x34*x53);
auto x64 = x103;
auto x65 = x4 + x57;
auto x66 = -(x110*(x109 - x47) + x37 - x38*(x0*x65 + x111*(x39 - x64) + x112*(x39 + x64) + x113 + x12*((x65)*(x65))));
auto x67 = -x117;
auto x68 = x14*x33;
auto x69 = x34*x68;
auto x118 = x113 + x69;
auto x70 = x43*x68;
auto x71 = x114 + x116*(x42*x67 + x55 + x70) + x51*(x0*x4 + x118 + x12*x98 + x32*x67);
auto x119 = x66 + x71;
auto x72 = 0;
auto x73 = V{2}*x27*x33;
auto x74 = x4 + x7;
auto x75 = x107 + x110*(x42*x72 + x43*x73 + x59) + x38*(x0*x74 + x113 + x12*((x74)*(x74)) + x32*x72 + x34*x73);
auto x76 = x117;
auto x77 = x114 + x116*(x42*x76 + x62 + x70) + x51*(x0*x7 + x118 + x12*x99 + x32*x76);
auto x120 = x48 + x77;
auto x78 = x119 + x120 + x24 + x56 + x60 + x63 + x75 + V{1};
auto x121 = V{1} / (x78);
auto x79 = -V{1}*x60;
auto x122 = x75 + x79;
cell[0] = x24;
cell[1] = x48;
cell[2] = x56;
cell[3] = x60;
cell[4] = x63;
cell[5] = x66;
cell[6] = x71;
cell[7] = x75;
cell[8] = x77;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x78;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x121*(x119 + x122 - V{1}*x48 - V{1}*x56);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x121*(x120 + x122 - V{1}*x63 - V{1}*x66);
return { x3, x8 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
