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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<0, 1>, momenta::BulkStress, momenta::DefineSeparately>, equilibria::ThirdOrder, collision::ParameterFromCell<collision::LES::SMAGORINSKY, collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB> >, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x79 = cell.template getFieldComponent<olb::collision::LES::SMAGORINSKY>(0);
auto x80 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x82 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1);
auto x84 = parameters.template get<descriptors::OMEGA>();
auto x81 = V{1} / (x80);
auto x83 = cell[4] + cell[8];
auto x86 = ((x82)*(x82));
auto x0 = V{3};
auto x87 = ((x0)*(x0));
auto x88 = V{1}/x0;
auto x1 = x81*(cell[0] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7] + x83 + V{1}) + V{-1};
auto x85 = ((x1)*(x1));
auto x2 = x85;
auto x3 = x2 + x86;
auto x4 = V{0.5}*x0*x3;
auto x97 = V{1} - x4;
auto x5 = V{0.444444444444444};
auto x6 = x87;
auto x7 = V{0.5}*x6;
auto x8 = x88;
auto x90 = V{1} - x8;
auto x9 = -x8;
auto x10 = x8*(x80 + V{-1});
auto x89 = cell[1] + cell[3] + cell[5] + cell[7] - x10;
auto x11 = -x80*x86 + x83 + x89;
auto x12 = x11*x9;
auto x13 = x1*x80;
auto x14 = cell[2] + cell[6] - x1*x13 + x89;
auto x15 = x14*x9;
auto x16 = V{1}/x84;
auto x17 = -V{1}*cell[1] + cell[3] - V{1}*cell[5] + cell[7] - x13*x82;
auto x91 = V{2}*x17;
auto x92 = x17*x6;
auto x18 = V{1} - V{2}*V{1} / (x16 + util::sqrt(V{2.828428}*x81*x87*((x79)*(x79))*util::sqrt(((x11)*(x11)) + ((x14)*(x14)) + V{2}*((x17)*(x17))) + ((x16)*(x16))));
auto x19 = x5*(x18*x7*(x12 + x15) - x80*(x4 + V{-1}) + V{-1});
auto x20 = x0*x6;
auto x21 = V{0.166667}*x20;
auto x103 = x21*x9;
auto x22 = x90;
auto x23 = -x88;
auto x24 = -x23 + V{-1};
auto x25 = x2*x82;
auto x26 = x1*x86;
auto x27 = x25 - x26;
auto x101 = x21*x27;
auto x28 = V{0.5}*x20;
auto x29 = x25 + x26;
auto x102 = x28*x29;
auto x30 = -V{1}*x1;
auto x31 = x30 + x82;
auto x32 = V{0.0277777777777778};
auto x100 = x18*x32;
auto x105 = -x32;
auto x33 = x32*x80;
auto x34 = -V{1}*x22;
auto x98 = -x34;
auto x35 = x1*x91 + x14*x82;
auto x36 = x1*x11 + x82*x91;
auto x37 = x35 - x36;
auto x99 = x21*x37*(x22 + x98);
auto x38 = x35 + x36;
auto x93 = x28*x38;
auto x39 = x11*x22;
auto x40 = x14*x22;
auto x41 = x7*(x39 + x40);
auto x42 = x41 - V{1}*x92 + x93*(x22 + x34);
auto x43 = -x23;
auto x44 = V{0.111111111111111};
auto x95 = -x44;
auto x96 = x18*x44;
auto x45 = x44*x80;
auto x46 = -V{1}*x9;
auto x94 = x21*x46;
auto x47 = x28*x46;
auto x48 = x38*x47;
auto x49 = x7*(x12 + x40);
auto x50 = -V{1}*x82;
auto x51 = x30 + x50;
auto x52 = x41 + x92;
auto x53 = x94;
auto x54 = x7*(x15 + x39);
auto x55 = x45*(x0*x50 + x27*x53 + x29*x47 + x7*((x50)*(x50)) + x97) + x95 + x96*(x37*x53 + x48 + x54);
auto x56 = x90;
auto x57 = x1 + x50;
auto x58 = -(x100*(-x42 + x99) + x32 - x33*(x0*x57 + x101*(x34 - x56) + x102*(x34 + x56) + x7*((x57)*(x57)) + x97));
auto x59 = -x103;
auto x60 = x28*x9;
auto x61 = x29*x60;
auto x104 = x61 + x97;
auto x62 = x38*x60;
auto x63 = x45*(x0*x1 + x104 + x27*x59 + x7*x85) + x95 + x96*(x37*x59 + x49 + x62);
auto x64 = 0;
auto x65 = V{2}*x22*x28;
auto x66 = x1 + x82;
auto x67 = x100*(x37*x64 + x38*x65 + x52) + x105 + x33*(x0*x66 + x27*x64 + x29*x65 + x7*((x66)*(x66)) + x97);
auto x68 = x103;
auto x69 = x45*(x0*x82 + x104 + x27*x68 + x7*x86) + x95 + x96*(x37*x68 + x54 + x62);
cell[0] = x19;
cell[1] = x100*(x42 + x99) + x105 + x33*(x0*x31 + x101*(x22 - x24) + x102*(x22 + x24) + x7*((x31)*(x31)) + x97);
cell[2] = x45*(x0*x30 - x101*x43 + x102*x43 + x7*((x30)*(x30)) + x97) + x95 + x96*(-x37*x94 + x48 + x49);
cell[3] = x100*(V{2}*x34*x93 + x52) + x105 + x33*(x0*x51 - x101*(x24 + x98) + x102*(x24 + x34) + x7*((x51)*(x51)) + x97);
cell[4] = x55;
cell[5] = x58;
cell[6] = x63;
cell[7] = x67;
cell[8] = x69;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x80;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x81*(x19 + x55 + V{2}*x58 + V{2}*x63 + V{2}*x67 + x69 + V{1}) + V{-1};
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x82;
return { x80, x3 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
