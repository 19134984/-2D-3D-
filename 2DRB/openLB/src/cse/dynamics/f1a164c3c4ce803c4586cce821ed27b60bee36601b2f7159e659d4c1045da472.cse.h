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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<0, 1>, momenta::BulkStress, momenta::DefineSeparately>, equilibria::ThirdOrder, collision::ParameterFromCell<collision::LES::SMAGORINSKY, collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB> >, dynamics::DefaultCombination>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x73 = cell.template getFieldComponent<olb::collision::LES::SMAGORINSKY>(0);
auto x74 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x76 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1);
auto x78 = parameters.template get<descriptors::OMEGA>();
auto x75 = V{1} / (x74);
auto x77 = cell[4] + cell[8];
auto x80 = ((x76)*(x76));
auto x0 = V{3};
auto x81 = ((x0)*(x0));
auto x82 = V{1}/x0;
auto x1 = x75*(cell[0] + V{2}*cell[5] + V{2}*cell[6] + V{2}*cell[7] + x77 + V{1}) + V{-1};
auto x79 = ((x1)*(x1));
auto x2 = x79;
auto x3 = x2 + x80;
auto x4 = V{0.5}*x0*x3;
auto x96 = V{1} - x4;
auto x5 = V{0.444444444444444};
auto x6 = x81;
auto x7 = V{0.5}*x6;
auto x8 = x82;
auto x84 = V{1} - x8;
auto x9 = -x8;
auto x10 = x8*(x74 + V{-1});
auto x83 = cell[1] + cell[3] + cell[5] + cell[7] - x10;
auto x11 = -x74*x80 + x77 + x83;
auto x12 = x11*x9;
auto x13 = x1*x74;
auto x14 = cell[2] + cell[6] - x1*x13 + x83;
auto x15 = x14*x9;
auto x16 = V{1}/x78;
auto x17 = -V{1}*cell[1] + cell[3] - V{1}*cell[5] + cell[7] - x13*x76;
auto x85 = V{2}*x17;
auto x86 = x17*x6;
auto x18 = V{1} - V{2}*V{1} / (x16 + util::sqrt(V{2.828428}*x75*x81*((x73)*(x73))*util::sqrt(((x11)*(x11)) + ((x14)*(x14)) + V{2}*((x17)*(x17))) + ((x16)*(x16))));
auto x19 = x0*x6;
auto x20 = V{0.166667}*x19;
auto x89 = x20*x9;
auto x21 = x84;
auto x22 = -x82;
auto x23 = -x22 + V{-1};
auto x24 = x2*x76;
auto x25 = x1*x80;
auto x26 = x24 - x25;
auto x94 = x20*x26;
auto x27 = V{0.5}*x19;
auto x28 = x24 + x25;
auto x95 = x27*x28;
auto x29 = -V{1}*x1;
auto x30 = x29 + x76;
auto x31 = V{0.0277777777777778};
auto x90 = -x31;
auto x93 = x18*x31;
auto x32 = x31*x74;
auto x33 = -V{1}*x21;
auto x91 = -x33;
auto x34 = x1*x85 + x14*x76;
auto x35 = x1*x11 + x76*x85;
auto x36 = x34 - x35;
auto x92 = x20*x36*(x21 + x91);
auto x37 = x34 + x35;
auto x87 = x27*x37;
auto x38 = x11*x21;
auto x39 = x14*x21;
auto x40 = x7*(x38 + x39);
auto x41 = x40 - V{1}*x86 + x87*(x21 + x33);
auto x42 = -x22;
auto x43 = V{0.111111111111111};
auto x97 = -x43;
auto x98 = x18*x43;
auto x44 = x43*x74;
auto x45 = -V{1}*x9;
auto x88 = x20*x45;
auto x46 = x27*x45;
auto x47 = x37*x46;
auto x48 = x7*(x12 + x39);
auto x49 = -V{1}*x76;
auto x50 = x29 + x49;
auto x51 = x40 + x86;
auto x52 = x88;
auto x53 = x7*(x15 + x38);
auto x54 = x84;
auto x55 = x1 + x49;
auto x56 = -x89;
auto x57 = x27*x9;
auto x58 = x28*x57;
auto x99 = x58 + x96;
auto x59 = x37*x57;
auto x60 = 0;
auto x61 = V{2}*x21*x27;
auto x62 = x1 + x76;
auto x63 = x89;
cell[0] = x5*(x18*x7*(x12 + x15) - x74*(x4 + V{-1}) + V{-1});
cell[1] = x32*(x0*x30 + x7*((x30)*(x30)) + x94*(x21 - x23) + x95*(x21 + x23) + x96) + x90 + x93*(x41 + x92);
cell[2] = x44*(x0*x29 - x42*x94 + x42*x95 + x7*((x29)*(x29)) + x96) + x97 + x98*(-x36*x88 + x47 + x48);
cell[3] = x32*(x0*x50 + x7*((x50)*(x50)) - x94*(x23 + x91) + x95*(x23 + x33) + x96) + x90 + x93*(V{2}*x33*x87 + x51);
cell[4] = x44*(x0*x49 + x26*x52 + x28*x46 + x7*((x49)*(x49)) + x96) + x97 + x98*(x36*x52 + x47 + x53);
cell[5] = -(x31 - x32*(x0*x55 + x7*((x55)*(x55)) + x94*(x33 - x54) + x95*(x33 + x54) + x96) + x93*(-x41 + x92));
cell[6] = x44*(x0*x1 + x26*x56 + x7*x79 + x99) + x97 + x98*(x36*x56 + x48 + x59);
cell[7] = x32*(x0*x62 + x26*x60 + x28*x61 + x7*((x62)*(x62)) + x96) + x90 + x93*(x36*x60 + x37*x61 + x51);
cell[8] = x44*(x0*x76 + x26*x63 + x7*x80 + x99) + x97 + x98*(x36*x63 + x53 + x59);
return { x74, x3 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
