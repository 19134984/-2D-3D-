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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FreeEnergyMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::FreeEnergy, collision::FreeEnergy, dynamics::DefaultCombination>, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x67 = cell.template getFieldComponent<olb::descriptors::CHEM_POTENTIAL>(0);
auto x68 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x69 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x70 = cell.template getFieldComponent<olb::descriptors::FORCE>(2);
auto x71 = parameters.template get<collision::FreeEnergy::GAMMA>();
auto x72 = parameters.template get<descriptors::OMEGA>();
auto x73 = ((x68)*(x68));
auto x74 = ((x69)*(x69));
auto x75 = ((x70)*(x70));
auto x0 = V{1} - x72;
auto x1 = V{3};
auto x2 = x73 + x74 + x75;
auto x3 = V{0.5}*x1*x2;
auto x77 = V{1} - x3;
auto x4 = V{0.333333333333333};
auto x5 = cell[0] + cell[10] + cell[11] + cell[12] + cell[13] + cell[14] + cell[15] + cell[16] + cell[17] + cell[18] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + cell[9] + V{1};
auto x6 = -x1*x67*x71 + x5;
auto x76 = x6*x72;
auto x7 = cell[0]*x0 - x4*x72*(x5*(x3 + V{-1}) + V{1}) - x76*(x4 + V{-1});
auto x8 = -V{1}*x68;
auto x9 = V{0.5}*((x1)*(x1));
auto x10 = V{0.0555555555555556};
auto x11 = x10*x5;
auto x12 = x10*x76;
auto x13 = -(-cell[1]*x0 + x12 + x72*(x10 - x11*(x1*x8 + x77 + x9*((x8)*(x8)))));
auto x14 = -V{1}*x69;
auto x15 = -(-cell[2]*x0 + x12 + x72*(x10 - x11*(x1*x14 + x77 + x9*((x14)*(x14)))));
auto x16 = -V{1}*x70;
auto x17 = -(-cell[3]*x0 + x12 + x72*(x10 - x11*(x1*x16 + x77 + x9*((x16)*(x16)))));
auto x18 = x14 + x8;
auto x19 = V{0.0277777777777778};
auto x20 = x19*x5;
auto x21 = x19*x76;
auto x22 = -(-cell[4]*x0 + x21 + x72*(x19 - x20*(x1*x18 + x77 + x9*((x18)*(x18)))));
auto x23 = x69 + x8;
auto x24 = -(-cell[5]*x0 + x21 + x72*(x19 - x20*(x1*x23 + x77 + x9*((x23)*(x23)))));
auto x25 = x16 + x8;
auto x26 = -(-cell[6]*x0 + x21 + x72*(x19 - x20*(x1*x25 + x77 + x9*((x25)*(x25)))));
auto x27 = x70 + x8;
auto x28 = -(-cell[7]*x0 + x21 + x72*(x19 - x20*(x1*x27 + x77 + x9*((x27)*(x27)))));
auto x29 = x14 + x16;
auto x30 = -(-cell[8]*x0 + x21 + x72*(x19 - x20*(x1*x29 + x77 + x9*((x29)*(x29)))));
auto x31 = x14 + x70;
auto x32 = -(-cell[9]*x0 + x21 + x72*(x19 - x20*(x1*x31 + x77 + x9*((x31)*(x31)))));
auto x33 = cell[10]*x0 - x12 - x72*(x10 - x11*(x1*x68 + x73*x9 + x77));
auto x34 = cell[11]*x0 - x12 - x72*(x10 - x11*(x1*x69 + x74*x9 + x77));
auto x35 = cell[12]*x0 - x12 - x72*(x10 - x11*(x1*x70 + x75*x9 + x77));
auto x36 = x68 + x69;
auto x37 = -(-cell[13]*x0 + x21 + x72*(x19 - x20*(x1*x36 + x77 + x9*((x36)*(x36)))));
auto x38 = x14 + x68;
auto x39 = -(-cell[14]*x0 + x21 + x72*(x19 - x20*(x1*x38 + x77 + x9*((x38)*(x38)))));
auto x40 = x68 + x70;
auto x41 = -(-cell[15]*x0 + x21 + x72*(x19 - x20*(x1*x40 + x77 + x9*((x40)*(x40)))));
auto x42 = x16 + x68;
auto x43 = -(-cell[16]*x0 + x21 + x72*(x19 - x20*(x1*x42 + x77 + x9*((x42)*(x42)))));
auto x44 = x69 + x70;
auto x45 = -(-cell[17]*x0 + x21 + x72*(x19 - x20*(x1*x44 + x77 + x9*((x44)*(x44)))));
auto x46 = x16 + x69;
auto x47 = -(-cell[18]*x0 + x21 + x72*(x19 - x20*(x1*x46 + x77 + x9*((x46)*(x46)))));
cell[0] = x7;
cell[10] = x33;
cell[11] = x34;
cell[12] = x35;
cell[13] = x37;
cell[14] = x39;
cell[15] = x41;
cell[16] = x43;
cell[17] = x45;
cell[18] = x47;
cell[1] = x13;
cell[2] = x15;
cell[3] = x17;
cell[4] = x22;
cell[5] = x24;
cell[6] = x26;
cell[7] = x28;
cell[8] = x30;
cell[9] = x32;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x13 + x15 + x17 + x22 + x24 + x26 + x28 + x30 + x32 + x33 + x34 + x35 + x37 + x39 + x41 + x43 + x45 + x47 + x7 + V{1};
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x68;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x69;
cell.template getFieldPointer<olb::fields::moments::U>()[2] = x70;
return { x5, x2 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
