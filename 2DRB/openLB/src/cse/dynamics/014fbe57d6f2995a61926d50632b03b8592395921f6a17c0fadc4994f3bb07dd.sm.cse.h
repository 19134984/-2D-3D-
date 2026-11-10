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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, forcing::Guo<momenta::ForcedWithStress> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x70 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x71 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x72 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x73 = parameters.template get<descriptors::OMEGA>();
auto x74 = cell[3] + cell[7];
auto x75 = cell[1] + cell[8];
auto x76 = cell[4] + x75;
auto x77 = cell[5] + cell[6];
auto x78 = cell[2] + x77;
auto x0 = V{1}/x73;
auto x1 = -V{1}*cell[3];
auto x80 = cell[7] + x1;
auto x2 = -V{1}*cell[5];
auto x3 = cell[0] + x74 + x76 + x78 + V{1};
auto x79 = V{1} / (x3);
auto x4 = V{0.5}*x71;
auto x5 = x4 + x79*(-V{1}*cell[4] + x2 + x75 + x80);
auto x85 = ((x5)*(x5));
auto x6 = -x4 + x5;
auto x7 = V{3};
auto x83 = ((x7)*(x7));
auto x87 = x5*x7;
auto x8 = V{1}*(x3 + V{-1})/x7;
auto x82 = x74 - x8;
auto x9 = V{0.5}*x3;
auto x81 = V{2}*x9;
auto x10 = x6*x71;
auto x11 = cell[5] + x10*x81 - x3*((x6)*(x6)) + x76 + x82;
auto x12 = -V{1}*cell[1];
auto x13 = V{0.5}*x70;
auto x14 = x13 + x79*(-V{1}*cell[2] + x12 + x77 + x80);
auto x84 = ((x14)*(x14));
auto x86 = x14*x7;
auto x15 = -x13 + x14;
auto x16 = x15*x3;
auto x17 = x12 - x16*x6 + x2 + x74 + x9*(x15*x71 + x6*x70);
auto x18 = x15*x70;
auto x19 = cell[1] - x15*x16 + x18*x81 + x78 + x82;
auto x20 = V{2}*V{1} / (x0 + util::sqrt(V{2.828428}*x79*x83*((x72)*(x72))*util::sqrt(((x11)*(x11)) + V{2}*((x17)*(x17)) + ((x19)*(x19))) + ((x0)*(x0))));
auto x21 = V{1} - x20;
auto x22 = x84 + x85;
auto x23 = V{0.5}*x22*x7;
auto x89 = V{1} - x23;
auto x24 = V{0.444444444444444};
auto x25 = -x70*x86;
auto x26 = -x71*x87;
auto x27 = -x3*(V{0.5}*x20 + V{-1});
auto x28 = cell[0]*x21 - x20*x24*(x3*(x23 + V{-1}) + V{1}) + x24*x27*(x25 + x26);
auto x29 = -V{1}*x14;
auto x30 = x29 + x5;
auto x31 = x83;
auto x90 = V{1}*x31;
auto x32 = V{0.5}*x31;
auto x33 = V{0.0277777777777778};
auto x34 = x3*x33;
auto x35 = x30*x31;
auto x36 = -x7*(x14 + V{1});
auto x88 = -x36;
auto x37 = -x7*(x5 + V{-1});
auto x38 = x27*x33;
auto x39 = -(-cell[1]*x21 + x20*(x33 - x34*(x30*x7 + x32*((x30)*(x30)) + x89)) + x38*(x70*(V{1}*x35 + x88) - x71*(x35 + x37)));
auto x40 = V{0.111111111111111};
auto x41 = x3*x40;
auto x42 = x27*x40;
auto x43 = cell[2]*x21 - x20*(x40 - x41*(x29*x7 + x32*((x29)*(x29)) + x89)) + x42*(x26 - x70*(x29*x90 + x88));
auto x44 = -V{1}*x5;
auto x45 = x29 + x44;
auto x46 = -x45*x90;
auto x47 = -x7*(x5 + V{1});
auto x48 = cell[3]*x21 - x20*(x33 - x34*(x32*((x45)*(x45)) + x45*x7 + x89)) + x38*(x70*(x36 + x46) + x71*(x46 + x47));
auto x49 = cell[4]*x21 - x20*(x40 - x41*(x32*((x44)*(x44)) + x44*x7 + x89)) + x42*(x25 - x71*(x44*x90 - x47));
auto x50 = x14 + x44;
auto x51 = x31*x50;
auto x52 = -x7*(x14 + V{-1});
auto x53 = cell[5]*x21 - x20*(x33 - x34*(x32*((x50)*(x50)) + x50*x7 + x89)) + x38*(x70*(x51 + x52) + x71*(x47 - V{1}*x51));
auto x54 = cell[6]*x21 - x20*(x40 - x41*(x32*x84 + x86 + x89)) + x42*(x26 + x70*(x14*x31 + x52));
auto x91 = x53 + x54;
auto x55 = x14 + x5;
auto x56 = x31*x55;
auto x57 = cell[7]*x21 - x20*(x33 - x34*(x32*((x55)*(x55)) + x55*x7 + x89)) + x38*(x70*(x52 + x56) + x71*(x37 + x56));
auto x58 = cell[8]*x21 - x20*(x40 - x41*(x32*x85 + x87 + x89)) + x42*(x25 + x71*(x31*x5 + x37));
auto x92 = x39 + x58;
auto x59 = x28 + x43 + x48 + x49 + x57 + x91 + x92 + V{1};
auto x93 = V{1} / (x59);
auto x60 = -V{1}*x48;
auto x94 = x57 + x60;
cell[0] = x28;
cell[1] = x39;
cell[2] = x43;
cell[3] = x48;
cell[4] = x49;
cell[5] = x53;
cell[6] = x54;
cell[7] = x57;
cell[8] = x58;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x59;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x13 + x93*(-V{1}*x39 - V{1}*x43 + x91 + x94);
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x4 + x93*(-V{1}*x49 - V{1}*x53 + x92 + x94);
return { x3, x22 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
