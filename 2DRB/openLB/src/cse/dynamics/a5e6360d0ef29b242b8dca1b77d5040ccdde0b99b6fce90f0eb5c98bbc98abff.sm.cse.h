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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::FixedVelocityMomentum, momenta::BulkStress, momenta::DefineUSeparately>, equilibria::SecondOrder, collision::BGK, forcing::MCGuo<momenta::Identity> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x49 = cell.template getFieldComponent<olb::descriptors::FORCE>(0);
auto x50 = cell.template getFieldComponent<olb::descriptors::FORCE>(1);
auto x51 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x52 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x53 = parameters.template get<descriptors::OMEGA>();
auto x54 = ((x51)*(x51));
auto x55 = ((x52)*(x52));
auto x0 = V{1} - x53;
auto x1 = V{3};
auto x56 = x1*x51;
auto x57 = x1*x52;
auto x2 = x54 + x55;
auto x3 = V{0.5}*x1*x2;
auto x59 = V{1} - x3;
auto x4 = V{0.444444444444444};
auto x5 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8] + V{1};
auto x6 = -x49*x56;
auto x7 = -x50*x57;
auto x8 = -x5*(V{0.5}*x53 + V{-1});
auto x9 = cell[0]*x0 - x4*x53*(x5*(x3 + V{-1}) + V{1}) + x4*x8*(x6 + x7);
auto x10 = -V{1}*x51;
auto x11 = x10 + x52;
auto x12 = ((x1)*(x1));
auto x60 = V{1}*x12;
auto x13 = V{0.5}*x12;
auto x14 = V{0.0277777777777778};
auto x15 = x14*x5;
auto x16 = x11*x12;
auto x17 = -x1*(x51 + V{1});
auto x58 = -x17;
auto x18 = -x1*(x52 + V{-1});
auto x19 = x14*x8;
auto x20 = -(-cell[1]*x0 + x19*(x49*(V{1}*x16 + x58) - x50*(x16 + x18)) + x53*(x14 - x15*(x1*x11 + x13*((x11)*(x11)) + x59)));
auto x21 = V{0.111111111111111};
auto x22 = x21*x5;
auto x23 = x21*x8;
auto x24 = -(-cell[2]*x0 + x23*(x49*(x10*x60 + x58) - x7) + x53*(x21 - x22*(x1*x10 + x13*((x10)*(x10)) + x59)));
auto x25 = -V{1}*x52;
auto x26 = x10 + x25;
auto x27 = -x26*x60;
auto x28 = -x1*(x52 + V{1});
auto x29 = cell[3]*x0 + x19*(x49*(x17 + x27) + x50*(x27 + x28)) - x53*(x14 - x15*(x1*x26 + x13*((x26)*(x26)) + x59));
auto x30 = -(-cell[4]*x0 + x23*(x50*(x25*x60 - x28) - x6) + x53*(x21 - x22*(x1*x25 + x13*((x25)*(x25)) + x59)));
auto x31 = x25 + x51;
auto x32 = x12*x31;
auto x33 = -x1*(x51 + V{-1});
auto x34 = cell[5]*x0 + x19*(x49*(x32 + x33) + x50*(x28 - V{1}*x32)) - x53*(x14 - x15*(x1*x31 + x13*((x31)*(x31)) + x59));
auto x35 = cell[6]*x0 + x23*(x49*(x12*x51 + x33) + x7) - x53*(x21 - x22*(x13*x54 + x56 + x59));
auto x36 = x51 + x52;
auto x37 = x12*x36;
auto x38 = cell[7]*x0 + x19*(x49*(x33 + x37) + x50*(x18 + x37)) - x53*(x14 - x15*(x1*x36 + x13*((x36)*(x36)) + x59));
auto x39 = cell[8]*x0 + x23*(x50*(x12*x52 + x18) + x6) - x53*(x21 - x22*(x13*x55 + x57 + x59));
cell[0] = x9;
cell[1] = x20;
cell[2] = x24;
cell[3] = x29;
cell[4] = x30;
cell[5] = x34;
cell[6] = x35;
cell[7] = x38;
cell[8] = x39;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x20 + x24 + x29 + x30 + x34 + x35 + x38 + x39 + x9 + V{1};
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x51;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x52;
return { x5, x2 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
