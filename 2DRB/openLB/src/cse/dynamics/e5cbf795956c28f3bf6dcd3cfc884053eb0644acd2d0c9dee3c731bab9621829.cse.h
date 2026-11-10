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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, guoZhao::GuoZhaoSecondOrder, collision::BGK, guoZhao::GuoZhaoForcing<momenta::GuoZhaoForced> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x47 = cell.template getFieldComponent<olb::descriptors::BODY_FORCE>(0);
auto x48 = cell.template getFieldComponent<olb::descriptors::BODY_FORCE>(1);
auto x49 = cell.template getFieldComponent<olb::descriptors::EPSILON>(0);
auto x52 = cell.template getFieldComponent<olb::descriptors::K>(0);
auto x53 = cell.template getFieldComponent<olb::descriptors::NU>(0);
auto x54 = parameters.template get<descriptors::OMEGA>();
auto x50 = cell[5] + cell[6];
auto x51 = cell[1] + cell[8];
auto x55 = x53/x52;
auto x57 = V{0.5}*x49;
auto x62 = V{1} / (x49);
auto x0 = V{1} - x54;
auto x1 = V{3};
auto x63 = x1*x62;
auto x2 = -V{1}*cell[3];
auto x59 = cell[7] + x2;
auto x3 = cell[0] + cell[2] + cell[3] + cell[4] + cell[7] + x50 + x51 + V{1};
auto x58 = V{1} / (x3);
auto x4 = V{0.25}*x49*x55 + V{0.5};
auto x5 = x4 + util::sqrt(((x4)*(x4)));
auto x56 = V{1} / (x5);
auto x6 = x56*(x47*x57 + x58*(-V{1}*cell[1] - V{1}*cell[2] + x50 + x59));
auto x60 = ((x6)*(x6));
auto x64 = x1*x6;
auto x67 = -x6;
auto x7 = x56*(x48*x57 + x58*(-V{1}*cell[4] - V{1}*cell[5] + x51 + x59));
auto x61 = ((x7)*(x7));
auto x65 = x1*x7;
auto x68 = -x7;
auto x8 = x60 + x61;
auto x9 = V{0.5}*x63*x8;
auto x74 = V{1} - x9;
auto x10 = V{0.444444444444444};
auto x11 = -V{1}*x6;
auto x12 = x49*(x11*x55 + x47);
auto x13 = -x12*x62*x64;
auto x14 = -V{1}*x7;
auto x15 = x49*(x14*x55 + x48);
auto x16 = -x15*x62*x65;
auto x17 = V{1} - V{0.5}*x54;
auto x18 = x11 + x7;
auto x19 = ((x1)*(x1));
auto x20 = V{0.5}*x19;
auto x73 = x20*x62;
auto x21 = V{0.0277777777777778};
auto x66 = x21*x3;
auto x72 = x17*x66;
auto x22 = x66;
auto x23 = x19*x62;
auto x70 = V{1}*x23;
auto x24 = x18*x23;
auto x25 = -V{1}*x49;
auto x26 = x63*(x25 + x67);
auto x71 = -x26;
auto x27 = x63*(x49 + x68);
auto x28 = V{0.111111111111111};
auto x69 = x28*x3;
auto x75 = x17*x69;
auto x29 = x69;
auto x30 = x11 + x14;
auto x31 = -x30*x70;
auto x32 = x63*(x25 + x68);
auto x33 = x14 + x6;
auto x34 = x23*x33;
auto x35 = x63*(x49 + x67);
auto x36 = x6 + x7;
auto x37 = x23*x36;
auto x38 = cell[0]*x0 + x10*x17*x3*(x13 + x16) - x10*x54*(x3*(x9 + V{-1}) + V{1});
auto x39 = -(-cell[1]*x0 + x54*(x21 - x22*(x1*x18 + x73*((x18)*(x18)) + x74)) + x72*(x12*(V{1}*x24 + x71) - x15*(x24 + x27)));
auto x40 = -(-cell[2]*x0 + x54*(x28 - x29*(x1*x11 + x73*((x11)*(x11)) + x74)) + x75*(x12*(x11*x70 + x71) - x16));
auto x41 = cell[3]*x0 - x54*(x21 - x22*(x1*x30 + x73*((x30)*(x30)) + x74)) + x72*(x12*(x26 + x31) + x15*(x31 + x32));
auto x42 = cell[4]*x0 - x54*(x28 - x29*(x1*x14 + x73*((x14)*(x14)) + x74)) + x75*(x13 - x15*(x14*x70 - x32));
auto x43 = cell[5]*x0 - x54*(x21 - x22*(x1*x33 + x73*((x33)*(x33)) + x74)) + x72*(x12*(x34 + x35) + x15*(x32 - V{1}*x34));
auto x44 = cell[6]*x0 - x54*(x28 - x29*(x60*x73 + x64 + x74)) + x75*(x12*(x23*x6 + x35) + x16);
auto x45 = cell[7]*x0 - x54*(x21 - x22*(x1*x36 + x73*((x36)*(x36)) + x74)) + x72*(x12*(x35 + x37) + x15*(x27 + x37));
auto x46 = cell[8]*x0 - x54*(x28 - x29*(x61*x73 + x65 + x74)) + x75*(x13 + x15*(x23*x7 + x27));
cell[0] = x38;
cell[1] = x39;
cell[2] = x40;
cell[3] = x41;
cell[4] = x42;
cell[5] = x43;
cell[6] = x44;
cell[7] = x45;
cell[8] = x46;
return { x3, x8 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
