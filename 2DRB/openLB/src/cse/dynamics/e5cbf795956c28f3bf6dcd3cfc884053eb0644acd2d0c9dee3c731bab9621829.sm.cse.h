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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, guoZhao::GuoZhaoSecondOrder, collision::BGK, guoZhao::GuoZhaoForcing<momenta::GuoZhaoForced> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x60 = cell.template getFieldComponent<olb::descriptors::BODY_FORCE>(0);
auto x61 = cell.template getFieldComponent<olb::descriptors::BODY_FORCE>(1);
auto x62 = cell.template getFieldComponent<olb::descriptors::EPSILON>(0);
auto x65 = cell.template getFieldComponent<olb::descriptors::K>(0);
auto x66 = cell.template getFieldComponent<olb::descriptors::NU>(0);
auto x67 = parameters.template get<descriptors::OMEGA>();
auto x63 = cell[5] + cell[6];
auto x64 = cell[1] + cell[8];
auto x68 = V{0.5}*x62;
auto x69 = x66/x65;
auto x75 = V{1} / (x62);
auto x0 = V{1} - x67;
auto x1 = V{3};
auto x76 = x1*x75;
auto x2 = -V{1}*cell[3];
auto x72 = cell[7] + x2;
auto x3 = cell[0] + cell[2] + cell[3] + cell[4] + cell[7] + x63 + x64 + V{1};
auto x71 = V{1} / (x3);
auto x4 = x60*x68;
auto x5 = V{0.25}*x62*x69 + V{0.5};
auto x6 = x5 + util::sqrt(((x5)*(x5)));
auto x70 = V{1} / (x6);
auto x7 = x70*(x4 + x71*(-V{1}*cell[1] - V{1}*cell[2] + x63 + x72));
auto x73 = ((x7)*(x7));
auto x77 = x1*x7;
auto x80 = -x7;
auto x8 = x61*x68;
auto x9 = x70*(x71*(-V{1}*cell[4] - V{1}*cell[5] + x64 + x72) + x8);
auto x74 = ((x9)*(x9));
auto x78 = x1*x9;
auto x81 = -x9;
auto x10 = x73 + x74;
auto x11 = V{0.5}*x10*x76;
auto x85 = V{1} - x11;
auto x12 = V{0.444444444444444};
auto x13 = -V{1}*x7;
auto x14 = x62*(x13*x69 + x60);
auto x15 = -x14*x75*x77;
auto x16 = -V{1}*x9;
auto x17 = x62*(x16*x69 + x61);
auto x18 = -x17*x75*x78;
auto x19 = V{1} - V{0.5}*x67;
auto x20 = cell[0]*x0 + x12*x19*x3*(x15 + x18) - x12*x67*(x3*(x11 + V{-1}) + V{1});
auto x21 = x13 + x9;
auto x22 = ((x1)*(x1));
auto x23 = V{0.5}*x22;
auto x84 = x23*x75;
auto x24 = V{0.0277777777777778};
auto x79 = x24*x3;
auto x83 = x19*x79;
auto x25 = x79;
auto x26 = x22*x75;
auto x87 = V{1}*x26;
auto x27 = x21*x26;
auto x28 = -V{1}*x62;
auto x29 = x76*(x28 + x80);
auto x82 = -x29;
auto x30 = x76*(x62 + x81);
auto x31 = -(-cell[1]*x0 + x67*(x24 - x25*(x1*x21 + x84*((x21)*(x21)) + x85)) + x83*(x14*(V{1}*x27 + x82) - x17*(x27 + x30)));
auto x32 = V{0.111111111111111};
auto x86 = x3*x32;
auto x88 = x19*x86;
auto x33 = x86;
auto x34 = -(-cell[2]*x0 + x67*(x32 - x33*(x1*x13 + x84*((x13)*(x13)) + x85)) + x88*(x14*(x13*x87 + x82) - x18));
auto x35 = x13 + x16;
auto x36 = -x35*x87;
auto x37 = x76*(x28 + x81);
auto x38 = cell[3]*x0 - x67*(x24 - x25*(x1*x35 + x84*((x35)*(x35)) + x85)) + x83*(x14*(x29 + x36) + x17*(x36 + x37));
auto x39 = cell[4]*x0 - x67*(x32 - x33*(x1*x16 + x84*((x16)*(x16)) + x85)) + x88*(x15 - x17*(x16*x87 - x37));
auto x40 = x16 + x7;
auto x41 = x26*x40;
auto x42 = x76*(x62 + x80);
auto x43 = cell[5]*x0 - x67*(x24 - x25*(x1*x40 + x84*((x40)*(x40)) + x85)) + x83*(x14*(x41 + x42) + x17*(x37 - V{1}*x41));
auto x44 = cell[6]*x0 - x67*(x32 - x33*(x73*x84 + x77 + x85)) + x88*(x14*(x26*x7 + x42) + x18);
auto x89 = x43 + x44;
auto x45 = x7 + x9;
auto x46 = x26*x45;
auto x47 = cell[7]*x0 - x67*(x24 - x25*(x1*x45 + x84*((x45)*(x45)) + x85)) + x83*(x14*(x42 + x46) + x17*(x30 + x46));
auto x48 = cell[8]*x0 - x67*(x32 - x33*(x74*x84 + x78 + x85)) + x88*(x15 + x17*(x26*x9 + x30));
auto x90 = x31 + x48;
auto x49 = x20 + x34 + x38 + x39 + x47 + x89 + x90 + V{1};
auto x91 = V{1} / (x49);
auto x50 = -V{1}*x38;
auto x92 = x47 + x50;
cell[0] = x20;
cell[1] = x31;
cell[2] = x34;
cell[3] = x38;
cell[4] = x39;
cell[5] = x43;
cell[6] = x44;
cell[7] = x47;
cell[8] = x48;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x49;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x70*(x4 + x91*(-V{1}*x31 - V{1}*x34 + x89 + x92));
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x70*(x8 + x91*(-V{1}*x39 - V{1}*x43 + x90 + x92));
return { x3, x10 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
