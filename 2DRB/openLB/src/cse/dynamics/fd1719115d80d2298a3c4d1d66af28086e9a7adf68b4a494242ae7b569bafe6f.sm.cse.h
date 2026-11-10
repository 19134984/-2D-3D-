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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SaveVelocity<collision::BGK>, dynamics::DefaultCombination>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<0, -1>, momenta::RegularizedBoundaryStress<0, -1>, momenta::DefineSeparately> >, olb::CodegenTarget::StoredMoments> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x78 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x80 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<0, -1>::VELOCITY>(1);
auto x81 = parameters.template get<descriptors::OMEGA>();
auto x76 = V{1} / (x78);
auto x77 = cell[4] + cell[8];
auto x79 = ((x80)*(x80));
auto x0 = V{3};
auto x83 = ((x0)*(x0));
auto x1 = V{0.5}*x0;
auto x2 = -V{1}*x76*(cell[0] + V{2}*cell[1] + V{2}*cell[2] + V{2}*cell[3] + x77 + V{1}) + V{1};
auto x82 = ((x2)*(x2));
auto x3 = x79;
auto x4 = x1*(x3 + x82);
auto x85 = -x4;
auto x86 = x85 + V{1};
auto x5 = V{0.444444444444444};
auto x6 = x5*x78;
auto x7 = V{1}/x0;
auto x8 = -x7;
auto x9 = V{-2};
auto x10 = -V{1}*x9;
auto x11 = -V{1}*x2;
auto x12 = x11 + x80;
auto x13 = V{0.5}*x83;
auto x14 = V{0.0277777777777778};
auto x84 = -x14;
auto x15 = x14*x78;
auto x16 = x15*(x0*x12 + x13*((x12)*(x12)) + x86) + x84;
auto x17 = cell[1] - x16;
auto x18 = V{0.111111111111111};
auto x87 = -x18;
auto x19 = x18*x78;
auto x20 = x19*(x0*x11 + x13*((x11)*(x11)) + x86) + x87;
auto x21 = -V{1}*x80;
auto x22 = x11 + x21;
auto x23 = x15*(x0*x22 + x13*((x22)*(x22)) + x86) + x84;
auto x24 = x10*(cell[3] - x23);
auto x25 = x10*x17 + x10*(cell[2] - x20) + x24;
auto x26 = x25*x8;
auto x27 = x0*x21 + x13*((x21)*(x21)) + V{1};
auto x28 = -x18 + x19*(x27 + x85);
auto x29 = x0*x80 + x13*x79 + V{1};
auto x30 = -x18 + x19*(x29 + x85);
auto x31 = V{2}*x17 + x24 - x28 - x30 + x77;
auto x32 = x31*x8;
auto x33 = V{0.5}*x5*x83*(x26 + x32) - x5 - x6*(x4 + V{-1});
auto x34 = V{1} - x81;
auto x35 = V{0.5}*x18*x83;
auto x36 = V{1} - x7;
auto x37 = x31*x36;
auto x38 = x35*(x26 + x37);
auto x39 = x30 + x38;
auto x40 = x28 + x38;
auto x41 = V{0.5}*x14*x83;
auto x42 = x25*x36;
auto x88 = x37 + x42;
auto x43 = x17*x9 + x24;
auto x44 = x41*(V{2}*x43 + x88);
auto x45 = x23 + x44;
auto x46 = x41*(x43*x9 + x88);
auto x47 = x16 + x46;
auto x48 = x35*(x32 + x42);
auto x49 = x20 + x48;
auto x50 = -V{1}*x76*(x33 + x39 + x40 + V{2}*x45 + V{2}*x47 + V{2}*x49 + V{1}) + V{1};
auto x89 = ((x50)*(x50));
auto x51 = x3 + x89;
auto x52 = x1*x51;
auto x90 = -x52;
auto x91 = x90 + V{1};
auto x53 = x33*x34 - x81*(x5 + x6*(x52 + V{-1}));
auto x54 = -V{1}*x50;
auto x55 = x54 + x80;
auto x56 = x34*x47 - x81*(x14 - x15*(x0*x55 + x13*((x55)*(x55)) + x91));
auto x57 = x34*x49 - x81*(x18 - x19*(x0*x54 + x13*((x54)*(x54)) + x91));
auto x58 = x21 + x54;
auto x59 = x34*x45 - x81*(x14 - x15*(x0*x58 + x13*((x58)*(x58)) + x91));
auto x60 = x34*x40 - x81*(x18 - x19*(x27 + x90));
auto x61 = x2 + x21;
auto x62 = x21 + x50;
auto x63 = x2 + x80;
auto x64 = x50 + x80;
auto x65 = x34*x39 - x81*(x18 - x19*(x29 + x90));
auto x66 = -V{1}*x76*(x53 + V{2}*x56 + V{2}*x57 + V{2}*x59 + x60 + x65 + V{1}) + V{1};
cell[0] = x53;
cell[1] = x56;
cell[2] = x57;
cell[3] = x59;
cell[4] = x60;
cell[5] = x34*(x15*(x0*x61 + x13*((x61)*(x61)) + x86) + x46 + x84) - x81*(x14 - x15*(x0*x62 + x13*((x62)*(x62)) + x91));
cell[6] = x34*(x19*(x0*x2 + x13*x82 + x86) + x48 + x87) - x81*(x18 - x19*(x0*x50 + x13*x89 + x91));
cell[7] = x34*(x15*(x0*x63 + x13*((x63)*(x63)) + x86) + x44 + x84) - x81*(x14 - x15*(x0*x64 + x13*((x64)*(x64)) + x91));
cell[8] = x65;
cell.template getFieldPointer<olb::descriptors::VELOCITY>()[0] = x66;
cell.template getFieldPointer<olb::descriptors::VELOCITY>()[1] = x80;
cell.template getFieldPointer<olb::fields::moments::RHO>()[0] = x78;
cell.template getFieldPointer<olb::fields::moments::U>()[0] = x66;
cell.template getFieldPointer<olb::fields::moments::U>()[1] = x80;
return { x78, x51 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
