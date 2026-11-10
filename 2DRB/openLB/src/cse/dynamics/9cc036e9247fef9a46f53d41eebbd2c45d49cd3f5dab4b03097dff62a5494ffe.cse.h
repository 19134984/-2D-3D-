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
struct CSE<CombinedRLBdynamics<T, descriptors::D2Q9<FIELDS...>, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::SaveVelocity<collision::BGK>, dynamics::DefaultCombination>, momenta::Tuple<momenta::FixedDensity, momenta::FixedPressureMomentum<0, 1>, momenta::RegularizedBoundaryStress<0, 1>, momenta::DefineSeparately> >, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x76 = cell.template getFieldComponent<olb::momenta::FixedDensity::RHO>(0);
auto x78 = cell.template getFieldComponent<olb::momenta::FixedPressureMomentum<0, 1>::VELOCITY>(1);
auto x79 = parameters.template get<descriptors::OMEGA>();
auto x74 = V{1} / (x76);
auto x75 = V{2}*cell[6];
auto x77 = cell[4] + cell[8];
auto x80 = ((x78)*(x78));
auto x0 = V{3};
auto x82 = ((x0)*(x0));
auto x1 = V{0.5}*x0;
auto x2 = x74*(cell[0] + V{2}*cell[5] + V{2}*cell[7] + x75 + x77 + V{1}) + V{-1};
auto x81 = ((x2)*(x2));
auto x3 = x80;
auto x4 = x1*(x3 + x81);
auto x84 = -x4;
auto x85 = x84 + V{1};
auto x5 = V{0.444444444444444};
auto x6 = x5*x76;
auto x7 = V{1}/x0;
auto x8 = -x7;
auto x9 = -V{1}*x78;
auto x10 = x2 + x9;
auto x11 = V{0.5}*x82;
auto x12 = V{0.0277777777777778};
auto x83 = -x12;
auto x13 = x12*x76;
auto x14 = x13*(x0*x10 + x11*((x10)*(x10)) + x85) + x83;
auto x15 = cell[5] - x14;
auto x16 = V{0.111111111111111};
auto x86 = -x16;
auto x17 = x16*x76;
auto x18 = x17*(x0*x2 + x11*x81 + x85) + x86;
auto x19 = x2 + x78;
auto x20 = x13*(x0*x19 + x11*((x19)*(x19)) + x85) + x83;
auto x21 = V{2}*cell[7] - V{2}*x20;
auto x22 = V{2}*x15 - V{2}*x18 + x21 + x75;
auto x23 = x22*x8;
auto x24 = x0*x9 + x11*((x9)*(x9)) + V{1};
auto x25 = -x16 + x17*(x24 + x84);
auto x26 = x0*x78 + x11*x80 + V{1};
auto x27 = -x16 + x17*(x26 + x84);
auto x28 = V{-2};
auto x87 = x15*x28;
auto x29 = x21 - x25 - x27 + x77 - V{1}*x87;
auto x30 = x29*x8;
auto x31 = V{0.5}*x5*x82*(x23 + x30) - x5 - x6*(x4 + V{-1});
auto x32 = V{1} - x79;
auto x33 = V{0.5}*x16*x82;
auto x34 = V{1} - x7;
auto x35 = x29*x34;
auto x36 = x33*(x23 + x35);
auto x37 = x27 + x36;
auto x38 = x25 + x36;
auto x39 = V{0.5}*x12*x82;
auto x40 = x22*x34;
auto x88 = x35 + x40;
auto x41 = x21 + x87;
auto x42 = x39*(V{2}*x41 + x88);
auto x43 = x20 + x42;
auto x44 = x39*(x28*x41 + x88);
auto x45 = x14 + x44;
auto x46 = x33*(x30 + x40);
auto x47 = x18 + x46;
auto x48 = x74*(x31 + x37 + x38 + V{2}*x43 + V{2}*x45 + V{2}*x47 + V{1}) + V{-1};
auto x89 = ((x48)*(x48));
auto x49 = x3 + x89;
auto x50 = x1*x49;
auto x90 = -x50;
auto x91 = x90 + V{1};
auto x51 = x31*x32 - x79*(x5 + x6*(x50 + V{-1}));
auto x52 = -V{1}*x2;
auto x53 = x52 + x78;
auto x54 = -V{1}*x48;
auto x55 = x54 + x78;
auto x56 = x52 + x9;
auto x57 = x54 + x9;
auto x58 = x32*x38 - x79*(x16 - x17*(x24 + x90));
auto x59 = x48 + x9;
auto x60 = x32*x45 - x79*(x12 - x13*(x0*x59 + x11*((x59)*(x59)) + x91));
auto x61 = x32*x47 - x79*(x16 - x17*(x0*x48 + x11*x89 + x91));
auto x62 = x48 + x78;
auto x63 = x32*x43 - x79*(x12 - x13*(x0*x62 + x11*((x62)*(x62)) + x91));
auto x64 = x32*x37 - x79*(x16 - x17*(x26 + x90));
cell[0] = x51;
cell[1] = x32*(x13*(x0*x53 + x11*((x53)*(x53)) + x85) + x44 + x83) - x79*(x12 - x13*(x0*x55 + x11*((x55)*(x55)) + x91));
cell[2] = x32*(x17*(x0*x52 + x11*((x52)*(x52)) + x85) + x46 + x86) - x79*(x16 - x17*(x0*x54 + x11*((x54)*(x54)) + x91));
cell[3] = x32*(x13*(x0*x56 + x11*((x56)*(x56)) + x85) + x42 + x83) - x79*(x12 - x13*(x0*x57 + x11*((x57)*(x57)) + x91));
cell[4] = x58;
cell[5] = x60;
cell[6] = x61;
cell[7] = x63;
cell[8] = x64;
cell.template getFieldPointer<olb::descriptors::VELOCITY>()[0] = x74*(x51 + x58 + V{2}*x60 + V{2}*x61 + V{2}*x63 + x64 + V{1}) + V{-1};
cell.template getFieldPointer<olb::descriptors::VELOCITY>()[1] = x78;
return { x76, x49 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
