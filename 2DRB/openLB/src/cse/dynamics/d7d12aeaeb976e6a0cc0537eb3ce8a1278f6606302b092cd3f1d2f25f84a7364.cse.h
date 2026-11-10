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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::TrackAveragePopulation<collision::BGK>, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x10 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(1);
auto x17 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(8);
auto x15 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(6);
auto x20 = parameters.template get<descriptors::LATTICE_TIME_T>();
auto x9 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(0);
auto x13 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(4);
auto x19 = parameters.template get<descriptors::OMEGA>();
auto x11 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(2);
auto x16 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(7);
auto x12 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(3);
auto x18 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x14 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(5);
auto x21 = x19 + V{-1};
auto x22 = cell[0] + cell[1] + cell[2] + cell[3] + cell[4] + cell[5] + cell[6] + cell[7] + cell[8];
auto x23 = x22 + V{1};
auto x24 = x22 + V{1};
auto x25 = V{1} / ((x24)*(x24));
auto x26 = V{1.5}*x25;
auto x27 = ((x18)*(x18));
auto x28 = cell[1] - cell[5];
auto x29 = -cell[4] + cell[8];
auto x30 = -cell[3] + cell[7] + x28 + x29;
auto x31 = ((x30)*(x30));
auto x32 = x27*x31;
auto x33 = x26*x32;
auto x34 = cell[2] - cell[6];
auto x35 = cell[3] - cell[7] + x28 + x34;
auto x36 = -x35;
auto x37 = x26*x27*((x36)*(x36)) + V{-1};
auto x38 = x33 + x37;
auto x39 = cell[0]*x21 + V{0.444444444444444}*x19*(x23*x38 + V{1});
auto x40 = V{0.0277777777777778}*x19;
auto x41 = V{3}*x18/x24;
auto x42 = x30*x41;
auto x43 = -x42;
auto x44 = V{4.5}*x25;
auto x45 = V{2}*cell[1] - V{2}*cell[5] + x29 + x34;
auto x46 = x36*x41;
auto x47 = x38 + x46;
auto x48 = cell[1]*x21 + x40*(x23*(-x27*x44*((x45)*(x45)) + x43 + x47) + V{1});
auto x49 = V{0.111111111111111}*x19;
auto x50 = x35*x41;
auto x51 = V{3}*x25;
auto x52 = ((x35)*(x35));
auto x53 = x27*x52;
auto x54 = V{1} - x33;
auto x55 = x51*x53 + x54;
auto x56 = -cell[2]*x21 + x49*(x23*(x50 + x55) + V{-1});
auto x57 = x27*x44*((-2*cell[3] - cell[4] + 2*cell[7] + cell[8] - x34)*(-2*cell[3] - cell[4] + 2*cell[7] + cell[8] - x34));
auto x58 = -x57;
auto x59 = cell[3]*x21 + x40*(x23*(x42 + x47 + x58) + V{1});
auto x60 = x32*x51;
auto x61 = cell[4]*x21 + x49*(x23*(x37 + x42 - x60) + V{1});
auto x62 = cell[5]*x21 + x40*(x23*(-x27*x44*((x45)*(x45)) + x38 + x42 - x46) + V{1});
auto x63 = -x50;
auto x64 = -cell[6]*x21 + x49*(x23*(x55 + x63) + V{-1});
auto x65 = cell[7]*x21;
auto x66 = x26*x53;
auto x67 = x42 - x66;
auto x68 = -cell[8]*x21 + x49*(x23*(x60 + x67 + V{1}) + V{-1});
auto x69 = V{1} / (x20);
auto x70 = x20 + V{-1};
cell[0] = -x39;
cell[1] = -x48;
cell[2] = x56;
cell[3] = -x59;
cell[4] = -x61;
cell[5] = -x62;
cell[6] = x64;
cell[7] = x40*(x23*(x54 + x57 + x63 + x67) + V{-1}) - x65;
cell[8] = x68;
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[0] = x69*(-x39 + x70*x9);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[1] = x69*(x10*x70 - x48);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[2] = x69*(x11*x70 + x56);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[3] = x69*(x12*x70 - x59);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[4] = x69*(x13*x70 - x61);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[5] = x69*(x14*x70 - x62);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[6] = x69*(x15*x70 + x64);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[7] = x69*(x16*x70 - x40*(x23*(x33 + x43 + x50 + x58 + x66 + V{-1}) + V{1}) - x65);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[8] = x69*(x17*x70 + x68);
return { x24, V{1}*x25*x27*(x31 + x52) };
}
};

}

}

// Generation Info: commit=b0dc2eff70ace778d8c09131cb5af79150b5dff9
