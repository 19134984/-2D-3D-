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

namespace operators {

template <typename T, typename...FIELDS>
struct CSE_O<operators::TrackAverageFunctorO<functors::DerivativeF<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, descriptors::POPULATION, dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination> >, opti::AVERAGE_DJDF>,descriptors::D2Q9<FIELDS...>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS>
void apply(CELL& cell, PARAMETERS& parameters) any_platform {
using V = typename CELL::value_t;
auto x92 = parameters.template get<parameters::CONVERSION_VELOCITY>();
auto x96 = parameters.template get<descriptors::LATTICE_TIME_T>();
auto x91 = parameters.template get<parameters::LATTICE_VISCOSITY>();
auto x93 = parameters.template get<parameters::PHYS_CHAR_VISCOSITY>();
auto x95 = parameters.template get<parameters::FACTOR>();
auto x94 = parameters.template get<parameters::DT>();
auto x90 = parameters.template get<parameters::DX>();
auto x89 = parameters.template get<parameters::OMEGA>();
auto x97 = cell[1] + cell[8];
auto x98 = cell[4] + x97;
auto x99 = cell[5] + cell[7];
auto x100 = cell[6] + x99;
auto x101 = cell[2] + x100;
auto x109 = V{1} / (x96);
auto x111 = ((x92)*(x92));
auto x0 = V{0.5};
auto x110 = x0*x93;
auto x1 = -V{1}*cell[3];
auto x2 = -V{1}*cell[5];
auto x102 = cell[7] + x2;
auto x3 = cell[0] + cell[3] + x101 + x98 + V{1};
auto x4 = V{1}/x3;
auto x106 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*x4;
auto x5 = x4*(-V{1}*cell[4] + x1 + x102 + x97);
auto x6 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*x5;
auto x7 = x3*x6;
auto x8 = V{3};
auto x9 = V{1}/x8;
auto x10 = x9*(x3 + V{-1});
auto x103 = cell[3] - x10;
auto x11 = x103 - x6*x7 + x98 + x99;
auto x12 = -V{1}*cell[1];
auto x13 = x4*(-V{1}*cell[2] + x1 + x100 + x12);
auto x14 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0)*x13;
auto x15 = x14*x3;
auto x16 = cell[3] + x102 + x12 - x15*x6;
auto x104 = V{2}*x16;
auto x17 = x104;
auto x18 = cell[1] + x101 + x103 - x14*x15;
auto x19 = x4*x8*x89;
auto x20 = V{1}/x94;
auto x105 = x19*x20;
auto x21 = x105;
auto x22 = -V{2}*x105*util::pow(x21, V{1})*x4*(x16*x17 + ((x11)*(x11)) + ((x18)*(x18)));
auto x23 = V{1}*x5;
auto x24 = -x106*x23;
auto x25 = V{1}*x6;
auto x26 = V{1}*x9;
auto x107 = x26 + V{-1};
auto x27 = -x11*(x24*x7 + x26 + x6*(x24*x3 + x25));
auto x28 = V{2}*x27;
auto x29 = x15*x24;
auto x30 = V{1}*x13;
auto x31 = -x106*x30;
auto x32 = V{1}*x14;
auto x33 = x3*x31 + x32;
auto x34 = x33*x6;
auto x35 = -x29 - x34;
auto x36 = -x18*(x14*x33 + x15*x31 + x26);
auto x37 = V{2}*x36;
auto x38 = util::pow(x21, V{2});
auto x39 = x14*x31;
auto x40 = V{2}*x39;
auto x41 = x24*x6;
auto x42 = V{2}*x41;
auto x43 = -V{1}*x89*x93*(cell.template getFieldComponent<olb::descriptors::POROSITY>(0) + V{-1})/(((x90)*(x90))*x91);
auto x44 = x96 + V{-1};
auto x45 = -x106*(x23 + V{-1});
auto x46 = -x11*(x107 + x45*x7 + x6*(x25 + x3*x45));
auto x47 = V{2}*x46;
auto x48 = V{-1};
auto x108 = -x48;
auto x49 = x15*x45;
auto x50 = -x106*(x108 + x30);
auto x51 = x3*x50 + x32;
auto x52 = x51*x6;
auto x53 = -x108 - x49 - x52;
auto x54 = -x18*(x107 + x14*x51 + x15*x50);
auto x55 = V{2}*x54;
auto x56 = x14*x50;
auto x57 = V{2}*x56;
auto x58 = x45*x6;
auto x59 = V{2}*x58;
auto x60 = -x29 - x52;
auto x61 = -x106*(x108 + x23);
auto x62 = -x11*(x107 + x6*(x25 + x3*x61) + x61*x7);
auto x63 = V{2}*x62;
auto x64 = x15*x61;
auto x65 = -x52 - x64 + V{1};
auto x66 = x6*x61;
auto x67 = V{2}*x66;
auto x68 = -x34 - x64;
auto x69 = -x106*(x30 + V{-1});
auto x70 = x3*x69 + x32;
auto x71 = x6*x70;
auto x72 = -x108 - x64 - x71;
auto x73 = -x18*(x107 + x14*x70 + x15*x69);
auto x74 = V{2}*x73;
auto x75 = x14*x69;
auto x76 = V{2}*x75;
auto x77 = -x29 - x71;
auto x78 = -x49 - x71 + V{1};
auto x79 = -x34 - x49;
auto x80 = cell[0];
auto x81 = cell[1];
auto x82 = cell[2];
auto x83 = cell[3];
auto x84 = cell[4];
auto x85 = cell[5];
auto x86 = cell[6];
auto x87 = cell[7];
auto x88 = cell[8];
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[0] = x109*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(0)*x44 + x95*(x110*(x22 + x38*(x104*x35 + x17*x35 + x28 + x37)) + x111*x43*(x40 + x42)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[1] = x109*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(1)*x44 + x95*(x110*(x22 + x38*(x104*x53 + x17*x53 + x47 + x55)) + x111*x43*(x57 + x59)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[2] = x109*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(2)*x44 + x95*(x110*(x22 + x38*(x104*x60 + x17*x60 + x28 + x55)) + x111*x43*(x42 + x57)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[3] = x109*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(3)*x44 + x95*(x110*(x22 + x38*(x104*x65 + x17*x65 + x55 + x63)) + x111*x43*(x57 + x67)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[4] = x109*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(4)*x44 + x95*(x110*(x22 + x38*(x104*x68 + x17*x68 + x37 + x63)) + x111*x43*(x40 + x67)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[5] = x109*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(5)*x44 + x95*(x110*(x22 + x38*(x104*x72 + x17*x72 + x63 + x74)) + x111*x43*(x67 + x76)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[6] = x109*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(6)*x44 + x95*(x110*(x22 + x38*(x104*x77 + x17*x77 + x28 + x74)) + x111*x43*(x42 + x76)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[7] = x109*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(7)*x44 + x95*(x110*(x22 + x38*(x104*x78 + x17*x78 + x47 + x74)) + x111*x43*(x59 + x76)));
cell.template getFieldPointer<olb::opti::AVERAGE_DJDF>()[8] = x109*(cell.template getFieldComponent<olb::opti::AVERAGE_DJDF>(8)*x44 + x95*(x110*(x22 + x38*(x104*x79 + x17*x79 + x37 + x47)) + x111*x43*(x40 + x59)));
cell[0] = x80;
cell[1] = x81;
cell[2] = x82;
cell[3] = x83;
cell[4] = x84;
cell[5] = x85;
cell[6] = x86;
cell[7] = x87;
cell[8] = x88;
}
};

}

}

// Generation Info: commit=21521319fcfc23b85b6fd87340ebc70ad50346d4
