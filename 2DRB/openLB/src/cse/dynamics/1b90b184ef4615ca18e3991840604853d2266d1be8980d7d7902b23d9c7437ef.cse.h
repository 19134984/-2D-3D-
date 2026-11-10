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
struct CSE<dynamics::Tuple<T, descriptors::D2Q9<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::MovingPorousMomentumCombination<momenta::BulkMomentum>, momenta::BulkStress, momenta::DefineToNEq>, equilibria::ThirdOrder, collision::ParameterFromCell<collision::HYBRID, collision::LocalSmagorinskyEffectiveOmega<collision::HRR> >, dynamics::DefaultCombination>, olb::CodegenTarget::Generic> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x86 = cell.template getFieldComponent<olb::collision::HYBRID>(0);
auto x88 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x89 = cell.template getFieldComponent<olb::descriptors::TENSOR>(0);
auto x90 = cell.template getFieldComponent<olb::descriptors::TENSOR>(1);
auto x91 = cell.template getFieldComponent<olb::descriptors::TENSOR>(2);
auto x92 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(0);
auto x93 = cell.template getFieldComponent<olb::descriptors::VELOCITY>(1);
auto x94 = cell.template getFieldComponent<olb::descriptors::WMPOROSITY>(0);
auto x95 = cell.template getFieldComponent<olb::descriptors::WMVELOCITY>(0);
auto x96 = cell.template getFieldComponent<olb::descriptors::WMVELOCITY>(1);
auto x98 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x99 = parameters.template get<descriptors::OMEGA>();
auto x87 = cell[1] + cell[7];
auto x97 = cell[8] + x87;
auto x100 = cell[4] + x97;
auto x101 = cell[5] + cell[6];
auto x102 = cell[2] + x101;
auto x0 = V{3};
auto x108 = ((x0)*(x0));
auto x109 = V{1} / (x0);
auto x110 = V{1}*x109;
auto x1 = V{1} - x94;
auto x2 = V{1} - x88;
auto x3 = -V{1}*cell[1];
auto x103 = cell[7] + x3;
auto x4 = -V{1}*cell[3];
auto x5 = cell[0] + cell[3] + x100 + x102 + V{1};
auto x104 = V{1} / (x5);
auto x105 = x104*x88;
auto x6 = x1*x95 + x94*(x105*(-V{1}*cell[2] + x101 + x103 + x4) + x2*x92);
auto x106 = ((x6)*(x6));
auto x7 = x106;
auto x8 = -V{1}*cell[5];
auto x9 = x1*x96 + x94*(x105*(-V{1}*cell[4] + x4 + x8 + x97) + x2*x93);
auto x107 = ((x9)*(x9));
auto x10 = x107 + x7;
auto x11 = V{0.5}*x0*x10;
auto x125 = V{1} - x11;
auto x12 = V{0.444444444444444};
auto x13 = x108;
auto x14 = V{0.5}*x13;
auto x15 = x110;
auto x113 = V{1} - x15;
auto x16 = -x15;
auto x17 = V{1}/x99;
auto x18 = x15*(x5 + V{-1});
auto x111 = cell[3] - x18;
auto x19 = cell[5] + x100 - x107*x5 + x111;
auto x20 = x5*x6;
auto x21 = cell[3] + x103 - x20*x9 + x8;
auto x22 = x102 + x111 - x20*x6 + x87;
auto x23 = V{2}*V{1} / (x17 + util::sqrt(V{2.828428}*x104*x108*((x98)*(x98))*util::sqrt(((x19)*(x19)) + V{2}*((x21)*(x21)) + ((x22)*(x22))) + ((x17)*(x17))));
auto x112 = V{1} / (x23);
auto x25 = V{2}*x109*x112*x5*(x86 + V{-1});
auto x26 = x19*x86 + x25*x91;
auto x27 = x16*x26;
auto x28 = x22*x86 + x25*x89;
auto x29 = x16*x28;
auto x30 = V{1} - x23;
auto x31 = x0*x13;
auto x32 = V{0.166667}*x31;
auto x118 = x16*x32;
auto x33 = x113;
auto x34 = -x110;
auto x35 = -x34 + V{-1};
auto x36 = x7*x9;
auto x37 = x107*x6;
auto x38 = x36 - x37;
auto x123 = x32*x38;
auto x39 = V{0.5}*x31;
auto x40 = x36 + x37;
auto x124 = x39*x40;
auto x41 = -V{1}*x6;
auto x42 = x41 + x9;
auto x43 = V{0.0277777777777778};
auto x119 = -x43;
auto x122 = x30*x43;
auto x44 = x43*x5;
auto x45 = -V{1}*x33;
auto x120 = -x45;
auto x46 = x21*x86 + x25*x90;
auto x114 = V{2}*x46;
auto x115 = x13*x46;
auto x47 = x114*x6 + x28*x9;
auto x48 = x114*x9 + x26*x6;
auto x49 = x47 - x48;
auto x121 = x32*x49*(x120 + x33);
auto x50 = x47 + x48;
auto x116 = x39*x50;
auto x51 = x26*x33;
auto x52 = x28*x33;
auto x53 = x14*(x51 + x52);
auto x54 = -V{1}*x115 + x116*(x33 + x45) + x53;
auto x55 = -x34;
auto x56 = V{0.111111111111111};
auto x126 = -x56;
auto x127 = x30*x56;
auto x57 = x5*x56;
auto x58 = -V{1}*x16;
auto x117 = x32*x58;
auto x59 = x39*x58;
auto x60 = x50*x59;
auto x61 = x14*(x27 + x52);
auto x62 = -V{1}*x9;
auto x63 = x41 + x62;
auto x64 = x115 + x53;
auto x65 = x117;
auto x66 = x14*(x29 + x51);
auto x67 = x113;
auto x68 = x6 + x62;
auto x69 = -x118;
auto x70 = x16*x39;
auto x71 = x40*x70;
auto x128 = x125 + x71;
auto x72 = x50*x70;
auto x73 = 0;
auto x74 = V{2}*x33*x39;
auto x75 = x6 + x9;
auto x76 = x118;
cell[0] = x12*(x14*x30*(x27 + x29) - x5*(x11 + V{-1}) + V{-1});
cell[1] = x119 + x122*(x121 + x54) + x44*(x0*x42 + x123*(x33 - x35) + x124*(x33 + x35) + x125 + x14*((x42)*(x42)));
cell[2] = x126 + x127*(-x117*x49 + x60 + x61) + x57*(x0*x41 - x123*x55 + x124*x55 + x125 + x14*((x41)*(x41)));
cell[3] = x119 + x122*(V{2}*x116*x45 + x64) + x44*(x0*x63 - x123*(x120 + x35) + x124*(x35 + x45) + x125 + x14*((x63)*(x63)));
cell[4] = x126 + x127*(x49*x65 + x60 + x66) + x57*(x0*x62 + x125 + x14*((x62)*(x62)) + x38*x65 + x40*x59);
cell[5] = -(x122*(x121 - x54) + x43 - x44*(x0*x68 + x123*(x45 - x67) + x124*(x45 + x67) + x125 + x14*((x68)*(x68))));
cell[6] = x126 + x127*(x49*x69 + x61 + x72) + x57*(x0*x6 + x106*x14 + x128 + x38*x69);
cell[7] = x119 + x122*(x49*x73 + x50*x74 + x64) + x44*(x0*x75 + x125 + x14*((x75)*(x75)) + x38*x73 + x40*x74);
cell[8] = x126 + x127*(x49*x76 + x66 + x72) + x57*(x0*x9 + x107*x14 + x128 + x38*x76);
cell.template getFieldPointer<olb::descriptors::OMEGA>()[0] = x23;
return { x5, x10 };
}
};

}

}

// Generation Info: commit=e62fb513c7fed8d9804dfb87eb8d98b331948423
