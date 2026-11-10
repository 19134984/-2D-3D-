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
struct CSE<dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::TrackAveragePopulation<collision::SmagorinskyEffectiveOmega<collision::TRT> >, dynamics::DefaultCombination>> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x34 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(15);
auto x25 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(6);
auto x39 = parameters.template get<descriptors::OMEGA>();
auto x29 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(10);
auto x37 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(18);
auto x36 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(17);
auto x40 = parameters.template get<collision::TRT::MAGIC>();
auto x30 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(11);
auto x20 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(1);
auto x22 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(3);
auto x38 = cell.template getFieldComponent<olb::descriptors::POROSITY>(0);
auto x28 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(9);
auto x41 = parameters.template get<collision::LES::SMAGORINSKY>();
auto x27 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(8);
auto x31 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(12);
auto x19 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(0);
auto x35 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(16);
auto x32 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(13);
auto x23 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(4);
auto x33 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(14);
auto x24 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(5);
auto x42 = parameters.template get<descriptors::LATTICE_TIME_T>();
auto x26 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(7);
auto x21 = cell.template getFieldComponent<olb::descriptors::AVERAGE_POPULATION>(2);
auto x43 = cell[13] + cell[15];
auto x44 = cell[10] + x43;
auto x45 = cell[11] + cell[5];
auto x46 = cell[12] + cell[7] + cell[9];
auto x47 = cell[2] + cell[8];
auto x48 = cell[18] + cell[3];
auto x49 = cell[16] + cell[6];
auto x50 = cell[0] + cell[14] + cell[17] + cell[1] + cell[4] + x44 + x45 + x46 + x47 + x48 + x49 + V{1};
auto x51 = V{1} / (x50);
auto x52 = -cell[13];
auto x53 = cell[13] + cell[17];
auto x54 = -cell[2];
auto x55 = -cell[9];
auto x56 = -cell[8];
auto x57 = cell[18] + x56;
auto x58 = x54 + x55 + x57;
auto x59 = -cell[14];
auto x60 = -cell[4];
auto x61 = x45 + x59 + x60;
auto x62 = x53 + x58 + x61;
auto x63 = -cell[1];
auto x64 = -cell[7];
auto x65 = -cell[6];
auto x66 = cell[16] + x65;
auto x67 = x63 + x64 + x66;
auto x68 = cell[14] + x60;
auto x69 = -cell[5] + x68;
auto x70 = x44 + x67 + x69;
auto x71 = x51*x70;
auto x72 = ((x38)*(x38));
auto x73 = -cell[15];
auto x74 = cell[15] + cell[17];
auto x75 = -cell[18];
auto x76 = -cell[3];
auto x77 = x56 + x75 + x76;
auto x78 = -cell[16];
auto x79 = x65 + x78;
auto x80 = x46 + x74 + x77 + x79;
auto x81 = x51*x62;
auto x82 = -cell[17];
auto x83 = cell[9] + x82;
auto x84 = V{1}*x51;
auto x85 = ((x80)*(x80));
auto x86 = x72*x85;
auto x87 = V{0.333333333333333}*cell[0];
auto x88 = V{0.333333333333333}*cell[10] - V{0.666666666666667}*cell[17] - V{0.666666666666667}*cell[18] + V{0.333333333333333}*cell[1] - V{0.666666666666667}*cell[8] - V{0.666666666666667}*cell[9] + x87;
auto x89 = V{0.333333333333333}*cell[11] - V{0.666666666666667}*cell[15] - V{0.666666666666667}*cell[16] + V{0.333333333333333}*cell[2] - V{0.666666666666667}*cell[6] - V{0.666666666666667}*cell[7];
auto x90 = ((x62)*(x62));
auto x91 = x72*x90;
auto x92 = V{0.333333333333333}*cell[12] - V{0.666666666666667}*cell[13] - V{0.666666666666667}*cell[14] + V{0.333333333333333}*cell[3] - V{0.666666666666667}*cell[4] - V{0.666666666666667}*cell[5];
auto x93 = ((x70)*(x70));
auto x94 = x72*x93;
auto x95 = util::sqrt(V{25.455852}*x51*((x41)*(x41))*util::sqrt(V{2}*((x57 + x72*x80*x81 + x83)*(x57 + x72*x80*x81 + x83)) + V{2}*((cell[5] + x52 + x62*x71*x72 + x68)*(cell[5] + x52 + x62*x71*x72 + x68)) + V{2}*((cell[7] + x66 + x71*x72*x80 + x73)*(cell[7] + x66 + x71*x72*x80 + x73)) + ((-V{0.666666666666667}*cell[11] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] - V{0.666666666666667}*cell[2] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + x84*x91 + x88 + x92)*(-V{0.666666666666667}*cell[11] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] - V{0.666666666666667}*cell[2] + V{0.333333333333333}*cell[6] + V{0.333333333333333}*cell[7] + x84*x91 + x88 + x92)) + ((-V{0.666666666666667}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] - V{0.666666666666667}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + x84*x86 + x88 + x89)*(-V{0.666666666666667}*cell[12] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] - V{0.666666666666667}*cell[3] + V{0.333333333333333}*cell[4] + V{0.333333333333333}*cell[5] + x84*x86 + x88 + x89)) + ((-V{0.666666666666667}*cell[10] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] - V{0.666666666666667}*cell[1] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] + x84*x94 + x87 + x89 + x92)*(-V{0.666666666666667}*cell[10] + V{0.333333333333333}*cell[17] + V{0.333333333333333}*cell[18] - V{0.666666666666667}*cell[1] + V{0.333333333333333}*cell[8] + V{0.333333333333333}*cell[9] + x84*x94 + x87 + x89 + x92))) + V{1}/((x39)*(x39))) + V{1} / (x39);
auto x96 = V{1} / (x95);
auto x97 = V{1} / ((x50)*(x50));
auto x98 = V{1.5}*x97;
auto x99 = x94*x98;
auto x100 = x91*x98;
auto x101 = x86*x98;
auto x102 = x100 + x101 + V{-1};
auto x103 = x102 + x99;
auto x104 = cell[0] - x96*(V{2}*cell[0] + x103*(V{0.666666666666667}*cell[0] + V{0.666666666666667}*cell[10] + V{0.666666666666667}*cell[11] + V{0.666666666666667}*cell[12] + V{0.666666666666667}*cell[13] + V{0.666666666666667}*cell[14] + V{0.666666666666667}*cell[15] + V{0.666666666666667}*cell[16] + V{0.666666666666667}*cell[17] + V{0.666666666666667}*cell[18] + V{0.666666666666667}*cell[1] + V{0.666666666666667}*cell[2] + V{0.666666666666667}*cell[3] + V{0.666666666666667}*cell[4] + V{0.666666666666667}*cell[5] + V{0.666666666666667}*cell[6] + V{0.666666666666667}*cell[7] + V{0.666666666666667}*cell[8] + V{0.666666666666667}*cell[9] + V{0.666666666666667}) + V{0.666666666666667});
auto x105 = V{1} / (V{2}*x40/(x95 + V{-1}) + V{0.5});
auto x106 = V{0.5}*cell[1];
auto x107 = V{0.0277777777777778}*cell[0] + V{0.0277777777777778}*cell[10] + V{0.0277777777777778}*cell[11] + V{0.0277777777777778}*cell[12] + V{0.0277777777777778}*cell[13] + V{0.0277777777777778}*cell[14] + V{0.0277777777777778}*cell[15] + V{0.0277777777777778}*cell[16] + V{0.0277777777777778}*cell[17] + V{0.0277777777777778}*cell[18] + V{0.0277777777777778}*cell[1] + V{0.0277777777777778}*cell[2] + V{0.0277777777777778}*cell[3] + V{0.0277777777777778}*cell[4] + V{0.0277777777777778}*cell[5] + V{0.0277777777777778}*cell[6] + V{0.0277777777777778}*cell[7] + V{0.0277777777777778}*cell[8] + V{0.0277777777777778}*cell[9] + V{0.0277777777777778};
auto x108 = V{3}*x38;
auto x109 = x108*x71;
auto x110 = x102 + x109 - V{3}*x72*x93*x97;
auto x111 = V{3}*x97;
auto x112 = -x100;
auto x113 = V{1} - x101;
auto x114 = x112 + x113;
auto x115 = x109 + x114;
auto x116 = x111*x94 + x115;
auto x117 = x107*x116;
auto x118 = V{0.0555555555555556}*cell[0] + V{0.0555555555555556}*cell[10] + V{0.0555555555555556}*cell[11] + V{0.0555555555555556}*cell[12] + V{0.0555555555555556}*cell[13] + V{0.0555555555555556}*cell[14] + V{0.0555555555555556}*cell[15] + V{0.0555555555555556}*cell[16] + V{0.0555555555555556}*cell[17] + V{0.0555555555555556}*cell[18] + V{0.0555555555555556}*cell[1] + V{0.0555555555555556}*cell[2] + V{0.0555555555555556}*cell[3] + V{0.0555555555555556}*cell[4] + V{0.0555555555555556}*cell[5] + V{0.0555555555555556}*cell[6] + V{0.0555555555555556}*cell[7] + V{0.0555555555555556}*cell[8] + V{0.0555555555555556}*cell[9] + V{0.0555555555555556};
auto x119 = x96*(V{1}*cell[10] + V{1}*cell[1] + x110*x118 - x116*x118 + V{0.111111111111111});
auto x120 = -x119;
auto x121 = cell[1] + x120;
auto x122 = V{0.5}*cell[2];
auto x123 = x108*x81;
auto x124 = x99 + V{-1};
auto x125 = x101 + x123 + x124 - V{3}*x72*x90*x97;
auto x126 = -x99;
auto x127 = x123 + x126;
auto x128 = x111*x91 + x113 + x127;
auto x129 = x107*x128;
auto x130 = x96*(V{1}*cell[11] + V{1}*cell[2] + x118*x125 - x118*x128 + V{0.111111111111111});
auto x131 = -x130;
auto x132 = cell[2] + x131;
auto x133 = V{0.5}*cell[3];
auto x134 = x108*x51*x80;
auto x135 = x100 + x124 + x134 - V{3}*x72*x85*x97;
auto x136 = x126 + x134;
auto x137 = x111*x86 + x112 + x136 + V{1};
auto x138 = x107*x137;
auto x139 = x96*(V{1}*cell[12] + V{1}*cell[3] + x118*x135 - x118*x137 + V{0.111111111111111});
auto x140 = -x139;
auto x141 = cell[3] + x140;
auto x142 = V{0.5}*cell[4];
auto x143 = V{0.0138888888888889}*cell[0] + V{0.0138888888888889}*cell[10] + V{0.0138888888888889}*cell[11] + V{0.0138888888888889}*cell[12] + V{0.0138888888888889}*cell[13] + V{0.0138888888888889}*cell[14] + V{0.0138888888888889}*cell[15] + V{0.0138888888888889}*cell[16] + V{0.0138888888888889}*cell[17] + V{0.0138888888888889}*cell[18] + V{0.0138888888888889}*cell[1] + V{0.0138888888888889}*cell[2] + V{0.0138888888888889}*cell[3] + V{0.0138888888888889}*cell[4] + V{0.0138888888888889}*cell[5] + V{0.0138888888888889}*cell[6] + V{0.0138888888888889}*cell[7] + V{0.0138888888888889}*cell[8] + V{0.0138888888888889}*cell[9] + V{0.0138888888888889};
auto x144 = cell[10] + x67;
auto x145 = ((cell[11] + 2*cell[13] - 2*cell[4] + x144 + x58 + x74)*(cell[11] + 2*cell[13] - 2*cell[4] + x144 + x58 + x74));
auto x146 = x103 + x109;
auto x147 = x123 - V{4.5}*x145*x72*x97 + x146;
auto x148 = V{4.5}*x97;
auto x149 = x115 + x127 + x145*x148*x72;
auto x150 = x143*x149;
auto x151 = x96*(V{1}*cell[13] + V{1}*cell[4] + x107*x147 - x107*x149 + V{0.0555555555555556});
auto x152 = -x151;
auto x153 = cell[4] + x152;
auto x154 = -cell[11];
auto x155 = V{2}*cell[14] + cell[15] - V{2}*cell[5] + x144 + x154 + x47 + x75 + x83;
auto x156 = -x123;
auto x157 = x146 + x156;
auto x158 = x157 - V{4.5}*x72*x97*((x155)*(x155));
auto x159 = -x148*x72*((x155)*(x155));
auto x160 = x103 - x109;
auto x161 = x123 + x159 + x160;
auto x162 = V{0.5}*cell[14] - V{0.5}*cell[5];
auto x163 = V{1}*cell[14] + V{1}*cell[5] + x107*x161 + V{0.0555555555555556};
auto x164 = x96*(x107*x158 + x163);
auto x165 = V{0.5}*cell[6];
auto x166 = cell[10] + x63 + x69;
auto x167 = ((cell[12] + 2*cell[15] - 2*cell[6] + cell[9] + x166 + x53 + x77)*(cell[12] + 2*cell[15] - 2*cell[6] + cell[9] + x166 + x53 + x77));
auto x168 = x134 + x146 - V{4.5}*x167*x72*x97;
auto x169 = x115 + x136 + x148*x167*x72;
auto x170 = x143*x169;
auto x171 = x96*(V{1}*cell[15] + V{1}*cell[6] + x107*x168 - x107*x169 + V{0.0555555555555556});
auto x172 = -x171;
auto x173 = cell[6] + x172;
auto x174 = -cell[12];
auto x175 = cell[13] + x174;
auto x176 = V{2}*cell[16] - V{2}*cell[7] + cell[8] + x166 + x175 + x48 + x55 + x82;
auto x177 = -x134;
auto x178 = x146 + x177;
auto x179 = x178 - V{4.5}*x72*x97*((x176)*(x176));
auto x180 = -x148*x72*((x176)*(x176));
auto x181 = x134 + x160 + x180;
auto x182 = V{0.5}*cell[16] - V{0.5}*cell[7];
auto x183 = V{1}*cell[16] + V{1}*cell[7] + x107*x181 + V{0.0555555555555556};
auto x184 = x96*(x107*x179 + x183);
auto x185 = V{0.5}*cell[8];
auto x186 = x54 + x61;
auto x187 = ((cell[12] + 2*cell[17] + cell[7] - 2*cell[8] + x186 + x43 + x76 + x79)*(cell[12] + 2*cell[17] + cell[7] - 2*cell[8] + x186 + x43 + x76 + x79));
auto x188 = x103 + x123;
auto x189 = x134 - V{4.5}*x187*x72*x97 + x188;
auto x190 = x114 + x127 + x134 + x148*x187*x72;
auto x191 = x143*x190;
auto x192 = x96*(V{1}*cell[17] + V{1}*cell[8] + x107*x189 - x107*x190 + V{0.0555555555555556});
auto x193 = -x192;
auto x194 = cell[8] + x193;
auto x195 = V{2}*cell[18] + cell[3] - V{2}*cell[9] + x175 + x186 + x49 + x64 + x73;
auto x196 = x177 + x188;
auto x197 = x196 - V{4.5}*x72*x97*((x195)*(x195));
auto x198 = -x148*x72*((x195)*(x195));
auto x199 = x103 + x134 + x156 + x198;
auto x200 = V{0.5}*cell[18] - V{0.5}*cell[9];
auto x201 = V{1}*cell[18] + V{1}*cell[9] + x107*x199 + V{0.0555555555555556};
auto x202 = x96*(x107*x197 + x201);
auto x203 = -V{0.5}*cell[10] + x106 + x107*x110 + x117;
auto x204 = -x105*x203;
auto x205 = -V{0.5}*cell[11] + x107*x125 + x122 + x129;
auto x206 = -x105*x205;
auto x207 = -V{0.5}*cell[12] + x107*x135 + x133 + x138;
auto x208 = -x105*x207;
auto x209 = -V{0.5}*cell[13] + x142 + x143*x147 + x150;
auto x210 = -x105*x209;
auto x211 = x143*x161 + x162;
auto x212 = x105*(-x143*x158 + x211);
auto x213 = -V{0.5}*cell[15] + x143*x168 + x165 + x170;
auto x214 = -x105*x213;
auto x215 = x143*x181 + x182;
auto x216 = x105*(-x143*x179 + x215);
auto x217 = -V{0.5}*cell[17] + x143*x189 + x185 + x191;
auto x218 = -x105*x217;
auto x219 = x143*x199 + x200;
auto x220 = x105*(-x143*x197 + x219);
auto x221 = V{1} / (x42);
auto x222 = x42 + V{-1};
auto x223 = x157 + x159;
auto x224 = -x96*(x107*x223 + x163);
auto x225 = x178 + x180;
auto x226 = -x96*(x107*x225 + x183);
auto x227 = x196 + x198;
auto x228 = -x96*(x107*x227 + x201);
cell[0] = x104;
cell[1] = x105*(V{0.5}*cell[10] - x106 - x107*x110 - x117) + x121;
cell[2] = x105*(V{0.5}*cell[11] - x107*x125 - x122 - x129) + x132;
cell[3] = x105*(V{0.5}*cell[12] - x107*x135 - x133 - x138) + x141;
cell[4] = x105*(V{0.5}*cell[13] - x142 - x143*x147 - x150) + x153;
cell[5] = cell[5] + x105*(-x143*x158 + x143*x161 + x162) - x164;
cell[6] = x105*(V{0.5}*cell[15] - x143*x168 - x165 - x170) + x173;
cell[7] = cell[7] + x105*(-x143*x179 + x143*x181 + x182) - x184;
cell[8] = x105*(V{0.5}*cell[17] - x143*x189 - x185 - x191) + x194;
cell[9] = cell[9] + x105*(-x143*x197 + x143*x199 + x200) - x202;
cell[10] = cell[10] - x119 - x204;
cell[11] = -x130 - x154 - x206;
cell[12] = -x139 - x174 - x208;
cell[13] = -x151 - x210 - x52;
cell[14] = -x164 - x212 - x59;
cell[15] = -x171 - x214 - x73;
cell[16] = -x184 - x216 - x78;
cell[17] = -x192 - x218 - x82;
cell[18] = -x202 - x220 - x75;
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[0] = x221*(x104 + x19*x222);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[1] = x221*(x121 + x20*x222 + x204);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[2] = x221*(x132 + x206 + x21*x222);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[3] = x221*(x141 + x208 + x22*x222);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[4] = x221*(x153 + x210 + x222*x23);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[5] = x221*(cell[5] + x212 + x222*x24 + x224);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[6] = x221*(x173 + x214 + x222*x25);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[7] = x221*(cell[7] + x216 + x222*x26 + x226);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[8] = x221*(x194 + x218 + x222*x27);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[9] = x221*(cell[9] + x220 + x222*x28 + x228);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[10] = x221*(cell[10] + x105*x203 + x120 + x222*x29);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[11] = x221*(cell[11] + x105*x205 + x131 + x222*x30);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[12] = x221*(cell[12] + x105*x207 + x140 + x222*x31);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[13] = x221*(cell[13] + x105*x209 + x152 + x222*x32);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[14] = x221*(cell[14] - x105*(-x143*x223 + x211) + x222*x33 + x224);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[15] = x221*(cell[15] + x105*x213 + x172 + x222*x34);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[16] = x221*(cell[16] - x105*(-x143*x225 + x215) + x222*x35 + x226);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[17] = x221*(cell[17] + x105*x217 + x193 + x222*x36);
cell.template getFieldPointer<olb::descriptors::AVERAGE_POPULATION>()[18] = x221*(cell[18] - x105*(-x143*x227 + x219) + x222*x37 + x228);
return { x50, V{1}*x72*x97*(x85 + x90 + x93) };
}
};

}

}

// Generation Info: commit=cc8b060277f106540151cd2376323351f6e6960e
