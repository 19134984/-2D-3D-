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
struct CSE<CombinedRLBdynamics<T, descriptors::D3Q19<FIELDS...>, dynamics::Tuple<T, descriptors::D3Q19<FIELDS...>, momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq>, equilibria::SecondOrder, collision::RLB, dynamics::DefaultCombination>, momenta::Tuple<momenta::VelocityBoundaryDensity<0, 1>, momenta::FixedVelocityMomentumGeneric, momenta::RegularizedBoundaryStress<0, 1>, momenta::DefineSeparately> >> {
template <concepts::Cell CELL, concepts::Parameters PARAMETERS, concepts::BaseType V=typename CELL::value_t>
CellStatistic<V> collide(CELL& cell, PARAMETERS& parameters) any_platform {
auto x21 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(2);
auto x22 = parameters.template get<descriptors::OMEGA>();
auto x19 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(0);
auto x20 = cell.template getFieldComponent<olb::momenta::FixedVelocityMomentumGeneric::VELOCITY>(1);
auto x23 = x22 + V{-1};
auto x24 = V{1} / (x19 + V{1});
auto x25 = V{0.0277777777777778}*x24;
auto x26 = V{3}*x21;
auto x27 = V{4.5}*((x20 + x21)*(x20 + x21));
auto x28 = V{3}*x20;
auto x29 = ((x21)*(x21));
auto x30 = V{1.5}*x29;
auto x31 = ((x19)*(x19));
auto x32 = V{1.5}*x31;
auto x33 = ((x20)*(x20));
auto x34 = V{1.5}*x33;
auto x35 = x32 + x34 + V{-1};
auto x36 = x30 + x35;
auto x37 = x28 + x36;
auto x38 = x26 - x27 + x37;
auto x39 = V{2}*cell[14] + V{2}*cell[16];
auto x40 = cell[0] + V{2}*cell[10] + cell[11] + cell[12] + V{2}*cell[13] + V{2}*cell[15] + cell[17] + cell[18] + cell[2] + cell[3] + cell[8] + cell[9] + x39 + V{1};
auto x41 = x24*x40;
auto x42 = V{0.0277777777777778}*x41;
auto x43 = -x21;
auto x44 = x20 + x43;
auto x45 = -x26;
auto x46 = x37 + x45;
auto x47 = x46 - V{4.5}*((x44)*(x44));
auto x48 = V{0.0555555555555555}*x41;
auto x49 = V{3}*x29;
auto x50 = x26 + x35 - x49;
auto x51 = V{3}*x19;
auto x52 = -x51;
auto x53 = x19 + x43;
auto x54 = x26 + x36;
auto x55 = x52 + x54 - V{4.5}*((x53)*(x53));
auto x56 = V{0.0555555555555555}*x41;
auto x57 = V{3}*x33;
auto x58 = x30 + V{-1};
auto x59 = x28 + x32 - x57 + x58;
auto x60 = V{0.333333333333333}*x24;
auto x61 = x40*x60;
auto x62 = V{0.0277777777777778}*x41;
auto x63 = V{0.0555555555555555}*x41;
auto x64 = x19 - x20;
auto x65 = x37 + x52 - V{4.5}*((x64)*(x64));
auto x66 = V{0.0277777777777778}*x41;
auto x67 = -x28;
auto x68 = -V{4.5}*((x44)*(x44));
auto x69 = x54 + x67 + x68;
auto x70 = V{4.5}*((x19 + x21)*(x19 + x21));
auto x71 = -x30;
auto x72 = V{1} - x32;
auto x73 = x71 + x72;
auto x74 = -x34;
auto x75 = x26 + x74;
auto x76 = x51 + x70 + x73 + x75;
auto x77 = V{0.0555555555555556}*x24;
auto x78 = x40*x77;
auto x79 = x28 + x73;
auto x80 = x57 + x79;
auto x81 = x49 + x72 + x75;
auto x82 = V{3}*x31;
auto x83 = x51 + x74;
auto x84 = x71 + x82 + x83 + V{1};
auto x85 = V{4.5}*((x19 + x20)*(x19 + x20));
auto x86 = x79 + x83 + x85;
auto x87 = x27 + x75 + x79;
auto x88 = V{4.44089209850063e-16}*cell[10] + V{1.66533453693773e-16}*cell[11] + V{1.11022302462516e-16}*cell[12] + V{6.66133814775094e-16}*cell[13] + V{6.66133814775094e-16}*cell[14] + V{8.88178419700125e-16}*cell[15] + V{4.44089209850063e-16}*cell[16] + V{4.71844785465692e-16}*cell[17] + V{5.27355936696949e-16}*cell[18] + V{1.66533453693773e-16}*cell[2] + V{1.11022302462516e-16}*cell[3] + V{4.71844785465692e-16}*cell[8] + V{5.27355936696949e-16}*cell[9] + V{0.111111111111111}*x41*x84 + x48*x76 + x63*x86 + x66*x87 + x78*x80 + x78*x81 + V{2.22044604925031e-16};
auto x89 = -x36*x61 - x38*x62 - x48*x50 - x48*x55 - x56*x59 - x63*x65 - x66*x69 + x88;
auto x90 = -x42*x47 + x89;
auto x91 = x38*x90;
auto x92 = x25*x91;
auto x93 = x47*x90;
auto x94 = x25*x93;
auto x95 = x69*x90;
auto x96 = x25*x95;
auto x97 = V{3.46944695195361e-18}*x41;
auto x98 = V{5.20417042793042e-18}*x41;
auto x99 = x50*x90;
auto x100 = x59*x90;
auto x101 = x65*x90;
auto x102 = x55*x90;
auto x103 = V{1.21430643318376e-17}*x41;
auto x104 = V{2.89120579329468e-19}*x41;
auto x105 = -x47;
auto x106 = -x38;
auto x107 = V{1.44560289664734e-19}*x41;
auto x108 = -x55;
auto x109 = V{0.0185185185185185}*x24;
auto x110 = -x50;
auto x111 = -x59;
auto x112 = -x65;
auto x113 = -x69;
auto x114 = x105*x42 + x106*x62 + x108*x48 + x110*x48 + x111*x56 + x112*x63 + x113*x66 - x36*x61 + x88;
auto x115 = x114*x84;
auto x116 = V{0.00462962962962963}*x24;
auto x117 = x114*x87;
auto x118 = x113*x114;
auto x119 = x105*x114;
auto x120 = x106*x114;
auto x121 = V{2.89120579329468e-19}*x41;
auto x122 = V{0.0833333333333333}*cell[11];
auto x123 = V{0.0833333333333333}*cell[2];
auto x124 = V{9.62964972193618e-35}*x41;
auto x125 = x111*x124;
auto x126 = x110*x124;
auto x127 = V{1.92592994438724e-34}*x41*x84;
auto x128 = V{0.00462962962962963}*x24;
auto x129 = x114*x76;
auto x130 = x108*x114;
auto x131 = x114*x80;
auto x132 = x111*x114;
auto x133 = -x116*x131 - x116*x132 + x122 + x123 + x125 + x126 + x127 + x128*x129 + x128*x130 + V{-6.93889390390723e-18};
auto x134 = V{0.0833333333333333}*cell[12];
auto x135 = V{0.0833333333333333}*cell[3];
auto x136 = V{0.166666666666667}*cell[14];
auto x137 = -x136;
auto x138 = V{1.44560289664734e-19}*x41;
auto x139 = x138*x86;
auto x140 = x114*x86;
auto x141 = x112*x114;
auto x142 = x114*x81;
auto x143 = x110*x114;
auto x144 = V{2.89120579329468e-19}*x41;
auto x145 = x144*x81;
auto x146 = -x145;
auto x147 = x112*x138 - x116*x142 - x116*x143 + x128*x140 + x128*x141 + x134 + x135 + x137 + x139 + x146;
auto x148 = -V{0.333333333333333}*cell[10] - V{0.166666666666667}*cell[13] - V{0.166666666666667}*cell[15] - V{0.166666666666667}*cell[16] + V{0.166666666666667}*cell[17] + V{0.166666666666667}*cell[18] + V{0.166666666666667}*cell[8] + V{0.166666666666667}*cell[9] + x104*x87 + x107*x76 - x144*x80;
auto x149 = V{0.0833333333333334}*cell[17];
auto x150 = V{0.0833333333333334}*cell[18];
auto x151 = V{0.0833333333333334}*cell[8];
auto x152 = V{0.0833333333333334}*cell[9];
auto x153 = V{0.166666666666667}*cell[11];
auto x154 = V{0.166666666666667}*cell[2];
auto x155 = V{0.166666666666667}*cell[13];
auto x156 = V{0.00925925925925926}*x24;
auto x157 = V{2.89120579329468e-19}*x41;
auto x158 = x157*x76;
auto x159 = V{0.00925925925925926}*x24;
auto x160 = V{1.44560289664734e-19}*x41;
auto x161 = x160*x87;
auto x162 = x127 + V{-6.93889390390723e-18};
auto x163 = V{1.44560289664734e-19}*x41;
auto x164 = V{0.00231481481481482}*x24;
auto x165 = V{0.166666666666667}*cell[10] + x113*x163 - x115*x159 + x117*x164 + x118*x164 + x119*x164 + x120*x164;
auto x166 = V{0.0833333333333334}*cell[17];
auto x167 = V{0.0833333333333334}*cell[18];
auto x168 = V{0.0833333333333334}*cell[8];
auto x169 = V{0.0833333333333334}*cell[9];
auto x170 = V{0.166666666666667}*cell[12];
auto x171 = V{0.166666666666667}*cell[16];
auto x172 = V{0.166666666666667}*cell[3];
auto x173 = V{0.166666666666667}*cell[15];
auto x174 = V{1.44560289664734e-19}*x41;
auto x175 = x157*x86;
auto x176 = V{1.44560289664734e-19}*x41;
auto x177 = x176*x87;
auto x178 = x121*x80;
auto x179 = V{0.833333333333333}*cell[13];
auto x180 = V{0.0231481481481481}*x24;
auto x181 = V{7.2280144832367e-19}*x41;
auto x182 = V{1.59016318631207e-18}*x41;
auto x183 = V{0.00462962962962963}*x24;
auto x184 = V{2.89120579329468e-19}*x41;
auto x185 = V{3.61400724161835e-19}*x41;
auto x186 = V{7.2280144832367e-19}*x41;
auto x187 = V{0.00231481481481481}*x24;
auto x188 = V{0.0833333333333334}*cell[15];
auto x189 = V{0.0833333333333334}*cell[16];
auto x190 = V{0.0833333333333334}*cell[11];
auto x191 = V{0.0833333333333334}*cell[2];
auto x192 = V{9.39641882820771e-19}*x41;
auto x193 = V{4.33680868994202e-19}*x41;
auto x194 = x111*x193;
auto x195 = x110*x193;
auto x196 = V{8.67361737988404e-19}*x41*x84;
auto x197 = V{0.00462962962962963}*x24;
auto x198 = x108*x192 + x129*x164 + x130*x164 + x131*x197 + x132*x197 - x188 - x189 - x190 - x191 + x192*x76 + x194 + x195 + x196 + V{-0.0416666666666667};
auto x199 = V{0.0416666666666667}*cell[17];
auto x200 = V{0.0416666666666667}*cell[18];
auto x201 = V{0.0416666666666667}*cell[8];
auto x202 = V{0.0416666666666667}*cell[9];
auto x203 = V{0.166666666666667}*cell[10];
auto x204 = V{0.00115740740740741}*x24;
auto x205 = V{5.05961013826569e-19}*x41;
auto x206 = x113*x205 + x115*x156 + x117*x204 + x118*x204 + x119*x204 + x120*x204 - x199 - x200 - x201 - x202 - x203;
auto x207 = V{0.0416666666666667}*cell[12] + V{0.0416666666666667}*cell[3] + x105*x185 + x106*x185 - x142*x187 - x143*x187 + x184*x81 + x185*x87 + x186*x80 + x198 + x206;
auto x208 = x36 + x51;
auto x209 = -x23*(V{0.166666666666666}*cell[13] - V{0.833333333333333}*cell[14] + x112*x181 - x140*x183 + x141*x180 + x182*x86 + x207) + V{0.0277777777777778};
auto x210 = V{0.833333333333333}*cell[15];
auto x211 = V{0.0833333333333334}*cell[13];
auto x212 = V{3.61400724161835e-19}*x41;
auto x213 = V{0.0833333333333334}*cell[14];
auto x214 = V{0.0833333333333334}*cell[12];
auto x215 = V{0.0833333333333334}*cell[3];
auto x216 = V{7.2280144832367e-19}*x41;
auto x217 = x112*x192 + x140*x164 + x141*x164 + x142*x197 + x143*x197 + x192*x86 - x213 - x214 - x215 + x216*x81;
auto x218 = V{0.0416666666666667}*cell[11] + V{0.0416666666666667}*cell[2] + x105*x212 + x106*x212 - x131*x187 - x132*x187 + x184*x80 + x194 + x195 + x196 + x206 - x211 + x212*x87 + x217 + V{-0.0416666666666667};
auto x219 = -x23*(V{0.166666666666666}*cell[15] - V{0.833333333333333}*cell[16] + x108*x181 - x129*x183 + x130*x180 + x182*x76 + x218) + V{0.0277777777777778};
auto x220 = V{0.416666666666667}*cell[17];
auto x221 = V{0.416666666666667}*cell[8];
auto x222 = V{3.61400724161835e-20}*x41;
auto x223 = V{0.0115740740740741}*x24;
auto x224 = V{0.00231481481481481}*x24;
auto x225 = V{0.0833333333333335}*cell[13];
auto x226 = V{1.44560289664734e-19}*x41;
auto x227 = V{0.0833333333333333}*cell[10] + x105*x226 + x113*x186 - x115*x116 + x198 + x216*x80 + x217 - x225;
auto x228 = V{2.52980506913284e-19}*x41;
auto x229 = -x23*(V{0.0833333333333333}*cell[17] - V{0.416666666666667}*cell[18] + V{0.0833333333333333}*cell[8] - V{0.416666666666667}*cell[9] + x106*x228 - x117*x224 + x118*x223 + x119*x223 - x120*x224 + x227 + x228*x87) + V{0.0277777777777778};
auto x230 = x86*x90;
auto x231 = x76*x90;
auto x232 = x84*x90;
auto x233 = x124*x50;
auto x234 = x124*x59;
auto x235 = x80*x90;
auto x236 = x116*x235;
auto x237 = x81*x90;
auto x238 = x116*x237;
auto x239 = x138*x65;
auto x240 = x101*x128;
auto x241 = x102*x128;
auto x242 = -x42*(x46 + x68) + x89;
auto x243 = -V{0.166666666666667}*cell[10] + x159*x232 + x163*x69 + x164*x91 + x164*x93 + x164*x95 + x233 + x234 - V{1.92592994438724e-34}*x24*x40*x84 - V{0.00231481481481482}*x24*x87*x90 + V{6.93889390390723e-18};
auto x244 = x193*x50;
auto x245 = x193*x59;
auto x246 = -V{8.67361737988404e-19}*x24*x40*x84;
auto x247 = x100*x197 + x102*x164 + x188 + x189 + x190 + x191 + x192*x55 - V{9.39641882820771e-19}*x24*x40*x76 - V{0.00231481481481482}*x24*x76*x90 - V{0.00462962962962963}*x24*x80*x90 + x244 + x245 + x246 + V{0.0416666666666667};
auto x248 = x199 + x200 + x201 + x202 + x203 + x204*x91 + x204*x93 + x204*x95 + x205*x69 - V{0.00925925925925926}*x24*x84*x90 - V{0.00115740740740741}*x24*x87*x90;
auto x249 = x101*x164 + x192*x65 + x197*x99 + x213 + x214 + x215 - V{7.2280144832367e-19}*x24*x40*x81 - V{9.39641882820771e-19}*x24*x40*x86 - V{0.00462962962962963}*x24*x81*x90 - V{0.00231481481481482}*x24*x86*x90;
auto x0 = -x23*(-V{1}*cell[10] - V{0.5}*cell[11] - V{0.5}*cell[12] - V{2}*cell[13] - V{2}*cell[15] - V{1}*cell[17] - V{1}*cell[18] - V{0.5}*cell[2] - V{0.5}*cell[3] - V{1}*cell[8] - V{1}*cell[9] - x100*x25 - x101*x77 - x102*x77 - x103*x55 - x103*x65 + V{1.21430643318376e-17}*x24*x40*x76 + V{6.93889390390723e-18}*x24*x40*x80 + V{6.93889390390723e-18}*x24*x40*x81 + V{1.04083408558608e-17}*x24*x40*x84 + V{1.21430643318376e-17}*x24*x40*x86 + V{3.46944695195361e-18}*x24*x40*x87 + V{0.0555555555555556}*x24*x76*x90 + V{0.0277777777777778}*x24*x80*x90 + V{0.0277777777777778}*x24*x81*x90 + V{0.0555555555555556}*x24*x84*x90 + V{0.0555555555555556}*x24*x86*x90 + V{0.0277777777777778}*x24*x87*x90 - x25*x99 - x38*x97 - x39 - V{6.93889390390723e-18}*x41*x69 - x47*x97 - x50*x98 - x59*x98 - x92 - x94 - x96 + V{-0.5}) - x36*x60*x90 + V{-0.333333333333333};
auto x1 = x23*(x104*x105 + x104*x106 + x107*x108 + x109*x115 - x113*x121 - x116*x117 - x116*x118 - x116*x119 - x116*x120 + x133 + x147 + x148) - x77*x90*(x34 + x51 + x58 - x82) + V{-0.0555555555555556};
auto x2 = -x100*x77 + x23*(V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] - x105*x160 - x106*x160 - x108*x157 + x125 + x126 - x129*x159 - x130*x159 + x131*x156 + x132*x156 + x147 - x149 - x150 - x151 - x152 - x153 - x154 - x155 - x158 - x161 + x162 + x165 + V{5.78241158658936e-19}*x41*x80) + V{-0.0555555555555556};
auto x3 = x23*(V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] - x105*x176 - x106*x176 + x108*x174 - x112*x157 + x133 - x140*x159 - x141*x159 + x142*x156 + x143*x156 + x165 - x166 - x167 - x168 - x169 - x170 - x171 - x172 - x173 + x174*x76 - x175 - x177 - x178 + V{5.78241158658936e-19}*x41*x81) - x77*x99 + V{-0.0555555555555556};
auto x4 = x23*(V{0.166666666666667}*cell[14] + x112*x182 + x140*x180 - x141*x183 - x179 + x181*x86 + x207) - x25*x90*(x37 + x51 - x85) + V{-0.0277777777777778};
auto x5 = -(x209 + x25*x90*(x208 + x67 - V{4.5}*((x64)*(x64))));
auto x6 = x23*(V{0.166666666666667}*cell[16] + x108*x182 + x129*x180 - x130*x183 + x181*x76 - x210 + x218) - x25*x90*(x51 + x54 - x70) + V{-0.0277777777777778};
auto x7 = -(x219 + x25*x90*(x208 + x45 - V{4.5}*((x53)*(x53))));
auto x8 = x23*(V{0.0833333333333333}*cell[18] + V{0.0833333333333333}*cell[9] + x106*x222 + x117*x223 - x118*x224 - x119*x224 + x120*x223 - x220 - x221 + x222*x87 + x227) - x92 + V{-0.0277777777777778};
auto x9 = -x229 - x94;
auto x10 = x23*(x100*x116 - x104*x38 - x104*x47 - x107*x55 + x109*x232 - x116*x87*x90 + x116*x91 + x116*x93 + x116*x95 + x116*x99 + x121*x69 + x122 + x123 + x128*x230 + x128*x231 + x134 + x135 + x137 + x139 + x146 + x148 + x162 - x233 - x234 - x236 - x238 - x239 - x240 - x241) + x242*x77*x84 + V{-0.0555555555555556};
auto x11 = x23*(V{0.0833333333333333}*cell[12] + V{0.333333333333333}*cell[15] + V{0.333333333333333}*cell[16] + V{0.0833333333333333}*cell[3] - x100*x156 - x136 - x145 - x149 - x150 - x151 - x152 - x153 - x154 - x155 - x158 - x159*x231 - x161 - x238 - x239 + V{1.44560289664734e-19}*x24*x38*x40 + V{1.44560289664734e-19}*x24*x40*x47 + V{2.89120579329468e-19}*x24*x40*x55 + V{5.78241158658936e-19}*x24*x40*x80 + V{1.44560289664734e-19}*x24*x40*x86 + V{0.00462962962962963}*x24*x50*x90 + V{0.00925925925925926}*x24*x55*x90 + V{0.00925925925925926}*x24*x80*x90 + V{0.00462962962962963}*x24*x86*x90 - x240 - x243) + x242*x77*x80 + V{-0.0555555555555556};
auto x12 = x23*(V{0.0833333333333333}*cell[11] + V{0.333333333333333}*cell[13] + V{0.333333333333333}*cell[14] + V{0.0833333333333333}*cell[2] - x156*x99 - x159*x230 - x166 - x167 - x168 - x169 - x170 - x171 - x172 - x173 - x174*x55 - x175 - x177 - x178 - x236 + V{1.44560289664734e-19}*x24*x38*x40 + V{1.44560289664734e-19}*x24*x40*x47 + V{2.89120579329468e-19}*x24*x40*x65 + V{1.44560289664734e-19}*x24*x40*x76 + V{5.78241158658936e-19}*x24*x40*x81 + V{0.00462962962962963}*x24*x59*x90 + V{0.00925925925925926}*x24*x65*x90 + V{0.00462962962962963}*x24*x76*x90 + V{0.00925925925925926}*x24*x81*x90 - x241 - x243) + x242*x77*x81 + V{-0.0555555555555556};
auto x13 = x23*(V{0.0416666666666667}*cell[12] + V{0.166666666666667}*cell[14] + V{0.0416666666666667}*cell[3] - x179 - x182*x65 - x185*x38 - x185*x47 - x187*x237 + V{7.2280144832367e-19}*x24*x40*x80 + V{2.89120579329468e-19}*x24*x40*x81 + V{7.2280144832367e-19}*x24*x40*x86 + V{3.61400724161835e-19}*x24*x40*x87 + V{0.00231481481481481}*x24*x50*x90 + V{0.00462962962962963}*x24*x65*x90 + V{0.0231481481481481}*x24*x86*x90 - x247 - x248) + x242*x25*x86 + V{-0.0277777777777778};
auto x14 = -x101*x25 - x209;
auto x15 = x23*(V{0.0416666666666667}*cell[11] + V{0.166666666666667}*cell[16] + V{0.0416666666666667}*cell[2] - x182*x55 - x187*x235 - x210 - x211 - x212*x38 - x212*x47 + V{7.2280144832367e-19}*x24*x40*x76 + V{2.89120579329468e-19}*x24*x40*x80 + V{3.61400724161835e-19}*x24*x40*x87 + V{0.00462962962962963}*x24*x55*x90 + V{0.00231481481481481}*x24*x59*x90 + V{0.0231481481481481}*x24*x76*x90 - x244 - x245 - x246 - x248 - x249 + V{-0.0416666666666667}) + x242*x25*x76 + V{-0.0277777777777778};
auto x16 = -x102*x25 - x219;
auto x17 = x23*(V{0.0833333333333333}*cell[10] + V{0.0833333333333333}*cell[18] + V{0.0833333333333333}*cell[9] - x116*x232 - x186*x69 - x220 - x221 - x222*x38 - x223*x91 - x225 - x226*x47 + V{7.2280144832367e-19}*x24*x40*x80 + V{3.61400724161835e-20}*x24*x40*x87 + V{0.00231481481481481}*x24*x47*x90 + V{0.00231481481481481}*x24*x69*x90 + V{0.0115740740740741}*x24*x87*x90 - x247 - x249) + x242*x25*x87 + V{-0.0277777777777778};
auto x18 = -x229 - x96;
cell[0] = x0;
cell[1] = x1;
cell[2] = x2;
cell[3] = x3;
cell[4] = x4;
cell[5] = x5;
cell[6] = x6;
cell[7] = x7;
cell[8] = x8;
cell[9] = x9;
cell[10] = x10;
cell[11] = x11;
cell[12] = x12;
cell[13] = x13;
cell[14] = x14;
cell[15] = x15;
cell[16] = x16;
cell[17] = x17;
cell[18] = x18;
return { V{1}*x24*x242, x29 + x31 + x33 };
}
};

}

}

// Generation Info: commit=4df8ae52d6cd7beb4f789c0d5fe28d16c9126065
