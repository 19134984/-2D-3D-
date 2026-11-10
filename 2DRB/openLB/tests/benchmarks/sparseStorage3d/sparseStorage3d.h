/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2026 Adrian Kummerlaender
 *
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

// Field and operator definitions shared with the mixed-mode device codegen
// (olbcuda.cu includes this via OLB_MIXED_MODE_INCLUDE_LOCAL_HEADERS).

#ifndef SPARSE_STORAGE_3D_H
#define SPARSE_STORAGE_3D_H

struct DENSEVEL  : public olb::descriptors::FIELD_BASE<3> { };
struct SPARSEVEL : public olb::descriptors::FIELD_BASE<3> {
  static constexpr olb::FieldStorage storage = olb::FieldStorage::Sparse;
};

using BENCH_DESCRIPTOR = olb::descriptors::D3Q19<DENSEVEL,SPARSEVEL>;

/// Isolated stage so only this operator runs when timed
struct TouchStage { };

/// Cheap read + arithmetic + write so storage access dominates the kernel
template <typename FIELD>
struct TouchO {
  static constexpr olb::OperatorScope scope = olb::OperatorScope::PerCell;
  int getPriority() const { return 0; }
  template <typename CELL>
  void apply(CELL& cell) any_platform {
    auto v = cell.template getField<FIELD>();
    v[0] = v[0] + v[1] + v[2] + typename CELL::value_t{1};
    v[1] = v[1] - v[2];
    v[2] = v[0] + typename CELL::value_t{0.5};
    cell.template setField<FIELD>(v);
  }
};

/// Same as TouchO but also performs the runtime 3D index decompose (2 divides +
/// 2 mods, exactly the getLatticeR a 3D-tile accessor would do) and feeds the
/// result into the output so it is not eliminated. Run on the dense path, the
/// delta vs TouchO is the decompose cost under bandwidth (Stage-1 go/no-go).
template <typename FIELD>
struct TouchDecomposeO {
  static constexpr olb::OperatorScope scope = olb::OperatorScope::PerCell;
  int getPriority() const { return 0; }
  template <typename CELL>
  void apply(CELL& cell) any_platform {
    const auto r = cell.getLatticeR();
    auto v = cell.template getField<FIELD>();
    v[0] = v[0] + v[1] + v[2] + typename CELL::value_t{1}
         + typename CELL::value_t((r[0] ^ r[1] ^ r[2]) & 1);
    v[1] = v[1] - v[2];
    v[2] = v[0] + typename CELL::value_t{0.5};
    cell.template setField<FIELD>(v);
  }
};

#endif
