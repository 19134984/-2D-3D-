/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2026 Adrian Kummerlaender
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

#ifndef FUNCTORS_PRIMITIVE_VOLUME_SAMPLERS_H
#define FUNCTORS_PRIMITIVE_VOLUME_SAMPLERS_H

#include "functors/primitive/interpolation3d.h"
#include "descriptor/cellPhysR.h"

namespace olb {

namespace functors {

namespace pointwise {

struct PhysVelocityF {
  using parameters = meta::list<fields::converter::PHYS_VELOCITY>;

  using result_t = descriptors::VELOCITY;

  template <typename CELL, typename PHYS_R, typename PARAMETERS>
  auto compute(CELL& cell, PHYS_R samplingR, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    FieldD<V,DESCRIPTOR,result_t> result{};
    const V physU = parameters.template get<fields::converter::PHYS_VELOCITY>();
    cell.computeU(result.data());
    result *= physU;
    return result;
  }
};

/// Porosity at the cell.
struct PorosityF {
  using parameters = meta::list<>;
  using result_t = descriptors::POROSITY;

  template <typename CELL, typename PHYSR, typename PARAMETERS>
  auto compute(CELL& cell, PHYSR& physR, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    FieldD<V,DESCRIPTOR,result_t> result{};
    result[0] = cell.template getField<descriptors::POROSITY>();
    return result;
  }
};

/// Physical pressure at the cell.
struct PhysPressureF {
  using parameters = meta::list<fields::converter::PHYS_PRESSURE>;
  using result_t = descriptors::PRESSURE_LEVEL;

  template <typename CELL, typename PHYSR, typename PARAMETERS>
  auto compute(CELL& cell, PHYSR& physR, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    const V physPressure = parameters.template get<fields::converter::PHYS_PRESSURE>();
    V rho{};
    Vector<V,DESCRIPTOR::d> u{};
    cell.computeRhoU(rho, u.data());
    FieldD<V,DESCRIPTOR,result_t> result{};
    result[0] = util::pressureFromDensity<V,DESCRIPTOR>(rho) * physPressure;
    return result;
  }
};

/// Physical vorticity magnitude at the cell (central-difference curl, per-level dx scaling).
struct VorticityMagnitudeF {
  using parameters = meta::list<fields::converter::PHYS_VELOCITY,
                                fields::converter::PHYS_DELTA_X>;
  using result_t = descriptors::VORTICITY;

  template <typename CELL, typename PHYSR, typename PARAMETERS>
  auto compute(CELL& cell, PHYSR& physR, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    const V physVelocity = parameters.template get<fields::converter::PHYS_VELOCITY>();
    const V physDeltaX   = parameters.template get<fields::converter::PHYS_DELTA_X>();
    Vector<V,DESCRIPTOR::d> gradU[DESCRIPTOR::d];
    for (int axis = 0; axis < DESCRIPTOR::d; ++axis) {
      Vector<int,DESCRIPTOR::d> e{};
      V rho{};
      Vector<V,DESCRIPTOR::d> uPlus{}, uMinus{};
      e[axis] = 1;
      cell.neighbor(e).computeRhoU(rho, uPlus.data());
      e[axis] = -1;
      cell.neighbor(e).computeRhoU(rho, uMinus.data());
      gradU[axis] = V{0.5} * (uPlus - uMinus);
    }
    const V curlScale = physVelocity / physDeltaX;
    Vector<V,DESCRIPTOR::d> vorticity{
      curlScale * (gradU[1][2] - gradU[2][1]),
      curlScale * (gradU[2][0] - gradU[0][2]),
      curlScale * (gradU[0][1] - gradU[1][0])
    };
    FieldD<V,DESCRIPTOR,result_t> result{};
    result[0] = util::norm<DESCRIPTOR::d>(vorticity);
    return result;
  }
};

namespace detail {

template <typename... Ls>
struct merge_all;
template <>
struct merge_all<> {
  using type = meta::list<>;
};
template <typename L>
struct merge_all<L> {
  using type = L;
};
template <typename L, typename... Rest>
struct merge_all<L,Rest...> {
  using type = meta::merge<L, typename merge_all<Rest...>::type>;
};

}

/// Samples several value-at-cell functors and concatenates their result fields.
/// The cell-functor analogue of AnalyticalComposed: the combined result_field is the
/// concatenation of the wrapped fields (in order), the parameters are their merged union.
template <typename... Fs>
struct Composed {
  using result_t = meta::list<typename Fs::result_t...>;

  using parameters = typename detail::merge_all<typename Fs::parameters...>::type;

  template <typename CELL, typename PHYSR, typename PARAMETERS>
  auto compute(CELL& cell, PHYSR& physR, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    FieldD<V,DESCRIPTOR,descriptors::packed_field_t<result_t>> result{};
    unsigned offset = 0;
    ([&] {
      const auto sub = Fs().compute(cell, physR, parameters);
      constexpr unsigned n = Fs::result_t::template size<DESCRIPTOR>();
      for (unsigned i = 0; i < n; ++i) {
        result[offset + i] = sub[i];
      }
      offset += n;
    }(), ...);
    return result;
  }
};

/// Trilinear interpolation wrapper: turns a value-at-cell functor (compute returning values
/// AT the containing cell) into an on-plane point functor by interpolating the wrapped values
/// over the 8 cells surrounding physR. Positions come only from getPhysR -- no hand offset.
/// Any functor already implementing compute(cell, physR, parameters) can be used directly
/// instead (the escape hatch), wrapped or not.
template <typename VALUE_F>
struct Trilinear {
  using parameters   = meta::merge<typename VALUE_F::parameters,
                                   meta::list<fields::converter::PHYS_DELTA_X, fields::BLOCK_LOWER>>;
  using result_t = typename VALUE_F::result_t;

  template <typename CELL, typename PHYSR, typename PARAMETERS>
  auto compute(CELL& cell, PHYSR& physR, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;
    constexpr unsigned N = descriptors::packed_field_t<result_t>::template size<DESCRIPTOR>();
    const V dx = parameters.template get<fields::converter::PHYS_DELTA_X>();
    const auto cellR = descriptors::getPhysR(cell, parameters);
    Vector<V,3> distance{};
    for (unsigned iD = 0; iD < 3; ++iD) {
      distance[iD] = (physR[iD] - cellR[iD]) / dx;
    }
    FieldD<V,DESCRIPTOR,descriptors::packed_field_t<result_t>> result{};
    auto f = [&parameters,&physR](auto nCell, V* res) -> void {
      const auto v = VALUE_F().compute(nCell, physR, parameters);
      for (unsigned i = 0; i < N; ++i) {
        res[i] = v[i];
      }
    };
    interpolate3d<N>(cell, result, distance, f);
    return result;
  }
};


}

}

}

#endif
