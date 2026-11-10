/*  This file is part of the OpenLB library
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

#ifndef DESCRIPTOR_DEFINITION_VECTOR_H
#define DESCRIPTOR_DEFINITION_VECTOR_H

namespace olb {

namespace concepts {

template <typename DESCRIPTOR>
concept VectorLatticeDescriptor = requires {
  { DESCRIPTOR::d } -> std::convertible_to<std::size_t>;
  requires std::is_base_of_v<meta::list_base, DESCRIPTOR>;
  /// Count of discrete velocitity components is defined
  { DESCRIPTOR::q } -> std::convertible_to<std::size_t>;
  /// Count of dimension per discrete velocity is defined
  { DESCRIPTOR::n } -> std::convertible_to<std::size_t>;
};

}

namespace descriptors {

template <typename DESCRIPTOR, unsigned N, typename... FIELDS>
struct VECTOR : public VECTOR_LATTICE_DESCRIPTOR<DESCRIPTOR::d,DESCRIPTOR::q,N,POPULATION,FIELDS...> {
  VECTOR() = delete;

  using base_descriptor_t = DESCRIPTOR;
};

template <concepts::VectorLatticeDescriptor DESCRIPTOR>
constexpr int d() any_platform
{
  return DESCRIPTOR::d;
}

template <concepts::VectorLatticeDescriptor DESCRIPTOR>
constexpr int q() any_platform
{
  return DESCRIPTOR::q;
}

/// Map vector index to linear population index
template <concepts::VectorLatticeDescriptor DESCRIPTOR>
constexpr int q(unsigned iPopV, unsigned iD) any_platform
{
  return iPopV*DESCRIPTOR::n + iD;
}

template <concepts::VectorLatticeDescriptor DESCRIPTOR>
constexpr int vicinity() any_platform
{
  return vicinity<DESCRIPTOR::d, DESCRIPTOR::q / DESCRIPTOR::n>();
}

template <concepts::VectorLatticeDescriptor DESCRIPTOR>
constexpr int c(unsigned iPop, unsigned iDim) any_platform
{
  return c<DESCRIPTOR::d, DESCRIPTOR::q / DESCRIPTOR::n>(iPop / DESCRIPTOR::n, iDim);
}

template <concepts::VectorLatticeDescriptor DESCRIPTOR>
constexpr Vector<int,DESCRIPTOR::d> c(unsigned iPop) any_platform
{
  return c<DESCRIPTOR::d, DESCRIPTOR::q / DESCRIPTOR::n>(iPop / DESCRIPTOR::n);
}

template <concepts::VectorLatticeDescriptor DESCRIPTOR>
constexpr Vector<int,DESCRIPTOR::d> cVector(unsigned iPopV) any_platform
{
  return c<DESCRIPTOR::d, DESCRIPTOR::q / DESCRIPTOR::n>(iPopV);
}

template <concepts::VectorLatticeDescriptor DESCRIPTOR>
constexpr int opposite(unsigned iPop) any_platform
{
  return opposite<DESCRIPTOR::d, DESCRIPTOR::q / DESCRIPTOR::n>(iPop / DESCRIPTOR::n);
}

template <concepts::VectorLatticeDescriptor DESCRIPTOR>
constexpr int oppositeVector(unsigned iPopV) any_platform
{
  return opposite<DESCRIPTOR::d, DESCRIPTOR::q / DESCRIPTOR::n>(iPopV);
}

template <typename T, concepts::VectorLatticeDescriptor DESCRIPTOR>
constexpr T invCs2() any_platform
{
  return invCs2<T, DESCRIPTOR::d, DESCRIPTOR::q / DESCRIPTOR::n>();
}

template <typename T, concepts::VectorLatticeDescriptor DESCRIPTOR>
constexpr T t(unsigned iPop) any_platform
{
  return t<T, DESCRIPTOR::d, DESCRIPTOR::q / DESCRIPTOR::n>(iPop, descriptors::tag::DEFAULT{});
}

}  // namespace descriptors

}  // namespace olb

#endif
