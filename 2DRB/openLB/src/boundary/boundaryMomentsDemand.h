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

#ifndef OLB_BOUNDARY_MOMENTS_DEMAND_H
#define OLB_BOUNDARY_MOMENTS_DEMAND_H

#include "core/momentsDemand.h"

namespace olb {

template <typename DESCRIPTOR, int direction, int orientation>
std::vector<MomentsDemand<DESCRIPTOR::d>> planeFdMomentsDemand()
{
  using O = Vector<std::int32_t,DESCRIPTOR::d>;
  O n1{}; n1[direction] =   -orientation;
  O n2{}; n2[direction] = -2*orientation;
  std::vector<O> velocityNeighbours{n1, n2};
  for (unsigned iD=0; iD < DESCRIPTOR::d; ++iD) {
    if (iD != static_cast<unsigned>(direction)) {
      O tp{}; tp[iD] =  1;
      O tm{}; tm[iD] = -1;
      velocityNeighbours.push_back(tp);
      velocityNeighbours.push_back(tm);
    }
  }
  return {
    {MomentSet::RhoU, Freshness::Fresh, {O{}}},
    {MomentSet::U,    Freshness::Fresh, velocityNeighbours},
  };
}

template <typename DESCRIPTOR, int plane, int normal1, int normal2>
std::vector<MomentsDemand<DESCRIPTOR::d>> outerVelocityEdgeMomentsDemand()
{
  using O = Vector<std::int32_t,DESCRIPTOR::d>;
  constexpr int direction1 = (plane+1)%3;
  constexpr int direction2 = (plane+2)%3;
  O d1a{}; d1a[direction1] =   -normal1;
  O d1b{}; d1b[direction1] = -2*normal1;
  O d2a{}; d2a[direction2] =   -normal2;
  O d2b{}; d2b[direction2] = -2*normal2;
  O ep{};  ep[plane] =  1;
  O em{};  em[plane] = -1;
  return {
    {MomentSet::U,    Freshness::Fresh, {O{}, ep, em}},
    {MomentSet::RhoU, Freshness::Fresh, {d1a, d1b, d2a, d2b}},
  };
}

template <unsigned D>
std::vector<MomentsDemand<D>> outerVelocityCornerMomentsDemand(Vector<std::int32_t,D> normal)
{
  using O = Vector<std::int32_t,D>;
  std::vector<O> rhoU;
  for (unsigned iD=0; iD < D; ++iD) {
    O a{}; a[iD] =   -normal[iD];
    O b{}; b[iD] = -2*normal[iD];
    rhoU.push_back(a);
    rhoU.push_back(b);
  }
  return {
    {MomentSet::U,    Freshness::Fresh, {O{}}},
    {MomentSet::RhoU, Freshness::Fresh, rhoU},
  };
}

}

#endif
