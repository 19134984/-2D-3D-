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

#ifndef DESCRIPTOR_CELL_PHYS_R_H
#define DESCRIPTOR_CELL_PHYS_R_H

#include "core/vector.h"
#include "fields.h"

namespace olb {

namespace fields {
namespace converter {
struct PHYS_DELTA_X;
}
}

namespace descriptors {

/// Physical position of a cell from its lattice coordinates and block origin
/**
 * Reconstructs cuboid.getPhysR(latticeR) = origin + latticeR * deltaR on host and
 * device, replacing a stored per-cell LOCATION field. Requires the operator's
 * PARAMETERS to carry BLOCK_LOWER and converter::PHYS_DELTA_X (both are supplied by
 * SuperLatticePointCoupling per block).
 **/
template <typename CELL, typename PARAMETERS>
auto getPhysR(CELL& cell, PARAMETERS& params) any_platform
{
  return params.template get<fields::BLOCK_LOWER>()
       + cell.getLatticeR() * params.template get<fields::converter::PHYS_DELTA_X>();
}

}

}

#endif
