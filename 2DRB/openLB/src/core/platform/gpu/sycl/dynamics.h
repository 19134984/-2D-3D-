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

#ifndef GPU_SYCL_DYNAMICS_H
#define GPU_SYCL_DYNAMICS_H

#include "descriptor/base.h"

#include <cstdint>

namespace olb {

namespace gpu {

namespace sycl {

/// Virtual interface for device-side dynamically-dispatched dynamics access
template <typename T, typename DESCRIPTOR>
struct Dynamics;

/// On-device per-cell index into the per-block dynamics dispatch table (0 == nullptr)
template <typename T, typename DESCRIPTOR>
struct DYNAMICS : public descriptors::TYPED_FIELD_BASE<std::uint16_t,1> { };

}

}

}

#endif
