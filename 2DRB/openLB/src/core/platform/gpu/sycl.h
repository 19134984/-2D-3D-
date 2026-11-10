/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Adrian Kummerlaender
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

#ifndef CORE_PLATFORM_GPU_SYCL_H
#define CORE_PLATFORM_GPU_SYCL_H

#include "sycl/column.h"
#include "sycl/communicator.h"
#include "sycl/dynamics.h"
#include "sycl/context.h"
#include "sycl/device.h"
#include "sycl/mask.h"
#include "sycl/operator.h"
#include "sycl/registry.h"
#include "sycl/statistics.h"

#ifdef OLB_PRECOMPILING_SYCL
#include "sycl/precompiled.hh"
#endif

#endif
