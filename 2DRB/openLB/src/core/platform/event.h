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

#ifndef CORE_PLATFORM_EVENT_H
#define CORE_PLATFORM_EVENT_H

#include "platform.h"

namespace olb {

/// Completion token for work dispatched to PLATFORM
/**
 * Events are reusable slots: record captures completion of all work
 * previously dispatched to an execution lane while wait / barrierOn
 * always refer to the most recent record. Never-recorded and already
 * completed events are satisfied immediately, allowing consumers to
 * over-wait freely without per-stage bookkeeping.
 *
 * Models concepts::Completable on all platforms. The device
 * specializations, declared alongside each platform's device::Stream,
 * additionally model concepts::DeviceOrderable to express producer /
 * consumer relationships without host-side blocking.
 *
 * This primary template covers the host-synchronous CPU platforms on
 * which any dispatched work is complete before a consumer can observe
 * the event.
 **/
template <Platform PLATFORM>
class Event {
public:
  /// Host-blocking join, no-op as CPU block execution is synchronous
  void wait() const { }
};

// Device-specific completion tokens are declared alongside each
// platform's device::Stream (see core/platform/gpu/*/device.h).
// Declared here so that a missing include fails loudly instead of
// silently instantiating the no-op primary template.
#ifdef PLATFORM_GPU_CUDA
template <> class Event<Platform::GPU_CUDA>;
#endif

#ifdef PLATFORM_GPU_HIP
template <> class Event<Platform::GPU_HIP>;
#endif

#ifdef PLATFORM_GPU_SYCL
template <> class Event<Platform::GPU_SYCL>;
#endif

}

#endif
