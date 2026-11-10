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

// Test-local second propagatable field for the two-field masked exchange proof.
// Shared between the host TU and the mixed-mode olbcuda.cu (pulled in via
// OLB_MIXED_MODE_INCLUDE_LOCAL_HEADERS) so the descriptor's GPU instantiations
// resolve the type. Shares the base descriptor stencil (same q, same c_i).

#ifndef TESTS_COMMUNICATION_TWO_POP_H
#define TESTS_COMMUNICATION_TWO_POP_H

#include <olb.h>

struct POPULATION2 : olb::descriptors::PROPAGATABLE_FIELD_BASE { };

#endif
