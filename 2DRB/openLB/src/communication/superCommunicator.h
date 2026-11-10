/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Adrian Kummerlaender
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

#ifndef SUPER_COMMUNICATOR_H
#define SUPER_COMMUNICATOR_H

#include <set>
#include <vector>
#include <algorithm>
#include <map>
#include <string_view>
#include <typeindex>

#include "core/meta.h"
#include "mpiManager.h"
#include "loadBalancer.h"
#include "blockCommunicator.h"
#include "blockCommunicationNeighborhood.h"
#include "superCommunicationTagCoordinator.h"
#include "utilities/functorPtr.h"

namespace olb {


/// Generic communicator for overlaps between blocks of SUPER
/**
 * This class provides inter-block communication for any
 * super structures containing overlapping blocks.
 *
 * SUPER must expose load balancer via getLoadBalancer
 * SUPER must expose cuboid geometry via getCuboidDecomposition
 *
 * Blocks must be exposed via SUPER::getBlock
 * Blocks must implement FieldsCommunicatable
 **/
template <typename T, typename SUPER>
class SuperCommunicator {
private:
  SUPER& _super;

#ifdef PARALLEL_MODE_MPI
  SuperCommunicationTagCoordinator<T> _tagCoordinator;
  MPI_Comm _neighborhoodComm; /// Communicator for neighborhood negotation
  MPI_Comm _communicatorComm; /// Communicator for overlap communication
#endif

  /// Neighborhood negotiation
  std::vector<std::unique_ptr<BlockCommunicationNeighborhood<T,SUPER::d>>> _blockNeighborhoods;
  /// Per-block communicators constructed to satify requested exchanges
  std::vector<std::unique_ptr<BlockCommunicator>> _blockCommunicators;

  /// List of requested FIELDS
  std::vector<std::type_index> _fieldsRequested;
  /// Demangled name of each requested field, for introspection
  std::map<std::type_index, std::string_view> _fieldNames;
  /// Set of non-local neighbor cuboids
  std::set<int> _remoteCuboidNeighborhood;

  /// True iff any fields or cells were requested
  bool _enabled = false;
  /// True iff requests are synced between processes
  bool _ready   = false;

public:
  SuperCommunicator(SUPER& super);
  ~SuperCommunicator();

  /// Request FIELD for communication
  /**
   * Communication of dynamic fields depends on the availability at
   * both ends of the exchange.
   **/
  template <typename FIELD>
  void requestField() {
    if (std::find(_fieldsRequested.begin(), _fieldsRequested.end(), typeid(FIELD)) == _fieldsRequested.end()) {
      _fieldsRequested.emplace_back(typeid(FIELD));
      _fieldNames.emplace(typeid(FIELD), meta::nice_name<FIELD>());
    }
    // Always forward (idempotent per block): a dense request must reach the
    // neighborhoods even if the field was previously requested with a
    // component selector, as dense overrides all selectors.
    for (auto& neighborhood : _blockNeighborhoods) {
      neighborhood->template requestField<FIELD>();
    }
  }

  /// Request FIELD restricted to a per-cell component selector
  /**
   * selector(latticeR, core) returns the component bitmask to exchange for
   * an overlap cell (block-local latticeR of a block with the given core
   * extent), e.g. overlap::streamingMask<DESCRIPTOR>() for the directional
   * streaming exchange. Selectors of the same FIELD merge by bitwise OR;
   * a dense requestField<FIELD>() overrides all its selectors.
   **/
  template <typename FIELD, typename SELECTOR>
  void requestField(SELECTOR&& selector) {
    if (std::find(_fieldsRequested.begin(), _fieldsRequested.end(), typeid(FIELD)) == _fieldsRequested.end()) {
      _fieldsRequested.emplace_back(typeid(FIELD));
      _fieldNames.emplace(typeid(FIELD), meta::nice_name<FIELD>());
    }
    for (auto& neighborhood : _blockNeighborhoods) {
      neighborhood->template requestField<FIELD>(selector);
    }
  }

  /// Convenience method for requesting multiple FIELDS in one call
  template <typename... FIELDS>
  void requestFields() {
    (requestField<FIELDS>(), ...);
  }

  /// Request single cell in the padding area for communication
  /**
   * latticeR[0] must be the local cuboid ID
   **/
  void requestCell(LatticeR<SUPER::d+1> latticeR);
  /// Request all cells in overlap of width for communication
  void requestOverlap(int width);
  /// Request all indicated cells in overlap of width for communication
  void requestOverlap(int width, FunctorPtr<SuperIndicatorF<T,SUPER::d>>&& indicatorF);

  /// Remove all requested cells
  void clearRequestedCells();

  /// Exchange requests between processes
  void exchangeRequests();

  void prepareCommunicate();
  void continueCommunicate();

  /// Perform communication
  void communicate();

  /// Returns set of non-local neighborhood cuboid indices
  const std::set<int>& getRemoteCuboids() const;

  /// Fields requested for this communication stage (for introspection)
  const std::vector<std::type_index>& getFieldsRequested() const {
    return _fieldsRequested;
  }

  /// Demangled name of a requested field (for introspection), raw name if unknown
  std::string_view getFieldName(std::type_index field) const {
    const auto it = _fieldNames.find(field);
    return it != _fieldNames.end() ? it->second : std::string_view(field.name());
  }

  /// Read-only access to a local block's communication neighborhood (for introspection)
  const BlockCommunicationNeighborhood<T,SUPER::d>& getBlockNeighborhood(int iBlock) const {
    return *_blockNeighborhoods[iBlock];
  }

};

template <typename SUPER>
SuperCommunicator(SUPER&) -> SuperCommunicator<typename SUPER::value_t,SUPER>;

}

#endif
