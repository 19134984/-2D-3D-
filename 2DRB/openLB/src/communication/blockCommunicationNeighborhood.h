/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Adrian Kummerlaender, Mathias J. Krause
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

#ifndef BLOCK_COMMUNICATION_NEIGHBORHOOD_H
#define BLOCK_COMMUNICATION_NEIGHBORHOOD_H

#include "core/platform/column.h"
#include "core/blockStructure.h"
#include "descriptor/fields.h"
#include "communication/mpiRequest.h"
#include "communication/loadBalancer.h"
#include "communication/communicatable.h"
#include "utilities/aliases.h"

#include <map>
#include <set>
#include <memory>
#include <functional>

namespace olb {

// *INDENT-OFF*

template<typename T> class SuperCommunicationTagCoordinator;

/// Per-field communication accounting for a single neighbor, for introspection
struct FieldCommSummary {
  /// Cells actually exchanged for the field (mask non-zero cells, or all for dense)
  std::size_t effectiveCells;
  /// Summed carried populations across the exchanged cells (0 marks a dense field)
  std::size_t carriedPopulations;
};

/// Configurable overlap communication neighborhood of a block
/**
 * Managed by SuperCommunicator
 **/
template<typename T, unsigned D>
class BlockCommunicationNeighborhood {
private:
  CuboidDecomposition<T,D>& _cuboidDecomposition;
  LoadBalancer<T>& _loadBalancer;

  const int _iC;
  const int _padding;

  std::vector<std::type_index> _fieldsRequested;
  cpu::sisd::Column<bool>      _fieldsAvailable;
  std::map<int, std::vector<std::type_index>> _fieldsCommonWith;

  /// Per-field component selectors (currently supported for POPULATION only)
  /**
   * Evaluated per requested overlap cell during maintain(); arguments are
   * the cell's block-local latticeR and the block's core extent. A field
   * without selector (or also requested densely) is exchanged densely.
   **/
  std::map<std::type_index,
           std::function<std::uint64_t(LatticeR<D>,LatticeR<D>)>> _fieldSelectors;
  /// Fields explicitly requested densely — dense overrides any selector
  std::set<std::type_index> _fieldsDenselyRequested;

  std::map<int, std::vector<CellID>> _cellsInboundFrom;
  std::map<int, std::vector<CellID>> _cellsOutboundTo;
  std::map<int, std::vector<CellID>> _cellsRequestedFrom;

  /// Per-field, per-cell component masks, index-aligned with _cellsInboundFrom
  /// resp. _cellsRequestedFrom (shared permutation, see maintain())
  std::map<int, overlap::FieldMasks> _masksInboundFrom;
  /// Masks negotiated for outbound cells (received counterpart of the
  /// neighbor's _masksInboundFrom), index-aligned with _cellsOutboundTo
  std::map<int, overlap::FieldMasks> _masksOutboundTo;

  /// Per-neighbor, per-field inbound exchange accounting, filled by maintain()
  std::map<int, std::map<std::type_index, FieldCommSummary>> _fieldCommSummary;

  /// Selector-active fields in field-request order (masked, not dense-overridden)
  std::vector<std::type_index> selectorActiveFields() const {
    std::vector<std::type_index> active;
    for (const auto& field : _fieldsRequested) {
      if (_fieldsDenselyRequested.find(field) == _fieldsDenselyRequested.end()
          && _fieldSelectors.find(field) != _fieldSelectors.end()) {
        active.push_back(field);
      }
    }
    return active;
  }

#ifdef PARALLEL_MODE_MPI
  MPI_Comm _neighborhoodComm;

  /// Flattened per-neighbor mask buffer backing the async mask send
  std::map<int, std::vector<std::uint64_t>> _maskSendBuffer;

  std::map<int, std::unique_ptr<MpiSendRequest>> _fieldRequests;
  std::map<int, std::unique_ptr<MpiSendRequest>> _cellsRequests;
  std::map<int, std::unique_ptr<MpiSendRequest>> _maskRequests;
#endif

public:
  BlockCommunicationNeighborhood(  CuboidDecomposition<T,D>& cuboidDecomposition
                                 , LoadBalancer<T>& loadBalancer
                                 , int iC
                                 , int padding
#ifdef PARALLEL_MODE_MPI
                                 , MPI_Comm comm
#endif
                                 );

  /// Request field and provides local availability
  template <typename FIELD>
  void requestField() {
    if (std::find(_fieldsRequested.begin(), _fieldsRequested.end(), typeid(FIELD)) == _fieldsRequested.end()) {
      _fieldsRequested.emplace_back(typeid(FIELD));
      _fieldsAvailable.resize(_fieldsRequested.size());
    }
    _fieldsDenselyRequested.insert(typeid(FIELD));
  }

  /// Request field restricted to a per-cell component selector
  /**
   * selector(latticeR, core) must return the component bitmask to exchange
   * for the overlap cell at block-local latticeR of a block with the given
   * core extent (e.g. overlap::streamingMask<DESCRIPTOR>() for directional
   * streaming exchange). A mask bit selects a streaming direction of a
   * propagatable field resp. one component of any other field (which must
   * fit the 64-bit mask, see overlap::maskedFieldTraits). Multiple selectors
   * for the same field merge by bitwise OR; a dense request of the same
   * field overrides all selectors.
   **/
  template <typename FIELD, typename SELECTOR>
  void requestField(SELECTOR&& selector) {
    if (std::find(_fieldsRequested.begin(), _fieldsRequested.end(), typeid(FIELD)) == _fieldsRequested.end()) {
      _fieldsRequested.emplace_back(typeid(FIELD));
      _fieldsAvailable.resize(_fieldsRequested.size());
    }
    auto& stored = _fieldSelectors[typeid(FIELD)];
    if (stored) {
      stored = [previous=stored, next=std::function<std::uint64_t(LatticeR<D>,LatticeR<D>)>(
                  std::forward<SELECTOR>(selector))](LatticeR<D> latticeR, LatticeR<D> core) {
        return previous(latticeR, core) | next(latticeR, core);
      };
    } else {
      stored = std::forward<SELECTOR>(selector);
    }
  }

  /// Returns true iff any field is exchanged with a component selector
  bool hasComponentSelectors() const {
    for (const auto& [field, selector] : _fieldSelectors) {
      if (_fieldsDenselyRequested.find(field) == _fieldsDenselyRequested.end()) {
        return true;
      }
    }
    return false;
  }
  /// Request individual cell for communication
  void requestCell(LatticeR<D> latticeR);
  void requestCell(LatticeR<D> latticeR, int* lastC);
  /// Request all cells in overlap of size width for communication
  void requestOverlap(int width);
  /// Request all indicated cells in overlap of width for communication
  void requestOverlap(int width, BlockIndicatorF<T,D>& indicatorF);

  /// Remove all requested cells
  void clearRequestedCells();

  /// Update local availability of previously requested field
  void setFieldAvailability(std::type_index field, bool available);
  /// Update outbound availabilities for locally available neighbor block
  /**
   * Used for local non-MPI exchange of field availabilities
   **/
  template <typename BLOCK>
  void setFieldsAvailability(int iC, BLOCK& block);

  void maintain();

#ifdef PARALLEL_MODE_MPI
  void send(SuperCommunicationTagCoordinator<T>&);
  void receive(SuperCommunicationTagCoordinator<T>&);
  void wait();
#endif

  /// Calls f(iC) for every neighboring cuboid ID iC
  template <typename F>
  void forNeighbors(F f) const {
    for (const auto& [iC, _] : _cellsRequestedFrom) {
      f(iC);
    }
  }

  /// All fields requested from this neighborhood (dense or selector)
  const std::vector<std::type_index>& getFieldsRequested() const {
    return _fieldsRequested;
  }

  /// Calls f(field, FieldCommSummary) for each field exchanged with neighbor iC
  template <typename F>
  void forFieldCommSummary(int iC, F f) const {
    const auto it = _fieldCommSummary.find(iC);
    if (it != _fieldCommSummary.end()) {
      for (const auto& [field, summary] : it->second) {
        f(field, summary);
      }
    }
  }

  const std::vector<std::type_index>& getFieldsCommonWith(int iC) const;

  const std::vector<CellID>& getCellsOutboundTo(int iC) const;
  const std::vector<CellID>& getCellsInboundFrom(int iC) const;
  const std::vector<CellID>& getCellsRequestedFrom(int iC) const;

  /// Per-field masks aligned with getCellsInboundFrom (empty if dense)
  const overlap::FieldMasks& getMasksInboundFrom(int iC) const;
  /// Per-field masks aligned with getCellsOutboundTo (empty if dense)
  /**
   * Negotiated: the neighbor's inbound masks for the cells it requested.
   * In the non-MPI/local case read the requesting neighborhood's
   * getMasksInboundFrom instead (mirroring getCellsRequestedFrom usage).
   **/
  const overlap::FieldMasks& getMasksOutboundTo(int iC) const;

};

// *INDENT-ON*

}

#endif
