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

#ifndef BLOCK_COMMUNICATION_NEIGHBORHOOD_HH
#define BLOCK_COMMUNICATION_NEIGHBORHOOD_HH

#include <stdexcept>
#include <algorithm>
#include <optional>
#include <bit>

#include "blockCommunicationNeighborhood.h"

#include "utilities/aliases.h"

namespace olb {

// *INDENT-OFF*

template <typename T, unsigned D>
BlockCommunicationNeighborhood<T,D>::BlockCommunicationNeighborhood(
    CuboidDecomposition<T,D>& cuboidDecomposition
  , LoadBalancer<T>& loadBalancer
  , int iC
  , int padding
#ifdef PARALLEL_MODE_MPI
  , MPI_Comm comm
#endif
):
    _cuboidDecomposition(cuboidDecomposition)
  , _loadBalancer(loadBalancer)
  , _iC(iC)
  , _padding(padding)
#ifdef PARALLEL_MODE_MPI
  , _neighborhoodComm(comm)
#endif
{
  // Ensure that any neighboring cuboids are interacted with during request negotiation
  requestOverlap(1);
  for (auto& [_, cells] : _cellsInboundFrom) {
    cells.clear();
  }
  for (auto& [_, cells] : _cellsRequestedFrom) {
    cells.clear();
  }
}

template <typename T, unsigned D>
void BlockCommunicationNeighborhood<T,D>::requestCell(LatticeR<D> latticeR)
{
  const auto& cuboid = _cuboidDecomposition.get(_iC);
  BlockStructureD<D> paddedBlock(cuboid.getExtent(), _padding);
  if (paddedBlock.isPadding(latticeR)) {
    // Read physR using global cuboid geometry instead of local cuboid
    // to resolve periodic boundaries
    auto physR = _cuboidDecomposition.getPhysR(latticeR.withPrefix(_iC));
    if (auto remoteLatticeR = _cuboidDecomposition.getSurroundingLatticeR(physR, _iC, nullptr)) {
      _cellsInboundFrom[(*remoteLatticeR)[0]].emplace_back(paddedBlock.getCellId(latticeR));

      const auto& remoteCuboid = _cuboidDecomposition.get((*remoteLatticeR)[0]);
      BlockStructureD<D> remotePaddedBlock(remoteCuboid.getExtent(), _padding);

      _cellsRequestedFrom[(*remoteLatticeR)[0]].emplace_back(
        remotePaddedBlock.getCellId(&(*remoteLatticeR)[1]));
    }
  } else {
    throw std::logic_error("Requested cell is outside of available padding");
  }
}

template <typename T, unsigned D>
void BlockCommunicationNeighborhood<T,D>::requestCell(LatticeR<D> latticeR, int* lastC)
{
  const auto& cuboid = _cuboidDecomposition.get(_iC);
  BlockStructureD<D> paddedBlock(cuboid.getExtent(), _padding);
  if (paddedBlock.isPadding(latticeR)) {
    // Read physR using global cuboid geometry instead of local cuboid
    // to resolve periodic boundaries
    auto physR = _cuboidDecomposition.getPhysR(latticeR.withPrefix(_iC));
    if (auto remoteLatticeR = _cuboidDecomposition.getSurroundingLatticeR(physR, _iC, lastC)) {
      _cellsInboundFrom[(*remoteLatticeR)[0]].emplace_back(paddedBlock.getCellId(latticeR));

      const auto& remoteCuboid = _cuboidDecomposition.get((*remoteLatticeR)[0]);
      BlockStructureD<D> remotePaddedBlock(remoteCuboid.getExtent(), _padding);

      _cellsRequestedFrom[(*remoteLatticeR)[0]].emplace_back(
        remotePaddedBlock.getCellId(&(*remoteLatticeR)[1]));
    }
  } else {
    throw std::logic_error("Requested cell is outside of available padding");
  }
}

template <typename T, unsigned D>
void BlockCommunicationNeighborhood<T,D>::requestOverlap(int width)
{
  if (width < 1) {
    throw std::logic_error("Overlap requests must have width >= 1");
  }
  if (width > _padding) {
    throw std::logic_error("Requested overlap exceeds available padding");
  }

  const auto& cuboid = _cuboidDecomposition.get(_iC);
  BlockStructureD<D> overlapBlock(cuboid.getExtent(), width);

  int lastC = -1;
  overlapBlock.forPaddingLocations(width, [&](LatticeR<D> latticeR) {
    requestCell(latticeR, &lastC);
  });
}

template <typename T, unsigned D>
void BlockCommunicationNeighborhood<T,D>::requestOverlap(int width, BlockIndicatorF<T,D>& indicatorF)
{
  if (width < 1) {
    throw std::logic_error("Overlap requests must have width >= 1");
  }
  if (width > _padding) {
    throw std::logic_error("Requested overlap exceeds available padding");
  }

  const auto& cuboid = _cuboidDecomposition.get(_iC);
  BlockStructureD<D> overlapBlock(cuboid.getExtent(), width);

  int lastC = -1;
  overlapBlock.forPaddingLocations(width, [&](LatticeR<D> latticeR) {
    if (indicatorF(latticeR)) {
      requestCell(latticeR, &lastC);
    }
  });
}

template <typename T, unsigned D>
void BlockCommunicationNeighborhood<T,D>::clearRequestedCells()
{
  _cellsInboundFrom.clear();
  _cellsOutboundTo.clear();
  _cellsRequestedFrom.clear();
  _masksInboundFrom.clear();
  _masksOutboundTo.clear();
}

template <typename T, unsigned D>
void BlockCommunicationNeighborhood<T,D>::setFieldAvailability(std::type_index field, bool available)
{
  auto iter = std::find(_fieldsRequested.begin(), _fieldsRequested.end(), field);
  if (iter != _fieldsRequested.end()) {
    _fieldsAvailable[iter - _fieldsRequested.begin()] = available;
  }
}

template <typename T, unsigned D>
template <typename BLOCK>
void BlockCommunicationNeighborhood<T,D>::setFieldsAvailability(int iC, BLOCK& block)
{
  auto& fieldsCommonWith = _fieldsCommonWith[iC];
  fieldsCommonWith.clear();
  for (unsigned iField=0; iField < _fieldsRequested.size(); ++iField) {
    auto field = _fieldsRequested[iField];
    if (block.hasCommunicatable(field) && _fieldsAvailable[iField]) {
      fieldsCommonWith.emplace_back(field);
    }
  }
}

template <typename T, unsigned D>
void BlockCommunicationNeighborhood<T,D>::maintain()
{
  for (auto& [iC, _] : _cellsInboundFrom) {
    auto& iCells = _cellsInboundFrom[iC];
    auto& oCells = _cellsRequestedFrom[iC];
    // Compute sorted indices w.r.t. cell IDs in iCells
    std::vector<std::ptrdiff_t> p(iCells.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(), [&iCells](auto i, auto j) { return iCells[i] < iCells[j]; });
    // Reduce indices to only include unique values of iCells
    auto pu = std::unique(p.begin(), p.end(), [&iCells](auto i, auto j) { return iCells[i] == iCells[j]; });
    p.erase(pu, p.end());
    // Copy cell IDs at computed indices to new in and out lists
    {
      std::vector<CellID> buffer(p.size());
      std::transform(p.begin(), p.end(), buffer.begin(), [&iCells](auto i) { return iCells[i]; });
      iCells = buffer;
    }
    {
      std::vector<CellID> buffer(p.size());
      std::transform(p.begin(), p.end(), buffer.begin(), [&oCells](auto i) { return oCells[i]; });
      oCells = buffer;
    }
  }

  // Evaluate component selectors into per-cell masks aligned with the
  // (deduplicated, sorted) inbound cell lists. Masks are receiver-computed:
  // only this block knows its own core geometry, the provider obeys the
  // negotiated masks (mirroring how the cell lists themselves work).
  _masksInboundFrom.clear();
  if (hasComponentSelectors()) {
    const auto& cuboid = _cuboidDecomposition.get(_iC);
    BlockStructureD<D> paddedBlock(cuboid.getExtent(), _padding);
    const LatticeR<D> core(cuboid.getExtent());
    const auto activeFields = selectorActiveFields();
    for (auto& [iC, iCells] : _cellsInboundFrom) {
      auto& fieldMasks = _masksInboundFrom[iC];
      for (const auto& field : activeFields) {
        const auto& selector = _fieldSelectors.at(field);
        std::vector<std::uint64_t> masks(iCells.size());
        for (std::size_t i=0; i < iCells.size(); ++i) {
          masks[i] = selector(paddedBlock.getLatticeR(iCells[i]), core);
        }
        fieldMasks.emplace_back(field, std::move(masks));
      }
    }
  }

  // Per-field inbound exchange accounting for introspection. Each field's own
  // selector is evaluated independently (dense fields exchange every cell), so
  // this reflects true per-field cell and population counts even where the
  // shared mask channel above collapses distinct selectors into one.
  _fieldCommSummary.clear();
  {
    std::optional<BlockStructureD<D>> paddedBlock;
    LatticeR<D> core{};
    if (hasComponentSelectors()) {
      const auto& cuboid = _cuboidDecomposition.get(_iC);
      paddedBlock.emplace(cuboid.getExtent(), _padding);
      core = LatticeR<D>(cuboid.getExtent());
    }
    for (const auto& [iC, iCells] : _cellsInboundFrom) {
      auto& summary = _fieldCommSummary[iC];
      for (const auto& field : _fieldsRequested) {
        const bool dense = _fieldsDenselyRequested.find(field) != _fieldsDenselyRequested.end();
        const auto selectorIt = _fieldSelectors.find(field);
        if (dense || selectorIt == _fieldSelectors.end()) {
          summary[field] = FieldCommSummary{iCells.size(), 0};
        } else {
          std::size_t cells = 0;
          std::size_t populations = 0;
          for (const auto cellId : iCells) {
            const std::uint64_t mask = selectorIt->second(paddedBlock->getLatticeR(cellId), core);
            if (mask) {
              ++cells;
              populations += std::popcount(mask);
            }
          }
          summary[field] = FieldCommSummary{cells, populations};
        }
      }
    }
  }
}

#ifdef PARALLEL_MODE_MPI

template <typename T, unsigned D>
void BlockCommunicationNeighborhood<T,D>::send(SuperCommunicationTagCoordinator<T>& coordinator)
{
  for (auto& [iC, cells] : _cellsRequestedFrom) {
    _fieldRequests[iC] = std::make_unique<MpiSendRequest>(
      _fieldsAvailable.data(), _fieldsRequested.size(),
      _loadBalancer.rank(iC), coordinator.get(_iC, iC, 0), _neighborhoodComm);
    _fieldRequests[iC]->start();

    _cellsRequests[iC] = std::make_unique<MpiSendRequest>(
      cells.data(), cells.size(),
      _loadBalancer.rank(iC), coordinator.get(_iC, iC, 1), _neighborhoodComm);
    _cellsRequests[iC]->start();

    // Component masks travel alongside the requested cells. Presence is
    // deducible on both sides: requests are stage-symmetric, so every rank
    // registered the same selectors for this communicator.
    if (hasComponentSelectors()) {
      // Flatten the per-field masks into one contiguous channel-2 message; the
      // receiver splits by the (stage-symmetric) selector-active field order.
      auto& buffer = _maskSendBuffer[iC];
      buffer.clear();
      for (const auto& [field, masks] : _masksInboundFrom[iC]) {
        buffer.insert(buffer.end(), masks.begin(), masks.end());
      }
      _maskRequests[iC] = std::make_unique<MpiSendRequest>(
        buffer.data(), buffer.size(),
        _loadBalancer.rank(iC), coordinator.get(_iC, iC, 2), _neighborhoodComm);
      _maskRequests[iC]->start();
    }
  }
}

template <typename T, unsigned D>
void BlockCommunicationNeighborhood<T,D>::receive(SuperCommunicationTagCoordinator<T>& coordinator)
{
  forNeighbors([&](int iC) {
    std::unique_ptr<bool[]> fieldsAvailable(new bool[_fieldsRequested.size()] { });
    auto& fieldsCommonWith = _fieldsCommonWith[iC];
    fieldsCommonWith.clear();

    singleton::mpi().receive(
      fieldsAvailable.get(), _fieldsRequested.size(),
      _loadBalancer.rank(iC),
      coordinator.get(iC, _iC, 0),
      _neighborhoodComm);

    for (unsigned iField=0; iField < _fieldsRequested.size(); ++iField) {
      if (fieldsAvailable[iField] && _fieldsAvailable[iField]) {
        fieldsCommonWith.emplace_back(_fieldsRequested[iField]);
      }
    }

    static_assert(sizeof(CellID) == 4 || sizeof(CellID) == 8,
                  "CellID must be either 4 or 8 byte unsigned integer");
    std::size_t newOutboundCount = singleton::mpi().probeReceiveSize<CellID>(_loadBalancer.rank(iC),
                                                                             coordinator.get(iC, _iC, 1),
                                                                             _neighborhoodComm);
    _cellsOutboundTo[iC].resize(newOutboundCount);

    singleton::mpi().receive(
      _cellsOutboundTo[iC].data(),
      _cellsOutboundTo[iC].size(),
      _loadBalancer.rank(iC),
      coordinator.get(iC, _iC, 1),
      _neighborhoodComm);

    if (hasComponentSelectors()) {
      const auto activeFields = selectorActiveFields();
      std::vector<std::uint64_t> flat(activeFields.size() * newOutboundCount);
      singleton::mpi().receive(
        flat.data(), flat.size(),
        _loadBalancer.rank(iC),
        coordinator.get(iC, _iC, 2),
        _neighborhoodComm);
      auto& fieldMasks = _masksOutboundTo[iC];
      fieldMasks.clear();
      std::size_t off = 0;
      for (const auto& field : activeFields) {
        fieldMasks.emplace_back(field,
          std::vector<std::uint64_t>(flat.begin() + off, flat.begin() + off + newOutboundCount));
        off += newOutboundCount;
      }
    }
  });
}

template <typename T, unsigned D>
void BlockCommunicationNeighborhood<T,D>::wait()
{
  forNeighbors([&](int iC) {
    _fieldRequests[iC]->wait();
    _cellsRequests[iC]->wait();
    if (auto iter = _maskRequests.find(iC); iter != _maskRequests.end()) {
      std::get<1>(*iter)->wait();
    }
  });
}

#endif // PARALLEL_MODE_MPI

template <typename T, unsigned D>
const std::vector<std::type_index>&
BlockCommunicationNeighborhood<T,D>::getFieldsCommonWith(int iC) const
{
  return _fieldsCommonWith.at(iC);
}

template <typename T, unsigned D>
const std::vector<CellID>&
BlockCommunicationNeighborhood<T,D>::getCellsOutboundTo(int iC) const
{
  return _cellsOutboundTo.at(iC);
}

template <typename T, unsigned D>
const std::vector<CellID>&
BlockCommunicationNeighborhood<T,D>::getCellsInboundFrom(int iC) const
{
  return _cellsInboundFrom.at(iC);
}

template <typename T, unsigned D>
const std::vector<CellID>&
BlockCommunicationNeighborhood<T,D>::getCellsRequestedFrom(int iC) const
{
  return _cellsRequestedFrom.at(iC);
}

template <typename T, unsigned D>
const overlap::FieldMasks&
BlockCommunicationNeighborhood<T,D>::getMasksInboundFrom(int iC) const
{
  static const overlap::FieldMasks dense{};
  auto iter = _masksInboundFrom.find(iC);
  return iter != _masksInboundFrom.end() ? std::get<1>(*iter) : dense;
}

template <typename T, unsigned D>
const overlap::FieldMasks&
BlockCommunicationNeighborhood<T,D>::getMasksOutboundTo(int iC) const
{
  static const overlap::FieldMasks dense{};
  auto iter = _masksOutboundTo.find(iC);
  return iter != _masksOutboundTo.end() ? std::get<1>(*iter) : dense;
}

// *INDENT-ON*

}

#endif
