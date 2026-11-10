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

#ifndef CORE_DATA_HH
#define CORE_DATA_HH

#include "data.h"

#include "platform/deviceMemory.h"
#include "io/textTable.h"

#include <ostream>
#include <sstream>
#include <string>

namespace olb {

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteData<T,DESCRIPTOR,PLATFORM>::writeFieldStorageCSV(std::ostream& out, std::size_t coreCells) const
{
  out << "name; type; components; bytesPerCell; storage; coreCells; allocatedCells; tile; "
         "tileCapacity; allocatedTiles; directorySlots; setCells; fillRatio; payloadBytes; dirBytes; bytes"
      << std::endl;
  forEach([&](const AnyFieldTypeD<T,DESCRIPTOR,PLATFORM>& any) {
    const auto casted = any.template tryAs<AbstractFieldArrayBase>();
    if (!casted) {
      return;
    }
    const FieldStorageInfo info = (*casted)->getStorageInfo();
    const std::size_t bytesPerCell = std::size_t(info.components) * info.componentBytes;
    out << any.getName() << "; " << info.valueType << "; " << info.components
        << "; " << bytesPerCell << "; ";
    if (info.storage == FieldStorage::Sparse) {
      const std::size_t setCells = info.allocatedTiles * std::size_t(info.tileCapacity);
      const std::size_t payloadBytes = (info.allocatedTiles + 1) * std::size_t(info.tileCapacity) * bytesPerCell;
      const std::size_t dirBytes = info.components * info.directorySlots * sizeof(std::int32_t);
      const double fillRatio = coreCells ? static_cast<double>(setCells) / coreCells : 0.0;
      out << "Sparse; " << coreCells << "; " << info.allocatedCells << "; "
          << info.tileShape.tx << "x" << info.tileShape.ty << "x" << info.tileShape.tz
          << "; " << info.tileCapacity << "; " << info.allocatedTiles << "; " << info.directorySlots
          << "; " << setCells << "; " << fillRatio
          << "; " << payloadBytes << "; " << dirBytes << "; " << (payloadBytes + dirBytes) << std::endl;
    } else {
      const std::size_t bytes = info.allocatedCells * bytesPerCell;
      out << "Dense; " << coreCells << "; " << info.allocatedCells
          << "; -1; -1; -1; -1; -1; -1; " << bytes << "; 0; " << bytes << std::endl;
    }
  });
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
void ConcreteData<T,DESCRIPTOR,PLATFORM>::describeFieldStorage(std::ostream& out, std::size_t coreCells) const
{
  using Align = io::text::Table::Align;
  std::size_t actualTotal = 0;
  std::size_t denseTotal  = 0;
  const std::size_t deviceBytes = getDeviceTotalMemory<PLATFORM>();

  io::text::Table table({
    {"Field",   Align::Left,  48},
    {"Type",    Align::Left,  24},
    {"Storage", Align::Left,   0},
    {"Cells",   Align::Left,   0},
    {"Size",    Align::Right,  0},
    {"Device",  Align::Right,  0},
  });

  forEach([&](const AnyFieldTypeD<T,DESCRIPTOR,PLATFORM>& any) {
    const auto casted = any.template tryAs<AbstractFieldArrayBase>();
    if (!casted) {
      return;
    }
    const FieldStorageInfo info = (*casted)->getStorageInfo();
    const std::size_t bytesPerCell = std::size_t(info.components) * info.componentBytes;
    const std::size_t denseBytes = info.allocatedCells * bytesPerCell;
    denseTotal += denseBytes;

    const std::string type = std::string(info.valueType) + "[" + std::to_string(info.components) + "]";
    std::string storage;
    std::string cells;
    std::size_t bytes = 0;
    if (info.storage == FieldStorage::Sparse) {
      const std::size_t setCells = info.allocatedTiles * std::size_t(info.tileCapacity);
      const std::size_t payloadBytes = (info.allocatedTiles + 1) * std::size_t(info.tileCapacity) * bytesPerCell;
      const std::size_t dirBytes = info.components * info.directorySlots * sizeof(std::int32_t);
      bytes = payloadBytes + dirBytes;
      std::ostringstream s;
      s << "Sparse " << info.tileShape.tx << "x" << info.tileShape.ty << "x" << info.tileShape.tz
        << " cap" << info.tileCapacity << " (" << info.allocatedTiles << "t)";
      storage = s.str();
      cells = std::to_string(setCells) + "/" + std::to_string(coreCells);
    } else {
      bytes = denseBytes;
      storage = "Dense";
      cells = std::to_string(info.allocatedCells);
    }
    actualTotal += bytes;

    const std::string device = deviceBytes
      ? io::text::percent(static_cast<double>(bytes) / deviceBytes, 1)
      : std::string("n/a");
    table.row(std::string(any.getName()), type, storage, cells, io::text::humanBytes(bytes), device);
  });
  table.rule();

  out << std::endl;
  out << "---Field Storage---" << std::endl;
  table.render(out);
  const double ratio = denseTotal ? static_cast<double>(actualTotal) / denseTotal : 1.0;
  out << "  total: " << io::text::humanBytes(actualTotal) << " actual, "
      << io::text::humanBytes(denseTotal) << " dense-equivalent (" << io::text::percent(ratio) << ")"
      << std::endl;
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
std::size_t ConcreteData<T,DESCRIPTOR,PLATFORM>::getNblock() const
{
  std::size_t nBlock = 0;
  for (const auto& [_, serializable] : _serializationNames) {
    nBlock += serializable->getNblock();
  }
  return nBlock;
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
std::size_t ConcreteData<T,DESCRIPTOR,PLATFORM>::getSerializableSize() const
{
  std::size_t size = 0;
  for (const auto& [_, serializable] : _serializationNames) {
    size += serializable->getSerializableSize();
  }
  return size;
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
bool* ConcreteData<T,DESCRIPTOR,PLATFORM>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  for (const auto& [_, serializable] : _serializationNames) {
    this->registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, *serializable, loadingMode);
  }

  return dataPtr;
}

}

#endif
