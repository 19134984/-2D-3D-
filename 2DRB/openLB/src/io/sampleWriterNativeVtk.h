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

#ifndef IO_SAMPLE_WRITER_NATIVE_VTK_H
#define IO_SAMPLE_WRITER_NATIVE_VTK_H

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "core/sampledData.h"
#include "core/vector.h"

namespace olb {

namespace io {

/// Axis-aligned physical grid describing a VTK ImageData block (origin, per-axis spacing, extent)
template <typename T>
struct RegularGrid {
  Vector<T,3>   origin{0,0,0};
  Vector<T,3>   spacing{1,1,1};
  Vector<int,3> extent{1,1,1};

  std::size_t nodeCount() const {
    return static_cast<std::size_t>(extent[0]) * extent[1] * extent[2];
  }
};

template <typename T>
void writeVti(const SampledData<T>& data, const RegularGrid<T>& grid,
              const std::string& path, std::size_t iT)
{
  const auto n = grid.extent;
  std::ofstream out(path + "_" + std::to_string(iT) + ".vti");
  out << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
      << "  <ImageData WholeExtent=\"0 " << (n[0]-1) << " 0 " << (n[1]-1) << " 0 " << (n[2]-1) << "\" "
      << "Origin=\"" << grid.origin[0] << " " << grid.origin[1] << " " << grid.origin[2] << "\" "
      << "Spacing=\"" << grid.spacing[0] << " " << grid.spacing[1] << " " << grid.spacing[2] << "\">\n"
      << "    <Piece Extent=\"0 " << (n[0]-1) << " 0 " << (n[1]-1) << " 0 " << (n[2]-1) << "\">\n"
      << "      <PointData>\n";
  for (const auto& f : data.fields) {
    out << "        <DataArray type=\"Float32\" Name=\"" << f.name
        << "\" NumberOfComponents=\"" << f.nComp << "\" format=\"ascii\">\n";
    for (std::size_t i=0; i < data.nLocal; ++i) {
      for (unsigned iComp=0; iComp < f.nComp; ++iComp) {
        out << ' ' << f.components[iComp][i];
      }
    }
    out << "\n        </DataArray>\n";
  }
  out << "      </PointData>\n"
      << "    </Piece>\n"
      << "  </ImageData>\n"
      << "</VTKFile>\n";
}

template <typename T>
void writeVthb(const std::string& path, std::size_t iT,
               const std::vector<RegularGrid<T>>& levelGrids,
               const std::vector<std::string>& levelVtiFiles)
{
  const Vector<T,3> origin0 = levelGrids.empty() ? Vector<T,3>{0,0,0} : levelGrids[0].origin;
  const T delta0 = levelGrids.empty() ? T{1} : levelGrids[0].spacing[0];

  std::ofstream out(path + "_" + std::to_string(iT) + ".vthb");
  out << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"vtkOverlappingAMR\" version=\"1.1\" byte_order=\"LittleEndian\">\n"
      << "  <vtkOverlappingAMR origin=\""
      << (origin0[0]-0.5*delta0) << " " << (origin0[1]-0.5*delta0) << " " << (origin0[2]-0.5*delta0)
      << "\" grid_description=\"XYZ\">\n";
  for (std::size_t iLevel=0; iLevel < levelGrids.size(); ++iLevel) {
    const auto& g = levelGrids[iLevel];
    const T delta_l = g.spacing[0];
    int refRatio = (iLevel == 0) ? 1 : 2;
    out << "    <Block level=\"" << iLevel << "\" refinement_ratio=\"" << refRatio << "\" spacing=\""
        << delta_l << " " << delta_l << " " << delta_l << "\">\n";
    int imin[3], imax[3];
    for (unsigned iD=0; iD<3; ++iD) {
      imin[iD] = static_cast<int>(std::round((g.origin[iD]-origin0[iD])/delta_l + 0.5*(delta0/delta_l - 1.0)));
      imax[iD] = imin[iD] + g.extent[iD] - 1;
    }
    out << "      <DataSet index=\"0\" amr_box=\""
        << imin[0] << " " << imax[0] << " " << imin[1] << " " << imax[1] << " " << imin[2] << " " << imax[2]
        << "\" file=\"" << levelVtiFiles[iLevel] << "\" />\n"
        << "    </Block>\n";
  }
  out << "  </vtkOverlappingAMR>\n"
      << "</VTKFile>\n";
}

template <typename T>
void writeVtu(const SampledData<T>& data, const std::string& path, std::size_t iT)
{
  const std::string base = path + "_" + std::to_string(iT);
  const std::size_t slash = base.find_last_of('/');
  const std::string stem = (slash == std::string::npos) ? base : base.substr(slash+1);
  const std::size_t n = data.nLocal;

  std::ofstream out(base + "_p" + std::to_string(data.rank) + ".vtu");
  out << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
      << "  <UnstructuredGrid>\n"
      << "    <Piece NumberOfPoints=\"" << n << "\" NumberOfCells=\"" << n << "\">\n"
      << "      <Points>\n"
      << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (std::size_t i=0; i<n; ++i) {
    for (unsigned iD=0; iD<3; ++iD) {
      out << ' ' << (iD < data.pointCoords.size() ? data.pointCoords[iD][i] : T{0});
    }
  }
  out << "\n        </DataArray>\n      </Points>\n      <Cells>\n"
      << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
  for (std::size_t i=0; i<n; ++i) {
    out << ' ' << i;
  }
  out << "\n        </DataArray>\n"
      << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
  for (std::size_t i=0; i<n; ++i) {
    out << ' ' << (i+1);
  }
  out << "\n        </DataArray>\n"
      << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (std::size_t i=0; i<n; ++i) {
    out << " 1";
  }
  out << "\n        </DataArray>\n      </Cells>\n      <PointData>\n";
  for (const auto& f : data.fields) {
    out << "        <DataArray type=\"Float32\" Name=\"" << f.name
        << "\" NumberOfComponents=\"" << f.nComp << "\" format=\"ascii\">\n";
    for (std::size_t i=0; i<n; ++i) {
      for (unsigned iComp=0; iComp<f.nComp; ++iComp) {
        out << ' ' << f.components[iComp][i];
      }
    }
    out << "\n        </DataArray>\n";
  }
  out << "      </PointData>\n    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
  out.close();

  if (data.rank == 0) {
    std::ofstream p(base + ".pvtu");
    p << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
      << "  <PUnstructuredGrid GhostLevel=\"0\">\n"
      << "    <PPoints>\n      <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n    </PPoints>\n"
      << "    <PPointData>\n";
    for (const auto& f : data.fields) {
      p << "      <PDataArray type=\"Float32\" Name=\"" << f.name
        << "\" NumberOfComponents=\"" << f.nComp << "\"/>\n";
    }
    p << "    </PPointData>\n";
    for (int iRank=0; iRank<data.nRanks; ++iRank) {
      p << "    <Piece Source=\"" << stem << "_p" << iRank << ".vtu\"/>\n";
    }
    p << "  </PUnstructuredGrid>\n</VTKFile>\n";
  }
}

}  // namespace io

}

#endif
