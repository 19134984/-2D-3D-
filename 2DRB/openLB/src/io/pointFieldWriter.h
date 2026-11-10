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

#ifndef IO_POINT_FIELD_WRITER_H
#define IO_POINT_FIELD_WRITER_H

#include <functional>
#include <string>
#include <vector>
#include <fstream>

#include "core/vector.h"
#include "core/concepts.h"
#include "utilities/samplingGenerators.h"
#include "communication/mpiManager.h"
#include "sampleWriterNativeVtk.h"
#include "core/olbInit.h"

namespace olb {

namespace io {

/// ImageData grid matching an axis-aligned plane lattice (Nx along U, Ny along V, 1 along the normal)
template <typename T>
RegularGrid<T> regularGrid(const HyperplaneLattice3D<T>& lattice)
{
  auto axisOf = [](Vector<T,3> d) {
    unsigned a = 0;
    for (unsigned k=1; k<3; ++k) {
      if ((d[k] < T{0} ? -d[k] : d[k]) > (d[a] < T{0} ? -d[a] : d[a])) { a = k; }
    }
    return a;
  };
  RegularGrid<T> grid;
  grid.origin  = lattice.getPhysOrigin();
  grid.spacing = {lattice.getPhysSpacing(), lattice.getPhysSpacing(), lattice.getPhysSpacing()};
  grid.extent  = {1, 1, 1};
  grid.extent[axisOf(lattice.getVectorU())] = lattice.getNx();
  grid.extent[axisOf(lattice.getVectorV())] = lattice.getNy();
  return grid;
}

/// Named data arrays addressed by (local point index, component), bound to a PointSampleSource
template <typename T>
struct PointFields {
  struct Field {
    std::string name;
    unsigned    nComp;
    std::function<T(std::size_t,unsigned)> at;
  };
  std::vector<Field> list;

  void add(std::string name, unsigned nComp, std::function<T(std::size_t,unsigned)> at) {
    list.push_back({std::move(name), nComp, std::move(at)});
  }
};

/// PointFields derived from a sampler's functor result_t: one array per output field, named by
/// the field type, reading the matching component range of getValue.
template <typename SAMPLER>
PointFields<typename SAMPLER::value_t> pointFieldsFor(SAMPLER& sampler) {
  using T          = typename SAMPLER::value_t;
  using DESCRIPTOR = typename SAMPLER::descriptor_t;
  PointFields<T> fields;
  auto* s = &sampler;
  unsigned offset = 0;
  descriptors::field_list_of<typename SAMPLER::functor_t::result_t>::for_each([&](auto id) {
    using FIELD = typename decltype(id)::type;
    const std::string_view qn = meta::nice_name<FIELD>();
    const auto sep = qn.rfind("::");
    std::string name{sep == std::string_view::npos ? qn : qn.substr(sep + 2)};
    const unsigned base = offset;
    const unsigned n    = FIELD::template size<DESCRIPTOR>();
    fields.add(std::move(name), n,
               [s,base](std::size_t i, unsigned c){ return s->getValue(i)[base + c]; });
    offset += n;
  });
  return fields;
}

/// Owned point samples collected on the root rank; itself a PointSampleSource
template <typename T, unsigned D>
struct GatheredSamples {
  static constexpr unsigned d = D;

  std::size_t                             n = 0;
  std::size_t                             nGlobal = 0;
  std::vector<T>                          coords;
  std::vector<std::size_t>                gIndex;
  std::vector<std::pair<std::string,unsigned>> fieldMeta;
  std::vector<std::vector<T>>             data;

  std::size_t getLocalSize()  const { return n; }
  std::size_t getGlobalSize() const { return nGlobal; }
  std::size_t getGlobalIndex(std::size_t i) const { return gIndex[i]; }

  Vector<T,D> getPhysR(std::size_t i) const {
    Vector<T,D> p;
    for (unsigned k=0; k<D; ++k) {
      p[k] = coords[i*D + k];
    }
    return p;
  }
};

namespace detail {

template <concepts::PointSampleSource SOURCE, typename T>
GatheredSamples<T,std::decay_t<SOURCE>::d>
gatherSamples(const SOURCE& source, const PointFields<T>& fields, int root=0)
{
  constexpr unsigned D = std::decay_t<SOURCE>::d;
  const std::size_t nLocal = source.getLocalSize();

  std::vector<T>           coordLocal(nLocal*D);
  std::vector<std::size_t> gIdxLocal(nLocal);
  for (std::size_t i=0; i<nLocal; ++i) {
    auto p = source.getPhysR(i);
    for (unsigned k=0; k<D; ++k) {
      coordLocal[i*D + k] = p[k];
    }
    gIdxLocal[i] = source.getGlobalIndex(i);
  }
  std::vector<std::vector<T>> fieldLocal(fields.list.size());
  for (std::size_t f=0; f<fields.list.size(); ++f) {
    const auto& fd = fields.list[f];
    fieldLocal[f].resize(nLocal*fd.nComp);
    for (std::size_t i=0; i<nLocal; ++i) {
      for (unsigned c=0; c<fd.nComp; ++c) {
        fieldLocal[f][i*fd.nComp + c] = fd.at(i,c);
      }
    }
  }

  GatheredSamples<T,D> g;
  g.nGlobal = source.getGlobalSize();
  for (const auto& fd : fields.list) {
    g.fieldMeta.emplace_back(fd.name, fd.nComp);
  }
  g.data.resize(fields.list.size());

#ifdef PARALLEL_MODE_MPI
  const int rank = singleton::mpi().getRank();
  const int size = singleton::mpi().getSize();

  int localCount = static_cast<int>(nLocal);
  std::vector<int> counts(size, 0);
  singleton::mpi().gather(&localCount, 1, counts.data(), 1, root);

  std::size_t total = 0;
  std::vector<int> displs(size, 0);
  if (rank == root) {
    for (int r=0; r<size; ++r) {
      displs[r] = static_cast<int>(total);
      total += counts[r];
    }
  }
  g.n = (rank == root) ? total : 0;

  auto scaledGatherv = [&](auto* sendBuf, auto* recvBuf, unsigned stride) {
    std::vector<int> c(size), d(size);
    for (int r=0; r<size; ++r) {
      c[r] = counts[r]*stride;
      d[r] = displs[r]*stride;
    }
    singleton::mpi().gatherv(sendBuf, localCount*static_cast<int>(stride),
                             recvBuf, c.data(), d.data(), root);
  };

  g.coords.resize(g.n*D);
  scaledGatherv(coordLocal.data(), g.coords.data(), D);
  g.gIndex.resize(g.n);
  scaledGatherv(gIdxLocal.data(), g.gIndex.data(), 1);
  for (std::size_t f=0; f<fields.list.size(); ++f) {
    g.data[f].resize(g.n*fields.list[f].nComp);
    scaledGatherv(fieldLocal[f].data(), g.data[f].data(), fields.list[f].nComp);
  }
#else
  g.n      = nLocal;
  g.coords = std::move(coordLocal);
  g.gIndex = std::move(gIdxLocal);
  g.data   = std::move(fieldLocal);
#endif

  return g;
}

template <concepts::PointSampleSource SOURCE, typename T>
void writeVtuPiece(std::ostream& out, const SOURCE& source, const PointFields<T>& fields)
{
  constexpr unsigned D = std::decay_t<SOURCE>::d;
  const std::size_t n = source.getLocalSize();

  out << "    <Piece NumberOfPoints=\"" << n << "\" NumberOfCells=\"" << n << "\">\n"
      << "      <Points>\n"
      << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for (std::size_t i=0; i<n; ++i) {
    auto p = source.getPhysR(i);
    for (unsigned k=0; k<3; ++k) {
      out << ' ' << (k < D ? p[k] : T{0});
    }
  }
  out << "\n        </DataArray>\n      </Points>\n      <Cells>\n"
      << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
  for (std::size_t i=0; i<n; ++i) { out << ' ' << i; }
  out << "\n        </DataArray>\n"
      << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
  for (std::size_t i=0; i<n; ++i) { out << ' ' << (i+1); }
  out << "\n        </DataArray>\n"
      << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (std::size_t i=0; i<n; ++i) { out << " 1"; }
  out << "\n        </DataArray>\n      </Cells>\n      <PointData>\n";
  for (const auto& f : fields.list) {
    out << "        <DataArray type=\"Float32\" Name=\"" << f.name
        << "\" NumberOfComponents=\"" << f.nComp << "\" format=\"ascii\">\n";
    for (std::size_t i=0; i<n; ++i) {
      for (unsigned c=0; c<f.nComp; ++c) {
        out << ' ' << f.at(i,c);
      }
    }
    out << "\n        </DataArray>\n";
  }
  out << "      </PointData>\n    </Piece>\n";
}

}  // namespace detail

/// Gather a source's result values into a dense rank-0 vector indexed by global index (empty elsewhere)
template <concepts::PointSampleSource SOURCE>
auto gatherValues(const SOURCE& source, std::size_t nGlobal)
{
  using Value = std::decay_t<decltype(source.getValue(std::size_t{0}))>;
  using T     = std::decay_t<decltype(source.getValue(std::size_t{0})[0])>;
  constexpr unsigned R = Value::d;

  PointFields<T> f;
  f.add("value", R, [&source](std::size_t i, unsigned c){ return source.getValue(i)[c]; });
  auto g = detail::gatherSamples(source, f);

  std::vector<Value> full(singleton::mpi().getRank() == 0 ? nGlobal : 0);
  for (std::size_t i = 0; i < g.n; ++i) {
    Value value;
    for (unsigned c = 0; c < R; ++c) {
      value[c] = g.data[0][i*R + c];
    }
    full[g.gIndex[i]] = value;
  }
  return full;
}

/// Distributed unstructured write: one .vtu piece per rank plus a root .pvtu index
template <concepts::PointSampleSource SOURCE, typename T>
void writeVtu(const SOURCE& source, const PointFields<T>& fields,
              const std::string& path, std::size_t iT)
{
  const int rank   = singleton::mpi().getRank();
  const int nRanks = singleton::mpi().getSize();
  const std::string base = path + "_" + std::to_string(iT);
  const std::size_t slash = base.find_last_of('/');
  const std::string stem = (slash == std::string::npos) ? base : base.substr(slash+1);

  std::ofstream out(base + "_p" + std::to_string(rank) + ".vtu");
  out << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
      << "  <UnstructuredGrid>\n";
  detail::writeVtuPiece(out, source, fields);
  out << "  </UnstructuredGrid>\n</VTKFile>\n";
  out.close();

  if (rank == 0) {
    std::ofstream p(base + ".pvtu");
    p << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
      << "  <PUnstructuredGrid GhostLevel=\"0\">\n"
      << "    <PPoints>\n      <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n    </PPoints>\n"
      << "    <PPointData>\n";
    for (const auto& f : fields.list) {
      p << "      <PDataArray type=\"Float32\" Name=\"" << f.name
        << "\" NumberOfComponents=\"" << f.nComp << "\"/>\n";
    }
    p << "    </PPointData>\n";
    for (int iRank=0; iRank<nRanks; ++iRank) {
      p << "    <Piece Source=\"" << stem << "_p" << iRank << ".vtu\"/>\n";
    }
    p << "  </PUnstructuredGrid>\n</VTKFile>\n";
  }
}

/// Reduced unstructured write: gather to root and emit a single consolidated .vtu
template <concepts::PointSampleSource SOURCE, typename T>
void writeVtuReduced(const SOURCE& source, const PointFields<T>& fields,
                     const std::string& path, std::size_t iT)
{
  auto g = detail::gatherSamples(source, fields);
  if (singleton::mpi().getRank() != 0) {
    return;
  }
  PointFields<T> gf;
  for (std::size_t f=0; f<g.fieldMeta.size(); ++f) {
    const unsigned nComp = g.fieldMeta[f].second;
    gf.add(g.fieldMeta[f].first, nComp,
           [&g,f,nComp](std::size_t i, unsigned c){ return g.data[f][i*nComp + c]; });
  }
  std::ofstream out(path + "_" + std::to_string(iT) + ".vtu");
  out << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
      << "  <UnstructuredGrid>\n";
  detail::writeVtuPiece(out, g, gf);
  out << "  </UnstructuredGrid>\n</VTKFile>\n";
}

namespace detail {

/// Write an already-gathered sample set into a structured .vti (root-only, self-contained)
template <typename T>
void writeVtiFromGathered(const GatheredSamples<T,3>& g, const RegularGrid<T>& grid,
                          const std::string& path, std::size_t iT)
{
  const auto nNode = grid.nodeCount();
  const auto ext = grid.extent;

  std::ofstream out(path + "_" + std::to_string(iT) + ".vti");
  out << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
      << "  <ImageData WholeExtent=\"0 " << (ext[0]-1) << " 0 " << (ext[1]-1) << " 0 " << (ext[2]-1) << "\" "
      << "Origin=\"" << grid.origin[0] << " " << grid.origin[1] << " " << grid.origin[2] << "\" "
      << "Spacing=\"" << grid.spacing[0] << " " << grid.spacing[1] << " " << grid.spacing[2] << "\">\n"
      << "    <Piece Extent=\"0 " << (ext[0]-1) << " 0 " << (ext[1]-1) << " 0 " << (ext[2]-1) << "\">\n"
      << "      <PointData>\n";
  std::vector<T> dense;
  for (std::size_t f=0; f<g.fieldMeta.size(); ++f) {
    const unsigned nComp = g.fieldMeta[f].second;
    dense.assign(nNode*nComp, T{0});
    for (std::size_t i=0; i<g.n; ++i) {
      const std::size_t node = g.gIndex[i];
      for (unsigned c=0; c<nComp; ++c) {
        dense[node*nComp + c] = g.data[f][i*nComp + c];
      }
    }
    out << "        <DataArray type=\"Float32\" Name=\"" << g.fieldMeta[f].first
        << "\" NumberOfComponents=\"" << nComp << "\" format=\"ascii\">\n";
    for (std::size_t k=0; k<dense.size(); ++k) {
      out << ' ' << dense[k];
    }
    out << "\n        </DataArray>\n";
  }
  out << "      </PointData>\n    </Piece>\n  </ImageData>\n</VTKFile>\n";
}

/// Write an already-gathered sample set into a single reduced .vtu (+ .pvtu index), root-only
template <typename T>
void writeVtuFromGathered(const GatheredSamples<T,3>& g, const std::string& path, std::size_t iT)
{
  PointFields<T> gf;
  for (std::size_t f=0; f<g.fieldMeta.size(); ++f) {
    const unsigned nComp = g.fieldMeta[f].second;
    gf.add(g.fieldMeta[f].first, nComp,
           [&g,f,nComp](std::size_t i, unsigned c){ return g.data[f][i*nComp + c]; });
  }
  const std::string base = path + "_" + std::to_string(iT);
  const std::size_t slash = base.find_last_of('/');
  const std::string stem = (slash == std::string::npos) ? base : base.substr(slash+1);

  std::ofstream out(base + "_p0.vtu");
  out << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
      << "  <UnstructuredGrid>\n";
  writeVtuPiece(out, g, gf);
  out << "  </UnstructuredGrid>\n</VTKFile>\n";
  out.close();

  std::ofstream p(base + ".pvtu");
  p << "<?xml version=\"1.0\"?>\n"
    << "<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n"
    << "  <PUnstructuredGrid GhostLevel=\"0\">\n"
    << "    <PPoints>\n      <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n    </PPoints>\n"
    << "    <PPointData>\n";
  for (const auto& f : gf.list) {
    p << "      <PDataArray type=\"Float32\" Name=\"" << f.name
      << "\" NumberOfComponents=\"" << f.nComp << "\"/>\n";
  }
  p << "    </PPointData>\n    <Piece Source=\"" << stem << "_p0.vtu\"/>\n"
    << "  </PUnstructuredGrid>\n</VTKFile>\n";
}

}  // namespace detail

/// Reduced structured write: gather to root and write the .vti
template <concepts::PointSampleSource SOURCE, typename T>
void writeVti(const SOURCE& source, const PointFields<T>& fields,
              const RegularGrid<T>& grid, const std::string& path, std::size_t iT)
{
  auto g = detail::gatherSamples(source, fields);
  if (singleton::mpi().getRank() == 0) {
    detail::writeVtiFromGathered(g, grid, path, iT);
  }
}

/// As writeVti, but the file write runs on the background ThreadPool. The sampled data is snapshot
/// before scheduling (so a later update() cannot race it), backpressure caps pending writes to bound
/// memory, and it falls back to a synchronous write when the pool has no workers.
template <concepts::PointSampleSource SOURCE, typename T>
void writeVtiAsync(const SOURCE& source, const PointFields<T>& fields,
                   const RegularGrid<T>& grid, const std::string& path, std::size_t iT)
{
  auto g = detail::gatherSamples(source, fields);
  if (singleton::mpi().getRank() != 0) {
    return;
  }
  if (singleton::pool().size() == 0) {
    detail::writeVtiFromGathered(g, grid, path, iT);
    return;
  }
  singleton::pool().waitUntilBelow(16);
  singleton::pool().scheduleAndForget(
    [g=std::move(g), grid, path, iT]() { detail::writeVtiFromGathered(g, grid, path, iT); });
}

/// As writeVtuReduced, but the file write runs on the background ThreadPool (see writeVtiAsync)
template <concepts::PointSampleSource SOURCE, typename T>
void writeVtuAsync(const SOURCE& source, const PointFields<T>& fields,
                   const std::string& path, std::size_t iT)
{
  auto g = detail::gatherSamples(source, fields);
  if (singleton::mpi().getRank() != 0) {
    return;
  }
  if (singleton::pool().size() == 0) {
    detail::writeVtuFromGathered(g, path, iT);
    return;
  }
  singleton::pool().waitUntilBelow(16);
  singleton::pool().scheduleAndForget(
    [g=std::move(g), path, iT]() { detail::writeVtuFromGathered(g, path, iT); });
}

}  // namespace io

}  // namespace olb

#endif
