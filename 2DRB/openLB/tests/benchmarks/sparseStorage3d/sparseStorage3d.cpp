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

/* Single-GPU sparse-vs-dense storage microbenchmark.
 *
 * A trivial per-cell operator reads and writes a 3-component field through the
 * genuine device access path (Cell::getField / setField -> DeviceSparseColumn /
 * DeviceColumnAccessor), so the measured cost is the real storage access, not a
 * synthetic reimplementation. Four variants isolate the two questions:
 *
 *   V1 plain-dense        : a dense Column field, whole domain   (baseline, no branch)
 *   V2 sparse-in-dense    : sparse-capable field, Dense mode      (Q-A: the dir==nullptr branch tax)
 *   V3 sparse-full        : sparse-capable field, all tiles       (tile-indirection at full occupancy)
 *   V4 sparse-surface     : sparse-capable field, shell tiles    (Q-B: the payoff -- memory + access)
 *
 * Timing: per-apply GPU time via device-synchronised wall clock, warmup dropped,
 * median + p5/p95 over many iterations. Memory: analytical field footprint.
 * Pin CUDA_VISIBLE_DEVICES=1; no VTK/IO.
 */

#include <olb.h>
#include "sparseStorage3d.h"

#include <chrono>
#include <algorithm>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdint>
#include <unordered_set>
#include <functional>

using namespace olb;
using namespace olb::descriptors;

using T = float;
using DESCRIPTOR = BENCH_DESCRIPTOR;

// ---------------------------------------------------------------------------
// Stage 0: host-only tile-occupancy analysis. For a given allocated cell set
// (data pattern) it counts the distinct tiles under several tile shapes and
// reports the predicted payload memory -- exactly, with no device code. Tiles
// all hold TILE=64 cells; the shapes are factorisations of 64, so they differ
// only in how a 3D region maps onto tiles. The current scheme is 1x1x64 along
// the fastest (z) axis: i = x*Ny*Nz + y*Nz + z, tile = i / 64.
// ---------------------------------------------------------------------------

struct AnalysisTileShape { int tx, ty, tz; const char* name; };

template <typename Pred>
void analyzePattern(const std::string& pname, int N,
                    Pred inside, const std::vector<AnalysisTileShape>& shapes)
{
  constexpr unsigned TILE = 64;
  constexpr unsigned dim  = 3;
  const std::size_t ncells = std::size_t(N)*N*N;
  const std::uint64_t Ny = N, Nz = N;

  std::size_t cells = 0;
  std::unordered_set<std::uint64_t> linearTiles;
  std::vector<std::unordered_set<std::uint64_t>> shapeTiles(shapes.size());

  for (int x=0; x < N; ++x)
  for (int y=0; y < N; ++y)
  for (int z=0; z < N; ++z) {
    if (!inside(x,y,z)) continue;
    ++cells;
    const std::uint64_t i = (std::uint64_t(x)*Ny + y)*Nz + z;
    linearTiles.insert(i / TILE);
    for (std::size_t s=0; s < shapes.size(); ++s) {
      const auto& sh = shapes[s];
      const std::uint64_t gx = x/sh.tx, gy = y/sh.ty, gz = z/sh.tz;
      const std::uint64_t Gy = (N+sh.ty-1)/sh.ty, Gz = (N+sh.tz-1)/sh.tz;
      shapeTiles[s].insert((gx*Gy + gy)*Gz + gz);
    }
  }

  auto mb  = [&](std::size_t tiles){ return double(tiles)*TILE*dim*sizeof(T)/(1024.0*1024.0); };
  auto occ = [&](std::size_t tiles){ return tiles ? 100.0*cells/double(tiles*TILE) : 0.0; };
  const double denseMB = double(ncells)*dim*sizeof(T)/(1024.0*1024.0);

  std::printf("%-11s %5.1f%% dense=%6.1f | 1x1x64 %6.1fMB %3.0f%%",
              pname.c_str(), 100.0*cells/ncells, denseMB,
              mb(linearTiles.size()), occ(linearTiles.size()));
  for (std::size_t s=0; s < shapes.size(); ++s) {
    std::printf("  %s %6.1f %3.0f%%", shapes[s].name, mb(shapeTiles[s].size()), occ(shapeTiles[s].size()));
  }
  std::printf("\n");
}

void runAnalysis(int N)
{
  const std::vector<AnalysisTileShape> shapes = {
    {4,4,4, "4x4x4"}, {8,8,1, "8x8x1"}, {2,2,16,"2x2x16"}, {64,1,1,"64x1x1"},
  };
  std::printf("N=%d, TILE=64, columns are payload MB and occupancy%% (higher occ = tighter)\n", N);
  std::printf("pattern      cells  dense   scheme/shape ->\n");
  const int t1=1, t3=3, t8=8, s=N/4;
  analyzePattern("shell t=1", N, [=](int x,int y,int z){ return x<t1||x>=N-t1||y<t1||y>=N-t1||z<t1||z>=N-t1; }, shapes);
  analyzePattern("shell t=3", N, [=](int x,int y,int z){ return x<t3||x>=N-t3||y<t3||y>=N-t3||z<t3||z>=N-t3; }, shapes);
  analyzePattern("shell t=8", N, [=](int x,int y,int z){ return x<t8||x>=N-t8||y<t8||y>=N-t8||z<t8||z>=N-t8; }, shapes);
  analyzePattern("slab_z<8",  N, [=](int x,int y,int z){ return z<t8; }, shapes);
  analyzePattern("slab_x<8",  N, [=](int x,int y,int z){ return x<t8; }, shapes);
  analyzePattern("strip_z",   N, [=](int x,int y,int z){ return x<4 && y<4; }, shapes);
  analyzePattern("scatter/8", N, [=](int x,int y,int z){ return x%8==0 && y%8==0 && z%8==0; }, shapes);
  analyzePattern("blob N/4",  N, [=](int x,int y,int z){ return x>=s&&x<3*s&&y>=s&&y<3*s&&z>=s&&z<3*s; }, shapes);
}

struct Result {
  std::string name;
  std::size_t cells;
  double medUs, p5Us, p95Us;
  double megabytes;
};

static double percentile(std::vector<double>& v, double p) {
  return v[std::min(v.size()-1, std::size_t(p*v.size()))];
}

template <typename FIELD, template<typename> class OP = TouchO>
Result runVariant(const std::string& name, int N, int shell,
                  FieldStorage mode, bool allocateWholeDomain, bool applyShellOnly,
                  int warmup, int iters, olb::sparse::TileShape shape = {1,1,64})
{
  const T dx = T{1};
  CuboidDecomposition3D<T> cuboidDecomposition(Vector<T,3>(0), dx, Vector<int,3>(N), 1);
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);
  for (int iC=0; iC < loadBalancer.size(); ++iC) {
    loadBalancer.setPlatform(iC, Platform::GPU_CUDA);
  }
  SuperGeometry<T,3> superGeometry(cuboidDecomposition, loadBalancer);
  // shell (material 2) around an interior (material 1)
  superGeometry.rename(0, 2);
  IndicatorCuboid3D<T> interior(Vector<T,3>((N-2*shell)*dx), Vector<T,3>(shell*dx));
  superGeometry.rename(2, 1, interior);
  superGeometry.communicate();

  UnitConverter<T,DESCRIPTOR> converter(T{1},T{1},T{1},T{1},T{1},T{1},T{0});
  SuperLattice<T,DESCRIPTOR> superLattice(converter, superGeometry);
  dynamics::set<NoDynamics>(superLattice, superGeometry.getMaterialIndicator({1,2}));

  constexpr bool isSparse = std::is_same_v<FIELD,SPARSEVEL>;
  if constexpr (isSparse) {
    superLattice.template setFieldStorage<FIELD>(mode);
    if (mode == FieldStorage::Sparse) {
      superLattice.template setTileShape<FIELD>(shape);
      superLattice.template allocateField<FIELD>(
        allocateWholeDomain ? superGeometry.getMaterialIndicator({1,2})
                            : superGeometry.getMaterialIndicator(2));
    }
  }

  auto region = applyShellOnly ? superGeometry.getMaterialIndicator(2)
                               : superGeometry.getMaterialIndicator({1,2});
  const std::size_t cells = applyShellOnly
    ? superGeometry.getStatistics().getNvoxel(2)
    : superGeometry.getStatistics().getNvoxel(1) + superGeometry.getStatistics().getNvoxel(2);

  superLattice.template addPostProcessor<TouchStage>(std::move(region), meta::id<OP<FIELD>>{});
  superLattice.initialize();
  superLattice.setProcessingContext(ProcessingContext::Simulation);

  for (int i=0; i < warmup; ++i) {
    superLattice.template executePostProcessors<TouchStage>();
  }
  gpu::cuda::device::synchronize();

  std::vector<double> us;
  us.reserve(iters);
  for (int i=0; i < iters; ++i) {
    gpu::cuda::device::synchronize();
    const auto t0 = std::chrono::high_resolution_clock::now();
    superLattice.template executePostProcessors<TouchStage>();
    gpu::cuda::device::synchronize();
    const auto t1 = std::chrono::high_resolution_clock::now();
    us.push_back(std::chrono::duration<double,std::micro>(t1-t0).count());
  }
  std::sort(us.begin(), us.end());

  // Analytical field footprint
  const std::size_t ncells = std::size_t(N)*N*N;
  constexpr unsigned dim = 3;
  double bytes = 0;
  if constexpr (isSparse) {
    if (mode == FieldStorage::Dense) {
      bytes = double(ncells) * dim * sizeof(T);
    } else {
      constexpr unsigned TILE = 64;
      const std::size_t tileGrid = (ncells + TILE - 1) / TILE;
      const std::size_t tiles = superLattice.getBlock(0).template getField<FIELD>().getAllocatedTiles();
      // per dimension: (tiles + scratch) payload + directory
      bytes = dim * (double(tiles+1)*TILE*sizeof(T) + double(tileGrid)*sizeof(std::int32_t));
    }
  } else {
    bytes = double(ncells) * dim * sizeof(T);
  }

  return Result{name, cells,
                us[us.size()/2], percentile(us,0.05), percentile(us,0.95),
                bytes/(1024.0*1024.0)};
}

int main(int argc, char* argv[])
{
  initialize(&argc, &argv);

  if (argc > 1 && std::string(argv[1]) == "analyze") {
    runAnalysis((argc > 2) ? std::atoi(argv[2]) : 160);
    return 0;
  }

  int N      = (argc > 1) ? std::atoi(argv[1]) : 160;
  int shell  = (argc > 2) ? std::atoi(argv[2]) : 3;
  int warmup = 50;
  int iters  = 500;

  OstreamManager clout(std::cout, "bench");
  clout << "N=" << N << " (" << std::size_t(N)*N*N << " cells), shell=" << shell
        << ", warmup=" << warmup << ", iters=" << iters << std::endl;

  std::vector<Result> results;
  results.push_back(runVariant<DENSEVEL> ("V1 plain-dense        ", N, shell, FieldStorage::Dense,  false, false, warmup, iters));
  results.push_back(runVariant<SPARSEVEL>("V2 sparse-in-dense    ", N, shell, FieldStorage::Dense,  false, false, warmup, iters));
  results.push_back(runVariant<SPARSEVEL>("V3 sparse-full        ", N, shell, FieldStorage::Sparse, true,  false, warmup, iters));
  results.push_back(runVariant<SPARSEVEL>("V4 sparse-surface     ", N, shell, FieldStorage::Sparse, false, true,  warmup, iters));
  results.push_back(runVariant<SPARSEVEL>("V4b surface 4x4x4      ", N, shell, FieldStorage::Sparse, false, true,  warmup, iters, olb::sparse::TileShape{4,4,4}));
  results.push_back(runVariant<DENSEVEL,TouchDecomposeO>("V5 dense+3D-decompose ", N, shell, FieldStorage::Dense, false, false, warmup, iters));

  clout << "variant                  cells        med(us)   p5(us)   p95(us)   ns/cell   MB" << std::endl;
  for (const auto& r : results) {
    char line[256];
    std::snprintf(line, sizeof(line), "%-23s %10zu  %8.1f %8.1f %8.1f  %8.3f %8.2f",
                  r.name.c_str(), r.cells, r.medUs, r.p5Us, r.p95Us,
                  1000.0*r.medUs/double(r.cells), r.megabytes);
    clout << line << std::endl;
  }
  return 0;
}
