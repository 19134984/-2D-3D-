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

#ifndef TESTS_COMMUNICATION_HARNESS_H
#define TESTS_COMMUNICATION_HARNESS_H

#include <olb.h>

#include <cmath>
#include <cstring>
#include <functional>
#include <iostream>
#include <map>
#include <sstream>

namespace test {

struct Context {
  int rank = 0;
  int ranks = 1;
  long failures = 0;
  long checks = 0;
  long reported = 0;

  static Context& get() {
    static Context ctx;
    return ctx;
  }
};

inline void fail(const std::string& message)
{
  auto& ctx = Context::get();
  ctx.failures++;
  if (ctx.reported++ < 32) {
    std::cout << "[rank " << ctx.rank << "] FAIL " << message << std::endl;
  }
}

inline void check(bool condition, const std::string& message)
{
  Context::get().checks++;
  if (!condition) {
    fail(message);
  }
}

inline void announce(const std::string& scenario)
{
  auto& ctx = Context::get();
  ctx.reported = 0;
  if (ctx.rank == 0) {
    std::cout << "=== " << scenario << std::endl;
  }
}

inline void skip(const std::string& scenario, const std::string& reason)
{
  if (Context::get().rank == 0) {
    std::cout << "=== " << scenario << " SKIP (" << reason << ")" << std::endl;
  }
}

/// Returns process-wide failure count, printing a rank 0 summary
inline long summarize()
{
  auto& ctx = Context::get();
  long failures = ctx.failures;
  long checks = ctx.checks;
#ifdef PARALLEL_MODE_MPI
  MPI_Allreduce(MPI_IN_PLACE, &failures, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &checks, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (ctx.rank == 0) {
    std::cout << (failures == 0 ? "PASSED" : "FAILED")
              << " (" << checks << " checks, " << failures << " failures)" << std::endl;
  }
  return failures;
}

inline std::vector<olb::Platform> availablePlatforms()
{
  std::vector<olb::Platform> platforms{olb::Platform::CPU_SISD};
#ifdef PLATFORM_CPU_SIMD
  platforms.emplace_back(olb::Platform::CPU_SIMD);
#endif
#ifdef PLATFORM_GPU_CUDA
  try {
    olb::checkPlatform<olb::Platform::GPU_CUDA>();
    platforms.emplace_back(olb::Platform::GPU_CUDA);
  } catch (std::exception&) { }
#endif
#ifdef PLATFORM_GPU_HIP
  try {
    olb::checkPlatform<olb::Platform::GPU_HIP>();
    platforms.emplace_back(olb::Platform::GPU_HIP);
  } catch (std::exception&) { }
#endif
#ifdef PLATFORM_GPU_SYCL
  try {
    olb::checkPlatform<olb::Platform::GPU_SYCL>();
    platforms.emplace_back(olb::Platform::GPU_SYCL);
  } catch (std::exception&) { }
#endif
  return platforms;
}

/// CPU_SISD and CPU_SIMD have no heterogeneous local copy task -- mixing SIMD
/// widths across blocks of one rank is not a supported configuration -- so the
/// homogeneous copy path would reinterpret the neighbour as the wrong column
/// type. Skip such pairs (homogeneous CPU_SIMD-CPU_SIMD is still exercised).
inline bool isHeterogeneousCpuPair(olb::Platform a, olb::Platform b)
{
  const auto isCpu = [](olb::Platform p) {
    return p == olb::Platform::CPU_SISD || p == olb::Platform::CPU_SIMD;
  };
  return a != b && isCpu(a) && isCpu(b);
}

inline int wrap(int g, int n)
{
  return ((g % n) + n) % n;
}

/// Injective float-exact stamp of a global cell's field component
inline float stamp(unsigned iField, olb::Vector<int,3> g, int n, unsigned iD)
{
  const std::size_t cell = (static_cast<std::size_t>(g[0])*n + g[1])*n + g[2];
  return static_cast<float>((std::size_t{1} << 22)*iField + (cell << 5) + iD);
}

constexpr float sentinel = -777.f;

/// Periodic unit-spacing mesh with per-cuboid platform assignment
template <typename T, unsigned D>
olb::Mesh<T,D> makeMesh(int n, int nCuboids, unsigned overlap,
                        std::function<olb::Platform(int)> platformOf)
{
  auto decomposition = std::make_unique<olb::CuboidDecomposition<T,D>>(
    olb::Vector<T,D>(0), T{1}, olb::Vector<int,D>(n), nCuboids);
  decomposition->setPeriodicity(olb::Vector<bool,D>(true));
  std::unique_ptr<olb::LoadBalancer<T>> balancer =
    std::make_unique<olb::BlockLoadBalancer<T>>(*decomposition);
  for (int iC=0; iC < balancer->size(); ++iC) {
    balancer->setPlatform(iC, platformOf(balancer->glob(iC)));
  }
  olb::Mesh<T,D> mesh(std::move(decomposition), std::move(balancer));
  mesh.setOverlap(overlap);
  return mesh;
}

/// Integer origin of a block's core in global lattice coordinates
template <typename T, unsigned D>
olb::Vector<int,D> blockOrigin(olb::Mesh<T,D>& mesh, int iC)
{
  auto physOrigin = mesh.getCuboidDecomposition().get(mesh.getLoadBalancer().glob(iC)).getOrigin();
  olb::Vector<int,D> origin;
  for (unsigned iDim=0; iDim < D; ++iDim) {
    origin[iDim] = static_cast<int>(std::lround(physOrigin[iDim]));
  }
  return origin;
}

}

#endif
