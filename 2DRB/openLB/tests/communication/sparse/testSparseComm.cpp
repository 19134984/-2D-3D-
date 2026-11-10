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

/* Sparse-field halo communication oracle.
 *
 * A sparse (tile-backed) field must communicate across cuboid halos identically
 * to a dense field carrying the same values: the type-erased any-fields kernels
 * gather/scatter/copy through the sparse device accessor rather than raw memory.
 *
 * The test builds a periodic multi-cuboid lattice with a DENSE reference field
 * (FORCE) and a SPARSE field of equal shape, stamps them identically on core
 * cells only, exchanges the overlap, and asserts the two fields agree bit-for-
 * bit on every cell -- core, transmitted halo (gather + scatter + copy_any_
 * fields), and untransmitted halo (absent tile -> default). Leaving the sparse
 * padding unstamped means the inbound halo tiles exist only if the communicator
 * pre-allocates them (ensureAllocated), which the device scatter requires since
 * the kernel never allocates.
 */

#include "../harness.h"

using namespace olb;
using namespace olb::names;

using T = float;

/// Dense reference field and its tile-backed sibling of identical shape.
/// WMVELOCITY (FieldStorage::Sparse, FIELD_BASE<0,1,0>) is a library field so the
/// mixed-mode device codegen sees it without folding in this test translation unit.
using SPARSE_FIELD = descriptors::WMVELOCITY;
using DENSE_FIELD  = descriptors::FORCE;

using DESCRIPTOR = descriptors::D3Q19<DENSE_FIELD,SPARSE_FIELD>;

constexpr int N = 12;

int paddingDepth(LatticeR<3> latticeR, Vector<int,3> extent)
{
  int depth = 0;
  for (unsigned iDim=0; iDim < 3; ++iDim) {
    depth = std::max(depth, -latticeR[iDim]);
    depth = std::max(depth, latticeR[iDim] - (extent[iDim]-1));
  }
  return depth;
}

using CaseD = Case<NavierStokes, Lattice<T,DESCRIPTOR>>;

/// Stamp core cells of FIELD with the injective global-coordinate stamp; leave
/// padding at the field's initial value (0). For the sparse field this keeps the
/// halo tiles unmaterialised until the communicator allocates the inbound ones.
template <typename FIELD>
void stampCore(CaseD& myCase)
{
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& load = lattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    const auto origin = test::blockOrigin(myCase.getMesh(), iC);
    auto& block = lattice.getBlock(iC);
    auto& array = block.template getField<FIELD>();
    block.forSpatialLocations([&](LatticeR<3> latticeR) {
      if (paddingDepth(latticeR, block.getExtent()) == 0) {
        const CellID cell = block.getCellId(latticeR);
        for (unsigned iD=0; iD < DESCRIPTOR::template size<FIELD>(); ++iD) {
          array[iD][cell] = test::stamp(0, origin + latticeR, N, iD);
        }
      }
    });
  }
}

/// Assert the sparse field agrees with the dense reference on every cell
void checkAgainstDense(CaseD& myCase, const std::string& name)
{
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& load = lattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    auto& block = lattice.getBlock(iC);
    auto& dense  = block.template getField<DENSE_FIELD>();
    auto& sparse = block.template getField<SPARSE_FIELD>();
    block.forSpatialLocations([&](LatticeR<3> latticeR) {
      const CellID cell = block.getCellId(latticeR);
      for (unsigned iD=0; iD < DESCRIPTOR::template size<SPARSE_FIELD>(); ++iD) {
        const float d = dense[iD][cell];
        const float s = sparse[iD][cell];
        if (std::memcmp(&d, &s, sizeof(float)) != 0) {
          std::ostringstream m;
          m << name << " block " << load.glob(iC)
            << " latticeR " << latticeR[0] << "," << latticeR[1] << "," << latticeR[2]
            << " iD " << iD << ": sparse " << s << " != dense " << d
            << " (padding " << paddingDepth(latticeR, block.getExtent()) << ")";
          test::fail(m.str());
        } else {
          test::Context::get().checks++;
        }
      }
    });
  }
}

void runOverlap(unsigned overlap, Platform a, Platform b, FieldStorage mode,
                sparse::TileShape shape = {1,1,64})
{
  struct SparseStage { };
  auto mesh = test::makeMesh<T,3>(N, 4, overlap, [&](int glob) { return glob % 2 == 0 ? a : b; });
  CaseD::ParametersD parameters;
  CaseD myCase(parameters, mesh);
  auto& lattice = myCase.getLattice(NavierStokes{});
  // Dense mode: the sparse-capable field is stored as a flat per-cell buffer and
  // must communicate identically to the dense reference with no explicit allocation.
  lattice.template setFieldStorage<SPARSE_FIELD>(mode);
  // A non-default tile shape must still communicate bit-identically to the dense reference
  if (mode == FieldStorage::Sparse) {
    lattice.template setTileShape<SPARSE_FIELD>(shape);
  }
  auto& communicator = lattice.getCommunicator(SparseStage{});
  communicator.requestField<DENSE_FIELD>();
  communicator.requestField<SPARSE_FIELD>();
  communicator.requestOverlap(overlap);
  communicator.exchangeRequests();

  stampCore<DENSE_FIELD>(myCase);
  stampCore<SPARSE_FIELD>(myCase);
  lattice.setProcessingContext(ProcessingContext::Simulation);

  // With FEATURES=SPARSE_SCRATCH_CANARY the scatter wrappers throw if any inbound halo
  // cell lands in the discard scratch tile (ensureAllocated missed a tile); a clean run
  // therefore proves every communicated sparse tile was pre-allocated.
  communicator.communicate();

  lattice.setProcessingContext(ProcessingContext::Evaluation);
  checkAgainstDense(myCase, "overlap" + std::to_string(overlap));
}

void runTests()
{
  const auto platforms = test::availablePlatforms();
  for (Platform a : platforms) {
    for (Platform b : platforms) {
      if (test::isHeterogeneousCpuPair(a, b)) {
        continue;
      }
      for (unsigned overlap : {1u, 2u}) {
        for (auto [mode,label] : {std::pair{FieldStorage::Sparse, "sparse"},
                                  std::pair{FieldStorage::Dense,  "dense"}}) {
          test::announce(std::string(label) + " overlap" + std::to_string(overlap) + " "
                         + std::string(getName(a)) + "-" + getName(b));
          runOverlap(overlap, a, b, mode);
        }
        // Non-default tile shape (4x4x4) on the sparse field must still match the dense reference
        test::announce(std::string("sparse-4x4x4 overlap") + std::to_string(overlap) + " "
                       + std::string(getName(a)) + "-" + getName(b));
        runOverlap(overlap, a, b, FieldStorage::Sparse, sparse::TileShape{4,4,4});
      }
    }
  }
}

int main(int argc, char* argv[])
{
  olb::initialize(&argc, &argv, false, false);
  auto& ctx = test::Context::get();
#ifdef PARALLEL_MODE_MPI
  ctx.rank = singleton::mpi().getRank();
  ctx.ranks = singleton::mpi().getSize();
#endif
  if (ctx.rank == 0) {
    std::cout << "platforms:";
    for (Platform p : test::availablePlatforms()) {
      std::cout << " " << getName(p);
    }
    std::cout << ", ranks: " << ctx.ranks << std::endl;
  }

  runTests();

  return test::summarize() > 0 ? 1 : 0;
}
