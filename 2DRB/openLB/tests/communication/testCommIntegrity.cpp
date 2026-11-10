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

#include "harness.h"
#include "twoPop.h"

using namespace olb;
using namespace olb::names;

using T = float;
using DESCRIPTOR = descriptors::D3Q19<descriptors::FORCE,descriptors::POROSITY>;

// f16 population storage is implemented on CPU_SISD/CPU_SIMD/GPU_CUDA only; skip
// the f16 section on HIP/SYCL builds so the store-precision static_assert holds.
#if !defined(PLATFORM_GPU_HIP) && !defined(PLATFORM_GPU_SYCL)
#define OLB_TEST_F16_STORE
#endif

constexpr int N = 12;

constexpr unsigned fieldComponents[3] = {DESCRIPTOR::q, DESCRIPTOR::d, 1};

constexpr std::uint64_t denseMask(unsigned iField)
{
  return (std::uint64_t{1} << fieldComponents[iField]) - 1;
}

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
using CommunicatorD = SuperCommunicator<T,SuperLattice<T,DESCRIPTOR>>;

/// One exchange micro-test: request recipe plus per-cell expected components
struct Recipe {
  std::string name;
  unsigned width;
  /// Bit iField set marks a selector-masked field (canary poisons its halo)
  std::uint64_t maskedFields;
  std::function<void(CaseD&, CommunicatorD&)> request;
  std::function<std::uint64_t(unsigned iField, LatticeR<3> latticeR,
                              Vector<int,3> extent, Vector<int,3> globalR)> expect;
};

template <typename FIELD>
void stampField(CaseD& myCase, unsigned iField)
{
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& load = lattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    const auto origin = test::blockOrigin(myCase.getMesh(), iC);
    auto& block = lattice.getBlock(iC);
    auto& array = block.template getField<FIELD>();
    block.forSpatialLocations([&](LatticeR<3> latticeR) {
      const CellID cell = block.getCellId(latticeR);
      const bool core = paddingDepth(latticeR, block.getExtent()) == 0;
      for (unsigned iD=0; iD < DESCRIPTOR::template size<FIELD>(); ++iD) {
        array[iD][cell] = core ? test::stamp(iField, origin + latticeR, N, iD)
                               : test::sentinel;
      }
    });
  }
}

template <typename FIELD>
void checkField(CaseD& myCase, const Recipe& recipe, unsigned iField)
{
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& load = lattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    const auto origin = test::blockOrigin(myCase.getMesh(), iC);
    auto& block = lattice.getBlock(iC);
    auto& array = block.template getField<FIELD>();
    block.forSpatialLocations([&](LatticeR<3> latticeR) {
      const CellID cell = block.getCellId(latticeR);
      const auto extent = block.getExtent();
      Vector<int,3> globalR = origin + latticeR;
      const bool core = paddingDepth(latticeR, extent) == 0;
      const std::uint64_t expected = core ? denseMask(iField)
                                          : recipe.expect(iField, latticeR, extent, globalR);
      Vector<int,3> wrapped;
      for (unsigned iDim=0; iDim < 3; ++iDim) {
        wrapped[iDim] = test::wrap(globalR[iDim], N);
      }
      for (unsigned iD=0; iD < DESCRIPTOR::template size<FIELD>(); ++iD) {
        const float value = array[iD][cell];
        const bool transmitted = expected & (std::uint64_t{1} << iD);
        const bool poisoned = overlap::canaryEnabled() && !transmitted && !core
                           && ((recipe.maskedFields >> iField) & 1)
                           && paddingDepth(latticeR, extent) <= static_cast<int>(recipe.width);
        float reference = transmitted ? test::stamp(iField, wrapped, N, iD)
                                      : test::sentinel;
        bool ok = poisoned ? std::isnan(value) : value == reference;
        if (poisoned) {
          reference = std::numeric_limits<float>::quiet_NaN();
        }
        if (!ok) {
          std::ostringstream s;
          s << recipe.name << " field " << iField
            << " block " << load.glob(iC)
            << " latticeR " << latticeR[0] << "," << latticeR[1] << "," << latticeR[2]
            << " iD " << iD << ": got " << value << " expected " << reference;
          test::fail(s.str());
        } else {
          test::Context::get().checks++;
        }
      }
    });
  }
}

void runRecipe(const Recipe& recipe, Platform a, Platform b)
{
  struct MicroStage { };
  auto mesh = test::makeMesh<T,3>(N, 4, 2, [&](int glob) { return glob % 2 == 0 ? a : b; });
  CaseD::ParametersD parameters;
  CaseD myCase(parameters, mesh);
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& communicator = lattice.getCommunicator(MicroStage{});
  recipe.request(myCase, communicator);
  communicator.exchangeRequests();

  stampField<descriptors::POPULATION>(myCase, 0);
  stampField<descriptors::FORCE>(myCase, 1);
  stampField<descriptors::POROSITY>(myCase, 2);
  lattice.setProcessingContext(ProcessingContext::Simulation);

  communicator.communicate();

  lattice.setProcessingContext(ProcessingContext::Evaluation, /*communicate=*/false);
  checkField<descriptors::POPULATION>(myCase, recipe, 0);
  checkField<descriptors::FORCE>(myCase, recipe, 1);
  checkField<descriptors::POROSITY>(myCase, recipe, 2);
}

std::vector<Recipe> makeRecipes()
{
  const auto streaming = [](unsigned iField, LatticeR<3> latticeR,
                            Vector<int,3> extent, Vector<int,3>) -> std::uint64_t {
    if (iField != 0 || paddingDepth(latticeR, extent) > 1) {
      return 0;
    }
    return overlap::requiredPopulations<DESCRIPTOR>(latticeR, extent);
  };

  return {
    {"dense-population", 1, 0,
     [](CaseD&, CommunicatorD& c) {
       c.requestField<descriptors::POPULATION>();
       c.requestOverlap(1);
     },
     [](unsigned iField, LatticeR<3> latticeR, Vector<int,3> extent, Vector<int,3>) -> std::uint64_t {
       return (iField == 0 && paddingDepth(latticeR, extent) <= 1) ? denseMask(0) : 0;
     }},
    {"streaming-population", 1, 0b001,
     [](CaseD&, CommunicatorD& c) {
       c.requestField<descriptors::POPULATION>(overlap::streamingMask<DESCRIPTOR>());
       c.requestOverlap(1);
     },
     streaming},
    {"multi-field-dense", 1, 0,
     [](CaseD&, CommunicatorD& c) {
       c.requestField<descriptors::FORCE>();
       c.requestField<descriptors::POPULATION>();
       c.requestField<descriptors::POROSITY>();
       c.requestOverlap(1);
     },
     [](unsigned, LatticeR<3> latticeR, Vector<int,3> extent, Vector<int,3>) -> std::uint64_t {
       return paddingDepth(latticeR, extent) <= 1 ? ~std::uint64_t{0} : 0;
     }},
    {"multi-field-streaming", 1, 0b001,
     [](CaseD&, CommunicatorD& c) {
       c.requestField<descriptors::FORCE>();
       c.requestField<descriptors::POPULATION>(overlap::streamingMask<DESCRIPTOR>());
       c.requestField<descriptors::POROSITY>();
       c.requestOverlap(1);
     },
     [streaming](unsigned iField, LatticeR<3> latticeR, Vector<int,3> extent, Vector<int,3> globalR) -> std::uint64_t {
       if (iField == 0) {
         return streaming(0, latticeR, extent, globalR);
       }
       return paddingDepth(latticeR, extent) <= 1 ? ~std::uint64_t{0} : 0;
     }},
    {"overlap2-dense", 2, 0,
     [](CaseD&, CommunicatorD& c) {
       c.requestField<descriptors::POPULATION>();
       c.requestOverlap(2);
     },
     [](unsigned iField, LatticeR<3> latticeR, Vector<int,3> extent, Vector<int,3>) -> std::uint64_t {
       return (iField == 0 && paddingDepth(latticeR, extent) <= 2) ? denseMask(0) : 0;
     }},
    {"overlap2-streaming", 2, 0b001,
     [](CaseD&, CommunicatorD& c) {
       c.requestField<descriptors::POPULATION>(overlap::streamingMask<DESCRIPTOR>());
       c.requestOverlap(2);
     },
     [](unsigned iField, LatticeR<3> latticeR, Vector<int,3> extent, Vector<int,3>) -> std::uint64_t {
       if (iField != 0 || paddingDepth(latticeR, extent) > 2) {
         return 0;
       }
       return overlap::requiredPopulations<DESCRIPTOR>(latticeR, extent);
     }},
    {"masked-force", 1, 0b010,
     [](CaseD&, CommunicatorD& c) {
       c.requestField<descriptors::POPULATION>();
       c.requestField<descriptors::FORCE>([](LatticeR<3> latticeR, LatticeR<3>) -> std::uint64_t {
         return latticeR[0] % 2 == 0 ? 0b101 : 0b010;
       });
       c.requestOverlap(1);
     },
     [](unsigned iField, LatticeR<3> latticeR, Vector<int,3> extent, Vector<int,3>) -> std::uint64_t {
       if (paddingDepth(latticeR, extent) > 1) {
         return 0;
       }
       switch (iField) {
       case 0:  return denseMask(0);
       case 1:  return latticeR[0] % 2 == 0 ? 0b101 : 0b010;
       default: return 0;
       }
     }},
    {"sparse-cells", 1, 0,
     [](CaseD& myCase, CommunicatorD& c) {
       c.requestField<descriptors::POPULATION>();
       auto& lattice = myCase.getLattice(NavierStokes{});
       auto& load = lattice.getLoadBalancer();
       for (int iC=0; iC < load.size(); ++iC) {
         const auto origin = test::blockOrigin(myCase.getMesh(), iC);
         auto& block = lattice.getBlock(iC);
         block.forSpatialLocations([&](LatticeR<3> latticeR) {
           const auto globalR = origin + latticeR;
           if (paddingDepth(latticeR, block.getExtent()) == 1
               && (globalR[0] + globalR[1] + globalR[2]) % 2 == 0) {
             c.requestCell(latticeR.withPrefix(iC));
           }
         });
       }
     },
     [](unsigned iField, LatticeR<3> latticeR, Vector<int,3> extent, Vector<int,3> globalR) -> std::uint64_t {
       if (iField != 0 || paddingDepth(latticeR, extent) != 1
           || (globalR[0] + globalR[1] + globalR[2]) % 2 != 0) {
         return 0;
       }
       return denseMask(0);
     }},
  };
}

void runMicroTests()
{
  const auto platforms = test::availablePlatforms();
  const auto recipes = makeRecipes();
  for (Platform a : platforms) {
    for (Platform b : platforms) {
      if (test::isHeterogeneousCpuPair(a, b)) {
        continue;
      }
      for (const auto& recipe : recipes) {
        test::announce("micro " + recipe.name + " " + std::string(getName(a)) + "-" + getName(b));
        runRecipe(recipe, a, b);
      }
    }
  }
}

/// Physics A/B: N steps of periodic BGK, core populations must be
/// bitwise invariant under the directional streaming exchange
std::map<int,std::vector<float>> runPhysics(std::function<Platform(int)> platformOf,
                                            bool directional)
{
  auto mesh = test::makeMesh<T,3>(2*N, 4, 2, std::move(platformOf));
  CaseD::ParametersD parameters;
  CaseD myCase(parameters, mesh);

  auto& geometry = myCase.getGeometry();
  geometry.rename(0, 1);
  auto& lattice = myCase.getLattice(NavierStokes{});
  // Directional streaming exchange is the default; force the dense reference
  // through the same runtime path the reversible-streaming fallback uses.
  if (!directional) {
    lattice.requireDenseStreamingExchange();
  }
  UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> converter(
    2*N, 0.55, 1., 1., 0.1, 1.);
  lattice.setUnitConverter<DESCRIPTOR>(converter);
  dynamics::set<BGKdynamics>(lattice, geometry, 1);
  lattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  auto& load = lattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    const auto origin = test::blockOrigin(myCase.getMesh(), iC);
    auto& block = lattice.getBlock(iC);
    auto& populations = block.template getField<descriptors::POPULATION>();
    block.forSpatialLocations([&](LatticeR<3> latticeR) {
      const CellID cell = block.getCellId(latticeR);
      const auto globalR = origin + latticeR;
      const float bump = 0.01f * std::sin(2*M_PI*globalR[0]/(2.f*N))
                               * std::cos(2*M_PI*globalR[1]/(2.f*N))
                               * std::cos(4*M_PI*globalR[2]/(2.f*N));
      for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        populations[iPop][cell] = bump * descriptors::t<T,DESCRIPTOR>(iPop);
      }
    });
  }
  lattice.setProcessingContext(ProcessingContext::Simulation);
  lattice.initialize();

  for (int iT=0; iT < 48; ++iT) {
    lattice.collideAndStream();
  }
  lattice.setProcessingContext(ProcessingContext::Evaluation, /*communicate=*/false);

  std::map<int,std::vector<float>> result;
  for (int iC=0; iC < load.size(); ++iC) {
    auto& block = lattice.getBlock(iC);
    auto& populations = block.template getField<descriptors::POPULATION>();
    auto& values = result[load.glob(iC)];
    block.forSpatialLocations([&](LatticeR<3> latticeR) {
      if (paddingDepth(latticeR, block.getExtent()) == 0) {
        const CellID cell = block.getCellId(latticeR);
        for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
          values.emplace_back(populations[iPop][cell]);
        }
      }
    });
  }
  return result;
}

void runPhysicsTests()
{
  const auto platforms = test::availablePlatforms();
  std::vector<std::pair<std::string,std::function<Platform(int)>>> patterns;
  for (Platform p : platforms) {
    patterns.emplace_back(getName(p), [p](int) { return p; });
  }
  for (Platform a : platforms) {
    for (Platform b : platforms) {
      if (a != b && !test::isHeterogeneousCpuPair(a, b)) {
        patterns.emplace_back(std::string(getName(a)) + "/" + getName(b),
                              [a,b](int glob) { return glob % 2 == 0 ? a : b; });
      }
    }
  }

  // CPU_SISD has no cyclic ring, so it is the shift-independent reference a
  // same-platform dense-vs-directional comparison cannot provide.
  const float cpuTol = 1e-4f;
  auto cpuReference = runPhysics([](int) { return olb::Platform::CPU_SISD; }, true);

  for (auto& [name, pattern] : patterns) {
    test::announce("physics dense-vs-directional " + name);
    auto reference = runPhysics(pattern, false);
    auto subject = runPhysics(pattern, true);

    float maxCpuDiff = 0.f;
    for (auto& [glob, values] : subject) {
      auto& cpuValues = cpuReference[glob];
      const std::size_t n = std::min(values.size(), cpuValues.size());
      for (std::size_t i=0; i < n; ++i) {
        test::Context::get().checks++;
        maxCpuDiff = std::max(maxCpuDiff, std::abs(values[i] - cpuValues[i]));
      }
    }
    test::announce("  vs CPU reference: max |diff| = " + std::to_string(maxCpuDiff));
    test::check(maxCpuDiff < cpuTol,
                name + ": core populations diverge from CPU reference by "
                + std::to_string(maxCpuDiff) + " (tol " + std::to_string(cpuTol) + ")");

    for (auto& [glob, values] : reference) {
      auto& subjectValues = subject[glob];
      test::check(values.size() == subjectValues.size(),
                  "core cell count mismatch in block " + std::to_string(glob));
      std::size_t mismatches = 0;
      for (std::size_t i=0; i < values.size(); ++i) {
        test::Context::get().checks++;
        if (std::memcmp(&values[i], &subjectValues[i], sizeof(float)) != 0) {
          mismatches++;
        }
      }
      if (mismatches > 0) {
        test::fail("block " + std::to_string(glob) + ": "
                   + std::to_string(mismatches) + " of " + std::to_string(values.size())
                   + " core values differ");
      }
      bool nanFree = std::none_of(subjectValues.begin(), subjectValues.end(),
                                  [](float v) { return std::isnan(v); });
      test::check(nanFree, "block " + std::to_string(glob) + ": NaN in core");
    }
  }
}

/// Stage 2: two propagatable fields (POPULATION + POPULATION2) sharing the one
/// descriptor stencil and the one negotiated streaming mask. Proves the masked
/// halo exchange generalised from hardcoded POPULATION to the propagatable set
/// end-to-end across every available platform pair, each field receiving R(p)
/// through its own masked segment.
namespace twofield {

using DESCRIPTOR2 = descriptors::D3Q19<POPULATION2>;
using Case2 = Case<NavierStokes, Lattice<T,DESCRIPTOR2>>;

template <typename FIELD>
void stamp(Case2& myCase, unsigned iField)
{
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& load = lattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    const auto origin = test::blockOrigin(myCase.getMesh(), iC);
    auto& block = lattice.getBlock(iC);
    auto& array = block.template getField<FIELD>();
    block.forSpatialLocations([&](LatticeR<3> latticeR) {
      const CellID cell = block.getCellId(latticeR);
      const bool core = paddingDepth(latticeR, block.getExtent()) == 0;
      for (unsigned iD=0; iD < DESCRIPTOR2::template size<FIELD>(); ++iD) {
        array[iD][cell] = core ? test::stamp(iField, origin + latticeR, N, iD)
                               : test::sentinel;
      }
    });
  }
}

template <typename FIELD>
void check(Case2& myCase, unsigned iField)
{
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& load = lattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    const auto origin = test::blockOrigin(myCase.getMesh(), iC);
    auto& block = lattice.getBlock(iC);
    auto& array = block.template getField<FIELD>();
    block.forSpatialLocations([&](LatticeR<3> latticeR) {
      const CellID cell = block.getCellId(latticeR);
      const auto extent = block.getExtent();
      Vector<int,3> globalR = origin + latticeR;
      const int depth = paddingDepth(latticeR, extent);
      const bool core = depth == 0;
      const std::uint64_t dense = (std::uint64_t{1} << DESCRIPTOR2::q) - 1;
      const std::uint64_t expected = core ? dense
                                    : (depth <= 1 ? overlap::requiredPopulations<DESCRIPTOR2>(latticeR, extent)
                                                  : 0);
      Vector<int,3> wrapped;
      for (unsigned iDim=0; iDim < 3; ++iDim) {
        wrapped[iDim] = test::wrap(globalR[iDim], N);
      }
      for (unsigned iD=0; iD < DESCRIPTOR2::template size<FIELD>(); ++iD) {
        const float value = array[iD][cell];
        const bool transmitted = expected & (std::uint64_t{1} << iD);
        const bool poisoned = overlap::canaryEnabled() && !transmitted && !core && depth <= 1;
        float reference = transmitted ? test::stamp(iField, wrapped, N, iD)
                                      : test::sentinel;
        bool ok = poisoned ? std::isnan(value) : value == reference;
        if (!ok) {
          std::ostringstream s;
          s << "two-field field " << iField << " block " << load.glob(iC)
            << " latticeR " << latticeR[0] << "," << latticeR[1] << "," << latticeR[2]
            << " iD " << iD << ": got " << value;
          test::fail(s.str());
        } else {
          test::Context::get().checks++;
        }
      }
    });
  }
}

void run(Platform a, Platform b)
{
  struct MicroStage { };
  auto mesh = test::makeMesh<T,3>(N, 4, 2, [&](int glob) { return glob % 2 == 0 ? a : b; });
  Case2::ParametersD parameters;
  Case2 myCase(parameters, mesh);
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& communicator = lattice.getCommunicator(MicroStage{});
  communicator.requestField<descriptors::POPULATION>(overlap::streamingMask<DESCRIPTOR2>());
  communicator.requestField<POPULATION2>(overlap::streamingMask<DESCRIPTOR2>());
  communicator.requestOverlap(1);
  communicator.exchangeRequests();

  stamp<descriptors::POPULATION>(myCase, 0);
  stamp<POPULATION2>(myCase, 1);
  lattice.setProcessingContext(ProcessingContext::Simulation);

  communicator.communicate();

  lattice.setProcessingContext(ProcessingContext::Evaluation, /*communicate=*/false);
  check<descriptors::POPULATION>(myCase, 0);
  check<POPULATION2>(myCase, 1);
}

/// Two propagatable fields carrying genuinely different selectors: POPULATION
/// streams as usual while POPULATION2 exchanges nothing. A single shared mask
/// channel would alias POPULATION2 onto POPULATION's mask and wrongly fill its
/// halo -- this pins the per-field channels.
void runDifferentSelectors(Platform a, Platform b)
{
  struct MicroStage { };
  auto mesh = test::makeMesh<T,3>(N, 4, 2, [&](int glob) { return glob % 2 == 0 ? a : b; });
  Case2::ParametersD parameters;
  Case2 myCase(parameters, mesh);
  auto& lattice = myCase.getLattice(NavierStokes{});
  auto& communicator = lattice.getCommunicator(MicroStage{});
  communicator.requestField<descriptors::POPULATION>(overlap::streamingMask<DESCRIPTOR2>());
  communicator.requestField<POPULATION2>([](LatticeR<3>, LatticeR<3>) -> std::uint64_t { return 0; });
  communicator.requestOverlap(1);
  communicator.exchangeRequests();

  stamp<descriptors::POPULATION>(myCase, 0);
  stamp<POPULATION2>(myCase, 1);
  lattice.setProcessingContext(ProcessingContext::Simulation);
  communicator.communicate();
  lattice.setProcessingContext(ProcessingContext::Evaluation, /*communicate=*/false);

  check<descriptors::POPULATION>(myCase, 0);

  // POPULATION2 exchanged nothing: overlap cells keep the sentinel (NaN under canary)
  auto& load = lattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    auto& block = lattice.getBlock(iC);
    auto& array = block.template getField<POPULATION2>();
    block.forSpatialLocations([&](LatticeR<3> latticeR) {
      const CellID cell = block.getCellId(latticeR);
      const int depth = paddingDepth(latticeR, block.getExtent());
      if (depth == 0) {
        return;
      }
      const bool poisoned = overlap::canaryEnabled() && depth <= 1;
      for (unsigned iD=0; iD < DESCRIPTOR2::q; ++iD) {
        const float value = array[iD][cell];
        if (poisoned ? std::isnan(value) : value == test::sentinel) {
          test::Context::get().checks++;
        } else {
          test::fail("POPULATION2 halo exchanged despite a zero selector");
        }
      }
    });
  }
}

}

void runTwoFieldTests()
{
  const auto platforms = test::availablePlatforms();
  for (Platform a : platforms) {
    for (Platform b : platforms) {
      if (test::isHeterogeneousCpuPair(a, b)) {
        continue;
      }
      test::announce("two-field-streaming " + std::string(getName(a)) + "-" + getName(b));
      twofield::run(a, b);
      test::announce("two-field-different-selectors " + std::string(getName(a)) + "-" + getName(b));
      twofield::runDifferentSelectors(a, b);
    }
  }
}

#ifdef OLB_TEST_F16_STORE
/// f16 population storage across the halo exchange. Guards the halved-traffic
/// serialize/deserialize (cpu/sisd, cpu/simd F16C, gpu/cuda __half) and the
/// storage round-trip: on each platform the masked directional exchange must
/// move the exact stored f16 bytes the dense path does, and every platform must
/// agree with the CPU_SISD f16 reference within f16 tolerance.
namespace f16store {

using DESCRIPTOR16 = descriptors::D3Q19<
  descriptors::tag::store_precision<descriptors::POPULATION, FieldPrecision::F16>>;
using Case16 = Case<NavierStokes, Lattice<T,DESCRIPTOR16>>;

std::map<int,std::vector<float>> run(std::function<Platform(int)> platformOf, bool directional)
{
  auto mesh = test::makeMesh<T,3>(2*N, 4, 2, std::move(platformOf));
  Case16::ParametersD parameters;
  Case16 myCase(parameters, mesh);

  auto& geometry = myCase.getGeometry();
  geometry.rename(0, 1);
  auto& lattice = myCase.getLattice(NavierStokes{});
  if (!directional) {
    lattice.requireDenseStreamingExchange();
  }
  UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR16> converter(
    2*N, 0.55, 1., 1., 0.1, 1.);
  lattice.setUnitConverter<DESCRIPTOR16>(converter);
  dynamics::set<BGKdynamics>(lattice, geometry, 1);
  lattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  auto& load = lattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    const auto origin = test::blockOrigin(myCase.getMesh(), iC);
    auto& block = lattice.getBlock(iC);
    auto& populations = block.template getField<descriptors::POPULATION>();
    block.forSpatialLocations([&](LatticeR<3> latticeR) {
      const CellID cell = block.getCellId(latticeR);
      const auto globalR = origin + latticeR;
      const float bump = 0.01f * std::sin(2*M_PI*globalR[0]/(2.f*N))
                               * std::cos(2*M_PI*globalR[1]/(2.f*N))
                               * std::cos(4*M_PI*globalR[2]/(2.f*N));
      for (unsigned iPop=0; iPop < DESCRIPTOR16::q; ++iPop) {
        populations[iPop][cell] = bump * descriptors::t<T,DESCRIPTOR16>(iPop);
      }
    });
  }
  lattice.setProcessingContext(ProcessingContext::Simulation);
  lattice.initialize();

  for (int iT=0; iT < 48; ++iT) {
    lattice.collideAndStream();
  }
  lattice.setProcessingContext(ProcessingContext::Evaluation, /*communicate=*/false);

  std::map<int,std::vector<float>> result;
  for (int iC=0; iC < load.size(); ++iC) {
    auto& block = lattice.getBlock(iC);
    auto& populations = block.template getField<descriptors::POPULATION>();
    auto& values = result[load.glob(iC)];
    block.forSpatialLocations([&](LatticeR<3> latticeR) {
      if (paddingDepth(latticeR, block.getExtent()) == 0) {
        const CellID cell = block.getCellId(latticeR);
        for (unsigned iPop=0; iPop < DESCRIPTOR16::q; ++iPop) {
          values.emplace_back(float(populations[iPop][cell]));
        }
      }
    });
  }
  return result;
}

void runTests()
{
  const auto platforms = test::availablePlatforms();
  std::vector<std::pair<std::string,std::function<Platform(int)>>> patterns;
  for (Platform p : platforms) {
    patterns.emplace_back(getName(p), [p](int) { return p; });
  }
  for (Platform a : platforms) {
    for (Platform b : platforms) {
      if (a != b && !test::isHeterogeneousCpuPair(a, b)) {
        patterns.emplace_back(std::string(getName(a)) + "/" + getName(b),
                              [a,b](int glob) { return glob % 2 == 0 ? a : b; });
      }
    }
  }

  // f16 quantises the stored populations, so cross-platform agreement is checked
  // within an f16 tolerance -- gross enough to catch structural breakage, loose
  // enough to survive CPU-F16C vs device-__half rounding divergence over 48 steps.
  const float f16Tol = 5e-3f;
  auto cpuReference = run([](int) { return olb::Platform::CPU_SISD; }, true);

  for (auto& [name, pattern] : patterns) {
    test::announce("f16-store dense-vs-directional " + name);
    auto reference = run(pattern, false);
    auto subject = run(pattern, true);

    // Directional (masked, halved-traffic) f16 comm must move the exact stored
    // bytes the dense reference does -> bitwise identity per platform.
    for (auto& [glob, values] : reference) {
      auto& subjectValues = subject[glob];
      test::check(values.size() == subjectValues.size(),
                  "f16 core cell count mismatch in block " + std::to_string(glob));
      std::size_t mismatches = 0;
      for (std::size_t i=0; i < values.size(); ++i) {
        test::Context::get().checks++;
        if (std::memcmp(&values[i], &subjectValues[i], sizeof(float)) != 0) {
          mismatches++;
        }
      }
      if (mismatches > 0) {
        test::fail("f16 block " + std::to_string(glob) + ": "
                   + std::to_string(mismatches) + " of " + std::to_string(values.size())
                   + " core values differ (dense vs directional)");
      }
      bool nanFree = std::none_of(subjectValues.begin(), subjectValues.end(),
                                  [](float v) { return std::isnan(v); });
      test::check(nanFree, "f16 block " + std::to_string(glob) + ": NaN in core");
    }

    // Every platform agrees with the CPU_SISD f16 reference within f16 tolerance.
    float maxCpuDiff = 0.f;
    for (auto& [glob, values] : subject) {
      auto& cpuValues = cpuReference[glob];
      const std::size_t n = std::min(values.size(), cpuValues.size());
      for (std::size_t i=0; i < n; ++i) {
        test::Context::get().checks++;
        maxCpuDiff = std::max(maxCpuDiff, std::abs(values[i] - cpuValues[i]));
      }
    }
    test::announce("  vs CPU f16 reference: max |diff| = " + std::to_string(maxCpuDiff));
    test::check(maxCpuDiff < f16Tol,
                name + ": f16 core populations diverge from CPU f16 reference by "
                + std::to_string(maxCpuDiff) + " (tol " + std::to_string(f16Tol) + ")");
  }
}

}
#endif  // OLB_TEST_F16_STORE

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
    std::cout << ", ranks: " << ctx.ranks
              << ", canary: " << (overlap::canaryEnabled() ? "on" : "off") << std::endl;
  }

  runMicroTests();
  runPhysicsTests();
  runTwoFieldTests();
#ifdef OLB_TEST_F16_STORE
  f16store::runTests();
#endif

  return test::summarize() > 0 ? 1 : 0;
}
