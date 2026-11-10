/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022-2026 Fedor Bukreev, Adrian Kummerlaender
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

#ifndef BOUNDARIES_VORTEX_METHOD_OPTIMIZED_H
#define BOUNDARIES_VORTEX_METHOD_OPTIMIZED_H

#include "utilities/functorPtr.h"
#include "functors/lattice/indicator/superIndicatorF3D.h"

#include "rngHelper.h"

namespace olb {
namespace optimized {

struct VortexMethodPreProcessor;

namespace stage {

struct VortexMethod { };

}

namespace VortexMethod {

struct U_PROFILE : public descriptors::FIELD_BASE<0,1> { };
struct U_PROFILE_SCALE : public descriptors::FIELD_BASE<1> { };
struct CONVERSION_FACTOR_LENGTH : public descriptors::FIELD_BASE<1> { };
struct CONVERSION_FACTOR_VELOCITY : public descriptors::FIELD_BASE<1> { };
struct SEEDS : public descriptors::TEMPLATE_FIELD_BASE<std::add_pointer_t,0,1> { };
struct SEEDS_VORTICITY : public descriptors::TEMPLATE_FIELD_BASE<std::add_pointer_t,1> { };
struct SEEDS_AXIS_DIRECTION : public descriptors::TEMPLATE_FIELD_BASE<std::add_pointer_t,0,1> { };
struct SEEDS_SIGMA : public descriptors::TEMPLATE_FIELD_BASE<std::add_pointer_t,1> { };
struct AXIS_DIRECTION : public descriptors::FIELD_BASE<0,1> { };

struct MAX_R : public descriptors::FIELD_BASE<0,1> { };
struct MIN_R : public descriptors::FIELD_BASE<0,1> { };
struct L_MIN : public descriptors::FIELD_BASE<1> { };
struct L_CHAR : public descriptors::FIELD_BASE<1> { };
struct U_CONVECTION : public descriptors::FIELD_BASE<1> { };
struct REFERENCE_HEIGHT : public descriptors::FIELD_BASE<1> { };
struct SEED_COUNT : public descriptors::TYPED_FIELD_BASE<std::size_t,1> { };
struct INLET_AREA : public descriptors::FIELD_BASE<1> { };
struct INTENSITY : public descriptors::FIELD_BASE<1> { };
struct CUTOFF_RATIO : public descriptors::FIELD_BASE<1> { };
struct SEED : public descriptors::TYPED_FIELD_BASE<std::uint64_t,1> { };

}

struct VortexMethodSeedConvectionO {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  using parameters = meta::list<VortexMethod::MAX_R,
                                VortexMethod::MIN_R,
                                VortexMethod::REFERENCE_HEIGHT,
                                VortexMethod::L_MIN,
                                VortexMethod::L_CHAR,
                                VortexMethod::U_CONVECTION,
                                VortexMethod::SEED_COUNT,
                                VortexMethod::INLET_AREA,
                                VortexMethod::INTENSITY,
                                VortexMethod::SEED,
                                fields::converter::PHYS_DELTA_T>;

  template <typename CELLS, typename PARAMETERS>
  void apply(CELLS& cells, PARAMETERS& parameters) any_platform
  {
    using V = typename CELLS::template value_t<names::Points>::value_t;
#if defined(OLB_PRECOMPILING_SYCL) || defined(__HIPCC__) || defined(__CUDACC__)
    const Vector<V,3> maxR = parameters.template get<VortexMethod::MAX_R>();
    const Vector<V,3> minR = parameters.template get<VortexMethod::MIN_R>();

    const auto height = maxR[2];
    const auto uConvection = parameters.template get<VortexMethod::U_CONVECTION>();
    const auto lMin = parameters.template get<VortexMethod::L_MIN>();
    const auto lChar = parameters.template get<VortexMethod::L_CHAR>();
    const auto deltaT = parameters.template get<fields::converter::PHYS_DELTA_T>();
    const auto seed = parameters.template get<VortexMethod::SEED>();
    const auto nSeeds = parameters.template get<VortexMethod::SEED_COUNT>();
    const auto inletArea = parameters.template get<VortexMethod::INLET_AREA>();
    auto intensity = parameters.template get<VortexMethod::INTENSITY>();

    auto point = cells.template get<names::Points>();
    auto physR = point.template getField<fields::PHYS_R>();
    auto oldPhysR = point.template getField<fields::PHYS_R>();

    physR[0] += uConvection * deltaT;

    if (physR[0] < maxR[0]) {
      point.template setField<fields::PHYS_R>(physR);
      return;

    }

    PlatformRNG<V> rng(seed, point.getCellId());

    physR[0] = minR[0];
    physR[1] = minR[1] + rng.uniform() * (maxR[1] - minR[1]);
    physR[2] = minR[2] + rng.uniform() * (maxR[2] - minR[2]);

    Vector<V,3> velocity;
    const V z0   = 0.001;
    const V zRef = parameters.template get<VortexMethod::REFERENCE_HEIGHT>();
    const V uRef = uConvection;
    velocity[0] = uRef * ((util::log(physR[2]/z0)/util::log(zRef/z0)));

    V normVel = util::norm<3>(velocity);
    V kinE = 3./2. * util::pow(normVel * intensity, 2.);
    V sigmaI = 0;

    if (normVel > 1e-6) {
      V mixingLength = 0.4 * physR[2];
      mixingLength = util::max(mixingLength, lMin);
      mixingLength = util::min(mixingLength, lChar);

      const V eps = util::pow(0.09, 3./4.) * util::pow(kinE, 3./2.) / mixingLength;
      sigmaI = 0.16*util::pow(kinE, 3./2.) / 2. / eps;
      sigmaI *= rng.uniform_size();
    }

    V vorticity = 4. * util::sqrt(std::numbers::pi_v<V> * inletArea * kinE / 3. / nSeeds / (2.*util::log(3.) - 3.*util::log(2.)));

    Vector<V,3> axis{};
    V phi = rng.uniform() * 2 * std::numbers::pi_v<V>;
    V z   = rng.uniform() * 2 - 1;
    V r   = util::sqrt(1 - z*z);
    axis[0] = r * util::cos(phi);
    axis[1] = r * util::sin(phi);
    axis[2] = z;

    //sycl::ext::oneapi::experimental::printf("physR [%f %f %f]; axis [%f %f %f]; vorticity %f; sigma %f\n",
    //                                        physR[0], physR[1], physR[2],
    //                                        axis[0], axis[1], axis[2],
    //                                        vorticity,
    //                                        sigmaI);

    point.template setField<fields::PHYS_R>(physR);
    point.template setField<descriptors::VORTICITY>(vorticity);
    point.template setField<descriptors::SIGMA>(sigmaI);
    point.template setField<VortexMethod::AXIS_DIRECTION>(axis);
#endif
  }

};

template <typename T, typename DESCRIPTOR>
class VortexMethodTurbulentVelocityBoundary final {
private:
  FunctorPtr<SuperIndicatorF3D<T>> _inletLatticeI;
  FunctorPtr<IndicatorF3D<T>> _inletPhysI;
  std::shared_ptr<AnalyticalF3D<T,T>> _velocityProfileF;
  std::shared_ptr<AnalyticalF3D<T,T>> _intensityProfileF;

  const UnitConverter<T,DESCRIPTOR>& _converter;
  SuperLattice<T, DESCRIPTOR>& _sLattice;

  const int _nSeeds;
  T _inletArea;
  T _sigmaConst;
  T _lMin;
  T _lChar;
  T _uConvection;

  Vector<T,3> _axisDirection;

  using SEED_DESCRIPTOR = descriptors::D3<fields::PHYS_R,
                                          descriptors::VORTICITY,
                                          descriptors::SIGMA,
                                          VortexMethod::AXIS_DIRECTION>;
  SuperD<T,SEED_DESCRIPTOR> _seedsD;

  SuperLatticePointCoupling<
    VortexMethodSeedConvectionO,
    meta::map<names::Lattice1, descriptors::VALUED_DESCRIPTOR<T,DESCRIPTOR>,
              names::Points,   descriptors::VALUED_DESCRIPTOR<T,SEED_DESCRIPTOR>>
  > _convectionO;

  std::mt19937 _generator;

  void generateSeeds();

public:
  VortexMethodTurbulentVelocityBoundary(
    FunctorPtr<SuperIndicatorF3D<T>>&& inletLatticeI,
    FunctorPtr<IndicatorF3D<T>>&& inletPhysI,
    const UnitConverter<T,DESCRIPTOR>& converter,
    SuperLattice<T, DESCRIPTOR>& sLattice,
    int nSeeds,
    T uConvection,
    T sigmaConst,
    T lMin,
    T lChar,
    T inletArea,
    Vector<T,3> axisDirection);

  void setVelocityProfile(std::shared_ptr<AnalyticalF3D<T,T>> velocityProfileF);
  void setIntensityProfile(std::shared_ptr<AnalyticalF3D<T,T>> intensityProfileF);
  void setConvectionU(T uConvection);

  void apply(std::size_t iT);

};

template <typename T, typename DESCRIPTOR>
VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR>::VortexMethodTurbulentVelocityBoundary(
  FunctorPtr<SuperIndicatorF3D<T>>&& inletLatticeI,
  FunctorPtr<IndicatorF3D<T>>&& inletPhysI,
  const UnitConverter<T,DESCRIPTOR>& converter,
  SuperLattice<T, DESCRIPTOR>& sLattice,
  int nSeeds,
  T uConvection,
  T sigmaConst,
  T lMin,
  T lChar,
  T inletArea,
  Vector<T,3> axisDirection)
  : _inletLatticeI(std::move(inletLatticeI)),
    _inletPhysI(std::move(inletPhysI)),
    _converter(converter),
    _sLattice(sLattice),
    _nSeeds(nSeeds),
    _uConvection(uConvection),
    _sigmaConst(sigmaConst),
    _lMin(lMin),
    _lChar(lChar),
    _inletArea(inletArea),
    _axisDirection(axisDirection),
    _seedsD(sLattice.getLoadBalancer()),
    _convectionO(VortexMethodSeedConvectionO{},
                 names::Lattice1{}, _sLattice,
                 names::Points{},   _seedsD)
{
  OstreamManager clout(std::cout, "VortexMethod");

  {
    std::seed_seq seq{{0x9e3779b9u, 0xedb7d151u, 0x8212c45du, 0x33ca8075u,
                       0x11122233u, 0x44556677u, 0x8899aabbu, 0xccddeeffu}};
    _generator.seed(seq);
    _generator.discard(1e5);
  }

  sLattice.template addPostProcessor<::olb::stage::PostStream>(*_inletLatticeI,
                                                        meta::id<VortexMethodPreProcessor>{});

  auto& cGeometry = sLattice.getCuboidDecomposition();

  auto& load = _sLattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    T blockArea = 0.;
    auto& cuboid = cGeometry.get(load.glob(iC));
    auto& block = _sLattice.getBlock(iC);
    for (int iX=0; iX < block.getNx(); ++iX) {
      for (int iY=0; iY < block.getNy(); ++iY) {
        for (int iZ=0; iZ < block.getNz(); ++iZ) {
          int latticeR[3] {iX, iY, iZ};
          auto physR = cuboid.getPhysR(latticeR);
          block.get(iX,iY,iZ).template setField<fields::PHYS_R>(physR);
        }
      }
    }
    block.template getData<Array<fields::PHYS_R>>()
         .setProcessingContext(ProcessingContext::Simulation);
  }

  #ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(_inletArea, MPI_SUM);
  #endif

  clout << "inletArea=" << _inletArea << std::endl;

  for (int iC=0; iC < load.size(); ++iC) {
    _seedsD.getBlock(iC).resize({_nSeeds,1,1});
    _seedsD.getBlock(iC).setProcessingContext(ProcessingContext::Simulation);
  }

  for (int iC=0; iC < load.size(); ++iC) {
    if (_nSeeds > 0) {
      auto& block = _sLattice.getBlock(iC);
      auto& seeds = _seedsD.getBlock(iC);
      block.template setParameter<VortexMethod::SEED_COUNT>(_nSeeds);
      block.template setParameter<VortexMethod::SEEDS>(seeds.template getField<fields::PHYS_R>());
      block.template setParameter<VortexMethod::SEEDS_VORTICITY>(seeds.template getField<descriptors::VORTICITY>());
      block.template setParameter<VortexMethod::SEEDS_SIGMA>(seeds.template getField<descriptors::SIGMA>());
      block.template setParameter<VortexMethod::SEEDS_AXIS_DIRECTION>(seeds.template getField<VortexMethod::AXIS_DIRECTION>());
      block.template getField<VortexMethod::U_PROFILE>();
    }
  }

  _sLattice.template setParameter<VortexMethod::CONVERSION_FACTOR_LENGTH>(
    _converter.getConversionFactorLength());
  _sLattice.template setParameter<VortexMethod::CONVERSION_FACTOR_VELOCITY>(
    _converter.getConversionFactorVelocity());
  _sLattice.template setParameter<VortexMethod::CUTOFF_RATIO>(-2.0 * util::log(1e-5));
}

template <typename T, typename DESCRIPTOR>
void VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR>::generateSeeds()
{
  std::uniform_real_distribution<T> distribution(0, 1);
  std::uniform_real_distribution<T> sizeDistribution(0.5, 1.5);
  auto randomness = [&]() -> T {
    return distribution(_generator);
  };

  std::vector<Vector<T,3>> initPositions(_nSeeds);
  std::vector<T> initVorticities(_nSeeds);
  std::vector<T> initSigmas(_nSeeds);
  std::vector<Vector<T,3>> initAxes(_nSeeds);

  for (int j = 0; j < _nSeeds; j++) {
    auto inlet = _inletPhysI->getSample(randomness);
    initPositions[j] = inlet;

    Vector<T,3> velocity;
    _velocityProfileF->operator()(velocity.data(), inlet.data());
    Vector<T,1> intensity;
    _intensityProfileF->operator()(intensity.data(), inlet.data());

    T normVel = olb::norm(velocity);
    T kinE = 3./2. * util::pow(normVel * intensity[0], 2. );
    T sigmaI = _sigmaConst;

    if (_sigmaConst == 0 && normVel > 1e-6) {
      T mixingLength = 0.4 * util::abs(inlet[2]);
      mixingLength = util::max(mixingLength, _lMin);
      mixingLength = util::min(mixingLength, _lChar);

      const T eps = util::pow(0.09, 3./4.) * util::pow(kinE, 3./2.) / mixingLength;
      sigmaI = 0.16*util::pow(kinE, 3./2.) / 2. / eps;
      sigmaI *= sizeDistribution(_generator);
    }

    T vorticity = 4. * util::sqrt(std::numbers::pi_v<T> * _inletArea * kinE / 3. / _nSeeds / (2.*util::log(3.) - 3.*util::log(2.)));

    initVorticities[j] = vorticity;
    initSigmas[j] = sigmaI;

    Vector<T,3> axis{};
    T phi = randomness() * 2 * std::numbers::pi_v<T>;
    T z   = randomness() * 2 - 1;
    T r   = util::sqrt(1 - z*z);
    axis[0] = r * util::cos(phi);
    axis[1] = r * util::sin(phi);
    axis[2] = z;
    initAxes[j] = axis;
  }

  auto& load = _sLattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    auto& block = _seedsD.getBlock(iC);
    auto& seeds = block.template getField<fields::PHYS_R>();
    auto& seedsVorticity = block.template getField<descriptors::VORTICITY>();
    auto& seedsSigma = block.template getField<descriptors::SIGMA>();
    auto& seedsAxisDirection = block.template getField<VortexMethod::AXIS_DIRECTION>();

    for (int j = 0; j < _nSeeds; j++) {
      seeds.set(j, initPositions[j]);
      seedsVorticity.set(j, initVorticities[j]);
      seedsSigma.set(j, initSigmas[j]);
      seedsAxisDirection.set(j, initAxes[j]);
    }

    block.setProcessingContext(ProcessingContext::Simulation);
  }
}

template <typename T, typename DESCRIPTOR>
void VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR>::setVelocityProfile(
  std::shared_ptr<AnalyticalF3D<T,T>> velocityProfileF)
{
  _velocityProfileF = velocityProfileF;
  //_sLattice.template defineField<VortexMethod::U_PROFILE>(*_inletLatticeI, *_velocityProfileF);
  //_sLattice.template setProcessingContext<Array<VortexMethod::U_PROFILE>>(ProcessingContext::Simulation);
}

template <typename T, typename DESCRIPTOR>
void VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR>::setIntensityProfile(std::shared_ptr<AnalyticalF3D<T,T>> intensityProfileF)
{
  _intensityProfileF = intensityProfileF;
}

template <typename T, typename DESCRIPTOR>
void VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR>::setConvectionU(T uConvection)
{
  _uConvection = uConvection;
}

template <typename T, typename DESCRIPTOR>
void VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR>::apply(std::size_t iT)
{
  if (iT == 0) {
    generateSeeds();

    _convectionO.template setParameter<VortexMethod::MAX_R>(_inletPhysI->getMax());
    _convectionO.template setParameter<VortexMethod::MIN_R>(_inletPhysI->getMin());
    _convectionO.template setParameter<VortexMethod::L_MIN>(_lMin);
    _convectionO.template setParameter<VortexMethod::L_CHAR>(_lChar);
    _convectionO.template setParameter<VortexMethod::U_CONVECTION>(_uConvection);

    _convectionO.template setParameter<VortexMethod::SEED_COUNT>(_nSeeds);
    _convectionO.template setParameter<VortexMethod::INLET_AREA>(_inletArea);

    _convectionO.template setParameter<fields::converter::PHYS_DELTA_T>(_converter.getPhysDeltaT());
    _convectionO.template setParameter<VortexMethod::INTENSITY>(0.05);
    _convectionO.template setParameter<VortexMethod::REFERENCE_HEIGHT>(68.5);
  } else {
    _convectionO.template setParameter<VortexMethod::U_CONVECTION>(_uConvection);
  }

  _convectionO.template setParameter<VortexMethod::SEED>(_generator());
  _convectionO.apply();
}


struct VortexMethodPreProcessor {
  using parameters = meta::list<
    VortexMethod::CONVERSION_FACTOR_LENGTH,
    VortexMethod::CONVERSION_FACTOR_VELOCITY,
    VortexMethod::SEED_COUNT,
    VortexMethod::SEEDS,
    VortexMethod::SEEDS_VORTICITY,
    VortexMethod::SEEDS_SIGMA,
    VortexMethod::SEEDS_AXIS_DIRECTION,
    VortexMethod::U_PROFILE_SCALE,
    VortexMethod::CUTOFF_RATIO
  >;

  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;

    const auto seeds = parameters.template get<VortexMethod::SEEDS>();
    const auto seedsVorticity = parameters.template get<VortexMethod::SEEDS_VORTICITY>();
    const auto seedsAxisDirection = parameters.template get<VortexMethod::SEEDS_AXIS_DIRECTION>();
    const auto seedsSigma = parameters.template get<VortexMethod::SEEDS_SIGMA>();
    const V cutoffRatio2 = parameters.template get<VortexMethod::CUTOFF_RATIO>();
    const V scale = parameters.template get<VortexMethod::U_PROFILE_SCALE>();

    Vector<V,DESCRIPTOR::d> x = cell.template getField<fields::PHYS_R>();
    Vector<V,DESCRIPTOR::d> output = cell.template getField<VortexMethod::U_PROFILE>();
    const V conversionFactorVelocity = parameters.template get<VortexMethod::CONVERSION_FACTOR_VELOCITY>();
    output /= conversionFactorVelocity;

    Vector<V,3> uVortex {0, 0, 0};
    for (std::size_t j=0; j < parameters.template get<VortexMethod::SEED_COUNT>(); ++j) {
      const V sig = seedsSigma[j];
      const V sig2 = sig * sig;
      const V maxR2 = cutoffRatio2 * sig2;

      V dx = seeds[0][j] - x[0];
      V dx2 = dx * dx;
      if (dx2 > maxR2) {
        continue;
      }

      V dy = seeds[1][j] - x[1];
      V dy2 = dy * dy;
      if (dx2 + dy2 > maxR2) {
        continue;
      }

      V dz = seeds[2][j] - x[2];
      V norm2 = dx2 + dy2 + dz * dz;
      if (norm2 > maxR2) {
        continue;
      }

      Vector<V,3> diffX{dx, dy, dz};
      Vector<V,3> axisDirection(seedsAxisDirection, j);
      Vector<V,3> crossP = crossProduct3D(diffX, axisDirection);

      V expTerm = util::exp(V{-0.5} * norm2 / sig2);
      uVortex += (V{0.5} / std::numbers::pi_v<V> * seedsVorticity[j] * crossP / norm2 * (V{1} - expTerm) * expTerm)
               / conversionFactorVelocity;
    }

    output += uVortex;
    output *= scale;
    cell.template setField<momenta::FixedVelocityMomentumGeneric::VELOCITY>(output);
  }

};

}

}  // namespace optimized
#endif
