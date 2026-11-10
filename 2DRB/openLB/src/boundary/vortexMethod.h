/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Fedor Bukreev, Adrian Kummerlaender
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

#ifndef BOUNDARIES_VORTEX_METHOD_H
#define BOUNDARIES_VORTEX_METHOD_H

#include "utilities/functorPtr.h"
#include "functors/lattice/indicator/superIndicatorF3D.h"

namespace olb {

struct VortexMethodPreProcessor;

struct U_PROFILE : public descriptors::FIELD_BASE<0,1> { };

struct CONVERSION_FACTOR_LENGTH : public descriptors::FIELD_BASE<1> { };
struct CONVERSION_FACTOR_VELOCITY : public descriptors::FIELD_BASE<1> { };

struct SEEDS_COUNT : public descriptors::TYPED_FIELD_BASE<std::size_t,1> { };
struct SEEDS : public descriptors::TEMPLATE_FIELD_BASE<std::add_pointer_t,0,1> { };
struct SEEDS_VORTICITY : public descriptors::TEMPLATE_FIELD_BASE<std::add_pointer_t,1> { };
struct SEEDS_AXIS_DIRECTION : public descriptors::TEMPLATE_FIELD_BASE<std::add_pointer_t,0,1> { };
struct SEEDS_SIGMA : public descriptors::TEMPLATE_FIELD_BASE<std::add_pointer_t,1> { };

struct AXIS_DIRECTION : public descriptors::FIELD_BASE<0,1> { };

namespace stage {

struct VortexMethod { };

}

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
  T _lChar;
  T _uConvection;

  Vector<T,3> _axisDirection;

  std::vector<T> _AiC;
  std::vector<int> _NiC;

  using SEED_DESCRIPTOR = descriptors::D3<descriptors::LOCATION,
                                          descriptors::VORTICITY,
                                          descriptors::SIGMA,
                                          AXIS_DIRECTION>;
  SuperD<T,SEED_DESCRIPTOR> _seedsD;

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
    T lChar,
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
  T lChar,
  Vector<T,3> axisDirection)
  : _inletLatticeI(std::move(inletLatticeI)),
    _inletPhysI(std::move(inletPhysI)),
    _converter(converter),
    _sLattice(sLattice),
    _nSeeds(nSeeds),
    _uConvection(uConvection),
    _sigmaConst(sigmaConst),
    _lChar(lChar),
    _axisDirection(axisDirection),
    _AiC(_sLattice.getLoadBalancer().size()),
    _NiC(_sLattice.getLoadBalancer().size()),
    _seedsD(sLattice.getLoadBalancer())
{
  OstreamManager clout(std::cout, "VortexMethod");

  {
    std::seed_seq seq{{0x9e3779b9u, 0xedb7d151u, 0x8212c45du, 0x33ca8075u,
                       0x11122233u, 0x44556677u, 0x8899aabbu, 0xccddeeffu}};
    _generator.seed(seq);
    _generator.discard(1e5);
  }

  sLattice.template addPostProcessor<stage::VortexMethod>(*_inletLatticeI,
                                                          meta::id<VortexMethodPreProcessor>{});
  _inletArea = 0.;

  auto& cGeometry = sLattice.getCuboidDecomposition();

  auto& load = _sLattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    _AiC[iC] = 0.;
    auto& cuboid = cGeometry.get(load.glob(iC));
    auto& block = _sLattice.getBlock(iC);
    bool partOfInlet = false;
    for (int iX=0; iX < block.getNx(); ++iX) {
      for (int iY=0; iY < block.getNy(); ++iY) {
        for (int iZ=0; iZ < block.getNz(); ++iZ) {
          int latticeR[3] {iX, iY, iZ};
          bool inInlet{};
          _inletLatticeI->getBlockIndicatorF(iC)(&inInlet, latticeR);
          if (inInlet) {
            partOfInlet = true;
            auto physR = cuboid.getPhysR(latticeR);
            block.get(iX,iY,iZ).template setField<descriptors::LOCATION>(physR);
            _AiC[iC] += util::pow(_converter.getPhysDeltaX(), 2.);
          }
        }
      }
    }
    if (partOfInlet) {
      block.template getData<Array<descriptors::LOCATION>>()
           .setProcessingContext(ProcessingContext::Simulation);
      _inletArea += _AiC[iC];
    }
  }

  #ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(_inletArea, MPI_SUM);
  #endif

  clout << "inletArea=" << _inletArea << std::endl;

  // initial seeding of random vortex points of number N
  for (int iC=0; iC < load.size(); ++iC) {
    _NiC[iC] = int(_nSeeds*_AiC[iC] / _inletArea);
    _seedsD.getBlock(iC).resize({_NiC[iC],1,1});
  }

  for (int iC=0; iC < load.size(); ++iC) {
    if (_NiC[iC] > 0) {
      auto& block = _sLattice.getBlock(iC);
      block.template setParameter<SEEDS_COUNT>(_NiC[iC]);
      auto& seedBlock = _seedsD.getBlock(iC);
      block.template setParameter<SEEDS>(seedBlock.template getField<descriptors::LOCATION>());
      block.template setParameter<SEEDS_VORTICITY>(seedBlock.template getField<descriptors::VORTICITY>());
      block.template setParameter<SEEDS_SIGMA>(seedBlock.template getField<descriptors::SIGMA>());
      block.template setParameter<SEEDS_AXIS_DIRECTION>(seedBlock.template getField<AXIS_DIRECTION>());
      block.template getField<U_PROFILE>();
    }
  }

  _sLattice.template setParameter<CONVERSION_FACTOR_LENGTH>(_converter.getConversionFactorLength());
  _sLattice.template setParameter<CONVERSION_FACTOR_VELOCITY>(_converter.getConversionFactorVelocity());
}

template <typename T, typename DESCRIPTOR>
void VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR>::generateSeeds()
{
  std::uniform_real_distribution<T> distribution(0, 1);
  std::uniform_real_distribution<T> sizeDistribution(0.8, 1.2);
  auto randomness = [&]() -> T {
    return distribution(_generator);
  };

  auto& load = _sLattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    auto& block = _seedsD.getBlock(iC);
    auto& seeds = block.template getField<descriptors::LOCATION>();
    auto& seedsVorticity = block.template getField<descriptors::VORTICITY>();
    auto& seedsSigma = block.template getField<descriptors::SIGMA>();
    auto& seedsAxisDirection = block.template getField<AXIS_DIRECTION>();

    if (_NiC[iC] > 0) {
      #ifdef PARALLEL_MODE_OMP
      #pragma omp parallel for schedule(static)
      #endif
      for (int j = 0; j < _NiC[iC]; j++) {
        auto inlet = _inletPhysI->getSample(randomness);
        seeds.set(j, inlet);

        Vector<T,3> velocity;
        _velocityProfileF->operator()(velocity.data(), inlet.data());
        Vector<T,1> intensity;
        _intensityProfileF->operator()(intensity.data(), inlet.data());

        T normVel = olb::norm(velocity);
        T kinE = 3./2. * util::pow(normVel * intensity[0], 2. );
        T sigmaI = _sigmaConst;

        if (_sigmaConst == 0 && normVel > 1e-6) {
          T mixingLength = 0.4 * util::abs(inlet[2]);
          mixingLength = util::max(mixingLength, _converter.getPhysDeltaX());
          mixingLength = util::min(mixingLength, _lChar);

          const T eps = util::pow(0.09, 3./4.) * util::pow(kinE, 3./2.) / mixingLength;
          sigmaI = 0.16*util::pow(kinE, 3./2.) / 2. / eps;
          sigmaI *= sizeDistribution(_generator);
        }

        T vorticity = 4. * util::sqrt(std::numbers::pi_v<T> * _AiC[iC] * kinE / 3. / _NiC[iC] / (2.*util::log(3.) - 3.*util::log(2.)) );

        seedsVorticity.set(j,vorticity);
        seedsSigma.set(j,sigmaI);

        Vector<T,3> axis{};
        T phi = randomness() * 2 * std::numbers::pi_v<T>;
        T z   = randomness() * 2 - 1;
        T r   = util::sqrt(1 - z*z);
        axis[0] = r * util::cos(phi);
        axis[1] = r * util::sin(phi);
        axis[2] = z;
        seedsAxisDirection.set(j, axis);
      }

      block.setProcessingContext(ProcessingContext::Simulation);
    }
  }
}

template <typename T, typename DESCRIPTOR>
void VortexMethodTurbulentVelocityBoundary<T,DESCRIPTOR>::setVelocityProfile(std::shared_ptr<AnalyticalF3D<T,T>> velocityProfileF)
{
  _velocityProfileF = velocityProfileF;
  _sLattice.template defineField<U_PROFILE>(*_inletLatticeI, *_velocityProfileF);
  _sLattice.template setProcessingContext<Array<U_PROFILE>>(ProcessingContext::Simulation);
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
  std::uniform_real_distribution<T> distribution(0, 1);
  std::uniform_real_distribution<T> sizeDistribution(0.8, 1.2);

  auto randomness = [&]() -> T {
    return distribution(_generator);
  };

  if (iT == 0) {
    generateSeeds();
  }

  auto& load = _sLattice.getLoadBalancer();
  for (int iC=0; iC < load.size(); ++iC) {
    if (_NiC[iC] > 0) {
      auto& block = _seedsD.getBlock(iC);
      auto& seeds = block.template getField<descriptors::LOCATION>();
      auto& seedsVorticity = block.template getField<descriptors::VORTICITY>();
      auto& seedsSigma = block.template getField<descriptors::SIGMA>();
      auto& seedsAxisDirection = block.template getField<AXIS_DIRECTION>();
      for (int i = 0; i < _NiC[iC]; ++i) {
        seeds[0][i] += _uConvection * _converter.getPhysDeltaT() * _axisDirection[0];
        if (seeds[0][i] > _inletPhysI->getMax()[0]) {
          Vector newPos = _inletPhysI->getSample(randomness);
          seeds[0][i] = _inletPhysI->getMin()[0];
          seeds[1][i] = newPos[1];
          seeds[2][i] = newPos[2];

          Vector<T,3> inlet = seeds.get(i);

          Vector<T,3> velocity;
          _velocityProfileF->operator()(velocity.data(), inlet.data());
          Vector<T,1> intensity;
          _intensityProfileF->operator()(intensity.data(), inlet.data());

          T normVel = olb::norm(velocity);
          T kinE = 3./2. * util::pow(normVel * intensity[0], 2. );
          T sigmaI = _sigmaConst;

          if (_sigmaConst == 0 && normVel > 1e-6) {
            T mixingLength = 0.4 * util::abs(inlet[2]);
            mixingLength = util::max(mixingLength, _converter.getPhysDeltaX());
            mixingLength = util::min(mixingLength, _lChar);

            const T eps = util::pow(0.09, 3./4.) * util::pow(kinE, 3./2.) / mixingLength;
            sigmaI = 0.16*util::pow(kinE, 3./2.) / 2. / eps;
            sigmaI *= sizeDistribution(_generator);
          }

          T vorticity = 4. * util::sqrt(std::numbers::pi_v<T> * _AiC[iC] * kinE / 3. / _NiC[iC] / (2.*util::log(3.) - 3.*util::log(2.)) );

          seedsVorticity.set(i,vorticity);
          seedsSigma.set(i,sigmaI);

          Vector<T,3> axis{};
          T phi = randomness() * 2 * std::numbers::pi_v<T>;
          T z   = randomness() * 2 - 1;
          T r   = util::sqrt(1 - z*z);
          axis[0] = r * util::cos(phi);
          axis[1] = r * util::sin(phi);
          axis[2] = z;
          seedsAxisDirection.set(i,axis);
        }
      }

      block.setProcessingContext(ProcessingContext::Simulation);
    }
  }

  {
    auto& load = _sLattice.getLoadBalancer();
    for (int iC=0; iC < load.size(); ++iC) {
      if (_NiC[iC] > 0) {
        auto& block = _sLattice.getBlock(iC);
        block.template setParameter<SEEDS_COUNT>(_NiC[iC]);
        auto& seedBlock = _seedsD.getBlock(iC);
        block.template setParameter<SEEDS>(seedBlock.template getField<descriptors::LOCATION>());
        block.template setParameter<SEEDS_VORTICITY>(seedBlock.template getField<descriptors::VORTICITY>());
        block.template setParameter<SEEDS_SIGMA>(seedBlock.template getField<descriptors::SIGMA>());
        block.template setParameter<SEEDS_AXIS_DIRECTION>(seedBlock.template getField<AXIS_DIRECTION>());
      }
    }
  }

  _sLattice.executePostProcessors(stage::VortexMethod{});
}


struct VortexMethodPreProcessor {
  using parameters = meta::list<
    CONVERSION_FACTOR_LENGTH,
    CONVERSION_FACTOR_VELOCITY,
    SEEDS_COUNT,
    SEEDS,
    SEEDS_VORTICITY,
    SEEDS_SIGMA,
    SEEDS_AXIS_DIRECTION
  >;

  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    using DESCRIPTOR = typename CELL::descriptor_t;

    const auto seeds = parameters.template get<SEEDS>();
    const auto seedsVorticity = parameters.template get<SEEDS_VORTICITY>();
    const auto seedsAxisDirection = parameters.template get<SEEDS_AXIS_DIRECTION>();
    const auto seedsSigma = parameters.template get<SEEDS_SIGMA>();

    Vector<V,DESCRIPTOR::d> x = cell.template getField<descriptors::LOCATION>();
    Vector<V,DESCRIPTOR::d> output = cell.template getField<U_PROFILE>();
    const V conversionFactorVelocity = parameters.template get<CONVERSION_FACTOR_VELOCITY>();
    output /= conversionFactorVelocity;

    // calculation of velocities from the placed vortexes
    Vector<V,3> uVortex {0, 0, 0};
    for (std::size_t j=0; j < parameters.template get<SEEDS_COUNT>(); ++j) {
      Vector<V,DESCRIPTOR::d> diffX = {0,0,0};
      for (int iD = 0; iD < DESCRIPTOR::d; ++iD) {
        diffX[iD] = seeds[iD][j] - x[iD];
      }
      Vector<V,3> axisDirection(seedsAxisDirection, j);
      Vector<V,3> crossP = crossProduct3D(diffX, axisDirection);
      V norm2 = diffX[0]*diffX[0] + diffX[1]*diffX[1] + diffX[2]*diffX[2];
      V expTerm = util::exp(V{-0.5}*norm2/(seedsSigma[j]*seedsSigma[j]));
      uVortex += (V{0.5}/std::numbers::pi_v<V> * seedsVorticity[j] * crossP / norm2 * (V{1} - expTerm)*expTerm)
               / conversionFactorVelocity;
    }

    output += uVortex;
    cell.defineU(output.data());
  }

};

}

#endif
