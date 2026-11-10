/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Adrian Kummerlaender, Fedor Bukreev
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

#include <olb.h>
#include "boundary/postprocessor/turbulentWallModelPostProcessorOptimized.h"

using namespace olb;

using T = float;
using DESCRIPTOR = descriptors::D3Q19<
  descriptors::POROSITY,
  descriptors::VELOCITY,

  fields::fsi::ELEMENT_TAG,
  fields::fsi::ELEMENT_TORQUE,

  collision::LES::SMAGORINSKY
>;

using FSI_DESCRIPTOR = descriptors::fsi::REFERENCE_POROSITY_ELEMENT_WITH_WALL_MODEL<DESCRIPTOR::d>;

using BulkDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
    momenta::BulkMomentum,
    momenta::BulkStress,
    momenta::DefineToNEq
  >,
  equilibria::ThirdOrder,
  collision::ParameterFromCell<collision::LES::SMAGORINSKY,
                               collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB>>
>;

namespace olb::stage {

struct Slice { };

};

void prepareGeometry(const UnitConverter<T,DESCRIPTOR>& converter,
                     SuperGeometry<T,3>& sGeometry,
                     Vector<T,3> extent,
                     Vector<T,3> rotorPos,
                     T rotorL,
                     T rotorR)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  sGeometry.setWriteIncrementalVTK(false);

  Vector<T,3> origin{};

  sGeometry.rename(0, 3);
  sGeometry.rename(3, 1, {1,1,1});

  {
    std::unique_ptr<IndicatorF3D<T>> fsiRegionI(new IndicatorSphere3D<T>(rotorPos, 1.5*rotorR));
    sGeometry.rename(1, 2, fsiRegionI);
  }

  // Set material number for inflow
  auto inflowExtent = extent;
  inflowExtent[2] = 2.*converter.getPhysDeltaX();
  origin[2] = -converter.getPhysDeltaX();
  IndicatorCuboid3D<T> inflow(inflowExtent, origin);
  sGeometry.rename(3, 4, 1, inflow);

  // Set material number for outflow
  origin[2] = extent[2]-converter.getPhysDeltaX();
  IndicatorCuboid3D<T> outflow(extent, origin);
  sGeometry.rename(3, 5, 1, outflow);

  sGeometry.checkForErrors();
  sGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(SuperLattice<T, DESCRIPTOR>& sLattice,
                    const UnitConverter<T, DESCRIPTOR>& converter,
                    SuperGeometry<T,3>& sGeometry,
                    Vector<T,3> extent)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  sLattice.defineDynamics<BulkDynamics>(sGeometry, 1);

  WallModelParameters<T> wallModelParameters;
  wallModelParameters.bodyForce = false;
  wallModelParameters.rhoMethod = 0;
  wallModelParameters.fNeqMethod = 0;
  wallModelParameters.samplingCellDistance = 2.5;
  wallModelParameters.interpolateSampleVelocity = true;
  wallModelParameters.useVanDriest = false;
  wallModelParameters.wallFunctionProfile = 1;
  wallModelParameters.movingWall = true;
  wallModelParameters.averageVelocity = false;

  //setTurbulentWallModelDynamics(sLattice, sGeometry.getMaterialIndicator(2), wallModelParameters);
  sLattice.template defineDynamics<WMHRRdynamics>(sGeometry.getMaterialIndicator(2));

  boundary::set<boundary::FullSlip>(sLattice, sGeometry, 3);
  boundary::set<boundary::LocalVelocity>(sLattice, sGeometry, 5);
  boundary::set<boundary::InterpolatedPressure>(sLattice, sGeometry, 4);

  sLattice.addPostProcessor(sGeometry.getMaterialIndicator(2),
                            meta::id<optimized::TurbulentWallModelPostProcessor<true,true,false,1,true>>{});
  sLattice.addPostProcessor(sGeometry.getMaterialIndicator(2),
                            meta::id<optimized::TurbulentWallModelPorousFneqFDMPostProcessor<0>>{});

  sLattice.setParameter<descriptors::SAMPLING_DISTANCE>(2.5);
  sLattice.setParameter<descriptors::OMEGA>(omega);
  sLattice.setParameter<collision::LES::SMAGORINSKY>(0.18);

  {
    AnalyticalConst3D<T,T> smagorinskyF(0.18);
    sLattice.defineField<collision::LES::SMAGORINSKY>(sGeometry.getMaterialIndicator({0,1,2,3,4,5}),
                                                      smagorinskyF);
  }
  {
    Vector<T,3> fringeOrigin(-converter.getPhysDeltaX(),
                             -converter.getPhysDeltaX(),
                             0);
    Vector<T,3> fringeExtent(extent[0] + 2*converter.getPhysDeltaX(),
                             extent[1] + 2*converter.getPhysDeltaX(),
                             30*converter.getPhysDeltaX());
    IndicatorCuboid3D<T> out(fringeExtent, fringeOrigin);
    SuperIndicatorFfromIndicatorF3D<T> outI(out, sGeometry);

    AnalyticalFfromCallableF<3,T,T> smagorinskyF([&](Vector<T,3> physR) -> Vector<T,1> {
      return 0.18 + ((30*converter.getPhysDeltaX() - physR[2]) / (30*converter.getPhysDeltaX())) * (3-0.18);
    });
    sLattice.defineField<collision::LES::SMAGORINSKY>(outI, smagorinskyF);
  }

  {
    AnalyticalConst3D<T,T> porosityF(1);
    sLattice.defineField<descriptors::POROSITY>(sGeometry.getMaterialIndicator({0,1,2,3,4,5}), porosityF);
    sLattice.defineField<descriptors::WMPOROSITY>(sGeometry.getMaterialIndicator({0,1,2,3,4,5}), porosityF);
  }

  AnalyticalConst3D<T,T> rhoF(1);
  AnalyticalConst3D<T,T> uF(0,0,0);
  sLattice.iniEquilibrium(sGeometry.getMaterialIndicator({1,2,3,4,5}), rhoF, uF);
  sLattice.defineRhoU(sGeometry.getMaterialIndicator({1,2,3,4,5}), rhoF, uF);

  sLattice.defineField<fields::moments::RHO>(sGeometry.getMaterialIndicator({1,2,3,4,5}), rhoF);
  sLattice.defineField<fields::moments::U>(sGeometry.getMaterialIndicator({1,2,3,4,5}), uF);
  sLattice.defineField<descriptors::Y1>(sGeometry.getMaterialIndicator(2), uF);

  sLattice.initialize();
  sLattice.writeSummary();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void getResults(SuperLattice<T,DESCRIPTOR>& sLattice,
                const UnitConverter<T,DESCRIPTOR>& converter,
                SuperGeometry<T,3>& sGeometry,
                util::Timer<T>& timer,
                std::size_t iT,
                std::size_t statT,
                std::size_t outputT,
                bool writeOutputVTK)
{
  OstreamManager clout(std::cout, "getResults");

  if (iT == 0) {
    SuperVTMwriter3D<T> vtmWriter("flow");
    vtmWriter.createMasterFile();
  }

  if (iT % statT == 0) {
    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    if (std::isnan(sLattice.getStatistics().getAverageRho())) {
      std::exit(-1);
    }
  }

  if (writeOutputVTK && iT % outputT == 0) {
    sLattice.executePostProcessors(stage::Evaluation{});
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    sLattice.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtkWriter("flow");
      SuperLatticePhysVelocity3D<T,DESCRIPTOR> velocityF(sLattice, converter);
      vtkWriter.addFunctor(velocityF);
      SuperLatticeField3D<T,DESCRIPTOR,descriptors::POROSITY> porosityF(sLattice);
      vtkWriter.addFunctor(porosityF);
      task(vtkWriter, iT);
    });
  }
}

Vector<T,3*3> calculateRotationMatrix(Vector<T,3> r, T a) {
  return Vector<T,3*3>{r[0]*r[0]*(1-util::cos(a)) +      util::cos(a), r[0]*r[1]*(1-util::cos(a)) - r[2]*util::sin(a), r[0]*r[2]*(1-util::cos(a)) + r[1]*util::sin(a),
                       r[0]*r[1]*(1-util::cos(a)) + r[2]*util::sin(a), r[1]*r[1]*(1-util::cos(a)) +      util::cos(a), r[1]*r[2]*(1-util::cos(a)) - r[0]*util::sin(a),
                       r[0]*r[2]*(1-util::cos(a)) - r[1]*util::sin(a), r[1]*r[2]*(1-util::cos(a)) + r[0]*util::sin(a), r[2]*r[2]*(1-util::cos(a)) +      util::cos(a)};
}

int main(int argc, char* argv[]) {
  initialize(&argc, &argv);
  OstreamManager clout(std::cout, "main");

  CLIreader args(argc, argv);

  const int resolution = args.getValueOrFallback("--resolution", 10);
  const int downsample = args.getValueOrFallback("--downsample", 1);

  const bool writeOutputVTK = !args.contains("--no-vtk");
  const bool writeOutputCDF = !args.contains("--no-cdf");

  const std::string rotorPath = args.getValueOrFallback<std::string>("--rotor", "rotor.stl");
  const T rotorL = args.getValueOrFallback("--rotorL", 0.08);
  const T rotorR = args.getValueOrFallback("--rotorR", 0.09);

  const T lengthX = (2*rotorR) * 3*args.getValueOrFallback<T>("--width", 1);
  const T lengthY = (2*rotorR) * 3*args.getValueOrFallback<T>("--height", 1);
  const T lengthZ = args.getValueOrFallback<T>("--length", 1.2);

  const Vector<T,3> rotorPos{lengthX/2,
                             lengthY/2,
                             lengthZ - T{0.3}};

  const T yaw = args.getValueOrFallback<T>("--yaw", 0.0);

  const T baseRotationX = args.getValueOrFallback<T>("--baseRotationX", 0.0);
  const T baseRotationY = args.getValueOrFallback<T>("--baseRotationY", yaw * std::numbers::pi_v<T> / 180);
  const T baseRotationZ = args.getValueOrFallback<T>("--baseRotationZ", 0.0);

  const T maxInflowU = args.getValueOrFallback("--maxInflowU", 0.56);

  const T tsr = args.getValueOrFallback<T>("--tsr", 3.03); // tip speed ratio
  const T rotationU = (tsr * maxInflowU) / rotorR;

  const T rampUpT = args.getValueOrFallback("--startT", 2.0);

  const int samplesPerRotation = args.getValueOrFallback("--samplesPerRotation", 3);
  const T rotationsToSample = args.getValueOrFallback<T>("--rotationsToSample", 10);
  const T maxPhysT = args.getValueOrFallback<T>(
    "--maxT", 2*rampUpT + rotationsToSample*(2*std::numbers::pi_v<T>/rotationU));



  singleton::directories().setOutputDir(
    args.getValueOrFallback<std::string>("--output", "./tmp/"));

  const UnitConverterFromResolutionAndLatticeVelocity<T,DESCRIPTOR> converter(
    int   {resolution},   // resolution: number of voxels per charPhysL
    (T)   0.1,            // charLatticeVelocity
    (T)   0.023,          // charPhysLength: reference length of simulation geometry
    (T)   tsr * maxInflowU, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1.0035e-6,      // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.23            // physDensity: physical density in __kg / m^3__
  );
  converter.print();

  const std::size_t startT = converter.getLatticeTime(rampUpT);

  const T rotationPeriod = 2.0 * std::numbers::pi_v<T> / rotationU;
  const T physOutputT = util::max(converter.getPhysDeltaT(),
                                  rotationPeriod / samplesPerRotation);
  const T physStatT = util::max(converter.getPhysDeltaT(),
                                rotationPeriod / (5*samplesPerRotation));

  const std::size_t statT = converter.getLatticeTime(args.getValueOrFallback("--statT", physStatT));
  const std::size_t outputT = converter.getLatticeTime(args.getValueOrFallback("--outputT", physOutputT));

  Vector<T,3> origin;
  Vector<T,3> extent(lengthX, lengthY, lengthZ);
  IndicatorCuboid3D<T> cuboid(extent, origin);

  CuboidDecomposition3D cuboidDecomposition(cuboid,
                                            converter.getPhysDeltaX(),
                                            1);
  //cuboidDecomposition.splitFractional(0, 2, {0.98, 0.02});
  cuboidDecomposition.split(0, singleton::mpi().getSize());
  cuboidDecomposition.buildNeighborCache(4);

  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  SuperGeometry<T,3> sGeometry(cuboidDecomposition, loadBalancer, 4);
  prepareGeometry(converter, sGeometry, extent, rotorPos, rotorL, rotorR);

  SuperLattice<T,DESCRIPTOR> sLattice(converter, sGeometry);
  prepareLattice(sLattice, converter, sGeometry, extent);

  SuperPorousElementEmbeddingO fsiEmbeddingO(sLattice);
  fsiEmbeddingO.registerElementType<ArbitrarilyRotateableReferenceLatticeWithWallModelPorosityF>(
    sGeometry.getMaterialIndicator(2));

  SuperPorousElementReductionO<T,DESCRIPTOR,fields::fsi::ELEMENT_TORQUE> fsiReductionO(
    sLattice,
    sGeometry.getMaterialIndicator(2));
  fsiReductionO.resize(1);
  fsiReductionO.addCollectionO(meta::id<CollectPorousBoundaryTorqueO>{});

  auto base = util::calculateRotationMatrix<T,DESCRIPTOR::d>({baseRotationX, baseRotationY, baseRotationZ});
  auto axis = util::matrixVectorProduct(base, Vector<T,3>{0.0,0.0,1.0});
  auto turn = calculateRotationMatrix(axis, 0.0);
  auto rotation = util::matrixMatrixProduct<T,DESCRIPTOR::d>(turn, base);

  std::unique_ptr<BlockD<T,FSI_DESCRIPTOR>> rotorD;
  if (fsiReductionO.rankDoesFSI()) {
    clout << "Set up rotor element" << std::endl;
    const T elementDeltaX = converter.getPhysDeltaX() / 2;
    std::shared_ptr<STLreader<T>> rotorI(
      new STLreader<T>(rotorPath, elementDeltaX, 0.001, RayMode::Robust));
    std::shared_ptr<IndicatorF3D<T>> rotorBoundingI(
      new IndicatorLayer3D<T>(*rotorI, converter.getPhysDeltaX()));
    Vector<std::size_t,DESCRIPTOR::d> rotorDim = (rotorBoundingI->getMax() - rotorBoundingI->getMin())
                                               / elementDeltaX + 1.5;
    clout << rotorI->getMin() << " " << rotorI->getMax() << std::endl;
    rotorD.reset(makeSharedBlockD<T,FSI_DESCRIPTOR>(
      loadBalancer, rotorDim).release());
    {
      rotorD->forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> latticeR) {
        bool inside{};
        Vector<T,DESCRIPTOR::d> physR = rotorBoundingI->getMin()
                                      + latticeR * elementDeltaX;
        rotorI->operator()(&inside, physR.data());
        if (inside) {
          rotorD->setField<descriptors::POROSITY>(latticeR, 0);
        } else {
          rotorD->setField<descriptors::POROSITY>(latticeR, 1);
        }
      });
      rotorD->forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> latticeR) {
        if (rotorD->getField<descriptors::POROSITY>(latticeR) == 1) {
          Vector<T,3> normal{};
          for(int jPop = 1; jPop < descriptors::D3Q27<>::q; ++jPop){
            auto latticeR2 = latticeR + descriptors::c<descriptors::D3Q27<>>(jPop);
            if (rotorD->isInside(latticeR2) && rotorD->getField<descriptors::POROSITY>(latticeR2) > 0) {
              normal += descriptors::c<descriptors::D3Q27<>>(jPop);
            }
          }
          if (util::norm<3>(normal) > 0) {
            normal /= util::norm<3>(normal);
            auto y1 = (elementDeltaX / converter.getPhysDeltaX()) * normal;
            rotorD->setField<descriptors::Y1>(latticeR, y1);
          } else {
            for(int jPop = 1; jPop < descriptors::D3Q27<>::q; ++jPop){
              auto latticeR2 = latticeR + 2*descriptors::c<descriptors::D3Q27<>>(jPop);
              if (rotorD->isInside(latticeR2) && rotorD->getField<descriptors::POROSITY>(latticeR2) > 0) {
                normal += descriptors::c<descriptors::D3Q27<>>(jPop);
              }
            }
            if (util::norm<3>(normal) > 0) {
              normal /= util::norm<3>(normal);
              auto y1 = ((2*elementDeltaX) / converter.getPhysDeltaX()) * normal;
              rotorD->setField<descriptors::Y1>(latticeR, y1);
            } else {
              for(int jPop = 1; jPop < descriptors::D3Q27<>::q; ++jPop){
                auto latticeR2 = latticeR + 3*descriptors::c<descriptors::D3Q27<>>(jPop);
                if (rotorD->isInside(latticeR2) && rotorD->getField<descriptors::POROSITY>(latticeR2) > 0) {
                  normal += descriptors::c<descriptors::D3Q27<>>(jPop);
                }
              }
              if (util::norm<3>(normal) > 0) {
                normal /= util::norm<3>(normal);
                auto y1 = ((3*elementDeltaX) / converter.getPhysDeltaX()) * normal;
                rotorD->setField<descriptors::Y1>(latticeR, y1);
              } else {
                for(int jPop = 1; jPop < descriptors::D3Q27<>::q; ++jPop){
                  auto latticeR2 = latticeR + 4*descriptors::c<descriptors::D3Q27<>>(jPop);
                  if (rotorD->isInside(latticeR2) && rotorD->getField<descriptors::POROSITY>(latticeR2) > 0) {
                    normal += descriptors::c<descriptors::D3Q27<>>(jPop);
                  }
                }
                if (util::norm<3>(normal) > 0) {
                  normal /= util::norm<3>(normal);
                  auto y1 = ((4*elementDeltaX) / converter.getPhysDeltaX()) * normal;
                  rotorD->setField<descriptors::Y1>(latticeR, y1);
                }
              }
            }
          }
        }
      });
    }
    rotorD->setProcessingContext(ProcessingContext::Simulation);

    clout << rotorBoundingI->getMin() << std::endl;

    {
      Vector lower{rotorBoundingI->getMin()[0]+rotorPos[0],
                   rotorBoundingI->getMin()[1]+rotorPos[1],
                   rotorPos[2]};
      auto params = makeParametersD<T,descriptors::SPATIAL_DESCRIPTOR<3>>(
        fields::fsi::ELEMENT_TAG{},   1,
        fields::fsi::ELEMENT_LOWER{}, lower,
        fields::fsi::ELEMENT_PIVOT{}, rotorPos,
        fields::fsi::ELEMENT_ROTATION_MATRIX{}, rotation,
        fields::fsi::ELEMENT_U_ROTATION{}, 0,
        fields::fsi::ELEMENT_REFERENCE_DELTA_X{},  elementDeltaX,
        fields::fsi::ELEMENT_REFERENCE_EXTENT{},   rotorDim,
        fields::fsi::ELEMENT_REFERENCE_POROSITY{}, rotorD->template getField<descriptors::POROSITY>(),
        fields::fsi::ELEMENT_REFERENCE_Y1{},       rotorD->template getField<descriptors::Y1>()
      );
      fsiEmbeddingO.add(params);
    }
  }

  clout << "Initialized porosity element" << std::endl;

  sLattice.setParameter<fields::converter::PHYS_LENGTH>(
    1/converter.getConversionFactorLength());
  sLattice.setParameter<fields::converter::PHYS_VELOCITY>(
    1/converter.getConversionFactorVelocity());
  sLattice.setParameter<fields::converter::PHYS_DELTA_X>(
    converter.getPhysDeltaX());

  fsiEmbeddingO.initialize();

  Vector<T,DESCRIPTOR::d> volumeSamplesOrigin = 2*converter.getPhysDeltaX();
  const T volumeSamplesDeltaX = converter.getPhysDeltaX();
  Vector<std::size_t,DESCRIPTOR::d> volumeSamplesExtent = (extent-4*converter.getPhysDeltaX())
                                                        / volumeSamplesDeltaX;
  //clout << "volumeSamplesDeltaX " << volumeSamplesDeltaX << std::endl;
  //clout << "volumeSamplesOrigin " << volumeSamplesOrigin << std::endl;
  //clout << "volumeSamplesExtent " << volumeSamplesExtent << std::endl;
  //std::vector<Vector<T,DESCRIPTOR::d>> volumeSamplesR;
  //for (std::size_t iX=0; iX < volumeSamplesExtent[0]; ++iX) {
  //  for (std::size_t iY=0; iY < volumeSamplesExtent[1]; ++iY) {
  //    for (std::size_t iZ=0; iZ < volumeSamplesExtent[2]; ++iZ) {
  //      volumeSamplesR.emplace_back(volumeSamplesOrigin + Vector{iX,iY,iZ}*volumeSamplesDeltaX);
  //    }
  //  }
  //}
  //SuperLatticePointExtraction<T,DESCRIPTOR,functors::pointwise::PhysVelocityF> pointExtractionO(
  //  sLattice, volumeSamplesR);
  //pointExtractionO.template setParameter<fields::converter::PHYS_VELOCITY>(
  //  converter.getConversionFactorVelocity());

  if (singleton::mpi().getRank() == 0) {
    std::ofstream meta("tmp/setup.meta");
    meta << "resolution=" << resolution << "\n"
         << "downsample=" << downsample << "\n"
         << "rotorPath=" << rotorPath << "\n"
         << "rotorL=" << rotorL << "\n"
         << "rotorR=" << rotorR << "\n"
         << "domainExtent=" << lengthX << " " << lengthY << " " << lengthZ << "\n"
         << "yaw=" << yaw << "\n"
         << "baseRotation=" << baseRotationX << " " << baseRotationY << " " << baseRotationZ << "\n"
         << "maxPhysT=" << maxPhysT << "\n"
         << "rampUpT=" << rampUpT << "\n"
         << "maxInflowU=" << maxInflowU << "\n"
         << "tsr=" << tsr << "\n"
         << "rotationU=" << rotationU << "\n"
         << "samplesPerRotation=" << samplesPerRotation << "\n"
         << "latticeCharU=" << converter.getCharLatticeVelocity() << "\n"
         << "physDeltaX=" << converter.getPhysDeltaX() << "\n"
         << "physDeltaT=" << converter.getPhysDeltaT() << "\n"
         << "physViscosity=" << converter.getPhysViscosity() << "\n";
    meta.close();
  }

  util::Timer<T> timer(converter.getLatticeTime(maxPhysT),
                       sGeometry.getStatistics().getNvoxel());
  timer.start();

  T angle = 0.0;

  for (std::size_t iT = 0; iT <= converter.getLatticeTime(maxPhysT); ++iT) {
    fsiEmbeddingO.apply();

    const T rampUpFrac = iT <= startT
                       ? 1 - util::cos(0.5*std::numbers::pi_v<T> * converter.getPhysTime(iT) / rampUpT)
                       : 1;
    // Ribeiro et al. reference: a uniform inflow ramped smoothly up to maxInflowU.
    // During ramp-up the inlet velocity (material 5, LocalVelocity boundary) is
    // updated each step; afterwards the boundary keeps imposing the final value.
    if (iT <= startT) {
      AnalyticalConst3D<T,T> uInflow(0, 0, converter.getLatticeVelocity(-rampUpFrac * maxInflowU));
      sLattice.defineU(sGeometry.getMaterialIndicator(5), uInflow);
      sLattice.setProcessingContext(ProcessingContext::Simulation);
    }

    sLattice.collideAndStream();
    getResults(sLattice, converter, sGeometry, timer, iT, statT, outputT, writeOutputVTK);

    //if (writeOutputCDF && iT % outputT == 0) {
    //  sLattice.waitForBackgroundTasks(stage::Evaluation{});
    //  pointExtractionO.update();
    //  if (singleton::mpi().getRank() == 0) {
    //  }
    //}

    fsiReductionO.apply();

    Vector angularU{0.0,0.0,-rotationU};
    angularU *= rampUpFrac;

    if (fsiReductionO.rankDoesFSI()) {
      for (unsigned iElement=1; iElement <= fsiReductionO.getElementCount(); ++iElement) {
        angle += angularU[2] * converter.getPhysDeltaT();
        angle = util::fmod(angle, 2*std::numbers::pi_v<T>);

        auto turn = calculateRotationMatrix(axis, angle);
        auto rotation = util::matrixMatrixProduct<T,DESCRIPTOR::d>(turn, base);

        fsiEmbeddingO.setField<fields::fsi::ELEMENT_ROTATION_MATRIX>(iElement-1, rotation);
        fsiEmbeddingO.setField<fields::fsi::ELEMENT_U_ROTATION>(iElement-1, angularU);
      }

      fsiEmbeddingO.setProcessingContext<Array<fields::fsi::ELEMENT_ROTATION_MATRIX>>(
        ProcessingContext::Simulation);
      fsiEmbeddingO.setProcessingContext<Array<fields::fsi::ELEMENT_U_ROTATION>>(
        ProcessingContext::Simulation);
    }
  }

  sLattice.setProcessingContext(ProcessingContext::Evaluation);
  timer.stop();
  timer.printSummary();
}
