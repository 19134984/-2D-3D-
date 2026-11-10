/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2024 Fedor Bukreev
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
 *  GNU General Public License for mor details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this proquotiegram; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

/*
This showcase represents electrochemical saturation reaction of phosphate on a resolved porous CSH microparticle. The detailed description of the model and results is given in the following two publications:

Bukreev, F., Kummerländer, A., Jeßberger, J. et al. A hybrid Lattice-Boltzmann model for hydro-electrochemical modeling and sensitivity analysis of crystallization potential in nanoporous media. Part I: simulation model. Engineering with Computers (2025). https://doi.org/10.1007/s00366-025-02216-x

Bukreev, F., Kummerländer, A., Jeßberger, J. et al. A hybrid Lattice-Boltzmann model for hydro-electrochemical modeling and sensitivity analysis of crystallization potential in nanoporous media. Part II: application to the identification and quantification of influencing factors of phosphate saturation . Engineering with Computers (2025). https://doi.org/10.1007/s00366-025-02217-w

The porous rock geometry can be downloaded under: https://bwsyncandshare.kit.edu/s/bfgMMXWzRre9oyf
This geometry is part of the publication:
DyMAS: A direct multi-scale pore-level simulation approach. vol. Day 3 Tue, April 25, 2017 of SPE Western Regional Meeting; 2017. Available from: https://doi.org/10.2118/185720-MS
*/

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;

using MyCase = Case<
  Poisson, Lattice<double, descriptors::D3Q19<VELOCITY,SOURCE,NORMAL,CRYSTLAYER,NUCL,SOLUBILITY,MOMENTA_DENSITY>>,
  Concentration<0>, Lattice<double, descriptors::D3Q19<VELOCITY,SOURCE,NORMAL,CRYSTLAYER,NUCL,SOLUBILITY,MOMENTA_DENSITY>>,
  Concentration<1>, Lattice<double, descriptors::D3Q19<VELOCITY,SOURCE,NORMAL,CRYSTLAYER,NUCL,SOLUBILITY,MOMENTA_DENSITY>>,
  Concentration<2>, Lattice<double, descriptors::D3Q19<VELOCITY,SOURCE,NORMAL,CRYSTLAYER,NUCL,SOLUBILITY,MOMENTA_DENSITY>>,
  NavierStokes, Lattice<double, descriptors::D3Q19<FORCE,NORMAL,CRYSTLAYER,MOMENTA_DENSITY,MOMENTA_VELOCITY>>
>;

namespace olb::parameters {

struct IONS_RELAXATION_TIME : public descriptors::FIELD_BASE<1> { };
struct POISSON_RELAXATION_TIME : public descriptors::FIELD_BASE<1> { };
struct NSE_RELAXATION_TIME : public descriptors::FIELD_BASE<1> { };
struct IONS_PHYS_CHAR_VELOCITY : public descriptors::FIELD_BASE<1> { };
struct POISSON_PHYS_CHAR_VELOCITY : public descriptors::FIELD_BASE<1> { };
struct NSE_PHYS_CHAR_VELOCITY : public descriptors::FIELD_BASE<1> { };
struct DIFFUSION_C : public descriptors::FIELD_BASE<1> { };
struct DIFFUSION_P : public descriptors::FIELD_BASE<1> { };
struct DIFFUSION_H : public descriptors::FIELD_BASE<1> { };
struct VALENCEC : public descriptors::FIELD_BASE<1> { };
struct VALENCEP : public descriptors::FIELD_BASE<1> { };
struct VALENCEH : public descriptors::FIELD_BASE<1> { };
struct C_0 : public descriptors::FIELD_BASE<1> { };
struct PSI_BC : public descriptors::FIELD_BASE<1> { };
struct U_INFLOW : public descriptors::FIELD_BASE<1> { };
struct MFP : public descriptors::FIELD_BASE<1> { };
struct CRYSTCOEFF : public descriptors::FIELD_BASE<1> { };
struct CRYSTORDER : public descriptors::FIELD_BASE<1> { };
struct EQCONST : public descriptors::FIELD_BASE<1> { };
struct CRYSTMOLARMASS : public descriptors::FIELD_BASE<1> { };
struct CRYSTDENSITY : public descriptors::FIELD_BASE<1> { };
struct DESORPCOEFF : public descriptors::FIELD_BASE<1> { };
struct DESORPEQ : public descriptors::FIELD_BASE<1> { };
struct C0 : public descriptors::FIELD_BASE<1> { };
struct CONCENTRATION_H : public descriptors::FIELD_BASE<1> { };
struct CONCENTRATION_P : public descriptors::FIELD_BASE<1> { };
struct CONCENTRATION_C : public descriptors::FIELD_BASE<1> { };
struct A_ACTIVITY : public descriptors::FIELD_BASE<1> { };
struct B_ACTIVITY : public descriptors::FIELD_BASE<1> { };
struct ADIST_ACTIVITY : public descriptors::FIELD_BASE<1> { };
struct A_NUCL : public descriptors::FIELD_BASE<1> { };
struct MOLECVOL : public descriptors::FIELD_BASE<1> { };
struct BOUNDARY_DISTANCE : public descriptors::FIELD_BASE<1> { };
struct SCALING : public descriptors::FIELD_BASE<1> { };
struct SIMULATION_CUT : public descriptors::FIELD_BASE<1> { };
struct INPUT_VTI : descriptors::TYPED_FIELD_BASE<std::string,1> { };

}

// Rock geometry
std::shared_ptr<BlockVTIreader3D<MyCase::value_t, MyCase::value_t>> vtiReader; // to ensure that the data persists

std::shared_ptr<IndicatorBlockData3D<MyCase::value_t>>
generateIndicatorFromVTI(const std::string vtiFile, const std::string arrayName, MyCase::ParametersD& parameters)
{
  const MyCase::value_t sourceScale = parameters.get<parameters::SCALING>();

  vtiReader = std::make_shared<BlockVTIreader3D<MyCase::value_t, MyCase::value_t>>(vtiFile, arrayName);

  auto           cuboidSample     = vtiReader->getCuboid();
  MyCase::value_t              deltaRsample     = cuboidSample.getDeltaR() * sourceScale;
  Vector<int, 3> extentSample     = cuboidSample.getExtent();
  Vector<MyCase::value_t, 3>   originSamplePhys = cuboidSample.getOrigin() * sourceScale;
  Vector<MyCase::value_t, 3>   extentSamplePhys = {deltaRsample * MyCase::value_t(extentSample[0] + 0.5),
                                     deltaRsample * MyCase::value_t(extentSample[1] + 0.5),
                                     deltaRsample * MyCase::value_t(extentSample[2] + 0.5)};
  parameters.set<parameters::DOMAIN_EXTENT>(extentSamplePhys);

  std::shared_ptr<IndicatorBlockData3D<MyCase::value_t>> ind(
      new IndicatorBlockData3D<MyCase::value_t>(vtiReader->getBlockData(), extentSamplePhys,
                                  originSamplePhys, deltaRsample, false));

  return ind;
}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  const T L = parameters.get<parameters::PHYS_DELTA_X>();
  const T simulationCut = parameters.get<parameters::SIMULATION_CUT>();
  const std::string name = parameters.get<parameters::INPUT_VTI>();
  const int boundaryDistance = parameters.get<parameters::BOUNDARY_DISTANCE>();
  std::shared_ptr<IndicatorBlockData3D<T>> rock =
      generateIndicatorFromVTI(name, "Tiff Scalars", parameters);

  Vector<T,3> extend( simulationCut*(*rock).getMax()[0] - simulationCut*(*rock).getMin()[0] + 2.*boundaryDistance*L, simulationCut*(*rock).getMax()[1] - simulationCut*(*rock).getMin()[1] + 2.*boundaryDistance*L, simulationCut*(*rock).getMax()[2] - simulationCut*(*rock).getMin()[2] + 2.*boundaryDistance*L );
  Vector<T,3> origin( simulationCut*(*rock).getMin()[0] - boundaryDistance*L, simulationCut*(*rock).getMin()[1] - boundaryDistance*L, simulationCut*(*rock).getMin()[2] - boundaryDistance*L );
  IndicatorCuboid3D<T> cuboid( extend, origin );

  Mesh<T,MyCase::d> mesh(cuboid, L, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  const T L = parameters.get<parameters::PHYS_DELTA_X>();
  const T simulationCut = parameters.get<parameters::SIMULATION_CUT>();
  const std::string name = parameters.get<parameters::INPUT_VTI>();
  std::shared_ptr<IndicatorBlockData3D<T>> rock =
      generateIndicatorFromVTI(name, "Tiff Scalars", parameters);

  geometry.rename( 0, 3 );
  geometry.rename( 3, 1, {1,1,1} );

  Vector<T,3> extend( simulationCut*(*rock).getMax()[0] - simulationCut*(*rock).getMin()[0] - 2.*L, simulationCut*(*rock).getMax()[1] - simulationCut*(*rock).getMin()[1] - 2.*L, simulationCut*(*rock).getMax()[2] - simulationCut*(*rock).getMin()[2] - 2.*L );
  Vector<T,3> origin( simulationCut*(*rock).getMin()[0] + L, simulationCut*(*rock).getMin()[1] + L, simulationCut*(*rock).getMin()[2] + L );
  IndicatorCuboid3D<T> cuboid( extend, origin );

  geometry.rename( 1, 10, cuboid );
  geometry.rename( 10, 1, (*rock) );
  geometry.rename( 10, 2 );

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  // Removes all not needed boundary voxels inside the surface
  geometry.innerClean();
  geometry.checkForErrors();
  geometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLatticePoisson(MyCase& myCase) {
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice Poisson..." << std::endl;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<Poisson>;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(Poisson{});

  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const T relaxationTimePoisson = parameters.get<parameters::POISSON_RELAXATION_TIME>();
  const T charPhysVelocityPoisson = parameters.get<parameters::POISSON_PHYS_CHAR_VELOCITY>();
  const T dX = parameters.get<parameters::PHYS_DELTA_X>();
  const T simulationCut = parameters.get<parameters::SIMULATION_CUT>();

  // Set up a unit converter with the characteristic physical units
  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR>>(
    int{ static_cast<int>(simulationCut * extent[1] / dX) },
    relaxationTimePoisson,
    simulationCut*extent[1],
    charPhysVelocityPoisson,
    1.,
    1.
  );
  lattice.getUnitConverter().print();

  const T omega = lattice.getUnitConverter().getLatticeRelaxationFrequency();

  dynamics::set<SourcedAdvectionDiffusionBGKdynamics>(lattice, geometry.getMaterialIndicator({1, 3}));
  dynamics::set<EquilibriumBoundaryFirstOrder>(lattice, geometry, 2);
  boundary::set<boundary::AdvectionDiffusionDirichlet>(lattice, geometry, 3);

  lattice.setParameter<descriptors::OMEGA>(omega);

  clout << "Prepare Lattice Poisson ... OK" << std::endl;
}

template<size_t ID>
void prepareLatticeNernstPlanck(MyCase& myCase) {
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice Nernst-Planck ..." << std::endl;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<Concentration<ID>>;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(Concentration<ID>{});

  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const T relaxationTimeNPE = parameters.get<parameters::IONS_RELAXATION_TIME>();
  const T charPhysVelocityNPE = parameters.get<parameters::IONS_PHYS_CHAR_VELOCITY>();
  const T dX = parameters.get<parameters::PHYS_DELTA_X>();
  const T simulationCut = parameters.get<parameters::SIMULATION_CUT>();
  const T diffusion = (ID==0)*parameters.get<parameters::DIFFUSION_H>()
                    + (ID==1)*parameters.get<parameters::DIFFUSION_P>()
                    + (ID==2)*parameters.get<parameters::DIFFUSION_C>();

  // Set up a unit converter with the characteristic physical units
  lattice.template setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR>>(
    int{ static_cast<int>(simulationCut * extent[1] / dX) },
    relaxationTimeNPE,
    simulationCut*extent[1],
    charPhysVelocityNPE,
    diffusion,
    1.
  );
  lattice.getUnitConverter().print();

  const T omega = lattice.getUnitConverter().getLatticeRelaxationFrequency();
  dynamics::set<CrystalSourcedAdvectionDiffusionBGKdynamics>(lattice, geometry, 1);
  dynamics::set<BounceBack>(lattice, geometry, 2);
  if (ID != 2) {
    dynamics::set<CrystalSourcedAdvectionDiffusionBGKdynamics>(lattice, geometry, 3);
    boundary::set<boundary::AdvectionDiffusionDirichlet>(lattice, geometry, 3);
  } else {
    dynamics::set<ZeroDistributionDynamics>(lattice, geometry, 3);
    boundary::set<boundary::ZeroDistribution>(lattice, geometry, 3);
  }

  lattice.template setParameter<descriptors::OMEGA>(omega);

  clout << "Prepare Lattice Nernst-Planck ... OK" << std::endl;
}

void prepareLatticeNSE(MyCase& myCase) {
  OstreamManager clout( std::cout,"prepareLatticeNSE" );
  clout << "Prepare Lattice NSE..." << std::endl;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<NavierStokes>;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});

  const Vector extent = parameters.get<parameters::DOMAIN_EXTENT>();
  const T dX = parameters.get<parameters::PHYS_DELTA_X>();
  const T simulationCut = parameters.get<parameters::SIMULATION_CUT>();
  const T relaxationTimeNSE = parameters.get<parameters::NSE_RELAXATION_TIME>();
  const T charPhysVelocityNSE = parameters.get<parameters::NSE_PHYS_CHAR_VELOCITY>();
  const T charPhysViscosity = parameters.get<parameters::PHYS_CHAR_VISCOSITY>();
  const T charPhysDensity = parameters.get<parameters::PHYS_CHAR_DENSITY>();

  // Set up a unit converter with the characteristic physical units
  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR>>(
    int{ static_cast<int>(simulationCut * extent[1] / dX) },
    relaxationTimeNSE,
    simulationCut*extent[1],
    charPhysVelocityNSE,
    charPhysViscosity,
    charPhysDensity
  );
  lattice.getUnitConverter().print();

  const T omega = lattice.getUnitConverter().getLatticeRelaxationFrequency();

  dynamics::set<CrystalForcedBGKdynamics>(lattice, geometry.getMaterialIndicator({1,3}));
  boundary::set<boundary::FullSlip>(lattice, geometry, 2);
  boundary::set<boundary::InterpolatedVelocity>(lattice, geometry, 3);

  lattice.setParameter<descriptors::OMEGA>(omega);

  clout << "Prepare Lattice NSE... OK" << std::endl;
}

void prepareLatticeCoupling(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& latticePoisson = myCase.getLattice(Poisson{});
  auto& latticeH = myCase.getLattice(Concentration<0>{});
  auto& latticePO4 = myCase.getLattice(Concentration<1>{});
  auto& latticeCa = myCase.getLattice(Concentration<2>{});
  auto& latticeNSE = myCase.getLattice(NavierStokes{});
  auto& converterPoisson = latticePoisson.getUnitConverter();
  auto& converterNernstPlanck = latticeH.getUnitConverter();
  auto& converterNSE = latticeNSE.getUnitConverter();
  const T valenceH = parameters.get<parameters::VALENCEH>();
  const T valencePO4 = parameters.get<parameters::VALENCEP>();
  const T valenceCa = parameters.get<parameters::VALENCEC>();
  const T dielectricC = parameters.get<parameters::DIELECTRIC_CONSTANT>();
  const T diffH = parameters.get<parameters::DIFFUSION_H>();
  const T diffPO4 = parameters.get<parameters::DIFFUSION_P>();
  const T diffCa = parameters.get<parameters::DIFFUSION_C>();
  const T temperature = parameters.get<parameters::TEMPERATURE>();
  const T crystCoeff = parameters.get<parameters::CRYSTCOEFF>();
  const T crystOrder = parameters.get<parameters::CRYSTORDER>();
  const T Kocp = parameters.get<parameters::EQCONST>();
  const T crystMolarMass = parameters.get<parameters::CRYSTMOLARMASS>();
  const T crystDensity = parameters.get<parameters::CRYSTDENSITY>();
  const T desorpCoeff = parameters.get<parameters::DESORPCOEFF>();
  const T Kcsh = parameters.get<parameters::DESORPEQ>();
  const T Aactivity = parameters.get<parameters::A_ACTIVITY>();
  const T Bactivity = parameters.get<parameters::B_ACTIVITY>();
  const T aDistActivity = parameters.get<parameters::ADIST_ACTIVITY>();
  const T A_nucl = parameters.get<parameters::A_NUCL>();
  const T molecVolume = parameters.get<parameters::MOLECVOL>();

  const T npVelCoeffH = physConstants::elementaryCharge<T>() * diffH / physConstants::boltzmannConstant<T>() / temperature;
  const T npVelCoeffPO4 = physConstants::elementaryCharge<T>() * diffPO4 / physConstants::boltzmannConstant<T>() / temperature;
  const T npVelCoeffCa = physConstants::elementaryCharge<T>() * diffCa / physConstants::boltzmannConstant<T>() / temperature;
  const T sourceCoeff = 1./dielectricC * physConstants::faradayConstant<T>();
  const T forceCoeff = physConstants::faradayConstant<T>()/converterNSE.getPhysDensity();

  auto& coupling = myCase.setCouplingOperator(
    "NSPNPCrystal",
    NSPNPCrystalDynamicCoupling{},
    names::Concentration0{}, latticeH,
    names::Temperature{}, latticePoisson,
    names::Concentration1{}, latticePO4,
    names::Concentration2{}, latticeCa,
    names::NavierStokes{}, latticeNSE);
  coupling.setParameter<NSPNPCrystalDynamicCoupling::DX>(converterPoisson.getConversionFactorLength());
  coupling.setParameter<NSPNPCrystalDynamicCoupling::VALENCEH>(valenceH);
  coupling.setParameter<NSPNPCrystalDynamicCoupling::VALENCEP>(valencePO4);
  coupling.setParameter<NSPNPCrystalDynamicCoupling::VALENCEC>(valenceCa);
  coupling.setParameter<NSPNPCrystalDynamicCoupling::NPVELCOEFFH>(npVelCoeffH / converterNernstPlanck.getConversionFactorVelocity());
  coupling.setParameter<NSPNPCrystalDynamicCoupling::NPVELCOEFFP>(npVelCoeffPO4 / converterNernstPlanck.getConversionFactorVelocity());
  coupling.setParameter<NSPNPCrystalDynamicCoupling::NPVELCOEFFC>(npVelCoeffCa / converterNernstPlanck.getConversionFactorVelocity());
  coupling.setParameter<NSPNPCrystalDynamicCoupling::POISSONCOEFF>(sourceCoeff * converterPoisson.getConversionFactorTime());
  coupling.setParameter<NSPNPCrystalDynamicCoupling::FORCECOEFF>(forceCoeff * converterNSE.getConversionFactorMass() / converterNSE.getConversionFactorForce());
  coupling.setParameter<NSPNPCrystalDynamicCoupling::DTADE>(converterNernstPlanck.getConversionFactorTime());
  coupling.setParameter<NSPNPCrystalDynamicCoupling::DTNSE>(converterNSE.getConversionFactorTime());
  coupling.setParameter<NSPNPCrystalDynamicCoupling::OMEGA>(converterPoisson.getLatticeRelaxationFrequency());
  coupling.setParameter<NSPNPCrystalDynamicCoupling::TAUNSE>(1. / converterNSE.getLatticeRelaxationFrequency());
  coupling.setParameter<NSPNPCrystalDynamicCoupling::CRYSTCOEFF>(crystCoeff * converterNernstPlanck.getConversionFactorTime());
  coupling.setParameter<NSPNPCrystalDynamicCoupling::CRYSTORDER>(crystOrder);
  coupling.setParameter<NSPNPCrystalDynamicCoupling::EQCONST>(Kocp);
  coupling.setParameter<NSPNPCrystalDynamicCoupling::CRYSTMOLARMASS>(crystMolarMass);
  coupling.setParameter<NSPNPCrystalDynamicCoupling::CRYSTDENSITY>(crystDensity);
  coupling.setParameter<NSPNPCrystalDynamicCoupling::DESORPCOEFF>(desorpCoeff * converterNernstPlanck.getConversionFactorTime());
  coupling.setParameter<NSPNPCrystalDynamicCoupling::DESORPEQ>(Kcsh);
  coupling.setParameter<NSPNPCrystalDynamicCoupling::CMASS>(converterNSE.getPhysDensity());
  coupling.setParameter<NSPNPCrystalDynamicCoupling::A_ACTIVITY>(Aactivity);
  coupling.setParameter<NSPNPCrystalDynamicCoupling::B_ACTIVITY>(Bactivity);
  coupling.setParameter<NSPNPCrystalDynamicCoupling::ADIST_ACTIVITY>(aDistActivity);
  coupling.setParameter<NSPNPCrystalDynamicCoupling::A_NUCL>(A_nucl);
  coupling.setParameter<NSPNPCrystalDynamicCoupling::MOLECVOL>(molecVolume);
  coupling.setParameter<NSPNPCrystalDynamicCoupling::TEMPERATURE>(temperature);
  coupling.restrictTo(geometry.getMaterialIndicator({1}));
}

void setInitialValuesPoisson(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(Poisson{});
  auto& parameters = myCase.getParameters();

  const T psi0 = parameters.get<parameters::PSI_BC>();
  momenta::setElectricPotential(lattice, geometry.getMaterialIndicator(2), psi0);
  momenta::setElectricPotential(lattice, geometry.getMaterialIndicator({1, 3}), T(0));

  // Make the lattice ready for simulation
  lattice.initialize();
}

template<size_t ID>
void setInitialValuesNernstPlanck(MyCase& myCase) {
  OstreamManager clout( std::cout,"setInitialValuesNernstPlanck" );
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(Concentration<ID>{});

  const T rhoOut = (ID==0)*parameters.get<parameters::CONCENTRATION_H>()
                 + (ID==1)*parameters.get<parameters::CONCENTRATION_P>()
                 + (ID==2)*parameters.get<parameters::CONCENTRATION_C>();
  if( ID != 2) {
    fields::set<descriptors::CRYSTLAYER>(lattice, geometry.getMaterialIndicator({1, 3}), T(0));
    fields::set<descriptors::CRYSTLAYER>(lattice, geometry.getMaterialIndicator(2), T(1));
    momenta::setConcentration(lattice, geometry.getMaterialIndicator(3), rhoOut);
    momenta::setConcentration(lattice, geometry.getMaterialIndicator({1, 2}), T(0));

    AnalyticalConst3D<T,T> rhoF( rhoOut ); //C0
    AnalyticalConst3D<T,T> rho0( 0 );
    AnalyticalConst3D<T,T> uF( 0,0,0);
    auto bulkIndicator = geometry.getMaterialIndicator({3});
    lattice.iniEquilibrium( bulkIndicator, rhoF, uF );

    auto bulkIndicator1 = geometry.getMaterialIndicator({2, 1});
    lattice.iniEquilibrium( bulkIndicator1, rho0, uF );
  } else {
    fields::set<descriptors::CRYSTLAYER>(lattice, geometry.getMaterialIndicator({1, 3}), T(0));
    fields::set<descriptors::CRYSTLAYER>(lattice, geometry.getMaterialIndicator(2), T(1));
    momenta::setConcentration(lattice, geometry.getMaterialIndicator({1, 2, 3}), T(0));

    AnalyticalConst3D<T,T> rho0( 0 );
    AnalyticalConst3D<T,T> uF( 0,0,0);
    auto bulkIndicator1 = geometry.getMaterialIndicator({1, 2, 3});
    lattice.iniEquilibrium( bulkIndicator1, rho0, uF );
  }

  // Make the lattice ready for simulation
  lattice.initialize();
}

void setInitialValuesNSE(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(NavierStokes{});

  fields::set<descriptors::CRYSTLAYER>(lattice, geometry.getMaterialIndicator({1, 3}), T(0));
  fields::set<descriptors::CRYSTLAYER>(lattice, geometry.getMaterialIndicator(2), T(1));

  // Make the lattice ready for simulation
  lattice.initialize();
}

void setTemporalValues(MyCase& myCase,
                       std::size_t iT)
{
  OstreamManager clout( std::cout,"setTemporalValues" );
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& latticeH = myCase.getLattice(Concentration<0>{});
  auto& latticePO4 = myCase.getLattice(Concentration<1>{});
  auto& latticeCa = myCase.getLattice(Concentration<2>{});
  auto& latticeNSE = myCase.getLattice(NavierStokes{});
  auto& converterNernstPlanck = latticeH.getUnitConverter();

  const T uinflow = parameters.get<parameters::U_INFLOW>();
  T maxVelocity = converterNernstPlanck.getLatticeVelocity(uinflow);
  int iTmaxStart = 2000;

  if (int(iT) <= iTmaxStart && int(iT)%10 == 0.) {
    SinusStartScale<T, int> nSinusStartScale(iTmaxStart, 1.);
    int iTvec[1] = { int(iT) };
    T frac[1] = { T() };
    nSinusStartScale(frac, iTvec);
    AnalyticalConst3D<T,T> wantedVelocity(frac[0]*maxVelocity,T(),T());
    AnalyticalConst3D<T,T> wantedVelocityNSE(frac[0]*uinflow,T(),T());

    fields::set<VELOCITY>(latticeH, geometry.getMaterialIndicator(3), wantedVelocity);
    fields::set<VELOCITY>(latticePO4, geometry.getMaterialIndicator(3), wantedVelocity);
    fields::set<VELOCITY>(latticeCa, geometry.getMaterialIndicator(3), wantedVelocity);
    momenta::setVelocity(latticeNSE, geometry.getMaterialIndicator(3), wantedVelocityNSE);
    latticeH.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
    latticePO4.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
    latticeCa.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
    latticeNSE.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(ProcessingContext::Simulation);
  }
}


void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  OstreamManager clout( std::cout,"getResults" );
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<Poisson>;
  using DESCRIPTORNSE = MyCase::descriptor_t_of<NavierStokes>;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& latticePoisson = myCase.getLattice(Poisson{});
  auto& latticeH = myCase.getLattice(Concentration<0>{});
  auto& latticePO4 = myCase.getLattice(Concentration<1>{});
  auto& latticeCa = myCase.getLattice(Concentration<2>{});
  auto& latticeNSE = myCase.getLattice(NavierStokes{});
  auto& converter = latticeH.getUnitConverter();
  auto& converterNSE = latticeNSE.getUnitConverter();
  const T maxPhysT = parameters.get<parameters::PHYS_MAX_T>();
  const int vtkIter = parameters.get<parameters::PHYS_VTK_ITER_T>();
  const int statIter = parameters.get<parameters::PHYS_STAT_ITER_T>();

  static Gnuplot<double> csvWriter("solubility");

  if ( iT==0 ) {
    SuperVTMwriter3D<T> vtmWriterPoisson( "poisson" );
    SuperGeometryF3D<T> geom( geometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( latticeH );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( latticeH );

    vtmWriterPoisson.write( geom );
    vtmWriterPoisson.write( cuboid );
    vtmWriterPoisson.write( rank );
    vtmWriterPoisson.createMasterFile();

    SuperVTMwriter3D<T> vtmWriterH( "hIon" );
    vtmWriterH.createMasterFile();

    SuperVTMwriter3D<T> vtmWriterPO4( "po4Ion" );
    vtmWriterPO4.createMasterFile();

    SuperVTMwriter3D<T> vtmWriterCa( "caIon" );
    vtmWriterCa.createMasterFile();

    SuperVTMwriter3D<T> vtmWriterNSE( "nse" );
    vtmWriterNSE.createMasterFile();
  }

  // Get statistics
  if ( iT%statIter == 0) {

    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    latticeH.getStatistics().print( iT,converter.getPhysTime( iT ) );
    latticePO4.getStatistics().print( iT,converter.getPhysTime( iT ) );
    latticeCa.getStatistics().print( iT,converter.getPhysTime( iT ) );
    latticeNSE.getStatistics().print( iT,converterNSE.getPhysTime( iT ) );
    latticePoisson.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }

  // Writes the VTK files
  if ( iT%vtkIter == 0 ) {
    latticePoisson.setProcessingContext(ProcessingContext::Evaluation);
    latticePoisson.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriterPoisson( "poisson" );
      SuperLatticeDensity3D psi( latticePoisson);
      psi.getName() = "psi";
      vtmWriterPoisson.addFunctor(psi);
      task(vtmWriterPoisson, iT);
    });
  }

  if ( iT%vtkIter == 0 ) {
    latticeH.setProcessingContext(ProcessingContext::Evaluation);
    latticeH.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriterH( "hIon" );
      SuperLatticeDensity3D concentrationH( latticeH);
      concentrationH.getName() = "concentrationH";
      SuperLatticePhysField3D<T,DESCRIPTOR,VELOCITY> velH( latticeH, converter.getConversionFactorVelocity() );
      velH.getName() = "velH";
      vtmWriterH.addFunctor(concentrationH);
      vtmWriterH.addFunctor(velH);
      task(vtmWriterH, iT);
    });
  }

  if ( iT%vtkIter == 0 ) {
    latticePO4.setProcessingContext(ProcessingContext::Evaluation);
    latticePO4.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriterPO4( "po4Ion" );
      SuperLatticeDensity3D concentrationPO4( latticePO4);
      concentrationPO4.getName() = "concentrationPO4";
      SuperLatticePhysField3D<T,DESCRIPTOR,VELOCITY> velPO4( latticePO4, converter.getConversionFactorVelocity() );
      velPO4.getName() = "velPO4";
      vtmWriterPO4.addFunctor(concentrationPO4);
      vtmWriterPO4.addFunctor(velPO4);
      task(vtmWriterPO4, iT);
    });
  }

  if ( iT%vtkIter == 0 ) {
    latticeCa.setProcessingContext(ProcessingContext::Evaluation);
    latticeCa.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriterCa( "caIon" );
      SuperLatticeDensity3D concentrationCa( latticeCa);
      concentrationCa.getName() = "concentrationCa";
      SuperLatticePhysField3D<T,DESCRIPTOR,VELOCITY> velCa( latticeCa, converter.getConversionFactorVelocity() );
      velCa.getName() = "velCa";
      SuperLatticeField3D<T,DESCRIPTOR,NUCL> nucl( latticeCa);
      nucl.getName() = "nucl";
      SuperLatticeField3D<T,DESCRIPTOR,SOLUBILITY> sol( latticeCa);
      sol.getName() = "sol";
      SuperLatticeField3D<T,DESCRIPTOR,NORMAL> normal( latticeCa);
      normal.getName() = "normal";
      SuperLatticeField3D<T,DESCRIPTOR,CRYSTLAYER> surfFraction( latticeCa);
      surfFraction.getName() = "surfFraction";
      vtmWriterCa.addFunctor(concentrationCa);
      vtmWriterCa.addFunctor(velCa);
      vtmWriterCa.addFunctor(nucl);
      vtmWriterCa.addFunctor(sol);
      vtmWriterCa.addFunctor(normal);
      vtmWriterCa.addFunctor(surfFraction);
      task(vtmWriterCa, iT);
    });
  }

  if ( iT%vtkIter == 0 ) {
    latticeNSE.setProcessingContext(ProcessingContext::Evaluation);
    latticeNSE.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriterNSE( "nse" );
      SuperLatticePhysVelocity3D physVel( latticeNSE, converterNSE );
      physVel.getName() = "physVelNSE";
      SuperLatticePhysPressure3D press( latticeNSE, converterNSE );
      SuperLatticeField3D<T,DESCRIPTORNSE,NORMAL> normal( latticeNSE);
      normal.getName() = "normal";
      vtmWriterNSE.addFunctor(physVel);
      vtmWriterNSE.addFunctor(press);
      vtmWriterNSE.addFunctor(normal);
      task(vtmWriterNSE, iT);
    });
  }

  if ( iT%statIter == 0 && iT > 0 ) {
    latticeCa.setProcessingContext(ProcessingContext::Evaluation);
    SuperLatticeField3D<T,DESCRIPTOR,SOLUBILITY> sol( latticeCa);
    int dummy[4] = {0};
    T result[4];
    SuperAverage3D<T> (sol, geometry, 1).operator()(result, dummy);
    csvWriter.setData(converter.getPhysTime(iT), result[0], "Solubility development", "bottom right");
    clout << "average solubility: " << result[0] << std::endl;
  }

  if (iT == (converter.getLatticeTime(maxPhysT) -1)) {
    latticeCa.setProcessingContext(ProcessingContext::Evaluation);
    csvWriter.writePDF();
  }
}

void simulatePoisson(MyCase& myCase) {
  using T = MyCase::value_t;
  util::Timer<T> timer( 1000, myCase.getGeometry().getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( 100, 1.e-7 );
  timer.start();

  for ( std::size_t iT = 0; iT < 1000; ++iT ) {

    myCase.getLattice(Poisson{}).collideAndStream();

    if(converge.hasConverged()) {
      break;
    }

    converge.takeValue( myCase.getLattice(Poisson{}).getStatistics().getAverageRho(), false );
  }

  timer.stop();
}

void simulate(MyCase& myCase) {
  OstreamManager clout( std::cout,"simulate" );
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  const std::size_t iTmax = myCase.getLattice(Concentration<0>{}).getUnitConverter().getLatticeTime(
    parameters.get<parameters::PHYS_MAX_T>());

  util::Timer<T> timer( iTmax, myCase.getGeometry().getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( 1000, 1.e-9 );
  timer.start();

  for (std::size_t iT=0; iT < iTmax; ++iT) {
    myCase.getLattice(Poisson{}).collideAndStream();
    //simulatePoisson(myCase);

    setTemporalValues(myCase, iT);

    myCase.getOperator("NSPNPCrystal").apply();

    // === 6th Step: Collide and Stream Execution ===
    myCase.getLattice(NavierStokes{}).collideAndStream();
    myCase.getLattice(Concentration<0>{}).collideAndStream();
    myCase.getLattice(Concentration<1>{}).collideAndStream();
    myCase.getLattice(Concentration<2>{}).collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults(myCase, timer, iT);

    if(iT > 1.6e6){
      if(converge.hasConverged()) break;
      converge.takeValue( myCase.getLattice(Concentration<2>{}).getStatistics().getAverageRho(), false );
    }
  }

  timer.stop();
  timer.printSummary();
  singleton::pool().wait();
}

int main(int argc, char* argv[]) {
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    //electric params
    ///constants
    myCaseParameters.set<DIELECTRIC_CONSTANT>(6.95e-10); //  C/V/m
    myCaseParameters.set<VALENCEH>(1.);
    myCaseParameters.set<VALENCEP>(-3.);
    myCaseParameters.set<VALENCEC>(2.);
    ///TUNABLE
    myCaseParameters.set<PSI_BC>(0.01); // V

    //hydrodynamic BC
    ///TUNABLE
    myCaseParameters.set<U_INFLOW>(0.2);
    ///constants
    myCaseParameters.set<PHYS_CHAR_VISCOSITY>(1.e-6);

    //chemical params
    ///constants
    myCaseParameters.set<MFP>(1.26e-10); // mean free path of H2O
    myCaseParameters.set<DIFFUSION_H>(9.31e-9);
    myCaseParameters.set<DIFFUSION_P>(0.612e-9);
    myCaseParameters.set<DIFFUSION_C>(0.792e-9);
    myCaseParameters.set<CRYSTCOEFF>(3.67e-8);
    myCaseParameters.set<CRYSTORDER>(4);
    myCaseParameters.set<EQCONST>(util::pow(10.,-48.86));
    myCaseParameters.set<CRYSTMOLARMASS>(0.50231);
    myCaseParameters.set<CRYSTDENSITY>(3050.);
    myCaseParameters.set<DESORPCOEFF>(3.1e-11);
    myCaseParameters.set<DESORPEQ>(util::pow(10.,11.15));
    myCaseParameters.set<PHYS_CHAR_DENSITY>(1000.);
    ///TUNABLE
    myCaseParameters.set<C0>(1); // mol/m³ in DFG between 0.1 and 3.0
    myCaseParameters.set<CONCENTRATION_H>([] (MyCase::ParametersD& p) {
      return 6.*p.get<C_0>();
    });
    myCaseParameters.set<CONCENTRATION_P>([] (MyCase::ParametersD& p) {
      return 1.*p.get<C_0>();
    });
    myCaseParameters.set<CONCENTRATION_C>(0.);
    myCaseParameters.set<parameters::TEMPERATURE>(298.15);

    //activity coeffs
    ///constants
    myCaseParameters.set<A_ACTIVITY>(0.016102); //depends on temperature!!!!!!!
    myCaseParameters.set<B_ACTIVITY>(1.0391e8); //depends on temperature!!!!!!!
    myCaseParameters.set<ADIST_ACTIVITY>(2.33e-10);

    //nucleation params
    ///constants
    myCaseParameters.set<A_NUCL>(util::pow(10.,35.5));
    myCaseParameters.set<MOLECVOL>(2.735e-28);

    //geometrical params
    ///TUNABLE
    myCaseParameters.set<KNUDSEN>(0.09);
    myCaseParameters.set<PHYS_DELTA_X>([] (MyCase::ParametersD& p) {
      return p.get<MFP>() / p.get<KNUDSEN>();
    });
    myCaseParameters.set<BOUNDARY_DISTANCE>(10);
    myCaseParameters.set<parameters::SCALING>(2.5e-9*0.5);
    myCaseParameters.set<SIMULATION_CUT>(0.45);

    //simulation params
    myCaseParameters.set<PHYS_MAX_T>(2.1e-4);
    myCaseParameters.set<PHYS_VTK_ITER_T>(10000);
    myCaseParameters.set<PHYS_STAT_ITER_T>(1000);
    myCaseParameters.set<POISSON_RELAXATION_TIME>(1);
    myCaseParameters.set<POISSON_PHYS_CHAR_VELOCITY>(0.001);
    myCaseParameters.set<IONS_RELAXATION_TIME>(0.56);
    myCaseParameters.set<IONS_PHYS_CHAR_VELOCITY>(0.5);
    myCaseParameters.set<NSE_RELAXATION_TIME>(0.9);
    myCaseParameters.set<NSE_PHYS_CHAR_VELOCITY>(0.5);
    myCaseParameters.set<INPUT_VTI>(std::string("stone.vti"));
  }

  /// === Step 3: Create Mesh ===
  Mesh mesh = createMesh(myCaseParameters);

  /// === Step 4: Create Case ===
  MyCase myCase(myCaseParameters, mesh);

  /// === Step 5: Prepare Geometry ===
  prepareGeometry(myCase);

  /// === Step 6: Prepare Lattice ===
  prepareLatticePoisson(myCase);
  prepareLatticeNernstPlanck<0>(myCase);
  prepareLatticeNernstPlanck<1>(myCase);
  prepareLatticeNernstPlanck<2>(myCase);
  prepareLatticeNSE(myCase);

  prepareLatticeCoupling(myCase);

  /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
  setInitialValuesPoisson(myCase);
  setInitialValuesNernstPlanck<0>(myCase);
  setInitialValuesNernstPlanck<1>(myCase);
  setInitialValuesNernstPlanck<2>(myCase);
  setInitialValuesNSE(myCase);

  /// === Step 8: Simulate ===
  simulate(myCase);
}
