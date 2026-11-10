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
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

 /*
  * the expressions used in this file are used for validational part of the paper
  * Bukreev, F., Kummerländer, A., Jeßberger, J. et al. A hybrid Lattice-Boltzmann model for hydro-electrochemical modeling
  * and sensitivity analysis of crystallization potential in nanoporous media. Part I: simulation model. Engineering with
  * Computers (2025). https://doi.org/10.1007/s00366-025-02216-x
  *
  * The current example solves the Poisson-Boltzmann equation describing Gouy-Chapman double electric layer model for one dimension.
  * Example contains 3 lattices:
  * 1 for Poisson equation,
  * 1 for Nernst-Planck equation of cation,
  * 1 for Nernst-Planck equation of anion,
  */

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::names;

// Analytical Profile for electric potential
template <typename T, typename S, typename DESCRIPTOR>
class PotentialProfile1D : public AnalyticalF3D<T, S>
{
private:
  T y0, y1, psi0, Debye;
  UnitConverter<T, DESCRIPTOR> const &converter;
public:
  PotentialProfile1D(T y0_, T y1_, T psi0_, T Debye_, UnitConverter<T, DESCRIPTOR> const &converter_) :
  AnalyticalF3D<T, S>(1), y0(y0_), y1(y1_), psi0(psi0_), Debye(Debye_), converter(converter_)
  {
    this->getName() = "PotentialProfile1D";
  };

  bool operator()(T output[1], const S x[3])
  {
    T distY = x[1] - y0;
    output[0] = psi0 * util::exp( -distY/Debye );
    return true;
  };
};

// Analytical Profile for ion concentrations
template <typename T, typename S, typename DESCRIPTOR>
class ConcentrationProfile1D : public AnalyticalF3D<T, S>
{
private:
  T C0, valence, temperature;
  PotentialProfile1D<T,T,DESCRIPTOR> &potential;

public:
  ConcentrationProfile1D(T C0_, T valence_, T temperature_, PotentialProfile1D<T,T,DESCRIPTOR> &potential_) :
  AnalyticalF3D<T, S>(1), C0(C0_), valence(valence_), temperature(temperature_), potential(potential_)
  {
    this->getName() = "ConcentrationProfile1D";
  };

  bool operator()(T output[1], const S x[3])
  {
    T psi[1];
    potential(psi,x);
    output[0] = C0 * util::exp(-physConstants::elementaryCharge<T>()*valence*psi[0]/physConstants::boltzmannConstant<T>()/temperature);
    return true;
  };
};

using MyCase = Case<
  Poisson, Lattice<double, descriptors::D3Q19<VELOCITY,SOURCE>>,
  Concentration<0>, Lattice<double, descriptors::D3Q19<VELOCITY,SOURCE>>,
  Concentration<1>, Lattice<double, descriptors::D3Q19<VELOCITY,SOURCE>>
>;

namespace olb::parameters {

struct IONS_RELAXATION_TIME : public descriptors::FIELD_BASE<1> { };
struct POISSON_RELAXATION_TIME : public descriptors::FIELD_BASE<1> { };
struct C_0 : public descriptors::FIELD_BASE<1> { };
struct PSI_BC : public descriptors::FIELD_BASE<1> { };
struct CB_CATION : public descriptors::FIELD_BASE<1> { };
struct CB_ANION : public descriptors::FIELD_BASE<1> { };
struct ERROR_NORMS_PSI : public descriptors::FIELD_BASE<3> { };
struct ERROR_NORMS_CATION : public descriptors::FIELD_BASE<3> { };
struct ERROR_NORMS_ANION : public descriptors::FIELD_BASE<3> { };
struct POISSON_LOOP : public descriptors::TYPED_FIELD_BASE<bool,1> { };
struct POISSON_LOOP_PER_TOTAL_TIME_STEP : public descriptors::TYPED_FIELD_BASE<int,1> { };

}

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const T physDeltaX = extend[1]/parameters.get<parameters::RESOLUTION>();

  extend[2] = 5*physDeltaX;
  Vector<T,3> origin( 0,0,0 );
  IndicatorCuboid3D<T> cuboid( extend, origin );

  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({true,false,true});
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();

  const Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const T physDeltaX = extend[1]/parameters.get<parameters::RESOLUTION>();

  geometry.rename( 0,2 );
  geometry.rename( 2,1,{0,1,0} );
  geometry.clean();

  Vector<T,3> extendI(extend[0], 2*physDeltaX, 5*physDeltaX);
  Vector<T,3> origin;
  IndicatorCuboid3D<T> inlet( extendI, origin );
  geometry.rename( 2,3,1,inlet );

  origin[1] = extend[1]-2*physDeltaX;
  IndicatorCuboid3D<T> outlet( extendI, origin );
  geometry.rename( 2,4,1,outlet );

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

  const int N = parameters.get<parameters::RESOLUTION>();
  const T relaxationTimePoisson = parameters.get<parameters::POISSON_RELAXATION_TIME>();
  Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const T charPhysVelocity = parameters.get<parameters::PHYS_CHAR_VELOCITY>();

  // Set up a unit converter with the characteristic physical units
  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR>>(
    N,
    relaxationTimePoisson,
    extend[1],
    charPhysVelocity,
    1.,
    1.
  );
  lattice.getUnitConverter().print();

  const T omega = lattice.getUnitConverter().getLatticeRelaxationFrequency();

  // Material=1,3,4 -->bulk dynamics
  dynamics::set<SourcedAdvectionDiffusionBGKdynamics>(lattice, geometry.getMaterialIndicator({1, 3, 4}));

  // Material=3,4 -->Dirichlet boundary for electric potential
  boundary::set<boundary::AdvectionDiffusionDirichlet>(lattice, geometry, 3);
  boundary::set<boundary::AdvectionDiffusionDirichlet>(lattice, geometry, 4);

  lattice.setParameter<descriptors::OMEGA>(omega);

  clout << "Prepare Lattice Poisson ... OK" << std::endl;
}

template<int ID>
void prepareLatticeNernstPlanck(MyCase& myCase) {
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice Nernst-Planck ..." << std::endl;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<Concentration<ID>>;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(Concentration<ID>{});

  const int N = parameters.get<parameters::RESOLUTION>();
  const T relaxationTime = parameters.get<parameters::IONS_RELAXATION_TIME>();
  Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const T charPhysVelocity = parameters.get<parameters::PHYS_CHAR_VELOCITY>();
  const T diffusion = parameters.get<parameters::DIFFUSION>();

  // Set up a unit converter with the characteristic physical units
  lattice.template setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR>>(
    N,
    relaxationTime,
    extend[1],
    charPhysVelocity,
    diffusion,
    1.
  );
  lattice.getUnitConverter().print();

  const T omega = lattice.getUnitConverter().getLatticeRelaxationFrequency();

  // Material=1, 4 -->bulk dynamics
  dynamics::set<AdvectionDiffusionBGKdynamics>(lattice, geometry.getMaterialIndicator({1, 4}));

  // Material=3 -->Bounce Back (Zero Gradient)
  dynamics::set<BounceBack>(lattice, geometry, 3);

  // Material=4 -->Dirichlet boundary for concentrations
  boundary::set<boundary::AdvectionDiffusionDirichlet>(lattice, geometry, 4);

  lattice.template setParameter<descriptors::OMEGA>(omega);

  clout << "Prepare Lattice Nernst-Planck ... OK" << std::endl;
}

void prepareLatticeCoupling(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& latticePoisson = myCase.getLattice(Poisson{});
  auto& latticeCation = myCase.getLattice(Concentration<0>{});
  auto& latticeAnion = myCase.getLattice(Concentration<1>{});
  const T diffusion = parameters.get<parameters::DIFFUSION>();
  const T valence = parameters.get<parameters::VALENCE>();
  const T temperature = parameters.get<parameters::TEMPERATURE>();
  const T dielectricC = parameters.get<parameters::DIELECTRIC_CONSTANT>();
  const T npVelCoeff = physConstants::elementaryCharge<T>() * valence * diffusion / physConstants::boltzmannConstant<T>() / temperature / latticeCation.getUnitConverter().getConversionFactorVelocity();
  const T sourceCoeff = 1./dielectricC * physConstants::faradayConstant<T>() * valence * latticePoisson.getUnitConverter().getConversionFactorTime();
  auto& coupling = myCase.setCouplingOperator(
    "PNP",
    PNPCoupling{},
    names::Concentration0{}, latticeCation,
    names::Concentration1{}, latticePoisson,
    names::Concentration2{}, latticeAnion);
  coupling.setParameter<PNPCoupling::DX>(latticePoisson.getUnitConverter().getPhysDeltaX());
  coupling.setParameter<PNPCoupling::NPVELCOEFF>(npVelCoeff);
  coupling.setParameter<PNPCoupling::POISSONCOEFF>(sourceCoeff);
  coupling.setParameter<PNPCoupling::OMEGA>(latticePoisson.getUnitConverter().getLatticeRelaxationFrequency());
  coupling.restrictTo(geometry.getMaterialIndicator({1}));
}

void setInitialValuesPoisson(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(Poisson{});
  auto& parameters = myCase.getParameters();

  const T psi0 = parameters.get<parameters::PSI_BC>();
  momenta::setElectricPotential(lattice, geometry.getMaterialIndicator(3), psi0);
  momenta::setElectricPotential(lattice, geometry.getMaterialIndicator({1, 4}), T(0));

  // Make the lattice ready for simulation
  lattice.initialize();
}

void setInitialValuesCation(MyCase& myCase) {
  OstreamManager clout( std::cout,"setInitialValuesCation" );
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(Concentration<0>{});

  const T C0 = parameters.get<parameters::C_0>();
  momenta::setConcentration(lattice, geometry.getMaterialIndicator(4), C0);
  momenta::setConcentration(lattice, geometry.getMaterialIndicator({1, 3}), T(0));

  // Make the lattice ready for simulation
  lattice.initialize();
}

void setInitialValuesAnion(MyCase& myCase) {
  OstreamManager clout( std::cout,"setInitialValuesAnion" );
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(Concentration<1>{});

  const T C0 = parameters.get<parameters::C_0>();
  momenta::setConcentration(lattice, geometry.getMaterialIndicator(4), C0);
  momenta::setConcentration(lattice, geometry.getMaterialIndicator({1, 3}), T(0));

  // Make the lattice ready for simulation
  lattice.initialize();
}

void evaluateError(MyCase& myCase) {
  OstreamManager clout( std::cout,"error" );
  using T = MyCase::value_t;
  using PDESCRIPTOR = MyCase::descriptor_t_of<Poisson>;
  using CDESCRIPTOR = MyCase::descriptor_t_of<Concentration<0>>;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& latticePoisson = myCase.getLattice(Poisson{});
  auto& latticeCation = myCase.getLattice(Concentration<0>{});
  auto& latticeAnion = myCase.getLattice(Concentration<1>{});
  latticePoisson.setProcessingContext(ProcessingContext::Evaluation);
  latticeCation.setProcessingContext(ProcessingContext::Evaluation);
  latticeAnion.setProcessingContext(ProcessingContext::Evaluation);
  int tmp[]= { };
  T result[2] = { };

  const Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const T valence = parameters.get<parameters::VALENCE>();
  const T psi0 = parameters.get<parameters::PSI_BC>();
  const T C0 = parameters.get<parameters::C_0>();
  const T Debye = parameters.get<parameters::DEBYE>();
  const T temperature = parameters.get<parameters::TEMPERATURE>();
  PotentialProfile1D<T,T,PDESCRIPTOR> psiSol(0., extend[1], psi0, Debye, latticePoisson.getUnitConverter());
  ConcentrationProfile1D<T,T,CDESCRIPTOR> concSolCation(C0, valence, temperature, psiSol);
  ConcentrationProfile1D<T,T,CDESCRIPTOR> concSolAnion(C0, -valence, temperature, psiSol);
  SuperLatticeDensity3D<T, PDESCRIPTOR> psi( latticePoisson );
  SuperLatticeDensity3D<T, CDESCRIPTOR> cation( latticeCation );
  SuperLatticeDensity3D<T, CDESCRIPTOR> anion( latticeAnion );

  auto material = geometry.getMaterialIndicator(1);
  SuperRelativeErrorL1Norm3D<T>   errorPsiL1Norm(psi, psiSol, *material);
  SuperRelativeErrorL2Norm3D<T>   errorPsiL2Norm(psi, psiSol, *material);
  SuperRelativeErrorLinfNorm3D<T> errorPsiLinfNorm(psi, psiSol, *material);

  SuperRelativeErrorL1Norm3D<T>   errorConcL1Norm(cation, concSolCation, *material);
  SuperRelativeErrorL2Norm3D<T>   errorConcL2Norm(cation, concSolCation, *material);
  SuperRelativeErrorLinfNorm3D<T> errorConcLinfNorm(cation, concSolCation, *material);

  SuperRelativeErrorL1Norm3D<T>   errorConc2L1Norm(anion, concSolAnion, *material);
  SuperRelativeErrorL2Norm3D<T>   errorConc2L2Norm(anion, concSolAnion, *material);
  SuperRelativeErrorLinfNorm3D<T> errorConc2LinfNorm(anion, concSolAnion, *material);

  Vector<T,3> errorsPsi;
  Vector<T,3> errorsCation;
  Vector<T,3> errorsAnion;
  errorPsiL1Norm(result,tmp);
  clout << "Relative Potential-L1-error: " << result[0] << std::endl;
  errorsPsi[0] = result[0];

  errorPsiL2Norm(result,tmp);
  clout << "Relative Potential-L2-error: " << result[0] << std::endl;
  errorsPsi[1] = result[0];

  errorPsiLinfNorm(result,tmp);
  clout << "Relative Potential-Linf-error: " << result[0] << std::endl;
  errorsPsi[2] = result[0];

  errorConcL1Norm(result,tmp);
  clout << "Relative Cation Conc-L1-error: " << result[0] << std::endl;
  errorsCation[0] = result[0];

  errorConcL2Norm(result,tmp);
  clout << "Relative Cation Conc-L2-error: " << result[0] << std::endl;
  errorsCation[1] = result[0];

  errorConcLinfNorm(result,tmp);
  clout << "Relative Cation Conc-Linf-error: " << result[0] << std::endl;
  errorsCation[2] = result[0];

  errorConc2L1Norm(result,tmp);
  clout << "Relative Anion Conc-L1-error: " << result[0] << std::endl;
  errorsAnion[0] = result[0];

  errorConc2L2Norm(result,tmp);
  clout << "Relative Anion Conc-L2-error: " << result[0] << std::endl;
  errorsAnion[1] = result[0];

  errorConc2LinfNorm(result,tmp);
  clout << "Relative Anion Conc-Linf-error: " << result[0] << std::endl;
  errorsAnion[2] = result[0];

  parameters.set<parameters::ERROR_NORMS_PSI>(errorsPsi);
  parameters.set<parameters::ERROR_NORMS_CATION>(errorsCation);
  parameters.set<parameters::ERROR_NORMS_ANION>(errorsAnion);
}

void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  OstreamManager clout( std::cout,"getResults" );
  using T = MyCase::value_t;
  using PDESCRIPTOR = MyCase::descriptor_t_of<Poisson>;
  using CDESCRIPTOR = MyCase::descriptor_t_of<Concentration<0>>;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& latticePoisson = myCase.getLattice(Poisson{});
  auto& latticeCation = myCase.getLattice(Concentration<0>{});
  auto& latticeAnion = myCase.getLattice(Concentration<1>{});
  const T maxPhysT = parameters.get<parameters::PHYS_MAX_T>();
  bool hasConverged = parameters.get<parameters::CONVERGED>();
  const bool lastTimeStep = ( hasConverged || (iT + 1 == latticeCation.getUnitConverter().getLatticeTime( maxPhysT )) );
  const Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const T valence = parameters.get<parameters::VALENCE>();
  const T psi0 = parameters.get<parameters::PSI_BC>();
  const T C0 = parameters.get<parameters::C_0>();
  const T Debye = parameters.get<parameters::DEBYE>();
  const T temperature = parameters.get<parameters::TEMPERATURE>();
  const T outputInterval = parameters.get<parameters::PHYS_STAT_ITER_T>();

  const int vtkIter  = latticeCation.getUnitConverter().getLatticeTime( outputInterval );
  const int statIter = latticeCation.getUnitConverter().getLatticeTime( outputInterval );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperVTMwriter3D<T> vtmWriterCation( "cation" );
    SuperLatticeCuboid3D<T, CDESCRIPTOR> cuboid( latticeCation );
    SuperLatticeRank3D<T, CDESCRIPTOR> rank( latticeCation );

    vtmWriterCation.write( cuboid );
    vtmWriterCation.write( rank );
    vtmWriterCation.createMasterFile();

    SuperVTMwriter3D<T> vtmWriterPoisson( "poisson" );
    vtmWriterPoisson.createMasterFile();

    SuperVTMwriter3D<T> vtmWriterAnion( "anion" );
    vtmWriterAnion.createMasterFile();
  }

  // Get statistics
  if ( iT%statIter == 0 || lastTimeStep ) {

    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    latticeCation.getStatistics().print( iT,latticeCation.getUnitConverter().getPhysTime( iT ) );
  }

  // Writes the VTK files
  if ( iT%vtkIter == 0 || lastTimeStep) {
    latticeCation.setProcessingContext(ProcessingContext::Evaluation);
    latticeAnion.setProcessingContext(ProcessingContext::Evaluation);
    latticePoisson.setProcessingContext(ProcessingContext::Evaluation);

    latticePoisson.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriterPoisson( "poisson" );
      SuperLatticeDensity3D psi( latticePoisson);
      psi.getName() = "psi";
      PotentialProfile1D<T,T,PDESCRIPTOR> psiSol(0., extend[1], psi0, Debye, latticePoisson.getUnitConverter());
      SuperLatticeFfromAnalyticalF3D<T,PDESCRIPTOR> analyticalPsi(psiSol, latticePoisson);
      analyticalPsi.getName() = "analytical potential";
      SuperGeometryF3D<T> geom(geometry);
      vtmWriterPoisson.addFunctor(geom);
      vtmWriterPoisson.addFunctor(psi);
      vtmWriterPoisson.addFunctor(analyticalPsi);
      task(vtmWriterPoisson, iT);
    });

    latticeCation.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriterCation( "cation" );
      SuperLatticeDensity3D<T,CDESCRIPTOR> cation( latticeCation );
      cation.getName() = "cation";
      SuperLatticePhysField3D<T,CDESCRIPTOR,VELOCITY> velCation( latticeCation, latticeCation.getUnitConverter().getConversionFactorVelocity() );
      velCation.getName() = "velCation";
      PotentialProfile1D<T,T,PDESCRIPTOR> psiSol(0., extend[1], psi0, Debye, latticePoisson.getUnitConverter());
      ConcentrationProfile1D<T,T,CDESCRIPTOR> concSolCation(C0, valence, temperature, psiSol);
      SuperLatticeFfromAnalyticalF3D<T,CDESCRIPTOR> analyticalConcCation(concSolCation, latticeCation);
      analyticalConcCation.getName() = "analytical concentration cation";
      SuperGeometryF3D<T> geom(geometry);
      vtmWriterCation.addFunctor(geom);
      vtmWriterCation.addFunctor(cation);
      vtmWriterCation.addFunctor(velCation);
      vtmWriterCation.addFunctor(analyticalConcCation);
      task(vtmWriterCation, iT);
    });

    latticeAnion.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriterAnion( "anion" );
      SuperLatticeDensity3D<T,CDESCRIPTOR> anion( latticeAnion );
      anion.getName() = "anion";
      SuperLatticePhysField3D<T,CDESCRIPTOR,VELOCITY> velAnion( latticeAnion, latticeCation.getUnitConverter().getConversionFactorVelocity() );
      velAnion.getName() = "velAnion";
      PotentialProfile1D<T,T,PDESCRIPTOR> psiSol(0., extend[1], psi0, Debye, latticePoisson.getUnitConverter());
      ConcentrationProfile1D<T,T,CDESCRIPTOR> concSolAnion(C0, -valence, temperature, psiSol);
      SuperLatticeFfromAnalyticalF3D<T,CDESCRIPTOR> analyticalConcAnion(concSolAnion, latticeAnion);
      analyticalConcAnion.getName() = "analytical concentration anion";
      SuperGeometryF3D<T> geom(geometry);
      vtmWriterAnion.addFunctor(geom);
      vtmWriterAnion.addFunctor(anion);
      vtmWriterAnion.addFunctor(velAnion);
      vtmWriterAnion.addFunctor(analyticalConcAnion);
      task(vtmWriterAnion, iT);
    });
   }

   if ( lastTimeStep ) {
     evaluateError(myCase);
   }
}

void simulatePoisson(MyCase& myCase) {
  OstreamManager clout( std::cout,"simulatePoisson" );
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(Poisson{});
  const T maxPhysT = parameters.get<parameters::PHYS_MAX_T>();
  const T residuum = parameters.get<parameters::RESIDUUM>();
  const bool poissonLoop = parameters.get<parameters::POISSON_LOOP>();

  util::Timer<T> timer( lattice.getUnitConverter().getLatticeTime( maxPhysT ), geometry.getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( 10, residuum );
  timer.start();

  for ( std::size_t iT = 0; iT < (!poissonLoop)*1 + (poissonLoop)*lattice.getUnitConverter().getLatticeTime(maxPhysT ); ++iT ) {

    lattice.collideAndStream();

    if(converge.hasConverged()) { break; }

    converge.takeValue( lattice.getStatistics().getAverageRho(), false );
  }
  timer.stop();
}

void simulate(MyCase& myCase) {
  OstreamManager clout( std::cout,"simulate" );
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  const bool poissonLoop = parameters.get<parameters::POISSON_LOOP>();
  const int poissonLoopPerTotalSteps = parameters.get<parameters::POISSON_LOOP_PER_TOTAL_TIME_STEP>();

  const std::size_t iTmax = myCase.getLattice(Concentration<0>{}).getUnitConverter().getLatticeTime(
    parameters.get<parameters::PHYS_MAX_T>());

  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  const T residuum = parameters.get<parameters::RESIDUUM>();
  util::ValueTracer<T> converge( 500, residuum );
  timer.start();

  parameters.set<parameters::CONVERGED>(false);
  for (std::size_t iT=0; iT < iTmax; ++iT) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << std::endl;
      parameters.set<parameters::CONVERGED>(true);
      getResults(myCase, timer, iT);

      break;
    }
    // internal Poisson loop. In terms of efficiency choose an appropriate frequency of Poisson equation updates
    if( iT == 0 || iT%((!poissonLoop)*1 + (poissonLoop)*poissonLoopPerTotalSteps) == 0){
      simulatePoisson(myCase);
    }

    myCase.getOperator("PNP").apply();

    myCase.getLattice(Concentration<0>{}).collideAndStream();
    myCase.getLattice(Concentration<1>{}).collideAndStream();

    getResults(myCase, timer, iT);
    converge.takeValue( myCase.getLattice(Poisson{}).getStatistics().getAverageRho(), false );
  }

  timer.stop();
  timer.printSummary();
  singleton::pool().wait();
}
