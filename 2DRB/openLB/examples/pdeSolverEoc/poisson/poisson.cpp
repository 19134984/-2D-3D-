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
  * In current example, Poisson equation with one-dimensional analytical solution is discretized with LBM.
  * This equation is often used for electric potential distribution. The analytical solution for this case is taken from
  * Z. Chai, B. Shi, A novel lattice boltzmann model for the poisson equation, Applied Mathematical Modelling 32 (10) (2008) 2050–2058
  * doi:https://doi.org/10.1016/j.apm.2007.06.033.
  */

#include <olb.h>

using namespace olb;
using namespace olb::names;
using namespace olb::descriptors;

// Anayltical solution for potential profile
template <typename T, typename S, typename DESCRIPTOR>
class PotentialProfile1D : public AnalyticalF3D<T, S>
{
private:
  T x0, x1, k;
  UnitConverter<T, DESCRIPTOR> const &converter;

public:
  PotentialProfile1D(T x0_, T x1_, T k_, UnitConverter<T, DESCRIPTOR> const &converter_) : AnalyticalF3D<T, S>(1),
                                                                                                                  x0(x0_), x1(x1_), k(k_), converter(converter_)
  {
      this->getName() = "PotentialProfile1D";
  };

  bool operator()(T output[1], const S x[3])
  {
      T distX = (x[0] - x0) / (x1 - x0);
      T L = converter.getPhysDeltaX();
      output[0] = (util::exp(k) - 1.)/(util::exp(k) - util::exp(-k))*util::exp(-k*distX) + (1. - util::exp(-k))/(util::exp(k) - util::exp(-k))*util::exp(k*distX);
      if(x[0] <= L/2. || x[0] >= x1 - L/2.) output[0] = 1.;

      return true;
  };
};

// PostProcessor for updating of source term in each time cell in each cell
struct ComputeSourceTerm {
  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  struct C_T : public descriptors::FIELD_BASE<1> { };
  struct CONSTK : public descriptors::FIELD_BASE<1> { };

  using parameters = meta::list<C_T,CONSTK>;

  int getPriority() const {
    return 0;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform {
    using V = typename CELL::value_t;
    V psi = cell.computeRho();
    V k = parameters.template get<CONSTK>();
    V source = -k*k*psi*parameters.template get<C_T>();
    cell.template setField<descriptors::SOURCE>(source);
  }
};

namespace olb::parameters {

struct CONSTK : public descriptors::FIELD_BASE<1> { };
struct ERROR_NORMS : public descriptors::FIELD_BASE<3> { };

}

using MyCase = Case<
  Poisson, Lattice<double, descriptors::D3Q7<VELOCITY,SOURCE>>
>;

Mesh<MyCase::value_t,MyCase::d> createMesh(MyCase::ParametersD& parameters) {
  using T = MyCase::value_t;
  Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const T physDeltaX = extend[1]/parameters.get<parameters::RESOLUTION>();

  extend[2] = 3*physDeltaX;
  Vector<T,3> origin( 0,0,0 );
  IndicatorCuboid3D<T> cuboid( extend, origin );

  Mesh<T,MyCase::d> mesh(cuboid, physDeltaX, singleton::mpi().getSize());
  mesh.setOverlap(parameters.get<parameters::OVERLAP>());
  mesh.getCuboidDecomposition().setPeriodicity({false,true,true});
  return mesh;
}

void prepareGeometry(MyCase& myCase) {
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  auto& geometry = myCase.getGeometry();

  geometry.rename( 0,2 );
  geometry.rename( 2,1,{1,0,0} );
  geometry.clean();

  // Removes all not needed boundary voxels outside the surface
  geometry.clean();
  // Removes all not needed boundary voxels inside the surface
  geometry.innerClean();
  geometry.checkForErrors();
  geometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(MyCase& myCase) {
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<Poisson>;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(Poisson{});

  const int N = parameters.get<parameters::RESOLUTION>();
  const T relaxationTime = parameters.get<parameters::LATTICE_RELAXATION_TIME>();
  Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();

  // Set up a unit converter with the characteristic physical units
  lattice.setUnitConverter<UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR>>(
    N,
    relaxationTime,
    extend[1],
    1.,
    1.,
    1.
  );
  lattice.getUnitConverter().print();

  const T omega = lattice.getUnitConverter().getLatticeRelaxationFrequency();

  // Material=1,2 -->bulk dynamics
  dynamics::set<SourcedAdvectionDiffusionBGKdynamics>(lattice, geometry.getMaterialIndicator({1, 2}));

  // Material=2 -->Dirichlet boundary
  boundary::set<boundary::AdvectionDiffusionDirichlet>(lattice, geometry, 2);

  lattice.setParameter<descriptors::OMEGA>(omega);

  const T constK = parameters.get<parameters::CONSTK>();
  lattice.addPostProcessor<stage::PreCollide>(meta::id<ComputeSourceTerm>{});
  lattice.setParameter<ComputeSourceTerm::C_T>(lattice.getUnitConverter().getConversionFactorTime());
  lattice.setParameter<ComputeSourceTerm::CONSTK>(constK);

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setInitialValues(MyCase& myCase) {
  using T = MyCase::value_t;
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(Poisson{});
  momenta::setElectricPotential(lattice, geometry.getMaterialIndicator(1), T(0));
  momenta::setElectricPotential(lattice, geometry.getMaterialIndicator(2), T(1));
  lattice.initialize();
}

void evaluateError(MyCase& myCase) {
  OstreamManager clout( std::cout,"error" );
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<Poisson>;
  auto& parameters = myCase.getParameters();
  auto& geometry = myCase.getGeometry();
  auto& lattice = myCase.getLattice(Poisson{});
  lattice.setProcessingContext(ProcessingContext::Evaluation);

  int tmp[]= { };
  T result[2] = { };

  const Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const T constK = parameters.get<parameters::CONSTK>();
  PotentialProfile1D<T,T,DESCRIPTOR> psiSol(0, extend[0], constK, lattice.getUnitConverter());
  SuperLatticeDensity3D<T, DESCRIPTOR> density( lattice );

  auto material = geometry.getMaterialIndicator(1);
  SuperRelativeErrorL1Norm3D<T>   errorPsiL1Norm(density, psiSol, *material);
  SuperRelativeErrorL2Norm3D<T>   errorPsiL2Norm(density, psiSol, *material);
  SuperRelativeErrorLinfNorm3D<T> errorPsiLinfNorm(density, psiSol, *material);

  Vector<T,3> errors;
  errorPsiL1Norm(result,tmp);
  clout << "Psi-L1-error: " << result[0] << std::endl;
  errors[0] = result[0];

  errorPsiL2Norm(result,tmp);
  clout << "Psi-L2-error: " << result[0] << std::endl;
  errors[1] = result[0];

  errorPsiLinfNorm(result,tmp);
  clout << "Psi-Linf-error: " << result[0] << std::endl;
  errors[2] = result[0];

  parameters.set<parameters::ERROR_NORMS>(errors);
}

void getResults(MyCase& myCase,
                util::Timer<MyCase::value_t>& timer,
                std::size_t iT)
{
  OstreamManager clout( std::cout,"getResults" );
  using T = MyCase::value_t;
  using DESCRIPTOR = MyCase::descriptor_t_of<Poisson>;
  auto& parameters = myCase.getParameters();
  auto& lattice = myCase.getLattice(Poisson{});
  const T maxPhysT = parameters.get<parameters::PHYS_MAX_T>();
  bool hasConverged = parameters.get<parameters::CONVERGED>();
  const bool lastTimeStep = ( hasConverged || (iT + 1 == lattice.getUnitConverter().getLatticeTime( maxPhysT )) );

  const int vtkIter  = lattice.getUnitConverter().getLatticeTime( .1 );
  const int statIter = lattice.getUnitConverter().getLatticeTime( .1 );

  const Vector extend = parameters.get<parameters::DOMAIN_EXTENT>();
  const T constK = parameters.get<parameters::CONSTK>();
  PotentialProfile1D<T,T,DESCRIPTOR> psiSol(0, extend[0], constK, lattice.getUnitConverter());
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> analyticalSol(psiSol, lattice);
  analyticalSol.getName() = "analytical solution";

  SuperVTMwriter3D<T> vtmWriter( "poisson1d" );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( lattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( lattice );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Get statistics
  if ( iT%statIter == 0 ) {

    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    lattice.getStatistics().print( iT,lattice.getUnitConverter().getPhysTime( iT ) );
  }

  // Writes the VTK files
  if ( iT%vtkIter == 0 || lastTimeStep) {
    lattice.setProcessingContext(ProcessingContext::Evaluation);

    SuperLatticeDensity3D<T,DESCRIPTOR> potential( lattice );
    potential.getName() = "potential";

    vtmWriter.addFunctor( potential );
    vtmWriter.addFunctor( analyticalSol );
    vtmWriter.write( iT );
  }

  if ( iT%statIter==0 || lastTimeStep ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    lattice.getStatistics().print( iT,lattice.getUnitConverter().getPhysTime( iT ) );

    // Error norms
    if ( lastTimeStep ) {
      evaluateError(myCase);
    }
  }
}

void simulate(MyCase& myCase) {
  OstreamManager clout( std::cout,"simulate" );
  using T = MyCase::value_t;
  auto& parameters = myCase.getParameters();
  const T maxPhysT = parameters.get<parameters::PHYS_MAX_T>();
  const T residuum = parameters.get<parameters::RESIDUUM>();

  const std::size_t iTmax = myCase.getLattice(Poisson{}).getUnitConverter().getLatticeTime(maxPhysT);
  util::Timer<T> timer(iTmax, myCase.getGeometry().getStatistics().getNvoxel());
  timer.start();
  util::ValueTracer<T> converge( 100, residuum );
  timer.start();

  parameters.set<parameters::CONVERGED>(false);
  for (std::size_t iT=0; iT < iTmax; ++iT) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << std::endl;
      parameters.set<parameters::CONVERGED>(true);
      getResults(myCase, timer, iT);

      break;
    }

    myCase.getLattice(Poisson{}).collideAndStream();

    getResults(myCase, timer, iT);
    converge.takeValue( myCase.getLattice(Poisson{}).getStatistics().getAverageRho(), false );
  }

  timer.stop();
  timer.printSummary();
  singleton::pool().wait();
}

static Gnuplot<MyCase::value_t> gplot(
  "Psi_eoc",
  false,
  "set terminal png size 720, 720 font 'Arial,10'",
  Gnuplot<MyCase::value_t>::LOGLOG,
  Gnuplot<MyCase::value_t>::LINREG);

int main(int argc, char* argv[]) {
  initialize(&argc, &argv);

  /// === Step 2: Set Parameters ===
  MyCase::ParametersD myCaseParameters;
  {
    using namespace olb::parameters;
    myCaseParameters.set<RESOLUTION_START>(20);
    myCaseParameters.set<RESOLUTION_END>(80);
    myCaseParameters.set<LATTICE_RELAXATION_TIME>(1.);
    myCaseParameters.set<DOMAIN_EXTENT>({1.,0.25,0.});
    myCaseParameters.set<CONSTK>(27.79);
    myCaseParameters.set<PHYS_MAX_T>(10.);
    myCaseParameters.set<RESIDUUM>(1.e-9);
    myCaseParameters.set<CONVERGED>(false);
    myCaseParameters.set<ERROR_NORMS>({0.,0.,0.});
  }

  gplot.setLabel("Resolution test", "average Error");

  for(int simuN =  myCaseParameters.get<parameters::RESOLUTION_START>();
      simuN <=  myCaseParameters.get<parameters::RESOLUTION_END>();
      simuN *=2) {
    myCaseParameters.set<parameters::RESOLUTION>(simuN);
    /// === Step 3: Create Mesh ===
    Mesh mesh = createMesh(myCaseParameters);

    /// === Step 4: Create Case ===
    MyCase myCase(myCaseParameters, mesh);

    /// === Step 5: Prepare Geometry ===
    prepareGeometry(myCase);

    /// === Step 6: Prepare Lattice ===
    prepareLattice(myCase);

    /// === Step 7: Definition of Initial, Boundary Values, and Fields ===
    setInitialValues(myCase);

    /// === Step 8: Simulate ===
    simulate(myCase);

    auto errors = myCaseParameters.get<parameters::ERROR_NORMS>();
    gplot.setData (
      MyCase::value_t(simuN),
      { errors[0], errors[1], errors[2]},
      { "psi L1 rel Error","psi L2 rel Error",
        "psi Linf rel error"},
      "top right",
      { 'p','p','p' }
    );
  }

  gplot.writePNG();
}
