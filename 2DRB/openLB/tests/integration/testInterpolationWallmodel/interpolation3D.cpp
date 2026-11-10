/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006, 2007, 2012 Jonas Latt, Mathias J. Krause,
 *  Louis Kronberg, Christian Vorwerk, Bastian Schäffauer
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

/* interpolation3D.cpp:
 * Test environment for the interpolation operator in turbulentWallModelPostProcessor3D.
 * 3D scalar field is defined on a cuboid, where a ball is placed inside as obstacle.
 * For every boundary cell on the surface of the ball, the scalar value is interpolated
 * at the chosen wall normal distance, called the exchange location, and compared to the analytical value.
 * The difference between interpolated value and analytical value - adjusted by the wall normal distance -
 * is set to an ERROR scalar field and its L1,L2 and Linf norms are plotted for increasing resolutions of the mesh.
 */

#include <olb.h>
#include <cmath>
#include <iostream>

using namespace std;

using namespace olb;
using namespace olb::descriptors;
using namespace olb::util;
using namespace olb::graphics;

struct COORDINATES : public FIELD_BASE<0,1,0> {}; // Field for Coordinates to determine analytical values of scalar fields.
struct ERROR : public FIELD_BASE<1> {};           // Field for error.
struct FUNC : public FIELD_BASE<1> {};            // Field for function values.

using T = FLOATING_POINT_TYPE;;
using DESCRIPTOR = D3Q19<SCALAR, FUNC, COORDINATES, ERROR>;

// Mathmatical constants
const T pi = util::acos(-1);

const T lx  = 4.;       // x-side of cuboid
const T ly  = 4.;       // y-side of cuboid
const T lz  = 4.;       // z-side of cuboid
const T r = 1.5;        // Radius of ball
int N = 10;             // resolution of the model
T dx = 0.;

// variables for eoc analysis
T valueAbsL1Error = 0;
T valueAbsL2Error = 0;
T valueAbsLinfError = 0;

// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR>& converter,
                      SuperGeometry<T,3>& superGeometry,
                      IndicatorSphere3D<T> ball )
{

  OstreamManager clout(std::cout, "prepareGeometry");

  clout << "Prepare Geometry ..." << std::endl;

// Rectangle from vector origin to vector extend
  superGeometry.rename(0, 1);

  superGeometry.rename( 1, 2, ball );

  superGeometry.checkForErrors();

  superGeometry.print();

  dx = converter.getPhysDeltaX();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Define field COORDINATES to adjust location of analytical values by wall normal distance.

// Functor -------------------------------------------------------
template <typename T, typename S, typename DESCRIPTOR>
class Coordinates : public AnalyticalF3D<T,S> {

protected:
  T s;
public:
  Coordinates(UnitConverter<S,DESCRIPTOR> const& converter) : AnalyticalF3D<T,S>(1)
  {
    this->getName() = "Coordinates";
  };

  bool operator()(T output[3], const S input[3])
  {

    T x = input[0];
    T y = input[1];
    T z = input[2];

    output[0] = x;
    output[1] = y;
    output[2] = z;

    return true;
  };
};

// Define scalar field.

// Functor -------------------------------------------------------
template <typename T, typename S, typename DESCRIPTOR>
class SinusField : public AnalyticalF3D<T,S> {

protected:
  T s;
public:
  SinusField(UnitConverter<S,DESCRIPTOR> const& converter) : AnalyticalF3D<T,S>(1)
  {
    this->getName() = "SinusField";
  };

  bool operator()(T output[3], const S input[3])
  {

    T x = input[0];
    T y = input[1];
    T z = input[2];

    //Nonlinear scalar field g1(x,y,z) = 2sin(2*pi*x) + 2sin(2*pi*y) + 2sin(2*pi*z).
    output[0] = util::sin((T) (2*pi*x)) *2 + util::sin((T) (2*pi*y)) *2 + util::sin((T) (2*pi*z)) *2;

    // Linear scalar field g2(x,y,z) = x + y + z.
    //output[0] = x + y + z; // For Interpolation Value

    return true;
  };
};

// Define interpolation operator as in turbulentWallModelPostProcessor3D.

/// @brief trilinear interpolation of function f
/// @tparam DataDim Dimension of function return value
/// @param distance Distance from cell
/// @param result Result is written here
/// @param f Function which is to interpolate
// Caution: Developer has to guarantee that the data are accessible.

template <unsigned DataDim, typename CELL, typename F, typename T = typename CELL::value_t>
void interpolate3d1(CELL cell, const Vector<T,3> distance, T result[DataDim], F f) any_platform{
  OstreamManager clout( std::cout,"interpolate3d1" );
  const Vector<int,3> floorV (util::floor(distance[0]), util::floor(distance[1]), util::floor(distance[2]));
  Vector<Vector<int,3>,8> surroundingPoints (floorV);
  surroundingPoints[1][0] += 1;
  surroundingPoints[2][1] += 1;
  surroundingPoints[3][2] += 1;
  surroundingPoints[4][0] += 1; surroundingPoints[4][1] += 1;
  surroundingPoints[5][0] += 1; surroundingPoints[5][2] += 1;
  surroundingPoints[6][1] += 1; surroundingPoints[6][2] += 1;
  surroundingPoints[7][0] += 1; surroundingPoints[7][1] += 1; surroundingPoints[7][2] += 1;

  T resNeighbor[DataDim];
  std::fill_n(result, DataDim, T(0));
  for (auto point : surroundingPoints) {
    const Vector<T,3> dist = distance - point;
    const T weight = (T(1) - util::abs(dist[0])) * (T(1) - util::abs(dist[1])) * (T(1) - util::abs(dist[2]));
    f(cell.neighbor(point), resNeighbor);
    for (unsigned i = 0; i < DataDim; ++i) {
      //clout << "i =" << i << std::endl;
      result[i] += weight * resNeighbor[i];
      //clout << "result[i] =" << result[i] << std::endl;
    }
  //clout << "result =" << result[0] << std::endl;
  }
}

// as above, but with Vector typed result
template <unsigned DataDim, typename CELL, typename F, typename T = typename CELL::value_t>
void interpolate3d1(CELL cell, const Vector<T,3> distance, Vector<T,DataDim>& result, F f) any_platform{
  interpolate3d1<DataDim,CELL,F,T>(cell, distance, result.data(), f);
}

/// Trilinear interpolation of FIELD
template <typename FIELD, typename CELL, typename T = typename CELL::value_t>
void interpolate3d1(CELL cell, const Vector<T,3> distance, T result[]) any_platform{
  OstreamManager clout( std::cout,"interpolate3d1" );
  FieldD<T, typename CELL::descriptor_t, FIELD> v;
  constexpr unsigned DataDim = CELL::descriptor_t::template size<FIELD>();
  auto f = [DataDim](CELL c, T res[DataDim]) -> void {
    OstreamManager clout( std::cout,"f" );
    for (unsigned i = 0; i < DataDim; ++i) {
      res[i] = c.template getFieldComponent<FIELD>(i);
      //clout << "res =" << res[i] << std::endl;
    }
  };
  interpolate3d1<DataDim,CELL,decltype(f)>(cell,distance,result,f);
}

/* Modified method "apply" from turbulentWallModelPostProcessor3D such that the difference of interpolated value
 * and analytical value at exchange location are set to ERROR field for every boundary cell.
 */

class TurbulentWallModelPostProcessor1 {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return -1;
  }

template <typename CELL, typename T = typename CELL::value_t>
void apply(CELL& cell) any_platform{
  OstreamManager clout( std::cout,"apply" );

  // Calculate distance between boundary cell and exchange location.
  auto y2 = cell.template getField<descriptors::Y2>();
  auto y1 = cell.template getField<descriptors::Y1>();
  auto y1y2 = y2 - y1;

  // Get coordinates of cell.
  Vector<T,3> Point = cell.template getField<COORDINATES>();

  // interpolate values at exchange location.
  T result[1];
  interpolate3d1<SCALAR>(cell, y1y2, result);

  // Define error as difference between interpolated value and analytical value adjusted by wall normal distance.
  Point += y1y2*dx;

  // Analytical value of scalar field at exchange locaiton.
  T valueSinus = util::sin((T) (2*pi*Point[0])) *2 + util::sin((T) (2*pi*Point[1])) *2 + util::sin((T) (2*pi*Point[2])) *2; // Nonlinear field.
  //T valuer = Point[0] + Point[1] + Point[2];  // Linear field.

  // Difference of interpolated value and analytical value.
  T sinerror = util::abs(result[0] - valueSinus); // Nonlinear field
  //T rerror = util::abs(result[0] - valuer);     // Linear field

  // Set error value to scalar field ERROR.
  cell.template setField<ERROR>(sinerror);  // Nonlinear field
  //cell.template setField<ERROR>(rerror);  // Linear field

}
};

// Set up the geometry of the simulation

void prepareLattice( SuperLattice<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T, DESCRIPTOR> converter,
                     SuperGeometry<T,3>& superGeometry)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  AnalyticalConst3D<T, T> y(0.,0.,0.);
  sLattice.defineField<descriptors::Y1>(superGeometry.getMaterialIndicator({0,1,2}), y);
  sLattice.defineField<descriptors::Y2>(superGeometry.getMaterialIndicator({0,1,2}), y);

  // Define scalar field:

  SinusField<T,T,DESCRIPTOR> sSol(converter);

  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> sSolLattice(sSol, sLattice);

  sLattice.defineField<SCALAR>( superGeometry.getMaterialIndicator({1, 2}), sSolLattice);
  sLattice.defineField<FUNC>( superGeometry.getMaterialIndicator({1, 2}), sSolLattice);

  // Coordinates:
  Coordinates<T,T,DESCRIPTOR> coord(converter);
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> coordLattice(coord, sLattice);
  sLattice.defineField<COORDINATES>( superGeometry.getMaterialIndicator({1, 2}), coordLattice);

  // Initialise modified wall model.
  Vector<T, 3> originB(0., 0., 0.);
  IndicatorSphere3D<T> ball(originB, r);

  setBouzidiBoundary<T, DESCRIPTOR, TurbulentWallModelPostProcessor1>(sLattice, superGeometry, 2, ball);
  setTurbulentWallModel<T, DESCRIPTOR>(sLattice, superGeometry, 2, ball, 1.2);

  // Error:
  AnalyticalConst3D<T,T> zero(0.);
  sLattice.defineField<ERROR>( superGeometry.getMaterialIndicator({1, 2}), zero);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Calculation of error norms.

void computeErrors( SuperLattice<T,DESCRIPTOR>& sLattice,
                    UnitConverter<T,DESCRIPTOR>& converter,
                    SuperGeometry<T,3>& superGeometry)
{

    OstreamManager clout(std::cout,"error");

    int tmp[1] = {0}; T result[1]; //T result1[1];

    SuperLatticeField3D<T, DESCRIPTOR, ERROR> s( sLattice ); // Calculated Value

    // L1 norm:
    SuperL1Norm3D<T> pL1Norm(s,superGeometry,1); // L1-Norm
    pL1Norm(result,tmp); //pSolL1Norm(result1,tmp);
    valueAbsL1Error = result[0];
    //clout << "valueAbsL1Error = " << valueAbsL1Error << endl;

    // L2 norm:
    SuperL2Norm3D<T> pL2Norm(s,superGeometry,1); // L2-Norm
    pL2Norm(result,tmp); //pSolL2Norm(result1,tmp);
    valueAbsL2Error = result[0];
    //clout << "valueAbsL2Error = " << valueAbsL2Error << endl;

    // Linf norm:
    SuperLinfNorm3D<T> pLinfNorm(s,superGeometry,1); // Linf-Norm
    pLinfNorm(result,tmp); //pSolLinfNorm(result1,tmp);
    valueAbsLinfError = result[0];
    //clout << "valueAbsLinfError = " << valueAbsLinfError << endl;

}

// Output to console and files
void getResults( int N,
                 SuperLattice<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR>& converter,
                 SuperGeometry<T,3>& superGeometry,
                 Gnuplot<T>& gplot)
{

  OstreamManager clout( std::cout,"getResults" );

  // VTK and image output only if no EOC analysis
  SuperVTMwriter3D<T> vtmWriter( "interpolation3D" );
  SuperLatticeField3D<T, DESCRIPTOR, ERROR> err( sLattice );
  err.getName() = "err";
  SuperLatticeField3D<T, DESCRIPTOR, SCALAR> value( sLattice );
  value.getName() = "value";
  SuperLatticeField3D<T, DESCRIPTOR, FUNC> function( sLattice );
  function.getName() = "function";

  vtmWriter.addFunctor( err );
  vtmWriter.addFunctor( value );
  vtmWriter.addFunctor( function );
  SuperLatticeField3D<T, DESCRIPTOR, Y1> y1( sLattice );
  y1.getName() = "y1";
  SuperLatticeField3D<T, DESCRIPTOR, Y2> y2( sLattice);
  y2.getName() = "y2";
  vtmWriter.addFunctor( y1 );
  vtmWriter.addFunctor( y2 );

  //vtmWriter.write( err );
  //vtmWriter.write( value );
  //vtmWriter.write( function );

  vtmWriter.createMasterFile();

  vtmWriter.write( err );
  vtmWriter.write( value );
  vtmWriter.write( function );
  vtmWriter.write( y1 );
  vtmWriter.write( y2 );

  // Error Calculation :
  computeErrors( sLattice, converter, superGeometry); // Compute norms of error
  //clout << "valueAbsLinfError = " << valueAbsLinfError << endl;

  // Gnuplot output
  T reso = converter.getResolution();
  gplot.setData( reso, { valueAbsL1Error, valueAbsL2Error, valueAbsLinfError }, { "L1 abs Error","L2 abs Error", "Linf abs error"}, "bottom right", { 'l','l','l' } );
  //gplot.setData( reso, { valueAbsL1Error }, { "value L1 abs Error"}, "bottom right", { 'l' } );

  clout << "Resolution = " << reso << endl;

}

// Method to interpolate the values and calculate norms of the error for resolution N.
void simulateScalar(int N, Gnuplot<T>& gplot)
{
  // === 1st Step: Initialization ===
  OstreamManager clout( std::cout,"scalar" );

  UnitConverterFromResolutionAndLatticeVelocity<T, DESCRIPTOR> converter(
    int {N},                  // resolution: number of voxels per charPhysL
    (T)   0.05,                // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   lx,           // charPhysLength: reference length of simulation geometry
    (T)   1.0,              // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1.0,                 // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0             // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("interpolation3D");

  // === 2nd Step: Prepare Geometry ===

  // Define cuboid with center in (0,0,0).
  Vector<T, 3> origin(-0.5*lx, -0.5*ly, -0.5*lz);
  Vector<T, 3> extend(lx, ly, lz);
  IndicatorCuboid3D<T> rec(extend, origin);

  // Define ball inside cuboid with center (0,0,0).
  Vector<T, 3> originB(0., 0., 0.);
  IndicatorSphere3D<T> ball(originB, r);

  // Instantiation of a cuboidDecomposition with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else // ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = 1;
#endif // ifdef PARALLEL_MODE_MPI
  CuboidDecomposition3D<T> cuboidDecomposition( rec, converter.getPhysDeltaX(), noOfCuboids);

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidDecomposition);

  // Instantiation of a superGeometry
  const int overlap = 2;
  SuperGeometry<T,3> superGeometry(cuboidDecomposition, loadBalancer, overlap);
  prepareGeometry(converter, superGeometry, ball);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( converter, superGeometry );

  //prepareLattice and setBoundaryConditions
  prepareLattice(sLattice, converter, superGeometry);

  // === 4th Step: Main ===
  clout << "starting simulation..." << std::endl;
  clout << "dx = " << converter.getPhysDeltaX() << std::endl;

  getResults(N, sLattice, converter, superGeometry, gplot);

}

// Main
int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  initialize( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  static Gnuplot<T> gplot( "plot", Gnuplot<T>::LINEAR );
  gplot.setLabel("Resolution", "Error");

  // Define maximus number of resolution N.
  T maxN = N+91;

  // loop over the different simulations
  for(int simuN = N; simuN < maxN; simuN = simuN + 20){

    // Run the simulations
    clout << "Starting next simulation with N = " << simuN << std::endl;
    simulateScalar(simuN, gplot); // Method
  }

  gplot.writePNG();
  clout << "Done " << std::endl;

}
