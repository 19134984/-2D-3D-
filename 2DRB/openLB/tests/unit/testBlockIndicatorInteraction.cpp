/*  This file is part of the OpenLB library
*
*  Copyright (C) 2021 Nicolas Hafen
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

#include <type_traits>

#include <gtest/gtest.h>
#include "helper/unitTestHelper.h"

#include <olb.h>

using namespace olb;
using namespace descriptors;
typedef double T;

//#define VERBOSE_TESTS

/// Create BlockGeometry environment
template<typename T, unsigned D>
BlockGeometry<T,D> createBlockGeometry( PhysR<T,D> origin, PhysR<T,D> extent, T delX,
                                        int padding = 2, int noCuboids=1 )
{
  int iCglob = 0;
  IndicatorCuboid<T,D> cuboidIndicator( extent, origin );
  CuboidDecomposition<T,D> cuboidDecomposition( cuboidIndicator, delX, noCuboids );
  auto& cuboid = cuboidDecomposition.get(iCglob);
  BlockGeometry<T,D> blockGeometry( cuboid, padding, iCglob);
  return blockGeometry;
}

template<typename T, typename DESCRIPTOR>
std::shared_ptr<BlockLattice<T,DESCRIPTOR>> createBlockLattice( PhysR<T,DESCRIPTOR::d> origin, PhysR<T,DESCRIPTOR::d> extent,
    T delX, int padding = 2, int noCuboids=1 )
{
  int iCglob = 0;
  IndicatorCuboid<T,DESCRIPTOR::d> cuboidIndicator( extent, origin );
  CuboidDecomposition<T,DESCRIPTOR::d> cuboidDecomposition( cuboidIndicator, delX, noCuboids );
  auto& cuboid = cuboidDecomposition.get(iCglob);
  auto& extendedCuboid = cuboid.getExtent();
  std::shared_ptr<BlockLattice<T,DESCRIPTOR>> blockLatticePtr( new ConcreteBlockLattice<T,DESCRIPTOR>( extendedCuboid, padding ));
  return blockLatticePtr;
}


/// 3D Versions

TEST(TestSetBlockParticleField, 3DsmallRadiusComprehensive)
{
  //Settings
  typedef D3Q19<VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR,POROSITY> DESCRIPTOR;
  constexpr unsigned D = DESCRIPTOR::d;
  PhysR<T,D> origin (0.5);
  T baseSize = 1.;
  T delX = 0.1;
  int padding = 2;
  PhysR<T,D> indicatorPosition( origin );
  T indicatorCircumRadius = 0.2*baseSize;
  T epsilon = .5*delX;

  //Create BlockGeometry
  PhysR<T,D> extent (baseSize);
  BlockGeometry<T,D> blockGeometry = createBlockGeometry<T,D>( origin, extent, delX, padding );

  //Create UnitConverter
  UnitConverter<T,DESCRIPTOR> converter( delX, delX, baseSize, 1., 1., 1.);

  //Create eccentric velocity field
  Vector<T,D> velocity( 5.0, 6.0, 7.0 );
  Vector<T,utilities::dimensions::convert<D>::rotation> angularVelocity( 0. );
  EccentricLatticeVelocityField<T,T,DESCRIPTOR> eccentricVelocity( indicatorPosition,
      velocity, angularVelocity, converter );

  //Create SmoothIndicator
  T radius = indicatorCircumRadius-.5*epsilon;
  Vector<T,utilities::dimensions::convert<D>::rotation> angle( 0. );
  IndicatorSphere3D<T> sphere(indicatorPosition, radius);
  std::shared_ptr<SmoothIndicatorF<T,T,D,true>> sIndicatorPtr(
        new SmoothIndicatorSphere3D<T,T,true>(sphere, epsilon ));

  //Initialize Fields
  BlockIndicatorMaterial3D<T> blockIndicator( blockGeometry, 0 );
  AnalyticalConst<D,T,T> one(1.);

//  Vector<T,utilities::dimensions::convert<D>::matrix> rotationMatrix =
//    util::calculateRotationMatrix<T,D>( angle );

//  //Create BlockLattice
//  BlockLattice<T,DESCRIPTOR>& blockLattice = *createBlockLattice<T,DESCRIPTOR>( origin, extent, delX, padding );

//  blockLattice.defineField<POROSITY>( blockIndicator, one ); //TODO: why does this not work anymore?

//  //Set block particle field
//  setBlockParticleField( blockGeometry, eccentricVelocity, *sIndicatorPtr, blockLattice,
//    indicatorPosition, rotationMatrix, indicatorCircumRadius );
//
//  //Setup Test values
//  Vector<T,D> fieldVelNumerator;
//  T fieldVelDenominator;
//  T fieldPorosity;
//
//  //Test inside
//  LatticeR<D> testPosLatticeIn(0);
//  fieldVelNumerator = blockLattice.get(testPosLatticeIn).template getField<descriptors::VELOCITY_NUMERATOR>();
//  fieldVelDenominator = blockLattice.get(testPosLatticeIn).template getField<descriptors::VELOCITY_DENOMINATOR>();
//  fieldPorosity = blockLattice.get(testPosLatticeIn).template getFieldPointer<descriptors::POROSITY>()[0];

//#ifdef VERBOSE_TESTS
//  std::cout << "VelocityNumeratorIn=   " << fieldVelNumerator << std::endl;
//  std::cout << "VelocityDenominatorIn= " << fieldVelDenominator << std::endl;
//  std::cout << "PorosityIn=            " << fieldPorosity << std::endl;
//#endif
//
//  EXPECT_EQ(fieldVelNumerator, velocity);
//  EXPECT_EQ(fieldVelDenominator, 1.);
//  EXPECT_EQ(fieldPorosity, 0.);
//
//
//  //Test outside
//  LatticeR<D> testPosLatticeOut(baseSize/delX);
//  fieldVelNumerator = blockLattice.get(testPosLatticeOut).template getField<descriptors::VELOCITY_NUMERATOR>();
//  fieldVelDenominator = blockLattice.get(testPosLatticeOut).template getField<descriptors::VELOCITY_DENOMINATOR>();
//  fieldPorosity = blockLattice.get(testPosLatticeOut).template getFieldPointer<descriptors::POROSITY>()[0];
//
//#ifdef VERBOSE_TESTS
//  std::cout << "VelocityNumeratorOut=   " << fieldVelNumerator << std::endl;
//  std::cout << "VelocityDenominatorOut= " << fieldVelDenominator << std::endl;
//  std::cout << "PorosityOut=            " << fieldPorosity << std::endl;
//#endif
//
//  Vector<T,D> velZero(0.);
//  EXPECT_EQ(fieldVelNumerator, velZero);
//  EXPECT_EQ(fieldVelDenominator, 0.);
//  EXPECT_EQ(fieldPorosity, 1.);
//
//  //Test border
//  LatticeR<D> testPosLatticeBorder(-1);
//  fieldVelNumerator = blockLattice.get(testPosLatticeBorder).template getField<descriptors::VELOCITY_NUMERATOR>();
//  fieldVelDenominator = blockLattice.get(testPosLatticeBorder).template getField<descriptors::VELOCITY_DENOMINATOR>();
//  fieldPorosity = blockLattice.get(testPosLatticeBorder).template getFieldPointer<descriptors::POROSITY>()[0];
//
//#ifdef VERBOSE_TESTS
//  std::cout << "VelocityNumeratorBorder=   " << fieldVelNumerator << std::endl;
//  std::cout << "VelocityDenominatorBorder= " << fieldVelDenominator << std::endl;
//  std::cout << "PorosityBorder=            " << fieldPorosity << std::endl;
//#endif
//
//  EXPECT_GT(fieldVelNumerator, velZero);
//  EXPECT_LT(fieldVelNumerator, velocity);
//  EXPECT_GT(fieldVelDenominator, 0.);
//  EXPECT_LT(fieldVelDenominator, 1.);
//  EXPECT_GT(fieldPorosity, 0.);
//  EXPECT_LT(fieldPorosity, 1.);
}



TEST(TestSetBlockParticleField, 2DsmallRadiusComprehensive)
{
  //Settings
  typedef D2Q9<VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR,POROSITY> DESCRIPTOR;
  constexpr unsigned D = DESCRIPTOR::d;
  PhysR<T,D> origin (0.5);
  T baseSize = 1.;
  T delX = 0.1;
  int padding = 2;
  PhysR<T,D> indicatorPosition( origin );
  T indicatorCircumRadius = 0.3*baseSize;
  T epsilon = .5*delX;

  //Create BlockGeometry
  PhysR<T,D> extent (baseSize);
  BlockGeometry<T,D> blockGeometry = createBlockGeometry<T,D>( origin, extent, delX, padding );

  //Create UnitConverter
  UnitConverter<T,DESCRIPTOR> converter( delX, delX, baseSize, 1., 1., 1.);

  //Create eccentric velocity field
  Vector<T,D> velocity( 8., 9. );
  Vector<T,utilities::dimensions::convert<D>::rotation> angularVelocity( 0. );
  EccentricLatticeVelocityField<T,T,DESCRIPTOR> eccentricVelocity( indicatorPosition,
      velocity, angularVelocity, converter );

  //Create SmoothIndicator
  T radius = indicatorCircumRadius-.5*epsilon;
  Vector<T,utilities::dimensions::convert<D>::rotation> angle( 0. );
  IndicatorCircle2D<T> circle(indicatorPosition, radius);
  std::shared_ptr<SmoothIndicatorF<T,T,D,true>> sIndicatorPtr(
        new SmoothIndicatorCircle2D<T,T,true>(circle, epsilon ));

  //Initialize fields
  BlockIndicatorMaterial2D<T> blockIndicator( blockGeometry, 0);
  AnalyticalConst<D,T,T> one(1.);

//  Vector<T,utilities::dimensions::convert<D>::matrix> rotationMatrix =
//    util::calculateRotationMatrix<T,D>( angle );

//  //Create BlockLattice
//  BlockLattice<T,DESCRIPTOR>& blockLattice = *createBlockLattice<T,DESCRIPTOR>( origin, extent, delX, padding );

//  blockLattice.defineField<POROSITY>( blockIndicator, one );

//  //Set block particle field
//  setBlockParticleField( blockGeometry, eccentricVelocity, *sIndicatorPtr, blockLattice,
//    indicatorPosition, rotationMatrix, indicatorCircumRadius );

//  //Setup Test values
//  Vector<T,D> fieldVelNumerator;
//  T fieldVelDenominator;
//  T fieldPorosity;

//  //Test inside
//  LatticeR<D> testPosLatticeIn(0);
//  fieldVelNumerator = blockLattice.get(testPosLatticeIn).template getField<descriptors::VELOCITY_NUMERATOR>();
//  fieldVelDenominator = blockLattice.get(testPosLatticeIn).template getField<descriptors::VELOCITY_DENOMINATOR>();
//  fieldPorosity = blockLattice.get(testPosLatticeIn).template getFieldPointer<descriptors::POROSITY>()[0];

//#ifdef VERBOSE_TESTS
//  std::cout << "VelocityNumeratorIn=   " << fieldVelNumerator << std::endl;
//  std::cout << "VelocityDenominatorIn= " << fieldVelDenominator << std::endl;
//  std::cout << "PorosityIn=            " << fieldPorosity << std::endl;
//#endif

//  EXPECT_EQ(fieldVelNumerator, velocity);
//  EXPECT_EQ(fieldVelDenominator, 1.);
//  EXPECT_EQ(fieldPorosity, 0.);

//  //Test outside
//  LatticeR<D> testPosLatticeOut(baseSize/delX);
//  fieldVelNumerator = blockLattice.get(testPosLatticeOut).template getField<descriptors::VELOCITY_NUMERATOR>();
//  fieldVelDenominator = blockLattice.get(testPosLatticeOut).template getField<descriptors::VELOCITY_DENOMINATOR>();
//  fieldPorosity = blockLattice.get(testPosLatticeOut).template getFieldPointer<descriptors::POROSITY>()[0];

//#ifdef VERBOSE_TESTS
//  std::cout << "VelocityNumeratorOut=   " << fieldVelNumerator << std::endl;
//  std::cout << "VelocityDenominatorOut= " << fieldVelDenominator << std::endl;
//  std::cout << "PorosityOut=            " << fieldPorosity << std::endl;
//#endif

//  Vector<T,D> velZero(0.);
//  EXPECT_EQ(fieldVelNumerator, velZero);
//  EXPECT_EQ(fieldVelDenominator, 0.);
//  EXPECT_EQ(fieldPorosity, 1.);

//  //Test border
//  LatticeR<D> testPosLatticeBorder(-2,-2);
//  fieldVelNumerator = blockLattice.get(testPosLatticeBorder).template getField<descriptors::VELOCITY_NUMERATOR>();
//  fieldVelDenominator = blockLattice.get(testPosLatticeBorder).template getField<descriptors::VELOCITY_DENOMINATOR>();
//  fieldPorosity = blockLattice.get(testPosLatticeBorder).template getFieldPointer<descriptors::POROSITY>()[0];

//#ifdef VERBOSE_TESTS
//  std::cout << "VelocityNumeratorBorder=   " << fieldVelNumerator << std::endl;
//  std::cout << "VelocityDenominatorBorder= " << fieldVelDenominator << std::endl;
//  std::cout << "PorosityBorder=            " << fieldPorosity << std::endl;
//#endif

//  EXPECT_GT(fieldVelNumerator, velZero);
//  EXPECT_LT(fieldVelNumerator, velocity);
//  EXPECT_GT(fieldVelDenominator, 0.);
//  EXPECT_LT(fieldVelDenominator, 1.);
//  EXPECT_GT(fieldPorosity, 0.);
//  EXPECT_LT(fieldPorosity, 1.);


}
