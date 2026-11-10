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
using namespace olb::descriptors;
using namespace olb::particles;
typedef double T;

//#define VERBOSE_TESTS

TEST(TestParticleFunctor, test2DPosition)
{
  const int D = 2;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>> PARTICLETYPE;
  typedef GENERAL GROUP;
  typedef POSITION FIELD;
  const int noP = 7;
  const int iP = 4;
  const int expectedTargetDim = 2;

  //Create particle system
  ParticleSystem<T,PARTICLETYPE> particleSystem(noP);

  //Write some data
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP){
    particleSystem.getFieldPointer<GROUP,FIELD>(iP)
      = Vector<T,expectedTargetDim>(2.*iP);
  }

  //Create particle position functor
  GroupedFieldF<T,PARTICLETYPE,GROUP,FIELD> groupedFieldF(
    particleSystem.get() );

  //Read Data
  const int Df = groupedFieldF.getTargetDim();
  T output[Df];
  for(int iDim=0; iDim<Df; ++iDim){ output[iDim] = 0.; }
  int input[1] = {iP};
  groupedFieldF( output, input );
  Vector<T,expectedTargetDim> value(output);

#ifdef VERBOSE_TESTS
  std::cout << "FunctorTargetDim= " << groupedFieldF.getTargetDim() << std::endl;
  std::cout << "ContainerSize=    " << groupedFieldF.getContainerSize() << std::endl;
  std::cout << "iP=" << iP << ", value=" << value << std::endl;
#endif

  //Test general
  EXPECT_EQ(groupedFieldF.getTargetDim(), expectedTargetDim);
  EXPECT_EQ(groupedFieldF.getContainerSize(), noP);

  //Test setting specifc
  Vector<T,expectedTargetDim> valZero(0.);
  EXPECT_GE(value,valZero);


}


TEST(TestParticleFunctor, test2DRotationMatrix)
{
  const int D = 2;
  typedef PARTICLE_DESCRIPTOR<D,SURFACE_RESOLVED<D>> PARTICLETYPE;
  typedef SURFACE GROUP;
  typedef ROT_MATRIX FIELD;
  const int noP = 10;
  const int iP = 5;
  const int expectedTargetDim = 4;

  //Create particle system
  ParticleSystem<T,PARTICLETYPE> particleSystem(noP);

  //Write some data
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP){
    particleSystem.getFieldPointer<GROUP,FIELD>(iP)
      = Vector<T,expectedTargetDim>(2.*iP);
  }

  //Create particle position functor
  GroupedFieldF<T,PARTICLETYPE,GROUP,FIELD> groupedFieldF(
    particleSystem.get() );

  //Read Data
  const int Df = groupedFieldF.getTargetDim();
  T output[Df];
  for(int iDim=0; iDim<Df; ++iDim){ output[iDim] = 0.; }
  int input[1] = {iP};
  groupedFieldF( output, input );
  Vector<T,expectedTargetDim> value(output);

#ifdef VERBOSE_TESTS
  std::cout << "FunctorTargetDim= " << groupedFieldF.getTargetDim() << std::endl;
  std::cout << "ContainerSize=    " << groupedFieldF.getContainerSize() << std::endl;
  std::cout << "iP=" << iP << ", value=" << value << std::endl;
#endif

  //Test general
  EXPECT_EQ(groupedFieldF.getTargetDim(), expectedTargetDim);
  EXPECT_EQ(groupedFieldF.getContainerSize(), noP);

  //Test setting specifc
  Vector<T,expectedTargetDim> valZero(0.);
  EXPECT_GE(value,valZero);
}

TEST(TestParticleFunctor, test3DPosition)
{
  const int D = 3;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>> PARTICLETYPE;
  typedef GENERAL GROUP;
  typedef POSITION FIELD;
  const int noP = 5;
  const int iP = 3;
  const int expectedTargetDim = 3;

  //Create particle system
  ParticleSystem<T,PARTICLETYPE> particleSystem(noP);

  //Write some data
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP){
    particleSystem.getFieldPointer<GROUP,FIELD>(iP)
      = Vector<T,expectedTargetDim>(2.*iP);
  }

  //Create particle position functor
  GroupedFieldF<T,PARTICLETYPE,GROUP,FIELD> groupedFieldF(
    particleSystem.get() );

  //Read Data
  const int Df = groupedFieldF.getTargetDim();
  T output[Df];
  for(int iDim=0; iDim<Df; ++iDim){ output[iDim] = 0.; }
  int input[1] = {iP};
  groupedFieldF( output, input );
  Vector<T,expectedTargetDim> value(output);

#ifdef VERBOSE_TESTS
  std::cout << "FunctorTargetDim= " << groupedFieldF.getTargetDim() << std::endl;
  std::cout << "ContainerSize=    " << groupedFieldF.getContainerSize() << std::endl;
  std::cout << "iP=" << iP << ", value=" << value << std::endl;
#endif

  //Test general
  EXPECT_EQ(groupedFieldF.getTargetDim(), expectedTargetDim);
  EXPECT_EQ(groupedFieldF.getContainerSize(), noP);

  //Test setting specifc
  Vector<T,expectedTargetDim> valZero(0.);
  EXPECT_GE(value,valZero);
}


TEST(TestParticleFunctor, test3DRotationMatrix)
{
  const int D = 3;
  typedef PARTICLE_DESCRIPTOR<D,SURFACE_RESOLVED<D>> PARTICLETYPE;
  typedef SURFACE GROUP;
  typedef ROT_MATRIX FIELD;
  const int noP = 3;
  const int iP = 2;
  const int expectedTargetDim = 9;

  //Create particle system
  ParticleSystem<T,PARTICLETYPE> particleSystem(noP);

  //Write some data
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP){
    particleSystem.getFieldPointer<GROUP,FIELD>(iP)
      = Vector<T,expectedTargetDim>(2.*iP);
  }

  //Create particle position functor
  GroupedFieldF<T,PARTICLETYPE,GROUP,FIELD> groupedFieldF(
    particleSystem.get() );

  //Read Data
  const int Df = groupedFieldF.getTargetDim();
  T output[Df];
  for(int iDim=0; iDim<Df; ++iDim){ output[iDim] = 0.; }
  int input[1] = {iP};
  groupedFieldF( output, input );
  Vector<T,expectedTargetDim> value(output);

#ifdef VERBOSE_TESTS
  std::cout << "FunctorTargetDim= " << groupedFieldF.getTargetDim() << std::endl;
  std::cout << "ContainerSize=    " << groupedFieldF.getContainerSize() << std::endl;
  std::cout << "iP=" << iP << ", value=" << value << std::endl;
#endif

  //Test general
  EXPECT_EQ(groupedFieldF.getTargetDim(), expectedTargetDim);
  EXPECT_EQ(groupedFieldF.getContainerSize(), noP);

  //Test setting specifc
  Vector<T,expectedTargetDim> valZero(0.);
  EXPECT_GE(value,valZero);
}


TEST(TestParticleFunctor, test2DCircumRadius)
{
  const int D = 2;
  typedef PARTICLE_DESCRIPTOR<D,SURFACE_RESOLVED<D>> PARTICLETYPE;
  const int noP = 1;
  const int expectedTargetDim = 1;
  const double radius = 0.5;
  const double epsilon = .05*radius;
  int iP = 0;

  //Create particle system
  ParticleSystem<T,PARTICLETYPE> particleSystem(noP);

  //Defining smooth indicator
  using namespace descriptors;
  typedef SmoothIndicatorF2D<T,T,true> SIndicatorBaseType;
  //Create SmoothIndicator
  std::unique_ptr<SIndicatorBaseType> sIndicatorPtr(
    new SmoothIndicatorCircle2D<T, T, true>(Vector<T,2>(0.), radius, epsilon ));
  //Pass smart pointer to particleSystem
  auto& vectorOfIndicators = particleSystem.template getAssociatedData<
                             std::vector<std::unique_ptr<SIndicatorBaseType>>>();
  //Add to field
  auto particle = particleSystem.get(iP);
  particle.template setField<SURFACE,SINDICATOR>( sIndicatorPtr.get() );
  vectorOfIndicators.push_back( std::move(sIndicatorPtr) );

  //Define functor
  ParticleCircumRadiusF<T,PARTICLETYPE> particleCircumRadiusF(
    particleSystem.get() );

  //Read Data
  const int Df = particleCircumRadiusF.getTargetDim();
  T output[Df];
  int input[1] = {iP};
  particleCircumRadiusF( output, input );
  T value = output[0];

#ifdef VERBOSE_TESTS
  std::cout << "FunctorTargetDim= " << particleCircumRadiusF.getTargetDim() << std::endl;
  std::cout << "ContainerSize=    " << particleCircumRadiusF.getContainerSize() << std::endl;
  std::cout << "iP=" << iP << ", value=" << value << std::endl;
#endif

  //Test general
  EXPECT_EQ(particleCircumRadiusF.getTargetDim(), expectedTargetDim);
  EXPECT_EQ(particleCircumRadiusF.getContainerSize(), noP);

  //Test setting specifc
  T circumRadiusExpected = radius+.5*epsilon;
  EXPECT_EQ(value,circumRadiusExpected);

}


TEST(TestParticleFunctor, test3DCircumRadius)
{
  const int D = 3;
  typedef PARTICLE_DESCRIPTOR<D,SURFACE_RESOLVED<D>> PARTICLETYPE;
  const int noP = 4;
  const int expectedTargetDim = 1;
  const double xLength = 1.6;
  const double yLength = 2.3;
  const double zLength = 1.0;
  const double epsilon = .2*zLength;
  int iP = 3;

  //Create particle system
  ParticleSystem<T,PARTICLETYPE> particleSystem(noP);

  //Defining smooth indicator
  using namespace descriptors;
  typedef SmoothIndicatorF3D<T,T,true> SIndicatorBaseType;
  //Create SmoothIndicator
  std::unique_ptr<SIndicatorBaseType> sIndicatorPtr(
    new SmoothIndicatorCuboid3D<T, T, true>( xLength,yLength,zLength,
      Vector<T,3>(0.), epsilon ));
  //Pass smart pointer to particleSystem
  auto& vectorOfIndicators = particleSystem.template getAssociatedData<
                             std::vector<std::unique_ptr<SIndicatorBaseType>>>();
  //Add to field
  auto particle = particleSystem.get(iP);
  particle.template setField<SURFACE,SINDICATOR>( sIndicatorPtr.get() );
  vectorOfIndicators.push_back( std::move(sIndicatorPtr) );

  //Define functor
  ParticleCircumRadiusF<T,PARTICLETYPE> particleCircumRadiusF(
    particleSystem.get() );

  //Read Data
  const int Df = particleCircumRadiusF.getTargetDim();
  T output[Df];
  int input[1] = {iP};
  particleCircumRadiusF( output, input );
  T value = output[0];

#ifdef VERBOSE_TESTS
  std::cout << "FunctorTargetDim= " << particleCircumRadiusF.getTargetDim() << std::endl;
  std::cout << "ContainerSize=    " << particleCircumRadiusF.getContainerSize() << std::endl;
  std::cout << "iP=" << iP << ", value=" << value << std::endl;
#endif

  //Test general
  EXPECT_EQ(particleCircumRadiusF.getTargetDim(), expectedTargetDim);
  EXPECT_EQ(particleCircumRadiusF.getContainerSize(), noP);

  //Test setting specifc
  T circumRadiusExpected = .5*(util::sqrt(xLength*xLength+yLength*yLength+zLength*zLength))+0.5*epsilon;
  EXPECT_EQ(value,circumRadiusExpected);

}
