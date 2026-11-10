/*  This file is part of the OpenLB library
*
*  Copyright (C) 2021 Jan E. Marquardt
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

#include <iostream>
#include <gtest/gtest.h>
#include "helper/unitTestHelper.h"

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::particles;

typedef double T;
typedef ResolvedParticle3D PARTICLETYPE;


//#define VERBOSE_TESTS

// test velocityVerletintegration 3D



//TEST(velocityVerletIntegration, test)
//{
//using namespace descriptors;
//constexpr unsigned D = PARTICLETYPE::d;
//  const unsigned Drot = utilities::dimensions::convert<D>::rotation;
//  const T timeStepSize = 1;


//}

// test velocityVerletintegration 3D

TEST(velocityVerletIntegration, test)
{
  //Settings
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  const unsigned Drot = utilities::dimensions::convert<D>::rotation;
  const T timeStepSize = 1;
  const T delX = 0.1;
  const T baseSize = 1.0;
  T particleCircumRadius = 0.2*baseSize;
  PhysR<T,D> particlePosition( .5*baseSize );

  //Create ParticleSystem with one particle
  ParticleSystem<T,PARTICLETYPE> particleSystem;

  //Create particle
  T epsilon = .5*delX;
  T radius = particleCircumRadius-.5*epsilon;
  T density = 1.0;
  particles::creators::addResolvedSphere3D( particleSystem, particlePosition, radius, epsilon, density );

  //Retrieve particle
  auto particle = particleSystem.get(0);

  //Define test quantities
  Vector<T,Drot> angularAcceleration(0.1, 0., 0.);
  Vector<T,Drot> acceleration(0.6, 0., 0.);

  //Perform velocity verlet integtation
  particles::dynamics::velocityVerletIntegration<T, PARTICLETYPE>(
    particle, timeStepSize, acceleration, angularAcceleration);

  //Retrieve new quantities
  Vector<T,D> position( particle.template getField<GENERAL,POSITION>() );
  Vector<T,D> velocity( particle.template getField<MOBILITY,VELOCITY>() );
  Vector<T,Drot> angle( particle.template getField<SURFACE,ANGLE>() );
  Vector<T,Drot> angularVelocity( particle.template getField<MOBILITY,ANG_VELOCITY>() );

#ifdef VERBOSE_TESTS
  std::cout <<  "Position:        " << position << std::endl;
  std::cout <<  "Velocity:        " << velocity << std::endl;
  std::cout <<  "Angle:           " << angle << std::endl;
  std::cout <<  "AngularVelocity: " << angularVelocity << std::endl;
#endif

  //Test quantities
  Vector<T,D> positionExpected(0.5, 0.5, 0.5);
  EXPECT_EQ( position, positionExpected);

  Vector<T,D> angularVelocityExpected(0.05, 0.0, 0.0);
  EXPECT_EQ( angularVelocity, angularVelocityExpected );

}
