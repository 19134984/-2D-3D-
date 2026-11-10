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


TEST(TestVTKwriter, test2D)
{
  const int D = 2;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>,PHYSPROPERTIES_RESOLVED<D>> PARTICLETYPE;
  const int noP = 5;

  //Create particle system
  ParticleSystem<T,PARTICLETYPE> particleSystem(noP);

  //Create particle writer
  GroupedFieldF<T,PARTICLETYPE,GENERAL,POSITION> particlePositionFunctor( particleSystem.get() );
  GroupedFieldF<T,PARTICLETYPE,PHYSPROPERTIES,MASS> particleMassF( particleSystem.get() );
  VTUwriter<T,PARTICLETYPE,false> particleVTKwriter("particles", false,false);
  particleVTKwriter.addFunctor(particlePositionFunctor, "Particles");
  particleVTKwriter.addFunctor(particleMassF, "Mass");
}


TEST(TestVTKwriter, test3D)
{
  const int D = 3;
  typedef PARTICLE_DESCRIPTOR<D,GENERAL_TMP<D>,PHYSPROPERTIES_RESOLVED<D>> PARTICLETYPE;
  const int noP = 12;

  //Create particle system
  ParticleSystem<T,PARTICLETYPE> particleSystem(noP);

  //Create particle writer
  GroupedFieldF<T,PARTICLETYPE,GENERAL,POSITION> particlePositionFunctor( particleSystem.get() );
  GroupedFieldF<T,PARTICLETYPE,PHYSPROPERTIES,MASS> particleMassF( particleSystem.get() );
  VTUwriter<T,PARTICLETYPE,false> particleVTKwriter("particles", false,false);
  particleVTKwriter.addFunctor(particlePositionFunctor, "Position");
  particleVTKwriter.addFunctor(particleMassF, "Mass");
}
