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
typedef double T;
const T eps = std::numeric_limits<T>::epsilon();

//#define VERBOSE_TESTS




///SCALARS

TEST(TestRandom, scalarRandom)
{
  util::Randomizer<T> randomizer;
  auto rndmOutput = randomizer.generate<T>();
#ifdef VERBOSE_TESTS
  std::cout << rndmOutput << std::endl;
#endif

  //Test general
  EXPECT_LE(rndmOutput, 1.+eps);
  EXPECT_GE(rndmOutput, 0.-eps);
}

TEST(TestRandom, scalarStored)
{
  std::vector<T> sequence = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8};
  util::Randomizer<T,true> randomizer( sequence );
  auto rndmOutput = randomizer.generate<T>();
#ifdef VERBOSE_TESTS
  std::cout << rndmOutput << std::endl;
#endif

  //Test general
  EXPECT_LE(rndmOutput, 1.+eps);
  EXPECT_GE(rndmOutput, 0.-eps);
  EXPECT_EQ(rndmOutput,sequence[0]);
}


TEST(TestRandom, scalarStoredFile)
{
  std::string filePathSequence = "./testSequence.dat";

  //Write sequence
  util::Randomizer<T,false> randomizerA;
  randomizerA.writeSequence(5,filePathSequence);

  //Create randomizer from written sequence
  util::Randomizer<T,true> randomizer(filePathSequence);
  std::remove("testSequence.dat");
  auto rndmOutput = randomizer.generate<T>();
#ifdef VERBOSE_TESTS
  std::cout << rndmOutput << std::endl;
#endif

  //Test general
  EXPECT_LE(rndmOutput, 1.+eps);
  EXPECT_GE(rndmOutput, 0.-eps);
}


TEST(TestRandom, scalarRandomNormal)
{
  T avg = 1.5;
  T stdDev = 0.3;
  util::Randomizer<T> randomizer;
  auto rndmOutput = randomizer.generateScalarNormal(avg, stdDev);
#ifdef VERBOSE_TESTS
  std::cout << rndmOutput << std::endl;
#endif
}








///VECTORS

TEST(TestRandom, vectorRandom)
{
  const int D = 3;
  util::Randomizer<T> randomizer;
  auto rndmOutput = randomizer.generate<Vector<T,D>>();
#ifdef VERBOSE_TESTS
  std::cout << rndmOutput << std::endl;
#endif

  Vector<T,D> maxExp(1.+eps);
  Vector<T,D> minExp(0.-eps);

  //Test general
  EXPECT_LE(rndmOutput, maxExp);
  EXPECT_GE(rndmOutput, minExp);
}

TEST(TestRandom, vectorStored)
{
  const int D = 3;
  std::vector<T> sequence = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8};
  util::Randomizer<T,true> randomizer( sequence );
  auto rndmOutput = randomizer.generate<Vector<T,D>>();
#ifdef VERBOSE_TESTS
  std::cout << rndmOutput << std::endl;
#endif

  Vector<T,D> maxExp(1.+eps);
  Vector<T,D> minExp(0.-eps);
  Vector<T,D> vecExp({sequence.begin(), sequence.begin() + D});

  //Test general
  EXPECT_LE(rndmOutput, maxExp);
  EXPECT_GE(rndmOutput, minExp);
  EXPECT_EQ(rndmOutput, vecExp);
}


///HANDLING


TEST(TestRandom, nonExistingFileSequence)
{
  const int D = 3;
  std::string fileName = "./testSequence.dat";
  bool enforceStored = false;
  util::Randomizer<T,true> randomizer(fileName,enforceStored);

  auto rndmOutput = randomizer.generate<Vector<T,D>>();
#ifdef VERBOSE_TESTS
  std::cout << rndmOutput << std::endl;
#endif

  std::remove("testSequence.dat");

  Vector<T,D> maxExp(1.+eps);
  Vector<T,D> minExp(0.-eps);

  //Test general
  EXPECT_LE(rndmOutput, maxExp);
  EXPECT_GE(rndmOutput, minExp);
}
