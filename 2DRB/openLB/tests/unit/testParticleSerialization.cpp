/*  This file is part of the OpenLB library
*
*  Copyright (C) 2022 Nicolas Hafen
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
using namespace olb::particles::dynamics;
using namespace olb::particles::communication;
typedef double T;

//#define VERBOSE_TESTS

TEST(TestParticleUtilities, testMultiFieldArrayCommunicatable)
{
  using DESCRIPTOR = descriptors::D2Q9<>;
  MultiFieldArrayD<T,DESCRIPTOR,Platform::CPU_SISD,descriptors::FORCE,descriptors::OMEGA> data(10);

  std::vector<CellID> indices;
  indices.emplace_back(5);
  indices.emplace_back(8);

  data.template get<descriptors::OMEGA>()[0][5] = 1.0;
  data.template get<descriptors::OMEGA>()[0][8] = 2.0;

  std::size_t size = ConcreteCommunicatable(data).size(indices);
  EXPECT_EQ(size, std::size_t{48});

  std::unique_ptr<std::uint8_t[]> buffer(new std::uint8_t[size] { });

  ConcreteCommunicatable(data).serialize(indices, buffer.get());

  MultiFieldArrayD<T,DESCRIPTOR,Platform::CPU_SISD,descriptors::FORCE,descriptors::OMEGA> data2(10);

  ConcreteCommunicatable(data2).deserialize(indices, buffer.get());

  EXPECT_EQ(data2.template get<descriptors::OMEGA>()[0][5], 1.0);
  EXPECT_EQ(data2.template get<descriptors::OMEGA>()[0][8], 2.0);
}


TEST(TestParticleUtilities, testDynamicFieldGroupsCommunicatable)
{
  using DESCRIPTOR = ResolvedParticle3D;
  using GROUP_A = GENERAL_TMP<3>;
  using FIELD_A = POSITION;
  using GROUP_B = PHYSPROPERTIES_RESOLVED<3>;
  using FIELD_B = MASS;
  using DATA = DynamicFieldGroupsD<T, typename DESCRIPTOR::fields_t>;

  //Create two empty DynamicFieldGroups
  DATA data(10);
  DATA data2(10);

  //Define indices to be communicated
  std::vector<CellID> indices;
  indices.emplace_back(5);
  indices.emplace_back(8);

  //Set some data at first index
  data.template get<GROUP_A>().template get<FIELD_A>()[0][indices[0]] = 1.0;
  data.template get<GROUP_A>().template get<FIELD_A>()[1][indices[0]] = 2.0;
  data.template get<GROUP_A>().template get<FIELD_A>()[2][indices[0]] = 3.0;
  data.template get<GROUP_B>().template get<FIELD_B>()[0][indices[0]] = 4.0;

  //Set some data at second index
  data.template get<GROUP_A>().template get<FIELD_A>()[0][indices[1]] = 31.0;
  data.template get<GROUP_A>().template get<FIELD_A>()[1][indices[1]] = 32.0;
  data.template get<GROUP_A>().template get<FIELD_A>()[2][indices[1]] = 33.0;
  data.template get<GROUP_B>().template get<FIELD_B>()[0][indices[1]] = 34.0;

  //Print data
#ifdef VERBOSE_TESTS
  particles::io::printGenericParticleInfo<T,DESCRIPTOR>( data, indices[0] );
  particles::io::printGenericParticleInfo<T,DESCRIPTOR>( data, indices[1] );
#endif

  //Evaluate size of GROUP_A
  std::size_t sizeA = ConcreteCommunicatable(data.template get<GROUP_A>()).size(indices);
  EXPECT_EQ(sizeA, std::size_t{48});

  //Evaluate size of GROUP_B
  std::size_t sizeB = ConcreteCommunicatable(data.template get<GROUP_B>()).size(indices);
  EXPECT_EQ(sizeB, std::size_t{64});

  //Define GroupedDataCommunicatable
  using Communicator = GroupedDataCommunicatableHelper<DATA,DESCRIPTOR>::type;

  //Evaluate size of data
  auto size = Communicator(data).size(indices);
  EXPECT_EQ(size, std::size_t{608});

#ifdef VERBOSE_TESTS
  std::cout << "Group_A size=" << sizeA << std::endl;
  std::cout << "Group_B size=" << sizeB<< std::endl;
  std::cout << "TOTAL size=  " << size << std::endl;
#endif

  //Create buffer with respective size
  std::unique_ptr<std::uint8_t[]> buffer(new std::uint8_t[size] { });

  //Serialize data and copy to buffer
  Communicator(data).serialize(indices, buffer.get());

  //Deserialize buffer and copy to data2
  Communicator(data2).deserialize(indices, buffer.get());

  //Print data2
#ifdef VERBOSE_TESTS
  particles::io::printGenericParticleInfo<T,DESCRIPTOR>( data2, indices[0] );
  particles::io::printGenericParticleInfo<T,DESCRIPTOR>( data2, indices[1] );
#endif

}

TEST(TestParticleUtilities, testCommunicatableFromParticleInterface)
{
  using PARTICLETYPE = ResolvedParticle3D;

  //Create particle systems
  ParticleSystem<T,PARTICLETYPE> particleSystem(10);
  ParticleSystem<T,PARTICLETYPE> particleSystem2(6);

  //Create particle at index idx
  std::size_t idx = 3;
  Vector<T,3> position(2., 3., 4.);
  T radius = 0.7;
  T epsilon = 0.1;
  creators::setResolvedCuboid3D<T,PARTICLETYPE>( particleSystem,
    idx, position, radius, epsilon );

  //Retrieve particle
  auto particle = particleSystem.get(idx);

#ifdef VERBOSE_TESTS
  particle.print();
#endif

  //Evaluate serialSize of data
  auto serialSize = particle.getSerialSize();
  EXPECT_EQ(serialSize, std::size_t{304});

#ifdef VERBOSE_TESTS
  std::cout << "SerialSize=" << serialSize << std::endl;
#endif

  //Create buffer
  std::unique_ptr<std::uint8_t[]> buffer(new std::uint8_t[serialSize] { });

  //Serialize particleA
  particle.serialize(buffer.get());

  //Deserialize into particleB
  auto particleB = particleSystem2.get(idx);
  particleB.deserialize(buffer.get());

#ifdef VERBOSE_TESTS
  particleB.print();
#endif

  Vector<T,3> positionB( particleB.getField<GENERAL,POSITION>() );
  EXPECT_EQ(position, positionB);
}


#define VERBOSE_TESTS


TEST(TestParticleUtilities, testSerialComposure)
{
  using PARTICLETYPE = ResolvedParticle3D;

  //Create particle systems
  ParticleSystem<T,PARTICLETYPE> particleSystem(10);
  ParticleSystem<T,PARTICLETYPE> particleSystem2(6);
  FieldD<T,PARTICLETYPE,IC> field(5);
  FieldD<T,PARTICLETYPE,IC> field2;

  //Create particle at index idx
  std::size_t idx = 3;
  Vector<T,3> position(2., 3., 4.);
  T radius = 0.7;
  T epsilon = 0.1;
  creators::setResolvedCuboid3D<T,PARTICLETYPE>( particleSystem,
    idx, position, radius, epsilon );

  //Retrieve particle
  auto particle = particleSystem.get(idx);

#ifdef VERBOSE_TESTS
  std::cout << "IC(at sending)=" << field[0] << std::endl;
  particle.print();
#endif

  //Evaluate serialSizes
  auto serialSizeParticle = particle.getSerialSize();
  EXPECT_EQ(serialSizeParticle, std::size_t{304});

  //Create communicatable for field
  ConcreteCommunicatable fieldCommunicatable(field);
  const std::vector<CellID> indicesTmp{0};
  auto serialSizeField = fieldCommunicatable.size(indicesTmp);

  //Calculate serial size of composure
  int serialSizeComposure = serialSizeParticle + serialSizeField;

#ifdef VERBOSE_TESTS
  std::cout << "Serialsize Particle=" << serialSizeParticle << std::endl;
  std::cout << "Serialsize Field=   " << serialSizeField << std::endl;
  std::cout << "SerialComposure=    " << serialSizeComposure << std::endl;
#endif

  //Create buffer
  std::unique_ptr<std::uint8_t[]> buffer(new std::uint8_t[serialSizeComposure] { });

  //Retrieve buffers pointers
  std::uint8_t* bufferRaw = buffer.get();
  std::uint8_t* bufferShifted = &bufferRaw[serialSizeField];

  //Serialize field and particle
  fieldCommunicatable.serialize(indicesTmp, bufferRaw);
  particle.serialize(bufferShifted);

  //Deserialize field and particle
  auto particle2 = particleSystem2.get(idx);
  particle2.deserialize(bufferShifted);
  ConcreteCommunicatable(field2).deserialize(indicesTmp, bufferRaw);

#ifdef VERBOSE_TESTS
  std::cout << "IC(at receiving)=" << field2[0] << std::endl;
  particle2.print();
#endif

  Vector<T,3> position2( particle2.getField<GENERAL,POSITION>() );
  EXPECT_EQ(position, position2);
  EXPECT_EQ(field, field2);

}



TEST(TestParticleUtilities, testSerialComposureFunctions)
{
  using PARTICLETYPE = ResolvedParticle3D;

  //Create particle systems
  ParticleSystem<T,PARTICLETYPE> particleSystem(10);
  ParticleSystem<T,PARTICLETYPE> particleSystem2(6);
  int globiCdest = 5;

  //Create particle at index idx
  std::size_t idx = 3;
  Vector<T,3> position(2., 3., 4.);
  T radius = 0.7;
  T epsilon = 0.1;
  creators::setResolvedCuboid3D<T,PARTICLETYPE>( particleSystem,
    idx, position, radius, epsilon );

  //Retrieve particle
  auto particle = particleSystem.get(idx);

#ifdef VERBOSE_TESTS
  std::cout << "IC(at sending)=" << globiCdest << std::endl;
  particle.print();
#endif

  //Evaluate serialSizes
  auto serialSizeParticle = particle.getSerialSize();
  EXPECT_EQ(serialSizeParticle, std::size_t{304});

  //Create communicatable for field
  FieldD<T,PARTICLETYPE,IC> field;
  ConcreteCommunicatable fieldCommunicatable(field);
  const std::vector<CellID> indicesTmp{0};
  auto serialSizeField = fieldCommunicatable.size(indicesTmp);

  //Calculate serial size of composure
  int serialSizeComposure = serialSizeParticle + serialSizeField;

#ifdef VERBOSE_TESTS
  std::cout << "Serialsize Particle=" << serialSizeParticle << std::endl;
  std::cout << "Serialsize Field=   " << serialSizeField << std::endl;
  std::cout << "SerialComposure=    " << serialSizeComposure << std::endl;
#endif

  //Create buffer
  std::unique_ptr<std::uint8_t[]> buffer(new std::uint8_t[serialSizeComposure] { });

  //Serialize field and particle
  serializeIcDestAndParticle(globiCdest,particle,buffer.get());

  //Deserialize field and particle
  auto particle2 = particleSystem2.get(idx);
  int globiCdest2;
  deserializeIcDestAndParticle(globiCdest2, particle2, buffer.get());


#ifdef VERBOSE_TESTS
  std::cout << "IC(at receiving)=" << globiCdest2 << std::endl;
  particle2.print();
#endif

  Vector<T,3> position2( particle2.getField<GENERAL,POSITION>() );
  EXPECT_EQ(position, position2);
  EXPECT_EQ(globiCdest, globiCdest2);

}


TEST(TestParticleUtilities, testSerializationParticleForce)
{
  using PARTICLETYPE = ResolvedParticle3D;
  const std::vector<unsigned int> indices{0};

  //Add fieldArray for ID
  int id = 7;
  FieldArrayD<T,PARTICLETYPE,Platform::CPU_SISD,ID> fieldArrayID(1);
  fieldArrayID.setField(0,id);
  auto communicatableID = ConcreteCommunicatable(fieldArrayID);
  auto serialSizeID = communicatableID.size(indices);
  auto dataID = fieldArrayID.getField(0);
  std::unique_ptr<std::uint8_t[]> bufferA(new std::uint8_t[serialSizeID] { });
  communicatableID.serialize( indices, bufferA.get() );

  FieldArrayD<T,PARTICLETYPE,Platform::CPU_SISD,ID> fieldArrayIDrec(1);
  auto communicatableIDrec = ConcreteCommunicatable(fieldArrayIDrec);
  communicatableIDrec.deserialize( indices, bufferA.get() );
  auto dataIDrec = fieldArrayIDrec.getField(0);

  //Add fieldArray for FORCE
  Vector<T,3> force(3.,5.,6.);
  FieldArrayD<T,PARTICLETYPE,Platform::CPU_SISD,FORCE> fieldArrayFORCE(1);
  fieldArrayFORCE.setField(0,force);
  auto communicatableFORCE = ConcreteCommunicatable(fieldArrayFORCE);
  auto serialSizeFORCE = communicatableFORCE.size(indices);
  auto dataFORCE = fieldArrayFORCE.getField(0);
  std::unique_ptr<std::uint8_t[]> bufferB(new std::uint8_t[serialSizeFORCE] { });
  communicatableFORCE.serialize( indices, bufferB.get() );

  FieldArrayD<T,PARTICLETYPE,Platform::CPU_SISD,FORCE> fieldArrayFORCErec(1);
  auto communicatableFORCErec = ConcreteCommunicatable(fieldArrayFORCErec);
  communicatableFORCErec.deserialize( indices, bufferB.get() );
  auto dataFORCErec = fieldArrayFORCErec.getField(0);

#ifdef VERBOSE_TESTS
  std::cout << "FieldID(size=" << serialSizeID << ")=" << dataID
            << ", deserialized=" << dataIDrec
            << std::endl;
  std::cout << "FieldFORCE(size=" << serialSizeFORCE << ")=" << dataFORCE
            << ", deserialized=" << dataFORCErec
            << std::endl;
#endif

  EXPECT_EQ(serialSizeID, 8);
  EXPECT_EQ(dataIDrec, id);
  EXPECT_EQ(serialSizeFORCE, 24);
  EXPECT_EQ(dataFORCErec, force);
}


TEST(TestParticleUtilities, testStdArrayCommunication)
{
  //Create data
  std::array<T,7> data = {-4.5,5,-6.33,9.6,4,1.002,4};

  //Create communuicatable
  auto arrayCommunicatable = ConcreteCommunicatable(data);
  std::size_t size = arrayCommunicatable.size();

  //Serialize data
  std::unique_ptr<std::uint8_t[]> buffer(new std::uint8_t[size] { });
  arrayCommunicatable.serialize(buffer.get());

  //Create data2 and deserialize buffer
  std::array<T,7> data2;
  auto arrayCommunicatable2 = ConcreteCommunicatable(data2);
  arrayCommunicatable2.deserialize(buffer.get());

  //Check bufer size and correct values
  EXPECT_EQ(size, std::size_t{56});
  EXPECT_EQ(data2, data);
}
