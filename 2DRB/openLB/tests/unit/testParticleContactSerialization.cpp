/*  This file is part of the OpenLB library
*
*  Copyright (C) 2022 Jan E. Marquardt
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

#include "helper/framework.h"
#include "helper/unitTestHelper.h"
#include <olb.h>
#include <gtest/gtest.h>
#include <iostream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::particles;
using namespace olb::particles::contact;

typedef double T;

TEST(ParticleContactSerialization, serializationAndDeserialization3d)
{
  constexpr unsigned                                            D = 3;
  typedef ParticleContactArbitraryFromOverlapVolume<T, D, true> CONTACTTYPE;

  std::array<std::size_t, 2> ids({5, 8});
  CONTACTTYPE                contact(ids);
  std::array<PhysR<T, D>, 2> positions;
  positions[0] = PhysR<T, D>(1, 2, 3);
  positions[1] = PhysR<T, D>(6, 5, 4);
  contact.setParticlePositions(positions);
  PhysR<T, D> min(-3);
  contact.updateMinMax(min);
  PhysR<T, D> max(7);
  contact.updateMinMax(max);
  contact.setParticlePositionUpdated(true);
  contact.setResponsibleRank(9);
  contact.isNew(false);

  std::size_t serialSize = contact.getSerialSize();
  EXPECT_EQ(126, serialSize);

  std::uint8_t* buffer = new std::uint8_t[serialSize] {};
  EXPECT_EQ(serialSize, contact.serialize(buffer));

  CONTACTTYPE newcontact;
  EXPECT_EQ(serialSize, newcontact.deserialize(buffer));

#ifdef VERBOSE
  std::cout << "IDs: " << newcontact.getIDs()[0] << std::endl;
  std::cout << "IDs: " << newcontact.getIDs()[1] << std::endl;
  std::cout << "positions: " << newcontact.getParticlePositions()[0]
            << std::endl;
  std::cout << "positions: " << newcontact.getParticlePositions()[1]
            << std::endl;
#endif

  EXPECT_TRUE(newcontact.isParticlePositionUpdated());
  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_EQ(positions[0][iD], newcontact.getParticlePositions()[0][iD]);
    EXPECT_EQ(positions[1][iD], newcontact.getParticlePositions()[1][iD]);
    EXPECT_EQ(min[iD], newcontact.getMin()[iD]);
    EXPECT_EQ(max[iD], newcontact.getMax()[iD]);
  }
  EXPECT_EQ(ids[0], newcontact.getIDs()[0]);
  EXPECT_EQ(ids[1], newcontact.getIDs()[1]);
  EXPECT_EQ(contact.getResponsibleRank(), newcontact.getResponsibleRank());
  EXPECT_EQ(contact.getDampingFactor(), newcontact.getDampingFactor());
  EXPECT_EQ(contact.isNew(), newcontact.isNew());
}

TEST(ParticleContactSerialization, serializationAndDeserialization2d)
{
  constexpr unsigned                                            D = 2;
  typedef ParticleContactArbitraryFromOverlapVolume<T, D, true> CONTACTTYPE;

  std::array<std::size_t, 2> ids({9, 2});
  CONTACTTYPE                contact(ids);
  std::array<PhysR<T, D>, 2> positions;
  positions[0] = PhysR<T, D>(1, 2);
  positions[1] = PhysR<T, D>(6, 5);
  contact.setParticlePositions(positions);
  PhysR<T, D> min(-13);
  contact.updateMinMax(min);
  PhysR<T, D> max(77);
  contact.updateMinMax(max);
  contact.setParticlePositionUpdated(false);
  contact.setResponsibleRank(2);
  contact.isNew(true);

  std::size_t serialSize = contact.getSerialSize();
  EXPECT_EQ(94, serialSize);

  std::uint8_t* buffer = new std::uint8_t[serialSize] {};
  EXPECT_EQ(serialSize, contact.serialize(buffer));

  CONTACTTYPE newcontact;
  EXPECT_EQ(serialSize, newcontact.deserialize(buffer));

#ifdef VERBOSE
  std::cout << "IDs: " << newcontact.getIDs()[0] << std::endl;
  std::cout << "IDs: " << newcontact.getIDs()[1] << std::endl;
  std::cout << "positions: " << newcontact.getParticlePositions()[0]
            << std::endl;
  std::cout << "positions: " << newcontact.getParticlePositions()[1]
            << std::endl;
#endif

  EXPECT_TRUE(!newcontact.isParticlePositionUpdated());
  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_EQ(positions[0][iD], newcontact.getParticlePositions()[0][iD]);
    EXPECT_EQ(positions[1][iD], newcontact.getParticlePositions()[1][iD]);
    EXPECT_EQ(min[iD], newcontact.getMin()[iD]);
    EXPECT_EQ(max[iD], newcontact.getMax()[iD]);
  }
  EXPECT_EQ(util::min(ids[0], ids[1]), newcontact.getIDs()[0]);
  EXPECT_EQ(util::max(ids[0], ids[1]), newcontact.getIDs()[1]);
  EXPECT_EQ(contact.getResponsibleRank(), newcontact.getResponsibleRank());
  EXPECT_EQ(contact.getDampingFactor(), newcontact.getDampingFactor());
  EXPECT_EQ(contact.isNew(), newcontact.isNew());
}

TEST(WallContactSerialization, serializationAndDeserialization3d)
{
  constexpr unsigned                                        D = 3;
  typedef WallContactArbitraryFromOverlapVolume<T, D, true> CONTACTTYPE;

  std::size_t id(8);
  unsigned    wallID = 2;
  CONTACTTYPE contact(id, wallID);
  PhysR<T, D> position;
  position = PhysR<T, D>(1, 2, 3);
  contact.setParticlePosition(position);
  PhysR<T, D> min(-3);
  contact.updateMinMax(min);
  PhysR<T, D> max(7);
  contact.updateMinMax(max);
  contact.setParticlePositionUpdated(true);
  contact.setResponsibleRank(13);
  contact.isNew(false);

  std::size_t serialSize = contact.getSerialSize();
  EXPECT_EQ(98, serialSize);

  std::uint8_t* buffer = new std::uint8_t[serialSize] {};
  EXPECT_EQ(serialSize, contact.serialize(buffer));

  CONTACTTYPE newcontact;
  EXPECT_EQ(serialSize, newcontact.deserialize(buffer));

  EXPECT_TRUE(newcontact.isParticlePositionUpdated());
  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_EQ(position[iD], newcontact.getParticlePosition()[iD]);
    EXPECT_EQ(min[iD], newcontact.getMin()[iD]);
    EXPECT_EQ(max[iD], newcontact.getMax()[iD]);
  }
  EXPECT_EQ(id, newcontact.getParticleID());
  EXPECT_EQ(wallID, newcontact.getWallID());
  EXPECT_EQ(contact.getResponsibleRank(), newcontact.getResponsibleRank());
  EXPECT_EQ(contact.getDampingFactor(), newcontact.getDampingFactor());
  EXPECT_EQ(contact.isNew(), newcontact.isNew());
}

TEST(WallContactSerialization, serializationAndDeserialization2d)
{
  constexpr unsigned                                        D = 2;
  typedef WallContactArbitraryFromOverlapVolume<T, D, true> CONTACTTYPE;

  std::size_t id(4);
  unsigned    wallID = 3;
  CONTACTTYPE contact(id, wallID);
  PhysR<T, D> position;
  position = PhysR<T, D>(5, 4);
  contact.setParticlePosition(position);
  PhysR<T, D> min(-5.5);
  contact.updateMinMax(min);
  PhysR<T, D> max(9.7);
  contact.updateMinMax(max);
  contact.setParticlePositionUpdated(false);
  contact.setResponsibleRank(36871);
  contact.isNew(true);

  std::size_t serialSize = contact.getSerialSize();
  EXPECT_EQ(74, serialSize);

  std::uint8_t* buffer = new std::uint8_t[serialSize] {};
  EXPECT_EQ(serialSize, contact.serialize(buffer));

  CONTACTTYPE newcontact;
  EXPECT_EQ(serialSize, newcontact.deserialize(buffer));

  EXPECT_TRUE(!newcontact.isParticlePositionUpdated());
  for (unsigned iD = 0; iD < D; ++iD) {
    EXPECT_EQ(position[iD], newcontact.getParticlePosition()[iD]);
    EXPECT_EQ(min[iD], newcontact.getMin()[iD]);
    EXPECT_EQ(max[iD], newcontact.getMax()[iD]);
  }
  EXPECT_EQ(id, newcontact.getParticleID());
  EXPECT_EQ(wallID, newcontact.getWallID());
  EXPECT_EQ(contact.getResponsibleRank(), newcontact.getResponsibleRank());
  EXPECT_EQ(contact.getDampingFactor(), newcontact.getDampingFactor());
  EXPECT_EQ(contact.isNew(), newcontact.isNew());
}
