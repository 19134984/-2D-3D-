/*  This file is part of the OpenLB library
*
*  Copyright (C) 2022 Nicolas Hafen, Martin Sadric
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
#include <type_traits>
#include <gtest/gtest.h>
#include "helper/unitTestHelper.h"
#include "helper/framework.h"

#include <olb.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::particles;
using namespace olb::particles::dynamics;
using namespace olb::particles::access;

typedef double T;


::testing::Environment* const env = ::testing::AddGlobalTestEnvironment(new OlbInitTestEnvironment);

//#define VERBOSE_TESTS

//test structs for meta test
struct A {
  static constexpr int id = 1;
  static constexpr bool select=true;
};

struct B {
  static constexpr int id = 2;
  static constexpr bool select=true;
};

struct C {
  static constexpr int id = 3;
  static constexpr bool select=false;
};

struct D {
  static constexpr int id = 4;
  static constexpr bool select=true;
};

struct E {
  static constexpr int id = 5;
  static constexpr bool select=false;
};

struct F {
  static constexpr int id = 6;
  static constexpr bool select=false;
};

template <typename TASK>
using is_selected = std::integral_constant<bool, TASK::select>;

template <typename TASK>
using is_not_selected = std::integral_constant<bool, !TASK::select>;

template <typename... TYPES>
using getFilteredTypes = typename meta::filter<is_selected,TYPES...>::type;

template <typename... TYPES>
using getNonSelectedTypes = typename meta::filter<is_not_selected,TYPES...>::type;

template <typename... TYPES>
using getAllTypes = typename meta::list<TYPES...>;

template <typename TASKLIST>
void printTasks(){
  TASKLIST::for_each([&](auto task){
      using TASK=typename decltype(task)::type;
      int id = TASK::id;
#ifdef VERBOSE_TESTS
      std::cout << "Task: " << id << std::endl;
#endif
  });
}

template <typename ...TASKLIST>
void execute(){
  using tasksFiltered = getFilteredTypes<TASKLIST...>;
  using tasksRev = typename tasksFiltered::template decompose_into<meta::reverse_t>;
  printTasks<tasksFiltered>();
#ifdef VERBOSE_TESTS
  std::cout << "-Reversed order:" << std::endl;
#endif
  printTasks<tasksRev>();
}

TEST(TestParticleManager, MetaTest)
{
  // TODO: Add expect?
  execute<A,B,C,D>();
}


TEST(TestParticleManager, MetaListSelectedIndices)
{
  //Define task list
  using taskList = meta::list<A,B,C,D,A,B,C,D>;
  //Define index list
  auto indexList = std::index_sequence<1,3,4,5,5>();
#ifdef VERBOSE_TESTS
  //Print id for all tasks specified
  meta::list_for_each_index<taskList>([](auto task){
    using TASK = typename decltype(task)::type;
    std::cout << TASK::id << std::endl;
  },indexList);
#endif

  EXPECT_EQ(taskList::size, 8);
  EXPECT_EQ(indexList.size(), 5);
}

TEST(TestParticleManager, MetaListFilterIndices)
{
  //Define task list
  using taskList = meta::list<A,B,C,D,A,B,C,D>;
  //Define index list
  auto indexList = std::index_sequence<0,1,2,3,4,5,5,5>();

  //Filter index list by evaluating respective tasks
  auto indexListFiltered =
    meta::filter_index_sequence<taskList>([](auto task){
      using TASK = typename decltype(task)::type;
      return TASK::select;
    },indexList);

#ifdef VERBOSE_TESTS
  //Print id for all tasks specified
  meta::list_for_each_index<taskList>([](auto task){
    using TASK = typename decltype(task)::type;
    std::cout << TASK::id << std::endl;
  },indexListFiltered);
#endif

  EXPECT_EQ(taskList::size, 8);
  EXPECT_EQ(indexList.size(), 8);
  EXPECT_EQ(indexListFiltered.size(), 7);
}

TEST(TestParticleManager, MetaListFilterIndicesByCondition)
{
  //Define task list
  using taskList = meta::list<A,B,C,D,A,B,C,D>;
  //Define index list
  auto indexList = std::index_sequence<0,1,2,3,4,5,5,5>();

  //Filter index list by evaluating COND
  auto indexListFiltered =
    meta::filter_index_sequence<taskList,is_selected>();

#ifdef VERBOSE_TESTS
  //Print id for all tasks specified
  meta::list_for_each_index<taskList>([](auto task){
    using TASK = typename decltype(task)::type;
    std::cout << TASK::id << std::endl;
  },indexListFiltered);
#endif

  EXPECT_EQ(taskList::size, 8);
  EXPECT_EQ(indexList.size(), 8);
  EXPECT_EQ(indexListFiltered.size(), 6);
}

TEST(TestParticleManager, shiftIndexSequence)
{
  //Define task list
  using taskList = meta::list<A,B,C,D,A,B,C,D>;
  //Define index list
  auto indexList = std::index_sequence<0,1,2,3,4,5,5,5>();

  //Create index list from range
  auto indexListFromRange =
    meta::make_index_sequence_in_range<2,6>();

#ifdef VERBOSE_TESTS
  //Print id for all tasks specified
  meta::list_for_each_index<taskList>([](auto task){
    using TASK = typename decltype(task)::type;
    std::cout << TASK::id << std::endl;
  },indexListFromRange);
#endif

  EXPECT_EQ(taskList::size, 8);
  EXPECT_EQ(indexList.size(), 8);
  EXPECT_EQ(indexListFromRange.size(), 4);
}



TEST(TestParticleManager, subIndexSequenceVerbose)
{
  //Define task list and intended number of tasks and subtasks
  // - defined by condition is_not_selected
  using taskList = meta::list<C,A,B,C,D,A,C,D,A,B,C,D,A,C>;
  int iTasks = 0;
  int iSubTasks = 0;

#ifdef VERBOSE_TESTS
  //Print id for all tasks
  std::cout << "All Tasks" << std::endl;
  meta::list_for_each_index<taskList>([](auto task, std::size_t idx){
    using TASK = typename decltype(task)::type;
    std::cout << idx
              << " name=" << typeid(TASK).name()
              << " id=" << TASK::id
              << " select=" << TASK::select
              << std::endl;
  },std::make_index_sequence<taskList::size>());
  std::cout << std::endl;
#endif

  //Filter index list by evaluating COND
  auto indexListMainTasks =
    meta::filter_index_sequence<taskList,is_not_selected>();

#ifdef VERBOSE_TESTS
  //Print main and sub index_sequences
  std::cout << "IndexListFiltered: " << indexListMainTasks << std::endl;
  meta::index_sequence_for_subsequence<taskList>([](auto task){},
  [](auto seq){
      std::cout << " SubIndexList: " << seq << std::endl;
  }, indexListMainTasks);
  std::cout << std::endl;
#endif

  //Call tasks
  meta::index_sequence_for_subsequence<taskList>(
  [&](auto task){
#ifdef VERBOSE_TESTS
    using TASK = typename decltype(task)::type;
    std::cout << ""
              << "name=" << typeid(TASK).name()
              << " select=" << TASK::select
              << std::endl;
#endif
    iTasks++;
  },
  [&](auto seq){
    //Iterate over subtasks
    meta::list_for_each_index<taskList>([&](auto task){
#ifdef VERBOSE_TESTS
      using TASK = typename decltype(task)::type;
      std::cout << "-"
                << " name=" << typeid(TASK).name()
                << " select=" << TASK::select
                << std::endl;
#endif
      iSubTasks++;
    },seq);
  }, indexListMainTasks);

  EXPECT_EQ(iTasks, 5);
  EXPECT_EQ(iSubTasks,9);
}


TEST(TestParticleManager, subIndexSequenceRepresentative)
{
  //Define task list and intended number of tasks and subtasks
  // - defined by condition is_not_selected
  using taskList = meta::list<C,E,A,B,C,D>;
  int iTasks = 0;
  int iSubTasks = 0;

  std::vector<int> sampledSequence;

  //Filter index list by evaluating COND
  auto indexListMainTasks =
    meta::filter_index_sequence<taskList,is_not_selected>();

  //Call tasks
  meta::index_sequence_for_subsequence<taskList>(
  [&](auto task){
    using TASK = typename decltype(task)::type;
    sampledSequence.push_back(TASK::id);
    iTasks++;
  },
  [&](auto seq){
    //Iterate over subtasks
    meta::list_for_each_index<taskList>([&](auto task){
      using TASK = typename decltype(task)::type;
      sampledSequence.push_back(TASK::id);
      iSubTasks++;
    },seq);
  }, indexListMainTasks);

  Vector<int,6> seqVec( sampledSequence );
  Vector<int,6> seqVecRef(3,5,1,2,3,4);

  EXPECT_EQ(seqVec, seqVecRef);
  EXPECT_EQ(iTasks, 3);
  EXPECT_EQ(iSubTasks,3);
}


TEST(TestParticleManager, subIndexSequenceLong)
{
  //Define task list and intended number of tasks and subtasks
  // - defined by condition is_not_selected
  using taskList = meta::list<C,E,A,B,C,D,E,A,B,C,D,E,A,A,A,B,C,D>;
  int iTasks = 0;
  int iSubTasks = 0;

  std::vector<int> sampledSequence;

  //Filter index list by evaluating COND
  auto indexListMainTasks =
    meta::filter_index_sequence<taskList,is_not_selected>();

  //Call tasks
  meta::index_sequence_for_subsequence<taskList>(
  [&](auto task){
    using TASK = typename decltype(task)::type;
    sampledSequence.push_back(TASK::id);
    iTasks++;
  },
  [&](auto seq){
    //Iterate over subtasks
    meta::list_for_each_index<taskList>([&](auto task){
      using TASK = typename decltype(task)::type;
      sampledSequence.push_back(TASK::id);
      iSubTasks++;
    },seq);
  }, indexListMainTasks);

  Vector<int,18> seqVec( sampledSequence );
  Vector<int,18> seqVecRef(3,5,1,2,3,4,5,1,2,3,4,5,1,1,1,2,3,4);

  EXPECT_EQ(seqVec, seqVecRef);
  EXPECT_EQ(iTasks, 7);
  EXPECT_EQ(iSubTasks,11);
}

TEST(TestParticleManager, TestFluidVelocityWithExternalAcceleration2D)
{
  //Settings
  typedef D2Q9<> DESCRIPTOR;
  constexpr unsigned D = DESCRIPTOR::d;
  typedef ResolvedCircle2D PARTICLETYPE;
  PhysR<T,D> origin (0.0);
  const T baseSize = 1.0;
  const T delX = 0.1;
  const T particleDensity = 0.57;
  Vector<T,D> externalAcceleration(-3.25);

  //Set up environment
  Vector<T,D> extend( baseSize );
  IndicatorCuboid<T,D> cuboid( extend, origin );
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 4;
#endif
  CuboidDecomposition<T,D> cuboidDecomposition( cuboid, delX, noOfCuboids );
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,D> sGeometry( cuboidDecomposition, loadBalancer, 2 );
  UnitConverter<T,DESCRIPTOR> converter( delX, delX, baseSize, 1., 1., 1.);
  SuperLattice<T, DESCRIPTOR> sLattice( converter, sGeometry );


  //Create particle system
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  creators::addResolvedCircle2D( particleSystem, Vector<T,D>(0.), 1., .5*delX, particleDensity );
  creators::addResolvedCircle2D( particleSystem, Vector<T,D>(0.), 1., .5*delX, particleDensity/2 );
  creators::addResolvedCircle2D( particleSystem, Vector<T,D>(0.), 1., .5*delX, particleDensity*2 );
  // Create and assign resolved particle dynamics
  particleSystem.defineDynamics<
    VerletParticleDynamics<T,PARTICLETYPE>>();

  //Create particle manager
  ParticleManager<T,DESCRIPTOR,PARTICLETYPE> particleManager(
      particleSystem, sGeometry, sLattice, converter, externalAcceleration );

  //Execute manager
  particleManager.execute<
    couple_lattice_to_particles_single_cuboid<T,DESCRIPTOR,PARTICLETYPE>,
    apply_gravity<T,PARTICLETYPE>,
    process_dynamics_single_cuboid<T,PARTICLETYPE>,
    couple_particles_to_lattice_single_cuboid<T,DESCRIPTOR,PARTICLETYPE>
  >();

  for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
    auto particle = particleSystem.get(iP);
    auto force = getForce<T,PARTICLETYPE>( particle );
    auto position = getPosition<T,PARTICLETYPE>( particle );
    auto velocity = getVelocity<T,PARTICLETYPE>( particle );
    auto avgAcc = particle.template getField<MOBILITY,ACCELERATION_STRD>();
    for(unsigned iD=0; iD<D; ++iD) {
      const T mass = particle.template getField<PHYSPROPERTIES,MASS>();
      EXPECT_EQ(force[iD], externalAcceleration[iD] * mass);
      EXPECT_NEAR(position[iD], 0., 1e-16);
      EXPECT_NEAR(velocity[iD], 0.05*externalAcceleration[iD], 1e-16);
      EXPECT_NEAR(avgAcc[iD], externalAcceleration[iD], 1e-16);
    }
  }
}

TEST(TestParticleManager, TestFluidVelocityWithExternalAcceleration3D)
{
  //Settings
  typedef D3Q19<> DESCRIPTOR;
  constexpr unsigned D = DESCRIPTOR::d;
  typedef ResolvedSphere3D PARTICLETYPE;
  PhysR<T,D> origin (0.0);
  const T baseSize = 1.0;
  const T delX = 0.1;
  const T particleDensity = 5700;
  Vector<T,D> externalAcceleration(-5.7);

  //Set up environment
  Vector<T,D> extend( baseSize );
  IndicatorCuboid<T,D> cuboid( extend, origin );
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 4;
#endif
  CuboidDecomposition<T,D> cuboidDecomposition( cuboid, delX, noOfCuboids );
  HeuristicLoadBalancer<T> loadBalancer( cuboidDecomposition );
  SuperGeometry<T,D> sGeometry( cuboidDecomposition, loadBalancer, 2 );
  UnitConverter<T,DESCRIPTOR> converter( delX, delX, baseSize, 1., 1., 1.);
  SuperLattice<T, DESCRIPTOR> sLattice( converter, sGeometry );


  //Create particle system
  ParticleSystem<T,PARTICLETYPE> particleSystem;
  creators::addResolvedSphere3D( particleSystem, Vector<T,D>(0.), 1., .5*delX, particleDensity );
  creators::addResolvedSphere3D( particleSystem, Vector<T,D>(0.), 1., .5*delX, particleDensity/2 );
  creators::addResolvedSphere3D( particleSystem, Vector<T,D>(0.), 1., .5*delX, particleDensity*2 );
  // Create and assign resolved particle dynamics
  particleSystem.defineDynamics<
    VerletParticleDynamics<T,PARTICLETYPE>>();

  //Create particle manager
  ParticleManager<T,DESCRIPTOR,PARTICLETYPE> particleManager(
      particleSystem, sGeometry, sLattice, converter, externalAcceleration );

  //Execute manager
  particleManager.execute<
    couple_lattice_to_particles_single_cuboid<T,DESCRIPTOR,PARTICLETYPE>,
    apply_gravity<T,PARTICLETYPE>,
    process_dynamics_single_cuboid<T,PARTICLETYPE>,
    couple_particles_to_lattice_single_cuboid<T,DESCRIPTOR,PARTICLETYPE>
  >();

  for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
    auto particle = particleSystem.get(iP);
    auto force = getForce<T,PARTICLETYPE>( particle );
    auto position = getPosition<T,PARTICLETYPE>( particle );
    auto velocity = getVelocity<T,PARTICLETYPE>( particle );
    auto avgAcc = particle.template getField<MOBILITY,ACCELERATION_STRD>();
    for(unsigned iD=0; iD<D; ++iD) {
      const T mass = particle.template getField<PHYSPROPERTIES,MASS>();
      EXPECT_EQ(force[iD], externalAcceleration[iD] * mass);
      EXPECT_NEAR(position[iD], 0., 1e-16);
      EXPECT_NEAR(velocity[iD], 0.05*externalAcceleration[iD], 1e-16);
      EXPECT_NEAR(avgAcc[iD], externalAcceleration[iD], 1e-16);
    }
  }
}


