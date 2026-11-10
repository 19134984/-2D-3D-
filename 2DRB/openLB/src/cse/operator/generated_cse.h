/*  ========================================================
 *  ==  WARNING: This is an automatically generated file, ==
 *  ==                  do not modify.                    ==
 *  ========================================================
 */

#ifndef DISABLE_CSE

//OuterVelocityCornerProcessor2D<Expr,descriptors::D2Q9<>,-1,-1>;descriptors::D2Q9<>
#include "e6228df2db446705f50cfff1fbceed5b29dc95632fb2c44fbaebb1d28b4eb24f.cse.h"

//OuterVelocityCornerProcessor2D<Expr,descriptors::D2Q9<>,-1,1>;descriptors::D2Q9<>
#include "7f1c38ec23060725a2e83e1352d92fa69aee53728a1efe50eaae064e7e19d059.cse.h"

//OuterVelocityCornerProcessor2D<Expr,descriptors::D2Q9<>,1,-1>;descriptors::D2Q9<>
#include "0b56821c54dfa17605cfe8bc40263dbcc04315079206bc85c268e0f0de107b45.cse.h"

//OuterVelocityCornerProcessor2D<Expr,descriptors::D2Q9<>,1,1>;descriptors::D2Q9<>
#include "1b923fea59b4ad57e90209b7c338c0ebcc8e0e55b2276358473565572ce89294.cse.h"

//PlaneFdBoundaryProcessor3D<Expr,descriptors::D3Q19<>,0,-1>;descriptors::D3Q19<>
#include "94e7d721b27f6f13b5b94ba700c55072e60252d715a596cb4f07cedd5aaa9011.cse.h"

//StraightConvectionBoundaryProcessor3D<Expr,descriptors::D3Q19<>,0,-1>;descriptors::D3Q19<>
#include "20d20494d9cf1941db1bae2473a60414ccb02b4518ad90812f541238500d9317.cse.h"

//StraightFdBoundaryProcessor2D<Expr,descriptors::D2Q9<>,0,-1>;descriptors::D2Q9<>
#include "c9d13dd4a4c09f42524778dd16abcbfcf972bd5ca3caafed51894b32c73ab522.cse.h"

//operators::TrackAverageFunctorO<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, opti::AVERAGE_J>;descriptors::D3Q19<>
#include "4d9fbc50b25871271408cc71e27d4181b96ffa53a34b830590d91398d036e3b1.cse.h"

//operators::TrackAverageFunctorO<functors::AddF<functors::DissipationF,functors::PorousDissipationF>,opti::AVERAGE_J>;descriptors::D2Q9<>
#include "f8f2ba62fc1f06142b5bdfd1ecbec28aa7200d23286211964c5d8b5aec16691f.cse.h"

//operators::TrackAverageFunctorO<functors::DerivativeF<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, descriptors::POPULATION, dynamics::Tuple<Expr, descriptors::D2Q9<>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination> >, opti::AVERAGE_DJDF>;descriptors::D2Q9<>
#include "8649ede6cd66d6b9a47ac8f526a12618edfe5333fa7cea8f6bab6df853756cd8.cse.h"

//operators::TrackAverageFunctorO<functors::DerivativeF<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, descriptors::POPULATION, dynamics::Tuple<Expr, descriptors::D3Q19<>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination> >, opti::AVERAGE_DJDF>;descriptors::D3Q19<>
#include "d70b17d794143381026ab1de30c5a54fc43f67f7cb13c410991ab1f6a98cfcb6.cse.h"

//operators::TrackAverageFunctorO<functors::DerivativeF<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, descriptors::POPULATION, dynamics::Tuple<Expr, descriptors::D3Q19<>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::ThirdOrder, collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB>, dynamics::DefaultCombination> >, opti::AVERAGE_DJDF>;descriptors::D3Q19<>
#include "abad0cf1e082712cc935c905fdb25dc1af0c3c13d1df42877d5b3ee02fccebfb.cse.h"

//operators::TrackAverageFunctorO<functors::DerivativeF<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, descriptors::POROSITY, dynamics::Tuple<Expr, descriptors::D2Q9<>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination> >, opti::AVERAGE_DJDALPHA<descriptors::POROSITY> >;descriptors::D2Q9<>
#include "947da5e42fc6b67b446a47d4be28cafb4d78afe17592412c8fe35256d36766f7.cse.h"

//operators::TrackAverageFunctorO<functors::DerivativeF<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, descriptors::POROSITY, dynamics::Tuple<Expr, descriptors::D3Q19<>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination> >, opti::AVERAGE_DJDALPHA<descriptors::POROSITY> >;descriptors::D3Q19<>
#include "aa508c3f95646648664473ac99b32634ac428c48ce297072e7d6db247bc19b4f.cse.h"

//operators::TrackAverageFunctorO<functors::DerivativeF<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, descriptors::POROSITY, dynamics::Tuple<Expr, descriptors::D3Q19<>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, dynamics::DefaultCombination> >, opti::AVERAGE_DJDALPHA<descriptors::POROSITY> >;descriptors::D3Q19<>
#include "7d5d640784f5fcac6b31580459a7df9139bdbe73fb9b2803428055bb59023388.cse.h"

//operators::TrackAverageFunctorO<functors::DerivativeF<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, descriptors::POROSITY, dynamics::Tuple<Expr, descriptors::D3Q19<>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::ThirdOrder, collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB>, dynamics::DefaultCombination> >, opti::AVERAGE_DJDALPHA<descriptors::POROSITY> >;descriptors::D3Q19<>
#include "8937eca9f25063ba4833fd6d96c62ff2a8c8ca5de12ffd409e5f822c461cc8ee.cse.h"

//operators::TrackAverageFunctorO<functors::DerivativeF<functors::DissipationF, descriptors::POROSITY, dynamics::Tuple<Expr, descriptors::D3Q19<>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, dynamics::DefaultCombination> >, opti::AVERAGE_DJDALPHA<descriptors::POROSITY> >;descriptors::D3Q19<>
#include "71010b6d7646be762858b74e8972a1cef3c5d2c8156ebf701e762e1309fe8640.cse.h"

//operators::TrackAverageFunctorO<functors::DerivativeF<functors::FarFieldDragF<0>, descriptors::POPULATION, dynamics::Tuple<Expr, descriptors::D3Q19<>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, dynamics::DefaultCombination> >, opti::AVERAGE_DJDF>;descriptors::D3Q19<>
#include "aa539399513d9863f2c60870b0addd3315e3a1140b3a4348996ab4e55097b458.cse.h"

//operators::TrackAverageFunctorO<functors::DerivativeF<functors::FarFieldDragF<0>, descriptors::POPULATION, dynamics::Tuple<Expr, descriptors::D3Q19<>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::ThirdOrder, collision::FrozenOmegaFromCell<collision::ThirdOrderRLB>, dynamics::DefaultCombination> >, opti::AVERAGE_DJDF>;descriptors::D3Q19<>
#include "35526bd6325cf51c4038eb01af9031c0f86bfeb91187f982b10a03c13cb9d4ed.cse.h"

//operators::TrackAverageFunctorO<functors::DerivativeF<functors::FarFieldDragF<0>, descriptors::POPULATION, dynamics::Tuple<Expr, descriptors::D3Q19<>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::ThirdOrder, collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB>, dynamics::DefaultCombination> >, opti::AVERAGE_DJDF>;descriptors::D3Q19<>
#include "901cdba618f30ab93391ab14ab9483cc2bc1223ecb045e161a626d4e1bd922e2.cse.h"

//operators::TrackAverageFunctorO<functors::DerivativeF<functors::FarFieldDragF<0>, descriptors::POROSITY, dynamics::Tuple<Expr, descriptors::D3Q19<>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, dynamics::DefaultCombination> >, opti::AVERAGE_DJDALPHA<descriptors::POROSITY> >;descriptors::D3Q19<>
#include "d2719d07d011db85576d723b5ec09b34b05c1ab447237a5cded4bb762e0faba8.cse.h"

//operators::TrackAverageFunctorO<functors::DerivativeF<functors::FarFieldDragF<0>, descriptors::POROSITY, dynamics::Tuple<Expr, descriptors::D3Q19<>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::ThirdOrder, collision::FrozenOmegaFromCell<collision::ThirdOrderRLB>, dynamics::DefaultCombination> >, opti::AVERAGE_DJDALPHA<descriptors::POROSITY> >;descriptors::D3Q19<>
#include "8ae9902692ea5245bb9314973f7305775178e0e9721412f962007868199e5119.cse.h"

//operators::TrackAverageFunctorO<functors::DerivativeF<functors::FarFieldDragF<0>, descriptors::POROSITY, dynamics::Tuple<Expr, descriptors::D3Q19<>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::ThirdOrder, collision::SmagorinskyEffectiveOmega<collision::ThirdOrderRLB>, dynamics::DefaultCombination> >, opti::AVERAGE_DJDALPHA<descriptors::POROSITY> >;descriptors::D3Q19<>
#include "108cb0a7d6e75ec5c36b52ec6b970c6bbd81ff8e1d76959e4794f1fc835ece77.cse.h"

//operators::TrackAverageFunctorO<functors::DissipationF, opti::AVERAGE_J>;descriptors::D3Q19<>
#include "bb7eead6857680addc048f9499f94bd1cf2625e73dd30ea7bf64c2f028bea1ff.cse.h"

//operators::TrackAverageFunctorO<functors::FarFieldDragF<0>, opti::AVERAGE_J>;descriptors::D3Q19<>
#include "28a700ca418027546e2d3d147bb6d1c533937458b3c538a03c5d76cfeee6729a.cse.h"

//operators::TrackAverageFunctorO<functors::ReverseDerivativeF<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, descriptors::POPULATION, dynamics::Tuple<Expr, descriptors::D2Q9<>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination> >, opti::AVERAGE_DJDF>;descriptors::D2Q9<>
#include "5a3f155467dd0c8b8d6b0fd1b02bdae9ab67d8084bf789b3b5e5582e1e20cb4f.cse.h"

//operators::TrackAverageFunctorO<functors::ReverseDerivativeF<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, descriptors::POPULATION, dynamics::Tuple<Expr, descriptors::D3Q19<>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::BGK, dynamics::DefaultCombination> >, opti::AVERAGE_DJDF>;descriptors::D3Q19<>
#include "b8e96a6d28e1879205aa9e54fd6d333225ecc36251adf7dea184cebe96d4f0f0.cse.h"

//operators::TrackAverageFunctorO<functors::ReverseDerivativeF<functors::AddF<functors::DissipationF, functors::PorousDissipationF>, descriptors::POPULATION, dynamics::Tuple<Expr, descriptors::D3Q19<>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, dynamics::DefaultCombination> >, opti::AVERAGE_DJDF>;descriptors::D3Q19<>
#include "ff0729e6ed399c07128d5600ed3da7851ab805ec5830e3e125e1fb6a651fa96c.cse.h"

//operators::TrackAverageFunctorO<functors::ReverseDerivativeF<functors::DissipationF, descriptors::POPULATION, dynamics::Tuple<Expr, descriptors::D3Q19<>, momenta::Porous<momenta::Tuple<momenta::BulkDensity, momenta::BulkMomentum, momenta::BulkStress, momenta::DefineToNEq> >, equilibria::SecondOrder, collision::SmagorinskyEffectiveOmega<collision::BGK>, dynamics::DefaultCombination> >, opti::AVERAGE_DJDF>;descriptors::D3Q19<>
#include "ca7e7fc972c1250244a5b8825d8b1cdbde4b70725585a8a0cbccf743a73b2714.cse.h"

#endif
