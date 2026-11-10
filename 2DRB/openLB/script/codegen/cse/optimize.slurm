#!/bin/bash
#SBATCH --partition=dev-g 
#SBATCH --nodes=1 
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=7 
#SBATCH --time=00:30:00 
#SBATCH --account=project_XXXX
#SBATCH --output=optimize_%j.out
#SBATCH --error=optimize_%j.err
#SBATCH --mem=0
module --force purge > /dev/null 2>&1
module load LUMI/25.03 partition/G PrgEnv-amd rocm
module load cray-python

source $HOME/roofline_pyenv/bin/activate

#./optimize_dynamics.sh --install "dynamics::Tuple<float, olb::descriptors::D3Q19<olb::descriptors::POROSITY, olb::descriptors::VELOCITY, olb::fields::fsi::ELEMENT_TAG, olb::fields::fsi::ELEMENT_TORQUE, olb::collision::LES::SMAGORINSKY>, olb::momenta::Tuple<olb::momenta::BulkDensity, olb::momenta::BulkMomentum, olb::momenta::BulkStress, olb::momenta::DefineToNEq>, olb::equilibria::ThirdOrder, olb::collision::ParameterFromCell<olb::collision::LES::SMAGORINSKY, olb::collision::SmagorinskyEffectiveOmega<olb::collision::ThirdOrderRLB>>>"

./optimize_dynamics.sh --install "dynamics::Tuple<float, olb::descriptors::D3Q19<olb::descriptors::POROSITY, olb::descriptors::VELOCITY, olb::fields::fsi::ELEMENT_TAG, olb::fields::fsi::ELEMENT_TORQUE, olb::collision::LES::SMAGORINSKY>, olb::momenta::Tuple<olb::momenta::BulkDensity, olb::momenta::MovingPorousMomentumCombination<olb::momenta::BulkMomentum>, olb::momenta::BulkStress, olb::momenta::DefineToNEq>, olb::equilibria::ThirdOrder, olb::collision::ParameterFromCell<olb::collision::HYBRID, olb::collision::LocalSmagorinskyEffectiveOmega<olb::collision::HRR>>>"

#./benchmark_dynamics.sh "dynamics::Tuple<double, olb::descriptors::D2Q9<>, olb::momenta::Tuple<olb::momenta::BulkDensity, olb::momenta::BulkMomentum, olb::momenta::BulkStress, olb::momenta::DefineToNEq>, olb::equilibria::SecondOrder, olb::collision::BGK, olb::dynamics::DefaultCombination>"
