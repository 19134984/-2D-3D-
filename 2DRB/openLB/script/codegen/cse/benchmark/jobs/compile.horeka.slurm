#!/bin/bash
#SBATCH --partition=cpuonly
#SBATCH --nodes=1 
#SBATCH --tasks-per-node=6
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00 
#SBATCH --output=compile_%j.out
#SBATCH --error=compile_%j.err

module purge
module use /software/easybuild/modules/all/
module load OpenMPI/5.0.7-GCC-14.2.0-CUDA-12.8.0

source $HOME/roofline_pyenv/bin/activate

./compile_multiple_dynamics.sh --precision float ../target/wind.dynamics

