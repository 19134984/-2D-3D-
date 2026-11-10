#!/bin/bash
#SBATCH --partition=dev-g
#SBATCH --nodes=1 
#SBATCH --tasks-per-node=7
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00 
#SBATCH --account=project_XXXX
#SBATCH --output=compile_%j.out
#SBATCH --error=compile_%j.err

module --force purge > /dev/null 2>&1
module load LUMI/25.03 partition/G PrgEnv-amd rocm cray-python

# IMPORTANT: set your python environment here
source $HOME/roofline_pyenv/bin/activate

./compile_multiple_dynamics.sh --precision float ../target/wind.dynamics

