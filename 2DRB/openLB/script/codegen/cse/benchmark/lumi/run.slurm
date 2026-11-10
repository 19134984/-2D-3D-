#!/bin/bash
#SBATCH --partition=dev-g 
#SBATCH --nodes=1 
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=7 
#SBATCH --gpus-per-node=1
#SBATCH --time=00:10:00 
#SBATCH --account=project_XXXX
#SBATCH --output=benchmark_%j.out
#SBATCH --error=benchmark_%j.err

module --force purge > /dev/null 2>&1
module load LUMI/25.03 partition/G PrgEnv-amd rocm cray-python

# IMPORTANT: set your python environment here
source ../.roofline_pyenv/bin/activate

./run_multiple_dynamics.sh --size 10 ../target/roofline.dynamics

